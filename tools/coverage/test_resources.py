#!/usr/bin/env python3
"""Replay the reference's own tests on both sides: criterion (1) of a bit-identical tool.

The inputs that apply to a tool are the ones Picard's own tests (`src/test/java/**/*Test.java`)
hand it, with the arguments they hand them with, out of `testdata/`. This reads those invocations
out of the test sources and runs each one against the pinned oracle and against the port, reading
and comparing what they write with `run_array.py`'s own readers and comparison.

The invocations are read STATICALLY: the test classes are parsed (tree-sitter) and their test
methods partially evaluated -- constants, String arrays and lists, helper methods, `@DataProvider`
rows, loops and branches whose condition is known, `runPicardCommandLine`, `instanceMain` and a
configured tool's `doWork()`. Running the Java tests instead would need their whole classpath and
would run the reference's assertions rather than a comparison. A test whose arguments cannot be
evaluated, or whose input is a file the test itself writes in code this does not run
(`SAMRecordSetBuilder`, `FileUtils.write`), is listed under `unextracted` with the reason, so what
the harness does not replay is counted rather than missed.

    # once per pinned release, on a host with the sources; the result is committed
    python3 tools/coverage/test_resources.py extract --picard ../picard-3.5.0 \
        --constants ../htsjdk-5.0.0/src/main/java
    # what CI runs (stdlib only), one shard of the replay
    python3 tools/coverage/test_resources.py run --picard picard-src --port-dir target/release --shard 0/8
    python3 tools/coverage/test_resources.py merge shard-*.json
    python3 tools/coverage/test_resources.py report   # the column in docs/VALIDATION.md

`extract` needs `pip install tree-sitter==0.26.1 tree-sitter-java==0.23.5`; nothing else does.
"""

import argparse
import hashlib
import concurrent.futures
import json
import os
import re
import shlex
import subprocess
import sys
import tempfile
import threading
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import run_array  # noqa: E402

REPO = HERE.parents[1]
DATA = HERE / "test_resources"
INVOCATIONS = DATA / "picard.json"
EXCLUSIONS = DATA / "exclusions.json"
MEASURED = DATA / "measured.json"
MANIFEST = HERE.parent / "conformance" / "manifest.json"
VALIDATION = REPO / "docs" / "VALIDATION.md"

# Where the reference's checkout is mounted on both sides; Picard's tests run with it as `user.dir`
# and name their inputs relative to it (`testdata/picard/...`).
ROOT = "/work/fixtures"
OUT = "/work/out"
# TestNG groups whose tests read from a bucket or need credentials, which neither side has.
CLOUD_GROUPS = {"cloud", "bucket", "gcs", "spark_cloud", "cloud_requester_pays"}
RELATIVE_RESOURCE = re.compile(r"(?<![\w/.])testdata/")


# ---------------------------------------------------------------------------------------------
# Values of the partial evaluator.


class Unresolved(Exception):
    """An expression whose value the sources alone do not determine."""


class Return(Exception):
    def __init__(self, value):
        self.value = value


class Thrown(Exception):
    """A `throw` reached: the test method ends there."""


class Poison:
    """A variable whose value could not be computed; reading it is what fails, not binding it."""

    def __init__(self, reason):
        self.reason = reason


class ClassRef:
    def __init__(self, name, info=None):
        self.name, self.info = name, info


class EnumValue:
    def __init__(self, name, info=None):
        self.name, self.info = name, info


class Tool:
    """A tool instance a test builds (`new SortSam()`) and configures field by field before it
    calls `doWork()` or `instanceMain(args)` on it. The fields it set are its command line."""

    def __init__(self, name):
        self.name, self.fields = name, {}

    def tokens(self):
        out = []
        for name, value in self.fields.items():
            for item in value if isinstance(value, list) else [value]:
                out.append(f"{name}={jstr(item)}")
        return out


class Optional:
    """`java.util.Optional`."""

    def __init__(self, value=None):
        self.value = value


class Lambda:
    def __init__(self, node, env, this):
        self.node, self.env, self.this = node, env, this


def jdouble(value):
    """`Double.toString`, for the doubles a test concatenates into an argument."""
    if value != value or value in (float("inf"), float("-inf")):
        return {True: "NaN"}.get(value != value, "Infinity" if value > 0 else "-Infinity")
    if value == 0:
        return "0.0"
    if 1e-3 <= abs(value) < 1e7:
        text = repr(value)
        return text if "." in text or "e" in text else text + ".0"
    mantissa, exponent = f"{value:.17e}".split("e")
    shortest = repr(value)
    digits = shortest.replace("-", "").split("e")[0].replace(".", "").lstrip("0").rstrip("0") or "0"
    sign = "-" if value < 0 else ""
    return f"{sign}{digits[0]}.{digits[1:] or '0'}E{int(exponent)}"


def jstr(value):
    """`String.valueOf`."""
    if isinstance(value, Poison):
        raise Unresolved(value.reason)
    if value is None:
        return "null"
    if value is True:
        return "true"
    if value is False:
        return "false"
    if isinstance(value, float):
        return jdouble(value)
    if isinstance(value, (int, str)):
        return str(value)
    if isinstance(value, EnumValue):
        return value.name
    if isinstance(value, ClassRef):
        # An unresolved capitalised name falls back to a class; its string form is never an argument.
        raise Unresolved(f"no string form for the class {value.name}")
    if isinstance(value, list):
        return "[" + ", ".join(jstr(v) for v in value) + "]"
    raise Unresolved(f"no string form for {type(value).__name__}")


def absolute(value):
    """`File.getAbsolutePath()`: Picard's tests run from the checkout's root."""
    text = jstr(value)
    if text.startswith(("/", "gs://", "gendb", "hdfs:", "file:", "http")) or not text:
        return text
    return f"{ROOT}/{text}"


def normalized(path):
    """What `new File(path)` and `Paths.get(path)` keep: no repeated or trailing separator."""
    if "://" in path:
        return path
    path = re.sub(r"/+", "/", path)
    return path.rstrip("/") if len(path) > 1 else path


def join_path(parent, child):
    parent, child = jstr(parent), jstr(child)
    if child.startswith("/") or not parent:
        return normalized(child)
    return normalized(parent + "/" + child)


def jformat(fmt, args):
    """`String.format` for the conversions tests use."""
    args = list(args)
    position = 0

    def one(match):
        nonlocal position
        index, flags, width, precision, conversion = match.groups()
        if conversion == "%":
            return "%"
        if conversion == "n":
            return "\n"
        if index:
            value = args[int(index[:-1]) - 1]
        else:
            if position >= len(args):
                raise Unresolved("String.format with too few arguments")
            value = args[position]
            position += 1
        if conversion in "sS":
            text = jstr(value)
        elif conversion == "d":
            text = str(int(value))
        elif conversion == "f":
            text = f"{float(value):.{int(precision[1:]) if precision else 6}f}"
        elif conversion == "b":
            text = jstr(bool(value))
        else:
            raise Unresolved(f"String.format conversion %{conversion}")
        if width:
            text = text.ljust(int(width)) if "-" in (flags or "") else text.rjust(int(width))
        return text.upper() if conversion == "S" else text

    return re.sub(r"%(\d+\$)?([-#+ 0,(<]*)?(\d+)?(\.\d+)?([a-zA-Z%])", one, fmt)


def java_string(node):
    """A string literal's value, text blocks included."""
    raw = node.text.decode()
    if raw.startswith('"""'):
        body = raw[3:-3].split("\n", 1)[1] if "\n" in raw[3:-3] else raw[3:-3]
        lines = body.split("\n")
        indent = min((len(l) - len(l.lstrip()) for l in lines if l.strip()), default=0)
        body = "\n".join(l[indent:].rstrip() for l in lines)
    else:
        body = raw[1:-1]
    return re.sub(
        r"\\(u+[0-9a-fA-F]{4}|[0-7]{1,3}|.)",
        lambda m: (
            chr(int(m.group(1).lstrip("u"), 16)) if m.group(1)[0] == "u"
            else chr(int(m.group(1), 8)) if m.group(1)[0] in "01234567"
            else {"n": "\n", "t": "\t", "r": "\r", "b": "\b", "f": "\f", "s": " ", "\n": ""}.get(m.group(1), m.group(1))
        ),
        body,
    )


# ---------------------------------------------------------------------------------------------
# The sources: every class of the checkout's main and test trees, by simple name.


class ClassInfo:
    def __init__(self, node, package, path, outer=None, static_imports=()):
        self.node, self.package, self.path, self.outer = node, package, path, outer
        self.static_imports = static_imports
        self.name = node.child_by_field_name("name").text.decode()
        superclass = node.child_by_field_name("superclass")
        self.superclass = None
        if superclass is not None:
            type_node = kids(superclass)[0]
            if type_node.type == "generic_type":
                type_node = kids(type_node)[0]
            self.superclass = type_node.text.decode().split(".")[-1]
        self.fields, self.methods, self.enum_constants, self.nested = {}, {}, [], {}
        self.constant_nodes = {}
        self.annotations = annotations(node)
        self.abstract = has_modifier(node, "abstract")
        body = node.child_by_field_name("body")
        for child in kids(body) if body is not None else ():
            if child.type == "field_declaration":
                for declarator in child.children_by_field_name("declarator"):
                    self.fields[declarator.child_by_field_name("name").text.decode()] = (child, declarator)
            elif child.type == "method_declaration":
                self.methods.setdefault(child.child_by_field_name("name").text.decode(), []).append(child)
            elif child.type == "enum_body_declarations":
                for inner in kids(child):
                    if inner.type == "method_declaration":
                        self.methods.setdefault(inner.child_by_field_name("name").text.decode(), []).append(inner)
            elif child.type == "enum_constant":
                self.enum_constants.append(child.child_by_field_name("name").text.decode())
                self.constant_nodes[self.enum_constants[-1]] = child


def kids(node):
    """Named children, comments left out: a comment inside an array initializer is a child too."""
    return [c for c in node.named_children if c.type not in ("line_comment", "block_comment")]


def has_modifier(node, word):
    for child in node.children:
        if child.type == "modifiers":
            return any(c.type == word for c in child.children)
    return False


def annotations(node):
    """{name: {key: value node}} of a declaration's annotations."""
    found = {}
    for child in node.children:
        if child.type != "modifiers":
            continue
        for annotation in kids(child):
            if annotation.type not in ("annotation", "marker_annotation"):
                continue
            name = annotation.child_by_field_name("name").text.decode().split(".")[-1]
            values = {}
            arguments = annotation.child_by_field_name("arguments")
            for pair in kids(arguments) if arguments is not None else ():
                if pair.type == "element_value_pair":
                    values[pair.child_by_field_name("key").text.decode()] = pair.child_by_field_name("value")
                else:
                    values["value"] = pair
            found[name] = values
    return found


class Sources:
    def __init__(self, roots):
        import tree_sitter
        import tree_sitter_java

        self.parser = tree_sitter.Parser(tree_sitter.Language(tree_sitter_java.language()))
        self.classes = {}
        for root in roots:
            for path in sorted(Path(root).rglob("*.java")):
                self.load(path)

    def load(self, path):
        tree = self.parser.parse(path.read_bytes())
        package, static_imports = "", []
        for child in kids(tree.root_node):
            if child.type == "package_declaration":
                package = kids(child)[0].text.decode()
            if child.type == "import_declaration" and any(c.type == "static" for c in child.children):
                parts = kids(child)[0].text.decode().split(".")
                star = any(c.type == "asterisk" for c in child.children)
                static_imports.append((parts[-1], "*") if star else (parts[-2], parts[-1]))

        def walk(node, outer):
            for child in kids(node):
                if child.type in ("class_declaration", "enum_declaration", "interface_declaration"):
                    info = ClassInfo(child, package, path, outer, static_imports)
                    self.classes.setdefault(info.name, info)
                    if outer is not None:
                        outer.nested[info.name] = info
                    body = child.child_by_field_name("body")
                    if body is not None:
                        walk(body, info)
                        for inner in kids(body):
                            if inner.type == "enum_body_declarations":
                                walk(inner, info)

        walk(tree.root_node, None)

    def chain(self, info):
        """The class and its superclasses that the sources hold, nearest first."""
        seen = []
        while info is not None and info not in seen:
            seen.append(info)
            info = self.classes.get(info.superclass) if info.superclass else None
        return seen


# ---------------------------------------------------------------------------------------------
# The evaluator.


class Run:
    """One test method's execution: the temporary files it made and the commands it ran."""

    def __init__(self, test, sources):
        self.test, self.sources = test, sources
        self.counter, self.touch, self.mkdir, self.invocations = 0, [], [], []
        self.steps = 0
        self.instance = {}
        # The methods being evaluated: an overload of the same arity (`f(String)` calling
        # `f(File)`) cannot be told apart without types, and must not be taken for itself.
        self.stack = []
        # Temporary paths handed to a statement this cannot evaluate before a tool reads them:
        # what the test wrote there (a SAM built in code, a copied input) is not on either side.
        self.tainted = {}
        self.output_dir = None
        # The names of the test class's methods that reach a runner (`reaching_runners`).
        self.reaching = set()

    def temp(self, prefix, suffix, directory=False, parent=None):
        self.counter += 1
        name = re.sub(r"[^\w.+-]", "_", f"{jstr(prefix)}{self.counter}{'' if suffix is None else jstr(suffix)}")
        path = f"{OUT if parent is None else jstr(parent)}/{name}"
        (self.mkdir if directory else self.touch).append(path)
        return path

    def temp_output_dir(self):
        """`CommandLineProgramTest.getTempOutputDir()`: one directory per test class."""
        if self.output_dir is None:
            self.output_dir = self.temp("tempOutputDir", "", directory=True)
        return self.output_dir


class Interpreter:
    def __init__(self, run, this):
        self.run, self.this, self.sources = run, this, run.sources

    # --- names --------------------------------------------------------------------------------

    def tool_name(self):
        """`getCommandLineProgramName()`, which every `CommandLineProgramTest` defines."""
        for info in self.sources.chain(self.this):
            if info.name in FRAMEWORK:
                break
            for method in info.methods.get("getCommandLineProgramName", []):
                if method.child_by_field_name("body") is not None:
                    return jstr(self.call(method, [], info))
        name = self.this.name
        return name[: -len("Test")] if name.endswith("Test") else name

    def overridden(self, name, args):
        """A test class's own definition of a framework helper, which wins as it does in Java."""
        for owner in self.sources.chain(self.this):
            if owner.name in FRAMEWORK:
                return NOT_NATIVE
            for method in owner.methods.get(name, []):
                if len(kids(method.child_by_field_name("parameters"))) == len(args) and method not in self.run.stack:
                    return self.call(method, args, owner)
        return NOT_NATIVE

    def field(self, info, name):
        """A field's value, looked up through the class chain and the enclosing classes."""
        for owner in self.sources.chain(info):
            key = (owner.name, name)
            if key in self.run.instance:
                return self.run.instance[key]
            if name in owner.enum_constants:
                return EnumValue(name, owner)
            if name in owner.fields:
                declaration, declarator = owner.fields[name]
                value_node = declarator.child_by_field_name("value")
                if value_node is None:
                    raise Unresolved(f"field {owner.name}.{name} has no initializer")
                try:
                    value = Interpreter(self.run, owner).eval(value_node, {"scopes": []})
                except (Unresolved, Return, Thrown, ArithmeticError, LookupError, TypeError, ValueError, AttributeError) as error:
                    value = Poison(str(error))
                self.run.instance[key] = value
                return value
        if info is not None and info.outer is not None:
            return self.field(info.outer, name)
        for owner, member in info.static_imports if info is not None else ():
            if member in (name, "*") and owner in self.sources.classes and self.sources.classes[owner] is not info:
                try:
                    return self.field(self.sources.classes[owner], name)
                except Unresolved:
                    continue
        raise Unresolved(f"unknown name {name}")

    def lookup(self, name, env):
        for scope in reversed(env["scopes"]):
            if name in scope:
                return scope[name]
        try:
            return self.field(self.this, name)
        except Unresolved:
            if name in self.sources.classes:
                return ClassRef(name, self.sources.classes[name])
            if name[:1].isupper() and not name.isupper():
                return ClassRef(name)
            raise

    # --- statements ---------------------------------------------------------------------------

    def call(self, method, args, owner, scopes=None):
        params = method.child_by_field_name("parameters")
        scope = {}
        names = []
        for param in kids(params) if params is not None else ():
            if param.type == "formal_parameter":
                names.append((param.child_by_field_name("name").text.decode(), False))
            elif param.type == "spread_parameter":
                names.append((kids(param)[-1].child_by_field_name("name").text.decode(), True))
        for i, (name, spread) in enumerate(names):
            if spread:
                rest = args[i:]
                scope[name] = rest[0] if len(rest) == 1 and isinstance(rest[0], list) else rest
            elif i < len(args):
                scope[name] = args[i]
        body = method.child_by_field_name("body")
        if body is None:
            raise Unresolved(f"abstract method {owner.name}.{method.child_by_field_name('name').text.decode()}")
        inner = Interpreter(self.run, owner if has_modifier(method, "static") else self.this)
        env = {"scopes": (scopes or []) + [scope]}
        self.run.stack.append(method)
        try:
            inner.block(body, env)
        except Return as done:
            return done.value
        finally:
            self.run.stack.pop()
        return None

    def block(self, node, env):
        env["scopes"].append({})
        try:
            for statement in kids(node):
                self.statement(statement, env)
        finally:
            env["scopes"].pop()

    def bind(self, env, name, value):
        env["scopes"][-1][name] = value

    def assign(self, env, name, value):
        for scope in reversed(env["scopes"]):
            if name in scope:
                scope[name] = value
                return
        for owner in self.sources.chain(self.this):
            if name in owner.fields:
                self.run.instance[(owner.name, name)] = value
                return
        self.bind(env, name, value)

    def taint(self, node, env, reason):
        """The temporary paths a statement we cannot evaluate is handed: it may write them."""
        for child in walk(node):
            if child.type != "identifier":
                continue
            try:
                value = self.lookup(child.text.decode(), env)
            except Unresolved:
                continue
            for item in value if isinstance(value, list) else [value]:
                if isinstance(item, str) and (item.startswith(OUT + "/") or missing_from_checkout(item)):
                    snippet = " ".join(node.text.decode().split())[:120]
                    self.run.tainted.setdefault(item, f"{reason}, in `{snippet}`")

    def poison(self, node, env, reason):
        """Every variable a statement we cannot evaluate might change."""
        self.taint(node, env, reason)
        for child in walk(node):
            target = None
            if child.type == "assignment_expression":
                target = child.child_by_field_name("left")
            elif child.type == "method_invocation":
                target = child.child_by_field_name("object")
            elif child.type == "variable_declarator":
                target = child.child_by_field_name("name")
            if target is not None and target.type == "identifier":
                self.assign(env, target.text.decode(), Poison(reason))
        if any(c.type == "method_invocation" and c.child_by_field_name("name").text.decode() in RUNNERS | self.run.reaching
               for c in walk(node)):
            raise Unresolved(reason)

    def statement(self, node, env):
        self.run.steps += 1
        if self.run.steps > 200_000:
            raise Unresolved("evaluation budget exhausted")
        kind = node.type
        if kind in ("line_comment", "block_comment", "empty_statement", "assert_statement",
                    "local_class_declaration", "class_declaration", "labeled_statement_"):
            return
        if kind == "block":
            return self.block(node, env)
        if kind == "local_variable_declaration":
            for declarator in node.children_by_field_name("declarator"):
                name = declarator.child_by_field_name("name").text.decode()
                value = declarator.child_by_field_name("value")
                try:
                    self.bind(env, name, None if value is None else self.eval(value, env))
                except Unresolved as error:
                    self.taint(value, env, str(error))
                    self.bind(env, name, Poison(str(error)))
                    if value is not None and any(
                        c.type == "method_invocation" and c.child_by_field_name("name").text.decode() in RUNNERS | self.run.reaching
                        for c in walk(value)
                    ):
                        raise
            return
        if kind == "expression_statement":
            try:
                self.eval(kids(node)[0], env)
            except Unresolved as error:
                self.poison(node, env, str(error))
            return
        if kind == "if_statement":
            try:
                condition = truthy(self.eval(node.child_by_field_name("condition"), env))
            except Unresolved as error:
                return self.poison(node, env, f"branch on {error}")
            branch = node.child_by_field_name("consequence" if condition else "alternative")
            if branch is not None:
                self.statement(branch, env)
            return
        if kind == "enhanced_for_statement":
            try:
                items = self.eval(node.child_by_field_name("value"), env)
                if isinstance(items, Poison):
                    raise Unresolved(items.reason)
                items = list(items)
            except (Unresolved, TypeError) as error:
                return self.poison(node, env, f"loop over {error}")
            name = node.child_by_field_name("name").text.decode()
            for item in items:
                env["scopes"].append({name: item})
                try:
                    self.statement(node.child_by_field_name("body"), env)
                finally:
                    env["scopes"].pop()
            return
        if kind == "for_statement":
            env["scopes"].append({})
            try:
                for init in node.children_by_field_name("init"):
                    self.statement(init, env) if init.type.endswith("declaration") else self.eval(init, env)
                for _ in range(10_000):
                    condition = node.child_by_field_name("condition")
                    if condition is not None and not truthy(self.eval(condition, env)):
                        break
                    self.statement(node.child_by_field_name("body"), env)
                    for update in node.children_by_field_name("update"):
                        self.eval(update, env)
            except Unresolved as error:
                self.poison(node, env, f"loop on {error}")
            finally:
                env["scopes"].pop()
            return
        if kind == "try_statement":
            try:
                self.statement(node.child_by_field_name("body"), env)
            finally:
                for child in kids(node):
                    if child.type == "finally_clause":
                        self.statement(kids(child)[-1], env)
            return
        if kind == "try_with_resources_statement":
            env["scopes"].append({})
            try:
                resources = node.child_by_field_name("resources")
                for resource in kids(resources):
                    if resource.type == "resource" and resource.child_by_field_name("name") is not None:
                        try:
                            value = self.eval(resource.child_by_field_name("value"), env)
                        except Unresolved as error:
                            value = Poison(str(error))
                        self.bind(env, resource.child_by_field_name("name").text.decode(), value)
                self.statement(node.child_by_field_name("body"), env)
            finally:
                env["scopes"].pop()
            return
        if kind == "return_statement":
            raise Return(self.eval(kids(node)[0], env) if kids(node) else None)
        if kind == "throw_statement":
            raise Thrown()
        if kind in ("while_statement", "do_statement", "switch_expression", "switch_statement",
                    "synchronized_statement", "yield_statement", "break_statement", "continue_statement"):
            return self.poison(node, env, f"{kind} is not evaluated")
        raise Unresolved(f"statement {kind}")

    # --- expressions --------------------------------------------------------------------------

    def eval(self, node, env):
        kind = node.type
        if kind == "string_literal":
            return java_string(node)
        if kind == "character_literal":
            return java_string(node)
        if kind in ("decimal_integer_literal", "hex_integer_literal", "octal_integer_literal", "binary_integer_literal"):
            return int(node.text.decode().rstrip("lL").replace("_", ""), 0 if kind != "decimal_integer_literal" else 10)
        if kind == "decimal_floating_point_literal":
            return float(node.text.decode().rstrip("fFdD").replace("_", ""))
        if kind == "true":
            return True
        if kind == "false":
            return False
        if kind == "null_literal":
            return None
        if kind == "this":
            return ClassRef("this", self.this)
        if kind == "identifier":
            value = self.lookup(node.text.decode(), env)
            if isinstance(value, Poison):
                raise Unresolved(value.reason)
            return value
        if kind in ("parenthesized_expression", "cast_expression"):
            return self.eval(kids(node)[-1] if kind == "cast_expression" else kids(node)[0], env)
        if kind == "class_literal":
            name = kids(node)[0].text.decode().split(".")[-1].split("<")[0]
            return ClassRef(name, self.sources.classes.get(name))
        if kind == "field_access":
            return self.field_access(node, env)
        if kind == "array_access":
            return self.eval(node.child_by_field_name("array"), env)[int(self.eval(node.child_by_field_name("index"), env))]
        if kind == "binary_expression":
            return self.binary(node, env)
        if kind == "unary_expression":
            operand = self.eval(node.child_by_field_name("operand"), env)
            operator = node.child_by_field_name("operator").text.decode()
            return {"!": lambda v: not v, "-": lambda v: -v, "+": lambda v: v}[operator](operand)
        if kind == "ternary_expression":
            branch = "consequence" if truthy(self.eval(node.child_by_field_name("condition"), env)) else "alternative"
            return self.eval(node.child_by_field_name(branch), env)
        if kind == "assignment_expression":
            left, right = node.child_by_field_name("left"), node.child_by_field_name("right")
            operator = node.child_by_field_name("operator").text.decode()
            value = self.eval(right, env)
            if operator == "+=":
                current = self.eval(left, env)
                value = current + value if isinstance(current, (int, float)) and isinstance(value, (int, float)) else jstr(current) + jstr(value)
            elif operator != "=":
                raise Unresolved(f"assignment {operator}")
            if left.type == "identifier":
                self.assign(env, left.text.decode(), value)
            elif left.type == "field_access" and left.child_by_field_name("object").type == "this":
                self.assign(env, left.child_by_field_name("field").text.decode(), value)
            elif left.type == "array_access":
                self.eval(left.child_by_field_name("array"), env)[int(self.eval(left.child_by_field_name("index"), env))] = value
            elif left.type == "field_access" and isinstance(self.eval(left.child_by_field_name("object"), env), Tool):
                tool = self.eval(left.child_by_field_name("object"), env)
                tool.fields[left.child_by_field_name("field").text.decode()] = value
            else:
                raise Unresolved("assignment to an expression")
            return value
        if kind == "update_expression":
            target = next(c for c in kids(node))
            delta = 1 if "++" in node.text.decode() else -1
            value = self.eval(target, env)
            self.assign(env, target.text.decode(), value + delta)
            return value
        if kind == "array_initializer":
            return [self.eval(child, env) for child in kids(node)]
        if kind == "array_creation_expression":
            value = node.child_by_field_name("value")
            if value is not None:
                return self.eval(value, env)
            dimensions = node.child_by_field_name("dimensions")
            sizes = [c for c in kids(node) if c.type == "dimensions_expr"]
            if sizes:
                return [None] * int(self.eval(sizes[0].named_children[0], env))
            return []
        if kind == "object_creation_expression":
            return self.construct(node, env)
        if kind == "method_invocation":
            return self.invoke(node, env)
        if kind == "lambda_expression":
            return Lambda(node, env, self.this)
        if kind == "instanceof_expression":
            raise Unresolved("instanceof")
        raise Unresolved(f"expression {kind}")

    def binary(self, node, env):
        operator = node.child_by_field_name("operator").text.decode()
        left = self.eval(node.child_by_field_name("left"), env)
        if operator == "&&":
            return truthy(left) and truthy(self.eval(node.child_by_field_name("right"), env))
        if operator == "||":
            return truthy(left) or truthy(self.eval(node.child_by_field_name("right"), env))
        right = self.eval(node.child_by_field_name("right"), env)
        if operator == "+":
            if isinstance(left, (int, float)) and isinstance(right, (int, float)) and not isinstance(left, bool):
                return left + right
            return jstr(left) + jstr(right)
        if operator in ("==", "!="):
            same = left is right if not isinstance(left, (int, float, str, type(None), bool)) else left == right
            if isinstance(left, EnumValue) and isinstance(right, EnumValue):
                same = left.name == right.name
            return same if operator == "==" else not same
        numbers = {"-": lambda a, b: a - b, "*": lambda a, b: a * b, "<": lambda a, b: a < b,
                   ">": lambda a, b: a > b, "<=": lambda a, b: a <= b, ">=": lambda a, b: a >= b,
                   "%": lambda a, b: a % b,
                   "/": lambda a, b: a // b if isinstance(a, int) and isinstance(b, int) else a / b}
        if operator in numbers:
            return numbers[operator](left, right)
        raise Unresolved(f"operator {operator}")

    def field_access(self, node, env):
        target, name = node.child_by_field_name("object"), node.child_by_field_name("field").text.decode()
        if target.type in ("this", "super"):
            return self.field(self.this, name)
        try:
            owner = self.eval(target, env)
        except Unresolved:
            # A qualified class name: `org.x.Tool.CONSTANT`.
            owner = ClassRef(target.text.decode().split(".")[-1], self.sources.classes.get(target.text.decode().split(".")[-1]))
        if isinstance(owner, ClassRef):
            if (owner.name, name) in (("File", "separator"), ("File", "pathSeparator")):
                return "/" if name == "separator" else ":"
            if owner.info is None:
                raise Unresolved(f"constant {owner.name}.{name} is outside the sources")
            if name in owner.info.nested:
                return ClassRef(name, owner.info.nested[name])
            return self.field(owner.info, name)
        if name == "length" and isinstance(owner, (list, str)):
            return len(owner)
        if isinstance(owner, EnumValue):
            return self.enum_field(owner, name)
        if isinstance(owner, Tool):
            if name in owner.fields:
                return owner.fields[name]
            raise Unresolved(f"{owner.name}.{name} read before the test set it")
        raise Unresolved(f"field {name} of a {type(owner).__name__}")

    def enum_field(self, constant, name):
        """A field of an enum constant, as its constructor sets it from the constant's arguments."""
        info = constant.info
        node = info.constant_nodes.get(constant.name) if info is not None else None
        if node is None:
            raise Unresolved(f"field {name} of enum constant {constant.name}")
        key = (info.name, constant.name)
        if key not in self.run.instance:
            args = self.arguments(node, {"scopes": []})
            fields = {}
            for child in walk(info.node):
                if child.type == "constructor_declaration":
                    params = kids(child.child_by_field_name("parameters"))
                    if len(params) == len(args):
                        scope = {p.child_by_field_name("name").text.decode(): a for p, a in zip(params, args)}
                        for statement in walk(child.child_by_field_name("body")):
                            if statement.type == "assignment_expression":
                                left = statement.child_by_field_name("left")
                                field = left.child_by_field_name("field") if left.type == "field_access" else left
                                try:
                                    fields[field.text.decode()] = Interpreter(self.run, info).eval(
                                        statement.child_by_field_name("right"), {"scopes": [scope]})
                                except Unresolved as error:
                                    fields[field.text.decode()] = Poison(str(error))
                        break
            self.run.instance[key] = fields
        fields = self.run.instance[key]
        if name not in fields:
            return self.field(info, name)
        if isinstance(fields[name], Poison):
            raise Unresolved(fields[name].reason)
        return fields[name]

    def arguments(self, node, env):
        arguments = node.child_by_field_name("arguments")
        return [self.eval(child, env) for child in kids(arguments)] if arguments is not None else []

    def construct(self, node, env):
        type_node = node.child_by_field_name("type")
        name = type_node.text.decode().split("<")[0].split(".")[-1]
        if kids(node)[-1].type == "class_body":
            raise Unresolved(f"anonymous subclass of {name}")
        args = self.arguments(node, env)
        if name == "File":
            return normalized(jstr(args[0])) if len(args) == 1 else join_path(args[0], args[1])
        if name in ("PicardHtsPath", "HtsPath", "URI"):
            return jstr(args[0])
        if name in ("ArrayList", "LinkedList", "HashSet", "LinkedHashSet", "TreeSet"):
            items = list(args[0]) if args and isinstance(args[0], (list, tuple)) else []
            return sorted(items) if name == "TreeSet" else items
        if name in ("String", "StringBuilder"):
            return jstr(args[0]) if args else ""
        if name == "PicardCommandLine":
            return ClassRef("PicardCommandLine")
        info = self.sources.classes.get(name)
        if info is not None and any(c.name == "CommandLineProgram" for c in self.sources.chain(info)):
            # `PicardCommandLine` names a tool by its class's simple name.
            return Tool(name)
        if info is not None and name.endswith("Test"):
            return ClassRef("this", info)
        if info is not None and any(c.name == "SamFileTester" for c in self.sources.chain(info)):
            raise Unresolved(f"new {name}: a SamFileTester, whose input SAMRecordSetBuilder builds in code")
        raise Unresolved(f"new {name}")

    def invoke(self, node, env):
        name = node.child_by_field_name("name").text.decode()
        target_node = node.child_by_field_name("object")
        arguments = node.child_by_field_name("arguments")
        arg_nodes = kids(arguments) if arguments is not None else []
        if target_node is None or target_node.type in ("this", "super"):
            return self.own_method(name, arg_nodes, env, super_call=target_node is not None and target_node.type == "super")
        try:
            target = self.eval(target_node, env)
        except Unresolved:
            if target_node.type in ("identifier", "field_access", "scoped_identifier") and target_node.text.decode()[:1].isupper():
                target = ClassRef(target_node.text.decode().split(".")[-1], self.sources.classes.get(target_node.text.decode().split(".")[-1]))
            else:
                self.lambdas(arg_nodes, env)
                raise
        if isinstance(target, ClassRef) and target.name == "this":
            return Interpreter(self.run, target.info).own_method(name, arg_nodes, env, interpreter=self)
        args = [self.eval(child, env) for child in arg_nodes] if name not in LAZY else None
        try:
            return self.member(target, name, args, arg_nodes, env)
        except Unresolved:
            self.lambdas(arg_nodes, env)
            raise

    def lambdas(self, arg_nodes, env):
        """Run the zero-argument lambdas of a call we cannot evaluate: `assertThrows(X.class, () -> runCommandLine(args))`."""
        for child in arg_nodes:
            if child.type == "lambda_expression":
                self.run_lambda(Lambda(child, env, self.this))

    def run_lambda(self, fn, args=()):
        body = fn.node.child_by_field_name("body")
        params = fn.node.child_by_field_name("parameters")
        names = [n.text.decode() for n in walk(params) if n.type == "identifier"] if params is not None else []
        names = names[-len(args):] if args else []
        env = {"scopes": fn.env["scopes"] + [dict(zip(names, args))]}
        inner = Interpreter(self.run, fn.this)
        try:
            if body.type == "block":
                inner.block(body, env)
                return None
            return inner.eval(body, env)
        except Return as done:
            return done.value
        except Thrown:
            return None

    def own_method(self, name, arg_nodes, env, super_call=False, interpreter=None):
        caller = interpreter or self
        args = [caller.eval(child, env) for child in arg_nodes]
        if name == "runPicardCommandLine":
            return self.runner(name, args)
        native = self.native(name, args)
        if native is not NOT_NATIVE:
            return native
        chain = self.sources.chain(self.this)
        if super_call:
            chain = chain[1:]
        for owner in chain + ([self.this.outer] if self.this.outer else []):
            for method in owner.methods.get(name, []):
                params = method.child_by_field_name("parameters")
                count = len([p for p in kids(params) if p.type in ("formal_parameter", "spread_parameter")])
                spread = any(p.type == "spread_parameter" for p in kids(params))
                if count == len(args) or (spread and len(args) >= count - 1):
                    if method.child_by_field_name("body") is None or method in self.run.stack:
                        continue
                    return Interpreter(self.run, self.this).call(method, args, owner)
        raise Unresolved(f"method {name}/{len(args)}")

    def native(self, name, args):
        """The test framework's own helpers, unless the test class overrides them."""
        if name not in ("getCommandLineProgramName", "getClass"):
            value = self.overridden(name, args)
            if value is not NOT_NATIVE:
                return value
        if name == "getTempOutputDir":
            return self.run.temp_output_dir()
        if name == "getTempOutputFile":
            return self.run.temp(args[0], args[1], parent=self.run.temp_output_dir())
        if name == "makePicardCommandLineArgs":
            if len(args) == 2:
                return [jstr(args[0])] + [jstr(a) for a in args[1]]
            if isinstance(args[0], list):
                return [self.tool_name()] + [jstr(a) for a in args[0]]
            raise Unresolved("makePicardCommandLineArgs(Map)")
        if name == "getCommandLineProgramName":
            return self.tool_name()
        if name == "getClass":
            return ClassRef(self.this.name, self.this)
        return NOT_NATIVE

    def runner(self, name, args):
        """`runPicardCommandLine` and its relatives: the invocation this harness exists to read."""
        if name == "runPicardCommandLine":
            if len(args) == 2:
                tool, tokens = jstr(args[0]), args[1]
            elif isinstance(args[0], list):
                tool, tokens = self.tool_name(), args[0]
            else:
                raise Unresolved("runPicardCommandLine(Map)")
            self.record(tool, [jstr(t) for t in tokens])
            return 0
        raise Unresolved(f"runner {name}")

    def record(self, tool, tokens, via=None):
        tokens = [RELATIVE_RESOURCE.sub(f"{ROOT}/testdata/", jstr(t)) for t in tokens]
        for path, reason in self.run.tainted.items():
            if any(path in t for t in tokens):
                raise Unresolved(f"{tool} reads {path}, which the test writes in code not evaluated ({reason})")
        # A `testdata/` path the checkout does not hold is one the test writes (an output beside its
        # inputs) or one it expects to be missing; the checkout is read-only on both sides, so it
        # moves to the output directory, under the same relative path, as `moved` records.
        moved = sorted({m for t in tokens for m in re.findall(rf"{ROOT}/(testdata/[^\s,;:=]+)", t)
                        if missing_from_checkout(m)})
        for path in moved:
            tokens = [t.replace(f"{ROOT}/{path}", f"{OUT}/{path}") for t in tokens]
            parent = f"{OUT}/{path}".rsplit("/", 1)[0]
            if parent not in self.run.mkdir:
                self.run.mkdir.append(parent)
        self.run.invocations.append({
            "tool": tool, "args": tokens, **({"via": via} if via else {}),
            **({"moved": moved} if moved else {}),
            "touch": list(self.run.touch), "mkdir": list(self.run.mkdir),
        })

    def member(self, target, name, args, arg_nodes, env):
        if isinstance(target, Poison):
            raise Unresolved(target.reason)
        if isinstance(target, Tool):
            if name == "instanceMain":
                # The fields the test set are the values the parser starts from: the command line.
                self.record(target.name, target.tokens() + [jstr(t) for t in args[0]])
                return 0
            if name == "doWork":
                self.record(target.name, target.tokens(), via="doWork")
                return 0
            raise Unresolved(f"{target.name}.{name}")
        if isinstance(target, ClassRef) and target.name == "PicardCommandLine" and name == "instanceMain":
            tokens = [jstr(t) for t in args[0]]
            if not tokens:
                raise Unresolved("PicardCommandLine.instanceMain with no tool")
            self.record(tokens[0], tokens[1:])
            return 0
        if isinstance(target, Optional):
            if name in ("isPresent", "isEmpty"):
                return (target.value is not None) == (name == "isPresent")
            if name == "get" and target.value is not None:
                return target.value
            if name == "orElse":
                return target.value if target.value is not None else args[0]
            raise Unresolved(f"Optional.{name}")
        if isinstance(target, Lambda):
            return self.run_lambda(target, args)
        if isinstance(target, ClassRef):
            return self.static(target, name, args)
        if isinstance(target, EnumValue):
            if name in ("name", "toString"):
                return target.name
            raise Unresolved(f"enum method {name}")
        if isinstance(target, str):
            return string_method(target, name, args)
        if isinstance(target, list) and name == "forEach" and isinstance(args[0], Lambda):
            for item in list(target):
                self.run_lambda(args[0], [item])
            return None
        if isinstance(target, list):
            return list_method(target, name, args)
        if isinstance(target, (int, float)) and name in ("toString", "intValue", "doubleValue", "longValue"):
            return jstr(target) if name == "toString" else target
        raise Unresolved(f"method {name} of a {type(target).__name__}")

    def static(self, target, name, args):
        cls = target.name
        if cls == "Class" and name == "forName":
            simple = jstr(args[0]).split(".")[-1]
            return ClassRef(simple, self.sources.classes.get(simple))
        if name == "getDeclaredConstructor":
            return target
        if name == "newInstance":
            info = target.info
            if info is not None and any(c.name == "CommandLineProgram" for c in self.sources.chain(info)):
                return Tool(target.name)
            raise Unresolved(f"{cls}.newInstance")
        if cls == "Optional" and name in ("empty", "of", "ofNullable"):
            return Optional(args[0] if args else None)
        if cls == "String" and name == "format":
            return jformat(jstr(args[0]), args[1:])
        if cls == "String" and name == "join":
            items = args[1] if len(args) == 2 and isinstance(args[1], list) else args[1:]
            return jstr(args[0]).join(jstr(v) for v in items)
        if cls == "String" and name == "valueOf" or cls in ("Integer", "Long", "Double", "Boolean") and name == "toString":
            return jstr(args[0])
        if cls in ("Integer", "Long") and name in ("parseInt", "parseLong", "valueOf"):
            return int(jstr(args[0]))
        if cls == "Double" and name in ("parseDouble", "valueOf"):
            return float(jstr(args[0]))
        if (cls, name) in (("Arrays", "asList"), ("List", "of"), ("ImmutableSet", "of"), ("Stream", "concat"), ("Lists", "newArrayList"), ("Set", "of"),
                           ("Collections", "singletonList"), ("Collections", "singleton"), ("ImmutableList", "of"),
                           ("Stream", "of"), ("Sets", "newHashSet")):
            if name == "concat":
                return list(args[0]) + list(args[1])
            return list(args[0]) if len(args) == 1 and isinstance(args[0], list) and cls in ("Arrays", "Lists", "Stream") else list(args)
        if (cls, name) in (("Collections", "emptyList"), ("Collections", "emptySet")):
            return []
        if (cls, name) in (("Paths", "get"), ("Path", "of"), ("IOUtils", "getPath")):
            path = normalized(jstr(args[0]))
            for part in args[1:]:
                path = join_path(path, part)
            return path
        if cls in ("File", "IOUtil") and name in ("createTempFile", "newTempFile"):
            return self.run.temp(args[0], args[1], parent=args[2] if len(args) > 2 else None)
        if cls == "Files" and name == "createTempFile":
            return self.run.temp(*(args[1:3] if len(args) > 2 else args[:2]), parent=args[0] if len(args) > 2 else None)
        if cls == "IOUtil" and name in ("createTempDir", "newTempPath"):
            return self.run.temp(args[0] if args else "tmp", args[1] if len(args) > 1 else "", directory=True)
        if cls == "Files" and name == "createTempDirectory":
            if len(args) > 1:
                return self.run.temp(args[1], "", directory=True, parent=args[0])
            return self.run.temp(args[0] if args else "tmp", "", directory=True)
        if cls in ("ArrayUtils", "ObjectArrays") and name in ("addAll", "concat"):
            return [v for part in args if isinstance(part, list) for v in part]
        if name == "values" and target.info is not None and target.info.enum_constants:
            return [EnumValue(c, target.info) for c in target.info.enum_constants]
        if cls == "Assert" or name in ("assertThrows", "expectThrows"):
            for arg in args:
                if isinstance(arg, Lambda):
                    self.run_lambda(arg)
            return None
        if cls == "System" and name == "getProperty":
            if jstr(args[0]) == "user.dir":
                return ROOT
            if len(args) > 1:
                return args[1]
            raise Unresolved(f"system property {jstr(args[0])}")
        if cls == "StringUtils" and name == "join":
            return jstr(args[1] if len(args) > 1 else " ").join(jstr(v) for v in args[0])
        if name == "getSimpleName":
            return cls
        if name == "getName" and target.info is not None:
            return f"{target.info.package}.{cls}"
        if cls in FRAMEWORK or (target.info and name in target.info.methods):
            native = self.native(name, args)
            if native is not NOT_NATIVE:
                return native
            if target.info is not None:
                for method in target.info.methods.get(name, []):
                    params = method.child_by_field_name("parameters")
                    if len([p for p in kids(params) if p.type == "formal_parameter"]) == len(args) and method not in self.run.stack:
                        return Interpreter(self.run, target.info).call(method, args, target.info)
        raise Unresolved(f"{cls}.{name}")


NOT_NATIVE = object()
FRAMEWORK = {"CommandLineProgramTest"}
RUNNERS = {"runPicardCommandLine", "instanceMain", "doWork"}
LAZY = set()


def truthy(value):
    if isinstance(value, Poison):
        raise Unresolved(value.reason)
    if not isinstance(value, bool):
        raise Unresolved(f"condition on a {type(value).__name__}")
    return value


def walk(node):
    yield node
    for child in node.children:
        yield from walk(child)


def string_method(text, name, args):
    simple = {
        "toString": lambda: text, "trim": lambda: text.strip(), "getAbsolutePath": lambda: absolute(text),
        "getPath": lambda: text, "toPath": lambda: text, "toFile": lambda: text, "toAbsolutePath": lambda: absolute(text),
        "getCanonicalPath": lambda: absolute(text), "getAbsoluteFile": lambda: absolute(text),
        "getURIString": lambda: text, "getRawInputString": lambda: text, "toUri": lambda: absolute(text),
        "toURI": lambda: absolute(text), "getName": lambda: text.rstrip("/").rsplit("/", 1)[-1],
        "getFileName": lambda: text.rstrip("/").rsplit("/", 1)[-1], "isEmpty": lambda: text == "",
        "getParent": lambda: text.rstrip("/").rsplit("/", 1)[0], "getParentFile": lambda: text.rstrip("/").rsplit("/", 1)[0],
        "toLowerCase": lambda: text.lower(), "toUpperCase": lambda: text.upper(), "length": lambda: len(text),
        "intern": lambda: text, "delete": lambda: True, "deleteOnExit": lambda: None, "mkdirs": lambda: True,
        "mkdir": lambda: True,
    }
    if name in simple and not args:
        return simple[name]()
    if name in ("resolve", "resolveSibling"):
        base = text if name == "resolve" else text.rstrip("/").rsplit("/", 1)[0]
        return join_path(base, args[0])
    if name == "concat":
        return text + jstr(args[0])
    if name == "replace":
        return text.replace(jstr(args[0]), jstr(args[1]))
    if name == "replaceAll":
        return re.sub(jstr(args[0]), jstr(args[1]).replace("$", "\\"), text)
    if name in ("equals", "contentEquals"):
        return text == (args[0] if isinstance(args[0], str) else None)
    if name == "equalsIgnoreCase":
        return isinstance(args[0], str) and text.lower() == args[0].lower()
    if name in ("startsWith", "endsWith", "contains"):
        return {"startsWith": text.startswith, "endsWith": text.endswith, "contains": text.__contains__}[name](jstr(args[0]))
    if name == "split":
        return [p for p in re.split(jstr(args[0]), text)]
    if name == "substring":
        return text[int(args[0]): int(args[1]) if len(args) > 1 else None]
    if name == "formatted":
        return jformat(text, args)
    raise Unresolved(f"String.{name}")


def list_method(items, name, args):
    if name in ("add", "addAll") and len(args) == 1:
        if name == "addAll":
            items.extend(args[0])
        else:
            items.append(args[0])
        return True
    if name == "add" and len(args) == 2:
        items.insert(int(args[0]), args[1])
        return None
    if name in ("toArray", "iterator", "stream", "toList", "collect", "clone", "asList"):
        return list(items)
    if name == "size":
        return len(items)
    if name == "get":
        return items[int(args[0])]
    if name == "isEmpty":
        return not items
    if name == "contains":
        return any(jstr(v) == jstr(args[0]) for v in items)
    if name == "toString":
        return jstr(items)
    if name == "subList":
        return items[int(args[0]):int(args[1])]
    raise Unresolved(f"List.{name}")


# ---------------------------------------------------------------------------------------------
# Extraction.


def element_strings(node):
    if node is None:
        return []
    if node.type == "string_literal":
        return [java_string(node)]
    return [java_string(n) for n in walk(node) if n.type == "string_literal"]


def test_methods(sources, info):
    """(owner, method, @Test values) of every TestNG test the class runs, inherited ones included."""
    seen, found = set(), []
    class_test = None
    for owner in sources.chain(info):
        if class_test is None and "Test" in owner.annotations:
            class_test = owner.annotations["Test"]
        for name, methods in owner.methods.items():
            for method in methods:
                key = (name, len(kids(method.child_by_field_name("parameters"))))
                if key in seen:
                    continue
                seen.add(key)
                marks = annotations(method)
                if "Test" in marks:
                    found.append((owner, method, marks["Test"]))
                elif class_test is not None and has_modifier(method, "public") and not has_modifier(method, "static") \
                        and not any(m.startswith(("Before", "After", "DataProvider")) for m in marks) \
                        and method.child_by_field_name("type").type == "void_type":
                    found.append((owner, method, {}))
    return found, class_test or {}


def provider_rows(sources, info, test, values, interpreter):
    """The rows a `@DataProvider` returns, or [None] for a test that takes none."""
    if "dataProvider" not in values:
        return [None]
    provider = element_strings(values["dataProvider"])[0]
    owner_info = info
    if "dataProviderClass" in values:
        name = values["dataProviderClass"].text.decode().removesuffix(".class").split(".")[-1]
        owner_info = sources.classes.get(name)
        if owner_info is None:
            raise Unresolved(f"data provider class {name} is outside the sources")
    for owner in sources.chain(owner_info):
        for methods in owner.methods.values():
            for method in methods:
                marks = annotations(method)
                if "DataProvider" not in marks:
                    continue
                named = element_strings(marks["DataProvider"].get("name"))
                if (named[0] if named else method.child_by_field_name("name").text.decode()) == provider:
                    rows = Interpreter(interpreter.run, owner_info).call(method, [], owner)
                    if isinstance(rows, Poison):
                        raise Unresolved(rows.reason)
                    return [row if isinstance(row, list) else [row] for row in rows]
    raise Unresolved(f"data provider {provider} not found")


def setup_methods(sources, info):
    """`@BeforeClass`/`@BeforeMethod`/... methods, superclasses first."""
    found = []
    for owner in reversed(sources.chain(info)):
        for methods in owner.methods.values():
            for method in methods:
                if any(m in annotations(method) for m in ("BeforeClass", "BeforeMethod", "BeforeTest", "BeforeSuite")):
                    found.append((owner, method))
    return found


def extract_class(sources, info, relative, called=frozenset()):
    invocations, unextracted = [], []
    tests, class_test = test_methods(sources, info)
    reaching = reaching_runners(sources, info, called)
    # What the `@Before*` methods run, recorded once as TestNG runs it. A class whose setup runs a
    # tool is replayed as one group, setup first, so a test reading what it wrote reads it.
    setup = Run(f"{info.name}#setup", sources)
    setup.reaching = reaching | called
    for setup_owner, method in setup_methods(sources, info):
        try:
            Interpreter(setup, info).call(method, [], setup_owner)
        except (Unresolved, Thrown, RecursionError, ArithmeticError, LookupError, TypeError, ValueError, AttributeError):
            pass
    group = {"group": info.name} if setup.invocations else {}
    for k, invocation in enumerate(setup.invocations):
        invocations.append({"id": f"{info.name}#setup.{k}", "source": relative, **group, **invocation})
    names = [m.child_by_field_name("name").text.decode() for _, m, _ in tests]
    for owner, method, values in tests:
        name = method.child_by_field_name("name").text.decode()
        if name not in reaching:
            continue  # a unit test: it runs no tool, so it has no invocation to replay
        base = f"{info.name}#{name}"
        if names.count(name) > 1:  # overloads: told apart by their number of parameters
            base += f"/{len(kids(method.child_by_field_name('parameters')))}"
        merged = {**class_test, **values}
        groups = set(element_strings(merged.get("groups")))
        if merged.get("enabled") is not None and merged["enabled"].type == "false":
            unextracted.append({"id": base, "reason": "disabled (enabled = false)"})
            continue
        if groups & CLOUD_GROUPS:
            unextracted.append({"id": base, "reason": f"reads a bucket (groups {sorted(groups & CLOUD_GROUPS)})"})
            continue
        try:
            probe = Run(base, sources)
            rows = provider_rows(sources, info, method, merged, Interpreter(probe, info))
        except (Unresolved, Return, Thrown, RecursionError, ArithmeticError, LookupError, TypeError, ValueError, AttributeError) as error:
            unextracted.append({"id": base, "reason": f"data provider: {error}"})
            continue
        for index, row in enumerate(rows):
            test_id = base if row is None else f"{base}[{index}]"
            run = Run(test_id, sources)
            run.reaching = reaching | called
            interpreter = Interpreter(run, info)
            start = 0
            try:
                # Run again for the fields they set, which the test reads; what they ran is the
                # setup's, recorded above.
                for setup_owner, method_ in setup_methods(sources, info):
                    try:
                        interpreter.call(method_, [], setup_owner)
                    except (Unresolved, Thrown):
                        pass
                start = len(run.invocations)
                interpreter.call(method, row or [], owner)
            except Thrown:
                pass
            except (Unresolved, RecursionError, ArithmeticError, LookupError, TypeError, ValueError, AttributeError) as error:
                # The invocations recorded before it are whole: what the test ran first is kept, and
                # the failure is marked partial (most often the test's assertions on the output).
                missed = {"id": test_id, "reason": str(error) or type(error).__name__}
                unextracted.append({**missed, "partial": True} if len(run.invocations) > start else missed)
            else:
                if len(run.invocations) == start:
                    unextracted.append({"id": test_id, "reason": "reaches a runner, but no invocation was evaluated"})
            for k, invocation in enumerate(run.invocations[start:]):
                invocations.append({"id": f"{test_id}.{k}", "source": relative, **group, **invocation})
    return invocations, unextracted


RUNNER_CALL = re.compile(r"\b(runPicardCommandLine|instanceMain|doWork)\s*\(")


def reaching_runners(sources, info, called=frozenset()):
    """The names of the methods of a test class (inherited and enclosing ones included) whose body
    calls a runner, directly, through another such method, or through a `called` method of another
    test class (`tester.runTest()`): a test among them runs a tool; any other test is a unit test."""
    owners = sources.chain(info) + ([info.outer] if info.outer else [])
    bodies = {}
    for owner in owners:
        for name, methods in owner.methods.items():
            bodies.setdefault(name, []).extend(m.text.decode() for m in methods)
    elsewhere = re.compile(r"\.(" + "|".join(map(re.escape, sorted(called))) + r")\s*\(") if called else None
    reach, changed = set(), True
    while changed:
        changed = False
        # Only an own-method call: `x.get(` is not a call of this class's `get`.
        calls = re.compile(r"(?<![.\w])(" + "|".join(map(re.escape, reach)) + r")\s*\(") if reach else None
        for name, texts in bodies.items():
            if name not in reach and any(RUNNER_CALL.search(t) or (calls and calls.search(t))
                                         or (elsewhere and elsewhere.search(t)) for t in texts):
                reach.add(name)
                changed = True
    return reach


def called_elsewhere(sources, tests):
    """The names of the test sources' methods that reach a runner, by name across every test class,
    without the names a main-source class also declares: `x.build(` is not known to be a tester's."""
    test_bodies, main_names = {}, set()
    for info in sources.classes.values():
        if info.path.is_relative_to(tests):
            for name, methods in info.methods.items():
                test_bodies.setdefault(name, []).extend(m.text.decode() for m in methods)
        else:
            main_names.update(info.methods)
    reach, changed = set(), True
    while changed:
        changed = False
        calls = re.compile(r"\b(" + "|".join(map(re.escape, reach)) + r")\s*\(") if reach else None
        for name, texts in test_bodies.items():
            if name not in reach and any(RUNNER_CALL.search(t) or (calls and calls.search(t)) for t in texts):
                reach.add(name)
                changed = True
    return frozenset(reach - main_names)


CHECKOUT = None


def missing_from_checkout(path):
    """A `testdata/` path the checkout does not hold: the test writes it in code before a tool
    reads it (CollectMultipleMetricsTest's setup builds CollectGcBias2.bam), so both sides lack it."""
    if CHECKOUT is None or not re.match(rf"({ROOT}/)?testdata/", path):
        return False
    return not (CHECKOUT / path.removeprefix(ROOT + "/")).exists()


def extract(options):
    global CHECKOUT
    picard = Path(options.picard).resolve()
    CHECKOUT = picard
    roots = [picard / "src" / "main" / "java", picard / "src" / "test" / "java"]
    sources = Sources(roots + [Path(c) for c in options.constants])
    invocations, unextracted = [], []
    files = sorted((picard / "src" / "test" / "java").rglob("*Test.java"))
    called = called_elsewhere(sources, picard / "src" / "test" / "java")
    for path in files:
        info = sources.classes.get(path.stem)
        if info is None or info.abstract or info.path != path:
            continue
        found, missed = extract_class(sources, info, path.relative_to(picard).as_posix(), called)
        invocations += found
        unextracted += missed
    sha = subprocess.run(["git", "-C", str(picard), "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    payload = {
        "$comment": "Generated by tools/coverage/test_resources.py extract from the pinned Picard sources; do not edit.",
        "source": {"repository": "broadinstitute/picard", "tag": options.tag, "sha": sha, "classes": len(files)},
        "invocations": invocations,
        "unextracted": unextracted,
    }
    DATA.mkdir(exist_ok=True)
    INVOCATIONS.write_text(json.dumps(payload, indent=1) + "\n")
    tools = {i["tool"] for i in invocations}
    print(f"classes={len(files)} invocations={len(invocations)} tools={len(tools)} unextracted={len(unextracted)}")
    return 0


# ---------------------------------------------------------------------------------------------
# Replay.


def coverage_entries():
    """The tools whose port lives here, with how `run_array.py` compares each one's answer."""
    return {e["tool"]: e for e in json.loads(MANIFEST.read_text())["coverage"]["tools"]}


def exclusions():
    return json.loads(EXCLUSIONS.read_text()).get("invocations", {}) if EXCLUSIONS.is_file() else {}


def test_of(invocation_id):
    return invocation_id.rsplit(".", 1)[0]


def prepare(invocation, out_dir):
    """What the test's `createTempFile`s and `getTempOutputDir` had made before the tool ran."""
    for path in invocation["mkdir"]:
        (out_dir / path[len(OUT) + 1:]).mkdir(parents=True, exist_ok=True)
    for path in invocation["touch"]:
        made = out_dir / path[len(OUT) + 1:]
        made.parent.mkdir(parents=True, exist_ok=True)
        if not made.exists():
            made.touch()


def snapshot(directory):
    return {p: (p.stat().st_size, p.stat().st_mtime_ns) for p in directory.rglob("*") if p.is_file()}


def written(directory, before, strip_pg, back=lambda text: text):
    """Every file the invocation wrote or rewrote, by its path under the output directory, read as
    `run_array.read_output` reads one; a test's outputs have the names it chose, so none is named
    in advance, and a file an earlier invocation of the same test wrote is not this one's answer."""
    names = sorted(str(p.relative_to(directory)) for p, stamp in snapshot(directory).items() if before.get(p) != stamp)
    out = ""
    for name in names:
        raw = (directory / name).read_bytes()
        if raw[:4] == b"TIDX":
            answer = tribble_index(raw, back)
        else:
            answer = run_array.read_output(directory, strip_pg, name)
        out += f"== {name}\n{answer}\n"
    return out


def tribble_index(raw, back):
    """A tribble `.idx`, digested with its header's timestamp masked: the indexed file's
    modification time, which no two runs share, so the reference itself has no stable byte there
    (htsjdk-rs's index suite masks the same eight bytes). The indexed path it records is mapped
    back to the container's, as every other path the port was handed is."""
    end = raw.index(b"\0", 12)
    path = back(raw[12:end].decode("utf-8", "replace")).encode()
    rest = raw[end + 1:]
    canonical = raw[:12] + path + b"\0" + rest[:8] + b"\0" * 8 + rest[16:]
    return f"TRIBBLE-INDEX(timestamp masked) sha256={hashlib.sha256(canonical).hexdigest()} bytes={len(canonical)}"


def oracle(invocation, picard, out_dir, entry):
    # The tokens go to the jar as the test handed them, from the checkout's root as `user.dir`:
    # on `NAME=value` words Picard 3.5.0 picks its legacy parser itself, as it does for a user.
    cli = " ".join(shlex.quote(t) for t in invocation["args"])
    container = f"picard-rs-test-{os.getpid()}-{os.urandom(4).hex()}"
    result = run_array.run_captured([
        "docker", "run", "--rm", "--name", container, "--platform", run_array.PLATFORM,
        "-v", f"{picard}:{ROOT}:ro", "-v", f"{out_dir}:{OUT}", "-w", ROOT, run_array.IMAGE,
        f"java -jar $PICARD_JAR {shlex.quote(invocation['tool'])} {cli}",
    ], side="reference", container=container)
    return result, lambda text: text


def port(invocation, picard, out_dir, entry, port_dir):
    def host(value):
        return value.replace(ROOT, str(picard)).replace(OUT, str(out_dir))

    def back(text):
        return text.replace(str(picard), ROOT).replace(str(out_dir), OUT)

    env = dict(os.environ)
    env["PICARD_RS_PATH_MAP"] = f"{picard}={ROOT};{out_dir}={OUT}"
    argv = [str(Path(port_dir).resolve() / entry["port"])] + [host(t) for t in invocation["args"]]
    return run_array.run_captured(argv, env=env, cwd=picard), back


ANSWER = threading.Lock()

# The notice Picard logs before running a `NAME=value` command line (`CommandLineProgram.instanceMain`):
# the syntax's transition, not the tool's answer, and a refusal printed without a stack trace is
# read from the last block of standard error, which the notice can be.
LEGACY_NOTE = re.compile(r"(?m)^INFO\t[^\n]*\t\n(\n|\*{10}[^\n]*\n)+")


def side_answer(run, invocation, out_dir, entry):
    """(exit code, the answer compared), read as `run_array.py` reads a row's."""
    before = snapshot(out_dir)
    try:
        result, back = run()
    except run_array.RowTimeout as hung:
        return None, f"TIMEOUT: the {hung} ran past {run_array.ROW_TIMEOUT} s"
    error = run_array.first_error(LEGACY_NOTE.sub("", back(result.stderr or result.stdout)))
    on_stdout = "stderr" if entry.get("output") == "stderr" else entry.get("output") == "stdout"
    if on_stdout:
        text = run_array.stream_answer(result, on_stdout)
    else:
        # A tool given no output file prints its answer (ValidateSamFile's report) on stdout.
        text = written(out_dir, before, bool(entry.get("strip_program_records")), back)
        text += f"== <stdout>\n{result.stdout}\n" if result.stdout else ""
    with ANSWER:  # `canonical` reads a module-wide list of prefixes
        run_array.EXTRA_PREFIXES[:] = entry.get("strip_line_prefixes", [])
        answer = run_array.outcome(result.returncode, back(text), error, invocation["tool"],
                                   bool(entry.get("exit_code_is_a_result")))
    return result.returncode, answer


def replay_test(group, options, entries):
    """One test's invocations in the order it ran them, each side in one output directory of its
    own, so an invocation that reads what an earlier one wrote reads that side's file."""
    picard = Path(options.picard).resolve()
    results = []
    # shortcut: the reference's files are root's on a Linux host, so they are left for the runner
    # to discard rather than removed.
    with tempfile.TemporaryDirectory(prefix="picard-test-resources-", ignore_cleanup_errors=True) as tmp:
        ours, theirs = Path(tmp) / "port", Path(tmp) / "reference"
        ours.mkdir()
        theirs.mkdir()
        for invocation in group:
            entry = entries[invocation["tool"]]
            prepare(invocation, theirs)
            prepare(invocation, ours)
            code, reference = side_answer(lambda: oracle(invocation, picard, theirs, entry), invocation, theirs, entry)
            _, answer = side_answer(lambda: port(invocation, picard, ours, entry, options.port_dir), invocation, ours, entry)
            result = {"id": invocation["id"], "tool": invocation["tool"], "reference_exit": code,
                      "matched": answer == reference, "reference": reference[:4000], "port": answer[:4000]}
            if answer != reference:
                result["first_difference"] = first_difference(reference, answer)
            results.append(result)
    return results


def first_difference(reference, port):
    """Where two answers part, since a head kept to 4000 characters may not reach it."""
    left, right = reference.split("\n"), port.split("\n")
    for line, (a, b) in enumerate(zip(left, right)):
        if a != b:
            return {"line": line, "reference": a[:2000], "port": b[:2000]}
    return {"line": min(len(left), len(right)), "reference": f"{len(left)} lines", "port": f"{len(right)} lines"}


def chosen(options):
    """The tests this shard replays: every `count`-th, with each one's invocations of a ported
    tool that is not excluded. A test is never split, so its invocations run in order."""
    data = json.loads(INVOCATIONS.read_text())
    entries, excluded = coverage_entries(), exclusions()
    groups = {}
    for invocation in data["invocations"]:
        if invocation["tool"] in entries and invocation["id"] not in excluded \
                and (not options.tool or invocation["tool"] in options.tool):
            groups.setdefault(invocation.get("group") or test_of(invocation["id"]), []).append(invocation)
    index, count = (int(x) for x in options.shard.split("/"))
    return list(groups.values())[index::count]


def run(options):
    entries = coverage_entries()
    picked = chosen(options)
    results = []
    with concurrent.futures.ThreadPoolExecutor(options.jobs) as pool:
        for group in pool.map(lambda g: replay_test(g, options, entries), picked):
            for result in group:
                print(f"{result['id']}: {'match' if result['matched'] else 'DIFFERS'} (reference exit {result['reference_exit']})", flush=True)
            results += group
    Path(options.out).write_text(json.dumps({"shard": options.shard, "results": results}, indent=1) + "\n")
    print(f"replayed={len(results)} matched={sum(r['matched'] for r in results)}")
    # A divergence is a finding; an oracle that never started is a broken shard.
    return 1 if results and all(r["reference_exit"] in (None, 125, 126, 127) for r in results) else 0


def merge(options):
    """Fold the shards' results into `measured.json`, the file VALIDATION.md's column reads. Every
    ported tool has an entry, so a tool no test runs reads as such rather than as unmeasured."""
    data = json.loads(INVOCATIONS.read_text())
    entries, excluded = coverage_entries(), exclusions()
    results = [r for path in options.shards for r in json.loads(Path(path).read_text())["results"]]
    tools = {tool: {"extracted": 0, "replayed": 0, "matched": 0, "excluded": 0, "unextracted": 0, "diverged": []}
             for tool in entries}
    by_class = {}
    for invocation in data["invocations"]:
        by_class.setdefault(invocation["id"].split("#")[0], invocation["tool"])
        if invocation["tool"] in tools:
            tools[invocation["tool"]]["extracted"] += 1
            tools[invocation["tool"]]["excluded"] += invocation["id"] in excluded
    for missed in data["unextracted"]:
        if missed.get("partial"):
            continue  # it ran an invocation, counted above; what failed came after it
        name = missed["id"].split("#")[0]
        tool = by_class.get(name, name.removesuffix("Test"))
        if tool in tools:
            tools[tool]["unextracted"] += 1
    for result in results:
        found = tools[result["tool"]]
        found["replayed"] += 1
        found["matched"] += result["matched"]
        if not result["matched"]:
            found["diverged"].append(result["id"])
    for found in tools.values():
        found["diverged"].sort()
    MEASURED.write_text(json.dumps({"source": data["source"], "tools": dict(sorted(tools.items()))}, indent=1) + "\n")
    print(f"tools={len(tools)} replayed={len(results)} matched={sum(r['matched'] for r in results)}")
    return 0


COLUMN = "Reference tests"


def cell(found):
    if found is None:
        return "not measured"
    if not found["extracted"] and not found["unextracted"]:
        return "none in Picard's tests"
    text = f"{found['matched']} of {found['replayed']} match"
    if found["matched"] == found["replayed"] and found["replayed"]:
        text = f"**{text}**"
    notes = [f"{found[k]} {k}" for k in ("excluded", "unextracted") if found[k]]
    return text + (f" ({', '.join(notes)})" if notes else "")


def report(options):
    """The `Reference tests` column of docs/VALIDATION.md's status table, from `measured.json`."""
    measured = json.loads(MEASURED.read_text())["tools"] if MEASURED.is_file() else {}
    lines = VALIDATION.read_text().split("\n")
    for i, line in enumerate(lines):
        if not line.startswith("| "):
            continue
        cells = line.strip("|").split(" | ")
        cells = [c.strip() for c in cells]
        if cells[0] == "Tool":
            has = COLUMN in cells
            cells = cells if has else cells + [COLUMN]
        elif set(cells[0]) <= set("-"):
            cells = cells if has else cells + ["---"]
        elif cells[0].startswith("`"):
            value = cell(measured.get(cells[0].strip("`")))
            cells = cells[:-1] + [value] if has else cells + [value]
        lines[i] = "| " + " | ".join(cells) + " |"
    VALIDATION.write_text("\n".join(lines))
    return 0


def self_check():
    assert jdouble(0.5) == "0.5" and jdouble(1e-5) == "1.0E-5" and jdouble(3.0) == "3.0" and jdouble(1.5e7) == "1.5E7"
    assert jformat("O=%s x=%d%%", ["/a", 3]) == "O=/a x=3%"
    assert first_difference("a\nb", "a\nc") == {"line": 1, "reference": "b", "port": "c"}
    assert join_path("/a/b/", "c/") == "/a/b/c" and normalized("x//y/") == "x/y"
    assert RELATIVE_RESOURCE.sub(f"{ROOT}/testdata/", "I=testdata/picard/x.sam") == f"I={ROOT}/testdata/picard/x.sam"
    tool = Tool("SortSam")
    tool.fields.update(INPUT=["/a", "/b"], CREATE_INDEX=True)
    assert tool.tokens() == ["INPUT=/a", "INPUT=/b", "CREATE_INDEX=true"]
    assert test_of("A#b[2].0") == "A#b[2]"
    one = b"TIDX" + bytes(8) + b"/h/x.vcf\0" + (5).to_bytes(8, "little") + (7).to_bytes(8, "little") + b"md5\0"
    two = one.replace((7).to_bytes(8, "little"), (9).to_bytes(8, "little"))
    assert tribble_index(one, lambda t: t.replace("/h", "/work/out")) == tribble_index(two, lambda t: t.replace("/h", "/work/out"))
    note = "INFO\t2026-11-24 20:46:46\tSortSam\t\n\n********** NOTE: x\n**********\n\n\nERROR: boom\n"
    assert LEGACY_NOTE.sub("", note) == "ERROR: boom\n"
    assert cell({"extracted": 2, "replayed": 2, "matched": 2, "excluded": 0, "unextracted": 1}) == "**2 of 2 match** (1 unextracted)"


def main(argv):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    commands = parser.add_subparsers(dest="command", required=True)
    one = commands.add_parser("extract")
    one.add_argument("--picard", required=True, help="the pinned Picard checkout")
    one.add_argument("--constants", action="append", default=[],
                     help="another source tree to read constants from (htsjdk's src/main/java)")
    one.add_argument("--tag", default="3.5.0", help="the release the checkout is at, which CI clones")
    two = commands.add_parser("run")
    two.add_argument("--picard", required=True)
    two.add_argument("--port-dir", required=True, help="where every tool's port binary is, named as the manifest names it")
    two.add_argument("--shard", default="0/1")
    two.add_argument("--jobs", type=int, default=2)
    two.add_argument("--tool", action="append")
    two.add_argument("--out", default="test-resources-shard.json")
    three = commands.add_parser("merge")
    three.add_argument("shards", nargs="+")
    commands.add_parser("report")
    options = parser.parse_args(argv)
    self_check()
    return {"extract": extract, "run": run, "merge": merge, "report": report}[options.command](options)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
