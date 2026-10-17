import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import org.broadinstitute.barclay.argparser.LegacyCommandLineArgumentParser;
import org.broadinstitute.barclay.argparser.NamedArgumentDefinition;

/**
 * The names each of Picard's two parsers accepts, per tool.
 *
 * A `NAME=value` command line is parsed by `LegacyCommandLineArgumentParser`, which looks every
 * name up, upper-cased, in its private `optionMap`. A key of that map is a name the parser
 * accepts and nothing else is, so the map is dumped rather than the documentation: the hidden
 * `DROP_MISSING_CONTIGS` of `BedToIntervalList` is in it and in no inventory, and
 * `MarkDuplicatesWithMateCigar` has no `TAGGING_POLICY` although `MarkDuplicates` does.
 *
 * A command line with a `-` word is parsed by Barclay's own parser, which declares a different
 * set: `help`, `version` and `arguments_file`, and the flow arguments of a collection the legacy
 * parser never sees. Its full and short names are the `barclay` rows.
 *
 * Rows: `legacy <tool> <key> <long name>`, `barclay <tool> <name>`, and
 * `positional <tool> <min> <max>` for a tool that takes positional words.
 */
public class CommandLineDeclarationsDump {
    public static void main(final String[] args) throws Exception {
        System.out.println("# CommandLineDeclarationsDump: the names Picard's two parsers accept, per tool");
        for (final Map.Entry<String, Class<?>> entry : Declarations.tools().entrySet()) {
            final String tool = entry.getKey();
            final List<NamedArgumentDefinition> named = Declarations.barclay(entry.getValue()).getNamedArgumentDefinitions();
            final LegacyCommandLineArgumentParser parser = Declarations.legacy(entry.getValue());
            final Map<String, Object> map = (Map<String, Object>) Declarations.read(parser, "optionMap");
            for (final Map.Entry<String, Object> option : new TreeMap<>(map).entrySet()) {
                System.out.println("legacy\t" + tool + "\t" + option.getKey() + "\t" + Declarations.read(option.getValue(), "name"));
            }
            for (final NamedArgumentDefinition definition : named) {
                System.out.println("barclay\t" + tool + "\t" + definition.getFullName());
                if (!definition.getShortName().isEmpty()) {
                    System.out.println("barclay\t" + tool + "\t" + definition.getShortName());
                }
            }
            if (Declarations.takesPositional(parser)) {
                System.out.println("positional\t" + tool + "\t" + Declarations.read(parser, "minPositionalArguments")
                    + "\t" + Declarations.read(parser, "maxPositionalArguments"));
            }
        }
    }
}
