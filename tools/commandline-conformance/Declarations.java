import java.lang.reflect.Field;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import org.broadinstitute.barclay.argparser.CommandLineArgumentParser;
import org.broadinstitute.barclay.argparser.LegacyCommandLineArgumentParser;
import org.broadinstitute.barclay.argparser.NamedArgumentDefinition;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.PicardCommandLine;

/** Shared by the two dumps: every tool, and the two parsers Picard builds for one. */
final class Declarations {
    private Declarations() {}

    /** Every command-line program in `picard`, by simple name, as `PicardCommandLine` lists them. */
    static Map<String, Class<?>> tools() {
        final Map<String, Class<?>> tools = new TreeMap<>();
        PicardCommandLine.processAllCommandLinePrograms(
            Collections.singletonList("picard"), (clazz, properties) -> tools.put(clazz.getSimpleName(), clazz));
        return tools;
    }

    /**
     * The parser a `NAME=value` command line gets. `getCommandLineParserForArgs` is the method
     * `instanceMain` goes through: it drops the Barclay special arguments and builds a
     * `LegacyCommandLineArgumentParser` on the tool, so asking it is asking the real thing.
     */
    static LegacyCommandLineArgumentParser legacy(final Class<?> clazz) throws Exception {
        final CommandLineProgram tool = (CommandLineProgram) clazz.getDeclaredConstructor().newInstance();
        return (LegacyCommandLineArgumentParser) tool.getCommandLineParserForArgs(new String[] {"X=1"});
    }

    /**
     * The parser a command line with a `-` word gets. `useLegacyParser` is a static that the
     * legacy call sets for every tool built after it, and a tool reads it when it is constructed,
     * so it is put back first.
     */
    static CommandLineArgumentParser barclay(final Class<?> clazz) throws Exception {
        final Field legacy = CommandLineProgram.class.getDeclaredField("useLegacyParser");
        legacy.setAccessible(true);
        legacy.set(null, Boolean.FALSE);
        final CommandLineProgram tool = (CommandLineProgram) clazz.getDeclaredConstructor().newInstance();
        return (CommandLineArgumentParser) tool.getCommandLineParser();
    }

    /** A private field of the legacy parser or of one of its definitions; nothing public reads them. */
    static Object read(final Object owner, final String name) throws Exception {
        final Field field = owner.getClass().getDeclaredField(name);
        field.setAccessible(true);
        return field.get(owner);
    }

    static boolean takesPositional(final LegacyCommandLineArgumentParser parser) throws Exception {
        return read(parser, "positionalArguments") != null;
    }
}
