/*
 * `FifoBuffer`, taken from the reference.
 *
 * The smallest tool in Picard: it copies its input to its output through a circular byte buffer,
 * so a pipeline can decouple a slow writer from a slow reader. What there is to measure is that
 * the bytes come out unchanged whatever the buffer's shape, and what it does when the shape is
 * impossible.
 *
 * Six behaviours this is built to catch.
 *
 *   - THE BYTES ARE UNCHANGED, whatever their length and whatever the buffer's;
 *   - A BUFFER SMALLER THAN THE INPUT still copies all of it, because the buffer is CIRCULAR and
 *     the two threads take turns;
 *   - AN IO SIZE LARGER THAN THE BUFFER is a refusal rather than a clamp;
 *   - AN EMPTY INPUT IS AN EMPTY OUTPUT and a zero status, not a hang;
 *   - BINARY INPUT SURVIVES, including the bytes a text tool would mangle: a NUL, a lone carriage
 *     return and a byte above 127;
 *   - AND THE TOOL IS QUIET BY DEFAULT when constructed on the streams, which is why a dump can
 *     read its output at all.
 *
 * Output:
 *
 *     out\t<case>\t<the output's bytes, base64>
 *     code\t<case>\t<the exit status>
 *     error\t<case>\t<exception class>:<message>
 *
 * Usage: FifoBufferDump
 */

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Base64;
import java.util.List;

public class FifoBufferDump {

    static final StringBuilder buf = new StringBuilder();

    static void emit(final String kind, final String name, final String payload) {
        buf.append(kind).append('\t').append(name).append('\t')
                .append(payload.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n"))
                .append('\n');
    }

    static void run(final String name, final byte[] input, final String... extra) {
        final ByteArrayOutputStream captured = new ByteArrayOutputStream();
        final picard.util.FifoBuffer tool = new picard.util.FifoBuffer(
                new ByteArrayInputStream(input),
                new PrintStream(captured, true, StandardCharsets.UTF_8));
        final List<String> argv = new ArrayList<>(Arrays.asList(extra));
        final PrintStream original = System.out;
        final PrintStream originalError = System.err;
        final ByteArrayOutputStream said = new ByteArrayOutputStream();
        try {
            System.setOut(new PrintStream(said, true, StandardCharsets.UTF_8));
            System.setErr(new PrintStream(said, true, StandardCharsets.UTF_8));
            final Object code = tool.instanceMain(argv.toArray(new String[0]));
            System.setOut(original);
            System.setErr(originalError);
            emit("code", name, String.valueOf(code));
            emit("out", name, Base64.getEncoder().encodeToString(captured.toByteArray()));
        } catch (final Exception | AssertionError e) {
            System.setOut(original);
            System.setErr(originalError);
            Throwable cause = e;
            while (cause.getCause() != null && cause.getCause() != cause) {
                cause = cause.getCause();
            }
            emit("error", name, cause.getClass().getName() + ":" + String.valueOf(cause.getMessage()));
        } finally {
            System.setOut(original);
            System.setErr(originalError);
        }
    }

    public static void main(final String[] args) {
        final byte[] text = "the quick brown fox\n".getBytes(StandardCharsets.UTF_8);
        run("a-line-of-text", text);
        run("nothing-at-all", new byte[0]);

        // A buffer smaller than the input, and an IO size larger than the buffer.
        run("a-small-buffer", text, "BUFFER_SIZE=8");
        run("a-buffer-of-one-byte", text, "BUFFER_SIZE=1");
        run("an-io-size-above-the-buffer", text, "BUFFER_SIZE=8", "IO_SIZE=64");

        // Bytes a text tool would mangle.
        final byte[] binary = new byte[]{0, 13, (byte) 200, 10, 9, (byte) 255};
        run("bytes-that-are-not-text", binary);

        // Something longer than any single read.
        final byte[] long_ = new byte[100_000];
        for (int index = 0; index < long_.length; index++) {
            long_[index] = (byte) (index % 251);
        }
        run("a-hundred-thousand-bytes", long_);
        run("a-hundred-thousand-bytes-through-a-small-buffer", long_, "BUFFER_SIZE=1024");

        System.out.print(buf);
    }
}
