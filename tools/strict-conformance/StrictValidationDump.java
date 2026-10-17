/*
 * Oracle dump harness for the record validation every htsjdk 5.0.0 BAM reader runs.
 *
 * BAMFileReader numbers each record it decodes and, unless the reader is SILENT, calls
 * SAMRecord.isValid on it: STRICT throws the first error as a SAMFormatException, LENIENT prints
 * every one as "Ignoring SAM validation error: ...". Picard reads at its VALIDATION_STRINGENCY,
 * STRICT by default, so a malformed record ends any tool that reads it. Each case below is a BAM
 * of three records: a valid one, the malformed one, and a valid unmapped one, written by
 * SamFormatConverter at SILENT (or, for an index past the dictionary, by hand). The rows:
 *
 *   bam    <case>                         the input, the bytes
 *   run    <case>/<tool>/<stringency>     what the run printed: each "Ignoring SAM validation
 *                                         error" line, then "exception: <class>: <message>" for
 *                                         the throwable that escaped instanceMain, or "exit <n>"
 *
 *   java -cp picard.jar:. StrictValidationDump
 */
import htsjdk.samtools.BAMRecordCodec;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.TextCigarCodec;
import htsjdk.samtools.util.BinaryCodec;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;

public class StrictValidationDump {
    static final File DIR = new File("/work/strict");
    static final String HEADER = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n@RG\tID:rg1\tSM:s1\tLB:lib1\n";
    static final String[] STRINGENCIES = {"STRICT", "LENIENT", "SILENT"};

    /** The malformed record of each case, as a SAM line without its QNAME, which is "bad". */
    static final String[][] CASES = {
        {"mate-unmapped-flag", "65\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII"},
        {"unpaired-flags", "234\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII"},
        {"unpaired-mate-reference", "0\tchr1\t10\t60\t4M\tchr1\t5\t0\tACGT\tIIII"},
        {"not-first-or-second", "9\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII"},
        {"unaligned-mate-start", "73\tchr1\t10\t60\t4M\t*\t50\t0\tACGT\tIIII"},
        {"mate-start-past-end", "97\tchr1\t10\t60\t4M\tchr1\t5000\t0\tACGT\tIIII"},
        {"mate-reference-without-start", "65\tchr1\t10\t60\t4M\tchr1\t0\t0\tACGT\tIIII"},
        {"insert-size-out-of-range", "0\tchr1\t10\t60\t4M\t*\t0\t-2147483648\tACGT\tIIII"},
        {"unmapped-secondary", "260\tchr1\t10\t5\t*\t*\t0\t0\tACGT\tIIII"},
        {"unmapped-supplementary", "2052\tchr1\t10\t0\t*\t*\t0\t0\tACGT\tIIII"},
        {"mapped-without-cigar", "0\tchr1\t10\t60\t*\t*\t0\t0\tACGT\tIIII"},
        {"mapped-without-reference", "0\t*\t0\t60\t4M\t*\t0\t0\tACGT\tIIII"},
        {"read-group-not-found", "0\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:nope"},
        {"start-past-end", "0\tchr1\t2000\t60\t4M\t*\t0\t0\tACGT\tIIII"},
        {"maps-off-reference", "0\tchr1\t998\t60\t4M\t*\t0\t0\tACGT\tIIII"},
        {"cigar-sequence-mismatch", "0\tchr1\t10\t60\t5M\t*\t0\t0\tACGT\tIIII"},
        {"adjacent-indel", "0\tchr1\t10\t60\t2M1I1I\t*\t0\t0\tACGT\tIIII"},
        {"hard-clip-inside", "0\tchr1\t10\t60\t2M1H2M\t*\t0\t0\tACGT\tIIII"},
        {"soft-clip-inside", "0\tchr1\t10\t60\t1M1S2M\t*\t0\t0\tACGT\tIIII"},
        {"padding-at-end", "0\tchr1\t10\t60\t4M1P\t*\t0\t0\tACGT\tIIII"},
        {"no-real-operator", "0\tchr1\t10\t60\t4S\t*\t0\t0\tACGT\tIIII"},
        {"empty-sequence-bad-cigar", "0\tchr1\t10\t60\t2M1H2M\t*\t0\t0\t*\t*"},
    };

    public static void main(final String[] args) throws Exception {
        DIR.mkdirs();
        for (final String[] c : CASES) {
            final File sam = new File(DIR, c[0] + ".sam");
            // Every record names the read group, because MarkDuplicates keys a mapped pair on
            // getReadGroup().getReadGroupId() and throws a NullPointerException without one.
            final String bad = c[1].contains("RG:Z:") ? c[1] : c[1] + "\tRG:Z:rg1";
            Files.write(sam.toPath(), (HEADER
                    + "first\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:rg1\n"
                    + "bad\t" + bad + "\n"
                    + "after\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\tRG:Z:rg1\n").getBytes(StandardCharsets.UTF_8));
            final File bam = new File(DIR, c[0] + ".bam");
            final int rc = new picard.sam.SamFormatConverter().instanceMain(new String[] {
                "INPUT=" + sam, "OUTPUT=" + bam, "VALIDATION_STRINGENCY=SILENT"});
            if (rc != 0) throw new IllegalStateException(c[0] + ": SamFormatConverter exited " + rc);
        }
        // A reference index past the dictionary cannot be written through a SAMFileWriter, which
        // resolves every record against its header, so these two are encoded by hand.
        handWritten("reference-index-past-dictionary", "chr2", "*");
        handWritten("mate-index-past-dictionary", "chr1", "chr2");

        final String[] names = new String[CASES.length + 2];
        for (int i = 0; i < CASES.length; i++) names[i] = CASES[i][0];
        names[CASES.length] = "reference-index-past-dictionary";
        names[CASES.length + 1] = "mate-index-past-dictionary";
        for (final String name : names) {
            final File bam = new File(DIR, name + ".bam");
            emit("bam", name, hex(bam));
            final String in = "INPUT=" + bam;
            for (final String s : STRINGENCIES) {
                final String vs = "VALIDATION_STRINGENCY=" + s;
                run(name, "CollectQualityYieldMetrics", s, new picard.analysis.CollectQualityYieldMetrics(),
                        in, "OUTPUT=/work/strict/out.cqym.txt", vs);
                run(name, "SortSam", s, new picard.sam.SortSam(),
                        in, "OUTPUT=/work/strict/out.sorted.bam", "SORT_ORDER=queryname", vs);
                run(name, "MarkDuplicates", s, new picard.sam.markduplicates.MarkDuplicates(),
                        in, "OUTPUT=/work/strict/out.md.bam", "METRICS_FILE=/work/strict/out.md.txt", vs);
                run(name, "CleanSam", s, new picard.sam.CleanSam(),
                        in, "OUTPUT=/work/strict/out.clean.bam", vs);
            }
        }
    }

    /** A BAM of chr1 alone whose middle record names `ref` and `mateRef` of a two-contig header. */
    static void handWritten(final String name, final String ref, final String mateRef) throws Exception {
        final SAMFileHeader one = new SAMFileHeader();
        one.addSequence(new SAMSequenceRecord("chr1", 1000));
        final SAMFileHeader two = new SAMFileHeader();
        two.addSequence(new SAMSequenceRecord("chr1", 1000));
        two.addSequence(new SAMSequenceRecord("chr2", 1000));
        final String text = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n";
        try (BlockCompressedOutputStream out = new BlockCompressedOutputStream(
                new FileOutputStream(new File(DIR, name + ".bam")), (File) null)) {
            final BinaryCodec codec = new BinaryCodec(out);
            codec.writeBytes("BAM\1".getBytes(StandardCharsets.US_ASCII));
            final byte[] t = text.getBytes(StandardCharsets.US_ASCII);
            codec.writeInt(t.length);
            codec.writeBytes(t);
            codec.writeInt(1);
            codec.writeInt(5);
            codec.writeBytes("chr1\0".getBytes(StandardCharsets.US_ASCII));
            codec.writeInt(1000);
            final BAMRecordCodec records = new BAMRecordCodec(two);
            records.setOutputStream(out);
            records.encode(record(two, "first", 0, "chr1", 1, "*", 0));
            records.encode(record(two, "bad", "*".equals(mateRef) ? 0 : 97, ref, 10, mateRef, "*".equals(mateRef) ? 0 : 20));
            records.encode(record(two, "after", 4, "*", 0, "*", 0));
        }
    }

    static SAMRecord record(final SAMFileHeader header, final String name, final int flags, final String ref,
                            final int start, final String mateRef, final int mateStart) {
        final SAMRecord r = new SAMRecord(header);
        r.setReadName(name);
        r.setFlags(flags);
        r.setReferenceName(ref);
        r.setAlignmentStart(start);
        r.setMappingQuality((flags & 4) != 0 ? 0 : 60);
        r.setCigar(TextCigarCodec.decode((flags & 4) != 0 ? "*" : "4M"));
        r.setMateReferenceName(mateRef);
        r.setMateAlignmentStart(mateStart);
        r.setReadBases("ACGT".getBytes(StandardCharsets.US_ASCII));
        r.setBaseQualities(new byte[] {40, 40, 40, 40});
        return r;
    }

    /** One tool run, with what reached stderr through System.err and what escaped instanceMain. */
    static void run(final String name, final String tool, final String stringency,
                    final picard.cmdline.CommandLineProgram program, final String... args) {
        final PrintStream err = System.err;
        final ByteArrayOutputStream captured = new ByteArrayOutputStream();
        System.setErr(new PrintStream(captured, true));
        String outcome;
        try {
            final int rc = program.instanceMain(args);
            outcome = rc == 0 ? "ok" : "exit " + rc;
        } catch (final Throwable t) {
            outcome = "exception: " + t.getClass().getName() + ": " + t.getMessage();
        } finally {
            System.setErr(err);
        }
        final StringBuilder payload = new StringBuilder();
        for (final String line : captured.toString().split("\n")) {
            if (line.startsWith("Ignoring SAM validation error")) payload.append(line).append('\n');
        }
        payload.append(outcome);
        emit("run", name + "/" + tool + "/" + stringency, esc(payload.toString()));
    }

    static String hex(final File f) throws Exception {
        final byte[] b = Files.readAllBytes(f.toPath());
        final StringBuilder sb = new StringBuilder(b.length * 2);
        for (final byte x : b) sb.append(String.format("%02x", x));
        return sb.toString();
    }

    static String esc(final String p) {
        return p.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n");
    }

    static void emit(final String kind, final String kase, final String payload) {
        System.out.println(kind + "\t" + kase + "\t" + payload);
    }
}
