/*
 * Oracle dump harness for CollectRnaSeqMetrics conformance in picard-rs.
 *
 * Emits the whole conformance corpus as an escaped TSV to stdout, one row per input or output:
 *
 *   sam       <case>  the exact SAM the oracle consumed
 *   refflat   <case>  the gene model
 *   ribosomal <case>  the ribosomal interval_list (only for cases that use one)
 *   metrics   <case>  the .rna_metrics file the oracle produced
 *
 * The committed corpus is `java ... RnaSeqMetricsDump | gzip > tests/data/rnaseq_metrics.txt.gz`,
 * and CI regenerates it in the pinned amd64 oracle and compares the metrics body (the command-line
 * and timestamp header lines are stripped, being JVM-chosen paths and a clock read). The Rust test
 * parses the sam/refflat/ribosomal rows as its inputs, so those are regenerated too.
 *
 * Two cases:
 *   basic     the reference's own CollectRnaSeqMetricsTest.testBasic. Its single 451bp transcript
 *             is below MINIMUM_LENGTH=500, so the histogram is empty and the MEDIAN_* metrics are 0;
 *             this pins the whole metrics row.
 *   coverage  three genes on chr1, each a single-exon 600bp coding transcript with deliberately
 *             different coverage depth, so the normalized_coverage histogram is non-empty and its
 *             floating-point fold over the transcripts is order-sensitive.
 *
 * Run under the pinned en_US locale so metrics format with '.' not ',' (htsjdk-rs decision 0011):
 *   java -Duser.language=en -Duser.country=US -cp picard-fat.jar:. RnaSeqMetricsDump
 */
import htsjdk.samtools.*;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import picard.analysis.CollectRnaSeqMetrics;

import java.io.File;
import java.io.PrintStream;
import java.nio.file.Files;

public class RnaSeqMetricsDump {
    public static void main(final String[] args) throws Exception {
        emitBasic();
        emitCoverage();
    }

    /** CollectRnaSeqMetricsTest.testBasic, verbatim. */
    private static void emitBasic() throws Exception {
        final String sequence = "chr1";
        final String ignoredSequence = "chrM";

        final SAMRecordSetBuilder builder = new SAMRecordSetBuilder(true, SAMFileHeader.SortOrder.coordinate);
        builder.setRandomSeed(0);
        final int idx = builder.getHeader().getSequenceIndex(sequence);
        builder.addPair("pair1", idx, 45, 475);
        builder.addPair("pair2", idx, 90, 225);
        builder.addPair("pair3", idx, 120, 600);
        builder.addFrag("frag1", idx, 150, true);
        builder.addFrag("frag2", idx, 450, true);
        builder.addFrag("frag3", idx, 225, false);
        builder.addPair("rrnaPair", idx, 400, 500);
        builder.addFrag("ignoredFrag", builder.getHeader().getSequenceIndex(ignoredSequence), 1, false);

        final File samFile = writeSam(builder);

        final IntervalList rRna = new IntervalList(builder.getHeader());
        rRna.add(new Interval(sequence, 300, 520, true, "rRNA"));
        final File rRnaFile = File.createTempFile("rnaseq-basic-rRNA-", ".interval_list");
        rRnaFile.deleteOnExit();
        rRna.write(rRnaFile);

        final File refFlat = writeRefFlat("myGene\tmyTranscript\tchr1\t+\t49\t500\t74\t400\t2\t49,249\t200,500");

        final File metrics = run(new String[] {
                "INPUT=" + samFile.getAbsolutePath(),
                "REF_FLAT=" + refFlat.getAbsolutePath(),
                "RIBOSOMAL_INTERVALS=" + rRnaFile.getAbsolutePath(),
                "STRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND",
                "IGNORE_SEQUENCE=" + ignoredSequence
        });

        emit("sam", "basic", samFile);
        emit("refflat", "basic", refFlat);
        emit("ribosomal", "basic", rRnaFile);
        emit("metrics", "basic", metrics);
    }

    /** Three long qualifying transcripts of different depth: exercises the histogram and MEDIAN_* metrics. */
    private static void emitCoverage() throws Exception {
        final String sequence = "chr1";

        final SAMRecordSetBuilder builder = new SAMRecordSetBuilder(true, SAMFileHeader.SortOrder.coordinate);
        builder.setRandomSeed(0);
        final int idx = builder.getHeader().getSequenceIndex(sequence);

        final int[] geneStarts = {1000, 3000, 5000};
        final int[] fragCounts = {45, 30, 22};
        final int[] fragSteps  = {12, 18, 26};
        int frag = 0;
        for (int g = 0; g < geneStarts.length; ++g) {
            for (int k = 0; k < fragCounts[g]; ++k) {
                final int start = geneStarts[g] + 1 + k * fragSteps[g];
                builder.addFrag("f" + (frag++), idx, start, false);
            }
        }

        final File samFile = writeSam(builder);
        final File refFlat = writeRefFlat(
                "covGeneA\tcovTxA\tchr1\t+\t1000\t1600\t1000\t1600\t1\t1000,\t1600,\n"
              + "covGeneB\tcovTxB\tchr1\t+\t3000\t3600\t3000\t3600\t1\t3000,\t3600,\n"
              + "covGeneC\tcovTxC\tchr1\t+\t5000\t5600\t5000\t5600\t1\t5000,\t5600,");

        final File metrics = run(new String[] {
                "INPUT=" + samFile.getAbsolutePath(),
                "REF_FLAT=" + refFlat.getAbsolutePath(),
                "STRAND_SPECIFICITY=NONE"
        });

        emit("sam", "coverage", samFile);
        emit("refflat", "coverage", refFlat);
        emit("metrics", "coverage", metrics);
    }

    private static File writeSam(final SAMRecordSetBuilder builder) throws Exception {
        final File samFile = File.createTempFile("rnaseq-", ".sam");
        samFile.deleteOnExit();
        final SAMFileWriter w = new SAMFileWriterFactory().makeSAMWriter(builder.getHeader(), false, samFile);
        for (final SAMRecord rec : builder.getRecords()) w.addAlignment(rec);
        w.close();
        return samFile;
    }

    private static File writeRefFlat(final String contents) throws Exception {
        final File f = File.createTempFile("rnaseq-", ".refFlat");
        f.deleteOnExit();
        final PrintStream ps = new PrintStream(f);
        ps.println(contents);
        ps.close();
        return f;
    }

    private static File run(final String[] toolArgs) throws Exception {
        final File metrics = File.createTempFile("rnaseq-", ".rna_metrics");
        metrics.deleteOnExit();
        final String[] full = new String[toolArgs.length + 1];
        full[0] = "OUTPUT=" + metrics.getAbsolutePath();
        System.arraycopy(toolArgs, 0, full, 1, toolArgs.length);
        final int rc = new CollectRnaSeqMetrics().instanceMain(full);
        if (rc != 0) {
            System.err.println("CollectRnaSeqMetrics exited " + rc);
            System.exit(rc);
        }
        return metrics;
    }

    private static void emit(final String kind, final String kase, final File file) throws Exception {
        final String payload = new String(Files.readAllBytes(file.toPath()))
                .replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n");
        System.out.println(kind + "\t" + kase + "\t" + payload);
    }
}
