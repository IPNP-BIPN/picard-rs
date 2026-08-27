/*
 * CheckDuplicateMarking's verdict, taken from the reference.
 *
 * The tool walks a queryname-sorted file and complains whenever two records of the same name
 * disagree about their duplicate flag. What is measured is which records are compared at all,
 * which disagreement is reported, and what the exit code and the bad-name file carry.
 *
 * Eleven behaviours this is built to catch.
 *
 *   - THE EXIT CODE IS THE COUNT'S SIGN: zero when every name agrees and one when any does not,
 *     whatever the count is;
 *   - THE COMPARISON IS AGAINST THE FIRST RECORD OF EACH NAME, not against the record before, so
 *     a name whose flags go true, false, false reports TWO bad records and not one;
 *   - A NAME THAT APPEARS ONCE IS NEVER BAD;
 *   - --OUTPUT HOLDS ONE LINE PER BAD RECORD, so a name that disagrees twice is written twice;
 *   - THE FILE IS SORTED BY QUERYNAME FIRST when it is not already, so records of one name that
 *     are far apart in a coordinate-sorted file are still compared;
 *   - --MODE ALL COMPARES EVERY RECORD, secondary and supplementary ones included;
 *   - --MODE PRIMARY_ONLY SKIPS THE SECONDARY AND SUPPLEMENTARY ONES, so a disagreement carried
 *     only by them disappears;
 *   - --MODE PRIMARY_MAPPED_ONLY SKIPS THE UNMAPPED ONES TOO;
 *   - --MODE PRIMARY_PROPER_PAIR_ONLY SKIPS EVERYTHING THAT IS NOT A PROPER PAIR, which includes
 *     every unpaired read;
 *   - THE MODE'S SKIP HAPPENS BEFORE THE TALLY, so a skipped record is neither compared nor
 *     remembered: the same file is clean under PRIMARY_ONLY and bad under ALL;
 *   - AND --OUTPUT IS OPTIONAL, the verdict being the same without it.
 *
 * Output:
 *
 *     sam\t<case>\t<that bam as sam, without its header, escaped>
 *     verdict\t<case>\t<exit code>
 *     bad\t<case>\t<the OUTPUT file, escaped>
 */

import htsjdk.samtools.*;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.zip.DeflaterFactory;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class CheckDuplicateMarkingDump {

    static final StringBuilder buf = new StringBuilder();

    static void emit(final String kind, final String name, final String payload) {
        buf.append(kind).append('\t').append(name).append('\t')
                .append(payload.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n"))
                .append('\n');
    }

    /** name, start, duplicate, secondary, supplementary, unmapped, properPair. */
    record Rec(String name, int start, boolean duplicate, boolean secondary,
               boolean supplementary, boolean unmapped, boolean properPair) { }

    static Rec plain(final String name, final int start, final boolean duplicate) {
        return new Rec(name, start, duplicate, false, false, false, true);
    }

    static SAMRecord record(final SAMFileHeader header, final Rec spec) {
        final SAMRecord read = new SAMRecord(header);
        read.setReadName(spec.name());
        read.setReadPairedFlag(true);
        read.setFirstOfPairFlag(true);
        read.setReadBases("ACGTACGTAC".getBytes());
        final byte[] quality = new byte[10];
        Arrays.fill(quality, (byte) 30);
        read.setBaseQualities(quality);
        if (spec.unmapped()) {
            read.setReadUnmappedFlag(true);
            read.setReferenceIndex(SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX);
            read.setMateUnmappedFlag(true);
            read.setMateReferenceIndex(SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX);
        } else {
            read.setReferenceIndex(0);
            read.setAlignmentStart(spec.start());
            read.setCigarString("10M");
            read.setMappingQuality(60);
            read.setMateReferenceIndex(0);
            read.setMateAlignmentStart(spec.start() + 100);
            read.setProperPairFlag(spec.properPair());
        }
        read.setDuplicateReadFlag(spec.duplicate());
        read.setSecondaryAlignment(spec.secondary());
        read.setSupplementaryAlignmentFlag(spec.supplementary());
        return read;
    }

    static void run(final String name, final List<Rec> specs,
                    final SAMFileHeader.SortOrder order, final String... extra)
            throws Exception {
        final Path dir = Files.createTempDirectory("cdm");
        final SAMFileHeader header = new SAMFileHeader();
        final SAMSequenceDictionary dictionary = new SAMSequenceDictionary();
        dictionary.addSequence(new SAMSequenceRecord("chr1", 100000));
        header.setSequenceDictionary(dictionary);
        header.setSortOrder(order);
        final File bam = new File(dir.toFile(), "in.bam");
        try (final SAMFileWriter writer = new SAMFileWriterFactory()
                .setUseAsyncIo(false)
                .makeBAMWriter(header, false, bam)) {
            for (final Rec spec : specs) {
                writer.addAlignment(record(header, spec));
            }
        }
        final StringBuilder sam = new StringBuilder();
        try (final SamReader reader = SamReaderFactory.makeDefault().open(bam)) {
            for (final SAMRecord read : reader) {
                sam.append(read.getSAMString());
            }
        }
        emit("sam", name, sam.toString());

        final File out = new File(dir.toFile(), "bad.txt");
        final List<String> argv = new ArrayList<>(List.of("I=" + bam.getAbsolutePath()));
        final List<String> extras = Arrays.asList(extra);
        if (!extras.contains("NO_OUTPUT")) {
            argv.add("O=" + out.getAbsolutePath());
        }
        for (final String value : extras) {
            if (!value.equals("NO_OUTPUT")) {
                argv.add(value);
            }
        }
        final int code = new picard.sam.markduplicates.CheckDuplicateMarking()
                .instanceMain(argv.toArray(new String[0]));
        emit("verdict", name, Integer.toString(code));
        if (out.exists()) {
            emit("bad", name, Files.readString(out.toPath()));
        } else {
            emit("bad", name, "");
        }
    }

    public static void main(final String[] args) throws Exception {
        BlockCompressedOutputStream.setDefaultDeflaterFactory(new DeflaterFactory());

        // Every name agrees.
        run("all-agree", List.of(
                plain("a", 100, true), plain("a", 200, true),
                plain("b", 300, false), plain("b", 400, false)),
                SAMFileHeader.SortOrder.queryname);

        // One name disagrees once.
        run("one-disagreement", List.of(
                plain("a", 100, true), plain("a", 200, false),
                plain("b", 300, false), plain("b", 400, false)),
                SAMFileHeader.SortOrder.queryname);

        // A name whose flags go true, false, false: the comparison is against the FIRST record,
        // so both of the later ones are bad.
        run("two-disagreements-one-name", List.of(
                plain("a", 100, true), plain("a", 200, false), plain("a", 300, false)),
                SAMFileHeader.SortOrder.queryname);

        // A name that appears once.
        run("single-record", List.of(plain("a", 100, true), plain("b", 200, false)),
                SAMFileHeader.SortOrder.queryname);

        // A coordinate-sorted file whose records of one name are far apart.
        run("coordinate-sorted", List.of(
                plain("a", 100, true), plain("b", 200, false),
                plain("a", 900, false), plain("b", 1000, false)),
                SAMFileHeader.SortOrder.coordinate);

        // The disagreement carried only by a secondary record.
        final List<Rec> secondary = List.of(
                plain("a", 100, true),
                new Rec("a", 200, false, true, false, false, true));
        run("secondary-disagrees-all", secondary, SAMFileHeader.SortOrder.queryname);
        run("secondary-disagrees-primary-only", secondary, SAMFileHeader.SortOrder.queryname,
                "MODE=PRIMARY_ONLY");

        // The disagreement carried only by a supplementary record.
        final List<Rec> supplementary = List.of(
                plain("a", 100, true),
                new Rec("a", 200, false, false, true, false, true));
        run("supplementary-disagrees-all", supplementary, SAMFileHeader.SortOrder.queryname);
        run("supplementary-disagrees-primary-only", supplementary,
                SAMFileHeader.SortOrder.queryname, "MODE=PRIMARY_ONLY");

        // The disagreement carried only by an unmapped record.
        final List<Rec> unmapped = List.of(
                plain("a", 100, true),
                new Rec("a", 0, false, false, false, true, false));
        run("unmapped-disagrees-primary-only", unmapped, SAMFileHeader.SortOrder.queryname,
                "MODE=PRIMARY_ONLY");
        run("unmapped-disagrees-mapped-only", unmapped, SAMFileHeader.SortOrder.queryname,
                "MODE=PRIMARY_MAPPED_ONLY");

        // The disagreement carried only by a record that is not a proper pair.
        final List<Rec> improper = List.of(
                plain("a", 100, true),
                new Rec("a", 200, false, false, false, false, false));
        run("improper-disagrees-mapped-only", improper, SAMFileHeader.SortOrder.queryname,
                "MODE=PRIMARY_MAPPED_ONLY");
        run("improper-disagrees-proper-only", improper, SAMFileHeader.SortOrder.queryname,
                "MODE=PRIMARY_PROPER_PAIR_ONLY");

        // A secondary record that disagrees with the name's primaries. The writer sorts by
        // query name, which puts it LAST, so this shows the skip and not the first-record rule.
        run("skipped-first-record", List.of(
                new Rec("a", 100, true, true, false, false, true),
                plain("a", 200, false), plain("a", 300, false)),
                SAMFileHeader.SortOrder.queryname, "MODE=PRIMARY_ONLY");

        // No output file at all.
        run("no-output-file", List.of(
                plain("a", 100, true), plain("a", 200, false)),
                SAMFileHeader.SortOrder.queryname, "NO_OUTPUT");

        System.out.print(buf);
    }
}
