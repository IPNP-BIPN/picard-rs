/*
 * Builds the fixture corpus the covering arrays run against.
 *
 * Usage: MakeFixtures <output directory>
 *
 * A covering array cannot invent a file path: the value has to be a file that exists and holds
 * content the tool accepts (gatk-rs tools/coverage/domains.py excludes every path-typed argument
 * for exactly that reason, which is most of what it excludes). This produces that corpus, small
 * and deterministic, at fixed paths, so a row of the array can be turned into a command line.
 *
 * Three properties matter and are the reason this is a program rather than a directory of
 * committed files:
 *
 *   1. Deterministic. Fixed seed, fixed content, no timestamps, no temp directories. Two runs
 *      produce the same bytes, so a divergence between the oracle and the port is about the tool
 *      and not about the input.
 *   2. Small. Every row of a covering array runs the tool once; HaplotypeCaller's t=2 array is 62
 *      rows and the whole inventory is 19,437. The corpus is sized for that, not for realism.
 *   3. Branchy. Uniform perfect reads exercise one path. These reads carry unmapped mates,
 *      duplicates, secondary and supplementary alignments, soft clips, an indel, no-calls, both
 *      strands, two read groups and two libraries, so a row that flips a filtering argument
 *      actually changes the answer. A corpus where every argument produces the same output would
 *      make a covering array look green while testing nothing.
 *
 * The reference is two short contigs, which is enough for a sequence dictionary, an interval
 * list, and a tool that needs REFERENCE_SEQUENCE.
 */

import htsjdk.samtools.*;
import htsjdk.samtools.util.SequenceUtil;
import java.io.File;
import java.io.PrintWriter;
import java.util.Random;

public class MakeFixtures {

    static final int CHR1 = 2_000;
    static final int CHR2 = 1_000;
    static final int READS = 400;
    static final int READ_LENGTH = 50;

    public static void main(String[] args) throws Exception {
        File dir = new File(args.length > 0 ? args[0] : "fixtures");
        dir.mkdirs();

        String chr1 = reference(CHR1, 20260729L);
        String chr2 = reference(CHR2, 20260730L);

        writeFasta(new File(dir, "ref.fasta"), chr1, chr2);
        writeFai(new File(dir, "ref.fasta.fai"), chr1, chr2);

        SAMFileHeader header = header(SAMFileHeader.SortOrder.coordinate);
        writeBam(new File(dir, "small.bam"), header, reads(header, true), true);
        writeSam(new File(dir, "small.sam"), header, reads(header, true));

        SAMFileHeader queryname = header(SAMFileHeader.SortOrder.queryname);
        writeBam(new File(dir, "queryname.bam"), queryname, reads(queryname, false), false);

        SAMFileHeader oneGroup = header(SAMFileHeader.SortOrder.coordinate);
        oneGroup.setReadGroups(java.util.Collections.singletonList(oneGroup.getReadGroup("rg1")));
        java.util.List<SAMRecord> oneGroupReads = reads(oneGroup, true);
        for (SAMRecord r : oneGroupReads) r.setAttribute("RG", "rg1");
        writeBam(new File(dir, "one_read_group.bam"), oneGroup, oneGroupReads, false);

        // Illumina-style read names, for the tools that read a physical location out of one.
        // `PositionBasedDownsampleSam` keeps reads by where they sit on the flowcell, and the
        // names everywhere else in this corpus (`read0314`) carry no location at all: the parser
        // refuses them, every read lands on the same defaulted tile, and the mask becomes
        // all-or-nothing. These names give it tiles and coordinates to work on, in a file of its
        // own so that no existing fixture's bytes move.
        SAMFileHeader tiled = header(SAMFileHeader.SortOrder.coordinate);
        java.util.List<SAMRecord> tiledReads = reads(tiled, true);
        for (SAMRecord r : tiledReads) {
            // The new name is derived from the old one, so the two ends of a pair still share it:
            // `read0314` gives 314, and both mates carry that number.
            int n = Integer.parseInt(r.getReadName().replaceAll("[^0-9]", ""));
            int tile = 1101 + (n % 4);
            int x = 1000 + ((n * 137) % 20000);
            int y = 1000 + ((n * 251) % 20000);
            r.setReadName(String.format("INST:1:FLOWCELL:1:%d:%d:%d", tile, x, y));
        }
        writeBam(new File(dir, "tiled.bam"), tiled, tiledReads, false);

        // Reads carrying `MC`, for the duplicate markers that read the mate's cigar instead of
        // waiting for the mate. Without the tag `SimpleMarkDuplicatesWithMateCigar` refuses every
        // file outright and `MarkDuplicatesWithMateCigar` skips every pair, so a corpus without it
        // measures the refusal and nothing else. `SamPairUtil.setMateInformation` writes the tag
        // the same way `FixMateInformation --ADD_MATE_CIGAR` does; the file is new, so no existing
        // fixture's bytes move.
        SAMFileHeader mateCigarHeader = header(SAMFileHeader.SortOrder.coordinate);
        java.util.List<SAMRecord> mateCigarReads = reads(mateCigarHeader, true);
        java.util.Map<String, java.util.List<SAMRecord>> byName = new java.util.LinkedHashMap<>();
        for (SAMRecord r : mateCigarReads) {
            byName.computeIfAbsent(r.getReadName(), k -> new java.util.ArrayList<>()).add(r);
        }
        for (java.util.List<SAMRecord> template : byName.values()) {
            if (template.size() == 2) {
                htsjdk.samtools.SamPairUtil.setMateInfo(template.get(0), template.get(1), true);
            }
        }
        writeBam(new File(dir, "mate_cigar.bam"), mateCigarHeader, mateCigarReads, true);

        // Records of one query name that DISAGREE about their duplicate flag, for
        // `CheckDuplicateMarking`. Everywhere else in this corpus both ends of a pair carry the
        // same flag, so the tool finds nothing whatever `MODE` is asked for: the array covers the
        // argument and observes none of it.
        //
        // Which record of the pair is flipped decides which modes see the disagreement, and that
        // is what makes the four values four answers. A flipped SECONDARY or SUPPLEMENTARY record
        // is seen by `ALL` alone; a flipped UNMAPPED one is also seen by `PRIMARY_ONLY`; a flipped
        // record of a pair that is not proper is seen by those and by `PRIMARY_MAPPED_ONLY`; and a
        // flipped ordinary mate is seen by all four.
        //
        // Coordinate-sorted, so the tool has to sort it into query-name order itself: the order it
        // sorts into decides which record of a name is the one the others are compared against.
        SAMFileHeader inconsistentHeader = header(SAMFileHeader.SortOrder.coordinate);
        java.util.List<SAMRecord> inconsistentReads = reads(inconsistentHeader, true);
        java.util.Map<String, java.util.List<SAMRecord>> inconsistentTemplates = new java.util.LinkedHashMap<>();
        for (SAMRecord r : inconsistentReads) {
            inconsistentTemplates.computeIfAbsent(r.getReadName(), k -> new java.util.ArrayList<>()).add(r);
        }
        int flipped = 0;
        for (java.util.List<SAMRecord> template : inconsistentTemplates.values()) {
            if (template.size() != 2) continue;
            SAMRecord one = template.get(0).getFirstOfPairFlag() ? template.get(0) : template.get(1);
            SAMRecord two = template.get(0).getFirstOfPairFlag() ? template.get(1) : template.get(0);
            int n = Integer.parseInt(one.getReadName().replaceAll("[^0-9]", ""));
            SAMRecord flip;
            if (n % 16 == 0) {
                flip = one;          // secondary
            } else if (n % 18 == 0) {
                flip = two;          // supplementary
            } else if (n % 20 == 0) {
                flip = two;          // unmapped
            } else if (n % 5 == 0) {
                flip = two;          // not a proper pair
            } else if (n % 6 == 0) {
                flip = two;          // an ordinary mate
            } else {
                continue;
            }
            flip.setDuplicateReadFlag(!flip.getDuplicateReadFlag());
            flipped++;
        }
        if (flipped == 0) throw new IllegalStateException("no duplicate flag was flipped");
        writeBam(new File(dir, "inconsistent_duplicates.bam"), inconsistentHeader, inconsistentReads, false);

        // Reads whose ends really are Illumina adapters, for `MarkIlluminaAdapters`. On the random
        // bases of the other fixtures no adapter is ever found, so every accepted row produces the
        // same output and the array covers the search without running it.
        //
        // Three families are planted, and which family a pair gets is what makes `--ADAPTERS`
        // observable. Picard's parser APPENDS to a list argument's default, so every row searches
        // INDEXED, DUAL_INDEXED and PAIRED_END whatever it asks for: a corpus carrying only a
        // PAIRED_END adapter answers the same thing for all nine values. NEXTERA_V2 and
        // TRUSEQ_SMALLRNA are not in that default, so the pairs carrying them are marked only by
        // the rows that name them.
        //
        // Two planting shapes, because the paired rule has two branches. A one-sided plant puts
        // the three prime adapter in read one and leaves read two alone, which is the branch that
        // re-checks the single match against twice the minimum and then marks BOTH reads. A
        // two-sided plant puts the three prime adapter in read one and the five prime adapter in
        // read-order in read two, at the same offset, which is the branch where the two indices
        // agree and the pair is marked immediately.
        SAMFileHeader adapterHeader = header(SAMFileHeader.SortOrder.queryname);
        java.util.List<SAMRecord> adapterReads = reads(adapterHeader, false);
        java.util.Map<String, java.util.List<SAMRecord>> adapterTemplates = new java.util.LinkedHashMap<>();
        for (SAMRecord r : adapterReads) {
            adapterTemplates.computeIfAbsent(r.getReadName(), k -> new java.util.ArrayList<>()).add(r);
        }
        int planted = 0;
        for (java.util.List<SAMRecord> template : adapterTemplates.values()) {
            if (template.size() != 2) continue;
            SAMRecord one = template.get(0).getFirstOfPairFlag() ? template.get(0) : template.get(1);
            SAMRecord two = template.get(0).getFirstOfPairFlag() ? template.get(1) : template.get(0);
            int n = Integer.parseInt(one.getReadName().replaceAll("[^0-9]", ""));
            String fivePrime, threePrime;
            boolean twoSided;
            if (n % 12 == 0) {
                fivePrime = "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT";
                threePrime = "AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG";
                twoSided = false;
            } else if (n % 12 == 4) {
                fivePrime = "AATGATACGGCGACCACCGAGATCTACACNNNNNNNNTCGTCGGCAGCGTCAGATGTGTATAAGAGACAG";
                threePrime = "CTGTCTCTTATACACATCTCCGAGCCCACGAGACNNNNNNNNATCTCGTATGCCGTCTTCTGCTTG";
                twoSided = true;
            } else if (n % 12 == 8) {
                fivePrime = "AATGATACGGCGACCACCGAGATCTACACGTTCAGAGTTCTACAGTCCGACGATC";
                threePrime = "TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG";
                twoSided = false;
            } else {
                continue;
            }
            // 24 bases, which is well over the paired minimum and over the single-end one too.
            plantAdapter(one, threePrime, 24);
            if (twoSided) {
                plantAdapter(two, SequenceUtil.reverseComplement(fivePrime), 24);
            }
            planted++;
        }
        if (planted == 0) throw new IllegalStateException("no adapter was planted");
        writeBam(new File(dir, "adapters.bam"), adapterHeader, adapterReads, false);

        // The same reads, coordinate-sorted and indexed, for the locus walkers. `CollectWgsMetrics`
        // and its siblings refuse a queryname-sorted input outright and small.bam carries no
        // adapter, so without this file PCT_EXC_ADAPTER is zero on every row the array can run.
        // The plants above sit at the 3' end, which `AdapterUtility.isAdapter` never reads: it
        // tests the read's first bases in sequencing order, and only on a read whose mapping
        // quality is zero (one that mapped with confidence is never an adapter read). So the
        // mapped, quality-zero ends of the n % 12 == 0 family get the five prime adapter at their
        // START here as well. Written after adapters.bam, from the records already planted, so
        // no other fixture's bytes move.
        SAMFileHeader adapterCoordinateHeader = header(SAMFileHeader.SortOrder.coordinate);
        for (SAMRecord r : adapterReads) {
            r.setHeader(adapterCoordinateHeader);
            int n = Integer.parseInt(r.getReadName().replaceAll("[^0-9]", ""));
            if (n % 12 == 0 && !r.getReadUnmappedFlag() && r.getMappingQuality() == 0) {
                byte[] bases = r.getReadBases();
                if (r.getReadNegativeStrandFlag()) SequenceUtil.reverseComplement(bases);
                byte[] adapter = "AATGATACGGCGACCACCGAGATCTACAC".getBytes("UTF-8");
                System.arraycopy(adapter, 0, bases, 0, 24);
                if (r.getReadNegativeStrandFlag()) SequenceUtil.reverseComplement(bases);
                r.setReadBases(bases);
            }
        }
        adapterReads.sort(new SAMRecordCoordinateComparator());
        writeBam(new File(dir, "adapters_coordinate.bam"), adapterCoordinateHeader, adapterReads, true);

        // Read pairs that really do repeat, for `EstimateLibraryComplexity`. The tool groups pairs
        // by the first MIN_IDENTICAL_BASES of both ends and counts how big each group of duplicates
        // is; on the random bases of the other fixtures every group holds one pair, every bin is
        // one, and MIN_GROUP_COUNT then drops the lot, so the metrics are zeros whatever the
        // arguments say.
        //
        // Twenty-four families of one to four identical pairs each, which gives bins of one, two,
        // three and four. Every fourth family mutates its later copies -- two bases, inside the
        // default MAX_DIFF_RATE's allowance, or eight, outside it -- so the rate decides whether
        // those copies join the group. Every sixth family is written at quality fifteen, below the
        // default MIN_MEAN_QUALITY, so the quality filter has something to drop.
        //
        // The names are Illumina-style because the optical duplicate finder reads a tile and a
        // position out of them: the first two copies of a family share a tile twenty pixels apart,
        // which is inside the default pixel distance and outside a smaller one, and the later
        // copies sit on tiles of their own.
        //
        // The two read groups alternate by family, never within one, so a family stays in one
        // library and the library split does not cut a group in half.
        SAMFileHeader complexityHeader = header(SAMFileHeader.SortOrder.queryname);
        java.util.List<SAMRecord> complexityReads = new java.util.ArrayList<>();
        Random complexityRng = new Random(20260908L);
        String complexityBases = "ACGT";
        for (int family = 0; family < 24; family++) {
            byte[] readOne = new byte[READ_LENGTH];
            byte[] readTwo = new byte[READ_LENGTH];
            for (int b = 0; b < READ_LENGTH; b++) {
                readOne[b] = (byte) complexityBases.charAt(complexityRng.nextInt(4));
                readTwo[b] = (byte) complexityBases.charAt(complexityRng.nextInt(4));
            }
            int copies = 1 + (family % 4);
            String group = (family % 2 == 0) ? "rg1" : "rg2";
            byte quality = (byte) (family % 6 == 5 ? 15 : 35);
            for (int copy = 0; copy < copies; copy++) {
                byte[] one = readOne.clone();
                if (copy > 0 && family % 4 == 3) {
                    int errors = (family % 8 == 3) ? 2 : 8;
                    for (int e = 0; e < errors; e++) {
                        one[20 + e] = mutateBase(one[20 + e]);
                    }
                }
                int tile = 1101 + (copy < 2 ? 0 : copy);
                int x = 1000 + family * 7 + (copy < 2 ? copy * 20 : copy * 5000);
                int y = 2000 + family * 11 + (copy < 2 ? copy * 20 : copy * 5000);
                String name = String.format("INST:1:FLOWCELL:1:%d:%d:%d", tile, x, y);

                SAMRecord first = new SAMRecord(complexityHeader);
                SAMRecord second = new SAMRecord(complexityHeader);
                for (SAMRecord r : new SAMRecord[] {first, second}) {
                    byte[] quals = new byte[READ_LENGTH];
                    java.util.Arrays.fill(quals, quality);
                    r.setReadName(name);
                    r.setBaseQualities(quals);
                    r.setReferenceIndex(0);
                    r.setCigarString(READ_LENGTH + "M");
                    r.setMappingQuality(60);
                    r.setAttribute("RG", group);
                    r.setReadPairedFlag(true);
                    r.setProperPairFlag(true);
                }
                first.setFirstOfPairFlag(true);
                second.setSecondOfPairFlag(true);
                first.setAlignmentStart(100 + family);
                second.setAlignmentStart(400 + family);
                first.setReadBases(one);
                // The second end is on the negative strand, so the file stores it reverse
                // complemented and the tool has to complement it back before it compares: a
                // corpus whose ends are all forward never runs that path.
                byte[] two = readTwo.clone();
                SequenceUtil.reverseComplement(two);
                second.setReadBases(two);
                second.setReadNegativeStrandFlag(true);
                first.setMateNegativeStrandFlag(true);
                SamPairUtil.setMateInfo(first, second, false);
                complexityReads.add(first);
                complexityReads.add(second);
            }
        }
        complexityReads.sort(new SAMRecordQueryNameComparator());
        writeBam(new File(dir, "complexity.bam"), complexityHeader, complexityReads, false);

        SAMFileHeader unmappedHeader = header(SAMFileHeader.SortOrder.unsorted);
        writeBam(new File(dir, "unmapped.bam"), unmappedHeader, unmapped(unmappedHeader), false);

        // A read-name list, for FilterSamReads' includeReadList / excludeReadList. Every fourth
        // pair, so both filters keep something and drop something.
        try (PrintWriter out = new PrintWriter(new File(dir, "read_names.txt"), "UTF-8")) {
            for (int i = 0; i < READS; i += 8) out.printf("read%04d%n", i);
        }

        // Three small VCFs, for the tools that read variants rather than reads. `variants.vcf` has
        // two samples and a mix of SNPs and an indel, half of its sites also in `dbsnp.vcf`, so a
        // tool that partitions by novelty sees both partitions; `single_sample.vcf` is the
        // one-sample file the tools that refuse more than one need. Written through htsjdk's own
        // writer, index and all, because a hand-written VCF is a fixture whose bugs become the
        // tool's answers.
        writeVcf(new File(dir, "variants.vcf"), chr1, chr2, true);
        writeVcf(new File(dir, "dbsnp.vcf"), chr1, chr2, false);
        writeVcf(new File(dir, "single_sample.vcf"), chr1, chr2, true, 1);

        // Two `CollectQualityYieldMetrics` outputs, for the tools that accumulate metrics files
        // rather than reads. The header comments are what that tool writes, command line and
        // timestamp included, because the accumulator reads past them to the table and a fixture
        // that dropped them would not be the file it is given in practice. The two differ in every
        // counter, so a row that reads one is a different answer from a row that reads the other.
        writeQualityYield(new File(dir, "quality_yield_one.metrics"),
                400, 400, 50, 20000, 20000, 10547, 10547, 5200, 5200, 20456, 20456);
        writeQualityYield(new File(dir, "quality_yield_two.metrics"),
                150, 120, 60, 9000, 7200, 4100, 3300, 2000, 1600, 8800, 7000);
        // Reads carrying a UMI, for `CollectUmiPrevalenceMetrics`. The tool groups the file into
        // duplicate sets and counts the DISTINCT barcodes in each, so a corpus needs reads that
        // duplicate one another and carry different barcodes when they do.
        //
        // The barcode quality filter is the reason for the `BQ` values here. It drops a read whose
        // barcode has NO base under the floor, which is the reverse of what its name says, so a
        // file of well-formed barcodes reports nothing at all. One base at twenty is what keeps a
        // read at the default floor of thirty and drops it at fifteen; a barcode with no low base
        // is dropped either way; and a read with no `BQ` tag at all is never dropped, because the
        // filter returns before it looks.
        //
        // Mapping qualities straddle the default floor of thirty, and every third family is
        // written unpaired, so FILTER_UNPAIRED_READS decides something too.
        SAMFileHeader umiHeader = header(SAMFileHeader.SortOrder.coordinate);
        java.util.List<SAMRecord> umiReads = new java.util.ArrayList<>();
        String[] umis = {"AACCGGTT", "TTGGCCAA", "ACACGTGT"};
        for (int family = 0; family < 18; family++) {
            int copies = 1 + (family % 3);
            int start = 100 + family * 40;
            boolean unpaired = family % 3 == 2;
            for (int copy = 0; copy < copies; copy++) {
                SAMRecord read = new SAMRecord(umiHeader);
                byte[] bases = new byte[READ_LENGTH];
                byte[] quals = new byte[READ_LENGTH];
                for (int b = 0; b < READ_LENGTH; b++) {
                    bases[b] = (byte) "ACGT".charAt((family + b) % 4);
                    quals[b] = 35;
                }
                read.setReadName(String.format("umi%04d_%d", family, copy));
                read.setReadBases(bases);
                read.setBaseQualities(quals);
                read.setReferenceIndex(0);
                read.setAlignmentStart(start);
                read.setCigarString(READ_LENGTH + "M");
                // Every fourth family is mapped under the default floor of thirty.
                read.setMappingQuality(family % 4 == 3 ? 25 : 60);
                read.setAttribute("RG", family % 2 == 0 ? "rg1" : "rg2");
                read.setAttribute("RX", umis[copy % umis.length]);
                if (copy == 1) {
                    // One base at twenty: kept at the default floor, dropped at a lower one.
                    read.setAttribute("BQ", "I5IIIIII");
                } else if (copy == 2) {
                    // No base under any floor the array uses, so this read is always dropped.
                    read.setAttribute("BQ", "IIIIIIII");
                }
                umiReads.add(read);
                if (!unpaired) {
                    // The mate is written too: the duplicate-set iterator reads the mate's cigar
                    // out of the `MC` tag, and refuses a paired read that has none.
                    SAMRecord mate = new SAMRecord(umiHeader);
                    mate.setReadName(read.getReadName());
                    mate.setReadBases(bases.clone());
                    mate.setBaseQualities(quals.clone());
                    mate.setReferenceIndex(0);
                    mate.setAlignmentStart(start + 200);
                    mate.setCigarString(READ_LENGTH + "M");
                    mate.setMappingQuality(read.getMappingQuality());
                    mate.setAttribute("RG", read.getStringAttribute("RG"));
                    mate.setAttribute("RX", read.getStringAttribute("RX"));
                    if (read.getStringAttribute("BQ") != null) {
                        mate.setAttribute("BQ", read.getStringAttribute("BQ"));
                    }
                    mate.setReadNegativeStrandFlag(true);
                    read.setReadPairedFlag(true);
                    mate.setReadPairedFlag(true);
                    read.setFirstOfPairFlag(true);
                    mate.setSecondOfPairFlag(true);
                    read.setProperPairFlag(true);
                    mate.setProperPairFlag(true);
                    SamPairUtil.setMateInfo(read, mate, true);
                    umiReads.add(mate);
                }
            }
        }
        umiReads.sort(new SAMRecordCoordinateComparator());
        writeBam(new File(dir, "umi.bam"), umiHeader, umiReads, false);

        // Reads carrying cell and molecular barcodes, for `SamToFastqWithTags`. That tool writes
        // the ordinary read FASTQ and, beside it, one FASTQ per SEQUENCE_TAG_GROUP whose reads are
        // built from TAG VALUES rather than from bases, so a corpus needs records that carry the
        // tags and records that do not: a group naming a tag a read is missing is a refusal, and
        // that is a row of the array.
        //
        // Every read carries every tag the array's groups name: a read missing one is a refusal
        // ("does have a value for tag"), and a corpus that refused most rows would measure the
        // check rather than the writing.
        SAMFileHeader taggedHeader = header(SAMFileHeader.SortOrder.queryname);
        java.util.List<SAMRecord> taggedReads = reads(taggedHeader, false);
        for (SAMRecord r : taggedReads) {
            int n = Integer.parseInt(r.getReadName().replaceAll("[^0-9]", ""));
            r.setAttribute("CB", "ACGTAC" + (char) ('A' + (n % 4)));
            r.setAttribute("CY", "IIIIIII");
            r.setAttribute("UR", "TTGCA");
            r.setAttribute("UY", "IIIII");
        }
        writeBam(new File(dir, "tagged.bam"), taggedHeader, taggedReads, false);

        // Reads carrying `OQ` and no `MC`, for `RevertOriginalBaseQualitiesAndAddMateCigar`. The
        // tool skips a file outright when its first informative record has no OQ and its mate
        // cigar is already there, so a corpus needs a file it cannot skip: these carry original
        // qualities that differ from the current ones, and no mate cigar at all.
        //
        // Every record is primary. The tool's mate-info pass has a separate path for secondary and
        // supplementary records, and a corpus that mixed them in would measure that path here
        // instead of this one.
        SAMFileHeader originalQualsHeader = header(SAMFileHeader.SortOrder.coordinate);
        java.util.List<SAMRecord> originalQualsReads = new java.util.ArrayList<>();
        for (int pair = 0; pair < 30; pair++) {
            SAMRecord first = new SAMRecord(originalQualsHeader);
            SAMRecord second = new SAMRecord(originalQualsHeader);
            int start = 100 + pair * 30;
            for (SAMRecord r : new SAMRecord[] {first, second}) {
                byte[] bases = new byte[READ_LENGTH];
                byte[] quals = new byte[READ_LENGTH];
                StringBuilder original = new StringBuilder();
                for (int b = 0; b < READ_LENGTH; b++) {
                    bases[b] = (byte) "ACGT".charAt((pair + b) % 4);
                    quals[b] = (byte) (20 + (b % 10));
                    // The original qualities are five higher, so a reverted read is a different
                    // read and the argument that reverts them is observable.
                    original.append((char) (33 + 25 + (b % 10)));
                }
                r.setReadName(String.format("oq%04d", pair));
                r.setReadBases(bases);
                r.setBaseQualities(quals);
                r.setReferenceIndex(0);
                r.setCigarString(READ_LENGTH + "M");
                r.setMappingQuality(60);
                r.setAttribute("RG", pair % 2 == 0 ? "rg1" : "rg2");
                // Every third pair keeps its current qualities, so RESTORE_ORIGINAL_QUALITIES is
                // not all-or-nothing over the file.
                if (pair % 3 != 2) r.setAttribute("OQ", original.toString());
                r.setReadPairedFlag(true);
                r.setProperPairFlag(true);
            }
            first.setFirstOfPairFlag(true);
            second.setSecondOfPairFlag(true);
            first.setAlignmentStart(start);
            second.setAlignmentStart(start + 150);
            second.setReadNegativeStrandFlag(true);
            first.setMateNegativeStrandFlag(true);
            first.setMateReferenceIndex(0);
            second.setMateReferenceIndex(0);
            first.setMateAlignmentStart(second.getAlignmentStart());
            second.setMateAlignmentStart(first.getAlignmentStart());
            first.setInferredInsertSize(150 + READ_LENGTH);
            second.setInferredInsertSize(-(150 + READ_LENGTH));
            originalQualsReads.add(first);
            originalQualsReads.add(second);
        }
        originalQualsReads.sort(new SAMRecordCoordinateComparator());
        writeBam(new File(dir, "original_quals.bam"), originalQualsHeader, originalQualsReads, false);

        // A trio and its pedigree, for `FindMendelianViolations`. The sites are chosen for the
        // classes the tool counts: a de novo heterozygote in a child of two reference parents, a
        // heterozygote from two homozygous-variant parents, a homozygote from a reference parent
        // and a variant one, a homozygous-variant child of a reference and a heterozygote -- and
        // inherited sites that are no violation at all, so the denominator is not the numerator.
        //
        // The genotype qualities straddle the default MIN_GQ of thirty and the depths the default
        // MIN_DP of zero, and one child heterozygote has an allele depth so lopsided (10 and 1)
        // that MIN_HET_FRACTION refuses to judge it.
        writeTrioVcf(new File(dir, "trio.vcf"), chr1);
        try (PrintWriter out = new PrintWriter(new File(dir, "trio.ped"), "UTF-8")) {
            // family, individual, father, mother, sex (1 male, 2 female), phenotype
            out.print("fam1\tchild\tfather\tmother\t1\t2\n");
            out.print("fam1\tfather\t0\t0\t1\t1\n");
            out.print("fam1\tmother\t0\t0\t2\t1\n");
        }

        // A UCSC chain, for the tools that lift coordinates from one build to another. It maps the
        // first thousand bases of chr1 five hundred bases to the right and the first five hundred
        // of chr2 onto themselves, and covers nothing else: an interval past chr1:1000 lifts only
        // in part, which is what MIN_LIFTOVER_PCT decides about, and one past chr1:1500 does not
        // lift at all.
        try (PrintWriter out = new PrintWriter(new File(dir, "lift.chain"), "UTF-8")) {
            out.print("chain 1000 chr1 " + CHR1 + " + 0 1000 chr1 " + CHR1 + " + 500 1500 1\n");
            out.print("1000\n");
            out.print("\n");
            out.print("chain 1000 chr2 " + CHR2 + " + 0 500 chr2 " + CHR2 + " + 0 500 2\n");
            out.print("500\n");
            out.print("\n");
        }

        // The same two contigs mapped onto themselves end to end, so that a tool given this chain
        // lifts everything and one given the other drops what falls outside it: the pair is what
        // makes CHAIN an argument that decides something.
        try (PrintWriter out = new PrintWriter(new File(dir, "lift_all.chain"), "UTF-8")) {
            out.print("chain 1000 chr1 " + CHR1 + " + 0 " + CHR1 + " chr1 " + CHR1 + " + 0 " + CHR1 + " 1\n");
            out.print(CHR1 + "\n");
            out.print("\n");
            out.print("chain 1000 chr2 " + CHR2 + " + 0 " + CHR2 + " chr2 " + CHR2 + " + 0 " + CHR2 + " 2\n");
            out.print(CHR2 + "\n");
            out.print("\n");
        }

        // A haplotype database, for the fingerprinting tools. Two blocks of two SNPs and one of
        // one, anchored the way the format wants: the first SNP of a block names no anchor and
        // every other names it. Two of the SNPs sit past chr1:1500, where lift.chain covers
        // nothing, so a liftover drops them and reports the failure in its exit code.
        try (PrintWriter out = new PrintWriter(new File(dir, "haplotypes.txt"), "UTF-8")) {
            out.print("@HD\tVN:1.0\tSO:coordinate\n");
            out.print("@SQ\tSN:chr1\tLN:" + CHR1 + "\n");
            out.print("@SQ\tSN:chr2\tLN:" + CHR2 + "\n");
            out.print("#CHROMOSOME\tPOSITION\tNAME\tMAJOR_ALLELE\tMINOR_ALLELE\tMAF\tANCHOR_SNP\tPANELS\n");
            out.print("chr1\t100\trs1\tA\tC\t0.30\t\t\n");
            out.print("chr1\t200\trs2\tG\tT\t0.25\trs1\t\n");
            out.print("chr1\t1600\trs3\tC\tG\t0.40\t\t\n");
            out.print("chr1\t1700\trs4\tT\tA\t0.15\trs3\t\n");
            out.print("chr2\t100\trs5\tG\tA\t0.20\t\t\n");
        }

        writeIntervals(new File(dir, "targets.interval_list"));
        writeBed(new File(dir, "targets.bed"));
        writeMixedBed(new File(dir, "targets_mixed.bed"));
        writeMixedIntervals(new File(dir, "targets_mixed.interval_list"));
        writeBaits(new File(dir, "baits.interval_list"));
        writeDescribedFasta(new File(dir, "described.fasta"), chr2);
        writeDict(new File(dir, "ref.dict"), chr1, chr2);

        // Hand-written VCFs for the picard.vcf manipulation tools (MakeSitesOnlyVcf,
        // RenameSampleInVcf, VcfToIntervalList, SortVcf, UpdateVcfSequenceDictionary). These are
        // text on purpose, unlike variants.vcf: what these tools write depends on bytes htsjdk's
        // own writer never produces. A genotype block whose FORMAT is not `GT` first and the rest
        // sorted, or whose missing trailing fields were not trimmed, is copied verbatim when the
        // file's sample names are in sorted order and re-encoded when they are not, so the same
        // records are written under both sample orders.
        writeVcfUtilityFixtures(dir);
        writeVcfMergeFixtures(dir);
        writeVcfGatherFixtures(dir);
        writeVcfSplitFixture(dir);
        writeVcfConverterFixture(dir);
        writeVcfFixHeaderFixtures(dir);
        writeFastq(new File(dir, "reads_1.fastq"), 1);
        writeFastq(new File(dir, "reads_2.fastq"), 2);

        // Detail tables in the shape `CollectSequencingArtifactMetrics` writes them, for
        // `ConvertSequencingArtifactToOxoG`, which reads two of them and nothing else. Written
        // through Picard's own beans and `MetricsFile`, with no header, so the bytes do not carry
        // the command line or the second they were made. `artifacts` is the whole table for two
        // libraries, where CCG, its reverse complement CGG, and GCA count nothing at all, which
        // is what puts an infinite rate and a NaN one in the output; `artifacts_subset` is one
        // library and the four contexts a CONTEXTS_TO_PRINT run leaves, each with its reverse
        // complement; and `artifacts_broken` keeps `ACA` without `TGT`, the reverse complement
        // the tool looks up.
        writeArtifactTables(dir, "artifacts", new String[] {"lib1", "lib2"}, "sample1,sample2", null);
        writeArtifactTables(dir, "artifacts_subset", new String[] {"lib1"}, "sample1",
                new String[] {"ACA", "CCC", "GGG", "TGT"});
        writeArtifactTables(dir, "artifacts_broken", new String[] {"lib1"}, "sample1",
                new String[] {"ACA", "CCC", "GGG"});

        // Bisulfite-converted reads, for `CollectRrbsMetrics`. The other fixtures' bases are
        // random, so every read fails the tool's mismatch filter and the metrics count nothing.
        // These copy the reference and convert it the way bisulfite does: a cytosine outside a
        // CpG is read as T nineteen times in twenty, one inside a CpG (methylated, so protected)
        // one time in four. A reverse-strand read converts the other strand, so its stored bases
        // show G read as A. Every eighth read also carries eight mismatches, which the default
        // MAX_MISMATCH_RATE rejects; every tenth is four bases long, under the default
        // MINIMUM_READ_LENGTH; one in thirty is placed but unmapped; and the qualities straddle
        // the two quality thresholds.
        SAMFileHeader rrbsHeader = header(SAMFileHeader.SortOrder.coordinate);
        writeBam(new File(dir, "rrbs.bam"), rrbsHeader, rrbsReads(rrbsHeader, chr1, chr2), false);

        // GFF3 for `SortGff`, as text: what the tool writes depends on lines htsjdk's writer never
        // produces (see writeGffFixtures).
        writeGffFixtures(dir);

        // A genotyping-array VCF and zCall's PED, MAP and thresholds, for `MergePedIntoVcf`
        // (see writeZCallFixtures).
        writeZCallFixtures(dir);

        // Unaligned pairs for `BamToBfq` (see bfqReads).
        SAMFileHeader bfqHeader = header(SAMFileHeader.SortOrder.queryname);
        writeBam(new File(dir, "bfq.bam"), bfqHeader, bfqReads(bfqHeader, false), false);
        SAMFileHeader bfqOddHeader = header(SAMFileHeader.SortOrder.unsorted);
        writeBam(new File(dir, "bfq_odd.bam"), bfqOddHeader, bfqReads(bfqOddHeader, true), false);

        // Insert-size metrics files that agree and disagree, for `CompareMetrics`.
        writeCompareMetricsFixtures(dir);

        // An unmapped bam and its alignments, for `MergeBamAlignment` (see writeMbaFixtures).
        writeMbaFixtures(dir, chr1, chr2);

        System.out.println("fixtures written to " + dir.getAbsolutePath());
        for (File f : dir.listFiles()) {
            System.out.printf("%s\t%d%n", f.getName(), f.length());
        }
    }

    /** A reference with a fixed but non-uniform base composition, so GC-dependent tools vary. */
    static String reference(int length, long seed) {
        Random rng = new Random(seed);
        char[] bases = new char[length];
        String alphabet = "ACGT";
        for (int i = 0; i < length; i++) {
            // A GC-rich stretch in the middle, and a run of Ns, so the reference is not featureless.
            if (i > length / 2 && i < length / 2 + 100) bases[i] = (i % 2 == 0) ? 'G' : 'C';
            else if (i > length - 60 && i < length - 40) bases[i] = 'N';
            else bases[i] = alphabet.charAt(rng.nextInt(4));
        }
        return new String(bases);
    }

    static void writeFasta(File f, String chr1, String chr2) throws Exception {
        try (PrintWriter p = new PrintWriter(f)) {
            writeContig(p, "chr1", chr1);
            writeContig(p, "chr2", chr2);
        }
    }

    static void writeContig(PrintWriter p, String name, String bases) {
        p.println(">" + name);
        for (int i = 0; i < bases.length(); i += 60) {
            p.println(bases.substring(i, Math.min(i + 60, bases.length())));
        }
    }

    /**
     * The FASTA index, whose offsets have to be the file's real ones.
     *
     * The byte length of a contig is its bases plus one newline per line, and the last line is
     * short: `ceil(len / 60) * 61` over-counts it by `60 - (len % 60)`. chr1 is 2000 bases, so the
     * old arithmetic put chr2's offset 40 bytes past where chr2 begins.
     *
     * Nothing noticed until a tool took the indexed path. `ReferenceSequenceFileFactory` opens
     * `IndexedFastaSequenceFile` only when the caller asks for names truncated at whitespace
     * ("Using faidx requires truncateNamesAtWhitespace"), so NormalizeFasta read this index only
     * with TRUNCATE_SEQUENCE_NAMES_AT_WHITESPACE=true, and then sliced chr2 from the wrong byte:
     * its output carried the file's own line terminators as if they were bases, one of them
     * producing an empty line. The covering array is what ran that combination.
     */
    static void writeFai(File f, String chr1, String chr2) throws Exception {
        try (PrintWriter p = new PrintWriter(f)) {
            int lineWidth = 60, lineBytes = 61;
            long offset1 = ">chr1\n".length();
            long chr1Bytes = bytesOnDisk(chr1.length(), lineWidth);
            long offset2 = offset1 + chr1Bytes + ">chr2\n".length();
            p.printf("chr1\t%d\t%d\t%d\t%d%n", chr1.length(), offset1, lineWidth, lineBytes);
            p.printf("chr2\t%d\t%d\t%d\t%d%n", chr2.length(), offset2, lineWidth, lineBytes);
        }
    }

    /**
     * The same corpus under a single read group.
     *
     * Every other fixture carries the same two read groups, so a tool whose whole output is a
     * function of the header's `@RG` records answers identically on all of them:
     * CalculateReadGroupChecksum's array was nine rows and one digest, which covers its arguments
     * without testing them. This file differs in exactly the thing that tool reads.
     */
    static SAMFileHeader header(SAMFileHeader.SortOrder order) {
        SAMFileHeader h = new SAMFileHeader();
        SAMSequenceDictionary d = new SAMSequenceDictionary();
        d.addSequence(new SAMSequenceRecord("chr1", CHR1));
        d.addSequence(new SAMSequenceRecord("chr2", CHR2));
        h.setSequenceDictionary(d);
        h.setSortOrder(order);
        // Two read groups in two libraries: the multi-level collectors have a LIBRARY and a
        // READ_GROUP accumulation level, and one read group would leave both untested.
        for (String[] rg : new String[][] {{"rg1", "lib1", "sample1"}, {"rg2", "lib2", "sample2"}}) {
            SAMReadGroupRecord r = new SAMReadGroupRecord(rg[0]);
            r.setLibrary(rg[1]);
            r.setSample(rg[2]);
            r.setPlatform("ILLUMINA");
            r.setPlatformUnit("unit-" + rg[0]);
            h.addReadGroup(r);
        }
        return h;
    }

    /** The next base in ACGT order, for planting a mismatch that is still a real base. */
    static byte mutateBase(byte base) {
        switch (base) {
            case 'A': return 'C';
            case 'C': return 'G';
            case 'G': return 'T';
            default: return 'A';
        }
    }

    /** One `QualityYieldMetrics` row, written the way `CollectQualityYieldMetrics` writes it. */
    static void writeQualityYield(File f, long totalReads, long pfReads, int readLength,
                                  long totalBases, long pfBases, long q20, long pfQ20, long q30,
                                  long pfQ30, long q20Yield, long pfQ20Yield) throws Exception {
        try (PrintWriter out = new PrintWriter(f, "UTF-8")) {
            out.print("## htsjdk.samtools.metrics.StringHeader\n");
            out.print("# CollectQualityYieldMetrics INPUT=/work/fixtures/small.bam OUTPUT=/work/out/output.txt\n");
            out.print("## htsjdk.samtools.metrics.StringHeader\n");
            out.print("# Started on: Mon Sep 07 00:00:00 UTC 2026\n");
            out.print("\n");
            out.print("## METRICS CLASS\tpicard.analysis.CollectQualityYieldMetrics$QualityYieldMetrics\n");
            out.print("TOTAL_READS\tPF_READS\tREAD_LENGTH\tTOTAL_BASES\tPF_BASES\tQ20_BASES\tPF_Q20_BASES\tQ30_BASES\tPF_Q30_BASES\tQ20_EQUIVALENT_YIELD\tPF_Q20_EQUIVALENT_YIELD\n");
            out.printf("%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d%n",
                    totalReads, pfReads, readLength, totalBases, pfBases, q20, pfQ20, q30, pfQ30,
                    q20Yield, pfQ20Yield);
            out.print("\n");
        }
    }

    /**
     * Put the first `length` bases of an adapter at the end of a read, IN READ ORDER.
     *
     * A record on the negative strand stores its bases reverse complemented, and the tool searches
     * what the sequencer read rather than what the file stores: it reverse complements a copy
     * before it looks. Planting into the stored bases would therefore put the adapter at the front
     * of half the reads, where the search never looks, so the plant is done on the read-order copy
     * and complemented back.
     */
    static void plantAdapter(SAMRecord read, String adapter, int length) throws Exception {
        byte[] bases = read.getReadBases();
        if (read.getReadNegativeStrandFlag()) SequenceUtil.reverseComplement(bases);
        byte[] planted = adapter.substring(0, length).getBytes("UTF-8");
        System.arraycopy(planted, 0, bases, bases.length - planted.length, planted.length);
        if (read.getReadNegativeStrandFlag()) SequenceUtil.reverseComplement(bases);
        read.setReadBases(bases);
    }

    static java.util.List<SAMRecord> reads(SAMFileHeader header, boolean coordinateSorted) {
        Random rng = new Random(20260729L);
        java.util.List<SAMRecord> out = new java.util.ArrayList<>();
        String alphabet = "ACGT";

        for (int i = 0; i < READS; i += 2) {
            String name = String.format("read%04d", i);
            boolean chr2 = i % 8 == 0;
            int contig = chr2 ? 1 : 0;
            int limit = (chr2 ? CHR2 : CHR1) - READ_LENGTH - 10;
            int start = 1 + rng.nextInt(limit);

            SAMRecord first = new SAMRecord(header);
            SAMRecord second = new SAMRecord(header);
            for (SAMRecord r : new SAMRecord[] {first, second}) {
                byte[] bases = new byte[READ_LENGTH];
                byte[] quals = new byte[READ_LENGTH];
                for (int b = 0; b < READ_LENGTH; b++) {
                    bases[b] = (byte) alphabet.charAt(rng.nextInt(4));
                    // Qualities span the CollectQualityYieldMetrics thresholds (Q20, Q30) rather
                    // than sitting above both, which would make those counters constant.
                    quals[b] = (byte) (2 + rng.nextInt(38));
                }
                // A no-call every so often, so the base-distribution and N-handling paths run.
                if (i % 10 == 0) bases[3] = 'N';
                r.setReadName(name);
                r.setReadBases(bases);
                r.setBaseQualities(quals);
                r.setReferenceIndex(contig);
                r.setMappingQuality(i % 12 == 0 ? 0 : 20 + rng.nextInt(40));
                r.setAttribute("RG", i % 4 == 0 ? "rg2" : "rg1");
            }

            first.setAlignmentStart(start);
            second.setAlignmentStart(start + 100 <= limit ? start + 100 : start);
            first.setCigarString(cigarFor(i));
            second.setCigarString(READ_LENGTH + "M");
            first.setReadPairedFlag(true);
            second.setReadPairedFlag(true);
            first.setFirstOfPairFlag(true);
            second.setSecondOfPairFlag(true);
            first.setReadNegativeStrandFlag(i % 3 == 0);
            second.setReadNegativeStrandFlag(!first.getReadNegativeStrandFlag());
            first.setProperPairFlag(i % 5 != 0);
            second.setProperPairFlag(i % 5 != 0);
            first.setDuplicateReadFlag(i % 14 == 0);
            second.setDuplicateReadFlag(i % 14 == 0);
            // Secondary and supplementary are properties of an *alignment*, so htsjdk's validation
            // rejects them on an unmapped read ("Supplementary alignment flag should not be set for
            // unaligned read"). The first covering-array run hit that: nine of eleven rows failed,
            // and one class of failure was this fixture rather than the argument under test. A
            // corpus that is invalid under STRICT tests the validator, not the tool.
            boolean secondUnmapped = i % 20 == 0;
            if (i % 16 == 0) first.setNotPrimaryAlignmentFlag(true);
            if (i % 18 == 0 && !secondUnmapped) second.setSupplementaryAlignmentFlag(true);

            // One pair in twenty has an unmapped mate: the paired-metrics and mate-info paths
            // behave differently there, and a corpus of clean pairs never reaches them.
            if (secondUnmapped) {
                second.setReadUnmappedFlag(true);
                second.setAlignmentStart(first.getAlignmentStart());
                second.setCigarString("*");
                second.setMappingQuality(0);
                first.setMateUnmappedFlag(true);
            }
            SamPairUtil.setMateInfo(first, second, false);
            out.add(first);
            out.add(second);
        }

        if (coordinateSorted) {
            out.sort(new SAMRecordCoordinateComparator());
        } else {
            out.sort(new SAMRecordQueryNameComparator());
        }
        return out;
    }

    /** Soft clips, an insertion and a deletion, so cigar-walking tools take more than one branch. */
    static String cigarFor(int i) {
        switch (i % 6) {
            case 0: return READ_LENGTH + "M";
            case 1: return "5S" + (READ_LENGTH - 5) + "M";
            case 2: return (READ_LENGTH - 8) + "M8S";
            case 3: return "20M2I" + (READ_LENGTH - 22) + "M";
            case 4: return "20M3D" + (READ_LENGTH - 20) + "M";
            default: return "10M5N" + (READ_LENGTH - 10) + "M";
        }
    }

    static java.util.List<SAMRecord> unmapped(SAMFileHeader header) {
        java.util.List<SAMRecord> out = new java.util.ArrayList<>();
        Random rng = new Random(20260731L);
        for (int i = 0; i < 40; i += 2) {
            SAMRecord first = new SAMRecord(header);
            SAMRecord second = new SAMRecord(header);
            for (SAMRecord r : new SAMRecord[] {first, second}) {
                byte[] bases = new byte[READ_LENGTH];
                byte[] quals = new byte[READ_LENGTH];
                for (int b = 0; b < READ_LENGTH; b++) {
                    bases[b] = (byte) "ACGT".charAt(rng.nextInt(4));
                    quals[b] = (byte) (2 + rng.nextInt(38));
                }
                r.setReadName(String.format("unmapped%04d", i));
                r.setReadBases(bases);
                r.setBaseQualities(quals);
                r.setReadUnmappedFlag(true);
                r.setReferenceIndex(SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX);
                r.setAlignmentStart(SAMRecord.NO_ALIGNMENT_START);
                r.setMappingQuality(0);
                r.setAttribute("RG", "rg1");
                r.setReadPairedFlag(true);
                r.setMateUnmappedFlag(true);
            }
            first.setFirstOfPairFlag(true);
            second.setSecondOfPairFlag(true);
            out.add(first);
            out.add(second);
        }
        return out;
    }

    /**
     * The caller passes -Dsamjdk.try_use_intel_deflater=false: the fixture must be
     * byte-reproducible, and the GKL deflater emits different bytes than zlib for the same input.
     * The oracle contract pins the JDK deflater for the same reason.
     *
     * Only the coordinate-sorted BAM is indexed; indexing a queryname-sorted or unsorted file is
     * an error, not an option.
     */
    static void writeBam(File f, SAMFileHeader header, java.util.List<SAMRecord> records,
                         boolean index) {
        SAMFileWriterFactory factory = new SAMFileWriterFactory().setUseAsyncIo(false);
        factory.setCreateIndex(index);
        try (SAMFileWriter w = factory.makeBAMWriter(header, true, f)) {
            for (SAMRecord r : records) w.addAlignment(r);
        }
    }

    static void writeSam(File f, SAMFileHeader header, java.util.List<SAMRecord> records) {
        try (SAMFileWriter w = new SAMFileWriterFactory().makeSAMWriter(header, true, f)) {
            for (SAMRecord r : records) w.addAlignment(r);
        }
    }

    static void writeIntervals(File f) throws Exception {
        try (PrintWriter p = new PrintWriter(f)) {
            p.println("@HD\tVN:1.6");
            p.printf("@SQ\tSN:chr1\tLN:%d%n", CHR1);
            p.printf("@SQ\tSN:chr2\tLN:%d%n", CHR2);
            p.println("chr1\t100\t400\t+\ttarget1");
            p.println("chr1\t900\t1200\t+\ttarget2");
            p.println("chr2\t50\t200\t-\ttarget3");
        }
    }

    /**
     * The same three targets as the interval list, in BED coordinates.
     *
     * A BED start is 0-based and its end is exclusive, where an interval list is 1-based and
     * inclusive, so the same target is one lower on the left here. Writing both from one set of
     * numbers is the point: a tool that converts between them can then be checked against a
     * fixture that says what the answer is, rather than against its own output.
     */
    static void writeBed(File f) throws Exception {
        try (PrintWriter p = new PrintWriter(f)) {
            p.println("chr1\t99\t400\ttarget1\t0\t+");
            p.println("chr1\t899\t1200\ttarget2\t0\t+");
            p.println("chr2\t49\t200\ttarget3\t0\t-");
        }
    }

    /**
     * A FASTA whose headers carry a description and whose lines are not the output length.
     *
     * ref.fasta has bare contig names and is already wrapped at the length NormalizeFasta writes,
     * so TRUNCATE_SEQUENCE_NAMES_AT_WHITESPACE has nothing to truncate and normalizing is the
     * identity: the array covers both arguments without observing either. Here each header is
     * "name description", so truncation changes the header line, and the bases are wrapped at 37
     * rather than 100, so normalizing rewraps them.
     *
     * It deliberately has no .fai beside it. With one, ReferenceSequenceFileFactory opens the
     * indexed reader, whose index would have to agree with the names; without one it opens
     * FastaSequenceFile, which is the path this port reproduces.
     */
    static void writeDescribedFasta(File f, String chr2) throws Exception {
        try (PrintWriter p = new PrintWriter(f)) {
            p.println(">seq1 first sequence, described");
            for (int i = 0; i < chr2.length(); i += 37) {
                p.println(chr2.substring(i, Math.min(i + 37, chr2.length())));
            }
            p.println(">seq2\ta tab-separated description");
            for (int i = 0; i < 120; i += 37) {
                p.println(chr2.substring(i, Math.min(i + 37, 120)));
            }
        }
    }

    /**
     * An interval list whose order is not the coordinate order.
     *
     * targets.interval_list is already sorted, so SORT produces the same file with it on or off
     * and an array over that argument covers it without testing it. Here chr2 leads, the chr1
     * entries are out of order, and both strands appear, so sorting moves lines and the
     * strand-then-name tiebreak of IntervalCoordinateComparator is reachable.
     */
    static void writeMixedIntervals(File f) throws Exception {
        try (PrintWriter p = new PrintWriter(f)) {
            p.println("@HD\tVN:1.6");
            p.printf("@SQ\tSN:chr1\tLN:%d%n", CHR1);
            p.printf("@SQ\tSN:chr2\tLN:%d%n", CHR2);
            p.println("chr2\t50\t200\t-\ttargetB");
            p.println("chr1\t300\t500\t+\ttargetC");
            p.println("chr1\t100\t400\t+\ttargetA");
            p.println("chr1\t600\t700\t-\ttargetD");
        }
    }

    /**
     * Baits for the hybrid-selection tools, against targets.interval_list's targets.
     *
     * Using the targets as their own baits makes every on-target base on-bait and leaves the
     * near-bait band empty, so the bait columns would only restate the target ones. These overhang
     * two targets, fall short of the third, and add one bait over no target at all, so ON_, NEAR_
     * and OFF_BAIT_BASES and BAIT_DESIGN_EFFICIENCY all move. The file name is what the tool
     * reports as BAIT_SET when none is given.
     */
    static void writeBaits(File f) throws Exception {
        try (PrintWriter p = new PrintWriter(f)) {
            p.println("@HD\tVN:1.6");
            p.printf("@SQ\tSN:chr1\tLN:%d%n", CHR1);
            p.printf("@SQ\tSN:chr2\tLN:%d%n", CHR2);
            p.println("chr1\t80\t420\t+\tbait1");
            p.println("chr1\t950\t1180\t+\tbait2");
            p.println("chr1\t1500\t1560\t+\tbait3");
            p.println("chr2\t60\t190\t-\tbait4");
        }
    }

    /**
     * A BED the interval tools' arguments can actually be observed on.
     *
     * targets.bed is already sorted, disjoint and length-nonzero, so SORT, UNIQUE and
     * KEEP_LENGTH_ZERO_INTERVALS all produce the same file on it: the array covers those
     * arguments without testing them, which the runner says out loud. This one is built so that
     * each of the three changes the output.
     *
     * Out of coordinate order, so SORT moves lines. Two overlapping features and two abutting
     * ones, so UNIQUE merges and concatenates names. One feature whose BED start equals its end,
     * which becomes `start == end + 1` and is dropped unless KEEP_LENGTH_ZERO_INTERVALS is set.
     */
    static void writeMixedBed(File f) throws Exception {
        try (PrintWriter p = new PrintWriter(f)) {
            p.println("chr2\t49\t200\ttargetB\t0\t-");
            p.println("chr1\t299\t500\ttargetC\t0\t+");
            p.println("chr1\t99\t400\ttargetA\t0\t+");
            p.println("chr1\t599\t700\ttargetD\t0\t+");
            p.println("chr1\t699\t800\ttargetE\t0\t+");
            p.println("chr1\t900\t900\tzeroLength\t0\t+");
        }
    }

    /**
     * The sequence dictionary, as its own file.
     *
     * `SAMSequenceDictionaryExtractor` reads a FASTA's dictionary through
     * `ReferenceSequenceFileFactory`, which does not derive one: it looks for the `.dict` beside
     * the reference and throws "Could not find dictionary next to reference file" when there is
     * none. Every tool taking a SEQUENCE_DICTIONARY therefore needed this file before it could be
     * given an array at all.
     */
    static void writeTrioVcf(File f, String chr1) throws Exception {
        htsjdk.samtools.SAMSequenceDictionary dict = new htsjdk.samtools.SAMSequenceDictionary();
        dict.addSequence(new SAMSequenceRecord("chr1", chr1.length()));
        dict.addSequence(new SAMSequenceRecord("chr2", CHR2));

        java.util.Set<htsjdk.variant.vcf.VCFHeaderLine> lines = new java.util.LinkedHashSet<>();
        lines.add(new htsjdk.variant.vcf.VCFFormatHeaderLine(
                "GT", 1, htsjdk.variant.vcf.VCFHeaderLineType.String, "Genotype"));
        lines.add(new htsjdk.variant.vcf.VCFFormatHeaderLine(
                "GQ", 1, htsjdk.variant.vcf.VCFHeaderLineType.Integer, "Genotype quality"));
        lines.add(new htsjdk.variant.vcf.VCFFormatHeaderLine(
                "DP", 1, htsjdk.variant.vcf.VCFHeaderLineType.Integer, "Depth"));
        lines.add(new htsjdk.variant.vcf.VCFFormatHeaderLine(
                "AD", htsjdk.variant.vcf.VCFHeaderLineCount.R,
                htsjdk.variant.vcf.VCFHeaderLineType.Integer, "Allele depths"));
        lines.add(new htsjdk.variant.vcf.VCFFormatHeaderLine(
                "PL", htsjdk.variant.vcf.VCFHeaderLineCount.G,
                htsjdk.variant.vcf.VCFHeaderLineType.Integer, "Phred-scaled likelihoods"));

        htsjdk.variant.vcf.VCFHeader header = new htsjdk.variant.vcf.VCFHeader(
                lines, java.util.Arrays.asList("father", "mother", "child"));
        header.setSequenceDictionary(dict);

        // position, father, mother, child, child's allele depths, genotype quality
        int[][] rows = {
                //   pos  fa fa  mo mo  ch ch   AD0 AD1  GQ
                {  100,  0, 0,   0, 0,   0, 1,   6,  6,  40 },  // de novo het
                {  200,  1, 1,   1, 1,   0, 1,   6,  6,  40 },  // het from two hom-var parents
                {  300,  0, 0,   1, 1,   0, 0,   9,  0,  40 },  // hom from ref x hom-var
                {  400,  0, 0,   0, 1,   1, 1,   0,  9,  40 },  // hom-var from ref x het
                {  500,  0, 1,   0, 1,   0, 1,   6,  6,  40 },  // inherited, no violation
                {  600,  0, 1,   0, 0,   0, 1,   6,  6,  40 },  // inherited, no violation
                {  700,  0, 0,   0, 0,   0, 1,   6,  6,  20 },  // de novo, under the default MIN_GQ
                {  800,  0, 0,   0, 0,   0, 1,  10,  1,  40 },  // de novo, too lopsided to judge
                {  900,  1, 1,   0, 0,   0, 1,   6,  6,  40 },  // inherited, no violation
                { 1000,  0, 0,   0, 0,   0, 0,   9,  0,  40 },  // not variant at all
        };

        try (htsjdk.variant.variantcontext.writer.VariantContextWriter writer =
                     new htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder()
                             .setOutputFile(f)
                             .setReferenceDictionary(dict)
                             .setOptions(java.util.EnumSet.of(
                                     htsjdk.variant.variantcontext.writer.Options.INDEX_ON_THE_FLY))
                             .build()) {
            writer.writeHeader(header);
            for (int[] row : rows) {
                int pos = row[0];
                htsjdk.variant.variantcontext.Allele ref =
                        htsjdk.variant.variantcontext.Allele.create(
                                chr1.substring(pos - 1, pos), true);
                String altBase = chr1.charAt(pos - 1) == 'A' ? "C" : "A";
                htsjdk.variant.variantcontext.Allele alt =
                        htsjdk.variant.variantcontext.Allele.create(altBase, false);
                java.util.List<htsjdk.variant.variantcontext.Allele> alleles =
                        java.util.Arrays.asList(ref, alt);

                java.util.List<htsjdk.variant.variantcontext.Genotype> genotypes =
                        new java.util.ArrayList<>();
                String[] names = {"father", "mother", "child"};
                for (int sample = 0; sample < 3; sample++) {
                    int first = row[1 + sample * 2];
                    int second = row[2 + sample * 2];
                    htsjdk.variant.variantcontext.GenotypeBuilder gb =
                            new htsjdk.variant.variantcontext.GenotypeBuilder(names[sample],
                                    java.util.Arrays.asList(
                                            first == 0 ? ref : alt, second == 0 ? ref : alt));
                    gb.GQ(sample == 2 ? row[9] : 40);
                    // Two sites are shallower than the others, so MIN_DP decides something short
                    // of deciding everything: at ten they drop out and the rest stay.
                    gb.DP(sample == 2 && (pos == 500 || pos == 900) ? 8 : 12);
                    gb.AD(sample == 2 ? new int[] {row[7], row[8]} : new int[] {6, 6});
                    // `MendelianViolationDetector.accumulate` reads PL without checking for it,
                    // so a genotype without one is a NullPointerException rather than a call the
                    // tool skips: every genotype here carries the three likelihoods of a diploid
                    // biallelic site, zero for the call that was made.
                    // The likelihoods carry the same contrast as the declared quality, because
                    // the detector reads the quality out of them: a PL that disagreed with GQ
                    // would make MIN_GQ decide nothing at all.
                    int called = first + second;
                    int contrast = sample == 2 ? row[9] : 40;
                    int[] pl = {called == 0 ? 0 : contrast,
                                called == 1 ? 0 : contrast,
                                called == 2 ? 0 : contrast};
                    gb.PL(pl);
                    genotypes.add(gb.make());
                }
                writer.add(new htsjdk.variant.variantcontext.VariantContextBuilder(
                        "fixture", "chr1", pos, pos, alleles).genotypes(genotypes).make());
            }
        }
    }

    static void writeVcf(File f, String chr1, String chr2, boolean withGenotypes) throws Exception {
        writeVcf(f, chr1, chr2, withGenotypes, 2);
    }

    static void writeVcf(File f, String chr1, String chr2, boolean withGenotypes, int sampleCount)
            throws Exception {
        htsjdk.samtools.SAMSequenceDictionary dict = new htsjdk.samtools.SAMSequenceDictionary();
        dict.addSequence(new SAMSequenceRecord("chr1", chr1.length()));
        dict.addSequence(new SAMSequenceRecord("chr2", chr2.length()));

        java.util.Set<htsjdk.variant.vcf.VCFHeaderLine> lines = new java.util.LinkedHashSet<>();
        lines.add(new htsjdk.variant.vcf.VCFFormatHeaderLine(
                "GT", 1, htsjdk.variant.vcf.VCFHeaderLineType.String, "Genotype"));
        lines.add(new htsjdk.variant.vcf.VCFFormatHeaderLine(
                "GQ", 1, htsjdk.variant.vcf.VCFHeaderLineType.Integer, "Genotype quality"));
        lines.add(new htsjdk.variant.vcf.VCFFormatHeaderLine(
                "DP", 1, htsjdk.variant.vcf.VCFHeaderLineType.Integer, "Depth"));
        lines.add(new htsjdk.variant.vcf.VCFInfoHeaderLine(
                "AC", 1, htsjdk.variant.vcf.VCFHeaderLineType.Integer, "Allele count"));
        lines.add(new htsjdk.variant.vcf.VCFFilterHeaderLine("LowQual", "Low quality"));

        java.util.List<String> samples = new java.util.ArrayList<>();
        if (withGenotypes) {
            for (int i = 1; i <= sampleCount; i++) samples.add("sample" + i);
        }
        htsjdk.variant.vcf.VCFHeader header = new htsjdk.variant.vcf.VCFHeader(lines, samples);
        header.setSequenceDictionary(dict);

        try (htsjdk.variant.variantcontext.writer.VariantContextWriter writer =
                     new htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder()
                             .setOutputFile(f)
                             .setReferenceDictionary(dict)
                             .setOption(htsjdk.variant.variantcontext.writer.Options.INDEX_ON_THE_FLY)
                             .build()) {
            writer.writeHeader(header);
            // The first four are on chr1 (2,000 bases), the last two on chr2 (1,000).
            int[] positions = {100, 300, 500, 700, 200, 600};
            for (int i = 0; i < positions.length; i++) {
                // The sites-only file keeps every other variant, so half of the full file is known.
                if (!withGenotypes && i % 2 == 1) continue;
                String contig = i < 4 ? "chr1" : "chr2";
                int position = positions[i];
                String reference = String.valueOf((i < 4 ? chr1 : chr2).charAt(position - 1));
                boolean indel = i == 3;
                htsjdk.variant.variantcontext.Allele ref = htsjdk.variant.variantcontext.Allele
                        .create(indel ? reference + "AT" : reference, true);
                htsjdk.variant.variantcontext.Allele alt = htsjdk.variant.variantcontext.Allele
                        .create(indel ? reference : (reference.equals("A") ? "G" : "A"), false);
                htsjdk.variant.variantcontext.VariantContextBuilder builder =
                        new htsjdk.variant.variantcontext.VariantContextBuilder()
                                .chr(contig)
                                .start(position)
                                .stop(position + ref.length() - 1)
                                .alleles(java.util.Arrays.asList(ref, alt))
                                .attribute("AC", 1 + (i % 2));
                if (i == 5) builder.filter("LowQual");
                if (withGenotypes) {
                    java.util.List<htsjdk.variant.variantcontext.Genotype> genotypes =
                            new java.util.ArrayList<>();
                    for (int g = 0; g < samples.size(); g++) {
                        java.util.List<htsjdk.variant.variantcontext.Allele> called = (i + g) % 3 == 0
                                ? java.util.Arrays.asList(ref, ref)
                                : ((i + g) % 3 == 1
                                        ? java.util.Arrays.asList(ref, alt)
                                        : java.util.Arrays.asList(alt, alt));
                        genotypes.add(new htsjdk.variant.variantcontext.GenotypeBuilder(
                                samples.get(g), called)
                                .GQ(20 + 7 * ((i + g) % 5))
                                .DP(10 + ((i + g) % 4))
                                .make());
                    }
                    builder.genotypes(genotypes);
                }
                writer.add(builder.make());
            }
        }
    }

    /**
     * The corpus of the picard.vcf manipulation tools.
     *
     * `vcf_sorted_samples.vcf` and `vcf_unsorted_samples.vcf` are the same nine records under the
     * two sample orders. The records are out of coordinate order and on both contigs, so SortVcf
     * has work to do; two of them are filtered, so VcfToIntervalList's INCLUDE_FILTERED decides
     * something; they overlap and abut, so its merge does; some have IDs (one of them two) and
     * some do not, so VARIANT_ID_METHOD does; and one is a symbolic deletion whose END is the
     * interval's end. `vcf_no_contigs.vcf` is the first file without its contig lines, which is
     * what the tools that need a sequence dictionary refuse. `vcf_one_sample.vcf` is the one
     * sample RenameSampleInVcf accepts and `vcf_sites_only.vcf` the none it also accepts.
     * The two `##source` lines are two lines to every tool but SortVcf, whose header goes
     * through `VCFUtils.smartMergeHeaders`: that keys an unstructured line by its key alone and
     * keeps the first in sorted order, so one of them is dropped.
     * `other.dict` disagrees with the corpus on chr2's length and adds a contig, and carries an
     * assembly, so a header rebuilt from it is visibly a different header.
     */
    static void writeVcfUtilityFixtures(File dir) throws Exception {
        String meta = String.join("\n",
                "##fileformat=VCFv4.2",
                "##FILTER=<ID=LowQual,Description=\"Low quality\">",
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
                "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype quality\">",
                "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">",
                "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count\">",
                "##INFO=<ID=DB,Number=0,Type=Flag,Description=\"dbSNP membership\">",
                "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position\">",
                "##ALT=<ID=DEL,Description=\"Deletion\">",
                "##source=handwritten",
                "##source=a second source line") + "\n";
        String contigs = "##contig=<ID=chr1,length=" + CHR1 + ">\n"
                + "##contig=<ID=chr2,length=" + CHR2 + ",assembly=test>\n";
        String[] records = {
            "chr2\t300\trs5\tG\tT\t12.5\tPASS\tAC=1\tGT:GQ:DP\t0/1:30:8\t0/0:.:.",
            "chr1\t700\t.\tC\tT\t.\tLowQual\tAC=2\tGT:GQ:DP\t1/1:12:3\t./.:.:.",
            "chr1\t100\trs1\tA\tG\t50\tPASS\tAC=1;DB\tGT:GQ:DP\t0/1:30:8\t0/0:.:.",
            "chr1\t101\t.\tT\tC\t40\t.\tAC=1\tGT:GQ:DP\t0/1:25:9\t0|1:20:7",
            "chr1\t200\trs2;rs3\tACGT\tA\t30\tPASS\tAC=1\tGT:GQ:DP\t0/1:.:4\t1/1:9:.",
            "chr1\t202\t.\tG\tA,C\t.\tPASS\tAC=1,1\tGT:GQ:DP\t1/2:5:5\t0/0:5:5",
            "chr1\t500\t.\tN\t<DEL>\t.\t.\tEND=560;AC=1\tGT\t0/1\t0/0",
            "chr1\t550\trs4\tA\tT\t.\t.\tAC=0\tGT:GQ:DP\t0/0:50:20\t0/0:40:20",
            "chr2\t100\t.\tA\tG\t.\tLowQual\tAC=1\tGT:GQ:DP\t0/1:30:8\t./.:.:.",
        };
        String columns = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";
        writeVcfText(new File(dir, "vcf_sorted_samples.vcf"), meta + contigs,
                columns + "\tFORMAT\tsampleA\tsampleB", records, 11);
        writeVcfText(new File(dir, "vcf_unsorted_samples.vcf"), meta + contigs,
                columns + "\tFORMAT\tsampleB\tsampleA", records, 11);
        writeVcfText(new File(dir, "vcf_no_contigs.vcf"), meta,
                columns + "\tFORMAT\tsampleA\tsampleB", records, 11);
        writeVcfText(new File(dir, "vcf_one_sample.vcf"), meta + contigs,
                columns + "\tFORMAT\tsampleA", records, 10);
        writeVcfText(new File(dir, "vcf_sites_only.vcf"), meta + contigs, columns, records, 8);

        try (PrintWriter p = new PrintWriter(new File(dir, "other.dict"), "UTF-8")) {
            p.print("@HD\tVN:1.6\n");
            p.print("@SQ\tSN:chr1\tLN:" + CHR1 + "\tAS:other\tM5:0123456789abcdef0123456789abcdef\n");
            p.print("@SQ\tSN:chr2\tLN:" + (CHR2 + 500) + "\n");
            p.print("@SQ\tSN:chr3\tLN:300\tSP:test\n");
        }
    }

    /**
     * The corpus of MergeVcfs: two files whose records are each in coordinate order and meet at
     * three positions (chr1:100, chr1:550, chr2:100), which is where the PriorityQueue inside
     * htsjdk's MergingIterator decides the order and the input order does not.
     *
     * `merge_a.vcf` has its samples in sorted order, so its genotype blocks are copied; `merge_b.vcf`
     * has them reversed, so its are decoded and re-encoded in the sorted order. The two disagree on
     * XC's Number, which smartMergeHeaders promotes to `.` (AC would not do: it is a standard line,
     * which the reader repairs to Number=A before the merge sees it), and b adds an INFO line.
     * `merge_b_no_contigs.vcf` is b without contig lines, which needs SEQUENCE_DICTIONARY, and
     * `merge_swapped_contigs.vcf` declares the same contigs in the other order, which the
     * comparator refuses as incompatible.
     */
    static void writeVcfMergeFixtures(File dir) throws Exception {
        String common = String.join("\n",
                "##fileformat=VCFv4.2",
                "##FILTER=<ID=LowQual,Description=\"Low quality\">",
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
                "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype quality\">",
                "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">",
                "##INFO=<ID=DB,Number=0,Type=Flag,Description=\"dbSNP membership\">") + "\n";
        String metaA = common
                + "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count\">\n"
                + "##INFO=<ID=XC,Number=2,Type=Integer,Description=\"Two counts\">\n"
                + "##source=merge_a\n";
        String metaB = common
                + "##INFO=<ID=AC,Number=1,Type=Integer,Description=\"Allele count in b\">\n"
                + "##INFO=<ID=XB,Number=1,Type=Float,Description=\"Only in b\">\n"
                + "##INFO=<ID=XC,Number=3,Type=Integer,Description=\"Three counts\">\n"
                + "##source=merge_b\n";
        String contigs = "##contig=<ID=chr1,length=" + CHR1 + ">\n"
                + "##contig=<ID=chr2,length=" + CHR2 + ",assembly=test>\n";
        String swapped = "##contig=<ID=chr2,length=" + CHR2 + ",assembly=test>\n"
                + "##contig=<ID=chr1,length=" + CHR1 + ">\n";
        String columns = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
        String[] recordsA = {
            "chr1\t100\trs1\tA\tG\t50\tPASS\tAC=1;DB\tGT:GQ:DP\t0/1:30:8\t0/0:.:.",
            "chr1\t200\t.\tACGT\tA\t30\tPASS\tAC=1\tGT:GQ:DP\t0/1:.:4\t1/1:9:.",
            "chr1\t550\trs4\tA\tT\t.\t.\tAC=0\tGT:GQ:DP\t0/0:50:20\t0/0:40:20",
            "chr2\t100\t.\tA\tG\t.\tLowQual\tAC=1\tGT:GQ:DP\t0/1:30:8\t./.:.:.",
        };
        // Columns in b's own order, sampleB first.
        String[] recordsB = {
            "chr1\t100\t.\tA\tC\t20\tPASS\tAC=1;XB=0.5\tGT:DP:GQ\t0/0:.:.\t0/1:7:22",
            "chr1\t150\trs9\tC\tCT\t.\tPASS\tAC=2\tGT:GQ\t1/1:15\t0/1:.",
            "chr1\t550\t.\tA\tC\t15\tLowQual\tXB=2\tGT:GQ:DP\t0/1:5:.\t0/0:.:3",
            "chr1\t700\t.\tC\tT\t.\t.\tAC=1\tGT:GQ:DP\t./.:.:.\t0/1:12:3",
            "chr2\t50\t.\tT\tA\t8\tPASS\t.\tGT\t0|1\t1|0",
            "chr2\t100\trs7\tA\tT\t.\tPASS\tAC=1\tGT:GQ:DP\t0/1:10:10\t0/0:10:10",
        };
        writeVcfText(new File(dir, "merge_a.vcf"), metaA + contigs,
                columns + "\tsampleA\tsampleB", recordsA, 11);
        writeVcfText(new File(dir, "merge_b.vcf"), metaB + contigs,
                columns + "\tsampleB\tsampleA", recordsB, 11);
        writeVcfText(new File(dir, "merge_b_no_contigs.vcf"), metaB,
                columns + "\tsampleB\tsampleA", recordsB, 11);
        writeVcfText(new File(dir, "merge_swapped_contigs.vcf"), metaB + swapped,
                columns + "\tsampleB\tsampleA", recordsB, 11);
    }

    /**
     * The corpus of GatherVcfs: `gather_1.vcf` and `gather_2.vcf` are two consecutive stretches of
     * one call set, `gather_empty.vcf` has the header and no records, `gather_overlap.vcf` starts
     * after gather_1 starts but before it ends (the first-record check passes and the gather's own
     * check does not), and `gather_2_swapped.vcf` is gather_2 with its sample columns in the
     * other order, which GatherVcfs compares as a list and refuses beside the others. Alone it is
     * accepted, and its output keeps its own column order with every genotype re-encoded, where
     * the other files' sorted columns are copied.
     */
    static void writeVcfGatherFixtures(File dir) throws Exception {
        String meta = String.join("\n",
                "##fileformat=VCFv4.2",
                "##FILTER=<ID=LowQual,Description=\"Low quality\">",
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
                "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype quality\">",
                "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">",
                "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count\">",
                "##source=gather",
                "##source=a second source line",
                "##contig=<ID=chr1,length=" + CHR1 + ">",
                "##contig=<ID=chr2,length=" + CHR2 + ",assembly=test>") + "\n";
        String columns = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
        String sorted = columns + "\tsampleA\tsampleB";
        String unsorted = columns + "\tsampleB\tsampleA";
        String[] first = {
            "chr1\t100\trs1\tA\tG\t50\tPASS\tAC=1\tGT:GQ:DP\t0/1:30:8\t0/0:.:.",
            "chr1\t200\t.\tACGT\tA\t30\tPASS\tAC=1\tGT:GQ:DP\t0/1:.:4\t1/1:9:.",
            "chr1\t550\trs4\tA\tT\t.\t.\tAC=0\tGT:DP:GQ\t0/0:20:50\t0/0:20:.",
        };
        String[] second = {
            "chr1\t700\t.\tC\tT\t.\tLowQual\tAC=1\tGT:GQ:DP\t./.:.:.\t0/1:12:3",
            "chr2\t50\t.\tT\tA\t8\tPASS\t.\tGT\t0|1\t1|0",
            "chr2\t100\trs7\tA\tT\t.\tPASS\tAC=1\tGT:GQ:DP\t0/1:10:10\t0/0:10:.",
        };
        String[] overlap = {
            "chr1\t300\t.\tG\tC\t9\tPASS\tAC=1\tGT:GQ\t0/1:9\t0/0:.",
            "chr1\t800\t.\tT\tG\t9\tPASS\tAC=1\tGT:GQ\t0/1:9\t0/0:.",
        };
        writeVcfText(new File(dir, "gather_1.vcf"), meta, sorted, first, 11);
        writeVcfText(new File(dir, "gather_2.vcf"), meta, sorted, second, 11);
        writeVcfText(new File(dir, "gather_empty.vcf"), meta, sorted, new String[0], 11);
        writeVcfText(new File(dir, "gather_overlap.vcf"), meta, sorted, overlap, 11);
        writeVcfText(new File(dir, "gather_2_swapped.vcf"), meta, unsorted, second, 11);
    }

    /**
     * The corpus of SplitVcfs beside vcf_sorted_samples.vcf (whose one odd record is SYMBOLIC):
     * every VariantContext type in one file, in coordinate order. The first record that is
     * neither a SNP nor an indel is a MIXED site, which is the type STRICT names; the multiallelic
     * SNP, the MNP, the site with no ALT (NO_VARIATION) and the insertion come after it.
     */
    static void writeVcfSplitFixture(File dir) throws Exception {
        String meta = String.join("\n",
                "##fileformat=VCFv4.2",
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
                "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">",
                "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total depth\">",
                "##contig=<ID=chr1,length=" + CHR1 + ">",
                "##contig=<ID=chr2,length=" + CHR2 + ">") + "\n";
        String[] records = {
            "chr1\t10\t.\tA\tG\t.\t.\tDP=5\tGT:DP\t0/1:5",
            "chr1\t20\t.\tAC\tA\t.\t.\t.\tGT:DP\t1/1:.",
            "chr1\t30\t.\tA\tG,AT\t.\t.\t.\tGT\t1/2",
            "chr1\t40\t.\tC\tA,T\t9\tPASS\tDP=3\tGT:DP\t1/2:3",
            "chr1\t50\t.\tAC\tGT\t.\t.\t.\tGT\t0/1",
            "chr1\t60\t.\tT\t.\t.\t.\t.\tGT\t0/0",
            "chr2\t5\trs5\tG\tGTT\t.\tPASS\t.\tGT:DP\t0|1:7",
            "chr2\t9\t.\tC\tT\t.\t.\t.\tGT:DP\t./.:.",
        };
        writeVcfText(new File(dir, "split_types.vcf"), meta,
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsampleA", records, 10);
    }

    /**
     * The corpus of VcfFormatConverter beside variants.vcf (which has its .idx): gather_1.vcf as
     * a block-compressed file with its tabix index, written by htsjdk's own writer, so the one
     * compressed input the tool reads with REQUIRE_INDEX has the index it looks for; and
     * `unsorted_indexed.vcf`, the out-of-order records beside a borrowed `.idx`.
     */
    static void writeVcfConverterFixture(File dir) throws Exception {
        // The unsorted records under an index the reader only needs to find and load, so
        // REQUIRE_INDEX lets them through to a compressed output, whose tabix index refuses them.
        java.nio.file.Files.copy(new File(dir, "vcf_unsorted_samples.vcf").toPath(),
                new File(dir, "unsorted_indexed.vcf").toPath());
        java.nio.file.Files.copy(new File(dir, "variants.vcf.idx").toPath(),
                new File(dir, "unsorted_indexed.vcf.idx").toPath());
        try (htsjdk.variant.vcf.VCFFileReader in =
                     new htsjdk.variant.vcf.VCFFileReader(new File(dir, "gather_1.vcf"), false);
             htsjdk.variant.variantcontext.writer.VariantContextWriter out =
                     new htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder()
                             .setOutputFile(new File(dir, "gather_1.vcf.gz"))
                             .setReferenceDictionary(in.getFileHeader().getSequenceDictionary())
                             .setOption(htsjdk.variant.variantcontext.writer.Options.INDEX_ON_THE_FLY)
                             .build()) {
            out.writeHeader(in.getFileHeader());
            for (htsjdk.variant.variantcontext.VariantContext vc : in) {
                out.add(vc);
            }
        }
    }

    /**
     * The corpus of FixVcfHeader. `fix_missing.vcf` uses a FILTER, three INFO keys (one a flag)
     * and three FORMAT keys (one of them GQ, a standard key) that its header does not define, at
     * most one undefined INFO key per record, so which one the writer names first never depends on
     * a HashMap's order. `fix_missing_unsorted.vcf` is the same under the other column order, so
     * its genotypes are decoded and their FORMAT keys checked. `fix_header.vcf` is a header
     * defining all of them, with no records, for HEADER, and `fix_header_swapped.vcf` the same
     * with its sample columns in the other order, which ENFORCE_SAME_SAMPLES lets through.
     */
    static void writeVcfFixHeaderFixtures(File dir) throws Exception {
        String contigs = "##contig=<ID=chr1,length=" + CHR1 + ">\n"
                + "##contig=<ID=chr2,length=" + CHR2 + ">\n";
        String meta = String.join("\n",
                "##fileformat=VCFv4.2",
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
                "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">",
                "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total depth\">",
                "##source=fix") + "\n" + contigs;
        String full = String.join("\n",
                "##fileformat=VCFv4.2",
                "##FILTER=<ID=q10,Description=\"Quality below 10\">",
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
                "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype quality\">",
                "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">",
                "##FORMAT=<ID=XF,Number=1,Type=String,Description=\"A string\">",
                "##FORMAT=<ID=XG,Number=1,Type=Integer,Description=\"An integer\">",
                "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total depth\">",
                "##INFO=<ID=XA,Number=1,Type=Integer,Description=\"One\">",
                "##INFO=<ID=XB,Number=2,Type=Integer,Description=\"Two\">",
                "##INFO=<ID=FLAGX,Number=0,Type=Flag,Description=\"A flag\">",
                "##source=the replacement header") + "\n" + contigs;
        String columns = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
        String[] records = {
            "chr1\t100\t.\tA\tG\t50\tPASS\tDP=10;XA=1\tGT:GQ:XF\t0/1:30:a\t0/0:.:.",
            "chr1\t200\t.\tC\tT\t.\tq10\tDP=5;FLAGX\tGT:DP\t1/1:3\t0/1:2",
            "chr2\t50\t.\tG\tA\t.\tPASS\tXB=2,3\tGT:XG\t0|1:7\t1|1:8",
        };
        // The same records with their two sample columns swapped, for the other column order.
        String[] swapped = new String[records.length];
        for (int i = 0; i < records.length; i++) {
            String[] f = records[i].split("\t");
            String t = f[9];
            f[9] = f[10];
            f[10] = t;
            swapped[i] = String.join("\t", f);
        }
        writeVcfText(new File(dir, "fix_missing.vcf"), meta,
                columns + "\tsampleA\tsampleB", records, 11);
        writeVcfText(new File(dir, "fix_missing_unsorted.vcf"), meta,
                columns + "\tsampleB\tsampleA", swapped, 11);
        writeVcfText(new File(dir, "fix_header.vcf"), full,
                columns + "\tsampleA\tsampleB", new String[0], 11);
        writeVcfText(new File(dir, "fix_header_swapped.vcf"), full,
                columns + "\tsampleB\tsampleA", new String[0], 11);
    }

    /** A VCF as text: the meta lines, the column line, and each record cut to its first columns. */
    static void writeVcfText(File f, String meta, String columns, String[] records, int keep)
            throws Exception {
        try (PrintWriter p = new PrintWriter(f, "UTF-8")) {
            p.print(meta);
            p.print(columns + "\n");
            for (String record : records) {
                String[] fields = record.split("\t");
                p.print(String.join("\t", java.util.Arrays.copyOf(fields, keep)) + "\n");
            }
        }
    }

    static void writeDict(File f, String chr1, String chr2) throws Exception {
        try (PrintWriter p = new PrintWriter(f)) {
            p.println("@HD\tVN:1.6\tSO:unsorted");
            p.printf("@SQ\tSN:chr1\tLN:%d%n", chr1.length());
            p.printf("@SQ\tSN:chr2\tLN:%d%n", chr2.length());
        }
    }

    /** A contig's bytes in the file: its bases, plus the newline that ends each line. */
    static long bytesOnDisk(int bases, int lineWidth) {
        long lines = (bases + lineWidth - 1) / lineWidth;
        return bases + lines;
    }

    static void writeFastq(File f, int end) throws Exception {
        Random rng = new Random(20260800L + end);
        try (PrintWriter p = new PrintWriter(f)) {
            for (int i = 0; i < 40; i++) {
                StringBuilder bases = new StringBuilder();
                StringBuilder quals = new StringBuilder();
                for (int b = 0; b < READ_LENGTH; b++) {
                    bases.append("ACGT".charAt(rng.nextInt(4)));
                    quals.append((char) (33 + 2 + rng.nextInt(38)));
                }
                p.printf("@fq%04d/%d%n%s%n+%n%s%n", i, end, bases, quals);
            }
        }
    }

    /** The pre-adapter and bait-bias detail tables of one artifact run; see the call site. */
    static void writeArtifactTables(File dir, String base, String[] libraries, String sample,
                                    String[] keep) throws Exception {
        Random rng = new Random(20260930L + base.length());
        htsjdk.samtools.metrics.MetricsFile<picard.analysis.artifacts.SequencingArtifactMetrics.PreAdapterDetailMetrics, Integer> pre =
                new htsjdk.samtools.metrics.MetricsFile<>();
        htsjdk.samtools.metrics.MetricsFile<picard.analysis.artifacts.SequencingArtifactMetrics.BaitBiasDetailMetrics, Integer> bait =
                new htsjdk.samtools.metrics.MetricsFile<>();
        String bases = "ACGT";
        java.util.List<String> contexts = new java.util.ArrayList<>();
        for (char a : bases.toCharArray())
            for (char b : bases.toCharArray())
                for (char c : bases.toCharArray()) contexts.add("" + a + b + c);
        for (String library : libraries) {
            for (char ref : bases.toCharArray()) {
                for (char alt : bases.toCharArray()) {
                    if (ref == alt) continue;
                    for (int i = 0; i < contexts.size(); i++) {
                        String context = contexts.get(i);
                        if (context.charAt(1) != ref) continue;
                        if (keep != null && !java.util.Arrays.asList(keep).contains(context)) continue;
                        boolean empty = context.equals("CCG") || context.equals("CGG") || context.equals("GCA");
                        picard.analysis.artifacts.SequencingArtifactMetrics.PreAdapterDetailMetrics p =
                                new picard.analysis.artifacts.SequencingArtifactMetrics.PreAdapterDetailMetrics();
                        p.SAMPLE_ALIAS = sample;
                        p.LIBRARY = library;
                        p.REF_BASE = ref;
                        p.ALT_BASE = alt;
                        p.CONTEXT = context;
                        p.PRO_REF_BASES = empty ? 0 : 50 + rng.nextInt(200);
                        p.PRO_ALT_BASES = empty ? 0 : rng.nextInt(30);
                        p.CON_REF_BASES = empty ? 0 : 50 + rng.nextInt(200);
                        p.CON_ALT_BASES = empty ? 0 : rng.nextInt(30);
                        p.calculateDerivedStatistics();
                        pre.addMetric(p);
                        picard.analysis.artifacts.SequencingArtifactMetrics.BaitBiasDetailMetrics b =
                                new picard.analysis.artifacts.SequencingArtifactMetrics.BaitBiasDetailMetrics();
                        b.SAMPLE_ALIAS = sample;
                        b.LIBRARY = library;
                        b.REF_BASE = ref;
                        b.ALT_BASE = alt;
                        b.CONTEXT = context;
                        b.FWD_CXT_REF_BASES = empty ? 0 : 50 + rng.nextInt(200);
                        b.FWD_CXT_ALT_BASES = empty ? 0 : rng.nextInt(30);
                        b.REV_CXT_REF_BASES = empty ? 0 : 50 + rng.nextInt(200);
                        b.REV_CXT_ALT_BASES = empty ? 0 : rng.nextInt(30);
                        b.calculateDerivedStatistics();
                        bait.addMetric(b);
                    }
                }
            }
        }
        pre.write(new File(dir, base + ".pre_adapter_detail_metrics"));
        bait.write(new File(dir, base + ".bait_bias_detail_metrics"));
    }

    /** The bisulfite-converted reads of `rrbs.bam`; see the call site. */
    static java.util.List<SAMRecord> rrbsReads(SAMFileHeader header, String chr1, String chr2) {
        Random rng = new Random(20261001L);
        java.util.List<SAMRecord> out = new java.util.ArrayList<>();
        for (int i = 0; i < 240; i++) {
            boolean onChr2 = i % 5 == 0;
            String contig = onChr2 ? chr2 : chr1;
            int length = i % 10 == 9 ? 4 : READ_LENGTH;
            // Clear of the run of Ns near each contig's end.
            int start = 1 + rng.nextInt(contig.length() - 80 - length);
            boolean negative = i % 3 == 1;
            byte[] bases = contig.substring(start - 1, start - 1 + length).getBytes();
            byte[] quals = new byte[length];
            for (int b = 0; b < length; b++) {
                quals[b] = (byte) (8 + rng.nextInt(33));
                int pos = start - 1 + b;
                if (!negative && bases[b] == 'C') {
                    boolean cpg = pos + 1 < contig.length() && contig.charAt(pos + 1) == 'G';
                    boolean convert = cpg ? rng.nextInt(4) == 0 : rng.nextInt(20) != 0;
                    if (convert) bases[b] = 'T';
                } else if (negative && bases[b] == 'G') {
                    boolean cpg = pos > 0 && contig.charAt(pos - 1) == 'C';
                    boolean convert = cpg ? rng.nextInt(4) == 0 : rng.nextInt(20) != 0;
                    if (convert) bases[b] = 'A';
                }
            }
            if (i % 8 == 3 && length == READ_LENGTH) {
                for (int e = 0; e < 8; e++) {
                    int b = 5 * e + 2;
                    bases[b] = mutateBase(bases[b]);
                }
            }
            SAMRecord r = new SAMRecord(header);
            r.setReadName(String.format("rrbs%04d", i));
            r.setReadBases(bases);
            r.setBaseQualities(quals);
            r.setReferenceIndex(onChr2 ? 1 : 0);
            r.setAlignmentStart(start);
            r.setCigarString(length + "M");
            r.setMappingQuality(60);
            r.setReadNegativeStrandFlag(negative);
            r.setAttribute("RG", i % 2 == 0 ? "rg1" : "rg2");
            if (i % 30 == 7) {
                // Placed but unmapped, which the tool skips before it asks the reference.
                r.setReadUnmappedFlag(true);
                r.setMappingQuality(0);
                r.setCigarString("*");
            }
            out.add(r);
        }
        out.sort(new SAMRecordCoordinateComparator());
        return out;
    }

    /**
     * `SortGff`'s inputs. `sort.gff3` has its features out of order across five contigs, two of
     * which no dictionary in the corpus names and one (chr3) only other.dict does; children before
     * their parents, a feature with two parents, an ID shared by two CDS lines far apart, and two
     * features starting together (the sort is stable). Its attribute values need escaping on the
     * way out (`:`, `,`, `;`, a space), its scores are absent, whole, fractional, negative and in
     * exponent form, and one strand is `?`. Comments sit before and between features, an unknown
     * directive is ignored, a flush directive is in the input, and a FASTA section ends it. The
     * other two are refusals: a feature outside its sequence-region, and a file with directives
     * and comments but no feature, which `canDecode` rejects.
     */
    static void writeGffFixtures(File dir) throws Exception {
        String t = "\t";
        try (PrintWriter p = new PrintWriter(new File(dir, "sort.gff3"), "UTF-8")) {
            p.print("##gff-version 3\n");
            p.print("# annotation for SortGff\n");
            p.print("##sequence-region chr2 1 1000\n");
            p.print("##sequence-region chr1 1 2000\n");
            p.print("#!genome-build test\n");
            String[] lines = {
                "chr2|src|gene|500|900|.|+|.|ID=gene:B;Name=geneB",
                "chr2|src|mRNA|500|900|.|+|.|ID=tx:B1;Parent=gene:B",
                "chr2|src|exon|700|900|12|+|.|ID=exon:B2;Parent=tx:B1",
                "chr2|src|exon|500|600|0.5|+|.|ID=exon:B1;Parent=tx:B1",
                "chr1|src|CDS|1500|1600|.|-|0|ID=cds:A;Parent=tx:A1,tx:A2",
                "chr1|src|mRNA|1000|1800|.|-|.|ID=tx:A1;Parent=gene:A",
                "chr1|src|gene|1000|1900|-3|-|.|ID=gene:A;Name=gene A;Note=a%2Cb%3Bc",
                "###",
                "chr10|other|region|1|300|1e8|.|.|.",
                "chr1|src|mRNA|1000|1700|.|-|.|ID=tx:A2;Parent=gene:A",
                "# a comment between features",
                "chr3|src|repeat|50|80|.|?|.|Alias=r1,r2;note= spaced value ",
                "chrM|src|gene|10|20|0|+|.|ID=gene:M",
                "chr1|src|exon|100|200|.|+|.|Parent=tx:Z",
                "chr1|src|gene|100|300|.|+|.|ID=gene:Z",
                "chr1|src|mRNA|100|250|.|+|.|ID=tx:Z;Parent=gene:Z",
                "chr1|src|CDS|120|180|.|+|2|ID=cds:Z;Parent=tx:Z",
                "chr1|src|CDS|400|500|.|+|1|ID=cds:Z;Parent=tx:Z",
                "##species https://example.org/species",
                "chr1|src|exon|260|290|.|+|.|Parent=tx:Z",
                "chr2|src|gene|100|200|.|+|.|ID=gene:C",
                "##FASTA",
                ">chr1",
                "ACGTACGT",
            };
            for (String line : lines) {
                p.print(line.replace("|", t) + "\n");
            }
        }
        try (PrintWriter p = new PrintWriter(new File(dir, "gff_bad_region.gff3"), "UTF-8")) {
            p.print("##gff-version 3\n##sequence-region chr1 1 100\n");
            p.print("chr1" + t + "src" + t + "gene" + t + "10" + t + "50" + t + "." + t + "+" + t + "." + t + "ID=g1\n");
            p.print("chr1" + t + "src" + t + "gene" + t + "50" + t + "150" + t + "." + t + "+" + t + "." + t + "ID=g2\n");
        }
        // A circular contig: its landmark feature spans the whole sequence-region and says
        // Is_circular, after which a feature that runs off the region's end only has to overlap
        // it. The version directive is not the writer's own, which the output does not keep.
        try (PrintWriter p = new PrintWriter(new File(dir, "gff_circular.gff3"), "UTF-8")) {
            p.print("##gff-version 3.1.26\n##sequence-region chrC 1 100\n");
            String[] lines = {
                "chrC|src|region|1|100|.|+|.|ID=chrC;Is_circular=true",
                "chrC|src|gene|90|120|.|+|.|ID=wrap;Name=wrapping gene",
                "chrC|src|gene|5|20|7.25|-|.|ID=early",
                "chrB|src|gene|5|20|.|.|.|ID=other;Dbxref=db:1,db:2",
            };
            for (String line : lines) {
                p.print(line.replace("|", t) + "\n");
            }
        }
        try (PrintWriter p = new PrintWriter(new File(dir, "gff_no_features.gff3"), "UTF-8")) {
            p.print("##gff-version 3\n# no features here\n##sequence-region chr1 1 2000\n");
        }
    }

    /**
     * `MergePedIntoVcf`'s inputs. `zcall.vcf` is one sample in the shape GtcToVcf writes: ALLELE_A
     * and ALLELE_B name the two array alleles, the reference one starred, and the genotypes carry
     * IGC, X and Y beside GT (one also GQ, which the merged genotype drops). Its records cover a
     * starred A and a starred B, a no-call, an indel, a missing IGC and a site with no ALT.
     * `zcall.ped` calls every SNP of `zcall.map` in order; `zcall_illegal.ped` calls rs5 with a
     * letter that is neither A, B nor 0 and `zcall_two_lines.ped` is two samples; `zcall_partial.map`
     * leaves rs3 out, so the PED pairs shift and rs3 has no call. The thresholds name some SNPs,
     * NA for both of one, and the second file has one NA alone, which the tool refuses.
     */
    static void writeZCallFixtures(File dir) throws Exception {
        String meta = String.join("\n",
                "##fileformat=VCFv4.2",
                "##FILTER=<ID=LowQual,Description=\"Low quality\">",
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
                "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype quality\">",
                "##FORMAT=<ID=IGC,Number=1,Type=Float,Description=\"Illumina GenCall Confidence Score\">",
                "##FORMAT=<ID=X,Number=1,Type=Integer,Description=\"Raw X intensity\">",
                "##FORMAT=<ID=Y,Number=1,Type=Integer,Description=\"Raw Y intensity\">",
                "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count\">",
                "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency\">",
                "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Allele number\">",
                "##INFO=<ID=ALLELE_A,Number=1,Type=String,Description=\"A allele\">",
                "##INFO=<ID=ALLELE_B,Number=1,Type=String,Description=\"B allele\">",
                "##autocallVersion=3.0.0",
                "##contig=<ID=chr1,length=" + CHR1 + ">",
                "##contig=<ID=chr2,length=" + CHR2 + ">") + "\n";
        String[] records = {
            "chr1|100|rs1|C|T|.|PASS|AC=1;AF=0.500;AN=2;ALLELE_A=C*;ALLELE_B=T|GT:IGC:X:Y|0/1:0.8100:1000:950",
            "chr1|200|rs2|G|A|.|PASS|AC=2;AF=1.00;AN=2;ALLELE_A=A;ALLELE_B=G*|GT:IGC:X:Y|1/1:0.9000:200:1500",
            "chr1|300|rs3|T|C|.|PASS|AC=0;AF=0.00;AN=2;ALLELE_A=T*;ALLELE_B=C|GT:IGC:X:Y|0/0:0.7000:1800:100",
            "chr1|900|rs4|A|G|50|LowQual|AC=0;AF=0.00;AN=0;ALLELE_A=A*;ALLELE_B=G|GT:IGC:X:Y|./.:.:300:310",
            "chr2|50|rs5|AT|A|.|PASS|AC=1;AF=0.500;AN=2;ALLELE_A=AT*;ALLELE_B=A|GT:IGC:X:Y|0/1:0.5500:700:720",
            "chr2|400|rs6|C|T|.|PASS|AC=1;AF=0.500;AN=2;ALLELE_A=C*;ALLELE_B=T|GT:GQ:IGC:X:Y|0/1:30:0.6000:640:660",
            "chr2|600|rs7|G|.|.|PASS|AN=2;ALLELE_A=G*;ALLELE_B=G|GT:IGC:X:Y|0/0:0.9500:900:50",
        };
        try (PrintWriter p = new PrintWriter(new File(dir, "zcall.vcf"), "UTF-8")) {
            p.print(meta);
            p.print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1\n");
            for (String record : records) {
                p.print(record.replace("|", "\t") + "\n");
            }
        }
        String[] calls = {"A\tB", "B\tB", "0\t0", "A\tA", "A\tB", "B\tA", "A\tA"};
        String prefix = "FAM1\tSAMPLE1\t0\t0\t2\t-9";
        try (PrintWriter p = new PrintWriter(new File(dir, "zcall.ped"), "UTF-8")) {
            p.print(prefix + "\t" + String.join("\t", calls) + "\n");
        }
        String[] illegal = calls.clone();
        illegal[4] = "A\tC";
        try (PrintWriter p = new PrintWriter(new File(dir, "zcall_illegal.ped"), "UTF-8")) {
            p.print(prefix + "\t" + String.join("\t", illegal) + "\n");
        }
        try (PrintWriter p = new PrintWriter(new File(dir, "zcall_two_lines.ped"), "UTF-8")) {
            p.print(prefix + "\t" + String.join("\t", calls) + "\n");
            p.print("FAM1\tSAMPLE2\t0\t0\t1\t-9\t" + String.join("\t", calls) + "\n");
        }
        int[] positions = {100, 200, 300, 900, 50, 400, 600};
        try (PrintWriter all = new PrintWriter(new File(dir, "zcall.map"), "UTF-8");
             PrintWriter partial = new PrintWriter(new File(dir, "zcall_partial.map"), "UTF-8")) {
            for (int i = 0; i < positions.length; i++) {
                String line = (i < 4 ? "1" : "2") + "\trs" + (i + 1) + "\t0\t" + positions[i] + "\n";
                all.print(line);
                if (i != 2) partial.print(line);
            }
        }
        try (PrintWriter p = new PrintWriter(new File(dir, "zcall_thresholds.txt"), "UTF-8")) {
            p.print("Name\tThr_X\tThr_Y\nrs1\t1.5\t2.5\nrs2\tNA\tNA\nrs3\t0.25\t3.0\nrs6\t10\t20\nrs7\t1e-3\t0\n");
        }
        try (PrintWriter p = new PrintWriter(new File(dir, "zcall_thresholds_half_na.txt"), "UTF-8")) {
            p.print("rs1\t1.5\tNA\nrs2\t0.5\t0.5\n");
        }
    }

    /**
     * `BamToBfq`'s reads: sixteen unaligned pairs, adjacent and in name order, most named under a
     * `RUN1:` prefix the tool can strip. Pair 1 is noise on both reads (XN=1) and pair 2 on one;
     * pair 3 carries XN=2, which the writer keeps and READS_TO_ALIGN's count does not; one read of
     * pair 4 fails vendor QC; pair 5 is wholly clipped (XT=1) on one read; pairs 6 and 7 mark an
     * adapter at XT=40 and XT=10; pair 8 has four no-calls, three of them in the seed region; pair
     * 9 has qualities above 63; pair 11 is thirty bases long; pair 15 is under another prefix. `odd` leaves pair 13 a single read and keeps the header
     * unsorted, which the paired path and READS_TO_ALIGN's sort check both refuse.
     */
    static java.util.List<SAMRecord> bfqReads(SAMFileHeader header, boolean odd) {
        Random rng = new Random(20261002L);
        java.util.List<SAMRecord> out = new java.util.ArrayList<>();
        for (int i = 0; i < 16; i++) {
            String name = (i == 15 ? "RUN2:" : "RUN1:") + String.format("frag%02d", i);
            int length = i == 11 ? 30 : READ_LENGTH;
            SAMRecord[] pair = new SAMRecord[2];
            for (int end = 0; end < 2; end++) {
                byte[] bases = new byte[length];
                byte[] quals = new byte[length];
                for (int b = 0; b < length; b++) {
                    bases[b] = (byte) "ACGT".charAt(rng.nextInt(4));
                    quals[b] = (byte) (i == 9 ? 50 + rng.nextInt(44) : 2 + rng.nextInt(39));
                }
                if (i == 8) {
                    for (int b : new int[] {0, 5, 10, 35}) bases[b] = 'N';
                }
                SAMRecord r = new SAMRecord(header);
                r.setReadName(name);
                r.setReadBases(bases);
                r.setBaseQualities(quals);
                r.setReadUnmappedFlag(true);
                r.setReferenceIndex(SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX);
                r.setAlignmentStart(SAMRecord.NO_ALIGNMENT_START);
                r.setMappingQuality(0);
                r.setReadPairedFlag(true);
                r.setMateUnmappedFlag(true);
                r.setFirstOfPairFlag(end == 0);
                r.setSecondOfPairFlag(end == 1);
                r.setAttribute("RG", "rg1");
                pair[end] = r;
            }
            if (i == 1) { pair[0].setAttribute("XN", 1); pair[1].setAttribute("XN", 1); }
            if (i == 2) pair[0].setAttribute("XN", 1);
            if (i == 3) { pair[0].setAttribute("XN", 2); pair[1].setAttribute("XN", 2); }
            if (i == 4) pair[1].setReadFailsVendorQualityCheckFlag(true);
            if (i == 5) pair[0].setAttribute("XT", 1);
            if (i == 6) { pair[0].setAttribute("XT", 40); pair[1].setAttribute("XT", 40); }
            if (i == 7) pair[1].setAttribute("XT", 10);
            if (odd && i == 13) {
                out.add(pair[0]);
            } else {
                out.add(pair[0]);
                out.add(pair[1]);
            }
        }
        return out;
    }

    /**
     * `CompareMetrics`' inputs: InsertSizeMetrics tables of three rows (all reads, a sample, a read
     * group) with a two-column histogram, written through Picard's own bean and `MetricsFile`
     * with no header. `ism_b` changes one value of each kind (a double by 1%, a long by 0.5%, an
     * int by 20%, the enum) and one histogram bin; `ism_c` is `ism_a` with its rows in another
     * order; `ism_d` lacks the WIDTH_OF_99_PERCENT column; `ism_e` keeps only the first row. The
     * changes are ratios a double holds exactly enough that `Double.toString` prints them short.
     */
    static void writeCompareMetricsFixtures(File dir) throws Exception {
        java.util.List<picard.analysis.InsertSizeMetrics> a = insertSizeRows(false);
        java.util.List<picard.analysis.InsertSizeMetrics> b = insertSizeRows(true);
        writeInsertSizeMetrics(new File(dir, "ism_a.txt"), a, false, null);
        writeInsertSizeMetrics(new File(dir, "ism_b.txt"), b, true, null);
        java.util.List<picard.analysis.InsertSizeMetrics> c = new java.util.ArrayList<>();
        c.add(a.get(2));
        c.add(a.get(0));
        c.add(a.get(1));
        writeInsertSizeMetrics(new File(dir, "ism_c.txt"), c, false, null);
        writeInsertSizeMetrics(new File(dir, "ism_d.txt"), a, false, "WIDTH_OF_99_PERCENT");
        writeInsertSizeMetrics(new File(dir, "ism_e.txt"), a.subList(0, 1), false, null);
    }

    static java.util.List<picard.analysis.InsertSizeMetrics> insertSizeRows(boolean changed) {
        java.util.List<picard.analysis.InsertSizeMetrics> rows = new java.util.ArrayList<>();
        String[][] levels = {{null, null, null}, {"sample1", null, null}, {"sample1", "lib1", "rg1"}};
        double[] means = {200.0, 210.0, 205.0};
        long[] pairs = {1000, 600, 400};
        for (int i = 0; i < 3; i++) {
            picard.analysis.InsertSizeMetrics m = new picard.analysis.InsertSizeMetrics();
            m.SAMPLE = levels[i][0];
            m.LIBRARY = levels[i][1];
            m.READ_GROUP = levels[i][2];
            m.MEDIAN_INSERT_SIZE = 200.0 + i;
            m.MODE_INSERT_SIZE = 198.0;
            m.MEDIAN_ABSOLUTE_DEVIATION = 20.0;
            m.MIN_INSERT_SIZE = 50;
            m.MAX_INSERT_SIZE = 400;
            m.MEAN_INSERT_SIZE = (changed && i == 0) ? 202.0 : means[i];
            m.STANDARD_DEVIATION = 30.5;
            m.READ_PAIRS = (changed && i == 1) ? 603 : pairs[i];
            m.PAIR_ORIENTATION = (changed && i == 2) ? SamPairUtil.PairOrientation.RF
                    : SamPairUtil.PairOrientation.FR;
            m.WIDTH_OF_10_PERCENT = (changed && i == 2) ? 12 : 10;
            m.WIDTH_OF_20_PERCENT = 20;
            m.WIDTH_OF_30_PERCENT = 30;
            m.WIDTH_OF_40_PERCENT = 40;
            m.WIDTH_OF_50_PERCENT = 50;
            m.WIDTH_OF_60_PERCENT = 60;
            m.WIDTH_OF_70_PERCENT = 70;
            m.WIDTH_OF_80_PERCENT = 80;
            m.WIDTH_OF_90_PERCENT = 90;
            m.WIDTH_OF_95_PERCENT = 95;
            m.WIDTH_OF_99_PERCENT = 99;
            rows.add(m);
        }
        return rows;
    }

    static void writeInsertSizeMetrics(File f, java.util.List<picard.analysis.InsertSizeMetrics> rows,
                                       boolean changed, String dropColumn) throws Exception {
        htsjdk.samtools.metrics.MetricsFile<picard.analysis.InsertSizeMetrics, Integer> file =
                new htsjdk.samtools.metrics.MetricsFile<>();
        for (picard.analysis.InsertSizeMetrics m : rows) file.addMetric(m);
        htsjdk.samtools.util.Histogram<Integer> all = new htsjdk.samtools.util.Histogram<>("insert_size", "All_Reads.fr_count");
        htsjdk.samtools.util.Histogram<Integer> sample = new htsjdk.samtools.util.Histogram<>("insert_size", "sample1.fr_count");
        for (int size = 100; size <= 300; size += 50) {
            all.increment(size, size == 150 && changed ? 41 : size / 5);
            sample.increment(size, size / 10);
        }
        file.addHistogram(all);
        file.addHistogram(sample);
        java.io.StringWriter text = new java.io.StringWriter();
        file.write(text);
        String out = text.toString();
        if (dropColumn != null) {
            StringBuilder kept = new StringBuilder();
            int drop = -1;
            boolean inTable = false;
            for (String line : out.split("\n", -1)) {
                if (line.startsWith("## METRICS CLASS")) {
                    inTable = true;
                } else if (inTable && line.isEmpty()) {
                    inTable = false;
                } else if (inTable) {
                    java.util.List<String> cells = new java.util.ArrayList<>(java.util.Arrays.asList(line.split("\t", -1)));
                    if (drop < 0) drop = cells.indexOf(dropColumn);
                    cells.remove(drop);
                    line = String.join("\t", cells);
                }
                kept.append(line).append("\n");
            }
            out = kept.substring(0, kept.length() - 1);
        }
        try (PrintWriter p = new PrintWriter(f, "UTF-8")) {
            p.print(out);
        }
    }

    /**
     * `MergeBamAlignment`'s inputs: `mba_unmapped.bam`, queryname-sorted reads as they came off the
     * sequencer (four fragments, then seventeen pairs), and `mba_aligned.bam`, what an aligner made
     * of them. Each template is one case of the merge: a plain proper pair (t00), an overlapping
     * one the merge clips (t01), a pair with one end unaligned (t02) and with neither aligned
     * (t03), two hits with the aligner's primary marked by HI (t04) and two with none marked
     * (t05), a supplementary alignment (t06), an end that runs off chr2 (t07), two indels against
     * MAX_INSERTIONS_OR_DELETIONS (t08), aligner hard clips (t09), a short doubly clipped
     * alignment the contamination filter takes (t10), adapter positions in XT (t11), an aligned
     * read shorter than the original (t12), aligner tags of every kind against the retain and
     * remove lists (t13), OQ and E2 on a negative-strand end (t14), ends on two contigs (t15) and
     * an RF pair (t16). The fragments have one hit, two unmarked hits, none, and a negative hit.
     * `mba_aligned_simple.bam` keeps one hit per paired end, has two @PG lines (so none is
     * adopted) and is out of queryname order (so the merge retries with the aligned reads sorted).
     * Read bases are the reference's at the true location with two substitutions, so NM, MD and
     * UQ have something to count.
     */
    static void writeMbaFixtures(File dir, String chr1, String chr2) throws Exception {
        String[] refs = {chr1, chr2};
        // name, end (0 fragment, 1 first, 2 second), contig, start, cigar, strand, mapq, flags, HI
        String[][] hits = {
            {"f00", "0", "0", "50", "50M", "+", "60", "", ""},
            {"f01", "0", "0", "150", "5S45M", "+", "20", "sec", ""},
            {"f01", "0", "1", "500", "50M", "+", "20", "sec", ""},
            {"f03", "0", "1", "700", "50M", "-", "50", "", ""},
            {"t00", "1", "0", "100", "50M", "+", "60", "", ""},
            {"t00", "2", "0", "300", "50M", "-", "60", "", ""},
            {"t01", "1", "0", "400", "50M", "+", "60", "", ""},
            {"t01", "2", "0", "390", "50M", "-", "55", "", ""},
            {"t02", "1", "1", "200", "50M", "+", "37", "", ""},
            {"t02", "2", "1", "200", "*", "+", "0", "unm", ""},
            {"t04", "1", "0", "600", "50M", "+", "30", "", "0"},
            {"t04", "2", "0", "750", "50M", "-", "30", "", "0"},
            {"t04", "1", "1", "600", "50M", "+", "5", "sec", "1"},
            {"t04", "2", "1", "750", "50M", "-", "5", "sec", "1"},
            {"t05", "1", "0", "900", "50M", "+", "20", "sec", "0"},
            {"t05", "2", "0", "1050", "50M", "-", "20", "sec", "0"},
            {"t05", "1", "0", "1200", "50M", "+", "40", "sec", "1"},
            {"t05", "2", "0", "1350", "50M", "-", "40", "sec", "1"},
            {"t06", "1", "0", "1500", "30M20S", "+", "60", "", ""},
            {"t06", "1", "1", "800", "30S20M", "+", "20", "sup", ""},
            {"t06", "2", "0", "1650", "50M", "-", "60", "", ""},
            {"t07", "1", "1", "960", "50M", "+", "60", "", ""},
            {"t07", "2", "1", "900", "50M", "-", "60", "", ""},
            {"t08", "1", "0", "1700", "20M1I10M1D19M", "+", "60", "", ""},
            {"t08", "2", "0", "1800", "50M", "-", "60", "", ""},
            {"t09", "1", "0", "200", "5H45M", "+", "60", "", ""},
            {"t09", "2", "0", "350", "45M5H", "-", "60", "", ""},
            {"t10", "1", "1", "300", "10S20M20S", "+", "25", "", ""},
            {"t10", "2", "1", "420", "50M", "-", "60", "", ""},
            {"t11", "1", "0", "1000", "50M", "+", "60", "", ""},
            {"t11", "2", "0", "1100", "50M", "-", "60", "", ""},
            {"t12", "1", "0", "500", "45M", "+", "60", "", ""},
            {"t12", "2", "0", "560", "50M", "-", "60", "", ""},
            {"t13", "1", "1", "100", "50M", "+", "60", "", ""},
            {"t13", "2", "1", "250", "50M", "-", "60", "", ""},
            {"t14", "1", "0", "1300", "50M", "+", "60", "", ""},
            {"t14", "2", "0", "1450", "50M", "-", "60", "", ""},
            {"t15", "1", "0", "1900", "50M", "+", "60", "", ""},
            {"t15", "2", "1", "50", "50M", "-", "60", "", ""},
            {"t16", "1", "0", "1600", "50M", "-", "60", "", ""},
            {"t16", "2", "0", "1700", "50M", "+", "60", "", ""},
        };
        String[] fragments = {"f00", "f01", "f02", "f03"};
        String[] pairs = new String[17];
        for (int i = 0; i < 17; i++) pairs[i] = String.format("t%02d", i);

        Random rng = new Random(20261003L);
        SAMFileHeader unmappedHeader = header(SAMFileHeader.SortOrder.queryname);
        java.util.List<SAMRecord> unmapped = new java.util.ArrayList<>();
        java.util.Map<String, byte[]> sequenced = new java.util.HashMap<>();
        java.util.List<String> names = new java.util.ArrayList<>();
        for (String f : fragments) names.add(f);
        for (String p : pairs) names.add(p);
        for (String name : names) {
            int ends = name.startsWith("f") ? 1 : 2;
            for (int end = (ends == 1 ? 0 : 1); end <= (ends == 1 ? 0 : 2); end++) {
                String[] truth = null;
                for (String[] h : hits) {
                    if (h[0].equals(name) && Integer.parseInt(h[1]) == end && !h[4].equals("*")) {
                        truth = h;
                        break;
                    }
                }
                byte[] bases = new byte[READ_LENGTH];
                if (truth == null) {
                    for (int b = 0; b < READ_LENGTH; b++) bases[b] = (byte) "ACGT".charAt(rng.nextInt(4));
                } else {
                    String contig = refs[Integer.parseInt(truth[2])];
                    int start = Integer.parseInt(truth[3]);
                    // Leading clips and hard clips are read bases before the aligned start.
                    int lead = 0;
                    java.util.regex.Matcher m = java.util.regex.Pattern.compile("^(\\d+)[SH]").matcher(truth[4]);
                    if (truth[5].equals("+") && m.find()) lead = Integer.parseInt(m.group(1));
                    java.util.regex.Matcher t = java.util.regex.Pattern.compile("(\\d+)[SH]$").matcher(truth[4]);
                    if (truth[5].equals("-") && t.find()) lead = Integer.parseInt(t.group(1));
                    for (int b = 0; b < READ_LENGTH; b++) {
                        int pos = start - 1 - (truth[5].equals("+") ? lead : 0) + b;
                        bases[b] = (pos >= 0 && pos < contig.length()) ? (byte) contig.charAt(pos) : (byte) 'A';
                        if (bases[b] == 'N') bases[b] = 'A';
                    }
                    bases[7] = mutateBase(bases[7]);
                    bases[33] = mutateBase(bases[33]);
                    if (truth[5].equals("-")) SequenceUtil.reverseComplement(bases);
                }
                byte[] quals = new byte[READ_LENGTH];
                for (int b = 0; b < READ_LENGTH; b++) quals[b] = (byte) (2 + rng.nextInt(39));
                SAMRecord r = new SAMRecord(unmappedHeader);
                r.setReadName(name);
                r.setReadBases(bases);
                r.setBaseQualities(quals);
                r.setReadUnmappedFlag(true);
                r.setReferenceIndex(SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX);
                r.setAlignmentStart(SAMRecord.NO_ALIGNMENT_START);
                r.setMappingQuality(0);
                if (end > 0) {
                    r.setReadPairedFlag(true);
                    r.setMateUnmappedFlag(true);
                    r.setFirstOfPairFlag(end == 1);
                    r.setSecondOfPairFlag(end == 2);
                }
                r.setAttribute("RG", name.compareTo("t08") < 0 ? "rg1" : "rg2");
                if (name.equals("t11")) r.setAttribute("XT", end == 1 ? 40 : 30);
                if ((name.equals("t14") && end == 2) || name.equals("f03")) {
                    StringBuilder oq = new StringBuilder();
                    for (int b = 0; b < READ_LENGTH; b++) oq.append((char) (33 + 10 + (b % 30)));
                    r.setAttribute("OQ", oq.toString());
                    r.setAttribute("E2", new String(bases, 0, READ_LENGTH).replace('A', 'C'));
                }
                if (name.equals("t13")) r.setAttribute("XU", "kept-" + end);
                sequenced.put(name + "/" + end, bases);
                unmapped.add(r);
            }
        }
        writeBam(new File(dir, "mba_unmapped.bam"), unmappedHeader, unmapped, false);

        for (int variant = 0; variant < 2; variant++) {
            boolean simple = variant == 1;
            SAMFileHeader alignedHeader = header(simple ? SAMFileHeader.SortOrder.unsorted : SAMFileHeader.SortOrder.queryname);
            SAMProgramRecord bwa = new SAMProgramRecord("bwa");
            bwa.setProgramName("bwa");
            bwa.setProgramVersion("0.7.17-r1188");
            bwa.setCommandLine("bwa mem ref.fasta reads.fq");
            alignedHeader.addProgramRecord(bwa);
            if (simple) {
                SAMProgramRecord other = new SAMProgramRecord("samtools");
                other.setProgramName("samtools");
                other.setPreviousProgramGroupId("bwa");
                alignedHeader.addProgramRecord(other);
            }
            java.util.List<SAMRecord> aligned = new java.util.ArrayList<>();
            java.util.Map<String, SAMRecord> byKey = new java.util.LinkedHashMap<>();
            for (String[] h : hits) {
                if (simple && (h[7].equals("sec") && !h[0].startsWith("f"))) continue;
                int end = Integer.parseInt(h[1]);
                byte[] read = sequenced.get(h[0] + "/" + end).clone();
                boolean negative = h[5].equals("-");
                if (negative) SequenceUtil.reverseComplement(read);
                SAMRecord r = new SAMRecord(alignedHeader);
                r.setReadName(h[0]);
                String cigar = h[4];
                int hardLead = 0, hardTail = 0;
                java.util.regex.Matcher lead = java.util.regex.Pattern.compile("^(\\d+)H").matcher(cigar);
                if (lead.find()) hardLead = Integer.parseInt(lead.group(1));
                java.util.regex.Matcher tail = java.util.regex.Pattern.compile("(\\d+)H$").matcher(cigar);
                if (tail.find()) hardTail = Integer.parseInt(tail.group(1));
                int length = cigar.equals("*") ? READ_LENGTH : TextCigarCodec.decode(cigar).getReadLength();
                if (cigar.equals("*")) {
                    r.setReadBases(read);
                } else {
                    r.setReadBases(java.util.Arrays.copyOfRange(read, hardLead, hardLead + length));
                }
                byte[] quals = new byte[r.getReadBases().length];
                for (int b = 0; b < quals.length; b++) quals[b] = (byte) (20 + (b % 15));
                r.setBaseQualities(quals);
                if (end > 0) {
                    r.setReadPairedFlag(true);
                    r.setFirstOfPairFlag(end == 1);
                    r.setSecondOfPairFlag(end == 2);
                }
                r.setReferenceIndex(Integer.parseInt(h[2]));
                r.setAlignmentStart(Integer.parseInt(h[3]));
                if (h[7].equals("unm")) {
                    r.setReadUnmappedFlag(true);
                    r.setMappingQuality(0);
                } else {
                    r.setCigarString(cigar);
                    r.setMappingQuality(Integer.parseInt(h[6]));
                    r.setReadNegativeStrandFlag(negative);
                    r.setAttribute("NM", 2);
                    r.setAttribute("AS", 40);
                    r.setAttribute("XS", 12);
                }
                r.setNotPrimaryAlignmentFlag(h[7].equals("sec"));
                r.setSupplementaryAlignmentFlag(h[7].equals("sup"));
                if (!h[8].isEmpty()) r.setAttribute("HI", Integer.parseInt(h[8]));
                r.setAttribute("RG", h[0].compareTo("t08") < 0 ? "rg1" : "rg2");
                if (h[0].equals("t13")) {
                    r.setAttribute("X0", 1);
                    r.setAttribute("ZZ", "zz-" + end);
                    r.setAttribute("YY", 7);
                    r.setAttribute("PG", "bwa");
                    r.setAttribute("MD", "50");
                    r.setAttribute("xa", "lower");
                }
                aligned.add(r);
                String key = h[0] + "/" + end + "/" + h[8] + "/" + h[7];
                byKey.put(key, r);
            }
            // Mate information as an aligner writes it: each hit with the other end's same hit.
            for (SAMRecord r : aligned) {
                if (!r.getReadPairedFlag() || !r.getFirstOfPairFlag() || r.getSupplementaryAlignmentFlag()) continue;
                Integer hi = r.getIntegerAttribute("HI");
                for (SAMRecord mate : aligned) {
                    if (mate.getReadName().equals(r.getReadName()) && mate.getSecondOfPairFlag()
                            && !mate.getSupplementaryAlignmentFlag()
                            && java.util.Objects.equals(hi, mate.getIntegerAttribute("HI"))
                            && mate.getNotPrimaryAlignmentFlag() == r.getNotPrimaryAlignmentFlag()) {
                        SamPairUtil.setMateInfo(r, mate, true);
                        SamPairUtil.setProperPairFlags(r, mate, java.util.Collections.singletonList(SamPairUtil.PairOrientation.FR));
                    }
                }
            }
            for (SAMRecord r : aligned) {
                if (!r.getSupplementaryAlignmentFlag()) continue;
                for (SAMRecord mate : aligned) {
                    if (mate.getReadName().equals(r.getReadName()) && mate.getSecondOfPairFlag()
                            && !mate.isSecondaryOrSupplementary()) {
                        SamPairUtil.setMateInformationOnSupplementalAlignment(r, mate, true);
                    }
                }
            }
            if (simple) {
                // Out of queryname order: the first template's records go last.
                java.util.List<SAMRecord> moved = new java.util.ArrayList<>();
                for (SAMRecord r : aligned) if (r.getReadName().equals("f00")) moved.add(r);
                aligned.removeAll(moved);
                aligned.addAll(moved);
                writeBam(new File(dir, "mba_aligned_simple.bam"), alignedHeader, aligned, false);
            } else {
                SAMFileWriterFactory factory = new SAMFileWriterFactory().setUseAsyncIo(false);
                try (SAMFileWriter w = factory.makeBAMWriter(alignedHeader, false, new File(dir, "mba_aligned.bam"))) {
                    for (SAMRecord r : aligned) w.addAlignment(r);
                }
            }
        }
    }
}
