/*
 * The VCF tools with a `.bcf` on either side, taken from the reference.
 *
 * Picard takes a BCF wherever it takes a VCF: VCFFileReader picks the BCF2 codec for a name ending
 * in `.bcf`, and VariantContextWriterBuilder picks BCF2Writer for an output named so. What comes
 * out is htsjdk's BCF, which is BCF 2.1 and uncompressed, and four of its behaviours decide the
 * bytes:
 *
 *   - INFO IS IN HASH ORDER, the bucket order of the HashMaps a record's attributes passed through;
 *   - A RECORD FROM A VCF IS FULLY DECODED FIRST, so a value the header disagrees with is the
 *     writer's refusal, and an integer list with a missing element is a NullPointerException;
 *   - GQ IS CAPPED AT 99;
 *   - A BCF RECORD NOBODY DECODED HAS ITS GENOTYPE BLOCK COPIED into the BCF written from it.
 *
 * Each case is a VCF, converted to BCF and back by VcfFormatConverter, and the BCF is then fed to
 * VcfFormatConverter again, to SortVcf and to MergeVcfs, while MakeSitesOnlyVcf writes a BCF from
 * the VCF. The inputs are in the output, so the port replays exactly what was measured.
 *
 * Output:
 *
 *     vcf\t<case>\t<hex of the input VCF>
 *     <step>\t<case>\t<hex of the output>
 *     error\t<case>/<step>\t<class>: <message>
 *
 * Usage: BcfDump
 */

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

public class BcfDump {

    static final StringBuilder buf = new StringBuilder();
    static Path dir;

    static void emit(final String kind, final String key, final String value) {
        buf.append(kind).append('\t').append(key).append('\t').append(value).append('\n');
    }

    static String hex(final byte[] bytes) {
        final StringBuilder b = new StringBuilder(bytes.length * 2);
        for (final byte x : bytes) {
            b.append(Character.forDigit((x >> 4) & 0xf, 16)).append(Character.forDigit(x & 0xf, 16));
        }
        return b.toString();
    }

    /** One tool run; the output's bytes, or the exception it threw. */
    static boolean run(final String name, final String step, final Path out, final picard.cmdline.CommandLineProgram tool,
                       final String... argv) {
        final PrintStream realErr = System.err;
        final PrintStream realOut = System.out;
        try {
            System.setErr(new PrintStream(new ByteArrayOutputStream(), true, StandardCharsets.UTF_8));
            System.setOut(new PrintStream(new ByteArrayOutputStream(), true, StandardCharsets.UTF_8));
            final int code = tool.instanceMain(argv);
            System.setErr(realErr);
            System.setOut(realOut);
            if (code != 0) {
                emit("error", name + "/" + step, "exit " + code);
                return false;
            }
            emit(step, name, hex(Files.readAllBytes(out)));
            return true;
        } catch (final Throwable e) {
            System.setErr(realErr);
            System.setOut(realOut);
            emit("error", name + "/" + step, e.getClass().getName() + ": "
                    + String.valueOf(e.getMessage()).replace(dir.toAbsolutePath().toString(), "<dir>")
                    .replace('\n', ' ').replace('\t', ' '));
            return false;
        }
    }

    public static void main(final String[] args) throws Exception {
        dir = Files.createTempDirectory("bcf-dump");
        final Map<String, String> cases = new LinkedHashMap<>();
        cases.put("edge", """
##fileformat=VCFv4.2
##FILTER=<ID=LowQual,Description="Low quality">
##FILTER=<ID=q10,Description="Quality below 10">
##FILTER=<ID=s50,Description="Less than 50% of samples">
##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral allele">
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total alleles">
##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP membership">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">
##INFO=<ID=H2,Number=0,Type=Flag,Description="HapMap2">
##INFO=<ID=MQ,Number=1,Type=Float,Description="Mapping quality">
##INFO=<ID=NS,Number=1,Type=Integer,Description="Samples with data">
##INFO=<ID=BIG,Number=1,Type=Integer,Description="Large integer">
##INFO=<ID=MID,Number=.,Type=Integer,Description="Unbounded ints">
##INFO=<ID=STRS,Number=.,Type=String,Description="Strings">
##INFO=<ID=CH,Number=1,Type=Character,Description="Char">
##INFO=<ID=RF,Number=R,Type=Float,Description="Per allele floats">
##INFO=<ID=GI,Number=G,Type=Integer,Description="Per genotype ints">
##INFO=<ID=F2,Number=2,Type=Float,Description="Two floats">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
##FORMAT=<ID=FT,Number=1,Type=String,Description="Genotype filter">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
##FORMAT=<ID=HQ,Number=2,Type=Integer,Description="Haplotype quality">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred likelihoods">
##FORMAT=<ID=GF,Number=1,Type=Float,Description="A float">
##FORMAT=<ID=GS,Number=1,Type=String,Description="A string">
##FORMAT=<ID=MI,Number=.,Type=Integer,Description="Unbounded ints">
##contig=<ID=20,length=62435964>
##contig=<ID=21,length=48129895>
##contig=<ID=X,length=155270560>
##reference=file:///ref.fa
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA1	NA2	NA3
20	14370	rs6054257	G	A	29	PASS	NS=3;DP=14;AF=0.5;DB;H2;AA=G;AN=6;AC=3;MQ=60.25;BIG=70000;MID=1,300,70000;STRS=a,bb,ccc;CH=Z;RF=0.1,0.2;GI=1,2,3;F2=1.5,2.5	GT:GQ:DP:HQ:AD:PL:GF:GS:MI	0|0:48:1:51,51:1,0:0,10,100:0.5:x:1	1|0:48:8:51,51:3,4:10,0,100:1.25:yy:1,2	1/1:43:5:.,.:.:.:.:.:.
20	17330	.	T	A	3	q10	NS=3;DP=11;AF=0.017	GT:GQ:DP:HQ:FT	0|0:49:3:58,50:PASS	0|1:3:5:65,3:q10	0/0:41:3:.:q10;s50
20	1110696	rs6040355	A	G,T	67	PASS	NS=2;DP=10;AF=0.333,0.667;AA=T;DB;AC=1,2;GI=1,2,3,4,5,6	GT:GQ:DP:HQ:AD:PL	1|2:21:6:23,27:1,2,3:1,2,3,4,5,6	2|1:2:0:18,2:.:.	2/2:35:4:.:4,5,6:7,8,9,10,11,12
20	1230237	.	T	.	47	.	NS=3;DP=13;AA=T	GT:GQ:DP:HQ	0|0:54:7:56,60	0|0:48:4:51,51	0/0:61:2:.
20	1234567	microsat1	GTC	G,GTCT	50	PASS	NS=3;DP=9;AA=G;MID=5	GT:GQ:DP	0/1:35:4	0/2:17:2	1/1:40:3
20	1234568	.	G	GA,GAC,GACT	.	q10;s50	DP=200;BIG=-2147483647	GT:GQ:PL	0/3:120:0,1,2,3,4,5,6,7,8,9	.	1|3:.:.
21	100	.	A	C	0	LowQual	.	GT:DP	./.:.	0:3	.|1:2
21	200	id1;id2	ACGT	A	1e-5	PASS	DP=0;MQ=0.001;AF=1e-30	GT	0/1	1/1	0/0
X	10	.	C	CAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA	123.456	PASS	STRS=onlyone;CH=.;DP=.	GT:GS:GF	1:hello:.	0:.:3.14159	.:z:-0.0
X	20	.	A	<DEL>	99	PASS	.	GT:AD	0/1:5,5	0/0:10,0	1/1:0,10
""");
        cases.put("edge-sorted", """
##fileformat=VCFv4.2
##FILTER=<ID=LowQual,Description="Low quality">
##FILTER=<ID=q10,Description="Quality below 10">
##FILTER=<ID=s50,Description="Less than 50% of samples">
##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral allele">
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total alleles">
##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP membership">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">
##INFO=<ID=H2,Number=0,Type=Flag,Description="HapMap2">
##INFO=<ID=MQ,Number=1,Type=Float,Description="Mapping quality">
##INFO=<ID=NS,Number=1,Type=Integer,Description="Samples with data">
##INFO=<ID=BIG,Number=1,Type=Integer,Description="Large integer">
##INFO=<ID=MID,Number=.,Type=Integer,Description="Unbounded ints">
##INFO=<ID=STRS,Number=.,Type=String,Description="Strings">
##INFO=<ID=CH,Number=1,Type=Character,Description="Char">
##INFO=<ID=RF,Number=R,Type=Float,Description="Per allele floats">
##INFO=<ID=GI,Number=G,Type=Integer,Description="Per genotype ints">
##INFO=<ID=F2,Number=2,Type=Float,Description="Two floats">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
##FORMAT=<ID=FT,Number=1,Type=String,Description="Genotype filter">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
##FORMAT=<ID=HQ,Number=2,Type=Integer,Description="Haplotype quality">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred likelihoods">
##FORMAT=<ID=GF,Number=1,Type=Float,Description="A float">
##FORMAT=<ID=GS,Number=1,Type=String,Description="A string">
##FORMAT=<ID=MI,Number=.,Type=Integer,Description="Unbounded ints">
##contig=<ID=20,length=62435964>
##contig=<ID=21,length=48129895>
##contig=<ID=X,length=155270560>
##reference=file:///ref.fa
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3
20	14370	rs6054257	G	A	29	PASS	NS=3;DP=14;AF=0.5;DB;H2;AA=G;AN=6;AC=3;MQ=60.25;BIG=70000;MID=1,300,70000;STRS=a,bb,ccc;CH=Z;RF=0.1,0.2;GI=1,2,3;F2=1.5,2.5	GT:GQ:DP:HQ:AD:PL:GF:GS:MI	0|0:48:1:51,51:1,0:0,10,100:0.5:x:1	1|0:48:8:51,51:3,4:10,0,100:1.25:yy:1,2	1/1:43:5:.,.:.:.:.:.:.
20	17330	.	T	A	3	q10	NS=3;DP=11;AF=0.017	GT:GQ:DP:HQ:FT	0|0:49:3:58,50:PASS	0|1:3:5:65,3:q10	0/0:41:3:.:q10;s50
20	1110696	rs6040355	A	G,T	67	PASS	NS=2;DP=10;AF=0.333,0.667;AA=T;DB;AC=1,2;GI=1,2,3,4,5,6	GT:GQ:DP:HQ:AD:PL	1|2:21:6:23,27:1,2,3:1,2,3,4,5,6	2|1:2:0:18,2:.:.	2/2:35:4:.:4,5,6:7,8,9,10,11,12
20	1230237	.	T	.	47	.	NS=3;DP=13;AA=T	GT:GQ:DP:HQ	0|0:54:7:56,60	0|0:48:4:51,51	0/0:61:2:.
20	1234567	microsat1	GTC	G,GTCT	50	PASS	NS=3;DP=9;AA=G;MID=5	GT:GQ:DP	0/1:35:4	0/2:17:2	1/1:40:3
20	1234568	.	G	GA,GAC,GACT	.	q10;s50	DP=200;BIG=-2147483647	GT:GQ:PL	0/3:120:0,1,2,3,4,5,6,7,8,9	.	1|3:.:.
21	100	.	A	C	0	LowQual	.	GT:DP	./.:.	0:3	.|1:2
21	200	id1;id2	ACGT	A	1e-5	PASS	DP=0;MQ=0.001;AF=1e-30	GT	0/1	1/1	0/0
X	10	.	C	CAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA	123.456	PASS	STRS=onlyone;CH=.;DP=.	GT:GS:GF	1:hello:.	0:.:3.14159	.:z:-0.0
X	20	.	A	<DEL>	99	PASS	.	GT:AD	0/1:5,5	0/0:10,0	1/1:0,10
""");
        cases.put("edge2", """
##fileformat=VCFv4.2
##contig=<ID=1,length=1000000>
##contig=<ID=2,length=1000000>
##INFO=<ID=K01,Number=1,Type=Integer,Description="k">
##INFO=<ID=K02,Number=1,Type=Integer,Description="k">
##INFO=<ID=K03,Number=1,Type=Integer,Description="k">
##INFO=<ID=K04,Number=1,Type=Integer,Description="k">
##INFO=<ID=K05,Number=1,Type=Integer,Description="k">
##INFO=<ID=K06,Number=1,Type=Integer,Description="k">
##INFO=<ID=K07,Number=1,Type=Integer,Description="k">
##INFO=<ID=K08,Number=1,Type=Integer,Description="k">
##INFO=<ID=K09,Number=1,Type=Integer,Description="k">
##INFO=<ID=K10,Number=1,Type=Integer,Description="k">
##INFO=<ID=K11,Number=1,Type=Integer,Description="k">
##INFO=<ID=K12,Number=1,Type=Integer,Description="k">
##INFO=<ID=K13,Number=1,Type=Integer,Description="k">
##INFO=<ID=K14,Number=1,Type=Integer,Description="k">
##INFO=<ID=K15,Number=1,Type=Integer,Description="k">
##INFO=<ID=K16,Number=1,Type=Integer,Description="k">
##INFO=<ID=K17,Number=1,Type=Integer,Description="k">
##INFO=<ID=K18,Number=1,Type=Integer,Description="k">
##INFO=<ID=K19,Number=1,Type=Integer,Description="k">
##INFO=<ID=K20,Number=1,Type=Integer,Description="k">
##INFO=<ID=K21,Number=1,Type=Integer,Description="k">
##INFO=<ID=K22,Number=1,Type=Integer,Description="k">
##INFO=<ID=K23,Number=1,Type=Integer,Description="k">
##INFO=<ID=K24,Number=1,Type=Integer,Description="k">
##INFO=<ID=K25,Number=1,Type=Integer,Description="k">
##INFO=<ID=Zeta,Number=1,Type=Integer,Description="k">
##INFO=<ID=alpha,Number=1,Type=Integer,Description="k">
##INFO=<ID=AC2,Number=1,Type=Integer,Description="k">
##INFO=<ID=BB,Number=1,Type=Integer,Description="k">
##INFO=<ID=LONG,Number=.,Type=Integer,Description="many">
##INFO=<ID=FLT,Number=.,Type=Float,Description="many floats">
##INFO=<ID=FT,Number=1,Type=String,Description="ft in info">
##FORMAT=<ID=GT,Number=1,Type=String,Description="g">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="d">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="q">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="a">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="p">
##FORMAT=<ID=BIGI,Number=1,Type=Integer,Description="b">
##FORMAT=<ID=VEC,Number=.,Type=Integer,Description="v">
##FORMAT=<ID=FV,Number=.,Type=Float,Description="fv">
##FORMAT=<ID=SV,Number=.,Type=String,Description="sv">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	b	a	c
1	100	.	A	C	50	PASS	K01=01;K02=02;K03=03;K04=04;K05=05;K06=06;K07=07;K08=08;K09=09;K10=10;K11=11;K12=12;K13=13;K14=14;K15=15;K16=16;K17=17;K18=18;K19=19;K20=20;K21=21;K22=22;K23=23;K24=24;K25=25;Zeta=1;alpha=2;AC2=3;BB=4;LONG=1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20;FLT=0.1,0.25,1e10,-3.5,1.0e-40;FT=x	GT:DP:GQ:AD:PL:BIGI:VEC:FV:SV	0/1:300:20:100,200:0,40000,90000:70000:1,2,3:0.5,1.5:a,b	1/1/1:5:150:1,2:.:-200:.:.:.	0|1:.:.:3,.:1,2,3:.:1:2.5:onlyone
1	200	.	A	C,G,T,AC,AG,AT,CC,CG,CT,GA,GC,GG,GT,TA,TC	.	.	K01=1	GT:AD	15/14:1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16	./.:.	0:.
1	300	.	AC	A	10	PASS	K02=-130;K03=40000;K04=-40000	GT	1	0/0/0/1	.
2	5	rsX	G	A	1.5	PASS	.	GT:DP	0/1:1	0/1:2	0/1:3
""");
        cases.put("sitesonly", """
##fileformat=VCFv4.2
##FILTER=<ID=LowQual,Description="Low quality">
##FILTER=<ID=q10,Description="Quality below 10">
##FILTER=<ID=s50,Description="Less than 50% of samples">
##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral allele">
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total alleles">
##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP membership">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">
##INFO=<ID=H2,Number=0,Type=Flag,Description="HapMap2">
##INFO=<ID=MQ,Number=1,Type=Float,Description="Mapping quality">
##INFO=<ID=NS,Number=1,Type=Integer,Description="Samples with data">
##INFO=<ID=BIG,Number=1,Type=Integer,Description="Large integer">
##INFO=<ID=MID,Number=.,Type=Integer,Description="Unbounded ints">
##INFO=<ID=STRS,Number=.,Type=String,Description="Strings">
##INFO=<ID=CH,Number=1,Type=Character,Description="Char">
##INFO=<ID=RF,Number=R,Type=Float,Description="Per allele floats">
##INFO=<ID=GI,Number=G,Type=Integer,Description="Per genotype ints">
##INFO=<ID=F2,Number=2,Type=Float,Description="Two floats">
##contig=<ID=20,length=62435964>
##contig=<ID=21,length=48129895>
##contig=<ID=X,length=155270560>
##reference=file:///ref.fa
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
20	14370	rs6054257	G	A	29	PASS	NS=3;DP=14;AF=0.5;DB;H2;AA=G;AN=6;AC=3;MQ=60.25;BIG=70000;MID=1,300,70000;STRS=a,bb,ccc;CH=Z;RF=0.1,0.2;GI=1,2,3;F2=1.5,2.5
20	17330	.	T	A	3	q10	NS=3;DP=11;AF=0.017
20	1110696	rs6040355	A	G,T	67	PASS	NS=2;DP=10;AF=0.333,0.667;AA=T;DB;AC=1,2;GI=1,2,3,4,5,6
20	1230237	.	T	.	47	.	NS=3;DP=13;AA=T
20	1234567	microsat1	GTC	G,GTCT	50	PASS	NS=3;DP=9;AA=G;MID=5
20	1234568	.	G	GA,GAC,GACT	.	q10;s50	DP=200;BIG=-2147483647
21	100	.	A	C	0	LowQual	.
21	200	id1;id2	ACGT	A	1e-5	PASS	DP=0;MQ=0.001;AF=1e-30
X	10	.	C	CAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA	123.456	PASS	STRS=onlyone;CH=.;DP=.
X	20	.	A	<DEL>	99	PASS	.
""");
        cases.put("trio", """
##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000000>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">
##INFO=<ID=AC,Number=A,Type=Integer,Description="AC">
##INFO=<ID=AF,Number=A,Type=Float,Description="AF">
##INFO=<ID=AN,Number=1,Type=Integer,Description="AN">
##FILTER=<ID=LowQual,Description="Low quality">
##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="AD">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="DP">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="GQ">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="PL">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	child	father	mother
chr1	10	.	A	C	2475	PASS	DP=26;AC=1;AF=0.167;AN=6	GT:AD:DP:GQ:PL	1/1:13,7:20:18:13,21,189	./.	0/1:7,0:7:12:27,7,129
chr1	11	.	A	T	1566	PASS	DP=44;AC=1;AF=0.167;AN=6	GT:AD:DP:GQ:PL	./.	0/1:25,17:42:1:20,22,108	1/1:18,4:22:28:43,3,23
chr1	12	.	A	G	2374	PASS	DP=69;AC=1;AF=0.167;AN=6	GT:AD:DP:GQ:PL	./.	0/1:32,7:39:49:10,17,75	./.
chr1	13	.	A	C	2505	PASS	DP=63;AC=1;AF=0.167;AN=6	GT:AD:DP:GQ:PL	0/1:36,18:54:11:29,27,25	0/1:8,14:22:82:46,5,94	1/1:10,8:18:90:87,20,18
chr1	14	.	A	C	2950	PASS	DP=100;AC=1;AF=0.167;AN=6	GT:AD:DP:GQ:PL	0/1:10,24:34:35:81,22,142	0/1:28,20:48:99:7,7,8	1/1:22,8:30:9:27,29,145
chr1	15	.	A	G	1492	PASS	DP=76;AC=1;AF=0.167;AN=6	GT:AD:DP:GQ:PL	0/1:17,29:46:19:33,4,63	./.	./.
chr1	16	.	AAA	A	2177	PASS	DP=84;AC=1;AF=0.167;AN=6	GT:AD:DP:GQ:PL	./.	0/0:53,7:60:20:80,5,174	0/1:39,4:43:50:48,19,119
chr1	17	.	AA	A	1232	PASS	DP=149;AC=1;AF=0.167;AN=6	GT:AD:DP:GQ:PL	0/0:14,34:48:97:34,24,164	1/1:8,4:12:56:20,14,0	1/1:26,11:37:65:13,27,160
chr1	18	.	A	C	2333	PASS	DP=40;AC=1;AF=0.167;AN=6	GT:AD:DP:GQ:PL	1/1:43,10:53:70:67,29,0	./.	0/0:7,5:12:40:30,1,61
chr1	19	.	A	G	1286	PASS	DP=122;AC=1;AF=0.167;AN=6	GT:AD:DP:GQ:PL	0/0:19,34:53:99:16,4,168	0/1:30,10:40:34:67,27,155	0/1:1,17:18:97:88,6,182
chr1	20	.	AAAAAAA	A	286	PASS	DP=28;AC=1;AF=0.167;AN=6	GT:AD:DP:GQ:PL	1/1:0,33:33:58:15,7,57	0/0:26,0:26:76:70,7,150	0/1:5,0:5:91:80,1,58
chr1	21	.	A	AC	2709	PASS	DP=130;AC=1;AF=0.167;AN=6	GT:AD:DP:GQ:PL	1/1:16,31:47:28:69,4,185	./.	0/1:25,30:55:53:24,3,24
chr1	22	.	A	G	761	PASS	DP=91;AC=1;AF=0.167;AN=6	GT:AD:DP:GQ:PL	0/1:14,46:60:7:86,20,165	0/0:2,6:8:94:43,25,27	0/1:11,6:17:69:57,4,108
chr1	23	.	A	C	885	PASS	DP=122;AC=1;AF=0.167;AN=6	GT:AD:DP:GQ:PL	0/1:1,55:56:71:12,1,166	./.	0/0:38,15:53:22:52,15,123
chr1	24	.	A	ATA	2382	PASS	DP=158;AC=1;AF=0.167;AN=6	GT:AD:DP:GQ:PL	0/1:7,14:21:37:54,22,187	./.	0/1:11,3:14:38:27,30,14
chr1	25	.	A	G	2343	PASS	DP=83;AC=1;AF=0.167;AN=6	GT:AD:DP:GQ:PL	0/0:1,7:8:65:67,5,14	./.	0/0:39,4:43:87:30,12,30
chr1	26	.	A	C	1631	PASS	DP=53;AC=1;AF=0.167;AN=6	GT:AD:DP:GQ:PL	./.	./.	1/1:15,6:21:86:40,7,67
chr1	27	.	A	G	1096	PASS	DP=53;AC=1;AF=0.167;AN=6	GT:AD:DP:GQ:PL	0/1:1,24:25:10:1,14,159	./.	./.
chr1	28	.	AA	A	2276	PASS	DP=59;AC=1;AF=0.167;AN=6	GT:AD:DP:GQ:PL	0/1:19,9:28:21:56,26,139	1/1:3,41:44:68:1,21,141	1/1:41,6:47:18:33,3,27
chr1	29	.	A	T	1376	PASS	DP=53;AC=1;AF=0.167;AN=6	GT:AD:DP:GQ:PL	0/1:29,21:50:27:87,20,67	./.	0/0:0,10:10:55:35,1,0
chr1	30	.	A	G	1272	PASS	DP=113;AC=1;AF=0.167;AN=6	GT:AD:DP:GQ:PL	0/1:24,28:52:71:90,13,143	0/0:11,1:12:89:19,17,9	1/1:7,35:42:19:55,4,10
chr1	31	.	A	AG	979	PASS	DP=61;AC=1;AF=0.167;AN=6	GT:AD:DP:GQ:PL	0/1:33,15:48:86:13,11,199	./.	./.
chr1	32	.	AAA	A	923	PASS	DP=26;AC=1;AF=0.167;AN=6	GT:AD:DP:GQ:PL	0/1:5,1:6:95:42,25,105	0/1:17,5:22:90:13,12,9	0/1:13,6:19:59:44,9,58
chr1	33	.	A	G	482	PASS	DP=86;AC=1;AF=0.167;AN=6	GT:AD:DP:GQ:PL	1/1:20,2:22:99:35,11,164	./.	./.
chr1	34	.	A	G	2215	PASS	DP=70;AC=1;AF=0.167;AN=6	GT:AD:DP:GQ:PL	0/0:2,9:11:56:44,23,80	0/1:11,32:43:15:49,28,147	0/1:20,1:21:91:55,0,133
chr1	35	.	A	C	795	PASS	DP=127;AC=1;AF=0.167;AN=6	GT:AD:DP:GQ:PL	1/1:24,20:44:85:15,23,76	./.	0/1:13,12:25:90:37,17,32
chr1	36	.	A	G	1948	PASS	DP=63;AC=1;AF=0.167;AN=6	GT:AD:DP:GQ:PL	0/1:8,36:44:39:51,17,0	1/1:17,6:23:56:74,19,167	1/1:6,28:34:57:86,6,130
chr1	37	.	A	G	199	PASS	DP=82;AC=1;AF=0.167;AN=6	GT:AD:DP:GQ:PL	./.	./.	0/1:29,19:48:29:25,4,6
chr1	38	.	AAAAAA	A	2860	PASS	DP=152;AC=1;AF=0.167;AN=6	GT:AD:DP:GQ:PL	0/0:8,26:34:81:73,6,183	0/1:11,25:36:32:18,20,176	0/0:4,49:53:14:54,7,45
chr1	39	.	A	T	1836	PASS	DP=60;AC=1;AF=0.167;AN=6	GT:AD:DP:GQ:PL	0/1:52,7:59:59:17,25,118	./.	1/1:25,28:53:79:64,13,140
chr1	40	.	A	G	1385	PASS	DP=101;AC=1;AF=0.167;AN=6	GT:AD:DP:GQ:PL	0/1:14,7:21:82:35,24,199	./.	1/1:29,4:33:92:36,7,69
chr1	41	.	A	AC	90	PASS	DP=167;AC=1;AF=0.167;AN=6	GT:AD:DP:GQ:PL	0/1:7,12:19:89:19,22,54	0/0:5,26:31:43:69,14,106	0/0:5,13:18:50:74,30,178
chr1	42	.	A	C	1602	PASS	DP=106;AC=1;AF=0.167;AN=6	GT:AD:DP:GQ:PL	1/1:0,24:24:50:53,17,191	./.	0/1:22,14:36:35:55,15,7
chr1	43	.	A	G	1562	PASS	DP=103;AC=1;AF=0.167;AN=6	GT:AD:DP:GQ:PL	0/1:29,29:58:17:79,17,6	0/1:6,36:42:85:3,2,164	0/1:0,13:13:60:23,1,66
chr1	44	.	A	G	1022	PASS	DP=71;AC=1;AF=0.167;AN=6	GT:AD:DP:GQ:PL	1/1:29,24:53:36:53,8,20	0/1:1,5:6:70:6,30,89	0/1:42,4:46:84:5,24,7
chr1	45	.	A	T	479	PASS	DP=61;AC=1;AF=0.167;AN=6	GT:AD:DP:GQ:PL	0/1:16,4:20:61:85,3,144	0/1:18,16:34:99:47,5,155	./.
chr1	46	.	AA	A	1004	PASS	DP=46;AC=1;AF=0.167;AN=6	GT:AD:DP:GQ:PL	./.	./.	0/1:38,12:50:10:75,22,160
chr1	47	.	A	G	2662	PASS	DP=107;AC=1;AF=0.167;AN=6	GT:AD:DP:GQ:PL	./.	./.	1/1:12,27:39:85:47,2,129
chr1	48	.	A	G	2436	PASS	DP=88;AC=1;AF=0.167;AN=6	GT:AD:DP:GQ:PL	0/1:5,6:11:47:81,28,117	0/1:21,11:32:94:66,30,166	1/1:10,34:44:62:59,13,187
chr1	49	.	A	C	2882	PASS	DP=90;AC=1;AF=0.167;AN=6	GT:AD:DP:GQ:PL	0/0:8,14:22:32:59,18,156	0/1:26,0:26:64:41,5,124	0/1:2,25:27:34:43,8,152
""");
        cases.put("refuse-okflag", """
##fileformat=VCFv4.2
##FILTER=<ID=LowQual,Description="Low quality">
##FILTER=<ID=q10,Description="Quality below 10">
##FILTER=<ID=s50,Description="Less than 50% of samples">
##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral allele">
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total alleles">
##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP membership">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">
##INFO=<ID=H2,Number=0,Type=Flag,Description="HapMap2">
##INFO=<ID=MQ,Number=1,Type=Float,Description="Mapping quality">
##INFO=<ID=NS,Number=1,Type=Integer,Description="Samples with data">
##INFO=<ID=BIG,Number=1,Type=Integer,Description="Large integer">
##INFO=<ID=MID,Number=.,Type=Integer,Description="Unbounded ints">
##INFO=<ID=STRS,Number=.,Type=String,Description="Strings">
##INFO=<ID=CH,Number=1,Type=Character,Description="Char">
##INFO=<ID=RF,Number=R,Type=Float,Description="Per allele floats">
##INFO=<ID=GI,Number=G,Type=Integer,Description="Per genotype ints">
##INFO=<ID=F2,Number=2,Type=Float,Description="Two floats">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
##FORMAT=<ID=FT,Number=1,Type=String,Description="Genotype filter">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
##FORMAT=<ID=HQ,Number=2,Type=Integer,Description="Haplotype quality">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred likelihoods">
##FORMAT=<ID=GF,Number=1,Type=Float,Description="A float">
##FORMAT=<ID=GS,Number=1,Type=String,Description="A string">
##FORMAT=<ID=MI,Number=.,Type=Integer,Description="Unbounded ints">
##contig=<ID=20,length=62435964>
##contig=<ID=21,length=48129895>
##contig=<ID=X,length=155270560>
##reference=file:///ref.fa
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA1
20	10	.	A	C	5	PASS	DB=1;H2=true;AA=.	GT	0/1
""");
        cases.put("refuse-nullint", """
##fileformat=VCFv4.2
##FILTER=<ID=LowQual,Description="Low quality">
##FILTER=<ID=q10,Description="Quality below 10">
##FILTER=<ID=s50,Description="Less than 50% of samples">
##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral allele">
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total alleles">
##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP membership">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">
##INFO=<ID=H2,Number=0,Type=Flag,Description="HapMap2">
##INFO=<ID=MQ,Number=1,Type=Float,Description="Mapping quality">
##INFO=<ID=NS,Number=1,Type=Integer,Description="Samples with data">
##INFO=<ID=BIG,Number=1,Type=Integer,Description="Large integer">
##INFO=<ID=MID,Number=.,Type=Integer,Description="Unbounded ints">
##INFO=<ID=STRS,Number=.,Type=String,Description="Strings">
##INFO=<ID=CH,Number=1,Type=Character,Description="Char">
##INFO=<ID=RF,Number=R,Type=Float,Description="Per allele floats">
##INFO=<ID=GI,Number=G,Type=Integer,Description="Per genotype ints">
##INFO=<ID=F2,Number=2,Type=Float,Description="Two floats">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
##FORMAT=<ID=FT,Number=1,Type=String,Description="Genotype filter">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
##FORMAT=<ID=HQ,Number=2,Type=Integer,Description="Haplotype quality">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred likelihoods">
##FORMAT=<ID=GF,Number=1,Type=Float,Description="A float">
##FORMAT=<ID=GS,Number=1,Type=String,Description="A string">
##FORMAT=<ID=MI,Number=.,Type=Integer,Description="Unbounded ints">
##contig=<ID=20,length=62435964>
##contig=<ID=21,length=48129895>
##contig=<ID=X,length=155270560>
##reference=file:///ref.fa
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA1
20	10	.	A	C,G	5	PASS	AC=1,.	GT	0/1
""");
        cases.put("refuse-discord", """
##fileformat=VCFv4.2
##FILTER=<ID=LowQual,Description="Low quality">
##FILTER=<ID=q10,Description="Quality below 10">
##FILTER=<ID=s50,Description="Less than 50% of samples">
##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral allele">
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total alleles">
##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP membership">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">
##INFO=<ID=H2,Number=0,Type=Flag,Description="HapMap2">
##INFO=<ID=MQ,Number=1,Type=Float,Description="Mapping quality">
##INFO=<ID=NS,Number=1,Type=Integer,Description="Samples with data">
##INFO=<ID=BIG,Number=1,Type=Integer,Description="Large integer">
##INFO=<ID=MID,Number=.,Type=Integer,Description="Unbounded ints">
##INFO=<ID=STRS,Number=.,Type=String,Description="Strings">
##INFO=<ID=CH,Number=1,Type=Character,Description="Char">
##INFO=<ID=RF,Number=R,Type=Float,Description="Per allele floats">
##INFO=<ID=GI,Number=G,Type=Integer,Description="Per genotype ints">
##INFO=<ID=F2,Number=2,Type=Float,Description="Two floats">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
##FORMAT=<ID=FT,Number=1,Type=String,Description="Genotype filter">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
##FORMAT=<ID=HQ,Number=2,Type=Integer,Description="Haplotype quality">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred likelihoods">
##FORMAT=<ID=GF,Number=1,Type=Float,Description="A float">
##FORMAT=<ID=GS,Number=1,Type=String,Description="A string">
##FORMAT=<ID=MI,Number=.,Type=Integer,Description="Unbounded ints">
##contig=<ID=20,length=62435964>
##contig=<ID=21,length=48129895>
##contig=<ID=X,length=155270560>
##reference=file:///ref.fa
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA1
20	10	.	A	C	5	PASS	AC=1,2	GT	0/1
""");
        cases.put("refuse-noinfo", """
##fileformat=VCFv4.2
##FILTER=<ID=LowQual,Description="Low quality">
##FILTER=<ID=q10,Description="Quality below 10">
##FILTER=<ID=s50,Description="Less than 50% of samples">
##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral allele">
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total alleles">
##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP membership">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">
##INFO=<ID=H2,Number=0,Type=Flag,Description="HapMap2">
##INFO=<ID=MQ,Number=1,Type=Float,Description="Mapping quality">
##INFO=<ID=NS,Number=1,Type=Integer,Description="Samples with data">
##INFO=<ID=BIG,Number=1,Type=Integer,Description="Large integer">
##INFO=<ID=MID,Number=.,Type=Integer,Description="Unbounded ints">
##INFO=<ID=STRS,Number=.,Type=String,Description="Strings">
##INFO=<ID=CH,Number=1,Type=Character,Description="Char">
##INFO=<ID=RF,Number=R,Type=Float,Description="Per allele floats">
##INFO=<ID=GI,Number=G,Type=Integer,Description="Per genotype ints">
##INFO=<ID=F2,Number=2,Type=Float,Description="Two floats">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
##FORMAT=<ID=FT,Number=1,Type=String,Description="Genotype filter">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
##FORMAT=<ID=HQ,Number=2,Type=Integer,Description="Haplotype quality">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred likelihoods">
##FORMAT=<ID=GF,Number=1,Type=Float,Description="A float">
##FORMAT=<ID=GS,Number=1,Type=String,Description="A string">
##FORMAT=<ID=MI,Number=.,Type=Integer,Description="Unbounded ints">
##contig=<ID=20,length=62435964>
##contig=<ID=21,length=48129895>
##contig=<ID=X,length=155270560>
##reference=file:///ref.fa
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA1
20	10	.	A	C	5	PASS	ZZ=1	GT	0/1
""");
        cases.put("refuse-noformat", """
##fileformat=VCFv4.2
##FILTER=<ID=LowQual,Description="Low quality">
##FILTER=<ID=q10,Description="Quality below 10">
##FILTER=<ID=s50,Description="Less than 50% of samples">
##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral allele">
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total alleles">
##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP membership">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">
##INFO=<ID=H2,Number=0,Type=Flag,Description="HapMap2">
##INFO=<ID=MQ,Number=1,Type=Float,Description="Mapping quality">
##INFO=<ID=NS,Number=1,Type=Integer,Description="Samples with data">
##INFO=<ID=BIG,Number=1,Type=Integer,Description="Large integer">
##INFO=<ID=MID,Number=.,Type=Integer,Description="Unbounded ints">
##INFO=<ID=STRS,Number=.,Type=String,Description="Strings">
##INFO=<ID=CH,Number=1,Type=Character,Description="Char">
##INFO=<ID=RF,Number=R,Type=Float,Description="Per allele floats">
##INFO=<ID=GI,Number=G,Type=Integer,Description="Per genotype ints">
##INFO=<ID=F2,Number=2,Type=Float,Description="Two floats">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
##FORMAT=<ID=FT,Number=1,Type=String,Description="Genotype filter">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
##FORMAT=<ID=HQ,Number=2,Type=Integer,Description="Haplotype quality">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred likelihoods">
##FORMAT=<ID=GF,Number=1,Type=Float,Description="A float">
##FORMAT=<ID=GS,Number=1,Type=String,Description="A string">
##FORMAT=<ID=MI,Number=.,Type=Integer,Description="Unbounded ints">
##contig=<ID=20,length=62435964>
##contig=<ID=21,length=48129895>
##contig=<ID=X,length=155270560>
##reference=file:///ref.fa
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA1
20	10	.	A	C	5	PASS	.	GT:ZZ	0/1:3
""");
        cases.put("refuse-badint", """
##fileformat=VCFv4.2
##FILTER=<ID=LowQual,Description="Low quality">
##FILTER=<ID=q10,Description="Quality below 10">
##FILTER=<ID=s50,Description="Less than 50% of samples">
##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral allele">
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total alleles">
##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP membership">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">
##INFO=<ID=H2,Number=0,Type=Flag,Description="HapMap2">
##INFO=<ID=MQ,Number=1,Type=Float,Description="Mapping quality">
##INFO=<ID=NS,Number=1,Type=Integer,Description="Samples with data">
##INFO=<ID=BIG,Number=1,Type=Integer,Description="Large integer">
##INFO=<ID=MID,Number=.,Type=Integer,Description="Unbounded ints">
##INFO=<ID=STRS,Number=.,Type=String,Description="Strings">
##INFO=<ID=CH,Number=1,Type=Character,Description="Char">
##INFO=<ID=RF,Number=R,Type=Float,Description="Per allele floats">
##INFO=<ID=GI,Number=G,Type=Integer,Description="Per genotype ints">
##INFO=<ID=F2,Number=2,Type=Float,Description="Two floats">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
##FORMAT=<ID=FT,Number=1,Type=String,Description="Genotype filter">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
##FORMAT=<ID=HQ,Number=2,Type=Integer,Description="Haplotype quality">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred likelihoods">
##FORMAT=<ID=GF,Number=1,Type=Float,Description="A float">
##FORMAT=<ID=GS,Number=1,Type=String,Description="A string">
##FORMAT=<ID=MI,Number=.,Type=Integer,Description="Unbounded ints">
##contig=<ID=20,length=62435964>
##contig=<ID=21,length=48129895>
##contig=<ID=X,length=155270560>
##reference=file:///ref.fa
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA1
20	10	.	A	C	5	PASS	DP=1.5	GT	0/1
""");
        cases.put("refuse-nocontig", """
##fileformat=VCFv4.2
##FILTER=<ID=LowQual,Description="Low quality">
##FILTER=<ID=q10,Description="Quality below 10">
##FILTER=<ID=s50,Description="Less than 50% of samples">
##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral allele">
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total alleles">
##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP membership">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">
##INFO=<ID=H2,Number=0,Type=Flag,Description="HapMap2">
##INFO=<ID=MQ,Number=1,Type=Float,Description="Mapping quality">
##INFO=<ID=NS,Number=1,Type=Integer,Description="Samples with data">
##INFO=<ID=BIG,Number=1,Type=Integer,Description="Large integer">
##INFO=<ID=MID,Number=.,Type=Integer,Description="Unbounded ints">
##INFO=<ID=STRS,Number=.,Type=String,Description="Strings">
##INFO=<ID=CH,Number=1,Type=Character,Description="Char">
##INFO=<ID=RF,Number=R,Type=Float,Description="Per allele floats">
##INFO=<ID=GI,Number=G,Type=Integer,Description="Per genotype ints">
##INFO=<ID=F2,Number=2,Type=Float,Description="Two floats">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
##FORMAT=<ID=FT,Number=1,Type=String,Description="Genotype filter">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
##FORMAT=<ID=HQ,Number=2,Type=Integer,Description="Haplotype quality">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred likelihoods">
##FORMAT=<ID=GF,Number=1,Type=Float,Description="A float">
##FORMAT=<ID=GS,Number=1,Type=String,Description="A string">
##FORMAT=<ID=MI,Number=.,Type=Integer,Description="Unbounded ints">
##contig=<ID=20,length=62435964>
##contig=<ID=21,length=48129895>
##contig=<ID=X,length=155270560>
##reference=file:///ref.fa
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA1
22	10	.	A	C	5	PASS	.	GT	0/1
""");
        cases.put("refuse-badfilter", """
##fileformat=VCFv4.2
##FILTER=<ID=LowQual,Description="Low quality">
##FILTER=<ID=q10,Description="Quality below 10">
##FILTER=<ID=s50,Description="Less than 50% of samples">
##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral allele">
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total alleles">
##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP membership">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">
##INFO=<ID=H2,Number=0,Type=Flag,Description="HapMap2">
##INFO=<ID=MQ,Number=1,Type=Float,Description="Mapping quality">
##INFO=<ID=NS,Number=1,Type=Integer,Description="Samples with data">
##INFO=<ID=BIG,Number=1,Type=Integer,Description="Large integer">
##INFO=<ID=MID,Number=.,Type=Integer,Description="Unbounded ints">
##INFO=<ID=STRS,Number=.,Type=String,Description="Strings">
##INFO=<ID=CH,Number=1,Type=Character,Description="Char">
##INFO=<ID=RF,Number=R,Type=Float,Description="Per allele floats">
##INFO=<ID=GI,Number=G,Type=Integer,Description="Per genotype ints">
##INFO=<ID=F2,Number=2,Type=Float,Description="Two floats">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
##FORMAT=<ID=FT,Number=1,Type=String,Description="Genotype filter">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
##FORMAT=<ID=HQ,Number=2,Type=Integer,Description="Haplotype quality">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred likelihoods">
##FORMAT=<ID=GF,Number=1,Type=Float,Description="A float">
##FORMAT=<ID=GS,Number=1,Type=String,Description="A string">
##FORMAT=<ID=MI,Number=.,Type=Integer,Description="Unbounded ints">
##contig=<ID=20,length=62435964>
##contig=<ID=21,length=48129895>
##contig=<ID=X,length=155270560>
##reference=file:///ref.fa
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA1
20	10	.	A	C	5	NOPE	.	GT	0/1
""");
        cases.put("refuse-flagfalse", """
##fileformat=VCFv4.2
##FILTER=<ID=LowQual,Description="Low quality">
##FILTER=<ID=q10,Description="Quality below 10">
##FILTER=<ID=s50,Description="Less than 50% of samples">
##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral allele">
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total alleles">
##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP membership">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">
##INFO=<ID=H2,Number=0,Type=Flag,Description="HapMap2">
##INFO=<ID=MQ,Number=1,Type=Float,Description="Mapping quality">
##INFO=<ID=NS,Number=1,Type=Integer,Description="Samples with data">
##INFO=<ID=BIG,Number=1,Type=Integer,Description="Large integer">
##INFO=<ID=MID,Number=.,Type=Integer,Description="Unbounded ints">
##INFO=<ID=STRS,Number=.,Type=String,Description="Strings">
##INFO=<ID=CH,Number=1,Type=Character,Description="Char">
##INFO=<ID=RF,Number=R,Type=Float,Description="Per allele floats">
##INFO=<ID=GI,Number=G,Type=Integer,Description="Per genotype ints">
##INFO=<ID=F2,Number=2,Type=Float,Description="Two floats">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
##FORMAT=<ID=FT,Number=1,Type=String,Description="Genotype filter">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
##FORMAT=<ID=HQ,Number=2,Type=Integer,Description="Haplotype quality">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred likelihoods">
##FORMAT=<ID=GF,Number=1,Type=Float,Description="A float">
##FORMAT=<ID=GS,Number=1,Type=String,Description="A string">
##FORMAT=<ID=MI,Number=.,Type=Integer,Description="Unbounded ints">
##contig=<ID=20,length=62435964>
##contig=<ID=21,length=48129895>
##contig=<ID=X,length=155270560>
##reference=file:///ref.fa
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA1
20	10	.	A	C	5	PASS	DB=false	GT	0/1
""");
        cases.put("refuse-nocontiglines", """
##fileformat=VCFv4.2
##FILTER=<ID=LowQual,Description="Low quality">
##FILTER=<ID=q10,Description="Quality below 10">
##FILTER=<ID=s50,Description="Less than 50% of samples">
##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral allele">
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total alleles">
##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP membership">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">
##INFO=<ID=H2,Number=0,Type=Flag,Description="HapMap2">
##INFO=<ID=MQ,Number=1,Type=Float,Description="Mapping quality">
##INFO=<ID=NS,Number=1,Type=Integer,Description="Samples with data">
##INFO=<ID=BIG,Number=1,Type=Integer,Description="Large integer">
##INFO=<ID=MID,Number=.,Type=Integer,Description="Unbounded ints">
##INFO=<ID=STRS,Number=.,Type=String,Description="Strings">
##INFO=<ID=CH,Number=1,Type=Character,Description="Char">
##INFO=<ID=RF,Number=R,Type=Float,Description="Per allele floats">
##INFO=<ID=GI,Number=G,Type=Integer,Description="Per genotype ints">
##INFO=<ID=F2,Number=2,Type=Float,Description="Two floats">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
##FORMAT=<ID=FT,Number=1,Type=String,Description="Genotype filter">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
##FORMAT=<ID=HQ,Number=2,Type=Integer,Description="Haplotype quality">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred likelihoods">
##FORMAT=<ID=GF,Number=1,Type=Float,Description="A float">
##FORMAT=<ID=GS,Number=1,Type=String,Description="A string">
##FORMAT=<ID=MI,Number=.,Type=Integer,Description="Unbounded ints">
##reference=file:///ref.fa
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA1	NA2	NA3
20	14370	rs6054257	G	A	29	PASS	NS=3;DP=14;AF=0.5;DB;H2;AA=G;AN=6;AC=3;MQ=60.25;BIG=70000;MID=1,300,70000;STRS=a,bb,ccc;CH=Z;RF=0.1,0.2;GI=1,2,3;F2=1.5,2.5	GT:GQ:DP:HQ:AD:PL:GF:GS:MI	0|0:48:1:51,51:1,0:0,10,100:0.5:x:1	1|0:48:8:51,51:3,4:10,0,100:1.25:yy:1,2	1/1:43:5:.,.:.:.:.:.:.
20	17330	.	T	A	3	q10	NS=3;DP=11;AF=0.017	GT:GQ:DP:HQ:FT	0|0:49:3:58,50:PASS	0|1:3:5:65,3:q10	0/0:41:3:.:q10;s50
20	1110696	rs6040355	A	G,T	67	PASS	NS=2;DP=10;AF=0.333,0.667;AA=T;DB;AC=1,2;GI=1,2,3,4,5,6	GT:GQ:DP:HQ:AD:PL	1|2:21:6:23,27:1,2,3:1,2,3,4,5,6	2|1:2:0:18,2:.:.	2/2:35:4:.:4,5,6:7,8,9,10,11,12
20	1230237	.	T	.	47	.	NS=3;DP=13;AA=T	GT:GQ:DP:HQ	0|0:54:7:56,60	0|0:48:4:51,51	0/0:61:2:.
20	1234567	microsat1	GTC	G,GTCT	50	PASS	NS=3;DP=9;AA=G;MID=5	GT:GQ:DP	0/1:35:4	0/2:17:2	1/1:40:3
20	1234568	.	G	GA,GAC,GACT	.	q10;s50	DP=200;BIG=-2147483647	GT:GQ:PL	0/3:120:0,1,2,3,4,5,6,7,8,9	.	1|3:.:.
21	100	.	A	C	0	LowQual	.	GT:DP	./.:.	0:3	.|1:2
21	200	id1;id2	ACGT	A	1e-5	PASS	DP=0;MQ=0.001;AF=1e-30	GT	0/1	1/1	0/0
X	10	.	C	CAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA	123.456	PASS	STRS=onlyone;CH=.;DP=.	GT:GS:GF	1:hello:.	0:.:3.14159	.:z:-0.0
X	20	.	A	<DEL>	99	PASS	.	GT:AD	0/1:5,5	0/0:10,0	1/1:0,10
""");

        for (final Map.Entry<String, String> c : cases.entrySet()) {
            final String name = c.getKey();
            final Path vcf = dir.resolve(name + ".vcf");
            Files.writeString(vcf, c.getValue(), StandardCharsets.UTF_8);
            emit("vcf", name, hex(Files.readAllBytes(vcf)));
            final Path bcf = dir.resolve(name + ".bcf");
            final String noIndex = "CREATE_INDEX=false";
            if (!run(name, "converter-bcf", bcf, new picard.vcf.VcfFormatConverter(),
                    "I=" + vcf, "O=" + bcf, "REQUIRE_INDEX=false", noIndex)) {
                continue;
            }
            final Path back = dir.resolve(name + ".back.vcf");
            run(name, "converter-vcf", back, new picard.vcf.VcfFormatConverter(),
                    "I=" + bcf, "O=" + back, "REQUIRE_INDEX=false", noIndex);
            final Path again = dir.resolve(name + ".again.bcf");
            run(name, "converter-bcf2", again, new picard.vcf.VcfFormatConverter(),
                    "I=" + bcf, "O=" + again, "REQUIRE_INDEX=false", noIndex);
            final Path sorted = dir.resolve(name + ".sorted.bcf");
            run(name, "sortvcf-bcf", sorted, new picard.vcf.SortVcf(),
                    "I=" + bcf, "O=" + sorted, noIndex);
            final Path merged = dir.resolve(name + ".merged.bcf");
            run(name, "mergevcfs-bcf", merged, new picard.vcf.MergeVcfs(),
                    "I=" + bcf, "O=" + merged, noIndex);
            final Path sites = dir.resolve(name + ".sites.bcf");
            run(name, "sitesonly-bcf", sites, new picard.vcf.MakeSitesOnlyVcf(),
                    "I=" + vcf, "O=" + sites, noIndex);
        }
        System.out.print(buf);
    }
}
