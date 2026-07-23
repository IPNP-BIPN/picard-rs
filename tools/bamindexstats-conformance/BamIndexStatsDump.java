import htsjdk.samtools.*;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.zip.DeflaterFactory;
import java.io.*; import java.nio.file.*; import java.util.*;
public class BamIndexStatsDump {
  static StringBuilder buf=new StringBuilder();
  static void emit(String k,String c,String p){buf.append(k).append('\t').append(c).append('\t')
    .append(p.replace("\\","\\\\").replace("\t","\\t").replace("\n","\\n")).append('\n');}
  static String hex(byte[] b){StringBuilder s=new StringBuilder();for(byte x:b)s.append(String.format("%02x",x));return s.toString();}
  static SAMRecord rd(SAMFileHeader h,String n,int ref,int start){
    SAMRecord r=new SAMRecord(h); r.setReadName(n);
    if(ref>=0){r.setReferenceIndex(ref);r.setAlignmentStart(start);r.setMappingQuality(60);r.setCigarString("10M");
      r.setReadBases("AAAAAAAAAA".getBytes());byte[] q=new byte[10];Arrays.fill(q,(byte)30);r.setBaseQualities(q);}
    else{r.setReadUnmappedFlag(true);r.setReadBases("AAAA".getBytes());byte[] q=new byte[4];Arrays.fill(q,(byte)30);r.setBaseQualities(q);}
    return r;
  }
  static void run(String c,int[] lens,List<SAMRecord> recs) throws Exception {
    File dir=Files.createTempDirectory("bis").toFile();
    SAMFileHeader h=new SAMFileHeader(); SAMSequenceDictionary d=new SAMSequenceDictionary();
    for(int i=0;i<lens.length;i++) d.addSequence(new SAMSequenceRecord("chr"+(i+1),lens[i]));
    h.setSequenceDictionary(d); h.setSortOrder(SAMFileHeader.SortOrder.coordinate);
    File bam=new File(dir,"in.bam");
    SAMFileWriter w=new SAMFileWriterFactory().setCreateIndex(true).setUseAsyncIo(false).makeBAMWriter(h,true,bam);
    for(SAMRecord r: recs) w.addAlignment(r); w.close();
    ByteArrayOutputStream bout=new ByteArrayOutputStream(); PrintStream old=System.out; System.setOut(new PrintStream(bout));
    new picard.sam.BamIndexStats().instanceMain(new String[]{"I="+bam.getAbsolutePath()});
    System.setOut(old);
    emit("bam",c,hex(Files.readAllBytes(bam.toPath())));
    emit("stats",c,bout.toString());
  }
  public static void main(String[] x) throws Exception {
    BlockCompressedOutputStream.setDefaultDeflaterFactory(new DeflaterFactory());
    SAMFileHeader h1=new SAMFileHeader(); // just to build records against
    // single ref, two mapped + two unplaced
    {SAMFileHeader h=new SAMFileHeader();SAMSequenceDictionary d=new SAMSequenceDictionary();
     d.addSequence(new SAMSequenceRecord("chr1",100000));h.setSequenceDictionary(d);h.setSortOrder(SAMFileHeader.SortOrder.coordinate);
     run("basic",new int[]{100000},List.of(rd(h,"a",0,10),rd(h,"b",0,500),rd(h,"u1",-1,0),rd(h,"u2",-1,0)));}
    // two refs, second has NO reads -> the null-metadata line
    {SAMFileHeader h=new SAMFileHeader();SAMSequenceDictionary d=new SAMSequenceDictionary();
     d.addSequence(new SAMSequenceRecord("chr1",100000));d.addSequence(new SAMSequenceRecord("chr2",50000));
     h.setSequenceDictionary(d);h.setSortOrder(SAMFileHeader.SortOrder.coordinate);
     run("empty_second_ref",new int[]{100000,50000},List.of(rd(h,"a",0,10),rd(h,"b",0,500)));}
    System.out.print(buf);
  }
}
