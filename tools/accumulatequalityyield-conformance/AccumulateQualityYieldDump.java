/*
 * Oracle dump for AccumulateQualityYieldMetrics against Picard 3.4.0. Writes two QualityYieldMetrics
 * files (built programmatically so their format is exactly Picard's), runs
 * AccumulateQualityYieldMetrics over them, and emits the two inputs and the combined output.
 *   java -cp picard-fat.jar:. AccumulateQualityYieldDump | gzip -n > accumulate_quality_yield.txt.gz
 */
import java.io.*; import java.nio.file.*;
import htsjdk.samtools.metrics.MetricsFile;
import picard.analysis.CollectQualityYieldMetrics.QualityYieldMetrics;
public class AccumulateQualityYieldDump {
  static StringBuilder buf = new StringBuilder();
  static void emit(String k, String c, String p) {
    buf.append(k).append('\t').append(c).append('\t')
       .append(p.replace("\\","\\\\").replace("\t","\\t").replace("\n","\\n")).append('\n');
  }
  static File dir;
  static File writeInput(String name, long tr, long pfr, long tb, long pfb,
                         long q20, long pfq20, long q30, long pfq30, long q20y, long pfq20y) {
    QualityYieldMetrics m = new QualityYieldMetrics(false);
    m.TOTAL_READS = tr; m.PF_READS = pfr; m.TOTAL_BASES = tb; m.PF_BASES = pfb;
    m.Q20_BASES = q20; m.PF_Q20_BASES = pfq20; m.Q30_BASES = q30; m.PF_Q30_BASES = pfq30;
    m.Q20_EQUIVALENT_YIELD = q20y; m.PF_Q20_EQUIVALENT_YIELD = pfq20y;
    m.calculateDerivedFields();
    MetricsFile<QualityYieldMetrics, ?> mf = new MetricsFile<>();
    mf.addMetric(m);
    File f = new File(dir, name);
    mf.write(f);
    return f;
  }
  public static void main(String[] x) throws Exception {
    dir = Files.createTempDirectory("aqy").toFile();
    File f1 = writeInput("shard1.quality_yield_metrics",
        100, 90, 15000, 13500, 14000, 12800, 12000, 11000, 16000, 15000);
    File f2 = writeInput("shard2.quality_yield_metrics",
        50, 48, 7600, 7300, 7100, 6900, 6000, 5800, 8000, 7700);
    File out = new File(dir, "combined.quality_yield_metrics");
    int rc = new picard.util.AccumulateQualityYieldMetrics().instanceMain(new String[]{
      "INPUT=" + f1.getAbsolutePath(), "INPUT=" + f2.getAbsolutePath(),
      "OUTPUT=" + out.getAbsolutePath()});
    emit("input1", "two_shards", new String(Files.readAllBytes(f1.toPath())));
    emit("input2", "two_shards", new String(Files.readAllBytes(f2.toPath())));
    emit("rc", "two_shards", String.valueOf(rc));
    emit("output", "two_shards", new String(Files.readAllBytes(out.toPath())));
    System.out.print(buf);
  }
}
