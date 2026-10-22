import com.intel.gkl.compression.IntelDeflaterFactory;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.security.MessageDigest;
import java.util.Locale;
import java.util.zip.Deflater;

/**
 * GKL's deflater, the C library Picard 3.5.0 loads, on BGZF-sized blocks of a payload: one thread,
 * level 5, raw deflate as BlockCompressedOutputStream asks for it. Prints the rate and a sha256 of
 * every block's bytes in order, which the port's side prints too.
 */
public class GklDeflateBench {
  public static void main(String[] args) throws Exception {
    byte[] data = Files.readAllBytes(Paths.get(args[0]));
    int reps = Integer.parseInt(args[1]);
    int block = 65498;
    Deflater deflater = new IntelDeflaterFactory().makeDeflater(5, true);
    System.out.println("deflater " + deflater.getClass().getName());
    byte[] out = new byte[block + 1024];
    for (int rep = 0; rep <= reps; rep++) {
      MessageDigest digest = MessageDigest.getInstance("SHA-256");
      long compressed = 0;
      long start = System.nanoTime();
      for (int at = 0; at < data.length; at += block) {
        int n = Math.min(block, data.length - at);
        deflater.reset();
        deflater.setInput(data, at, n);
        deflater.finish();
        int got = 0;
        while (!deflater.finished()) got += deflater.deflate(out, got, out.length - got);
        digest.update(out, 0, got);
        compressed += got;
      }
      double seconds = (System.nanoTime() - start) / 1e9;
      StringBuilder hex = new StringBuilder();
      for (byte b : digest.digest()) hex.append(String.format("%02x", b));
      // The first pass warms the JIT and is not reported.
      if (rep > 0)
        System.out.printf(Locale.ROOT, "gkl_c %.1f MB/s ratio %.4f sha256 %s%n",
            data.length / 1e6 / seconds, (double) compressed / data.length, hex);
    }
  }
}
