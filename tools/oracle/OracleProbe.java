/*
 * Asserts that this container satisfies picard-rs's oracle contract, and records what it found.
 *
 * The principle, inherited from htsjdk-rs's probe: a degraded environment must FAIL rather than
 * produce a golden that looks exactly like a good one. Every check here exists because its
 * failure mode is silent.
 *
 * Exits 2 on violation, so the Docker build itself cannot complete in a bad environment.
 */

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Locale;
import java.util.TreeSet;

public class OracleProbe {

    private static final String EXPECTED_ARCH = "amd64";
    private static final String EXPECTED_JAVA_MAJOR = "17";
    private static final String EXPECTED_PICARD_VERSION = "3.4.0";

    /**
     * htsjdk-rs decision 0011: FormatUtil reaches NumberFormat.getNumberInstance(), which takes
     * the default locale, and nothing in Picard pins it. Under fr-FR a metrics file has commas
     * for decimal points; under ar-EG it has Eastern Arabic-Indic digits. So the locale is part
     * of the contract.
     */
    private static final String EXPECTED_LOCALE = "en_US";

    private static final String[] REQUIRED_CPU_FLAGS = {"avx", "avx2", "sse4_2"};

    public static void main(final String[] args) throws Exception {
        final List<String> failures = new ArrayList<>();

        final String arch = System.getProperty("os.arch");
        final String javaVersion = System.getProperty("java.version");
        final String javaVendor = System.getProperty("java.vendor");
        final String javaMajor = javaVersion.split("\\.")[0];
        final Locale locale = Locale.getDefault();
        final String decimalSample = new java.text.DecimalFormat("0.0#####").format(1.0 / 3.0);
        final TreeSet<String> cpuFlags = readCpuFlags();

        // The version Picard reports about itself, not the one we asked the build for. A jar
        // fetched from the wrong tag would otherwise be invisible.
        String picardVersion = "unknown";
        try {
            picardVersion = Class.forName("picard.cmdline.CommandLineProgram")
                    .getPackage().getImplementationVersion();
            if (picardVersion == null) picardVersion = "absent";
        } catch (final Throwable t) {
            failures.add("picard classes are not on the classpath: " + t);
        }

        // Picard installs Intel GKL as the default deflater factory unless USE_JDK_DEFLATER is
        // set. The oracle contract pins the JDK deflater, so what matters here is that GKL is
        // *present and working*, since its silent absence would mean the pin is doing nothing
        // and a future unpinned golden would be wrong without warning.
        boolean gklPresent = false;
        try {
            final Class<?> f = Class.forName("com.intel.gkl.compression.IntelDeflaterFactory");
            gklPresent = (Boolean) f.getMethod("usingIntelDeflater").invoke(f.getConstructor().newInstance());
        } catch (final Throwable t) {
            failures.add("Intel GKL is not usable: " + t
                    + ". It degrades silently, which is why this is checked.");
        }

        if (!EXPECTED_ARCH.equals(arch)) {
            failures.add("os.arch is '" + arch + "', expected '" + EXPECTED_ARCH + "'");
        }
        if (!EXPECTED_JAVA_MAJOR.equals(javaMajor)) {
            failures.add("java major is '" + javaMajor + "', expected '" + EXPECTED_JAVA_MAJOR + "'");
        }
        if (!EXPECTED_PICARD_VERSION.equals(picardVersion)) {
            failures.add("picard reports version '" + picardVersion + "', expected '"
                    + EXPECTED_PICARD_VERSION + "'");
        }
        if (!EXPECTED_LOCALE.equals(locale.toString())) {
            failures.add("default locale is '" + locale + "', expected '" + EXPECTED_LOCALE
                    + "'. Metrics number formatting is locale-dependent; see htsjdk-rs"
                    + " decision 0011.");
        }
        if (!"0.333333".equals(decimalSample)) {
            failures.add("a decimal formats as '" + decimalSample + "', expected '0.333333'.");
        }
        if (!gklPresent) {
            failures.add("usingIntelDeflater is false. The oracle pins the JDK deflater, but a"
                    + " GKL that cannot load means that pin is untested.");
        }
        // Ten of the 44 metrics tools require a chart argument and shell out to Rscript. A
        // missing R means those tools refuse to run, which is loud; a *different* R could in
        // principle change a chart, but never the metrics file, so only presence is checked.
        String rVersion = "absent";
        try {
            final Process p = new ProcessBuilder("Rscript", "--version").redirectErrorStream(true).start();
            rVersion = new String(p.getInputStream().readAllBytes()).trim();
            p.waitFor();
        } catch (final Exception e) {
            failures.add("Rscript is not available: " + e
                    + ". Tools with a required chart argument cannot produce a golden.");
        }

        for (final String flag : REQUIRED_CPU_FLAGS) {
            if (!cpuFlags.isEmpty() && !cpuFlags.contains(flag)) {
                failures.add("CPU flag '" + flag + "' is absent");
            }
        }

        final StringBuilder json = new StringBuilder();
        json.append("{\n");
        json.append("  \"os_arch\": \"").append(arch).append("\",\n");
        json.append("  \"java_version\": \"").append(javaVersion).append("\",\n");
        json.append("  \"java_vendor\": \"").append(javaVendor).append("\",\n");
        json.append("  \"picard_version\": \"").append(picardVersion).append("\",\n");
        json.append("  \"default_locale\": \"").append(locale).append("\",\n");
        json.append("  \"decimal_sample\": \"").append(decimalSample).append("\",\n");
        json.append("  \"rscript\": \"").append(rVersion.replace("\"", "'").replace("\n", " ")).append("\",\n");
        json.append("  \"using_intel_deflater\": ").append(gklPresent).append(",\n");
        json.append("  \"avx\": ").append(cpuFlags.contains("avx")).append(",\n");
        json.append("  \"avx2\": ").append(cpuFlags.contains("avx2")).append(",\n");
        json.append("  \"avx512f\": ").append(cpuFlags.contains("avx512f")).append(",\n");
        json.append("  \"contract_satisfied\": ").append(failures.isEmpty()).append("\n");
        json.append("}");

        if (!failures.isEmpty()) {
            System.err.println("ORACLE CONTRACT VIOLATED. No golden produced here may be trusted.");
            System.err.println(json);
            for (final String f : failures) System.err.println("  - " + f);
            System.exit(2);
        }
        System.out.println(json);
    }

    private static TreeSet<String> readCpuFlags() {
        final TreeSet<String> flags = new TreeSet<>();
        try (BufferedReader r = new BufferedReader(new FileReader("/proc/cpuinfo"))) {
            String line;
            while ((line = r.readLine()) != null) {
                if (line.startsWith("flags")) {
                    flags.addAll(Arrays.asList(line.split(":", 2)[1].trim().split("\\s+")));
                    break;
                }
            }
        } catch (final Exception ignored) {
            // Absent /proc/cpuinfo is not itself a violation; the flag checks then no-op.
        }
        return flags;
    }
}
