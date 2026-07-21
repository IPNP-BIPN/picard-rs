/*
 * Probe: does GcBiasUtils.calculateGc count a lowercase 'n' differently depending on which branch
 * of the sliding window it takes?
 *
 * The initialising branch counts no-calls with SequenceUtil.basesEqual(base, 'N'), which is
 * case-insensitive. The incremental branch counts the *incoming* base with a raw byte comparison:
 *
 *     else if (newBase == 'N') ++state.nCount;
 *
 * while still *decrementing* on the outgoing base with the case-insensitive test:
 *
 *     else if (SequenceUtil.basesEqual(state.priorBase, (byte)'N')) --state.nCount;
 *
 * So over a sequence containing lowercase 'n', nCount is decremented for bases that were never
 * counted, and can go negative. A window that should be rejected for having too many no-calls is
 * then accepted, and its GC value enters the histogram.
 *
 * calculateRefWindowsByGc uppercases the reference first, so it is safe. calculateAllGcs does
 * not, and takes whatever the caller passes.
 *
 * Prints the GC bin per window for an uppercase and a lowercase version of the same sequence.
 * If the two agree, the asymmetry is unreachable and this record is wrong.
 */

import htsjdk.samtools.util.StringUtil;
import picard.analysis.GcBiasUtils;
import java.lang.reflect.*;
import java.util.Arrays;

public class GcAsymmetryProbe {

    static byte[] gcs(String seq, int windowSize) throws Exception {
        byte[] bases = seq.getBytes();
        // calculateAllGcs(refBases, lastWindowStart, windowSize)
        Method m = GcBiasUtils.class.getMethod("calculateAllGcs", byte[].class, int.class, int.class);
        return (byte[]) m.invoke(null, bases, bases.length - windowSize, windowSize);
    }

    public static void main(String[] args) throws Exception {
        final int w = 10;
        // The no-calls must *enter* the window, not leave it: the buggy test is on the incoming
        // base. So the first window is clean and the no-calls sit after it.
        //
        // A first attempt put them at the start, where they only ever left the window, and the
        // two versions agreed. That probe proved nothing and is corrected here.
        String upper = "GCGCGCGCGC" + "GNNNNNNNGC" + "GCGCGCGCGC" + "ATATATATAT";
        String lower = upper.replace('N', 'n');

        byte[] u = gcs(upper, w);
        byte[] l = gcs(lower, w);

        System.out.println("upper: " + Arrays.toString(Arrays.copyOf(u, 22)));
        System.out.println("lower: " + Arrays.toString(Arrays.copyOf(l, 22)));
        System.out.println("CASE_ASYMMETRY=" + !Arrays.equals(u, l));

        // And the negative-count consequence, stated as a value: how many windows are rejected
        // (-1) in each version.
        int ru = 0, rl = 0;
        for (int i = 1; i < upper.length() - w; i++) { if (u[i] == -1) ru++; if (l[i] == -1) rl++; }
        System.out.println("rejected_windows_upper=" + ru);
        System.out.println("rejected_windows_lower=" + rl);
    }
}
