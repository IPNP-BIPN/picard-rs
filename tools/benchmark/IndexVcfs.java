/** Tribble linear indexes (.idx) beside plain VCFs, as the corpus writes them: `IndexVcfs <vcf> ...`. */
public class IndexVcfs {
    public static void main(final String[] args) throws Exception {
        for (final String vcf : args) {
            htsjdk.tribble.index.IndexFactory.writeIndex(
                    htsjdk.tribble.index.IndexFactory.createLinearIndex(
                            new java.io.File(vcf), new htsjdk.variant.vcf.VCFCodec()),
                    new java.io.File(vcf + ".idx"));
        }
    }
}
