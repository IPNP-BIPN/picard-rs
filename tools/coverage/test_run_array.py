"""Checks of how run_array.py reads an answer, on texts the reference printed. Run directly."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import run_array  # noqa: E402

# CheckFingerprint returning EXIT_CODE_WHEN_EXPECTED_SAMPLE_NOT_FOUND ran to completion: the
# answer is its log line, not the banners around it, which carry the clock.
RAN = """14:11:28.063 INFO  NativeLibraryLoader - Loading libgkl_compression.so from jar:file:/x
[Tue Nov 24 14:11:28 UTC 2026] CheckFingerprint --INPUT /work/fixtures/fp_a.bam
[Tue Nov  4 14:11:28 UTC 2026] Executing as root@3ace3351deec on Linux
WARNING\t2026-11-24 14:11:28\tCheckFingerprint\tSample sampleB where fingerprint was expected not found in x.vcf
[Tue Nov 24 14:11:28 UTC 2026] picard.fingerprint.CheckFingerprint done. Elapsed time: 0.01 minutes.
Runtime.totalMemory()=132120576
To get help, see http://broadinstitute.github.io/picard/index.html#GettingHelp
"""
assert run_array.first_error(RAN) == (
    "WARNING\tCheckFingerprint\tSample sampleB where fingerprint was expected not found in x.vcf"
), run_array.first_error(RAN)

# A refusal that prints only its message is untouched.
assert run_array.first_error("INPUT= is required\n") == "INPUT= is required"
# The data types IlluminaDataProviderFactory has no file for are a HashSet's, compared as a set.
assert run_array.outcome(1, "", "picard.PicardException: Could not find a format with available "
                         "files for the following data types: QualityScores, PF, BaseCalls", "T") == (
    "EXIT=1 picard.PicardException: Could not find a format with available files for the "
    "following data types: BaseCalls, PF, QualityScores")
print("run_array answers read as expected")
