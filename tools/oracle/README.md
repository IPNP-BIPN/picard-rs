# The picard-rs oracle

A digest-pinned `linux/amd64` container holding Picard 3.4.0's release fat jar, used to
generate every golden in this repository.

## The contract

`OracleProbe.java` runs **during the image build**, so a degraded environment cannot produce an
image, let alone a golden. It exits 2 on any violation. Each check exists because its failure
mode is silent:

| checked | why |
|---|---|
| `os.arch = amd64` | Intel GKL ships x86-only natives |
| Java 17 | the pinned runtime |
| Picard reports `3.4.0` | a jar from the wrong tag is otherwise invisible |
| default locale is `en_US` | metrics number formatting is locale-dependent; see htsjdk-rs decision 0011 |
| a decimal formats as `0.333333` | the locale check passing is not the same as the formatting being right |
| `usingIntelDeflater` is true | the oracle pins the JDK deflater, and a GKL that cannot load means that pin is untested |
| `avx`, `avx2`, `sse4_2` | recorded, since PairHMM selection depends on them |

The locale check is not theoretical. Under `-Duser.language=fr` the same probe reports
`"decimal_sample": "0,333333"` and refuses to run.

## Why the fat jar

Picard's release artefact bundles htsjdk 4.2.0 and GKL exactly as `build.gradle` at tag 3.4.0
pins them. Assembling a classpath by hand would invite version skew between Picard and its
htsjdk, which is precisely the kind of thing that changes output without any error. The fat jar
is also what users actually run.

## Note on the entry point

The program's bit-identity claim is defined against `gatk <Tool>`, because that fixes the
argument parser to Barclay (`--INPUT`) rather than Picard's legacy syntax (`INPUT=`). This image
runs `picard.jar` directly, with legacy syntax.

For **metrics computation** the two are the same code path; what differs is parsing and the
command-line string recorded in the file header, which is canonicalized. That claim is stated
here so it can be checked rather than assumed, and confirming it against a GATK-entry oracle is
tracked work, not a settled question.

## Usage

```sh
docker build --platform linux/amd64 -t picard-rs-oracle:3.4.0 tools/oracle
docker run --rm --platform linux/amd64 picard-rs-oracle:3.4.0     # prints the provenance record
```
