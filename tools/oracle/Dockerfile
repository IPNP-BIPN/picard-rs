# Oracle image for picard-rs golden generation.
#
# Same contract as htsjdk-rs's oracle: the base image is pinned by digest, the jar by sha256,
# the platform is linux/amd64, and the probe fails the build rather than letting a degraded
# environment produce a golden.
#
# The jar is Picard's release fat jar, which bundles htsjdk 4.2.0 and GKL as
# `build.gradle` at tag 3.4.0 pins them. Taking the fat jar rather than assembling a classpath
# is deliberate: the version skew between Picard and its htsjdk is exactly the kind of thing
# that would silently change output, and the release artefact is what users actually run.

FROM eclipse-temurin@sha256:068a8f9ae4b74d9a20de3ecca771ba1a6437f1d7a8a8ad6deaf9dbdd2274397a

ARG PICARD_VERSION=3.4.0

# R is not optional for a large part of the archetype. Ten of the 44 metrics tools declare
# a chart argument as REQUIRED (CollectInsertSizeMetrics's HISTOGRAM_FILE among them) and shell
# out to Rscript, so without R they refuse to run at all and no golden can be produced. The
# stratification in docs/decisions/0001 predicted this as an oracle-environment cost paid once
# for the repository; this is that payment.
# Retry the apt step: r-base-core is a large download from the Debian mirrors and the fetch
# fails intermittently (exit 100) with no bearing on the build's correctness, which turned a
# transient mirror hiccup into a red oracle job. Three attempts with a short backoff turns the
# flake back into the non-event it is.
RUN for attempt in 1 2 3; do \
      apt-get update \
      && apt-get install -y --no-install-recommends curl ca-certificates r-base-core \
      && break; \
      echo "apt attempt $attempt failed; retrying"; sleep 10; \
    done \
 && command -v Rscript >/dev/null \
 && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/oracle/jars

RUN set -eux; \
    curl -fsSL -o picard.jar \
      https://github.com/broadinstitute/picard/releases/download/${PICARD_VERSION}/picard.jar; \
    echo "e76128c283889fc583c9dea33a3b7448974c067d102c9e35be152642d4d5f901  picard.jar" \
      | sha256sum -c -

ENV ORACLE_CP=/opt/oracle/jars/picard.jar:/opt/oracle/classes
ENV PICARD_JAR=/opt/oracle/jars/picard.jar

WORKDIR /opt/oracle
COPY OracleProbe.java .
RUN mkdir -p classes && javac -cp "$ORACLE_CP" -d classes OracleProbe.java

# Fail the build itself if the environment does not satisfy the contract.
RUN java -cp "$ORACLE_CP:classes" OracleProbe

WORKDIR /work
ENTRYPOINT ["/bin/sh", "-c"]
CMD ["java -cp \"$ORACLE_CP:/opt/oracle/classes\" OracleProbe"]
