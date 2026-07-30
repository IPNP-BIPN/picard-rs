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
    curl -fsSL --retry 8 --retry-all-errors --retry-delay 5 -o picard.jar \
      https://github.com/broadinstitute/picard/releases/download/${PICARD_VERSION}/picard.jar; \
    echo "e76128c283889fc583c9dea33a3b7448974c067d102c9e35be152642d4d5f901  picard.jar" \
      | sha256sum -c -

# JaCoCo, for the coverage-guided differential fuzzer. Pinned by sha256 like the Picard jar: the
# fuzzer's stopping condition is a branch-coverage threshold on the reference, so a silently
# different JaCoCo would silently change what "covered" means.
#
# The zip is taken rather than the individual jars because it carries the agent and the CLI (which
# bundles jacoco-core, needed to analyse the exec data in process) as one digest.
RUN set -eux; \
    curl -fsSL --retry 8 --retry-all-errors --retry-delay 5 -o /tmp/jacoco.zip \
      https://repo1.maven.org/maven2/org/jacoco/jacoco/0.8.12/jacoco-0.8.12.zip; \
    echo "35cf2a2b8937659ecbc8d0d3902dfa7f55b7ed2250e82424036a3e1d12462cd7  /tmp/jacoco.zip" \
      | sha256sum -c -; \
    mkdir -p /opt/jacoco; \
    cd /opt/jacoco && jar xf /tmp/jacoco.zip lib/jacocoagent.jar lib/jacococli.jar; \
    rm /tmp/jacoco.zip

# jacococli bundles only the CLI, so the analysis classes and their asm dependency are taken
# separately, each pinned by sha256. They are what turns the agent's execution data into a set of
# covered branches inside the fuzzer's own JVM, without shelling out per iteration.
RUN set -eux; \
    cd /opt/jacoco/lib; \
    for spec in \
      "org/jacoco/org.jacoco.core/0.8.12/org.jacoco.core-0.8.12.jar fca26db37c0c5fbd5dc4985237eb82866df9799d5082af899475a73f91f5b035" \
      "org/ow2/asm/asm/9.7/asm-9.7.jar adf46d5e34940bdf148ecdd26a9ee8eea94496a72034ff7141066b3eea5c4e9d" \
      "org/ow2/asm/asm-commons/9.7/asm-commons-9.7.jar 389bc247958e049fc9a0408d398c92c6d370c18035120395d4cba1d9d9304b7a" \
      "org/ow2/asm/asm-tree/9.7/asm-tree-9.7.jar 62f4b3bc436045c1acb5c3ba2d8ec556ec3369093d7f5d06c747eb04b56d52b1" \
    ; do \
      set -- $spec; \
      curl -fsSL --retry 8 --retry-all-errors --retry-delay 5 -O "https://repo1.maven.org/maven2/$1"; \
      echo "$2  $(basename $1)" | sha256sum -c -; \
    done

ENV JACOCO_AGENT=/opt/jacoco/lib/jacocoagent.jar
ENV JACOCO_CLI=/opt/jacoco/lib/jacococli.jar
ENV JACOCO_ANALYSIS=/opt/jacoco/lib/org.jacoco.core-0.8.12.jar:/opt/jacoco/lib/asm-9.7.jar:/opt/jacoco/lib/asm-commons-9.7.jar:/opt/jacoco/lib/asm-tree-9.7.jar
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
