# The speed image: Picard 3.5.0's release jar on the oracle's own JDK base, for the host's
# architecture. tools/oracle pins linux/amd64 because goldens must come from x86-64 (GKL); timing
# must not, because an emulated JVM measured against a native binary says nothing (decision 0004).
# On arm64 GKL has no native library and Picard falls back to the JDK's Deflater/Inflater, which is
# what a Picard user on that machine gets. Same base digest, same jar digest as tools/oracle.
FROM eclipse-temurin@sha256:068a8f9ae4b74d9a20de3ecca771ba1a6437f1d7a8a8ad6deaf9dbdd2274397a

ARG PICARD_VERSION=3.5.0

# R for the chart arguments some tools require; python3 runs timer.py.
RUN for attempt in 1 2 3; do \
      apt-get update \
      && apt-get install -y --no-install-recommends curl ca-certificates r-base-core python3 \
      && break; \
      echo "apt attempt $attempt failed; retrying"; sleep 10; \
    done \
 && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/oracle/jars
RUN set -eux; \
    curl -fsSL --retry 8 --retry-all-errors --retry-delay 5 -o picard.jar \
      https://github.com/broadinstitute/picard/releases/download/${PICARD_VERSION}/picard.jar; \
    echo "b7d97861c3a54ba5a421f5a317f38382f955803862d30ef4aca2bcdc54943631  picard.jar" \
      | sha256sum -c -

ENV PICARD_JAR=/opt/oracle/jars/picard.jar
WORKDIR /work
ENTRYPOINT ["/bin/sh", "-c"]
