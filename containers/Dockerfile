# Dockerfile for llmap container image. Used by build_llmap_container.sh to
# build a docker image first (where apt has working DNS), then convert that
# image to a Singularity/Apptainer .sif via `singularity build` against the
# docker-daemon. This works around the fact that singularity 3.11's build
# `%post` runs in an isolated net namespace that cannot resolve DNS.
#
# Output: localhost/llmap:1.0.0 (docker image), then llmap.sif (apptainer).

FROM ubuntu:24.04

LABEL org.opencontainers.image.title="LLmap" \
      org.opencontainers.image.version="1.0.0" \
      org.opencontainers.image.description="Lossless. LLM-Augmented. Wave-Particle. Sequence mapper." \
      org.opencontainers.image.authors="schlein.lab@gmail.com"

ENV DEBIAN_FRONTEND=noninteractive

# gcc-13 / libstdc++ from Ubuntu 24.04 = GLIBCXX_3.4.32 + GLIBC 2.39.
RUN apt-get update && apt-get install -y --no-install-recommends \
        gcc-13 g++-13 \
        cmake ninja-build make \
        git ca-certificates \
        pkg-config \
        python3 python3-pip python3-venv \
        libcurl4-openssl-dev libssl-dev libzstd-dev \
        zlib1g-dev libbz2-dev liblzma-dev \
        libdeflate-dev libomp-dev \
        wget curl less file procps \
        samtools minimap2 bwa tabix \
    && rm -rf /var/lib/apt/lists/* \
    && update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-13 100 \
         --slave /usr/bin/g++ g++ /usr/bin/g++-13

# Bundle ONNX Runtime libs/headers + the prebuilt llmap binary.
COPY third_party/onnxruntime/lib     /opt/onnxruntime/lib
COPY third_party/onnxruntime/include /opt/onnxruntime/include
COPY build/src/llmap                 /opt/llmap/bin/llmap

# Wire onnxruntime into the dynamic linker path.
RUN echo "/opt/onnxruntime/lib" > /etc/ld.so.conf.d/onnxruntime.conf \
    && ldconfig

ENV LC_ALL=C.UTF-8 \
    LANG=C.UTF-8 \
    LD_LIBRARY_PATH=/opt/onnxruntime/lib \
    PATH=/opt/llmap/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin

# Smoke test at build time.
RUN /opt/llmap/bin/llmap --version \
    && minimap2 --version \
    && samtools --version | head -1

ENTRYPOINT ["/opt/llmap/bin/llmap"]
