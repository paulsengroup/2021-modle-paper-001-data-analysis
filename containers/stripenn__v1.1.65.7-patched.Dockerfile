# Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT


FROM ubuntu:22.04 AS downloader

ARG CONTAINER_VERSION
ARG STRIPENN_VER=${CONTAINER_VERSION}


ARG URL="https://files.pythonhosted.org/packages/2f/e0/aaa97b2c84b24f8340582c1acb93e475188216dd726244fc68e8d87a922a/stripenn-1.1.65.7.tar.gz"

RUN apt-get update \
&& apt-get install -y curl diffutils patch tar \
&& cd /tmp \
&& curl -L "$URL" | tar -xzf -

COPY containers/patches/stripenn*.patch /tmp

RUN cd /tmp \
&& patch -p1 \
         --ignore-whitespace \
         --fuzz 3 < stripenn-*.fix_prng_seed.patch \
&& touch stripenn-1.1.65.7/README.md


FROM fedora:36 AS base

ARG CONTAINER_VERSION
ARG CONTAINER_TITLE
ARG STRIPENN_VER=${CONTAINER_VERSION}

ARG PIP_NO_CACHE_DIR=0

RUN if [ -z "$CONTAINER_VERSION" ]; then echo "Missing CONTAINER_VERSION --build-arg" && exit 1; fi

RUN echo 'max_parallel_downloads=10' >> /etc/dnf/dnf.conf \
&&  echo 'fastestmirror=True' >> /etc/dnf/dnf.conf

COPY --from=downloader /tmp/stripenn-*/README.md /tmp/stripenn/
COPY --from=downloader /tmp/stripenn-*/setup.* /tmp/stripenn/
COPY --from=downloader /tmp/stripenn-*/src /tmp/stripenn/src

RUN dnf update -y \
&&  dnf install -y gcc \
                   python3-joblib \
                   python3-matplotlib \
                   python3-numpy \
                   python3-opencv \
                   python3-pip \
                   python3-pandas \
                   python3-psutil \
                   python3-scikit-image \
                   python3-scipy \
                   python3-tqdm \
                   python3-typer \
                   zlib-devel \
&&  strip --remove-section=.note.ABI-tag /usr/lib64/libQt5Core.so.5 \
&&  ldconfig \
&&  pip3 install cooler \
&&  pip3 install /tmp/stripenn/ --no-deps \
&&  dnf remove -y gcc \
                  python3-pip \
                  zlib-devel \
&&  dnf autoremove -y \
&&  dnf clean all

# The strip command serve as a workaround for libQt5Core.so.5 failing to load when importing cv2
# This seems to occur only with Singularity, and only if the image was obtained with singularity pull

WORKDIR /data
ENTRYPOINT ["/usr/local/bin/stripenn"]

RUN python3 -c "import cv2"
RUN stripenn --help
RUN stripenn compute --help

LABEL org.opencontainers.image.authors='Roberto Rossini <roberros@uio.no>'
LABEL org.opencontainers.image.url='https://github.com/paulsengroup/2021-modle-paper-001-data-analysis'
LABEL org.opencontainers.image.documentation='https://github.com/2021-modle-paper-001-data-analysis'
LABEL org.opencontainers.image.source='https://github.com/paulsengroup/2021-modle-paper-001-data-analysis'
LABEL org.opencontainers.image.licenses='MIT'
LABEL org.opencontainers.image.title="${CONTAINER_TITLE:-stripenn}"
LABEL org.opencontainers.image.version="${CONTAINER_VERSION:-latest}"
