# Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

FROM fedora:35 AS base

ARG CONTAINER_VERSION
ARG CONTAINER_TITLE
ARG STRIPENN_VER=${CONTAINER_VERSION}

ARG PIP_NO_CACHE_DIR=0

RUN if [ -z "$CONTAINER_VERSION" ]; then echo "Missing CONTAINER_VERSION --build-arg" && exit 1; fi

RUN echo 'max_parallel_downloads=10' >> /etc/dnf/dnf.conf \
&&  echo 'fastestmirror=True' >> /etc/dnf/dnf.conf

RUN dnf update -y \
&&  dnf install -y gcc \
                   python3-joblib \
                   python3-matplotlib \
                   python3-numpy \
                   python3-opencv \
                   python3-pip \
                   python3-pandas \
                   python3-scikit-image \
                   python3-scipy \
                   python3-tqdm \
                   python3-typer \
                   zlib-devel \
&&  strip --remove-section=.note.ABI-tag /usr/lib64/libQt5Core.so.5 \
&&  ldconfig \
&&  pip3 install cooler \
&&  pip3 install "stripenn==${STRIPENN_VER}" --no-deps \
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

LABEL org.opencontainers.image.authors='Roberto Rossini <roberros@uio.no>'
LABEL org.opencontainers.image.url='https://github.com/paulsengroup/2021-modle-paper-001-data-analysis'
LABEL org.opencontainers.image.documentation='https://github.com/2021-modle-paper-001-data-analysis'
LABEL org.opencontainers.image.source='https://github.com/paulsengroup/2021-modle-paper-001-data-analysis'
LABEL org.opencontainers.image.licenses='MIT'
LABEL org.opencontainers.image.title="${CONTAINER_TITLE:-stripenn}"
LABEL org.opencontainers.image.version="${CONTAINER_VERSION:-latest}"
