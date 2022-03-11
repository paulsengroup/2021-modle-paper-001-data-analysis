# Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

FROM ubuntu:20.04 AS base

ARG CONTAINER_VERSION
ARG CONTAINER_TITLE
ARG STRIPENN_VER=${CONTAINER_VERSION}
ARG PIP_NO_CACHE_DIR=0


RUN apt-get update \
&&  apt-get install -y gcc \
                   python3 \
                   python3-dev \
                   python3-pip \
                   python3-setuptools \
                   libglvnd0 \
                   libglvnd-dev \
                   zlib1g \
                   zlib1g-dev \
&&  pip3 install stripenn==${STRIPENN_VER} --no-deps \
&&  pip3 check | cut -d ' ' -f 4 | sed 's/,//' | grep -v 'pathlib' \
    > /tmp/fix-deps.txt \
&&  echo psutil >> /tmp/fix-deps.txt \
&&  pip install -r /tmp/fix-deps.txt \
&&  pip install psutil \
&&  apt-get remove -y gcc \
                   python3-dev \
                   python3-pip \
                   python3-setuptools \
                   libglvnd-dev \
                   zlib1g-dev \
&&  apt-get autoremove -y \
&&  rm -rf /var/lib/apt/lists/*

WORKDIR /data
ENTRYPOINT ["/usr/local/bin/stripenn"]

RUN stripenn --help

LABEL org.opencontainers.image.authors='Roberto Rossini <roberros@uio.no>'
LABEL org.opencontainers.image.url='https://github.com/paulsengroup/2021-modle-paper-001-data-analysis'
LABEL org.opencontainers.image.documentation='https://github.com/2021-modle-paper-001-data-analysis'
LABEL org.opencontainers.image.source='https://github.com/paulsengroup/2021-modle-paper-001-data-analysis'
LABEL org.opencontainers.image.licenses='MIT'
LABEL org.opencontainers.image.title="${CONTAINER_TITLE:-stripenn}"
LABEL org.opencontainers.image.version="${CONTAINER_VERSION:-latest}"
