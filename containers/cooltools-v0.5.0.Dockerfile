# Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

FROM fedora:35 AS base

ARG CONTAINER_VERSION=0.5.0
ARG COOLTOOLS_VER=${CONTAINER_VERSION}
ARG PIP_NO_CACHE_DIR=0

# Update system repo and install required tools
RUN dnf update -y \
    && dnf install -y --setopt=install_weak_deps=False --best  \
                      gcc python3 python3-Cython python3-numpy \
                      python3-pip python3-devel zlib-devel     \
    && python3 -m pip install cooltools==${COOLTOOLS_VER}      \
    && dnf remove -y gcc python3-devel python3-pip zlib-devel  \
    && dnf clean all
     
ENV SHELL=/usr/bin/bash
ENV PATH='/usr/local/bin:/usr/bin'

LABEL maintainer='Roberto Rossini <roberros@uio.no>'
LABEL version=${CONTAINER_VERSION}
WORKDIR /data
ENTRYPOINT ["/usr/local/bin/cooltools"]

RUN cooltools --help
