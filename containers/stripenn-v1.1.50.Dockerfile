# Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

FROM fedora:35 AS base

ARG CONTAINER_VERSION=1.1.50
ARG STRIPENN_VER=${CONTAINER_VERSION}
ARG PSUTIL_MIN_VER='5.*'
ARG PIP_NO_CACHE_DIR=0

# Update system repo and install required tools
RUN dnf update -y \
    && dnf install -y --setopt=install_weak_deps=False --best               \
                      gcc python3-devel python3-pip python3-setuptools      \
                      libglvnd-devel zlib-devel                             \
    && pip3 install "stripenn==${STRIPENN_VER}" "psutil==${PSUTIL_MIN_VER}" \
    && dnf remove -y gcc python3-devel python3-pip python3-setuptools       \
                     zlib-devel                                             \
    && dnf clean all
     
ENV SHELL=/usr/bin/bash

LABEL maintainer='Roberto Rossini <roberros@uio.no>'
LABEL version=${CONTAINER_VERSION}
WORKDIR /data
ENTRYPOINT ["/usr/local/bin/stripenn"]

RUN stripenn --help