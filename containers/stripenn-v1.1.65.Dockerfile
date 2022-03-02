# Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

FROM fedora:34 AS base

ARG CONTAINER_VERSION=1.1.65
ARG STRIPENN_VER=${CONTAINER_VERSION}.2
ARG PIP_NO_CACHE_DIR=0

RUN dnf update -y \
    && dnf install -y --setopt=install_weak_deps=False --best                    \
                      gcc python3-devel python3-pip python3-setuptools           \
                      libglvnd-devel zlib-devel                                  \
    && pip install stripenn==${STRIPENN_VER} --no-deps                           \
    && pip check | cut -d $' ' -f 4 | sed 's/,//' | grep -v 'pathlib'            \
       > /tmp/fix-deps.txt                                                       \
    && pip install -r /tmp/fix-deps.txt                                          \
    && pip install psutil                                                        \
    && dnf remove -y gcc python3-devel python3-pip python3-setuptools zlib-devel \
    && dnf clean all

ENV SHELL=/usr/bin/bash

LABEL maintainer='Roberto Rossini <roberros@uio.no>'
LABEL version=${CONTAINER_VERSION}
WORKDIR /data
ENTRYPOINT ["/usr/local/bin/stripenn"]

RUN stripenn --help
