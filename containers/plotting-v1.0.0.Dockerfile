# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

FROM fedora:34 AS base

ARG CONTAINER_VERSION=1.0.0

ARG COOLER_VER=0.8.11
ARG PILLOW_VER=9.0.0
ARG MATPLOTLIB_VER=3.5.1

ARG PIP_NO_CACHE_DIR=0

# Update system repo and install required tools
RUN dnf update -y \
    && dnf install -y --setopt=install_weak_deps=False --best \
                      gcc                \
                      perl-Digest-SHA    \
                      poppler-utils      \
                      python3            \
                      python3-Cython     \
                      python3-numpy      \
                      python3-pip        \
                      python3-devel      \
                      python3-setuptools \
                      python3-wheel      \
                      zlib-devel \
    && pip3 install cooler>=${COOLER_VER}         \
                    Pillow>=${PILLOW_VER}         \
                    matplotlib>=${MATPLOTLIB_VER} \
                    numpy                         \
    && dnf remove -y gcc python3-devel python3-Cython python3-pip \
                     python3-wheel zlib-devel                     \
    && dnf clean all

COPY scripts /tmp/scripts

RUN install -d /usr/local/bin \
    && install -D -m 755 /tmp/scripts/convert_image_to_matrix_sanborn.py /usr/local/bin/ \
    && install -D -m 755 /tmp/scripts/make_modle_vs_md_plots.py /usr/local/bin/ \
    && rm -r /tmp/scripts

ENV SHELL=/usr/bin/bash
ENV PATH='/usr/local/bin:/usr/bin'

LABEL maintainer='Roberto Rossini <roberros@uio.no>'
LABEL version=${CONTAINER_VERSION}
WORKDIR /data

RUN convert_image_to_matrix_sanborn.py --help
