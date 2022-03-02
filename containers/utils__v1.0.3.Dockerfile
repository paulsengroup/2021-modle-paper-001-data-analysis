# Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

FROM fedora:35 AS base


ARG CONTAINER_VERSION

RUN dnf update -y \
    && dnf install -y --setopt=install_weak_deps=False --best \
                      findutils \
                      gawk \
                      perl-Digest-SHA \
                      pigz \
                      poppler-utils \
                      sed \
                      unzip \
                      zstd \
    && dnf clean all

ENV SHELL=/usr/bin/bash

LABEL maintainer='Roberto Rossini <roberros@uio.no>'
LABEL version=${CONTAINER_VERSION}
WORKDIR /data
