# Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

FROM fedora:35 AS base


ARG CONTAINER_VERSION=1.0.0

RUN dnf update -y \
    && dnf install -y --setopt=install_weak_deps=False --best        \
                      zstd pigz gawk sed                             \
    && dnf clean all

ENV SHELL=/usr/bin/bash
ENV PATH='/usr/local/bin:/usr/bin'

LABEL maintainer='Roberto Rossini <roberros@uio.no>'
LABEL version=${CONTAINER_VERSION}
WORKDIR /data
ENTRYPOINT ["/usr/bin/bash"]
