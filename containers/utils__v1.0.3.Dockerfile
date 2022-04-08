# Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

FROM fedora:35 AS base


ARG CONTAINER_VERSION
ARG CONTAINER_TITLE

RUN if [ -z "$CONTAINER_VERSION" ]; then echo "Missing CONTAINER_VERSION --build-arg" && exit 1; fi

RUN dnf update -y \
&&  dnf install -y findutils \
                   gawk \
                   perl-Digest-SHA \
                   pigz \
                   poppler-utils \
                   sed \
                   unzip \
                   zstd \
&&  dnf clean all

ENV SHELL=/usr/bin/bash

WORKDIR /data

LABEL org.opencontainers.image.authors='Roberto Rossini <roberros@uio.no>'
LABEL org.opencontainers.image.url='https://github.com/paulsengroup/2021-modle-paper-001-data-analysis'
LABEL org.opencontainers.image.documentation='https://github.com/2021-modle-paper-001-data-analysis'
LABEL org.opencontainers.image.source='https://github.com/paulsengroup/2021-modle-paper-001-data-analysis'
LABEL org.opencontainers.image.licenses='MIT'
LABEL org.opencontainers.image.title="${CONTAINER_TITLE:-utils}"
LABEL org.opencontainers.image.version="${CONTAINER_VERSION:-latest}"
