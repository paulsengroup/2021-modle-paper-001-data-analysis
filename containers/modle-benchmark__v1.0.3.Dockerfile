# Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

FROM ghcr.io/paulsengroup/modle:1.0.0-rc.2 AS base

ARG CONTAINER_VERSION
ARG CONTAINER_TITLE

RUN if [ -z "$CONTAINER_VERSION" ]; then echo "Missing CONTAINER_VERSION --build-arg" && exit 1; fi

WORKDIR /data

RUN apt-get update \
    && apt-get install -y --no-install-recommends time \
    && rm -rf /var/lib/apt/lists/*

RUN command time -v modle --help
RUN command time -v modle_tools --help

LABEL org.opencontainers.image.authors='Roberto Rossini <roberros@uio.no>'
LABEL org.opencontainers.image.url='https://github.com/paulsengroup/2021-modle-paper-001-data-analysis'
LABEL org.opencontainers.image.documentation='https://github.com/2021-modle-paper-001-data-analysis'
LABEL org.opencontainers.image.source='https://github.com/paulsengroup/2021-modle-paper-001-data-analysis'
LABEL org.opencontainers.image.licenses='MIT'
LABEL org.opencontainers.image.title="${CONTAINER_TITLE:-modle-benchmark}"
LABEL org.opencontainers.image.version="${CONTAINER_VERSION:-latest}"
