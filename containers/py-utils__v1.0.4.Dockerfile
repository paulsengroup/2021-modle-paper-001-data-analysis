# Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

FROM python:3.10-bullseye AS base

ARG CONTAINER_VERSION
ARG CONTAINER_TITLE

ARG BIOFRAME_VER='0.3.*'
ARG COOLER_VER='0.8.11'
ARG MATPLOTLIB_VER='3.5.*'
ARG NUMPY_VER='1.22.*'
ARG PANDAS_VER='1.3.*'
ARG PILLOW_VER='9.0.*'
ARG PYBIGWIG_VER='0.3.18'
ARG SCIPY_VER='1.8.*'

RUN if [ -z "$CONTAINER_VERSION" ]; then echo "Missing CONTAINER_VERSION --build-arg" && exit 1; fi

RUN pip install --no-cache-dir        \
        bioframe=="$BIOFRAME_VER"     \
        cooler=="$COOLER_VER"         \
        matplotlib=="$MATPLOTLIB_VER" \
        natsort                       \
        numpy=="$NUMPY_VER"           \
        pandas=="$PANDAS_VER"         \
        Pillow=="$PILLOW_VER"         \
        pyBigWig=="$BIOFRAME_VER"     \
        scipy=="$SCIPY_VER"

LABEL org.opencontainers.image.authors='Roberto Rossini <roberros@uio.no>'
LABEL org.opencontainers.image.url='https://github.com/paulsengroup/2021-modle-paper-001-data-analysis'
LABEL org.opencontainers.image.documentation='https://github.com/2021-modle-paper-001-data-analysis'
LABEL org.opencontainers.image.source='https://github.com/paulsengroup/2021-modle-paper-001-data-analysis'
LABEL org.opencontainers.image.licenses='MIT'
LABEL org.opencontainers.image.title="${CONTAINER_TITLE:-py-utils}"
LABEL org.opencontainers.image.version="${CONTAINER_VERSION:-latest}"
