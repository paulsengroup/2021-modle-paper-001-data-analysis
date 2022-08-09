# Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

FROM ghcr.io/paulsengroup/modle:sha-ed185da AS base

ARG CONTAINER_VERSION
ARG CONTAINER_TITLE

ARG BIOFRAME_VER='0.3.*'
ARG COOLER_VER='0.8.*'
ARG CLOUDPICKLE_VER='2.0.*'
ARG DEAP_VER='1.3.*'
ARG PANDAS_VER='1.4.*'
ARG PYBIGWIG_VER='0.3.*'
ARG SKIMAGE_VER='0.19.*'
ARG SKOPT_VER='0.9.*'
ARG SCIPY_VER='1.8.*'
ARG PIP_NO_CACHE_DIR=0

ENV SHELL=/usr/bin/bash

RUN if [ -z "$CONTAINER_VERSION" ]; then echo "Missing CONTAINER_VERSION --build-arg" && exit 1; fi

RUN apt-get update \
    && apt-get install -y                    \
                       bash                  \
                       gcc                   \
                       python3               \
                       python3-dev           \
                       python3-pkg-resources \
                       python3-pip           \
    && pip install bioframe==${BIOFRAME_VER}         \
                   cooler==${COOLER_VER}             \
                   cloudpickle==${CLOUDPICKLE_VER}   \
                   deap==${DEAP_VER}                 \
                   networkx                          \
                   pandas==${PANDAS_VER}             \
                   pyBigWig==${PYBIGWIG_VER}         \
                   scikit-image==${SKIMAGE_VER}      \
                   scikit-optimize==${SKOPT_VER}     \
                   scipy==${SCIPY_VER}               \
    && apt-get remove -y gcc python3-pip \
    && apt-get autoremove -y \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /data
ENTRYPOINT ["/usr/bin/python3"]

RUN modle --help
RUN modle_tools --help

LABEL org.opencontainers.image.authors='Roberto Rossini <roberros@uio.no>'
LABEL org.opencontainers.image.url='https://github.com/paulsengroup/2021-modle-paper-001-data-analysis'
LABEL org.opencontainers.image.documentation='https://github.com/2021-modle-paper-001-data-analysis'
LABEL org.opencontainers.image.source='https://github.com/paulsengroup/2021-modle-paper-001-data-analysis'
LABEL org.opencontainers.image.licenses='MIT'
LABEL org.opencontainers.image.title="${CONTAINER_TITLE:-param-optimization}"
LABEL org.opencontainers.image.version="${CONTAINER_VERSION:-latest}"
