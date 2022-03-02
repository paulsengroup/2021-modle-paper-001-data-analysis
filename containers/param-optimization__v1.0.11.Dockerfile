# Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

FROM ghcr.io/robomics/modle:sha-7080e51 AS base

ARG CONTAINER_VERSION
ARG BIOFRAME_VER='0.3.*'
ARG PANDAS_VER='1.4.*'
ARG PYBIGWIG_VER='0.3.*'
ARG SKOPT_VER='0.9.*'
ARG PIP_NO_CACHE_DIR=0

ENV SHELL=/usr/bin/bash

LABEL maintainer='Roberto Rossini <roberros@uio.no>'
LABEL version=${CONTAINER_VERSION}
WORKDIR /data

# Update system repo and install required tools
RUN apt-get update \
    && apt-get install -y --no-install-recommends  \
                       bash        \
                       gcc         \
                       python3     \
                       python3-dev \
                       python3-pip \
    && pip install bioframe==${BIOFRAME_VER}         \
                   pandas==${PANDAS_VER}             \
                   pyBigWig==${PYBIGWIG_VER}         \
                   scikit-optimize==${SKOPT_VER}     \
    && apt-get remove -y gcc python3-dev python3-pip \
    && rm -rf /var/lib/apt/lists/*

RUN modle --help
RUN modle_tools --help
