# Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

FROM ghcr.io/robomics/modle:sha-b76b0a6 AS base

ARG CONTAINER_VERSION=1.0.5
ARG PANDAS_VER='1.3.*'
ARG PYBIGWIG_VER='0.3.*'
ARG SKOPT_VER='0.9.*'
ARG PIP_NO_CACHE_DIR=0

ENV SHELL=/usr/bin/bash
ENV PATH='/usr/local/bin:/usr/bin'

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
    && pip install pandas==${PANDAS_VER}             \
                   pyBigWig==${PYBIGWIG_VER}         \
                   scikit-optimize==${SKOPT_VER}     \
    && apt-get remove -y gcc python3-dev python3-pip \
    && rm -rf /var/lib/apt/lists/*

COPY scripts /tmp/scripts

RUN install -d /usr/local/bin \
    && install -D -m 755 /tmp/scripts/optimize_modle_sim_params.py /usr/local/bin/ \
    && rm -r /tmp/scripts

RUN optimize_modle_sim_params.py --help
RUN modle --help
RUN modle_tools --help

ENTRYPOINT ["/usr/local/bin/optimize_modle_sim_params.py"]
