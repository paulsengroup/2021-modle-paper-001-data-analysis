# Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

FROM mambaorg/micromamba:0.22.0 AS base

ARG CONTAINER_VERSION
ARG CONTAINER_TITLE
ARG STRIPENN_VER=${CONTAINER_VERSION}

ARG MAMBA_DOCKERFILE_ACTIVATE=1
ARG PIP_NO_CACHE_DIR=0

RUN micromamba install -y \
    -c conda-forge \
    -c fastai \
    -c bioconda \
    cooler \
    joblib \
    matplotlib \
    numpy \
    opencv-python-headless \
    pandas \
    scikit-image \
    scipy \
    tqdm \
    typer \
&&  pip3 install "stripenn==${STRIPENN_VER}" --no-deps

ENV PATH="/opt/conda/bin:$PATH"

WORKDIR /data
ENTRYPOINT ["/opt/conda/bin/stripenn"]

RUN stripenn --help

LABEL org.opencontainers.image.authors='Roberto Rossini <roberros@uio.no>'
LABEL org.opencontainers.image.url='https://github.com/paulsengroup/2021-modle-paper-001-data-analysis'
LABEL org.opencontainers.image.documentation='https://github.com/2021-modle-paper-001-data-analysis'
LABEL org.opencontainers.image.source='https://github.com/paulsengroup/2021-modle-paper-001-data-analysis'
LABEL org.opencontainers.image.licenses='MIT'
LABEL org.opencontainers.image.title="${CONTAINER_TITLE:-stripenn}"
LABEL org.opencontainers.image.version="${CONTAINER_VERSION:-latest}"
