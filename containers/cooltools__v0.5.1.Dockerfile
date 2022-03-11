# Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

FROM mambaorg/micromamba:0.22.0 AS base

ARG MAMBA_DOCKERFILE_ACTIVATE=1

ARG CONTAINER_VERSION
ARG CONTAINER_TITLE
ARG COOLTOOLS_VER=${CONTAINER_VERSION}
ARG PIP_NO_CACHE_DIR=0

RUN micromamba install -y \
        -c conda-forge \
        -c bioconda \
        "cooltools=$COOLTOOLS_VER" \
&& micromamba clean --all -y

ENV PATH="/opt/conda/bin:$PATH"

ENTRYPOINT ["/opt/conda/bin/cooltools"]
WORKDIR /data


RUN cooltools --help
RUN cooler --help
RUN python3 -c "import cooler; import cooltools"

LABEL org.opencontainers.image.authors='Roberto Rossini <roberros@uio.no>'
LABEL org.opencontainers.image.url='https://github.com/paulsengroup/2021-modle-paper-001-data-analysis'
LABEL org.opencontainers.image.documentation='https://github.com/2021-modle-paper-001-data-analysis'
LABEL org.opencontainers.image.source='https://github.com/paulsengroup/2021-modle-paper-001-data-analysis'
LABEL org.opencontainers.image.licenses='MIT'
LABEL org.opencontainers.image.title="${CONTAINER_TITLE:-cooltools}"
LABEL org.opencontainers.image.version="${CONTAINER_VERSION:-latest}"

