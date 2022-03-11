# Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

FROM glotzerlab/software:2022.01.12-nompi AS base

ARG CONTAINER_VERSION
ARG CONTAINER_TITLE

USER root
RUN apt-get update \
    && apt-get install -y --no-install-recommends time \
    && rm -rf /var/lib/apt/lists/* \
    && pip3 install --no-cache-dir git+git://github.com/robomics/PolymerCpp.git@01486f15aafdb2ee40efec389816ccc2d41c8ad9

WORKDIR /data

RUN command time -v python -c "import hoomd; import PolymerCpp"

LABEL org.opencontainers.image.authors='Roberto Rossini <roberros@uio.no>'
LABEL org.opencontainers.image.url='https://github.com/paulsengroup/2021-modle-paper-001-data-analysis'
LABEL org.opencontainers.image.documentation='https://github.com/2021-modle-paper-001-data-analysis'
LABEL org.opencontainers.image.source='https://github.com/paulsengroup/2021-modle-paper-001-data-analysis'
LABEL org.opencontainers.image.licenses='MIT'
LABEL org.opencontainers.image.title="${CONTAINER_TITLE:-hoomd-blue-benchmark}"
LABEL org.opencontainers.image.version="${CONTAINER_VERSION:-latest}"
