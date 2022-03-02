# Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

FROM glotzerlab/software:2022.01.12-nompi AS base

ARG CONTAINER_VERSION
LABEL hoomd_version=2.9.7
LABEL maintainer='Roberto Rossini <roberros@uio.no>'
LABEL version=${CONTAINER_VERSION}
WORKDIR /data

USER root
RUN apt-get update \
    && apt-get install -y --no-install-recommends time \
    && rm -rf /var/lib/apt/lists/* \
    && pip3 install --no-cache-dir git+git://github.com/robomics/PolymerCpp.git@01486f15aafdb2ee40efec389816ccc2d41c8ad9

RUN command time -v python -c "import hoomd; import PolymerCpp"
