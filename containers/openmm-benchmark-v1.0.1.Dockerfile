# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

FROM ghcr.io/robomics/2021-modle-paper-001-data-analysis/openmm-helpers:7.7.0 AS base

ARG CONTAINER_VERSION=1.0.1
LABEL openmm_version=7.7.0
LABEL maintainer='Roberto Rossini <roberros@uio.no>'
LABEL version=${CONTAINER_VERSION}
WORKDIR /data

USER root
RUN apt-get update \
    && apt-get install -y --no-install-recommends time \
    && rm -rf /var/lib/apt/lists/*
