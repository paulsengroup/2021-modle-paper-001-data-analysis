# Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

FROM python:3.10-bullseye AS base

LABEL maintainer='Roberto Rossini <roberros@uio.no>'
LABEL version=${CONTAINER_VERSION}

ARG CONTAINER_VERSION=1.0.0

ARG NUMPY_VER='1.22.*'
ARG PANDAS_VER='1.3.*'
ARG MATPLOTLIB_VER='3.5.*'

RUN pip install --no-cache-dir \
        numpy=="$NUMPY_VER" \
        pandas=="$PANDAS_VER" \
        matplotlib=="$MATPLOTLIB_VER"
