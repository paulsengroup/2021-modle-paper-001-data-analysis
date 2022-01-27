# Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

FROM python:3.10-bullseye AS base

ARG CONTAINER_VERSION=1.0.3
LABEL maintainer='Roberto Rossini <roberros@uio.no>'
LABEL version=${CONTAINER_VERSION}

ARG BIOFRAME_VER='0.3.*'
ARG COOLER_VER='0.8.11'
ARG MATPLOTLIB_VER='3.5.*'
ARG NUMPY_VER='1.22.*'
ARG PANDAS_VER='1.3.*'
ARG PILLOW_VER='9.0.*'
ARG PYBIGWIG_VER='0.3.18'

RUN pip install --no-cache-dir        \
        bioframe=="$BIOFRAME_VER"     \
        cooler=="$COOLER_VER"         \
        matplotlib=="$MATPLOTLIB_VER" \
        natsort                       \
        numpy=="$NUMPY_VER"           \
        pandas=="$PANDAS_VER"         \
        Pillow=="$PILLOW_VER"          \
        pyBigWig=="$BIOFRAME_VER"
