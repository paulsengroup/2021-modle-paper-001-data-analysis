# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

FROM ubuntu:20.04 as downloader

ARG GET_QC_URL='https://raw.githubusercontent.com/dovetail-genomics/Micro-C/85d5d4aaf830b076658205f11a9413b24434789a/get_qc.py'
ARG GET_QC_SHA256='37b34225118c8c46a5abe05dfa79a304a637ffe1a393d9dc200621b0090c7da5'

RUN apt-get update \
&& apt-get install -y --no-install-recommends ca-certificates curl perl \
&& curl -L "$GET_QC_URL" > get_qc.py \
&& echo "$GET_QC_SHA256  get_qc.py" > sum.sha256 \
&& shasum -c sum.sha256 \
&& install -Dm755 get_qc.py /usr/local/bin/get_qc.py

FROM mambaorg/micromamba:0.22.0 AS base

ARG BEDTOOLS_VER='2.30.*'
ARG BWA_MEM2_VER='2.2.*'
ARG COOLTOOLS_VER='0.8.*'
ARG DEEPTOOLS_VER='3.4.*'
ARG MATPLOTLIB_VER='3.5.*'
ARG PAIRTOOLS_VER='0.3.*'
ARG PAIRIX_VER='0.3.*'
ARG PANDAS_VER='1.4.*'
ARG PRESEQ_VER='3.1.*'
ARG PYSAM_VER='0.17.*'
ARG SAMTOOLS_VER='1.14.*'
ARG TABULATE_VER='0.8.*'

ARG MAMBA_DOCKERFILE_ACTIVATE=1

ENV SHELL=/bin/bash

COPY --from=downloader '/usr/local/bin/get_qc.py' '/usr/local/bin/get_qc.py'

RUN micromamba install -y \
    -c conda-forge \
    -c bioconda \
    bedtools=$BEDTOOLS_VER \
    bwa-mem2=$BWA_MEM2_VER \
    cooler=$COOLER_VER \
    deeptools=$DEEPTOOLS_VER \
    matplotlib=$MATPLOTLIB_VER \
    pairtools=$PAIRTOOLS_VER \
    pairix=$PAIRIX_VER \
    pandas=$PANDAS_VER \
    preseq=$PRESEQ_VER \
    pysam=$PYSAM_VER \
    samtools=$SAMTOOLS_VER \
    tabulate=$TABULATE_VER \
&& micromamba clean --all -y

ENV PATH="/opt/conda/bin:$PATH"

ARG CONTAINER_VERSION
ARG CONTAINER_TITLE
WORKDIR /data

RUN get_qc.py --help

LABEL org.opencontainers.image.authors='Roberto Rossini <roberros@uio.no>'
LABEL org.opencontainers.image.url='https://github.com/paulsengroup/2021-modle-paper-001-data-analysis'
LABEL org.opencontainers.image.documentation='https://github.com/2021-modle-paper-001-data-analysis'
LABEL org.opencontainers.image.source='https://github.com/paulsengroup/2021-modle-paper-001-data-analysis'
LABEL org.opencontainers.image.licenses='MIT'
LABEL org.opencontainers.image.title="${CONTAINER_TITLE:-microc-processing}"
LABEL org.opencontainers.image.version="${CONTAINER_VERSION:-latest}"
