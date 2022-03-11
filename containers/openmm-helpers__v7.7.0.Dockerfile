# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

FROM mambaorg/micromamba:0.22.0 AS base

ARG CONTAINER_VERSION
ARG OPENMM_VERSION=${CONTAINER_VERSION}
ARG CUDA_VERSION=11.1
ARG PYTHON_VERSION=3.9

ARG MAMBA_DOCKERFILE_ACTIVATE=1
ARG PIP_NO_CACHE_DIR=0

ENV SHELL=/bin/bash

RUN if [ -z "$CONTAINER_VERSION" ]; then echo "Missing CONTAINER_VERSION --build-arg" && exit 1; fi

RUN micromamba install -y \
    -c conda-forge \
    -c bioconda \
    'biopython>=1.79' \
    click \
    'cooler>=0.8.11' \
    cudatoolkit=$CUDA_VERSION \
    Cython \
    gcc \
    git \
    gxx \
    'h5py>=2.5' \
    joblib \
    'numpy>=1.9' \
    openmm=$OPENMM_VERSION \
    'pandas>=0.19' \
    python=$PYTHON_VERSION \
    'scipy>=0.16' \
    seaborn \
    six \
    && micromamba clean --all -y

ARG POLYCHROM_COMMIT_HASH=0ebb3e0c66dabda3c27b855eff2bdcecaae7d2c9
RUN git clone https://github.com/open2c/polychrom.git /tmp/polychrom \
    && cd /tmp/polychrom \
    && git checkout $POLYCHROM_COMMIT_HASH \
    && find . -type f -exec sed -i 's/import simtk.openmm as openmm/import openmm/g' {} + \
    && echo "click" >> requirements.txt \
    && sed -i '/ext_modules=ext_modules,/a\    include_dirs=[numpy.get_include()],' setup.py \
    && sed -i '/"polychrom\/__polymer_math.cpp"],/a\            include_dirs=[numpy.get_include()]' setup.py \
    && pip install . \
    && cd /tmp && rm -rf polychrom*

USER root
RUN git clone https://github.com/open2c/polychrom.git /tmp/polychrom \
    && cd /tmp/polychrom \
    && git checkout $POLYCHROM_COMMIT_HASH \
    && find . -type f -exec sed -i 's/import simtk.openmm as openmm/import openmm/g' {} + \
    && install -Dm755 polychrom/cli/traj_convert.py /usr/local/bin/polychrom_traj_convert \
    && install -Dm755 polychrom/cli/dif /usr/local/bin/polychrom_dif \
    && install -Dm755 polychrom/cli/show /usr/local/bin/polychrom_show \
    && install -Dm755 polychrom/cli/xyz /usr/local/bin/polychrom_xyz \
    && cd /tmp && rm -rf polychrom*
USER mambauser

ARG OPENMM_POLYMER_LEGAGY_HASH=f35dafb47ac005b63a30211eaf6e8a1650fb8b97
RUN git clone https://github.com/mirnylab/openmm-polymer-legacy /tmp/openmm-polymer-legacy \
    && cd /tmp/openmm-polymer-legacy \
    && git checkout $OPENMM_POLYMER_LEGAGY_HASH \
    && find . -type f -exec sed -i 's/import simtk.openmm as openmm/import openmm/g' {} + \
    && pip install . \
    && cd /tmp && rm -rf openmm-polymer-legacy*

ARG MIRNYLIB_LEGACY_HASH=2d27793111b63ac7d5f82ba3369f927c6c9745da
RUN git clone https://github.com/mirnylab/mirnylib-legacy.git /tmp/mirnylib-legacy \
    && cd /tmp/mirnylib-legacy/ \
    && git checkout $MIRNYLIB_LEGACY_HASH \
    && find . -type f -exec sed -i 's/"-march=native" , //g' {} + \
    && pip install . \
    && cd /tmp && rm -rf mirnylib-legacy*

ARG LOOPLIB_HASH=bf39727a323e6aec383b8cee008c92cc3f0ff59a
RUN git clone https://github.com/golobor/looplib.git /tmp/looplib \
    && cd /tmp/looplib \
    && git checkout $LOOPLIB_HASH \
    && pip install . \
    && cd /tmp && rm -rf looplib*

ENV PATH="/opt/conda/bin:$PATH"

ENV CC=/opt/conda/bin/gcc
ENV CXX=/opt/conda/bin/g++

ARG CONTAINER_VERSION
ARG CONTAINER_TITLE
WORKDIR /data

LABEL org.opencontainers.image.authors='Roberto Rossini <roberros@uio.no>'
LABEL org.opencontainers.image.url='https://github.com/paulsengroup/2021-modle-paper-001-data-analysis'
LABEL org.opencontainers.image.documentation='https://github.com/2021-modle-paper-001-data-analysis'
LABEL org.opencontainers.image.source='https://github.com/paulsengroup/2021-modle-paper-001-data-analysis'
LABEL org.opencontainers.image.licenses='MIT'
LABEL org.opencontainers.image.title="${CONTAINER_TITLE:-openmm-helpers}"
LABEL org.opencontainers.image.version="${CONTAINER_VERSION:-latest}"
