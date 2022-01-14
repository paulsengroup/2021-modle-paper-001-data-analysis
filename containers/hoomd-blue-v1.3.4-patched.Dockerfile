# Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

FROM ubuntu:20.04 AS patch-hoomd

ARG HOOMD_VERSION=1.3.4

ARG HOOMD_SOURCE_URL="https://github.com/glotzerlab/hoomd-blue/archive/refs/tags/v${HOOMD_VERSION}.tar.gz"
ARG HOOMD_SOURCE_SHA256='e4b13310ad9c60813f0fc3d3cb887b8a60092161c94f9b3fd7471b88c4ce0be8'

# Download and patch hoomd
RUN apt-get update \
    && apt-get install -y --no-install-recommends  \
                       ca-certificates \
                       curl \
                       findutils \
                       perl \
    && cd /tmp \
    && curl -LO "${HOOMD_SOURCE_URL}" \
    && echo "${HOOMD_SOURCE_SHA256}  v${HOOMD_VERSION}.tar.gz" > hoomd-blue.sha256 \
    && shasum -a256 -c hoomd-blue.sha256 \
    && tar -xf "v${HOOMD_VERSION}.tar.gz" \
    && cd "hoomd-blue-${HOOMD_VERSION}" \
    && sed -i '320s/make_tuple/boost::python::make_tuple/' 'libhoomd/python/hoomd_module.cc' \
    && find . -type f -name "*.py" -exec sed -i 's/time.clock()/time.perf_counter()/g' {} + \
    && find . -type f -name "*.py" -exec sed -i 's/\.push_back(/.append(/g' {} +


FROM ubuntu:20.04 AS builder

ARG HOOMD_VERSION=1.3.4
ARG CONTAINER_VERSION=1.3.4-patched

ENV SHELL=/bin/bash

LABEL maintainer='Roberto Rossini <roberros@uio.no>'
LABEL version=${CONTAINER_VERSION}
WORKDIR /data

ARG DEBIAN_FRONTEND=noninteractive
ENV TZ=Etc/UTC

ENV OMPI_ALLOW_RUN_AS_ROOT=1
ENV OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
# https://github.com/open-mpi/ompi/issues/4948
ENV OMPI_MCA_btl_vader_single_copy_mechanism='none'

COPY --from=patch-hoomd "/tmp/hoomd-blue-${HOOMD_VERSION}" "/tmp/hoomd-blue-${HOOMD_VERSION}"

RUN apt-get update \
    && apt-get install -y --no-install-recommends  \
                       build-essential \
                       cmake \
                       libboost-chrono1.67-dev \
                       libboost-chrono1.67.0 \
                       libboost-filesystem1.67-dev \
                       libboost-filesystem1.67.0 \
                       libboost-iostreams1.67-dev \
                       libboost-iostreams1.67.0 \
                       libboost-numpy1.67-dev \
                       libboost-numpy1.67.0 \
                       libboost-program-options1.67-dev \
                       libboost-program-options1.67.0 \
                       libboost-python1.67-dev \
                       libboost-python1.67.0 \
                       libboost-serialization1.67-dev \
                       libboost-serialization1.67.0 \
                       libboost-signals1.67-dev \
                       libboost-signals1.67.0 \
                       libboost-system1.67-dev \
                       libboost-system1.67.0 \
                       libboost-test1.67-dev \
                       libboost-test1.67.0 \
                       libboost-thread1.67-dev \
                       libboost-thread1.67.0 \
                       libboost-timer1.67-dev \
                       libboost-timer1.67.0 \
                       libopenmpi-dev \
                       libpython3-dev \
                       openmpi-bin \
                       python3 \
                       python3-numpy \
    && cd "/tmp/hoomd-blue-${HOOMD_VERSION}" \
    && mkdir "build" \
    && cd "build" \
    && env GCC_ARCH='x86-64' cmake -DENABLE_MPI=ON .. \
    && cmake --build . -- -j $(nproc) \
    && ctest -E 'test_communication-mpi-cpu' -- -j $(nproc) \
    && cmake --install . \
    && find /usr/local/lib/hoomd -type f -name "*.py" -exec \
       sed -i 's/hoomd\.Bond(\(.*\),[[:space:]]\(.*\),[[:space:]]\(.*\))/hoomd.Bond(\1, int(\2), int(\3))/' {} + \
    && apt-get remove -y build-essential \
                         cmake \
                         libboost-chrono1.67-dev \
                         libboost-filesystem1.67-dev \
                         libboost-iostreams1.67-dev \
                         libboost-numpy1.67-dev \
                         libboost-program-options1.67-dev \
                         libboost-python1.67-dev \
                         libboost-serialization1.67-dev \
                         libboost-signals1.67-dev \
                         libboost-system1.67-dev \
                         libboost-test1.67-dev \
                         libboost-thread1.67-dev \
                         libboost-timer1.67-dev \
                         libopenmpi-dev \
    && apt-get autoremove -y \
    && cd / && rm -rf /var/lib/apt/lists/* /tmp/hoomd*

RUN which hoomd

ENTRYPOINT ["/usr/local/bin/hoomd"]
