# Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

FROM fedora:35 AS builder


ARG CONTAINER_VERSION
ARG BEDTOOLS_VER=${CONTAINER_VERSION}
ARG BEDTOOLS_URL="https://github.com/arq5x/bedtools2/releases/download/v$BEDTOOLS_VER/bedtools-$BEDTOOLS_VER.tar.gz"

# Update system repo and install required tools
RUN dnf update -y \
    && dnf install -y --setopt=install_weak_deps=False --best        \
    diffutils gcc-c++ make python xz bzip2-devel xz-devel zlib-devel

RUN cd /tmp && curl -L "$BEDTOOLS_URL" | tar -xzf -

RUN cd "/tmp/bedtools2"  \
    && make -j $(nproc)  \
    && make test 2>&1 | grep -q 'Tools failing:  negativecontrol'

FROM fedora:35 AS base

COPY --from=builder "/tmp/bedtools2/bin" "/usr/local/bin"
COPY --from=builder "/tmp/bedtools2/LICENSE" "/usr/local/share/licenses/bedtools/LICENSE"

RUN dnf update -y \
     && dnf install -y --setopt=install_weak_deps=False --best        \
     xz bzip2 zlib && dnf clean all

RUN chown -R root:root /usr/local/share/licenses/

ENV SHELL=/usr/bin/bash
ENV PATH='/usr/local/bin:/usr/bin'

LABEL maintainer='Roberto Rossini <roberros@uio.no>'
LABEL version=2.30.0
WORKDIR /data
ENTRYPOINT ["/usr/local/bin/bedtools"]

RUN bedtools --help
