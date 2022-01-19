# Copyright (C) 2021 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

FROM ubuntu:20.04 AS downloader

ARG CONTAINER_VERSION=1.0.1

ARG HYPERFINE_VER=1.12.0
ARG HYPERFINE_DEB_URL="https://github.com/sharkdp/hyperfine/releases/download/v${HYPERFINE_VER}/hyperfine_${HYPERFINE_VER}_amd64.deb"
ARG HYPERFINE_DEB_SHA256='2869d59cec4179a05c55fc9959597642307cde17e7110a8105269fb407fa762a'

RUN apt-get update \
    && apt-get install -y --no-install-recommends  \
                       ca-certificates \
                       curl \
                       perl \
    && cd / && curl -LO "$HYPERFINE_DEB_URL" \
    && echo "${HYPERFINE_DEB_SHA256}  $(basename "${HYPERFINE_DEB_URL}")" > hyperfine.sha256 \
    && shasum -a256 -c hyperfine.sha256

FROM glotzerlab/software:2022.01.12-nompi AS base

LABEL hoomd_version=2.9.7
LABEL maintainer='Roberto Rossini <roberros@uio.no>'
LABEL version=${CONTAINER_VERSION}
WORKDIR /data

COPY --from=downloader "/hyperfine_*.deb" '/tmp/'

USER root
RUN dpkg -i /tmp/hyperfine_*.deb \
    && rm /tmp/hyperfine_*.deb

USER glotzerlab-software
RUN python -c "import hoomd"
RUN hyperfine --help

ENTRYPOINT ["hyperfine"]
