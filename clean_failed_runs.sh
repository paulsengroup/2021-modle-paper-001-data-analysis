#!/usr/bin/env bash

# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

set -e
set -o pipefail
set -u

nextflow clean -f $(nextflow log | awk -F '[[:space:]]+' '$5=="ERR"{print $4}' | tr '\n' ' ')
