#!/usr/bin/env bash

set -e
set -o pipefail
set -u

nextflow clean -f $(nextflow log | awk -F '[[:space:]]+' '$5=="ERR"{print $4}' | tr '\n' ' ')
