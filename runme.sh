#!/usr/bin/env bash

set -e
set -u
set -o pipefail
set -x

nextflow run -c configs/fetch_data.config \
    -c configs/base.config \
    workflows/fetch_data.nf \
    -resume

nextflow run -c configs/preprocess_data.config \
    -c configs/base.config \
    workflows/preprocess_data.nf \
    -resume

nextflow run -c configs/optimize_modle_sim_params.config \
    -c configs/base.config \
    workflows/optimize_modle_sim_params.nf \
    -resume

nextflow run -c configs/genome_wide_analyses.config \
    -c configs/base.config \
    workflows/genome_wide_analyses.nf \
    -resume
