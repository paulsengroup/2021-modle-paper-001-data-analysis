#!/usr/bin/env bash

set -e
set -u
set -o pipefail
set -x


function run_workflow() {
    name="$1"
    if [[ $HOSTNAME == *.saga ]]; then
        base_config='configs/base_saga.config'
        args=("${@:2}"
              --max_memory=3000.G
              --max_cpus=64
              --max_time=336.h
              --project="${SLURM_PROJECT_ID-changeme}")
    else
        base_config='configs/base_hovig.config'
        args=("${@:2}")
    fi

    nextflow run \
        "${args[@]}" \
        -c "configs/$name.config" \
        -c "$base_config" \
        "workflows/$name.nf" \
        -resume
}


steps=(fetch_data
       preprocess_data
       optimize_modle_sim_params
#       genome_wide_analyses
)

for step in "${steps[@]}"; do
    run_workflow "$step"
done
