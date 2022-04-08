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
              --max_memory=400.G
              --max_cpus=52
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
       make_heatmaps_comparison
       optimize_modle_sim_params
)

for step in "${steps[@]}"; do
    run_workflow "$step"
done

if [[ $HOSTNAME == *.saga ]]; then
    mkdir -p .nextflow-benchmarks-wd
    (cd .nextflow-benchmarks-wd && ln -sf ../configs/ configs)
    (cd .nextflow-benchmarks-wd && ln -sf ../containers/ containers)
    (cd .nextflow-benchmarks-wd && ln -sf ../data/ data)
    (cd .nextflow-benchmarks-wd && ln -sf ../scripts/ scripts)
    (cd .nextflow-benchmarks-wd && ln -sf ../workflows/ workflows)

    (cd .nextflow-benchmarks-wd && \
        nextflow run \
            --max_memory=400.G \
            --max_cpus=52 \
            --max_time=336.h \
            --project="${SLURM_PROJECT_ID-changeme}" \
            -profile benchmark \
            -c "configs/run_benchmarks.config" \
            -c "configs/base_saga.config" \
            "workflows/run_benchmarks.nf" \
            -resume)
fi
