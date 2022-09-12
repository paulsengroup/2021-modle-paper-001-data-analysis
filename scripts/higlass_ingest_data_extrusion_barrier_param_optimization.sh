#!/usr/bin/env bash

# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

set -e
set -u
set -o pipefail

git_root="$(git rev-parse --show-toplevel)"
name="$(basename "$git_root")"

data_dir="$git_root/data/"
higlass_dir="$data_dir/higlass/"

# rm -rf "$higlass_dir"

"$(dirname "$0")/higlass_start_instance.sh"

function uuidgen {
    "$(dirname "$0")/higlass_uuid_generator.py" "$1"
}

export uuidgen

# Ingest extrusion barrier param optimization data
for bw in "$data_dir/output/extrusion_barrier_param_optimization/optimized_barriers/"*.bw; do
    bw_name="$(basename "$bw" .bw)"
    sudo -E higlass-manage ingest \
        --project-name "extrusion_barrier_param_optimization" \
        --name "$bw_name" \
        --assembly hg38 \
        --uid "$(uuidgen "$bw_name")" \
        --hg-name "$name" \
        --no-upload \
        "output/extrusion_barrier_param_optimization/optimized_barriers/$bw_name.bw" &
done
wait

for mcool in "$data_dir/output/extrusion_barrier_param_optimization/mcools/"*.mcool; do
    mcool_name="$(basename "$mcool" .mcool)"
    sudo -E higlass-manage ingest \
        --project-name "extrusion_barrier_param_optimization" \
        --name "$mcool_name" \
        --uid "$(uuidgen "$mcool_name")" \
        --hg-name "$name" \
        --no-upload \
        "output/extrusion_barrier_param_optimization/mcools/$mcool_name.mcool" &
done
wait

bw_name="GRCh38_H1_optimized_barriers_microc_003_extrusion_barriers_lef_1d_occupancy"
sudo -E higlass-manage ingest \
    --project-name "extrusion_barrier_param_optimization" \
    --name "$bw_name" \
    --assembly hg38 \
    --uid "$(uuidgen "$bw_name")" \
    --hg-name "$name" \
    --no-upload \
    "output/extrusion_barrier_param_optimization/simulations/$bw_name.bw"
