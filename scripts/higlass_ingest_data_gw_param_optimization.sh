#!/usr/bin/env bash

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

# Ingest matrices for genome-wide param optimization
for mcool in "$data_dir/output/gw_param_optimization/mcools/"*.mcool; do
    mcool_name="$(basename "$mcool" .mcool)"
    sudo -E higlass-manage ingest \
        --project-name "gw_param_optimization" \
        --name "$mcool_name" \
        --uid "$(uuidgen "$mcool_name")" \
        --hg-name "$name" \
        --no-upload \
        "output/gw_param_optimization/mcools/$mcool_name.mcool"
done
