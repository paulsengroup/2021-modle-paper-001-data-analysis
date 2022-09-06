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

# Ingest matrices for heatmap-comparison pt2
for mcool in "$data_dir/output/heatmap_comparison_pt2/prob_decay/"*.mcool; do
    mcool_name="$(basename "$mcool" .mcool)"
    sudo -E higlass-manage ingest \
        --project-name "heatmap_comparison_pt2_prob_decay" \
        --name "$mcool_name" \
        --uid "$(uuidgen "$mcool_name")" \
        --hg-name "$name" \
        --no-upload \
        "output/heatmap_comparison_pt2/prob_decay/$mcool_name.mcool" &
done
wait

for mcool in "$data_dir/output/heatmap_comparison_pt2/lef_param_exploration/mcools/"*.mcool; do
    mcool_name="$(basename "$mcool" .mcool)"
    sudo -E higlass-manage ingest \
        --project-name "heatmap_comparison_pt2_lef_param_exploration" \
        --name "$mcool_name" \
        --uid "$(uuidgen "$mcool_name")" \
        --hg-name "$name" \
        --no-upload \
        "output/heatmap_comparison_pt2/lef_param_exploration/mcools/$mcool_name.mcool" &
done
wait
