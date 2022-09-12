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

# Ingest chrom sizes
sudo -E higlass-manage ingest \
    --filetype chromsizes-tsv \
    --project-name "Chromosomes" \
    --datatype chromsizes \
    --assembly hg38 \
    --name GRCh38 \
    --uid "$(uuidgen "GRCh38.chrom.sizes")" \
    --no-upload \
    --hg-name "$name" \
    output/preprocessing/chrom_sizes/GRCh38.chrom.sizes

sudo -E higlass-manage ingest \
    --filetype chromsizes-tsv \
    --project-name "Chromosomes" \
    --datatype chromsizes \
    --assembly hg19 \
    --name GRCh37 \
    --uid "$(uuidgen "GRCh37.chrom.sizes")" \
    --no-upload \
    --hg-name "$name" \
    output/preprocessing/chrom_sizes/GRCh37.chrom.sizes

sudo -E higlass-manage ingest \
    --filetype chromsizes-tsv \
    --project-name "Chromosomes" \
    --datatype chromsizes \
    --assembly mm10 \
    --name GRCm38 \
    --uid "$(uuidgen "GRCm38.chrom.sizes")" \
    --no-upload \
    --hg-name "$name" \
    output/preprocessing/chrom_sizes/GRCm38.chrom.sizes

# Ingest reference matrices
for mcool in "$data_dir/output/preprocessing/mcools/"*.mcool; do
    mcool_name="$(basename "$mcool" .mcool)"
    sudo -E higlass-manage ingest \
        --project-name "Reference" \
        --name "$mcool_name" \
        --uid "$(uuidgen "$mcool_name")" \
        --hg-name "$name" \
        --no-upload \
        "output/preprocessing/mcools/$mcool_name.mcool"
done

# Ingest ChIP bigwig files
for bw in "$data_dir/input/"GRCh37*.bigwig; do
    bw_name="$(basename "$bw" .bigwig)"
    sudo -E higlass-manage ingest \
        --project-name "ChIP" \
        --name "$bw_name" \
        --assembly hg19 \
        --uid "$(uuidgen "$bw_name")" \
        --hg-name "$name" \
        --no-upload \
        "input/$bw_name.bigwig"
done

for bw in "$data_dir/input/"GRCh38*.bigwig; do
    bw_name="$(basename "$bw" .bigwig)"
    sudo -E higlass-manage ingest \
        --project-name "ChIP" \
        --name "$bw_name" \
        --assembly hg38 \
        --uid "$(uuidgen "$bw_name")" \
        --hg-name "$name" \
        --no-upload \
        "input/$bw_name.bigwig"
done

for bw in "$data_dir/input/"GRCm38*.bigwig; do
    bw_name="$(basename "$bw" .bigwig)"
    sudo -E higlass-manage ingest \
        --project-name "ChIP" \
        --name "$bw_name" \
        --assembly mm10 \
        --uid "$(uuidgen "$bw_name")" \
        --hg-name "$name" \
        --no-upload \
        "input/$bw_name.bigwig"
done

# Ingest matrices for heatmap-comparison pt1
for mcool in "$data_dir/output/heatmap_comparison_pt1/matrices/"*.mcool; do
    mcool_name="$(basename "$mcool" .mcool)"
    sudo -E higlass-manage ingest \
        --project-name "heatmap_comparison_pt1" \
        --name "$mcool_name" \
        --uid "$(uuidgen "$mcool_name")" \
        --hg-name "$name" \
        --no-upload \
        "output/heatmap_comparison_pt1/matrices/$mcool_name.mcool"
done

exit 0

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

exit 0

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
        "output/extrusion_barrier_param_optimization/optimized_barriers/$bw_name.bw"
done

# Ingest extrusion barrier param optimization data
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
        "output/extrusion_barrier_param_optimization/optimized_barriers/$bw_name.bw"
done

exit 0

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

# Ingest matrices for heatmap-comparison pt2
for mcool in "$data_dir/output/heatmap_comparison_pt2/prob_decay/"*.mcool; do
    mcool_name="$(basename "$mcool" .mcool)"
    sudo -E higlass-manage ingest \
        --project-name "heatmap_comparison_pt2_prob_decay" \
        --name "$mcool_name" \
        --uid "$(uuidgen "$mcool_name")" \
        --hg-name "$name" \
        --no-upload \
        "output/heatmap_comparison_pt2/prob_decay/$mcool_name.mcool"
done

for mcool in "$data_dir/output/heatmap_comparison_pt2/lef_param_exploration/mcools/"*.mcool; do
    mcool_name="$(basename "$mcool" .mcool)"
    sudo -E higlass-manage ingest \
        --project-name "heatmap_comparison_pt2_lef_param_exploration" \
        --name "$mcool_name" \
        --uid "$(uuidgen "$mcool_name")" \
        --hg-name "$name" \
        --no-upload \
        "output/heatmap_comparison_pt2/lef_param_exploration/mcools/$mcool_name.mcool"
done
