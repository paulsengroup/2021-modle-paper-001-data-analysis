#!/usr/bin/env bash

set -e
set -u
set -o pipefail

git_root="$(git rev-parse --show-toplevel)"
name="$(basename "$git_root")"

data_dir="$git_root/data/"
higlass_dir="$data_dir/higlass/"

rm -rf "$higlass_dir"

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

# matrices=("$data_dir/$name/matrices/"*.mcool "$data_dir/preprocessing/mcools/"GRCh38_H1_hESC*.mcool "$wd/data/media/uploads/"*.mcool)
# for m in "${matrices[@]}"; do
#     sudo -E higlass-manage ingest \
#         --assembly hg38 \
#         --name "$(basename "$m" .mcool)" \
#         "$(basename "$m")" \
#         --hg-name "$name" \
#         --no-upload
# done
