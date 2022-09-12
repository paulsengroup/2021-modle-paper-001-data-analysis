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

for bw in "$data_dir/output/hoxd_encode_chip/"*.bigwig; do
    bw_name="$(basename "$bw" .bigwig)"
    sudo -E higlass-manage ingest \
        --project-name "ChIP" \
        --name "$bw_name" \
        --assembly mm10 \
        --uid "$(uuidgen "$bw_name")" \
        --hg-name "$name" \
        --no-upload \
        "output/hoxd_encode_chip/$bw_name.bigwig"
done

bg="$data_dir/input/GRCm38_E12.5PL_H3K27ac.bedgraph.gz"
if [ ! -f "${bg%.bedgraph.gz}.bw" ]; then
    tmpdir="$(mktemp -d)"
    trap "rm -rf '$tmpdir'" EXIT

    zcat "$bg" |
        grep -E 'chr[[:digit:]XY]+[[:space:]]' |
        sort -k1,1V -k2,2n |
        tee "$tmpdir/tmp.bg" >/dev/null

    bedGraphToBigWig "$tmpdir/tmp.bg" \
        "$data_dir/output/preprocessing/chrom_sizes/GRCm38.chrom.sizes" \
        "$data_dir/input/GRCm38_E12.5PL_H3K27ac.bw"
fi

bw="$data_dir/input/GRCm38_E12.5PL_H3K27ac.bw"
bw_name="$(basename "$bw" .bw)"
sudo -E higlass-manage ingest \
    --project-name "ChIP" \
    --name "$bw_name" \
    --assembly mm10 \
    --uid "$(uuidgen "$bw_name")" \
    --hg-name "$name" \
    --no-upload \
    "input/$bw_name.bw"
