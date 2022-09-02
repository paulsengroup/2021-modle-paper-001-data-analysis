#!/usr/bin/env bash

# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

set -e
set -u
set -o pipefail
set -x

### IMPORTANT: make sure caper and cromwell are properly configured before attempting to run this script

function symlink_folders() {
    ln -sf ../configs/ configs
    ln -sf ../containers/ containers
    ln -sf ../data/ data
}

function fetch_workflow() {
    git_tag='v2.2.0'
    workflow_dir="chip-seq-pipeline2-$git_tag"

    if [ ! -d "$workflow_dir" ]; then
        git clone https://github.com/ENCODE-DCC/chip-seq-pipeline2.git \
            -b "$git_tag" \
            --depth=1 \
            "$workflow_dir"
    fi

    readlink -f "$workflow_dir/chip.wdl"
}

function setup_workflow() {
    json="$1"

    # caper really doesn't play well with relative paths, so we have to update the .json config with
    # the absolute paths on the machine where the script is being run
    git_root="$(git rev-parse --show-toplevel)"
    datadir="$(readlink -f "$git_root/data/")"
    config="$(basename "$json")"
    name="${config%.json}"

    mkdir -p "$name"
    config="$name/$config"
    sed "s|\*\*\*datadir\*\*\*|$datadir|g" "$json" | tee "$config" >/dev/null

    config="$(readlink -f "$config")"

    echo "$config"
}

export symlink_folders
export fetch_workflow
export setup_workflow

step='chipseq2_encode'
echo 1>&2 "Running step $step..."

wd=".cromwell-$step-wd"
mkdir -p "$wd"

workflow_file="$(cd "$wd" && fetch_workflow)"

(cd "$wd" && symlink_folders)

configs=()
for json in "$wd/configs/"*cromwell.json; do
    json="$(readlink -f "$json")"
    configs+=($(cd "$wd" && setup_workflow "$json"))
done

for config in "${configs[@]}"; do
    wd1="$(dirname "$config")"
    (cd "$wd1" && caper run "$workflow_file" -i "$config" --singularity |& tee "${config%.json}.log" >/dev/null) &
done

wait

for config in "${configs[@]}"; do
    wd="$(dirname "$config")"
    hash="$(grep -o 'workflow [[:alnum:]-]\+ submitted' "$wd/cromwell.out"* |
        sed -E 's/ *(workflow|submitted) *//g' | tail -n 1)"
    metadata="$wd/chip/$hash/metadata.json"
    croo --method=copy --out-dir "$wd" "$metadata"

    bwig="$(grep -F 'fold-enrichment' "$wd/croo.filetable.$hash.tsv" | rev | cut -f 2 | rev)"
    narrowpeak="$(grep 'Blacklist-filtered.*narrowpeak.*merged.*regionPeak.gz' "$wd/croo.filetable.$hash.tsv" | grep -v "Pseudoreplicate" | rev | cut -f 2 | rev)"

    mkdir -p data/hoxd_encode_chip
    cp "$bwig" "data/output/hoxd_encode_chip/$(basename "$config" .json)_fold_change.bigwig"
    cp "$narrowpeak" "data/output/hoxd_encode_chip/$(basename "$config" .json)_narrow_peaks.bed.gz"
done

ls -lah "data/hoxd_encode_chip/"
echo 1>&2 "Workflow \"$step\" successfully completed!"
