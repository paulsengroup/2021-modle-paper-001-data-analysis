#!/usr/bin/env bash

set -e
set -o pipefail
set -u
set -x

if [ $# -ne 14 ]; then
    echo 2>&1 "Usage: $0 script_dir nthreads chrom_sizes chrom_ranges extr_barriers bin_size target_contact_density ref_cooler sigma_ref sigma_mult_ref discr_thresh_ref sigma_tgt sigma_mult_tgt discr_thresh_tgt"
    exit 1
fi

stage=0
script_dir="$1"
nthreads="$2"
chrom_sizes="$3"
chrom_ranges="$4"
extr_barriers="$5"
bin_size="$6"
target_contact_density="$7"
ref_matrix="$8"
sigma_ref="$9"
sigma_mult_ref="${10}"
discr_thresh_ref="${11}"
sigma_tgt="${12}"
sigma_mult_tgt="${13}"
discr_thresh_tgt="${14}"

function compute_seeds {

    pyscript="$(mktemp /tmp/script.XXXXXX)"

    trap 'rm -rf "$pyscript"' EXIT

    # Hash this script and use that seed Python PRNG. Then generate 10 randints
    (
        printf '%s\n' 'import hashlib; import random'
        printf 'with open("%s", "rb") as f:\n' "$0"
        printf '\t%s\n' 'm = hashlib.sha256()'
        printf '\t%s\n' 'm.update(f.read())'
        printf '\t%s\n' 'random.seed(int(m.hexdigest(), base=16) % ((2**32) - 1))'
        printf '%s\n' 'print(" ".join([str(random.randint(0, (2**32) - 1)) for _ in range(10)]))'
    ) >>"$pyscript"

    python3 "$pyscript"
}

read -r -d '' -a seeds < <(compute_seeds && printf '\0')
echo 1>&2 "Seeds: ${seeds[*]}"

function run_optimization {
    echo 1>&2 "Running stage #$stage..."
    out_prefix="$(printf 'stage_%03d/stage_%03d' "$stage" "$stage")"
    mkdir -p "$(dirname "$out_prefix")"

    # shellcheck disable=SC2206
    args=(${@})
    "$script_dir/optimize_extrusion_barrier_params.py" \
        --chrom-sizes "$chrom_sizes" \
        --chrom-subranges "$chrom_ranges" \
        --bin-size "$bin_size" \
        --target-contact-density "$target_contact_density" \
        --reference-matrix "$ref_matrix" \
        --output-prefix "$out_prefix" \
        --nthreads "$nthreads" \
        --seed "${seeds[$stage]}" \
        --gaussian-blur-sigma-ref "$sigma_ref" \
        --gaussian-blur-sigma-multiplier-ref "$sigma_mult_ref" \
        --discretization-thresh-ref "$discr_thresh_ref" \
        --gaussian-blur-sigma-tgt "$sigma_tgt" \
        --gaussian-blur-sigma-multiplier-tgt "$sigma_mult_tgt" \
        --discretization-thresh-tgt "$discr_thresh_tgt" \
        "${args[@]}" |& tee "$out_prefix.log"

    echo 1>&2 "Done running stage #$stage..."
    stage=$((stage + 1))
}

function polish_annotation {
    threshold="$1"
    # shellcheck disable=SC2206
    input_annotations=(${@:2})

    "$script_dir/polish_optimized_barrier_annotation.py" \
        --threshold="$threshold" \
        "${input_annotations[@]}"
}

# Try to rapidly improve the initial (random) population by mutating a lot and fast
run_optimization --extrusion-barriers "$extr_barriers" \
    --pop-size 128 \
    --num-generations 25 \
    --lambda 256 \
    --cxpb 0.5 \
    --mutpb-individual 0.5 \
    --mutpb-locus 0.1 \
    --mut-sigma 0.1 \
    --tournament-size 3

path_to_initial_pop="$(printf 'stage_%03d/stage_%03d_population.pickle' $((stage - 1)) $((stage - 1)))"
extr_barriers="$(printf 'stage_%03d/stage_%03d_extrusion_barriers.bed.gz' $((stage - 1)) $((stage - 1)))"
run_optimization --extrusion-barriers "$extr_barriers" \
    --initial-population "$path_to_initial_pop" \
    --pop-size 128 \
    --num-generations 50 \
    --lambda 256 \
    --cxpb 0.3 \
    --mutpb-individual 0.7 \
    --mutpb-locus 0.03 \
    --tournament-size 5

path_to_initial_pop="$(printf 'stage_%03d/stage_%03d_population.pickle' $((stage - 1)) $((stage - 1)))"
extr_barriers="$(printf 'stage_%03d/stage_%03d_extrusion_barriers.bed.gz' $((stage - 1)) $((stage - 1)))"
run_optimization --extrusion-barriers "$extr_barriers" \
    --initial-population "$path_to_initial_pop" \
    --pop-size 256 \
    --num-generations 200 \
    --weak-barrier-purge-interval 15 \
    --lambda 256 \
    --cxpb 0.3 \
    --mutpb-individual 0.7 \
    --mutpb-locus 0.03 \
    --tournament-size 5

path_to_initial_pop="$(printf 'stage_%03d/stage_%03d_population.pickle' $((stage - 1)) $((stage - 1)))"
extr_barriers="$(printf 'stage_%03d/stage_%03d_extrusion_barriers.bed.gz' $((stage - 1)) $((stage - 1)))"
run_optimization --extrusion-barriers "$extr_barriers" \
    --initial-population "$path_to_initial_pop" \
    --pop-size 256 \
    --num-generations 50 \
    --lambda 375 \
    --cxpb 0.1 \
    --mutpb-individual 0.5 \
    --mutpb-locus 0.01 \
    --mut-sigma 0.025 \
    --tournament-size 5

cp "$(printf 'stage_%03d/stage_%03d' $((stage - 1)) $((stage - 1)))_extrusion_barriers.bed.gz" extrusion_barriers_final.bed.gz
