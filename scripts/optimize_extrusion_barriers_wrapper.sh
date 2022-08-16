#!/usr/bin/env bash

set -e
set -o pipefail
set -u
set -x

if [ $# -ne 15 ]; then
    echo 2>&1 "Usage: $0 script_dir nthreads chrom_sizes chrom_ranges extr_barriers bin_size target_contact_density ref_cooler polish_cutoff sigma_ref sigma_mult_reg sigma_tgt sigma_mult_tgt"
    exit 1
fi

stage=0
script_dir="$1"
nthreads="$2"
chrom_sizes="$3"
chrom_ranges="$4"
bin_size="$5"
target_contact_density="$6"
extr_barriers="$7"
ref_matrix="$8"
polish_cutoff="$9"
sigma_ref="${10}"
sigma_mult_ref="${11}"
discr_thresh_ref="${12}"
sigma_tgt="${13}"
sigma_mult_tgt="${14}"
discr_thresh_tgt="${15}"

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
    --pop-size 256 \
    --num-generations 15 \
    --lambda 512 \
    --cxpb 0.5 \
    --mutpb-individual 0.5 \
    --mutpb-locus 0.1 \
    --tournament-size 3

# Keep exploring with a more conservative mutation rate
path_to_initial_pop="$(printf 'stage_%03d/stage_%03d_population.pickle' $((stage - 1)) $((stage - 1)))"
run_optimization --extrusion-barriers "$extr_barriers" \
    --initial-population "$path_to_initial_pop" \
    --pop-size 256 \
    --num-generations 150 \
    --lambda 512 \
    --cxpb 0.3 \
    --mutpb-individual 0.7 \
    --mutpb-locus 0.03 \
    --tournament-size 5

# Hopefully at this point we have a pretty good set of solutions.
# Cut back even more on the mutation rate
path_to_initial_pop="$(printf 'stage_%03d/stage_%03d_population.pickle' $((stage - 1)) $((stage - 1)))"
run_optimization --extrusion-barriers "$extr_barriers" \
    --initial-population "$path_to_initial_pop" \
    --pop-size 256 \
    --num-generations 100 \
    --lambda 256 \
    --cxpb 0.3 \
    --mutpb-individual 0.7 \
    --mutpb-locus 0.01 \
    --tournament-size 5

# At this point we (usually) have a large number of barriers with low occupancy < 0.6.
# These barriers tend to have no/minor effects in the contact matrices produced by MoDLE.
polished_annotation="$(printf 'stage_%03d/stage_%03d_barrier_annotation_polished.bed' $((stage - 1)) $((stage - 1)))"
barrier_annotations=("$(printf 'stage_%03d/stage_%03d' $((stage - 1)) $((stage - 1)))"*.bed.gz)
polish_annotation "$polish_cutoff" "${barrier_annotations[@]}" >"$polished_annotation"

run_optimization --extrusion-barriers "$polished_annotation" \
    --pop-size 256 \
    --num-generations 200 \
    --lambda 256 \
    --cxpb 0.3 \
    --mutpb-individual 0.7 \
    --mutpb-locus 0.03 \
    --tournament-size 5 \
    --hof-size 256 \
    --ncells 8

path_to_initial_pop="$(printf 'stage_%03d/stage_%03d_population.pickle' $((stage - 1)) $((stage - 1)))"
run_optimization --extrusion-barriers "$polished_annotation" \
    --initial-population "$path_to_initial_pop" \
    --pop-size 256 \
    --num-generations 20 \
    --lambda 256 \
    --cxpb 0.7 \
    --mutpb-individual 0.3 \
    --mutpb-locus 0.01 \
    --tournament-size 5 \
    --hof-size 256 \
    --ncells 16

cp "$(printf 'stage_%03d/stage_%03d' $((stage - 1)) $((stage - 1)))"*hof_000.bed.gz final_annotation.bed.gz
