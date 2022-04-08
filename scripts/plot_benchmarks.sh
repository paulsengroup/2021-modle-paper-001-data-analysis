#!/usr/bin/env bash

# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

set -e
set -u
set -o pipefail

argc=$#

if [ $argc -ne 3 ]; then
    2>&1 echo "Usage: $(basename "$0") plotting_script_dir benchmark_summary_dir output_dir"
    exit 1
fi

script_dir="$1"

indir="$2"
outdir="$3"

modle_st="$indir/modle_st_weak_scaling_benchmark_summary.tsv"
modle_mt="$indir/modle_mt_weak_scaling_benchmark_summary.tsv"
md_cpu="$indir/openmm_cpu_weak_scaling_benchmark_summary.tsv"
md_gpu="$indir/openmm_gpu_weak_scaling_benchmark_summary.tsv"

modle_scaling="$indir/modle_strong_scaling_benchmark_summary.tsv"
modle_scaling_no_smt="$indir/modle_strong_scaling_no_smt_benchmark_summary.tsv"

mkdir -p "$outdir"

"$script_dir/plot_benchmark_wallclock_md_comparison.py"   \
                        --modle-mt-report-tsv "$modle_mt" \
                        --md-cpu-report-tsv "$md_cpu"     \
                        --md-gpu-report-tsv "$md_gpu"     \
                        --extrapolate                     \
                        --output-prefix "$outdir/benchmark_artificial_chrom_wc"


"$script_dir/plot_benchmark_memory_md_comparison.py"      \
                        --modle-mt-report-tsv "$modle_mt" \
                        --md-cpu-report-tsv "$md_cpu"     \
                        --md-gpu-report-tsv "$md_gpu"     \
                        --output-prefix "$outdir/benchmark_artificial_chrom_memory"

"$script_dir/plot_benchmark_wallclock_modle.py"           \
                        --modle-st-report-tsv "$modle_st" \
                        --modle-mt-report-tsv "$modle_mt" \
                        --output-prefix "$outdir/benchmark_artificial_chrom_wc_modle_only"

"$script_dir/plot_benchmark_scaling_modle.py"               \
                        --modle-report-tsv "$modle_scaling" \
                        --title "MoDLE strong scaling"      \
                        --output-prefix "$outdir/benchmark_hg38_modle_scaling"

"$script_dir/plot_benchmark_scaling_modle.py"                      \
                        --modle-report-tsv "$modle_scaling_no_smt" \
                        --title "MoDLE strong scaling wo/ SMT"     \
                        --yaxis-left-interval='-0.075,3.2'         \
                        --yaxis-right-interval='0,2750'            \
                        --output-prefix "$outdir/benchmark_hg38_modle_scaling_no_smt"


