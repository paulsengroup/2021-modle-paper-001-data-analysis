#!/usr/bin/env bash

set -e
set -u
set -o pipefail

argc=$#

if [ $argc -ne 3 ]; then
    2>&1 echo "Usage: $(basename "$0") path/to/summarize_benchmark_run.py benchmark_report_dirbenchmark_summary_dir output_folder"
    exit 1
fi

script="$1"

indir="$2"
outdir="$3"

files=("$indir/modle_mt_weak_scaling_benchmark_report.tsv"
       "$indir/modle_ncells_benchmark_report.tsv"
       "$indir/modle_strong_scaling_benchmark_report.tsv"
       "$indir/modle_strong_scaling_no_smt_benchmark_report.tsv"
       "$indir/modle_st_weak_scaling_benchmark_report.tsv"
       "$indir/openmm_cpu_weak_scaling_benchmark_report.tsv"
       "$indir/openmm_gpu_weak_scaling_benchmark_report.tsv")

keys=("chrom_size"
      "ncells"
      "ncores"
      "ncores"
      "chrom_size"
      "chrom_size"
      "chrom_size")

mkdir -p "$outdir"

for i in "${!files[@]}"; do
    f="${files[$i]}"
    outname="$outdir/$(basename "$f" _report.tsv)_summary.tsv"
    "$script" --groupby-keys "${keys[$i]}" "$f" | tee "$outname" > /dev/null
done
