#!/usr/bin/env bash

set -e
set -u
set -o pipefail

argc=$#

if [ $argc -ne 3 ]; then
    2>&1 echo "Usage: $(basename "$0") script_dir pickle_dir output_dir"
    exit 1
fi

script_dir="$1"
pickle_dir="$2"
output_dir="$3"

pickle1="$pickle_dir/modle_sim_param_optimization_tad_plus_loop_bayesian.pickle"
pickle2="$pickle_dir/modle_sim_param_optimization_loop_only_bayesian.pickle"

mkdir -p "$output_dir"

"$script_dir/plot_param_optimization_results.py"    \
    -o "$output_dir/$(basename "$pickle1" .pickle)" \
    "$pickle1"

"$script_dir/plot_param_optimization_results.py"    \
    -o "$output_dir/$(basename "$pickle2" .pickle)" \
    "$pickle2"
