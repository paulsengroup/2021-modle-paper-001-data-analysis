#!/usr/bin/env bash

# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

set -e
set -u
set -o pipefail
set -x

./run_fetch_data.sh
./run_chip_seq2_encode.sh

./run_preprocessing.sh
./run_heatmap_comparison_pt1.sh               # &
./run_heatmap_comparison_pt2.sh               # &
./run_gw_param_optimization.sh                # &
./run_extrusion_barrier_param_optimization.sh # &
./run_benchmarks.sh                           # &

# wait
