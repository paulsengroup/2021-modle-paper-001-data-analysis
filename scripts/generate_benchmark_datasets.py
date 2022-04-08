#!/usr/bin/env python3

# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import argparse
import os
import re

import numpy as np
import pandas as pd


def make_cli():
    cli = argparse.ArgumentParser()
    cli.add_argument("--hoomd-input-parameters",
                     type=str,
                     help="Path to parameters.txt.")
    cli.add_argument("--hoomd-input-probabilities",
                     required=True,
                     type=str,
                     help="Path to in.probs.")
    cli.add_argument("--output-prefix",
                     required=True,
                     type=str,
                     help="Prefix name to use for output.")
    cli.add_argument("--min-chrom-size-mbp",
                     required=True,
                     type=float,
                     help="Size in Mbp of the smallest region to be simulated.")
    cli.add_argument("--monomer-size-bp",
                     required=True,
                     type=int,
                     help="Size of one monomer in bp.")
    cli.add_argument("--max-chrom-size-mbp",
                     required=True,
                     type=float,
                     help="Size in Mbp of the largest region to be simulated.")
    cli.add_argument("--step-size-mbp",
                     type=float,
                     default=3.0,
                     help="Step size in Mbp.")
    cli.add_argument("--mode",
                     choices=["modle", "hoomd"],
                     required=True,
                     help="Kind of dataset to generate.")
    cli.add_argument("--chrom-name",
                     default="chrN",
                     help="Chromosome name to use when generating datasets for MoDLE.")
    return cli


def convert_probs_to_bed(df, chrom_name, monomer_size):
    rev_barriers = df[df["rev_prob"] != 1.0]
    fwd_barriers = df[df["fwd_prob"] != 1.0]

    rev_barriers.rename(columns={"rev_prob": "score", "fwd_prob": "foo"}, inplace=True)
    fwd_barriers.rename(columns={"fwd_prob": "score", "rev_prob": "foo"}, inplace=True)

    rev_barriers["strand"] = "-"
    fwd_barriers["strand"] = "+"

    rev_barriers["start"] = rev_barriers.index * monomer_size
    fwd_barriers["start"] = fwd_barriers.index * monomer_size

    barriers = pd.concat([rev_barriers, fwd_barriers])

    barriers["chrom"] = chrom_name
    barriers["end"] = barriers["start"] + monomer_size
    barriers["name"] = "."

    barriers = barriers[["chrom", "start", "end", "name", "score", "strand"]]

    return barriers.sort_values(["start", "strand"])


def generate_modle_barriers(barriers, chrom_name, bin_size, num_bins, chrom_size):
    iterations = int(np.ceil(chrom_size / (bin_size * num_bins)))
    offsets = np.repeat(np.arange(iterations), barriers.shape[0]) * bin_size * num_bins

    starts = np.tile(barriers["start"].to_numpy(), iterations) + offsets
    scores = np.tile(barriers["score"].to_numpy(), iterations)
    strands = np.tile(barriers["strand"].to_numpy(), iterations)

    return pd.DataFrame({"chrom": np.repeat(chrom_name, len(starts)),
                         "start": starts,
                         "end": starts + bin_size,
                         "name": np.repeat(".", len(starts)),
                         "score": scores,
                         "strand": strands})


def generate_hoomd_params(probs, input_param_file, output_prefix, min_size, max_size, step_size, monomer_size):
    base_num_monomers = probs.shape[0]
    fwd_probs = np.tile(probs["fwd_prob"].to_numpy(), int(np.ceil(max_size / (probs.shape[0] * monomer_size))))
    rev_probs = np.tile(probs["rev_prob"].to_numpy(), int(np.ceil(max_size / (probs.shape[0] * monomer_size))))

    probs = pd.DataFrame({"fwd_prob": fwd_probs, "rev_prob": rev_probs})

    with open(input_param_file, "r") as f:
        params = f.read()
    base_num_bonds = int(re.search(r"nbonds\s*=\s*(.*)", params, re.MULTILINE)[1])
    params = re.sub(r"N\s*=.*", "", params, re.MULTILINE)
    params = re.sub(r"nbonds\s*=.*", "", params, re.MULTILINE)

    padding_len = len(str(max_size))
    for chrom_size in range(min_size, max_size + 1, step_size):
        num_monomers = round(chrom_size / monomer_size)
        out_suffix = f"{chrom_size:0{padding_len}d}_bp"
        probs.head(num_monomers).to_csv(f"{output_prefix}_probabilities_{out_suffix}.tsv",
                                        sep="\t",
                                        header=False,
                                        index=False)

        with open(f"{output_prefix}_parameters_{out_suffix}.txt", "w") as f:
            num_bonds = int(np.round(base_num_bonds * (chrom_size / (base_num_monomers * monomer_size))))
            f.write(f"{params}\nN = {num_monomers}\nnbonds = {num_bonds}\n")


def generate_modle_datasets(hoomd_probs, output_prefix, chrom_name, min_size, max_size, step_size, bin_size):
    barriers = convert_probs_to_bed(hoomd_probs, chrom_name, bin_size)
    num_bins = hoomd_probs.shape[0]
    barriers = generate_modle_barriers(barriers, chrom_name, bin_size, num_bins, max_size)

    padding_len = len(str(max_size))
    for chrom_size in range(min_size, max_size + 1, step_size):
        out_suffix = f"{chrom_size:0{padding_len}d}_bp"

        with open(f"{output_prefix}_{out_suffix}.chrom.sizes", "w") as f:
            f.write(f"{chrom_name}\t{chrom_size}\n")

        barriers[barriers["end"] <= chrom_size].to_csv(f"{output_prefix}_extrusion_barriers_{out_suffix}.bed",
                                                       header=False,
                                                       index=False,
                                                       sep="\t")


if __name__ == "__main__":
    cli = make_cli()
    args = cli.parse_args()

    if args.mode == "hoomd" and not args.hoomd_input_parameters:
        cli.error("--hoomd-input-parameters is a mandatory argument when --mode=hoomd")

    out_dir = os.path.dirname(args.output_prefix)
    if out_dir != "":
        os.makedirs(out_dir, exist_ok=True)

    min_size = int(args.min_chrom_size_mbp * 1.0e6)
    max_size = int(args.max_chrom_size_mbp * 1.0e6)
    step_size = int(args.step_size_mbp * 1.0e6)
    monomer_size = args.monomer_size_bp

    assert min_size <= max_size
    assert step_size <= max_size

    hoomd_probs = pd.read_csv(args.hoomd_input_probabilities,
                              sep=r"\s+",
                              engine="python",
                              names=["fwd_prob", "rev_prob"])

    assert hoomd_probs.shape[1] == 2

    if args.mode == "modle":
        print("Generating datasets for MoDLE...")
        generate_modle_datasets(hoomd_probs,
                                args.output_prefix,
                                args.chrom_name,
                                min_size,
                                max_size,
                                step_size,
                                monomer_size)

    else:
        assert args.mode == "hoomd"
        print("Generating parameters for HOOMD-blue...")
        generate_hoomd_params(hoomd_probs,
                              args.hoomd_input_parameters,
                              args.output_prefix,
                              min_size,
                              max_size,
                              step_size,
                              monomer_size)
