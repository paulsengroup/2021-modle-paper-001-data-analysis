#!/usr/bin/env python3

# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import argparse
import sys

import pandas as pd


def make_cli():
    cli = argparse.ArgumentParser()
    cli.add_argument("input-bed",
                     type=str,
                     help="Path to BED file with the intervals to be rearranged.")
    cli.add_argument("rearrangement-bed",
                     type=str,
                     help="Path to BED file with the list of rearrangements to perform.")
    cli.add_argument("--shift-coords",
                     action="store_true",
                     default=False,
                     help="Shift coordinates when deleting/inserting.")
    return cli


def rearrange(df, chrom_name, start_pos, end_pos, rearr_type, shift):
    assert rearr_type in {"deletion", "duplication", "inversion"}
    assert start_pos < end_pos

    other_chroms = df[~(df["chrom"] == chrom_name)]
    upstream = df[(df["chrom"] == chrom_name) &
                  (df["start"] < start_pos)]

    mut = df[(df["chrom"] == chrom_name) &
             (df["start"] >= start_pos) &
             (df["end"] < end_pos)]

    downstream = df[(df["chrom"] == chrom_name) &
                    (df["start"] >= end_pos)]

    mut_span = end_pos - start_pos

    if rearr_type == "deletion":
        if shift:
            downstream["start"] = downstream["start"] - mut_span
            downstream["end"] = downstream["end"] - mut_span

        return pd.concat([other_chroms, upstream, downstream]).sort_values(by=["chrom", "start"])

    if rearr_type == "inversion":
        mut["start"] = end_pos - (mut["start"] - start_pos)
        mut["end"] = end_pos - (mut["end"] - start_pos)
        mut["strand"] = mut["strand"].map({"+": "-", "-": "+"})
        mut["name"] = [f"{name}_inv" for name in mut["name"]]
    else:
        dup = mut.copy()
        dup["start"] = dup["start"] + mut_span
        dup["end"] = dup["end"] + mut_span
        dup["name"] = [f"{name}_dup" for name in dup["name"]]

        mut = pd.concat([mut, dup])

        downstream["start"] = downstream["start"] + mut_span
        downstream["end"] = downstream["end"] + mut_span

    return pd.concat([other_chroms, upstream, mut, downstream]).sort_values(by=["chrom", "start"])


if __name__ == "__main__":
    args = vars(make_cli().parse_args())

    input_bed = pd.read_table(args["input-bed"],
                              names=["chrom", "start", "end", "name", "score", "strand"],
                              index_col=False,
                              comment="#")
    rearr_bed = pd.read_table(args["rearrangement-bed"],
                              names=["chrom", "start", "end", "rearrangement"],
                              index_col=False)

    for (_, (chrom, start, end, mut_type)) in rearr_bed.iterrows():
        input_bed = rearrange(input_bed, chrom, start, end, mut_type, args["shift_coords"])

    input_bed.to_csv(sys.stdout, sep="\t", header=None, index=False)
