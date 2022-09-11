#!/usr/bin/env python3

# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import argparse

import matplotlib.pyplot as plt
import pandas as pd
import pathlib


def make_cli():
    cli = argparse.ArgumentParser()
    cli.add_argument("--tsv", type=str, nargs="+", required=True)
    cli.add_argument("--name", type=str, nargs="+")
    cli.add_argument("--max-distance",
                     type=int,
                     default=int(2.0e6))
    cli.add_argument("--min-distance",
                     type=int,
                     default=int(100e3))
    cli.add_argument("--output-prefix",
                     required=True,
                     type=str)
    cli.add_argument("--chrom",
                     type=str,
                     default="all")
    return cli


def import_data(tsv, min_distance, max_distance):
    df = pd.read_table(tsv)
    return df[(df["start"] >= min_distance) & (df["end"] < max_distance)]


def plot(data, chrom):
    fig, ax = plt.subplots(1, 1, figsize=(6.4*2, 6.4))
    for name, df in data.items():
        df1 = df[df["chrom"] == chrom]
        if len(df1) != 0:
            ax.plot(df1["start"], df1["correlation"], label=name)

    ax.set(title=f"Correlation by diagonal ({chrom})",
           xlabel="Distance from diagonal (bp)",
           ylabel="Correlation")

    ax.legend()
    return fig, ax


if __name__ == "__main__":
    args = vars(make_cli().parse_args())

    assert args["min_distance"] < args["max_distance"]
    assert args["name"] is None or len(args["tsv"]) == len(args["name"])

    files = args["tsv"]
    names = args["name"]
    if names is None:
        names = [pathlib.Path(f).stem for f in files]

    data = {}
    chroms = set()
    for name, path in zip(names, files):
        df = import_data(path, args["min_distance"], args["max_distance"])
        if args["chrom"] != "all":
            df = df[df["chrom"] == args["chrom"]]

        chroms.update(df["chrom"].unique().tolist())
        data[name] = df
    assert len(data) != 0

    output_prefix = args["output_prefix"]
    for chrom in chroms:
        fig, _ = plot(data, chrom)
        fig.savefig(f"{output_prefix}_{chrom}.png", dpi=600)
        fig.savefig(f"{output_prefix}_{chrom}.svg")
