#!/usr/bin/env python3

# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import argparse

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def make_cli():
    cli = argparse.ArgumentParser()
    cli.add_argument("hic_cis_expected", type=str)
    cli.add_argument("microc_cis_expected", type=str)
    cli.add_argument("low_proc_cis_expected", type=str)
    cli.add_argument("normal_proc_cis_expected", type=str)
    cli.add_argument("high_proc_cis_expected", type=str)
    cli.add_argument("--diagonal-width",
                     type=int,
                     default=int(5.0e6))
    cli.add_argument("--output-prefix",
                     required=True,
                     type=str)
    cli.add_argument("--bin-size",
                     required=True,
                     type=int)
    cli.add_argument("--chrom",
                     type=str,
                     default="all")
    return cli


def import_data(hic, microc, low, normal, high):
    data = None

    for k, v in zip(["hic", "microc", "low_proc", "normal_proc", "high_proc"],
                    [hic, microc, low, normal, high]):
        df = pd.read_table(v)
        if "balanced.sum" in df.columns:
            df["count.sum"] = df["balanced.sum"]
            df["count.avg"] = df["balanced.avg"]

            df.drop(columns=["balanced.sum", "balanced.avg"], inplace=True)

        if data is None:
            data = df
        else:
            data = pd.merge(data, df, on=["region1", "region2", "dist"], suffixes=["", ""])

        data.rename(columns={"count.sum": f"{k}_sum",
                             "count.avg": f"{k}_avg",
                             "n_valid": f"{k}_n_valid"},
                    inplace=True)

    return data.dropna()


def plot(df, bin_size, chrom):
    fig, ax = plt.subplots(1, 1, figsize=(6.4*2, 6.4))
    x = df["dist"] * bin_size
    for pretty_cond, cond in zip(["Hi-C", "Micro-C", "MoDLE - low processivity", "MoDLE - normal processivity", "MoDLE - high processivity"],
                                 ["hic", "microc", "low_proc", "normal_proc", "high_proc"]):

        y = df[f"{cond}_sum"]
        min_ = np.min(y[y != 0])
        max_ = np.max(y)
        y = (y - min_) / (max_ - min_)
        ax.plot(x, y, label=pretty_cond)
        ax.set(title=f"Interaction decay curve ({chrom})",
               xlabel="Distance from diagonal (bp)",
               ylabel="Scaled sum (MinMax)",
               yscale="log")

    ax.legend()
    return fig, ax


if __name__ == "__main__":
    args = vars(make_cli().parse_args())

    files = [args[k] for k in ["hic_cis_expected",
                               "microc_cis_expected",
                               "low_proc_cis_expected",
                               "normal_proc_cis_expected",
                               "high_proc_cis_expected"]]
    data = import_data(*files)

    if args["chrom"] != "all":
        data = data[data["region1"] == args["chrom"]]

    data = data[data["dist"] < (args["diagonal_width"] / args["bin_size"])]

    output_prefix = args["output_prefix"]
    for chrom in data["region1"].unique():
        fig, _ = plot(data[data["region1"] == chrom], args["bin_size"], chrom)
        fig.savefig(f"{output_prefix}_{chrom}.png", dpi=600)
        fig.savefig(f"{output_prefix}_{chrom}.svg")
