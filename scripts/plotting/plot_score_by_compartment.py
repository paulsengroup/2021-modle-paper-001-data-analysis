#!/usr/bin/env python3

# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import argparse

import matplotlib.pyplot as plt
import pandas as pd


def make_cli():
    cli = argparse.ArgumentParser()
    cli.add_argument("report_modle", type=str)
    cli.add_argument("report_openmm", type=str)
    cli.add_argument("--output-prefix",
                     required=True,
                     type=str)
    return cli


def plot_histogram(ax, df1, df2, label1, label2, compartments, bins=20, density=True, ylim=[0.0, 0.025], xlim=[0, 425], num_points=None):
    if num_points is None:
        num_points = min(len(df1), len(df2))

    ax.hist(df1["score"], label=label1, alpha=0.5, bins=bins, density=density)
    ax.hist(df2["score"], label=label2, alpha=0.5, bins=bins, density=density)

    title = f"Score (selected regions, {compartments}"
    if num_points is not None:
        title += f", {num_points} points)"
    else:
        title += ")"

    ax.set(title=title,
           xlabel="score (lower is better)",
           ylabel="relative freq",
           xlim=xlim,
           ylim=ylim)
    ax.legend()


def import_data(path):
    return pd.read_table(path, names=["chrom", "start", "end", "name", "score"])


if __name__ == "__main__":
    args = vars(make_cli().parse_args())

    df1 = import_data(args["report_openmm"])
    df2 = import_data(args["report_modle"])

    fig, axs = plt.subplots(2, 3, figsize=(5*6.4, 2*6.4))

    num_a_points = min(df1["name"].value_counts()["A"], df2["name"].value_counts()["A"])
    num_b_points = min(df1["name"].value_counts()["B"], df2["name"].value_counts()["B"])

    plot_histogram(axs[0][0], df1, df2, "OpenMM", "MoDLE", "A + B")
    plot_histogram(axs[0][1], df1.loc[df1["name"] == "A"], df2.loc[df2["name"] == "A"], "OpenMM", "MoDLE", "A only", num_points=num_a_points)
    plot_histogram(axs[0][2], df1.loc[df1["name"] == "B"], df2.loc[df2["name"] == "B"], "OpenMM", "MoDLE", "B only", num_points=num_b_points)

    plot_histogram(axs[1][0], df1.loc[df1["name"] == "A"], df1.loc[df1["name"] == "B"], "OpenMM (A)", "OpenMM (B)", "OpenMM", num_points=None)
    plot_histogram(axs[1][1], df2.loc[df2["name"] == "A"], df2.loc[df2["name"] == "B"], "MoDLE (A)", "MoDLE (B)", "MoDLE", num_points=None)

    avg_scores = {}

    avg_scores["OpenMM A"] = df1.loc[df1["name"] == "A", "score"].mean()
    avg_scores["OpenMM B"] = df1.loc[df1["name"] == "B", "score"].mean()
    avg_scores["MoDLE A"] = df2.loc[df2["name"] == "A", "score"].mean()
    avg_scores["MoDLE B"] = df2.loc[df2["name"] == "B", "score"].mean()

    axs[-1][-1].bar(avg_scores.keys(), avg_scores.values())
    axs[-1][-1].set(title="Avg. score", ylabel="Score (lower is better)")

    output_prefix = args["output_prefix"]
    fig.savefig(f"{output_prefix}.png", dpi=600)
    fig.savefig(f"{output_prefix}.svg")
