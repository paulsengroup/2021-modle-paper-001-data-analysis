#!/usr/bin/env python3

# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import argparse
from collections import namedtuple

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from numpy.polynomial import Polynomial
from numpy.polynomial.polynomial import polyfit
from scipy.interpolate import UnivariateSpline


def make_cli():
    class SplitIntArgs(argparse.Action):
        def __call__(self, parser, namespace, values, option_string=None):
            setattr(namespace, self.dest, [int(v) for v in values.split(",")])

    cli = argparse.ArgumentParser()
    cli.add_argument("--modle-mt-report-tsv",
                     required=True,
                     type=str,
                     help="Path to MoDLE multi-threaded benchmark report in TSV format.")
    cli.add_argument("--md-cpu-report-tsv",
                     required=True,
                     type=str,
                     help="Path to molecular dynamics CPU benchmark report in TSV format.")
    cli.add_argument("--md-gpu-report-tsv",
                     required=True,
                     type=str,
                     help="Path to molecular dynamics GPU benchmark report in TSV format.")

    cli.add_argument("--title",
                     type=str,
                     default="MoDLE vs MD - Memory resident size of simulations on artificial chromosomes",
                     help="Plot title.")
    cli.add_argument("--xaxis-label",
                     type=str,
                     default="Chromosome size (Mbp)",
                     help="x axis label.")
    cli.add_argument("--yaxis-label",
                     type=str,
                     default="Peak memory usage (MB)",
                     help="Y axis label.")

    cli.add_argument("--xaxis-interval-mbp",
                     default=[-2, 252],
                     action=SplitIntArgs,
                     help="Plotting interval for the x axis in Mbp (e.g. '-5,255').")
    cli.add_argument("--yaxis-interval-mb",
                     default=[0, 1100],
                     action=SplitIntArgs,
                     help="Plotting interval for the y axis in MB (e.g. '0,3500').")

    cli.add_argument("--chrom-size-range-mbp",
                     default=[1, 250],
                     action=SplitIntArgs,
                     help="Minimum and maximum chromosome size used for plotting (e.g. '1,250').\n"
                          "This is mostly useful to limit the extrapolation extent.")

    cli.add_argument("--color-style",
                     type=str,
                     default="seaborn-paper",
                     help="Matplotlib style to use for plotting.")

    cli.add_argument("--extrapolate",
                     default=False,
                     action="store_true")
    cli.add_argument("--output-prefix",
                     required=True,
                     type=str,
                     help="Prefix name to use for output.")
    return cli


def plot_line(ax, x, y, xmax, label, color, extrapolate, alpha=0.5):
    l1 = ax.plot(x[x <= xmax],
                 y[x <= xmax],
                 label=label,
                 color=color,
                 alpha=alpha,
                 linestyle="-")
    if extrapolate:
        l2 = ax.plot(x[x > xmax],
                     y[x > xmax],
                     label=f"{label} (extrapolated)",
                     color=color,
                     alpha=alpha,
                     linestyle="--")
        return l1 + l2

    return l1


def plot_scatter(ax, x, y, label, color, marker):
    return ax.scatter(x[~np.isnan(x)],
                      y[~np.isnan(y)],
                      label=label,
                      color=color,
                      marker=marker)


def gen_interpolation_fx(x, y, mode="spline", **kwargs):
    if mode == "spline":
        bbox = [x[0], x[-1]]
        k = 3
        s = None
        if "bbox" in kwargs:
            bbox = kwargs["bbox"]
        if "k" in kwargs:
            k = kwargs["k"]
        if "s" in kwargs:
            s = kwargs["s"]
        return UnivariateSpline(x[x <= bbox[1]], y[x <= bbox[1]], k=k, s=s, bbox=bbox, ext="extrapolate")
    else:
        assert mode == "polynomial"
        deg = kwargs["deg"]
        return Polynomial(polyfit(x, y, deg=deg))


def import_report(path_to_report, min_chrom_size, max_chrom_size):
    df = pd.read_table(path_to_report)
    df["chrom_size"] /= 1.0e6
    return df[(df["chrom_size"] >= min_chrom_size) &
              (df["chrom_size"] <= max_chrom_size)].set_index("chrom_size")


def compute_data_extent(modle_mt_report, md_cpu_report, md_gpu_report, min_chrom_size, max_chrom_size):
    assert min_chrom_size <= max_chrom_size

    ExtentT = namedtuple("DataExtent", ["modle_x0", "md_cpu_x0", "md_gpu_x0",
                                        "modle_x1", "md_cpu_x1", "md_gpu_x1"])
    extent = ExtentT
    extent.modle_x0 = int(max(modle_mt_report.index.min(), min_chrom_size))
    extent.md_cpu_x0 = int(max(md_cpu_report.index.min(), min_chrom_size))
    extent.md_gpu_x0 = int(max(md_gpu_report.index.min(), min_chrom_size))

    extent.modle_x1 = int(min(modle_mt_report.index.max(), max_chrom_size))
    extent.md_cpu_x1 = int(min(md_cpu_report.index.max(), max_chrom_size))
    extent.md_gpu_x1 = int(min(md_gpu_report.index.max(), max_chrom_size))

    return extent


def generate_x_tracks(extent, extrapolate, points=2500):
    XTracksT = namedtuple("XTracks", ["modle", "md_cpu", "md_gpu"])

    if extrapolate:
        lb = min(extent.modle_x0, extent.md_cpu_x0, extent.md_gpu_x0)
        ub = max(extent.modle_x1, extent.md_cpu_x1, extent.md_gpu_x1)
        x = np.linspace(lb, ub, points)
        return XTracksT(x, x, x)

    tracks = XTracksT
    tracks.modle = np.linspace(extent.modle_x0, extent.modle_x1, points)
    tracks.md_cpu = np.linspace(extent.md_cpu_x0, extent.md_cpu_x1, points)
    tracks.md_gpu = np.linspace(extent.md_gpu_x0, extent.md_gpu_x1, points)

    return tracks


if __name__ == "__main__":
    args = make_cli().parse_args()

    min_chrom_size, max_chrom_size = args.chrom_size_range_mbp
    modle_mt_report = import_report(args.modle_mt_report_tsv, min_chrom_size, max_chrom_size)
    md_cpu_report = import_report(args.md_cpu_report_tsv, min_chrom_size, max_chrom_size)
    md_gpu_report = import_report(args.md_gpu_report_tsv, min_chrom_size, max_chrom_size)

    modle_mt_zoom_report = modle_mt_report.loc[md_cpu_report.index, :]

    x0, x1 = args.xaxis_interval_mbp
    y0, y1 = args.yaxis_interval_mb

    assert x0 < x1
    assert y0 < y1
    assert args.color_style in plt.style.library, "Unknown color style!"

    mpl.style.use(args.color_style)

    fig, ax1 = plt.subplots(1, 1)
    ax2 = ax1.inset_axes([0.65, 0.02, 0.34, 0.3])

    modle_mt_fx = gen_interpolation_fx(modle_mt_report.index,
                                       modle_mt_report["max_resident_mem_kb_median"] / 1000,
                                       mode="polynomial", deg=1)
    md_cpu_fx = gen_interpolation_fx(md_cpu_report.index,
                                     md_cpu_report["max_resident_mem_kb_median"] / 1000,
                                     mode="polynomial",
                                     deg=1)
    md_gpu_fx = gen_interpolation_fx(md_gpu_report.index,
                                     md_gpu_report["max_resident_mem_kb_median"] / 1000,
                                     mode="polynomial",
                                     deg=1)
    modle_mt_zoom_fx = gen_interpolation_fx(modle_mt_zoom_report.index,
                                            modle_mt_zoom_report["max_resident_mem_kb_median"] / 1000,
                                            mode="polynomial", deg=1)

    tracks = []

    tracks.append(
        plot_scatter(ax1,
                     modle_mt_report.index,
                     modle_mt_report["max_resident_mem_kb_median"] / 1000,
                     "MoDLE (MT)", color="C0",
                     marker="o"))

    tracks.append(
        plot_scatter(ax1,
                     md_gpu_report.index,
                     md_gpu_report["max_resident_mem_kb_median"] / 1000,
                     "OpenMM (GPU)",
                     color="C1",
                     marker="s"))

    tracks.append(
        plot_scatter(ax2,
                     md_cpu_report.index,
                     md_cpu_report["max_resident_mem_kb_median"] / 1000,
                     "OpenMM (CPU)",
                     color="C7",
                     marker="P"))
    plot_scatter(ax2,
                 modle_mt_zoom_report.index,
                 modle_mt_zoom_report["max_resident_mem_kb_median"] / 1000,
                 "MoDLE (MT) zoomed", color="C0",
                 marker="o")

    extent = compute_data_extent(modle_mt_report, md_cpu_report, md_gpu_report,
                                 min_chrom_size, max_chrom_size)
    x_tracks = generate_x_tracks(extent, args.extrapolate)

    if args.extrapolate:
        x1 = int(min(modle_mt_report.index.max(),
                     md_gpu_report.index.max(), max_chrom_size))

    modle_mt_y = modle_mt_fx(x_tracks.modle)
    md_cpu_y = md_cpu_fx(x_tracks.md_cpu)
    md_gpu_y = md_gpu_fx(x_tracks.md_gpu)
    modle_mt_zoom_cpu_y = modle_mt_zoom_fx(x_tracks.md_cpu)

    plot_line(ax1, x_tracks.modle, modle_mt_y,
              extent.modle_x1,
              "MoDLE (MT)",
              color="C0",
              extrapolate=args.extrapolate)
    plot_line(ax1, x_tracks.md_gpu, md_gpu_y,
              extent.md_gpu_x1,
              "OpenMM (GPU)",
              color="C1",
              extrapolate=args.extrapolate)

    plot_line(ax2, x_tracks.md_cpu, modle_mt_zoom_cpu_y,
              extent.md_cpu_x1,
              "MoDLE (MT)",
              color="C0",
              extrapolate=args.extrapolate)
    plot_line(ax2, x_tracks.md_cpu, md_cpu_y,
              extent.md_cpu_x1,
              "OpenMM (CPU)",
              color="C7",
              extrapolate=args.extrapolate)

    tracks = [t for t in tracks if not str(t).endswith("(extrapolated))")]
    if args.extrapolate:
        tracks += [Line2D([0], [0], color="tab:grey", linestyle="--", label="Extrapolated")]
    labels = [t.get_label() for t in tracks]

    ax1.set(title=args.title,
            xlabel=args.xaxis_label,
            ylabel=args.yaxis_label,
            xlim=[x0, x1 + 2],
            ylim=[y0, y1])

    y0, y1 = ax2.get_ylim()
    y_pad = y0 * 0.15
    ax2.set(ylim=[y0 - y_pad, y1 + y_pad])

    ax1.legend(tracks, labels)
    ax2.tick_params(axis="x", direction="in", pad=-12)

    fig.tight_layout()
    assert args.output_prefix
    fig.savefig(f"{args.output_prefix}.png", dpi=1200)
    fig.savefig(f"{args.output_prefix}.svg")
