#!/usr/bin/env python3

import argparse

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from numpy.polynomial import Polynomial
from numpy.polynomial.polynomial import polyfit
from scipy.interpolate import UnivariateSpline


def make_cli():
    class SplitIntArgs(argparse.Action):
        def __call__(self, parser, namespace, values, option_string=None):
            setattr(namespace, self.dest, [int(v) for v in values.split(",")])

    cli = argparse.ArgumentParser()
    cli.add_argument("--modle-st-report-tsv",
                     required=True,
                     type=str,
                     help="Path to MoDLE single-threaded benchmark report in TSV format.")
    cli.add_argument("--modle-mt-report-tsv",
                     required=True,
                     type=str,
                     help="Path to MoDLE multi-threaded benchmark report in TSV format.")

    cli.add_argument("--title",
                     type=str,
                     default="MoDLE - Wall clock of simulations on artificial chromosomes",
                     help="Plot title.")
    cli.add_argument("--xaxis-label",
                     type=str,
                     default="Chromosome size (Mbp)",
                     help="x axis label.")
    cli.add_argument("--yaxis-label",
                     type=str,
                     default="Wall clock (s)",
                     help="Y axis label.")

    cli.add_argument("--xaxis-interval-mbp",
                     default=[-5, 505],
                     action=SplitIntArgs,
                     help="Plotting interval for the x axis in Mbp (e.g. '-5,255').")
    cli.add_argument("--yaxis-interval-sec",
                     default=[-10, 800],
                     action=SplitIntArgs,
                     help="Plotting interval for the y axis in seconds (e.g. '0,1250').")

    cli.add_argument("--color-style",
                     type=str,
                     default="seaborn-paper",
                     help="Matplotlib style to use for plotting.")

    cli.add_argument("--output-prefix",
                     required=True,
                     type=str,
                     help="Prefix name to use for output.")
    return cli


def plot_line(ax, x, y, label, color, alpha=0.5):
    return ax.plot(x,
                   y,
                   label=label,
                   color=color,
                   alpha=alpha,
                   linestyle="-")


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


def import_report(path_to_report):
    df = pd.read_table(path_to_report)
    df["chrom_size"] /= 1.0e6
    return df.set_index("chrom_size")


def compute_data_extent(modle_mt_report, min_chrom_size, max_chrom_size):
    assert min_chrom_size <= max_chrom_size

    return (int(max(modle_mt_report.index.min(), min_chrom_size)),
            int(min(modle_mt_report.index.max(), max_chrom_size)))


def generate_x_track(extent, extrapolate, points=2500):
    if extrapolate:
        lb = min(extent.modle_x0, extent.md_cpu_x0, extent.md_gpu_x0)
        ub = max(extent.modle_x1, extent.md_cpu_x1, extent.md_gpu_x1)
        return np.linspace(lb, ub, points)

    return np.linspace(extent.modle_x0, extent.modle_x1, points)


if __name__ == "__main__":
    args = make_cli().parse_args()

    modle_st_report = import_report(args.modle_st_report_tsv)
    modle_mt_report = import_report(args.modle_mt_report_tsv)

    x0, x1 = args.xaxis_interval_mbp
    y0, y1 = args.yaxis_interval_sec

    assert x0 < x1
    assert y0 < y1
    assert args.color_style in plt.style.library, "Unknown color style!"

    mpl.style.use(args.color_style)

    fig, ax = plt.subplots(1, 1)

    modle_st_fx = gen_interpolation_fx(modle_st_report.index,
                                       modle_st_report["wall_clock_median"],
                                       mode="polynomial", deg=1)
    modle_mt_fx = gen_interpolation_fx(modle_mt_report.index,
                                       modle_mt_report["wall_clock_median"],
                                       mode="polynomial", deg=1)

    tracks = []

    x_track = np.linspace(modle_st_report.index.min(), modle_st_report.index.max(), 2500)

    modle_st_y = modle_st_fx(x_track)
    modle_mt_y = modle_mt_fx(x_track)

    plot_line(ax, x_track, modle_st_y,
              "MoDLE (ST)",
              color="C6")
    plot_line(ax, x_track, modle_mt_y,
              "MoDLE (MT)",
              color="C0")

    tracks.append(
        plot_scatter(ax,
                     modle_mt_report.index,
                     modle_mt_report["wall_clock_median"],
                     "MoDLE (MT)", color="C0",
                     marker="o"))
    tracks.append(
        plot_scatter(ax,
                     modle_st_report.index,
                     modle_st_report["wall_clock_median"],
                     "MoDLE (ST)", color="C6",
                     marker="X"))

    labels = [t.get_label() for t in tracks]

    ax.set(title=args.title,
           xlabel=args.xaxis_label,
           ylabel=args.yaxis_label,
           xlim=[x0, x1],
           ylim=[y0, y1])

    ax.legend(tracks, labels)

    fig.tight_layout()
    assert args.output_prefix
    fig.savefig(f"{args.output_prefix}.png", dpi=1200)
    fig.savefig(f"{args.output_prefix}.svg")
