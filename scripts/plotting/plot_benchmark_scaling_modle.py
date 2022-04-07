#!/usr/bin/env python3

import argparse

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from numpy.polynomial import Polynomial
from numpy.polynomial.polynomial import polyfit
from scipy.interpolate import UnivariateSpline


def make_cli():
    class SplitFloatArgs(argparse.Action):
        def __call__(self, parser, namespace, values, option_string=None):
            setattr(namespace, self.dest, [float(v) for v in values.split(",")])

    cli = argparse.ArgumentParser()
    cli.add_argument("--modle-report-tsv",
                     required=True,
                     type=str,
                     help="Path to MoDLE benchmark report in TSV format.")
    cli.add_argument("--output-prefix",
                     required=True,
                     type=str,
                     help="Prefix name to use for output")
    cli.add_argument("--color-style",
                     type=str,
                     default="seaborn-paper",
                     help="Matplotlib style to use for plotting.")
    cli.add_argument("--title",
                     required=True,
                     type=str,
                     help="Plot title.")
    cli.add_argument("--xaxis-label",
                     type=str,
                     default="# CPU cores",
                     help="x axis label.")
    cli.add_argument("--yaxis-left-label",
                     type=str,
                     default="Wall clock (h)",
                     help="Left y axis label.")
    cli.add_argument("--yaxis-right-label",
                     type=str,
                     default="Peak memory usage (MB)",
                     help="Right y axis label.")
    cli.add_argument("--xaxis-interval",
                     default=[-0.5, 55],
                     action=SplitFloatArgs,
                     help="Plotting interval for the x axis (e.g. '-2,54').")
    cli.add_argument("--yaxis-left-interval",
                     default=[-0.1, 5.5],
                     action=SplitFloatArgs,
                     help="Plotting interval for the left y axis (e.g. '-0.1,6').")
    cli.add_argument("--yaxis-right-interval",
                     default=[0, 6500],
                     action=SplitFloatArgs,
                     help="Plotting interval for the right y axis (ratio) (e.g. '0,5000').")
    return cli


def plot_line(ax, x, y, label, **kwargs):
    kwargs1 = kwargs.copy()
    if "linestyle" not in kwargs:
        kwargs["linestyle"] = "-"

    return ax.plot(x,
                   y,
                   label=label,
                   **kwargs1)


def plot_scatter(ax, x, y, label, **kwargs):
    return ax.scatter(x[~np.isnan(x)],
                      y[~np.isnan(y)],
                      label=label,
                      **kwargs)


def gen_interpolation_fx(x, y, mode="spline", **kwargs):
    if mode == "spline":
        bbox = [x[0], x[-1]]
        return UnivariateSpline(x, y, bbox=bbox, **kwargs)

    assert mode == "polynomial"
    return Polynomial(polyfit(x, y, **kwargs))


if __name__ == "__main__":
    args = make_cli().parse_args()

    assert args.color_style in plt.style.library, "Unknown color style!"

    mpl.style.use(args.color_style)

    modle_report = pd.read_table(args.modle_report_tsv).set_index("ncores")
    modle_report["cpu_time_median"] = modle_report["system_time_median"] + modle_report["user_time_median"]

    fig, ax1 = plt.subplots(1, 1)

    modle_wall_clock = modle_report["wall_clock_median"] / 3600
    modle_cputime = modle_report["cpu_time_median"] / 3600
    modle_mem = modle_report["max_resident_mem_kb_median"] / 1024

    modle_wc_fx = gen_interpolation_fx(modle_report.index, modle_wall_clock, mode="polynomial",
                                       deg=modle_report.shape[0] - 3)
    modle_cput_fx = gen_interpolation_fx(modle_report.index, modle_cputime, mode="polynomial",
                                         deg=2)


    def modle_perfect_scaling_time_fx(nthreads):
        return (modle_report.loc[1, "wall_clock_median"] / 3600) / nthreads


    modle_mem_fx = gen_interpolation_fx(modle_report.index, modle_mem, mode="polynomial", deg=2)

    x = np.linspace(1, modle_report.index.max(), 2500)

    ax2 = ax1.twinx()
    tracks = []

    tracks.append(
        plot_scatter(ax1,
                     modle_report.index,
                     modle_wall_clock,
                     label="Wall clock",
                     color="C0",
                     marker="o"))
    tracks.append(
        plot_scatter(ax1,
                     modle_report.index,
                     modle_cputime,
                     label="CPU time",
                     color="C1",
                     marker="X"))
    tracks.append(
        plot_scatter(ax2,
                     modle_report.index,
                     modle_mem,
                     label="Memory",
                     color="C2",
                     marker="D"))

    plot_line(ax1,
              x,
              modle_wc_fx(x),
              "Wall clock",
              color="C0")
    plot_line(ax1,
              x,
              modle_cput_fx(x),
              "CPU time",
              color="C1")

    plot_line(ax1,
              x,
              modle_perfect_scaling_time_fx(x),
              "Wall clock (perfect scaling)",
              color="C0",
              linestyle="dotted")

    plot_line(ax2,
              x,
              modle_mem_fx(x),
              "Memory",
              color="C2")

    tracks += [Line2D([0], [0], color="C0", linestyle="dotted", label="Perfect scaling")]
    labels = [t.get_label() for t in tracks]

    ax1.set(title=args.title,
            xlabel="# CPU cores",
            ylabel="Wall clock (h)",
            xlim=args.xaxis_interval,
            ylim=args.yaxis_left_interval)

    ax2.set(ylabel="Peak memory usage (MB)",
            ylim=args.yaxis_right_interval)

    ax1.legend(tracks, labels)

    for (i, (xp, yp)) in enumerate(zip(modle_report.index, modle_wall_clock)):
        if i == 0:
            ax1.annotate(f"{yp * 60:.1f}m", (xp + 1, yp - 0.025), size=7)
            continue

        ax1.annotate(f"{yp * 60:.1f}m", (xp - 0.25, yp + 0.075), size=7)

    fig.tight_layout()
    assert args.output_prefix
    fig.savefig(f"{args.output_prefix}.png", dpi=1200)
    fig.savefig(f"{args.output_prefix}.svg")
