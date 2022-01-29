#!/usr/bin/env python3

import argparse
import itertools
from collections import namedtuple

import cooler
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.cm import register_cmap
from matplotlib.patches import Rectangle

# Source: https://github.com/open2c/cooltools/blob/master/cooltools/lib/plotting.py
PALETTES = {
    "fall": np.array(
        (
            (255, 255, 255),
            (255, 255, 204),
            (255, 237, 160),
            (254, 217, 118),
            (254, 178, 76),
            (253, 141, 60),
            (252, 78, 42),
            (227, 26, 28),
            (189, 0, 38),
            (128, 0, 38),
            (0, 0, 0),
        )
    )
            / 255,
    "fall_mod": np.array(
        (
            (255, 255, 255),
            (255, 255, 204),
            (255, 237, 160),
            (254, 217, 118),
            (254, 178, 76),
            (253, 141, 60),
            (252, 78, 42),
            (227, 26, 28),
            (227, 26, 28),
            (227, 26, 28),
            (189, 0, 38),
        )
    )
                / 255,
    "blues": np.array(
        (
            (255, 255, 255),
            (180, 204, 225),
            (116, 169, 207),
            (54, 144, 192),
            (5, 112, 176),
            (4, 87, 135),
            (3, 65, 100),
            (2, 40, 66),
            (1, 20, 30),
            (0, 0, 0),
        )
    )
             / 255,
    "acidblues": np.array(
        (
            (255, 255, 255),
            (162, 192, 222),
            (140, 137, 187),
            (140, 87, 167),
            (140, 45, 143),
            (120, 20, 120),
            (90, 15, 90),
            (60, 10, 60),
            (30, 5, 30),
            (0, 0, 0),
        )
    )
                 / 255,
    "nmeth": np.array(
        (
            (236, 250, 255),
            (148, 189, 217),
            (118, 169, 68),
            (131, 111, 43),
            (122, 47, 25),
            (41, 0, 20),
        )
    )
             / 255,
}


def list_to_colormap(color_list, name=None):
    color_list = np.array(color_list)
    if color_list.min() < 0:
        raise ValueError("Colors should be 0 to 1, or 0 to 255")
    if color_list.max() > 1.0:
        if color_list.max() > 255:
            raise ValueError("Colors should be 0 to 1 or 0 to 255")
        else:
            color_list = color_list / 255.0
    return mpl.colors.LinearSegmentedColormap.from_list(name, color_list, 256)


def get_cmap(name):
    is_reversed = name.endswith("_r")
    try:
        if is_reversed:
            pal = PALETTES[name[:-2]][::-1]
        else:
            pal = PALETTES[name]
    except KeyError:
        raise ValueError('Palette not found "{}"'.format(name))
    return list_to_colormap(pal)


def _register_cmaps():
    for name, pal in PALETTES.items():
        register_cmap(name, list_to_colormap(pal))
        register_cmap(name + "_r", list_to_colormap(pal[::-1]))


def make_cli():
    # https://stackoverflow.com/questions/52132076/argparse-action-or-type-for-comma-separated-list
    class SplitArgs(argparse.Action):
        def __call__(self, parser, namespace, values, option_string=None):
            setattr(namespace, self.dest, values.split(","))

    class SplitIntervals(argparse.Action):
        def __call__(self, parser, namespace, values, option_string=None):
            values = values.split(",")
            values = [[float(n) for n in v.split("-")] for v in values]
            setattr(namespace, self.dest, values)

    class ParseUCSCString(argparse.Action):
        def __call__(self, parser, namespace, values, option_string=None):
            chrom, _, coord = values.partition(":")
            start, _, end = coord.partition("-")
            setattr(namespace, self.dest, [chrom, int(start), int(end)])

    cli = argparse.ArgumentParser()
    cli.add_argument("--path-to-hic-cool",
                     required=True,
                     type=str,
                     help="Path to HiC matrix in cool format (mcool URI syntax supported).")
    cli.add_argument("--path-to-md-cool",
                     required=True,
                     type=str,
                     help="Path to MD matrix in cool format (mcool URI syntax supported).")
    cli.add_argument("--path-to-modle-cool",
                     required=True,
                     type=str,
                     help="Path to MoDLE matrix in cool format (mcool URI syntax supported).")
    cli.add_argument("--path-to-microc-cool",
                     required=True,
                     type=str,
                     help="Path to microC matrix in cool format (mcool URI syntax supported).")
    cli.add_argument("--path-to-chrom-subranges-hic-comparison-bed",
                     required=True,
                     type=str,
                     help="Path to a BED file with the chromosome subranges to use for plotting the comparison with HiC data.")
    cli.add_argument("--color-scale-intervals",
                     action=SplitIntervals,
                     type=str,
                     default=[(0.0, 1.0),
                              (0.0, 1.0),
                              (0.0, 0.3),
                              (0.0, 0.01)],
                     help="Comma separated list of color scale intervals (e.g. 0.2-0.8,0.1-0.7...). Intervals refer to HiC, MD, MoDLE and microC dataset respectively.")
    cli.add_argument("--color-maps",
                     action=SplitArgs,
                     type=str,
                     default=["fall", "fall", "fall_mod", "fall"],
                     help="Comma separated list of colormaps to use for HiC, MD, MoDLE and microC dataset respectively.")
    cli.add_argument("--chrom-subrange-microc-comparison-ucsc",
                     type=str,
                     default=["chr7", 15000000, 16500000],
                     action=ParseUCSCString,
                     help="Genomic coordinates to use for the comparison with microC data (e.g. chr7:15000000-16500000).")
    cli.add_argument("--output-name",
                     required=True,
                     type=str)
    cli.add_argument("--labels",
                     action=SplitArgs,
                     default=["HiC", "MD", "MoDLE", "microC"],
                     type=str,
                     help="Comma separated list of plot labels (one per contact matrix).")
    cli.add_argument("--prefer-balanced-contacts",
                     type=bool,
                     default=True,
                     action=argparse.BooleanOptionalAction,
                     help="Prefer using balanced contacts whenever available.")
    return cli


def import_contact_matrix_cool(path_to_cool, chrom, start, end, prefer_balanced):
    assert start < end
    c = cooler.Cooler(path_to_cool)
    coord = f"{chrom}:{start}-{end}"
    if prefer_balanced:
        try:
            return np.nan_to_num(c.matrix(balance=prefer_balanced, as_pixels=False).fetch(coord).astype(float), nan=0)
        except ValueError:
            pass
    return c.matrix(balance=False, as_pixels=False).fetch(coord).astype(int)


def import_contacts(chrom_ranges, path_to_matrix, normalize=True, log_transform=False, prefer_balanced=True):
    matrices = {}
    if isinstance(chrom_ranges, list) or isinstance(chrom_ranges, tuple):
        assert len(chrom_ranges) == 3
        chrom_ranges = pd.DataFrame({"chrom": [chrom_ranges[0]],
                                     "chrom_start": [chrom_ranges[1]],
                                     "chrom_end": [chrom_ranges[2]]})

    for _, coords in chrom_ranges.iterrows():
        chrom = coords.chrom
        start = coords.chrom_start
        end = coords.chrom_end

        matrix = import_contact_matrix_cool(path_to_matrix, chrom, start, end, prefer_balanced)

        if normalize:
            matrix = matrix.astype(float) / np.max(matrix.ravel())
        if log_transform:
            matrix = np.log(matrix)

        matrices[tuple([chrom, start, end])] = matrix

    return matrices


def init_figure(labels_left_side, labels_right_side, chrom_ranges_left, chrom_range_right):
    """
    This function will initialize and empty figure with 12 subplots organized in two groups:
     - On the left side there's 3x3 group. This group is meant to hold the hic/md/modle
       comparison for 3 genomic regions (1 region per column)
     - On the right side the layout is 2 + 1: 2 smaller plots for modle/md data and a
       larger plot for microC

    This function returns the figure object as well as the axes as a dictionary.
    Axes are indexed with a tuple consisting of: [data_type, chrom, start, end]
    Example axs["hic", "chr1", 0, 1_000_000]
    """
    assert len(labels_left_side) == 3
    assert len(labels_right_side) == 3
    assert len(chrom_ranges_left) == 3

    fig = plt.figure(figsize=(6.4 * 5, 6.4 * 3), constrained_layout=True)
    # Create a gridspec such that the first three columns take up as much space as the last two
    gs = fig.add_gridspec(nrows=3, ncols=5, width_ratios=[1, 1, 1, 3 / 2, 3 / 2])
    T = namedtuple("Axes", "left right")
    axs = T({}, {})
    # Create the 3x3 subplots in the occupying the left half of the plot
    i = 0
    for label in labels_left_side:
        for (_, (name, start, end)) in chrom_ranges_left.iterrows():
            key = tuple([label, name, start, end])
            assert key not in axs.left
            axs.left[key] = fig.add_subplot(gs[int(i // 3), i % 3])
            i += 1

    # Create the 2+1 subplots in the occupying the right half of the plot
    key = tuple([labels_right_side[0], *chrom_range_right])
    assert key not in axs.right
    axs.right[key] = fig.add_subplot(gs[0, 3])

    key = tuple([labels_right_side[1], *chrom_range_right])
    assert key not in axs.right
    axs.right[key] = fig.add_subplot(gs[0, 4])

    key = tuple([labels_right_side[2], *chrom_range_right])
    assert key not in axs.right
    axs.right[key] = fig.add_subplot(gs[1:, 3:])

    return fig, axs


def draw_custom_legend(ax, cmap, vmin=0.0, vmax=1.0):
    if isinstance(cmap, str):
        cmap = get_cmap(cmap)
    custom_lines = [Rectangle((0.0, 0.0), 1, 1, facecolor=cmap(0.0), edgecolor="black", fill=True),
                    Rectangle((0.0, 0.0), 1, 1, facecolor=cmap(1.0), edgecolor="black", fill=True)]

    ax.legend(custom_lines, [str(vmin), str(vmax)], loc="lower left")


def plot_hic_comparison(axs, matrices, labels, color_maps, color_scale_intervals):
    keys = tuple(["hic", "md", "modle"])
    for coord, data in matrices.items():
        start = coord[1] / 1.0e6
        end = coord[2] / 1.0e6
        for key in keys:
            vmin, vmax = color_scale_intervals[key]
            matrix = getattr(matrices[tuple([*coord])], key)
            cmap = color_maps[key]
            key = tuple([key, *coord])
            axs[key].imshow(matrix,
                            cmap=cmap,
                            extent=(start, end,
                                    end, start),
                            vmin=vmin,
                            vmax=vmax,
                            origin="upper")
            draw_custom_legend(axs[key], cmap, vmin, vmax)

    # Add column titles
    for name, start, end in matrices.keys():
        key = tuple([keys[0], name, start, end])
        chrom = name.removeprefix("chrom").removeprefix("chr")
        start = f"{start / 1.0e6:.2f}".rstrip("0").rstrip(".")
        end = f"{end / 1.0e6:.2f}".rstrip("0").rstrip(".")
        axs[key].set_title(f"Chrom. {chrom} {start}-{end} Mbp", fontsize=18)

    # Add row titles
    coord = list(matrices.keys())[0]
    for key in keys:
        axs[tuple([key, *coord])].set_ylabel(labels[key], rotation=90, fontsize=24)


def plot_microc_comparison(axs, matrices, coord, labels, color_maps, color_scale_intervals):
    keys = tuple(["microc", "md", "modle"])

    matrices = matrices[tuple(coord)]
    chrom = coord[0]
    start = coord[1] / 1.0e6
    end = coord[2] / 1.0e6
    for key in keys:
        vmin, vmax = color_scale_intervals[key]
        matrix = getattr(matrices, key)
        cmap = color_maps[key]
        key = tuple([key, *coord])
        axs[key].imshow(matrix,
                        cmap=cmap,
                        extent=(start, end,
                                end, start),
                        vmin=vmin,
                        vmax=vmax,
                        origin="upper")

        draw_custom_legend(axs[key], cmap, vmin, vmax)

    # Add subplot titles
    chrom = chrom.removeprefix("chrom").removeprefix("chr")
    start = f"{start:.2f}".rstrip("0").rstrip(".")
    end = f"{end:.2f}".rstrip("0").rstrip(".")
    for key in keys:
        label = labels[key]
        key = tuple([key, *coord])

        axs[key].set_title(f"{label} - Chrom. {chrom} ({start}-{end})", fontsize=18)


def add_diagonal_signal(matrix, n=15, exp=1.33):
    max_ = np.inf
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            if abs(i - j) < n:
                matrix[i, j] += (n - abs(i - j)) ** exp

            if i == j:
                max_ = min(max_, matrix[i, j])
    return np.clip(matrix, 0, max_)


def import_matrices(chrom_ranges, path_to_hic_cool, path_to_md_cool, path_to_modle_cool, path_to_microc_cool,
                    prefer_balanced):
    hic_matrices = import_contacts(chrom_ranges, path_to_hic_cool, prefer_balanced=prefer_balanced)
    md_matrices = import_contacts(chrom_ranges, path_to_md_cool, prefer_balanced=prefer_balanced)
    modle_matrices = import_contacts(chrom_ranges, path_to_modle_cool, prefer_balanced=prefer_balanced)
    microc_matrices = import_contacts(chrom_ranges, path_to_microc_cool, prefer_balanced=prefer_balanced)

    MatrixTuple = namedtuple("MatrixTuple", "hic md modle microc")
    return {k: MatrixTuple(hic_matrices[k],
                           md_matrices[k],
                           modle_matrices[k],
                           microc_matrices[k]) for k in hic_matrices.keys()}


if __name__ == "__main__":
    _register_cmaps()

    args = make_cli().parse_args()

    assert len(args.labels) == 4
    assert len(args.color_scale_intervals) == 4
    assert len(args.color_maps) == 4
    # Assert color scale intervals are between 0 and 1
    assert sum((0 <= n <= 1 for n in itertools.chain.from_iterable(args.color_scale_intervals))) == len(
        args.labels) * 2

    # Index color scale intervals and labels by their respective datatype
    data_types = ["hic", "md", "modle", "microc"]
    color_scale_intervals = {k: interval for k, interval in zip(data_types, args.color_scale_intervals)}
    color_maps = {k: cmap for k, cmap in zip(data_types, args.color_maps)}
    labels = {k: label for k, label in zip(data_types, args.labels)}
    assert len(color_scale_intervals) == len(args.color_scale_intervals)
    assert len(labels) == len(args.labels)

    chrom_ranges_left = pd.read_table(args.path_to_chrom_subranges_hic_comparison_bed,
                                      names=["chrom", "chrom_start", "chrom_end"],
                                      usecols=range(0, 3))

    matrices = import_matrices(chrom_ranges_left,
                               args.path_to_hic_cool,
                               args.path_to_md_cool,
                               args.path_to_modle_cool,
                               args.path_to_microc_cool,
                               prefer_balanced=args.prefer_balanced_contacts)

    chrom_ranges_right = args.chrom_subrange_microc_comparison_ucsc
    fig, axs = init_figure(["hic", "md", "modle"],
                           ["modle", "md", "microc"],
                           chrom_ranges_left,
                           chrom_ranges_right)
    plot_hic_comparison(axs.left,
                        matrices,
                        labels,
                        color_maps,
                        color_scale_intervals)

    matrices = import_matrices(chrom_ranges_right,
                               args.path_to_hic_cool,
                               args.path_to_md_cool,
                               args.path_to_modle_cool,
                               args.path_to_microc_cool,
                               prefer_balanced=args.prefer_balanced_contacts)

    plot_microc_comparison(axs.right,
                           matrices,
                           chrom_ranges_right,
                           labels,
                           color_maps,
                           color_scale_intervals)

    fig.savefig(args.output_name)
