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


_register_cmaps()


def gridspec_inches(wcols, hrows, fig_kwargs={}):
    fig_height_inches = sum(hrows)

    fig_width_inches = sum(wcols)

    fig = plt.figure(
        figsize=(fig_width_inches, fig_height_inches),
        subplotpars=mpl.figure.SubplotParams(
            left=0, right=1, bottom=0, top=1, wspace=0, hspace=0.0
        ),
        # frameon=False,
        **fig_kwargs
    )
    fig.set_size_inches(fig_width_inches, fig_height_inches, forward=True)

    gs = mpl.gridspec.GridSpec(
        len(hrows),
        len(wcols),
        left=0,
        right=1,
        top=1,
        bottom=0,
        wspace=0,
        hspace=0,
        width_ratios=wcols,
        height_ratios=hrows,
    )

    return fig, gs


# https://stackoverflow.com/questions/52132076/argparse-action-or-type-for-comma-separated-list
class SplitArgs(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, values.split(','))


def make_cli():
    cli = argparse.ArgumentParser()
    cli.add_argument("--path-to-hic-tsv",
                     required=True,
                     type=str,
                     action=SplitArgs,
                     help="List of comma-separated paths to HiC matrices in TSV format.")
    cli.add_argument("--path-to-md-tsv",
                     required=True,
                     type=str,
                     action=SplitArgs,
                     help="List of comma-separated paths to MD matrices in TSV format.")
    cli.add_argument("--path-to-modle-cool",
                     required=True,
                     type=str,
                     help="Path to MoDLE's contact matrix in .cool format.")
    cli.add_argument("--path-to-chrom-subranges-bed",
                     required=True,
                     type=str,
                     help="Path to a BED file with the chromosome subranges corresponding to the TSV matrices "
                          "specified through --path-to-hic-tsv and --path-to-md-tsv")
    cli.add_argument("--output-name",
                     required=True,
                     type=str)
    return cli


def import_contact_matrix_tsv(path_to_tsv):
    data = np.genfromtxt(path_to_tsv, delimiter="\t", dtype=int)
    assert data.shape[0] == data.shape[1]  # Assert matrix is square
    return data


def import_contact_matrix_cool(path_to_cool, chrom, start, end):
    assert start < end
    c = cooler.Cooler(path_to_cool)
    coord = f"{chrom}:{start}-{end}"
    return c.matrix(balance=False, as_pixels=False).fetch(coord).astype(int)


def import_contacts(path_to_chrom_ranges, path_to_hic_tsv, path_to_md_tsv, path_to_modle_cool, normalize=True):
    data = {}
    for ((_, coords), hic, md, modle) in zip(path_to_chrom_ranges.iterrows(),
                                             path_to_hic_tsv,
                                             path_to_md_tsv,
                                             itertools.repeat(path_to_modle_cool)):
        chrom = coords.chrom
        start = coords.chrom_start
        end = coords.chrom_end

        hic_matrix = import_contact_matrix_tsv(hic)
        md_matrix = import_contact_matrix_tsv(md)
        modle_matrix = add_diagonal_signal(import_contact_matrix_cool(modle, chrom, start, end))

        if normalize:
            hic_matrix = hic_matrix.astype(float) / np.max(hic_matrix.ravel())
            md_matrix = md_matrix.astype(float) / np.max(md_matrix.ravel())
            modle_matrix = modle_matrix.astype(float) / np.max(modle_matrix.ravel())

        T = namedtuple("Contacts", ["hic", "md", "modle"])
        data[tuple([chrom, start, end])] = T(hic=hic_matrix, md=md_matrix, modle=modle_matrix)

    return data


def make_plot(matrices):
    nrows = len(matrices)
    ncols = len(tuple(matrices.values())[0])

    fig, axs = plt.subplots(nrows, ncols, figsize=((6.4 * ncols) + 1, 6.4 * nrows), constrained_layout=True)
    for (col, (coord, data)) in enumerate(matrices.items()):
        start = coord[1] / 1.0e6
        end = coord[2] / 1.0e6
        for row, matrix in enumerate(data):
            img = axs[row][col].imshow(matrix,
                                       cmap="fall",
                                       extent=(start, end,
                                               end, start),
                                       origin="upper")

    fig.colorbar(img, ax=axs.ravel().tolist())

    for (ax, (chrom, start, end)) in zip(axs[0], matrices.keys()):
        chrom_number = chrom.removeprefix("chrom").removeprefix("chr")
        start = f"{start / 1.0e6:.2f}".rstrip("0").rstrip(".")
        end = f"{end / 1.0e6:.2f}".rstrip("0").rstrip(".")
        ax.set_title(f"Chromosome {chrom_number} {start}-{end} Mbp", fontsize=24)

    for ax, row in zip(axs[:, 0], ["HiC", "MD", "MoDLE"]):
        ax.set_ylabel(row, rotation=90, fontsize=24)

    return fig


def add_diagonal_signal(matrix, n=15, exp=1.33):
    max_ = np.inf
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            if abs(i - j) < n:
                matrix[i, j] += (n - abs(i - j))**exp

            if i == j:
                max_ = min(max_, matrix[i, j])
    return np.clip(matrix, 0, max_)


if __name__ == "__main__":
    args = make_cli().parse_args()

    hic_tsv = tuple(args.path_to_hic_tsv)
    md_tsv = tuple(args.path_to_md_tsv)
    assert len(hic_tsv) == len(md_tsv)

    chrom_ranges = pd.read_table(args.path_to_chrom_subranges_bed,
                                 names=["chrom", "chrom_start", "chrom_end"],
                                 usecols=range(0, 3))
    assert chrom_ranges.shape[0] == len(hic_tsv)

    matrices = import_contacts(chrom_ranges, hic_tsv, md_tsv, args.path_to_modle_cool)

    fig = make_plot(matrices)
    fig.savefig(args.output_name)
