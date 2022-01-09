#!/usr/bin/env python3

import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.cm import register_cmap
import numpy as np
from PIL import Image, ImageOps
from pathlib import Path
import argparse

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


def make_cli():
    cli = argparse.ArgumentParser()
    cli.add_argument("--output-dir",
                     required=True,
                     type=str)
    cli.add_argument("inputs",
                     metavar="input-pics",
                     nargs="+",
                     type=str,
                     help="Path to one or more pictures to be processed")
    return cli


def import_image(path):
    img = ImageOps.invert(Image.open(path).convert('L'))
    return np.array(img, dtype=int)


if __name__ == "__main__":
    args = make_cli().parse_args()
    matrix_out_prefix = f"{args.output_dir}/matrices"
    plot_out_prefix = f"{args.output_dir}/plots"

    Path(matrix_out_prefix).mkdir(parents=True, exist_ok=True)
    Path(plot_out_prefix).mkdir(parents=True, exist_ok=True)

    for path in args.inputs:
        bname = Path(path).name
        stem = Path(bname).with_suffix("")

        print(f"Processing {bname}...")

        matrix = import_image(path)
        matrix_out_name = f"{matrix_out_prefix}/{stem}.tsv"
        np.savetxt(matrix_out_name, matrix, delimiter="\t", fmt="%d")

        fig, ax = plt.subplots(1)
        cmap = ax.imshow(matrix, cmap="fall", interpolation="nearest")
        fig.colorbar(cmap)
        fig.tight_layout()

        plot_out_name = f"{plot_out_prefix}/{stem}.png"
        fig.savefig(plot_out_name)
