#!/usr/bin/env python3
import argparse
import multiprocessing as mp
import os
import sys

import cooler
import numpy as np
import pandas as pd
from polychrom.contactmaps import monomerResolutionContactMapSubchains
from polychrom.hdf5_format import list_URIs, load_URI


def make_cli():
    cli = argparse.ArgumentParser()
    cli.add_argument("--input-folders",
                     nargs="+",
                     required=True,
                     type=str)
    cli.add_argument("--output-name",
                     required=True,
                     type=str)
    cli.add_argument("--threads",
                     default=mp.cpu_count(),
                     type=int)
    cli.add_argument("--bin-size",
                     required=True,
                     type=int)
    cli.add_argument("--chrom-sizes",
                     required=True,
                     type=str,
                     help="Path to chrom.sizes file.")
    cli.add_argument("--chrom-names",
                     nargs="+",
                     type=str,
                     required=True,
                     help="Name of the chromosomes to which files inside --input-folder refer to.")
    cli.add_argument("--offsets",
                     nargs="+",
                     type=int,
                     required=True,
                     help="Offset to apply to contacts before building the contact matrix.")

    return cli


def matrix_to_df(matrix, bins, chrom_name, offset, bin_size):
    assert matrix.shape[0] == matrix.shape[1]
    assert chrom_name in bins["chrom"].unique()
    first_bin, *_, last_bin = bins.index[bins["chrom"] == chrom_name].to_list()

    bin_offset = first_bin + (offset // bin_size)
    assert bin_offset <= last_bin

    num_bins = matrix.shape[0]
    pixels = pd.DataFrame({"bin1_id": np.repeat(np.arange(num_bins), num_bins) + bin_offset,
                           "bin2_id": np.tile(np.arange(num_bins), num_bins) + bin_offset,
                           "count": matrix.flatten()})
    pixels = pixels[pixels["count"] != 0]
    return pixels[pixels["bin1_id"] <= pixels["bin2_id"]].sort_values(["bin1_id", "bin2_id"])


def make_bins(path_to_chrom_sizes, bin_size):
    chrom_sizes = pd.read_table(path_to_chrom_sizes, names=["name", "length"]).set_index("name")["length"]
    assert all((name in chrom_sizes for name in chrom_names))

    return cooler.binnify(chrom_sizes, bin_size)


if __name__ == "__main__":
    args = make_cli().parse_args()

    assert len(args.input_folders) == len(args.chrom_names)
    assert len(args.input_folders) == len(args.offsets)

    assert all((offset >= 0 for offset in args.offsets))

    bin_size = args.bin_size
    chrom_names = args.chrom_names
    offsets = args.offsets

    bins = make_bins(args.chrom_sizes, bin_size)

    pixels = []
    for folder, chrom_name, offset in zip(args.input_folders, chrom_names, offsets):
        assert chrom_name in bins["chrom"].unique()

        URIs = list_URIs(folder)
        size = load_URI(URIs[0])["pos"].shape[0]

        print(f"Reading {len(URIs)} chunks for \"{chrom_name}\" from folder {folder} using size={size}...")
        matrix = monomerResolutionContactMapSubchains(filenames=URIs,
                                                      mapStarts=range(0, size, size),
                                                      mapN=size,
                                                      cutoff=6,
                                                      n=args.threads,
                                                      loadFunction=lambda x: load_URI(x)["pos"])

        pixels.append(matrix_to_df(matrix, bins, chrom_name, offset, bin_size))
        print(f"Read {sum(matrix.flatten())} contacts into a {matrix.shape[0]}x{matrix.shape[1]} matrix.")

    tot_contacts = sum((sum(m["count"]) for m in pixels))
    pixels = pd.concat(pixels, sort=True)
    tot_contacts_after_concat = sum(pixels["count"])

    if tot_contacts != tot_contacts_after_concat:
        sys.exit(
            f"It appears that some contacts have been lost after merging contacts (expected {tot_contacts}, got {pixels}).\n"
            f"This usually means that the input folders contained simulations for overlapping ranges")

    if os.path.dirname(args.output_name) not in ["", "."]:
        os.makedirs(os.path.dirname(args.output_name), exist_ok=True)

    print(f"Writing contacts to file {args.output_name}...")
    cooler.create_cooler(args.output_name, bins, pixels)
    print("DONE!")
