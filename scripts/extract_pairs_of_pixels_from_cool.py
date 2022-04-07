#!/usr/bin/env python3

# Important! This script is not meant to be used on large regions!

import argparse
import sys

import cooler
import numpy as np
import pandas as pd


def make_cli():
    cli = argparse.ArgumentParser()
    cli.add_argument("--ref-matrix",
                     type=str,
                     required=True,
                     help="Path to the reference matrix in cool or mcool format.")
    cli.add_argument("--tgt-matrix",
                     type=str,
                     required=True,
                     help="Path to the target matrix in cool or mcool format.")
    cli.add_argument("--chrom-ranges-bed",
                     type=str,
                     required=True,
                     help="Path to a BED file with the genomic coordinates of the regions to correlate.")
    cli.add_argument("--bin-size",
                     type=int,
                     required=True,
                     help="Matrix bin size.")
    cli.add_argument("--diagonal-width",
                     type=int,
                     required=True,
                     help="Diagonal width in bp.")
    return cli


def import_contacts(cooler, chrom, start, end, diagonal_width):
    coord = f"{chrom}:{start}-{end}"
    bins = cooler.bins().fetch(coord)[["chrom", "start", "end"]]
    pixels = bins.merge(bins, how="cross", suffixes=("1", "2"))
    pixels = pixels[(pixels["start2"] >= pixels["start1"]) &
                    (pixels["end2"] - pixels["end1"] <= diagonal_width)]

    try:
        counts = cooler.matrix(as_pixels=True,
                               balance=True,
                               join=True).fetch(coord)
    except ValueError:
        counts = cooler.matrix(as_pixels=True,
                               balance=False,
                               join=True).fetch(coord)
        counts["balanced"] = np.nan

    pixels = pixels.merge(counts,
                          on=tuple(pixels.columns),
                          how="left",
                          suffixes=("", ""),
                          indicator=True)

    pixels.loc[pixels["_merge"] != "both", "count"] = 0
    pixels.loc[pixels["_merge"] != "both", "balanced"] = 0

    return pixels.drop("_merge", axis="columns")


def path_to_uri(path, bin_size):
    if path.endswith(".mcool"):
        return path + f"::/resolutions/{bin_size}"
    return path


if __name__ == "__main__":
    args = make_cli().parse_args()

    bin_size = args.bin_size
    diagonal_width = args.diagonal_width

    ref = cooler.Cooler(path_to_uri(args.ref_matrix, bin_size))

    tgt = cooler.Cooler(path_to_uri(args.tgt_matrix, bin_size))

    assert ref.binsize == bin_size
    assert tgt.binsize == bin_size

    ranges = pd.read_table(args.chrom_ranges_bed,
                           usecols=range(0, 3),
                           names=["chrom", "start", "end"])
    ranges.drop_duplicates(inplace=True)

    buff = []

    header = ["chrom1", "start1", "end1",
              "chrom2", "start2", "end2",
              "count1", "count2",
              "balanced1", "balanced2"]

    cols_for_join = header[:6]

    print("\t".join(header))
    for chrom, start, end in zip(ranges["chrom"].to_numpy(),
                                 ranges["start"].to_numpy(),
                                 ranges["end"].to_numpy()):
        assert start <= end

        m1 = import_contacts(ref,
                             chrom, start, end,
                             diagonal_width)
        m2 = import_contacts(tgt,
                             chrom, start, end,
                             diagonal_width)

        assert m1.shape == m2.shape

        m = m1.merge(m2,
                    how="inner",
                    on=cols_for_join,
                    sort=True,
                    suffixes=("1", "2"))[header]

        m["count1"] = m["count1"].astype(int)
        m["count2"] = m["count2"].astype(int)

        m.to_csv(sys.stdout, sep="\t", na_rep="nan", header=False, index=False)
