#!/usr/bin/env python3

# Important! This script is not meant to be used on large regions!

import argparse
from collections import namedtuple
from math import isnan
# import matplotlib.pyplot as plt

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
                     help="Path to a BED file with the genomic coordinates of the regions to score.")
    cli.add_argument("--sites-of-interest-bed",
                     type=str,
                     help="Path to a BED file in BED6+ format with the set of genomic coordinates of interest.")
    cli.add_argument("--bin-size",
                     type=int,
                     required=True,
                     help="Matrix bin size.")
    cli.add_argument("--diagonals-to-exclude",
                     type=int,
                     default=5,
                     help="Number of diagonals to exclude.")
    cli.add_argument("--direction",
                     type=str,
                     choices={"horizontal", "vertical"},
                     required=True,
                     help="Stripe direction.")
    cli.add_argument("--diagonal-width",
                     type=int,
                     required=True,
                     help="Diagonal width in bp.")
    return cli


# Source: https://stackoverflow.com/a/18081653
def kth_diag_indices(a, k):
    rows, cols = np.diag_indices_from(a)
    if k < 0:
        return rows[-k:], cols[:k]
    elif k > 0:
        return rows[:-k], cols[k:]
    else:
        return rows, cols


def import_contacts(cooler, chrom, start, end, diagonals_to_exclude):
    assert diagonals_to_exclude >= 0
    coord = f"{chrom}:{start}-{end}"
    data = cooler.matrix(balance=False,
                         as_pixels=False).fetch(coord).astype(float)

    for k in range(diagonals_to_exclude):
        rows_i, cols_i = kth_diag_indices(data, k)
        data[rows_i, cols_i] = np.nan

    # Copy upper triangle to lower triangle
    i_lower = np.tril_indices(data.shape[0], -1)
    data[i_lower] = data.T[i_lower]
    return data


def compute_score(v1, v2):
    ResultT = namedtuple("ResultT", ["oo", "oi", "io", "ii", "nans"])
    oo = 0
    oi = 0
    io = 0
    ii = 0
    nans = 0
    for p1, p2 in zip(v1, v2, strict=True):
        if isnan(p1) or isnan(p2):
            nans += 1
        elif p1 == 0 and p2 == 0:
            oo += 1
        elif p1 == 0:
            oi += 1
        elif p1 != 0 and p2 != 0:
            ii += 1
        else:
            io += 1
    return ResultT(oo, oi, io, ii, nans)


def path_to_uri(path, bin_size):
    if path.endswith(".mcool"):
        return path + f"::/resolutions/{bin_size}"
    return path


if __name__ == "__main__":
    args = make_cli().parse_args()

    bin_size = args.bin_size
    diagonal_width = args.diagonal_width
    num_diagonals_to_exclude = args.diagonals_to_exclude

    ref = cooler.Cooler(path_to_uri(args.ref_matrix, bin_size))
    tgt = cooler.Cooler(path_to_uri(args.tgt_matrix, bin_size))

    assert ref.binsize == bin_size
    assert tgt.binsize == bin_size

    ranges = pd.read_table(args.chrom_ranges_bed,
                           usecols=range(0, 3),
                           names=["chrom", "start", "end"])
    ranges.drop_duplicates(inplace=True)

    print(f"chrom1\tstart1\tend1\tchrom2\tstart2\tend2\too\toi\tio\tii\tnans")
    for chrom, start, end in zip(ranges["chrom"].to_numpy(),
                                 ranges["start"].to_numpy(),
                                 ranges["end"].to_numpy()):
        assert start <= end

        num_diagonals = int(min(np.ceil(diagonal_width / bin_size),
                                np.ceil((end - start) / bin_size)))

        m1 = import_contacts(ref, chrom, start, end, num_diagonals_to_exclude)
        m2 = import_contacts(tgt, chrom, start, end, num_diagonals_to_exclude)
        # fig, axs = plt.subplots(2, 1)
        # axs[0].imshow(m1)
        # axs[1].imshow(m2)
        # plt.show()

        for i in range(m1.shape[0]):
            if args.direction == "vertical":
                # Take a column from the upper triangle
                i1 = min(i + 1, m1.shape[1])
                i0 = max(0, i1 - num_diagonals)
            else:
                # Take a column from the lower triangle (which is equivalent to taking a row from the upper triangle)
                assert args.direction == "horizontal"
                i0 = i
                i1 = min(i + num_diagonals, m1.shape[1])

            score = compute_score(m1[i0:i1, i], m2[i0:i1, i])
            bin_start1 = start + (i * bin_size)
            bin_end1 = min(bin_start1 + bin_size, end)

            bin_start2 = start + (i0 * bin_size)
            bin_end2 = min(start + (i1 * bin_size), end)

            print("\t".join((str(v) for v in (chrom, bin_start1, bin_end1, chrom, bin_start2, bin_end2,
                                              score.oo, score.oi, score.io, score.ii, score.nans))))
            assert sum([score.oo, score.oi, score.io, score.ii + score.nans]) == (bin_end2 - bin_start2) // bin_size
            assert bin_end2 - bin_start2 <= diagonal_width
