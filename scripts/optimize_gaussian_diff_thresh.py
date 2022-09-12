#!/usr/bin/env python3

# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import argparse
import os

import pandas as pd
import bioframe as bf
import cooler
import matplotlib.pyplot as plt
import numpy as np
from scipy import optimize


def make_cli():
    cli = argparse.ArgumentParser()
    cli.add_argument("--ref-cooler",
                     type=str,
                     required=True,
                     help="Path to the reference matrix after gaussian transformation and discretization.")
    cli.add_argument("--tgt-cooler",
                     type=str,
                     required=True,
                     help="Path to the target matrix after gaussian transformation (wo/ discretization).")
    cli.add_argument("--regions-of-interest-bed",
                     type=str,
                     help="Path to a BED file with the genomic regions to consider during the optimization.")
    cli.add_argument("--bin-size",
                     type=int,
                     required=True,
                     help="Contact matrix bin size.")
    cli.add_argument("--output-prefix",
                     required=True,
                     type=str)
    cli.add_argument("--seed",
                     type=int,
                     default=3173852078,
                     help="Random seed.")
    cli.add_argument("--diagonal-width",
                     type=int,
                     required=True,
                     help="Diagonal width in bp.")
    cli.add_argument("--quiet",
                     default=False,
                     action="store_true",
                     help="Don't print status updates to stdout.")
    cli.add_argument("--diagonals-to-exclude",
                     type=int,
                     default=10,
                     help="Number of diagonals to exclude.")
    return cli


def cool_path_to_uri(path, bin_size):
    if path.endswith(".mcool"):
        return path + f"::/resolutions/{bin_size}"
    return path


def import_bed(path_to_bed):
    if path_to_bed is None:
        return None
    return bf.read_table(path_to_bed,
                         usecols=list(range(3)),
                         names=["chrom", "start", "end"]).sort_values(by=["chrom", "start"])


def compute_regions_of_interest(path_to_regions_of_interest, chromsizes1, chromsizes2):
    chrom_sizes = {}
    chrom_names = {name for name, _ in chromsizes1.iteritems()}
    chrom_names = chrom_names.intersection({name for name, _ in chromsizes2.iteritems()})
    assert len(chrom_names) > 0, "Chromosome name intersection is empty!"
    for name in chrom_names:
        assert chromsizes1[name] == chromsizes2[name], f"Chromosome size mismatch for \"{name}\""
        chrom_sizes[name] = chromsizes1[name]

    chrom_sizes = bf.from_dict(chrom_sizes)
    bed = import_bed(path_to_regions_of_interest)
    if bed is None:
        return chrom_sizes

    overlaps = bf.overlap(chrom_sizes, bed, return_overlap=True)[["chrom", "overlap_start", "overlap_end"]].dropna()
    overlaps.columns = ["chrom", "start", "end"]
    return overlaps


# Source: https://stackoverflow.com/a/18081653
def kth_diag_indices(a, k):
    rows, cols = np.diag_indices_from(a)
    if k < 0:
        return rows[-k:], cols[:k]
    elif k > 0:
        return rows[:-k], cols[k:]
    else:
        return rows, cols


def import_contacts(cooler, chrom, start, end, diagonals_to_exclude=0, max_diagonals=np.inf, fill_value=0):
    assert diagonals_to_exclude >= 0
    assert max_diagonals > 0

    coord = f"{chrom}:{start}-{end}"
    data = cooler.matrix(balance=False,
                         as_pixels=False).fetch(coord)

    idx1, idx2 = np.meshgrid(np.arange(data.shape[0]), np.arange(data.shape[1]))
    mask = np.logical_or(np.abs(idx1 - idx2) <= num_diagonals_to_exclude,
                         np.abs(idx1 - idx2) >= max_diagonals)
    data[mask] = fill_value
    return data


def plot_matrices(m1, m2, chrom, start, end, thresh=None):
    fig, axs = plt.subplots(1, 2, figsize=(6.4 * 2, 6.4))
    axs[0].imshow(m1, interpolation="none")
    axs[1].imshow(m2, interpolation="none")

    if thresh is None:
        fig.suptitle(f"{chrom}:{start}-{end}")
    else:
        fig.suptitle(f"{chrom}:{start}-{end} ({thresh})")
    return fig


def obj_fx(threshold, objective, data):
    return abs(objective - np.sum(data >= threshold))


def numpy_pixels_to_dict(chrom, start, matrix, bins):
    bin_size = int(bins.iloc[0]["end"] - bins.iloc[0]["start"])
    pixels = {"bin1_id": [], "bin2_id": [], "count": []}

    bin_offset = bins.index[(bins["chrom"] == chrom) &
                            (bins["start"] >= start) &
                            (bins["end"] <= start + bin_size)]
    assert len(bin_offset) == 1

    bin_offset = int(bin_offset[0])
    for i in range(matrix.shape[0]):
        for j in range(i, matrix.shape[1]):
            if matrix[i][j] != 0:
                pixels["bin1_id"].append(bin_offset + i)
                pixels["bin2_id"].append(bin_offset + j)
                pixels["count"].append(matrix[i][j])

    return pixels


if __name__ == "__main__":
    args = make_cli().parse_args()

    bin_size = args.bin_size
    diagonal_width = args.diagonal_width
    num_diagonals_to_exclude = args.diagonals_to_exclude

    ref_cool = cooler.Cooler(cool_path_to_uri(args.ref_cooler, args.bin_size))
    tgt_cool = cooler.Cooler(cool_path_to_uri(args.tgt_cooler, args.bin_size))

    assert tgt_cool.binsize == ref_cool.binsize, "Bin size mismatch!"

    regions_of_interest = compute_regions_of_interest(args.regions_of_interest_bed,
                                                      ref_cool.chromsizes,
                                                      tgt_cool.chromsizes)
    assert regions_of_interest.shape[0] > 0

    outdir = os.path.dirname(args.output_prefix)
    if outdir not in ["", "."]:
        os.makedirs(outdir, exist_ok=True)

    bins = cooler.util.binnify(ref_cool.chromsizes, bin_size)
    transformed_pixels = {"bin1_id": [], "bin2_id": [], "count": []}
    with open(f"{args.output_prefix}_discretization_thresholds.tsv", "w") as out:
        out.write("chrom\tstart\tend\tthreshold\tabs_difference\n")

        for _, (chrom, start, end) in regions_of_interest.iterrows():
            max_diagonals = int(min(np.ceil(diagonal_width / bin_size),
                                    np.ceil((end - start) / bin_size)))

            m1 = import_contacts(ref_cool, chrom, start, end, num_diagonals_to_exclude, max_diagonals)
            m2 = import_contacts(tgt_cool, chrom, start, end, num_diagonals_to_exclude, max_diagonals)

            fig = plot_matrices(m1, m2, chrom, start, end)
            fig.savefig(f"{args.output_prefix}_{chrom}_{start}_{end}.png", dpi=300)

            if not args.quiet:
                print(f"Processing {chrom}:{start}-{end}...")

            bounds = optimize.Bounds([0.0], [float(np.max(m2.ravel()))])
            objective = np.sum(m1.ravel())


            def obj_fx_converged(xk, convergence):
                return obj_fx(xk, objective, m2) == 0


            result = optimize.differential_evolution(obj_fx,
                                                     args=(objective, m2.ravel()),
                                                     bounds=bounds,
                                                     disp=not args.quiet,
                                                     popsize=100,
                                                     workers=-1,
                                                     atol=0,
                                                     tol=0,
                                                     updating="deferred",
                                                     callback=obj_fx_converged,
                                                     seed=args.seed,
                                                     maxiter=int(100e6))

            solution = result.x[0]
            obj = obj_fx(solution, objective, m2)
            record = "\t".join((str(x) for x in (chrom, start, end, solution, obj)))
            out.write(record + "\n")

            m2[m2 < solution] = 0
            m2[m2 >= solution] = 1
            fig = plot_matrices(m1, m2, chrom, start, end, solution)
            fig.savefig(f"{args.output_prefix}_{chrom}_{start}_{end}_optimized.png", dpi=300)

            pixels = numpy_pixels_to_dict(chrom, start, m2, bins)
            for k in pixels.keys():
                transformed_pixels[k].extend(pixels[k])

        transformed_pixels = pd.DataFrame(transformed_pixels)
        cooler.create_cooler(f"{args.output_prefix}_optimal_threshold.cool", bins, transformed_pixels)
