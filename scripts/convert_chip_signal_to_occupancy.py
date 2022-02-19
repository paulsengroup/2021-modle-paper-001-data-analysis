#!/usr/bin/env python3

import argparse
from sys import stdout

import bioframe as bf
import cooler
import numpy as np
import pandas as pd
import pyBigWig
from natsort import natsort_keygen


def make_cli():
    def check_uint(value):
        x = int(value)
        if x <= 0:
            raise argparse.ArgumentTypeError(f"{x} is not a valid positive integer value")
        return x

    def check_is_valid_score(value):
        x = float(value)
        if not 0 <= x <= 1:
            raise argparse.ArgumentTypeError(f"{x} is not a number between 0.0 and 1.0.")
        return x

    cli = argparse.ArgumentParser()
    cli.add_argument("--chip-seq-bigwig",
                     type=str,
                     required=True,
                     help="Path to a bigwig file with CTCF fold-change over control.")
    cli.add_argument("--motifs-bed",
                     type=str,
                     required=True,
                     help="BED file with the genomic coordinates of a motif of interest.")
    cli.add_argument("--regions-of-interest-bed",
                     type=str,
                     nargs="+",
                     help="One or more BED files representing a set of regions of interests."
                          "The results returned by this script correspond to the intersection between --motifs-bed and the BED files specified though this option.")
    cli.add_argument("--bin-size-bp",
                     type=check_uint,
                     default=1000,
                     help="Bin size used to compute extr. barrier occupancies.")
    cli.add_argument("--mu",
                     type=check_uint,
                     default=3,
                     help="Mu parameter, see \"https://www.cell.com/cell-reports/pdf/S2211-1247(16)30530-7.pdf\".")
    cli.add_argument("--occupancy-lower-bound",
                     type=check_is_valid_score,
                     default=0.6,
                     help="Lower bound for the occupancy. Bins with an occupancy smaller than --occupancy-lower-bound will be discarded.")
    cli.add_argument("--occupancy-upper-bound",
                     type=check_is_valid_score,
                     default=1.0,
                     help="Upper bound for the occupancy. Bins with an occupancy larger than --occupancy-upper-bound will be discarded.")
    return cli


def log_fx(x, mu):
    return 1 / (1 + np.exp(-x / (20 - mu)))


def import_bed_file(path_to_bed, columns=("chrom", "start", "end")):
    return bf.read_table(path_to_bed, names=columns, usecols=range(0, len(columns)))


def import_motifs(path_to_bed):
    return import_bed_file(path_to_bed, columns=["chrom", "start", "end", "name", "score", "strand"])


def fetch_peak_summits(motifs, path_to_bwig):
    with pyBigWig.open(path_to_bwig) as bwfp:
        avail_chroms = bwfp.chroms()

        scores = np.full(motifs.shape[0], np.nan, dtype=float)
        for i, (_, row) in enumerate(motifs.iterrows()):
            name = row["chrom"]
            if name in avail_chroms:
                start = row["start"]
                end = row["end"]
                assert start < end
                assert start <= avail_chroms[name]
                assert end <= avail_chroms[name]

                scores[i] = bwfp.stats(name, start, end, type="max")[0]

    return scores


def import_chromosomes(path_to_bwig):
    with pyBigWig.open(path_to_bwig) as bwfp:
        chroms = bwfp.chroms()

    return pd.DataFrame({"chrom": chroms.keys(), "length": chroms.values()}).set_index("chrom")["length"]


def compute_occupancy(motifs, bins):
    occupancies = bf.overlap(bins, motifs, how="inner").drop(["start_", "end_", "name_"], axis="columns")

    occupancies = occupancies.groupby(by=["chrom", "start", "end"], observed=True).sum()
    occupancies.columns = [label.rstrip("_") for label in occupancies.columns]
    occupancies["score"] = log_fx(occupancies["score"], args.mu)

    return occupancies.reset_index()


if __name__ == "__main__":
    cli = make_cli()
    args = cli.parse_args()
    bin_size = args.bin_size_bp
    occupancy_lower_bound = args.occupancy_lower_bound
    occupancy_upper_bound = args.occupancy_upper_bound
    assert bin_size > 0
    assert occupancy_lower_bound <= occupancy_upper_bound, "--occupancy-lower-bound should be smaller than --occupancy-upper-bound"

    motifs = import_motifs(args.motifs_bed)
    if args.regions_of_interest_bed is not None:
        for bed in args.regions_of_interest_bed:
            motifs = bf.overlap(motifs, import_bed_file(bed), how="inner").drop(["chrom_", "start_", "end_"],
                                                                                axis="columns")

    motifs["score"] = fetch_peak_summits(motifs, args.chip_seq_bigwig)

    chroms = import_chromosomes(args.chip_seq_bigwig)
    bins = bf.from_any(cooler.util.binnify(chroms, bin_size))

    fwd_occupancies = compute_occupancy(motifs[motifs["strand"] == "+"], bins)
    rev_occupancies = compute_occupancy(motifs[motifs["strand"] == "-"], bins)

    # This workaround is required for the output of the script to be acceptetd by MoDLE
    rev_occupancies["start"] += 1
    rev_occupancies["end"] += 1

    fwd_occupancies["strand"] = "+"
    rev_occupancies["strand"] = "-"

    occ = pd.concat([rev_occupancies, fwd_occupancies]).sort_values(by=["chrom", "start", "strand"],
                                                                    key=natsort_keygen())

    occ = occ[occ["score"].between(occupancy_lower_bound, occupancy_upper_bound, inclusive="both")]

    if "name" not in occ:
        occ["name"] = "."

    cols = ["chrom", "start", "end", "name", "score", "strand"]
    occ[cols].to_csv(stdout, sep="\t", header=None, index=False, float_format="%.4g")
