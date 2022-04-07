#!/usr/bin/env python3

import argparse
import multiprocessing as mp
import subprocess as sp

import cooler
# noinspection PyUnresolvedReferences
import cooltools
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
    cli.add_argument("--output",
                     type=str,
                     required=True,
                     help="Path to output cool.")
    cli.add_argument("--chrom-ranges-bed",
                     type=str,
                     required=True,
                     help="Path to a BED file with the genomic coordinates of the regions to correlate.")
    cli.add_argument("--bin-size",
                     type=int,
                     required=True,
                     help="Matrix bin size.")
    cli.add_argument("--threads",
                     type=int,
                     default=mp.cpu_count(),
                     help="Threads.")
    cli.add_argument("--diagonal-width",
                     type=int,
                     required=True,
                     help="Diagonal width in bp.")
    cli.add_argument("--diagonal-width-to-exclude",
                     type=int,
                     default=50000,
                     help="Width in base-pairs of the band along the diagonal to ignore.")
    return cli


def import_contacts(cooler, chrom, start, end, diagonal_width, diagonal_width_to_exclude):
    coord = f"{chrom}:{start}-{end}"
    pixels = cooler.matrix(as_pixels=True,
                           balance=False,
                           join=True).fetch(coord)

    return pixels[(pixels["end2"] - pixels["end1"] >= diagonal_width_to_exclude) &
                  (pixels["end2"] - pixels["end1"] <= diagonal_width)]


def path_to_uri(path, bin_size):
    if path.endswith(".mcool"):
        return path + f"::/resolutions/{bin_size}"
    return path


def run_cooltools_random_sample(path_to_input, path_to_output, fraction, threads):
    sp.run(["cooltools", "random-sample",
            "--frac", str(fraction),
            "--nproc", str(threads),
            path_to_input,
            path_to_output],
           check=True)


if __name__ == "__main__":
    args = make_cli().parse_args()

    bin_size = args.bin_size
    diagonal_width = args.diagonal_width
    diagonal_width_to_exclude = args.diagonal_width_to_exclude

    ref = cooler.Cooler(path_to_uri(args.ref_matrix, bin_size))
    tgt = cooler.Cooler(path_to_uri(args.tgt_matrix, bin_size))

    assert ref.binsize == bin_size
    assert tgt.binsize == bin_size

    ranges = pd.read_table(args.chrom_ranges_bed,
                           usecols=range(0, 3),
                           names=["chrom", "start", "end"])
    ranges.drop_duplicates(inplace=True)

    ref_contacts = 0
    tgt_contacts = 0
    for chrom, start, end in zip(ranges["chrom"].to_numpy(),
                                 ranges["start"].to_numpy(),
                                 ranges["end"].to_numpy()):
        assert start <= end

        ref_contacts += import_contacts(ref,
                                        chrom, start, end,
                                        diagonal_width,
                                        diagonal_width_to_exclude)["count"].sum()

        tgt_contacts += import_contacts(tgt,
                                        chrom, start, end,
                                        diagonal_width,
                                        diagonal_width_to_exclude)["count"].sum()

        assert tgt_contacts >= ref_contacts

    fraction = ref_contacts / tgt_contacts

    print(f"Subsampling using fraction={fraction:.6G}")
    run_cooltools_random_sample(path_to_uri(args.tgt_matrix, bin_size), args.output, fraction, args.threads)
