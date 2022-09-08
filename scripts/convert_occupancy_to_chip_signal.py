#!/usr/bin/env python3

# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import argparse
import numpy as np
import bioframe as bf
import pyBigWig
import warnings


def make_cli():
    def check_uint(value):
        x = int(value)
        if x <= 0:
            raise argparse.ArgumentTypeError(f"{x} is not a valid positive integer value")
        return x

    cli = argparse.ArgumentParser()
    cli.add_argument("chrom_sizes",
                     type=str)
    cli.add_argument("extrusion_barrier_annotation",
                     type=str)
    cli.add_argument("-o", "--output-prefix", type=str, required=True)
    cli.add_argument("--mu",
                     type=check_uint,
                     default=3,
                     help="Mu parameter, see \"https://www.cell.com/cell-reports/pdf/S2211-1247(16)30530-7.pdf\".")
    cli.add_argument("--bin-signal", default=False, action="store_true")
    cli.add_argument("--bin-size", type=check_uint, default=5000)
    cli.add_argument("--clamp", default=False, action="store_true")
    cli.add_argument("--lb", type=float, default=0.0)
    cli.add_argument("--ub", type=float, default=999.0)
    return cli


def log_fx_inv(y, mu):
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=RuntimeWarning)
        if isinstance(y, (np.ndarray, np.generic)):
            y = np.clip(y, np.nextafter(0.0, 1.0), np.nextafter(1.0, 0.0))
        else:
            y = min(max(np.nextafter(0.0, 1.0), y), np.nextafter(1.0, 0.0))
        return (20 - mu) * (np.log(-y / (y - 1)))


def import_bed_file(path_to_bed, columns=("chrom", "start", "end", "name", "score", "strand")):
    return bf.read_table(path_to_bed, names=columns, usecols=range(0, len(columns)))


def import_chromosomes(path_to_chrom_sizes):
    return bf.read_table(path_to_chrom_sizes, names=["chrom", "length"]).set_index("chrom")["length"]


def bed_to_bigwig(outname, chrom_sizes, bed):
    with pyBigWig.open(outname, "w") as bw:
        bw.addHeader([tuple([k, v]) for k, v in chrom_sizes.to_dict().items()])

        for chrom in chrom_sizes.to_dict().keys():
            df = bed[bed["chrom"] == chrom]

            if len(df) != 0:
                chroms = [str(c) for c in df["chrom"].to_list()]
                starts = [int(s) for s in df["start"].to_list()]
                ends = [int(e) for e in df["end"].to_list()]
                values = [float(x) for x in np.nan_to_num(df["score"].to_list())]
                bw.addEntries(chroms, starts, ends=ends, values=values)


def clamp_scores(scores, lb, ub):
    return np.clip(np.nan_to_num(scores, nan=0.0, neginf=0.0, posinf=ub), lb, ub)


def bin_occupancies(chrom_sizes, occupancies, bin_size, strand=None):
    bins = bf.binnify(chrom_sizes, bin_size)
    df = bf.overlap(bins, occupancies).dropna() \
        .drop(columns=["chrom_", "start_", "end_"])

    df.columns = [k.rstrip("_") for k in df.columns]

    def aggregate_occupancies(df):
        return 1.0 - np.prod(1.0 - df["score"])

    if strand == "both":
        df = df.groupby(by=["chrom", "start", "end"]).apply(aggregate_occupancies).reset_index()
        df["strand"] = "."
    elif strand in {"+", "-"}:
        df = df[df["strand"] == strand].groupby(by=["chrom", "start", "end", "strand"]).apply(
            aggregate_occupancies).reset_index()
    else:
        df = df.groupby(by=["chrom", "start", "end", "strand"]).apply(aggregate_occupancies).reset_index()

    df["score"] = df[0]
    df["name"] = "."
    return df[["chrom", "start", "end", "name", "score", "strand"]]


def make_bed_df(chrom_sizes, barriers, bin_size):
    if bin_size is None:
        chip_df = barriers.copy()
    else:
        chip_df = bin_occupancies(chrom_sizes, barriers, bin_size)

    chip_df["score"] = log_fx_inv(chip_df["score"].to_numpy(), float(args["mu"]))
    if args["clamp"]:
        chip_df["score"] = clamp_scores(chip_df["score"], args["lb"], args["ub"])

    return chip_df


if __name__ == "__main__":
    args = vars(make_cli().parse_args())

    chrom_sizes = import_chromosomes(args["chrom_sizes"])
    extrusion_barriers = import_bed_file(args["extrusion_barrier_annotation"]).dropna()

    out_prefix = args["output_prefix"]

    df = make_bed_df(chrom_sizes, extrusion_barriers, int(args["bin_size"]) if args["bin_signal"] else None)
    df.to_csv(f"{out_prefix}.bed.gz", sep="\t", header=False, index=False, float_format="%.4g", na_rep="nan")

    strand = "both"

    bin_size = int(args["bin_size"])
    df = bin_occupancies(chrom_sizes, extrusion_barriers, bin_size, strand)
    df["score"] = log_fx_inv(df["score"].to_numpy(), float(args["mu"]))
    if args["clamp"]:
        df["score"] = clamp_scores(df["score"], args["lb"], args["ub"])

    df.to_csv(f"{out_prefix}_{bin_size}.bed.gz", sep="\t", header=False, index=False, float_format="%.4g", na_rep="nan")

    for strand in {"-", "+", "both"}:
        df = bin_occupancies(chrom_sizes, extrusion_barriers, bin_size, strand)
        df["score"] = log_fx_inv(df["score"].to_numpy(), float(args["mu"]))
        if args["clamp"]:
            df["score"] = clamp_scores(df["score"], args["lb"], args["ub"])

        if strand == "-":
            strand = "rev"
        elif strand == "+":
            strand = "fwd"

        bed_to_bigwig(f"{out_prefix}_{strand}_{bin_size}.bw", chrom_sizes, df)
