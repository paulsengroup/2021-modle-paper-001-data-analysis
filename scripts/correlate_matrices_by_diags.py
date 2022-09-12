#!/usr/bin/env python3

# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

# Important! This script is not meant to be used on large regions!

import argparse

import cooler
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr, kendalltau
from collections import namedtuple
import sys
import warnings


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
    cli.add_argument("--method",
                     type=str,
                     default="pearson",
                     choices=allowed_correlation_methods,
                     help="Correlation method to use. Allowed values are: " + ", ".join(allowed_correlation_methods))
    return cli


def path_to_uri(path, bin_size):
    if path.endswith(".mcool"):
        return path + f"::/resolutions/{bin_size}"
    return path


def fetch_region(cooler, chrom, start, end, balanced=False):
    query = f"{chrom}:{start}-{end}"
    try:
        pixels = np.nan_to_num(cooler.matrix(balance=balanced, as_pixels=False).fetch(query))
    except ValueError:
        pixels = np.nan_to_num(cooler.matrix(balance=False, as_pixels=False).fetch(query))
    return pixels / pixels.max()


def correlate(m1, m2, correlation_fx=pearsonr):
    assert m1.shape[0] == m1.shape[1]
    assert m1.shape == m2.shape

    corr = []
    pval = []
    for i in range(m1.shape[0]):
        v1 = np.diagonal(m1, offset=i)
        v2 = np.diagonal(m2, offset=i)

        assert len(v1) == len(v2)

        if len(v1) >= 2:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                c, pv = correlation_fx(v1, v2)
            corr.append(c)
            pval.append(pv)

    ResultT = namedtuple("ResultT", ["corr", "pvalue"])
    return ResultT(np.nan_to_num(corr), np.nan_to_num(pval))


if __name__ == "__main__":
    allowed_correlation_methods = {"pearson", "spearman", "kendall"}

    args = vars(make_cli().parse_args())

    bin_size = args["bin_size"]
    diagonal_width = args["diagonal_width"]

    ranges = pd.read_table(args["chrom_ranges_bed"], usecols=range(0, 3),
                           names=["chrom", "start", "end"]).drop_duplicates()

    ref = cooler.Cooler(path_to_uri(args["ref_matrix"], bin_size))
    tgt = cooler.Cooler(path_to_uri(args["tgt_matrix"], bin_size))

    assert ref.binsize == bin_size
    assert tgt.binsize == bin_size

    if args["method"] == "pearson":
        def correlation_fx(x, y): return tuple(pearsonr(x, y))
    elif args["method"] == "spearman":
        def correlation_fx(x, y): return tuple(spearmanr(x, y))
    elif args["method"] == "kendall":
        def correlation_fx(x, y): return tuple(kendalltau(x, y))
    else:
        assert args["method"] in allowed_correlation_methods

    df = None
    for chrom, start, end in zip(ranges["chrom"].to_numpy(),
                                 ranges["start"].to_numpy(),
                                 ranges["end"].to_numpy()):
        assert start <= end

        m1 = fetch_region(ref, chrom, start, end, correlation_fx)
        m2 = fetch_region(tgt, chrom, start, end, correlation_fx)

        corrs, pvals = correlate(m1, m2)
        starts = np.arange(len(corrs)) * bin_size

        df1 = pd.DataFrame({"chrom": [chrom] * len(corrs),
                            "start": starts,
                            "end": starts + bin_size,
                            "correlation": corrs,
                            "significance": pvals})

        if df is None:
            df = df1
        else:
            df = pd.concat([df, df1])

    df.to_csv(sys.stdout, index=False, sep="\t")
