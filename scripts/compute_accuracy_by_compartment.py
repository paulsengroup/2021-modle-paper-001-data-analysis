#!/usr/bin/env python3

# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import argparse
import pyBigWig
import numpy as np
import bioframe as bf


def make_cli():
    cli = argparse.ArgumentParser()

    cli.add_argument("stripe_annotation")
    cli.add_argument("bigwig_compartment")
    cli.add_argument("bigwig_horizontal")
    cli.add_argument("bigwig_vertical")
    cli.add_argument("--regions-of-interest",
                     type=str)
    cli.add_argument("--stripe-annotation-format",
                     choices={"stripenn", "bedpe"},
                     default="stripenn")

    return cli


def stripe_is_vertical(s1, e1, s2, e2):
    return ((s1 + e1) / 2) >= ((s2 + e2) / 2)


def import_stripes(stripe_path, regions_of_interest_path, format):
    if format == "stripenn":
        return import_stripes_stripenn(stripe_path, regions_of_interest_path)

    return import_stripes_bedpe(stripe_path, regions_of_interest_path)


def import_stripes_bedpe(stripe_path, regions_of_interest_path):
    df = bf.read_table(stripe_path, schema="bed6")
    if regions_of_interest_path is not None:
        df = bf.overlap(df, bf.read_table(regions_of_interest_path, schema="bed3"))
        df = df.dropna("start_").drop(columns=["chrom_", "start_", "end_"])
    direction = []
    for _, row in df.iterrows():
        if row["strand"] == "-":
            direction.append("vertical")
        elif row["strand"] == "+":
            direction.append("horizontal")
        else:
            direction.append(".")

    df["direction"] = direction
    return df


def import_stripes_stripenn(stripe_path, regions_of_interest_path):
    df = bf.read_table(stripe_path, header=0)
    if regions_of_interest_path is not None:
        df = bf.overlap(df,
                        bf.read_table(regions_of_interest_path, schema="bed3"),
                        cols1=["chr", "pos1", "pos2"])
        df = df.dropna(subset=["start_"]).drop(columns=["chrom_", "start_", "end_"])
    direction = []
    for _, row in df.iterrows():
        if stripe_is_vertical(row["pos1"], row["pos2"], row["pos3"], row["pos4"]):
            direction.append("vertical")
        else:
            direction.append("horizontal")

    df["direction"] = direction
    df = df[["chr", "pos1", "pos2", "direction"]]
    return df.rename(columns={"chr": "chrom",
                              "pos1": "start",
                              "pos2": "end"})


def import_scores_from_bwigs(horizontal_bwig, vertical_bwig, genomic_coords):
    df = genomic_coords.copy()
    scores = []
    with pyBigWig.open(horizontal_bwig) as h_bw, pyBigWig.open(vertical_bwig) as v_bw:
        # Fill scores vector
        for _, row in df.iterrows():
            if row["direction"] == "vertical":
                scores.append(v_bw.stats(row["chrom"],
                                         int(row["start"]),
                                         int(row["end"]),
                                         exact=True)[0])
            elif row["direction"] == "horizontal":
                scores.append(h_bw.stats(row["chrom"],
                                         int(row["start"]),
                                         int(row["end"]),
                                         exact=True)[0])
            else:
                scores.append(np.nan)

    df["score"] = np.array(scores, dtype=float)
    return df


if __name__ == "__main__":
    args = vars(make_cli().parse_args())

    stripes = import_stripes(args["stripe_annotation"],
                             args["regions_of_interest"],
                             args["stripe_annotation_format"])
    stripes = import_scores_from_bwigs(args["bigwig_horizontal"],
                                       args["bigwig_vertical"],
                                       stripes)

    stripes = bf.sort_bedframe(stripes, cols=["chrom", "start", "end"])

    with pyBigWig.open(args["bigwig_compartment"]) as bw:
        # print("chrom\tstart\tend\tname\tscore")
        for _, stripe in stripes.iterrows():
            score = bw.stats(stripe["chrom"], stripe["start"], stripe["end"], exact=True)[0]
            if score is None:
                compartment = "."
            elif score > 0.05:
                compartment = "A"
            elif score < -0.05:
                compartment = "B"
            else:
                compartment = "."

            print("\t".join(str(x) for x in [stripe["chrom"],
                                             stripe["start"],
                                             stripe["end"],
                                             compartment,
                                             stripe["score"]]))
