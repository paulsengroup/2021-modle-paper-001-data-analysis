#!/usr/bin/env python3

# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import argparse
import numpy as np
import pandas as pd
import sys


def make_cli():
    cli = argparse.ArgumentParser()

    cli.add_argument("barrier_annotation_bed", type=str)
    cli.add_argument("--puu", type=float, default=0.7)

    cli.add_argument("--remove-weak-barriers",
                     action="store_true",
                     default=False)
    cli.add_argument("--weak-barrier-occ-threshold",
                     type=float,
                     default=0.5)
    return cli


def import_barriers(path_to_barriers):
    return pd.read_table(path_to_barriers, names=["chrom", "start", "end", "name", "score", "strand"],
                         index_col=False)


def compute_pbb(puus, occs):
    tp_inactive_to_active = 1.0 - puus
    tp_active_to_inactive = (tp_inactive_to_active - (occs * tp_inactive_to_active)) / occs
    pbbs = np.clip(1.0 - tp_active_to_inactive, 0.0, 1.0)

    mask = (puus != 1.0) & (occs != 0.0)
    pbbs[~mask] = 0.0
    return pbbs


def compute_occs(puus, pbbs):
    tp_inactive_to_active = 1.0 - puus
    tp_active_to_inactive = 1.0 - pbbs

    occs = np.clip(tp_inactive_to_active / (tp_inactive_to_active + tp_active_to_inactive), 0.0, 1.0)

    mask = (puus != 1.0) & (pbbs != 0.0)
    occs[~mask] = 0.0

    return occs


def normalize_occupancies(annotation, puu):
    pbbs = compute_pbb(annotation["name"], annotation["score"])
    annotation["name"] = puu
    annotation["score"] = compute_occs(puu, pbbs)

    return annotation


if __name__ == "__main__":
    args = vars(make_cli().parse_args())

    assert 0.0 <= float(args["puu"]) <= 1.0
    assert 0.0 <= float(args["weak_barrier_occ_threshold"]) <= 1.0

    df = import_barriers(str(args["barrier_annotation_bed"]))
    df = normalize_occupancies(df, float(args["puu"]))

    if args["remove_weak_barriers"]:
        df = df[df["score"] >= float(args["weak_barrier_occ_threshold"])]

    df.to_csv(sys.stdout, sep="\t", header=None, index=False)
