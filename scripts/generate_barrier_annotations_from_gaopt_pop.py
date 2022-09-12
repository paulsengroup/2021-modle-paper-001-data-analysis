#!/usr/bin/env python3

# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import argparse
import numpy as np
import pandas as pd
from deap import base, creator
from collections import namedtuple
import cloudpickle
import os


def make_cli():
    cli = argparse.ArgumentParser()

    cli.add_argument("barrier_annotation_bed", type=str)
    cli.add_argument("pickled_pop", type=str)
    cli.add_argument("-o", "--output-prefix", type=str, required=True)
    cli.add_argument("--remove-weak-barriers",
                     action="store_true",
                     default=False)
    cli.add_argument("--weak-barrier-occ-threshold",
                     type=float,
                     default=0.5)
    cli.add_argument("--weak-barrier-puu-threshold",
                     type=float,
                     default=1.0)

    return cli


def import_population(path_to_pickled_pop):
    with open(path_to_pickled_pop, "rb") as fp:
        # Workaround for pickle not working reliably
        old_pop = cloudpickle.load(fp)
        new_pop = []
        for individual in old_pop:
            new_pop.append(creator.Individual([BarrierT(puu, occ) for puu, occ in individual]))
            new_pop[-1].fitness.values = individual.fitness.values
        return new_pop


def import_barriers(path_to_barriers):
    return pd.read_table(path_to_barriers, names=["chrom", "start", "end", "name", "score", "strand"],
                         index_col=False)


def generate_barrier_annotation(barrier_annotation, barrier_params):
    assert len(barrier_annotation) == len(barrier_params), f"{len(barrier_annotation)} != {len(barrier_params)}"
    df1 = barrier_annotation.copy()
    df1["name"] = np.array([puu for puu, _ in ind], dtype=float)
    df1["score"] = np.array([occ for _, occ in ind], dtype=float)

    return df1


if __name__ == "__main__":
    args = vars(make_cli().parse_args())
    assert 0.0 <= float(args["weak_barrier_occ_threshold"]) <= 1.0
    assert 0.0 <= float(args["weak_barrier_puu_threshold"]) <= 1.0

    BarrierT = namedtuple("BarrierT", ["puu", "occ"])
    creator.create("FitnessMulti", base.Fitness, weights=(-1.0,))
    creator.create("Individual", list, fitness=creator.FitnessMulti)

    pop = import_population(str(args["pickled_pop"]))
    base_annotation = import_barriers(str(args["barrier_annotation_bed"]))
    out_prefix = str(args["output_prefix"])

    if os.path.dirname(out_prefix) != "." and os.path.dirname(out_prefix) != "":
        os.makedirs(os.path.dirname(out_prefix), exist_ok=True)

    padding = len(str(len(pop)))
    for i, ind in enumerate(pop):
        df = generate_barrier_annotation(base_annotation, ind)
        if args["remove_weak_barriers"]:
            df = df[(df["score"] >= float(args["weak_barrier_occ_threshold"])) &
                    (df["name"] < float(args["weak_barrier_puu_threshold"]))]

        df.to_csv(f"{out_prefix}_{i:0{padding}}.bed.gz", sep="\t", header=None, index=False)
