#!/usr/bin/env python3

# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import numpy as np
import matplotlib.pyplot as plt
import argparse


def make_cli():
    cli = argparse.ArgumentParser()
    cli.add_argument("output_name",
                     type=str)
    cli.add_argument("--t0", type=float, default=0.5)
    cli.add_argument("--t1", type=float, default=0.975)
    cli.add_argument("--exp", type=float, default=10)

    return cli


def penalty_fx(scores, t0, t1, exp):
    penalties = np.zeros_like(scores)
    mask1 = scores < t0
    mask2 = (scores >= t0) & (scores < t1)
    mask3 = scores > t1

    penalties[mask1] = scores[mask1] / t0
    penalties[mask2] = (t1 - scores[mask2]) / (t1 - t0)
    penalties[mask3] = (scores[mask3] - t1) / (1.0 - t1)

    return 1 + (np.power(1 + penalties, exp) - 1) / (2**exp)


if __name__ == "__main__":
    args = vars(make_cli().parse_args())

    x = np.linspace(0.0, 1.0, 1000)
    y = penalty_fx(x,
                   t0=float(args["t0"]),
                   t1=float(args["t1"]),
                   exp=float(args["exp"]))

    fig, ax = plt.subplots(1, 1, figsize=(10, 6.4))

    ax.plot(x, y)
    ax.set(title="Penalty function", xlabel="Occupancy", ylabel="Penalty")

    fig.savefig(args["output_name"])
