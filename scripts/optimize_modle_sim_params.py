#!/usr/bin/env -S python3 -u

import argparse
import glob
import os
import pickle
import subprocess as sp
import sys
import uuid
import tempfile

import numpy as np
import pandas as pd
import pyBigWig
from skopt import dummy_minimize, forest_minimize, gbrt_minimize, gp_minimize
from skopt.space.space import Categorical


def try_convert_to_numeric(x):
    if x.isdigit():
        return int(x)
    try:
        return float(x)
    except Exception:
        return x


def eval(bed_file, bwig1, bwig2):
    bed = pd.read_csv(bed_file,
                      sep="\t",
                      header=1,  # None,
                      names=("chrom", "start", "end"),
                      usecols=list(range(3))).drop_duplicates(ignore_index=True)
    bed = bed[~bed["chrom"].isin(excluded_chroms)]

    scores = np.empty(bed.shape[0] * 2, dtype=float)

    with pyBigWig.open(bwig1) as bw1, pyBigWig.open(bwig2) as bw2:
        # Fill scores vector
        for (i, (chrom, start, end)) in bed.iterrows():
            i *= 2
            scores[i] = bw1.stats(chrom, int(start), int(end))[0]
            scores[i + 1] = bw2.stats(chrom, int(start), int(end))[0]

    # Drop nan and inf values
    scores = scores[(~np.isnan(scores)) & (~np.isinf(scores))]

    if len(scores) == 0:
        return diagonal_width // bin_size

    return np.average(scores)


def run_subprocess(cmd):
    """
    Run an external command as a subprocess, capture its output (stderr as well as stdout)
    and raise an exception if the subprocess has a non-zero return code
    """

    status = sp.run(cmd, encoding="utf-8", stdout=sp.PIPE, stderr=sp.PIPE)
    if status.returncode != 0:
        print(status.stderr, file=sys.stderr)
        print(status.stdout, file=sys.stdout)
        raise RuntimeError(
            f"Subprocess for {cmd[0]} {cmd[1]} exited with code {status.returncode}.\nCMD: " + " ".join(cmd))


def run_modle(params):
    cmd = ["modle", "sim",
           "-c", chrom_sizes,
           "--extrusion-barrier-file", extr_barriers,
           "-o", modle_out_prefix,
           "--randomize-contacts",
           "--threads", str(nthreads),
           "--diagonal-width", str(diagonal_width)]

    assert len(param_df.index) == len(params)
    for param, val in zip(param_df.index, params):
        cmd.extend((param, str(val)))

    run_subprocess(cmd)


def run_modle_tools_transform(input_name, output_name, sigma, sigma_mult, cutoff):
    cmd = ["modle_tools", "transform",
           "-i", input_name,
           "--bin-size", str(bin_size),
           "-m", "difference_of_gaussians",
           "--gaussian-blur-sigma", str(sigma),
           "--gaussian-blur-multiplier", str(sigma_mult),
           "--binary-discretization-value", str(cutoff),
           "-o", output_name,
           "-w", str(diagonal_width)]

    run_subprocess(cmd)


def run_modle_tools_eval(input_name, output_name):
    cmd = ["modle_tools", "eval",
           "-i", input_name,
           "--bin-size", str(bin_size),
           "--reference-matrix", reference_matrix_transformed,
           "-o", output_name,
           "--metric", modle_tools_eval_metric,
           "-w", str(diagonal_width)]

    run_subprocess(cmd)


def fx(params):
    global epoch
    msg = [epoch]
    msg.extend(params)
    msg = "\t".join([str(x) for x in msg])
    print(msg, end="")
    with open(f"{out_prefix}.tsv", "a") as fp:
        print(msg, end="", file=fp)

        run_modle(params)
        run_modle_tools_transform(f"{modle_out_prefix}.cool", f"{modle_tools_out_prefix}_transformed.cool",
                                  gaussian_blur_sigma_tgt, gaussian_blur_sigma_multiplier_tgt,
                                  discretization_thresh_tgt)
        run_modle_tools_eval(
            f"{modle_tools_out_prefix}_transformed.cool", modle_tools_out_prefix)
        n = eval(sites_for_eval,
                 glob.glob(f"{modle_tools_out_prefix}_*_horizontal.bw")[0],
                 glob.glob(f"{modle_tools_out_prefix}_*_vertical.bw")[0])

        tmpfiles = list(glob.glob(f"{modle_out_prefix}*"))
        tmpfiles.extend(list(glob.glob(f"{modle_tools_out_prefix}*")))
        for file in tmpfiles:
            os.remove(file)

        print(f"\t{n:.8G}")
        print(f"\t{n:.8G}", file=fp)

        epoch += 1
        return n


def make_cli():
    cli = argparse.ArgumentParser()
    cli.add_argument("--param-space-tsv", help="Path to a TSV file with three columns names param, start and end.",
                     required=True)
    cli.add_argument("--output-prefix", help="Output prefix", required=True)
    cli.add_argument("--chrom-sizes", help="Path to chrom.sizes file", required=True)
    cli.add_argument("--extrusion-barriers", help="Path to extrusion barrier BED file", required=True)
    cli.add_argument("--reference-matrix", help="Path to multi-res cooler to use as reference", required=True)
    cli.add_argument("--evaluation-sites", help="Path to a BED file with the genomic coordinates to use during evaluation", required=True)

    cli.add_argument("--x0", help="Comma-separated string of starting points")

    cli.add_argument("--excluded-chroms", default=["chrY", "chrM"], nargs="*")
    cli.add_argument("--gaussian-blur-sigma-ref", default=2.0, type=float)
    cli.add_argument("--gaussian-blur-sigma-multiplier_ref",
                     default=1.6, type=float)
    cli.add_argument("--discretization-thresh-ref", default=1.5, type=float)
    cli.add_argument("--gaussian-blur-sigma-tgt", default=1.0, type=float)
    cli.add_argument("--gaussian-blur-sigma-multiplier-tgt",
                     default=1.6, type=float)
    cli.add_argument("--discretization-thresh-tgt", default=1.0, type=float)
    cli.add_argument("--diagonal-width", default=int(3e6), type=int)

    cli.add_argument("--ncalls", default=200, type=int)
    cli.add_argument("--nrandom-starts", default=10, type=int)
    cli.add_argument("--optimization-method", default="bayesian", choices=["dummy", "forest", "gbrt", "bayesian"],
                     help="See https://scikit-optimize.github.io/stable/modules/minimize_functions.html#skopt.forest_minimize")
    cli.add_argument("--seed", default=1630986062, type=int)
    cli.add_argument("--modle-tools-eval-metric", default="custom",
                     choices=["custom", "eucl_dist", "pearson", "rmse", "spearman"])

    return cli


if __name__ == "__main__":
    args = make_cli().parse_args()

    data_dir = args.data_dir
    tmpdir = args.tmp_dir
    out_prefix = args.output_prefix

    uuid_ = uuid.uuid1()

    param_df = pd.read_csv(args.param_space_tsv, sep="\t")
    param_df.set_index("param", inplace=True)
    assert param_df.loc["--bin-size", "start"] == param_df.loc["--bin-size", "end"]

    bin_size = int(param_df.loc["--bin-size", "start"])
    diagonal_width = int(args.diagonal_width)
    chrom_sizes = args.chrom_sizes
    extr_barriers = args.extrusion_barriers
    reference_matrix = args.reference_matrix
    reference_matrix_transformed = reference_matrix.replace(".mcool", "_transformed.cool")
    excluded_chroms = set(args.excluded_chroms)
    nthreads = args.nthreads
    seed = args.seed

    gaussian_blur_sigma_ref = args.gaussian_blur_sigma_ref
    gaussian_blur_sigma_multiplier_ref = args.gaussian_blur_sigma_multiplier_ref
    discretization_thresh_ref = args.discretization_thresh_ref
    gaussian_blur_sigma_tgt = args.gaussian_blur_sigma_tgt
    gaussian_blur_sigma_multiplier_tgt = args.gaussian_blur_sigma_multiplier_tgt
    discretization_thresh_tgt = args.discretization_thresh_tgt

    modle_tools_eval_metric = args.modle_tools_eval_metric

    sites_for_eval = args.eval_sites_bed

    x0 = args.x0
    if x0 is not None:
        x0 = [try_convert_to_numeric(x) for x in x0.split(",")]

    ncalls = args.ncalls
    nrandom_starts = args.nrandom_starts

    optimizer = None
    method = args.optimization_method
    if method == "dummy":
        optimizer = dummy_minimize
    elif method == "forest":
        optimizer = forest_minimize
    elif method == "gbrt":
        optimizer = gbrt_minimize
    else:
        assert method == "bayesian"
        optimizer = gp_minimize

    with tempfile.TemporaryDirectory(suffix=f"optimize_modle_sim_params_{method}") as tmpdir:
        os.makedirs(os.path.dirname(out_prefix), exist_ok=True)
        modle_out_prefix = f"{tmpdir}/modle_"
        modle_tools_out_prefix = f"{tmpdir}/modle_tools_"

        run_modle_tools_transform(reference_matrix, reference_matrix_transformed, gaussian_blur_sigma_ref,
                                  gaussian_blur_sigma_multiplier_ref, discretization_thresh_ref)

        header = "epoch\t" + "\t".join(param_df.index) + "\tscore"
        print(header)
        with open(f"{out_prefix}.tsv", "w") as fp:
            print(header, file=fp)

        dimensions = []
        for param, vals in param_df.iterrows():
            assert len(vals) == 2
            start, end = vals

            try:
                if int(start) == start:
                    start = int(start)
                    end = int(end)
            except Exception:
                pass

            if start == end:
                dimensions.append(Categorical([start]))
            else:
                dimensions.append((start, end))

        assert x0 is None or len(x0) == len(dimensions)

        epoch = 0
        res = optimizer(fx,
                        dimensions=dimensions,
                        x0=x0,
                        n_jobs=nthreads,
                        n_calls=ncalls,
                        n_random_starts=nrandom_starts,
                        random_state=seed)

        with open(f"{out_prefix}_{method}.pickle", "wb") as fp:
            pickle.dump(res, fp)
