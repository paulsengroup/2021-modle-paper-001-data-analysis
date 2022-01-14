#!/usr/bin/env -S python3 -u

import argparse
import bioframe as bf
import glob
import itertools
import numpy as np
import os
import pandas as pd
import pickle
import pyBigWig
import shutil
import subprocess as sp
import sys
import tempfile
import uuid
from multiprocessing import cpu_count
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
    cols = ["chrom", "start", "end"]
    bed = bf.read_table(bed_file,
                        header=0,
                        names=cols,
                        usecols=list(range(3))).drop_duplicates(ignore_index=True)
    bed = bed[~bed["chrom"].isin(excluded_chroms)]

    if chrom_subranges is not None:
        intervals = bf.read_table(chrom_subranges,
                                  header=None,
                                  names=cols,
                                  usecols=list(range(3))).drop_duplicates(ignore_index=True)

        bed = bf.overlap(bed, intervals, how="left").dropna()[cols]

    scores = np.empty(bed.shape[0] * 2, dtype=float)
    with pyBigWig.open(bwig1) as bw1, pyBigWig.open(bwig2) as bw2:
        # Fill scores vector
        for (i, (chrom, start, end)) in bed.iterrows():
            scores[i] = bw1.stats(chrom, int(start), int(end))[0]
            scores[bed.shape[0] + i] = bw2.stats(chrom, int(start), int(end))[0]

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


def run_modle_sim(config, extr_barrier_file):
    tmp_config = f"{tmpdir}/{os.path.basename(config)}"
    shutil.copyfile(config, tmp_config)
    with open(tmp_config, "a") as f:
        f.write(f"output-prefix=\"{modle_out_prefix}\"\n")
        f.write(f"extrusion-barrier-file=\"{extr_barrier_file}\"\n")
    cmd = ["modle", "--config", tmp_config]
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
           "-w", str(diagonal_width),
           "--threads", str(nthreads)]

    run_subprocess(cmd)


def run_modle_tools_eval(input_name, output_name):
    cmd = ["modle_tools", "eval",
           "-i", input_name,
           "--bin-size", str(bin_size),
           "--reference-matrix", reference_matrix,
           "-o", output_name,
           "--metric", modle_tools_eval_metric,
           "-w", str(diagonal_width),
           "--threads", str(nthreads)]

    run_subprocess(cmd)


def fx(params, keep_files=False, print_score=True):
    global epoch

    local_barriers = extr_barriers
    local_barriers["occupancy"] = params

    local_barrier_file = f"{tmpdir}/extr_barriers_{epoch:06}.bed"
    local_barriers.to_csv(local_barrier_file, sep="\t", index=False, header=False)

    run_modle_sim(config, local_barrier_file)
    run_modle_tools_transform(f"{modle_out_prefix}.cool", f"{modle_tools_out_prefix}_transformed.cool",
                              gaussian_blur_sigma_tgt, gaussian_blur_sigma_multiplier_tgt,
                              discretization_thresh_tgt)
    run_modle_tools_eval(
        f"{modle_tools_out_prefix}_transformed.cool", modle_tools_out_prefix)
    n = eval(sites_for_eval,
             glob.glob(f"{modle_tools_out_prefix}_*_horizontal.bw")[0],
             glob.glob(f"{modle_tools_out_prefix}_*_vertical.bw")[0])

    if not keep_files:
        tmpfiles = list(glob.glob(f"{modle_out_prefix}*"))
        tmpfiles.extend(list(glob.glob(f"{modle_tools_out_prefix}*")))
        tmpfiles.append(local_barrier_file)
        for file in set(tmpfiles):
            os.remove(file)

    if print_score:
        assert method is not None
        msg = [epoch]
        msg.extend(params)
        msg = "\t".join([str(x) for x in msg])
        print(f"{msg}\t{n:.8G}")
        with open(f"{out_prefix}_{method}.tsv", "a") as fp:
            print(f"{msg}\t{n:.8G}", file=fp)
        epoch += 1

    return n


def make_cli():
    cli = argparse.ArgumentParser()

    cli.add_argument("--config",
                     help="Path to MoDLE's config file", required=True)
    cli.add_argument("--bin-size",
                     default=5000,
                     type=int)
    cli.add_argument("--diagonal-width",
                     default=int(3e6),
                     type=int)
    cli.add_argument("--chrom-subranges",
                     help="Path to BED file with chrom. subranges")
    cli.add_argument("--discretization-thresh-ref",
                     default=1.5,
                     type=float)
    cli.add_argument("--discretization-thresh-tgt",
                     default=1.0,
                     type=float)
    cli.add_argument("--excluded-chroms",
                     default=["chrY", "chrM"],
                     nargs="*")
    cli.add_argument("--extrusion-barriers",
                     help="Path to extrusion barrier BED file",
                     required=True)
    cli.add_argument("--gaussian-blur-sigma-multiplier-ref",
                     default=1.6,
                     type=float)
    cli.add_argument("--gaussian-blur-sigma-multiplier-tgt",
                     default=1.6,
                     type=float)
    cli.add_argument("--gaussian-blur-sigma-ref",
                     default=2.0,
                     type=float)
    cli.add_argument("--gaussian-blur-sigma-tgt",
                     default=1.0,
                     type=float)
    cli.add_argument("--modle-tools-eval-metric",
                     default="custom",
                     choices=["custom", "eucl_dist", "pearson", "rmse", "spearman"])
    cli.add_argument("--nthreads",
                     default=cpu_count(),
                     type=int)
    cli.add_argument("--output-prefix",
                     help="Output prefix",
                     required=True)
    cli.add_argument("--transformed-reference-matrix",
                     help="Path to the reference matrix after computing the difference of gaussians",
                     required=True)

    cli.add_argument("--ncalls",
                     default=5,
                     type=int)
    cli.add_argument("--nrandom-starts",
                     default=1,
                     type=int)
    cli.add_argument(
        "--optimization-method",
        default="bayesian",
        choices=["dummy", "forest", "gbrt", "bayesian"],
        help="See https://scikit-optimize.github.io/stable/modules/minimize_functions.html")

    cli.add_argument("--seed",
                     default=1630986062,
                     type=int)
    return cli


def import_barriers(path_to_barriers):
    return pd.read_table(path_to_barriers, names=["chrom", "start", "end", "name", "occupancy", "strand"], header=None)


def run_optimize():
    global modle_out_prefix
    global modle_tools_out_prefix
    global epoch

    ncalls = args.ncalls
    nrandom_starts = args.nrandom_starts

    optimizer = None
    if method == "dummy":
        optimizer = dummy_minimize
    elif method == "forest":
        optimizer = forest_minimize
    elif method == "gbrt":
        optimizer = gbrt_minimize
    else:
        assert method == "bayesian"
        optimizer = gp_minimize

    if os.path.dirname(out_prefix) != "":
        os.makedirs(os.path.dirname(out_prefix), exist_ok=True)

    modle_out_prefix = f"{tmpdir}/modle"
    modle_tools_out_prefix = f"{tmpdir}/modle_tools"

    dimensions = [(0.0, 1.0)] * len(extr_barriers)

    seed = args.seed
    res = optimizer(fx,
                    dimensions=dimensions,
                    x0=extr_barriers["occupancy"].tolist(),
                    n_jobs=nthreads,
                    n_calls=ncalls,
                    n_random_starts=nrandom_starts,
                    random_state=seed)

    with open(f"{out_prefix}_{method}.pickle", "wb") as fp:
        pickle.dump(res, fp)


if __name__ == "__main__":
    args = make_cli().parse_args()

    out_prefix = args.output_prefix

    uuid_ = uuid.uuid1()

    config = args.config
    chrom_subranges = args.chrom_subranges
    extr_barriers = args.extrusion_barriers
    reference_matrix = args.transformed_reference_matrix
    excluded_chroms = set(itertools.chain.from_iterable([str(tok).split(",") for tok in args.excluded_chroms]))
    nthreads = args.nthreads
    bin_size = args.bin_size
    diagonal_width = args.diagonal_width

    gaussian_blur_sigma_ref = args.gaussian_blur_sigma_ref
    gaussian_blur_sigma_multiplier_ref = args.gaussian_blur_sigma_multiplier_ref
    discretization_thresh_ref = args.discretization_thresh_ref
    gaussian_blur_sigma_tgt = args.gaussian_blur_sigma_tgt
    gaussian_blur_sigma_multiplier_tgt = args.gaussian_blur_sigma_multiplier_tgt
    discretization_thresh_tgt = args.discretization_thresh_tgt

    modle_tools_eval_metric = args.modle_tools_eval_metric

    sites_for_eval = args.extrusion_barriers

    extr_barriers = import_barriers(extr_barriers)

    epoch = 0

    method = args.optimization_method
    with tempfile.TemporaryDirectory(suffix=f"_optimize_modle_sim_params_{method}") as tmpdir:
        run_optimize()
