#!/usr/bin/env -S python3 -u

import argparse
import glob
import itertools
import os
import pickle
import subprocess as sp
import sys
import tempfile
import uuid
from multiprocessing import cpu_count

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


def run_modle_sim(params):
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
           "--reference-matrix", reference_matrix,
           "-o", output_name,
           "--metric", modle_tools_eval_metric,
           "-w", str(diagonal_width)]

    run_subprocess(cmd)


def fx(params, keep_files=False, print_score=True):
    if print_score:
        global epoch
        msg = [epoch]
        msg.extend(params)
        msg = "\t".join([str(x) for x in msg])
        print(msg, end="")
        with open(f"{out_prefix}_{method}.tsv", "a") as fp:
            print(msg, end="", file=fp)

    run_modle_sim(params)
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
        for file in set(tmpfiles):
            os.remove(file)

    if print_score:
        print(f"\t{n:.8G}")
        with open(f"{out_prefix}_{method}.tsv", "a") as fp:
            print(f"\t{n:.8G}", file=fp)
        epoch += 1

    return n


def make_cli():
    cli = argparse.ArgumentParser()

    parser = cli.add_subparsers(dest="cmd", required=True)
    test = parser.add_parser("test",
                             help=f"Test the set of \"optimal\" parameters identified by {os.path.basename(sys.argv[0])} optimize")
    optimize = parser.add_parser("optimize",
                                 help="Run the optimization procedure to find a set of \"optimal\" parameters for modle sim")

    def check_positive_arg(value):
        if int(value) <= 0:
            raise argparse.ArgumentTypeError(f"{value} is not a positive integer value")
        return int(value)

    def define_common_arguments(p):
        p.add_argument("--chrom-sizes",
                       help="Path to chrom.sizes file", required=True)
        p.add_argument("--diagonal-width",
                       default=int(3e6),
                       type=int)
        p.add_argument("--discretization-thresh-ref",
                       default=1.5,
                       type=float)
        p.add_argument("--discretization-thresh-tgt",
                       default=1.0,
                       type=float)
        p.add_argument("--evaluation-sites",
                       help="Path to a BED file with the genomic coordinates to use during evaluation",
                       required=True)
        p.add_argument("--excluded-chroms",
                       default=["chrY", "chrM"],
                       nargs="*")
        p.add_argument("--extrusion-barriers",
                       help="Path to extrusion barrier BED file",
                       required=True)
        p.add_argument("--gaussian-blur-sigma-multiplier-ref",
                       default=1.6,
                       type=float)
        p.add_argument("--gaussian-blur-sigma-multiplier-tgt",
                       default=1.6,
                       type=float)
        p.add_argument("--gaussian-blur-sigma-ref",
                       default=2.0,
                       type=float)
        p.add_argument("--gaussian-blur-sigma-tgt",
                       default=1.0,
                       type=float)
        p.add_argument("--modle-tools-eval-metric",
                       default="custom",
                       choices=["custom", "eucl_dist", "pearson", "rmse", "spearman"])
        p.add_argument("--nthreads",
                       default=cpu_count(),
                       type=int)
        p.add_argument("--output-prefix",
                       help="Output prefix",
                       required=True)
        p.add_argument("--transformed-reference-matrix",
                       help="Path to the reference matrix after computing the difference of gaussians",
                       required=True)

    define_common_arguments(optimize)
    define_common_arguments(test)

    optimize.add_argument("--ncalls",
                          default=5,
                          type=int)
    optimize.add_argument("--nrandom-starts",
                          default=1,
                          type=int)
    optimize.add_argument(
        "--optimization-method",
        default="bayesian",
        choices=["dummy", "forest", "gbrt", "bayesian"],
        help="See https://scikit-optimize.github.io/stable/modules/minimize_functions.html")

    optimize.add_argument("--param-space-tsv",
                          help="Path to a TSV file with three columns names param, start and end.",
                          required=True)
    optimize.add_argument("--seed",
                          default=1630986062,
                          type=int)
    optimize.add_argument("--x0",
                          help="Comma-separated string of starting points")

    test.add_argument("--optimization-result-tsv",
                      help=f"TSV containing the result of the optimization conducted by {os.path.basename(sys.argv[0])} optimize",
                      required=True,
                      type=str)
    test.add_argument(
        "--num-params-to-test",
        help="Number of parameter combinations to test. Parameters are first sorted in ascending order based on their score, then the top N parameter combinations are tested",
        default=1,
        type=check_positive_arg)

    return cli


def run_testing():
    global modle_out_prefix
    global modle_tools_out_prefix
    global param_df
    global bin_size

    result_df = pd.read_csv(args.optimization_result_tsv, sep="\t").sort_values("score", ascending=True)
    params = list(result_df.columns)
    params.remove("epoch")
    params.remove("score")

    param_df = pd.DataFrame(index=params)  # This df is used by run_modle_sim

    assert (result_df["--bin-size"] == result_df["--bin-size"][0]).all(0)
    bin_size = int(result_df.loc[0, "--bin-size"])

    with open(f"{out_prefix}_test_result.tsv", "w") as fp:
        header = "\t".join(param_df.index) + "\tscore"
        print(f"{header}")
        print(f"{header}", file=fp)

        for i in range(args.num_params_to_test):
            modle_out_prefix = f"{i:03d}_{out_prefix}_modle"
            modle_tools_out_prefix = f"{i:03d}_{out_prefix}_modle_tools"

            params = result_df.values.tolist()[i][1:-1]
            score = fx(params, keep_files=True, print_score=False)

            msg = "\t".join([str(x) for x in params])
            print(f"{msg}\t{score}")
            print(f"{msg}\t{score}", file=fp)


def run_optimize():
    global modle_out_prefix
    global modle_tools_out_prefix
    global epoch
    global param_df
    global bin_size

    param_df = pd.read_csv(args.param_space_tsv, sep="\t", dtype=str)
    param_df.set_index("param", inplace=True)
    assert param_df.loc["--bin-size", "start"] == param_df.loc["--bin-size", "end"]
    bin_size = int(param_df.loc["--bin-size", "start"])

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

    with tempfile.TemporaryDirectory(suffix=f"_optimize_modle_sim_params_{method}") as tmpdir:
        if os.path.dirname(out_prefix) != "":
            os.makedirs(os.path.dirname(out_prefix), exist_ok=True)

        modle_out_prefix = f"{tmpdir}/modle"
        modle_tools_out_prefix = f"{tmpdir}/modle_tools"

        header = "epoch\t" + "\t".join(param_df.index) + "\tscore"
        print(header)
        with open(f"{out_prefix}_{method}.tsv", "w") as fp:
            print(header, file=fp)

        dimensions = []
        for i, (param, vals) in enumerate(param_df.iterrows()):
            assert len(vals) == 2
            start, end = (try_convert_to_numeric(vals[0]), try_convert_to_numeric(vals[1]))
            assert isinstance(start, type(end))

            if start == end:
                dimensions.append(Categorical([start]))
                assert dimensions[i].rvs(1)[0] == x0[i]
            else:
                dimensions.append((start, end))

        assert x0 is None or len(x0) == len(dimensions)

        seed = args.seed
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


if __name__ == "__main__":
    args = make_cli().parse_args()

    out_prefix = args.output_prefix

    uuid_ = uuid.uuid1()

    diagonal_width = int(args.diagonal_width)
    chrom_sizes = args.chrom_sizes
    extr_barriers = args.extrusion_barriers
    reference_matrix = args.transformed_reference_matrix
    excluded_chroms = set(itertools.chain.from_iterable([str(tok).split(",") for tok in args.excluded_chroms]))
    nthreads = args.nthreads

    gaussian_blur_sigma_ref = args.gaussian_blur_sigma_ref
    gaussian_blur_sigma_multiplier_ref = args.gaussian_blur_sigma_multiplier_ref
    discretization_thresh_ref = args.discretization_thresh_ref
    gaussian_blur_sigma_tgt = args.gaussian_blur_sigma_tgt
    gaussian_blur_sigma_multiplier_tgt = args.gaussian_blur_sigma_multiplier_tgt
    discretization_thresh_tgt = args.discretization_thresh_tgt

    modle_tools_eval_metric = args.modle_tools_eval_metric

    sites_for_eval = args.evaluation_sites

    if args.cmd == "optimize":
        run_optimize()
    else:
        assert args.cmd == "test"
        run_testing()
