#!/usr/bin/env -S python3 -u

# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import argparse
import glob
import itertools
import pathlib
import shutil
import subprocess as sp
import sys
import tempfile
from collections import namedtuple
from multiprocessing import cpu_count

import bioframe as bf
import cloudpickle
import numpy as np
import pandas as pd
import pyBigWig
from skopt import dummy_minimize, forest_minimize, gbrt_minimize, gp_minimize
from skopt.space.space import Categorical


def try_convert_to_numeric(s):
    """
    Attempt to convert string s to a numeric type.
    Return the string itself if no conversion is possible
    """
    if s.isdigit():
        return int(s)
    try:
        return float(s)
    except Exception:
        return s


def tee(message, file1, quiet=False, mode="a"):
    """
    Append message to file1 and also print it to stdout if quiet is False
    """
    if isinstance(file1, (str, pathlib.Path)):
        with open(file1, mode) as fp:
            print(message, file=fp)
    else:
        print(message, file=file1)

    if not quiet:
        print(message, file=sys.stdout)


def generate_optimization_dimensions(params, x0):
    """
    Convert the params dataframe into a list of dimensions that can be readily used by skopt
    """
    dimensions = []
    for i, (_, vals) in enumerate(params.iterrows()):
        assert len(vals) == 2
        start = try_convert_to_numeric(vals[0])
        end = try_convert_to_numeric(vals[1])
        assert isinstance(start, type(end))

        if start == end:
            dimensions.append(Categorical([start]))
            if x0 is not None:
                # If dimension is categorical with cardinality == 1,
                # sampling a random variate using rvs should always
                # return the same value
                assert dimensions[i].rvs(1)[0] == x0[i]
        else:
            dimensions.append((start, end))

    if x0 is not None:
        assert len(x0) == len(dimensions)

    return tuple(dimensions)


def extract_bin_size_from_params(param_df, key="--resolution"):
    return int(param_df.loc[key, "start"])


def import_evaluation_sites(eval_sites_bed, excluded_chroms=None, sites_of_interest_bed=None):
    if excluded_chroms is None:
        excluded_chroms = set()

    cols = ["chrom1", "start1", "end1",
            "chrom2", "start2", "end2"]
    ncols = len(cols)

    # Import eval sites
    bed = bf.read_table(eval_sites_bed,
                        header=0,
                        names=cols,
                        usecols=list(range(ncols))).drop_duplicates(ignore_index=True)

    # Drop excluded chroms
    if excluded_chroms is not None:
        bed = bed[~bed["chrom1"].isin(excluded_chroms)]

    # Filter BED records using the sites of interests
    if sites_of_interest_bed is not None:
        intervals = bf.read_table(sites_of_interest_bed,
                                  header=None,
                                  names=cols[:3],
                                  usecols=list(range(3))).drop_duplicates(ignore_index=True)

        bed = bf.overlap(bed, intervals, how="left").dropna()[cols]

    return bed


def import_scores_from_bwigs(horizontal_bwig, vertical_bwig, genomic_coords):
    def stripe_is_vertical(s1, e1, s2, e2):
        return ((s1 + e1) / 2) >= ((s2 + e2) / 2)

    scores = []
    with pyBigWig.open(horizontal_bwig) as h_bw, pyBigWig.open(vertical_bwig) as v_bw:
        # Fill scores vector
        for (_, (chrom, start1, end1, _, start2, end2)) in genomic_coords.iterrows():
            if stripe_is_vertical(start1, end1, start2, end2):
                scores.append(v_bw.stats(chrom, int(start1), int(end1))[0])
            else:
                scores.append(h_bw.stats(chrom, int(start1), int(end1))[0])

    # Drop nan and inf values
    scores = np.array(scores, dtype=float)
    return scores[(~np.isnan(scores)) & (~np.isinf(scores))]


def objective_fx(eval_sites, horizontal_bwig, vertical_bwig, bin_size, **kwargs):
    scores = import_scores_from_bwigs(horizontal_bwig, vertical_bwig, eval_sites)

    # Deal with rows/columns of nans
    if len(scores) == 0:
        return int(args.diagonal_width) // bin_size

    return np.average(scores)


def run_subprocess(cmd, timeout=900, attempts=3):
    try:
        sp.check_output(cmd,
                        encoding="utf-8",
                        stdin=sp.DEVNULL,
                        stderr=sp.STDOUT,
                        timeout=timeout)
    except sp.TimeoutExpired:
        if attempts == 0:
            raise
        run_subprocess(cmd, timeout, attempts - 1)
    except sp.CalledProcessError as e:
        cmd_ = " ".join(cmd)
        msg = f"Subprocess for {cmd[0]} {cmd[1]} exited with code {e.returncode}.\nCMD: {cmd_}\n\n{e.output}"
        raise RuntimeError(msg)


def run_modle_sim(params, **kwargs):
    cmd = ["modle", "sim",
           "-c", str(kwargs.get("chrom_sizes")),
           "--extrusion-barrier-file", str(kwargs.get("extrusion_barriers")),
           "--force",
           "-o", str(kwargs.get("tmp_output_prefix")),
           "--threads", str(kwargs.get("threads")),
           "--diagonal-width", str(kwargs.get("diagonal_width"))]

    chrom_subranges = kwargs.get("chrom_subranges")
    if chrom_subranges is not None:
        cmd.extend(["--chrom-subranges", chrom_subranges])

    # Add optional arguments
    cli_options = kwargs.get("params_df").index
    assert len(cli_options) == len(params)
    cmd.extend(itertools.chain.from_iterable(([param, str(val)] for param, val in zip(cli_options, params))))
    run_subprocess(cmd)

    ModleSimOutputs = namedtuple("ModleSimOutputs", ["cooler", "log", "config"])
    out_prefix = kwargs.get("tmp_output_prefix")

    return ModleSimOutputs(out_prefix.with_suffix(".cool"),
                           out_prefix.with_suffix(".log"),
                           pathlib.Path(f"{out_prefix}_config.toml"))


def run_modle_tools_transform(input_name, sigma, sigma_multiplier, cutoff, **kwargs):
    output_name = input_name.with_suffix("")
    output_name = f"{output_name}_transformed.cool"
    cmd = ["modle_tools", "transform",
           "-i", str(input_name),
           "--resolution", str(kwargs.get("bin_size")),
           "-m", "difference_of_gaussians",
           "--gaussian-blur-sigma", str(sigma),
           "--gaussian-blur-multiplier", str(sigma_multiplier),
           "--force",
           "--binary-discretization-value", str(cutoff),
           "-o", str(output_name),
           "-w", str(kwargs.get("diagonal_width")),
           "--threads", str(kwargs.get("threads"))]

    run_subprocess(cmd)

    ModleToolsTransformOutputs = namedtuple("ModleToolsTransformOutputs", ["cooler"])
    return ModleToolsTransformOutputs(output_name)


def run_modle_tools_eval(input_name, **kwargs):
    cmd = ["modle_tools", "eval",
           "-i", str(input_name),
           "--resolution", str(kwargs.get("bin_size")),
           "--reference-matrix", str(kwargs.get("transformed_reference_matrix")),
           "-o", str(kwargs.get("tmp_output_prefix")),
           "--metric", kwargs.get("modle_tools_eval_metric"),
           "-w", str(kwargs.get("diagonal_width")),
           "--threads", str(kwargs.get("threads")),
           "--force"]

    run_subprocess(cmd)

    ModleToolsEvalOutputs = namedtuple("ModleToolsEvalOutputs",
                                       ["bigwig_horizontal", "bigwig_vertical", "tsv_horizontal", "tsv_vertical"])

    out_prefix = kwargs.get("tmp_output_prefix")
    bwig_horizontal = glob.glob(f"{out_prefix}_*_horizontal.bw")
    bwig_vertical = glob.glob(f"{out_prefix}_*_vertical.bw")
    tsv_horizontal = glob.glob(f"{out_prefix}_*_horizontal.tsv*")
    tsv_vertical = glob.glob(f"{out_prefix}_*_vertical.tsv*")

    assert len(bwig_horizontal) == 1
    assert len(bwig_vertical) == 1
    assert len(tsv_horizontal) == 1
    assert len(tsv_vertical) == 1

    return ModleToolsEvalOutputs(bwig_horizontal[0], bwig_vertical[0],
                                 tsv_horizontal[0], tsv_vertical[0])


def fx(params, **kwargs):
    global epoch

    kwargs["bin_size"] = extract_bin_size_from_params(kwargs.get("params_df"))

    padding_length = len(str(kwargs["ncalls"]))
    out_prefix = kwargs.get("output_prefix")
    out_name = f"{str(epoch).zfill(padding_length)}_{out_prefix.stem}"
    method = kwargs.get("optimization_method")

    with tempfile.TemporaryDirectory(suffix=f"_modle_sim_param_optimization_{method}") as tmpdir:
        tmp_output_dir = f"{tmpdir}/{out_name}"
        kwargs["tmp_output_prefix"] = pathlib.Path(f"{tmp_output_dir}/{out_name}")

        modle_cooler = run_modle_sim(params, **kwargs).cooler
        modle_cooler_transformed = run_modle_tools_transform(modle_cooler,
                                                             kwargs.get("gaussian_blur_sigma_tgt"),
                                                             kwargs.get("gaussian_blur_sigma_multiplier_tgt"),
                                                             kwargs.get("discretization_thresh_tgt"),
                                                             **kwargs).cooler
        bwig_horizontal, bwig_vertical, _, _ = run_modle_tools_eval(modle_cooler_transformed, **kwargs)

        training_score = objective_fx(kwargs.get("sites_for_eval_training"),
                                      bwig_horizontal,
                                      bwig_vertical,
                                      **kwargs)

        validation_score = objective_fx(kwargs.get("sites_for_eval_validation"),
                                        bwig_horizontal,
                                        bwig_vertical,
                                        **kwargs)

        # Save output files in a tar archive
        shutil.make_archive(f"{out_prefix}_{out_name}",
                            "tar",
                            root_dir=tmpdir,
                            base_dir=out_name)

        # Log fx evaluation to report file
        path_to_report = f"{out_prefix}_{method}.tsv"
        toks = [str(tok) for tok in [epoch] + params + [f"{training_score:.8G}", f"{validation_score:.8G}"]]
        tee("\t".join(toks), path_to_report)

        epoch += 1
        return training_score


def run_optimize(**kwargs):
    # Read param space from file
    param_df = pd.read_table(kwargs.get("param_space_tsv"), dtype=str)
    param_df.set_index("param", inplace=True)
    kwargs["params_df"] = param_df

    # Read bin size from param space
    kwargs["bin_size"] = extract_bin_size_from_params(kwargs.get("params_df"))

    # Process starting search point when appropriate
    x0 = kwargs.get("x0")
    if x0 is not None:
        x0 = [try_convert_to_numeric(x) for x in x0.split(",")]

    # Select the appropriate optimizer
    method = kwargs.get("optimization_method")
    if method == "dummy":
        optimizer = dummy_minimize
    elif method == "forest":
        optimizer = forest_minimize
    elif method == "gbrt":
        optimizer = gbrt_minimize
    else:
        assert method == "bayesian"
        optimizer = gp_minimize

    # Create output folder
    out_prefix = kwargs.get("output_prefix")
    out_prefix.mkdir(parents=True, exist_ok=True)

    # Write report header to file
    optimization_report_tsv = pathlib.Path(f"{out_prefix}_{method}.tsv")
    optimization_report_tsv.unlink(missing_ok=True)

    header = "epoch\t" + "\t".join(param_df.index) + "\ttraining_score\tvalidation_score"
    tee(header, optimization_report_tsv)

    # Prepare dimensions to use for parameter search
    dimensions = generate_optimization_dimensions(param_df, x0)

    # Run the optimization
    res = optimizer(lambda params: fx(params, **kwargs),
                    dimensions=dimensions,
                    x0=x0,
                    n_jobs=kwargs.get("threads"),
                    n_calls=kwargs.get("ncalls"),
                    n_initial_points=kwargs.get("nrandom_starts"),
                    random_state=kwargs.get("seed"))

    # Serialize optimization result to file
    with open(f"{out_prefix}_{method}.pickle", "wb") as fp:
        cloudpickle.dump(res, fp)


def make_cli():
    cli = argparse.ArgumentParser()

    parser = cli.add_subparsers(dest="cmd", required=True)
    optimize = parser.add_parser("optimize",
                                 help="Run the optimization procedure to find a set of \"optimal\" parameters for modle sim")

    optimize.add_argument("--chrom-sizes",
                          type=pathlib.Path,
                          help="Path to chrom.sizes file", required=True)
    optimize.add_argument("--chrom-subranges",
                          type=pathlib.Path,
                          help="Path to BED file with chrom. subranges")
    optimize.add_argument("--diagonal-width",
                          default=int(3e6),
                          type=int)
    optimize.add_argument("--evaluation-sites-training",
                          help="Path to a BED file with the genomic coordinates to use during evaluation (training).",
                          required=True)
    optimize.add_argument("--evaluation-sites-validation",
                          help="Path to a BED file with the genomic coordinates to use during evaluation (validation).",
                          required=True)
    optimize.add_argument("--excluded-chroms",
                          default={"chrY", "chrM"},
                          nargs="*")
    optimize.add_argument("--extrusion-barriers",
                          type=pathlib.Path,
                          help="Path to extrusion barrier BED file",
                          required=True)
    optimize.add_argument("--gaussian-blur-sigma-multiplier-tgt",
                          default=1.6,
                          type=float)
    optimize.add_argument("--gaussian-blur-sigma-tgt",
                          default=1.0,
                          type=float)
    optimize.add_argument("--discretization-thresh-tgt",
                          default=1.0,
                          type=float)
    optimize.add_argument("--modle-tools-eval-metric",
                          default="custom",
                          choices={"custom", "eucl_dist", "pearson", "rmse", "spearman"})
    optimize.add_argument("--threads",
                          default=cpu_count(),
                          type=int)
    optimize.add_argument("--output-prefix",
                          type=pathlib.Path,
                          help="Output prefix",
                          required=True)
    optimize.add_argument("--transformed-reference-matrix",
                          type=pathlib.Path,
                          help="Path to the reference matrix after computing the difference of gaussians",
                          required=True)
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
    return cli


if __name__ == "__main__":
    args = make_cli().parse_args()

    # Generate kwargs from the Namespace object returned by argparse
    kwargs = vars(args)

    # Import and add eval sites to kwargs
    kwargs["sites_for_eval_training"] = import_evaluation_sites(args.evaluation_sites_training,
                                                                args.excluded_chroms,
                                                                args.chrom_subranges)
    kwargs["sites_for_eval_validation"] = import_evaluation_sites(args.evaluation_sites_validation,
                                                                  args.excluded_chroms,
                                                                  args.chrom_subranges)

    epoch = 0
    run_optimize(**kwargs)
