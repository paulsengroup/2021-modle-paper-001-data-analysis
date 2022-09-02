#!/usr/bin/env -S python3 -u

# fmt: off
import os

os.environ["OMP_NUM_THREADS"] = "1"  # nopep8
os.environ["OPENBLAS_NUM_THREADS"] = "1"  # nopep8
os.environ["MKL_NUM_THREADS"] = "1"  # nopep8
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"  # nopep8
os.environ["NUMEXPR_NUM_THREADS"] = "1"  # nopep8
# fmt: on

import time
from deap import algorithms, base, creator, tools
import tempfile
import random
import shutil
import pyBigWig
import skimage
import pandas as pd
import cloudpickle
import threading
import json
import logging as log
import glob
import numpy as np
import cooler
import pathlib
import warnings
import multiprocessing as mp
import subprocess as sp
import bioframe as bf
import argparse
from collections import namedtuple


def make_cli():
    cli = argparse.ArgumentParser()
    cli.add_argument("--bin-size",
                     default=5000,
                     type=int)
    cli.add_argument("--diagonal-width",
                     default=int(3e6),
                     type=int)
    cli.add_argument("--chrom-sizes",
                     help="Path to .chrom.sizes file",
                     type=str,
                     required=True)
    cli.add_argument("--chrom-subranges",
                     type=str,
                     required=True,
                     help="Path to BED file with chrom. subranges")
    cli.add_argument("--extrusion-barriers",
                     help="Path to extrusion barrier BED file",
                     required=True)
    cli.add_argument("--target-contact-density",
                     type=float,
                     default=1.0)
    cli.add_argument("--modle-sigma",
                     type=float,
                     default=2000.0)
    cli.add_argument("--gaussian-blur-sigma-ref",
                     default=2.2,
                     type=float)
    cli.add_argument("--gaussian-blur-sigma-multiplier-ref",
                     default=1.05,
                     type=float)
    cli.add_argument("--discretization-thresh-ref",
                     default=0.035,
                     type=float)
    cli.add_argument("--gaussian-blur-sigma-tgt",
                     default=1.0,
                     type=float)
    cli.add_argument("--gaussian-blur-sigma-multiplier-tgt",
                     default=1.25,
                     type=float)
    cli.add_argument("--discretization-thresh-tgt",
                     default=0.015,
                     type=float)
    cli.add_argument("--num-diagonals-to-mask",
                     default=5,
                     type=int)
    cli.add_argument("--nthreads",
                     default=mp.cpu_count(),
                     type=int)
    cli.add_argument("--output-prefix",
                     help="Output prefix",
                     required=True)
    cli.add_argument("--reference-matrix",
                     help="Path to the raw (as in not transformed) reference matrix.",
                     required=True)
    cli.add_argument("--max-num-generations",
                     default=1000,
                     type=int,
                     help="Maximum number of generations to simulate on mainland.")
    cli.add_argument("--pop-size",
                     default=256,
                     type=int,
                     help="Mainland population size.")
    cli.add_argument("--lambda",
                     default=512,
                     type=int,
                     help="Number of offsprings to generate each generation (mu,lambda evolution of mainland).")
    cli.add_argument("--pop-size-island",
                     default=128,
                     type=int,
                     help="Population size for island(s).")
    cli.add_argument("--lambda-island",
                     default=256,
                     type=int,
                     help="Number of offsprings to generate each generation (mu,lambda evolution of island(s)).")
    cli.add_argument("--cxpb",
                     default=0.15,
                     type=float,
                     help="Crossover probability.")
    cli.add_argument("--mutpb-individual",
                     default=0.85,
                     type=float,
                     help="Individual mutation probability.")
    cli.add_argument("--mutpb-locus",
                     default=0.1,
                     type=float,
                     help="Locus mutation probability.")
    cli.add_argument("--mut-sigma",
                     type=float,
                     default=0.15,
                     help="Standard deviation of the normal distribution used for random mutation.")
    cli.add_argument("--hof-size",
                     default=256,
                     type=int,
                     help="Hall of fame size.")
    cli.add_argument("--seed",
                     default=1630986062,
                     type=int)
    cli.add_argument("--modle-exec",
                     default=shutil.which("modle"),
                     type=str,
                     help="Path to MoDLE executable.")
    cli.add_argument("--modle-tools-exec",
                     default=shutil.which("modle_tools"),
                     type=str,
                     help="Path to MoDLE tools executable.")
    cli.add_argument("--initial-population",
                     type=str,
                     help="Path to a .pickle of the population/hall_of_fame to use as initial mainland population")
    cli.add_argument("--ncells",
                     default=1,
                     type=int,
                     help="Number of MoDLE simulation instances or cells.")
    cli.add_argument("--occupancy-lb",
                     type=float,
                     default=0.5,
                     help="Extrusion occupancy threshold used to classify weak barriers.\n"
                          "Barriers with occupancy below this threshold are considered weak.")
    cli.add_argument("--occupancy-ub",
                     type=float,
                     default=0.975,
                     help="Extrusion occupancy threshold used to classify very strong barriers.")
    cli.add_argument("--early-stopping-window",
                     default=25,
                     type=int,
                     help="Number of generations to consider for early-stopping.")
    cli.add_argument("--early-stopping-pct",
                     type=float,
                     default=0.01,
                     help="Percentange improvement used as eatly stopping criterion.\n"
                          "Optimization is stopped early when the optimizer fails ro improve "
                          "the avg. population score by at least --early-stopping-pct over the "
                          "last --early-stopping-window generations.")
    cli.add_argument("--square-stripe-score-before-avg",
                     default=False,
                     action="store_true",
                     help="Square locus scores before computing the avg. score for an individual.")
    cli.add_argument("--num-islands",
                     type=int,
                     default=6,
                     help="Number of islands to simulate.\n"
                          "Set to 0 to only simulate mainland.\n"
                          "When this param is set to 1 the optimizer will simulate one island. Initial population of\n"
                          "this island is generated by copying mainland pop, and setting PUU and occupancy of weak\n"
                          "barriers to 1.0 and 0.0 respectively.\n"
                          "When param is set to 2 or more, one island is simulated as described above, while the rest\n"
                          "are simulated as follows. The initial pop is based on mainland population. Before starting\n"
                          "the island optimization, a random chunk of barriers (avg. size of 25) is disabled, that is\n"
                          "PUU and occupancy for these barriers is set to 1.0 and 0.0 respectively, and is not\n"
                          "allowed to change throughout the island optimization.")
    cli.add_argument("--max-generations-per-island",
                     type=int,
                     default=150)
    cli.add_argument("--selection-algorithm",
                     choices={"NSGA2", "best"},
                     default="NSGA2")
    return cli


def import_barriers(path_to_barriers, path_to_chrom_subranges):
    df1 = bf.read_table(path_to_barriers, names=["chrom", "start", "end", "name", "score", "strand"],
                        index_col=False)
    df2 = bf.read_table(path_to_chrom_subranges, names=["chrom", "start", "end"], index_col=False)

    return bf.overlap(df1, df2, how="left").dropna().drop(columns=["chrom_", "start_", "end_"]).reset_index(drop=True)


def kth_diag_indices(a, k):
    """
    https://stackoverflow.com/a/18081653
    """
    rows, cols = np.diag_indices_from(a)
    if k < 0:
        return rows[-k:], cols[:k]
    elif k > 0:
        return rows[:-k], cols[k:]
    else:
        return rows, cols


def apply_gaussian_diff(matrix, low_sigma, high_sigma, diagonals_to_mask=None):
    assert low_sigma >= 0
    assert low_sigma <= high_sigma

    matrix = np.clip(skimage.filters.difference_of_gaussians(matrix, low_sigma, high_sigma), 0, np.inf)
    if diagonals_to_mask is not None:
        assert diagonals_to_mask >= 0
        for k in range(-diagonals_to_mask, diagonals_to_mask + 1):
            row_idx, col_idx = kth_diag_indices(matrix, k)
            matrix[row_idx, col_idx] = 0

    return matrix / np.max(matrix)


def discretize_matrix(matrix, threshold):
    mask = matrix < threshold

    matrix[mask] = 0
    matrix[~mask] = 1
    return matrix


def read_pixels_from_cooler(cooler_file, chrom_name, chrom_start, chrom_end, balance=False):
    coord = f"{chrom_name}:{chrom_start}-{chrom_end}"
    m = np.nan_to_num(cooler_file.matrix(balance=balance, as_pixels=False).fetch(coord).astype(float))
    if not balance:
        return m

    mm = np.nan_to_num(cooler_file.matrix(balance=False, as_pixels=False).fetch(coord).astype(float))
    return (m / m.max()) * (mm.max() - mm.min())


def convert_dense_pixels_to_sparse(dense_pixels, chrom, start_pos, end_pos, bin_size):
    pixels = {"chrom1": [], "start1": [], "end1": [],
              "chrom2": [], "start2": [], "end2": [],
              "count": []}
    for i in range(dense_pixels.shape[0]):
        for j in range(i, dense_pixels.shape[1]):
            if dense_pixels[i][j] == 0:
                continue
            start1 = start_pos + (bin_size * i)
            end1 = min(start1 + bin_size, end_pos)

            start2 = start_pos + (bin_size * j)
            end2 = min(start2 + bin_size, end_pos)

            pixels["start1"].append(start1)
            pixels["end1"].append(end1)
            pixels["start2"].append(start2)
            pixels["end2"].append(end2)
            pixels["count"].append(dense_pixels[i][j])

    pixels["chrom1"] = [chrom] * len(pixels["start1"])
    pixels["chrom2"] = [chrom] * len(pixels["start1"])

    return pd.DataFrame(pixels)


def run_subprocess(cmd, timeout, attempts=5):
    try:
        sp.check_output(cmd,
                        encoding="utf-8",
                        stdin=sp.DEVNULL,
                        stderr=sp.STDOUT,
                        timeout=timeout)
    except sp.TimeoutExpired:
        if attempts == 0:
            raise
        log.warning(f"Subprocess {cmd[0]} {cmd[1]} timed out. Relaunching subprocess ({attempts} attempts left)")
        run_subprocess(cmd, timeout, attempts - 1)
    except sp.CalledProcessError as e:
        cmd_ = " ".join([str(x) for x in cmd])
        msg = f"Subprocess for {cmd[0]} {cmd[1]} exited with code {e.returncode}.\nCMD: {cmd_}\n\n{e.output}"
        log.error(msg)
        raise RuntimeError(msg)


def drain_named_pipe(path_to_pipe):
    with open(path_to_pipe, "rb") as f:
        # This shoudl be fine, as named pipes generated by this script should
        # never hold more than a few MBs of data
        f.read()


def run_modle_sim(barrier_annotation,
                  barrier_params,
                  out_prefix,
                  bin_size,
                  chrom_sizes,
                  chrom_subranges,
                  diagonal_width,
                  ncells,
                  target_contact_density,
                  sigma,
                  modle_exec,
                  contact_sampling_strategy="loop-only-with-noise"):
    path_to_barriers = pathlib.Path(f"{out_prefix}_barriers.bed.fifo")
    barrier_writer = write_barriers_background(path_to_barriers, barrier_annotation, barrier_params)
    cmd = [modle_exec, "sim",
           "-r", str(bin_size),
           "-c", str(chrom_sizes),
           "--chrom-subranges", str(chrom_subranges),
           "--contact-sampling-strategy", str(contact_sampling_strategy),
           "--extrusion-barrier-file", str(path_to_barriers),
           "--force",
           "--interpret-extrusion-barrier-name-as-not-bound-stp",
           "--ncells", str(ncells),
           "-o", str(out_prefix),
           "--sigma", str(sigma),
           "--target-contact-density", str(target_contact_density),
           "--threads", str(ncells),
           "-w", str(diagonal_width)]

    attempts_left = 5
    while True:
        try:
            run_subprocess(cmd, timeout=30, attempts=0)
            barrier_writer.join(timeout=0.1)
            assert not barrier_writer.is_alive()
            break
        except sp.TimeoutExpired:
            if attempts_left == 0:
                raise

            attempts_left -= 1
            barrier_writer.join(timeout=0.1)
            if barrier_writer.is_alive():
                drain_named_pipe(path_to_barriers)
            barrier_writer.join(timeout=0.1)
            assert not barrier_writer.is_alive()

            barrier_writer = write_barriers_background(path_to_barriers, barrier_annotation, barrier_params)
            log.warning(
                f"Subprocess {cmd[0]} {cmd[1]} timed out. Relaunching subprocess ({attempts_left} attempts left)")
            continue

    ModleSimOutputs = namedtuple("ModleSimOutputs", ["cooler", "log", "config"])

    return ModleSimOutputs(out_prefix.with_suffix(".cool"),
                           out_prefix.with_suffix(".log"),
                           pathlib.Path(f"{out_prefix}_config.toml"))


def transform_pixels(input_name, output_name, chrom_subranges, sigma, sigma_multiplier, num_diagonals_to_mask, cutoff,
                     balance=False):
    input_cool = cooler.Cooler(input_name)
    chroms = input_cool.chromsizes
    bin_size = input_cool.binsize
    bins = cooler.util.binnify(chroms, bin_size)

    pixels = None
    for (_, (chrom, start, end)) in pd.read_table(chrom_subranges,
                                                  names=["chrom", "start", "end"]).iterrows():
        pixels_dense = read_pixels_from_cooler(input_cool, chrom, start, end, balance)
        transformed_pixels = discretize_matrix(
            apply_gaussian_diff(pixels_dense, sigma, sigma * sigma_multiplier, num_diagonals_to_mask), cutoff)

        if pixels is None:
            pixels = convert_dense_pixels_to_sparse(transformed_pixels, chrom, start, end, bin_size)
        else:
            new_pixels = convert_dense_pixels_to_sparse(transformed_pixels, chrom, start, end, bin_size)
            pixels = pd.concat([pixels, new_pixels])

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sanitizer = cooler.create.sanitize_records(bins, schema="bg2")
        pixels = sanitizer(pixels)

    pixels = pixels[["bin1_id", "bin2_id", "count"]]
    pixels["count"] = pixels["count"].astype(int)

    log.getLogger().setLevel(log.ERROR)
    cooler.create_cooler(output_name, bins, pixels)
    log.getLogger().setLevel(log.INFO)

    return output_name


def run_modle_tools_eval(path_to_reference_matrix, path_to_target_matrix, tmp_out_prefix, chrom_subranges, bin_size,
                         diagonal_width, modle_tools_exec):
    cmd = [modle_tools_exec, "eval",
           "-i", str(path_to_target_matrix),
           "--resolution", str(bin_size),
           "--chrom-subranges", str(chrom_subranges),
           "--force",
           "--reference-matrix", str(path_to_reference_matrix),
           "-o", str(tmp_out_prefix),
           "--metric", "custom",
           "-w", str(diagonal_width),
           "--threads", "1"]

    run_subprocess(cmd, timeout=5)

    ModleToolsEvalOutputs = namedtuple("ModleToolsEvalOutputs",
                                       ["bigwig_horizontal", "bigwig_vertical", "tsv_horizontal", "tsv_vertical"])

    bwig_horizontal = glob.glob(f"{tmp_out_prefix}_*_horizontal.bw")
    bwig_vertical = glob.glob(f"{tmp_out_prefix}_*_vertical.bw")
    tsv_horizontal = glob.glob(f"{tmp_out_prefix}_*_horizontal.tsv*")
    tsv_vertical = glob.glob(f"{tmp_out_prefix}_*_vertical.tsv*")

    assert len(bwig_horizontal) == 1
    assert len(bwig_vertical) == 1
    assert len(tsv_horizontal) == 1
    assert len(tsv_vertical) == 1

    return ModleToolsEvalOutputs(bwig_horizontal[0], bwig_vertical[0],
                                 tsv_horizontal[0], tsv_vertical[0])


def import_scores_from_bwigs(horizontal_bwig, vertical_bwig, genomic_coords):
    scores = []
    with pyBigWig.open(horizontal_bwig) as h_bw, pyBigWig.open(vertical_bwig) as v_bw:
        # Fill scores vector
        for (_, (chrom, start, end, _, _, strand)) in genomic_coords.iterrows():
            if strand == "-":
                scores.append(v_bw.stats(chrom, int(start), int(end), exact=True)[0])
            else:
                assert strand == "+"
                scores.append(h_bw.stats(chrom, int(start), int(end), exact=True)[0])

    # Drop nan and inf values
    scores = np.array(scores, dtype=float)
    return scores[(~np.isnan(scores)) & (~np.isinf(scores))]


def compute_stripe_similarity(eval_sites, barrier_params, horizontal_bwig, vertical_bwig, bin_size, diagonal_width,
                              square_scores):
    scores = import_scores_from_bwigs(horizontal_bwig, vertical_bwig, eval_sites)

    if len(scores) == 0:
        return 2 * (diagonal_width // bin_size)

    scores = scores * (1.0 + barrier_params)
    if square_scores:
        return np.sqrt(np.average(scores ** 2))

    return np.average(scores)


def compute_penalties(barrier_params, occ_lb, occ_ub, exp=10):
    """
    Penalize occupancies around 0.5 as well as close to 1.0.
    The rationale is that we want the optimization do decide whether or not a given barrier should be active.
    Usually, barriers with occupancy < 0.55 do not produce visible stripes (regardless of whether the occupancy
    is e.g. 0.1 or 0.5), and so we want to push (truly) inactive barriers towards occupancy = 0.0.
    If a barrier is active, then we want to push the occupancy above 0.5 (but not too close to 1.0).
    This penalty score only makes sense when use toghether with our other similarity metrics.

    Example #1:
    A barrier is inactive according to the reference HiC matrix (i.e. there's no visible stripe/dot).
    High occupancy will be penalized by the similarity metrics, but low occupancy (i.e. <= 0.55) won't be,
    as low and very low occupancy lead to the same contact matrix in practice. The penalty score takes care of this.

    Example #2
    A barrier is active according to the reference HiC matrix (i.e. there's a visible stripe/dot).
    Low occupancies will be penalized by the similarity metrics. That being said, very high occupancies won't always be
    penalized. The penalty score takes care of this scenario: if very high and high occupancy produce similar
    accuracies, prefer the less extreme param combination.
    """
    assert occ_lb <= occ_ub

    def fx(x, exp=exp):
        return (np.power(1 + x, exp) - 1) / (2 ** exp)

    occs = get_occupancies_from_individual(barrier_params)
    occ_penalties = np.zeros_like(occs)
    mask1 = occs < occ_lb
    mask2 = (occs >= occ_lb) & (occs < occ_ub)
    mask3 = occs > occ_ub

    occ_penalties[mask1] = occs[mask1] / occ_lb
    occ_penalties[mask2] = (occ_ub - occs[mask2]) / (occ_ub - occ_lb)
    occ_penalties[mask3] = (occs[mask3] - occ_ub) / (1.0 - occ_ub)

    return fx(occ_penalties)


def generate_barrier_annotation(barrier_annotation, barrier_params):
    assert len(barrier_annotation) == len(barrier_params), f"{len(barrier_annotation)} != {len(barrier_params)}"
    df1 = barrier_annotation.copy()
    df1["name"] = get_puus_from_individual(barrier_params)
    df1["score"] = get_occupancies_from_individual(barrier_params)

    return df1


def write_barriers_background(path, barrier_annotation, barrier_params):
    if path.exists():
        path.unlink()
    os.mkfifo(path)

    def fx(path, barrier_annotation, barrier_params):
        new_barriers = generate_barrier_annotation(barrier_annotation, barrier_params)
        new_barriers.to_csv(path, sep="\t", header=None, index=False)

    writer = threading.Thread(target=fx,
                              args=(path, barrier_annotation, barrier_params))
    writer.start()

    return writer


def transform_reference_matrix(out_path, **kwargs):
    return transform_pixels(
        str(kwargs.get("reference_matrix_uri")),
        str(out_path),
        str(kwargs.get("chrom_subranges")),
        float(kwargs.get("gaussian_blur_sigma_ref")),
        float(kwargs.get("gaussian_blur_sigma_multiplier_ref")),
        int(kwargs.get("num_diagonals_to_mask")),
        float(kwargs.get("discretization_thresh_ref")),
        balance=True)


def evaluate(barrier_params,
             barrier_annotation,
             bin_size,
             diagonal_width,
             chrom_sizes,
             chrom_subranges,
             reference_cooler_transformed,
             gb_num_diagonals_to_mask,
             gb_sigma_tgt,
             gb_sigma_mult_tgt,
             gb_cutoff_tgt,
             occ_lb,
             occ_ub,
             square_stripe_score):
    with tempfile.TemporaryDirectory(suffix="_modle_sim_extr_barrier_opt") as tmpdir:
        tmp_out_prefix = pathlib.Path(f"{tmpdir}/out")

        modle_cooler = run_modle_sim(barrier_annotation,
                                     barrier_params,
                                     tmp_out_prefix,
                                     bin_size,
                                     chrom_sizes,
                                     chrom_subranges,
                                     diagonal_width,
                                     args.get("ncells"),
                                     args.get("target_contact_density"),
                                     args.get("modle_sigma"),
                                     args.get("modle_exec")).cooler

        modle_cooler_transformed = transform_pixels(str(modle_cooler),
                                                    str(modle_cooler.with_suffix("")) + "_transformed.cool",
                                                    chrom_subranges,
                                                    gb_sigma_tgt,
                                                    gb_sigma_mult_tgt,
                                                    gb_num_diagonals_to_mask,
                                                    gb_cutoff_tgt,
                                                    balance=False)

        bwig_horizontal, bwig_vertical, _, _ = run_modle_tools_eval(reference_cooler_transformed,
                                                                    modle_cooler_transformed,
                                                                    tmp_out_prefix,
                                                                    chrom_subranges,
                                                                    bin_size,
                                                                    diagonal_width,
                                                                    args.get("modle_tools_exec"))

        penalties = compute_penalties(barrier_params,
                                      occ_lb=occ_lb,
                                      occ_ub=occ_ub)

        score = compute_stripe_similarity(barrier_annotation,
                                          penalties,
                                          bwig_horizontal,
                                          bwig_vertical,
                                          bin_size,
                                          diagonal_width,
                                          square_stripe_score)

        return score,


def init_barriers(barriers, num_barriers, prng):
    return barriers([BarrierT(puu, occ) for puu, occ in zip(prng.uniform(size=num_barriers),
                                                            prng.uniform(size=num_barriers))])


def mutate_barriers(barriers, prng, indpb, sigma, mask=None, mu=0.0):
    num_barriers = len(barriers)

    if mask is None:
        mask = np.full(num_barriers, True)

    mask = mask & (prng.uniform(size=num_barriers) <= indpb)
    num_mutations = mask.sum()

    gauss = np.zeros_like(mask, dtype=float)
    gauss[mask] = prng.normal(loc=mu, scale=sigma, size=num_mutations)
    new_puu = np.clip(get_puus_from_individual(barriers) + gauss, 0.0, 1.0)

    gauss = np.zeros_like(mask, dtype=float)
    gauss[mask] = prng.normal(loc=mu, scale=sigma, size=num_mutations)
    new_occ = np.clip(get_occupancies_from_individual(barriers) + gauss, 0.0, 1.0)

    return creator.Individual([BarrierT(puu, occ) for puu, occ in zip(new_puu, new_occ)]),


def mask_barriers_random(population, prng, avg_num_barriers_to_disable=25, std=5.0):
    num_barriers = len(population[0])
    i = prng.integers(0, num_barriers)
    num_barriers_to_disable = int(min(prng.normal(loc=avg_num_barriers_to_disable, scale=std), num_barriers))
    i0 = max(0, int(i - (num_barriers_to_disable / 2)))
    i1 = min(num_barriers, i0 + num_barriers_to_disable)

    if i0 == i1:
        return np.full(num_barriers, True), population.copy()

    new_pop = [creator.Individual([x for x in ind]) for ind in population]
    for ind in new_pop:
        for i in range(i0, i1):
            ind[i] = BarrierT(1.0, 0.0)

    mask = np.full(num_barriers, True)
    mask[i0:i1 - 1] = False
    return mask, new_pop


def mask_weak_barriers(population, occ_threshold=0.5, puu_threshold=1.0):
    new_pop = []

    for ind in population:
        new_ind = [x for x in ind]
        for i in range(len(new_ind)):
            if new_ind[i].occ <= occ_threshold or new_ind[i].puu >= puu_threshold:
                new_ind[i] = BarrierT(1.0, 0.0)
        new_pop.append(creator.Individual(new_ind))

    return new_pop


def init_population(toolbox, pop_size):
    log.info(f"Initializing a population of size {pop_size}")
    return toolbox.population(n=pop_size)


def import_population(path_to_pickled_pop):
    log.info(f"Importing population from file {path_to_pickled_pop}")
    with open(path_to_pickled_pop, "rb") as fp:
        # Workaround for pickle not working reliably
        old_pop = cloudpickle.load(fp)
        new_pop = []
        for individual in old_pop:
            new_pop.append(creator.Individual([BarrierT(puu, occ) for puu, occ in individual]))
            new_pop[-1].fitness.values = individual.fitness.values
        log.info(f"Imported a population of {len(new_pop)} individuals")
        return new_pop


def init_or_import_population(toolbox, **kwargs):
    pop_size = int(kwargs.get("pop_size"))
    if kwargs.get("initial_population") is None:
        return init_population(toolbox, pop_size)

    return import_population(kwargs.get("initial_population"))


def simulate_mainland(barrier_annotation,
                      toolbox,
                      pop,
                      hof,
                      max_num_gens,
                      history,
                      stats,
                      early_stopping_offset,
                      logbook=None,
                      **kwargs):
    output = simulate_island(island_id=0,
                             barrier_annotation=barrier_annotation,
                             toolbox=toolbox,
                             pop=pop,
                             hof=hof,
                             max_num_gens=max_num_gens,
                             history=history,
                             stats=stats,
                             early_stopping_offset=early_stopping_offset,
                             logbook=logbook,
                             **kwargs)

    _, _, barrier_annotation, logbook, hof, _ = output
    out_prefix = str(kwargs.get("output_prefix")) + f"_mainland_gen{len(logbook):04}"
    run_modle_sim(barrier_annotation,
                  barrier_params=hof[0],
                  out_prefix=pathlib.Path(out_prefix),
                  bin_size=int(kwargs.get("bin_size")),
                  chrom_sizes=str(kwargs.get("chrom_sizes")),
                  chrom_subranges=str(kwargs.get("chrom_subranges")),
                  diagonal_width=int(kwargs.get("diagonal_width")),
                  ncells=int(kwargs.get("nthreads")),
                  sigma=args.get("modle_sigma"),
                  target_contact_density=max(100.0, float(kwargs.get("target_contact_density"))),
                  modle_exec=kwargs.get("modle_exec"),
                  contact_sampling_strategy="tad-plus-loop-with-noise")
    return output


def simulate_island(island_id,
                    barrier_annotation,
                    toolbox,
                    pop,
                    hof,
                    max_num_gens,
                    mu=None,
                    lambda_=None,
                    history=None,
                    stats=None,
                    logbook=None,
                    early_stopping_offset=0,
                    **kwargs):
    crossover_prob = float(kwargs.get("cxpb"))
    mutation_prob = float(kwargs.get("mutpb_individual"))

    if mu is None:
        mu = len(pop)

    if lambda_ is None:
        lambda_ = max(mu, int(kwargs.get("lambda")))

    if history is None:
        history = tools.History()
        history.update(pop)

    if island_id == 0:
        logfile = kwargs.get("output_prefix") + "_mainland.log"
    else:
        logfile = kwargs.get("output_prefix") + f"_isl{island_id:03d}.log"

    def stopping_fx(avg_scores, score_stds):
        return run_early_stopping_check(avg_scores=avg_scores,
                                        score_std=score_stds[-1],
                                        offset=early_stopping_offset,
                                        window=kwargs.get("early_stopping_window"),
                                        improvement_thresh=kwargs.get("early_stopping_pct"))

    t0 = time.time()
    gen, pop, logbook = ea_mu_comma_lambda(pop, toolbox,
                                           mu=mu,
                                           lambda_=lambda_,
                                           cxpb=crossover_prob,
                                           mutpb=mutation_prob,
                                           max_ngen=max_num_gens,
                                           stats=stats,
                                           halloffame=hof,
                                           logbook=logbook,
                                           stopping_criterion=stopping_fx,
                                           logfile=logfile)

    t1 = time.time()
    elapsed_time = time.strftime("%Hh%Mm%Ss", time.gmtime(t1 - t0))
    if island_id == 0:
        log.info(f"Mainland optimization took {elapsed_time} ({gen} generations)")
    else:
        log.info(f"Optimization of island {island_id:03d} took {elapsed_time} ({gen} generations)")
    barrier_annotation = generate_barrier_annotation(barrier_annotation, hof[0])
    return gen, pop, barrier_annotation, logbook, hof, history


def run_optimization(barrier_annotation, **kwargs):
    assert int(kwargs.get("num_islands")) >= 0

    max_num_generations = int(kwargs.get("max_num_generations"))
    hof_size = int(kwargs.get("hof_size"))
    tpool_size = int(np.ceil(int(kwargs.get("nthreads")) / int(kwargs.get("ncells"))))

    num_islands = int(kwargs.get("num_islands"))
    score_improvement_thresh = float(kwargs.get("early_stopping_pct"))

    island_id = 1
    num_generations = 0
    with tempfile.NamedTemporaryFile() as tmpcool, mp.Pool(tpool_size) as pool:
        transform_reference_matrix(tmpcool.name, **kwargs)
        kwargs["reference_cooler_transformed"] = str(tmpcool.name)

        toolbox = init_toolbox(barrier_annotation, pool.map, **kwargs)
        stats = init_stats("mainland", 0)
        logbook = None
        pop = init_or_import_population(toolbox, **kwargs)
        hof = tools.HallOfFame(hof_size)
        history = tools.History()
        history.update(pop)
        islands = {}
        latest_score = np.inf

        while num_generations < max_num_generations:
            num_generations_left = max_num_generations - num_generations
            if len(islands) != 0:
                pop = migrate(pop, islands, kwargs.get("prng"))

            gen, pop, barrier_annotation, logbook, hof, history = simulate_mainland(
                barrier_annotation,
                toolbox,
                pop,
                hof=hof,
                max_num_gens=num_generations_left,
                history=history,
                stats=stats,
                early_stopping_offset=num_generations,
                logbook=logbook,
                **kwargs)

            current_score = logbook.select("avg_score")[-1]
            num_generations = len(logbook)
            if num_generations == max_num_generations:
                log.info(
                    f"Terminating optimization as the maximum number of generations ({max_num_generations}) "
                    "has been reached")
                break

            score_improvement = 1.0 - (current_score / latest_score)
            if score_improvement < score_improvement_thresh:
                log.info(
                    "Terminating optimization as the last batch of generations failed to significantly improve the avg "
                    f"score ({score_improvement:.4f} < {score_improvement_thresh:.4f})")
                break
            latest_score = current_score

            if num_islands == 0:
                break

            islands = {}
            if num_islands > 1:
                for island_id in tuple(range(island_id, island_id + num_islands - 1)):
                    mask, isl_pop = mask_barriers_random(pop, kwargs.get("prng"))

                    isl_toolbox = init_toolbox(barrier_annotation, pool.map, mask=mask, **kwargs)
                    isl_hof = tools.HallOfFame(hof_size)

                    _, isl_pop, _, isl_logbook, isl_hof, _ = simulate_island(island_id,
                                                                             barrier_annotation,
                                                                             isl_toolbox,
                                                                             isl_pop,
                                                                             mu=int(kwargs["pop_size_island"]),
                                                                             lambda_=int(kwargs["lambda_island"]),
                                                                             hof=isl_hof,
                                                                             max_num_gens=int(
                                                                                 kwargs["max_generations_per_island"]),
                                                                             stats=init_stats(f"island_{island_id:03d}",
                                                                                              np.sum(~mask)),
                                                                             **kwargs)
                    islands[island_id] = isl_pop

                    write_optimization_state(args.get("output_prefix") + f"_island_{island_id:03d}",
                                             isl_pop,
                                             isl_logbook,
                                             isl_hof)

            isl_pop = mask_weak_barriers(pop, occ_threshold=kwargs.get("occupancy_lb"))
            isl_toolbox = init_toolbox(barrier_annotation, pool.map, mask=None, **kwargs)
            isl_hof = tools.HallOfFame(hof_size)

            island_id += 1
            _, isl_pop, _, isl_logbook, isl_hof, _ = simulate_island(island_id,
                                                                     barrier_annotation,
                                                                     isl_toolbox,
                                                                     isl_pop,
                                                                     mu=int(kwargs["pop_size_island"]),
                                                                     lambda_=int(kwargs["lambda_island"]),
                                                                     hof=isl_hof,
                                                                     max_num_gens=int(kwargs["max_generations_per_island"]),
                                                                     stats=init_stats(f"island_{island_id:03d}", 0),
                                                                     **kwargs)
            islands[island_id] = isl_pop
            write_optimization_state(args.get("output_prefix") + f"_island_{island_id:03d}", isl_pop, isl_logbook,
                                     isl_hof)
            island_id += 1

    write_optimization_state(args.get("output_prefix") + "_mainland", pop, logbook, hof)
    return gen, pop, barrier_annotation, logbook, hof, history


def migrate(mainland_pop, islands, prng, fraction=0.5):
    island_pops = []
    weights = []

    for pop in islands.values():
        island_pops.extend([ind for ind in pop])
        weights.extend([ind.fitness.values[0] for ind in pop])

    old_pop = mainland_pop.copy()
    old_pop.sort()

    weights = np.finfo(float).eps + np.max(weights) - np.array(weights)
    weights = weights / np.sum(weights)
    sample_size = int(round(len(mainland_pop) * fraction))
    idx = prng.choice(len(island_pops), size=sample_size, replace=False, p=weights)
    return old_pop[:-sample_size] + [island_pops[i] for i in idx]


def write_optimization_state(out_prefix, pop, logbook, hof):
    with open(f"{out_prefix}_population.pickle", "wb") as fp:
        cloudpickle.dump(pop, fp)

    with open(f"{out_prefix}_logbook.pickle", "wb") as fp:
        cloudpickle.dump(logbook, fp)

    with open(f"{out_prefix}_hall_of_fame.pickle", "wb") as fp:
        cloudpickle.dump(hof, fp)


def init_toolbox(barrier_annotation, tpool_map, mask=None, **kwargs):
    toolbox = base.Toolbox()
    toolbox.register("individual",
                     init_barriers,
                     creator.Individual,
                     num_barriers=len(barrier_annotation),
                     prng=kwargs.get("prng"))

    toolbox.register("population", tools.initRepeat, list, toolbox.individual)

    toolbox.register("mate", tools.cxTwoPoint)
    toolbox.register("mutate", mutate_barriers,
                     prng=kwargs.get("prng"),
                     indpb=float(kwargs.get("mutpb_locus")),
                     sigma=float(kwargs.get("mut_sigma")),
                     mask=mask)

    if kwargs["selection_algorithm"] == "NSGA2":
        toolbox.register("select", tools.selNSGA2)
    else:
        assert kwargs["selection_algorithm"] == "best"
        toolbox.register("select", tools.selBest)

    toolbox.register("map", tpool_map)

    toolbox.register("evaluate", evaluate,
                     barrier_annotation=barrier_annotation,
                     bin_size=int(kwargs.get("bin_size")),
                     diagonal_width=int(kwargs.get("diagonal_width")),
                     chrom_sizes=str(kwargs.get("chrom_sizes")),
                     chrom_subranges=str(kwargs.get("chrom_subranges")),
                     reference_cooler_transformed=str(kwargs.get("reference_cooler_transformed")),
                     gb_num_diagonals_to_mask=int(kwargs.get("num_diagonals_to_mask")),
                     gb_sigma_tgt=float(kwargs.get("gaussian_blur_sigma_tgt")),
                     gb_sigma_mult_tgt=float(kwargs.get("gaussian_blur_sigma_multiplier_tgt")),
                     gb_cutoff_tgt=float(kwargs.get("discretization_thresh_tgt")),
                     occ_lb=float(kwargs.get("occupancy_lb")),
                     occ_ub=float(kwargs.get("occupancy_ub")),
                     square_stripe_score=kwargs.get("square_stripe_score_before_avg"))

    return toolbox


def init_stats(label, num_barriers_masked):
    def compute_fraction_of_barriers_with_occ_in_range(pop, lb, ub):
        n1 = len(pop) * len(pop[0])
        n2 = 0
        for ind in pop:
            n2 += np.sum([lb <= occ < ub for _, occ in ind])
        return n2 / n1

    stats = tools.Statistics()
    stats.register("label", lambda _: label)
    stats.register("num_barriers_masked", lambda _: num_barriers_masked)
    stats.register("avg_score", lambda pop: np.mean([ind.fitness.values[0] for ind in pop]))
    stats.register("std_score", lambda pop: np.std([ind.fitness.values[0] for ind in pop]))
    stats.register("min_score", lambda pop: np.min([ind.fitness.values[0] for ind in pop]))
    stats.register("max_score", lambda pop: np.max([ind.fitness.values[0] for ind in pop]))

    stats.register("occ_lt01", lambda pop: compute_fraction_of_barriers_with_occ_in_range(pop, 0.0, 0.01))
    stats.register("occ_0_50", lambda pop: compute_fraction_of_barriers_with_occ_in_range(pop, 0.0, 0.5))
    stats.register("occ_50_60", lambda pop: compute_fraction_of_barriers_with_occ_in_range(pop, 0.5, 0.6))
    stats.register("occ_60_70", lambda pop: compute_fraction_of_barriers_with_occ_in_range(pop, 0.6, 0.7))
    stats.register("occ_70_80", lambda pop: compute_fraction_of_barriers_with_occ_in_range(pop, 0.7, 0.8))
    stats.register("occ_80_90", lambda pop: compute_fraction_of_barriers_with_occ_in_range(pop, 0.8, 0.9))
    stats.register("occ_90_100", lambda pop: compute_fraction_of_barriers_with_occ_in_range(pop, 0.9, 1.01))
    stats.register("occ_gt99", lambda pop: compute_fraction_of_barriers_with_occ_in_range(pop, 0.99, 1.01))

    return stats


def init_logbook(stats):
    logbook = tools.Logbook()
    logbook.header = ["gen", "nevals"] + (stats.fields if stats else [])
    return logbook


def run_early_stopping_check(avg_scores, score_std, window, improvement_thresh, offset=0, std_threshold=0.05):
    assert offset <= len(avg_scores), f"{offset} > {len(avg_scores)}"
    if score_std < std_threshold:
        log.info(f"Stopping as score variability is very low (std={score_std:.4f})")
        return True

    if len(avg_scores) - offset <= window:
        return False

    def compute_improvement_pct(vect, how="minimize", offset=offset, window=window):
        if how == "minimize":
            current_best = np.min(vect[-window:])
            previous_best = np.min(vect[offset:-window])

            return 1.0 - (current_best / previous_best)

        assert how == "maximize"
        current_best = np.max(vect[-window:])
        previous_best = np.max(vect[:-window])

        return (current_best / previous_best) - 1.0

    score_improvement = compute_improvement_pct(avg_scores, how="minimize")
    if score_improvement < improvement_thresh:
        log.info(
            "Stopping early as optimization failed to significantly improve any of the scoring metrics over the "
            f"last {window} generations: "
            f"score_improvement={score_improvement};")
        return True

    return False


def ea_mu_comma_lambda(population, toolbox, mu, lambda_, cxpb, mutpb, max_ngen, stopping_criterion,
                       stats=None, halloffame=None, logfile=None, logbook=None):
    # Based on eaMuCommaLambda from https://github.com/DEAP/deap/blob/master/deap/algorithms.py
    assert lambda_ >= mu, "lambda must be greater or equal to mu."

    # the individuals with an invalid fitness
    invalid_ind = [ind for ind in population if not ind.fitness.valid]
    fitnesses = toolbox.map(toolbox.evaluate, invalid_ind)
    for ind, fit in zip(invalid_ind, fitnesses):
        ind.fitness.values = fit

    if halloffame is not None:
        halloffame.update(population)

    if logbook is None:
        logbook = init_logbook(stats)

    gen = len(logbook)
    record = stats.compile(population) if stats is not None else {}
    logbook.record(gen=gen, nevals=len(invalid_ind), **record)
    if logfile:
        with open(logfile, "a") as f:
            print(logbook.stream, file=f)

    # Begin the generational process
    for gen in range(gen + 1, gen + max_ngen + 1):
        # Vary the population
        offspring = algorithms.varOr(population, toolbox, lambda_, cxpb, mutpb)

        #  the individuals with an invalid fitness
        invalid_ind = [ind for ind in offspring if not ind.fitness.valid]

        fitnesses = toolbox.map(toolbox.evaluate, invalid_ind)
        for ind, fit in zip(invalid_ind, fitnesses):
            ind.fitness.values = fit

        # Update the hall of fame with the generated individuals
        if halloffame is not None:
            halloffame.update(offspring)

        # Select the next generation population
        population[:] = toolbox.select(offspring, mu)

        # Update the statistics with the new population
        record = stats.compile(population) if stats is not None else {}
        logbook.record(gen=gen, nevals=len(invalid_ind), **record)
        if logfile:
            with open(logfile, "a") as f:
                print(logbook.stream, file=f)

        if stopping_criterion(logbook.select("avg_score"),
                              logbook.select("std_score")):
            break

    return gen, population, logbook


def get_occupancies_from_individual(ind):
    return np.array([occ for _, occ in ind], dtype=float)


def get_puus_from_individual(ind):
    return np.array([puu for puu, _ in ind], dtype=float)


if __name__ == "__main__":
    log.getLogger().setLevel(log.INFO)

    args = vars(make_cli().parse_args())
    if args.get("reference_matrix").endswith(".mcool"):
        args["reference_matrix_uri"] = args["reference_matrix"] + "::/resolutions/" + str(args.get("bin_size"))

    barrier_annotation = import_barriers(args.get("extrusion_barriers"), args.get("chrom_subranges"))
    num_barriers = len(barrier_annotation)
    log.info(f"Optimizing {num_barriers} extrusion barriers")

    args["prng"] = np.random.default_rng(int(args.get("seed")))
    random.seed(int(args.get("seed")))

    BarrierT = namedtuple("BarrierT", ["puu", "occ"])
    creator.create("FitnessMulti", base.Fitness, weights=(-1.0, ))
    creator.create("Individual", list, fitness=creator.FitnessMulti)

    out_prefix = args.get("output_prefix")
    os.makedirs(os.path.dirname(out_prefix), exist_ok=True)
    gen, pop, barrier_annotation, logbook, hof, history = run_optimization(barrier_annotation, **args)

    write_optimization_state(out_prefix, pop, logbook, hof)
    with open(f"{out_prefix}_settings.json", "w") as fp:
        args.pop("prng")
        print(json.dumps(args, indent=2), file=fp)

    barrier_annotation.to_csv(f"{out_prefix}_extrusion_barriers.bed.gz",
                              sep="\t",
                              header=None,
                              index=False,
                              compression="gzip")
