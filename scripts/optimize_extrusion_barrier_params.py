#!/usr/bin/env -S python3 -u
import os
from collections import namedtuple
import argparse
import random
import bioframe as bf
import subprocess as sp
import multiprocessing as mp
import warnings
import pathlib
import cooler
import numpy as np
import glob
import logging as log
import json
import threading
import cloudpickle
import pandas as pd
import skimage
import pyBigWig
import shutil
import tempfile
from deap import algorithms, base, creator, tools
import time


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
                     help="Path to the raw reference matrix.",
                     required=True)
    cli.add_argument("--num-generations",
                     default=25,
                     type=int)
    cli.add_argument("--pop-size",
                     default=10,
                     type=int)
    cli.add_argument("--lambda",
                     default=30,
                     type=int)
    cli.add_argument("--cxpb",
                     default=0.3,
                     type=float,
                     help="Crossover probability")
    cli.add_argument("--mutpb-individual",
                     default=0.7,
                     type=float,
                     help="Mutation probability")
    cli.add_argument("--mutpb-locus",
                     default=0.075,
                     type=float,
                     help="Mutation probability")
    cli.add_argument("--mut-sigma",
                     type=float,
                     default=0.05)
    cli.add_argument("--hof-size",
                     default=10,
                     type=int,
                     help="Hall of fame size")
    cli.add_argument("--tournament-size",
                     default=5,
                     type=int,
                     help="Tournament size")
    cli.add_argument("--weak-barrier-occ-thresh",
                     type=float,
                     default=0.65,
                     help="Barrier occupancy threshold used to discriminate weak barriers.")
    cli.add_argument("--weak-barrier-fraction-thresh",
                     type=float,
                     default=0.85,
                     help="Threshold to use when purging weak barriers. Barriers that have a low occupancy in e.g. >=85\% of the population are purged")
    cli.add_argument("--weak-barrier-purge-interval",
                     type=int,
                     default=-1,
                     help="Set to 0 or negative to not purge weak barriers.")
    cli.add_argument("--seed",
                     default=1630986062,
                     type=int)
    cli.add_argument("--modle-exec",
                     default=shutil.which("modle"),
                     type=str,
                     help="Path to modle executable")
    cli.add_argument("--modle-tools-exec",
                     default=shutil.which("modle_tools"),
                     type=str,
                     help="Path to modle_tools executable")
    cli.add_argument("--initial-population",
                     type=str,
                     help="Path to a .pickle file with the initial population.")
    cli.add_argument("--ncells",
                     default=1,
                     type=int)
    cli.add_argument("--occupancy-lb",
                     type=float,
                     default=0.5)
    cli.add_argument("--occupancy-ub",
                     type=float,
                     default=0.99)
    cli.add_argument("--puu-lb",
                     type=float,
                     default=0.0)
    cli.add_argument("--puu-ub",
                     type=float,
                     default=1.0)
    return cli


def import_barriers(path_to_barriers, path_to_chrom_subranges):
    df1 = bf.read_table(path_to_barriers, names=["chrom", "start", "end", "name", "score", "strand"],
                        index_col=False)
    df2 = bf.read_table(path_to_chrom_subranges, names=["chrom", "start", "end"], index_col=False)

    return bf.overlap(df1, df2, how="left").dropna().drop(columns=["chrom_", "start_", "end_"]).reset_index(drop=True)


# https://stackoverflow.com/a/18081653
def kth_diag_indices(a, k):
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


def run_subprocess(cmd, timeout=300, attempts=3):
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


def run_modle_sim(path_to_barriers, tmp_out_prefix, bin_size, chrom_sizes, chrom_subranges, diagonal_width, ncells,
                  target_contact_density, modle_exec):
    cmd = [modle_exec, "sim",
           "-r", str(bin_size),
           "-c", str(chrom_sizes),
           "--chrom-subranges", str(chrom_subranges),
           "--contact-sampling-strategy", "loop-only-with-noise",
           "--extrusion-barrier-file", str(path_to_barriers),
           "--force",
           "--ncells", str(ncells),
           "-o", str(tmp_out_prefix),
           "--target-contact-density", str(target_contact_density),
           "--threads", str(ncells),
           "-w", str(diagonal_width)]

    run_subprocess(cmd)

    ModleSimOutputs = namedtuple("ModleSimOutputs", ["cooler", "log", "config"])

    return ModleSimOutputs(tmp_out_prefix.with_suffix(".cool"),
                           tmp_out_prefix.with_suffix(".log"),
                           pathlib.Path(f"{tmp_out_prefix}_config.toml"))


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

    run_subprocess(cmd)

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


def objective_fx(eval_sites, horizontal_bwig, vertical_bwig, bin_size, diagonal_width):
    scores = import_scores_from_bwigs(horizontal_bwig, vertical_bwig, eval_sites)

    # Deal with rows/columns of nans
    if len(scores) == 0:
        return diagonal_width // bin_size

    return np.average(scores)


def generate_barrier_annotation(barrier_annotation, barrier_params):
    assert len(barrier_annotation) == len(barrier_params), f"{len(barrier_annotation)} != {len(barrier_params)}"
    df1 = barrier_annotation.copy()
    df1["name"] = [puu for puu, _ in barrier_params]
    df1["score"] = [occ for _, occ in barrier_params]

    return df1


def write_barriers_background(path, barrier_annotation, barrier_params):

    def fx(path, barrier_annotation, barrier_params):
        new_barriers = generate_barrier_annotation(barrier_annotation, barrier_params)
        new_barriers.to_csv(path, sep="\t", header=None, index=False)

    t = threading.Thread(name="write_barriers_background",
                         target=fx,
                         args=(path, barrier_annotation, barrier_params))
    t.start()

    return t


def transform_reference_matrix(reference_matrix_uri, out_path, chrom_subranges, sigma, sigma_mult, num_diagonals_to_mask, cutoff):
    return transform_pixels(str(reference_matrix_uri),
                            str(out_path),
                            chrom_subranges,
                            sigma,
                            sigma_mult,
                            num_diagonals_to_mask,
                            cutoff,
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
             gb_cutoff_tgt):

    with tempfile.TemporaryDirectory(suffix="_modle_sim_extr_barrier_opt") as tmpdir:
        path_to_barriers = pathlib.Path(f"{tmpdir}/barriers.bed.fifo")
        tmp_out_prefix = pathlib.Path(f"{tmpdir}/out")

        os.mkfifo(path_to_barriers)

        t1 = write_barriers_background(path_to_barriers, barrier_annotation, barrier_params)

        modle_cooler = run_modle_sim(path_to_barriers,
                                     tmp_out_prefix,
                                     bin_size,
                                     chrom_sizes,
                                     chrom_subranges,
                                     diagonal_width,
                                     args.get("ncells"),
                                     args.get("target_contact_density"),
                                     args.get("modle_exec")).cooler
        t1.join()

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

        return objective_fx(barrier_annotation, bwig_horizontal, bwig_vertical, bin_size, diagonal_width),


def init_barriers(barriers, num_barriers, occ_lb, occ_ub, puu_lb, puu_ub):
    return barriers([BarrierT(random.uniform(puu_lb, puu_ub),
                              random.uniform(occ_lb, occ_ub)) for _ in range(num_barriers)])


def mutate_barriers(barriers, indpb, sigma, occ_lb, occ_ub, puu_lb, puu_ub, mu=0.0):
    new_puu = np.clip(tools.mutGaussian([puu for puu, _ in barriers], mu, sigma, indpb)[0],
                      puu_lb,
                      puu_ub)
    new_occ = np.clip(tools.mutGaussian([occ for _, occ in barriers], mu, sigma, indpb)[0],
                      occ_lb,
                      occ_ub)

    for i in range(len(barriers)):
        # noinspection PyChainedComparisons
        if barriers[i].occ == 0 and new_occ[i] != 0 and new_occ[i] < occ_lb:
            new_occ[i] = min(new_occ[i] + occ_lb, occ_ub)
        elif barriers[i].occ >= occ_lb and new_occ[i] < occ_lb:
            new_occ[i] = 0.0

    return creator.Individual([BarrierT(puu, occ) for puu, occ in zip(new_puu, new_occ)]),


def identify_weak_barriers(pop, occupancy_thresh, fraction_thresh):
    sum1 = np.zeros(len(pop[0]), dtype=int)
    sum2 = np.zeros(len(pop[0]), dtype=int)

    for individual in pop:
        sum1 += np.array([occ < occupancy_thresh for _, occ in individual], dtype=int)
        sum2 += np.array([puu == 1.0 for puu, _ in individual], dtype=int)

    mask1 = (sum1.astype(float) / len(pop)) >= fraction_thresh
    mask2 = (sum2.astype(float) / len(pop)) >= fraction_thresh

    return (mask1 | mask2)


def purge_weak_barriers_from_pop(pop, weak_barrier_mask):
    new_pop = []
    for individual in pop:
        new_pop.append(creator.Individual(
            [barrier for drop, barrier in zip(weak_barrier_mask, individual) if not drop])
        )

    return new_pop


def purge_weak_barriers_from_annotation(annotation, weak_barrier_mask):
    return annotation.drop(np.where(weak_barrier_mask)[0], axis="rows").reset_index(drop=True)


def init_or_import_population(toolbox, **kwargs):
    if kwargs.get("initial_population") is None:
        pop_size = int(kwargs.get("pop_size"))
        log.info(f"Initializing a population of size {pop_size}")
        return toolbox.population(n=pop_size)

    path_to_pickled_pop = kwargs.get("initial_population")
    log.info(f"Importing population from file {path_to_pickled_pop}")
    with open(path_to_pickled_pop, "rb") as fp:
        # Workaround for pickle not working reliably
        old_pop = cloudpickle.load(fp)
        new_pop = []
        for individual in old_pop:
            new_pop.append(creator.Individual([BarrierT(puu, occ) for puu, occ in individual]))
            # new_pop[-1].fitness.values = individual.fitness.values
        log.info(f"Imported a population of {len(new_pop)} individuals")
        return new_pop


def run_optimization(barrier_annotation, **kwargs):
    num_generations = int(kwargs.get("num_generations"))
    generation_step = int(kwargs.get("weak_barrier_purge_interval"))

    purge_weak_barriers = generation_step > 0

    if purge_weak_barriers:
        gen_per_stage = [generation_step] * int(((num_generations + generation_step - 1) / generation_step))
        weak_barrier_occ_thresh = float(kwargs.get("weak_barrier_occ_thresh"))
        weak_barrier_frac_thresh = float(kwargs.get("weak_barrier_fraction_thresh"))
    else:
        gen_per_stage = [num_generations]

    crossover_prob = float(kwargs.get("cxpb"))
    mutation_prob = float(kwargs.get("mutpb_individual"))
    hof_size = int(kwargs.get("hof_size"))

    tpool_size = int(np.ceil(int(kwargs.get("nthreads")) / int(kwargs.get("ncells"))))

    fittest_individual = None
    barrier_annotation = barrier_annotation.copy()

    with tempfile.NamedTemporaryFile() as tmpcool:
        transform_reference_matrix(
            str(kwargs.get("reference_matrix_uri")),
            tmpcool.name,
            str(kwargs.get("chrom_subranges")),
            float(kwargs.get("gaussian_blur_sigma_ref")),
            float(kwargs.get("gaussian_blur_sigma_multiplier_ref")),
            int(kwargs.get("num_diagonals_to_mask")),
            float(kwargs.get("discretization_thresh_ref")))

        kwargs["reference_cooler_transformed"] = str(tmpcool.name)

        with mp.Pool(tpool_size) as pool:
            toolbox = init_toolbox(barrier_annotation, pool.map, **kwargs)
            pop = init_or_import_population(toolbox, **kwargs)

            mu = int(kwargs.get("pop_size"))
            lambda_ = int(kwargs.get("lambda"))
            hof = tools.HallOfFame(hof_size)
            history = tools.History()
            history.update(pop)

            t0 = time.time()
            current_gen = 0
            for n in gen_per_stage:
                pop, logbook = algorithms.eaMuCommaLambda(pop, toolbox,
                                                          mu=mu,
                                                          lambda_=lambda_,
                                                          cxpb=crossover_prob,
                                                          mutpb=mutation_prob,
                                                          ngen=n,
                                                          stats=stats,
                                                          halloffame=hof,
                                                          verbose=True)
                current_gen += n
                if fittest_individual is None or fittest_individual != hof[0]:
                    fittest_individual = hof[0]
                    best_barrier_annotation = generate_barrier_annotation(barrier_annotation, fittest_individual)
                    padding = len(str(num_generations))
                    best_barrier_annotation.to_csv(f"{out_prefix}_extrusion_barriers_gen_{current_gen:0>{padding}}.bed.gz",
                                                   sep="\t",
                                                   header=None,
                                                   index=False,
                                                   compression="gzip")

                if purge_weak_barriers:
                    log.info("Purging weak barriers...")
                    num_barriers = len(barrier_annotation)
                    weak_barrier_mask = identify_weak_barriers(pop,
                                                               weak_barrier_occ_thresh,
                                                               weak_barrier_frac_thresh)
                    if np.sum(weak_barrier_mask) == 0:
                        log.info(f"No barriers were purged ({num_barriers} barriers left)")
                        continue

                    barrier_annotation = purge_weak_barriers_from_annotation(barrier_annotation,
                                                                             weak_barrier_mask)

                    toolbox = init_toolbox(barrier_annotation, pool.map, **kwargs)
                    pop = purge_weak_barriers_from_pop(pop, weak_barrier_mask)

                    history.update(pop)
                    num_barriers2 = len(barrier_annotation)
                    log.info(f"Purged {num_barriers - num_barriers2} barriers ({num_barriers2} barriers left)")

        t1 = time.time()
        elapsed_time = time.strftime("%Hh%Mm%Ss", time.gmtime(t1 - t0))
        log.info(f"Optimization took {elapsed_time} ({num_generations} generations)")
        return pop, barrier_annotation, logbook, hof, history


def init_toolbox(barrier_annotation, tpool_map,  **kwargs):
    toolbox = base.Toolbox()
    toolbox.register("individual",
                     init_barriers,
                     creator.Individual,
                     num_barriers=len(barrier_annotation),
                     occ_lb=float(kwargs.get("occupancy_lb")),
                     occ_ub=float(kwargs.get("occupancy_ub")),
                     puu_lb=float(kwargs.get("puu_lb")),
                     puu_ub=float(kwargs.get("puu_ub")))

    toolbox.register("population", tools.initRepeat, list, toolbox.individual)

    toolbox.register("mate", tools.cxTwoPoint)
    toolbox.register("mutate", mutate_barriers,
                     indpb=float(kwargs.get("mutpb_locus")),
                     sigma=float(kwargs.get("mut_sigma")),
                     occ_lb=float(kwargs.get("occupancy_lb")),
                     occ_ub=float(kwargs.get("occupancy_ub")),
                     puu_lb=float(kwargs.get("puu_lb")),
                     puu_ub=float(kwargs.get("puu_ub")))

    toolbox.register("select", tools.selTournament, tournsize=int(kwargs.get("tournament_size")))
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
                     gb_cutoff_tgt=float(kwargs.get("discretization_thresh_tgt")))

    return toolbox


if __name__ == "__main__":
    log.getLogger().setLevel(log.INFO)

    args = vars(make_cli().parse_args())
    if args.get("reference_matrix").endswith(".mcool"):
        args["reference_matrix_uri"] = args["reference_matrix"] + "::/resolutions/" + str(args.get("bin_size"))

    barrier_annotation = import_barriers(args.get("extrusion_barriers"), args.get("chrom_subranges"))
    num_barriers = len(barrier_annotation)
    log.info(f"Optimizing {num_barriers} extrusion barriers")

    random.seed(int(args.get("seed")))

    BarrierT = namedtuple("BarrierT", ["puu", "occ"])
    creator.create("FitnessMin", base.Fitness, weights=(-1.0,))
    creator.create("Individual", list, fitness=creator.FitnessMin)

    stats = tools.Statistics(lambda ind: ind.fitness.values)
    stats.register("avg", np.mean)
    stats.register("std", np.std)
    stats.register("min", np.min)
    stats.register("max", np.max)

    out_prefix = args.get("output_prefix")
    os.makedirs(os.path.dirname(out_prefix), exist_ok=True)
    pop, barrier_annotation, logbook, hof, history = run_optimization(barrier_annotation, **args)

    with open(f"{out_prefix}_population.pickle", "wb") as fp:
        cloudpickle.dump(pop, fp)

    with open(f"{out_prefix}_logbook.pickle", "wb") as fp:
        cloudpickle.dump(logbook, fp)

    with open(f"{out_prefix}_hall_of_fame.pickle", "wb") as fp:
        cloudpickle.dump(hof, fp)

    with open(f"{out_prefix}_settings.json", "w") as fp:
        print(json.dumps(args, indent=2), file=fp)

    barrier_annotation.to_csv(f"{out_prefix}_extrusion_barriers.bed.gz",
                              sep="\t",
                              header=None,
                              index=False,
                              compression="gzip")
