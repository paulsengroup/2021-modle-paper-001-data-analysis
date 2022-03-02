#!/usr/bin/env python
import argparse
import json
import os

import numpy as np
from hoomd import *
from hoomd import md, deprecated
from hoomd.deprecated.init import create_random_polymers


def make_cli():
    cli = argparse.ArgumentParser()
    cli.add_argument("--output-dir",
                     required=True,
                     type=str)
    cli.add_argument("--parameters",
                     type=str,
                     required=True,
                     help="Path to parameters.txt")
    cli.add_argument("--probabilities",
                     type=str,
                     required=True,
                     help="Path to in.probs")
    cli.add_argument("--initial-structure",
                     type=str,
                     help="Path to initial structure")
    cli.add_argument("--skip-lj-collapse",
                     type=bool,
                     default=False)
    return cli


def try_convert_to_numeric(x):
    if x.isdigit():
        return int(x)
    try:
        return float(x)
    except Exception:
        return x


# Read parameter file in the simulation folder
def read_params(param_file):
    print("*** Loading parameters")
    with open(param_file, "r") as f:
        params = {}

        for line in f:
            if "#" not in line and "=" in line:
                tokens = line.split("=")
                assert len(tokens) == 2
                params[tokens[0].strip()] = try_convert_to_numeric(tokens[1].strip())

        if "seed" not in params:
            params["seed"] = np.random.randint(1e5)

        print(json.dumps(params, sort_keys=True, indent=4))
        return params


# Generates self-avoiding walk of length N using HOOMD functions
def get_poly_coords(N, seed, max_attempts=10):
    rng = np.random.default_rng(seed)
    Lmax = N / 10

    for attempt in range(max_attempts):
        print(f"*** Generating polymer coordinates, try {attempt + 1} of 10")

        # Generate polymer coordinates using HOOMD functions.
        # Use bond lengths 1/10th the size to save memory usage.
        polymer = dict(bond_len=0.1, type=["A"] * N, bond="linear", count=1)
        system = create_random_polymers(box=data.boxdim(L=Lmax),
                                        polymers=[polymer],
                                        separation=dict(A=0.0499),
                                        seed=rng.integers(0, 1e5))

        # read out coordinates, shift each axis so that no wrapping occurs
        coords = [p.position for p in system.particles]
        coords = np.mod(np.array(coords), Lmax * 2)
        for D in range(3):
            while np.min(coords[:, D]) < 10:
                shift = np.array([0, 0, 0])
                shift[D] = -Lmax / 10
                coords = np.mod(coords + shift, Lmax * 2)

        # rescale and shift coordinates
        coords *= 10
        coords -= np.mean(coords, axis=0)

        del system
        context.initialize()

        # verify that polymer is not wrapped around the box boundaries
        if np.max(np.abs(coords[:-1, :] - coords[1:, :])) < 2.0:
            return coords

    raise NotImplementedError


# Bonds become fixed with probability according to in.probs.
# Loops reset upon hitting each other, but fixed bonds remain.
# Fixed bonds reset with probability resetprob.

def get_extrusion_bonds(probs, seed, steps=1000, nbonds=10, resetprob=0.0):
    print("Generating extrusion bond list...")
    rng = np.random.default_rng(seed + 1)

    # load sliding probabilities. If two columns, CTCF is oriented.
    with open(probs, "r") as f:
        fwd_probs, rev_probs = [], []
        for line in f:
            tokens = list(map(float, line.split()))
            assert len(tokens) == 1 or len(tokens) == 2
            if len(tokens) == 2:
                fwd_probs.append(tokens[0])
                rev_probs.append(tokens[1])
            else:
                fwd_probs.append(tokens[0])
                rev_probs.append(tokens[0])

    polylen = len(fwd_probs)

    # list of fixed anchors - locations on the polymer
    # list of fixed anchors, bool for each anchor
    fixed = np.zeros((nbonds, 2)).astype(bool)

    # list of bonds: timesteps x num_bonds x 2
    bonds = np.zeros((steps + 1, nbonds, 2)).astype(int)

    #### Initialize the bonds! ####
    # Initialize length=2 bonds at random positions, independent of boundaries.
    # Bond positions run from 1 to polylen
    start_pos = rng.integers(3, polylen - 2, nbonds)
    for i in range(nbonds):
        while np.sum(np.abs(start_pos - start_pos[i]) <= 2) >= 2:
            start_pos[i] = rng.integers(3, polylen - 2)
        bonds[0, i, 0] = start_pos[i] - 1
        bonds[0, i, 1] = start_pos[i] + 1

    #### slide the bonds! ####
    for step in range(steps):
        if (step + 1) % 1000 == 0:
            print(f"step {step + 1}...")

        # shift all the extruding bonds
        for ndx in range(nbonds):
            bond = bonds[step, ndx]

            # reset bond?
            reset_bond = False
            # fix bonds?
            if rng.random() > fwd_probs[bond[0] - 1]:
                fixed[ndx, 0] = True
            if rng.random() > rev_probs[bond[1] - 1]:
                fixed[ndx, 1] = True

            # slide bonds
            nextbond = [max(bond[0] - 1, 2),
                        min(bond[1] + 1, polylen - 2)]

            # check for fixed loop anchors
            if fixed[ndx, 0]:
                nextbond[0] = bond[0]
            if fixed[ndx, 1]:
                nextbond[1] = bond[1]

            # RESET BOND?
            # list other bond anchors
            other_anchors = np.hstack([bonds[step + 1, :ndx].flatten(),
                                       bonds[step, ndx + 1:].flatten()])
            # check if nextbond conflicts with existing bonds
            if np.any(np.abs(other_anchors - nextbond[0]) == 0):
                reset_bond = True
            if np.any(np.abs(other_anchors - nextbond[1]) == 0):
                reset_bond = True
            # check endpoints
            if nextbond[0] == 2 or nextbond[1] == polylen - 2:
                reset_bond = True

            # reset bond randomly?
            if rng.random() < resetprob:
                reset_bond = True

            if reset_bond:
                # unfix anchors
                fixed[ndx, :] = False
                retry = True
                while retry:
                    start_pos = rng.integers(3, polylen - 2)
                    retry = np.any(np.abs(other_anchors - start_pos) <= 1)
                nextbond[0] = start_pos - 1
                nextbond[1] = start_pos + 1

            # record shifted bond
            bonds[step + 1, ndx, 0] = nextbond[0]
            bonds[step + 1, ndx, 1] = nextbond[1]

    return bonds


def print_bonds(bonds, outfile):
    steps, nbonds, _ = bonds.shape
    with open(outfile, "w") as f:
        for step in range(steps - 1):
            f.write("----- Sliding bonds -----\n")

            for b in range(nbonds):
                f.write("Shift:\t(%i, %i) --> (%i, %i)\n" % \
                        tuple(bonds[step, b].tolist() + bonds[step + 1, b].tolist()))


def generate_or_load_initial_structure(length, seed, path_to_structure=None):
    if path_to_structure:
        print(f"*** Loading initial coordinates from file \"{path_to_structure}\"...")
        s = np.loadtxt(path_to_structure)
        assert len(s) == length
        return s

    # generate self-avoiding walk for initial conditions
    return get_poly_coords(length, seed)


def init_system(coordinates):
    boxsize = np.max(np.abs(coordinates)) * 2.3
    print("Boxsize:", boxsize)
    snapshot = data.make_snapshot(N=len(coordinates),
                                  box=data.boxdim(L=boxsize),
                                  particle_types=["A"],
                                  bond_types=["bb", "loop"])
    return init.read_snapshot(snapshot)


def init_beads(system):
    for ndx in range(len(system.particles)):
        p = system.particles[ndx]
        p.position = tuple(coords[ndx])

        if ndx > 0:
            p0 = system.particles[ndx - 1]
            system.bonds.add("bb", p0.tag, p.tag)


def init_bonds_and_forces(k_bb, r0_bb, k_loop, r0_loop, lj_cut, epsilon, sigma):
    harmonic = md.bond.harmonic()
    harmonic.bond_coeff.set("bb", k=k_bb, r0=r0_bb)
    harmonic.bond_coeff.set("loop", k=k_loop, r0=r0_loop)

    nl_c = md.nlist.cell(deterministic=True)
    lj = md.pair.lj(r_cut=lj_cut, nlist=nl_c)
    lj.pair_coeff.set("A", "A", epsilon=epsilon, sigma=sigma)


def init_dynamics(dt, kT, gamma, seed=5):
    md.integrate.mode_standard(dt=dt)
    lang = md.integrate.langevin(group=group.all(), kT=kT, seed=seed)
    lang.set_gamma("A", gamma=gamma)


def init_output_files(out_dir, period):
    deprecated.dump.xml(group.all(), filename=f"{out_dir}/sim_topology.xml", vis=True)  # system topology
    dump.dcd(filename=f"{out_dir}/sim_trajectory.dcd", period=period)  # trajectory


def run_simulation_loop(num_extr, steps_per_extr, nbonds):
    for sndx in range(num_extr):
        # delete existing loop bonds.
        map(lambda tag: system.bonds.remove(tag),
            (bnd.tag for bnd in system.bonds if bnd.type == "loop"))

        # set new loop bonds
        map(lambda b1, b2: system.bonds.add("loop",
                                            b1 - 1,
                                            b2 - 1),
            extr_bonds[sndx, :nbonds])

        run(steps_per_extr)


if __name__ == "__main__":
    args, hoomd_args = make_cli().parse_known_args()
    out_dir = args.output_dir
    os.makedirs(out_dir, exist_ok=True)

    context.initialize(" ".join(hoomd_args))

    params = read_params(args.parameters)
    if "seed" not in params:
        params["seed"] = np.random.randint(1e5)

    ############
    # Initialize
    coords = generate_or_load_initial_structure(params["N"],
                                                params["seed"],
                                                args.initial_structure)

    system = init_system(coords)

    init_beads(system)
    init_bonds_and_forces(params["k_bb"],
                          params["r0_bb"],
                          params["k_loop"],
                          params["r0_loop"],
                          params["lj_cut"],
                          params["epsilon"],
                          params["sigma"])

    init_dynamics(params["dt"],
                  params["T"],
                  params["gamma"])

    init_output_files(out_dir, params["print_every"])

    ###############
    # Run LJ collapse
    if not args.skip_lj_collapse:
        run(params["steps"] + 1)

    ###############
    # Extrusion
    extr_bonds = get_extrusion_bonds(probs=args.probabilities,
                                     seed=params["seed"] + 1,
                                     steps=params["nextrude"],
                                     nbonds=params["nbonds"])

    print_bonds(extr_bonds, f"{out_dir}/bond_history.txt")

    run_simulation_loop(params["nextrude"], params["steps_per_extrude"], params["nbonds"])
