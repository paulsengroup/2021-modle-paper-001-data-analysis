#!/usr/bin/env python3

# Copyright (C) 2020 Edward Banigan, Aafke van den Berg, and Hugo Brandão
# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import argparse
import csv
import os
import pickle
import shutil
import sys
import time
import warnings
from collections import namedtuple

import numpy as np
import pyximport
from openmmlib import openmmlib
from openmmlib import polymerutils
from openmmlib.openmmlib import Simulation

pyximport.install(setup_args={'include_dirs': np.get_include()})
from smc_translocator_diffusive_mixed import smcTranslocatorDirectional

'''
Make sure that on average we have 1 smcstep/200 3D simulations
eg by setting
smcStepsPerBLock=1
steps=200

if smcStepsPerBlock<1, reduce the number of 3D simulations per simulation block accordingly
The speed of the smc is 1-PAUSEP.

One consideration is that the number of smc steps between 3D simulation blocks should never exceed 1 too much, as it will result in jerky motion of the polymer.
'''


def make_cli():
    class BoundCheckProbability(argparse.Action):
        def __call__(self, parser, namespace, value, option_string=None):
            if not 0.0 <= float(value) <= 1.0:
                parser.error(f"Invalid probability '{value}'")
            setattr(namespace, self.dest, float(value))

    class BoundCheckDistance(argparse.Action):
        def __call__(self, parser, namespace, value, option_string=None):
            if int(value) <= 0:
                parser.error(f"Invalid distance '{value}'")
            setattr(namespace, self.dest, int(value))

    cli = argparse.ArgumentParser()
    cli.add_argument("--gpu",
                     dest="gpu",
                     action="store_true")
    cli.add_argument("--cpu",
                     dest="gpu",
                     action="store_false")
    cli.add_argument("--monomer-size-bp",
                     type=int,
                     action=BoundCheckDistance,
                     required=True)
    cli.add_argument("--lef-processivity-bp",
                     type=int,
                     action=BoundCheckDistance,
                     required=True)
    cli.add_argument("--lef-separation-bp",
                     type=int,
                     action=BoundCheckDistance,
                     required=True)
    cli.add_argument("--pause-prob",
                     type=float,
                     action=BoundCheckProbability,
                     default=0.0)
    cli.add_argument("--barrier-occupancy",
                     type=float,
                     action=BoundCheckProbability)
    cli.add_argument("--pause-prob-passive-arm",
                     type=float,
                     action=BoundCheckProbability,
                     default=0.9)
    cli.add_argument("--fraction-asymmetric-extrusion",
                     type=float,
                     action=BoundCheckProbability,
                     default=0.0)
    cli.add_argument("--output-folder",
                     required=True,
                     type=str,
                     help="Prefix name to use for output")
    cli.add_argument("--extrusion-barrier-bed",
                     required=True,
                     type=str,
                     help="Prefix name to use for output")
    cli.add_argument("--simulation-size-mbp",
                     required=True,
                     type=float)
    cli.add_argument("--saved-simulation-blocks",
                     type=int,
                     default=5000)
    cli.add_argument("--skip-burnin",
                     dest="run_burnin",
                     action="store_false",
                     default=True)
    cli.add_argument("--force",
                     action="store_true",
                     default=False)

    cli.set_defaults(gpu=True)
    validate_args(cli, cli.parse_args())

    return cli


def validate_args(cli, args):
    if args.lef_processivity_bp % args.monomer_size_bp != 0:
        cli.error(f"--lef-processivity-bp is not a multiple of --monomer-size-bp")

    if args.lef_separation_bp % args.monomer_size_bp != 0:
        cli.error(f"--lef-processivity-bp is not a multiple of --monomer-size-bp")

    if int(args.simulation_size_mbp * 1.0e6) % args.monomer_size_bp != 0:
        cli.error(f"--simulation-size-mbp is not a multiple of --monomer-size-bp")


def save_Es_ts_Rg():
    with open(time_fname, "a+") as time_file:
        time_file.write('%f\n' % (a.state.getTime() / openmmlib.ps))
    with open(Ekin_fname, "a+") as Ekin_file:
        Ekin_file.write('%f\n' % ((a.state.getKineticEnergy()) / a.N / a.kT))
    with open(Epot_fname, "a+") as Epot_file:
        Epot_file.write('%f\n' % ((a.state.getPotentialEnergy()) / a.N / a.kT))
    # with open(Rg_fname, "a+") as Rg_file:
    #      Not sure where's analysis_plot_lib definition
    #      Rg_file.write('%f\n'%(analysis_plot_lib.rg(a.getData())))


class smcTranslocatorMilker(object):

    def __init__(self, smcTransObject):
        """
        :param smcTransObject: smc translocator object to work with
        """
        self.smcObject = smcTransObject
        self.allBonds = []

    def setParams(self, activeParamDict, inactiveParamDict):
        """
        A method to set parameters for bonds.
        It is a separate method because you may want to have a Simulation object already existing

        :param activeParamDict: a dict (argument:value) of addBond arguments for active bonds
        :param inactiveParamDict:  a dict (argument:value) of addBond arguments for inactive bonds

        """
        self.activeParamDict = activeParamDict
        self.inactiveParamDict = inactiveParamDict

    def setup(self, bondForce, blocks=100, smcStepsPerBlock=1):
        """
        A method that milks smcTranslocator object
        and creates a set of unique bonds, etc.

        :param bondForce: a bondforce object (new after simulation restart!)
        :param blocks: number of blocks to precalculate
        :param smcStepsPerBlock: number of smcTranslocator steps per block
        :return:
        """

        if len(self.allBonds) != 0:
            raise ValueError("Not all bonds were used; {0} sets left".format(len(self.allBonds)))

        self.bondForce = bondForce

        # precalculating all bonds
        allBonds = []
        for dummy in range(blocks):
            self.smcObject.steps(smcStepsPerBlock)
            left, right = self.smcObject.getSMCs()
            bonds = [(int(i), int(j)) for i, j in zip(left, right)]
            allBonds.append(bonds)

        self.allBonds = allBonds
        # 'sum' preserves order and makes one long list with bonds, 'set' creates a set with left bonds from different time points ordered from small to large and eliminates two equal bonds (also if they were created by different LEFs at different times). List turns set into a list with unique bonds at different time points.
        self.uniqueBonds = list(set(sum(allBonds, [])))

        # adding forces and getting bond indices
        self.bondInds = []
        self.curBonds = allBonds.pop(0)  # pop(0) removes and returns first list of bonds

        for bond in self.uniqueBonds:
            paramset = self.activeParamDict if (bond in self.curBonds) else self.inactiveParamDict
            ind = bondForce.addBond(bond[0], bond[1], **paramset)
            self.bondInds.append(ind)
        self.bondToInd = {i: j for i, j in zip(self.uniqueBonds, self.bondInds)}
        return self.curBonds, []

    def step(self, context, verbose=False):
        """
        Update the bonds to the next step.
        It sets bonds for you automatically!
        :param context:  context
        :return: (current bonds, previous step bonds); just for reference
        """
        if len(self.allBonds) == 0:
            raise ValueError("No bonds left to run; you should restart simulation and run setup  again")

        pastBonds = self.curBonds
        self.curBonds = self.allBonds.pop(0)  # getting current bonds
        bondsRemove = [i for i in pastBonds if i not in self.curBonds]
        bondsAdd = [i for i in self.curBonds if i not in pastBonds]
        bondsStay = [i for i in pastBonds if i in self.curBonds]
        if verbose:
            print("{0} bonds stay, {1} new bonds, {2} bonds removed".format(len(bondsStay),
                                                                            len(bondsAdd), len(bondsRemove)))
        bondsToChange = bondsAdd + bondsRemove
        bondsIsAdd = [True] * len(bondsAdd) + [False] * len(bondsRemove)
        for bond, isAdd in zip(bondsToChange, bondsIsAdd):
            ind = self.bondToInd[bond]
            paramset = self.activeParamDict if isAdd else self.inactiveParamDict
            self.bondForce.setBondParameters(ind, bond[0], bond[1], **paramset)  # actually updating bonds
        self.bondForce.updateParametersInContext(context)  # now run this to update things in the context
        return self.curBonds, pastBonds


def initModel(barriers):
    birthArray = np.ones(N, dtype=np.double) * 0.1
    deathArray = np.zeros(N, dtype=np.double) + 1. / (LIFETIME / dt)
    # Probability of stalling. Choose np.zeros(N, dtype=np.double) if no stalling
    stallLeftArray = np.zeros(N, dtype=np.double)

    stallRightArray = np.zeros(N, dtype=np.double)

    pauseArray = np.zeros(N, dtype=np.double) + PAUSEP
    slidePauseArray = np.zeros(N, dtype=np.double) + SLIDE_PAUSEP

    stallDeathArray = np.zeros(N, dtype=np.double) + 1 / (LIFETIME / dt)
    smcNum = N // SEPARATION

    for idx, direction, occupancy in barriers:
        if direction == "+":
            stallLeftArray[idx] = occupancy
        else:
            assert direction == "-"
            stallRightArray[idx] = occupancy

    oneSidedArray = np.ones(smcNum, dtype=np.int64)
    for i in range(int((1. - FRACTION_ONESIDED) * smcNum)):
        oneSidedArray[i] = 0

    belt_on_array = np.zeros(smcNum, dtype=np.double) + BELT_ON
    belt_off_array = np.zeros(smcNum, dtype=np.double) + BELT_OFF

    spf = slidePauseArray * (1. - (1. - SLIDE_PAUSEP) * np.exp(-1. * loop_prefactor))
    spb = slidePauseArray * (1. - (1. - SLIDE_PAUSEP) * np.exp(loop_prefactor))
    # sps=spf+spb
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=DeprecationWarning)
        SMCTran = smcTranslocatorDirectional(birthArray, deathArray, stallLeftArray, stallRightArray, pauseArray,
                                             stallDeathArray, smcNum, oneSidedArray, paired=PAIRED, slide=SLIDE,
                                             slidepauseForward=spf, slidepauseBackward=spb,
                                             # slidepauseSum=sps,
                                             switch=SWITCH, pushing=PUSH, belt_on=belt_on_array,
                                             belt_off=belt_off_array,
                                             SLIDE_PAUSEPROB=SLIDE_PAUSEP)

    return SMCTran


def import_barriers(path_to_bed, monomer_size_bp, default_barrier_occupancy):
    barriers = []
    BarrierT = namedtuple("BarrierT", ["idx", "direction", "occupancy"])
    with open(path_to_bed, "r", newline="") as f:
        for row in csv.reader(f, delimiter="\t"):
            assert len(row) >= 6
            idx = ((int(row[1]) + int(row[2])) // 2) // monomer_size_bp
            direction = row[5]
            if default_barrier_occupancy is not None:
                occupancy = default_barrier_occupancy
            else:
                occupancy = float(row[4])
                assert 0.0 <= occupancy <= 1.0

            if idx < N and direction in ["-", "+"]:
                barriers.append(BarrierT(idx, direction, occupancy))

    if os.path.isfile(path_to_bed) and os.path.getsize(path_to_bed) > 0 and len(barriers) == 0:
        raise RuntimeError(
            f"Unable to import any barrier from file {path_to_bed}. Is this intended?\n"
            "If you really want to run a simulation without extrusion barriers, "
            "you can do so by passing an empty file to --extrusion-barrier-bed")
    print(f"Imported {len(barriers)} from file {path_to_bed}...")
    return barriers


if __name__ == "__main__":
    args = make_cli().parse_args()

    # -------defining parameters----------
    # -- basic loop extrusion parameters--
    LIFETIME = args.lef_processivity_bp // args.monomer_size_bp  # Processivity
    SEPARATION = args.lef_separation_bp // args.monomer_size_bp  # Separation LEFs in number of monomers

    PAUSEP = args.pause_prob  # pause prob active arm, set to 0 if not using diffusive simulations
    SLIDE_PAUSEP = args.pause_prob_passive_arm  # pause prob passive arm
    FRACTION_ONESIDED = args.fraction_asymmetric_extrusion

    totalSavedBlocks = args.saved_simulation_blocks  # how many blocks to save (number of blocks done is totalSavedBlocks * saveEveryBlocks)

    N = int(
        (args.simulation_size_mbp * 1.0e6) + args.monomer_size_bp - 1) // args.monomer_size_bp  # number of monomers
    barriers = import_barriers(args.extrusion_barrier_bed, args.monomer_size_bp, args.barrier_occupancy)

    smcStepsPerBlock = 1  # I take something like 1/steppingprobability, because stepping is not determistic. I usually choose the probability of stepping to be max 0.1.
    stiff = 0  # Polymer siffness in unit of bead size
    dens = 0.2  # density in beads / volume. The density can roughly be estimated by looking at the amount of DNA in a nuclear volume.
    box = (N / dens) ** 0.33  # Define size of the bounding box for Periodic Boundary Conditions
    data = polymerutils.grow_rw(N, int(box) - 2)  # creates a compact conformation
    block = 0  # starting block

    # Time step of each smcTranslocatorPol iteration, the time step is chosen such that MaxRate*dt=0.1,
    # as this should give proper time step distributions.
    dt = 1 * (1 - PAUSEP)

    # nr of 3D simulation blocks btw advancing LEFs. For deterministic stepping choose 200-250 steps per block, otherwise, rescale with stepping probability.
    steps = int(200 * dt)

    stg = 0.8  # Probability of stalling at TAD boundary

    BELT_ON = 0
    BELT_OFF = 1
    SWITCH = 0
    PUSH = 0
    PAIRED = 0
    SLIDE = 1  # diffusive motion
    loop_prefactor = 1.5  # Should be the same as in smcTranslocator
    FULL_LOOP_ENTROPY = 1  # Should be the same as in smcTranslocator

    if dt == 0:
        sys.exit('WARNING: dt=0')

    '!!!!!!!Uncomment if entropy is considered'
    if (SLIDE_PAUSEP < 1 - np.exp(0.5 * loop_prefactor)) or (SLIDE_PAUSEP < 1 - 0.5 / np.cosh(0.5 * loop_prefactor)):
        sys.exit('WARNING: slidepause prob is too small! simulations will bias toward forward sliding!')

    # if ((PAUSEP>0.9) and (smcStepsPerBlock<10)):
    #    sys.exit('WARNING: Make sure that smc steps per block is effectively 1')

    #  10 save every 10 blocks (saving every block is now too much almost)
    saveEveryBlocks = int(100 / (smcStepsPerBlock * dt))
    # saveEveryBlocks = int(400/(smcStepsPerBlock*dt))

    # how many blocks (saved) to skip after you restart LEF positions
    skipSavedBlocksBeginning = int(20 / (smcStepsPerBlock * dt))

    restartMilkerEveryBlocks = int(100 / (smcStepsPerBlock * dt))  # int(100/(smcStepsPerBlock*dt))
    # Only one Hamiltonian can be loaded at a time to the simkt, but we want to update the bonds every time a LEF steps.
    # Loading a new Hamiltonian costs a lot of time.
    # Instead we precalculate bonds and load all positions at once as one big Hamiltonian and just change the prefactors.

    # parameters for smc bonds
    smcBondWiggleDist = 0.2
    smcBondDist = 0.5

    FullFileName = args.output_folder

    if not args.force and os.path.exists(FullFileName):
        sys.exit("Folder already exists, pass --force to overwrite")
    else:
        if os.path.exists(FullFileName):
            shutil.rmtree(FullFileName)  # Remove folder and content

    os.makedirs(FullFileName, exist_ok=True)

    # assertions for easy managing code below

    assert restartMilkerEveryBlocks % saveEveryBlocks == 0
    assert (skipSavedBlocksBeginning * saveEveryBlocks) % restartMilkerEveryBlocks == 0
    assert (totalSavedBlocks * saveEveryBlocks) % restartMilkerEveryBlocks == 0
    # max number of steps per smc block should not be too large to prevent 'jerky' polymer motion
    assert smcStepsPerBlock * (1 - PAUSEP) < 3

    savesPerMilker = restartMilkerEveryBlocks // saveEveryBlocks
    milkerInitsSkip = saveEveryBlocks * skipSavedBlocksBeginning // restartMilkerEveryBlocks
    milkerInitsTotal = (totalSavedBlocks + skipSavedBlocksBeginning) * saveEveryBlocks // restartMilkerEveryBlocks
    print(
        "Time step = {0}, Milker will be initialized {1} times, first {2} will be skipped".format(dt, milkerInitsTotal,
                                                                                                  milkerInitsSkip))

    # create filenames for Ekin, Epot and time

    Ekin_fname = os.path.join(FullFileName, 'Ekin.txt')
    Epot_fname = os.path.join(FullFileName, 'Epot.txt')
    time_fname = os.path.join(FullFileName, 'time.txt')
    Rg_fname = os.path.join(FullFileName, 'Rg.txt')
    Par_fname = os.path.join(FullFileName, 'Pars.txt')

    SMCTran = initModel(barriers)  # defining actual smc translocator object

    print(
        f"Simulating a {((N * args.monomer_size_bp) / 1.0e6):.2f}Mbp region with {len(barriers)} extrusion barriers...")

    # now polymer simulation code starts

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=DeprecationWarning)
        if args.run_burnin:
            # ------------feed smcTran to the milker---
            print("Running burn-in...")
            SMCTran.steps(1000000)  # first steps to "equilibrate" SMC dynamics. If desired of course.

    milker = smcTranslocatorMilker(SMCTran)  # now feed this thing to milker (do it once!)
    # --------- end new code ------------

    for milkerCount in range(milkerInitsTotal):
        doSave = milkerCount >= milkerInitsSkip

        # simulation parameters are defined below
        a = Simulation(timestep=80, thermostat=0.01)
        # Collision rate in inverse picoseconds, low collistion rate means ballistic like motion, default in openmmpolymer is 0.001. Motion polymer is not diffusive, this is ok for statistical average,
        # but not for dynamics of the polymer
        if args.gpu:
            # set up GPU here, PBC=Periodic Boundary Conditions. Default integrator is langevin with 300 K, friction coefficient of 1/ps, step size 0.002ps
            a.setup(platform="CUDA", PBC=True, PBCbox=[box, box, box], GPU=0, precision="mixed")
        else:
            # set up GPU here, PBC=Periodic Boundary Conditions. Default integrator is langevin with 300 K, friction coefficient of 1/ps, step size 0.002ps
            a.setup(platform="CPU", PBC=True, PBCbox=[box, box, box], precision="mixed")

        a.saveFolder(FullFileName)
        a.load(data)
        a.addHarmonicPolymerBonds(wiggleDist=0.1)  # WiggleDist controls distance at which energy of bond equals kT
        if stiff > 0:
            # Chain stiffness is introduced by an angular potential U(theta)=stiff(1-cos(theta-Pi))
            a.addGrosbergStiffness(stiff)

        # Polynomial repulsive potential between particles. Has value trunc=3.0 at zero, stays flat until 0.6-0.7 and then drops to zero.
        # For attraction between a selective set of particles, use LeonardJones or addSelectiveSSWForce (see blocks.py or ask Johannes)
        a.addPolynomialRepulsiveForce(trunc=1.5, radiusMult=1.05)

        a.step = block

        # ------------ initializing milker; adding bonds ---------
        # copied from addBond
        kbond = a.kbondScalingFactor / (smcBondWiggleDist ** 2)
        bondDist = smcBondDist * a.length_scale

        activeParams = {"length": bondDist, "k": kbond}
        inactiveParams = {"length": bondDist, "k": 0}
        milker.setParams(activeParams, inactiveParams)

        # this step actually puts all bonds in and sets first bonds to be what they should be
        milker.setup(bondForce=a.forceDict["HarmonicBondForce"],
                     blocks=restartMilkerEveryBlocks,  # default value; milk for 100 blocks
                     smcStepsPerBlock=smcStepsPerBlock)  #
        print("Restarting milker")

        a.doBlock(steps=steps, increment=False)  # do block for the first time with first set of bonds in
        for i in range(restartMilkerEveryBlocks - 1):
            curBonds, pastBonds = milker.step(a.context)  # this updates bonds. You can do something with bonds here
            if i % saveEveryBlocks == (saveEveryBlocks - 2):
                a.doBlock(steps=steps, increment=doSave)
                if doSave:
                    a.save()
                    pickle.dump(curBonds, open(os.path.join(a.folder, "SMC{0}.dat".format(a.step)), 'wb'))
                    save_Es_ts_Rg()  # save energies and time
            else:
                a.integrator.step(steps)  # do steps without getting the positions from the GPU (faster)

        data = a.getData()  # save data and step, and delete the simulation
        block = a.step
        del a

        time.sleep(0.2)  # wait 200ms for sanity (to let garbage collector do its magic)

    with open(Par_fname, "a+") as Parfile:
        Parfile.write(
            " tau=" + str(LIFETIME) + "\n Separation=" + str(SEPARATION) + "\n N=" + str(
                N) + "\n smcStepsPerBlock=" + str(
                smcStepsPerBlock) + "\n stiff=" + str(stiff) + "\n dens=" + str(dens) + "\n block=" + str(
                block) + "\n SaveEveryBlocks=" + str(saveEveryBlocks) + "\n skipSavedBlocksBeginning=" + str(
                skipSavedBlocksBeginning) + "\n totalSavedBlocks=" + str(
                totalSavedBlocks) + "\n restartMilkerEveryBlocks=" + str(
                restartMilkerEveryBlocks) + "\n smcBondWiggleDist=" + str(smcBondWiggleDist) + "\n smcBondDist=" + str(
                smcBondDist) + "\n SmcTimestep=" + str(dt) + "\n Number of barriers" + str(len(barriers)))
