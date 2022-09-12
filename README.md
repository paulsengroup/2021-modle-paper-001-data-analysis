<!--
Copyright (C) 2022 Roberto Rossini <roberros@uio.no>

SPDX-License-Identifier: MIT
-->

# Synopsis

This repository contains the source code and input data used to generate results presented in "MoDLE: High-performance stochastic modeling of DNA loop extrusion interactions" (preprint available on [bioRxiv](https://doi.org/10.1101/2022.04.13.488157)).

Input data download and subsequent analyses are automated using Nextflow and Singularity/Apptainer.

## Docker images availability

Docker images are hosted on GHCR and can be found in the [Packages](https://github.com/orgs/paulsengroup/packages?repo_name=2021-modle-paper-001-data-analysis) page of this repository.

Images were generated using the `build-dockerfiles.yml` GHA workflow using the Dockerfiles from the `containers` folder.

## Nextflow workflows

Nextflow workflows under `workflows` were developed and tested using Nextflow v21.10.6, and should in principle work with any version supporting Nextflow DSL2.

Each workflow is paired with a config file (see `configs` folder). As an example, `workflows/fetch_data.nf` is paired with config `configs/fetch_data.config`.

## Requirements

- Access to an internet connection (required to download input files and Docker images)
- Nextflow v20.07.1 or newer
- Apptainer/Singularity (tested with Singularity v3.7.2)
- One or more NVIDIA GPUs (required to run MD simulations using OpenMM)
- Caper and Croo (required to run the ENCODE TF ChIP-seq pipeline, tested with v2.2.2 and v0.6.0 respectively)

## Running workflows

The `fetch_data.nf`, ENCODE `chip-seq-pipeline2` and `preprocess_data.nf` workflows must be executed first (and in this order), as they are responsible for preparing data used by the other workflows.

Remaining workflows can be executed in any order.

Inside the `config` folder there are two base configs (`base_hovig.config` and `base_saga.config`). The first config can be used to run workflows on a single node/machine without using a job scheduler, while the second config can be used to run workflows on a compute cluster using the SLURM scheduler.

Both configs are specific to the machine and cluster we used during workflow development and data analysis and will most likely need to be updated in order to run on other machines/clusters.

Refer to `run_*.sh` scripts for examples on how to run workflows individually, or `runme.sh` to run every workflow in the correct order.

Please make sure Nextflow and/or Caper are properly installed and configured before running any of the workflows.
