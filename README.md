<!--
Copyright (C) 2022 Roberto Rossini <roberros@uio.no>

SPDX-License-Identifier: MIT
-->

# Synopsis

This repository contains the source code and input data used to generate data presented in "MoDLE: High-performance stochastic modeling of DNA loop extrusion interactions" (preprint available soon).

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

## Running workflows

The `fetch_data.nf` and `preprocess_data.nf` workflows must be executed first (and in this order), as they are responsible for preparing data used by the other workflows.

Remaining workflows can be executed in any order.

Inside the `config` folder there are two base configs (`base_hovig.config` and `base_saga.config`). The first config can be used to run workflows on a single node/machine without using a job scheduler, while the second config can be used to run workflows on a compute cluster using the SLURM scheduler.

Both configs are specific to the machine and cluster we used during workflow development and data analysis and will most likely need to be updated in order to run on other machines/clusters.

### Running workflows on a single node

```bash
nextflow run -c configs/fetch_data.config -c configs/base_hovig.config workflows/fetch_data.nf -resume
nextflow run -c configs/preprocess_data.config -c configs/base_hovig.config workflows/preprocess_data.nf -resume
nextflow run -c configs/heatmap_comparison_pt1.config -c configs/base_hovig.config workflows/heatmap_comparison_pt1.nf -resume
nextflow run -c configs/comparison_with_mut.config -c configs/base_hovig.config workflows/comparison_with_mut.nf -resume
nextflow run -c configs/gw_param_optimization.config -c configs/base_hovig.config workflows/gw_param_optimization.nf -resume
nextflow run --max_memory=400.G --max_cpus=52 --max_time=336.h --project=na -c configs/run_benchmarks.config -c configs/base_hovig.config workflows/run_benchmarks.nf -resume
```

### Running workflows on a compute cluster using SLURM

```bash
nextflow run --max_cpus=52 --max_memory=400.G --max_time=336.h --project="${SLURM_PROJECT_ID-changeme}" \
             -c configs/fetch_data.config -c configs/base_saga.config workflows/fetch_data.nf -resume
nextflow run --max_cpus=52 --max_memory=400.G --max_time=336.h --project="${SLURM_PROJECT_ID-changeme}" \
             -c configs/preprocess_data.config -c configs/base_saga.config workflows/preprocess_data.nf -resume
nextflow run --max_cpus=52 --max_memory=400.G --max_time=336.h --project="${SLURM_PROJECT_ID-changeme}" \
             -c configs/heatmap_comparison_pt1.config -c configs/base_saga.config workflows/heatmap_comparison_pt1.nf -resume
nextflow run --max_cpus=52 --max_memory=400.G --max_time=336.h --project="${SLURM_PROJECT_ID-changeme}" \
             -c configs/comparison_with_mut.config -c configs/base_saga.config workflows/comparison_with_mut.nf -resume
nextflow run --max_cpus=52 --max_memory=400.G --max_time=336.h --project="${SLURM_PROJECT_ID-changeme}" \
             -c configs/gw_param_optimization.config -c configs/base_saga.config workflows/gw_param_optimization.nf -resume
nextflow run --max_cpus=52 --max_memory=400.G --max_time=336.h --project="${SLURM_PROJECT_ID-changeme}" \
             -c configs/run_benchmarks.config -c configs/base_saga.config workflows/run_benchmarks.nf -resume
```
