#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { make_heatmap_comparison } from './modules/make_plots/heatmap_comparison.nfm' params(params)
//include { plot_benchmarks } from './modules/make_plots/benchmarks.nfm'

workflow {
    make_heatmap_comparison()
}
