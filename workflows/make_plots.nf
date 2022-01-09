#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { make_plots_md_comparison } from './modules/make_plots/modle_md_comparison.nfm' addParams(output_dir: params.output_dir)

workflow {
    make_plots_md_comparison(file(params.sanborn2015_suppl_pdf),
                             file(params.sanborn2015_plot_digests))
}