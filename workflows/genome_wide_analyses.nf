#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { compute_modle_eval_custom_metric } from './modules/modle.nfm' addParams(output_dir: params.output_dir,
                                                               ncpus: params.ncpus)
include { compare_file_names;
          dirname
        } from './modules/utils.nfm'


workflow {
    compute_modle_eval_custom_metric(file(params.modle_all_ctcf_sites_config),
                                     file(params.grch38_chrom_sizes),
                                     file(params.grch38_hela_extr_barriers_all),
                                     params.bin_size,
                                     params.diagonal_width,
                                     params.grch38_hela_reference_matrix,
                                     params.gaussian_blur_sigma_ref,
                                     params.gaussian_blur_sigma_tgt,
                                     params.gaussian_blur_multiplier_ref,
                                     params.gaussian_blur_multiplier_tgt,
                                     params.binary_discretization_value_ref,
                                     params.binary_discretization_value_tgt)
}
