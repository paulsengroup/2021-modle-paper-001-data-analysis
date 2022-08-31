#!/usr/bin/env nextflow
// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

nextflow.enable.dsl=2

workflow {
    basenames = Channel.of("GRCh38_H1_optimized_barriers_hic", "GRCh38_H1_optimized_barriers_microc")
    sigmas = Channel.of(params.hic_low_sigma, params.microc_low_sigma)
    sigma_mults = Channel.of(params.hic_sigma_mult, params.microc_sigma_mult)
    discr_thresholds = Channel.of(params.hic_discr_thresh, params.microc_discr_thresh)

    ref_matrices = Channel.of(file(params.reference_hic), file(params.reference_microc))

    optimize_extrusion_barriers(basenames,
                                file(params.chrom_sizes),
                                file(params.candidate_barriers),
                                file(params.regions_of_interest),
                                params.bin_size,
                                params.target_contact_density,
                                sigmas,
                                sigma_mults,
                                discr_thresholds,
                                params.modle_low_sigma,
                                params.modle_sigma_mult,
                                params.modle_discr_thresh,
                                ref_matrices)


    optimize_extrusion_barriers.out.optimized_barriers
                               .branch {
                                            hic: it =~ /optimized_barriers_hic/
                                            microc: true
                                       }
                               .set {extrusion_barriers}


    modle_configs = Channel.of(file(params.hic_modle_config), file(params.microc_modle_config))
    run_modle_sim(modle_configs,
                  file(params.chrom_sizes),
                  Channel.empty().concat(extrusion_barriers.hic, extrusion_barriers.microc),
                  file(params.regions_of_interest))

    cooler_zoomify(run_modle_sim.out.cool)
}

process optimize_extrusion_barriers {
    publishDir "${params.output_dir}/optimized_barriers", mode: 'copy'

    label 'process_very_high'
    label 'process_very_long'

    input:
        val basename
        path chrom_sizes
        path extr_barriers
        path regions_of_interest
        val bin_size
        val target_contact_density

        val gaussian_blur_sigma_ref
        val gaussian_blur_mult_ref
        val discretization_thresh_ref

        val gaussian_blur_sigma_tgt
        val gaussian_blur_mult_tgt
        val discretization_thresh_tgt

        path reference_matrix

    output:
        path "*_extrusion_barriers.bed.gz", emit: optimized_barriers
        path "*.tar.gz", emit: result
        path "*.log", emit: log


    shell:
        '''
        '!{params.script_dir}/optimize_extrusion_barrier_params_islands.py' \
            --bin-size '!{bin_size}' \
            --extrusion-barriers '!{extr_barriers}' \
            --output-prefix '!{basename}/!{basename}' \
            --chrom-sizes '!{chrom_sizes}' \
            --chrom-subranges '!{regions_of_interest}' \
            --target-contact-density '!{target_contact_density}' \
            --nthreads '!{task.cpus}' \
            --reference-matrix '!{reference_matrix}' \
            --gaussian-blur-sigma-ref '!{gaussian_blur_sigma_ref}' \
            --gaussian-blur-sigma-multiplier-ref '!{gaussian_blur_mult_ref}' \
            --discretization-thresh-ref '!{discretization_thresh_ref}' \
            --gaussian-blur-sigma-tgt '!{gaussian_blur_sigma_tgt}' \
            --gaussian-blur-sigma-multiplier-tgt '!{gaussian_blur_mult_tgt}' \
            --discretization-thresh-tgt '!{discretization_thresh_tgt}' |& tee '!{basename}.log'

        cp '!{basename}/!{basename}_extrusion_barriers.bed.gz' .
        tar -czf '!{basename}.tar.gz' '!{basename}'
        rm -r '!{basename}'
        '''
}


process run_modle_sim {

    label 'process_medium'

    input:
        path config
        path chrom_sizes
        path extr_barriers
        path regions_of_interest


    output:
        path "*.cool", emit: cool

    shell:
        out_prefix=extr_barriers.getSimpleName() - /_final_annotation/
        '''
        modle sim --config '!{config}'                       \
                  -c '!{chrom_sizes}'                        \
                  -b '!{extr_barriers}'                      \
                  --chrom-subranges '!{regions_of_interest}' \
                  -t !{task.cpus}                            \
                  -o '!{out_prefix}'
        '''


}

process cooler_zoomify {
    publishDir "${params.output_dir}/mcools", mode: 'copy'

    label 'process_medium'

    input:
        path cool

    output:
        path "*.mcool", emit: mcool

    shell:
        '''
        cooler zoomify -p !{task.cpus} -r N '!{cool}'
        '''
}
