#!/usr/bin/env nextflow
// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

nextflow.enable.dsl=2

workflow {
    optimize_extrusion_barriers_001("GRCh38_H1_optimized_barriers_microc_001",
                                    file(params.chrom_sizes),
                                    file(params.candidate_barriers),
                                    file(params.regions_of_interest),
                                    params.bin_size,
                                    params.target_contact_density,
                                    params.microc_low_sigma,
                                    params.microc_sigma_mult,
                                    params.microc_discr_thresh,
                                    params.modle_low_sigma,
                                    params.modle_sigma_mult,
                                    params.modle_discr_thresh,
                                    file(params.reference_microc),
                                    params.num_islands_001,
                                    params.cxpb_001,
                                    params.mutpb_ind_001,
                                    params.mutpb_locus_001,
                                    params.mutsigma_001)

    optimize_extrusion_barriers_002("GRCh38_H1_optimized_barriers_microc_002",
                                    file(params.chrom_sizes),
                                    file(params.candidate_barriers),
                                    file(params.regions_of_interest),
                                    params.bin_size,
                                    params.target_contact_density,
                                    params.microc_low_sigma,
                                    params.microc_sigma_mult,
                                    params.microc_discr_thresh,
                                    params.modle_low_sigma,
                                    params.modle_sigma_mult,
                                    params.modle_discr_thresh,
                                    file(params.reference_microc),
                                    optimize_extrusion_barriers_001.out.result,
                                    params.num_islands_002,
                                    params.cxpb_002,
                                    params.mutpb_ind_002,
                                    params.mutpb_locus_002,
                                    params.mutsigma_002)

    optimize_extrusion_barriers_003("GRCh38_H1_optimized_barriers_microc_003",
                                    file(params.chrom_sizes),
                                    file(params.candidate_barriers),
                                    file(params.regions_of_interest),
                                    params.bin_size,
                                    params.target_contact_density,
                                    params.microc_low_sigma,
                                    params.microc_sigma_mult,
                                    params.microc_discr_thresh,
                                    params.modle_low_sigma,
                                    params.modle_sigma_mult,
                                    params.modle_discr_thresh,
                                    file(params.reference_microc),
                                    optimize_extrusion_barriers_002.out.result,
                                    params.num_islands_003,
                                    params.cxpb_003,
                                    params.mutpb_ind_003,
                                    params.mutpb_locus_003,
                                    params.mutsigma_003)

    run_modle_sim(file(params.microc_modle_config),
                  file(params.chrom_sizes),
                  optimize_extrusion_barriers_003.out.optimized_barriers,
                  file(params.regions_of_interest))

    cooler_zoomify(run_modle_sim.out.cool)
}

process optimize_extrusion_barriers_001 {
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

        val num_islands
        val cxpb
        val mutpb_individual
        val mutpb_locus
        val mutsigma

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
            --discretization-thresh-tgt '!{discretization_thresh_tgt}' \
            --num-islands=!{num_islands} \
            --cxpb !{cxpb} \
            --mutpb-individual !{mutpb_individual} \
            --mutpb-locus !{mutpb_locus} \
            --early-stopping-pct=0.025 \
            --mut-sigma !{mutsigma} |& tee '!{basename}.log'

        find . -type f -name "*.fifo" -delete
        cp '!{basename}/!{basename}_extrusion_barriers.bed.gz' .
        tar -czf '!{basename}.tar.gz' '!{basename}'
        rm -r '!{basename}'
        '''
}

process optimize_extrusion_barriers_002 {
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
        path previous_simulation_tar

        val num_islands
        val cxpb
        val mutpb_individual
        val mutpb_locus
        val mutsigma

    output:
        path "*_extrusion_barriers.bed.gz", emit: optimized_barriers
        path "*.tar.gz", emit: result
        path "*.log", emit: log

    shell:
        '''
        tar -xvf '!{previous_simulation_tar}' --wildcards \
                                              --no-anchored \
                                              --strip-components=1 \
                                              '*mainland_hall_of_fame.pickle'
        ln -s *.pickle initial_pop.pickle

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
            --discretization-thresh-tgt '!{discretization_thresh_tgt}' \
            --initial-population initial_pop.pickle \
            --num-islands=!{num_islands} \
            --cxpb !{cxpb} \
            --mutpb-individual !{mutpb_individual} \
            --mutpb-locus !{mutpb_locus} \
            --mut-sigma !{mutsigma} |& tee '!{basename}.log'

        find . -type f -name "*.fifo" -delete
        cp '!{basename}/!{basename}_extrusion_barriers.bed.gz' .
        tar -czf '!{basename}.tar.gz' '!{basename}'
        rm -r '!{basename}'
        '''
}

process optimize_extrusion_barriers_003 {
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
        path previous_simulation_tar

        val num_islands
        val cxpb
        val mutpb_individual
        val mutpb_locus
        val mutsigma

    output:
        path "*_extrusion_barriers.bed.gz", emit: optimized_barriers
        path "*.tar.gz", emit: result
        path "*.log", emit: log

    shell:
        '''
        tar -xvf '!{previous_simulation_tar}' --wildcards \
                                              --no-anchored \
                                              --strip-components=1 \
                                              '*mainland_hall_of_fame.pickle'
        ln -s *.pickle initial_pop.pickle

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
            --discretization-thresh-tgt '!{discretization_thresh_tgt}' \
            --initial-population initial_pop.pickle \
            --num-islands=!{num_islands} \
            --cxpb !{cxpb} \
            --mutpb-individual !{mutpb_individual} \
            --mutpb-locus !{mutpb_locus} \
            --mut-sigma !{mutsigma} |& tee '!{basename}.log'

        find . -type f -name "*.fifo" -delete
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
        path "*.log", emit: log

    shell:
        out_prefix=extr_barriers.getSimpleName()
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
