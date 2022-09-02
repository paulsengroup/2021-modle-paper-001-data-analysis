#!/usr/bin/env nextflow
// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

nextflow.enable.dsl=2

workflow {
    optimization_params = \
        Channel.of(params.hoxd_wt_cell_line,
                   params.idh_wt_cell_line)
               .merge(
        Channel.of(file(params.chrom_sizes_hoxd),
                   file(params.chrom_sizes_idh)),
        Channel.of(file(params.regions_of_interest_hoxd),
                   file(params.regions_of_interest_idh)),
        Channel.of(file(params.candidate_barriers_hoxd),
                   file(params.candidate_barriers_idh)),
        Channel.of(params.bin_size_hoxd,
                   params.bin_size_idh),
        Channel.of(params.target_contact_density_hoxd,
                   params.target_contact_density_idh),
        Channel.of(params.hoxd_low_sigma,
                   params.idh_low_sigma),
        Channel.of(params.hoxd_sigma_mult,
                   params.idh_sigma_mult),
        Channel.of(params.hoxd_discr_thresh,
                   params.idh_discr_thresh),
        Channel.of(file(params.hoxd_wt_hic),
                   file(params.idh_wt_hic)))


    optimize_extrusion_barriers_001(optimization_params,
                                    params.modle_low_sigma,
                                    params.modle_sigma_mult,
                                    params.modle_discr_thresh,
                                    params.cxpb_001,
                                    params.mutpb_ind_001,
                                    params.mutpb_locus_001,
                                    params.mutsigma_001)

    optimize_extrusion_barriers_002(optimization_params,
                                    params.modle_low_sigma,
                                    params.modle_sigma_mult,
                                    params.modle_discr_thresh,
                                    optimize_extrusion_barriers_001.out.result,
                                    params.cxpb_002,
                                    params.mutpb_ind_002,
                                    params.mutpb_locus_002,
                                    params.mutsigma_002)

    optimize_extrusion_barriers_003(optimization_params,
                                    params.modle_low_sigma,
                                    params.modle_sigma_mult,
                                    params.modle_discr_thresh,
                                    optimize_extrusion_barriers_002.out.result,
                                    params.cxpb_003,
                                    params.mutpb_ind_003,
                                    params.mutpb_locus_003,
                                    params.mutsigma_003)


    optimize_extrusion_barriers_003.out.optimized_barriers
                                   .branch {
                                       hoxd: it.contains(params.hoxd_wt_cell_line)
                                       idh: it.contains(params.idh_wt_cell_line)
                                           }
                                   .set { optimized_barriers_wt }

    generate_mut_barrier_annotation(Channel.empty()
                                           .concat(optimized_barriers_wt.hoxd,
                                                   optimized_barriers_wt.hoxd,
                                                   optimized_barriers_wt.idh),
                                    Channel.of(file(params.hoxd_rearr1),
                                               file(params.hoxd_rearr2),
                                               file(params.idh_mut_rearr)))

    extr_barriers = Channel.empty()
                           .mix(optimize_extrusion_barriers_003.out.optimized_barriers,
                                generate_mut_barrier_annotation.out.bed)

    extr_barriers.map {
                        if (it.contains(params.hoxd_wt_cell_line)) {
                            return file(params.chrom_sizes_hoxd)
                        }
                        return file(params.chrom_sizes_idh)
                      }
                 .set { chrom_sizes }

    extr_barriers.map {
                        if (it.contains(params.hoxd_wt_cell_line)) {
                            return file(params.regions_of_interest_hoxd)
                        }
                        return file(params.regions_of_interest_idh)
                      }
                 .set { regions_of_interest }

    run_modle_sim(file(params.modle_config),
                  chrom_sizes,
                  extr_barriers,
                  regions_of_interest)

    cooler_zoomify(run_modle_sim.out.cool)

}

process optimize_extrusion_barriers_001 {
    publishDir "${params.output_dir}/optimized_barriers/001", mode: 'copy'

    label 'process_very_high'
    label 'process_very_long'

    input:
        tuple val(basename),
              path(chrom_sizes),
              path(regions_of_interest),
              path(extr_barriers),
              val(bin_size),
              val(target_contact_density),
              val(gaussian_blur_sigma_ref),
              val(gaussian_blur_mult_ref),
              val(discretization_thresh_ref),
              path(reference_matrix)

        val gaussian_blur_sigma_tgt
        val gaussian_blur_mult_tgt
        val discretization_thresh_tgt

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
            --early-stopping-pct=0.025 \
            --num-islands=1 \
            --cxpb !{cxpb} \
            --mutpb-individual !{mutpb_individual} \
            --mutpb-locus !{mutpb_locus} \
            --mut-sigma !{mutsigma} |& tee '!{basename}.log'

        cp '!{basename}/!{basename}_extrusion_barriers.bed.gz' .
        tar -czf '!{basename}.tar.gz' '!{basename}'
        rm -r '!{basename}'
        '''
}

process optimize_extrusion_barriers_002 {
    publishDir "${params.output_dir}/optimized_barriers/002", mode: 'copy'

    label 'process_very_high'
    label 'process_very_long'

    input:
        tuple val(basename),
              path(chrom_sizes),
              path(regions_of_interest),
              path(extr_barriers),
              val(bin_size),
              val(target_contact_density),
              val(gaussian_blur_sigma_ref),
              val(gaussian_blur_mult_ref),
              val(discretization_thresh_ref),
              path(reference_matrix)

        val gaussian_blur_sigma_tgt
        val gaussian_blur_mult_tgt
        val discretization_thresh_tgt

        path previous_simulation_tar

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
        mkdir previous
        tar -xf '!{previous_simulation_tar}' -C previous --strip-components 1

        barriers='previous/!{basename}_extrusion_barriers.bed.gz'
        initial_pop='previous/!{basename}_mainland_hall_of_fame.pickle'

        '!{params.script_dir}/optimize_extrusion_barrier_params_islands.py' \
            --bin-size '!{bin_size}' \
            --extrusion-barriers "$barriers" \
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
            --initial-population "$initial_pop" \
            --num-islands=1 \
            --cxpb !{cxpb} \
            --mutpb-individual !{mutpb_individual} \
            --mutpb-locus !{mutpb_locus} \
            --mut-sigma !{mutsigma} |& tee '!{basename}.log'

        cp '!{basename}/!{basename}_extrusion_barriers.bed.gz' .
        tar -czf '!{basename}.tar.gz' '!{basename}'
        rm -r '!{basename}'
        '''
}

process optimize_extrusion_barriers_003 {
    publishDir "${params.output_dir}/optimized_barriers/003", mode: 'copy'

    label 'process_very_high'
    label 'process_very_long'

    input:
        tuple val(basename),
              path(chrom_sizes),
              path(regions_of_interest),
              path(extr_barriers),
              val(bin_size),
              val(target_contact_density),
              val(gaussian_blur_sigma_ref),
              val(gaussian_blur_mult_ref),
              val(discretization_thresh_ref),
              path(reference_matrix)

        val gaussian_blur_sigma_tgt
        val gaussian_blur_mult_tgt
        val discretization_thresh_tgt

        path previous_simulation_tar

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
        mkdir previous
        tar -xf '!{previous_simulation_tar}' -C previous --strip-components 1

        barriers='previous/!{basename}_extrusion_barriers.bed.gz'
        initial_pop='previous/!{basename}_mainland_hall_of_fame.pickle'

        '!{params.script_dir}/optimize_extrusion_barrier_params_islands.py' \
            --bin-size '!{bin_size}' \
            --extrusion-barriers "$barriers" \
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
            --initial-population "$initial_pop" \
            --num-islands=1 \
            --cxpb !{cxpb} \
            --mutpb-individual !{mutpb_individual} \
            --mutpb-locus !{mutpb_locus} \
            --mut-sigma !{mutsigma} |& tee '!{basename}.log'

        cp '!{basename}/!{basename}_extrusion_barriers.bed.gz' .
        tar -czf '!{basename}.tar.gz' '!{basename}'
        rm -r '!{basename}'
        '''
}

process generate_mut_barrier_annotation {
    publishDir "${params.output_dir}/optimized_barriers", mode: 'copy'

    label 'process_very_short'

    input:
        path annotation
        path rearrangements

    output:
        path "*.bed.gz", emit: bed

    shell:
        outname="${annotation.simpleName}_${rearrangements.simpleName}.bed.gz"
        '''
        set -o pipefail


        '!{params.script_dir}/rearrange_bed_intervals.py' \
            '!{annotation}'                               \
            '!{rearrangements}'                           |
            gzip -9 > "!{outname}"
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
