#!/usr/bin/env nextflow
// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

nextflow.enable.dsl=2

workflow {
    optimization_params = \
        Channel.of(file(params.chrom_sizes_hoxd),
                   file(params.chrom_sizes_idh),
                   file(params.chrom_sizes_idh))
        .merge(
        Channel.of(file(params.regions_of_interest_hoxd),
                   file(params.regions_of_interest_idh),
                   file(params.regions_of_interest_idh)),
        Channel.of(file(params.candidate_barriers_grcm38),
                   file(params.candidate_barriers_grch37),
                   file(params.candidate_barriers_grch37)),
        Channel.of(params.bin_size_hoxd,
                   params.bin_size_idh,
                   params.bin_size_idh),
        Channel.of(params.target_contact_density_hoxd,
                   params.target_contact_density_idh,
                   params.target_contact_density_idh),
        Channel.of(params.hoxd_low_sigma,
                   params.gm12878_low_sigma,
                   params.imr90_low_sigma),
        Channel.of(params.hoxd_sigma_mult,
                   params.gm12878_sigma_mult,
                   params.imr90_sigma_mult),
        Channel.of(params.hoxd_discr_thresh,
                   params.gm12878_discr_thresh,
                   params.imr90_discr_thresh),
        Channel.of(file(params.hoxd_wt_hic),
                   file(params.gm12878_hic),
                   file(params.imr90_hic)))

    optimize_extrusion_barriers_001(Channel.of("GRCm38_${params.hoxd_wt_cell_line}_001",
                                               "GRCh37_${params.gm12878_cell_line}_001",
                                               "GRCh37_${params.imr90_cell_line}_001"),
                                    optimization_params,
                                    params.modle_low_sigma,
                                    params.modle_sigma_mult,
                                    params.modle_discr_thresh,
                                    params.num_islands_001,
                                    params.cxpb_001,
                                    params.mutpb_ind_001,
                                    params.mutpb_locus_001,
                                    params.mutsigma_001)

    optimize_extrusion_barriers_002(Channel.of("GRCm38_${params.hoxd_wt_cell_line}_002",
                                               "GRCh37_${params.gm12878_cell_line}_002",
                                               "GRCh37_${params.imr90_cell_line}_002"),
                                    optimization_params,
                                    params.modle_low_sigma,
                                    params.modle_sigma_mult,
                                    params.modle_discr_thresh,
                                    optimize_extrusion_barriers_001.out.result,
                                    params.num_islands_002,
                                    params.cxpb_002,
                                    params.mutpb_ind_002,
                                    params.mutpb_locus_002,
                                    params.mutsigma_002)

    optimize_extrusion_barriers_003(Channel.of("GRCm38_${params.hoxd_wt_cell_line}_003",
                                               "GRCh37_${params.gm12878_cell_line}_003",
                                               "GRCh37_${params.imr90_cell_line}_003"),
                                    optimization_params,
                                    params.modle_low_sigma,
                                    params.modle_sigma_mult,
                                    params.modle_discr_thresh,
                                    optimize_extrusion_barriers_002.out.result,
                                    params.num_islands_003,
                                    params.cxpb_003,
                                    params.mutpb_ind_003,
                                    params.mutpb_locus_003,
                                    params.mutsigma_003)

    optimize_extrusion_barriers_003.out.optimized_barriers
                                   .branch {
                                       hoxd: it.getName().contains(params.hoxd_wt_cell_line)
                                       gm12878: it.getName().contains(params.gm12878_cell_line)
                                       imr90: it.getName().contains(params.imr90_cell_line)
                                           }
                                   .set { optimized_barriers_wt }

    generate_mut_barrier_annotation(Channel.empty()
                                           .concat(Channel.of(
                                                        file(params.barriers_jm8n4),
                                                        file(params.barriers_jm8n4),
                                                        file(params.barriers_gm12878),
                                                        file(params.barriers_imr90)),
                                                   optimized_barriers_wt.hoxd,
                                                   optimized_barriers_wt.hoxd,
                                                   optimized_barriers_wt.gm12878,
                                                   optimized_barriers_wt.imr90),
                                    Channel.of(file(params.hoxd_rearr1),
                                               file(params.hoxd_rearr2),
                                               file(params.idh_mut_rearr),
                                               file(params.idh_mut_rearr),
                                               file(params.hoxd_rearr1),
                                               file(params.hoxd_rearr2),
                                               file(params.idh_mut_rearr),
                                               file(params.idh_mut_rearr)))

    extr_barriers = Channel.empty()
                           .mix(optimize_extrusion_barriers_003.out.optimized_barriers,
                                generate_mut_barrier_annotation.out.bed,
                                Channel.of(file(params.barriers_jm8n4),
                                           file(params.barriers_gm12878),
                                           file(params.barriers_imr90)))


    chrom_sizes = extr_barriers.map {
        if (it.getName().startsWith("GRCm38")) {
            return file(params.chrom_sizes_hoxd)
        }
        assert it.getName().startsWith("GRCh37")
        return file(params.chrom_sizes_idh)
    }

    regions_of_interest = extr_barriers.map {
        if (it.getName().startsWith("GRCm38")) {
            return file(params.regions_of_interest_hoxd)
        }

        assert it.getName().startsWith("GRCh37")
        return file(params.regions_of_interest_idh)
    }

    run_modle_sim(file(params.modle_config),
                  chrom_sizes,
                  extr_barriers,
                  regions_of_interest)

    subsample_contacts(run_modle_sim.out.cool,
                       0,
                       params.contact_subsampling_fract)

    cooler_zoomify(run_modle_sim.out.cool.mix(subsample_contacts.out.cool))

}

process optimize_extrusion_barriers_001 {
    publishDir "${params.output_dir}/optimized_barriers", mode: 'copy'

    label 'process_very_high'
    label 'process_very_long'
    input:
        val basename

        tuple path(chrom_sizes),
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
            --early-stopping-pct=0.025 \
            --num-islands=!{num_islands} \
            --cxpb !{cxpb} \
            --mutpb-individual !{mutpb_individual} \
            --mutpb-locus !{mutpb_locus} \
            --mut-sigma !{mutsigma} |& tee '!{basename}.log'

        find . -name "*.fifo" -delete
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

        tuple path(chrom_sizes),
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

        find . -name "*.fifo" -delete
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

        tuple path(chrom_sizes),
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

        find . -name "*.fifo" -delete
        cp '!{basename}/!{basename}_extrusion_barriers.bed.gz' .
        tar -czf '!{basename}.tar.gz' '!{basename}'
        rm -r '!{basename}'
        '''
}

process generate_mut_barrier_annotation {
    publishDir "${params.output_dir}/mutants", mode: 'copy'

    label 'process_very_short'

    input:
        path annotation
        path rearrangements

    output:
        path "*.bed.gz", emit: bed

    shell:
        '''
        set -o pipefail

        function strip_ext {
            f="$1"

            f="${f%.gz}"
            echo "${f%.bed}"
        }

        annotation_name="$(strip_ext '!{annotation}')"
        rearrangement_name="$(strip_ext '!{rearrangements}')"

        '!{params.script_dir}/rearrange_bed_intervals.py' \
            '!{annotation}'                               \
            '!{rearrangements}'                           |
            gzip -9 > "${annotation_name}_${rearrangement_name}.bed.gz"
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
        '''
        out_prefix="$(echo '!{extr_barriers}' | sed -E 's/_[[:digit:]]{3}_extrusion_barriers//')"

        chroms=($(cut -f 1 '!{regions_of_interest}' | sort -u))

        for chrom in "${chroms[@]}"; do
            zcat '!{extr_barriers}' | grep -P "^${chrom}\\b" >> barriers.bed
        done

        modle sim --config '!{config}'                       \
                  -c '!{chrom_sizes}'                        \
                  -b barriers.bed                            \
                  --chrom-subranges '!{regions_of_interest}' \
                  -t !{task.cpus}                            \
                  -o "${out_prefix%.bed.gz}"
        '''
}

process subsample_contacts {
    input:
        path cool
        val resolution
        val fraction

    output:
        path "*.cool", emit: cool

    shell:
        outname="${cool.baseName}_subsampled.cool"
        '''
        cooler='!{cool}'

        if [ !{resolution} -ne 0 ]; then
            cooler+='::/resolutions/!{resolution}'
        fi

        cooltools random-sample -f !{fraction} "$cooler" '!{outname}'
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
