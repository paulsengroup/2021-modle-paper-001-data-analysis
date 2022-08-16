#!/usr/bin/env nextflow
// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

nextflow.enable.dsl=2

workflow {

    basenames = Channel.of("hoxd_dl_wt", "hoxd_pl_wt", "idh_wt")

    chrom_sizes1 = Channel.of(file(params.chrom_sizes_hoxd),
                              file(params.chrom_sizes_hoxd),
                              file(params.chrom_sizes_idh))

    regions_of_interest1 = Channel.of(file(params.regions_of_interest_hoxd),
                                      file(params.regions_of_interest_hoxd),
                                      file(params.regions_of_interest_idh))

    candidate_barriers = Channel.of(file(params.candidate_barriers_hoxd),
                                    file(params.candidate_barriers_hoxd),
                                    file(params.candidate_barriers_idh))

    bin_sizes = Channel.of(params.bin_size_hoxd,
                           params.bin_size_hoxd,
                           params.bin_size_idh)

    target_contact_densities = Channel.of(params.target_contact_density_hoxd,
                                          params.target_contact_density_hoxd,
                                          params.target_contact_density_idh)

    low_sigmas = Channel.of(params.dl_low_sigma,
                            params.pl_low_sigma,
                            params.idh_low_sigma)

    sigma_mults = Channel.of(params.dl_sigma_mult,
                             params.pl_sigma_mult,
                             params.idh_sigma_mult)

    discr_thresholds = Channel.of(params.dl_discr_thresh,
                                  params.pl_discr_thresh,
                                  params.idh_discr_thresh)

    optimization_ref_matrices = Channel.of(file(params.dl_wt_hic),
                                           file(params.pl_wt_hic),
                                           file(params.idh_hic))

    optimize_extrusion_barriers(basenames,
                                chrom_sizes1,
                                candidate_barriers,
                                regions_of_interest1,
                                bin_sizes,
                                target_contact_densities,
                                low_sigmas,
                                sigma_mults,
                                discr_thresholds,
                                params.modle_low_sigma,
                                params.modle_sigma_mult,
                                params.modle_discr_thresh,
                                optimization_ref_matrices,
                                params.annotation_polish_cutoff)

    optimize_extrusion_barriers.out.optimized_barriers
                                   .branch {
                                       dl: it =~ /hoxd_dl/
                                       pl: it =~ /hoxd_pl/
                                       idh: it =~ /idh/
                                       other: true
                                   }.set { optimized_barriers_wt }

    generate_mut_barrier_annotation(Channel.empty()
                                           .concat(optimized_barriers_wt.dl,
                                                   optimized_barriers_wt.dl,
                                                   optimized_barriers_wt.pl,
                                                   optimized_barriers_wt.pl,
                                                   optimized_barriers_wt.idh),
                                    Channel.of(file(params.hoxd_rearr1),
                                               file(params.hoxd_rearr2),
                                               file(params.hoxd_rearr1),
                                               file(params.hoxd_rearr2),
                                               file(params.idh_mut_rearr)))

    extr_barriers = optimize_extrusion_barriers.out.optimized_barriers
                                               .concat(generate_mut_barrier_annotation.out.bed)

    extr_barriers.map {
                        if (it =~ /hoxd/) {
                            return file(params.chrom_sizes_hoxd)
                        }
                        return file(params.chrom_sizes_idh)
                      }
                 .set { chrom_sizes2 }

    extr_barriers.map {
                        if (it =~ /hoxd/) {
                            return file(params.regions_of_interest_hoxd)
                        }
                        return file(params.regions_of_interest_idh)
                      }
                 .set { region_of_interest2 }

    run_modle_sim(file(params.modle_config),
                  chrom_sizes2,
                  extr_barriers,
                  region_of_interest2)

    cooler_zoomify(run_modle_sim.out.cool)

}

process optimize_extrusion_barriers {
    publishDir "${params.output_dir}/optimized_barriers", mode: 'copy'

    label 'process_very_high'
    label 'process_long'

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
        val annotation_polish_cutoff

    output:
        path "*final_annotation.bed.gz", emit: optimized_barriers
        path "optimization_stages/*.tar.gz", emit: result


    shell:
        '''
        '!{params.script_dir}/optimize_extrusion_barriers_wrapper.sh' \
            '!{params.script_dir}' \
            '!{task.cpus}' \
            '!{chrom_sizes}' \
            '!{regions_of_interest}' \
            '!{bin_size}' \
            '!{target_contact_density}' \
            '!{extr_barriers}' \
            '!{reference_matrix}' \
            '!{annotation_polish_cutoff}' \
            '!{gaussian_blur_sigma_ref}' \
            '!{gaussian_blur_mult_ref}' \
            '!{discretization_thresh_ref}' \
            '!{gaussian_blur_sigma_tgt}' \
            '!{gaussian_blur_mult_tgt}' \
            '!{discretization_thresh_tgt}'

        mkdir optimization_stages
        for dir in *stage_???; do
            tar -czf "optimization_stages/!{basename}_$dir.tar.gz" "$dir/"
            rm -r "$dir/"
        done

        mv final_annotation.bed.gz '!{basename}_final_annotation.bed.gz'
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
        
        outname="$(echo '!{outname}' | sed 's/_wt//')"

        '!{params.script_dir}/rearrange_bed_intervals.py' '!{annotation}' '!{rearrangements}' |
            gzip -9 > "$outname"
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
