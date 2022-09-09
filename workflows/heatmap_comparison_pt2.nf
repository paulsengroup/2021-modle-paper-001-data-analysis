#!/usr/bin/env nextflow
// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

nextflow.enable.dsl=2

workflow {
    lef_params = Channel.of(params.lef_processivity.tokenize(",")).flatten()
        .combine(Channel.of(params.lef_density.tokenize(",")).flatten())
        .combine(Channel.of(params.lef_prob_bypass.tokenize(",")).flatten())

    run_modle_lef_param_exploration(file(params.modle_config_lef_param_exploration),
                                    file(params.chrom_sizes),
                                    file(params.extr_barriers),
                                    params.chroms_of_interest_regex,
                                    lef_params)

    run_modle_contact_prob_decay(Channel.fromPath(file(params.modle_configs_prob_decay)),
                                 file(params.chrom_sizes),
                                 file(params.extr_barriers))

    cooltools_expected_cis(Channel.empty()
                                  .mix(run_modle_contact_prob_decay.out.cool
                                  .map { tuple("", 0, it) },
                           Channel.of(file(params.hic_matrix), file(params.microc_matrix))
                                  .map { tuple("weight", 5000, it)}))

    cooler_zoomify(Channel.empty()
                          .mix(run_modle_lef_param_exploration.out.cool
                          .map { tuple("lef_param_exploration/mcools/", it) },
                               run_modle_contact_prob_decay.out.cool
                          .map { tuple("prob_decay/", it) }))

    matrices_for_plotting = \
        cooltools_expected_cis.out.tsv
                              .branch {
                              hic: it.getBaseName().toLowerCase().contains("_hic_")
                              microc: it.getBaseName().toLowerCase().contains("_microc_")
                              low: it.getBaseName().toLowerCase().contains("_low_proc_tad_plus_loop")
                              normal: it.getBaseName().toLowerCase().contains("_normal_proc_tad_plus_loop")
                              high: it.getBaseName().toLowerCase().contains("_high_proc_tad_plus_loop")
                              }

    plot_prob_decay(matrices_for_plotting.hic,
                    matrices_for_plotting.microc,
                    matrices_for_plotting.low,
                    matrices_for_plotting.normal,
                    matrices_for_plotting.high,
                    params.plot_bin_size,
                    params.plot_diag_width)
}

process run_modle_lef_param_exploration {
    label 'process_medium'

    input:
        path config
        path chrom_sizes
        path extr_barriers
        val chroms_regex
        tuple val(lef_processivity), val(lef_density), val(prob_lef_bypass)

    output:
        path "*.cool", emit: cool

    shell:
        '''
        set -o pipefail

        outprefix="$(basename '!{chrom_sizes}' .chrom.sizes)_"
        outprefix+='!{lef_density}_dens_!{lef_processivity}_proc_!{prob_lef_bypass}_bypass'

        modle sim --config '!{config}'                                   \
                  -b <(zcat '!{extr_barriers}' | grep '!{chroms_regex}') \
                  -c '!{chrom_sizes}'                                    \
                  -t !{task.cpus}                                        \
                  -o "$outprefix"                                        \
                  --lef-density !{lef_density}                           \
                  --avg-lef-processivity !{lef_processivity}             \
                  --probability-of-lef-bypass !{prob_lef_bypass}
        '''
}

process run_modle_contact_prob_decay {
    label 'process_medium'

    input:
        path config
        path chrom_sizes
        path extr_barriers

    output:
        path "*.cool", emit: cool

    shell:
        '''
        modle sim --config '!{config}' \
                  -c '!{chrom_sizes}'  \
                  -b '!{extr_barriers}'
        '''
}

process cooler_zoomify {
    publishDir "${params.outdir}/", mode: 'copy',
                                    saveAs: { "${outdir}/${it}" }

    label 'process_medium'

    input:
        tuple val(outdir), path(cool)

    output:
        path "*.mcool", emit: mcool

    shell:
        '''
        cooler zoomify -p !{task.cpus} -r N '!{cool}'
        '''
}

process cooltools_expected_cis {
    publishDir "${params.outdir}/prob_decay", mode: 'copy'

    input:
        tuple val(weight_dset_name), val(resolution), path(cool)

    output:
        path "*.tsv.gz", emit: tsv

    shell:
        outname="${cool.baseName}_expected_cis.tsv.gz"
        '''
        set -o pipefail

        cooler='!{cool}'

        if [ !{resolution} -ne 0 ]; then
            cooler+='::/resolutions/!{resolution}'
        fi

        cooltools expected-cis --clr-weight-name='!{weight_dset_name}' \
                               "$cooler" |
                  gzip -9 > '!{outname}'
        '''
}

process plot_prob_decay {
    publishDir "${params.outdir}/prob_decay", mode: 'copy'

    input:
        path hic_tsv
        path microc_tsv
        path low_proc_tsv
        path normal_proc_tsv
        path high_proc_tsv
        val bin_size
        val diag_width

    output:
        path "*.png", emit: png
        path "*.svg", emit: svg

    shell:
        '''
        '!{params.script_dir}/plotting/plot_interaction_decay_curve.py' \
            '!{hic_tsv}' '!{microc_tsv}' \
            '!{low_proc_tsv}' '!{normal_proc_tsv}' '!{high_proc_tsv}' \
            --bin-size '!{bin_size}' \
            --diagonal-width '!{diag_width}' \
            --output-prefix "interaction_decay_plot"
        '''
}
