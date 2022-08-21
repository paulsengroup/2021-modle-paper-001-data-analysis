#!/usr/bin/env nextflow
// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

nextflow.enable.dsl=2

workflow {
    main:
        Channel.fromPath(params.regions_of_interest)
                .splitText()
                .multiMap {
                    // File is in BED3 format
                    toks = it.trim().split('\t')

                    assert toks.size() >= 3
                    assert Long.parseLong(toks[2]) >= Long.parseLong(toks[1])

                    name: toks[0]
                    start: Long.parseLong(toks[1])
                    end: Long.parseLong(toks[2])
                    size: ((Float.parseFloat(toks[2]) - Float.parseFloat(toks[1])) / 1.0e6).round(5)
                    }
                .set { chrom_ranges }

        run_openmm(chrom_ranges.name,
                   chrom_ranges.size,
                   chrom_ranges.start,
                   file(params.extr_barriers),
                   params.monomer_size,
                   params.md_lef_processivity,
                   params.md_lef_separation)

        openmm_to_cool(run_openmm.out.tar.collect(),
                       file(params.chrom_sizes),
                       run_openmm.out.chrom_name.collect(),
                       run_openmm.out.offset.collect(),
                       params.monomer_size,
                       params.bin_size)

        run_modle(Channel.fromPath(params.modle_configs),
                                   file(params.chrom_sizes),
                                   file(params.extr_barriers))

        coolers0 = run_modle.out.cool.mix(openmm_to_cool.out.cool)
        subsample_contacts(file(params.microc_cool),
                           coolers0,
                           file(params.regions_of_interest),
                           params.bin_size,
                           params.contact_subsampling_diagonal_width)

        coolers1 = subsample_contacts.out.cool.mix(run_modle.out.cool,
                                                   openmm_to_cool.out.cool)
        cooler_zoomify(coolers1)

        chrom_ranges.name.reduce("") { it, buff ->
                                       if (buff == "") {
                                           return it;
                                       }
                                           return "${buff},${it}";
                                      }
                         .set { chroms_str }
        chroms_str = chroms_str - ~/,$/  // Strip trailing comma
        run_stripenn(file(params.microc_cool),
                     params.bin_size,
                     chroms_str)

        cooler_zoomify.out.mcool.mix(Channel.of(file(params.microc_cool))).branch {
                                         microc: it.getBaseName().startsWith(file(params.microc_cool).getBaseName())
                                         openmm: it.getBaseName().contains("openmm") && it.getBaseName().contains("subsampled")
                                         modle: it.getBaseName().contains("tad_plus_loop") && it.getBaseName().contains("subsampled")
                                         others: true
                                        }
                                .set { coolers2 }

        coolers3 = Channel.empty().concat(coolers2.microc,
                                          coolers2.openmm,
                                          coolers2.modle)
        mt_transform_cutoffs = Channel.of(params.diff_of_gaussians_cutoff_microc,
                                          params.diff_of_gaussians_cutoff_openmm,
                                          params.diff_of_gaussians_cutoff_modle)

        compute_diff_of_gaussian(coolers3,
                                 params.bin_size,
                                 mt_transform_cutoffs,
                                 params.diff_of_gaussians_diagonal_width)

        compute_diff_of_gaussian.out.transformed_cool.branch {
                                                              microc: it.getBaseName().startsWith(file(params.microc_cool).getBaseName())
                                                              openmm: it.getBaseName().contains("openmm")
                                                              modle: it.getBaseName().contains("tad_plus_loop")
                                                              others: true
                                                             }
                                                     .set { coolers4 }

        reference_coolers =
        target_coolers = Channel.empty().concat(coolers4.openmm,
                                                coolers4.modle,
                                                coolers4.openmm)
        run_modle_tools_eval(Channel.empty().concat(coolers4.openmm,
                                                    coolers4.modle,
                                                    coolers4.openmm),
                             Channel.empty().concat(coolers4.microc,
                                                    coolers4.microc,
                                                    coolers4.modle),
                             params.bin_size,
                             params.diff_of_gaussians_diagonal_width)

        filter_modle_tools_eval(run_modle_tools_eval.out.vertical_tsv_gz,
                                run_modle_tools_eval.out.horizontal_tsv_gz,
                                file(params.regions_of_interest),
                                run_stripenn.out.filtered)

        dump_pixels(Channel.empty().concat(coolers2.openmm,
                                           coolers2.modle,
                                           coolers2.openmm),
                    Channel.empty().concat(coolers2.microc,
                                           coolers2.microc,
                                           coolers2.modle),
                    file(params.regions_of_interest),
                    params.bin_size,
                    params.diagonal_width)

        compute_custom_scores_py(Channel.empty().concat(coolers2.openmm,
                                                        coolers2.modle,
                                                        coolers2.openmm),
                                 Channel.empty().concat(coolers2.microc,
                                                        coolers2.microc,
                                                        coolers2.modle),
                                 file(params.regions_of_interest),
                                 params.bin_size,
                                 params.diagonal_width)

        filter_custom_scores(Channel.empty().mix(compute_custom_scores_py.out.vertical_scores,
                                                 compute_custom_scores_py.out.horizontal_scores),
                             file(params.extr_barriers))
}

process liftover_cool {
    label 'process_short'

    input:
        path input_cool
        path chrom_sizes
        path chrom_ranges
        val bin_size
        val assembly_name

    output:
        path "*.cool", emit: cool

    shell:
        '''
        outname="!{assembly_name}_$(basename '!{input_cool}' .mcool)_!{bin_size}.cool"
        '!{params.script_dir}/liftover_matrix_cool.py'                                   \
                                  --input-cool '!{input_cool}::/resolutions/!{bin_size}' \
                                  --output-cool "$outname"                               \
                                  --chrom-sizes '!{chrom_sizes}'                         \
                                  --chrom-ranges '!{chrom_ranges}'                       \
                                  --assembly-name '!{assembly_name}'
        '''
}

process run_modle {
    publishDir params.matrix_outdir, mode: 'copy'
    label 'process_high'

    input:
        path config
        path chrom_sizes
        path extr_barriers

    output:
        path "*.cool", emit: cool

    shell:
        '''
        mkdir out

        sed 's|^chrom-sizes=.*|chrom-sizes="!{chrom_sizes}"|' "!{config}" |
        sed 's|^extrusion-barrier-file=.*|extrusion-barrier-file="!{extr_barriers}"|' |
        sed "s|^threads=.*|threads=$(( 2 * !{task.cpus} ))|" > out/config.tmp.toml

        modle --config out/config.tmp.toml
        '''
}

process run_openmm {
    cpus 1
    label 'process_very_long'

    input:
        val chrom_name
        val chrom_size
        val offset
        path extr_barriers
        val monomer_size
        val lef_processivity
        val lef_separation

    output:
        path "*.tar.gz", emit: tar
        val chrom_name, emit: chrom_name
        val chrom_size, emit: chrom_size
        val offset, emit: offset

    shell:
        '''
        mkdir tmp
        awk -F '\\t' -v offset=!{offset}   \
                     -v name=!{chrom_name} \
            'BEGIN{ OFS=FS } { if ($1==name && $2>=offset) { print $1,$2-offset,$3-offset,$4,$5,$6 } }' \
            <(zcat !{extr_barriers}) > tmp/barriers.bed

        '!{params.script_dir}/openmm/polymer_simulation_diffusive_mixed.py' \
            --gpu                                     \
            --output-folder !{chrom_name}             \
            --extrusion-barrier-bed tmp/barriers.bed  \
            --monomer-size-bp !{monomer_size}         \
            --lef-processivity-bp !{lef_processivity} \
            --lef-separation-bp !{lef_separation}     \
            --simulation-size-mbp !{chrom_size}

        tar -czf '!{chrom_name}.tar.gz' '!{chrom_name}'
        rm -rf '!{chrom_name}/' tmp/
        '''
}

process openmm_to_cool {
    publishDir params.matrix_outdir, mode: 'copy'

    label 'process_medium'

    input:
        path openmm_tars
        path chrom_sizes
        val chrom_names
        val offsets
        val monomer_size
        val bin_size

    output:
        path "*.cool", emit: cool

    shell:
        '''
        set -x

        input_folders=()
        for f in !{openmm_tars}; do
            tar -xf "$f"

            input_dir="$(basename "$f" .tar.gz)"
            output_dir="${input_dir}_hdf5"
            polychrom_traj_convert --allow-nonconsecutive "$input_dir" "$output_dir"
            input_folders+=("$output_dir")
        done

        chrom_names=( $(echo '!{chrom_names}' | tr '[],' ' ') )
        offsets=( $(echo '!{offsets}' | tr '[],' ' ') )
        out_name='GRCh38_H1_openmm_heatmaps_comparison.cool'

        '!{params.script_dir}/openmm/generate_contact_matrix.py' \
            --input-folders ${input_folders[@]}                  \
            --output-name "${out_name}.tmp"                      \
            --chrom-sizes '!{chrom_sizes}'                       \
            --bin-size !{monomer_size}                           \
            --chrom-names ${chrom_names[@]}                      \
            --threads $(( 2 * !{task.cpus} ))                    \
            --offsets ${offsets[@]}

        cooler zoomify -r !{bin_size}              \
                       -p $(( 2 * !{task.cpus} ))  \
                       -o "${out_name}.mcool"      \
                       "${out_name}.tmp"

        cooler cp "${out_name}.mcool::/resolutions/!{bin_size}" "${out_name}"

        rm -rf *.tmp *_hdf5 *.mcool
        '''
}

process subsample_contacts {
    publishDir params.matrix_outdir, mode: 'copy'

    label 'process_medium'

    input:
        path reference
        path target
        path chrom_ranges
        val bin_size
        val diagonal_width

    output:
        path "*.cool", emit: cool

    shell:
        '''
        out="$(basename '!{target}' .mcool)"
        out="$(basename "$out" .cool)_subsampled.cool"

        '!{params.script_dir}/subsample_contact_matrix.py' \
                --ref-matrix '!{reference}'                \
                --tgt-matrix '!{target}'                   \
                --bin-size '!{bin_size}'                   \
                --diagonal-width '!{diagonal_width}'       \
                --output "$out"                            \
                --chrom-ranges-bed '!{chrom_ranges}'       \
                --threads $(( 2 * !{task.cpus} ))
        '''
}

process cooler_zoomify {
    publishDir params.matrix_outdir, mode: 'copy'
    memory '16.G'

    input:
        path cool

    output:
        path "*.mcool", emit: mcool

    shell:
        '''
        cooler zoomify -r 5000N        \
                       -p !{task.cpus} \
                       '!{cool}'
        '''
}

process compute_diff_of_gaussian {
    publishDir params.matrix_outdir, mode: 'copy'

    label 'process_medium'

    input:
        path cool
        val bin_size
        val cutoff
        val diagonal_width

    output:
        path "*.cool", emit: transformed_cool

    shell:
        '''
        outname='!{cool}'
        outname="${outname%.mcool}"
        outname="${outname%.cool}_transformed.cool"

        modle_tools transform                               \
                    -i '!{cool}'                            \
                    -o "$outname"                           \
                    --resolution '!{bin_size}'              \
                    -m difference_of_gaussians              \
                    --binary-discretization-value !{cutoff} \
                    -w !{diagonal_width}                    \
                    --threads $(( 2 * !{task.cpus} ))
        '''
}

process run_modle_tools_eval {
    publishDir "${params.outdir}/eval", mode: 'copy'

    label 'process_medium'

    when:
        target_cool != reference_cool

    input:
        path target_cool
        path reference_cool
        val bin_size
        val diagonal_width

    output:
        path "*horizontal.bw", emit: horizontal_bwig
        path "*vertical.bw", emit: vertical_bwig
        path "*horizontal.tsv.gz", emit: horizontal_tsv_gz
        path "*vertical.tsv.gz", emit: vertical_tsv_gz

    shell:
        '''
        prefix1='!{target_cool}'
        prefix1="${prefix1%.mcool}"
        prefix1="${prefix1%.cool}"
        prefix2='!{reference_cool}'
        prefix2="${prefix2%.mcool}"
        prefix2="${prefix2%.cool}"

        out_prefix="${prefix1}_${prefix2}"

        modle_tools eval                                   \
                    -i '!{target_cool}'                    \
                    --reference-matrix '!{reference_cool}' \
                    -o "$out_prefix"                       \
                    --resolution '!{bin_size}'             \
                    -w !{diagonal_width}                   \
                    --threads $(( 2 * !{task.cpus} ))
        '''
}

process filter_modle_tools_eval {
    publishDir "${params.outdir}/eval", mode: 'copy'

    label 'process_medium'

    input:
        path vertical_scores_tsv
        path horizontal_scores_tsv
        path regions_of_interest
        path stripes_bed

    output:
        path "*vertical*.tsv.gz", emit: vertical_tsv_gz
        path "*horizontal*.tsv.gz", emit: horizontal_tsv_gz

    shell:
        '''
        bedtools intersect -a <(tail -n +2 '!{stripes_bed}') \
                           -b '!{regions_of_interest}'       \
           > regions_of_interest.bed.tmp

        outname_vertical="$(basename '!{vertical_scores_tsv}' .tsv.gz)_filtered.tsv.gz"
        outname_horizontal="$(basename '!{horizontal_scores_tsv}' .tsv.gz)_filtered.tsv.gz"

        header="$(head -n 1 '!{vertical_scores_tsv}')"


        bedtools intersect -a <(zcat '!{vertical_scores_tsv}') \
                           -b regions_of_interest.bed.tmp > bed1.tmp
        cat <(echo $header) bed1.tmp | gzip -9 > "$outname_vertical"

        bedtools intersect -a <(zcat '!{horizontal_scores_tsv}') \
                           -b regions_of_interest.bed.tmp > bed2.tmp

        cat <(echo $header) bed2.tmp | gzip -9 > "$outname_horizontal"

        rm *.tmp
        '''
}

process run_stripenn {
    publishDir "${params.outdir}/stripenn", mode: 'copy'

    label 'process_medium'
    label 'process_long'

    input:
        path cool
        val resolution
        val chroms

    output:
        path "*_filtered.tsv", emit: filtered
        path "*_unfiltered.tsv", emit: unfiltered
        path "*_stripenn.log", emit: log

    shell:
        '''
        input_name="!{cool}"
        if [[ $input_name == *.mcool ]]; then
            input_name="${input_name}::/resolutions/!{resolution}"
        fi

        # Apparently --norm and --chrom are required in order to get stripenn
        # to work reliably
        stripenn compute --cool "$input_name"               \
                         --out tmpout/                      \
                         --chrom "!{chroms}"                \
                         --norm weight                      \
                         --numcores $(( 2 * !{task.cpus} ))

        out_prefix="$(basename "!{cool}")"
        out_prefix="${out_prefix%.*}"

        mv tmpout/result_filtered.tsv "${out_prefix}_filtered.tsv"
        mv tmpout/result_unfiltered.tsv "${out_prefix}_unfiltered.tsv"
        mv tmpout/stripenn.log "${out_prefix}_stripenn.log"
        '''
}

process dump_pixels {
    publishDir "${params.outdir}/dumps", mode: 'copy'
    label 'process_short'

    input:
        path ref_cool
        path tgt_cool
        path regions_of_interest
        val resolution
        val diagonal_width

    output:
        path "*.tsv.gz", emit: pixels

    shell:
        out="${ref_cool.simpleName}_vs_${tgt_cool.simpleName}.tsv.gz"
        '''
        set -o pipefail

        '!{params.script_dir}/extract_pairs_of_pixels_from_cool.py' \
            --ref-matrix '!{ref_cool}'                              \
            --tgt-matrix '!{tgt_cool}'                              \
            --chrom-ranges-bed '!{regions_of_interest}'             \
            --bin-size '!{resolution}'                              \
            --diagonal-width '!{diagonal_width}'                    |
        gzip -9c > '!{out}'

        '''

}

process compute_custom_scores_py {
    publishDir "${params.outdir}/eval", mode: 'copy'
    label 'process_short'

    input:
        path ref_cool
        path tgt_cool
        path regions_of_interest
        val resolution
        val diagonal_width

    output:
        path "*_vertical.bedpe.gz", emit: vertical_scores
        path "*_horizontal.bedpe.gz", emit: horizontal_scores

    shell:
        outprefix="${ref_cool.simpleName}_vs_${tgt_cool.simpleName}"
        '''
        set -o pipefail

        out_vertical="$(echo '!{outprefix}' | sed 's/_transformed//g')_vertical.bedpe.gz"
        out_horizontal="$(echo '!{outprefix}' | sed 's/_transformed//g')_horizontal.bedpe.gz"

        '!{params.script_dir}/compute_custom_score.py'  \
            --ref-matrix '!{ref_cool}'                  \
            --tgt-matrix '!{tgt_cool}'                  \
            --chrom-ranges-bed '!{regions_of_interest}' \
            --bin-size '!{resolution}'                  \
            --direction 'horizontal'                    \
            --diagonal-width '!{diagonal_width}'        |
        gzip -9c > "$out_horizontal"

        '!{params.script_dir}/compute_custom_score.py'  \
            --ref-matrix '!{ref_cool}'                  \
            --tgt-matrix '!{tgt_cool}'                  \
            --chrom-ranges-bed '!{regions_of_interest}' \
            --bin-size '!{resolution}'                  \
            --direction 'vertical'                      \
            --diagonal-width '!{diagonal_width}'        |
        gzip -9c > "$out_vertical"
        '''
}

process filter_custom_scores {
    publishDir "${params.outdir}/eval", mode: 'copy'
    label 'process_short'

    input:
        path input_scores
        path filtering_bed

    output:
        path "*_filtered.bedpe.gz", emit: scores

    shell:
        out="${input_scores.simpleName}_filtered.bedpe.gz"
        '''
        if [[ '!{input_scores}' == *_vertical.bedpe.* ]]; then
            pattern='[[:space:]]-[[:space:]]*$'
        else
            pattern='[[:space:]]+[[:space:]]*$'
        fi

        header="$(zcat '!{input_scores}' | head -n 1 | tr -d '\\n')"
        printf "%s\\tbarrier_strength\\tbarrier_direction\\n" "$header" > '!{out}.tmp'

        bedtools intersect -a <(zcat '!{input_scores}'| tail -n +2 )       \
                           -b <(zcat '!{filtering_bed}' | grep "$pattern") \
                           -loj              |
                           grep -v '\\-1'    |
                           cut -f 1-11,16-17 >> '!{out}.tmp'

        gzip -9 < '!{out}.tmp' > '!{out}'
        rm '!{out}.tmp'
        '''
}
