#!/usr/bin/env nextflow

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

        run_openmm(file(params.run_openmm_script),
                   chrom_ranges.name,
                   chrom_ranges.size,
                   chrom_ranges.start,
                   file(params.extr_barriers),
                   params.monomer_size,
                   params.md_lef_processivity,
                   params.md_lef_separation)

        openmm_to_cool(file(params.openmm_to_cool_script),
                       run_openmm.out.tar.collect(),
                       file(params.chrom_sizes),
                       run_openmm.out.chrom_name.collect(),
                       run_openmm.out.offset.collect(),
                       params.monomer_size)

        modle_configs = Channel.of(file(params.modle_config_tad_plus_loop),
                                   file(params.modle_config_tad_only),
                                   file(params.modle_config_loop_only),
                                   file(params.modle_config_ctcf_kd),
                                   file(params.modle_config_wapl_kd))
        run_modle(modle_configs,
                  file(params.chrom_sizes),
                  file(params.extr_barriers))

        subsample_contacts(run_modle.out.cool,
                           0.15)

        coolers = subsample_contacts.out.cool.concat(openmm_to_cool.out.cool)
        cooler_zoomify(coolers)

        // plot_comparison(file(params.plot_comparison_script),
        //                 file(params.hic_cool),
        //                 openmm_to_cool.out.cool,
        //                 run_modle.out.cool,
        //                 file(params.microc_cool),
        //                 file(params.regions_of_interest),
        //                 params.chrom_ranges_microc_cmp,
        //                 params.bin_size,
        //                 params.assembly_name,
        //                 params.color_maps,
        //                 params.color_scale_intervals)
}

process liftover_cool {
    label 'process_short'

    input:
        path script
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
        ./liftover_matrix_cool.py --input-cool '!{input_cool}::/resolutions/!{bin_size}'   \
                                  --output-cool "$outname"                                 \
                                  --chrom-sizes '!{chrom_sizes}'                           \
                                  --chrom-ranges '!{chrom_ranges}'                         \
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
        sed 's|^threads=.*|threads=!{task.cpus}|' > out/config.tmp.toml

        modle --config out/config.tmp.toml
        '''
}

process run_openmm {
    cpus 1
    label 'process_long'

    input:
        path main_script
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

        python '!{main_script}' --gpu                 \
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
        path main_script
        path openmm_tars
        path chrom_sizes
        val chrom_names
        val offsets
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

        python '!{main_script}'                 \
            --input-folders ${input_folders[@]} \
            --output-name "$out_name"           \
            --chrom-sizes '!{chrom_sizes}'      \
            --bin-size !{bin_size}              \
            --chrom-names ${chrom_names[@]}     \
            --threads !{task.cpus}              \
            --offsets ${offsets[@]}

        rm -rf *_hdf5
        '''
}

process plot_comparison {
    publishDir params.plot_outdir, mode: 'copy',
                                   saveAs: { "GRCh38_H1_heatmaps_comparison.svg" }
    label 'process_very_short'

    input:
        path script
        path hic_cool
        path md_cool
        path modle_cool
        path microc_cool
        path chrom_ranges_hic_cmp
        val chrom_ranges_microc_cmp
        val bin_size
        val assembly_name
        val color_maps
        val color_scale_intervals

    output:
        path "*.svg", emit: plot

    shell:
        '''

        make_cool_uri () {
            path="$1"
            res="$2"

            if [[ $path == *.mcool ]]; then
                echo "$path::/resolutions/$res"
            else
                echo "$path"
            fi
        }

        hic_cool="$(make_cool_uri '!{hic_cool}' !{bin_size})"
        md_cool="$(make_cool_uri '!{md_cool}' !{bin_size})"
        modle_cool="$(make_cool_uri '!{modle_cool}' !{bin_size})"
        microc_cool="$(make_cool_uri '!{microc_cool}' !{bin_size})"

        ./make_heatmap_comparison_plot.py                                        \
            --path-to-hic-cool "$hic_cool"                                       \
            --path-to-md-cool "$md_cool"                                         \
            --path-to-modle-cool "$modle_cool"                                   \
            --path-to-microc-cool "$microc_cool"                                 \
            --color-maps '!{color_maps}'                                         \
            --color-scale-intervals '!{color_scale_intervals}'                   \
            --path-to-chrom-subranges-hic-comparison '!{chrom_ranges_hic_cmp}'   \
            --chrom-subrange-microc-comparison-ucsc '!{chrom_ranges_microc_cmp}' \
            --output-name '!{assembly_name}_comparison.svg'
        '''
}

process subsample_contacts {
    publishDir params.matrix_outdir, mode: 'copy'
    label 'process_short'

    input:
        path cool
        val fraction

    output:
        path "*.cool", emit: cool

    shell:
        '''
        out="$(basename '!{cool}' .cool)_subsampled.cool"
        cooltools random-sample --frac !{fraction} \
                                -p !{task.cpus}    \
                                '!{cool}' "$out"
        '''
}

process cooler_zoomify {
    publishDir params.matrix_outdir, mode: 'copy'

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
