#!/usr/bin/env nextflow
// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

nextflow.enable.dsl=2

workflow {
    run_stripenn(file(params.reference_matrix_file),
                 params.bin_size,
                 "all")

    generate_training_and_test_sites(run_stripenn.out.filtered,
                                     params.seed,
                                     params.eval_sites_fraction_for_training)

    transform_reference_matrix(file(params.reference_matrix_file),
                               params.bin_size,
                               params.diagonal_width,
                               params.gaussian_sigma_ref,
                               params.gaussian_sigma_multiplier_ref,
                               params.discretization_thresh_ref)

    param_files = Channel.of(file(params.param_space_file1),
                             file(params.param_space_file2))
    output_prefixes = Channel.of(file(params.output_prefix1),
                                 file(params.output_prefix2))
    starting_points = Channel.of(params.starting_point1,
                                 params.starting_point2)
    run_gw_optimization(param_files,
                        output_prefixes,
                        file(params.chrom_sizes_file),
                        file(params.extr_barrier_file),
                        generate_training_and_test_sites.out.training_set,
                        generate_training_and_test_sites.out.validation_set,
                        transform_reference_matrix.out.transformed_matrix.first(),
                        params.excluded_chroms,
                        starting_points,
                        params.gaussian_sigma_tgt,
                        params.gaussian_sigma_multiplier_tgt,
                        params.discretization_thresh_tgt,
                        params.diagonal_width,
                        params.num_calls,
                        params.num_random_starts,
                        params.optimization_method,
                        params.seed,
                        params.scoring_method)

    cooler_zoomify_tar(run_gw_optimization.out.tar.flatten())


    generate_modle_configs(file(params.chrom_sizes_file),
                           file(params.extr_barrier_file),
                           run_gw_optimization.out.summary_tsv,
                           file(params.param_intervals_of_interest),
                           params.target_contact_density_param_of_interest)

    run_modle(file(params.chrom_sizes_file),
              file(params.extr_barrier_file),
              generate_modle_configs.out.configs.flatten())

    cooler_zoomify(run_modle.out.cool)

    plot_optimization_results(run_gw_optimization.out.stats_pickle)

}

process generate_training_and_test_sites {
    publishDir "${params.output_dir}", mode: 'copy'

    label 'process_very_short'

    input:
        path input_tsv
        val seed
        val fraction_for_training

    output:
        path "*_training.tsv", emit: training_set
        path "*_validation.tsv", emit: validation_set

    shell:
        '''
        #!/usr/bin/env python3

        from numpy.random import RandomState
        import pandas as pd
        from pathlib import Path

        df = pd.read_csv("!{input_tsv}", sep="\\t") \
               .sort_values(["chr", "pos1", "pos2", "chr2", "pos3", "pos4"])
        bname = Path("!{input_tsv}").stem

        training_set_size = round(df.shape[0] * !{fraction_for_training})

        df1 = df.sample(n=training_set_size, replace=False, random_state=RandomState(!{seed}))
        df2 = df[~df.index.isin(df1.index)]
        assert df1.shape[0] + df2.shape[0] == df.shape[0]

        df1.to_csv(f"{bname}_eval_sites_for_training.tsv", sep="\\t", index=False)
        df2.to_csv(f"{bname}_eval_sites_for_validation.tsv", sep="\\t", index=False)
        '''
}

process transform_reference_matrix {
    publishDir "${params.output_dir}", mode: 'copy'

    label 'process_medium'

    input:
        path reference_matrix
        val bin_size
        val diagonal_width
        val gaussian_blur_sigma
        val gaussian_blur_multiplier
        val discretization_thresh

    output:
        path "${reference_matrix.baseName}_${bin_size}_transformed.cool", emit: transformed_matrix

    shell:
        out = "${reference_matrix.baseName}_${bin_size}_transformed.cool"
        """
        modle_tools transform -i "$reference_matrix"                               \
                              --resolution $bin_size                               \
                              -w $diagonal_width                                   \
                              --method difference_of_gaussians                     \
                              --gaussian-blur-sigma $gaussian_blur_sigma           \
                              --gaussian-blur-multiplier $gaussian_blur_multiplier \
                              --binary-discretization-value $discretization_thresh \
                              -t $task.cpus                                        \
                              -o "$out"
        """
}

process run_gw_optimization {
    publishDir "${params.output_dir}/optimization", mode: 'copy'

    label 'process_very_high'
    label 'process_very_long'

    input:
        path param_space_file
        val output_prefix
        path chrom_sizes
        path extrusion_barriers
        path evaluation_sites_training
        path evaluation_sites_validation
        path reference_matrix
        val excluded_chroms
        val starting_point
        val blur_sigma_tgt
        val blur_sigma_mult_tgt
        val discret_thresh_tgt
        val diagonal_width
        val ncalls
        val num_random_starts
        val optimization_method
        val seed
        val scoring_method

    output:
        path "*.pickle", emit: stats_pickle
        path "*.tsv", emit: summary_tsv
        path "*.tar", emit: tar
        val "${output_prefix.fileName}_${optimization_method}", emit: output_prefix

    shell:
        out="${output_prefix.fileName}"
        '''
        python3 '!{params.script_dir}/optimize_modle_sim_params.py' optimize \
             --param-space-tsv="!{param_space_file}"                         \
             --output-prefix="!{out}"                                        \
             --chrom-sizes="!{chrom_sizes}"                                  \
             --extrusion-barriers="!{extrusion_barriers}"                    \
             --evaluation-sites-training="!{evaluation_sites_training}"      \
             --evaluation-sites-validation="!{evaluation_sites_validation}"  \
             --transformed-reference-matrix="!{reference_matrix}"            \
             --excluded-chroms="!{excluded_chroms}"                          \
             --x0="!{starting_point}"                                        \
             --gaussian-blur-sigma-tgt="!{blur_sigma_tgt}"                   \
             --gaussian-blur-sigma-multiplier-tgt="!{blur_sigma_mult_tgt}"   \
             --discretization-thresh-tgt="!{discret_thresh_tgt}"             \
             --diagonal-width="!{diagonal_width}"                            \
             --ncalls="!{ncalls}"                                            \
             --nrandom-starts="!{num_random_starts}"                         \
             --optimization-method="!{optimization_method}"                  \
             --seed="!{seed}"                                                \
             --threads=!{task.cpus}                                          \
             --modle-tools-eval-metric="!{scoring_method}"
        '''
}

process run_stripenn {
    publishDir "${params.output_dir}/stripenn", mode: 'copy'

    label 'process_very_high'
    label 'process_very_long'

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
        stripenn compute --cool "$input_name"   \
                         --out tmpout/          \
                         --chrom "!{chroms}"    \
                         --norm weight          \
                         --numcores !{task.cpus}

        out_prefix="$(basename "!{cool}")"
        out_prefix="${out_prefix%.*}"

        mv tmpout/result_filtered.tsv "${out_prefix}_filtered.tsv"
        mv tmpout/result_unfiltered.tsv "${out_prefix}_unfiltered.tsv"
        mv tmpout/stripenn.log "${out_prefix}_stripenn.log"
        '''
}


process cooler_zoomify_tar {
    publishDir "${params.output_dir}/optimization", mode: 'copy',
                                                    saveAs: { file(it).getBaseName() } // Trim .new

    label 'process_medium'
    label 'process_short'

    input:
        path tar

    output:
        path "*tar.new", emit: tar

    shell:
        '''
        prefix="$(tar tf '!{tar}' | head -n 1)"
        prefix="${prefix%/}"

        tar -xf '!{tar}'

        matrix="$prefix/$prefix.cool"
        cooler zoomify -p !{task.cpus} -r N -o "$prefix/$prefix.mcool" "$matrix"

        rm "$matrix"

        tar -cf '!{tar}.new' "$prefix/"

        if ! grep -q '\\.mcool$' <(tar tf '!{tar}.new'); then
            echo "Tar archive does not contain any .mcool file"
            exit 1
        fi
        '''
}

process generate_modle_configs {

    label 'process_very_short'

    input:
        path chrom_sizes
        path barriers
        path optimization_report_tsv
        path param_intervals
        val target_contact_density

    output:
        path "*.toml", emit: configs

    shell:
        '''
        #!/usr/bin/env python3

        import pandas as pd

        report = pd.read_table("!{optimization_report_tsv}")
        param_intervals = pd.read_table("!{param_intervals}")
        outprefix = str("!{optimization_report_tsv}").removesuffix(".tsv")

        rows = []


        for i, row in param_intervals.iterrows():
            # Select param combinations falling in an interval of interest
            df = report[(report["--extrusion-barrier-occupancy"].between(row["occ_lb"], row["occ_ub"])) & \\
                        (report["--extrusion-barrier-not-bound-stp"].between(row["puu_lb"], row["puu_ub"]))]

            # Take the param combination leading to the worst score
            rows.append(df.sort_values(by=["training_score", "validation_score"]).tail(1))

        # Add params for the overal best and worst scores
        df = report.sort_values(by=["training_score", "validation_score"], ascending=False)
        rows.append(df.head(1))
        rows.append(df.tail(1))

        df = pd.concat(rows).sort_values(by=["training_score", "validation_score"], ascending=False)

        # Generate 1 MoDLE config for each param combination
        for i, row in df.iterrows():
            occ = row["--extrusion-barrier-occupancy"]
            puu = row["--extrusion-barrier-not-bound-stp"]
            score = row["training_score"]

            config = ["[simulate]"]
            config.append("chrom-sizes=\\"!{chrom_sizes}\\"")
            config.append("extrusion-barrier-file=\\"!{barriers}\\"")
            config.append(f"output-prefix=\\"{outprefix}_{score:.4f}_occ_{occ:.4f}_puu_{puu:.4f}\\"")
            config.append("target-contact-density=!{target_contact_density}")
            config.append(f"extrusion-barrier-occupancy={occ:.4f}")
            config.append(f"extrusion-barrier-not-bound-stp={puu:.4f}")

            with open(f"{outprefix}_{score:.4f}_occ_{occ:.4f}_puu_{puu:.4f}.toml", "w") as f:
                print("\\n".join(config), file=f)
        '''
}

process run_modle {
    label 'process_medium'

    input:
        path chrom_sizes
        path extrusion_barriers
        path config

    output:
        path "*.cool", emit: cool
        path config, emit: config
        path "*.log", emit: log


    shell:
        '''
        modle sim --config '!{config}' -t !{task.cpus}
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


process plot_optimization_results {
    publishDir "${params.output_dir}/optimization", mode: 'copy'

    input:
        path pickle

    output:
        path "*.png", emit: png
        path "*.svg", emit: svg
    shell:
        outprefix="${pickle.baseName}"
        '''
        '!{params.script_dir}/plotting/plot_param_optimization_results.py'  \
             -o '!{outprefix}'    \
             --vmax 65            \
             --gradient-levels 14 \
             '!{pickle}'
        '''
}
