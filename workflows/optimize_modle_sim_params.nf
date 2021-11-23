#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { run_stripenn } from './modules/stripenn.nfm' addParams(stripenn_output_prefix: "${params.output_dir}/stripenn",
                                                                 ncpus: params.ncpus)

workflow {
    run_stripenn(file(params.reference_matrix_file),
                 params.bin_size,
                 "all")

    generate_training_and_test_sites(run_stripenn.out.stripenn_filtered,
                                     params.seed,
                                     params.eval_sites_fraction_for_training)

    transform_reference_matrix(file(params.reference_matrix_file),
                               params.bin_size,
                               params.diagonal_width,
                               params.gaussian_sigma_ref,
                               params.gaussian_sigma_multiplier_ref,
                               params.discretization_thresh_ref)

    run_optimization(file(params.param_space_file),
                     file(params.output_prefix),
                     file(params.chrom_sizes_file),
                     file(params.extr_barrier_file),
                     generate_training_and_test_sites.out.training_set,
                     transform_reference_matrix.out.transformed_matrix,
                     params.excluded_chroms,
                     params.starting_point,
                     params.gaussian_sigma_ref,
                     params.gaussian_sigma_multiplier_ref,
                     params.discretization_thresh_ref,
                     params.gaussian_sigma_tgt,
                     params.gaussian_sigma_multiplier_tgt,
                     params.discretization_thresh_tgt,
                     params.diagonal_width,
                     params.num_calls,
                     params.num_random_starts,
                     params.optimization_method,
                     params.seed,
                     params.scoring_method)

    test_optimal_params(run_optimization.out.summary_tsv,
                        params.output_prefix,
                        file(params.chrom_sizes_file),
                        file(params.extr_barrier_file),
                        generate_training_and_test_sites.out.testing_set,
                        transform_reference_matrix.out.transformed_matrix,
                        params.excluded_chroms,
                        params.gaussian_sigma_ref,
                        params.gaussian_sigma_multiplier_ref,
                        params.discretization_thresh_ref,
                        params.gaussian_sigma_tgt,
                        params.gaussian_sigma_multiplier_tgt,
                        params.discretization_thresh_tgt,
                        params.diagonal_width,
                        params.scoring_method,
                        params.num_params_to_test)
}

process generate_training_and_test_sites {
    publishDir "${params.output_dir}", mode: 'copy'

    input:
        path input_tsv
        val seed
        val fraction_for_training

    output:
        path "*_training.tsv", emit: training_set
        path "*_testing.tsv", emit: testing_set

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
        df2.to_csv(f"{bname}_eval_sites_for_testing.tsv", sep="\\t", index=False)
        '''
}

process transform_reference_matrix {
    publishDir "${params.output_dir}", mode: 'copy'

    cpus params.ncpus

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
                              --bin-size $bin_size                                 \
                              -w $diagonal_width                                   \
                              --method difference_of_gaussians                     \
                              --gaussian-blur-sigma $gaussian_blur_sigma           \
                              --gaussian-blur-multiplier $gaussian_blur_multiplier \
                              --binary-discretization-value $discretization_thresh \
                              -t $task.cpus                                        \
                              -o "$out"
        """
}

process run_optimization {
    publishDir "${params.output_dir}/optimization", mode: 'copy'

    cpus params.ncpus

    input:
        path param_space_file
        val output_prefix
        path chrom_sizes
        path extrusion_barriers
        path evaluation_sites
        path reference_matrix
        val excluded_chroms
        val starting_point
        val blur_sigma_ref
        val blur_sigma_mult_ref
        val discret_thresh_ref
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
        path "${output_prefix.fileName}_${optimization_method}.pickle", emit: stats_pickle
        path "${output_prefix.fileName}_${optimization_method}.tsv", emit: summary_tsv

    shell:
        out="${output_prefix.fileName}"
        '''
        optimize_modle_sim_params.py optimize \
                                     --param-space-tsv="!{param_space_file}"                       \
                                     --output-prefix="!{out}"                            \
                                     --chrom-sizes="!{chrom_sizes}"                                \
                                     --extrusion-barriers="!{extrusion_barriers}"                  \
                                     --evaluation-sites="!{evaluation_sites}"                      \
                                     --transformed-reference-matrix="!{reference_matrix}"                      \
                                     --excluded-chroms="!{excluded_chroms}"                        \
                                     --x0="!{starting_point}"                                      \
                                     --gaussian-blur-sigma-ref="!{blur_sigma_ref}"                 \
                                     --gaussian-blur-sigma-multiplier-ref="!{blur_sigma_mult_ref}" \
                                     --discretization-thresh-ref="!{discret_thresh_ref}"           \
                                     --gaussian-blur-sigma-tgt="!{blur_sigma_tgt}"                 \
                                     --gaussian-blur-sigma-multiplier-tgt="!{blur_sigma_mult_tgt}" \
                                     --discretization-thresh-tgt="!{discret_thresh_tgt}"           \
                                     --diagonal-width="!{diagonal_width}"                          \
                                     --ncalls="!{ncalls}"                                          \
                                     --nrandom-starts="!{num_random_starts}"                       \
                                     --optimization-method="!{optimization_method}"                \
                                     --seed="!{seed}"                                              \
                                     --nthreads=!{task.cpus}                                       \
                                     --modle-tools-eval-metric="!{scoring_method}"
        '''
}

process test_optimal_params {
    publishDir "${params.output_dir}/test", mode: 'copy'

    cpus params.ncpus

    input:
        path optimization_result
        val output_prefix
        path chrom_sizes
        path extrusion_barriers
        path evaluation_sites
        path reference_matrix
        val excluded_chroms
        val blur_sigma_ref
        val blur_sigma_mult_ref
        val discret_thresh_ref
        val blur_sigma_tgt
        val blur_sigma_mult_tgt
        val discret_thresh_tgt
        val diagonal_width
        val scoring_method
        val num_params_to_test

    output:
        path "*.cool", emit: cool
        path "*.bw", emit: bigwig
        path "*.tsv", emit: report

    shell:
        '''
        optimize_modle_sim_params.py test \
                                     --optimization-result-tsv="!{optimization_result}"            \
                                     --num-params-to-test=!{num_params_to_test}                    \
                                     --output-prefix="!{output_prefix}"                            \
                                     --chrom-sizes="!{chrom_sizes}"                                \
                                     --extrusion-barriers="!{extrusion_barriers}"                  \
                                     --evaluation-sites="!{evaluation_sites}"                      \
                                     --transformed-reference-matrix="!{reference_matrix}"                      \
                                     --excluded-chroms="!{excluded_chroms}"                        \
                                     --gaussian-blur-sigma-ref="!{blur_sigma_ref}"                 \
                                     --gaussian-blur-sigma-multiplier-ref="!{blur_sigma_mult_ref}" \
                                     --discretization-thresh-ref="!{discret_thresh_ref}"           \
                                     --gaussian-blur-sigma-tgt="!{blur_sigma_tgt}"                 \
                                     --gaussian-blur-sigma-multiplier-tgt="!{blur_sigma_mult_tgt}" \
                                     --discretization-thresh-tgt="!{discret_thresh_tgt}"           \
                                     --diagonal-width="!{diagonal_width}"                          \
                                     --nthreads=!{task.cpus}                                       \
                                     --modle-tools-eval-metric="!{scoring_method}"
        '''
}
