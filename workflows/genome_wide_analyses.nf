#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { compute_modle_eval_custom_metric } from './modules/modle.nfm' addParams(output_dir: params.output_dir,
                                                               ncpus: params.ncpus)
include { compare_file_names;
          dirname
        } from './modules/utils.nfm'


workflow {
    compute_modle_eval_custom_metric(file(params.modle_all_ctcf_sites_config),
                                     file(params.grch38_chrom_sizes),
                                     file(params.grch38_hela_extr_barriers_all),
                                     params.bin_size,
                                     params.diagonal_width,
                                     params.grch38_hela_reference_matrix,
                                     params.gaussian_blur_sigma_ref,
                                     params.gaussian_blur_sigma_tgt,
                                     params.gaussian_blur_multiplier_ref,
                                     params.gaussian_blur_multiplier_tgt,
                                     params.binary_discretization_value_ref,
                                     params.binary_discretization_value_tgt)

    compute_custom_score_accuracy(compute_modle_eval_custom_metric.out.mt_eval_custom_tsv.flatten(),
                                  file(params.grch38_chrom_sizes),
                                  file(params.grch38_hela_extr_barriers_all))

    select_sites_of_interest(compute_custom_score_accuracy.out.tsv.flatten(),
                             file(params.grch38_hela_extr_barriers_all))
}


process compute_custom_score_accuracy {
    publishDir params.output_dir, mode: 'copy'

    label 'process_short'

    input:
        path mt_eval_tsv
        path chrom_sizes
        path sites_of_interest

    output:
        path "*.bw", emit: bw
        path "*.tsv.gz", emit: tsv

    shell:
        '''
        #!/usr/bin/env python3

        import pandas as pd
        import pyBigWig

        tsv_in = "!{mt_eval_tsv}"
        chrom_sizes = "!{chrom_sizes}"
        assert tsv_in.endswith(".tsv.gz")
        bw = tsv_in[: -len(".tsv.gz")] + "_accuracy.bw"
        tsv_out = tsv_in[: -len(".tsv.gz")] + "_accuracy.tsv.gz"

        df = pd.read_csv(tsv_in, sep="\\t")
        df["accuracy"] = df["correctly_classified_pixels"] / (df["correctly_classified_pixels"] + df["incorrectly_classified_pixels"])

        chroms = [(row[0], row[1]) for _, row in pd.read_csv(chrom_sizes, sep="\\t", header=None).iterrows()]

        with pyBigWig.open(bw, "w") as bwf:
            bwf.addHeader(chroms)
            bwf.addEntries(df["chrom"].tolist(),
                           df["chrom_start"].tolist(),
                           ends=df["chrom_end"].tolist(),
                           values=df["accuracy"].tolist())

        df.to_csv(tsv_out, sep="\\t", na_rep="nan", index=False)
        '''
}

process select_sites_of_interest {
    publishDir params.output_dir, mode: 'copy'

    label 'process_short'

    input:
        path data_tsv
        path sites_of_interest

    output:
        path "*.bed.gz", emit: bed
        path "*.tsv.gz", emit: tsv

    shell:
        baseout= "${data_tsv.simpleName}_sites_of_interest"
        '''
        # Write TSV header
        zcat '!{data_tsv}' | head -n 1 | gzip -9 > '!{baseout}.tsv.gz'

        # Intersect regions of interests with TSV.
        # In case multiple sites of interest overlap on entry in the TSV file,
        # bedtools' output will contain a single entry from the TSV file
        bedtools intersect                         \
                 -wa                               \
                 -u                                \
                 -a <(zcat '!{data_tsv}')          \
                 -b <(zcat '!{sites_of_interest}') |
        gzip -9 >> '!{baseout}.tsv.gz'

        # Convert TSV to BED
        zcat '!{baseout}.tsv.gz' |
            tail -n +1           |
            awk -F '\\t' 'BEGIN{ OFS=FS } NR>1 { print $1,$2,$3,"none",$6,"."; }' |
        gzip -9 > '!{baseout}.bed.gz'
        '''
}
