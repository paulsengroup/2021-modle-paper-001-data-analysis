#!/usr/bin/env nextflow
// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

nextflow.enable.dsl=2

def compare_basename(f1, f2) { file(f1).getBaseName() <=> file(f2).getBaseName() }

workflow {

    genome_assemblies = Channel.fromPath(params.genome_assemblies)
    assembly_reports = Channel.fromPath(params.assembly_reports)

    bnames = Channel.fromPath(params.assembly_reports)
                    .map { it = it.getBaseName(); it.takeWhile { it != '_' } }

    generate_chrom_sizes(bnames,
                         Channel.fromPath(params.assembly_reports))

    chrom_sizes = generate_chrom_sizes.out.chrom_sizes.collect(sort: { f1, f2 -> compare_basename(f1, f2) })
    chrom_sizes_bed = generate_chrom_sizes.out.bed.collect(sort: { f1, f2 -> compare_basename(f1, f2) })

    grch37_chrom_sizes = generate_chrom_sizes.out.chrom_sizes
                                             .filter { it.getBaseName().startsWith(params.grch37_assembly_name_short) }
                                             .first()
    grch38_chrom_sizes = generate_chrom_sizes.out.chrom_sizes
                                             .filter { it.getBaseName().startsWith(params.grch38_assembly_name_short) }
                                             .first()
    grcm38_chrom_sizes = generate_chrom_sizes.out.chrom_sizes
                                             .filter { it.getBaseName().startsWith(params.grcm38_assembly_name_short) }
                                             .first()

    grch38_chrom_sizes_bed = generate_chrom_sizes.out.bed
                                                 .filter { it.getBaseName().startsWith(params.grch38_assembly_name_short) }
                                                 .first()
    grcm38_chrom_sizes_bed = generate_chrom_sizes.out.bed
                                                 .filter { it.getBaseName().startsWith(params.grcm38_assembly_name_short) }
                                                 .first()

    extract_meme_matrix_from_zip(file(params.jaspar_2022_core_zip),
                                 params.motif_id)

    run_mast(params.motif_name,
             genome_assemblies,
             extract_meme_matrix_from_zip.out.meme.first())

    convert_mast_to_bed(run_mast.out.txt_gz.collect(sort: { f1, f2 -> compare_basename(f1, f2) }).flatten(),
                        chrom_sizes_bed.flatten())

    grch38_mast_bed = convert_mast_to_bed.out.bed_gz
                                         .filter { it.getBaseName().startsWith(params.grch38_assembly_name_short) }
                                         .first()

    generate_extr_barriers_bed(grch38_mast_bed,
                               file(params.h1_rad21_chip_fold_change),
                               file(params.h1_ctcf_chip_peaks),
                               file(params.h1_rad21_chip_peaks),
                               "${params.grch38_assembly_name_short}_${params.cell_line_name}_barriers_RAD21_occupancy")

    fixed_mcools = fix_mcool(Channel.fromPath(params.broken_mcools))

    convert_hic_to_mcool(Channel.fromPath(params.hic_files))
    rename_chromosomes_mcool(convert_hic_to_mcool.out.mcool)

    geo_tar_matrix_to_mcool(Channel.fromPath(params.geo_tar_matrices),
                            grch37_chrom_sizes,
                            20000,
                            'chr2')

    grch38_genome_assembly = genome_assemblies.filter { it.getBaseName().startsWith(params.grch38_assembly_name_short) }.first()
    call_compartments(fixed_mcools,
                      grch38_genome_assembly,
                      grch38_chrom_sizes_bed,
                      250000)
}

process generate_chrom_sizes {
    publishDir "${params.output_dir}/chrom_sizes", mode: 'copy'

    label 'process_short'

    input:
        val assembly_name
        path assembly_report

    output:
        path "${assembly_name}.bed", emit: bed
        path "${assembly_name}.chrom.sizes", emit: chrom_sizes

    shell:
        out_bed = "${assembly_name}.bed"
        out = "${assembly_name}.chrom.sizes"
        '''
        set -e
        set -u
        set -o pipefail

         # Extract chromosome sizes from assembly report
         gzip -dc "!{assembly_report}" |
         awk -F $'\t' 'BEGIN { OFS=FS } $2 == "assembled-molecule" { print "chr"$1,0,$9,$7 }' |
         grep -v 'chrMT' | sort -V > "!{out_bed}"

         # Convert chromosome sizes from bed to chrom.sizes format
         cut -f 1,3 "!{out_bed}" | sort -V > "!{out}"
        '''
}

process run_mast {
    publishDir "${params.output_dir}/extrusion_barriers", mode: 'copy'

    input:
        val motif_name
        path genome_assembly_fa
        path motif_meme

    output:
        path "*_mast_hits.txt.gz", emit: txt_gz

    shell:
        '''
        set -u
        set -o pipefail

        outname="$(echo '!{genome_assembly_fa}' | cut -d '_' -f 1)_!{motif_name}_mast_hits.txt.gz"

        trap 'rm -f ga.tmp.fa' EXIT

        # MAST does not like gzipped files nor FIFOs, so we have to actually write the inflated reference to disk
        gzip -dc '!{genome_assembly_fa}' > ga.tmp.fa

        mast -hit_list       \
             '!{motif_meme}' \
             ga.tmp.fa       |
        gzip -9 > "$outname"
        '''
}

process convert_mast_to_bed {
    publishDir "${params.output_dir}/extrusion_barriers", mode: 'copy'

    label 'process_short'

    input:
        path mast_output_txt
        path chrom_sizes_bed

    output:
        path "*_mast_hits.bed.gz", emit: bed_gz

    shell:
        '''
        set -e
        set -u
        set -o pipefail

        outname="$(basename "!{mast_output_txt}" .txt.gz).bed.gz"

        mkfifo tmp.fifo
        # Converts MAST output to BED format
        gzip -dc "!{mast_output_txt}" |
        grep -v '^#'                  |
        awk 'BEGIN { FS="[[:space:]]+"; OFS="\t"; } { print $1,$5,$6,$3";"$4,$7,substr($2, 1, 1) }' \
        > tmp.fifo &

        # - Map chromosome ids to chromosome names (e.g. NC_000001.11 -> chr1)
        awk -F $'\t' 'BEGIN { OFS=FS } NR==FNR{chroms[$4]=$1;next} { if ($1 in chroms) { print chroms[$1],$2,$3,$4,0,$6 } else { print $1,$2,$3,$4,0,$6 } }' \
            "!{chrom_sizes_bed}" \
            tmp.fifo             |
        gzip -9 > "$outname"

        rm -f tmp.fifo
        '''
}

process extract_meme_matrix_from_zip {
    publishDir "${params.output_dir}", mode: 'copy'

    label 'process_short'

    input:
        path zip
        val name

    output:
        path "*.meme", emit: meme

    shell:
        '''
        set -e
        set -u
        set -o pipefail

        unzip "!{zip}" -d tmp/

        mv "tmp/!{name}.meme" .
        rm -rf tmp
        '''
}

process generate_extr_barriers_bed {
    publishDir "${params.output_dir}/extrusion_barriers", mode: 'copy'

    label 'process_short'

    input:
        path ctcf_motifs_bed
        path rad21_fold_change_bwig
        path ctcf_chip_peaks_bed
        path rad21_chip_peaks_bed
        val bname

    output:
        path "${bname}.bed.gz", emit: bed_gz

    shell:
        out="${bname}.bed.gz"
        '''
        '!{params.script_dir}/convert_chip_signal_to_occupancy.py'  \
                --chip-seq-bigwig '!{rad21_fold_change_bwig}'       \
                --motifs-bed '!{ctcf_motifs_bed}'                   \
                --regions-of-interest-bed '!{ctcf_chip_peaks_bed}'  \
                                          '!{rad21_chip_peaks_bed}' |
                gzip -9 > '!{out}'
        '''
}

process convert_hic_to_mcool {
    // publishDir "${params.output_dir}/mcools", mode: 'copy'

    label 'process_very_long'
    label 'process_medium_memory'

    input:
        path hic

    output:
        path "*.mcool", emit: mcool

    shell:
        out = "${hic.baseName}.mcool"
        '''
        set -e
        set -u
        set -o pipefail

        hicConvertFormat -m "!{hic}"       \
                         -o "!{out}"       \
                         --inputFormat hic \
                         --outputFormat cool
        '''
}

process rename_chromosomes_mcool {
    publishDir "${params.output_dir}/mcools", mode: 'copy',
                                       saveAs: { fname ->
                                                 file(fname).getBaseName() // Trim .new
                                               }
    label 'process_short'

    input:
        path mcool

    output:
        path "*.mcool.new", emit: mcool

    shell:
        '''
        '!{params.script_dir}/normalize_cooler_chrom_names.py' '!{mcool}'
        '''
}

process fix_mcool {
    publishDir "${params.output_dir}/mcools", mode: 'copy'
    label 'process_high'
    label 'process_memory_high'

    input:
        path mcool

    output:
        path "*_fixed.mcool", emit: mcool

    shell:
        outprefix="${mcool.baseName}"
        '''
        # Extract mcool datasets. Datasets are sorted by resolution (ascending)
        mapfile -t dsets < \
                <(cooler info '!{mcool}' |&
                  grep 'KeyError' |
                  grep -o '/resolutions/[[:digit:]]\\+')

        cooler zoomify -p !{task.cpus}                     \
                       -r 5000N                            \
                       --balance                           \
                       --balance-args '-p !{task.cpus}'    \
                       -o '!{outprefix}_fixed.mcool'       \
                       "!{mcool}::${dsets[0]}"
        '''

}

process geo_tar_matrix_to_mcool {
    publishDir "${params.output_dir}/mcools", mode: 'copy'
    input:
        path tar
        path chrom_sizes
        val bin_size
        val chrom_name

    output:
        path "*.mcool", emit: mcool

    shell:
        '''
        tar -xf '!{tar}'

        input_matrices=( $(find . -type f -name '*__20kb__raw.matrix.gz') )

        for matrix in "${input_matrices[@]}"; do
            outprefix="${matrix%__20kb__raw.matrix.gz}"
            echo "Processing \\"$outprefix\\"..."

            '!{params.script_dir}/convert_2d_text_matrix_to_cool.py'  \
                    --input-matrices "$matrix"                        \
                    --output-matrices "$outprefix.cool"               \
                    --chrom-sizes '!{chrom_sizes}'                    \
                    --resolution '!{bin_size}'                        \
                    --chrom-names '!{chrom_name}'

            cooler zoomify -p '!{task.cpus}'                           \
                           --balance                                   \
                           --balance-args='-p !{task.cpus} --cis-only' \
                           "$outprefix.cool"
            mv "$outprefix.mcool" .
        done
        '''
}

process call_compartments {
    publishDir "${params.output_dir}/compartments", mode: 'copy'

    input:
        path mcool
        path ref_genome
        path chrom_sizes_bed
        val resolution

    output:
        path "*.vecs.tsv", emit: eigvect_txt
        path "*.lam.txt", emit: eigval_txt
        path "*.bw", emit: bw

    shell:
        '''
        trap 'rm -f ref.fna gc.bed' EXIT

        cooler='!{mcool}::/resolutions/!{resolution}'

        zcat '!{ref_genome}' > ref.fna
        cooler dump -t bins "$cooler" > bins.bed

        mkfifo tmp1.fifo tmp2.fifo

        # - Map chromosome names to chromosome ids (e.g. chr1 -> NC_000001.11)
        awk -F $'\t' 'BEGIN { OFS=FS } NR==FNR{chroms[$1]=$4;next} { print chroms[$1],$2,$3 }' \
            '!{chrom_sizes_bed}'             \
            <(cooler dump -t bins "$cooler") > tmp1.fifo &

        '!{params.script_dir}/compute_binned_gc_content.py' tmp1.fifo \
                                              ref.fna > tmp2.fifo &

        # - Map chromosome ids to chromosome names (e.g. NC_000001.11 -> chr1)
        awk -F $'\t' 'BEGIN { OFS=FS } NR==FNR{chroms[$4]=$1;next} { print chroms[$1],$2,$3,$4,0,$6 }' \
            '!{chrom_sizes_bed}'  \
            tmp2.fifo > gc.bed

        cooltools eigs-cis --phasing-track gc.bed             \
                           --bigwig                           \
                           -o "$(basename '!{mcool}' .mcool)" \
                           "$cooler"
        '''
}
