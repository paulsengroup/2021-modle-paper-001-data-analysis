// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#!/usr/bin/env nextflow

nextflow.enable.dsl=2

workflow {
    grch37_bname = "${params.grch37_assembly_name_short}"
    grch38_bname = "${params.grch38_assembly_name_short}"

    bnames = channel.of(grch37_bname, grch38_bname)

    assembly_reports = channel.of(file(params.grch37_assembly_report),
                                  file(params.grch38_assembly_report))

    generate_chrom_sizes(bnames,
                         assembly_reports)

    chrom_sizes = generate_chrom_sizes.out.chrom_sizes.toSortedList(
        { f1, f2 -> file(f1).getBaseName() <=> file(f2).getBaseName() }
    ).flatten()
    chrom_sizes_bed = generate_chrom_sizes.out.bed.toSortedList(
        { f1, f2 -> file(f1).getBaseName() <=> file(f2).getBaseName() }
    ).flatten()

    grch37_chrom_sizes = chrom_sizes.first()
    grch38_chrom_sizes = chrom_sizes.last()

    grch38_chrom_sizes_bed = chrom_sizes_bed.last()

    extract_meme_motif_from_zip(file(params.jaspar_2022_core_zip),
                                params.motif_name)

    motif = extract_meme_motif_from_zip.out.meme
    run_mast("${grch38_bname}_CTCF",
             file(params.grch38_genome_assembly),
             motif)

    convert_mast_to_bed(run_mast.out.txt_gz,
                        grch38_chrom_sizes_bed)

    generate_extr_barriers_bed(file(params.convert_chip_to_occupancy_script),
                               convert_mast_to_bed.out.bed_gz,
                               file(params.h1_rad21_chip_fold_change),
                               file(params.h1_ctcf_chip_peaks),
                               file(params.h1_rad21_chip_peaks),
                               "${grch38_bname}_${params.cell_line_name}_barriers_RAD21_occupancy")

   fix_mcool(file(params.microc_mcool),
             params.microc_base_bin_size)
}

process generate_chrom_sizes {
    publishDir "${params.output_dir}", mode: 'copy'

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
    publishDir "${params.output_dir}", mode: 'copy'

    input:
        val bname
        path genome_assembly_fa
        path motif_meme

    output:
        path "${bname}_mast_hits.txt.gz", emit: txt_gz

    shell:
        out = "${bname}_mast_hits.txt.gz"
        """
        set -e
        set -u
        set -o pipefail

        # MAST does not like gzipped files nor FIFOs, so we have to actually write the inflated reference to disk
        gzip -dc "$genome_assembly_fa" > ga.tmp.fa

        mast -hit_list     \
             "$motif_meme" \
             ga.tmp.fa     |
        gzip -9 > "$out"

        rm -f ga.tmp.fa
        """
}

process convert_mast_to_bed {
    publishDir "${params.output_dir}", mode: 'copy'

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

process extract_meme_motif_from_zip {
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
    publishDir "${params.output_dir}", mode: 'copy'

    label 'process_short'

    input:
        path main_script
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
        python3 '!{main_script}' \
                --chip-seq-bigwig '!{rad21_fold_change_bwig}'       \
                --motifs-bed '!{ctcf_motifs_bed}'                   \
                --regions-of-interest-bed '!{ctcf_chip_peaks_bed}'  \
                                          '!{rad21_chip_peaks_bed}' |
                gzip -9 > '!{bname}.bed.gz'
        '''
}

process convert_hic_to_mcool {
    publishDir "${params.output_dir}", mode: 'copy'

    label 'process_very_long'

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
    label 'process_short'

    input:
        path mcool

    output:
        path "*.mcool.tmp", emit: mcool

    shell:
        '''
        #!/usr/bin/env python3

        import os
        import re
        import shutil
        import sys

        import cooler
        pattern = re.compile(r"^chrom|chr", re.IGNORECASE)

        input_mcool = "!{mcool}"
        output_mcool = os.path.basename(input_mcool.strip()) + ".tmp"
        shutil.copyfile(input_mcool, output_mcool)

        for path in cooler.fileops.list_coolers(input_mcool):
            c = cooler.Cooler(f"{output_mcool}::{path}")
            mappings = {chrom: "chr" + pattern.sub("", chrom, 1) for chrom in c.chromnames}
            cooler.rename_chroms(c, mappings)
        '''
}

process balance_mcool {
    publishDir "${params.output_dir}", mode: 'copy',
                                       saveAs: { fname ->
                                                 file(fname).getBaseName() // Trim .new
                                               }
    label 'process_long'
    label 'process_high'

    memory {
        // 750 MB/core
        750e6 * task.cpus * task.attempt as Long
    }

    input:
        path mcool

    output:
        path "*.mcool.new", emit: mcool

    shell:
        '''
        mapfile -t dsets < \
             <(cooler info '!{mcool}' |&
               grep 'KeyError' |
               grep -o '/resolutions/[[:digit:]]\\+')

        mcool_out="$(basename '!{mcool}' .tmp).new"
        cp '!{mcool}' "$mcool_out"
        for dset in "${dsets[@]}"; do
            cooler balance -p !{task.cpus} "${mcool_out}::$dset"
        done
        '''
}

process fix_mcool {
    publishDir "${params.output_dir}", mode: 'copy'
    label 'process_high'
    label 'process_memory_high'

    input:
        path mcool
        val bin_size

    output:
        path "*_fixed.mcool", emit: mcool

    shell:
        outprefix="${mcool.baseName}"
        '''
        cooler zoomify -p !{task.cpus}                     \
                       -r 5000N                            \
                       --balance                           \
                       --balance-args '-p !{task.cpus}'    \
                       -o '!{outprefix}_fixed.mcool'       \
                       '!{mcool}::/resolutions/!{bin_size}'
        '''

}
