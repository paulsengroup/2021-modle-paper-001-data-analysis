#!/usr/bin/env nextflow

nextflow.enable.dsl=2

avail_cpus = Runtime.runtime.availableProcessors()

workflow {

    generate_chrom_sizes(file(params.assembly_report))

    
    run_mast(file(params.genome_assembly),
             file(params.ctcf_motif))
    
    convert_mast_to_bed(run_mast.out.txt_gz,
                        generate_chrom_sizes.out.bed)
    
    
    generate_extr_barriers_bed(convert_mast_to_bed.out.bed_gz,
                               file(params.hela_ctcf_chip),
                               file(params.hela_rad21_chip))

        
}

process generate_chrom_sizes {
    publishDir "${params.output_dir}", mode: 'copy'
    
    label 'process_short'
    
    input:
        path assembly_report

    output:
        path "${params.assembly_name_short}.bed", emit: bed
        path "${params.assembly_name_short}.chrom.sizes", emit: chrom_sizes

    shell:
        out_bed = "${params.assembly_name_short}.bed"
        out = "${params.assembly_name_short}.chrom.sizes"
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
        path genome_assembly_fa
        path ctcf_motif_meme

    output:
        path "${params.assembly_name_short}_${params.cell_line_name}_CTCF_mast_hits.txt.gz", emit: txt_gz

    shell:
        out = "${params.assembly_name_short}_${params.cell_line_name}_CTCF_mast_hits.txt.gz"
        """
        set -e
        set -u
        set -o pipefail
        
        # MAST does not like gzipped files nor FIFOs, so we have to actually write the inflated reference to disk
        gzip -dc "$genome_assembly_fa" > ga.tmp.fa &
        gzip -dc "$ctcf_motif_meme" > motif.meme &
        wait

        mast -hit_list     \
             motif.meme    \
             ga.tmp.fa     |
        gzip -9 > "$out"
        
        rm -f ga.tmp.fa motif.meme
        """
}

process convert_mast_to_bed {
    publishDir "${params.output_dir}", mode: 'copy'
    
    label 'process_short'
    
    input:
        path mast_output_txt
        path chrom_sizes_bed

    output:
        path "${params.assembly_name_short}_${params.cell_line_name}_CTCF_mast_hits.bed.gz", emit: bed_gz

    shell: 
        out = "${params.assembly_name_short}_${params.cell_line_name}_CTCF_mast_hits.bed.gz"
        '''
        set -e
        set -u
        set -o pipefail
        
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
        gzip -9 > "!{out}"
        
        rm -f tmp.fifo
        '''
}


process generate_extr_barriers_bed {
    publishDir "${params.output_dir}", mode: 'copy'
    
    label 'process_short'
    
    input:
        path ctcf_sites_bed
        path ctcf_chip_peaks
        path rad21_chip_peaks

    output:
        path "${params.assembly_name_short}_${params.cell_line_name}_CTCF_sites_filtered.bed.gz", emit: bed_gz

    shell:
        out = "${params.assembly_name_short}_${params.cell_line_name}_CTCF_sites_filtered.bed.gz"
        """
        set -e
        set -u
        set -o pipefail
        
        mkfifo tmp1.fifo
        mkfifo tmp2.fifo

        # The following bedtools intersect commands output records from file -a whenever
        # a record from -a overlaps with one or more records from b

        # Intersect CTCF sites with CTCF ChIP peaks
        bedtools intersect              \
                 -wa                    \
                 -u                     \
                 -a "$ctcf_sites_bed"   \
                 -b "$ctcf_chip_peaks" > tmp1.fifo &
                    
        # Intersect CTCF sites with RAD21 ChIP peaks
        bedtools intersect              \
                 -wa                    \
                 -u                     \
                 -a "$ctcf_sites_bed"   \
                 -b "$rad21_chip_peaks" > tmp2.fifo &

        # Equivalent of intersecting CTCF sites, CTCF ChIP peaks and RAD21 ChIP peaks 
        bedtools intersect        \
                 -wa              \
                 -u               \
                 -a tmp1.fifo     \
                 -b tmp2.fifo     |
            sort -V -k 1,1 -k 2,2 |
        gzip -9 > "$out"
        
        rm -f tmp?.fifo
        """
}

