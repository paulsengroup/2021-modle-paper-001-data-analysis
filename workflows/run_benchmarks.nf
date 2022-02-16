#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// For some reason it is not possible to use params.output_dir directly in the publishDir directive
output_dir = params.output_dir

def repeat_csv(csv, n) {
    toks = csv.tokenize(",")
    toks.multiply(n).flatten()
}

def generate_task_ids(tasks) {
    num_tasks = tasks.size()
    (0..num_tasks).toList()
}

workflow {

    if (params.run_weak_scaling_comparison) {
        openmm_gpu_chrom_sizes = repeat_csv(params.md_gpu_chrom_sizes, params.md_replicates)
        openmm_gpu_task_ids = generate_task_ids(openmm_gpu_chrom_sizes)

        openmm_cpu_chrom_sizes = repeat_csv(params.md_cpu_chrom_sizes, params.md_replicates)
        openmm_cpu_task_ids = generate_task_ids(openmm_cpu_chrom_sizes)

         modle_chrom_sizes = repeat_csv(params.modle_chrom_sizes, params.modle_replicates)
         modle_task_ids = generate_task_ids(modle_chrom_sizes)

         benchmark_openmm_gpu_weak_scaling(Channel.of(openmm_gpu_task_ids).flatten(),
                                           file(params.openmm_main_script),
                                           params.chrom_name,
                                           Channel.of(openmm_gpu_chrom_sizes).flatten(),
                                           file(params.extr_barriers_benchmark),
                                           params.md_monomer_size,
                                           params.md_lef_processivity,
                                           params.md_lef_separation)

         benchmark_openmm_cpu_weak_scaling(Channel.of(openmm_cpu_task_ids).flatten(),
                                           file(params.openmm_main_script),
                                           params.chrom_name,
                                           Channel.of(openmm_cpu_chrom_sizes).flatten(),
                                           file(params.extr_barriers_benchmark),
                                           params.md_monomer_size,
                                           params.md_lef_processivity,
                                           params.md_lef_separation)

         benchmark_modle_st_weak_scaling(Channel.of(modle_task_ids).flatten(),
                                         file(params.modle_template_config),
                                         params.chrom_name,
                                         Channel.of(modle_chrom_sizes).flatten(),
                                         params.md_monomer_size,
                                         file(params.extr_barriers_benchmark))

         benchmark_modle_mt_weak_scaling(Channel.of(modle_task_ids).flatten(),
                                         file(params.modle_template_config),
                                         params.chrom_name,
                                         Channel.of(modle_chrom_sizes).flatten(),
                                         params.md_monomer_size,
                                         file(params.extr_barriers_benchmark))

         out_tars = Channel.empty()
                           .concat(benchmark_openmm_gpu_weak_scaling.out.tar,
                                   benchmark_openmm_cpu_weak_scaling.out.tar)
                           .flatten()

         benchmark_openmm_gpu_weak_scaling.out.report.collectFile(keepHeader: true,
                                                                  name: "${output_dir}/openmm_gpu_weak_scaling_benchmark_report.tsv",
                                                                  sort: false,
                                                                  skip: 1)

         benchmark_openmm_cpu_weak_scaling.out.report.collectFile(keepHeader: true,
                                                                  name: "${output_dir}/openmm_cpu_weak_scaling_benchmark_report.tsv",
                                                                  sort: false,
                                                                  skip: 1)

         benchmark_modle_st_weak_scaling.out.report.collectFile(keepHeader: true,
                                                                name: "${output_dir}/modle_st_weak_scaling_benchmark_report.tsv",
                                                                sort: false,
                                                                skip: 1)

         benchmark_modle_mt_weak_scaling.out.report.collectFile(keepHeader: true,
                                                                name: "${output_dir}/modle_mt_weak_scaling_benchmark_report.tsv",
                                                                sort: false,
                                                                skip: 1)

        if (params.run_openmm_to_cool) {

            benchmark_openmm_to_cool(file(params.openmm_to_cool_script),
                                     out_tars,
                                     params.chrom_name,
                                     params.md_monomer_size)

            benchmark_openmm_to_cool.out.report.collectFile(keepHeader: true,
                                                            name: "${output_dir}/openmm_to_cool_benchmark_report.tsv",
                                                            sort: false,
                                                            skip: 1)
        }
    }

    if (params.run_strong_scaling_comparison) {
        nthreads = (0..params.max_cpus).step(4).collect { Math.max(1, it) }
        modle_nthreads = repeat_csv(params.modle_strong_scaling_nthreads, params.modle_replicates)
        modle_task_ids = generate_task_ids(modle_nthreads)

        benchmark_modle_strong_scaling(Channel.of(modle_task_ids).flatten(),
                                       file(params.modle_template_config),
                                       file(params.hg38_chrom_sizes),
                                       file(params.hg38_extr_barriers),
                                       params.md_monomer_size,
                                       Channel.of(modle_nthreads).flatten())

        benchmark_modle_strong_scaling.out.report.collectFile(keepHeader: true,
                                                              name: "${output_dir}/modle_strong_scaling_benchmark_report.tsv",
                                                              sort: false,
                                                              skip: 1)
    }

    if (params.run_modle_ncells) {
        modle_ncells = repeat_csv(params.modle_ncells_benchmark, params.modle_replicates)
        modle_task_ids = generate_task_ids(modle_ncells)

        benchmark_modle_ncells(Channel.of(modle_task_ids).flatten(),
                               file(params.modle_template_config),
                               file(params.hg38_chrom_sizes),
                               file(params.hg38_extr_barriers),
                               Channel.of(modle_ncells).flatten())

        benchmark_modle_ncells.out.report.collectFile(keepHeader: true,
                                                      name: "${output_dir}/modle_ncells_benchmark_report.tsv",
                                                      sort: false,
                                                      skip: 1)
    }
}

process benchmark_openmm_gpu_weak_scaling {
    cpus 1
    label 'error_retry'
    label 'process_very_long'

    input:
        val id
        path main_script
        val chrom_name
        val chrom_size
        path extr_barriers
        val monomer_size
        val lef_processivity
        val lef_separation

    output:
        stdout emit: report
        path '*.tar', emit: tar

    shell:
        '''
        # set -x

        CC="$(which gcc)"
        CXX="$(which g++)"

        chrom_size_bp=$(( !{chrom_size} * 1000000 ))
        out_prefix="$(printf '!{id}_openmm_gpu_%d' $chrom_size_bp)"

        mkdir "$out_prefix"

        run_simulation () {
            outdir=$1
            python '!{main_script}' --gpu                 \
                --output-folder "$outdir"                 \
                --extrusion-barrier-bed !{extr_barriers}  \
                --monomer-size-bp !{monomer_size}         \
                --lef-processivity-bp !{lef_processivity} \
                --lef-separation-bp !{lef_separation}     \
                --simulation-size-mbp !{chrom_size}       \
                --force
        }

        measure_peak_gpu_memory_usage () {
            log_file="$1"
            max=0
            sleep 10
            while true; do
                v=$(nvidia-smi --query-gpu=memory.used --format=csv | tail -n +2 | cut -d ' ' -f 1)
                if [ $v -gt $max ]; then
                    echo "$v" | tee "$log_file" > /dev/null
                    max=$v
                fi
                sleep .1
            done
        }

        export -f run_simulation
        export -f measure_peak_gpu_memory_usage

        measure_peak_gpu_memory_usage "$out_prefix/peak_memory_usage.log" &

        # Here the tee call is used to echo errors to stderr
        report="$(command time -f '%e %S %U %M' bash -c "run_simulation $out_prefix" |&
                    tee >(cat 1>&2) |
                    tail -n 1 |
                    sed 's/[[:space:]]\\+/\\t/g')"

        kill $!

        max_gpu_mem=$(<"$out_prefix/peak_memory_usage.log")
        max_gpu_mem="${max_gpu_mem//$'\\n'/ }"

        record="!{id}\topenmm_gpu_weak_scaling\t$chrom_size_bp\t$(hostname)\t$(date +"%T.%N")\t$report\t$((max_gpu_mem * 1000))"
        printf 'id\\tname\\tchrom_size\\thostname\\tdate\\twall_clock\\tsystem_time\\tuser_time\\tmax_resident_mem_kb\\tmax_resident_mem_kb_gpu\\n%s\\n' "${record//$'\\n'/}" |
            tee "$out_prefix/report.txt"

        tar -cf "$out_prefix.tar" "$out_prefix"
        rm -rf "$out_prefix"
        '''
}

process benchmark_openmm_cpu_weak_scaling {
     time {
        // Request 12h for every Mbp of DNA
        cs = Long.parseLong(chrom_size)
        time = 12 * 3600 * 1000 * cs as Long
        time * task.attempt
     }
    label 'error_retry'
    cpus = 10

    input:
        val id
        path main_script
        val chrom_name
        val chrom_size
        path extr_barriers
        val monomer_size
        val lef_processivity
        val lef_separation

    output:
        stdout emit: report
        path '*.tar', emit: tar

    shell:
        '''
        # set -x

        CC="$(which gcc)"
        CXX="$(which g++)"

        chrom_size_bp=$(( !{chrom_size} * 1000000 ))
        out_prefix="$(printf '!{id}_openmm_cpu_%d' $chrom_size_bp)"

        run_simulation () {
            outdir="$1"
            OPENMM_CPU_THREADS=16                         \
            python '!{main_script}' --cpu                 \
                --output-folder="$outdir"                 \
                --extrusion-barrier-bed !{extr_barriers}  \
                --monomer-size-bp !{monomer_size}         \
                --lef-processivity-bp !{lef_processivity} \
                --lef-separation-bp !{lef_separation}     \
                --simulation-size-mbp !{chrom_size}
        }

        export -f run_simulation

        # Here the tee call is used to echo errors to stderr
        report="$(command time -f '%e %S %U %M' bash -c "run_simulation $out_prefix" |&
                    tee >(cat 1>&2) |
                    tail -n 1 |
                    sed 's/[[:space:]]\\+/\\t/g')"

        record="!{id}\t!{task.cpus}\topenmm_cpu_weak_scaling\t$chrom_size_bp\t$(hostname)\t$(date +"%T.%N")\t$report"
        printf 'id\\tnthreads\\tname\\tchrom_size\\thostname\\tdate\\twall_clock\\tsystem_time\\tuser_time\\tmax_resident_mem_kb\\n%s\\n' "${record//$'\\n'/}" |
            tee "$out_prefix/report.txt"

        tar -cf "$out_prefix.tar" "$out_prefix"
        rm -rf "$out_prefix"
        '''
}

process benchmark_openmm_to_cool {
    label 'process_very_high'
    label 'error_retry'

    time {
        toks = "${openmm_tar}" =~ /.*_(\d+).tar$/
        sim_size = Float.parseFloat(toks[0][1])

        // Request 1h every 100 Mbp of DNA, and not less than 1h
        time = Math.ceil(3600 * (sim_size / 100.0e6)) as Long
        Math.max(3600, time * task.attempt) * 1000
    }

    input:
        path main_script
        path openmm_tar
        val chrom_name
        val bin_size

    output:
        stdout emit: report
        path "*.cool", emit: cool

    shell:
        '''
        # set -x

        tar -xf '!{openmm_tar}'
        input_dir="$(basename '!{openmm_tar}' .tar)"

        # input_dir has the following format (id)_openmm_[cg]pu_(chrom_size)
        id=$(echo "$input_dir" | grep -o '^[[:digit:]]\\+')
        name=$(echo "$input_dir" | grep -o 'openmm_[cg]pu')
        chrom_size_bp=$(echo "$input_dir" | grep -o '[[:digit:]]\\+$')

        mkdir tmp

        printf "%s\\t%d\\n" !{chrom_name} "$chrom_size_bp" > tmp/chrom.sizes

        make_cool_file () {
            input_dir="$1"
            chrom_size_bp=$2
            cool_name="$3.cool"

            polychrom_traj_convert --allow-nonconsecutive "$input_dir" tmp/
            python '!{main_script}'           \
                --input-folders tmp           \
                --output-name "$cool_name"    \
                --chrom-sizes tmp/chrom.sizes \
                --bin-size !{bin_size}        \
                --chrom-names !{chrom_name}   \
                --threads !{task.cpus}        \
                --offsets 0
        }

        export -f make_cool_file
        # Here the tee call is used to echo errors to stderr
        report="$(command time -f '%e %S %U %M' bash -c \
                    "make_cool_file $input_dir $chrom_size_bp $input_dir" |&
                    tee >(cat 1>&2) |
                    tail -n 1 |
                    sed 's/[[:space:]]\\+/\\t/g')"

        rm -rf tmp/ "$input_dir"

        record="$id\t!{task.cpus}\tmake_cool_file\t$chrom_size_bp\t$(hostname)\t$(date +"%T.%N")\t$report"
        printf 'id\\tnthreads\\tname\\tchrom_size\\thostname\\tdate\\twall_clock\\tsystem_time\\tuser_time\\tmax_resident_mem_kb\\n%s\\n' "${record//$'\\n'/}"
        '''
}

process benchmark_modle_st_weak_scaling {
    label 'error_retry'
    cpus 1

    input:
        val id
        path base_config
        val chrom_name
        val chrom_size
        val bin_size
        path extr_barriers

    output:
        stdout emit: report
        path "*.log", emit: log
        path "*.cool", emit: cool

    shell:
        '''
        set -o pipefail

        mkdir tmp

        chrom_size_bp=$(( !{chrom_size} * 1000000 ))
        out_prefix="$(printf '!{id}_modle_st_%d' $chrom_size_bp)"
        printf '%s\\t%d\\n' '!{chrom_name}' $chrom_size_bp > tmp/tmp.chrom.sizes

        sed 's|^chrom-sizes=.*|chrom-sizes="tmp/tmp.chrom.sizes"|' "!{base_config}" |
        sed 's|^extrusion-barrier-file=.*|extrusion-barrier-file="!{extr_barriers}"|' |
        sed 's|^bin-size=.*|bin-size=!{bin_size}|' |
        sed 's|^threads=.*|threads=!{task.cpus}|' |
        sed "s|^output-prefix=.*|output-prefix=\"$out_prefix\"|" > tmp/config.tmp.toml

        report="$(command time -f '%e %S %U %M' modle --config tmp/config.tmp.toml |&
                    tee >(cat 1>&2) |
                    tail -n 1 |
                    sed 's/[[:space:]]\\+/\\t/g')"

        rm -rf tmp/

        record="!{id}\t!{task.cpus}\tmodle_st_weak_scaling\t$chrom_size_bp\t$(hostname)\t$(date +"%T.%N")\t$report"
        printf 'id\\tnthreads\\tname\\tchrom_size\\thostname\\tdate\\twall_clock\\tsystem_time\\tuser_time\\tmax_resident_mem_kb\\n%s\\n' "${record//$'\\n'/}" |
            tee -a "$out_prefix.log"
        '''
}

process benchmark_modle_mt_weak_scaling {
    label 'process_very_high'
    label 'error_retry'
    label 'process_short'

    input:
        val id
        path base_config
        val chrom_name
        val chrom_size
        val bin_size
        path extr_barriers

    output:
        stdout emit: report
        path "*.log", emit: log
        path "*.cool", emit: cool

    shell:
        '''
        set -o pipefail

        mkdir tmp

        chrom_size_bp=$(( !{chrom_size} * 1000000 ))
        out_prefix="$(printf '!{id}_modle_mt_%d_%d' !{task.cpus} $chrom_size_bp)"
        printf '%s\\t%d\\n' '!{chrom_name}' $chrom_size_bp > tmp/tmp.chrom.sizes

        sed 's|^chrom-sizes=.*|chrom-sizes="tmp/tmp.chrom.sizes"|' "!{base_config}" |
        sed 's|^extrusion-barrier-file=.*|extrusion-barrier-file="!{extr_barriers}"|' |
        sed 's|^bin-size=.*|bin-size=!{bin_size}|' |
        sed 's|^threads=.*|threads=!{task.cpus}|' |
        sed "s|^output-prefix=.*|output-prefix=\"$out_prefix\"|" > tmp/config.tmp.toml

        report="$(command time -f '%e %S %U %M' modle --config tmp/config.tmp.toml |&
                    tee >(cat 1>&2) |
                    tail -n 1 |
                    sed 's/[[:space:]]\\+/\\t/g')"

        rm -rf tmp/

        record="!{id}\t!{task.cpus}\tmodle_mt_weak_scaling\t$chrom_size_bp\t$(hostname)\t$(date +"%T.%N")\t$report"
        printf 'id\\tnthreadst\\tname\\tchrom_size\\thostname\\tdate\\twall_clock\\tsystem_time\\tuser_time\\tmax_resident_mem_kb\\n%s\\n' "${record//$'\\n'/}" |
            tee -a "$out_prefix.log"
        '''
}

process benchmark_modle_strong_scaling {
    label 'error_retry'

    cpus { nthreads }
    memory { task.attempt * 10e9 as Long }
    time {
          thr = Float.parseFloat(nthreads)
          time = Math.ceil((3600 * 8) / thr) as Long
          1000 * time * task.attempt
         }

    input:
        val id
        path base_config
        path chrom_sizes
        path extr_barriers
        val bin_size
        val nthreads

    output:
        stdout emit: report
        path "*.log", emit: log
        path "*.cool", emit: cool

    shell:
        '''
        set -o pipefail

        mkdir tmp

        out_prefix="$(printf '!{id}_modle_strong_scaling_%d' !{task.cpus})"

        sed 's|^chrom-sizes=.*|chrom-sizes="!{chrom_sizes}"|' "!{base_config}" |
        sed 's|^extrusion-barrier-file=.*|extrusion-barrier-file="!{extr_barriers}"|' |
        sed 's|^bin-size=.*|bin-size=!{bin_size}|' |
        sed 's|^threads=.*|threads=!{task.cpus}|' |
        sed "s|^output-prefix=.*|output-prefix=\"$out_prefix\"|" |
        tee "tmp/config.tmp.toml" > /dev/null

        report="$(command time -f '%e %S %U %M' modle --config tmp/config.tmp.toml |&
                    tee >(cat 1>&2) |
                    tail -n 1 |
                    sed 's/[[:space:]]\\+/\\t/g')"

        rm -rf tmp/

        record="!{id}\t!{task.cpus}\tmodle_strong_scaling\t$(hostname)\t$(date +"%T.%N")\t$report"
        printf 'id\\tnthreads\\tname\\thostname\\tdate\\twall_clock\\tsystem_time\\tuser_time\\tmax_resident_mem_kb\\n%s\\n' "${record//$'\\n'/}" |
            tee -a "$out_prefix.log"
        '''
}

process benchmark_modle_ncells {
    label 'error_retry'

    cpus { params.max_cpus}
    memory { task.attempt * 10e9 as Long }
    time {
          cells = Float.parseFloat(ncells)
          time = Math.ceil(0.2 * cells) as Long
          Math.max(600, time) * 1000 * task.attempt
         }

    input:
        val id
        path base_config
        path chrom_sizes
        path extr_barriers
        val ncells

    output:
        stdout emit: report
        path "*.log", emit: log
        path "*.cool", emit: cool

    shell:
        '''
        set -o pipefail

        out_prefix='!{id}_modle_ncells_!{ncells}'

        mkdir tmp

        sed 's|^chrom-sizes=.*|chrom-sizes="!{chrom_sizes}"|' "!{base_config}" |
        sed 's|^extrusion-barrier-file=.*|extrusion-barrier-file="!{extr_barriers}"|' |
        sed 's|^threads=.*|threads=!{task.cpus}|' |
        sed 's|^ncells=.*|ncells=!{ncells}|' |
        sed "s|^output-prefix=.*|output-prefix=\"$out_prefix\"|" |
        tee tmp/config.tmp.toml > /dev/null

        report="$(command time -f '%e %S %U %M' modle --config tmp/config.tmp.toml |&
                    tee >(cat 1>&2) |
                    tail -n 1 |
                    sed 's/[[:space:]]\\+/\\t/g')"
        rm -r tmp/

        record="!{id}\tmodle_ncells\t!{ncells}\t$(hostname)\t$(date +"%T.%N")\t$report"
        printf 'id\\tname\\tncells\\thostname\\tdate\\twall_clock\\tsystem_time\\tuser_time\\tmax_resident_mem_kb\\n%s\\n' "${record//$'\\n'/}" |
            tee "$out_prefix.log"
        '''
}
