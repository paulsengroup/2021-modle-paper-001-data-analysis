#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// For some reason it is not possible to use params.output_dir directly in the publishDir directive
output_dir = params.output_dir

def repeat_channel_of_paths(input_ch, n) {
    input_ch.toList()
            .multiply(n)
            .flatten()
            .toSortedList({ f1, f2 -> file(f1).getBaseName() <=> file(f2).getBaseName() })
            .flatten()
}

workflow {
    gen_dsets_py = file(params.generate_test_datasets_script)
    hoomd_base_params = file(params.hoomd_base_params_txt)
    hoomd_base_probs = file(params.hoomd_base_probabilities_tsv)

    seed = 0
    for (line: hoomd_base_params.readLines()) {
        if ( line ==~ /^seed\s*=.*\d+$/) {
            seed = line.replaceFirst(/seed\s*=\s*/, "") as Long
        }
    }



    generate_test_datasets_hoomd(gen_dsets_py,
                                 hoomd_base_params,
                                 hoomd_base_probs,
                                 params.hoomd_bin_size,
                                 params.hoomd_min_size_mbp,
                                 params.hoomd_max_size_mbp,
                                 params.hoomd_step_mbp)

    generate_test_datasets_modle(gen_dsets_py,
                                 hoomd_base_probs,
                                 params.chrom_name,
                                 params.hoomd_bin_size,
                                 params.modle_min_size_mbp,
                                 params.modle_max_size_mbp,
                                 params.modle_step_mbp)

    max_num_monomers = {
        chrom_size = params.hoomd_max_size_mbp as Float * 1.0e6
        bin_size = params.hoomd_bin_size as Float

        Math.ceil(chrom_size / bin_size) as Long
    }

    generate_random_polymer_hoomd(file(params.generate_random_polymer_script),
                                  max_num_monomers,
                                  seed)


    hoomd_params = repeat_channel_of_paths(generate_test_datasets_hoomd.out.params,
                                          params.hoomd_runs)
    hoomd_probs = repeat_channel_of_paths(generate_test_datasets_hoomd.out.probs,
                                          params.hoomd_runs)

    modle_chrom_sizes = repeat_channel_of_paths(generate_test_datasets_modle.out.chrom_sizes,
                                                params.modle_runs)
    modle_barriers = repeat_channel_of_paths(generate_test_datasets_modle.out.extr_barriers_bed,
                                                params.modle_runs)

    benchmark_hoomd_cpu_weak_scaling(file(params.hoomd_main_script),
                                     generate_random_polymer_hoomd.out.structure,
                                     hoomd_params,
                                     hoomd_probs)

    benchmark_modle_cpu_weak_scaling(file(params.modle_template_config),
                                     modle_chrom_sizes,
                                     modle_barriers)


    nthreads = Channel.of(params.modle_nthreads_strong_scaling_bench
                                .tokenize(",")
                                .multiply(params.modle_runs))
                                .flatten()

    benchmark_modle_cpu_strong_scaling(file(params.modle_template_config),
                                       modle_chrom_sizes.last(),
                                       modle_barriers.last(),
                                       nthreads)

    benchmark_hoomd_cpu_weak_scaling.out.report.collectFile(keepHeader: true,
                                                            name: "${output_dir}/hoomd_weak_scaling_benchmark_report.tsv",
                                                            sort: false,
                                                            skip: 1)

    benchmark_modle_cpu_weak_scaling.out.report.collectFile(keepHeader: true,
                                                            name: "${output_dir}/modle_weak_scaling_benchmark_report.tsv",
                                                            sort: false,
                                                            skip: 1)

    benchmark_modle_cpu_strong_scaling.out.report.collectFile(keepHeader: true,
                                                              name: "${output_dir}/modle_strong_scaling_benchmark_report.tsv",
                                                              sort: false,
                                                              skip: 1)
}

process generate_test_datasets_hoomd {
    publishDir "${output_dir}/datasets/hoomd", mode: 'copy'
    label 'process_very_short'

    input:
        path script
        path param_file
        path prob_file
        val monomer_size
        val min_chrom_size
        val max_chrom_size
        val step_size

    output:
        path "*.txt", emit: params
        path "*.tsv", emit: probs

    shell:
        '''
        python "!{script}" \
                --mode hoomd                              \
                --hoomd-input-parameters "!{param_file}"   \
                --hoomd-input-probabilities "!{prob_file}" \
                --monomer-size-bp !{monomer_size}         \
                --output-prefix hoomd                      \
                --min-chrom-size-mbp !{min_chrom_size}    \
                --max-chrom-size-mbp !{max_chrom_size}    \
                --step-size-mbp !{step_size}
        '''
}

process generate_random_polymer_hoomd {
    publishDir "${output_dir}/datasets/hoomd", mode: 'copy'
    label 'process_very_short'

    input:
        path script
        val monomers
        val seed

    output:
        path "*.xyz", emit: structure

    shell:
        '''
            python '!{script}'                  \
               --number-of-monomers !{monomers} \
               --seed !{seed} > random_polymer_!{monomers}_beads.xyz
        '''
}

process generate_test_datasets_modle {
    publishDir "${output_dir}/datasets/modle", mode: 'copy'
    label 'process_very_short'

    input:
        path script
        path prob_file
        val chrom_name
        val monomer_size
        val min_chrom_size
        val max_chrom_size
        val step_size

    output:
        path "*.chrom.sizes", emit: chrom_sizes
        path "*.bed", emit: extr_barriers_bed

    shell:
        '''
        python "!{script}" \
                --mode modle                              \
                --chrom-name "!{chrom_name}"              \
                --hoomd-input-probabilities "!{prob_file}" \
                --monomer-size-bp !{monomer_size}         \
                --output-prefix modle                      \
                --min-chrom-size-mbp !{min_chrom_size}    \
                --max-chrom-size-mbp !{max_chrom_size}    \
                --step-size-mbp !{step_size}
        '''
}

process benchmark_hoomd_cpu_weak_scaling {
    cpus 1
    memory {
        toks = "${param_file}" =~ /.*_(\d+)_bp.txt$/
        sim_size = Float.parseFloat(toks[0][1])

        // Request 8.5 GB every Mbp of DNA
        mem = Math.ceil(8.5e9 * (sim_size / 1.0e6)) as Long
        mem * task.attempt
    }

     time {
        toks = "${param_file}" =~ /.*_(\d+)_bp.txt$/
        sim_size = Float.parseFloat(toks[0][1])

        // Request 4h for every Mbp of DNA
        time = Math.ceil(4 * 3600 * 1000 * (sim_size / 1.0e6)) as Long
        time * task.attempt
     }
    label 'error_retry'

    input:
        path main_script
        path initial_structure
        path param_file
        path prob_file

    output:
        stdout emit: report

    shell:
        '''
        set -x

        num_monomers=$(grep 'N[[:space:]]*=' '!{param_file}' | grep -o '[[:digit:]]\\+$')

        mkdir out/
        awk -F '\\t' 'NF == 4{ print $2,$3,$4 }' '!{initial_structure}' |
            head -n "$num_monomers" > out/polymer.txt

        set -o pipefail

        run_hoomd () {
            python '!{main_script}'                 \
                --output-dir='out/'                 \
                --parameters='!{param_file}'         \
                --probabilities='!{prob_file}'       \
                --initial-structure out/polymer.txt \
                --mode=cpu                          \
                --nthreads=1
        }

        export -f run_hoomd

        run_hoomd

        # Here the tee call is used to echo errors to stderr
        report="$(command time -f '%e %S %U %M' bash -c "run_hoomd" |&
                    tee >(cat 1>&2) |
                    tail -n 1 |
                    sed 's/[[:space:]]\\+/\\t/g')"

        rm -rf out/

        name="$(printf "!{param_file}" | grep -o '[[:digit:]]\\+_bp.txt$')"
        name="${name%_bp.txt}"

        printf 'name\twall_clock\tsystem_time\tuser_time\tmax_resident_mem_kb\\n%s\\t%s\\n' "$name" "$report"
        '''
}

process benchmark_modle_cpu_weak_scaling {
    cpus 1

    input:
        path base_config
        path chrom_sizes
        path extr_barriers

    output:
        stdout emit: report

    shell:
        '''
        set -o pipefail

        mkdir out

        sed 's|^chrom-sizes=.*|chrom-sizes="!{chrom_sizes}"|' "!{base_config}" |
        sed 's|^extrusion-barrier-file=.*|extrusion-barrier-file="!{extr_barriers}"|' |
        sed 's|^threads=.*|threads=!{task.cpus}|' |
        sed 's|^output-prefix=.*|output-prefix=out/modle|' > out/config.tmp.toml

        report="$(command time -f '%e %S %U %M' modle --config out/config.tmp.toml |&
                    tee >(cat 1>&2) |
                    tail -n 1 |
                    sed 's/[[:space:]]\\+/\\t/g')"
        rm -r out/

        name="$(printf "!{chrom_sizes}" | grep -o '[[:digit:]]\\+_bp.chrom.sizes$')"
        name="${name%_bp.chrom.sizes}"

        printf 'name\twall_clock\tsystem_time\tuser_time\tmax_resident_mem_kb\\n%s\\t%s\\n' "$name" "$report"
        '''
}

process benchmark_modle_cpu_strong_scaling {
    cpus { nthreads }

    input:
        path base_config
        path chrom_sizes
        path extr_barriers
        val nthreads

    output:
        stdout emit: report

    shell:
        '''
        set -o pipefail

        mkdir out

        sed 's|^chrom-sizes=.*|chrom-sizes="!{chrom_sizes}"|' "!{base_config}" |
        sed 's|^extrusion-barrier-file=.*|extrusion-barrier-file="!{extr_barriers}"|' |
        sed 's|^threads=.*|threads=!{task.cpus}|' |
        sed 's|^output-prefix=.*|output-prefix="out/benchmark_"|' |
        tee out/config.tmp.toml > /dev/null

        report="$(command time -f '%e %S %U %M' modle --config out/config.tmp.toml |&
                    tee >(cat 1>&2) |
                    tail -n 1 |
                    sed 's/[[:space:]]\\+/\\t/g')"
        rm -r out/

        name="$(printf '!{chrom_sizes}' | grep -o '[[:digit:]]\\+_bp.chrom.sizes$')"
        name="${name%_bp.chrom.sizes}_!{task.cpus}"

        printf 'name\twall_clock\tsystem_time\tuser_time\tmax_resident_mem_kb\\n%s\\t%s\\n' "$name" "$report"
        '''
}


process plot_weak_scaling {
    publishDir "${output_dir}/plots", mode: 'copy'

    label 'process_short'

    input:
        path modle_report
        val num_threads
        val num_cells

    output:
        path "*.png", emit: png
        path "*.svg", emit: svg
    shell:
        '''
#!/usr/bin/env python3

import matplotlib.pyplot as plt
import pandas as pd

modle_report = pd.read_table("!{modle_report}")
hoomd_report = pd.read_table("!{modle_report}") # CHANGEME

nthreads = int("!{num_threads}")
ncells = int("!{num_cells}")

modle_report["chrom_size"] = pd.to_numeric(modle_report["name"].replace(
                                          {r".*_(\\d+)_bp" : r"\\1"},
                                          regex=True),
                                          downcast="signed") / 1.0e6


fig, axs = plt.subplots(2, 1, figsize=(6.4, 6.4*2)

ax.errorbar(modle_report["chrom_size"], modle_report["mean"], modle_report["stddev"])
ax.plot(modle_report["chrom_size"], modle_report["user"] + modle_report["system"])

ax.set(title=f"MoDLE weak scaling ({nthreads} threads, {ncells} cells)",
       xlabel="Chromosome size (Mbp)",
       ylabel="Average wall clock (s)"


ax.set(title=f"MoDLE weak scaling ({nthreads} threads, {ncells} cells)",
       xlabel="Chromosome size (Mbp)",
       ylabel="Wall clock (s)"
fig.tight_layout()

        '''
}
