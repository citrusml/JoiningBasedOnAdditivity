nextflow.enable.dsl=2

/*
 * Nextflow DSL2 Pipeline for NELSI Autocorrelated Variance Benchmark
 *
 * JoiningBasedOnAdditivity の benchmark_bioperl_nelsi.py を Nextflow で並列実行し、
 * 各シャード結果を集約して benchmark_results.csv / summary_by_condition.csv を生成する。
 *
 * Execution:
 *   nextflow run main.nf -c next_configs/nextflow_nelsi_v06_v10.config -profile supercomputer -resume
 */

process RUN_NELSI_CHUNK {
    tag "v=${nelsi_v}_ir=${indel_rate}_N=${n_taxa}_D=${divergence}_c=${chunk_idx}/${num_chunks}"
    publishDir "${params.outdir}/shards", mode: 'copy', overwrite: true

    input:
    tuple val(nelsi_v), val(indel_rate), val(n_taxa), val(divergence), val(chunk_idx), val(num_chunks)

    output:
    path "shards/*.csv", emit: shard_csv

    script:
    def shard_idx = chunk_idx - 1
    """
    ln -sfn "${projectDir}/src" src
    export PYTHONPATH="${projectDir}:\${PYTHONPATH:-}"
    export NELSI_RSCRIPT="\${NELSI_RSCRIPT:-\${MAMBA_ROOT_PREFIX:-\$HOME/.micromamba}/envs/nelsi_env/bin/Rscript}"

    python -u "${projectDir}/benchmark_bioperl_nelsi.py" \\
        --tree-models "${params.tree_model}" \\
        --nelsi-v ${nelsi_v} \\
        --taxa ${n_taxa} \\
        --divergences ${divergence} \\
        --lba-ratios ${params.lba_ratio} \\
        --length ${params.length} \\
        --indel-rate ${indel_rate} \\
        --replicates ${params.replicates} \\
        --methods-set "${params.methods_set}" \\
        --threads ${task.cpus} \\
        --max-iter-search ${params.max_iter_search} \\
        --shard-idx ${shard_idx} \\
        --num-shards ${num_chunks} \\
        --seed ${params.seed} \\
        --outdir . \\
        --skip-existing
    """
}

process MERGE_RESULTS {
    tag "merge_all"
    publishDir "${params.outdir}", mode: 'copy', overwrite: true

    input:
    path shard_files, stageAs: "shards/*"

    output:
    path "benchmark_results.csv", emit: results_csv
    path "summary_by_condition.csv", emit: summary_csv

    script:
    """
    ln -sfn "${projectDir}/src" src
    export PYTHONPATH="${projectDir}:\${PYTHONPATH:-}"

    python -u "${projectDir}/benchmark_bioperl_nelsi.py" \\
        --merge-only \\
        --outdir .
    """
}

workflow {
    def v_list     = (params.nelsi_v_list instanceof Collection) ? params.nelsi_v_list.flatten() : [params.nelsi_v_list]
    def indel_src  = params.containsKey('indel_rates') ? params.indel_rates : params.indel_rate
    def indel_list = (indel_src instanceof Collection) ? indel_src.flatten() : [indel_src]
    def taxa_list  = (params.taxa_list    instanceof Collection) ? params.taxa_list.flatten()    : [params.taxa_list]
    def dist_list  = (params.distances    instanceof Collection) ? params.distances.flatten()    : [params.distances]

    def reps       = params.replicates as int
    def chk_size   = params.chunk_size as int
    def num_chunks = Math.max(1, Math.ceil(reps / chk_size) as int)

    ch_v      = Channel.fromList(v_list)
    ch_indel  = Channel.fromList(indel_list)
    ch_taxa   = Channel.fromList(taxa_list)
    ch_dist   = Channel.fromList(dist_list)
    ch_chunks = Channel.fromList((1..num_chunks).collect { c -> [c, num_chunks] })

    ch_tasks = ch_v
        .combine(ch_indel)
        .combine(ch_taxa)
        .combine(ch_dist)
        .combine(ch_chunks)
        .map { nv, ir, nt, d, c_idx, n_chk -> tuple(nv, ir, nt, d, c_idx, n_chk) }

    RUN_NELSI_CHUNK(ch_tasks)
    MERGE_RESULTS(RUN_NELSI_CHUNK.out.shard_csv.collect())
}
