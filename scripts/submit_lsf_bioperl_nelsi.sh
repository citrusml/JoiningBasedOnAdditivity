#!/bin/bash
# ==============================================================================
# High-Parallelism LSF (bsub) Job Array Script: BioPerl vs. NELSI Benchmark
# (Includes Automatic Post-Array CSV Merge Job & PhyloTree_Methods_Comparison Log Format)
#
# Evaluates:
#   - Tree Models: bioperl, nelsi (2 models)
#   - Taxa N:      20, 32, 64 (3 levels)
#   - Divergence D: 0.1, 0.2, 0.5, 1.0, 1.5, 2.0, 3.0 (7 levels)
#   - Length L:    500 aa
#   - Indel Rate:  0.1 (INDELible Zipfian POW 1.7 50)
#   - Replicates:  50 per condition (split into 10 shards -> 5 reps/job)
#
# Total LSF Array Jobs: 2 x 3 x 7 x 10 = 420 parallel jobs + 1 merge job
#   - Jobs   1..210 : tree_model = bioperl
#   - Jobs 211..420 : tree_model = nelsi
#   - Merge job     : runs automatically after all 420 jobs complete
#
# Recommended Usage (submits 420-job array + dependent merge job):
#   bash scripts/submit_lsf_bioperl_nelsi.sh
#
# Alternative (manual merge):
#   bash scripts/submit_lsf_bioperl_nelsi.sh --merge
# ==============================================================================

#BSUB -J bioperl_nelsi[1-420]
#BSUB -q mafft
#BSUB -n 1
#BSUB -M 4000
#BSUB -W 04:00
#BSUB -o logs/bioperl_nelsi_bsub_%I.out  # BSUBジョブ自体のstdoutを記録（OOMなど診断に必須）
#BSUB -e logs/bioperl_nelsi_bsub_%I.err  # BSUBジョブ自体のstderrを記録

set -uo pipefail

OUTDIR=${OUTDIR:-results_bioperl_nelsi}
TOTAL_REPS=${REPLICATES:-50}
SEQ_LEN=${SEQ_LENGTH:-500}
INDEL_RATE=${INDEL_RATE:-0.1}

_activate_env() {
    export PATH="$HOME/bin:$PATH"
    eval "$(micromamba shell hook --shell bash)" 2>/dev/null || true
    micromamba activate phylomethod_env
    export PATH="${CONDA_PREFIX:-}:$CONDA_PREFIX/bin:$HOME/bin:$PATH"
    export NELSI_RSCRIPT="${NELSI_RSCRIPT:-${MAMBA_ROOT_PREFIX:-$HOME/.micromamba}/envs/nelsi_env/bin/Rscript}"
}

# ------------------------------------------------------------------------------
# Mode 1: Merge-only execution (called by dependent LSF job or manually)
# ------------------------------------------------------------------------------
if [ "${1:-}" = "--merge" ]; then
    _activate_env
    cd "${LS_SUBCWD:-$(pwd)}"
    mkdir -p logs
    {
        echo "=== Merging all shard CSVs in ${OUTDIR}/shards/ | Start: $(date) ==="
        python -u benchmark_bioperl_nelsi.py --merge-only --outdir "${OUTDIR}"
        echo "=== Merge completed at $(date) ==="
    } > logs/bioperl_nelsi_merge.log 2>&1
    exit $?
fi

# ------------------------------------------------------------------------------
# Mode 2: Login-node launcher (submits array job + dependent merge job)
# ------------------------------------------------------------------------------
if [ -z "${LSB_JOBID:-}" ]; then
    mkdir -p logs "${OUTDIR}/shards"
    SCRIPT_PATH="$(cd "$(dirname "$0")" && pwd)/$(basename "$0")"

    echo "Submitting 420-job LSF array via ${SCRIPT_PATH}..."
    SUBMIT_OUT=$(AUTO_MERGE_SUBMITTED=1 OUTDIR="${OUTDIR}" REPLICATES="${TOTAL_REPS}" SEQ_LENGTH="${SEQ_LEN}" INDEL_RATE="${INDEL_RATE}" bsub < "${SCRIPT_PATH}")
    echo "${SUBMIT_OUT}"

    ARRAY_JOBID=$(echo "${SUBMIT_OUT}" | sed -n 's/.*Job <\([0-9]*\)>.*/\1/p')
    if [ -n "${ARRAY_JOBID}" ]; then
        echo "Submitting dependent merge job waiting on done(${ARRAY_JOBID})..."
        OUTDIR="${OUTDIR}" bsub \
            -w "done(${ARRAY_JOBID})" \
            -J "bioperl_nelsi_merge" \
            -q mafft \
            -n 1 -M 4000 -W 01:00 \
            -o "logs/bioperl_nelsi_merge_bsub.out" \
            -e "logs/bioperl_nelsi_merge_bsub.err" \
            "bash \"${SCRIPT_PATH}\" --merge"
    else
        echo "Warning: Could not parse array Job ID. Run 'bash ${SCRIPT_PATH} --merge' after jobs finish."
    fi
    exit 0
fi

# ------------------------------------------------------------------------------
# Mode 3: LSF compute-node worker (executes one shard of one condition)
# ------------------------------------------------------------------------------
_activate_env

cd "${LS_SUBCWD:-$(pwd)}"
mkdir -p logs "${OUTDIR}/shards"

# Fallback auto-schedule merge job if user ran `bsub < scripts/submit_lsf_bioperl_nelsi.sh` directly
if [ "${LSB_JOBINDEX:-0}" -eq 1 ] && [ "${AUTO_MERGE_SUBMITTED:-0}" != "1" ]; then
    OUTDIR="${OUTDIR}" bsub \
        -w "done(${LSB_JOBID})" \
        -J "bioperl_nelsi_merge" \
        -q mafft \
        -n 1 -M 4000 -W 01:00 \
        -o "logs/bioperl_nelsi_merge_bsub.out" \
        -e "logs/bioperl_nelsi_merge_bsub.err" \
        "bash scripts/submit_lsf_bioperl_nelsi.sh --merge" >/dev/null 2>&1 || true
fi

# Ensure INDELible v1.03 binary is available and compiled for this host OS/arch
bash scripts/install_indelible.sh

MODEL_LIST=(bioperl nelsi)
TAXA_LIST=(20 32 64)
DIV_LIST=(0.1 0.2 0.5 1.0 1.5 2.0 3.0)

N_TAXA_OPTS=${#TAXA_LIST[@]}     # 3
N_DIV_OPTS=${#DIV_LIST[@]}       # 7
N_COND_PER_MODEL=$(( N_TAXA_OPTS * N_DIV_OPTS ))  # 21
N_SHARDS=10
JOBS_PER_MODEL=$(( N_COND_PER_MODEL * N_SHARDS )) # 210

JOB_IDX=$(( ${LSB_JOBINDEX:-1} - 1 ))
MODEL_IDX=$(( JOB_IDX / JOBS_PER_MODEL ))
REM_MODEL=$(( JOB_IDX % JOBS_PER_MODEL ))
COND_IDX=$(( REM_MODEL / N_SHARDS ))
SHARD_IDX=$(( REM_MODEL % N_SHARDS ))

TAXA_IDX=$(( COND_IDX / N_DIV_OPTS ))
DIV_IDX=$(( COND_IDX % N_DIV_OPTS ))

TREE_MODEL=${MODEL_LIST[$MODEL_IDX]}
N_TAXA=${TAXA_LIST[$TAXA_IDX]}
DIV=${DIV_LIST[$DIV_IDX]}

export OMP_NUM_THREADS=1
LOG_FILE="logs/bioperl_nelsi_${TREE_MODEL}_N${N_TAXA}_D${DIV}_s${SHARD_IDX}.log"

# リアルタイムに進捗ログを出力して実行 (いつでも tail -f logs/bioperl_nelsi_*.log で確認可能)
{
    echo "=== Running ${TREE_MODEL} N=${N_TAXA}, D=${DIV}, L=${SEQ_LEN}, Indel=${INDEL_RATE}, Shard=${SHARD_IDX}/${N_SHARDS} | Host: $(hostname) | Start: $(date) ==="
    python -u benchmark_bioperl_nelsi.py \
      --tree-models "${TREE_MODEL}" \
      --taxa "${N_TAXA}" \
      --divergences "${DIV}" \
      --length "${SEQ_LEN}" \
      --indel-rate "${INDEL_RATE}" \
      --replicates "${TOTAL_REPS}" \
      --shard-idx "${SHARD_IDX}" \
      --num-shards "${N_SHARDS}" \
      --threads 1 \
      --skip-existing \
      --outdir "${OUTDIR}"
    echo "=== Finished ${TREE_MODEL} N=${N_TAXA}, D=${DIV}, Shard=${SHARD_IDX} | End: $(date) ==="
} > "${LOG_FILE}" 2>&1
