#!/bin/bash
# ==============================================================================
# High-Parallelism LSF (bsub) Job Array Script for Zipfian Indel (INDELible)
#
# Evaluates both Indel Rates (0.05 and 0.1) under Zipfian Indel model (POW 1.7 50)
# across all 12 (N, D) conditions:
#   - Indel Rate in {0.05, 0.1}
#   - Taxa N in {20, 32, 64, 128}
#   - Divergence D in {0.5, 1.0, 2.0}
#   - 10 replicate shards per condition -> 2 x 12 x 10 = 240 parallel jobs
#
# By default, runs REPLICATES=30 (3 replicates per shard across 10 shards):
#   - Jobs   1..120 : indel_rate = 0.05 -> results_hpc_indel_0.05/ (reps 1..30)
#   - Jobs 121..240 : indel_rate = 0.1  -> results_hpc_indel_0.1/  (reps 1..30)
#
# To override replicates:
#   REPLICATES=50 bsub < scripts/submit_lsf.sh
#
# Usage:
#   mkdir -p logs
#   bsub < scripts/submit_lsf.sh
#
# After all jobs complete, merge shards and generate all plots (nRF, KF/BSD,
# Patristic RMSE, 4PC Penalty, SSS vs nRF curves) with:
#   python analysis/merge_and_plot_hpc.py --outdir results_hpc_indel_0.05
#   python analysis/merge_and_plot_hpc.py --outdir results_hpc_indel_0.1
# ==============================================================================

#BSUB -J add_phylo_indel[1-240]
#BSUB -n 1
#BSUB -M 4000
#BSUB -W 04:00
#BSUB -o logs/add_phylo_indel_%J_%I.out
#BSUB -e logs/add_phylo_indel_%J_%I.err

set -uo pipefail

# Activate micromamba environment
source ~/.bashrc 2>/dev/null || true
eval "$(micromamba shell hook --shell bash)" 2>/dev/null || true
micromamba activate phylomethod_env

# Move to repository root directory
cd "${LS_SUBCWD:-$(pwd)}"

mkdir -p logs

# Ensure INDELible v1.03 binary is available and compiled for this host OS/arch
bash scripts/install_indelible.sh

# 2 Indel rates x 12 (N, D) combinations = 24 conditions
INDEL_LIST=(0.05 0.1)
TAXA_LIST=(20 20 20 32 32 32 64 64 64 128 128 128)
DIV_LIST=(0.5 1.0 2.0 0.5 1.0 2.0 0.5 1.0 2.0 0.5 1.0 2.0)

N_COND_PER_INDEL=12
N_SHARDS=10
JOBS_PER_INDEL=$(( N_COND_PER_INDEL * N_SHARDS ))  # 120
TOTAL_REPS=${REPLICATES:-30}

JOB_IDX=$(( ${LSB_JOBINDEX:-1} - 1 ))
INDEL_IDX=$(( JOB_IDX / JOBS_PER_INDEL ))
REM_IDX=$(( JOB_IDX % JOBS_PER_INDEL ))
COND_IDX=$(( REM_IDX / N_SHARDS ))
SHARD_IDX=$(( REM_IDX % N_SHARDS ))

INDEL_RATE=${INDEL_LIST[$INDEL_IDX]}
N_TAXA=${TAXA_LIST[$COND_IDX]}
DIV=${DIV_LIST[$COND_IDX]}
OUTDIR="results_hpc_indel_${INDEL_RATE}"

mkdir -p "${OUTDIR}"

# Calculate replicate range [REP_START, REP_END] for this shard
CHUNK_SIZE=$(( (TOTAL_REPS + N_SHARDS - 1) / N_SHARDS ))
REP_START=$(( SHARD_IDX * CHUNK_SIZE + 1 ))
REP_END=$(( (SHARD_IDX + 1) * CHUNK_SIZE ))
if [ "${REP_END}" -gt "${TOTAL_REPS}" ]; then
    REP_END=${TOTAL_REPS}
fi

if [ "${REP_START}" -gt "${TOTAL_REPS}" ]; then
    echo "Shard ${SHARD_IDX} has no replicates to run (REP_START=${REP_START} > ${TOTAL_REPS}). Exiting."
    exit 0
fi

export OMP_NUM_THREADS=1

echo "============================================================"
echo "LSF Job Array Index ${LSB_JOBINDEX:-1}/240 | Host: $(hostname)"
echo "Condition: Indel=${INDEL_RATE} (POW 1.7 50), N=${N_TAXA}, D=${DIV}"
echo "Replicates: ${REP_START}..${REP_END} (of ${TOTAL_REPS}) | Outdir: ${OUTDIR}"
echo "Start Time: $(date)"
echo "============================================================"

python -u benchmark_hpc.py \
  --taxa "${N_TAXA}" \
  --divergence "${DIV}" \
  --indel-rate "${INDEL_RATE}" \
  --rep-start "${REP_START}" \
  --rep-end "${REP_END}" \
  --outdir "${OUTDIR}"

echo "============================================================"
echo "Finished LSF Job ${LSB_JOBINDEX:-1}: Indel=${INDEL_RATE}, N=${N_TAXA}, D=${DIV}, Reps=${REP_START}..${REP_END}"
echo "End Time: $(date)"
echo "============================================================"
