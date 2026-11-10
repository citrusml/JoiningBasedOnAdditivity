#!/bin/bash
# ==============================================================================
# High-Parallelism LSF (bsub) Job Array Script (120 Parallel Jobs)
#
# Cluster capacity: ~232 free slots across mafftsrv1..mafftsrv5 & mafftweb.
# This script splits the 12 (N, D) conditions into 10 replicate shards each:
#   12 conditions x 10 shards = 120 parallel jobs (#BSUB -J add_phylo[1-120]).
#
# By default, runs REPLICATES=100 (10 replicates per job -> finishes in ~5-8 mins!).
# To run 30 replicates instead (3 replicates per job -> finishes in ~2 mins):
#   REPLICATES=30 bsub < scripts/submit_lsf.sh
#
# Usage:
#   mkdir -p logs
#   bsub < scripts/submit_lsf.sh
#
# After all jobs complete, merge shards and generate all plots with:
#   python analysis/merge_and_plot_hpc.py --outdir results_hpc
# ==============================================================================

#BSUB -J add_phylo[1-120]
#BSUB -n 1
#BSUB -M 4000
#BSUB -W 04:00
#BSUB -o logs/add_phylo_%J_%I.out
#BSUB -e logs/add_phylo_%J_%I.err

set -uo pipefail

# Activate micromamba environment
source ~/.bashrc 2>/dev/null || true
eval "$(micromamba shell hook --shell bash)" 2>/dev/null || true
micromamba activate phylomethod_env

# Move to repository root directory
cd "${LS_SUBCWD:-$(pwd)}"

mkdir -p logs results_hpc

# 12 (N, D) combinations
TAXA_LIST=(20 20 20 32 32 32 64 64 64 128 128 128)
DIV_LIST=(0.5 1.0 2.0 0.5 1.0 2.0 0.5 1.0 2.0 0.5 1.0 2.0)

# Number of shards per (N, D) condition = 10 (total 12 x 10 = 120 jobs)
N_SHARDS=10
TOTAL_REPS=${REPLICATES:-100}

JOB_IDX=$(( ${LSB_JOBINDEX:-1} - 1 ))
COND_IDX=$(( JOB_IDX / N_SHARDS ))
SHARD_IDX=$(( JOB_IDX % N_SHARDS ))

N_TAXA=${TAXA_LIST[$COND_IDX]}
DIV=${DIV_LIST[$COND_IDX]}

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
echo "LSF Job Array Index ${LSB_JOBINDEX:-1}/120 | Host: $(hostname)"
echo "Condition: N=${N_TAXA}, D=${DIV} | Replicates: ${REP_START}..${REP_END} (of ${TOTAL_REPS})"
echo "Start Time: $(date)"
echo "============================================================"

python -u benchmark_hpc.py \
  --taxa "${N_TAXA}" \
  --divergence "${DIV}" \
  --rep-start "${REP_START}" \
  --rep-end "${REP_END}" \
  --outdir results_hpc

echo "============================================================"
echo "Finished LSF Job ${LSB_JOBINDEX:-1}: N=${N_TAXA}, D=${DIV}, Reps=${REP_START}..${REP_END}"
echo "End Time: $(date)"
echo "============================================================"
