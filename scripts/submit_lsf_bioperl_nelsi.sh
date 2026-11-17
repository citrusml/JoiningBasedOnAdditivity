#!/bin/bash
# ==============================================================================
# High-Parallelism LSF (bsub) Job Array Script: BioPerl vs. NELSI Benchmark
#
# Evaluates:
#   - Tree Models: bioperl, nelsi (2 models)
#   - Taxa N:      20, 32, 64 (3 levels)
#   - Divergence D: 0.1, 0.2, 0.5, 1.0, 1.5, 2.0, 3.0 (7 levels)
#   - Length L:    500 aa
#   - Indel Rate:  0.1 (INDELible Zipfian POW 1.7 50)
#   - Replicates:  50 per condition (split into 10 shards -> 5 reps/job)
#
# Total LSF Array Jobs: 2 x 3 x 7 x 10 = 420 parallel jobs
#   - Jobs   1..210 : tree_model = bioperl
#   - Jobs 211..420 : tree_model = nelsi
#
# Usage:
#   mkdir -p logs
#   bsub < scripts/submit_lsf_bioperl_nelsi.sh
#
# To override replicates (e.g., 10):
#   REPLICATES=10 bsub < scripts/submit_lsf_bioperl_nelsi.sh
#
# After all 420 jobs complete, merge all shards into final CSVs:
#   python benchmark_bioperl_nelsi.py --merge-only --outdir results_bioperl_nelsi
# ==============================================================================

#BSUB -J bioperl_nelsi[1-420]
#BSUB -n 1
#BSUB -M 4000
#BSUB -W 04:00
#BSUB -o logs/bioperl_nelsi_%J_%I.out
#BSUB -e logs/bioperl_nelsi_%J_%I.err

set -uo pipefail

# Activate micromamba environment
source ~/.bashrc 2>/dev/null || true
eval "$(micromamba shell hook --shell bash)" 2>/dev/null || true
micromamba activate phylomethod_env

# Move to repository root directory
cd "${LS_SUBCWD:-$(pwd)}"

mkdir -p logs

# Ensure INDELible v1.03 binary is compiled for this host OS/arch
bash scripts/install_indelible.sh

MODEL_LIST=(bioperl nelsi)
TAXA_LIST=(20 32 64)
DIV_LIST=(0.1 0.2 0.5 1.0 1.5 2.0 3.0)

N_TAXA_OPTS=${#TAXA_LIST[@]}     # 3
N_DIV_OPTS=${#DIV_LIST[@]}       # 7
N_COND_PER_MODEL=$(( N_TAXA_OPTS * N_DIV_OPTS ))  # 21
N_SHARDS=10
JOBS_PER_MODEL=$(( N_COND_PER_MODEL * N_SHARDS )) # 210

TOTAL_REPS=${REPLICATES:-50}
SEQ_LEN=${SEQ_LENGTH:-500}
INDEL_RATE=${INDEL_RATE:-0.1}
OUTDIR=${OUTDIR:-results_bioperl_nelsi}

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

mkdir -p "${OUTDIR}/shards"

export OMP_NUM_THREADS=1

echo "============================================================"
echo "LSF Job Array Index ${LSB_JOBINDEX:-1}/420 | Host: $(hostname)"
echo "Model: ${TREE_MODEL} | N=${N_TAXA} | D=${DIV} | L=${SEQ_LEN} | Indel=${INDEL_RATE}"
echo "Shard: ${SHARD_IDX}/${N_SHARDS} (of ${TOTAL_REPS} reps) | Outdir: ${OUTDIR}"
echo "Start Time: $(date)"
echo "============================================================"

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
  --outdir "${OUTDIR}"

echo "============================================================"
echo "Finished LSF Job ${LSB_JOBINDEX:-1}: ${TREE_MODEL}, N=${N_TAXA}, D=${DIV}, Shard=${SHARD_IDX}"
echo "End Time: $(date)"
echo "============================================================"
