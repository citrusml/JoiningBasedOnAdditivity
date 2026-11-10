#!/bin/bash
# ==============================================================================
# LSF (bsub) Job Array Script for Additivity-Based Phylogenetic Benchmark
#
# Submits 12 parallel jobs corresponding to:
#   Taxa N in {20, 32, 64, 128} x Divergence D in {0.5, 1.0, 2.0}
#   Each job runs 30 replicates across 3 distance matrix types:
#     1. BLOSUM62 + Poisson distance
#     2. WAG + Gamma Maximum Likelihood distance
#     3. True Patristic distance
#   plus MSA + ML (MAFFT + FastTree), saving all Newick trees, SSS statistics,
#   and 4-point condition violation scores.
#
# Usage:
#   mkdir -p logs
#   bsub < scripts/submit_lsf.sh
#
# After all 12 jobs complete, merge and plot results with:
#   python analysis/merge_and_plot_hpc.py --outdir results_hpc
# ==============================================================================

#BSUB -J add_phylo[1-12]
#BSUB -n 4
#BSUB -M 16000
#BSUB -W 12:00
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

# 12 (N, D) combinations mapped to LSB_JOBINDEX = 1..12
TAXA_LIST=(20 20 20 32 32 32 64 64 64 128 128 128)
DIV_LIST=(0.5 1.0 2.0 0.5 1.0 2.0 0.5 1.0 2.0 0.5 1.0 2.0)

ARRAY_IDX=$(( ${LSB_JOBINDEX:-1} - 1 ))
N_TAXA=${TAXA_LIST[$ARRAY_IDX]}
DIV=${DIV_LIST[$ARRAY_IDX]}

REPLICATES=${REPLICATES:-30}
export OMP_NUM_THREADS=4

echo "============================================================"
echo "Starting LSF Job Array Index ${LSB_JOBINDEX:-1}: N=${N_TAXA}, D=${DIV}, Replicates=${REPLICATES}"
echo "Host: $(hostname) | Date: $(date)"
echo "============================================================"

python -u benchmark_hpc.py \
  --taxa "${N_TAXA}" \
  --divergence "${DIV}" \
  --replicates "${REPLICATES}" \
  --outdir results_hpc

echo "============================================================"
echo "Finished LSF Job Array Index ${LSB_JOBINDEX:-1}: N=${N_TAXA}, D=${DIV}"
echo "Date: $(date)"
echo "============================================================"
