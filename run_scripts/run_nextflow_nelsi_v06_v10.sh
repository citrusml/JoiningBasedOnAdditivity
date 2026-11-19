#!/bin/bash
#BSUB -J nextflow_nelsi_v06_v10
#BSUB -q mafft
#BSUB -n 1
#BSUB -M 16GB
#BSUB -W 120:00
#BSUB -o logs/nextflow_nelsi_v06_v10_bsub.out
#BSUB -e logs/nextflow_nelsi_v06_v10_bsub.err

export PATH="$HOME/bin:$PATH"
eval "$(micromamba shell hook --shell bash)" 2>/dev/null || true
micromamba activate phylomethod_env
export PATH="${CONDA_PREFIX:-}:$CONDA_PREFIX/bin:$HOME/bin:$PATH"

# NELSI 用 Rscript のパスを指定
export NELSI_RSCRIPT="${NELSI_RSCRIPT:-${MAMBA_ROOT_PREFIX:-$HOME/.micromamba}/envs/nelsi_env/bin/Rscript}"

cd "${LS_SUBCWD:-$(cd "$(dirname "$0")/.." && pwd)}"
mkdir -p logs results_nelsi_v06_v10

bash scripts/install_indelible.sh
export NXF_OPTS="-Xms2g -Xmx12g"

nextflow run main.nf \
    -c next_configs/nextflow_nelsi_v06_v10.config \
    -profile supercomputer \
    -resume \
    > logs/nextflow_nelsi_v06_v10.log 2>&1
