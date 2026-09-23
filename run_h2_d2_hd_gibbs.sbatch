#!/bin/bash
#SBATCH --job-name=h2d2hd-gibbs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --nodelist=node4
#SBATCH --time=01:00:00
#SBATCH --output=h2d2hd-gibbs-%j.out
#SBATCH --error=h2d2hd-gibbs-%j.err

set -euo pipefail
cd /home/acc-2026-27422/practice
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export MAESTRO_ALLOW_RAW=1
maestro-python h2_d2_hd_gibbs.py
