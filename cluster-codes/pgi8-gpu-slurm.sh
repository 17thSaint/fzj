#!/bin/bash
#SBATCH --job-name=hh_singlethreadcpu_ongpu
#SBATCH --output=log_file
#SBATCH --error=log_file
#SBATCH -p pgi-8-gpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --gres=gpu:a100:1

srun julia gpu-benchmarking.jl "benchmark_type" "gpu" "model" "HH" "min_mdim" 50 "max_mdim" 1000 "count_mdim" 50 "lx" 8
