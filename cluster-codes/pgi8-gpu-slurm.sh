#!/bin/bash
#SBATCH --job-name=j1j2_gputest
#SBATCH --output=log_file
#SBATCH --error=log_file
#SBATCH -p pgi-8-gpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --gres=gpu:a100:1

srun julia gpu-benchmarking.jl "benchmark_type" "gpu" "model" "HH" "min_mdim" 20 "max_mdim" 50 "count_mdim" 2 "lx" 8
