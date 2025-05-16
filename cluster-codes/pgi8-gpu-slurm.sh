#!/bin/bash
#SBATCH --job-name=gpu_10x5_pinning_ulr
#SBATCH --output=log_file_10x5_pinning_ulr
#SBATCH --error=log_file_10x5_pinning_ulr
#SBATCH -p pgi-8-gpu-h100
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=128G
#SBATCH --gres=gpu:h100:1

srun julia long-range-ttn.jl "Lx" 10 "Ly" 5 "particles" 5 "onsite_strength" 300.0 "mdim" 300 "pinning_strength" 0.001
