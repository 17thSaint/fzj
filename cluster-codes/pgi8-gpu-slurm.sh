#!/bin/bash
#SBATCH --job-name=pinning_14x7_pinning_ulr_star
#SBATCH --output=logs/log_file_14x7_pinning_ulr_star
#SBATCH --error=logs/log_file_14x7_pinning_ulr_star
#SBATCH -p pgi-8-gpu-h100
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=50G
#SBATCH --gres=gpu:h100:1

export PATH=/opt/nsight-2025.6.1/bin:$PATH

srun julia long-range-ttn.jl "Lx" 14 "Ly" 7 "particles" 7 "onsite_strength" 300.0 "pinning_strength" 0.0001
