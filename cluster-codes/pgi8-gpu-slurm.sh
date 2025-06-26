#!/bin/bash
#SBATCH --job-name=ulr1p0_14x4_gpu
#SBATCH --output=log_file_14x4_ulr_1p0
#SBATCH --error=log_file_14x4_ulr_1p0
#SBATCH -p pgi-8-gpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=64G
#SBATCH --gres=gpu:a100:1

srun julia long-range-ttn.jl "Lx" 14 "Ly" 4 "particles" 7 "onsite_strength" 1.0
