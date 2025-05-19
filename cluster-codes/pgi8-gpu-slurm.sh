#!/bin/bash
#SBATCH --job-name=gpu_12x6_pinning
#SBATCH --output=log_file_12x6_pinning
#SBATCH --error=log_file_12x6_pinning
#SBATCH -p pgi-8-gpu-h100
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=128G
#SBATCH --gres=gpu:h100:1

srun julia long-range-ttn.jl "Lx" 12 "Ly" 6 "particles" 6 "onsite_strength" 0.0 "mdim" 400 "pinning_strength" 0.1
