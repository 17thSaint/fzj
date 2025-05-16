#!/bin/bash
#SBATCH --job-name=gpu_withexpander_16x8
#SBATCH --output=log_file
#SBATCH --error=log_file
#SBATCH -p pgi-8-gpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=128G
#SBATCH --gres=gpu:h100:1

srun julia long-range-ttn.jl "Lx" 10 "Ly" 5 "particles" 5 "onsite_strength" 0.0 "mdim" 300 "pinning_strength" 0.1 
