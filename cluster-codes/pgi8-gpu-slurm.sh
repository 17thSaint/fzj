#!/bin/bash
#SBATCH --job-name=gpu_withexpander_16x8
#SBATCH --output=log_file
#SBATCH --error=log_file
#SBATCH -p pgi-8-gpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=128G
#SBATCH --gres=gpu:h100:1

srun julia long-range-ttn.jl "layers" 7 "particles" 8 "onsite_strength" 2.0 "mdim" 500 
