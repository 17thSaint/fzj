#!/bin/bash
#SBATCH --job-name=j1_8x8_gpu
#SBATCH --output=log_file_j1_8x8_gpu
#SBATCH --error=log_file_j1_8x8_gpu
#SBATCH -p pgi-8-gpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=80G
#SBATCH --gres=gpu:a100:1

srun julia j1j2.jl "Lx" 8 "Ly" 8 "j2" 0.0 "mdim" 300 "if_gpu" true
