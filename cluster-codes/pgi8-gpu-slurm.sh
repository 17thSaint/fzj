#!/bin/bash
#SBATCH --job-name=j1_4x4_gpu
#SBATCH --output=log_file_j1_4x4_gpu
#SBATCH --error=log_file_j1_4x4_gpu
#SBATCH -p pgi-8-gpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=64G
#SBATCH --gres=gpu:a100:1

srun julia j1j2.jl "Lx" 4 "Ly" 4 "j2" 0.0 "mdim" 300 "if_gpu" true
