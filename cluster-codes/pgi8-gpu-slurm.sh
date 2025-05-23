#!/bin/bash
#SBATCH --job-name=gpu_16x8_ulr_range
#SBATCH --output=log_file_16x8_ulr_range
#SBATCH --error=log_file_16x8_ulr_range
#SBATCH -p pgi-8-gpu-h100
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=128G
#SBATCH --gres=gpu:h100:1

srun julia long-range-ttn.jl
