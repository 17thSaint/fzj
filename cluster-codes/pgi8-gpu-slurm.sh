#!/bin/bash
#SBATCH --job-name=14x7_pinned_ulr
#SBATCH --output=logs/log_file_14x7_pinned_ulr
#SBATCH --error=logs/log_file_14x7_pinned_ulr
#SBATCH -p pgi-8-gpu-h100
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=80G
#SBATCH --gres=gpu:h100:1

srun julia long-range-ttn.jl
