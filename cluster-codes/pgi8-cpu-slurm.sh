#!/bin/bash
#SBATCH --job-name=16x8_flat
#SBATCH --output=logs/log_file_16x8_flat_50
#SBATCH --error=logs/log_file_16x8_flat_50
#SBATCH --partition=pgi-8
#SBATCH --ntasks=1             
#SBATCH --cpus-per-task=72
#SBATCH --mem=100G


srun julia execute-fpeps.jl
