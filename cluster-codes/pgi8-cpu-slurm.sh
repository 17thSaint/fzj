#!/bin/bash
#SBATCH --job-name=int_0p0_xiSD
#SBATCH --partition=pgi-8
#SBATCH --ntasks=1             
#SBATCH --cpus-per-task=5
#SBATCH --mem=20G
#SBATCH --output=logs/log_file_12x6_laugh
#SBATCH --error=logs/log_file_12x6_laugh


srun julia long-range-ttn.jl "Lx" 12 "Ly" 6 "particles" 6 "onsite_strength" 0.0 "pinning_strength" 0.001