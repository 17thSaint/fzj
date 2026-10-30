#!/bin/bash
#SBATCH --job-name=16x8_dd_0p55
#SBATCH --output=logs/log_file_16x8_dd_0p55
#SBATCH --error=logs/log_file_16x8_dd_0p55
#SBATCH -p pgi-8-gpu-h100
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=100G
#SBATCH --gres=gpu:h100:1

srun julia long-range-ttn.jl "Lx" 16 "Ly" 8 "particles" 8 "onsite_strength" 300.0 "magnetic_spacing" 0.55
