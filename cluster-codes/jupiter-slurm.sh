#!/bin/bash
#SBATCH --account=syqma2tens
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --job-name=test_gpu_8x8_laughlin
#SBATCH --output=logs/log_file_test_8x8_laughlin
#SBATCH --error=logs/log_file_test_8x8_laughlin
#SBATCH --time=00:00:59
#SBATCH --partition=all


module load Julia

srun julia long-range-ttn.jl "Lx" 8 "Ly" 8 "particles" 4 "onsite_strength" 0.0
