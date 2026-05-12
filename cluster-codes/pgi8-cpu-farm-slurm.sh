#!/bin/bash
#SBATCH --job-name=dtbenchmark_8x4
#SBATCH --partition=pgi-8
#SBATCH --array=0-17
#SBATCH --ntasks=1           
#SBATCH --cpus-per-task=5
#SBATCH --mem=10G
#SBATCH --output=logs/log_file_8x4_dtbenchmark_smaller_%a
#SBATCH --error=logs/log_file_8x4_dtbenchmark_smaller_%a

farm_values=(0.001 0.0015 0.002 0.0025 0.003 0.0035 0.004 0.0045 0.005 0.0055 0.006 0.0065 0.007 0.0075 0.008 0.0085 0.009 0.0095)
val=${farm_values[$SLURM_ARRAY_TASK_ID]}


julia execute-ed.jl "Lx" 8 "Ly" 4 "particles" 4 "dt" $val
