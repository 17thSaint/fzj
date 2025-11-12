#!/bin/bash
#SBATCH --job-name=int_0p0_xiSD
#SBATCH --partition=pgi-8
#SBATCH --array=0-10                
#SBATCH --ntasks=1             
#SBATCH --cpus-per-task=5
#SBATCH --mem=10G
#SBATCH --output=logs/log_file_8x4_0p0_%A_%a.out
#SBATCH --error=logs/log_file_8x4_0p0_%A_%a.err

# 11 farmed values (index = SLURM_ARRAY_TASK_ID)
farm_values=(0.0 0.4 0.8 1.2 1.6 2.0 2.4 2.8 3.2 3.6 4.0)
val=${farm_values[$SLURM_ARRAY_TASK_ID]}


julia execute-ed.jl "xi" "$val" "onsite_strength" 0.0
