#!/bin/bash
#SBATCH --job-name=8x8_dd_fpeps
#SBATCH --partition=pgi-8
#SBATCH --array=0-9
#SBATCH --ntasks=1           
#SBATCH --cpus-per-task=72
#SBATCH --mem=30G
#SBATCH --output=logs/log_file_8x8_dd_fpeps_%a
#SBATCH --error=logs/log_file_8x8_dd_fpeps_%a

farm_values=(1.126 0.968 0.832 0.715 0.614 0.528 0.453 0.390 0.335 0.288)
val=${farm_values[$SLURM_ARRAY_TASK_ID]}


julia execute-fpeps.jl "Lx" 8 "Ly" 8 "particles" 4 "magspac" $val
