#!/bin/bash
#SBATCH --job-name=10x8_ulrxi_farm
#SBATCH --output=logs/log_file_10x8_ulrxi_4p0_%a
#SBATCH --error=logs/log_file_10x8_ulrxi_4p0_%a
#SBATCH -p pgi-8-gpu-h100
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=60G
#SBATCH --gres=gpu:h100:1
#SBATCH --array=5-10  # Adjust based on number of values


# Map array indices to onsite_strength values
farm_values=(0.0 0.8 1.6 2.4 3.2 4.0 4.8 5.6 6.4 7.2 8.0)
local_farm_value=${farm_values[$SLURM_ARRAY_TASK_ID]}

# Run the Julia script with the desired parameters
srun julia long-range-ttn.jl "Lx" 10 "Ly" 8 "particles" 5 "onsite_strength" 4.0 "corr_length" "$local_farm_value"
