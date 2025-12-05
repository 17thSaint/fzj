#!/bin/bash
#SBATCH --job-name=ulrxi_16x8_farm
#SBATCH --output=logs/log_file_16x8_ulrxi_0p4_%a
#SBATCH --error=logs/log_file_16x8_ulrxi_0p4_%a
#SBATCH -p pgi-8-gpu-h100
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=128G
#SBATCH --gres=gpu:h100:1
#SBATCH --array=0-1  # Adjust based on number of values


# Map array indices to onsite_strength values
farm_values=(0.8 1.6)
local_farm_value=${farm_values[$SLURM_ARRAY_TASK_ID]}

# Run the Julia script with the desired parameters
srun julia long-range-ttn.jl "Lx" 16 "Ly" 8 "particles" 8 "onsite_strength" 0.4 "corr_length" "$local_farm_value"
