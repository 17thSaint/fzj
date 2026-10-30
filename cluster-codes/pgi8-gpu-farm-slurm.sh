#!/bin/bash
#SBATCH --job-name=14x7_dd
#SBATCH --output=logs/log_file_14x7_dd_small_%a
#SBATCH --error=logs/log_file_14x7_dd_small_%a
#SBATCH -p pgi-8-gpu-h100
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=100G
#SBATCH --gres=gpu:h100:1
#SBATCH --array=0-1  # Adjust based on number of values


# Map array indices to onsite_strength values
farm_values=(0.7 0.8 1.1 1.53 2.1 2.9 4.0 5.5 7.6 10.5 14.5)
local_farm_value=${farm_values[$SLURM_ARRAY_TASK_ID]}

# Run the Julia script with the desired parameters
srun julia long-range-ttn.jl "Lx" 14 "Ly" 7 "particles" 7 "onsite_strength" 300.0 "magnetic_spacing" "$local_farm_value"
