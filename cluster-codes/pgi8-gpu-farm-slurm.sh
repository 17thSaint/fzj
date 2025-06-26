#!/bin/bash
#SBATCH --job-name=12x6ulr_farm
#SBATCH --output=logs/log_file_12x6_ulr_%A_%a
#SBATCH --error=logs/log_file_12x6_ulr_%A_%a
#SBATCH -p pgi-8-gpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=64G
#SBATCH --gres=gpu:a100:1
#SBATCH --array=0-7  # Adjust based on number of values


# Map array indices to onsite_strength values
farm_values=(1.0 5.0 6.0 7.0 8.0 9.0 10.0 300.0)
local_farm_value=${farm_values[$SLURM_ARRAY_TASK_ID]}

# Set threads for Julia (optional)
export JULIA_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Run the Julia script with the desired parameters
srun julia long-range-ttn.jl "Lx" 12 "Ly" 6 "particles" 6 "onsite_strength" "$local_farm_value"
