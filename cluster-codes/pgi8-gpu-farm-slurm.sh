#!/bin/bash
#SBATCH --job-name=12x4anis_farm
#SBATCH --output=logs/log_file_12x4_hanis_%A_%a
#SBATCH --error=logs/log_file_12x4_hanis_%A_%a
#SBATCH -p pgi-8-gpu-h100
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=128G
#SBATCH --gres=gpu:h100:1
#SBATCH --array=0-8  # Adjust based on number of values


# Map array indices to onsite_strength values
farm_values=(1e-3 1e-2 1e-1 0.2 0.3 0.5 0.7 0.9 1.0)
local_farm_value=${farm_values[$SLURM_ARRAY_TASK_ID]}

# Set threads for Julia (optional)
export JULIA_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Run the Julia script with the desired parameters
srun julia long-range-ttn.jl "Lx" 12 "Ly" 4 "particles" 6 "onsite_strength" 300.0 "tx" "$local_farm_value"
