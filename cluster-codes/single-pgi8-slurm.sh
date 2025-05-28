#!/bin/bash
#SBATCH --job-name=ttn_14x7_laughlin_pinned
#SBATCH --output=log_file_ttn_14x7_laughlin_pinned
#SBATCH --error=log_file_ttn_14x7_laughlin_pinned
#SBATCH --nodes=1
#SBATCH --cpus-per-task=5
#SBATCH --ntasks-per-node=1

srun julia long-range-ttn.jl "Lx" 14 "Ly" 7 "particles" 7 "onsite_strength" 0.0 "mdim" 400 "pinning_strength" 0.1