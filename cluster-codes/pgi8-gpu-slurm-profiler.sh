#!/bin/bash
#SBATCH --job-name=nsighttest-8x4-ulr
#SBATCH --output=logs/log_file_8x4_nsighttest_ulr_h100
#SBATCH --error=logs/log_file_8x4_nsighttest_ulr_h100
#SBATCH -p pgi-8-gpu-h100
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=40G
#SBATCH --gres=gpu:h100:1

export PATH=/opt/nsight-2025.6.1/bin:$PATH

# (Optional but often needed for NVTX-triggered capture)
export NSYS_NVTX_PROFILER_REGISTER_ONLY=0

out="$(nvidia-smi)"
printf '%s\n' "$out"

srun nsys profile \
  -t cuda,nvtx,osrt \
  --capture-range=nvtx \
  --capture-range-end=stop \
  --force-overwrite=true \
  -o "logs/nsighttest_8x4_h100_ulr_${SLURM_JOB_ID}" \
  julia --project=. daily-things.jl

