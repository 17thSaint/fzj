#!/bin/bash -x

#SBATCH --account=netenesyquma
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --ntasks-per-node=10
#SBATCH --time=00:10:00
#SBATCH --partition=batch

echo "Date:"
date
echo "Nodes used:"

export SRUN_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK

#module load CUDA/11.5 cuDNN/8.3.1.22-CUDA-11.5 Intel/2021.4.0 ParaStationMPI mpi4py h5py

module load Julia

#source # !! HIER DEIN ENVIRONMENT AKTIVIEREN 
#export XLA_PYTHON_CLIENT_PREALLOCATE=false
#export XLA_FLAGS="--xla_gpu_cuda_data_dir=$CUDA_HOME"
#export CUDA_VISIBLE_DEVICES=0,1,2,3

# Get the script name, parameter name, start value, end value, and step size
script_name=$1
param=$2
START_VALUE=$3
END_VALUE=$4
STEP_SIZE=$(bc -l <<< "($END_VALUE - $START_VALUE) / ($SLURM_NTASKS_PER_NODE * $SLURM_JOB_NUM_NODES - 1)")
additional_params=("${@:5}")

# Calculate the number of iterations
#num_iters=$(echo "scale=2; ($END_VALUE - $START_VALUE) / $STEP_SIZE + 1" | bc)

# Get the current date and time
#alpha=$(date +"%Y-%m-%d_%H-%M")
datafolder="/p/project/netenesyquma/geraghty1/data/data-$alpha"

# Create the folder
#mkdir -p "$datafolder"


#value=$(( ( $SLURM_ARRAY_TASK_ID - 1 ) * $STEP_SIZE + $START_VALUE ))
#echo "$value"

srun intermediate-jsc.sh "$script_name" "$param" "$START_VALUE" "$STEP_SIZE" "$datafolder" "${additional_params[@]}"

wait
