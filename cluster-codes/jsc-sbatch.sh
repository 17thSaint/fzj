#!/bin/bash -x

echo "Date:"
date
echo "Nodes used:"
echo $SLURM_NODELIST

#module load CUDA/11.5 cuDNN/8.3.1.22-CUDA-11.5 Intel/2021.4.0 ParaStationMPI mpi4py h5py
#module load Julia

#source # !! HIER DEIN ENVIRONMENT AKTIVIEREN 
#export XLA_PYTHON_CLIENT_PREALLOCATE=false
#export XLA_FLAGS="--xla_gpu_cuda_data_dir=$CUDA_HOME"
#export CUDA_VISIBLE_DEVICES=0,1,2,3


# Get the script name, parameter name, start value, end value, and step size
script_name=$1
param=$2
START_VALUE=$3
END_VALUE=$4
STEP_SIZE=$5
additional_params=("${@:6}")

# Calculate the number of iterations
num_iters=$(echo "scale=2; ($END_VALUE - $START_VALUE) / $STEP_SIZE + 1" | bc)

# Get the current date and time
alpha=$(date +"%Y-%m-%d_%H-%M")
datafolder="/p/project/netenesyquma/geraghty1/data/data-$alpha"

# Create the folder
mkdir -p "$datafolder"

#SBATCH --array=1-$num_iters
#SBATCH --time=00:15:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
##SBATCH --mem=1GB
#SBATCH --account=netenesyquma


# Loop through the range of input parameter values
value=$START_VALUE
while (( $(bc <<< "$value <= $END_VALUE") )); do
    # Submit the specified script as a SLURM job
    sbatch -A netenesyquma run-script-jsc.sh "$script_name" "open_cores" 4 "dataloc" "$datafolder" "$param" "$value" "${additional_params[@]}"

    value=$(bc <<< "$value + $STEP_SIZE")
done

# Wait for all background jobs to finish
wait


