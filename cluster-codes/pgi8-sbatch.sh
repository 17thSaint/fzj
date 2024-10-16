#!/bin/bash

#SBATCH --job-name=testing-ttn-cluster
#SBATCH --nodes=5
#SBATCH --cpus-per-task=5
#SBATCH --ntasks-per-node=1
#SBATCH --time=00:10:00

export SRUN_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK

# Get the script name, parameter name, start value, end value, and step size
script_name=$1
param=$2
START_VALUE=$3
END_VALUE=$4

if (( $(bc -l <<< "$START_VALUE == $END_VALUE") )); then
    STEP_SIZE=0
else
    STEP_SIZE=$(bc -l <<< "($END_VALUE - $START_VALUE) / ($SLURM_NTASKS_PER_NODE * $SLURM_JOB_NUM_NODES - 1)")
fi

additional_params=("${@:5}")

# Calculate the number of iterations
#num_iters=$(echo "scale=2; ($END_VALUE - $START_VALUE) / $STEP_SIZE + 1" | bc)

# Get the current date and time
#alpha=$(date +"%Y-%m-%d_%H-%M")
datafolder="/local/geraghty/cluster-data/synth-dims/excited-states"


#value=$(( ( $SLURM_ARRAY_TASK_ID - 1 ) * $STEP_SIZE + $START_VALUE ))
#echo "$value"

srun intermediate-jsc.sh "$script_name" "$param" "$START_VALUE" "$STEP_SIZE" "$datafolder" "${additional_params[@]}"

wait