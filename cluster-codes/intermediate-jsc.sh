#!/bin/bash

param=$1
START_VALUE=$2
STEP_SIZE=$3
datafolder=$4
additional_params=("${@:5}")

value=$(( ( $SLURM_ARRAY_TASK_ID - 1 ) * $STEP_SIZE + $START_VALUE ))
echo "$value"

./run-script-jsc.sh "$script_name" "open_cores" 4 "dataloc" "$datafolder" "$param" "$value" "${additional_params[@]}"

