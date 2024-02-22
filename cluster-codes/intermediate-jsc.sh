#!/bin/bash

script_name=$1
param=$2
START_VALUE=$3
STEP_SIZE=$4
datafolder=$5
additional_params=("${@:6}")

value=$(bc -l <<< "($SLURM_PROCID) * $STEP_SIZE + $START_VALUE" )
echo "$value"

./run-script-jsc.sh "$script_name" "dataloc" "$datafolder" "$param" "$value" "${additional_params[@]}"

