#!/usr/bin/env bash

script_name=$1
param=$2
START_VALUE=$3
END_VALUE=$4
farm_iters=$5
STEP_SIZE=$(bc -l <<< "($END_VALUE - $START_VALUE) / ($farm_iters - 1)")
exp_name=$6
additional_params=("${@:7}")

encap run "$script_name" -args "$param","$START_VALUE","$STEP_SIZE","${additional_params[@]}" -i $farm_iters -n "$exp_name"

wait

