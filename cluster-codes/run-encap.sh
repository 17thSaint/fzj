#!/usr/bin/env bash

script_name=$1
param=$2
START_VALUE=$3
END_VALUE=$4
num_tasks=$5
STEP_SIZE=$(bc -l <<< "($END_VALUE - $START_VALUE) / ($farm_iters - 1)")
exp_name=$6
additional_params=("${@:7}")

#encap run "$script_name" -args "$param","$START_VALUE","$STEP_SIZE","${additional_params[@]}" -i $farm_iters -n "$exp_name"

num_cores_pertask=1
num_cores_pernode=48

# Check if all required variables are provided
if [ -z "$num_cores_pertask" ] || [ -z "$num_tasks" ] || [ -z "$num_cores_pernode" ]; then
  echo "Please provide values for num_cores_pertask, num_tasks, and num_cores_pernode."
  exit 1
fi

# Calculate the total number of cores needed to run all tasks
total_cores_needed=$((num_cores_pertask * num_tasks))

# Calculate the number of nodes required to evenly distribute tasks
nodes_needed=$((total_cores_needed / num_cores_pernode))

# If there are any remaining tasks, we need an extra node
if [ $((total_cores_needed % num_cores_pernode)) -ne 0 ]; then
  nodes_needed=$((nodes_needed + 1))
fi

tasks_per_node=$((num_tasks / nodes_needed))

echo "To run $num_tasks tasks with $num_cores_pertask cores per task on nodes with $num_cores_pernode cores each, you need $nodes_needed nodes with $tasks_per_node tasks per node."

encap run "$script_name" -args "$param","$START_VALUE","$STEP_SIZE","${additional_params[@]}" -n "$exp_name" -sl_i $nodes_needed -sl_ntpn $tasks_per_node

wait

