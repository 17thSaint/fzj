#!/bin/bash

script_name=$1
num_vers=$2
param=$3
START_VALUE=$4
END_VALUE=$5
STEP_SIZE=$6

num_iters=$(echo "scale=2; ($END_VALUE - $START_VALUE) / $STEP_SIZE + 1" | bc)
num_runs=$(echo "$num_iters * $num_vers" | bc)
num_cores_available=$(nproc)
open_cores=$(echo "scale=0; ($num_cores_available - 3) / $num_runs " | bc)
echo "Using $open_cores for each of the $num_runs Runs"

shift 6
additional_params=("$@")

# Loop through the range of input parameter values
value=$START_VALUE
while (( $(bc <<< "$value <= $END_VALUE") )); do
    # Call the specified script with the specific input parameter value and additional parameters in the background
    ./run-script.sh "$script_name" "open_cores" "$open_cores" "$param" "$value" "${additional_params[@]}" &
    value=$(bc <<< "$value + $STEP_SIZE")
done


# Wait for all background processes to finish
wait
