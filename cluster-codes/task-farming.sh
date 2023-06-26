#!/bin/bash

script_name=$1
param=$2
START_VALUE=$3
END_VALUE=$4
STEP_SIZE=$5

shift 5
additional_params=("$@")

# Loop through the range of input parameter values
value=$START_VALUE
while (( $(bc <<< "$value <= $END_VALUE") )); do
    # Call the specified script with the specific input parameter value and additional parameters in the background
    ./run-script.sh "$script_name" "$param" "$value" "${additional_params[@]}" &
    value=$(bc <<< "$value + $STEP_SIZE")
done


# Wait for all background processes to finish
wait
