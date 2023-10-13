#!/bin/bash

SCRIPTNAME="$1"

declare -A dictionary

for ((i=4; i<=$#; i+=2)); do
    key="${!i}"
    eval "value=\${$((i+1))}"
    dictionary["$key"]="$value"
done

JOBNAME=""
for key in "${!dictionary[@]}"; do
    value="${dictionary[$key]}"
    JOBNAME+="$key-$value-"
done

JOBNAME=${JOBNAME%-}

OUTPUTFILENAME="../cluster-data/output-$JOBNAME.txt"

# written by ChatGPT 15.06.2023
extension="${OUTPUTFILENAME##*.}"
basename="${OUTPUTFILENAME%.*}"
iter=1
while [ -e "$OUTPUTFILENAME" ]
do
        new_basename="${basename}-mk-$iter"
        OUTPUTFILENAME="${new_basename}.${extension}"
        iter=$((iter + 1))
done

echo "$@"

julia "$@" #> "${OUTPUTFILENAME}"

