#!/bin/bash

node_name="$1"

for folder in "review" "synthdims" "otherfuncs" "clustercodes" "clusterdata" "ttnkit"
do
	sh make-new-prf.sh "$node_name" "$folder"
done
