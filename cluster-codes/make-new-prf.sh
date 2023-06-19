#!/bin/bash

# written by ChatGPT 19.06.2023

og_loc=$PWD

cd $HOME/.unison

# Check if the node_name argument is provided
if [ -z "$1" ]; then
  echo "Please provide a node_name argument."
  exit 1
fi

node_name="$1"
folder_name="$2"
filename="iff1500-$folder_name.prf"

# Check if the iff1500-$folder.prf file exists
if [ ! -f "$filename" ]; then
  echo "File $filename does not exist."
  exit 1
fi

# Create a copy of the original file with the modified name
new_filename="${node_name}-$folder_name.prf"
if [ -f "$new_filename" ]; then
  echo "File $new_filename already exists."
  exit 1
fi
cp "$filename" "$new_filename"

# Replace all instances of "iff1500" with "node_name" in the copied file
sed -i "s/iff1500/$node_name/g" "$new_filename"

echo "Replacement completed. Modified file: $new_filename"

cd "$og_loc"
