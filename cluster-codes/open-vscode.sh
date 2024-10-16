#!/bin/bash

# input gives the node name
NODE_NAME="$1"

# input helps make the remote workspace path
REMOTE_WORKSPACE_PATH="/local/geraghty/${NODE_NAME}-main-git.code-workspace"

# Open VS Code with the remote workspace
code --new-window --remote "ssh-remote+${NODE_NAME}" "${REMOTE_WORKSPACE_PATH}"

# Check if VS Code opened successfully
if [ $? -eq 0 ]; then
    echo "VS Code opened successfully with the remote workspace on ${NODE_NAME}."
else
    echo "Failed to open VS Code with the remote workspace on ${NODE_NAME}."
    exit 1
fi

