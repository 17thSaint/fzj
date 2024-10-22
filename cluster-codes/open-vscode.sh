#!/bin/bash

# code written by Claude October 22, 2024

# Check if a node name was provided
if [ -z "$1" ]; then
    echo "Error: Node name not provided"
    echo "Usage: $0 <node_name>"
    exit 1
fi

# Store the node name from input
NODE_NAME="$1"

# Validate node name format (assuming alphanumeric and hyphen only)
if ! [[ $NODE_NAME =~ ^[a-zA-Z0-9-]+$ ]]; then
    echo "Error: Invalid node name format. Use only letters, numbers, and hyphens."
    exit 1
fi

# Check if VS Code is installed
if ! command -v code &> /dev/null; then
    echo "Error: VS Code (code) command not found. Please install VS Code."
    exit 1
fi

# Set workspace path and command based on node name
if [ "$NODE_NAME" = "laptop" ]; then
    WORKSPACE_PATH="$HOME/fzj/main-git/laptop-main-git.code-workspace"
    echo "Opening local workspace..."
    code --new-window "${WORKSPACE_PATH}"
else
    WORKSPACE_PATH="/local/geraghty/${NODE_NAME}-main-git.code-workspace"
    echo "Attempting to connect to ${NODE_NAME}..."
    code --new-window --remote "ssh-remote+${NODE_NAME}" "${WORKSPACE_PATH}"
fi

# Check if VS Code opened successfully
if [ $? -eq 0 ]; then
    echo "✓ VS Code opened successfully"
    echo "  Workspace: ${WORKSPACE_PATH}"
    echo "  Mode: $([ "$NODE_NAME" = "laptop" ] && echo "Local" || echo "Remote on ${NODE_NAME}")"
else
    echo "✗ Failed to open VS Code"
    if [ "$NODE_NAME" = "laptop" ]; then
        echo "  Please check:"
        echo "  - Workspace file exists at ${WORKSPACE_PATH}"
    else
        echo "  Please check:"
        echo "  - SSH connection to ${NODE_NAME} is configured"
        echo "  - Remote workspace file exists at ${WORKSPACE_PATH}"
        echo "  - VS Code Remote-SSH extension is installed"
    fi
    exit 1
fi

