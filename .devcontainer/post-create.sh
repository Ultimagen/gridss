#!/bin/bash
set -e

echo "Setting up GRIDSS development environment..."

# Install development tools not in production image
echo "Installing development tools (Maven, editors)..."
apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -y \
    maven \
    vim \
    nano \
    git \
    && rm -rf /var/lib/apt/lists/*

# Configure git
echo "Configuring git..."
git config --global --add safe.directory /workspace

# Initialize submodules
echo "Initializing git submodules..."
git submodule update --init --recursive

echo "Development environment setup complete!"
