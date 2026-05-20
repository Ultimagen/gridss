#!/bin/bash
# Download required dependencies for Docker build

set -e

echo "Downloading RepeatMasker..."
gsutil cp gs://concordance-data/scripts/RepeatMasker-4.1.2-p1.tar.gz .

echo "Downloading rmblast..."
gsutil cp gs://concordance-data/scripts/rmblast-2.14.1+-x64-linux.tar.gz .

echo "Dependencies downloaded successfully!"
echo "You can now run: docker build ..."
