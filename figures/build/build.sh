#!/bin/bash

VERSION="0.0.1"

# Build the Docker image
docker build --rm -t r_hkgi_flagship_figures:${VERSION} ./figures/build

# Convert the Docker image to a Singularity image file
apptainer build r_hkgi_flagship_figures-v${VERSION}.sif docker-daemon:r_hkgi_flagship_figures:${VERSION}

echo "Finished building r_hkgi_flagship_figures-v${VERSION}.sif"
echo
echo "To run the container, use:"
echo "apptainer run r_hkgi_flagship_figures-v${VERSION}.sif"
echo
echo "To shell into the container, use:"
echo "apptainer shell r_hkgi_flagship_figures-v${VERSION}.sif"