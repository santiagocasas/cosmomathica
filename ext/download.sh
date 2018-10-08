#!/bin/bash

echo "Downloading Eisenstein & Hu's transfer function code..."
mkdir -p tf
wget http://background.uchicago.edu/~whu/transfer/tf_fit.c -q -O tf/tf_fit.c
echo "Downloading Eisenstein & Hu's power spectra code..."
wget http://background.uchicago.edu/~whu/transfer/power.c -q -O tf/power.c
echo "Downloading Robert E. Smith's et al. halofit code..."
mkdir -p halofit
wget http://www.roe.ac.uk/~jap/haloes/halofit+.tar -q -O - | tar x -C halofit
# echo "Downloading CosmicEmulator..."
# mkdir -p CosmicEmulator
# wget http://www.lanl.gov/projects/cosmology/CosmicEmu/CosmicEmu_v1.1.tar.gz -q -O - | tar x -C CosmicEmulator
echo "Please download CAMB yourself at http://camb.info/"

