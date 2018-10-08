#!/bin/bash

# Switch the default compiler (ifort) in the CAMB Makefile to gfortran

makefile=halofit/Makefile
echo "Fixing gcc compiler options for halofit"
sed -i~ -e 's/^\(CC.*gcc.*\)/CC = gcc -Wl,--no-as-needed/'      \
    $makefile
