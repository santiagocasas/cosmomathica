#!/bin/bash

echo "Fixing makefile for halofit"
makefile=halofit/Makefile
echo "Fixing gcc compiler options for halofit"
sed -i~ -e 's/^\(CC.*gcc.*\)/CC = gcc -Wl,--no-as-needed/'      \
    $makefile

