#!/bin/bash


echo "Fixing makefile for camb"
makefilecamb=camb/Makefile_main
echo "Changing EQUATIONS variable to equations_ppf"
sed -i~ -e 's/^\(EQUATIONS.*equations.*\)/EQUATIONS   ?=   equations_ppf/'      \
    $makefilecamb
