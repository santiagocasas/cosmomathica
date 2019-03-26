#!/bin/bash


tfdir="tf"
halofitdir="halofit"
frankendir="FrankenEmu"
cambdir="camb"
classdir="class"
mgcambdir="mgcamb"

echo "Downloading of cosmological codes started"


if [ ! -d "$tfdir" ]; then
mkdir -p $tfdir
echo "Downloading Eisenstein & Hu's transfer function code..."
wget http://background.uchicago.edu/~whu/transfer/tf_fit.c -q -O $tfdir/tf_fit.c
echo "Downloading Eisenstein & Hu's power spectra code..."
wget http://background.uchicago.edu/~whu/transfer/power.c -q -O $tfdir/power.c
else
	echo "Eisenstein & Hu's ($tfdir) folder already exists"
fi

if [ ! -d "$halofitdir" ]; then
echo "Downloading Robert E. Smith's et al. halofit code..."
mkdir -p $halofitdir
wget http://www.roe.ac.uk/~jap/haloes/halofit+.tar -q -O - | tar x -C $halofitdir
else
	echo "Robert E. Smith's et al. ($halofitdir) folder already exists"
fi

if [ ! -d "$frankendir" ]; then
echo "Downloading FrankenEmulator..."
mkdir -p $frankendir
wget http://www.hep.anl.gov/cosmology/CosmicEmu/CosmicEmu_v2.tar.gz -q -O - | tar xz -C "$frankendir" --strip-components=1
else
	echo "FrankenEmulator ($frankendir) folder already exists."
fi

# wget http://www.lanl.gov/projects/cosmology/CosmicEmu/CosmicEmu_v1.1.tar.gz -q -O - | tar x -C CosmicEmulator

read -p "Download and clone CAMB git repository (y/n)? " answer
case ${answer:0:1} in
    y|Y )
	echo "removing camb/ folder if existent"
        if [ ! -d "$cambdir" ]; then
		rm -rv $cambdir
	fi
        echo "Cloning CAMB git repo"
	git clone https://github.com/cmbant/CAMB.git
	sleep 2s
	echo "Renaming CAMB to $cambdir"
	mv CAMB/ "$cambdir/"
    ;;
    * )
        echo No CAMB repo downloaded. You can add manually your own. Make sure the wrapper is still valid.
    ;;
esac

read -p "Download and clone CLASS git repository (y/n)? " answer
case ${answer:0:1} in
    y|Y )
	echo "removing class/ folder if existent"
        if [ ! -d "$classdir" ]; then
		rm -rv $classdir
	fi
        echo "Cloning CLASS git repo"
	git clone https://github.com/lesgourg/class_public.git
	sleep 2s
	echo "Renaming class_public to $classdir"
	mv class_public/ "$classdir/"
    ;;
    * )
        echo No CLASS repo downloaded. You can add manually your own. Make sure the wrapper is still valid.
    ;;
esac

read -p "Download and clone MGCAMB git repository (y/n)? " answer
case ${answer:0:1} in
    y|Y )
	echo "removing mgcamb/ folder if existent"
        if [ ! -d "$mgcambdir" ]; then
		rm -rv $mgcambdir
	fi
        echo "Cloning MGCAMB git repo"
	git clone https://github.com/sfu-cosmo/MGCAMB.git
	sleep 2s
	echo "Renaming MGCAMB to $mgcambdir"
	mv MGCAMB/ "$mgcambdir/"
    ;;
    * )
        echo No MGCAMB repo downloaded. You can add manually your own. Make sure the wrapper is still valid.
    ;;
esac

echo "Downloading of cosmological codes finished"
