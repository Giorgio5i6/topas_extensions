#!/bin/bash

export GAPSIZE=1 #mm
export BIN=1

for BIN in {80..400..10}
do
	mkdir -p Bin_${BIN}
	envsubst '${GAPSIZE} ${BIN}' < VoxelBasedEkin.txt > Bin_${BIN}/VoxelBasedEkin_Gap${GAPSIZE}_${BIN}.txt
	cd Bin_${BIN}
	topas VoxelBasedEkin_Gap${GAPSIZE}_${BIN}.txt
	cd ../
done
