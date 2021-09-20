#!/bin/bash
set -e

filename="../mymac_WScFi.mac"
num_threads=15

num_events=100
energies=(0.1 0.2 0.5 1 2 5 10 20 40 60)

for (( k=0; k<10; k++ ))
do
	sed -i "s/\/gps\/ene\/mono .*/\/gps\/ene\/mono ${energies[k]} GeV/" $filename
	sed -i "s/\/run\/beamOn .*/\/run\/beamOn ${num_events}/" $filename
	make
	echo "Working on energy ${energies[k]}"
	outfile="gev${energies[k]}.log"
	echo "./exampleB4c -m mymac_WScFi.mac -t ${num_threads} > ${outfile}"
	time ./exampleB4c -m mymac_WScFi.mac -t ${num_threads} > ${outfile}
done
