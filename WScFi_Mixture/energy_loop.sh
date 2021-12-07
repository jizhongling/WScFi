#!/bin/bash
set -o errexit

filename="mymac_WScFi.mac"
num_threads=1

particle="e-"

num_events=10000
energies=(0.1 0.2 0.5 1 2 5 10 20 40 60)
k=$1

mkdir -p ${particle}_FTFP
mkdir -p proc$1
cd proc$1

cp ../$filename .
sed -i "s/\/analysis\/setFileName .*/\/analysis\/setFileName out/" $filename
sed -i "s/\/gps\/particle .*/\/gps\/particle ${particle}/" $filename
sed -i "s/\/gps\/ene\/mono .*/\/gps\/ene\/mono ${energies[k]} GeV/" $filename
sed -i "s/\/run\/beamOn .*/\/run\/beamOn ${num_events}/" $filename
../build/exampleB4c -m mymac_WScFi.mac -t ${num_threads}
mv out.root ../${particle}_FTFP/${particle}_${energies[k]}GeV.root

cd ..
rm -rf proc$1
