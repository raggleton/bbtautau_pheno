#!/bin/bash

# Make a condor script to run over user-defined masses, number of events, etc

MASS_START=15
MASS_END=60
MASS_STEP=5

NEVT=100000

newfile="generateMasses_${MASS_START}_${MASS_END}_${NEVT}.condor"

echo "Writing new HTCondor script to $newfile"

cp generateMasses.condor $newfile

for (( m=$MASS_START; m <= $MASS_END; m = m + $MASS_STEP ));
do
    echo "arguments = -n ${NEVT} --write --mass ${m} --seed ${m} --name bbtautau_ma1_${m}_${NEVT}.hepmc" >> $newfile
    echo "queue" >> $newfile

done
