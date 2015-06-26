#!/bin/bash

# Make a condor script to run over user-defined masses, number of events, etc

MASS_START=4
MASS_END=12
MASS_STEP=4

NEVT=10000
NJOBS=10

newfile="generateMasses_${MASS_START}_${MASS_END}_${NEVT}.condor"

echo "Writing new HTCondor script to $newfile"

cp generateMasses.condor $newfile

for (( m=$MASS_START; m <= $MASS_END; m = m + $MASS_STEP ));
do
    echo "arguments = -n ${NEVT} --hepmc --mass ${m} --seed \$(process) --tautau\n" >> $newfile
    echo "queue $NJOBS" >> $newfile

done
