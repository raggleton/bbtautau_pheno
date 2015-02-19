#!/bin/bash

nevt=10000

for i in {10..10}
do
    echo "qsub -N bbtt$i -v args=\"-n ${nevt} --write --mass ${i} --seed ${i} --name /storage/ra12451/bbtautau/bbtautau_ma1_${i}_${nevt}.hepmc\" runMassBatchJobPBS.sh"
    qsub -N bbtt$i -v args="-n ${nevt} --write --mass ${i} --seed ${i} --name /storage/ra12451/bbtautau/bbtautau_ma1_${i}_${nevt}.hepmc" runMassBatchJobPBS.sh
    sleep 10
done
