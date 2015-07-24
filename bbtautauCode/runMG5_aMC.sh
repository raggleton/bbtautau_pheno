#!/bin/bash -e

# This script setups and run MG5_aMC on the worker node on HTCondor

# User must pass in name of process to generate - must correspond to a card in
# MG5_aMC_inputs dir
process="$1"
echo "Doing process $process"
BASE=$PWD


# First we need to setup CMSSW_7_4_4 to get a decent GCC version that actually works.
export SCRAM_ARCH=slc6_amd64_gcc491
VO_CMS_SW_DIR=/cvmfs/cms.cern.ch
. $VO_CMS_SW_DIR/cmsset_default.sh
scramv1 project CMSSW CMSSW_7_4_4_ROOT5
cd CMSSW_*/src/
eval `scramv1 runtime -sh`
cd ../..

# Setup HepMC2
# Use cvmfs one for now...
HEPMC_PATH=/cvmfs/sft.cern.ch/lcg/external/HepMC/2.06.08/x86_64-slc6-gcc48-opt/

# Now need to setup Pythia8 for showering
tar xvzf pythia8*.tgz
cd pythia8*
echo $PWD
./configure --with-hepmc2=$HEPMC_PATH
# ./configure --with-gzip --with-hepmc2=$HEPMC_PATH
make -j8
ls bin
cd ..

# Need to make sure it picks up correct libs and XML file since it uses PYTHIA8DATA var preferentially
export LD_LIBRARY_PATH=$HEPMC_PATH/lib:$LD_LIBRARY_PATH
export PYTHIA8DATA=$BASE/pythia8209/share/Pythia8/xmldoc


# Now setup MG5_aMC
tar xvzf MG5_aMC_v*.tar.gz
ls
cd MG5_aMC_v2_*
echo $PWD
# Hack to add in dl library to shower_card.dat - can't set it in the script file
# as case matters, and if you read in a script file it makes everything lowercase
sed -i "s/Fmcfio/Fmcfio dl/" Template/NLO/Cards/shower_card.dat

# printenv

# Now generate some MC!
# There should be a corresponding input card ${process}.txt,
# and an output dir Background/${process}
card=../MG5_aMC_inputs/"$process".txt
if [ ! -e "$card" ]
then
    echo "Process card file $process.txt does not exist!"
    exit 1
fi

# setup pythia8 location (change to relative?)
sed -i s@/users/ra12451@${BASE}@ $card

./bin/mg5_aMC "$card"

echo "mcatnlo.log"
cat pp_bbx/mcatnlo.log

if [ -e pp_bbx/MCatNLO/RUN_PYTHIA8_1/shower_card.dat ]; then
    echo "pp_bbx/MCatNLO/RUN_PYTHIA8_1/shower_card.dat"
    cat pp_bbx/MCatNLO/RUN_PYTHIA8_1/shower_card.dat
fi

if [ -e pp_bx/Cards/shower_card.dat ]; then
    echo "pp_bx/Cards/shower_card.dat"
    cat pp_bx/Cards/shower_card.dat
fi

if [ -e pp_bx/Cards/shower_card_default.dat ]; then
    echo "pp_bx/Cards/shower_card_default.dat"
    cat pp_bx/Cards/shower_card_default.dat
fi

ls pp_bbx/MCatNLO/lib/

# Move to HDFS
outdir=/hdfs/user/ra12451/bbtautau/Background/"$process"/
if [ ! -d "$outdir" ]
then
    echo "Making output dir $outdir"
    mkdir "$outdir"
fi
echo "Moving output to $outdir"
# N=`ls "$outdir" | wc -l`
# N=$(( $N + 1 ))

# make a unique run folder
dt=$(date '+%d_%b_%y_%H%M')
mv "$process"/Events/run_* "$outdir"/run_$dt