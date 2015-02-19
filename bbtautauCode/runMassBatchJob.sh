#!/bin/bash
#
#PBS -l walltime=11:59:00
#
cd ${HOME}/bbtautau_pheno/bbtautauCode
${HOME}/bbtautau_pheno/bbtautauCode/bbtautau.exe ${args}
