#!/bin/bash
echo "Starting job on host: $(hostname)"
echo "Current dir: $(pwd)"
echo "Listing initialdir content:"
ls -l

echo "Argument received: $1"
echo "Listing the file list argument:"
ls -l $1

source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARC=el8_amd64_gcc12
eval `scramv1 project CMSSW CMSSW_14_0_15`
cd CMSSW_14_0_15/src/
eval `scramv1 runtime -sh`
ls
pwd
cd ../..
ls
pwd
python3 analyzer_zmc.py $1
xrdcp test/histos.root /home/yide/fitting/output_histo/histos_$(basename $1 .txt).root
