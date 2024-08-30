#!/bin/bash
#cmsRun reco.py inFile=TOTEM21/270000/FABAB63C-8F50-E911-BDCA-0090FAA57730.root outDir=./kjljljlj.root

#cd /afs/cern.ch/user/y/yelberke/scratch1/CMSSW_10_1_7/src/UserCode/ProtonReco/scripts

# 1st argument: input file location
# 2nd argument: output file location
# 3rd argument: path to proxy file (see https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookStartingGrid).
# 4th argument: directory of executable, reco.py.

path="$4/scripts"
cd $path

mkdir -p $(dirname "$2")

# For file permissions
export X509_USER_PROXY=$3

cmsenv

cmsRun reco.py inFile=$1 outDir=$2
