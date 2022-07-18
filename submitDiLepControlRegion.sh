#!/bin/bash
ERA=$1
SKFlat.py -a diLepControlRegion -i DoubleMuon -n 10 -e ${ERA} &
SKFlat.py -a diLepControlRegion -i WJets_MG -n 30 -e ${ERA} &
SKFlat.py -a diLepControlRegion -i TTLJ_powheg -n 30 -e ${ERA} &
SKFlat.py -a diLepControlRegion -i DYJets -n 30 -e ${ERA} &
SKFlat.py -a diLepControlRegion -i TTLL_powheg -n 30 -e ${ERA} &
SKFlat.py -a diLepControlRegion -i SingleTop_sch_Lep -n 10 -e ${ERA} &
SKFlat.py -a diLepControlRegion -i SingleTop_tW_top_NoFullyHad -n 10 -e ${ERA} &
SKFlat.py -a diLepControlRegion -i SingleTop_tW_antitop_NoFullyHad -n 10 -e ${ERA} &
SKFlat.py -a diLepControlRegion -i WW_pythia -n 10 -e ${ERA} &
SKFlat.py -a diLepControlRegion -i WZ_pythia -n 10 -e ${ERA} &
SKFlat.py -a diLepControlRegion -i ZZ_pythia -n 10 -e ${ERA} &
