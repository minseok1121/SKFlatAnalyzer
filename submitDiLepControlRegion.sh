#!/bin/bash
ERA=$1
SKFlat.py -a diLepRegion -i DoubleMuon -n 10 -e ${ERA} &
SKFlat.py -a diLepRegion -i WJets_MG -n 30 -e ${ERA} &
SKFlat.py -a diLepRegion -i TTLJ_powheg -n 30 -e ${ERA} &
SKFlat.py -a diLepRegion -i DYJets -n 30 -e ${ERA} &
SKFlat.py -a diLepRegion -i DYJets10to50_MG -n 20 -e ${ERA} &
SKFlat.py -a diLepRegion -i TTLL_powheg -n 30 -e ${ERA} &
SKFlat.py -a diLepRegion -i SingleTop_sch_Lep -n 10 -e ${ERA} &
SKFlat.py -a diLepRegion -i SingleTop_tch_top_Incl -n 10 -e ${ERA} &
SKFlat.py -a diLepRegion -i SingleTop_tch_antitop_Incl -n 10 -e ${ERA} &
SKFlat.py -a diLepRegion -i SingleTop_tW_top_NoFullyHad -n 10 -e ${ERA} &
SKFlat.py -a diLepRegion -i SingleTop_tW_antitop_NoFullyHad -n 10 -e ${ERA} &
SKFlat.py -a diLepRegion -i WW_pythia -n 10 -e ${ERA} &
SKFlat.py -a diLepRegion -i WZ_pythia -n 10 -e ${ERA} &
SKFlat.py -a diLepRegion -i ZZ_pythia -n 10 -e ${ERA} &
