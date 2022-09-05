#!/bin/sh
flag=$1

SKFlat.py -a Selector -i DYJets -n 20 -e 2016preVFP  --userflag $flag & 
SKFlat.py -a Selector -i DYJets -n 20 -e 2016postVFP  --userflag $flag &
SKFlat.py -a Selector -i DYJets -n 30 -e 2017  --userflag $flag &
SKFlat.py -a Selector -i DYJets -n 40 -e 2018  --userflag $flag &

SKFlat.py -a Selector -i TTLL_powheg -n 20 -e 2016preVFP  --userflag $flag &
SKFlat.py -a Selector -i TTLL_powheg -n 20 -e 2016postVFP  --userflag $flag &
SKFlat.py -a Selector -i TTLL_powheg -n 30 -e 2017  --userflag $flag &
SKFlat.py -a Selector -i TTLL_powheg -n 40 -e 2018  --userflag $flag &

SKFlat.py -a Selector -i DoubleMuon -n 10 -e 2016preVFP  --userflag $flag &
SKFlat.py -a Selector -i DoubleMuon -n 10 -e 2016postVFP  --userflag $flag &
SKFlat.py -a Selector -i DoubleMuon -n 10 -e 2017  --userflag $flag &
SKFlat.py -a Selector -i DoubleMuon -n 10 -e 2018  --userflag $flag &

SKFlat.py -a Selector -l triLepSamples.txt -n 10 -e 2016preVFP --userflag $flag &
SKFlat.py -a Selector -l triLepSamples.txt -n 10 -e 2016postVFP --userflag $flag &
SKFlat.py -a Selector -l triLepSamples.txt -n 10 -e 2017 --userflag $flag &
SKFlat.py -a Selector -l triLepSamples.txt -n 10 -e 2018 --userflag $flag &
SKFlat.py -a Selector -l signalSamples2016a.txt -n 10 -e 2016preVFP --userflag $flag &
SKFlat.py -a Selector -l signalSamples2016b.txt -n 10 -e 2016postVFP --userflag $flag &
SKFlat.py -a Selector -l signalSamples2017.txt -n 10 -e 2017 --userflag $flag &
SKFlat.py -a Selector -l signalSamples2018.txt -n 10 -e 2018 --userflag $flag &
