#/bin/sh

SKFlat.py -a TutorialRun -i SingleMuon -n 10 -e 2017 --reduction 10 &
SKFlat.py -a TutorialRun -l MCList_tutorial.txt -n 10 -e 2017 --reduction 10 &
