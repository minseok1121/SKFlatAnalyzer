#!/bin/sh
SKFlat.py -a SampleValidation -l signalSamples2016a.txt -n 10 -e 2016preVFP &
SKFlat.py -a SampleValidation -l signalSamples2016b.txt -n 10 -e 2016postVFP &
SKFlat.py -a SampleValidation -l signalSamples2017.txt -n 10 -e 2017 &
SKFlat.py -a SampleValidation -l signalSamples2018.txt -n 10 -e 2018 &
