#!/bin/sh

python3 plot_DiMu.py --object muon --index 1 --variable pt
python3 plot_DiMu.py --object muon --index 2 --variable pt
python3 plot_DiMu.py --object jets --index 1 --variable pt
python3 plot_DiMu.py --object jets --index 2 --variable pt
python3 plot_DiMu.py --object jets --index 3 --variable pt
python3 plot_DiMu.py --object bjets --index 1 --variable pt
python3 plot_DiMu.py --object bjets --index 2 --variable pt
python3 plot_DiMu.py --object bjets --index 3 --variable pt
python3 plot_DiMu.py --object METv --index 0 --variable pt
python3 plot_DiMu.py --object pair --index 0 --variable pt

python3 plot_DiMu.py --object muon --index 1 --variable eta
python3 plot_DiMu.py --object muon --index 2 --variable eta
python3 plot_DiMu.py --object jets --index 1 --variable eta
python3 plot_DiMu.py --object jets --index 2 --variable eta
python3 plot_DiMu.py --object jets --index 3 --variable eta
python3 plot_DiMu.py --object bjets --index 1 --variable eta
python3 plot_DiMu.py --object bjets --index 2 --variable eta
python3 plot_DiMu.py --object bjets --index 3 --variable eta
python3 plot_DiMu.py --object pair --index 0 --variable eta

python3 plot_DiMu.py --object muon --index 1 --variable phi
python3 plot_DiMu.py --object muon --index 2 --variable phi
python3 plot_DiMu.py --object jets --index 1 --variable phi
python3 plot_DiMu.py --object jets --index 2 --variable phi
python3 plot_DiMu.py --object jets --index 3 --variable phi
python3 plot_DiMu.py --object bjets --index 1 --variable phi
python3 plot_DiMu.py --object bjets --index 2 --variable phi
python3 plot_DiMu.py --object bjets --index 3 --variable phi
python3 plot_DiMu.py --object METv --index 0 --variable phi
python3 plot_DiMu.py --object pair --index 0 --variable phi

python3 plot_DiMu.py --object pair --index 0 --variable mass

python3 plot_DiMu.py --object bjets --index 0 --variable size
python3 plot_DiMu.py --object jets --index 0 --variable size