import os
import sys
import argparse
#sys.path.remove('/opt/ohpc/pub/apps/root_6_12_06/lib')

from ROOT import TFile
from Plotter.DataAndMC import DataAndMC

# NOTE: This script use pyROOT with python version >= 3.6 because of f-strings
# If you want to use lower python version, change f"...{var}" to "...{}".format(var)
# In the tamsa server, setting an environment with CMSSW_11_* will be sufficient

# set global variables
parser = argparse.ArgumentParser()
parser.add_argument("--object", required=True, type=str, help="histkey")
parser.add_argument("--variable", required=True, type=str, help="histkey")
parser.add_argument("--index", required=True, type=str, help="histkey")
args = parser.parse_args()

if "0" in args.index: index = "/"
if "1" in args.index: index = "/1/"
if "2" in args.index: index = "/2/"
if "3" in args.index: index = "/3/"

if "muon" in args.object: histkey = "TTDiMu/Central/muons"
if "jet" in args.object: histkey = "TTDiMu/Central/jets"
if "bjet" in args.object: histkey = "TTDiMu/Central/bjets"
if "METv" in args.object: histkey = "TTDiMu/Central/METv"
if "pair" in args.object: histkey = "TTDiMu/Central/pair"

if "pt" in args.variable: 
    hist_params = {"x_title": "Pt(#pt)",
               "x_range": [0., 300.],
               "y_title": "Events",
               "ratio_title": "Data/MC",
               "ratio_range": [0.5, 1.5],
               "rebin" : 5,
               }
if "phi" in args.variable: 
    hist_params = {"x_title": "Phi(#phi)",
               "x_range": [-5., 5.],
               "y_title": "Events",
               "ratio_title": "Data/MC",
               "ratio_range": [0.5, 1.5],
               #"rebin" : 5,
               }
if "eta" in args.variable: 
    hist_params = {"x_title": "Eta(#eta)",
               "x_range": [-5., 5.],
               "y_title": "Events",
               "ratio_title": "Data/MC",
               "ratio_range": [0.5, 1.5],
               #"rebin" : 5,
               }
if "size" in args.variable: 
    hist_params = {"x_title": "N(#)",
               "x_range": [0., 10.],
               "y_title": "Events",
               "ratio_title": "Data/MC",
               "ratio_range": [0.5, 1.5],
               #"rebin" : 5,
               }
if "mass" in args.variable: 
    hist_params = {"x_title": "M(#M)",
               "x_range": [0., 300.],
               "y_title": "Events",
               "ratio_title": "Data/MC",
               "ratio_range": [0.5, 1.5],
               #"rebin" : 5,
               }







#histkey = "TTDiMu/Central/pair/mass"
DATASTREAM = "DoubleMuon"
MCs = ["TTLL", "VV", "ST", "DY", "TTLJ", "DYMG"]

# set plotting parameters
cvs_params = {"logy": False,
              "grid": False
             }

info_params = {"info": "L^{int} = 41.5 fb^{-1}",
               "cms_text": "CMS",
               "extra_text": "Work in progress"
              }

# get histograms
def get_hist(sample, histkey):
	if sample == "data":
		fkey = f"DATA/CR_DiLepton_{DATASTREAM}.root"
	else:
		fkey = f"MCSamples/CR_DiLepton_{sample}.root"
	#print(fkey)
	#print(histkey)
	f = TFile.Open(fkey)
	h = f.Get(f"{histkey}{index}{args.variable}")
	h.SetDirectory(0)
	return h

h_data = get_hist("data", histkey)
#h_data = 10*h_data

mc_size = dict()
for mc in MCs:
    print(mc)
    #directory = print('MCSamples/CR_TTbarDiLepton_',mc,'.root')
    print(f"MCSamples/CR_DiLepton_{mc}.root")
    mc_size[mc] = os.path.getsize(f"MCSamples/CR_DiLepton_{mc}.root") 

print(mc_size)
print(sorted(mc_size,key=lambda x:mc_size[x]))
MCs = sorted(mc_size,key=lambda x:mc_size[x])

h_mc = dict()
for mc in MCs:
	h_mc[mc] = get_hist(mc, histkey)

plotter = DataAndMC(cvs_params, hist_params, info_params)
plotter.get_hists(h_data, h_mc)
plotter.combine()
plotter.save(f"DataAndMC_{histkey.replace('/', '_')}{index.replace('/','_')}{args.variable}.png")
