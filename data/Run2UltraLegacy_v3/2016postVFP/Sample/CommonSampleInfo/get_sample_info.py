import os
from ROOT import TFile

sample_dir = "/gv0/Users/choij/SKFlat/2016postVFP"
prefix = "TTToHcToWAToMuMu"
suffix = "MultiLepFilter_TuneCP5_13TeV-madgraph-pythia8"
#masspoints = [
#        "MHc-70_MA-15",
#        "MHc-70_MA-40",
#        "MHc-70_MA-65",
#        "MHc-100_MA-15",
#        "MHc-100_MA-60",
#        "MHc-100_MA-95",
#        "MHc-130_MA-15",
#        "MHc-130_MA-55",
#        "MHc-130_MA-90",
#        "MHc-130_MA-125",
#        "MHc-160_MA-15",
#        "MHc-160_MA-80",
#        "MHc-160_MA-85",
#        "MHc-160_MA-155"]
masspoints = ["MHc-160_MA-120"]

for mp in masspoints:
		try:
				# get crab outputs
				parent_key = f"{sample_dir}/{prefix}_{mp}_{suffix}/SKFlat_Run2UltraLegacy_v3"
				assert os.path.exists(parent_key)
				assert len(os.listdir(parent_key)) == 1
				child_key = f"{parent_key}/{os.listdir(parent_key)[0]}/0000"
				assert os.path.exists(child_key)
				crab_outputs = list(filter(lambda x: "root" in x, os.listdir(child_key)))

				# make information
				alias = f"{prefix}_{mp}"
				PD = f"{alias}_MultiLepFilter_TuneCP5_13TeV-madgraph-pythia8"
				xsec = 0.015
				nmc = 0
				sumw = 0
				for output in crab_outputs:
						f = TFile.Open(f"{child_key}/{output}")
						tree = f.Get("recoTree/SKFlat")
						nmc += tree.GetEntries()
						for evt in tree:
								sumw += evt.gen_weight
						f.Close()

				# write outfile
				with open(f"{alias}.txt", "w") as f:
						f.write("# alias PD xsec nmc sumsign sumw\n")
						f.write(f"{alias}\t{PD}\t{xsec}\t{int(nmc)}\t{float(nmc)}\t{sumw:.2f}\n")
		except Exception as e:
				print(e)
				continue

