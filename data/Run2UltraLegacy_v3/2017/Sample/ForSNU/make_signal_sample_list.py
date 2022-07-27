import os

sample_dir = "/gv0/Users/choij/SKFlat/2017"
prefix = "TTToHcToWAToMuMu"
suffix = "MultiLepFilter_TuneCP5_13TeV-madgraph-pythia8"
masspoints = [
				"MHc-70_MA-15",
				"MHc-70_MA-40",
				"MHc-70_MA-65",
				"MHc-100_MA-15",
				"MHc-100_MA-60",
				"MHc-100_MA-95",
				"MHc-130_MA-15",
				"MHc-130_MA-55",
				"MHc-130_MA-90",
				"MHc-130_MA-125",
				"MHc-160_MA-15",
				"MHc-160_MA-85",
				"MHc-160_MA-120",
				"MHc-160_MA-155"]

for mp in masspoints:
		try:
				parent_key = f"{sample_dir}/{prefix}_{mp}_{suffix}/SKFlat_Run2UltraLegacy_v3"
				assert os.path.exists(parent_key)
				assert len(os.listdir(parent_key)) == 1 
				child_key = f"{sample_dir}/{prefix}_{mp}_{suffix}/SKFlat_Run2UltraLegacy_v3/{os.listdir(parent_key)[0]}/0000"
				assert os.path.exists(child_key)
				crab_outputs = os.listdir(child_key)
				crab_outputs = list(filter(lambda x: "root" in x, crab_outputs))

				# now write the file
				with open(f"TTToHcToWAToMuMu_{mp}.txt", "w") as f:
						for output in crab_outputs:
								f.write(f"{child_key}/{output}\n")
		except Exception as e:
				print(e)
				continue
		
