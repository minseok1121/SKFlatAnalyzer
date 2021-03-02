#include "RegionSelector.h"

vector<TString> RegionSelector::getCuts(TString region) {
	if (region == "SR_3mu")
		return {"noCut", "METFilter", "3mu", "trigger", "safePtCut", "ExistOSDiMu", "MOSDiMu_ge12", "Nj_ge2", "Nb_ge1"};
	else if (region == "SR_1e2mu")
		return {"noCut", "METFilter", "1e2mu", "trigger", "safePtCut", "OSDiMu", "MOSDiMu_ge12","Nj_ge2", "Nb_ge1"};
	else if (region == "WZ_3mu")
		return {"noCut", "METFilter", "3mu", "trigger", "safePtCut", "ExistOSDiMu", "MOsDiMu_ge12","Nj_le1", "NoBjet"};
	else if (region == "WZ_1e2mu")
		return {"noCut", "METFilter", "1e2mu", "trigger", "safePtCut", "OSDiMu", "MOsDiMu_ge12", "Nj_le1", "NoBjet"};
	else if (region == "DY_OSdimu")
		return {"noCut", "METFilter", "2mu", "trigger", "safePtCut", "OSDiMu", "OnshellZ", "NoBjet"};
	else if (region == "TT_OSdimu")
		return {"noCut", "METFilter", "2mu", "trigger", "safePtCut", "OSDiMu", "OffshellZ", "MOsDiMu_ge12", "MET_ge40", "dRll_ge04", "Nj_ge2", "Nb_ge2"};
	else if (region == "TT_OSemu")
		return {"noCut", "METFilter", "emu", "trigger", "safePtCut", "OSemu", "dRll_ge04", "Nj_ge2", "Nb_ge1"};
	else {
		cerr << "[RegionSelector::getCuts] wrong region " << region << endl;
		exit(EXIT_FAILURE);
	}
}

// Define each cuts
TString RegionSelector::Selector(
		Event& ev,
		vector<Muon>& muons_tight, vector<Electron>& electrons_tight,
		vector<Muon>& muons_loose, vector<Electron>& electrons_loose,
		vector<Jet>& jets, vector<Jet>& bjets, Particle& METv) const {
	return "";
}
