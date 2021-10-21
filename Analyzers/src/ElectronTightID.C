#include "ElectronTightID.h"

void ElectronTightID::initializeAnalyzer(){

}

void ElectronTightID::executeEvent(){
	
	if (!PassMETFilter()) return;
	Event ev = GetEvent();
	vector<Gen> gens = GetGens();
	vector<Muon> muons = GetAllMuons();
	vector<Electron> electrons = GetAllElectrons();
	
	sort(muons.begin(), muons.end(), PtComparing);
	sort(electrons.begin(), electrons.end(), PtComparing);

	// select leptons
	vector<Muon> muons_tight = SelectMuons(muons, "HcToWATight", 10., 2.4);
	vector<Muon> muons_loose = SelectMuons(muons, "HcToWALoose", 10., 2.4);
	vector<Electron> electrons_tight = SelectElectrons(electrons, "HcToWATight", 10., 2.5);
	vector<Electron> electrons_loose = SelectElectrons(electrons, "HcToWALoose", 10., 2.5);

	if (muons_tight.size() != 0) return;
	if (muons_loose.size() != 0) return;
	if (electrons_tight.size() != 2) return;
	if (electrons_loose.size() != 2) return;

	if (electrons_tight.at(0).Pt() < 30.) return;

	double weight = 1.;
	weight *= ev.MCweight() * weight_norm_1invpb;
	weight *= ev.GetTriggerLumi("Full");

	// TODO: check electron efficiency for
	// TODO: prompt & nonprompt leptons
	// TODO: with WP90 and WP80
	TString key = "";
	for (unsigned int i = 0; i < electrons_tight.size(); i++) {
		key = "pte" + TString::Itoa(i+1, 10);
		const Electron& ele = electrons_tight.at(i);
		const int lepType = GetLeptonType(ele, gens);
		if (lepType < 0)
			key += "_nonprompt";
		else if (lepType == 1)
			key += "_prompt";
		else if (lepType == 3)
			key += "_ewDaughter";
		else if (lepType == 4 || lepType == 5)
			key += "_conversion";
		else
			continue;

		FillHist(key, ele.Pt(), weight, 300, 0., 300.);
		if (ele.passMVAID_noIso_WP80()) {
			key += "_passWP80";
			FillHist(key, ele.Pt(), weight, 300, 0., 300.);
		}
	}
}

ElectronTightID::ElectronTightID(){

}

ElectronTightID::~ElectronTightID(){

}


