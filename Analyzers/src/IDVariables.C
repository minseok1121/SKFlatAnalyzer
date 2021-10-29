#include "IDVariables.h"

void IDVariables::initializeAnalyzer(){
		MuonID = "POGMedium";
		ElectronID = "passMVAID_noIso_WP90";

		cout << "[IDVariables::initializeAnalyzer] initialization finished" << endl;
}

void IDVariables::executeEvent(){
		// Get All objects
		AllMuons = GetAllMuons();
		AllElectrons = GetAllElectrons();
		GenParts = GetGens();

		if (!PassMETFilter()) return;

		Event ev = GetEvent();
		Particle METv = ev.GetMETVector();

		vector<Muon> muons = SelectMuons(AllMuons, MuonID, 10., 2.4);
		vector<Electron> electrons = SelectElectrons(AllElectrons, ElectronID, 10., 2.5);

		sort(muons.begin(), muons.end(), PtComparing);
		sort(electrons.begin(), electrons.end(), PtComparing);
		
		//==== Match leptons
		vector<Muon> muons_prompt, muons_nonprompt;
		vector<Electron> electrons_prompt, electrons_nonprompt;
		
		for (const auto &mu : muons) {
				const int lepType = GetLeptonType(mu, GenParts);
				if (lepType == 1)
						muons_prompt.emplace_back(mu);
				else if (lepType == -2 || lepType == -3 || lepType == -4)
						muons_nonprompt.emplace_back(mu);
				else
						continue;
		}
		for (const auto &ele : electrons) {
				const int lepType = GetLeptonType(ele, GenParts);
				if (lepType == 1)
						electrons_prompt.emplace_back(ele);
				else if (lepType == -2 || lepType == -3 || lepType == -4)
						electrons_nonprompt.emplace_back(ele);
				else
						continue;
		}

		// Fill histograms
		FillMuons("muons_prompt", muons_prompt, 1., true);
		FillMuons("muons_nonprompt", muons_nonprompt, 1., true);
		FillElectrons("electrons_prompt", electrons_prompt, 1., true);
		FillElectrons("electrons_nonprompt", electrons_nonprompt, 1., true);

}

IDVariables::IDVariables(){

}

IDVariables::~IDVariables(){

}


