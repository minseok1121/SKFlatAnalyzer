#include "IDVariables.h"

IDVariables::IDVariables() {}
IDVariables::~IDVariables() {}

void IDVariables::initializeAnalyzer(){
		MuonIDs = {"NOCUT", "POGLoose", "POGMedium"};
		ElectronIDs = {"NOCUT", "passLooseID", "passMVAID_noIso_WP90"};
		cout << "[IDVariables::initializeAnalyzer] initialization finished" << endl;
}

void IDVariables::executeEvent(){
		// Get All objects
		AllMuons = GetAllMuons();
		AllElectrons = GetAllElectrons();
		GenParts = GetGens();

		// Loop over IDs
		for (unsigned int it_MuonID = 0; it_MuonID < MuonIDs.size(); it_MuonID++) {
				TString MuonID = MuonIDs.at(it_MuonID);
				TString ElectronID = ElectronIDs.at(it_MuonID);
	
				if (!PassMETFilter()) return;


				Event ev = GetEvent();
				Particle METv = ev.GetMETVector();

				vector<Muon>       muons     = SelectMuons(AllMuons, MuonID, 10., 2.4);
				vector<Electron>   electrons = SelectElectrons(AllElectrons, ElectronID, 10., 2.5);
				sort(muons.begin(), muons.end(), PtComparing);
				sort(electrons.begin(), electrons.end(), PtComparing);

				// Match leptons based on loose IDs
				vector<Muon> muons_prompt;
				vector<Electron> electrons_prompt;
				for (const auto &mu: muons)
						if (GetLeptonType(mu, GenParts) == 1) muons_prompt.emplace_back(mu);
				for (const auto &ele: electrons)
						if (GetLeptonType(ele, GenParts) == 1) electrons_prompt.emplace_back(ele);

				// make channels
				TString channel = "";
				if      (muons_prompt.size() == 0 && electrons_prompt.size() == 1) channel = "ejj";
				else if (muons_prompt.size() == 1 && electrons_prompt.size() == 0) channel = "mujj";
				else                                                               channel = "unknown";
		
				// Fill histograms
				// Fill with the same binning in muon / electron efficiency measurements
				TString histkey;
				const vector<double> pt_range = {10., 15., 20., 25., 30., 40., 50., 60., 120., 200., 300., 500., 700., 1200.};
				const vector<double> eta_range = {0., 0.9, 1.2, 2.1, 2.4};
				for (const auto &mu: muons) {
						// Fill inclusive
						const int lepType = GetLeptonType(mu, GenParts);
						histkey = channel+"/muon/"+MuonID+"/Incl/LeptonType-("+TString::Itoa(lepType, 10)+")";
						FillMuon(histkey, mu, 1.);
						// Fill w.r.t. pt & eta
						// Find pt range
						TString this_pt_range = "", this_eta_range = "";
						for (unsigned int i = 0; i < pt_range.size() ; i++) {
								if (i == pt_range.size() - 1) {	// last bin
										this_pt_range = TString::Itoa(static_cast<int>(pt_range.at(i)), 10) + "_to_Inf";
								}
								else if ((pt_range.at(i) < mu.Pt()) && (mu.Pt() < pt_range.at(i+1))) {
										this_pt_range = TString::Itoa(static_cast<int>(pt_range.at(i)), 10) 
																		+ "_to_" + TString::Itoa(static_cast<int>(pt_range.at(i+1)), 10);
										break;
								}
								else {
										continue;
								}
						}
						histkey = channel+"/muon/"+MuonID+"/"+this_pt_range+"/LeptonType-("+TString::Itoa(lepType, 10)+")";
						FillMuon(histkey, mu, 1.);

						for (unsigned int i = 0; i < eta_range.size(); i++) {
								if ((eta_range.at(i) < fabs(mu.Eta())) && (fabs(mu.Eta()) < eta_range.at(i+1))) {
										// segment violation should not occur
										int d1 = static_cast<int>(eta_range.at(i)+1e-6);		// to avoid underflow
										int d2 = static_cast<int>((eta_range.at(i)-d1)*10.+1e-6); 
										TString first_eta= TString::Itoa(d1, 10) + "." + TString::Itoa(d2, 10);
										int d3 = static_cast<int>(eta_range.at(i+1)+1e-6);
										int d4 = static_cast<int>((eta_range.at(i+1)-d3)*10.+1e-6); 
										TString second_eta = TString::Itoa(d3, 10) + "." + TString::Itoa(d4, 10);
										this_eta_range = first_eta + "_to_" + second_eta;
										break;
								}
								else {
										continue;
								}
						}
						histkey = channel+"/muon/"+MuonID+"/"+this_eta_range+"/LeptonType-("+TString::Itoa(lepType, 10)+")";
						FillMuon(histkey, mu, 1.);
				}
				for (const auto &ele: electrons) {
						histkey = channel+"/electron/"+ElectronID+"/LeptonType-("+TString::Itoa(GetLeptonType(ele, GenParts), 10)+")";
						FillElectron(histkey, ele, 1.);
				}
		}
}

