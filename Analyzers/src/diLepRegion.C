#include "diLepRegion.h"

diLepRegion::diLepRegion() {}
diLepRegion::~diLepRegion() {
		delete h_muon_sf;
		delete h_ele_sf;
}

void diLepRegion::initializeAnalyzer(){
		// Userflags
		RunSysts = HasFlag("RunSysts");
		RunLowPt = HasFlag("RunLowPt");

		// Set systematic
		if (RunSysts) Systs = {"Central", "NoIDSF", "MuonIDSFUp", "MuonIDSFDown", "EleIdSFUp", "EleIDSFDown"};
		else Systs = {"Central"};

		// scale factors to use
		TFile* f_muon = new TFile("$SKFlat_WD/data/Run2UltraLegacy_v2/"+DataEra+"/ID/Muon/eff_muon.root");
		h_muon_sf = (TH2D*)f_muon->Get("SF_eta_pt");
		h_muon_sf->SetDirectory(0);
		f_muon->Close();

		TFile* f_electron = new TFile("$SKFlat_WD/data/Run2UltraLegacy_v2/"+DataEra+"/ID/Electron/eff_electron.root");
		h_ele_sf = (TH2D*)f_electron->Get("sf");
		h_ele_sf->SetDirectory(0);
		f_electron->Close();

		// trigger && ID settings
		if (DataEra == "2016preVFP") {
				DblMuTrigs = {
						"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v",
						"HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v",
				};
				EMuTrigs = {
						"HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v",
						"HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",
						"HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v",
						"HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v"
				};
				MuonIDs = {"HcToWATight", "HcToWALoose"};
				ElectronIDs = {"HcToWATight16", "HcToWALoose16"};
		}
		else if (DataEra == "2016postVFP") {
				DblMuTrigs = {
						"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v",
						"HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v",
				};
				EMuTrigs = {
						"HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v",
						"HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",
						"HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v",
						"HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v"
				};
				MuonIDs = {"HcToWATight", "HcToWALoose"};
				ElectronIDs = {"HcToWATight16", "HcToWALoose16"};
		}
		else if (DataEra == "2017") {
				DblMuTrigs = {
						"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v",
						"HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v",
						"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v"
				};
				EMuTrigs = {
						"HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v",
						"HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v"
				};
				MuonIDs = {"HcToWATight", "HcToWALoose"};
				ElectronIDs = {"HcToWATight17", "HcToWALoose17"};
		}
		else if (DataEra == "2018") {
				DblMuTrigs = {
						"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v",
						"HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v",
						"HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_Mass3p8_v"
				};
				EMuTrigs = {
						"HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v",
						"HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",
						"HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v"
				};
				MuonIDs = {"HcToWATight", "HcToWALoose"};
				ElectronIDs = {"HcToWATight18", "HcToWALoose18"};
		}

		// B-tagging
		// Copy and paste the code from PreLegacy, might buggy in UL
		vector<JetTagging::Parameters> jtps;
		jtps.emplace_back(
						JetTagging::Parameters(JetTagging::DeepJet, JetTagging::Medium, JetTagging::incl, JetTagging::mujets));
		mcCorr->SetJetTaggingParameters(jtps);

}


void diLepRegion::executeEvent(){

		if (!PassMETFilter()) return;
		Event ev = GetEvent();
		vector<Gen>				gens = GetGens();
		vector<Muon>			muons_all = GetAllMuons();
		vector<Electron>	electrons_all = GetAllElectrons();
		vector<Jet>				jets_all = GetAllJets();
		Particle					METv = ev.GetMETVector();

		// sort objects
		sort(muons_all.begin(), muons_all.end(), PtComparing);
		sort(electrons_all.begin(), electrons_all.end(), PtComparing);
		sort(jets_all.begin(), jets_all.end(), PtComparing);

		// select objects
		vector<Muon> muons_tight = SelectMuons(muons_all, MuonIDs.at(0), 10., 2.4);
		vector<Muon> muons_loose = SelectMuons(muons_all, MuonIDs.at(1), 10., 2.4);
		vector<Electron> electrons_tight = SelectElectrons(electrons_all, ElectronIDs.at(0), 10., 2.5);
		vector<Electron> electrons_loose = SelectElectrons(electrons_all, ElectronIDs.at(1), 10., 2.5);
		vector<Jet> jets_tight = SelectJets(jets_all, "tight", 20., 2.4);
		vector<Jet> jets = JetsVetoLeptonInside(jets_tight, electrons_loose, muons_loose, 0.4);
		vector<Jet> bjets;
		for (const auto &jet: jets) {
				const double score = jet.GetTaggerResult(JetTagging::DeepJet);
				if (score > mcCorr->GetJetTaggingCutValue(JetTagging::DeepJet, JetTagging::Medium))
						bjets.emplace_back(jet);
		}

		// event Selection
		TString channel = selectEvent(ev, muons_tight, electrons_tight,
																				muons_loose, electrons_loose,
																				jets, bjets, METv);
		if (channel == "") return;
		
		// set event weight and fill histograms
		if (IsDATA) {
				double weight = 1.;
				FillMuons(channel+"/Central/muons_tight", muons_tight, weight, true);
				FillElectrons(channel+"/Central/electrons_tight", electrons_tight, weight, true);
				FillJets(channel+"/Central/jets", jets, weight);
				FillJets(channel+"/Central/bjets", bjets, weight);
				FillObject(channel+"/Central/METv", METv, weight);
				if (channel.Contains("DblMu")) {
						Particle diMuon = muons_tight.at(0) + muons_tight.at(1);
						FillObject(channel+"/Central/diMuon", diMuon, weight);
				}
		}
		else {
				for (const auto &syst: Systs) {
						double weight = getWeight(ev, muons_tight, electrons_tight, jets, syst);
						FillMuons(channel+"/"+syst+"/muons_tight", muons_tight, weight, true);
						FillElectrons(channel+"/"+syst+"/electrons_tight", electrons_tight, weight, true);
						FillJets(channel+"/"+syst+"/jets", jets, weight);
						FillJets(channel+"/"+syst+"/bjets", bjets, weight);
						FillObject(channel+"/"+syst+"/METv", METv, weight);
						if (channel.Contains("DblMu")) {
								Particle diMuon = muons_tight.at(0) + muons_tight.at(1);
								FillObject(channel+"/"+syst+"/diMuon", diMuon, weight);
						}
				}
		}
}

TString diLepRegion::selectEvent(
				Event &ev, vector<Muon> &muons_tight, vector<Electron> &electrons_tight,
				vector<Muon> &muons_loose, vector<Electron> &electrons_loose,
				vector<Jet> &jets, vector<Jet> &bjets, Particle &METv) {

		// MuMu
		if (muons_tight.size() == 2 && muons_loose.size() == 2 && electrons_loose.size() == 0) {
				if (! ev.PassTrigger(DblMuTrigs)) return "";
				if (! (muons_tight.at(0).Charge() + muons_tight.at(1).Charge() == 0)) return "";
				if (! (muons_tight.at(0).Pt() > 20.)) return "";
				if (! (muons_tight.at(1).Pt() > 10.)) return "";

				const Particle ZCand = muons_tight.at(0) + muons_tight.at(1);
				if (fabs(ZCand.M() - 91.2) < 15.) {
						// DY candidate, jet conditions
						if (bjets.size() != 0) return "";
						
						return "DY_DblMu";
				}
				else {
						// avoid contributions from low mass resonance
						if (! (ZCand.M() > 12.)) return "";
						// TT candidate, jet  conditions
						if (! (jets.size() >= 2)) return "";
						if (! (bjets.size() >= 1)) return "";

						// MET condition
						if (! (METv.Pt() > 40.)) return "";

						return "TT_DblMu";
				}
		}
		// EMu
		else if (muons_tight.size() == 1 && muons_loose.size() == 1 && electrons_tight.size() == 1 && electrons_loose.size() == 1) {
				if (! ev.PassTrigger(EMuTrigs)) return "";
				const Muon &mu = muons_tight.at(0);
				const Electron &ele = electrons_tight.at(0);
				if (! (mu.Charge() + ele.Charge() == 0)) return "";
				if (! ((mu.Pt() > 10. && ele.Pt() > 25.) || (mu.Pt() > 25. && ele.Pt() > 15.))) return "";
				if (! (mu.DeltaR(ele) > 0.4)) return "";
				if (! (jets.size() >= 2)) return "";
				if (! (bjets.size() >= 1)) return "";
				return "TT_EMu";
		}
		else 
				return "";
}

double diLepRegion::getWeight(Event &ev, vector<Muon> &muons, vector<Electron>& electrons,
														  vector<Jet> &jets, TString syst) {

		if (IsDATA) return 1.;

		double weight = 1.;
		const double w_prefire = GetPrefireWeight(0);
		const double w_gen = ev.MCweight()*weight_norm_1invpb;
		const double w_lumi = ev.GetTriggerLumi("Full");
		const double w_pileup = GetPileUpWeight(nPileUp, 0);
		//cout << "w_prefire: " << w_prefire << endl;
		//cout << "w_gen: " << w_gen << endl;
		//cout << "w_lumi: " << w_lumi << endl;
		//cout << "w_pileup: " << w_pileup << endl;
		weight *= (w_prefire*w_gen*w_lumi*w_pileup);
		
		// Set ID scale factor
		double w_idsf = 1.;
		for (const auto &mu: muons) {
				if (syst == "NoIDSF") { w_idsf *= 1.; }
				else if (syst == "MuonIDSFUp") { w_idsf *= getMuonIDSF(mu, 1); }
				else if (syst == "MuonIDSFDown") { w_idsf *= getMuonIDSF(mu, -1); }
				else { w_idsf *= getMuonIDSF(mu, 0); } // Central or something something...
		}
		for (const auto &ele: electrons) {
				if (syst == "NoIDSF") { w_idsf *= 1.; }
				else if (syst == "EleIdSFUp") { w_idsf *= getElectronIDSF(ele, 1); }
				else if (syst == "EleIDSFDown") { w_idsf *= getElectronIDSF(ele, -1); }
				else { w_idsf *= getElectronIDSF(ele, 0); }
		}
		//cout << "[diLepRegion::getWeight - " << syst << "] w_idsf:" << w_idsf << endl;
		weight *= w_idsf;

		return weight;
}

double diLepRegion::getMuonIDSF(const Muon &mu, int sys) {
		if (IsDATA) return 1.;

		double value = 1., error = 0.;
		// Use 10 < pt < 15 GeV SFs directly
		if (!RunLowPt) {
				const double pt = min(mu.Pt(), 200.);
				const double eta = fabs(mu.Eta());
				value = h_muon_sf->GetBinContent(h_muon_sf->FindBin(eta, pt));
				error = h_muon_sf->GetBinError(h_muon_sf->FindBin(eta, pt));
		}
		else {
				const double pt = max(15., min(mu.Pt(), 200.));
				const double eta = fabs(mu.Eta());
				value = h_muon_sf->GetBinContent(h_muon_sf->FindBin(eta, pt));
				error = h_muon_sf->GetBinError(h_muon_sf->FindBin(eta, pt));
				if (mu.Pt() < 15.) error *= 2;
		}
		return value+sys*error;
}

double diLepRegion::getElectronIDSF(const Electron &ele, int sys) {
		if (IsDATA) return 1.;

		const double pt = min(ele.Pt(), 500.);
		const double sceta = fabs(ele.scEta());
		const double value = h_ele_sf->GetBinContent(h_ele_sf->FindBin(sceta, pt));
		const double error = h_ele_sf->GetBinError(h_ele_sf->FindBin(sceta, pt));

		return value+sys*error;
}
