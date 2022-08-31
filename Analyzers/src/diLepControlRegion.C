#include "diLepControlRegion.h"

diLepControlRegion::diLepControlRegion() {}
diLepControlRegion::~diLepControlRegion() {}

void diLepControlRegion::initializeAnalyzer(){
		// trigger & ID settings
		if (DataEra == "2016preVFP") {
				DblMuTriggers = {
            "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v",
            "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v",
        };
        EMuTriggers = {
            "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v",
            "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",
            "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v",
            "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v"
        };
				MuonIDs = {"HcToWATight", "HcToWALoose", "HcToWAVeto"};
				ElectronIDs = {"HcToWATight16a", "HcToWALoose16a", "HcToWAVeto16a"};
		}
		else if (DataEra == "2016postVFP") {
				DblMuTriggers = {
            "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v",
            "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v",
        };
        EMuTriggers = {
            "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v",
            "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",
            "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v",
            "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v"
        };
				MuonIDs = {"HcToWATight", "HcToWALoose", "HcToWAVeto"};
				ElectronIDs = {"HcToWATight16b", "HcToWALoose16b", "HcToWAVeto16b"};
		}
		else if (DataEra == "2017") {
				DblMuTriggers = {
            "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v",
            "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v",
            "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v"
        };
        EMuTriggers = {
            "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v",
            "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v"
        };
				MuonIDs = {"HcToWATight", "HcToWALoose", "HcToWAVeto"};
				ElectronIDs = {"HcToWATight17", "HcToWALoose17", "HcToWAVeto17"};
		}
		else if (DataEra == "2018") {
				DblMuTriggers = {
            "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v",
            "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v",
            "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v"
        };
        EMuTriggers = {
            "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v",
            "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",
            "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v"
        };
				MuonIDs = {"HcToWATight", "HcToWALoose", "HcToWAVeto"};
				ElectronIDs = {"HcToWATight18", "HcToWALoose18", "HcToWAVeto18"};
		}
		else {
				cerr << "[diLepControlRegion::initializeAnalyzer] Wrong era " << DataEra << endl;
				exit(EXIT_FAILURE);
		}

		// B-tagging
		vector<JetTagging::Parameters> jtps;
		jtps.emplace_back(
				JetTagging::Parameters(JetTagging::DeepJet, JetTagging::Medium, JetTagging::incl, JetTagging::mujets));
		mcCorr->SetJetTaggingParameters(jtps);

		// Systematics
		Systematics = {// to see how the correction changes
									 "Default",
									 "MuonIDSFOnly", "MuonIDSFOnlyUp", "MuonIDSFOnlyDown",
									 "Central",
									 "PileUpCorrUp", "PileUpCorrDown",
									 "MuonEnUp", "MuonEnDown",
									 "MuonIDSFUp", "MuonIDSFDown", 
									 "MuonTrigSFUp", "MuonTrigSFDown",
									 //"EleIDSFUp", "EleIDSFDown",
									 //"EleTrigSFUp", "EleTrigSFDown"
		              };

}

void diLepControlRegion::executeEvent(){
		
		if (! PassMETFilter()) return;

		Event ev = GetEvent();
		vector<Gen>	truth_coll = GetGens();
		vector<Muon> AllMuons = GetAllMuons();
		vector<Electron> AllElectrons = GetAllElectrons();
		vector<Jet>	AllJets = GetAllJets();
		Particle METv = ev.GetMETVector();

		for (const auto &syst: Systematics) {
				vector<Muon> muon_coll = AllMuons;
				vector<Electron> ele_coll = AllElectrons;
				vector<Jet> jet_coll = AllJets;
				if (syst == "MuonEnUp") { muon_coll = ScaleMuons(muon_coll, +1); }
				if (syst == "MuonEnDown") {muon_coll = ScaleMuons(muon_coll, -1); }

				// sort objects
				sort(muon_coll.begin(), muon_coll.end(), PtComparing);
				sort(ele_coll.begin(), ele_coll.end(), PtComparing);
				sort(jet_coll.begin(), jet_coll.end(), PtComparing);

				// select objects
				vector<Muon> muonT_coll = SelectMuons(muon_coll, MuonIDs.at(0), 10., 2.4);
				vector<Muon> muonL_coll = SelectMuons(muon_coll, MuonIDs.at(1), 10., 2.4);
				vector<Muon> muonV_coll = SelectMuons(muon_coll, MuonIDs.at(2), 10., 2.4);
				vector<Electron> eleT_coll = SelectElectrons(ele_coll, ElectronIDs.at(0), 10., 2.5);
				vector<Electron> eleL_coll = SelectElectrons(ele_coll, ElectronIDs.at(1), 10., 2.5);
				vector<Electron> eleV_coll = SelectElectrons(ele_coll, ElectronIDs.at(2), 10., 2.5);
				vector<Jet> jetT_coll = SelectJets(jet_coll, "tight", 20., 2.4);
				jetT_coll = JetsVetoLeptonInside(jetT_coll, eleV_coll, muonV_coll, 0.4);
				vector<Jet> bjet_coll;
				for (const auto &j: jetT_coll) {
						const double score = j.GetTaggerResult(JetTagging::DeepJet);
						if (score > mcCorr->GetJetTaggingCutValue(JetTagging::DeepJet, JetTagging::Medium))
								bjet_coll.emplace_back(j);
				}

				// event selection
				channel_t channel = selectEvent(ev, muonT_coll, muonL_coll, muonV_coll,
																				   eleT_coll, eleL_coll, eleV_coll,
																					 jetT_coll, bjet_coll, METv);
				if (channel == "" || channel == "emu_ttbar_candidate") return;

				// loop over systematics and make histograms
				double weight = 1.;
				weight = getWeight(ev, muonT_coll, muonL_coll, muonV_coll,
														   eleT_coll, eleL_coll, eleV_coll,
														   jetT_coll, bjet_coll, METv,
														   syst);
				FillMuons(channel+"/"+syst+"/muons", muonT_coll, weight, false);
				FillElectrons(channel+"/"+syst+"/electrons", eleT_coll, weight, false);
				FillJets(channel+"/"+syst+"/jets", jetT_coll, weight);
				FillJets(channel+"/"+syst+"/bjets", bjet_coll, weight);
				FillObject(channel+"/"+syst+"/METv", METv, weight);
				if (channel == "dimuon_drellyan_candidate") {
						Particle diMuon = muonT_coll.at(0)+muonT_coll.at(1);
						FillObject(channel+"/"+syst+"/diMuon", diMuon, weight);
				}																																        
		}
		return;
}


diLepControlRegion::channel_t diLepControlRegion::selectEvent(Event &ev,
																													   vector<Muon> &muonT_coll,
																														 vector<Muon> &muonL_coll,
																														 vector<Muon> &muonV_coll,
																														 vector<Electron> &eleT_coll,
																														 vector<Electron> &eleL_coll,
																														 vector<Electron> &eleV_coll,
																														 vector<Jet> &jetT_coll,
																														 vector<Jet> &bjet_coll,
																														 Particle &METv) {
		// Separate Regions by DblMuon and EMu first
		if (muonT_coll.size() == 2 && muonV_coll.size() == 2 && eleT_coll.size() == 0 && eleV_coll.size() == 0) {
				if (! ev.PassTrigger(DblMuTriggers)) return "";
				if (! (muonT_coll.at(0).Charge()+muonT_coll.at(1).Charge() == 0)) return "";
				if (! (muonT_coll.at(0).Pt() > 20.)) return "";
				if (! (muonT_coll.at(1).Pt() > 10.)) return "";

				const Particle ZCand = muonT_coll.at(0) + muonT_coll.at(1);
				if (fabs(ZCand.M() - 91.2) < 15.) {
						// DY
						if (bjet_coll.size() != 0) return "";
						return "dimuon_drellyan_candidate";
				}
				else {
						// TT
						if (! (ZCand.M() > 12.)) return "";
						if (! (jetT_coll.size() >= 2)) return "";
						if (! (bjet_coll.size() >= 1)) return "";
						if (! (METv.Pt() > 40.)) return "";
						return "dimuon_ttbar_candidate";
				}
		}
		else if (muonT_coll.size() == 1 && muonV_coll.size() == 1 && eleT_coll.size() == 1 && eleV_coll.size() == 1) {
				if (! ev.PassTrigger(EMuTriggers)) return "";
				const Muon &mu = muonT_coll.at(0);
				const Electron &ele = eleT_coll.at(0);
				if (! (mu.Charge() + ele.Charge() == 0)) return "";
				if (! ((mu.Pt() > 10. && ele.Pt() > 25.) || (mu.Pt() > 25. && ele.Pt() > 15.))) return "";
				if (! (jetT_coll.size() >= 2 )) return "";
				if (! (bjet_coll.size() >= 1 )) return "";
				return "emu_ttbar_candidate";
		}
		else {
				return "";
		}
}

double diLepControlRegion::getWeight(Event &ev,
																		 vector<Muon> &muonT_coll,
																		 vector<Muon> &muonL_coll,
																		 vector<Muon> &muonV_coll,
																		 vector<Electron> &eleT_coll,
																		 vector<Electron> &eleL_coll,
																		 vector<Electron> &eleV_coll,
																		 vector<Jet> &jetT_coll,
																		 vector<Jet> &bjet_coll,
																		 Particle &METv,
																		 diLepControlRegion::syst_t syst) {
		if (IsDATA) return 1.;
		double weight = 1.;
		double w_prefire = GetPrefireWeight(0);
		double w_lumi = ev.GetTriggerLumi("Full");
		double w_gen = MCweight();
		double w_pileup = GetPileUpWeight(nPileUp, 0);

		double w_idsf = mcCorr->GetMuonIDSF(MuonIDs.at(0), muonT_coll, 0);
		double w_trigsf = mcCorr->GetDoubleMuonTriggerSF(MuonIDs.at(0), muonT_coll, 0);

		// Set ID / Trigger scale factor
		if (syst == "Default") {
				w_idsf = 1.;
				w_trigsf = 1.;
		}
		else if (syst == "MuonIDSFOnly") {
				w_idsf = mcCorr->GetMuonIDSF(MuonIDs.at(0), muonT_coll, 0);
		}
		else if (syst == "MuonIDSFOnlyUp") {
				w_idsf = mcCorr->GetMuonIDSF(MuonIDs.at(0), muonT_coll, 1);
				w_trigsf = 1.;
		}
		else if (syst == "MuonIDSFOnlyDown") {
				w_idsf = mcCorr->GetMuonIDSF(MuonIDs.at(0), muonT_coll, -1);
				w_trigsf = -1.;
		}
		else if (syst == "Central" || syst == "MuonEnUp" || syst == "MuonEnDown") {
		}
		else if (syst == "PileUpCorrUp") {
				w_pileup = GetPileUpWeight(nPileUp, 1);
		}
		else if (syst == "PileUpCorrDown") {
				w_pileup = GetPileUpWeight(nPileUp, -1);
		}
		else if (syst == "MuonIDSFUp") {
				w_idsf = mcCorr->GetMuonIDSF(MuonIDs.at(0), muonT_coll, 1);
		}
		else if (syst == "MuonIDSFDown") {
				w_idsf = mcCorr->GetMuonIDSF(MuonIDs.at(0), muonT_coll, -1);
		}
		else if (syst == "MuonTrigSFUp") {
				w_trigsf = mcCorr->GetDoubleMuonTriggerSF(MuonIDs.at(0), muonT_coll, 1);
		}
		else if (syst == "MuonTrigSFDown") {
				w_trigsf = mcCorr->GetDoubleMuonTriggerSF(MuonIDs.at(0), muonT_coll, -1);
		}
		else {
				cerr << "[diLepControlRegion::getWeight] Wrong systematic " << syst << endl;
				exit(EXIT_FAILURE);
		}
		// electrons... not yet
		//
		weight *= (w_prefire*w_gen*w_lumi*w_pileup);
		weight *= (w_idsf*w_trigsf);

		// b-tagging weight
		JetTagging::Parameters jtp_DeepJet_Medium 
						= JetTagging::Parameters(JetTagging::DeepJet,
																		  JetTagging::Medium,
																			JetTagging::incl,
																			JetTagging::mujets);
		const double w_btag = mcCorr->GetBTaggingReweight_1a(jetT_coll, jtp_DeepJet_Medium);
		weight *= w_btag;
		
		//cout << "syst: " << syst << endl;
		//cout << "w_prefire: " << w_prefire << endl;
		//cout << "w_lumi: " << w_lumi << endl;
		//cout << "w_gen: " << w_gen << endl;
		//cout << "w_pileup: " << w_pileup << endl;
		//cout << "w_idsf: " << w_idsf << endl;
		//cout << "w_trigsf: " << w_trigsf << endl;
		//cout << "w_btag: " << w_btag << endl;
		//cout << "weight: " << weight << endl;
		
		return weight;
}

