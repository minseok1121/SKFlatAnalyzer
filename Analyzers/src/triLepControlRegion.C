#include "triLepControlRegion.h"

triLepControlRegion::triLepControlRegion() {}
triLepControlRegion::~triLepControlRegion() {}

void triLepControlRegion::initializeAnalyzer(){
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

		// get histogram for fake estimation
		TString datapath = getenv("DATA_DIR");
		TString fakepath = datapath+"/"+GetEra()+"/FakeRate/Muon/fakerate.root";
		TFile* f = new TFile(fakepath);
		h_fake = (TH2D*)f->Get("FR_cent_TopHNT_TopHNL");
		h_fake->SetDirectory(0);
		f->Close();
		if (!h_fake) {
				cerr << "[triLepControlRegion::InitializeAnalyzer] No histogram for fake rate" << endl;
				cerr << "[triLepControlRegion::InitializeAnalyzer] fakepath = " << fakepath << endl;
				exit(EXIT_FAILURE);
		}

		// B-tagging
		vector<JetTagging::Parameters> jtps;
		jtps.emplace_back(
						JetTagging::Parameters(JetTagging::DeepJet, JetTagging::Medium, JetTagging::incl, JetTagging::mujets)
						);
		mcCorr->SetJetTaggingParameters(jtps);

		// Systemtatics
		Systematics = {"Central",
								   "PileUpCorrUp", "PileUpCorrDown",
									 "MuonEnUp", "MuonEnDown",
									 "MuonIDSFUp", "MuonIDSFDown",
									 "MuonTrigSFUp", "MuonTrigSFDown",
									 "SystUpHtagUncorr", "SystDownHtagUncorr"
								  };
		FakeSystematics = {"Central", "MuonFakeUp", "MuonFakeDown"};
}

void triLepControlRegion::executeEvent(){

		if (! PassMETFilter()) return;

		Event ev = GetEvent();
		vector<Gen> truth_coll = GetGens();
		vector<Muon> RawMuons = GetAllMuons();
		vector<Electron> RawElectrons = GetAllElectrons();
		vector<Jet> RawJets = GetAllJets();
		Particle METv = ev.GetMETVector();

		for (const auto &syst: Systematics) {
				if (IsDATA && syst != "Central") continue;

				vector<Muon> muon_coll = RawMuons;
				vector<Electron> ele_coll = RawElectrons;
				vector<Jet> jet_coll = RawJets;
				if (syst == "MuonEnUp") { muon_coll = ScaleMuons(muon_coll, +1); }
				if (syst == "MuonEnDown") { muon_coll = ScaleMuons(muon_coll, -1); }

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
				if (channel == "") return;

				// loop over systematics and make histograms
				double weight = getWeight(ev, muonT_coll, muonL_coll, muonV_coll,
																		  eleT_coll, eleL_coll, eleV_coll,
																			jetT_coll, bjet_coll, METv, syst);
				
				// For MC && Data with tight leptons
				if (muonT_coll.size() == muonL_coll.size() && eleT_coll.size() == eleL_coll.size()) {
						FillMuons(channel+"/"+syst+"/muons", muonT_coll, weight, false);
						FillElectrons(channel+"/"+syst+"/electrons", eleT_coll, weight, false);
						FillJets(channel+"/"+syst+"/jets", jetT_coll, weight);
						FillJets(channel+"/"+syst+"/bjets", bjet_coll, weight);
						FillObject(channel+"/"+syst+"/METv", METv, weight);
				}
				else {
						for (auto &fsyst: FakeSystematics) {
								double evt_weight = getFakeWeight(muonL_coll, fsyst);
								FillMuons(channel+"/"+syst+"/"+fsyst+"/muons", muonT_coll, evt_weight, false);
								FillElectrons(channel+"/"+syst+"/"+fsyst+"/electrons", eleT_coll, evt_weight, false);
								FillJets(channel+"/"+syst+"/"+fsyst+"/jets", jetT_coll, evt_weight);
								FillJets(channel+"/"+syst+"/"+fsyst+".bjets", bjet_coll, evt_weight);
								FillObject(channel+"/"+syst+"/"+fsyst+"/METv", METv, evt_weight);	
						}
				}

				// MCsample depedent variables
				if ((!IsDATA) && (MCSample.Contains("TTToHcToWA"))) {
						if (channel == "signal_region_3mu") {
								const Particle mu1 = muonT_coll.at(0);
								const Particle mu2 = muonT_coll.at(1);
								const Particle mu3 = muonT_coll.at(2);
								Particle pair1, pair2;
								Particle MuonFromOffshellW;
								if (mu1.Charge() == mu2.Charge()) {
										pair1 = mu1 + mu3;
										pair2 = mu2 + mu3;
								}
								else if (mu1.Charge() == mu3.Charge()) {
										pair1 = mu1 + mu2;
										pair2 = mu2 + mu3;
								}
								else {  // mu2.Charge () == mu3.Charge()
										pair1 = mu1 + mu2;
										pair2 = mu1 + mu3;
								}
								// match pair to A
								Particle AMatched, NonAMatched;
								for (auto &gen: truth_coll) {
										if (! (gen.PID() == 36)) continue;
										FillObject(channel+"/"+syst+"/AGen", gen, weight);
										if (gen.DeltaR(pair1) < 0.1) {
												AMatched = pair1;
												NonAMatched = pair2;
										}
										else if (gen.DeltaR(pair2) < 0.1) {
												AMatched = pair2;
												NonAMatched = pair1;
										}
										else {
												continue;
										}
								}
								FillObject(channel+"/"+syst+"/AMatched", AMatched, weight);
								FillObject(channel+"/"+syst+"/NonAMatched", NonAMatched, weight);

								// Find muon from offshell W -> which is matched with Hc == 37
								for (auto &gen: truth_coll) {
										if (! (fabs(gen.PID()) == 37)) continue;
										FillObject(channel+"/"+syst+"/HcGen", gen, weight);
										if (gen.DeltaR(mu1) < 0.1) {
												MuonFromOffshellW = mu1;
										}
										else if (gen.DeltaR(mu2) < 0.1) {
												MuonFromOffshellW = mu2;
										}
										else if (gen.DeltaR(mu3) < 0.1) {
												MuonFromOffshellW = mu3;
										}
										else {
												continue;
										}
								}
								FillObject(channel+"/"+syst+"/MuonFromOffshellW", MuonFromOffshellW, weight);
						}
				}					
		}				
}

triLepControlRegion::channel_t triLepControlRegion::selectEvent(Event &ev,
																																vector<Muon> &muonT_coll,
																																vector<Muon> &muonL_coll,
																																vector<Muon> &muonV_coll,
																																vector<Electron> &eleT_coll,
																																vector<Electron> &eleL_coll,
																																vector<Electron> &eleV_coll,
																																vector<Jet> &jetT_coll,
																																vector<Jet> &bjet_coll,
																																Particle &METv) {
		// Separate regions by lepton multiplicity first
		const bool is3Mu = (muonL_coll.size() == 3 
										 && muonV_coll.size() == 3 
										 && eleL_coll.size() == 0 
										 && eleV_coll.size() == 0);
		const bool is1E2Mu = (muonL_coll.size() == 2
										   && muonV_coll.size() == 2
											 && eleL_coll.size() == 1
											 && eleV_coll.size() == 1);

		// At this point, only tight leptons remain in MC while data contains loose leptons
		if (!IsDATA) {
				if (muonT_coll.size() != muonL_coll.size()) return "";
				if (eleT_coll.size() != eleT_coll.size()) return "";
		}

		if (is1E2Mu) {
				// will be set later
				return "";
		}
		else if (is3Mu) {
				// baseline selection
				// 1. Pass DoubleMuon Triggers
				// 2. Exact 3 tight muons, no additional lepton
				// 3. Exist OS muon pair (|charge sum| == 1)
				// 4. All OS muon pair mass > 12GeV
				// For MC, Lcoll == Tcoll
				if (! ev.PassTrigger(DblMuTriggers)) return "";
				const Muon &mu1 = muonL_coll.at(0);
				const Muon &mu2 = muonL_coll.at(1);
				const Muon &mu3 = muonL_coll.at(2);
				if (! (mu1.Pt() > 20.)) return "";
				if (! (mu2.Pt() > 10.)) return "";
				if (! (mu3.Pt() > 10.)) return "";
				if (! (abs(mu1.Charge() + mu2.Charge() + mu3.Charge()) == 1)) return "";
				// make dimuon pairs
				Particle pair1, pair2;
				if (mu1.Charge() == mu2.Charge()) {
						pair1 = mu1 + mu3;
						pair2 = mu2 + mu3;
				}
				else if (mu1.Charge() == mu3.Charge()) {
						pair1 = mu1 + mu2;
						pair2 = mu2 + mu3;
				}
				else {	// mu2.Charge () == mu3.Charge()
						pair1 = mu1 + mu2;
						pair2 = mu1 + mu3;
				}
				if (! (pair1.M() > 12.)) return "";
				if (! (pair2.M() > 12.)) return "";

				// Orthogonality of Signal Region & Z+fake/Z+gamma region by bjet multiplicity
				// SR: Nj >= 2, Nb >= 1
				if (bjet_coll.size() >= 1) {
						if (jetT_coll.size() >= 2) return "signal_region_3mu";
						else return "";
				}
				else {
						// Z+fake: exist on-shell Z mass pair, Nb = 0
						// Z+gamma: no on-shell Z mass pair, M(3mu) on Z mass, Nb = 0
						// check whether pairs are on Z mass
						const bool isOnZ = ((fabs(pair1.M() - 91.2) < 10. || fabs(pair2.M() - 91.2) < 10.));
						if (isOnZ) return "z_fake_region_3mu";
						else {
								if (fabs((mu1+mu2+mu3).M() - 91.2) < 10.) return "z_gamma_region_3mu";
								else return "";
						}
				}
		}
		else {
				return "";
		}
}

double triLepControlRegion::getWeight(Event &ev,
														  vector<Muon> &muonT_coll,
                              vector<Muon> &muonL_coll,
                              vector<Muon> &muonV_coll,
                              vector<Electron> &eleT_coll,
                              vector<Electron> &eleL_coll,
                              vector<Electron> &eleV_coll,
                              vector<Jet> &jetT_coll,
                              vector<Jet> &bjet_coll,
                              Particle &METv,
                              triLepControlRegion::syst_t syst) {
		if (IsDATA) return 1.;
		double weight = 1.;
		double w_prefire = GetPrefireWeight(0);
		double w_lumi = ev.GetTriggerLumi("Full");
		double w_gen = MCweight();
		double w_pileup = GetPileUpWeight(nPileUp, 0);

		double w_idsf = mcCorr->GetMuonIDSF(MuonIDs.at(0), muonT_coll, 0);
		double w_trigsf = mcCorr->GetDoubleMuonTriggerSF(MuonIDs.at(0), muonT_coll, 0);
		
		JetTagging::Parameters jtp_DeepJet_Medium
				    = JetTagging::Parameters(JetTagging::DeepJet,
														         JetTagging::Medium,
																		 JetTagging::incl,
																		 JetTagging::mujets);
		double w_btag = mcCorr->GetBTaggingReweight_1a(jetT_coll, jtp_DeepJet_Medium);

		if (syst == "PileUpCorrUp")
				w_pileup = GetPileUpWeight(nPileUp, +1);
		if (syst == "PileUpCorrDown")
				w_pileup = GetPileUpWeight(nPileUp, -1);
		if (syst == "MuonIDSFUp")
				w_idsf = mcCorr->GetMuonIDSF(MuonIDs.at(0), muonT_coll, +1);
		if (syst == "MuonIDSFDown")
				w_idsf = mcCorr->GetMuonIDSF(MuonIDs.at(0), muonT_coll, -1);
		if (syst == "MuonTrigSFUp")
				w_trigsf = mcCorr->GetDoubleMuonTriggerSF(MuonIDs.at(0), muonT_coll, +1);
		if (syst == "MuonTrigSFDown")
				w_trigsf = mcCorr->GetDoubleMuonTriggerSF(MuonIDs.at(0), muonT_coll, -1);
		if (syst.Contains("HTag")) {
				string systString(syst.Data());
				w_btag = mcCorr->GetBTaggingReweight_1a(jetT_coll, jtp_DeepJet_Medium, systString);
		}

		weight *= (w_prefire*w_gen*w_lumi*w_pileup);
		weight *= (w_idsf*w_trigsf);
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

double triLepControlRegion::getFakeWeight(vector<Muon> &muonL_coll, triLepControlRegion::syst_t fsyst) {
		double out = -1.;
		for (auto &mu: muonL_coll) {
				if (mu.PassID(MuonIDs.at(0))) {
						continue;
				}
				else {
						double PtCorr = mu.Pt() * (1.+max(0., mu.MiniRelIso() - 0.1));
						double fr = h_fake->GetBinContent(h_fake->FindBin(PtCorr, fabs(mu.Eta())));
						if (fsyst == "MuonFakeUp") fr *= 1.3;
						if (fsyst == "MuonFakeDown") fr *= 0.7;
						out *= -1.*(fr / (1.+fr));
				}
		}

		return out;
}

















