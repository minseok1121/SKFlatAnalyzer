#include "FakeMeasurement.h"

FakeMeasurement::FakeMeasurement() {}
FakeMeasurement::~FakeMeasurement() {}

void FakeMeasurement::initializeAnalyzer(){
	// flags
	MeasMuon = HasFlag("MeasMuon");
	MeasElectron = HasFlag("MeasElectron");
	cout << "[FakeMeasurement::initializeAnalyzer] MeasMuon = " << MeasMuon << endl;
	cout << "[FakeMeasurement::initializeAnalyzer] MeasElectron = " << MeasElectron << endl;
	if (! (MeasMuon || MeasElectron)) {
		cerr << "At least one of the flag should be specified" << endl;
		exit(EXIT_FAILURE);
	}

	// Triggers
    if (DataYear == 2017) {
        trigs_sglmu.clear();
        trigs_sglmu.emplace_back("HLT_Mu8_TrkIsoVVL_v");    // 2.605*1.33461
        trigs_sglmu.emplace_back("HLT_Mu17_TrkIsoVVL_v"); // 70.039
        trigs_sglele.clear();
        trigs_sglele.emplace_back("HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_v");
        trigs_sglele.emplace_back("HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v");
        trigs_sglele.emplace_back("HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v");
    }
    else {
        cerr << "No trigger set for DataYear " << DataYear << endl;
        exit(EXIT_FAILURE);
    }

    // ID
    HcToWA_MuID = {"HcToWATight", "HcToWALoose"};
    HcToWA_EleID = {"HcToWATight", "HcToWALoose"};

    // B-tagging (for systematics)
    vector<JetTagging::Parameters> jtps;
    jtps.emplace_back(JetTagging::Parameters(JetTagging::DeepJet, JetTagging::Medium, JetTagging::incl, JetTagging::mujets));
    mcCorr->SetJetTaggingParameters(jtps);

	// Systematics
    Systematics = {
			"Central", "JetPtCutUp", "JetPtCutDown", 
			"JetResUp", "JetResDown", "JetEnUp", "JetEnDown",
			"MuonEnUp", "MuonEnDown", "ElectronResUp", "ElectronResDown",
			"ElectronEnUp", "ElectronEnDown", "HasBjet"
	};
}

void FakeMeasurement::executeEvent(){
	// Get All Objects
	muons_all = GetAllMuons();
	electrons_all = GetAllElectrons();
	jets_all = GetAllJets();

	for (const auto& syst: Systematics)
		executeEventWithSystematics(syst);
}

void FakeMeasurement::executeEventWithSystematics(const TString& syst) {
	if (!PassMETFilter()) return;

	Event ev = GetEvent();
	Particle METv = ev.GetMETVector();
	vector<Muon> muons = muons_all;
	vector<Electron> electrons = electrons_all;
	vector<Jet> jets = jets_all;

	// Jet Resolution, Energy Scale
    if (syst == "JetResUp")			jets = SmearJets(jets, +1);
    if (syst == "JetResDown")		jets = SmearJets(jets, -1);
    if (syst == "JetEnUp")			jets = ScaleJets(jets, +1);
    if (syst == "JetEnDown")		jets = ScaleJets(jets, -1);
    if (syst == "MuonEnUp")			muons = ScaleMuons(muons, +1);
    if (syst == "MuonEnDown")		muons = ScaleMuons(muons, -1);
    if (syst == "ElectronResUp")	electrons = SmearElectrons(electrons, +1);
    if (syst == "ElectronResDown")	electrons = SmearElectrons(electrons, -1);
    if (syst == "ElectronEnUp")		electrons = ScaleElectrons(electrons, +1);
    if (syst == "ElectronEnDown")	electrons = ScaleElectrons(electrons, -1);

	sort(muons.begin(), muons.end(), PtComparing);
	sort(electrons.begin(), electrons.end(), PtComparing);
	sort(jets.begin(), jets.end(), PtComparing);

	vector<Muon> muons_loose = SelectMuons(muons, "HcToWALoose", 10., 2.4);
	vector<Electron> electrons_loose = SelectElectrons(electrons, "HcToWALoose", 10., 2.5);
	vector<Jet> jets_tight = SelectJets(jets, "tight", 20., 2.4);
	vector<Jet> jets_lepVeto = JetsVetoLeptonInside(jets_tight, electrons_loose, muons_loose, 0.4);
	if (syst == "HasBjet") {
		vector<Jet> temp;
		for (const auto& jet: jets_lepVeto) {
			const double this_discr = jet.GetTaggerResult(JetTagging::DeepJet);
			if (this_discr > mcCorr->GetJetTaggingCutValue(JetTagging::DeepJet, JetTagging::Medium))
				temp.emplace_back(jet);
		}
		jets_lepVeto = temp;
	}
	// Event Selection
	// Fire the trigger for the given flavor of the loose lepton
	// Exactly one loose lepton present in the event with p_T > 10 GeV
	// At least one jet of p_T > 40 GeV with |\eta| < 2.4
	// deltaR(j, l) > 1.0
	TString channel = "";
	
	// JetPtCuts
	if (! (jets_lepVeto.size() > 0)) return;
	if (MeasMuon) {
		if (syst == "JetPtCutUp") {
			if (! (jets_lepVeto.at(0).Pt() > 60.)) return;
		}
		else if (syst == "JetPtCutDown") {
			if (! (jets_lepVeto.at(0).Pt() > 20.)) return;
		}
		else {
			if (! (jets_lepVeto.at(0).Pt() > 40.)) return;
		}
	}
	if (MeasElectron) {
		if (syst == "JetPtCutUp") {
			if (! (jets_lepVeto.at(0).Pt() > 60.)) return;
		}
		else if (syst == "JetPtCutDown") {
			if (! (jets_lepVeto.at(0).Pt() > 30.)) return;
		}
		else {
			if (! (jets_lepVeto.at(0).Pt() > 40.)) return;
		}
	}

	// Leptons
	if (muons_loose.size() == 1 && electrons_loose.size() == 0) channel = "SglMu";
    if (muons_loose.size() == 0 && electrons_loose.size() == 1) channel = "SglEle";
    if (muons_loose.size() == 2 && electrons_loose.size() == 0) channel = "DblMu";
    if (muons_loose.size() == 0 && electrons_loose.size() == 2) channel = "DblEle";
    if (channel == "") return;		

	if (channel == "SglMu")
		if (! (muons_loose.at(0).DeltaR(jets_lepVeto.at(0)) > 1.0)) return;
	if (channel == "SglEle")
		if (! (electrons_loose.at(0).DeltaR(jets_lepVeto.at(0)) > 1.0)) return;

	// Devide regions
	if (MeasMuon) {
		if (channel.Contains("Ele")) return;
		
		bool passMu8Path = ev.PassTrigger("HLT_Mu8_TrkIsoVVL_v");
		bool passMu17Path = ev.PassTrigger("HLT_Mu17_TrkIsoVVL_v");
		if (! (passMu8Path || passMu17Path)) return;

		if (channel == "SglMu") {
			// only need Central, JetPtCut, HasBjet
			if (! (syst == "Central" || syst.Contains("JetPtCut") || syst == "HasBjet")) return;
			const Muon& muon = muons_loose.at(0);
			const double ptCorr = muon.Pt()*(1.+max(muon.MiniRelIso()-0.1, 0.));
			const double absEta = fabs(muon.Eta());
			const double Mt = MT(muon, METv);
			const double MET = METv.Pt();
			double ptCorrBin[] = {10., 15., 20., 30., 50., 70., 100.};
			double absEtaBin[] = {0., 0.9, 1.6, 2.4};

			TString region = "";
			if (Mt > 80. && MET > 50.) region = "PromptEnriched";
			if (Mt < 25. && MET < 25.) region = "QCDEnriched";

			TString trig = "";
			TString ID = "";
			TString histkey = "";
			if (passMu8Path) {
				trig = "HLT_Mu8_TrkIsoVVL";
				ID = "PassLoose";
				
				// set weight
				double weight = 1.;
				if (!IsDATA) {
					double genWeight = ev.MCweight()*weight_norm_1invpb;
					double L1PrefireWeight = GetPrefireWeight(0);
					double pileupWeight = GetPileUpWeight(nPileUp, 0);
					weight = genWeight*L1PrefireWeight*pileupWeight*2.8977;
				}
				// Fill Inclusive
				histkey = channel+"/"+trig+"/"+ID+"/"+syst+"/Incl/";
				FillObjects(histkey+"muons", muons_loose, weight);
				FillObjects(histkey+"jets_lepVeto", jets_lepVeto, weight);
				FillHist(histkey+"MET", MET, weight, 200, 0., 200.);
				FillHist(histkey+"MT", Mt, weight, 200, 0., 200.);

				// Fill PromptEnriched
				if (region == "PromptEnriched") {
					histkey = channel+"/"+trig+"/"+ID+"/"+syst+"/"+region+"/";
					FillObjects(histkey+"muons", muons_loose, weight);
					FillObjects(histkey+"jets_lepVeto", jets_lepVeto, weight);
					FillHist(histkey+"MET", MET, weight, 200, 0., 200.);
					FillHist(histkey+"MT", Mt, weight, 200, 0., 200.);
				}
				if (region == "QCDEnriched") {
					histkey = channel+"/"+trig+"/"+ID+"/"+syst+"/"+region+"/";
					FillObjects(histkey+"muons", muons_loose, weight);
					FillObjects(histkey+"jets_lepVeto", jets_lepVeto, weight);
					FillHist(histkey+"MET", MET, weight, 200, 0., 200.);
					FillHist(histkey+"MT", Mt, weight, 200, 0., 200.);
					if (ptCorr < 30.) {
						histkey = channel+"/"+ID+"/"+syst;
						FillHist(histkey, ptCorr, absEta, weight, 6, ptCorrBin, 3, absEtaBin);
					}
				}

				if (muon.PassID("HcToWATight")) {
					ID = "PassTight";
					// Fill Inclusive
                	histkey = channel+"/"+trig+"/"+ID+"/"+syst+"/Incl/";
                	FillObjects(histkey+"muons", muons_loose, weight);
                	FillObjects(histkey+"jets_lepVeto", jets_lepVeto, weight);
                	FillHist(histkey+"MET", MET, weight, 200, 0., 200.);
                	FillHist(histkey+"MT", Mt, weight, 200, 0., 200.);

                	// Fill PromptEnriched
                	if (region == "PromptEnriched") {
                    	histkey = channel+"/"+trig+"/"+ID+"/"+syst+"/"+region+"/";
                    	FillObjects(histkey+"muons", muons_loose, weight);
                   	 	FillObjects(histkey+"jets_lepVeto", jets_lepVeto, weight);
                    	FillHist(histkey+"MET", MET, weight, 200, 0., 200.);
                    	FillHist(histkey+"MT", Mt, weight, 200, 0., 200.);
                	}
                	if (region == "QCDEnriched") {
                    	histkey = channel+"/"+trig+"/"+ID+"/"+syst+"/"+region+"/";
                    	FillObjects(histkey+"muons", muons_loose, weight);
                    	FillObjects(histkey+"jets_lepVeto", jets_lepVeto, weight);
                    	FillHist(histkey+"MET", MET, weight, 200, 0., 200.);
                    	FillHist(histkey+"MT", Mt, weight, 200, 0., 200.);
						if (ptCorr < 30.) {
							histkey = channel+"/"+ID+"/"+syst;
							FillHist(histkey, ptCorr, absEta, weight, 6, ptCorrBin, 3, absEtaBin);
						}
					}
				}
			}
			if (passMu17Path) {
                trig = "HLT_Mu17_TrkIsoVVL";
                ID = "PassLoose";

                // set weight
                double weight = 1.;
                if (!IsDATA) {
					double genWeight = ev.MCweight()*weight_norm_1invpb;
					double L1PrefireWeight = GetPrefireWeight(0);
					double pileupWeight = GetPileUpWeight(nPileUp, 0);
					weight = genWeight*L1PrefireWeight*pileupWeight*65.8989;
                }
                // Fill Inclusive
                histkey = channel+"/"+trig+"/"+ID+"/"+syst+"/Incl/";
                FillObjects(histkey+"muons", muons_loose, weight);
                FillObjects(histkey+"jets_lepVeto", jets_lepVeto, weight);
                FillHist(histkey+"MET", MET, weight, 200, 0., 200.);
                FillHist(histkey+"MT", Mt, weight, 200, 0., 200.);

                // Fill PromptEnriched
                if (region == "PromptEnriched") {
                    histkey = channel+"/"+trig+"/"+ID+"/"+syst+"/"+region+"/";
                    FillObjects(histkey+"muons", muons_loose, weight);
                    FillObjects(histkey+"jets_lepVeto", jets_lepVeto, weight);
                    FillHist(histkey+"MET", MET, weight, 200, 0., 200.);
                    FillHist(histkey+"MT", Mt, weight, 200, 0., 200.);
                }
                if (region == "QCDEnriched") {
                    histkey = channel+"/"+trig+"/"+ID+"/"+syst+"/"+region+"/";
                    FillObjects(histkey+"muons", muons_loose, weight);
                    FillObjects(histkey+"jets_lepVeto", jets_lepVeto, weight);
                    FillHist(histkey+"MET", MET, weight, 200, 0., 200.);
                    FillHist(histkey+"MT", Mt, weight, 200, 0., 200.);
					if (ptCorr > 30.) {
						histkey = channel+"/"+ID+"/"+syst;
						FillHist(histkey, ptCorr, absEta, weight, 6, ptCorrBin, 3, absEtaBin);
					}
				}
				
				if (muon.PassID("HcToWATight")) {
                    ID = "PassTight";
                    // Fill Inclusive
                    histkey = channel+"/"+trig+"/"+ID+"/"+syst+"/Incl/";
                    FillObjects(histkey+"muons", muons_loose, weight);
                    FillObjects(histkey+"jets_lepVeto", jets_lepVeto, weight);
                    FillHist(histkey+"MET", MET, weight, 200, 0., 200.);
                    FillHist(histkey+"MT", Mt, weight, 200, 0., 200.);

                    if (region == "PromptEnriched") {
                        histkey = channel+"/"+trig+"/"+ID+"/"+syst+"/"+region+"/";
                        FillObjects(histkey+"muons", muons_loose, weight);
                        FillObjects(histkey+"jets_lepVeto", jets_lepVeto, weight);
                        FillHist(histkey+"MET", MET, weight, 200, 0., 200.);
                        FillHist(histkey+"MT", Mt, weight, 200, 0., 200.);
                    }
                    if (region == "QCDEnriched") {
                        histkey = channel+"/"+trig+"/"+ID+"/"+syst+"/"+region+"/";
                        FillObjects(histkey+"muons", muons_loose, weight);
                        FillObjects(histkey+"jets_lepVeto", jets_lepVeto, weight);
                        FillHist(histkey+"MET", MET, weight, 200, 0., 200.);
                        FillHist(histkey+"MT", Mt, weight, 200, 0., 200.);
						if (ptCorr > 30.) {
							histkey = channel+"/"+ID+"/"+syst;
							FillHist(histkey, ptCorr, absEta, weight, 6, ptCorrBin, 3, absEtaBin);
						}
					}
                }
            }
		}
		if (channel == "DblMu") {
			if (syst.Contains("JetPtCut") || syst == "HasBjet") return;
			const Muon mu1 = muons_loose.at(0);
			const Muon mu2 = muons_loose.at(1);
			const Particle ZCand = mu1 + mu2;
			if (fabs(ZCand.M() - 91.2) > 15.) return;
			if (mu1.Charge() + mu2.Charge() != 0) return;

			TString trig = "";
			TString ID = "";
			TString histkey = "";
			if (passMu8Path) {
				trig = "HLT_Mu8_TrkIsoVVL";
				ID = "PassLoose";

				// set weight
				double weight = 1.;
				if (!IsData) {
					double genWeight = ev.MCweight()*weight_norm_1invpb;
					double L1PrefireWeight = GetPrefireWeight(0);
					double pileupWeight = GetPileUpWeight(nPileUp, 0);
					weight = genWeight*L1PrefireWeight*pileupWeight*2.8977;
				}

				histkey = channel+"/"+trig+"/"+ID+"/"+syst+"/";
				FillObjects(histkey+"muons", muons_loose, weight);
				FillObjects(histkey+"jets_lepVeto", jets_lepVeto, weight);
				FillObject(histkey+"ZCand", ZCand, weight);

				if (mu1.PassID("HcToWATight") && mu2.PassID("HcToWATight")) {
					ID = "PassTight";
					histkey = channel+"/"+trig+"/"+ID+"/"+syst+"/";
					FillObjects(histkey+"muons", muons_loose, weight);
					FillObjects(histkey+"jets_lepVeto", jets_lepVeto, weight);
					FillObject(histkey+"ZCand", ZCand, weight);
				}
			}
			if (passMu17Path) {
				trig = "HLT_Mu17_TrkIsoVVL";
                ID = "PassLoose";

                // set weight
                double weight = 1.;
                if (!IsData) {
					double genWeight = ev.MCweight()*weight_norm_1invpb;
					double L1PrefireWeight = GetPrefireWeight(0);
					double pileupWeight = GetPileUpWeight(nPileUp, 0);
					weight = genWeight*L1PrefireWeight*pileupWeight*65.8989;
                }

                histkey = channel+"/"+trig+"/"+ID+"/"+syst+"/";
                FillObjects(histkey+"muons", muons_loose, weight);
                FillObjects(histkey+"jets_lepVeto", jets_lepVeto, weight);
                FillObject(histkey+"ZCand", ZCand, weight);

                if (mu1.PassID("HcToWATight") && mu2.PassID("HcToWATight")) {
                    ID = "PassTight";
                    histkey = channel+"/"+trig+"/"+ID+"/"+syst+"/";
                    FillObjects(histkey+"muons", muons_loose, weight);
                    FillObjects(histkey+"jets_lepVeto", jets_lepVeto, weight);
                    FillObject(histkey+"ZCand", ZCand, weight);
                }
            }
		}
	}
	if (MeasElectron) {
		if (channel.Contains("Mu")) return;

		bool passEle8Path = ev.PassTrigger("HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_v");
		bool passEle12Path = ev.PassTrigger("HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v");
		bool passEle23Path = ev.PassTrigger("HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v");
		if (! (passEle8Path || passEle12Path || passEle23Path)) return;

		if (channel == "SglEle") {
			// only need Central, JetPtCut, HasBjet
			if (! (syst == "Central" || syst.Contains("JetPtCut")  || syst == "HasBjet")) return;
			const Electron& ele = electrons_loose.at(0);
			const double ptCorr = ele.Pt()*(1.+max(ele.MiniRelIso()-0.1, 0.));
			const double absEta = fabs(ele.Eta());
			const double Mt = MT(ele, METv);
			const double MET = METv.Pt();
			//double ptCorrBin[] = {10., 15., 20., 25., 35., 50., 70., 100.};
			double ptCorrBin[] = {10., 20., 35., 50., 70};
			double absEtaBin[] = {0., 0.8, 1.479, 2.5};

			TString region = "";
			if (Mt > 80. && MET > 50.) region = "PromptEnriched";
			if (Mt < 25. && MET < 25.) region = "QCDEnriched";

			TString trig = "";
			TString ID = "";
			TString histkey = "";

			if (passEle8Path) {
				trig = "HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30";
				ID = "PassLoose";

				// set weight
				// 3.973, 27.699, 43.468
				double weight = 1.;
				if (!IsDATA) {
					double genWeight = ev.MCweight()*weight_norm_1invpb;
					double L1PrefireWeight = GetPrefireWeight(0);
					double pileupWeight = GetPileUpWeight(nPileUp, 0);
					weight = genWeight*L1PrefireWeight*pileupWeight*3.973;
				}

				// Fill Inclusive
                histkey = channel+"/"+trig+"/"+ID+"/"+syst+"/Incl/";
                FillObjects(histkey+"electrons", electrons_loose, weight);
                FillObjects(histkey+"jets_lepVeto", jets_lepVeto, weight);
                FillHist(histkey+"MET", MET, weight, 200, 0., 200.);
                FillHist(histkey+"MT", Mt, weight, 200, 0., 200.);

                // Fill PromptEnriched
                if (region == "PromptEnriched") {
                    histkey = channel+"/"+trig+"/"+ID+"/"+syst+"/"+region+"/";
                    FillObjects(histkey+"electrons", electrons_loose, weight);
                    FillObjects(histkey+"jets_lepVeto", jets_lepVeto, weight);
                    FillHist(histkey+"MET", MET, weight, 200, 0., 200.);
                    FillHist(histkey+"MT", Mt, weight, 200, 0., 200.);
                }
                if (region == "QCDEnriched") {
                    histkey = channel+"/"+trig+"/"+ID+"/"+syst+"/"+region+"/";
                    FillObjects(histkey+"electrons", electrons_loose, weight);
                    FillObjects(histkey+"jets_lepVeto", jets_lepVeto, weight);
                    FillHist(histkey+"MET", MET, weight, 200, 0., 200.);
                    FillHist(histkey+"MT", Mt, weight, 200, 0., 200.);
                    if (ptCorr < 20.) {
						histkey = channel+"/"+ID+"/"+syst;
                        FillHist(histkey, ptCorr, absEta, weight, 4, ptCorrBin, 3, absEtaBin);
					}
				}
				if (ele.PassID("HcToWATight")) {
                    ID = "PassTight";
                    // Fill Inclusive
                    histkey = channel+"/"+trig+"/"+ID+"/"+syst+"/Incl/";
                    FillObjects(histkey+"electrons", electrons_loose, weight);
                    FillObjects(histkey+"jets_lepVeto", jets_lepVeto, weight);
                    FillHist(histkey+"MET", MET, weight, 200, 0., 200.);
                    FillHist(histkey+"MT", Mt, weight, 200, 0., 200.);

                    // Fill PromptEnriched
                    if (region == "PromptEnriched") {
                        histkey = channel+"/"+trig+"/"+ID+"/"+syst+"/"+region+"/";
                        FillObjects(histkey+"electrons", electrons_loose, weight);
                        FillObjects(histkey+"jets_lepVeto", jets_lepVeto, weight);
                        FillHist(histkey+"MET", MET, weight, 200, 0., 200.);
                        FillHist(histkey+"MT", Mt, weight, 200, 0., 200.);
                    }
                    if (region == "QCDEnriched") {
                        histkey = channel+"/"+trig+"/"+ID+"/"+syst+"/"+region+"/";
                        FillObjects(histkey+"electrons", electrons_loose, weight);
                        FillObjects(histkey+"jets_lepVeto", jets_lepVeto, weight);
                        FillHist(histkey+"MET", MET, weight, 200, 0., 200.);
                        FillHist(histkey+"MT", Mt, weight, 200, 0., 200.);
						if (ptCorr < 20.) {
							histkey = channel+"/"+ID+"/"+syst;
							FillHist(histkey, ptCorr, absEta, weight, 4, ptCorrBin, 3, absEtaBin);
						}
					}
                }
            }
			if (passEle12Path) {
                trig = "HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30";
                ID = "PassLoose";

                // set weight
                // 3.973, 27.699, 43.468
                double weight = 1.; 
                if (!IsDATA) {
                    double genWeight = ev.MCweight()*weight_norm_1invpb;
                    double L1PrefireWeight = GetPrefireWeight(0);
                    double pileupWeight = GetPileUpWeight(nPileUp, 0); 
                    weight = genWeight*L1PrefireWeight*pileupWeight*27.699;
                }

                // Fill Inclusive
                histkey = channel+"/"+trig+"/"+ID+"/"+syst+"/Incl/";
                FillObjects(histkey+"electrons", electrons_loose, weight);
                FillObjects(histkey+"jets_lepVeto", jets_lepVeto, weight);
                FillHist(histkey+"MET", MET, weight, 200, 0., 200.);
                FillHist(histkey+"MT", Mt, weight, 200, 0., 200.);

                // Fill PromptEnriched
                if (region == "PromptEnriched") {
                    histkey = channel+"/"+trig+"/"+ID+"/"+syst+"/"+region+"/";
                    FillObjects(histkey+"electrons", electrons_loose, weight);
                    FillObjects(histkey+"jets_lepVeto", jets_lepVeto, weight);
                    FillHist(histkey+"MET", MET, weight, 200, 0., 200.);
                    FillHist(histkey+"MT", Mt, weight, 200, 0., 200.);
                }
                if (region == "QCDEnriched") {
                    histkey = channel+"/"+trig+"/"+ID+"/"+syst+"/"+region+"/";
                    FillObjects(histkey+"electrons", electrons_loose, weight);
                    FillObjects(histkey+"jets_lepVeto", jets_lepVeto, weight);
                    FillHist(histkey+"MET", MET, weight, 200, 0., 200.);
                    FillHist(histkey+"MT", Mt, weight, 200, 0., 200.);
                    if (20. < ptCorr && ptCorr < 35.) {
						histkey = channel+"/"+ID+"/"+syst;
                        FillHist(histkey, ptCorr, absEta, weight, 4, ptCorrBin, 3, absEtaBin);
					}
				}
				
				if (ele.PassID("HcToWATight")) {
                    ID = "PassTight";
                    // Fill Inclusive
                    histkey = channel+"/"+trig+"/"+ID+"/"+syst+"/Incl/";
                    FillObjects(histkey+"electrons", electrons_loose, weight);
                    FillObjects(histkey+"jets_lepVeto", jets_lepVeto, weight);
                    FillHist(histkey+"MET", MET, weight, 200, 0., 200.);
                    FillHist(histkey+"MT", Mt, weight, 200, 0., 200.);

                    // Fill PromptEnriched
                    if (region == "PromptEnriched") {
                        histkey = channel+"/"+trig+"/"+ID+"/"+syst+"/"+region+"/";
                        FillObjects(histkey+"electrons", electrons_loose, weight);
                        FillObjects(histkey+"jets_lepVeto", jets_lepVeto, weight);
                        FillHist(histkey+"MET", MET, weight, 200, 0., 200.);
                        FillHist(histkey+"MT", Mt, weight, 200, 0., 200.);
                    }
                    if (region == "QCDEnriched") {
                        histkey = channel+"/"+trig+"/"+ID+"/"+syst+"/"+region+"/";
                        FillObjects(histkey+"electrons", electrons_loose, weight);
                        FillObjects(histkey+"jets_lepVeto", jets_lepVeto, weight);
                        FillHist(histkey+"MET", MET, weight, 200, 0., 200.);
                        FillHist(histkey+"MT", Mt, weight, 200, 0., 200.);
						if (20. < ptCorr && ptCorr < 35.) {
							histkey = channel+"/"+ID+"/"+syst;
							FillHist(histkey, ptCorr, absEta, weight, 4, ptCorrBin, 3, absEtaBin);
						}
					}
                }
            }
			if (passEle23Path) {
                trig = "HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30";
                ID = "PassLoose";

                // set weight
                // 3.973, 27.699, 43.468
                double weight = 1.;
                if (!IsDATA) {
                    double genWeight = ev.MCweight()*weight_norm_1invpb;
                    double L1PrefireWeight = GetPrefireWeight(0);
                    double pileupWeight = GetPileUpWeight(nPileUp, 0);
                    weight = genWeight*L1PrefireWeight*pileupWeight*43.468;
                }

                // Fill Inclusive
                histkey = channel+"/"+trig+"/"+ID+"/"+syst+"/Incl/";
                FillObjects(histkey+"electrons", electrons_loose, weight);
                FillObjects(histkey+"jets_lepVeto", jets_lepVeto, weight);
                FillHist(histkey+"MET", MET, weight, 200, 0., 200.);
                FillHist(histkey+"MT", Mt, weight, 200, 0., 200.);

                // Fill PromptEnriched
                if (region == "PromptEnriched") {
                    histkey = channel+"/"+trig+"/"+ID+"/"+syst+"/"+region+"/";
                    FillObjects(histkey+"electrons", electrons_loose, weight);
                    FillObjects(histkey+"jets_lepVeto", jets_lepVeto, weight);
                    FillHist(histkey+"MET", MET, weight, 200, 0., 200.);
                    FillHist(histkey+"MT", Mt, weight, 200, 0., 200.);
                }
                if (region == "QCDEnriched") {
                    histkey = channel+"/"+trig+"/"+ID+"/"+syst+"/"+region+"/";
                    FillObjects(histkey+"electrons", electrons_loose, weight);
                    FillObjects(histkey+"jets_lepVeto", jets_lepVeto, weight);
                    FillHist(histkey+"MET", MET, weight, 200, 0., 200.);
                    FillHist(histkey+"MT", Mt, weight, 200, 0., 200.);
                    if (ptCorr > 35.) {
						histkey = channel+"/"+ID+"/"+syst;
                        FillHist(histkey, ptCorr, absEta, weight, 4, ptCorrBin, 3, absEtaBin);
					}
				}
				if (ele.PassID("HcToWATight")) {
                    ID = "PassTight";
                    // Fill Inclusive
                    histkey = channel+"/"+trig+"/"+ID+"/"+syst+"/Incl/";
                    FillObjects(histkey+"electrons", electrons_loose, weight);
                    FillObjects(histkey+"jets_lepVeto", jets_lepVeto, weight);
                    FillHist(histkey+"MET", MET, weight, 200, 0., 200.);
                    FillHist(histkey+"MT", Mt, weight, 200, 0., 200.);

                    // Fill PromptEnriched
                    if (region == "PromptEnriched") {
                        histkey = channel+"/"+trig+"/"+ID+"/"+syst+"/"+region+"/";
                        FillObjects(histkey+"electrons", electrons_loose, weight);
                        FillObjects(histkey+"jets_lepVeto", jets_lepVeto, weight);
                        FillHist(histkey+"MET", MET, weight, 200, 0., 200.);
                        FillHist(histkey+"MT", Mt, weight, 200, 0., 200.);
                    }
                    if (region == "QCDEnriched") {
                        histkey = channel+"/"+trig+"/"+ID+"/"+syst+"/"+region+"/";
                        FillObjects(histkey+"electrons", electrons_loose, weight);
                        FillObjects(histkey+"jets_lepVeto", jets_lepVeto, weight);
                        FillHist(histkey+"MET", MET, weight, 200, 0., 200.);
                        FillHist(histkey+"MT", Mt, weight, 200, 0., 200.);
						if (ptCorr > 35.) {
							histkey = channel+"/"+ID+"/"+syst;
							FillHist(histkey, ptCorr, absEta, weight, 4, ptCorrBin, 3, absEtaBin);
						}
					}
                }
            }		
		}
		if (channel == "DblEle") {
			if (syst.Contains("JetPtCut") || syst == "HasBjet") return;
            const Electron& ele1 = electrons_loose.at(0);
            const Electron& ele2 = electrons_loose.at(1);
            const Particle ZCand = ele1 + ele2;
            if (fabs(ZCand.M() - 91.2) > 15.) return;
            if (ele1.Charge() + ele2.Charge() != 0) return;

            TString trig = "";
            TString ID = "";
            TString histkey = "";
            if (passEle8Path) {
                trig = "HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30";
                ID = "PassLoose";

                // set weight
                double weight = 1.;
                if (!IsData) {
                    double genWeight = ev.MCweight()*weight_norm_1invpb;
                    double L1PrefireWeight = GetPrefireWeight(0);
                    double pileupWeight = GetPileUpWeight(nPileUp, 0);
                    weight = genWeight*L1PrefireWeight*pileupWeight*3.973;
                }

                histkey = channel+"/"+trig+"/"+ID+"/"+syst+"/";
                FillObjects(histkey+"electrons", electrons_loose, weight);
                FillObjects(histkey+"jets_lepVeto", jets_lepVeto, weight);
                FillObject(histkey+"ZCand", ZCand, weight);

                if (ele1.PassID("HcToWATight") && ele2.PassID("HcToWATight")) {
                    ID = "PassTight";
                    histkey = channel+"/"+trig+"/"+ID+"/"+syst+"/";
                    FillObjects(histkey+"electrons", electrons_loose, weight);
                    FillObjects(histkey+"jets_lepVeto", jets_lepVeto, weight);
                    FillObject(histkey+"ZCand", ZCand, weight);
                }
            }
			if (passEle12Path) {
                trig = "HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30";
                ID = "PassLoose";

                // set weight
                double weight = 1.;
                if (!IsData) {
                    double genWeight = ev.MCweight()*weight_norm_1invpb;
                    double L1PrefireWeight = GetPrefireWeight(0);
                    double pileupWeight = GetPileUpWeight(nPileUp, 0);
                    weight = genWeight*L1PrefireWeight*pileupWeight*27.699;
                }

                histkey = channel+"/"+trig+"/"+ID+"/"+syst+"/";
                FillObjects(histkey+"electrons", electrons_loose, weight);
                FillObjects(histkey+"jets_lepVeto", jets_lepVeto, weight);
                FillObject(histkey+"ZCand", ZCand, weight);

                if (ele1.PassID("HcToWATight") && ele2.PassID("HcToWATight")) {
                    ID = "PassTight";
                    histkey = channel+"/"+trig+"/"+ID+"/"+syst+"/";
                    FillObjects(histkey+"electrons", electrons_loose, weight);
                    FillObjects(histkey+"jets_lepVeto", jets_lepVeto, weight);
                    FillObject(histkey+"ZCand", ZCand, weight);
                }
            }
			if (passEle23Path) {
                trig = "HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30";
                ID = "PassLoose";

                // set weight
                double weight = 1.;
                if (!IsData) {
                    double genWeight = ev.MCweight()*weight_norm_1invpb;
                    double L1PrefireWeight = GetPrefireWeight(0);
                    double pileupWeight = GetPileUpWeight(nPileUp, 0);
                    weight = genWeight*L1PrefireWeight*pileupWeight*43.468;
                }

                histkey = channel+"/"+trig+"/"+ID+"/"+syst+"/";
                FillObjects(histkey+"electrons", electrons_loose, weight);
                FillObjects(histkey+"jets_lepVeto", jets_lepVeto, weight);
                FillObject(histkey+"ZCand", ZCand, weight);

                if (ele1.PassID("HcToWATight") && ele2.PassID("HcToWATight")) {
                    ID = "PassTight";
                    histkey = channel+"/"+trig+"/"+ID+"/"+syst+"/";
                    FillObjects(histkey+"electrons", electrons_loose, weight);
                    FillObjects(histkey+"jets_lepVeto", jets_lepVeto, weight);
                    FillObject(histkey+"ZCand", ZCand, weight);
                }
            }	
		}
	}	
}
