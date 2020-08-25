#include "HNSignal_OS.h"


// To do::
// 1. trigger SF for POGTight - 2016
// 2. need check of trigger SF
// 3. complete signal selection

void HNSignal_OS::initializeAnalyzer(){
	// flags
	RunPOGTight = HasFlag("RunPOGTight");
	RunHighPt = HasFlag("RunHighPt");

	// get files
	// for POGTight
	f_trig_sf_lead_2016_BtoF = new TFile("/home/choij/SKFlat/data/Run2Legacy_v4/2016/ID/Muon/MuonTriggerEfficiency_ISR_DoubleIsoMu17Mu8_IsoMu17leg_RunBCDEF.root");
    f_trig_sf_lead_2016_GtoH = new TFile("/home/choij/SKFlat/data/Run2Legacy_v4/2016/ID/Muon/MuonTriggerEfficiency_ISR_DoubleIsoMu17Mu8_IsoMu17leg_RunGH.root");
	f_trig_sf_tail_2016_BtoF = new TFile("/home/choij/SKFlat/data/Run2Legacy_v4/2016/ID/Muon/MuonTriggerEfficiency_ISR_Mu8_OR_TkMu8_RunBCDEF.root");
	f_trig_sf_tail_2016_GtoH = new TFile("/home/choij/SKFlat/data/Run2Legacy_v4/2016/ID/Muon/MuonTriggerEfficiency_ISR_Mu8_OR_TkMu8_RunGH.root");
	f_trig_sf_lead_2017 = new TFile("/home/choij/SKFlat/data/Run2Legacy_v4/2017/ID/Muon/Mu_TrigEff_IsoMu17Mu8_Iso17_Run2017BCDEF_v1.root");
	f_trig_sf_tail_2017 = new TFile("/home/choij/SKFlat/data/Run2Legacy_v4/2017/ID/Muon/Mu_TrigEff_IsoMu17Mu8_Mu8_OR_TkMu8_Run2017BCDEF_v1.root");
	f_trig_sf_lead_2018 = new TFile("/home/choij/SKFlat/data/Run2Legacy_v4/2018/ID/Muon/2018_Mu17syst.root");
    f_trig_sf_tail_2018 = new TFile("/home/choij/SKFlat/data/Run2Legacy_v4/2018/ID/Muon/2018_Mu8syst.root");

	// for HighPt
	f_trig_sf_highpt_2016 = new TFile("/home/choij/SKFlat/data/Run2Legacy_v4/2016/ID/Muon/HighPtMuonTriggerSF_Run2016.root");
	f_trig_sf_highpt_2017 = new TFile("/home/choij/SKFlat/data/Run2Legacy_v4/2017/ID/Muon/HighPtMuonTriggerSF_Run2017.root");
	f_trig_sf_highpt_2018 = new TFile("/home/choij/SKFlat/data/Run2Legacy_v4/2018/ID/Muon/HighPtMuonTriggerSF_Run2018.root");
    f_fake_2016 = new TFile("/home/choij/SKFlat/FakeRate/Muon_HNtypeI_FakeRates_muonV1FR_2016.root");
    f_fake_2017 = new TFile("/home/choij/SKFlat/FakeRate/Muon_HNtypeI_FakeRates_muonV1FR_2017.root");
    f_fake_2018 = new TFile("/home/choij/SKFlat/FakeRate/Muon_HNtypeI_FakeRates_muonV1FR_2018.root");

	// triggers
	// only need doublemuon triggers
	if (DataYear == 2016) {
		trigs_POGTight = {
			"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v",		// BCDEF
      		"HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v",		// BCDEF
      		"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v",	// GH
      		"HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v"	// GH
		};
		trigs_HighPt = {
			"HLT_Mu50_v",
      		"HLT_TkMu50_v"
		};
	}
	else if (DataYear == 2017) {
		trigs_POGTight = {
			"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v"
    	};
		trigs_HighPt = {
			"HLT_Mu50_v",
      		"HLT_OldMu100_v",
      		"HLT_TkMu100_v"
    	};
	}
	else if (DataYear == 2018) {
		trigs_POGTight = {
			"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v"
		};
		trigs_HighPt = {
			"HLT_Mu50_v",
      		"HLT_OldMu100_v",
      		"HLT_TkMu100_v"
		};
	}
	else {
		cerr << "[HNSignal_OS::initializeAnalyzer] DataYear is wrong : " << DataYear << endl;
		exit(EXIT_FAILURE);
	}
	
	// ID settings
	HNType1_POGTight = {"HNType1_POGTight", "HNType1_POGLoose", "HNType1_POGVeto"};
	HNType1_HighPt = {"HNType1_HighPtTight", "HNType1_HighPtLoose", "HNType1_HighPtVeto"};
	HNType1_Electron = {"HNType1_CutBasedVeto"};
	HNType1_Jet = {"tight"};
	HNType1_FatJet = {"tightWithSDMass"};
	HNType1_IDSets = {"POGTight", "HighPt"};
	
	// B-tagging
	vector<JetTagging::Parameters> jtps;
	jtps.push_back(
			JetTagging::Parameters(JetTagging::DeepCSV, JetTagging::Medium, JetTagging::incl, JetTagging::comb));
	mcCorr->SetJetTaggingParameters(jtps);
}

void HNSignal_OS::executeEvent(){

	FillHist("preselection/cutflow", 0., 1., 10, 0., 10.); 
	// MET filter && trigger
	if (!PassMETFilter()) return;
	FillHist("preselection/cutflow", 1., 1., 10, 0., 10.);

	Event ev = GetEvent();
	METv = ev.GetMETVector();
	
	bool passTrigs;
	// trigger setting
	if (RunPOGTight) passTrigs = ev.PassTrigger(trigs_POGTight);
	else if (RunHighPt) passTrigs = ev.PassTrigger(trigs_HighPt);
	else {
		cerr << "[HNSignal_OS::executeEvent] Wrong flag" << endl;
		exit(EXIT_FAILURE);
	}
	
	// event should pass such trigger
	if (!passTrigs) return;
	FillHist("preselection/cutflow", 2., 1., 10, 0., 10.);

	// objects definition
	ClearCollections();
	gens = GetGens();
	muons = GetAllMuons();
	electrons = GetAllElectrons();
	fatjets = GetAllFatJets();
	jets = GetAllJets();

	sort(muons.begin(), muons.end(), PtComparing);
	sort(electrons.begin(), electrons.end(), PtComparing);
	sort(fatjets.begin(), fatjets.end(), PtComparing);
	sort(jets.begin(), jets.end(), PtComparing);

	muons = MuonPromptOnlyHNtype1(muons); // gen matching applied for MC
	if (RunPOGTight) {
		muons_veto = SelectMuons(muons, "HNType1_POGVeto", 5., 2.4);
		muons_loose = SelectMuons(muons_veto, "HNType1_POGLoose", 10., 2.4);
		muons_tight = SelectMuons(muons_loose, "HNType1_POGTight", 10, 2.4);
	}
	if (RunHighPt) {
		muons_veto = SelectMuons(muons, "HNType1_HighPtVeto", 5., 2.4);
		muons_loose = SelectMuons(muons_veto, "HNType1_HighPtLoose", 10., 2.4);
		muons_tight = SelectMuons(muons_loose, "HNType1_HighPtTight", 10., 2.4);
	}

	electrons_veto = SelectElectrons(electrons, "HNType1_CutBasedVeto", 10, 2.5);

	fatjets_tight = SelectFatJets(fatjets, "tightWithSDMass", 200., 2.7);
	fatjets_dR10 = FatJetsVetoLeptonInside(fatjets_tight, electrons_veto, muons_veto, 1.0);

	jets_tight = SelectJets(jets, "tight", 20, 2.7);
	jets_dR04 = JetsVetoLeptonInside(jets_tight, electrons_veto, muons_veto, 0.4);
	jets_awayFatJet = JetsAwayFromFatJet(jets_dR04, fatjets_dR10, 0.8);

	// can't use jets outside |eta| > 2.4 for b-taggings
	// no lepton veto for bjets
	bjets = SelectJets(jets_tight, "tight", 20., 2.4);
	for (const auto& j: bjets) {
		double this_discr = j.GetTaggerResult(JetTagging::DeepCSV);
		if (this_discr >  mcCorr->GetJetTaggingCutValue(JetTagging::DeepCSV, JetTagging::Medium)) 
			bjets_tight.push_back(j);
	}
	//bjets_dR04 = JetsVetoLeptonInside(bjets_tight, electrons_veto, muons_veto, 0.4);
	
	// preselection
	// 1. pass trigger. Done above.
	// 2. 2 os tight leptons. no 3rd lepton
	// 3. truth matching(?)
	// 4. muons should above pt threshold
	// 5. m(ll) > 10 GeV
	// 6. at least one jet or fatjet
	
	// for fake estimation, loose id first
	if (! (muons_loose.size() == 2)) return;
	if (! (muons_loose.at(0).Charge()*muons_loose.at(1).Charge() < 0)) return;
	if (! (muons_veto.size() == 2)) return;
	if (! (electrons_veto.size() == 0)) return;
	bool tightFlag = (muons_tight.size() == 2);
	//if (muons_tight.size() !=2) cout << tightFlag << endl; // should be false
	if (!IsDATA && !tightFlag) return; // for MC, no need for loose id
	if (tightFlag) FillHist("preselection/cutflow", 3., 1., 10, 0., 10.);

	if (RunPOGTight) {
		if (! (muons_loose.at(0).Pt() > 20.)) return;
		if (! (muons_loose.at(1).Pt() > 10.)) return;
	}
	if (RunHighPt) {
		if (! (muons_loose.at(0).Pt() > 50.)) return;
		if (! (muons_loose.at(1).Pt() > 50.)) return;
	}
	if (tightFlag) FillHist("preselection/cutflow", 4., 1., 10, 0., 10.);

	Particle reso = muons_loose.at(0) + muons_loose.at(1);
	if (! (reso.M() > 10)) return;
	if (tightFlag) FillHist("preselection/cutflow", 5., 1., 10, 0., 10.);
	
	if (jets_awayFatJet.size()==0 && fatjets_dR10.size()==0) return;
	if (tightFlag) FillHist("preselection/cutflow", 6., 1., 10, 0., 10.);

	// weight settings
	ClearWeights();
	double weight = 1.;
	if (!IsDATA) {
		w_prefire = GetPrefireWeight(0);
		w_gen = ev.MCweight()*weight_norm_1invpb;
		w_lumi = ev.GetTriggerLumi("Full");
		w_pileup = GetPileUpWeight(nPileUp, 0);

		TString id_key, iso_key;
		if (RunPOGTight) {
			id_key = "NUM_TightID_DEN_genTracks";
			iso_key = "NUM_TightRelIso_DEN_TightIDandIPCut";
			for (const auto& mu: muons_loose) {
				sf_muid *= mcCorr->MuonID_SF(id_key, mu.Eta(), mu.MiniAODPt());
				sf_muiso *= mcCorr->MuonISO_SF(iso_key, mu.Eta(), mu.MiniAODPt());
			}
			//sf_trig *= mcCorr->MuonTrigger_SF($ID, $ISO, muons_tight, 0);
		}
		else if (RunHighPt) {
			id_key = "NUM_HighPtID_DEN_genTracks";
			iso_key = "NUM_LooseRelTkIso_DEN_HighPtIDandIPCut";
			for (const auto& mu:muons_loose) {
				sf_muid *= mcCorr->MuonID_SF(id_key, mu.Eta(), mu.MiniAODPt());
				sf_muiso *= mcCorr->MuonISO_SF(iso_key, mu.Eta(), mu.MiniAODPt());
			}
		}
		sf_trig = GetTrigSF(muons_loose);
		
		JetTagging::Parameters jtp_DeepCSV_Medium = JetTagging::Parameters(JetTagging::DeepCSV, JetTagging::Medium, JetTagging::incl, JetTagging::comb);
		// no lepton veto for bjets in this analysis
		// bjets == jets_tight with |eta| < 2.4
		// sf_btag = mcCorr->GetBTaggingReweight_1a(bjets, jtp_DeepCSV_Medium); 

		//if (w_prefire == 1.) cerr << "Warning: w_prefire = 1." << endl;
		if (w_gen == 1.) cerr << "Warning: w_gen = 1."  << endl;
		if (w_lumi == 1.) cerr << "Warning:w_lumi = 1." << endl;
		if (w_pileup == 1.) cerr << "Warning: w_pileup = 1." << endl;
		if (sf_muid == 1. || sf_muid > 10. || sf_muid < 0.5) cerr << "Warning: sf_muid = " << sf_muid << endl;
		if (sf_muiso == 1. || sf_muiso > 10. || sf_muiso < 0.5) cerr << "Warning: sf_muiso = " << sf_muiso << endl;
		if (sf_trig == 1. || sf_trig > 10. || sf_trig < 0.5) cerr << "Warning:: sf_trig = " << sf_trig << endl;;
		// if (bjets.size()!=0 && sf_btag == 1.) cerr << "Warning::sf_btag = 1.";
	}
	weight *= w_prefire*w_gen*w_lumi*w_pileup;
	weight *= sf_muid*sf_muiso*sf_trig*sf_btag;

	// measure fake contribution from data
	// ID setting for fake contribution
    TString TightID, LooseID;
    if (RunPOGTight) {
        TightID = "HNType1_POGTight";
        LooseID = "HNType1_POGLoose";
    }
    else if (RunHighPt) {
        TightID = "HNType1_HighPtTight";
        LooseID = "HNType1_HighPtLoose";
    }	
	// apply weight as w_fake
	double w_fake = weight;
	if (IsDATA) {
	    if (tightFlag) w_fake *= 0.;
		else w_fake *= GetFakeWeight(TightID);
	}
	// now draw histograms for each region
	// IMPORTANT:: set path for fake contributions
	TString path;
	if (tightFlag) {
		// mc and data with tightFlag
		path = "preselection/";
		DrawHists(path + "muons_tight/", muons_loose, weight);
		DrawHists(path + "jets_awayFatJet/", jets_awayFatJet, weight);
		DrawHists(path + "bjets_tight/", bjets_tight, weight);
		DrawHists(path + "fatjets_dR10/", fatjets_dR10, weight);
		DrawHists(path + "METv/", METv, weight);
	}

	if (IsDATA) {
		// fake contribution
		path = "preselection/fake/";
		// loose muons? Get fake fraction tha pass tight ID from loose ones
		DrawHists(path + "muons_tight/", muons_loose, w_fake);
    	DrawHists(path + "jets_awayFatJet/", jets_awayFatJet, w_fake);
    	DrawHists(path + "bjets_tight/", bjets_tight, w_fake);
    	DrawHists(path + "fatjets_dR10/", fatjets_dR10, w_fake);
    	DrawHists(path + "METv/", METv, w_fake);
	}
}

// ==== other functions ====
void HNSignal_OS::ClearCollections() {
	muons.clear(); muons_tight.clear(); muons_loose.clear(), muons_veto.clear();
	electrons.clear(); electrons_veto.clear();
	jets.clear(); jets_tight.clear(); jets_dR04.clear(); jets_awayFatJet.clear();
	bjets.clear(); bjets_tight.clear(); bjets_dR04.clear();
	fatjets.clear(); fatjets_tight.clear(); fatjets_dR10.clear();
	gens.clear();
}

void HNSignal_OS::ClearWeights() {
	w_prefire = 1.; w_gen = 1.; w_lumi = 1.;  w_pileup = 1.;
	sf_muid = 1.; sf_muiso = 1.; sf_trig = 1.; sf_btag = 1.;
}

bool HNSignal_OS::isSignal(const SR& sig) const {
	switch (sig) {
	case SR1:
		if (! (jets_dR04.size() >= 2)) return false;
		if (! (fatjets_dR10.size() == 0)) return false;
		return true;
		break;
	case SR2:
		if (! (fatjets_dR10.size() > 0)) return false;
		return true;
		break;
	default:
		cerr << "[HNSignal_OS::isSignal] Wrong Signal region... sig = " << sig << endl;
		exit(EXIT_FAILURE);
	}
}

bool HNSignal_OS::isControl(const CR& con) const {
	cout << "[HNSignal_OS::isControl] Control region has not been set yet" << endl;
	exit(EXIT_FAILURE);

	switch (con) {
	case CR1:
		if (! (jets_dR04.size() >= 2)) return false;
		if (! (fatjets_dR10.size() == 0)) return false;
		return true;
		break;
	case CR2:
		if (! (fatjets_dR10.size() > 0)) return false;
		return true;
		break;
	default:
		cerr << "[HNSignal_OS::isControl] Wrong control region... con = " << con << endl;
		exit(EXIT_FAILURE);
	}
}


double HNSignal_OS::GetFakeWeight(const TString& TightID) {
	// get file to use
	TFile* f = nullptr;
	double TightIso = -999.;
	if (RunPOGTight) {
		TightIso = 0.15;
		if (DataYear==2016) f = f_fake_2016;
		else if (DataYear==2017) f = f_fake_2017;
		else if (DataYear==2018) f = f_fake_2018;
		else {
			cerr << "[HNSignal_OS::GetFakeWeight] Wrong year... DataYear = " << DataYear << endl;
			exit(EXIT_FAILURE);
		}
		if (!f) {
			cerr << "[HNSignal_OS::GetFakeWeight] null TFile" << endl;
			exit(EXIT_FAILURE);
		}
	}
	if (RunHighPt) {
		TightIso = 0.1;
		cout << "[HNSignal_OS::GetFakeWeight] TFile for high pt is not set yet" << endl;
		exit(EXIT_FAILURE);
	}

	TH2D* h = (TH2D*)f->Get("AwayJetPt40");

	double w_fake = 1.;
	for (const auto& mu: muons_loose) {
		if (mu.PassID(TightID)) continue;
		else {
			double pt = mu.Pt();
			double eta = fabs(mu.Eta());
			double RelIso = mu.RelIso();
			
			// pt cone
			double pt_cone = pt;
			if (RelIso > TightIso) pt_cone *= (1 + RelIso - TightIso);
			
			// get fake rate
			if (eta > 2.4) eta = 2.39;
			if (pt_cone > 60) pt_cone = 59.;

			int this_bin = h->FindBin(pt_cone, eta);
			double this_fake = h->GetBinContent(this_bin);

			// evaluate fake weight
			w_fake *= (-1. * this_fake / (1. - this_fake));
		}
	}
	w_fake *= -1.;
	
	delete h;
	return w_fake;
}

double HNSignal_OS::GetTrigSF(const vector<Muon>& muons) {
	// can use for factorized one only
	double sf = -999.;
	// for 2016
	double sf_BtoF = -999.;
	double sf_GtoH = -999.; 
	if (RunPOGTight) {
		TH2D* h_lead = nullptr;
    	TH2D* h_tail = nullptr;
		TH2D* h_lead_BtoF = nullptr;
		TH2D* h_tail_BtoF = nullptr;
		TH2D* h_lead_GtoH = nullptr;
		TH2D* h_tail_GtoH = nullptr;
		if (DataYear==2016) {
			h_lead_BtoF = (TH2D*)f_trig_sf_lead_2016_BtoF->Get("TriggerEff_ISR_BCDEF_SF");
			h_tail_BtoF = (TH2D*)f_trig_sf_tail_2016_BtoF->Get("TriggerEff_ISR_BCDEF_SF");
			h_lead_GtoH = (TH2D*)f_trig_sf_lead_2016_GtoH->Get("TriggerEff_ISR_GH_SF");
			h_tail_GtoH = (TH2D*)f_trig_sf_tail_2016_GtoH->Get("TriggerEff_ISR_GH_SF");
		}
		else if (DataYear==2017) {
			h_lead = (TH2D*)f_trig_sf_lead_2017->Get("TriggerEff_AFB_2017BCDEF_SF");
			h_tail = (TH2D*)f_trig_sf_tail_2017->Get("TriggerEff_AFB_2017BCDEF_SF");
		}
		else if (DataYear==2018) {
			h_lead = (TH2D*)f_trig_sf_lead_2018->Get("SF_abseta_pt");
			h_tail = (TH2D*)f_trig_sf_tail_2018->Get("SF_abseta_pt");
		}
		else {
			cerr << "[HNSignal_OS::GetTrigSF] Wrong year" << endl;
			exit(EXIT_FAILURE);
		}

		Muon lead = muons.at(0); Muon tail = muons.at(1);
    	double pt_lead = lead.Pt(), eta_lead = fabs(lead.Eta());
    	double pt_tail = tail.Pt(), eta_tail = fabs(tail.Eta());
    	if (pt_lead >= 120.) pt_lead = 119.;
    	if (pt_tail >= 120.) pt_tail = 119.;
    	if (eta_lead >= 2.4) eta_lead = 2.39;
    	if (eta_tail >= 2.4) eta_tail = 2.39;

		if (DataYear==2016) {
			if (!h_lead_BtoF || !h_lead_GtoH || !h_tail_BtoF || !h_tail_GtoH) {
				cerr << "[HNSignal_OS::GetTrigSF] null hist" << endl;
				exit(EXIT_FAILURE);
			}

			int bin_lead_BtoF = h_lead_BtoF->FindBin(eta_lead, pt_lead);
			int bin_lead_GtoH = h_lead_GtoH->FindBin(eta_lead, pt_lead);
			int bin_tail_BtoF = h_tail_BtoF->FindBin(eta_tail, pt_tail);
			int bin_tail_GtoH = h_tail_GtoH->FindBin(eta_tail, pt_tail);
			sf_BtoF = h_lead_BtoF->GetBinContent(bin_lead_BtoF) * h_tail_BtoF->GetBinContent(bin_tail_BtoF);
			sf_GtoH = h_lead_GtoH->GetBinContent(bin_lead_GtoH) * h_tail_GtoH->GetBinContent(bin_tail_GtoH);
			// using average weight
			double lumi_BtoF = 27.3, lumi_BtoH = 35.9;
			double lumi_GtoH = lumi_BtoH - lumi_BtoF;
			double w_BtoF = lumi_BtoF / lumi_BtoH;
			double w_GtoH = lumi_GtoH / lumi_BtoH;
			sf =  w_BtoF * sf_BtoF + w_GtoH * sf_GtoH;

			if (sf < 0.8 || sf > 1.2) {
				cerr << "weired sf" << endl;
				cerr << "sf_BtoF = " << sf_BtoF << endl;
				cerr << "sf_GtoH = " << sf_GtoH << endl;
			}
			delete h_lead_BtoF;
			delete h_tail_BtoF;
			delete h_lead_GtoH;
			delete h_tail_GtoH;
		}
		else {
			if (!h_lead || !h_tail) {
				cerr << "[HNSignal_OS::GetTrigSF] null hist" << endl;
				exit(EXIT_FAILURE);
			}

			int bin_lead = h_lead->FindBin(eta_lead, pt_lead);
			int bin_tail = h_tail->FindBin(eta_tail, pt_tail);
			sf = h_lead->GetBinContent(bin_lead) * h_tail->GetBinContent(bin_tail);

			if (sf < 0.8 || sf > 1.2) {
				cerr << "WARNING: sf too small" << endl;
				cerr << "sf = " << sf << endl;
				cerr << "pt_lead = " << pt_lead << endl;
				cerr << "eta_lead = " << eta_lead << endl;
				cerr << "sf_lead = " << h_lead->GetBinContent(bin_lead) << endl;
				cerr << "pt_tail = " << pt_tail << endl;
				cerr << "eta_tail = " << eta_tail << endl;
				cerr << "sf_tail = " << h_tail->GetBinContent(bin_tail) << endl << endl;
			}
			delete h_lead;
			delete h_tail;
		}
	}
	else if (RunHighPt) {
		// need remodeling
		TH2D* h = nullptr;
		if (DataYear==2016) {
			h = (TH2D*)f_trig_sf_highpt_2016->Get("SF");
		}
		else if (DataYear==2017) {
			h = (TH2D*)f_trig_sf_highpt_2017->Get("SF");
		}
		else if (DataYear==2018) {
			h = (TH2D*)f_trig_sf_highpt_2018->Get("SF");
		}

		Muon mu = muons.at(0);
		double pt = mu.Pt(); double eta = fabs(mu.Eta());
		if (pt > 1000.) pt = 999.;
		if (eta > 2.4) eta = 2.39;

		int bin = h->FindBin(pt, eta);
		sf = h->GetBinContent(bin);
	}

	return sf;
}

// need to set path before executing functions
void HNSignal_OS::DrawHists(TString path, const vector<Muon>& muons, const double& weight) {
	TString obj_path;
	for (unsigned int i = 0; i < muons.size(); i++) {
		obj_path = path + TString::Itoa(i+1, 10) + "/";
		FillHist(obj_path + "pt", muons.at(i).Pt(), weight, 300, 0., 300.);
		FillHist(obj_path + "eta", muons.at(i).Eta(), weight, 48, -2.4, 2.4);
		FillHist(obj_path + "phi", muons.at(i).Phi(), weight, 70, -3.5, 3.5);
		FillHist(obj_path + "RelIso", muons.at(i).RelIso(), weight, 50, 0., 0.5);
		FillHist(obj_path + "TrkIso", muons.at(i).TrkIso(), weight, 80, 0., 0.8);
	}
}

void HNSignal_OS::DrawHists(TString path, const vector<Electron>& electrons, const double& weight) {
	TString obj_path;
	for (unsigned int i = 0; i < electrons.size(); i++) {
		obj_path = path + TString::Itoa(i+1, 10) + "/";
		FillHist(obj_path + "pt", electrons.at(i).Pt(), weight, 300, 0., 300.);
		FillHist(obj_path + "eta", electrons.at(i).Eta(), weight, 50, -2.5, 2.5);
		FillHist(obj_path + "phi", electrons.at(i).Phi(), weight, 70, -3.5, 3.5);
		FillHist(obj_path + "RelIso", electrons.at(i).RelIso(), weight, 50, 0., 0.5);
		FillHist(obj_path + "dXY", fabs(electrons.at(i).dXY()), weight, 50, 0., 0.5);
		FillHist(obj_path + "dZ", fabs(electrons.at(i).dZ()), weight, 80, 0., 0.8);
	}
}
void HNSignal_OS::DrawHists(TString path, const vector<Jet>& jets, const double& weight) {
	TString obj_path;
	FillHist(path + "size", jets.size(), weight, 14, 0., 14.);
	for (unsigned int i = 0; i < jets.size(); i++) {
		obj_path = path + TString::Itoa(i+1, 10) + "/";
		FillHist(obj_path + "pt", jets.at(i).Pt(), weight, 300, 0., 300.);
		FillHist(obj_path + "eta", jets.at(i).Eta(), weight, 48, -2.4, 2.4);
		FillHist(obj_path + "phi", jets.at(i).Phi(), weight, 70, -3.5, 3.5);
	}	
}
void HNSignal_OS::DrawHists(TString path, const vector<FatJet>& jets, const double& weight) {
	TString obj_path;
	FillHist(path + "size", jets.size(), weight, 10, 0., 10.);
	for (unsigned int i = 0; i < jets.size(); i++) {
		obj_path = path + TString::Itoa(i+1, 10) + "/";
		FillHist(obj_path + "pt", jets.at(i).Pt(), weight, 500, 0., 500.);
		FillHist(obj_path + "eta", jets.at(i).Eta(), weight, 48, -2.4, 2.4);
		FillHist(obj_path + "phi", jets.at(i).Phi(), weight, 70, -3.5, 3.5);
		FillHist(obj_path + "SDMass", jets.at(i).SDMass(), weight, 300, 0., 300.);
		if (jets.at(i).PuppiTau1() != 0.) {
			double tau21 = jets.at(i).PuppiTau2() / jets.at(i).PuppiTau1();
			FillHist(obj_path + "tau21", tau21, weight, 100, 0., 1.);
		}
	}
}
void HNSignal_OS::DrawHists(TString path, const Particle& METv, const double& weight) {
	FillHist(path + "pt", METv.Pt(), weight, 300, 0., 300.);
	FillHist(path + "eta", METv.Eta(), weight, 48, -2.4, 2.4); // of course, 0.
	FillHist(path + "phi", METv.Phi(), weight, 70, -3.5, 3.5);
}

vector<Muon> HNSignal_OS::MuonPromptOnlyHNtype1(const vector<Muon>& muons) {
	if (IsDATA) return muons;

	vector<Muon> out;
	for (const auto& mu : muons) {
		if (-5 < GetLeptonType(mu, gens) && GetLeptonType(mu, gens) <= 0) continue;
		out.push_back(mu);
	}

	return out;
}

HNSignal_OS::HNSignal_OS(){}

HNSignal_OS::~HNSignal_OS(){}


