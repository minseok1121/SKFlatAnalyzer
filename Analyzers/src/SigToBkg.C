#include "SigToBkg.h"

SigToBkg::SigToBkg() {}
SigToBkg::~SigToBkg() {}
void SigToBkg::initializeAnalyzer(){
	// This is for 2017
	if (DataYear == 2017) {
		trigs_dimu.emplace_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v");
		trigs_dimu.emplace_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v");
		// This trigger might not exist
		trigs_dimu.emplace_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v");
		
		trigs_emu.emplace_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v");
		trigs_emu.emplace_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");
	}	

	// IDs
	HcToWA_MuID = {"HcToWATight", "HcToWALoose"};
	HcToWA_EleID = {"HcToWATight", "HcToWALoose"};

	// B-tagging
	vector<JetTagging::Parameters> jtps;
	jtps.emplace_back(
		JetTagging::Parameters(JetTagging::DeepCSV, JetTagging::Medium, 
			JetTagging::incl, JetTagging::mujets));
	mcCorr->SetJetTaggingParameters(jtps);

	// initiate regions and cutflows
	regions = {"3mu", "1e2mu"};
	cuts = {"nocut", "metfilter", "nleptons", "trigger", "passSafePtCut", "exist_osmu", "Nj_ge1", "Nj_ge2", "Nb_ge1"};
	InitiateCutflow("3mu", cuts);
	InitiateCutflow("1e2mu", cuts);
}

void SigToBkg::executeEvent(){
	//FillCutflow
	for (const auto& region: regions)
		FillCutflow(region, "nocut");

	if (!PassMETFilter()) return;
	for (const auto& region: regions)
		FillCutflow(region, "metfilter");
	
	Event ev = GetEvent();
	vector<Gen> gens = GetGens();
	vector<Muon> muons= GetAllMuons();
	vector<Electron> electrons = GetAllElectrons();
	vector<Jet> jets = GetAllJets();
	Particle METv = ev.GetMETVector();
	
	// sort at the first time, I don't want to get confused
	sort(muons.begin(), muons.end(), PtComparing);
	sort(electrons.begin(), electrons.end(), PtComparing);
	sort(jets.begin(), jets.end(), PtComparing);

	// select objects
	vector<Muon> muons_tight = SelectMuons(muons, "HcToWATight", 10., 2.4);
	vector<Muon> muons_loose = SelectMuons(muons, "HcToWALoose", 10., 2.4);
	vector<Electron> electrons_tight = SelectElectrons(electrons, "HcToWATight", 10., 2.5);
	vector<Electron> electrons_loose = SelectElectrons(electrons, "HcToWALoose", 10., 2.5);
	vector<Jet> jets_tight = SelectJets(jets, "tight", 20., 2.7);
	vector<Jet> jets_cleaned = JetsVetoLeptonInside(jets_tight, electrons_loose, muons_loose, 0.4);
	vector<Jet> bjets_cleaned;
	for (const auto& jet: jets_cleaned) {
		if (jet.Eta() > 2.4)
			continue;
		double this_discr = jet.GetTaggerResult(JetTagging::DeepCSV);
		if (this_discr > mcCorr->GetJetTaggingCutValue(JetTagging::DeepJet, JetTagging::Medium))
			bjets_cleaned.emplace_back(jet);
	}

	// NOTE: Cutflow will automatically generated inside the SignalSelector
	TString signal = SignalRegionSelector(
			ev, muons_tight, electrons_tight, 
			muons_loose, electrons_loose, 
			jets_cleaned, bjets_cleaned);
	
	if (signal=="") return;

	// Let's categorize leptons
	vector<Muon> muons_signal, muons_ewprompt, muons_fromOffW, muons_fromtau, muons_fake;
	vector<Electron> electrons_signal, electrons_ewprompt, electrons_fromOffW, electrons_fromtau, electrons_fake;

	for (const auto& mu: muons_tight) {
		const int LepType = GetLeptonType(mu, gens);
		if (LepType == 1)
			muons_ewprompt.emplace_back(mu);
		else if (LepType == 2)
			muons_signal.emplace_back(mu);
		else if (LepType == 6)
			muons_fromOffW.emplace_back(mu);
		else if (LepType == 3)
			muons_fromtau.emplace_back(mu);
		else if (LepType < 0)
			muons_fake.emplace_back(mu);
		else
			continue;
	}
	for (const auto& ele: electrons_tight) {
		const int LepType = GetLeptonType(ele, gens);
		if (LepType == 1)
			electrons_ewprompt.emplace_back(ele);
		else if (LepType == 2)
			electrons_signal.emplace_back(ele); // should be empty
		else if (LepType == 6)
			electrons_fromOffW.emplace_back(ele);
		else if (LepType == 3)
			electrons_fromtau.emplace_back(ele);
		else if (LepType < 0)
			electrons_fake.emplace_back(ele);
		else
			continue;
	}

	// set weight
	double weight = 1.;
	double w_prefire, w_gen, w_lumi;
	w_prefire = GetPrefireWeight(0);
	w_gen = ev.MCweight() * weight_norm_1invpb;
	w_lumi = ev.GetTriggerLumi("Full");
	weight *= w_prefire*w_gen*w_lumi;

	// Finally, fill histograms
	// TODO: ACand?
	// Let's define Event level variables
	// Theses variables might not use in Preselection, but can be used afterward!
	double HT = 0.;
	for (const auto& jet: jets_cleaned)
		HT += jet.Pt();
	double ST = HT;
	for (const auto& muon: muons_tight)
		ST += muon.Pt();
	for (const auto& ele: electrons_tight)
		ST += ele.Pt();
	const double HToverST = HT/ST;

	vector<Lepton> leptons;
	for (const auto& mu: muons_tight)
		leptons.emplace_back(mu);
	for (const auto& ele: electrons_tight) 
		leptons.emplace_back(ele);
	sort(leptons.begin(), leptons.end(), PtComparing);
	//cout << "size l: " << leptons.size() << endl;
	//cout << "size j: " << jets_cleaned.size() << endl;
	//cout << "size b: " << bjets_cleaned.size() << endl;
	const double dRjj = jets_cleaned.at(0).DeltaR(jets_cleaned.at(1));
	const double dRbl1 = bjets_cleaned.at(0).DeltaR(leptons.at(0));
	const double dRbl2 = bjets_cleaned.at(0).DeltaR(leptons.at(1));
	const double dRbl3 = bjets_cleaned.at(0).DeltaR(leptons.at(2));
	const double dRl1l2 = leptons.at(0).DeltaR(leptons.at(1));
	const double dRl2l3 = leptons.at(1).DeltaR(leptons.at(2));
	const double dRl1l3 = leptons.at(0).DeltaR(leptons.at(3));
	
	const double ZMass = 91.2;
	Particle ZCand;
	if (signal == "1e2mu")
		ZCand = muons_tight.at(0) + muons_tight.at(1);
	if (signal == "3mu") {
		// select ZCand which mass is closest to Z
		Particle ZCand1, ZCand2;
		vector<Muon> muons_plus, muons_minus;
		for (const auto& mu: muons_tight) {
			if (mu.Charge() > 0)
				muons_plus.emplace_back(mu);
			else
				muons_minus.emplace_back(mu);
		}
		if (muons_plus.size() == 2) {
			ZCand1 = muons_plus.at(0) + muons_minus.at(0);
			ZCand2 = muons_plus.at(1) + muons_minus.at(0);
		}
		else if (muons_minus.size() == 2) {
			ZCand1 = muons_plus.at(0) + muons_minus.at(0);
			ZCand2 = muons_plus.at(0) + muons_minus.at(1);
		}
		else {
			cerr << "[SigToBkg::ExecuteEvent] wrong charge balance" << endl;
			cerr << "[SigToBkg::ExecuteEvent] muons plus: " << muons_plus.size() << endl;
			cerr << "[SigToBkg::ExecuteEvent] muons minus: " << muons_minus.size() << endl;
			exit(EXIT_FAILURE);
		}
		ZCand = fabs(ZCand1.M() - ZMass) < fabs(ZCand2.M() - ZMass) ? ZCand1 : ZCand2;
	}

	FillObjects(signal + "/muons_tight", muons_tight, weight);
	FillObjects(signal + "/muons_ewprompt", muons_ewprompt, weight);
	FillObjects(signal + "/muons_signal", muons_signal, weight);
	FillObjects(signal + "/muons_fromOffW", muons_fromOffW, weight);
	FillObjects(signal + "/muons_fromtau", muons_fromtau, weight);
	FillObjects(signal + "/muons_fake", muons_fake, weight);
	FillObjects(signal + "/electrons_tight", electrons_tight, weight);
	FillObjects(signal + "/electrons_ewprompt", electrons_ewprompt, weight);
	FillObjects(signal + "/electrons_signal", electrons_signal, weight);
	FillObjects(signal + "/electrons_fromOffW", electrons_fromOffW, weight);
	FillObjects(signal + "/electrons_fromtau", electrons_fromtau, weight);
	FillObjects(signal + "/electrons_fake", electrons_fake, weight);
	FillObjects(signal + "/jets_cleaned", jets_cleaned, weight);
	FillObjects(signal + "/bjets_cleaned", bjets_cleaned, weight);
	FillObject(signal + "/METv", METv, weight);
	FillObject(signal + "/ZCand", ZCand, weight);
	FillHist(signal + "/Event/HT", HT, weight, 300, 0., 300.);
	FillHist(signal + "/Event/ST", ST, weight, 300, 0., 300.);
	FillHist(signal + "/Event/HToverST", HToverST, weight, 20, 0., 1.);
	FillHist(signal + "/Event/dRjj", dRjj, weight, 40, 0., 4.);
	FillHist(signal + "/Event/dRbl1", dRbl1, weight, 40, 0., 4.);
	FillHist(signal + "/Event/dRbl2", dRbl2, weight, 40, 0., 4.);
	FillHist(signal + "/Event/dRbl3", dRbl3, weight, 40, 0., 4.);
	FillHist(signal + "/Event/dRl1l2", dRl1l2, weight, 40, 0., 4.);
	FillHist(signal + "/Event/dRl2l3", dRl2l3, weight, 40, 0., 4.);
	FillHist(signal + "/Event/dRl1l3", dRl1l3, weight, 40, 0., 4.);
}

//==== This is actual selector!
//==== TODO: might divide in more functions?
TString SigToBkg::SignalRegionSelector(
		Event &ev,
		vector<Muon> &muons_tight, vector<Electron> &electrons_tight,
		vector<Muon> &muons_loose, vector<Electron> &electrons_loose,
		vector<Jet> &jets, vector<Jet> &bjets) {
	// let's count leptons first
	TString signal = "";
	if (muons_tight.size() == 3 && muons_loose.size() == 3) 
		signal = "3mu";
	if (muons_tight.size() == 2 && muons_loose.size() == 2 && electrons_tight.size() == 1 && electrons_loose.size() == 1)
		signal = "1e2mu";
	
	if (signal != "") {
		for (const auto& region: regions)
			FillCutflow(region, "nleptons");
	}

	if (signal == "3mu") {
		if (! ev.PassTrigger(trigs_dimu)) return "";
		FillCutflow(signal, "trigger");

		if (muons_tight.at(0).Pt() < 20.)
			return "";
		if (muons_tight.at(1).Pt() < 10.)
			return "";
		if (muons_tight.at(2).Pt() < 10.)
			return "";
		FillCutflow(signal, "passSafePtCut");
	  
		// jets
		if (jets.size() < 1) return "";
		FillCutflow(signal, "Nj_ge1");

		if (jets.size() < 2) return "";
		FillCutflow(signal, "Nj_ge2");

		if (bjets.size() < 1) return "";
		FillCutflow(signal, "Nb_ge1");

		// leptons
	    int chargeSum
			= muons_tight.at(0).Charge() + muons_tight.at(1).Charge() + muons_tight.at(2).Charge();
        if (abs(chargeSum) != 1)
            return "";
        FillCutflow(signal, "exist_osmu");

		/*
        Particle ACand1 = muons_tight.at(0) + muons_tight.at(1);
        Particle ACand2 = muons_tight.at(1) + muons_tight.at(2);
        Particle ACand3 = muons_tight.at(2) + muons_tight.at(0);
        if (fabs(ACand1.M() - ZMass) < 10.) return "";
        if (fabs(ACand2.M() - ZMass) < 10.) return "";
        if (fabs(ACand3.M() - ZMass) < 10.) return "";
        FillCutflow(signal, "off_Zmass");
		*/
		
		return signal;
	}


	else if (signal == "1e2mu") {
		if (! ev.PassTrigger(trigs_emu)) return "";
		FillCutflow(signal, "trigger");

		bool passSafeCut = false;
		if (muons_tight.at(0).Pt() > 10. && electrons_tight.at(0).Pt() > 25.)
			passSafeCut = true;
		if (muons_tight.at(0).Pt() > 25. && electrons_tight.at(0).Pt() > 15.)
			passSafeCut = true;
		if (!passSafeCut) return "";
		FillCutflow(signal, "passSafePtCut");

		// jets
		if (jets.size() < 1)
			return "";
		FillCutflow(signal, "Nj_ge1");

		if (jets.size() < 2) 
			return "";
		FillCutflow(signal, "Nj_ge2");
		
		if (bjets.size() < 1)
			return "";
		FillCutflow(signal, "Nb_ge1");

		// leptons
		int chargeSum = muons_tight.at(0).Charge() + muons_tight.at(1).Charge();
        if (chargeSum != 0)
            return "";
        FillCutflow(signal, "exist_osmu");
		/*

        Particle ACand = muons_tight.at(0) + muons_tight.at(1);
        if (fabs(ACand.M() - ZMass) < 10.)
            return "";
        FillCutflow(signal, "off_Zmass");
		*/

		return signal;
	}
	else
		return "";
}

