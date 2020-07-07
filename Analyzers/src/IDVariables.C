#include "IDVariables.h"

//==== This analyzer is made to ====
//==== see the ID Variable distributions ====
//==== for electrons and muons ====

void IDVariables::initializeAnalyzer(){

	//=== ID settings ====
	ElectronIDs = {"passVetoID", "passLooseID", "passMediumID", "passTightID"};
	ElectronMVAIDs = {"passMVAID_noIso_WP90", "passMVAID_noIso_WP80", "passMVAID_iso_WP90", "passMVAID_iso_WP80"};
	MuonIDs = {"POGLoose", "POGMedium", "POGTight"};

	//==== Btagging ====
	//vector<JetTagging::Parameters> jtps;
	//jtps.push_back( JetTagging::Parameters(JetTagging::DeepCSV, JetTagging::Medium, JetTagging::incl, JetTagging::comb));
	//mcCorr->SetJetTaggingParameters(jtps);

	//==== finish initialization ====
	cout << "[IDVariables::initializeAnalyzer] initialization finished" << endl;
}

void IDVariables::executeEvent(){

	AllMuons = GetAllMuons();
	AllJets = GetAllJets();
	AllElectrons = GetAllElectrons();
	AllGens = GetGens();

	//==== for reweight ====
	weight_Prefire = GetPrefireWeight(0);
	weight_PileUp = GetPileUpWeight(nPileUp, 0);
	weight_TopPt = mcCorr->GetTopPtReweight(AllGens);

	AnalyzerParameter param;

	//==== Loop over muonIDs and ElectronIDs ====
	for (unsigned int i = 0; i < MuonIDs.size(); i++) {
		for (unsigned int j = 0; j < ElectronIDs.size(); j++) {
			MuonVetoID = "POGLoose";
			ElectronVetoID = "passVetoID";
			MuonID = MuonIDs.at(i);
			ElectronID = ElectronIDs.at(j);
			JetID = "tight";


			param.Clear();
			param.syst_ = AnalyzerParameter::Central;
			param.Name = MuonID + "_" + ElectronID + "_Central";

			executeEventFromParameter(param);
		}

		for (unsigned int j = 0; j < ElectronMVAIDs.size(); j++) {
            MuonID = MuonIDs.at(i);
            ElectronID = ElectronMVAIDs.at(j);
			JetID = "tight";

            param.Clear();
            param.syst_ = AnalyzerParameter::Central;
            param.Name = MuonID + "_" + ElectronID + "_Central";

            executeEventFromParameter(param);
        }	
	}
}

void IDVariables::executeEventFromParameter(AnalyzerParameter param){

	//==== No cut ====
	TString path = param.Name + "/";
	FillHist(path + "NoCut", 0., 1., 1, 0., 1.);

	//=== METFilter ====
	if(!PassMETFilter()) return;

	Event ev = GetEvent();
	Particle METv = ev.GetMETVector(); // No usage for MET currently

	//==== Copy All Objects ====
	vector<Muon> this_AllMuons = AllMuons;
	vector<Electron> this_AllElectrons = AllElectrons;
	vector<Jet> this_AllJets = AllJets;

	vector<Muon> muons = SelectMuons(this_AllMuons, MuonID, 25., 2.4);
	vector<Muon> muons_veto = SelectMuons(this_AllMuons, MuonVetoID, 20., 2.4);
	vector<Electron> electrons = SelectElectrons(this_AllElectrons, ElectronID, 25., 2.5);
	vector<Electron> electrons_veto = SelectElectrons(this_AllElectrons, ElectronVetoID, 10., 2.5);
	vector<Jet> jets = SelectJets(this_AllJets, JetID, 30., 2.4);

	sort(muons.begin(), muons.end(), PtComparing);
	sort(muons_veto.begin(), muons_veto.end(), PtComparing);
	sort(electrons.begin(), electrons.end(), PtComparing);
	sort(electrons_veto.begin(), electrons_veto.end(), PtComparing);
	sort(jets.begin(), jets.end(), PtComparing);

	vector<Jet> jets_dR04 = JetsVetoLeptonInside(jets, electrons_veto, muons_veto, 0.4);

	//==== Btagging ====
	//JetTagging::Parameters jtp_DeepCSV_Medium = JetTagging::Parameters(JetTagging::DeepCSV, JetTagging::Medium, JetTagging::incl, JetTagging::comb);
	unsigned int Nb = 0;
	for (unsigned int i = 0; i < jets_dR04.size(); i++) {
		double this_discr = jets_dR04.at(i).GetTaggerResult(JetTagging::DeepCSV);
		if (this_discr > mcCorr->GetJetTaggingCutValue(JetTagging::DeepCSV, JetTagging::Medium)) Nb++;
	}
	//double btagWeight_DeepCSV = mcCorr->GetBTaggingReweight_1a(jets, jtp_DeepCSV_Medium);

	//==== Event selection ====
	//==== SRDY1
	//==== DrellYan ee channel
	//==== use SingleElectron dataset
	//==== 1. Exactly 2 SFOS electrons, no additional muons
	//==== 2. |M(ll) - 91.2| < 15GeV
	path = param.Name + "/SRDY1/";
	TString TriggerSRDY1 = "HLT_Ele27_WPTight_Gsf_v";
	double TriggerSafeCutSRDY1 = 30.;
	bool SRDY1 = true;


	if (! (ev.PassTrigger(TriggerSRDY1))) SRDY1 = false;
	if (! (electrons.size() == 2)) SRDY1 = false;
	if (! (electrons_veto.size() ==2)) SRDY1 = false;
	if (! (muons_veto.size() == 0)) SRDY1 = false;
	if (SRDY1) {
		if (! (electrons.at(0).Pt() > TriggerSafeCutSRDY1)) SRDY1 = false;
		if (! (electrons.at(0).Charge()*electrons.at(1).Charge() < 0)) SRDY1 = false;
		Particle ZCand1 = electrons.at(0) + electrons.at(1);
		if (! (IsOnZ(ZCand1.M(), 15.))) SRDY1 = false;
	}

	//==== Draw histograms ====
	double weight_SRDY1 = 1.;
	if (SRDY1) {
		FillHist(path + "LeadingElectronPt", electrons.at(0).Pt(), weight_SRDY1, 200, 0., 200.);
		FillHist(path + "LeadingElectronEta", electrons.at(0).Eta(), weight_SRDY1, 50, -2.5, 2.5);
		FillHist(path + "SubLeadingElectronPt", electrons.at(1).Pt(), weight_SRDY1, 200, 0., 200.);
		FillHist(path + "SubLeadingElectronEta", electrons.at(1).Eta(), weight_SRDY1, 50, -2.5, 2.5);

		DrawIDVariables(path, electrons.at(0), 0, weight_SRDY1);
		DrawIDVariables(path, electrons.at(1), 1, weight_SRDY1);
	}

	//==== SRDY2
	//==== DrellYan mumu channel
	//==== use SingleMuon dataset
	//==== 1. Exactly 2 SFOS muons, no additional electrons
	//==== 2. |M(ll) - 91.2| < 15 GeV
	path = param.Name + "/SRDY2/";
	TString TriggerSRDY2 = "HLT_IsoMu24_v";
	double TriggerSafeCutSRDY2 = 27.;
	bool SRDY2 = true;

	if (! (ev.PassTrigger(TriggerSRDY2))) SRDY2 = false;
	if (! (muons.size() == 2)) SRDY2 = false;
	if (! (muons_veto.size() == 2)) SRDY2 = false;
	if (! (electrons_veto.size() == 0)) SRDY2 = false;
	if (SRDY2) {
		if (! (muons.at(0).Pt() > TriggerSafeCutSRDY2)) SRDY2 = false;
		if (! (muons.at(0).Charge()*muons.at(1).Charge() < 0)) SRDY2 = false;
		Particle ZCand2 = muons.at(0) + muons.at(1);
		if (! (IsOnZ(ZCand2.M(), 15.))) SRDY2 = false;
	}

	//==== Draw histograms ====
	double weight_SRDY2 = 1.;
	if (SRDY2) {
		FillHist(path + "LeadingMuonPt", muons.at(0).Pt(), weight_SRDY2, 200, 0., 200.);
        FillHist(path + "LeadingMuonEta", muons.at(0).Eta(), weight_SRDY2, 50, -2.5, 2.5);
        FillHist(path + "SubLeadingMuonPt", muons.at(1).Pt(), weight_SRDY2, 200, 0., 200.);
        FillHist(path + "SubLeadingMuonEta", muons.at(1).Eta(), weight_SRDY2, 50, -2.5, 2.5);

		DrawIDVariables(path, muons.at(0), 0, weight_SRDY2);
		DrawIDVariables(path, muons.at(1), 1, weight_SRDY2);
    }

	//==== SRTT1
	//==== tt to WWee channel
	//==== use SingleElectron dataset
	//==== 1. at least 2 jets, at least 1 bjets (for cleaned jets)
	//==== 2. exactly two OS electrons, no additional veto leptons
	//==== 3. dR(ll) > 0.4
	path = param.Name + "/SRTT1/";
	TString TriggerSRTT1 = "HLT_Ele27_WPTight_Gsf_v";
	double TriggerSafeCutSRTT1 = 30.;
	bool SRTT1 = true;

	if (! (jets_dR04.size() >= 2)) SRTT1 = false;
	if (! (Nb >= 1)) SRTT1 = false;
	if (! (ev.PassTrigger(TriggerSRTT1))) SRTT1 = false;
	if (! (electrons.size() == 2)) SRTT1 = false;
	if (! (electrons_veto.size() == 2)) SRTT1 = false;
	if (! (muons_veto.size() == 0)) SRTT1 = false;
	if (SRTT1) {
		if (! (electrons.at(0).Pt() > TriggerSafeCutSRTT1)) SRTT1 = false;
		if (! (electrons.at(0).Charge()*electrons.at(1).Charge() < 0)) SRTT1 = false;
		if (! (electrons.at(0).DeltaR(electrons.at(1)) > 0.4)) SRTT1 = false;
	}

	//==== Draw histograms ====
	double weight_SRTT1 = 1.;
	if (SRTT1) {
		FillHist(path + "LeadingElectronPt", electrons.at(0).Pt(), weight_SRTT1, 200, 0., 200.);
        FillHist(path + "LeadingElectronEta", electrons.at(0).Eta(), weight_SRTT1, 50, -2.5, 2.5);
        FillHist(path + "SubLeadingElectronPt", electrons.at(1).Pt(), weight_SRTT1, 200, 0., 200.);
        FillHist(path + "SubLeadingElectronEta", electrons.at(1).Eta(), weight_SRTT1, 50, -2.5, 2.5);
		FillHist(path + "Njets", jets_dR04.size(), weight_SRTT1, 6, 0, 6);
		FillHist(path + "NBjets", Nb, weight_SRTT1, 6, 0, 6);

		DrawIDVariables(path, electrons.at(0), 0, weight_SRTT1);
		DrawIDVariables(path, electrons.at(1), 1, weight_SRTT1);
	}

	//==== SRTT2
    //==== tt to WWmumu channel
    //==== use SingleMuon dataset
    //==== 1. at least 2 jets, at least 1 bjets (for cleaned jets)
    //==== 2. exactly two OS muons, no additional veto leptons
    //==== 3. dR(ll) > 0.4
    path = param.Name + "/SRTT2/";
    TString TriggerSRTT2 = "HLT_IsoMu24_v";
    double TriggerSafeCutSRTT2 = 27.;
    bool SRTT2 = true;

	if (! (jets_dR04.size() >= 2)) SRTT2 = false;
    if (! (Nb >= 1)) SRTT2 = false;
    if (! (ev.PassTrigger(TriggerSRTT2))) SRTT2 = false;
    if (! (muons.size() == 2)) SRTT2 = false;
    if (! (muons_veto.size() == 2)) SRTT2 = false;
    if (! (electrons_veto.size() == 0)) SRTT2 = false;
    if (SRTT2) {
        if (! (muons.at(0).Pt() > TriggerSafeCutSRTT2)) SRTT2 = false;
        if (! (muons.at(0).Charge()*muons.at(1).Charge() < 0)) SRTT2 = false;
        if (! (muons.at(0).DeltaR(muons.at(1)) > 0.4)) SRTT2 = false;
    }

	//==== Draw histograms ====
	double weight_SRTT2 = 1.;
    if (SRTT2) {
        FillHist(path + "LeadingMuonPt", muons.at(0).Pt(), weight_SRTT2, 200, 0., 200.);
        FillHist(path + "LeadingMuonEta", muons.at(0).Eta(), weight_SRTT2, 50, -2.5, 2.5);
        FillHist(path + "SubLeadingMuonPt", muons.at(1).Pt(), weight_SRTT2, 200, 0., 200.);
        FillHist(path + "SubLeadingMuonEta", muons.at(1).Eta(), weight_SRTT2, 50, -2.5, 2.5);   
        FillHist(path + "Njets", jets_dR04.size(), weight_SRTT2, 6, 0, 6);
        FillHist(path + "NBjets", Nb, weight_SRTT2, 6, 0, 6);

		DrawIDVariables(path, muons.at(0), 0, weight_SRTT2);
		DrawIDVariables(path, muons.at(1), 1, weight_SRTT2);
    }

	//==== SRTT3
    //==== tt to WWemu channel
    //==== use SingleMuon dataset
    //==== 1. at least 2 jets, at least 1 bjets (for cleaned jets)
    //==== 2. exactly one OS electron and muon pair, no additional veto leptons
    //==== 3. dR(ll) > 0.4
    path = param.Name + "/SRTT3/";
    TString TriggerSRTT3 = "HLT_IsoMu24_v";
    double TriggerSafeCutSRTT3 = 27.;
    bool SRTT3 = true;

    if (! (jets_dR04.size() >= 2)) SRTT3 = false;
    if (! (Nb >= 1)) SRTT3 = false;
    if (! (ev.PassTrigger(TriggerSRTT3))) SRTT3 = false;
	if (! (electrons.size() == 1)) SRTT3 = false;
	if (! (electrons_veto.size() == 1)) SRTT3 = false;
	if (! (muons.size() == 1)) SRTT3 = false;
	if (! (muons_veto.size() == 1)) SRTT3 = false;
    if (SRTT3) {
        if (! (muons.at(0).Pt() > TriggerSafeCutSRTT3)) SRTT3 = false;
        if (! (electrons.at(0).Charge()*muons.at(0).Charge() < 0)) SRTT3 = false;
        if (! (electrons.at(0).DeltaR(muons.at(0)) > 0.4)) SRTT3 = false;
    }
	
	//==== Draw histograms ====
	double weight_SRTT3 = 1.;	
	if (SRTT3) {
        FillHist(path + "MuonPt", muons.at(0).Pt(), weight_SRTT3, 200, 0., 200.);
        FillHist(path + "MuonEta", muons.at(0).Eta(), weight_SRTT3, 50, -2.5, 2.5);
        FillHist(path + "ElectronPt", electrons.at(0).Pt(), weight_SRTT3, 200, 0., 200.);
        FillHist(path + "ElectronEta", electrons.at(0).Eta(), weight_SRTT3, 50, -2.5, 2.5);   
        FillHist(path + "Njets", jets_dR04.size(), weight_SRTT3, 6, 0, 6); 
        FillHist(path + "NBjets", Nb, weight_SRTT3, 6, 0, 6);

		DrawIDVariables(path, electrons.at(0), 0, weight_SRTT3);
		DrawIDVariables(path, muons.at(0), 0, weight_SRTT3);
    }	
	

}
//==== member functions ====
void IDVariables::DrawIDVariables(TString path, const Electron &e, unsigned int order, const double weight) {
	// order guard
	if (order > 3) {
		cout << "[IDVariables::DrawIDVariables] order = " << order << endl;
		cout << "[IDVariables::DrawIDVariables] Wrong order" << endl;
		exit(EXIT_FAILURE);
	}
	
	TString this_path = path + "e" + TString::Itoa(order+1, 10) + "/";
	FillHist(this_path + "RelIso", e.RelIso(), weight, 100, 0., 1.);
	FillHist(this_path + "dXY", e.dXY(), weight, 100, -0.5, 0.5);
	FillHist(this_path + "dXYerr", e.dXYerr(), weight, 100, -0.5, 0.5);
	if (e.dXYerr() != 0) FillHist(this_path + "SIP2D", e.dXY() / e.dXYerr(), weight, 200, -10, 10);
	FillHist(this_path + "dZ", e.dZ(), weight, 100, 0., 1.);

	FillHist(this_path + "Full5x5_sigmaIetaIeta", e.Full5x5_sigmaIetaIeta(), weight, 100, 0., 1.);
	FillHist(this_path + "dEtaSeed", e.dEtaSeed(), weight, 100, 0., 1.);
	FillHist(this_path + "dPhiIn", e.dPhiIn(), weight, 100, 0., 1.);
	FillHist(this_path + "HoverE", e.HoverE(), weight, 100, 0., 1.);
	FillHist(this_path + "InvEminusInvP", e.InvEminusInvP(), weight, 100, 0., 1.);
	FillHist(this_path + "e2x5OverE5x5", e.e2x5OverE5x5(), weight, 100, 0., 1.);
	FillHist(this_path + "e1x5OverE5x5", e.e1x5OverE5x5(), weight, 100, 0., 1.);
	FillHist(this_path + "trackIso", e.TrkIso(), weight, 100, 0., 1.);
	FillHist(this_path + "dr03EcalRecHitSumEt", e.dr03EcalRecHitSumEt(), weight, 100, 0., 1.);
	FillHist(this_path + "dr03HcalDepth1TowerSumEt", e.dr03HcalDepth1TowerSumEt(), weight, 100, 0., 1.);
	FillHist(this_path + "dr03HcalTowerSumEt", e.dr03HcalTowerSumEt(), weight, 100, 0., 1.);
	FillHist(this_path + "dr03TkSumPt", e.dr03TkSumPt(), weight, 100, 0., 1.);
	FillHist(this_path + "ecalPFClusterIso", e.ecalPFClusterIso(), weight, 100, 0., 1.);
	FillHist(this_path + "hcalPFClusterIso", e.hcalPFClusterIso(), weight, 100, 0., 1.);

}

void IDVariables::DrawIDVariables(TString path, const Muon &mu, unsigned int order, const double weight) {
	// odrder gaurd
	if (order > 3) {
		cout << "[IDVariables::DrawIDVariables] order = " << order << endl;
		cout << "[IDVariables::DrawIDVariables] Wrong order" << endl;
		exit(EXIT_FAILURE);
	}

	TString this_path = path + "mu" + TString::Itoa(order+1, 10) + "/";
	FillHist(this_path + "RelIso", mu.RelIso(), weight, 100, 0, 1.);
	FillHist(this_path + "dXY" , mu.dXY(), weight, 100, -0.5, 0.5);
	FillHist(this_path + "dXYerr", mu.dXYerr(), weight, 100, -0.5, 0.5);
	if (mu.dXYerr() != 0) {
		FillHist(this_path + "SIP2D", mu.dXY()/mu.dXYerr(), weight, 200, -10, 10);
	}
	FillHist(this_path + "dZ", mu.dZ(), weight, 100, -0.5, 0.5);
	FillHist(this_path + "Chi2", mu.Chi2(), weight, 200, -10, 10);
}

IDVariables::IDVariables(){

}

IDVariables::~IDVariables(){

}


