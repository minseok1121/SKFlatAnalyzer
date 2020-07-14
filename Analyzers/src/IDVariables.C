#include "IDVariables.h"

//==== This analyzer is made to ====
//==== see the ID Variable distributions ====
//==== for electrons and muons ====

void IDVariables::initializeAnalyzer(){

	//=== ID settings ====
	MuonIDs = {"POGLoose", "POGMedium", "POGTight"};

	//==== Btagging ====
	vector<JetTagging::Parameters> jtps;
	jtps.push_back( JetTagging::Parameters(JetTagging::DeepCSV, JetTagging::Medium, JetTagging::incl, JetTagging::comb));
	mcCorr->SetJetTaggingParameters(jtps);

	if (DataYear == 2016) {
		TrigNames.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v");
		TrigNames.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v");

		if (IsDATA && run > 280385) {
			TrigNames.clear();
			TrigNames.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v");
			TrigNames.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZv");
		}

		if (!IsDATA) {
			TrigNames.clear();
			TrigNames.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v");
			TrigNames.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v");
			TrigNames.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v");
			TrigNames.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZv");
		}
		TriggerSafePtCut1 = 20.;
		TriggerSafePtCut2 = 10.;
	}
	else if (DataYear == 2017) {
		TrigNames.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v");
		if (IsDATA && run >= 299368) 
			TrigNames.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v");
		if (!IsDATA) 
			TrigNames.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v");
		TriggerSafePtCut1 = 20.;
		TriggerSafePtCut2 = 10.;
	}
	else if (DataYear == 2018) {
		TrigNames.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v");
		TriggerSafePtCut1 = 20.;
		TriggerSafePtCut2 = 10.;
	}
	else {
		cout << "[IDVariables::initializeAnalyzer] Wrong Year" << endl;
		exit(EXIT_FAILURE);
	}

	//==== finish initialization ====
	cout << "[IDVariables::initializeAnalyzer] initialization finished" << endl;
}

void IDVariables::executeEvent(){

	AllMuons = GetAllMuons();
	AllJets = GetAllJets();
	AllElectrons = GetAllElectrons();
	AllGens = GetGens();

	AnalyzerParameter param;

	//==== Loopt over muonIDs ====
	for(const auto &mu : MuonIDs) {
		MuonVetoID = "POGLoose";
		ElecVetoID = "passLooseID";
		MuonID = mu;
		JetID = "tight";

		param.Clear();
		param.syst_ = AnalyzerParameter::Central;
		param.Name = MuonID + "_Central";

		executeEventFromParameter(param);
	}
}

void IDVariables::executeEventFromParameter(AnalyzerParameter param){

	//==== No cut ====
	TString path = param.Name + "/";
	FillHist(path + "PreSelection", 0., 1., 5, 0., 5.);

	//=== METFilter ====
	if(!PassMETFilter()) return;

	//==== Copy All Objects ====
	vector<Muon> this_AllMuons = AllMuons;
	vector<Electron> this_AllElectrons = AllElectrons;
	vector<Jet> this_AllJets = AllJets;

	Event ev = GetEvent();
    Particle METv = ev.GetMETVector();
	Particle METv_xyCorr(pfMET_Type1_PhiCor_pt*TMath::Cos(pfMET_Type1_PhiCor_phi), pfMET_Type1_PhiCor_pt*TMath::Sin(pfMET_Type1_PhiCor_phi), 0., pfMET_Type1_PhiCor_pt);

	vector<Muon> muons = SelectMuons(this_AllMuons, MuonID, 25., 2.4);
	vector<Muon> muons_veto = SelectMuons(this_AllMuons, MuonVetoID, 20., 2.4);
	vector<Electron> electrons_veto = SelectElectrons(this_AllElectrons, ElecVetoID, 20., 2.5);
	vector<Jet> jets = SelectJets(this_AllJets, JetID, 30., 2.4);
	vector<Jet> jets_dR04 = JetsVetoLeptonInside(jets, electrons_veto, muons_veto, 0.4);
	vector<Jet> bjets_dR04;

	JetTagging::Parameters jtp_DeepCSV_Medium 
		= JetTagging::Parameters(JetTagging::DeepCSV, JetTagging::Medium, JetTagging::incl, JetTagging::comb);


	for (const auto j : jets_dR04) {
		double this_discr = j.GetTaggerResult(JetTagging::DeepCSV);
		if (this_discr > mcCorr->GetJetTaggingCutValue(JetTagging::DeepCSV, JetTagging::Medium))
			bjets_dR04.push_back(j);
	}



	sort(muons.begin(), muons.end(), PtComparing);
	sort(muons_veto.begin(), muons_veto.end(), PtComparing);
	sort(electrons_veto.begin(), electrons_veto.end(), PtComparing);
	sort(jets.begin(), jets.end(), PtComparing);
	sort(jets_dR04.begin(), jets_dR04.end(), PtComparing);
	sort(bjets_dR04.begin(), bjets_dR04.end(), PtComparing);


	//==== Preselection ====
	if (! (muons.size() == 2)) return;
	if (! (muons_veto.size() == 2)) return;
	if (! (electrons_veto.size() == 0)) return;

	if (! ev.PassTrigger(TrigNames)) return;
	if (! (muons.at(0).Pt() > TriggerSafePtCut1)) return;
	if (! (muons.at(1).Pt() > TriggerSafePtCut2)) return;
	//cout << "path = " << path << endl;

	FillHist(path + "PreSelection", 1., 1., 5, 0., 5.);

	//==== event weight ====
	weight = 1.;
	w_gen = 1.;
	w_filter = 1.;
	w_toppt = 1.;
	w_lumi = 1.;
	w_pileup = 1.;
	w_prefire = 1.;
	sf_trig = 1.;
	sf_mutk = 1.;
	sf_muid = 1.;
	sf_muiso = 1.;
	sf_elreco = 1.;
	sf_elid = 1.;
	sf_btag = 1.;
	if (!IsDATA) {
		w_gen = ev.MCweight();
		w_lumi = weight_norm_1invpb*ev.GetTriggerLumi("Full");
		if (MCSample.Contains("TT") && MCSample.Contains("powheg"))
			w_toppt = mcCorr->GetTopPtReweight(AllGens);
		w_pileup = GetPileUpWeight(nPileUp, 0);
		w_prefire = GetPrefireWeight(0);
		sf_btag = mcCorr->GetBTaggingReweight_1a(jets_dR04, jtp_DeepCSV_Medium);
	}
	weight *= w_gen * w_filter * w_lumi * w_toppt * w_pileup * w_prefire;
	weight *= sf_trig * sf_mutk * sf_muid * sf_muiso * sf_elreco * sf_elid * sf_btag;

	//==== Event selection ====
	bool isDY = IsDY(muons, electrons_veto);
	bool isTTbar = IsTTbar(muons, electrons_veto, jets_dR04, bjets_dR04);

	if (isDY) {
		TString pathDY = path + "DY/";
		Particle ZCand = muons.at(0) + muons.at(1);
		double dRDY = muons.at(0).DeltaR(muons.at(1));
		
		FillHist(pathDY+"M(mumu)", ZCand.M(), weight, 60, 60, 120);
		FillHist(pathDY+"dR(mumu)", dRDY, weight, 80, 0., 2.);
		FillHist(pathDY+"nJets", jets_dR04.size(), weight, 10, 0., 10.);
		FillHist(pathDY+"nBJets", bjets_dR04.size(), weight, 10, 0., 10.);
		FillHist(pathDY+"METv", METv.Pt(), weight, 240, 0, 240);
		FillHist(pathDY+"METv_phi", METv.Phi(), weight, 70, -3.5, 3.5);
		FillHist(pathDY+"METv_xyCorr_pt", METv_xyCorr.Pt(), weight, 70, -3.5, 3.5);
		FillHist(pathDY+"METv_xyCorr_phi", METv_xyCorr.Phi(), weight, 70, -3.5, 3.5);
		DrawHists(pathDY, muons.at(0), 0, weight);
		DrawHists(pathDY, muons.at(1), 1, weight);
		if (jets_dR04.size() > 0) DrawHists(pathDY, jets_dR04.at(0), 0, weight);
		if (jets_dR04.size() > 1) DrawHists(pathDY, jets_dR04.at(1), 1, weight);
	}

	if (isTTbar) {
		TString pathTT = path + "TT/";
		double dRTT = muons.at(0).DeltaR(muons.at(1));

		FillHist(pathTT+"dR(mumu)", dRTT, weight, 80, 0., 2.);
		FillHist(pathTT+"nJets", jets_dR04.size(), weight, 10, 0., 10.);
		FillHist(pathTT+"nBJets", bjets_dR04.size(), weight, 10, 0., 10.);
		FillHist(pathTT+"METv", METv.Pt(), weight, 240, 0, 240);
		FillHist(pathTT+"METv_phi", METv.Phi(), weight, 70, -3.5, 3.5);
		FillHist(pathTT+"METv_xyCorr_pt", METv_xyCorr.Pt(), weight, 70, -3.5, 3.5);
		FillHist(pathTT+"METv_xyCorr_phi", METv_xyCorr.Phi(), weight, 70, -3.5, 3.5);
		DrawHists(pathTT, muons.at(0), 0, weight);
		DrawHists(pathTT, muons.at(1), 1, weight);
		DrawHists(pathTT, jets_dR04.at(0), 0, weight);
		DrawHists(pathTT, jets_dR04.at(1), 1, weight);
		DrawHists(pathTT, bjets_dR04.at(0), 0, weight);
		if (bjets_dR04.size() > 1) DrawHists(pathTT, bjets_dR04.at(1), 1, weight);
		
	}

	return;


}
//==== member functions ====
bool IDVariables::IsDY(const vector<Muon> &muons, const vector<Electron> &electrons_veto) {
	TString path = MuonID + "_Central/cutflow_DY";
	FillHist(path, 0., 1., 5, 0., 5.);	
	if (! (muons.at(0).Charge() * muons.at(1).Charge() < 0)) return false;
	FillHist(path, 1., 1., 5, 0., 5.);
	Particle ZCand = muons.at(0) + muons.at(1);
	if (! (IsOnZ(ZCand.M(), 10))) return false;
	FillHist(path, 2., 1., 5, 0., 5.);

	return true;
}

bool IDVariables::IsTTbar(const vector<Muon> &muons, const vector<Electron> &electrons_veto, const vector<Jet> &jets, const vector<Jet> &bjets) {
	TString path = MuonID + "_Central/cutflow_TT";
	FillHist(path, 0., 1., 5, 0., 5.);
	if (! (jets.size() >= 2)) return false;
	FillHist(path, 1., 1., 5, 0., 5.);
	if (! (bjets.size() >= 1)) return false;
	FillHist(path, 2., 1., 5, 0., 5.);
	if (! (muons.at(0).Charge() * muons.at(1).Charge() < 0)) return false;
	FillHist(path, 3., 1., 5, 0., 5.);
	if (! (muons.at(0).DeltaR(muons.at(1)) > 0.4)) return false;
	FillHist(path, 4., 1., 5, 0., 5.);

	return true;
}

void IDVariables::DrawHists(TString path, const Muon &mu, unsigned int order, const double weight) {

	TString this_path = path + "mu" + TString::Itoa(order+1, 10) + "/";
	FillHist(this_path + "pt", mu.Pt(), weight, 240, 0, 240);
	FillHist(this_path + "eta", mu.Eta(), weight, 50, -2.5, 2.5);
	FillHist(this_path + "phi", mu.Phi(), weight, 70, -3.5, 3.5);
	
	this_path += "IDvar/";
	FillHist(this_path + "RelIso", mu.RelIso(), weight, 100, 0, 1.);
    FillHist(this_path + "dXY" , mu.dXY(), weight, 100, -0.1, 0.1);
    FillHist(this_path + "dXYerr", mu.dXYerr(), weight, 100, 0., 0.1);
    if (mu.dXYerr() != 0) {
        FillHist(this_path + "SIP2D", fabs(mu.dXY())/mu.dXYerr(), weight, 100, 0., 10.);
    }
    FillHist(this_path + "dZ", mu.dZ(), weight, 100, -0.1, 0.1);
    FillHist(this_path + "Chi2", mu.Chi2(), weight, 100, 0., 10.);
}

void IDVariables::DrawHists(TString path, const Jet &j, unsigned int order, const double weight) {

	TString this_path = path + "j" + TString::Itoa(order+1, 10) + "/";
	FillHist(this_path + "pt", j.Pt(), weight, 240, 0, 240);
	FillHist(this_path + "eta", j.Eta(), weight, 50, -2.5, 2.5);
	FillHist(this_path + "phi", j.Phi(), weight, 70, -3.5, 3.5);
}

IDVariables::IDVariables(){

}

IDVariables::~IDVariables(){

}


