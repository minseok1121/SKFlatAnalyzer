#include "TestAnalyzer.h"

TestAnalyzer::TestAnalyzer(): __called(0) {}
TestAnalyzer::~TestAnalyzer() {}

void TestAnalyzer::initializeAnalyzer(){
	// flags
	RunSyst = HasFlag("RunSyst");

	// triggers
	// Double muon triggers for mumumu channel
	// TODO: Set a trigger for emumu channel
	// TODO: Set trigger safe pt cut
	if (DataYear == 2016) {
		DblMuTriggers = {
			"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v",			// BCDEF
			"HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v",			// BCDEF
			"HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v",		// BCDEF
			"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v",		// GH
			"HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v",		// GH
			"HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v"		// GH
		};
		EMuTriggers = {
		};
	}
	else if (DataYear == 2017) {
		DblMuTriggers = {
			"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DA_Mass3p8_v",	// whole year
			"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v"		// run >= 299368
		};
		EMuTriggers = {
		};
	}
	else if (DataYear == 2018) {
		DblMuTriggers = {
			"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v"		// whole year
		};
		EMuTriggers = {
		};
	}
	else {
		cerr << "[TestAnalyzer::initializeAnalyzer] DataYear is wrong" << endl;
		exit(EXIT_FAILURE);
	}
	TriggerSafePtCuts = {20, 10};

	// ID settings
	//MuonIDs = {"POGMedium", "POGLoose"}		// Medium for the siganl, Loose for veto
	//Electron = {}
	
	// B-tagging
	// Use 'mujets' in ttbar topology since it's statistically indepedent
	vector<JetTagging::Parameters> jtps;
	jtps.emplace_back(
		JetTagging::Parameters(JetTagging::DeepCSV, JetTagging::Medium, JetTagging::incl, JetTagging::mujets));
	mcCorr->SetJetTaggingParameters(jtps);
}

void TestAnalyzer::executeEvent(){
	// initiate cutflow
	MyCutflowMaker();
	
	// MET Filter
	if (!PassMETFilter()) return;
	MyCutflowMaker();

	Event ev = GetEvent();
	Particle METv = ev.GetMETVector();

	// object definition
	// No ID criteria, classify leptons using GetLeptonType_Public()
	vector<Gen> gens = GetGens();
	vector<Muon> muons = GetAllMuons();
	vector<Electron> electrons = GetAllElectrons();
	vector<Jet> jets = GetAllJets();

	// pt sort
	sort(muons.begin(), muons.end(), PtComparing);
	sort(electrons.begin(), electrons.end(), PtComparing);
	sort(jets.begin(), jets.end(), PtComparing);
	
	// Categorize muons
	vector<Muon> muons_prompt, muons_fake;
	vector<Electron> electrons_prompt, electrons_fake;

	for (const auto &ele: electrons) {
		for (unsigned int idx = 2; idx < gens.size(); idx++) {
			auto this_gen = gens.at(idx);
			if (this_gen.Status() != 1) continue;
			if (abs(this_gen.PID()) != 11) continue;
			if (this_gen.DeltaR(ele) > 0.1) continue;
			int LepType = GetLeptonType_Public(idx, gens);
			if (LepType > 0) electrons_prompt.emplace_back(ele);
			if (LepType < 0) electrons_fake.emplace_back(ele);
		}
	}

	for (const auto &mu: muons) {
		for (unsigned int idx = 2; idx < gens.size(); idx++) {
			auto this_gen = gens.at(idx);
			if (this_gen.Status() != 1) continue;
			if (abs(this_gen.PID()) != 13) continue;
			if (this_gen.DeltaR(mu) > 0.1) continue;
			int LepType = GetLeptonType_Public(idx, gens);
			if (LepType > 0) muons_prompt.emplace_back(mu);
			if (LepType < 0) muons_fake.emplace_back(mu);
		}
	}
	
	// set weight
	double weight = 1.;
	double w_prefire, w_gen, w_lumi, w_pileup;
	w_prefire = GetPrefireWeight(0);
	w_gen = ev.MCweight()*weight_norm_1invpb;
	w_lumi = ev.GetTriggerLumi("Full");
	w_pileup = GetPileUpWeight(nPileUp, 0);
	weight *= w_prefire*w_gen*w_lumi*w_pileup;

	MyHistoMaker("muons_prompt/", muons_prompt, weight);
	MyHistoMaker("muons_fake/", muons_fake, weight);
	MyHistoMaker("electrons_prompt/", electrons_prompt, weight);
	MyHistoMaker("electorns_fake/", electrons_fake, weight);
	MyHistoMaker("METv/", METv, weight);
}

// Other functions
void TestAnalyzer::MyCutflowMaker() {
	FillHist("cutflow", __called, 1., 20, 0., 20.);
	__called++;
}
void TestAnalyzer::MyHistoMaker(TString path, const vector<Muon> &muons, const double &weight) {
	TString obj_path;
	FillHist(path+"size", muons.size(), weight, 14, 0., 14.);
	for (unsigned int i = 0; i < muons.size(); i++) {
		obj_path = path + TString::Itoa(i+1, 10) + "/";
		FillHist(obj_path+"pt", muons.at(i).Pt(), weight, 300, 0., 300.);
	    FillHist(obj_path + "eta", muons.at(i).Eta(), weight, 50, -2.5, 2.5);
		FillHist(obj_path + "phi", muons.at(i).Phi(), weight, 70, -3.5, 3.5);
		FillHist(obj_path + "RelIso", muons.at(i).RelIso(), weight, 50, 0., 0.5);
		FillHist(obj_path + "TrkIso", muons.at(i).TrkIso()/muons.at(i).Pt(), weight, 80, 0., 0.8);
		FillHist(obj_path + "dXY", fabs(muons.at(i).dXY()), weight, 50, 0., 0.5);
		FillHist(obj_path + "dZ", fabs(muons.at(i).dZ()), weight, 80, 0., 0.8);
		FillHist(obj_path + "MVA", muons.at(i).MVA(), weight, 100, -1., 1);
	}
}
void TestAnalyzer::MyHistoMaker(TString path, const vector<Electron>& electrons, const double& weight) {
	TString obj_path;
	FillHist(path+"size", electrons.size(), weight, 14, 0., 14.);
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
void TestAnalyzer::MyHistoMaker(TString path, const vector<Jet>& jets, const double& weight) {
	TString obj_path;
	FillHist(path + "size", jets.size(), weight, 14, 0., 14.);
	for (unsigned int i = 0; i < jets.size(); i++) {
		obj_path = path + TString::Itoa(i+1, 10) + "/";
		FillHist(obj_path + "pt", jets.at(i).Pt(), weight, 300, 0., 300.);
		FillHist(obj_path + "eta", jets.at(i).Eta(), weight, 48, -2.4, 2.4);
		FillHist(obj_path + "phi", jets.at(i).Phi(), weight, 70, -3.5, 3.5);
	}	
}

void TestAnalyzer::MyHistoMaker(TString path, const Particle& METv, const double& weight) {
	FillHist(path + "pt", METv.Pt(), weight, 300, 0., 300.);
	FillHist(path + "eta", METv.Eta(), weight, 48, -2.4, 2.4); // of course, 0.
	FillHist(path + "phi", METv.Phi(), weight, 70, -3.5, 3.5);
}



































