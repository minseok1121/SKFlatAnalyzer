#include "MakeNPV.h"

/////////////////////////////////////////////////////////////////////////////
//==== MakeNPV                                                   ====//
//==== This analyzer is for evaluation of the fake rate of electorns   ====//
//==== from jet activity                                               ====// 
//==== using data                                                      ====//
//==== Contact: choij@cern.ch                                          ====//
//==== Based on SKFlatAnalyzer - Run2Legacy_v4                         ====// 
/////////////////////////////////////////////////////////////////////////////

// FakeEstiamtor: Central estimation of Fake Rate

void MakeNPV::initializeAnalyzer(){

	//==== ID setting for electrons
	ElectronIDs = {"passLooseID", "passTightID", "FakeLooseID", "FakeTightID"};
	
	//==== Trigger Setting
	if (DataYear == 2016) {
		HLTElecTriggerName = "HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v";
		TriggerSafePtCut = 25.;
	}
	else if (DataYear == 2017) {
		HLTElecTriggerName = "HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v";
		TriggerSafePtCut = 25.;
	}
	else if (DataYear == 2018) {
		cout << "[MakeNPV::initializeAnalyzer] Trigger is not set for 2018" << endl;
		exit(EXIT_FAILURE);
	}
	else {
		cout << "[MakeNPV::initializeAnalyzer] Wrong Year" << endl;
		exit(EXIT_FAILURE);
	}

	cout << "[MakeNPV::initializeAnalyzer] HLTElecTriggerName = " << HLTElecTriggerName << endl;
	cout << "[MakeNPV::initializeAnalyzer] TriggerSafePtCut = " << TriggerSafePtCut << endl;

	//==== B-Tagging
	std::vector<JetTagging::Parameters> jtps;
	jtps.push_back( JetTagging::Parameters(JetTagging::DeepCSV, JetTagging::Medium, JetTagging::incl, JetTagging::comb) );
	mcCorr->SetJetTaggingParameters(jtps);
	
	cout << "[FakeEstiamtor::initializeAnalyer] Finish initialization" << endl;
}

void MakeNPV::executeEvent(){

	//==== Copy all objects ====
	AllMuons = GetAllMuons();
	AllElectrons = GetAllElectrons();
	AllJets = GetAllJets();
	AllGens = GetGens();

	//==== Get L1Prefire. pileup reweight ====
	weight_Prefire = GetPrefireWeight(0);
	weight_PileUp = GetPileUpWeight(nPileUp, 0);

	AnalyzerParameter param;

	//==== Loop over electron IDs ====
	for (unsigned int i = 0; i < ElectronIDs.size(); i++) {
		MuonID = "POGLoose"; // to veto loose muons
		ElectronID = ElectronIDs.at(i);
		JetID = "tight";

		param.Clear();
		param.syst_ = AnalyzerParameter::Central;
		param.Name = ElectronID + "_Central";
		param.Jet_ID = JetID;
		

		executeEventFromParameter(param);
	}			
}

void MakeNPV::executeEventFromParameter(AnalyzerParameter param){

	if(!PassMETFilter()) return;

    Event ev = GetEvent();
    Particle METv = ev.GetMETVector();

    //==== Trigger ====
    if (! (ev.PassTrigger(HLTElecTriggerName) )) return;
	
	//==== Copy all objects ====
	vector<Muon> this_AllMuons = AllMuons;
	vector<Electron> this_AllElectrons = AllElectrons;
	vector<Jet> this_AllJets = AllJets;

	//==== ID Selection ====
	muons = SelectMuons(this_AllMuons, MuonID, 10., 2.4);
	electrons = SelectElectrons(this_AllElectrons, ElectronID, 25, 2.5);
	jets = SelectJets(this_AllJets, param.Jet_ID, 20., 2.4);

	if (ElectronID.Contains("pass")) electrons_loose = SelectElectrons(this_AllElectrons, "passLooseID", 25, 2.5);
	else if (ElectronID.Contains("Fake")) electrons_loose = SelectElectrons(this_AllElectrons, "FakeLooseID", 25, 2.5);
	
	std::sort(muons.begin(), muons.end(), PtComparing);
	std::sort(electrons.begin(), electrons.end(), PtComparing);
	std::sort(electrons_loose.begin(), electrons_loose.end(), PtComparing);
	std::sort(jets.begin(), jets.end(), PtComparing);

	int NBjets_NoSF = 0;
	int NBjets_WithSF_2a = 0;
	JetTagging::Parameters jtp_DeepCSV_Medium = JetTagging::Parameters(JetTagging::DeepCSV, JetTagging::Medium, JetTagging::incl, JetTagging::comb);
	double btagWeight = mcCorr->GetBTaggingReweight_1a(jets, jtp_DeepCSV_Medium);

	//==== leading electron should pass trigger-safe pt cut
	if (electrons.size() == 0) return;
	if (electrons.at(0).Pt() <= TriggerSafePtCut) return;

	///////////////////////////////////////////////////////////////////////
	//==== measure nPV reweight ====
	//==== 1. no loose muon
	//==== 2. exactly 1 electron
	//==== 3. jet - with systematic source, but use clean_jets04
	///////////////////////////////////////////////////////////////////////
	
	if (muons.size() != 0) return;
	double weight = 1.;

	if (electrons.size() == 1) {
		clean_jets04 = JetsVetoLeptonInside(jets, electrons_loose, muons, 0.4);

		if (clean_jets04.size() == 0) return;
		if (clean_jets04.at(0).Pt() < 30) return;

		//==== Get weights
		if (!IsDATA) {
			weight *= weight_norm_1invpb*ev.GetTriggerLumi(HLTElecTriggerName);
			weight *= ev.MCweight();
			weight *= weight_Prefire;
		}
		
		// Measure nPV reweights
		if (nPV > 100) nPV = 100;

		if (clean_jets04.at(0).Pt() > 30) {
			JSFillHist(param.Name, "nPV_" + param.Name + "_JetPtCut30", nPV, weight, 100, 0, 100);
		}
		if (clean_jets04.at(0).Pt() > 35) {
			JSFillHist(param.Name, "nPV_" + param.Name + "_JetPtCut35", nPV, weight, 100, 0, 100);
		}
		if (clean_jets04.at(0).Pt() > 40) {
			JSFillHist(param.Name, "nPV_" + param.Name + "_JetPtCut40", nPV, weight, 100, 0, 100);
			if (NBjets_WithSF_2a > 1) {
				JSFillHist(param.Name, "nPV_" + param.Name + "_BtagDep", nPV, weight, 100, 0, 100);
			}
		}
		if (clean_jets04.at(0).Pt() > 45) {
			JSFillHist(param.Name, "nPV_" + param.Name + "_JetPtCut45", nPV, weight, 100, 0, 100);
		}
		if (clean_jets04.at(0).Pt() > 50) {
			JSFillHist(param.Name, "nPV_" + param.Name + "_JetPtCut50", nPV, weight, 100, 0, 100);
		}
		if (clean_jets04.at(0).Pt() > 55) {
			JSFillHist(param.Name, "nPV_" + param.Name + "_JetPtCut55", nPV, weight, 100, 0, 100);
		}
		if (clean_jets04.at(0).Pt() > 60) {
			JSFillHist(param.Name, "nPV_" + param.Name + "_JetPtCut60", nPV, weight, 100, 0, 100);
		}
		return;
	}
	else return;
}

MakeNPV::MakeNPV(){

}

MakeNPV::~MakeNPV(){

}


