#include "FakeEstimator.h"

/////////////////////////////////////////////////////////////////////////////
//==== FakeEstimator                                                   ====//
//==== This analyzer is for evaluation of the fake rate of electorns   ====//
//==== from jet activity                                               ====// 
//==== using data                                                      ====//
//==== Contact: choij@cern.ch                                          ====//
//==== Based on SKFlatAnalyzer - Run2Legacy_v4                         ====// 
/////////////////////////////////////////////////////////////////////////////

void FakeEstimator::initializeAnalyzer(){

	//==== Flag for systematic sources
	RunSysts = HasFlag("RunSysts");
	RunNormSysts = HasFlag("RunNormSysts");
	RunXsecSyst = HasFlag("RunXsecSyst");

	Systs = {"Central", "FlavorDep", "JetPtCut"};

	cout << "[FakeEstimator::initializeAnalyzer] Systematic Sources" << endl;
	cout << "[FakeEstimator::initializeAnalyzer] RunSysts = " << RunSysts << endl;
	cout << "[FakeEstimator::initializeAnalyzer] RunNormSysts = " << RunNormSysts << endl;
	cout << "[FakeEstimator::initializeAnalyzer] RunXsecSyst = " << RunXsecSyst << endl;
	if (RunSysts == true) {
		cout << "[FakeEstimator::initializeAnalyzer] Systematic sources: ";
		for (unsigned int i = 0; i < Systs.size(); i++) {
			cout << Systs.at(i) << " ";
		}
		cout << endl;
	}

	//==== ID setting for electrons
	ElectronIDs = {"PassLooseID", "PassTightID", "FakeLooseID", "FakeTightID"};
	
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
		cout << "[FakeEstimator::initializeAnalyzer] Trigger is not set for 2018" << endl;
		exit(EXIT_FAILURE);
	}
	else {
		cout << "[FakeEstimator::initializeAnalyzer] Wrong Year" << endl;
		exit(EXIT_FAILURE);
	}

	cout << "[FakeEstimator::initializeAnalyzer] HLTElecTriggerName = " << HLTElecTriggerName << endl;
	cout << "[FakeEstimator::initializeAnalyzer] TriggerSafePtCut = " << TriggerSafePtCut << endl;

	//==== B-Tagging
	std::vector<JetTagging::Parameters> jtps;
	jtps.push_back( JetTagging::Parameters(JetTagging::DeepCSV, JetTagging::Medium, JetTagging::incl, JetTagging::comb) );
	mcCorr->SetJetTaggingParameters(jtps);
	
	cout << "[FakeEstiamtor::initializeAnalyer] Finish initialization" << endl;
}

void FakeEstimator::executeEvent(){

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

		//==== Systematic sources ====
		//==== 1. B-tagging
		//==== 2. Jet Cut Up & Down

		if (RunSysts) {
			for (unsigned int i = 1; i < Systs.size(); i++) {
				param.syst_ = AnalyzerParameter::Syst(i);
				param.Name = ElectronID + "_Syst_" + Systs.at(i);
				cout << "[FakeEstimator::executeEvent] execute " << param.Name << endl;
				executeEventFromParameter(param);
			}
		}

		//==== Normalization systematic sources ====
		if (RunNormSysts) {
			for (unsigned int i = 3; i < AnalyzerParameter::NSyst; i++) {
				param.syst_ = AnalyzerParameter::Syst(i);
				param.Name = ElectronID + "_Syst_" + Systs.at(i);
				cout << "[FakeEstimator::executeEvent] execute " << param.Name << endl;
				executeEventFromParameter(param);
			}
		}

		//==== Xsec syst ====
		if (RunXsecSyst) {
			cout << "[FakeEstimator::executeEvent] Xsec systematics is not set yet" << endl;
			exit(EXIT_FAILURE);
		}
	}			
}

void FakeEstimator::executeEventFromParameter(AnalyzerParameter param){

	if(!PassMETFilter()) return;

    Event ev = GetEvent();
    Particle METv = ev.GetMETVector();

    //==== Trigger ====
    if (! (ev.PassTrigger(HLTElecTriggerName) )) return;
	
	//==== Copy all objects ====
	vector<Muon> this_AllMuons = AllMuons;
	vector<Electron> this_AllElectrons = AllElectrons;
	vector<Jet> this_AllJets = AllJets;

	//==== Normalization Systematic Sources ====
	if (param.syst_ == AnalyzerParameter::Central) {
	}
	else if (param.syst_ == AnalyzerParameter::JetResUp) {
		this_AllJets = SmearJets( this_AllJets, +1 );
	}
	else if(param.syst_ == AnalyzerParameter::JetResDown) {
		this_AllJets = SmearJets( this_AllJets, -1 );
	}
	else if(param.syst_ == AnalyzerParameter::JetEnUp) {
		this_AllJets = ScaleJets( this_AllJets, +1 );
	}
	else if(param.syst_ == AnalyzerParameter::JetEnDown) {
		this_AllJets = ScaleJets( this_AllJets, -1 );
	}
	else if(param.syst_ == AnalyzerParameter::ElectronResUp) {
		this_AllElectrons = SmearElectrons( this_AllElectrons, +1 );
	}
	else if(param.syst_ == AnalyzerParameter::ElectronResDown) {
		this_AllElectrons = SmearElectrons( this_AllElectrons, -1 );
	}
	else if(param.syst_ == AnalyzerParameter::ElectronEnUp) {
		this_AllElectrons = ScaleElectrons( this_AllElectrons, +1 );
	}
	else if(param.syst_ == AnalyzerParameter::ElectronEnDown) {
		this_AllElectrons = ScaleElectrons( this_AllElectrons, -1 );
	}

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
}

FakeEstimator::FakeEstimator(){

}

FakeEstimator::~FakeEstimator(){

}


