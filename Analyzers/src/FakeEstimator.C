#include "FakeEstimator.h"

/////////////////////////////////////////////////////////////////////////////
//==== FakeEstimator                                                   ====//
//==== This analyzer is for evaluation of the fake rate of electorns   ====//
//==== from jet activity                                               ====// 
//==== using data                                                      ====//
//==== Contact: choij@cern.ch                                          ====//
//==== Based on SKFlatAnalyzer - Run2Legacy_v4                         ====// 
/////////////////////////////////////////////////////////////////////////////

// FakeEstiamtor: Central estimation of Fake Rate

void FakeEstimator::initializeAnalyzer(){

	//==== Flag for systematic sources
	RunSysts = HasFlag("RunSysts");
	RunXsecSyst = HasFlag("RunXsecSyst");

	cout << "[FakeEstimator::initializeAnalyzer] RunSysts = " << RunSysts << endl;
	cout << "[FakeEstimator::initializeAnalyzer] RunXsecSyst = " << RunXsecSyst << endl;

	//==== Systematic sources
	Systs = {"BtagDep", "JetPtCut30", "JetPtCut40", "JetPtCut50", "JetPtCut60"};
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

	f_nPV = new TFile("/home/choij/SKFlat/data/Run2Legacy_v4/2016/nPV/nPV_reweight.root");
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
		//==== need update for MET sources
		if (RunSysts) {
			for (unsigned int i = 1; i < AnalyzerParameter::NSyst; i++) {
				param.syst_ = AnalyzerParameter::Syst(i);
				param.Name = ElectronID + "_Syst_" + param.GetSystType();
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

	/////////////////////////////////////////////////////////////////////
	//==== Event Selection
	/////////////////////////////////////////////////////////////////////
	//==== QCD dominated region
	//==== 1. Exactly 1 electron
	//==== 2. At least 1 jet with Pt > 40 GeV
	//==== 3. Delta(e, j) > 1.0
	//==== 4. MET < 25 GeV, MT(e, METv) < 25 GeV
	/////////////////////////////////////////////////////////////////////
	//==== W boson dominated region
	//==== 1. Exactly 1 electron
	//==== 2. At least 1 jet with Pt > 40 GeV
	//==== 3. Delta(e, j) > 0.4
	//==== 4. MET > 50 GeV & MT(e, METv) > 70 GeV
	//////////////////////////////////////////////////////////////////////
	//==== Z boson dominated region
	//==== 1. Exactly 2 electrons
	//==== 2, At least 1 jet with Pt > 40 GeV
	//==== 3. Delta(e, j ) > 0.4
	//==== 4. |M(ee) - 91.2| < 15
	///////////////////////////////////////////////////////////////////////
	if (muons.size() != 0) return;
	double jetPtCut = 40;
	double weight = 1.;

	//==== QCD dominated & W boson dominated region
	if (electrons.size() == 1) {
		Particle elec = Particle(electrons.at(0));
		double Mt = MT(elec, METv);
		clean_jets04 = JetsVetoLeptonInside(jets, electrons_loose, muons, 0.4);
		clean_jets10 = JetsVetoLeptonInside(jets, electrons_loose, muons, 1.0);

		if (clean_jets04.size() == 0) return;
		if (clean_jets04.at(0).Pt() < jetPtCut) return;

		if (!IsDATA) {
			weight *= weight_norm_1invpb*ev.GetTriggerLumi(HLTElecTriggerName);
			weight *= ev.MCweight();
			weight *= weight_Prefire;
			cout << GetNPVReweight(ElectronID, Systs.at(2)) << endl;
			weight *= GetNPVReweight(ElectronID, Systs.at(2));
		}	
	}
		
}

FakeEstimator::FakeEstimator(){

}

FakeEstimator::~FakeEstimator(){

}

// Private Functions
double FakeEstimator::GetNPVReweight(TString id, TString syst) {
	TDirectory* temp_dir = (TDirectory*)f_nPV->GetDirectory(id + "_Central");
	TH1D* h = (TH1D*)temp_dir->Get("nPV_reweight_" + id + "_" + syst);

	if (nPV > 100) nPV = 100;
	int this_bin = nPV;

	return h->GetBinContent(this_bin);
}
