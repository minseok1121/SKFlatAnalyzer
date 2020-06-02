#include "FakeValidator.h"

void FakeValidator::initializeAnalyzer(){

	ElectronIDs = {"passLooseID", "passTightID"};
	
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

	//==== nPV reweight
    if (DataYear == 2016) {
        f_nPV = new TFile("/home/choij/SKFlat/data/Run2Legacy_v4/2016/nPV/nPV_reweight.root");
    }
    else if (DataYear == 2017 || DataYear == 2018) {
        cout << "[FakeEstimator::initializeAnalyzer] nPV_reweight is not set yet" << endl;
        exit(EXIT_FAILURE);
    }
    else {
        cout << "[FakeEstiamtor::initializeAnalyzer] Wrong Year" << endl;
        exit(EXIT_FAILURE);
    }
}

void FakeValidator::executeEvent(){

	//==== Copy all objects ====
    AllMuons = GetAllMuons();
    AllElectrons = GetAllElectrons();
    AllJets = GetAllJets();

	//==== Get L1Prefire. pileup reweight ====
    weight_Prefire = GetPrefireWeight(0);	

    AnalyzerParameter param;

	//==== Loop over electron IDs ====
	for (unsigned int i = 0; i < ElectronIDs.size(); i++) {
		MuonID = "POGLoose"; // to veto loose muons
        ElectronID = ElectronIDs.at(i);
        JetID = "tight";
        
		param.Clear();
		param.syst_ = AnalyzerParameter::Central;
		param.Name = ElectronID + "_Central";
	

  		executeEventFromParameter(param);
	}
}

void FakeValidator::executeEventFromParameter(AnalyzerParameter param){

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
	electrons_loose = SelectElectrons(this_AllElectrons, "passLooseID", 25, 2.5);
    jets = SelectJets(this_AllJets, param.Jet_ID, 20., 2.4);
	clean_jets04 = JetsVetoLeptonInside(jets, electrons_loose, muons, 0.4);

	std::sort(muons.begin(), muons.end(), PtComparing);
    std::sort(electrons.begin(), electrons.end(), PtComparing);
    std::sort(electrons_loose.begin(), electrons_loose.end(), PtComparing);
    std::sort(jets.begin(), jets.end(), PtComparing);

	//==== leading electron should pass trigger-safe pt cut ====
    if (electrons.size() == 0) return;
    if (electrons.at(0).Pt() <= TriggerSafePtCut) return;

	//==== Event Selection ====
	if (electrons.size() != 3) return; // trielectron
	if (clean_jets04.size() == 0) return; // at least one clean jet
	if (muons.size() != 0) return; // no loose muon

	// distributions to draw?
}

FakeValidator::FakeValidator(){

}

FakeValidator::~FakeValidator(){

}


