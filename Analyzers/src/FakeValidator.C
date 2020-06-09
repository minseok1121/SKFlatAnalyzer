#include "FakeValidator.h"

void FakeValidator::initializeAnalyzer(){
	//==== No systemtaic setup yet =====
	RunFake = HasFlag("RunFake");
	cout << "[FakeValidator::initializeAnalyzer] RunFake = " << RunFake << endl;
	
	//==== Electron ID setting ====
	ElectronIDs = {"passLooseID", "passTightID", "FakeLooseID", "FakeTightID"};
	idsets = {"POG", "Fake"};

	//==== Trigger setting
	if (DataYear == 2016) {
		HLTElecTriggerName = "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v";
		TriggerSafePtCut1 = 25.;
		TriggerSafePtCut2 = 15.;
	}
	else if (DataYear == 2017) {
		cout << "[FakeValidator::initializeAnalyzer] Trigger is not set for 2017" << endl;
		exit(EXIT_FAILURE);
	}
	else if (DataYear == 2018) {
		cout << "[FakeValidator::initializeAnalyzer] Trigger is not set for 2018" << endl;
		exit(EXIT_FAILURE);
	}
	else {
		cout << "[FakeValidator::initializeAnalyzer] Wrong Year" << endl;
		exit(EXIT_FAILURE);
	}
	
	cout << "[FakeValidator::initializeAnalyzer] HLTElecTriggerName  = " << HLTElecTriggerName << endl;
	cout << "[FakeValidator::initializeAnalyzer] TriggerSafePtCut1 = " << TriggerSafePtCut1 << endl;
	cout << "[FakeValidator::initializeAnalyzer] TriggerSafePtCut2 = " << TriggerSafePtCut2 << endl;

	//==== No b-tagging ====
	
}

void FakeValidator::executeEvent(){

	//==== Copy all objects ====
	AllMuons = GetAllMuons();
	AllElectrons = GetAllElectrons();
	AllJets = GetAllJets();

	//==== Get L1Prefire, pileup reweight ====
	weight_Prefire = GetPrefireWeight(0);
	weight_PileUp = GetPileUpWeight(nPileUp, 0);

	AnalyzerParameter param;

	//==== Loop over ElectornID sets =====
	for (unsigned int i = 0; i < idsets.size(); i++) {
		MuonID = "POGLoose";
		if (idsets.at(i).Contains("POG")) {
			ElectronID = ElectronIDs.at(i);
			ElectronTightID = ElectronIDs.at(i+1);
		}
		else if (idsets.at(i).Contains("Fake")) {
			ElectronID = ElectronIDs.at(2*i);
			ElectronTightID = ElectronIDs.at(2*i + 1);
		}
		JetID = "tight";
	
		param.Clear();
		param.syst_ = AnalyzerParameter::Central;
		param.Name = idsets.at(i) + "_Central";
		param.Jet_ID = JetID;

		executeEventFromParameter(param);
	}

}

void FakeValidator::executeEventFromParameter(AnalyzerParameter param){

	if(!PassMETFilter()) return;

	Event ev = GetEvent();
	Particle METv = ev.GetMETVector();

	//==== Trigger ====
	if (! ev.PassTrigger(HLTElecTriggerName)) return;

	//==== Copy all objects ====
	vector<Muon> this_AllMuons = AllMuons;
	vector<Electron> this_AllElectrons = AllElectrons;
	vector<Jet> this_AllJets = AllJets;

	//==== ID selection ====
	muons = SelectMuons(this_AllMuons, MuonID, 10, 2.4);
	electrons = SelectElectrons(this_AllElectrons, ElectronID, 15, 2.5);
	jets = SelectJets(this_AllJets, param.Jet_ID, 20., 2.4);
	
	std::sort(muons.begin(), muons.end(), PtComparing);
	std::sort(electrons.begin(), electrons.end(), PtComparing);
	std::sort(jets.begin(), jets.end(), PtComparing);

	//==== Event Selectrion ====
	//==== 1. 3 loose electrons, no loose muon ====
	//==== 2. pt(1) > 25, pt(2) > 15, pt(3) > 15 ====
	//==== 3. At least 1 opposite charge electron pair with on-Z
	if (muons.size() != 0) return;
	if (electrons.size() != 3) return;
	if (electrons.at(0).Pt() <= TriggerSafePtCut1) return;
	if (electrons.at(1).Pt() <= TriggerSafePtCut2) return;
	if (electrons.at(2).Pt() < 15) return;

	double charge[3];
	Particle pair[3];
	bool isOnZ = false;
	for (int i = 0; i < 3; i++) {
		charge[i] = electrons.at(i%3).Charge() * electrons.at((i+1)%3).Charge();
		pair[i] = electrons.at(i%3) + electrons.at((i+1)%3);

		if (charge[i] < 0 && IsOnZ(pair[i].M(), 10.)) isOnZ = true;
	}

	if (!isOnZ) return;

	double weight = 1.;
	if (!IsDATA) {
		weight *= weight_norm_1invpb*ev.GetTriggerLumi("Full");
		weight *= ev.MCweight();
		weight *= weight_Prefire;

		//ID scale factor
	}

	//==== Fill hist ====
	JSFillHist(param.Name, "1st_electron_pt_" + param.Name, electrons.at(0).Pt(), weight, 140, 0., 140.);
	JSFillHist(param.Name, "2nd_electron_pt_" + param.Name, electrons.at(1).Pt(), weight, 140, 0., 140.);
	JSFillHist(param.Name, "3rd_electron_pt_" + param.Name, electrons.at(2).Pt(), weight, 140, 0., 140.);
	JSFillHist(param.Name, "1st_electron_eta_" + param.Name, electrons.at(0).Eta(), weight, 50, -2.5, 2.5);
	JSFillHist(param.Name, "2nd_electron_eta_" + param.Name, electrons.at(1).Eta(), weight, 50, -2.5, 2.5);
	JSFillHist(param.Name, "3rd_electron_eta_" + param.Name, electrons.at(2).Eta(), weight, 50, -2.5, 2.5);

}

FakeValidator::FakeValidator(){

}

FakeValidator::~FakeValidator(){

}


