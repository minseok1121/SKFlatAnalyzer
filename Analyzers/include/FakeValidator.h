#ifndef FakeValidator_h
#define FakeValidator_h

#include "AnalyzerCore.h"

class FakeValidator : public AnalyzerCore {

public:
	//==== initializeAnalyzer ====
	bool RunFakeSystUp, RunFakeSystDown;
	void initializeAnalyzer();
	vector <TString> ElectronIDs;
	vector <TString> idsets;
	TString HLTElecTriggerName;
	double TriggerSafePtCut1, TriggerSafePtCut2;
	TFile* f = new TFile("/home/choij/SKFlat/data/Run2Legacy_v4/2016/FakeRate/Electron/Electron_fake_rate.root");
	//==== executeEvent ====
	void executeEvent();
	vector<Muon> AllMuons;
	vector<Jet> AllJets;
	vector<Electron> AllElectrons;
	vector<Gen> AllGens;

	TString MuonID, ElectronID, ElectronTightID, JetID;
	double weight_Prefire, weight_PileUp;

	//==== executeEventFromParameter ====
	void executeEventFromParameter(AnalyzerParameter param);
	vector<Electron> electrons;
	vector<Muon> muons;
	vector<Jet> jets;

	//==== Member functions ====
	double GetCorrPt(Electron &e);
	double GetFakeRate(Electron &e, TString id, int sys);
	
	FakeValidator();
	~FakeValidator();

};



#endif

