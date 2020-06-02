#ifndef FakeValidator_h
#define FakeValidator_h

#include "AnalyzerCore.h"

class FakeValidator : public AnalyzerCore {

public:
	
	//==== initializeAnalyzer
	void initializeAnalyzer();
	vector<TString> ElectronIDs;
	TFile* f_nPV;
	TString HLTElecTriggerName;
	double TriggerSafePtCut;

	//==== executeEvent
	void executeEvent();
	vector<Muon> AllMuons;
	vector<Jet> AllJets;
	vector<Electron> AllElectrons;
	vector<Gen> AllGens;

	TString MuonID, ElectronID, JetID;
	double weight_Prefire;

	//==== executeEventFromParameter
    void executeEventFromParameter(AnalyzerParameter param);
	vector<Electron> electrons, electrons_loose;
	vector<Muon> muons;
	vector<Jet> jets, clean_jets04;

	//==== Member functions
	double GetNPVReweight(TString id, TString syst);
	double GetCorrPt(Electron e);

    FakeValidator();
    ~FakeValidator();

};



#endif

