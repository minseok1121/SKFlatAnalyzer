#ifndef fr_electron_h
#define fr_electron_h

#include "AnalyzerCore.h"

class fr_electron : public AnalyzerCore {

public:
	vector<TString> ElectronIDs;
	vector<TString> Systs;
	vector<TString> Prompts;
	vector<TString> Regions;

	//==== initializeAnalyzer() ====
	void initializeAnalyzer();
	bool RunSysts, RunPrompts, RunXSecSyst, RunNPV;
	TFile* f_nPV = NULL;
	TString HLTElecTriggerName;
	double TriggerSafePtCut;
	double JetPtCut;

	//==== executeEvent() ====
	void executeEvent();
	vector<Muon> AllMuons;
	vector<Jet> AllJets;
	vector<Electron> AllElectrons;
	vector<Gen> AllGens;

	TString MuonLooseID, ElectronLooseID, ElectronID, JetID, syst, prompt;
	double weight_Prefire, weight_PileUp;

	//==== executeEventFromParameter() ====
	void executeEventFromParameter(AnalyzerParameter param);
	vector<Electron> electrons, electrons_loose;
	vector<Muon> muons, muons_loose;
	vector<Jet> jets, jets_dR04, jets_dR10;

	//==== Member Functions ====
	double GetJetPtCut(const TString syst);
	double GetNPVReweight(const TString id);
	double GetCorrPt(const Electron &e);

	fr_electron();
	~fr_electron();

};



#endif

