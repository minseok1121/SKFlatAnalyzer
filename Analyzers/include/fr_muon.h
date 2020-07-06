#ifndef fr_muon_h
#define fr_muon_h

#include "AnalyzerCore.h"

class fr_muon : public AnalyzerCore {

public:
	vector<TString> MuonIDs;
    vector<TString> Systs;
    vector<TString> Prompts;
    vector<TString> Regions;

	//==== initializeAnalyzer() ====
    void initializeAnalyzer();
	bool RunSysts, RunPrompts, RunXSecSyst, RunNPV;
    TFile* f_nPV = NULL;
    TString HLTMuonTriggerName[2];
    double TriggerSafePtCut[2];
    double JetPtCut;

	//==== executeEvent() ====
    void executeEvent();
    vector<Muon> AllMuons;
    vector<Jet> AllJets;
    vector<Electron> AllElectrons;
    vector<Gen> AllGens;

	TString MuonLooseID, ElectronLooseID, MuonID, JetID, syst, prompt;
    double weight_Prefire, weight_PileUp;

	//==== executeEventFromParameter() ====
    void executeEventFromParameter(AnalyzerParameter param);
    vector<Electron> electrons, electrons_loose;
    vector<Muon> muons, muons_loose;
    vector<Jet> jets, jets_dR04, jets_dR10;

	//==== Member Functions ====
    double GetJetPtCut(const TString syst);
    //double GetNPVReweight(const TString id);
    double GetCorrPt(const Muon &mu);

	fr_muon();
  	~fr_muon();

};



#endif

