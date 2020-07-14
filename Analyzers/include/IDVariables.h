#ifndef IDVaribales_h
#define IDVaribales_h

#include "AnalyzerCore.h"

class IDVariables : public AnalyzerCore {

public:
	//==== initializeAnalyzer() ====
	void initializeAnalyzer();
	vector<TString> MuonIDs;

	vector<TString> TrigNames;
	double TriggerSafePtCut1, TriggerSafePtCut2;

	//==== executeEvent() ====
	void executeEvent();
	vector<Muon> AllMuons;
	vector<Electron> AllElectrons;
	vector<Jet> AllJets;
	vector<Gen> AllGens;
	
	TString MuonID, JetID;
	TString MuonVetoID, ElecVetoID;

	//==== ExecuteEventFromParameter() ====
	void executeEventFromParameter(AnalyzerParameter param);
	double weight;
	double w_gen, w_filter, w_toppt, w_lumi, w_pileup, w_prefire;
	double sf_trig, sf_muid, sf_mutk, sf_muiso, sf_elreco, sf_elid, sf_btag;

	IDVariables();
	~IDVariables();

	//==== member functions ====
	bool IsDY(const vector<Muon> &muons, const vector<Electron> &electrons_veto);
	bool IsTTbar(const vector<Muon> &muons, const vector<Electron> &electrons_veto, const vector<Jet> &jets, const vector<Jet> &bjets);
	void DrawHists(TString path, const Muon &mu, unsigned int order, const double weight);
	void DrawHists(TString path, const Jet &j, unsigned int order, const double weight);
};



#endif

