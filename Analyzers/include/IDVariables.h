#ifndef IDVaribales_h
#define IDVaribales_h

#include "AnalyzerCore.h"

class IDVariables : public AnalyzerCore {

public:
	//==== initializeAnalyzer() ====
	void initializeAnalyzer();
	vector<TString> ElectronIDs, ElectronMVAIDs;
	vector<TString> MuonIDs;

	//==== executeEvent() ====
	void executeEvent();
	vector<Muon> AllMuons;
	vector<Electron> AllElectrons;
	vector<Jet> AllJets;
	vector<Gen> AllGens;
	
	double weight_PileUp, weight_Prefire, weight_TopPt;
	TString MuonID, ElectronID, JetID;
	TString MuonVetoID, ElectronVetoID;
	void executeEventFromParameter(AnalyzerParameter param);

	IDVariables();
	~IDVariables();

	//==== member functions ====
	void DrawIDVariables(TString path, const Electron &e, unsigned int order, const double weight);
	void DrawIDVariables(TString path, const Muon &mu, unsigned int order, const double weight);
};



#endif

