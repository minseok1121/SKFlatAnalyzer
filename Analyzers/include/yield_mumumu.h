#ifndef yield_mumumu_h
#define yield_mumumu_h

#include "AnalyzerCore.h"

class yield_mumumu : public AnalyzerCore {

public:

	//==== initializeAnalyzer() ====
	void initializeAnalyzer();
	vector<TString> MuonIDs, ElectronIDs;
	TString HLTMuonTriggerName;
	double TriggerSafePtCut1, TriggerSafePtCut2;

	//==== executeEvent() ====
	void executeEvent();
	vector<Muon> AllMuons;
	vector<Electron> AllElectrons;
	vector<Jet> AllJets;
	vector<Gen> AllGens;

	TString MuonLooseID, ElecLooseID, MuonTightID, ElecTightID, JetID;
	
	//==== executeEventFromParameter(param) ====
	void executeEventFromParameter(AnalyzerParameter param);
	vector<Electron> electrons, electrons_loose;
	vector<Muon> muons, muons_loose;
	vector<Jet> jets, jets_dR04;

	typedef struct _Pair {
		double Mass;
		bool isOS;
		bool isAbove12;
		bool isOnZ;
		bool isOffZ;
	} Pair;

	yield_mumumu();
	~yield_mumumu();

};



#endif

