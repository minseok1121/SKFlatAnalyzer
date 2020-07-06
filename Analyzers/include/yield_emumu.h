#ifndef yield_emumu_h
#define yield_emumu_h

#include "AnalyzerCore.h"

class yield_emumu : public AnalyzerCore {

public:

	//==== initializeAnalyzer() ====
	void initializeAnalyzer();
	vector<TString> MuonIDs, ElectronIDs;
	TString HLTEMuTriggerName;
	double TriggerSafeElecPtCut, TriggerSafeMuPtCut;

	//==== executeEvent() ====
	void executeEvent();
	vector<Muon> AllMuons;
    vector<Electron> AllElectrons;
    vector<Jet> AllJets;
    vector<Gen> AllGens;

    TString MuonLooseID, ElecLooseID, MuonTightID, ElecTightID, JetID;

	//==== executeEventFromParamter() ====
  	void executeEventFromParameter(AnalyzerParameter param);
	vector<Electron> electrons, electrons_loose;
    vector<Muon> muons, muons_loose;
    vector<Jet> jets, jets_dR04;

 	yield_emumu();
    ~yield_emumu();

};



#endif

