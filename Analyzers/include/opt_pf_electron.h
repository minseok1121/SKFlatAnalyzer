#ifndef opt_pf_electron_h
#define opt_pf_electron_h

#include "AnalyzerCore.h"

class opt_pf_electron : public AnalyzerCore {

public:
	//==== initializeAnalyzer() ====
	void initializeAnalyzer();
	bool RunSysts;
	vector<TString> ElectronIDs;
	vector<TString> TrigNames;
	double TriggerSafePtCut;
	
	//==== executeEvent() ====
	void executeEvent();
	vector<Muon> AllMuons;
	vector<Jet> AllJets;
	vector<Electron> AllElectrons;
	vector<Gen> AllGens;

	TString MuonID, ElectronID, JetID, MuonVetoID, ElecVetoID;
	double w_prefire, w_pileup;

	//==== executeEventFromParameter() ====
	void executeEventFromParameter(AnalyzerParameter param);
	vector<Electron> electrons, electrons_loose;
	vector<Muon> muons, muons_loose;
	vector<Jet> jets, jets_dR04, jets_dR10;

	//==== Member Functions ====
	bool IsPromptEvent(const vector<Electron> &electrons, const vector<Electron> &electrons_veto, const vector<Jet> &jets, const Particle &METv);
	bool IsFakeEvent(const vector<Electron> &electrons, const vector<Electron> &electrons_veto, const vector<Jet> &jets, const Particle &METv);
	bool IsOnZEvent(const vector<Electron> &electrons, const vector<Jet> &jets);

	void DrawHists(TString path, const Electron &e, unsigned int order, const double weight);
	void DrawHists(TString path, const Jet &j, unsigned int order, const double weight);
	void DrawHists(TString path, const Particle &METv, const Electron &e, const double weight);


	opt_pf_electron();
	~opt_pf_electron();
};



#endif

