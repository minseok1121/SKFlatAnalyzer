#ifndef FakeEstimator_h
#define FakeEstimator_h

#include "AnalyzerCore.h"

class FakeEstimator : public AnalyzerCore {

public:
	// Constructor and destructor
	FakeEstimator();
	~FakeEstimator();
	
	// Initialize
	void initializeAnalyzer();
	bool RunSysts;
	bool SkimSglMu, SkimSglEle, SkimDblMu, SkimDblEle;
	vector<TString> trigs_sglmu, trigs_sglele;
	vector<TString> HcToWA_MuID, HcToWA_EleID;
	vector<TString> Systematics;

	TTree* tree;
	double genWeight, L1PrefireWeight;
	int passMu8Path, passMu17Path;
	int passEle8Path, passEle12Path, passEle23Path;
	unsigned int nMuons;
	double muons_pt[2], muons_eta[2], muons_phi[2], muons_mass[2];
	double muons_miniIso[2], muons_relIso[2];
	int muons_charge[2];
	int muons_lepType[2];
	int muons_isTight[2];

	unsigned int nElectrons;
	double electrons_pt[2], electrons_eta[2], electrons_phi[2], electrons_mass[2];
	double electrons_miniIso[2], electrons_relIso[2];
	int electrons_charge[2];
	int electrons_lepType[2];
	int electrons_isTight[2];

	unsigned int nJets;
	double jets_pt[28], jets_eta[28], jets_phi[28];
	unsigned int jets_isBtagged[28];
	double METv_pt, METv_phi;

	// ExecuteEvent
	void executeEvent();
	vector<Muon> muons_all;
	vector<Electron> electrons_all;
	vector<Jet> jets_all;
	vector<Gen> gens;
	
	// Systematics
	void executeEventWithSystematics(const TString& syst);
};



#endif

