#ifndef FakeEstimator_h
#define FakeEstimator_h

#include "AnalyzerCore.h"

class FakeEstimator : public AnalyzerCore {

public:
	// Constructor and destructor
	FakeEstimator();
	~FakeEstimator();

	// Initialize
	void initializeAnalyzer();;
	bool SkimSglMu, SkimSglEle, SkimDblMu, SkimDblEle;
	vector<TString> trigs_sglmu, trigs_sglele;
	vector<TString> HcToWA_MuID, HcToWA_EleID;

	TTree* tree;
	// Events
	bool passMu8Path, passMu17Path;
	bool passEle8Path, passEle12Path, passEle23Path;
	float METv_pt, METv_eta, METv_phi;

	// Muons
	unsigned int nMuons;
	float muons_pt[2], muons_eta[2], muons_phi[2], muons_mass[2], muons_miniIso[2];
	float muons_scaleUp[2], muons_scaleDown[2];
	int muons_charge[2], muons_lepType[2];
	bool muons_isTight[2];

	// Electrons
	unsigned int nElectrons;
	float electrons_pt[2], electrons_eta[2], electrons_phi[2], electrons_mass[2], electrons_miniIso[2];
	float electrons_scaleUp[2], electrons_scaleDown[2], electrons_smearUp[2], electrons_smearDown[2];
	int electrons_charge[2], electrons_lepType[2];
	bool electrons_isTight[2];

	// Jets
	unsigned int nJets;
	float jets_pt[24], jets_eta[24], jets_phi[24], jets_mass[24];
	float jets_scaleUp[24], jets_scaleDown[24], jets_smearUp[24], jets_smearDown[24];
	bool jets_isBtagged[24];
	
	// Weights
	float genWeight, trigLumi;
	// 0 for central, 1 for up, 2 for down
	float weights_L1Prefire[3], weights_btag[3];

	// ExecuteEvent
	void executeEvent();
};



#endif

