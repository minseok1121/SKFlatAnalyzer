#ifndef Selector_h
#define Selector_h

#include "AnalyzerCore.h"

class Selector : public AnalyzerCore {
public:
	// Constructor and destructor
	Selector();
	~Selector();

	// Initialize
	void initializeAnalyzer();
	bool SkimDiMu, SkimEMu, Skim3Mu, Skim1E2Mu;
	bool SkimSignal;
	vector<TString> trigs_dblmu, trigs_emu;
	vector<TString> HcToWA_MuID, HcToWA_EleID;

	TTree* tree;
	unsigned int nMuons;
	bool passDblMuTrigs, passEMuTrigs;
	float muons_pt[3];
	float muons_eta[3];
	float muons_phi[3];
	float muons_mass[3];
	float muons_miniIso[3];
	float muons_scaleUp[3];
	float muons_scaleDown[3];
	int muons_charge[3];
	int muons_lepType[3];
	bool muons_isTight[3];

	unsigned int nElectrons;
	float electrons_pt[1];
	float electrons_eta[1];
	float electrons_phi[1];
	float electrons_mass[1];
	float electrons_miniIso[1];
	float electrons_scaleUp[1];
	float electrons_scaleDown[1];
	float electrons_smearUp[1];
	float electrons_smearDown[1];
	int electrons_charge[1];
	int electrons_lepType[1];
	bool electrons_isTight[1];

	unsigned int nJets;
	float jets_pt[20];
	float jets_eta[20];
	float jets_phi[20];
	float jets_mass[20];
	float jets_scaleUp[20];
	float jets_scaleDown[20];
	float jets_smearUp[20];
	float jets_smearDown[20];
	float jets_btagScore[20];
	bool jets_isBtagged[20];

	float METv_pt, METv_eta, METv_phi;
	float A_pt, A_eta, A_phi, A_mass;
	float Hc_pt, Hc_eta, Hc_phi, Hc_mass;

	float genWeight;
	float trigLumi;
	// 0 for central, 1 for up, 2 for down
	float weights_L1Prefire[3];
	float weights_id[3];
	float weights_trigs_dblmu[3];
	float weights_trigs_emu[3];
	float weights_btag[3];

	void executeEvent();
};



#endif

