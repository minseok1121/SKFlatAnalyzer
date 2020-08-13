#ifndef PromptSelector_h
#define PromptSelector_h

#include "AnalyzerCore.h"

class PromptSelector : public AnalyzerCore {

public:

    void initializeAnalyzer();
    bool RunTTLL, RunTTLJ;

    void executeEvent();
    vector<Muon> AllMuons;
    vector<Electron> AllElectrons;
    vector<Jet> AllJets;
    vector<Gen> gens;
	enum LEPTYPE {ALL, FAKE, EWPROMPT, SIGDAUGHTER, EWTAUDAUGHTER, OTHERS};
	const bool JET = false, BJET = true;
	double wegiht_Prefire;
    double sf_trig, sf_mutk, sf_muid, sf_muiso, sf_elreco, sf_elid, sf_btag;
    PromptSelector();
    ~PromptSelector();

    // member functions
	// for RunTTLL
    bool isDiMu(const vector<Muon>& muons, const vector<Muon>& muons_veto, const vector<Electron>& electrons_veto);
    bool isDiElec(const vector<Electron>& electrons, const vector<Electron>& electrons_veto, const vector<Muon>& muons_veto);
    bool isEMu(const vector<Electron>& electrons, const vector<Muon>& muons, const vector<Electron>& electrons_veto, const vector<Muon>& muons_veto);
	bool isMuJJ(const vector<Muon>& muons, const vector<Muon>& muons_veto, const vector<Electron>& electrons_veto);
	bool isElecJJ(const vector<Electron>& electrons, const vector<Electron>& electrons_veto, const vector<Muon>& muons_veto);

	// DrawHist
	void DrawHists(TString path, const vector<Electron>& electrons, const vector<Gen>& gens, const double& weight, LEPTYPE TYPE);
	void DrawHists(TString path, const vector<Muon>& muons, const vector<Gen>& gens, const double& weight, LEPTYPE TYPE);
	void DrawHists(TString path, const vector<Jet>& jets, const double& weight, const bool& isBJET);
	void DrawHists(TString path, const Particle& METv, const double& weight);
	int GetLeptonMotherType(const Lepton& lep);
};



#endif

