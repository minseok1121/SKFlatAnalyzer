#ifndef Preselection_h
#define Preselection_h

#include "AnalyzerCore.h"

class Preselection : public AnalyzerCore {
public:
	// constructor and destructor
	Preselection();
	~Preselection();

	// Initialization
	void initializeAnalyzer();
	bool RunDeepCSV, RunPUVeto;
	bool SkimBaseline, SkipPUweight;
	TString path_output;
	TTree* tree;
	double ptMu1, ptMu2, ptMu3, etaMu1, etaMu2, etaMu3, phiMu1, phiMu2, phiMu3;
	double ptJ1, ptJ2, ptB1, etaJ1, etaJ2, etaJ3, phiJ1, phiJ2, phiJ3;
	double dRl1l2, dRl1l3, dRl2l3;
	double dRj1l1, dRj1l2, dRj1l3;
	double dRj2l1, dRj2l2, dRj2l3;
	double dRj1j2;
	double HT, ST, HToverST;
	double MET;
	//double mMuMu, mMuMu2;
	vector<TString> trigs_dimu, trigs_emu;
	vector<TString> HcToWA_MuID, HcToWA_EleID;
	vector<TString> regions;
	vector<TString> getCuts(TString region);
    
	// ExecuteEvent
	void executeEvent();
	TString RegionSelector(
			Event& ev,
			vector<Muon>& muons_tight, vector<Electron>& electrons_tight,
			vector<Muon>& muons_loose, vector<Electron>& electrons_loose,
			vector<Jet>& jets, vector<Jet>& bjets);

	void finalizeAnalyzer();
};

#endif

