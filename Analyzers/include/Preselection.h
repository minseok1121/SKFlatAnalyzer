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
	bool RunLowPtJet, RunDeepCSV;
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
};

#endif

