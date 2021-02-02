#ifndef SigToBkg_h
#define SigToBkg_h

#include "AnalyzerCore.h"

class SigToBkg : public AnalyzerCore {
public:
	// constructor and destructor
	SigToBkg();
	~SigToBkg();
	
	// Initialization
	void initializeAnalyzer();
	vector<TString> trigs_dimu, trigs_emu;
	vector<TString> HcToWA_MuID, HcToWA_EleID;
	vector<TString> regions;
	vector<TString> cuts;
	
	// ExectueEvent
	void executeEvent();
	TString SignalRegionSelector(
			Event &ev,
			vector<Muon> &muons_tight, vector<Electron> &electrons_tight,
			vector<Muon> &muon_loose, vector<Electron> &electrons_loose,
			vector<Jet> &jets, vector<Jet> &bjets);

};



#endif

