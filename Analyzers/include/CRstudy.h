#ifndef CRstudy_h
#define CRstudy_h

#include "AnalyzerCore.h"

class CRstudy : public AnalyzerCore {
public:
	// constructor and destructor
	CRstudy();
	~CRstudy();

	// Initialization
	void initializeAnalyzer();
	vector<TString> trigs_dimu, trigs_emu;
	vector<TString> HcToWA_MuID, HcToWA_EleID;
	vector<TString> regions;
	vector<TString> getCuts(TString region);

	// ExecuteEvent
	void executeEvent();
	TString ControlRegionSelector(
			Event &ev,
			vector<Muon> &muons_tight, vector<Electron> &electrons_tight,
			vector<Muon> &muons_loose, vector<Electron> &electrons_loose,
			vector<Jet> &jets, vector<Jet> &bjets);

};



#endif

