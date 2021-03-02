#ifndef RegionSelector_h
#define RegionSelector_h

#include "AnalyzerCore.h"

class RegionSelector : public AnalyzerCore {
public:
	// Constructor and destructor
	RegionSelector();
	~RegionSelector();

	// Initialize
	void initializeAnalyzer();
	bool RunDeepCSV; // PU veto will not be applied
	bool EMuTrigOnly;
	bool SkimBaseline; // No need to set for a moment, after WZ trilep control region is finished
	double weight;
	vector<TString> trigs_dblmu, trigs_emu;
	vector<TString> HcToWA_MuID, HcToWA_EleID;
	vector<TString> regions;
	vector<TString> getCuts(TString region);

	// Execute Event
    void executeEvent();
	TString Selector(
			Event& ev,
			vector<Muon>& muons_tight, vector<Electron>& electrons_tight,
			vector<Muon>& muons_loose, vector<Electron>& electrons_loose,
			vector<Jet>& jets, vector<Jet>& bjets, Particle& METv) const;
};


#endif

