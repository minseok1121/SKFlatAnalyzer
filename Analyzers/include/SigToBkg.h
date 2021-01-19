#ifndef SigToBkg_h
#define SigToBkg_h

#include "AnalyzerCore.h"
#include "HistoMaker.h"

class SigToBkg : public AnalyzerCore {
public:
	//==== Aliases ====
	typedef enum { NONE, SR3MU, SR1E2MU } SIGNAL;

	// constructor and destructor
	SigToBkg();
	~SigToBkg();
	
	// Initialization
	void initializeAnalyzer();
	vector<TString> HcToWA_MuID, HcToWA_EleID;
	vector<TString> cuts;
	map<SIGNAL, HistoMaker*> makers;
	
	// ExectueEvent
	void executeEvent();
	SIGNAL SignalRegionSelector(
			vector<Muon> &muons_tight, vector<Electron> &electrons_tight,
			vector<Muon> &muon_loose, vector<Electron> &electrons_loose,
			vector<Jet> &jets, vector<Jet> &bjets);

	// other member functions
	void myCutflowMaker(TString region);	// region == emumu, mumumu this case
	void myHistoMaker(TString path, const vector<Muon> &muons, const double &weight);
	void myHistoMaker(TString path, const vector<Electron> &electrons, const double &weight);
	void myHistoMaker(TString path, const vector<Jet> &jets, const double &weight);
	void myHistoMaker(TString path, const Particle &part, const double &weight);
};



#endif

