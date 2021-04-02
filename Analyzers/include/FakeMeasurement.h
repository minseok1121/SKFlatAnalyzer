#ifndef FakeMeasurement_h
#define FakeMeasurement_h

#include "AnalyzerCore.h"

class FakeMeasurement : public AnalyzerCore {

public:
	// Constructor and Destructor
	FakeMeasurement();
	~FakeMeasurement();

	// Initialize
	void initializeAnalyzer();
	bool RunSysts;
	bool MeasMuon, MeasElectron;
	vector<TString> trigs_sglmu, trigs_sglele;
	vector<TString> HcToWA_MuID, HcToWA_EleID;
	vector<TString> Systematics;

	// ExecuteEvent
	void executeEvent();
	vector<Muon> muons_all;
	vector<Electron> electrons_all;
	vector<Jet> jets_all;
	
	// Systematics
	void executeEventWithSystematics(const TString& syst);
};



#endif

