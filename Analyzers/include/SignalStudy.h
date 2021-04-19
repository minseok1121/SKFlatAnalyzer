#ifndef SignalStudy_h
#define SignalStudy_h

#include "AnalyzerCore.h"

class SignalStudy : public AnalyzerCore {

public:
	// Constructor and destructor
	SignalStudy();
	~SignalStudy();

	// Initialize
    void initializeAnalyzer();
	bool Skim3Mu, Skim1E2Mu;
	vector<TString> trigs_dblmu, trigs_emu;
	vector<TString> HcToWA_MuID, HcToWA_EleID;
	
	void executeEvent();
};



#endif

