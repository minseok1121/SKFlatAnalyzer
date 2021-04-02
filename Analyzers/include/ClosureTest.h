#ifndef ClosureTest_h
#define ClosureTest_h

#include "AnalyzerCore.h"

class ClosureTest : public AnalyzerCore {

public:
	// Constructor and destructor
	ClosureTest();
	~ClosureTest();

	// Initialize
	void initializeAnalyzer();
	bool RunTrigger, RunFakeRate;
	vector<TString> trigs_dblmu, trigs_emu;

	TTree* tree;
	double genWeight, L1PrefireWeight;
	double trigDblMuEff, trigDblMuEffUp, trigDblMuEffDown;
	double trigEMuEff, trigEMuEffUp, trigEMuEffDown;
	bool passDblMuTrigs, passEMuTrigs;
	unsigned int nMuons;
	double muons_pt[3], muons_eta[3], muons_phi[3], muons_mass[3];	
	unsigned int nElectrons;
	double electrons_pt[1], electrons_eta[1], electrons_phi[1], electrons_mass[1];

	// ExecuteEvent
	void executeEvent();
};



#endif

