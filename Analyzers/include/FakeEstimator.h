#ifndef FakeEstimator_h
#define FakeEstimator_h

#include "AnalyzerCore.h"

class FakeEstimator : public AnalyzerCore {

public:
	typedef enum{QCDEnriched, WEnriched, ZEnriched, RegionSize} Regions;
	bool RunNPV;

	//==== initialization ====
	void initializeAnalyzer();
	vector<TString> MuonIDs;
	TString ElectronVetoID;
	// TFile* f_nPV;
	vector<TString> trigs_dlmu;

	//==== executeEvent ====
	void executeEvent();
	vector<Muon> muons_precoll;
	vector<Electron> electrons_precoll;
	vector<Jet> jets_precoll;
	vector<Gen> truth_precoll;

	//==== executeEventFromParameter====
	void executeEventFromParameter(AnalyzerParameter param);
	vector<Muon> muons_tight, muons_loose;
	vector<Electron> electrons_veto;
	vector<Jet> jets_dR04, jets_dR10;

	FakeEstimator();
	~FakeEstimator();

};



#endif

