#ifndef BaselineSelector_h
#define BaselineSelector_h

#include "AnalyzerCore.h"

class BaselineSelector : public AnalyzerCore {

public:
	typedef enum { None, SR_3mu, SR_1e2mu, CR_DYdimu, CR_TTdimu, CR_TTemu, CR_SSdimu, CR_SSemu, CR_Conv3mu, CR_Conv1e2mu, CR_RevBjet } REGION;
	//==== constructor and destructor
	BaselineSelector();
	~BaselineSelector();

	//==== Initialization
	void initializeAnalyzer();
	bool RunDeepCSV;
	bool SkimBaseline;
	vector<TString> trigs_dimu, trigs_emu;
	vector<TString> HcToWA_MuID, HcToWA_EleID, POG_MuID, POG_EleID;
	vector<TString> regions;
	vector<TString> getCuts(const TString& region);
	TTree* Events = nullptr;
	double ptl1, ptl2, ptl3, etal1, etal2, etal3, phil1, phil2, phil3;
	double ptj1, ptj2, ptb1, etaj1, etaj2, etab1, phij1, phij2, phib3;
	double dRl1l2, dRl1l3, dRl2l3;
	double dRj1l1, dRj1l2, dRj1l3;
	double dRj2l1, dRj2l2, dRj2l3;
	double dRj1j2;
	double HT, LT, MET, ST, HToverST, LToverST, METoverST;
	// mass?
	
	//==== ExecuteEvent
	void executeEventFromParameter(AnalyzerParameter param);
	void executeEvent();
	Event ev; vector<Gen> gens;
	vector<Muon> muons; vector<Electron> electrons;
	vector<Jet> jets; Particle METv;

	REGION RegionSelector(
			const Event& ev,
			const vector<Muon>& muons_tight, const vector<Electron>& electrons_tight,
			const vector<Muon>& muons_loose, const vector<Electron>& electrons_loose,
			const vector<Jet>& jets, const vector<Jet>& bjets);
};



#endif

