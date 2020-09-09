#ifndef WWW_SS_h
#define WWW_SS_h

#include "AnalyzerCore.h"

class WWW_SS : public AnalyzerCore {

public:
	//==== aliases ====
	typedef enum {SR1, SR2, SRSize} SR;
	typedef enum {CR1, CR2, CRSize} CR;

	//==== initialization ====
	void initializeAnalyzer();
	vector<TString> trigs_POGTight; 
	vector<TString> WWW_SSs;
	//==== execute events ====
	void executeEvent();
	// objects
	void ClearCollections();
	vector<Muon> muons, muons_tight, muons_loose, muons_veto;
	vector<Electron> electrons, electrons_veto;
	vector<Jet> jets, jets_tight, jets_dR04;
	vector<Jet> bjets, bjets_dR04;
	vector<Gen> gens;
	Particle METv;
	
	// To do:: trigger sf
	// To do:: measure fake contribution

	//==== other functions ====
	bool isSignal(const SR& sig) const;
	bool isControl(const CR& con) const;

	// update muon collection
	//vector<Muon> UpdateMuColl(const vector<Muon>& muons, 

	// DrawHists
	void DrawHists(TString path, const vector<Muon>& muons, const double& weight);
	void DrawHists(TString path, const vector<Electron>& electrons, const double& weight);
	void DrawHists(TString path, const vector<Jet>& jets, const double& weight);
	void DrawHists(TString path, const vector<FatJet>& jets, const double& weight);
	void DrawHists(TString path, const Particle& METv, const double& weight);

	//==== constructor, destructor ====
	WWW_SS();
	~WWW_SS();

};



#endif

