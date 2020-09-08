#ifndef MuonID_OS_h
#define MuonID_OS_h

#include "AnalyzerCore.h"

class MuonID_OS : public AnalyzerCore {

public:
	//==== aliases ====
	typedef enum {SR1, SR2, SRSize} SR;
	typedef enum {CR1, CR2, CRSize} CR;
	typedef enum {LL, LT, TL, TT} IDFlag;

	//==== initialization ====
	void initializeAnalyzer();
	vector<TString> trigs_POGTight; 
	vector<TString> MuonID_OSs;
	//==== execute events ====
	void executeEvent();
	// objects
	void ClearCollections();
	vector<Muon> muons, muons_tight, muons_loose, muons_veto;
	vector<Electron> electrons, electrons_veto;
	vector<Jet> jets, jets_tight, jets_dR04, jets_awayFatJet;
	vector<Jet> bjets, bjets_tight, bjets_dR04;
	vector<FatJet> fatjets, fatjets_tight, fatjets_dR10;
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
	MuonID_OS();
	~MuonID_OS();

};



#endif

