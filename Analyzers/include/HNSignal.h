#ifndef HNSignal_h
#define HNSignal_h

#include "AnalyzerCore.h"

class HNSignal : public AnalyzerCore {

public:
	//==== aliases ====
	typedef enum {SR1, SR2, SRSize} SR;
	typedef enum {CR1, CR2, CRSize} CR;

	//==== initialization ====
	bool RunPOGTight, RunHighPt;
	void initializeAnalyzer();
	vector<TString> trigs_POGTight, trigs_HighPt; 
	vector<TString> HNType1_POGTight, HNType1_HighPt, HNType1_Electron, HNType1_Jet, HNType1_FatJet;
	vector<TString> HNType1_IDSets;

	TFile* f_trig_sf_lead_2016_BtoF = nullptr;
	TFile* f_trig_sf_tail_2016_BtoF = nullptr;
	TFile* f_trig_sf_lead_2016_GtoH = nullptr;
	TFile* f_trig_sf_tail_2016_GtoH = nullptr;
	TFile* f_trig_sf_lead_2017 = nullptr;
	TFile* f_trig_sf_tail_2017 = nullptr;
	TFile* f_trig_sf_lead_2018 = nullptr;
	TFile* f_trig_sf_tail_2018 = nullptr;
	TFile* f_trig_sf_highpt_2016 = nullptr;
	TFile* f_trig_sf_highpt_2017 = nullptr;
	TFile* f_trig_sf_highpt_2018 = nullptr;
	TFile* f_fake_2016 = nullptr;
	TFile* f_fake_2017 = nullptr;
	TFile* f_fake_2018 = nullptr;

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
	
	// weights
	void ClearWeights();
	double w_prefire,  w_gen, w_lumi, w_pileup;  
	double sf_trig, sf_muid, sf_muiso, sf_btag;
	
	// To do:: trigger sf
	// To do:: measure fake contribution

	//==== other functions ====
	bool isSignal(const SR& sig) const;
	bool isControl(const CR& con) const;
	double GetFakeWeight(const TString& TightID);
	double GetTrigSF(const vector<Muon>& muons);

	// DrawHists
	void DrawHists(TString path, const vector<Muon>& muons, const double& weight);
	void DrawHists(TString path, const vector<Electron>& electrons, const double& weight);
	void DrawHists(TString path, const vector<Jet>& jets, const double& weight);
	void DrawHists(TString path, const vector<FatJet>& jets, const double& weight);
	void DrawHists(TString path, const Particle& METv, const double& weight);

	vector<Muon> MuonPromptOnlyHNtype1(const vector<Muon>& muons);
	//==== constructor, destructor ====
	HNSignal();
	~HNSignal();

};



#endif

