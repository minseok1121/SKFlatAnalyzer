#ifndef id_var_prompt_h
#define id_var_prompt_h

#include "AnalyzerCore.h"

class id_var_prompt : public AnalyzerCore {

public:

	//==== initializeAnalyze() ====
	void initializeAnalyzer();
	bool RunMuon, RunElectron;

	//==== executeEvent() ====
	void executeEvent();
	vector<Muon> AllMuons;
	vector<Electron> AllElectrons;
	vector<Jet> AllJets;
	vector<Gen> Gens;

	//==== executeEventFromParameter() ====
	void executeEventFromParameter(AnalyzerParameter param);
	double weight;
	double w_gen, w_kfactor, w_toppt, w_lumi, w_pileup, w_prefire;
	double sf_trig, sf_muid, sf_mutk, sf_muiso, sf_elreco, sf_elid, sf_btag;

	id_var_prompt();
	~id_var_prompt();

	//==== member function ====
	bool IsDY(const vector<Muon> &muons, const vector<Electron> &electrons_veto);
	bool IsDY(const vector<Electron> &electrons, const vector<Muon> &muons_veto);
	bool IsTT(const vector<Muon> &muons, const vector<Electron> &electrons_veto, const vector<Jet> &jets, const vector<Jet> &bjets, const Particle &METv);
	bool IsTT(const vector<Electron> &electrons, const vector<Muon> &muons_veto, const vector<Jet> &jets, const vector<Jet> &bjets, const Particle &METv);
	void DrawHists(TString path, const Muon &mu, unsigned int order, int type, const double weight);
	void DrawHists(TString path, const Electron &e, unsigned int order, int type, const double weight);
	void DrawHists(TString path, const Jet &j, unsigned int order, const double weight);
};



#endif

