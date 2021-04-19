#ifndef Selector_h
#define Selector_h

#include "AnalyzerCore.h"

class Selector : public AnalyzerCore {
public:
	// Constructor and destructor
	Selector();
	~Selector();

	// aliases

	// Initialize
	void initializeAnalyzer();
	bool RunSysts;
	bool RunDeepCSV; // PU veto will not be applied
	bool EMuTrigOnly;
	bool SkimSR3mu, SkimSR1e2mu, SkimWZ3mu, SkimWZ1e2mu;
	bool AllRegions, DiLepOnly, TriLepOnly; // on/off manually
	double weight;
	vector<TString> trigs_dblmu, trigs_emu;
	vector<TString> HcToWA_MuID, HcToWA_EleID;
	vector<TString> regions;
	vector<TString> getCuts(TString region);
	vector<TString> Systematics;

	TTree* tree;
	unsigned int nMuons; 
	double muons_pt[3];
	double muons_eta[3];
	double muons_phi[3];
	double muons_mass[3];
	int muons_charge[3];
	int muons_lepType[3];
	unsigned int nElectrons; 
	double electrons_pt[1];
	double electrons_eta[1];
	double electrons_phi[1];
	double electrons_mass[1];
	int electrons_charge[1];
	int electrons_lepType[1];
	unsigned int nJets;
	double jets_pt[14];
	double jets_eta[14];
	double jets_phi[14];
	unsigned int jets_isBtagged[14]; // 0 for not tagged, 1 for tagged
	double METv_pt, METv_phi;
	
	

	// Execute Event
    void executeEvent();
	vector<Muon> muons_all;
	vector<Electron> electrons_all;
	vector<Jet> jets_all;
	vector<Gen> gens;

	// With systematics
	void executeEventWithSystematics(const TString& syst);
	TString RegionSelector(
			Event& ev,
			vector<Muon>& muons_tight, vector<Electron>& electrons_tight,
			vector<Muon>& muons_loose, vector<Electron>& electrons_loose,
			vector<Jet>& jets, vector<Jet>& bjets, Particle& METv, const TString syst);
	double getWeight(
			const TString channel, Event& ev,
			vector<Muon>& muons, vector<Electron>& electrons,
			vector<Jet>& jets, const TString syst);
};


#endif

