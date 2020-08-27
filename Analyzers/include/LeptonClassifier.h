#ifndef LeptonClassifier_h
#define LeptonClassifier_h

#include "AnalyzerCore.h"

class LeptonClassifier : public AnalyzerCore {
public:
	// aliases
	typedef enum {DIMU, DIELEC, EMU, EJJ, MUJJ, SRSize} SR_t;
	typedef enum {CRSize} CR_t;
	typedef enum {EWPrompt, Signal, CHDecay, MisID, Conv, Others} Lepton_t;
	typedef enum {Lep, CH, Gamma} FinalState_t;

	// Gen level particle
	// not only record mother, also daughter
	typedef struct _Gen {
		Gen me;
		Gen mom;
		vector<Gen> sons;
	} Gen_t;

	// initializeAnalyzer
	void initializeAnalyzer();
	unsigned int nEvent;
	vector<TString> MuonIDs, ElectronIDs, JetIDs;

	// executeEvent
	void executeEvent();
	vector<Muon> muons, muons_tight, muons_veto;
	vector<Electron> electrons, electrons_tight, electrons_veto;
	vector<Jet> jets, jets_tight, jets_dR04;
	vector<Jet> bjets, bjets_tight, bjets_dR04;
	vector<Gen> gens; vector<Gen_t> Gens;
	Particle METv;

	double w_prefire, w_gen, w_lumi;

	// signal region selection
	bool isSignal(const SR_t& sig);
	bool isControl(const CR_t& con);

	// DrawHists
	void DrawHists(TString path, const vector<Muon>& muons, const double& weight);
	void DrawHists(TString path, const vector<Electron>& electrons, const double& weight);
	void DrawHists(TString path, const vector<Jet>& jets, const double& weight);
	void DrawHists(TString path, const Particle& METv, const double& weight);

	// other functions
	void ClearCollections();
	void InitializeWeights();
	
	// to record history
	vector<Gen_t> RewriteGens();
	vector<vector<int>> MakeHistory();
	void PrintHistory(const FinalState_t& state);

	LeptonClassifier();
	~LeptonClassifier();

};



#endif

