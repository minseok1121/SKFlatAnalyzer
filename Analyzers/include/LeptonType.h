#ifndef LeptonType_h
#define LeptonType_h

#include "AnalyzerCore.h"

class LeptonType : public AnalyzerCore {
public:
	// type as enum, cut as vector
	typedef enum {JET, BJET, JetTypeSize} JetType; 
	typedef enum {ALL, FAKE, EWPROMPT, SIGNAL, EWTAU, OTHERS, LepTypeSize} LepType;
	typedef enum {DIMU, DIELEC, EMU, MUJJ, EJJ, SignalSize} Signal; 
	typedef enum {ControlSize} Control;

	typedef struct _Gen {
		const Gen* me;
		unsigned int mother_idx;
		vector<unsigned int> daughters_idx;
	} Gen_t;

    void initializeAnalyzer();
	vector<TString> MuonIDs, ElectronIDs, JetIDs;
	int nEvent;

	void executeEvent();
	vector<Muon> muons, muons_tight, muons_veto;
	vector<Electron> electrons, electrons_tight, electrons_veto;
	vector<Jet> jets, jets_tight, jets_dR04, bjets_tight, bjets_dR04;
	vector<Gen> gens;
	Particle METv;
	double w_prefire;
	
	void executeEventFromParameter(AnalyzerParameter param);
	LeptonType();
	~LeptonType();

	// clearing vector containers for each event
	void ClearPreColl();
	void ClearSubColl();

	// history for Gen
	vector<Gen_t> RecordHistory();
	vector<vector<unsigned int>> ReturnHistory(const vector<Gen_t>& Gens) const;
	void PrintHistory();

	// singal region selection
	bool isSignal(const Signal& sig);
	bool isControl(const Control& con);

	// DrawHists
	void DrawHists(TString path, const vector<Muon>& muons, const double& weight, const TString& muonid);
	void DrawHists(TString path, const vector<Electron>& electrons, const double& weight, const TString& elecid);
	void DrawHists(TString path, const vector<Jet>& jets, const double& weight, const JetType& type, const TString& jetid);
	void DrawHists(TString path, const Particle& METv, const double& weight);
	void PrintGens();

};



#endif

