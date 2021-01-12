#ifndef TestAnalyzer_h
#define TestAnalyzer_h

#include "AnalyzerCore.h"

class TestAnalyzer : public AnalyzerCore {
private:
	int __called;
public:
	//==== CONSTURCTOR & DESTRUCTOR ====
	TestAnalyzer();
	~TestAnalyzer();

	//==== ALIASES ====
	//==== INITIALIZE ====
	void initializeAnalyzer();
	bool RunSyst;
	vector<TString> DblMuTriggers, EMuTriggers;
	vector<double> TriggerSafePtCuts;

	//==== EXECUTE ====
	void executeEventFromParameter(AnalyzerParameter param);
	void executeEvent();

	//==== OTHER MEMBER FUNCTIONS ====
	void MyCutflowMaker();
	void MyHistoMaker(TString path, const vector<Muon> &muons, const double &weight);
	void MyHistoMaker(TString path, const vector<Electron> &electrons, const double &weight);
	void MyHistoMaker(TString path, const vector<Jet> &jets, const double &weight);
	void MyHistoMaker(TString path, const Particle &METv, const double &weight);

};



#endif

