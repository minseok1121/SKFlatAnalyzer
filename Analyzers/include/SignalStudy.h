#ifndef SignalStudy_h
#define SignalStudy_h

#include "AnalyzerCore.h"

class SignalStudy : public AnalyzerCore {
private:
	int __called;
public:
	// Aliases
	typedef enum { NONE, SR3MU, SR1E2MU } SIGNAL;
	// Constructor and Destructor
	SignalStudy();
	~SignalStudy();

	// Aliases
	// Initialization
	void initializeAnalyzer();

	// execute
	void executeEvent();
	SIGNAL SignalSelector(
			vector<Muon> &muons_passID, vector<Electron> &electrons_passID, 
			vector<Muon> &muons_veto, vector<Electron> &electrons_veto);

	// other member functions
	void myCutflowMaker();
	void myHistoMaker(TString path, const vector<Muon> &muons, const double &weight);
	void myHistoMaker(TString path, const vector<Electron> &electrons, const double &weight);
	void myHistoMaker(TString path, const vector<Jet> &jets, const double &weight);
	void myHistoMaker(TString path, const Particle &part, const double &weight);
};



#endif

