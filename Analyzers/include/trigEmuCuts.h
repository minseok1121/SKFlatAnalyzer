#ifndef trigEmuCuts_h
#define trigEmuCuts_h

#include "AnalyzerCore.h"

class trigEmuCuts : public AnalyzerCore {

public:
		// constructor & destructor
		trigEmuCuts();
		~trigEmuCuts();

		// initialization
		void initializeAnalyzer();
		vector<TString> Ele12Trigs, Ele23Trigs;
		 vector<TString> AllEMuTrigs;

		// execution
		void executeEvent();
		TString selectEvent(
						Event &ev, vector<Muon> &muons_tight, vector<Muon> &muons_loose, 
						vector<Electron> &electrons, vector<Jet> &jets, vector<Jet> &bjets,
						Particle &METv);

};



#endif

