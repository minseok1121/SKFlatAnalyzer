#ifndef IDVariables_h
#define IDVariables_h

#include "AnalyzerCore.h"

class IDVariables : public AnalyzerCore {

public:
  void initializeAnalyzer();
	vector<TString> MuonIDs;
	vector<TString> ElectronIDs;
	
  void executeEvent();
	vector<Muon> AllMuons;
	vector<Electron> AllElectrons;
	vector<Gen> GenParts;

  IDVariables();
  ~IDVariables();
};



#endif

