#ifndef IDVariables_h
#define IDVariables_h

#include "AnalyzerCore.h"

class IDVariables : public AnalyzerCore {

public:
  void initializeAnalyzer();
	
  void executeEvent();
	vector<Muon> AllMuons;
	vector<Electron> AllElectrons;
	vector<Gen> GenParts;

	TString MuonID, ElectronID;

  IDVariables();
  ~IDVariables();
};



#endif

