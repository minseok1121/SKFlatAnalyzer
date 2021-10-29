#ifndef TutorialRun_h
#define TutorialRun_h

#include "AnalyzerCore.h"

class TutorialRun : public AnalyzerCore {

public:

  void initializeAnalyzer();
  void executeEvent();

  TString IsoMuTriggerName;
	double TriggerSafePtCut;

	vector<TString> MuonIDs, MuonIDSFKeys;
	TString MuonVetoID, ElectronVetoID;
	vector<Muon> AllMuons;
	vector<Electron> AllElectrons;
	vector<Jet> AllJets;

  TutorialRun();
  ~TutorialRun();

};



#endif

