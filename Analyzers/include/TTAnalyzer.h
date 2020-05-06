#ifndef TTAnalyzer_h
#define TTAnalyzer_h

#include "AnalyzerCore.h"

class TTAnalyzer : public AnalyzerCore {

public:

  void initializeAnalyzer();
  void executeEventFromParameter(AnalyzerParameter param);
  void executeEvent();

  bool RunSyst, RunNewPDF, RunXSecSyst;

  vector<TString> MuonIDs, MuonIDSFKeys, MuonIsoSFKeys;
  vector<TString> ElectronIDs, ElectronIDSFKeys;

  TString IsoMuTriggerName;
  double TriggerSafePtCut;

  vector<Muon> AllMuons;
  vector<Jet> AllJets;
  vector<Electron> AllElectrons;
  vector<Gen> AllGens;

  double weight_Prefire, weight_PileUp, weight_TopPt;

  TTAnalyzer();
  ~TTAnalyzer();

};



#endif

