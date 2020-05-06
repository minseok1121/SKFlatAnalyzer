#ifndef FakeEstimator_h
#define FakeEstimator_h

#include "AnalyzerCore.h"

class FakeEstimator : public AnalyzerCore {

public:

  void initializeAnalyzer();
  void executeEventFromParameter(AnalyzerParameter param);
  void executeEvent();

  FakeEstimator();
  ~FakeEstimator();

  //==== initializeAnalyzer
  bool RunSysts, RunNormSysts, RunXsecSyst;
  
  vector<TString> Systs;
  vector<TString> ElectronIDs;

  TString HLTElecTriggerName;
  double TriggerSafePtCut;

  //==== executeEvent
  vector<Muon> AllMuons;
  vector<Jet> AllJets;
  vector<Electron> AllElectrons;
  vector<Gen> AllGens;

  TString MuonID, ElectronID, JetID;

  double weight_Prefire, weight_PileUp;
  
  //==== executeEventFromParameter
  vector<Electron> electrons, electrons_loose;
  vector<Muon> muons;
  vector<Jet> jets, clean_jets04, clean_jets10;
};



#endif

