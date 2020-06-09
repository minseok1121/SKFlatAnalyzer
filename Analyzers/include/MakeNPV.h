#ifndef MakeNPV_h
#define MakeNPV_h

#include "AnalyzerCore.h"

class MakeNPV : public AnalyzerCore {

public:

  void initializeAnalyzer();
  void executeEventFromParameter(AnalyzerParameter param);
  void executeEvent();

  MakeNPV();
  ~MakeNPV();

  //==== initializeAnalyzer
  bool RunSysts, RunXsecSyst;
  
  vector<TString> ElectronIDs;

  TString HLTElecTriggerName, HLTElecTriggerName1, HLTElecTriggerName2;
  double TriggerSafePtCut, TriggerSafePtCut1, TriggerSafePtCut2;

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
  vector<Jet> jets, clean_jets04;

};



#endif

