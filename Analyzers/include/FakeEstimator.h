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
  bool RunSysts, RunXsecSyst;
  vector<TString> Systs; 
  vector<TString> ElectronIDs;
  TFile* f_nPV;

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

  //==== private functions
  map< TString, map< TString, TH1D*> > maphist_TH1D;
  map< TString, map< TString, TH2D*> > maphist_TH2D;
  map< TString, map< TString, TDirectory*> > mapDirectory;
  double GetNPVReweight(TString id, TString syst);
  double GetCorrPt(Electron e);
  TH1D* GetHist1D(TString suffix, TString histname);
  TH2D* GetHist2D(TString suffix, TString histname);
  void FillHist(TString syst, TString histname, double value, double weight, int n_bin, double x_min, double x_max);
  void FillHist(TString syst, TString histname, double value_x, double value_y, double weight,
												double n_binx, double* xbins, 
												double n_biny, double* ybins);
  void WriteHist();
};



#endif

