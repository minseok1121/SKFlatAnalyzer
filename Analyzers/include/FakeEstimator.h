#ifndef FakeEstimator_h
#define FakeEstimator_h

#include "AnalyzerCore.h"

class FakeEstimator : public AnalyzerCore {

public:
  vector<TString> IDs;
  vector<TString> Systs;
  vector<TString> Prompts;
  vector<TString> Regions;

  //==== initializeAnalyzer
  void initializeAnalyzer();
  bool RunSysts, RunXsecSyst;
  vector<TString> ElectronIDs = IDs;
  TFile* f_nPV;
  TString HLTElecTriggerName, HLTElecTriggerName1, HLTElecTriggerName2;
  double TriggerSafePtCut, TriggerSafePtCut1, TriggerSafePtCut2;
  double JetPtCut;

  //==== executeEvent
  void executeEvent();
  vector<Muon> AllMuons;
  vector<Jet> AllJets;
  vector<Electron> AllElectrons;
  vector<Gen> AllGens;

  TString MuonID, ElectronID, JetID, syst, prompt;
  double weight_Prefire, weight_PileUp;

  //==== executeEventFromParameter
  void executeEventFromParameter(AnalyzerParameter param);
  vector<Electron> electrons, electrons_loose;
  vector<Muon> muons;
  vector<Jet> jets, clean_jets04, clean_jets10;

  //==== Member functions
  double GetJetPtCut(TString syst);
  map<TString, TH1D*> maphist_TH1D;
  map<TString, TH2D*> maphist_Th2D;
  map<TString,
	  map<TString,
	      map<TString,
			map<TString, TDirectory*>>>> mapDirectory;
  double GetNPVReweight(TString id, TString syst);
  double GetCorrPt(Electron e);
  TH1D* GetHist1D(TString histname);
  TH2D* GetHist2D(TString histname);
  void FillHist(TString histname, double value, double weight, int n_bin, double x_min, double x_max);
  void FillHist(TString histname, double value_x, double value_y, double weight,
				double n_binx, double* xbins,
				double n_biny, double* ybins);
  void WriteHist();

  FakeEstimator();
  ~FakeEstimator();

};



#endif

