#ifndef CR_DiLepton_h
#define CR_DiLepton_h

#include "AnalyzerCore.h"

class CR_DiLepton : public AnalyzerCore {

public:
    bool RunSyst;     // flags
    vector<TString> ElectronIDs, MuonIDs;
    vector<TString> DblMuTriggers, EMuTriggers;
    vector<JetTagging::Parameters> jtps;
    vector<TString> systematics;
    
    TH2D *hMuonIDSF;
    TH2D *hMu17Leg1_DATA, *hMu17Leg1_MC;
    TH2D *hMu8Leg2_DATA, *hMu8Leg2_MC;

    void initializeAnalyzer();
    void executeEvent();

    double getMuonIDSF(const Muon &mu, int sys);
    double getTriggerEff(const Muon &mu, TString histkey, bool isDataEff, int sys);
    double getDblMuTriggerEff(vector<Muon> &muons, bool isDATA, int sys);
    double getDblMuTriggerSF(vector<Muon> &muons, int sys);

    CR_DiLepton();
    ~CR_DiLepton();
};



#endif

