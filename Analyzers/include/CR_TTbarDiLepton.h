#ifndef CR_TTbarDiLepton_h
#define CR_TTbarDiLepton_h

#include "AnalyzerCore.h"

class CR_TTbarDiLepton : public AnalyzerCore {

public:
    bool TTDiMu, TTEMu;     // flags
    vector<TString> ElectronIDs, MuonIDs;
    vector<TString> DblMuTriggers, EMuTriggers;
    vector<JetTagging::Parameters> jtps;

    void initializeAnalyzer();
    void executeEvent();

    CR_TTbarDiLepton();
    ~CR_TTbarDiLepton();
};



#endif

