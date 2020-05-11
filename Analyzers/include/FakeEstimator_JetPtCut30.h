#ifndef FakeEstimator_JetPtCut30_h
#define FakeEstimator_JetPtCut30_h

#include "AnalyzerCore.h"

class FakeEstimator_JetPtCut30 : public AnalyzerCore {

public:

  void initializeAnalyzer();
  void executeEventFromParameter(AnalyzerParameter param);
  void executeEvent();

  FakeEstimator_JetPtCut30();
  ~FakeEstimator_JetPtCut30();

};



#endif

