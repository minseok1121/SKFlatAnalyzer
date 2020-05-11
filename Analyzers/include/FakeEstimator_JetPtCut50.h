#ifndef FakeEstimator_JetPtCut50_h
#define FakeEstimator_JetPtCut50_h

#include "AnalyzerCore.h"

class FakeEstimator_JetPtCut50 : public AnalyzerCore {

public:

  void initializeAnalyzer();
  void executeEventFromParameter(AnalyzerParameter param);
  void executeEvent();

  FakeEstimator_JetPtCut50();
  ~FakeEstimator_JetPtCut50();

};



#endif

