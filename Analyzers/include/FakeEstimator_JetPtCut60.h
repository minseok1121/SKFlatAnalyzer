#ifndef FakeEstimator_JetPtCut60_h
#define FakeEstimator_JetPtCut60_h

#include "AnalyzerCore.h"

class FakeEstimator_JetPtCut60 : public AnalyzerCore {

public:

  void initializeAnalyzer();
  void executeEventFromParameter(AnalyzerParameter param);
  void executeEvent();

  FakeEstimator_JetPtCut60();
  ~FakeEstimator_JetPtCut60();

};



#endif

