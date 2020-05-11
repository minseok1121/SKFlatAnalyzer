#ifndef FakeEstimator_BtagDep_h
#define FakeEstimator_BtagDep_h

#include "AnalyzerCore.h"

class FakeEstimator_BtagDep : public AnalyzerCore {

public:

  void initializeAnalyzer();
  void executeEventFromParameter(AnalyzerParameter param);
  void executeEvent();

  FakeEstimator_BtagDep();
  ~FakeEstimator_BtagDep();

};



#endif

