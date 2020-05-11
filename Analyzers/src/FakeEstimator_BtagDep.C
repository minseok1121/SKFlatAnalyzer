#include "FakeEstimator_BtagDep.h"

void FakeEstimator_BtagDep::initializeAnalyzer(){

}

void FakeEstimator_BtagDep::executeEvent(){


  AnalyzerParameter param;

  executeEventFromParameter(param);

}

void FakeEstimator_BtagDep::executeEventFromParameter(AnalyzerParameter param){

  if(!PassMETFilter()) return;

  Event ev = GetEvent();

}

FakeEstimator_BtagDep::FakeEstimator_BtagDep(){

}

FakeEstimator_BtagDep::~FakeEstimator_BtagDep(){

}


