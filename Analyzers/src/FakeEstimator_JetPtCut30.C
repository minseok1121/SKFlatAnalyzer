#include "FakeEstimator_JetPtCut30.h"

void FakeEstimator_JetPtCut30::initializeAnalyzer(){

}

void FakeEstimator_JetPtCut30::executeEvent(){


  AnalyzerParameter param;

  executeEventFromParameter(param);

}

void FakeEstimator_JetPtCut30::executeEventFromParameter(AnalyzerParameter param){

  if(!PassMETFilter()) return;

  Event ev = GetEvent();

}

FakeEstimator_JetPtCut30::FakeEstimator_JetPtCut30(){

}

FakeEstimator_JetPtCut30::~FakeEstimator_JetPtCut30(){

}


