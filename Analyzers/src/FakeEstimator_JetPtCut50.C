#include "FakeEstimator_JetPtCut50.h"

void FakeEstimator_JetPtCut50::initializeAnalyzer(){

}

void FakeEstimator_JetPtCut50::executeEvent(){


  AnalyzerParameter param;

  executeEventFromParameter(param);

}

void FakeEstimator_JetPtCut50::executeEventFromParameter(AnalyzerParameter param){

  if(!PassMETFilter()) return;

  Event ev = GetEvent();

}

FakeEstimator_JetPtCut50::FakeEstimator_JetPtCut50(){

}

FakeEstimator_JetPtCut50::~FakeEstimator_JetPtCut50(){

}


