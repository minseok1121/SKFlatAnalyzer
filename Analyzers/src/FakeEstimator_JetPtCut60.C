#include "FakeEstimator_JetPtCut60.h"

void FakeEstimator_JetPtCut60::initializeAnalyzer(){

}

void FakeEstimator_JetPtCut60::executeEvent(){


  AnalyzerParameter param;

  executeEventFromParameter(param);

}

void FakeEstimator_JetPtCut60::executeEventFromParameter(AnalyzerParameter param){

  if(!PassMETFilter()) return;

  Event ev = GetEvent();

}

FakeEstimator_JetPtCut60::FakeEstimator_JetPtCut60(){

}

FakeEstimator_JetPtCut60::~FakeEstimator_JetPtCut60(){

}


