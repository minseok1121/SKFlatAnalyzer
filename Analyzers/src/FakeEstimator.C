#include "FakeEstimator.h"

void FakeEstimator::initializeAnalyzer(){

}

void FakeEstimator::executeEvent(){


  AnalyzerParameter param;

  executeEventFromParameter(param);

}

void FakeEstimator::executeEventFromParameter(AnalyzerParameter param){

  if(!PassMETFilter()) return;

  Event ev = GetEvent();

}

FakeEstimator::FakeEstimator(){

}

FakeEstimator::~FakeEstimator(){

}


