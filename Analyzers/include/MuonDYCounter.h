#ifndef MuonDYCounter_h
#define MuonDYCounter_h

#include "AnalyzerCore.h"

class MuonDYCounter : public AnalyzerCore {

public:

  void initializeAnalyzer();
  void executeEvent();

  MuonDYCounter();
  ~MuonDYCounter();

};



#endif

