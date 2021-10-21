#ifndef ElectronTightID_h
#define ElectronTightID_h

#include "AnalyzerCore.h"

class ElectronTightID : public AnalyzerCore {

public:

  void initializeAnalyzer();
  void executeEvent();

  ElectronTightID();
  ~ElectronTightID();

};



#endif

