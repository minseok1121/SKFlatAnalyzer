#ifndef SampleValidation_h
#define SampleValidation_h

#include "AnalyzerCore.h"

class SampleValidation : public AnalyzerCore {

public:
		SampleValidation();
		~SampleValidation();

		void initializeAnalyzer();
		void executeEvent();

};

#endif

