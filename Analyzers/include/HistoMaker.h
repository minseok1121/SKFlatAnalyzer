#ifndef HistoMaker_h__
#define HistoMaker_h__

#include "AnalyzerCore.h"

class HistoMaker: public AnalyzerCore {
private:
	TString _region;
	vector<TString> _cuts;
public:
	HistoMaker();
	HistoMaker(const TString& region, const vector<TString> &cuts);
	~HistoMaker();
	void SetRegion(const TString &region);
	TString GetRegion() const;
	// basically it will be stored like region/path
	void FillCutflow(const TString &cut);
	void FillObjects(const TString path, const vector<Muon> &muons, const double &weight);
	void FillObjects(const TString path, const vector<Electron> &electrons, const double &weight);
	void FillObjects(const TString path, const vector<Jet> &jets, const double &weight);
	void FillObjects(const TString path, const Particle &part, const double &weight);
};

#endif
