#ifndef diLepRegion_h
#define diLepRegion_h

#include "AnalyzerCore.h"

class diLepRegion : public AnalyzerCore {

public:
		// constructor & destructor
		diLepRegion();
		~diLepRegion();

		// initialization
		bool RunLowPt, RunSysts;
		void initializeAnalyzer();
		vector<TString> Systs;
		TH2D* h_muon_sf;
		TH2D* h_ele_sf;

		vector<TString> DblMuTrigs, EMuTrigs;
		vector<TString> MuonIDs, ElectronIDs;

		// execution
		void executeEvent();
		TString selectEvent(
				Event &ev, vector<Muon> &muons_tight, vector<Electron> &electrons_tight,
				vector<Muon> &muons_loose, vector<Electron> &electrons_loose,
				vector<Jet> &jets, vector<Jet> &bjets, Particle &METv);
		double getWeight(Event &ev, 
										 vector<Muon> &muons, vector<Electron> &electrons, vector<Jet> &jets,
										 TString syst);
		double getMuonIDSF(const Muon &mu, int sys);
		double getElectronIDSF(const Electron &ele, int sys);
};



#endif

