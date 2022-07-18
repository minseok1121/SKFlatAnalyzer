#ifndef diLepControlRegion_h
#define diLepControlRegion_h

#include "AnalyzerCore.h"

class diLepControlRegion : public AnalyzerCore {

public:
		typedef TString channel_t;
		typedef TString syst_t;
		// constructor & destructor
		diLepControlRegion();
		~diLepControlRegion();

		// initialization
		void initializeAnalyzer();
		vector<TString> DblMuTriggers, EMuTriggers;
		vector<TString> MuonIDs, ElectronIDs;
		vector<syst_t> Systematics;

		// execution
		void executeEvent();
		channel_t selectEvent(Event &ev,
												 vector<Muon> &muonT_coll,
												 vector<Muon> &muonL_coll,
												 vector<Muon> &muonV_coll,
												 vector<Electron> &eleT_coll,
												 vector<Electron> &eleL_coll,
												 vector<Electron> &eleV_coll,
												 vector<Jet> &jetT_coll,
												 vector<Jet> &bjet_coll,
												 Particle &METv);
		double getWeight(Event &ev,
										 vector<Muon> &muonT_coll,
										 vector<Muon> &muonL_coll,
										 vector<Muon> &muonV_coll,
										 vector<Electron> &eleT_coll,
										 vector<Electron> &eleL_coll,
										 vector<Electron> &eleV_coll,
										 vector<Jet> &jetT_coll,
										 vector<Jet> &bjet_coll,
										 Particle &METv,
										 syst_t syst);
};



#endif

