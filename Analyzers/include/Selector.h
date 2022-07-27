#ifndef Selector_h
#define Selector_h

#include "AnalyzerCore.h"

class Selector : public AnalyzerCore {
private:
		bool Skim1E2Mu, Skim3Mu;
		vector<TString> ElectronIDs, MuonIDs;
		vector<TString> DblMuTriggers, EMuTriggers;
		TTree* Events;

		// events
		bool passDblMuTrigs, passEMuTrigs;
		float trigLumi;
		float genWeight;
		float METv_pt, METv_phi;
		// muons
		unsigned int nMuons;
		float muons_pt[3], muons_eta[3], muons_phi[3], muons_mass[3], muons_relIso[3], muons_miniIso[3];
		int muons_charge[3], muons_lepType[3];
		bool muons_passTight[3], muons_passLoose[3];
		// electrons
		unsigned int nElectrons;
		float electrons_pt[1], electrons_eta[1], electrons_phi[1], electrons_mass[1], electrons_relIso[3], electrons_miniIso[3];
		int electrons_charge[1], electrons_lepType[1];
		bool electrons_passTight[1], electrons_passLoose[1];
		// jets
		unsigned int nJets;
		float jets_pt[30], jets_eta[30], jets_phi[30], jets_mass[30], jets_btagScore[30];
		bool jets_isBtagged[30];

public:
		Selector();
		~Selector();
		void initializeAnalyzer();
		void executeEvent();

};



#endif

