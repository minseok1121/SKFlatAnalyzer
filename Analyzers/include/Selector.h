#ifndef Selector_h
#define Selector_h

#include "AnalyzerCore.h"

class Selector : public AnalyzerCore {
private:
		bool Skim1E2Mu, Skim3Mu;
		TString TightElecID, LooseElecID;
		TString TightMuonID, LooseMuonID;
		vector<TString> trigs_dblmu, trigs_emu;
		TTree* tree;
		// events
		bool passDblMuTrigs, passEMuTrigs;
		float trigLumi;
		float genWeight;
		float METv_pt, METv_phi;
		// muons
		unsigned int nMuons;
		float muons_pt[3], muons_eta[3], muons_phi[3], muons_mass[3], muons_miniIso[3];
		int muons_charge[3], muons_lepType[3];
		bool muons_isTight[3];
		// electrons
		unsigned int nElectrons;
		float electrons_pt[1], electrons_eta[1], electrons_phi[1], electrons_mass[1], electrons_miniIso[1];
		int electrons_charge[1], electrons_lepType[1];
		bool electrons_isTight[1];

		// jets
		unsigned int nJets;
		float jets_pt[30], jets_eta[30], jets_phi[30], jets_mass[30], jets_btagScore[30];
		float jets_cH[30], jets_nH[30], jets_cEM[30], jets_nEM[30], jets_muE[30];
		float jets_cM[30], jets_nM[30];
		bool jets_isBtagged[30];

public:
		Selector();
		~Selector();
		void SetMuonIDs(TString TightID, TString LooseID);
		void SetElectronIDs(TString TightID, TString LooseID);
		void SetTriggers(vector<TString> emuTrigs, vector<TString> dblmuTrigs);
		void LinkTreeContents();

		void initializeAnalyzer();
		void executeEvent();

};



#endif

