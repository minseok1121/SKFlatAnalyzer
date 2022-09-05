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
		bool PassDblMuTrigs, PassEMuTrigs;
		float TrigLumi;
		float GenWeight;
		float L1PrefireWeight, L1PrefireWeightUp, L1PrefireWeightDown;

		// METVector
		float METvPt, METvPhi;

		// muons
		unsigned int nMuons;
		float MuonPtColl[3];
		float MuonPtColl_MomentumShiftUp[3];
		float MuonPtColl_MomentumShiftDown[3];
		float MuonEtaColl[3];
		float MuonPhiColl[3];
		float MuonMassColl[3];
		float MuonRelIsoColl[3];
		float MuonMiniRelIsoColl[3];
		int   MuonChargeColl[3];
		int   MuonLepTypeColl[3];
		bool  MuonPassTightColl[3];
		bool  MuonPassLooseColl[3];

		// electrons
		unsigned int nElectrons;
		float ElectronPtColl[1];
		float ElectronPtColl_EnShiftUp[1];
		float ElectronPtColl_EnShiftDown[1];
		float ElectronPtColl_ResShiftUp[1];
		float ElectronPtColl_ResShiftDown[1];
		float ElectronEtaColl[1];
		float ElectronPhiColl[1];
		float ElectronMassColl[1];
		float ElectronRelIsoColl[1];
		float ElectronMiniRelIsoColl[1];
		int   ElectronChargeColl[1];
		int   ElectronLepTypeColl[1];
		bool  ElectronPassTightColl[1];
		bool  ElectronPassLooseColl[1];

		// jets
		unsigned int nJets;
		float JetPtColl[20];
		float JetPtColl_EnShiftUp[20];
		float JetPtColl_EnShiftDown[20];
		float JetPtColl_ResShiftUp[20];
		float JetPtColl_ResShiftDown[20];
		float JetEtaColl[20];
		float JetPhiColl[20];
		float JetMassColl[20];
		int   JetChargeColl[20];
		int   JetPartonFlavourColl[20];
		int   JetHadronFlavourColl[20];
		float JetBtagScoreColl[20];
		bool  JetIsBtaggedColl[20];

public:
		Selector();
		~Selector();
		void initializeAnalyzer();
		void executeEvent();

};



#endif

