#include "Selector.h"

Selector::Selector() {
		// Link Tree Contents
		Events = new TTree("Events", "Events");
		// events
		Events->Branch("IsDATA", &IsDATA);
		Events->Branch("DataStream", &DataStream);
    Events->Branch("run", &run); 
		Events->Branch("event", &event); 
		Events->Branch("lumi", &lumi);
		Events->Branch("TrigLumi", &TrigLumi);
		Events->Branch("GenWeight", &GenWeight);
		Events->Branch("L1PrefireWeight", &L1PrefireWeight);
		Events->Branch("L1PrefireWeightUp", &L1PrefireWeightUp);
		Events->Branch("L1PrefireWeightDown", &L1PrefireWeightDown);
    Events->Branch("nPV", &nPV); 
		Events->Branch("nPileUp", &nPileUp);
    Events->Branch("PassDblMuTrigs", &PassDblMuTrigs); 
		Events->Branch("PassEMuTrigs", &PassEMuTrigs);
    Events->Branch("METvPt", &METvPt);
    Events->Branch("METvPhi", &METvPhi);
    
		// muons
    Events->Branch("nMuons", &nMuons);
		Events->Branch("MuonPtColl", MuonPtColl, "MuonPtColl[nMuons]/F");
		Events->Branch("MuonPtColl_MomentumShiftUp", MuonPtColl_MomentumShiftUp, "MuonPtColl_MomentumShiftUp[nMuons]/F");
		Events->Branch("MuonPtColl_MomentumShiftDown", MuonPtColl_MomentumShiftDown, "MuonPtColl_MomentumShiftDown[nMuons]/F");
		Events->Branch("MuonEtaColl", MuonEtaColl, "MuonEtaColl[nMuons]/F");
    Events->Branch("MuonPhiColl", MuonPhiColl, "MuonPhiColl[nMuons]/F");
    Events->Branch("MuonMassColl", MuonMassColl, "MuonMassColl[nMuons]/F");
		Events->Branch("MuonRelIsoColl", MuonRelIsoColl, "MuonRelIsoColl[nMuons]/F");
    Events->Branch("MuonMiniRelIsoColl", MuonMiniRelIsoColl, "MuonMiniRelIsoColl[nMuons]/F");
    Events->Branch("MuonChargeColl", MuonChargeColl, "MuonChargeColl[nMuons]/I");
    Events->Branch("MuonLepTypeColl", MuonLepTypeColl, "MuonLepTypeColl[nMuons]/I");
    Events->Branch("MuonPassTightColl", MuonPassTightColl, "MuonPassTightColl[nMuons]/O");
		Events->Branch("MuonPassLooseColl", MuonPassLooseColl, "MuonPassLooseColl[nMuons]/O");

    // electrons
		Events->Branch("nElectrons", &nElectrons);
    Events->Branch("ElectronPtColl", ElectronPtColl, "ElectronPtColl[nElectrons]/F");
    Events->Branch("ElectronPtColl_EnShiftUp", ElectronPtColl_EnShiftUp, "ElectronPtColl_EnShiftUp[nElectrons]/F");
		Events->Branch("ElectronPtColl_EnShiftDown", ElectronPtColl_EnShiftDown, "ElectronPtColl_EnShiftDown[nElectrons]/F");
		Events->Branch("ElectronPtColl_ResShiftUp", ElectronPtColl_ResShiftUp, "ElectronPtColl_ResShiftUp[nElectrons]/F");
		Events->Branch("ElectronPtColl_ResShiftDown", ElectronPtColl_ResShiftDown, "ElectronPtColl_ResShiftDown[nElectrons]/F");
		Events->Branch("ElectronEtaColl", ElectronEtaColl, "ElectronEtaColl[nElectrons]/F");
    Events->Branch("ElectronPhiColl", ElectronPhiColl, "ElectronPhiColl[nElectrons]/F");
    Events->Branch("ElectronMassColl", ElectronMassColl, "ElectronMassColl[nElectrons]/F");
    Events->Branch("ElectronRelIsoColl", ElectronRelIsoColl, "ElectronRelIsoColl[nElectrons]/F");
    Events->Branch("ElectronMiniRelIsoColl", ElectronMiniRelIsoColl, "ElectronMiniRelIsoColl[nElectrons]/F");
    Events->Branch("ElectronChargeColl", ElectronChargeColl, "ElectronChargeColl[nElectrons]/I");
    Events->Branch("ElectronLepTypeColl", ElectronLepTypeColl, "ElectronLepTypeColl[nElectrons]/I");
    Events->Branch("ElectronPassTightColl", ElectronPassTightColl, "ElectronPassTightColl[nElectrons]/O");
    Events->Branch("ElectronPassLooseColl", ElectronPassLooseColl, "ElectronPassLooseColl[nElectrons]/O");
    
		// jets
		Events->Branch("nJets", &nJets);
    Events->Branch("JetPtColl", JetPtColl, "JetPtColl[nJets]/F");
    Events->Branch("JetPtColl_EnShiftUp", JetPtColl_EnShiftUp, "JetPtColl_EnShiftUp[nJets]/F");
    Events->Branch("JetPtColl_EnShiftDown", JetPtColl_EnShiftDown, "JetPtColl_EnShiftDown[nJets]/F");
    Events->Branch("JetPtColl_ResShiftUp", JetPtColl_ResShiftUp, "JetPtColl_ResShiftUp[nJets]/F");
    Events->Branch("JetPtColl_ResShiftDown", JetPtColl_ResShiftDown, "JetPtColl_ResShiftDown[nJets]/F");
    Events->Branch("JetEtaColl", JetEtaColl, "JetEtaColl[nJets]/F");
    Events->Branch("JetPhiColl", JetPhiColl, "JetPhiColl[nJets]/F");
    Events->Branch("JetMassColl", JetMassColl, "JetMassColl[nJets]/F");
		Events->Branch("JetChargeColl", JetChargeColl, "JetChargeColl[nJets]/I");
		Events->Branch("JetPartonFlavourColl", JetPartonFlavourColl, "JetPartonFlavourColl[nJets]/I");
		Events->Branch("JetHadronFlavourColl", JetHadronFlavourColl, "JetHadronFlavourColl[nJets]/I");
		Events->Branch("JetBtagScoreColl", JetBtagScoreColl, "JetBtagScoreColl[nJets]/F");
		Events->Branch("JetIsBtaggedColl", JetIsBtaggedColl, "JetIsBtaggedColl[nJets]/O");

}

Selector::~Selector() {
		outfile->cd();
		Events->Write();
}

void Selector::initializeAnalyzer(){
		// flags
		Skim1E2Mu = HasFlag("Skim1E2Mu");
		Skim3Mu = HasFlag("Skim3Mu");

		// trigger & ID settings
		if (DataEra == "2016preVFP") {
        DblMuTriggers = {
            "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v",
            "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v",
        };
        EMuTriggers = {
            "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v",
            "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",
            "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v",
            "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v"
        };
        MuonIDs = {"HcToWATight", "HcToWALoose", "HcToWAVeto"};
        ElectronIDs = {"HcToWATight16a", "HcToWALoose16a", "HcToWAVeto16a"};
    }
    else if (DataEra == "2016postVFP") {
        DblMuTriggers = {
            "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v",
            "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v",
        };
        EMuTriggers = {
            "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v",
            "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",
            "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v",
            "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v"
        };
        MuonIDs = {"HcToWATight", "HcToWALoose", "HcToWAVeto"};
        ElectronIDs = {"HcToWATight16b", "HcToWALoose16b", "HcToWAVeto16b"};
    }
		else if (DataEra == "2017") {
        DblMuTriggers = {
            "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v",
            "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v",
            "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v"
        };
        EMuTriggers = {
            "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v",
            "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v"
        };
        MuonIDs = {"HcToWATight", "HcToWALoose", "HcToWAVeto"};
        ElectronIDs = {"HcToWATight17", "HcToWALoose17", "HcToWAVeto17"};
    }
    else if (DataEra == "2018") {
        DblMuTriggers = {
            "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v",
            "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v",
            "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v"
        };
        EMuTriggers = {
            "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v",
            "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",
            "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v"
        };
        MuonIDs = {"HcToWATight", "HcToWALoose", "HcToWAVeto"};
        ElectronIDs = {"HcToWATight18", "HcToWALoose18", "HcToWAVeto18"};
    }
    else {
        cerr << "[diLepControlRegion::initializeAnalyzer] Wrong era " << DataEra << endl;
        exit(EXIT_FAILURE);
    }

		// Jet tagger
		vector<JetTagging::Parameters> jtps;
    jtps.emplace_back(
            JetTagging::Parameters(JetTagging::DeepJet, JetTagging::Medium, JetTagging::incl, JetTagging::mujets));
    mcCorr->SetJetTaggingParameters(jtps);
}

void Selector::executeEvent(){
		
		if (!PassMETFilter()) return;
		Event ev = GetEvent();
		const Particle METv = ev.GetMETVector();

		// Object Definitions
		vector<Gen> truth_coll = GetGens();
		vector<Muon> muon_coll = GetAllMuons();
		vector<Electron> ele_coll = GetAllElectrons();
		vector<Jet> jet_coll = GetAllJets();

		// sort the objects in the first time
		sort(muon_coll.begin(), muon_coll.end(), PtComparing);
		sort(ele_coll.begin(), ele_coll.end(), PtComparing);
		sort(jet_coll.begin(), jet_coll.end(), PtComparing);

		vector<Muon> muonV_coll = SelectMuons(muon_coll, MuonIDs.at(2), 10., 2.4);
		vector<Electron> eleV_coll = SelectElectrons(ele_coll, ElectronIDs.at(2), 10., 2.4);
		vector<Jet> jetT_coll = SelectJets(jet_coll, "tight", 20., 2.4);
		jetT_coll = JetsVetoLeptonInside(jetT_coll, eleV_coll, muonV_coll, 0.4);

		// Channel
		if (Skim3Mu){
				if (! (muonV_coll.size() == 3 && eleV_coll.size() == 0)) return;
		}
		else if (Skim1E2Mu) {
				if (! (muonV_coll.size() == 2 && eleV_coll.size() == 1)) return;
		}
		else {
				cerr << "Wrong flag" << endl;
				exit(EXIT_FAILURE);
		}
		
		// Initialize contents
		// events
		PassDblMuTrigs = ev.PassTrigger(DblMuTriggers);
		PassEMuTrigs = ev.PassTrigger(EMuTriggers);
		METvPt = METv.Pt(); METvPhi = METv.Phi();
		GenWeight = 1.; TrigLumi = 1.;
		if (!IsDATA) {
				GenWeight = MCweight();
				TrigLumi = ev.GetTriggerLumi("Full");
		}
		L1PrefireWeight = GetPrefireWeight(0);
		L1PrefireWeightUp = GetPrefireWeight(1);
		L1PrefireWeightDown = GetPrefireWeight(-1);

		// muons
		nMuons = muonV_coll.size();
		for (unsigned int i = 0; i < nMuons; i++) {
				const Muon &mu = muonV_coll.at(i);
				MuonPtColl[i] = mu.Pt(); 
				MuonPtColl_MomentumShiftUp[i] = mu.MomentumShift(1);
				MuonPtColl_MomentumShiftDown[i] = mu.MomentumShift(-1);
				MuonEtaColl[i] = mu.Eta(); 
				MuonPhiColl[i] = mu.Phi();
				MuonMassColl[i] = mu.M();
				MuonRelIsoColl[i] = mu.RelIso();
				MuonMiniRelIsoColl[i] = mu.MiniRelIso();
				MuonChargeColl[i] = mu.Charge();
				MuonLepTypeColl[i] = GetLeptonType(mu, truth_coll);
				MuonPassTightColl[i] = mu.PassID(MuonIDs.at(0));
				MuonPassLooseColl[i] = mu.PassID(MuonIDs.at(1));
		}
		// electrons
		nElectrons = eleV_coll.size();
		for (unsigned int i = 0; i < nElectrons; i++) {
        const Electron &ele = eleV_coll.at(i);
        ElectronPtColl[i] = ele.Pt();
				ElectronPtColl_EnShiftUp[i] = ele.Pt() * ele.EnShift(1);
				ElectronPtColl_EnShiftDown[i] = ele.Pt() * ele.EnShift(-1);
				ElectronPtColl_ResShiftUp[i] = ele.Pt() * ele.ResShift(1);
				ElectronPtColl_ResShiftDown[i] = ele.Pt() * ele.ResShift(-1);
        ElectronEtaColl[i] = ele.Eta();
        ElectronPhiColl[i] = ele.Phi();
        ElectronMassColl[i] = ele.M();
				ElectronRelIsoColl[i] = ele.RelIso();
        ElectronMiniRelIsoColl[i] = ele.MiniRelIso();
				ElectronChargeColl[i] = ele.Charge();
        ElectronLepTypeColl[i] = GetLeptonType(ele, truth_coll);
        ElectronPassTightColl[i] = ele.PassID(ElectronIDs.at(0));
				ElectronPassLooseColl[i] = ele.PassID(ElectronIDs.at(1));
    }
		// jets
		nJets = jetT_coll.size();
		const float bcut = mcCorr->GetJetTaggingCutValue(JetTagging::DeepJet, JetTagging::Medium);
		for (unsigned int i = 0; i < nJets; i++) {
        const Jet& j = jetT_coll.at(i);
        JetPtColl[i] = j.Pt();
				JetPtColl_EnShiftUp[i] = j.Pt() * j.EnShift(1);
				JetPtColl_EnShiftDown[i] = j.Pt() * j.EnShift(-1);
				JetPtColl_ResShiftUp[i] = j.Pt() * j.ResShift(1);
				JetPtColl_ResShiftDown[i] = j.Pt() * j.ResShift(-1);
        JetEtaColl[i] = j.Eta();
        JetPhiColl[i] = j.Phi();
        JetMassColl[i] = j.M();
				JetChargeColl[i] = j.Charge();
				JetPartonFlavourColl[i] = j.partonFlavour();
				JetHadronFlavourColl[i] = j.hadronFlavour();
        JetBtagScoreColl[i] = j.GetTaggerResult(JetTagging::DeepJet);
        JetIsBtaggedColl[i] = j.GetTaggerResult(JetTagging::DeepJet) > bcut;
    }

    Events->Fill();
}

