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
    Events->Branch("nPV", &nPV); Events->Branch("nPileUp", &nPileUp);
    Events->Branch("passDblMuTrigs", &passDblMuTrigs); 
		Events->Branch("passEMuTrigs", &passEMuTrigs);
    Events->Branch("trigLumi", &trigLumi);
    Events->Branch("METv_pt", &METv_pt);
    Events->Branch("METv_phi", &METv_phi);
    Events->Branch("genWeight", &genWeight);
    // muons
    Events->Branch("nMuons", &nMuons);
    Events->Branch("muons_pt", muons_pt, "muons_pt[nMuons]/F");
    Events->Branch("muons_eta", muons_eta, "muons_eta[nMuons]/F");
    Events->Branch("muons_phi", muons_phi, "muons_phi[nMuons]/F");
    Events->Branch("muons_mass", muons_mass, "muons_mass[nMuons]/F");
    Events->Branch("muons_miniIso", muons_miniIso, "muons_miniIso[nMuons]/F");
    Events->Branch("muons_charge", muons_charge, "muons_charge[nMuons]/I");
    Events->Branch("muons_lepType", muons_lepType, "muons_lepType[nMuons]/I");
    Events->Branch("muons_passTight", muons_passTight, "muons_passTight[nMuons]/O");
		Events->Branch("muons_passLoose", muons_passLoose, "muons_passLoose[nMuons]/O");
    // electrons
    Events->Branch("nElectrons", &nElectrons);
    Events->Branch("electrons_pt", electrons_pt, "electrons_pt[nElectrons]/F");
    Events->Branch("electrons_eta", electrons_eta, "electrons_eta[nElectrons]/F");
    Events->Branch("electrons_phi", electrons_phi, "electrons_phi[nElectrons]/F");
    Events->Branch("electrons_mass", electrons_mass, "electrons_mass[nElectrons]/F");
    Events->Branch("electrons_miniIso", electrons_miniIso, "electrons_miniIso[nElectrons]/F");
    Events->Branch("electrons_charge", electrons_charge, "electrons_charge[nElectrons]/I");
    Events->Branch("electrons_lepType", electrons_lepType, "electrons_lepType[nElectrons]/I");
    Events->Branch("electrons_passTight", electrons_passTight, "electrons_passTight[nElectrons]/O");
		Events->Branch("electrons_passLoose", electrons_passLoose, "electrons_passLoose[nElectrons]/O");
    // jets
    Events->Branch("nJets", &nJets);
    Events->Branch("jets_pt", jets_pt, "jets_pt[nJets]/F");
    Events->Branch("jets_eta", jets_eta, "jets_eta[nJets]/F");
    Events->Branch("jets_phi", jets_phi, "jets_phi[nJets]/F");
    Events->Branch("jets_mass", jets_mass, "jets_mass[nJets]/F");
    Events->Branch("jets_btagScore", jets_btagScore, "jets_btagScore[nJets]/F");
    Events->Branch("jets_isBtagged", jets_isBtagged, "jets_isBtagged[nJets]/O");
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
		passDblMuTrigs = ev.PassTrigger(DblMuTriggers);
		passEMuTrigs = ev.PassTrigger(EMuTriggers);
		METv_pt = METv.Pt(); METv_phi = METv.Phi();
		genWeight = 1.; trigLumi = -1.;
		if (!IsDATA) {
				genWeight = MCweight();
				trigLumi = ev.GetTriggerLumi("Full");
		}

		// muons
		nMuons = muonV_coll.size();
		for (unsigned int i = 0; i < nMuons; i++) {
				const Muon &mu = muonV_coll.at(i);
				muons_pt[i] = mu.Pt(); 
				muons_eta[i] = mu.Eta(); 
				muons_phi[i] = mu.Phi();
				muons_mass[i] = mu.M();
				muons_miniIso[i] = mu.MiniRelIso();
				muons_charge[i] = mu.Charge();
				muons_lepType[i] = GetLeptonType(mu, truth_coll);
				muons_passTight[i] = mu.PassID(MuonIDs.at(0));
				muons_passLoose[i] = mu.PassID(MuonIDs.at(1));
		}
		// electrons
		nElectrons = eleV_coll.size();
		for (unsigned int i = 0; i < nElectrons; i++) {
        const Electron &ele = eleV_coll.at(i);
        electrons_pt[i] = ele.Pt();
        electrons_eta[i] = ele.Eta();
        electrons_phi[i] = ele.Phi();
        electrons_mass[i] = ele.M();
        electrons_miniIso[i] = ele.MiniRelIso();
        electrons_charge[i] = ele.Charge();
        electrons_lepType[i] = GetLeptonType(ele, truth_coll);
        electrons_passTight[i] = ele.PassID(ElectronIDs.at(0));
				electrons_passLoose[i] = ele.PassID(ElectronIDs.at(1));
    }
		// jets
		nJets = jetT_coll.size();
		const float bcut = mcCorr->GetJetTaggingCutValue(JetTagging::DeepJet, JetTagging::Medium);
		for (unsigned int i = 0; i < nJets; i++) {
        const Jet& j = jetT_coll.at(i);
        jets_pt[i] = j.Pt();
        jets_eta[i] = j.Eta();
        jets_phi[i] = j.Phi();
        jets_mass[i] = j.M();
        jets_btagScore[i] = j.GetTaggerResult(JetTagging::DeepJet);
        jets_isBtagged[i] = j.GetTaggerResult(JetTagging::DeepJet) > bcut;
    }
    Events->Fill();
}

