#include "Selector.h"

Selector::Selector() {}
Selector::~Selector() {
	outfile->cd();
	tree->Write();
}

void Selector::initializeAnalyzer(){
	// flags
	SkimDiMu = HasFlag("SkimDiMu");
	SkimEMu = HasFlag("SkimEMu");
	Skim3Mu = HasFlag("Skim3Mu");
	Skim1E2Mu = HasFlag("Skim1E2Mu");
	cout << "[Selector::initializeAnalyzer] SkimDiMu = " << SkimDiMu << endl;
	cout << "[Selector::initializeAnalyzer] SkimEMu = " << SkimEMu << endl;
	cout << "[Selector::initializeAnalyzer] Skim3Mu = " << Skim3Mu << endl;
	cout << "[Selector::initializeAnalyzer] Skim1E2Mu = " << Skim1E2Mu << endl;

	// Triggers
	if (DataYear == 2017) {
        trigs_dblmu.clear(); trigs_emu.clear();
        trigs_dblmu.emplace_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v");
        trigs_dblmu.emplace_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v");
        trigs_dblmu.emplace_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v");
        trigs_emu.emplace_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v");
        trigs_emu.emplace_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");
    }
    else {
        cerr << "DataYear " << DataYear << " is not set" << endl;
        exit(EXIT_FAILURE);
    }

	// ID
    HcToWA_MuID = {"HcToWATight", "HcToWALoose"};
    HcToWA_EleID = {"HcToWATight", "HcToWALoose"};

	// Jet tagger
	vector<JetTagging::Parameters> jtps;
	jtps.emplace_back(
			JetTagging::Parameters(JetTagging::DeepJet, JetTagging::Medium, JetTagging::incl, JetTagging::mujets));
	//jtps.emplace_back(
	//		JetTagging::Parameters(JetTagging::DeepCSV, JetTagging::Medium, JetTagging::incl, JetTagging::mujets));
	mcCorr->SetJetTaggingParameters(jtps);

	// Set TTree
	tree = new TTree("Events", "Events");
	tree->Branch("IsDATA", &IsDATA); tree->Branch("DataStream", &DataStream);
	tree->Branch("run", &run); tree->Branch("event", &event); tree->Branch("lumi", &lumi);
	tree->Branch("nPV", &nPV); tree->Branch("nPileUp", &nPileUp);
	tree->Branch("passDblMuTrigs", &passDblMuTrigs);
	tree->Branch("passEMuTrigs", &passEMuTrigs);
    tree->Branch("nMuons", &nMuons);
    tree->Branch("muons_pt", muons_pt, "muons_pt[nMuons]/F");
    tree->Branch("muons_eta", muons_eta, "muons_eta[nMuons]/F");
    tree->Branch("muons_phi", muons_phi, "muons_phi[nMuons]/F");
    tree->Branch("muons_mass", muons_mass, "muons_mass[nMuons]/F");
	tree->Branch("muons_miniIso", muons_miniIso, "muons_miniIso[nMuons]/F");
    tree->Branch("muons_charge", muons_charge, "muons_charge[nMuons]/I");
    tree->Branch("muons_lepType", muons_lepType, "muons_lepType[nMuons]/I");
	tree->Branch("muons_isTight", muons_isTight, "muons_isTight[nMuons]/O");
    tree->Branch("nElectrons", &nElectrons);
    tree->Branch("electrons_pt", electrons_pt, "electrons_pt[nElectrons]/F");
    tree->Branch("electrons_eta", electrons_eta, "electrons_eta[nElectrons]/F");
    tree->Branch("electrons_phi", electrons_phi, "electrons_phi[nElectrons]/F");
    tree->Branch("electrons_mass", electrons_mass, "electrons_mass[nElectrons]/F");
	tree->Branch("electrons_miniIso", electrons_miniIso, "electrons_miniIso[nElectrons]/F");
    tree->Branch("electrons_charge", electrons_charge, "electrons_charge[nElectrons]/I");
    tree->Branch("electrons_lepType", electrons_lepType, "electrons_lepType[nElectrons]/I");
	tree->Branch("electrons_isTight", electrons_isTight, "electrons_isTight[nElectrons]/O");
    tree->Branch("nJets", &nJets);
    tree->Branch("jets_pt", jets_pt, "jets_pt[nJets]/F");
    tree->Branch("jets_eta", jets_eta, "jets_eta[nJets]/F");
    tree->Branch("jets_phi", jets_phi, "jets_phi[nJets]/F");
	tree->Branch("jets_mass", jets_mass, "jets_mass[nJets]/F");
	tree->Branch("jets_isBtagged", jets_isBtagged, "jets_isBtagged[nJets]/O");
	tree->Branch("METv_Pt", &METv_pt);
	tree->Branch("METv_eta", &METv_eta);
	tree->Branch("METv_phi", &METv_phi);
	tree->Branch("genWeight", &genWeight);
	tree->Branch("trigLumi", &trigLumi);
	tree->Branch("weights_L1Prefire", weights_L1Prefire, "weights_L1Prefire[3]/F");
	tree->Branch("weights_id", weights_id, "weights_id[3]/F");
	tree->Branch("weights_trigs_dblmu", weights_trigs_dblmu, "weights_trigs_dblmu[3]/F");
	tree->Branch("weights_trigs_emu", weights_trigs_emu, "weights_trigs_emu[3]/F");
	tree->Branch("weights_btag", weights_btag, "weights_btag[3]/F");
}

void Selector::executeEvent(){

	if (!PassMETFilter()) return;
	Event ev = GetEvent();
	Particle METv = ev.GetMETVector();

	// Object definitions
	vector<Gen> gens = GetGens();
	vector<Muon> muons = GetAllMuons();
	vector<Electron> electrons = GetAllElectrons();
	vector<Jet> jets = GetAllJets();
	
	// sort at the first time, don't want to be confused
	sort(muons.begin(), muons.end(), PtComparing);
	sort(electrons.begin(), electrons.end(), PtComparing);
	sort(jets.begin(), jets.end(), PtComparing);

	vector<Muon> muons_loose = SelectMuons(muons, "HcToWALoose", 10., 2.4);
	vector<Electron> electrons_loose = SelectElectrons(electrons, "HcToWALoose", 10., 2.5);
	vector<Jet> jets_tight = SelectJets(jets, "tight", 20., 2.4);
	vector<Jet> jets_lepVeto = JetsVetoLeptonInside(jets_tight, electrons_loose, muons_loose, 0.4);

	// Events should pass triggers & trigger safe pt cut
	if (SkimDiMu) {
		if (! (muons_loose.size() == 2 && electrons_loose.size() == 0)) return;
	}
	else if (SkimEMu) {
		if (! (muons_loose.size() == 1 && electrons_loose.size() == 1)) return;
	}
	else if (Skim3Mu) {
		if (! (muons_loose.size() == 3 && electrons_loose.size() == 0)) return;
	}
	else if (Skim1E2Mu) {
		if (! (muons_loose.size() == 2 && electrons_loose.size() == 1)) return;
	}
	else {
		cerr << "[Selector::executeEvent] Wrong flag" << endl;
		exit(EXIT_FAILURE);
	}

	// Initialize contents
	passDblMuTrigs = ev.PassTrigger(trigs_dblmu);
	passEMuTrigs = ev.PassTrigger(trigs_emu);
	genWeight = 1.;
	trigLumi = 1.;
	for (unsigned int i = 0; i < 3; i++) {
		weights_L1Prefire[i] = 1.;
		weights_id[i] = 1.;
		weights_trigs_dblmu[i] = 1.;
		weights_trigs_emu[i] = 1.;
		weights_btag[i] = 1.;
	}

	if (!IsDATA) {
		genWeight = ev.MCweight()*weight_norm_1invpb;
		trigLumi = ev.GetTriggerLumi("Full");
		weights_L1Prefire[0] = GetPrefireWeight(0); // central
		weights_L1Prefire[1] = GetPrefireWeight(1); // up
		weights_L1Prefire[2] = GetPrefireWeight(-1); // down
		const TString ID = "HcToWATight";

		if (SkimDiMu) {
			const Muon& mu1 = muons_loose.at(0);
			const Muon& mu2 = muons_loose.at(1);
			const float mu1_idsf = mcCorr->MuonID_SF(ID, mu1.Eta(), mu1.MiniAODPt(), 0);
			const float mu1_idsf_up = mcCorr->MuonID_SF(ID, mu1.Eta(), mu1.MiniAODPt(), 1);
			const float mu1_idsf_down = mcCorr->MuonID_SF(ID, mu1.Eta(), mu1.MiniAODPt(), -1);
			const float mu2_idsf = mcCorr->MuonID_SF(ID, mu2.Eta(), mu2.MiniAODPt(), 0);
			const float mu2_idsf_up = mcCorr->MuonID_SF(ID, mu2.Eta(), mu2.MiniAODPt(), 1);
			const float mu2_idsf_down = mcCorr->MuonID_SF(ID, mu2.Eta(), mu2.MiniAODPt(), -1);
			weights_id[0] = mu1_idsf*mu2_idsf;
			weights_id[1] = mu1_idsf_up*mu2_idsf_up;
			weights_id[2] = mu1_idsf_down*mu2_idsf_down;
			weights_trigs_dblmu[0] = mcCorr->GetTriggerSF(electrons_loose, muons_loose, "DiMuIso_HNTopID", "");
			weights_trigs_dblmu[1] = mcCorr->GetTriggerSF(electrons_loose, muons_loose, "DiMuIso_HNTopID", "Syst_Up");
			weights_trigs_dblmu[2] = mcCorr->GetTriggerSF(electrons_loose, muons_loose, "DiMuIso_HNTopID", "Syst_Down");
			weights_trigs_emu[0] = 0.;
			weights_trigs_emu[1] = 0.;
			weights_trigs_emu[2] = 0.;
		}
		else if (SkimEMu) {
			const Muon& mu = muons_loose.at(0);
        	const Electron& ele = electrons_loose.at(0);
        	const float mu_idsf = mcCorr->MuonID_SF(ID, mu.Eta(), mu.MiniAODPt(), 0);
			const float mu_idsf_up = mcCorr->MuonID_SF(ID, mu.Eta(), mu.MiniAODPt(), +1);
			const float mu_idsf_down = mcCorr->MuonID_SF(ID, mu.Eta(), mu.MiniAODPt(), -1);
        	const float ele_idsf = mcCorr->ElectronID_SF(ID, ele.scEta(), ele.Pt(), 0);
			const float ele_idsf_up = mcCorr->ElectronID_SF(ID, ele.scEta(), ele.Pt(), +1);
			const float ele_idsf_down = mcCorr->ElectronID_SF(ID, ele.scEta(), ele.Pt(), -1);
			weights_id[0] = mu_idsf*ele_idsf;
			weights_id[1] = mu_idsf_up*ele_idsf_up;
			weights_id[2] = mu_idsf_down*ele_idsf_down;
			weights_trigs_emu[0] = mcCorr->GetTriggerSF(electrons_loose, muons_loose, "EMuIso_HNTopID", "");
			weights_trigs_emu[1] = mcCorr->GetTriggerSF(electrons_loose, muons_loose, "EMuIso_HNTopID", "Syst_Up");
			weights_trigs_emu[2] = mcCorr->GetTriggerSF(electrons_loose, muons_loose, "EMuIsoHNTopID", "Syst_Down");
			weights_trigs_dblmu[0] = 0.;
			weights_trigs_dblmu[1] = 0.;
			weights_trigs_dblmu[2] = 0.;
		}
		else if (Skim3Mu) {
			const Muon& mu1 = muons_loose.at(0);
			const Muon& mu2 = muons_loose.at(1);
			const Muon& mu3 = muons_loose.at(2);
			const float mu1_idsf = mcCorr->MuonID_SF(ID, mu1.Eta(), mu1.MiniAODPt(), 0);
			const float mu1_idsf_up = mcCorr->MuonID_SF(ID, mu1.Eta(), mu1.MiniAODPt(), +1);
			const float mu1_idsf_down = mcCorr->MuonID_SF(ID, mu1.Eta(), mu1.MiniAODPt(), -1);
			const float mu2_idsf = mcCorr->MuonID_SF(ID, mu2.Eta(), mu2.MiniAODPt(), 0);
			const float mu2_idsf_up = mcCorr->MuonID_SF(ID, mu2.Eta(), mu2.MiniAODPt(), +1);
			const float mu2_idsf_down = mcCorr->MuonID_SF(ID, mu2.Eta(), mu2.MiniAODPt(), -1);
			const float mu3_idsf = mcCorr->MuonID_SF(ID, mu3.Eta(), mu3.MiniAODPt(), 0);
			const float mu3_idsf_up = mcCorr->MuonID_SF(ID, mu3.Eta(), mu3.MiniAODPt(), +1);
			const float mu3_idsf_down = mcCorr->MuonID_SF(ID, mu3.Eta(), mu3.MiniAODPt(), -1);
			weights_id[0] = mu1_idsf*mu2_idsf*mu3_idsf;
			weights_id[1] = mu1_idsf_up*mu2_idsf_up*mu3_idsf_up;
			weights_id[2] = mu1_idsf_down*mu2_idsf_down*mu3_idsf_down;
			weights_trigs_dblmu[0] = mcCorr->GetTriggerSF(electrons_loose, muons_loose, "DiMuIso_HNTopID", "");
			weights_trigs_dblmu[1] = mcCorr->GetTriggerSF(electrons_loose, muons_loose, "DiMuIso_HNTopID", "Syst_Up");
			weights_trigs_dblmu[2] = mcCorr->GetTriggerSF(electrons_loose, muons_loose, "DiMuIso_HNTopID", "Syst_Down");
			weights_trigs_emu[0] = 0.;
			weights_trigs_emu[1] = 0.;
			weights_trigs_emu[2] = 0.;
		}
		else if (Skim1E2Mu) {
			const Muon& mu1 = muons_loose.at(0);
        	const Muon& mu2 = muons_loose.at(1);
        	const Electron& ele = electrons_loose.at(0);
			const float mu1_idsf = mcCorr->MuonID_SF(ID, mu1.Eta(), mu1.MiniAODPt(), 0);
			const float mu1_idsf_up = mcCorr->MuonID_SF(ID, mu1.Eta(), mu1.MiniAODPt(), +1);
			const float mu1_idsf_down = mcCorr->MuonID_SF(ID, mu1.Eta(), mu1.MiniAODPt(), -1);
			const float mu2_idsf = mcCorr->MuonID_SF(ID, mu2.Eta(), mu2.MiniAODPt(), 0);
			const float mu2_idsf_up = mcCorr->MuonID_SF(ID, mu2.Eta(), mu2.MiniAODPt(), +1);
			const float mu2_idsf_down = mcCorr->MuonID_SF(ID, mu2.Eta(), mu2.MiniAODPt(), -1);
			const float ele_idsf = mcCorr->ElectronID_SF(ID, ele.scEta(), ele.Pt(), 0);
			const float ele_idsf_up = mcCorr->ElectronID_SF(ID, ele.scEta(), ele.Pt(), +1);
			const float ele_idsf_down = mcCorr->ElectronID_SF(ID, ele.scEta(), ele.Pt(), -1);
			weights_id[0] = mu1_idsf*mu2_idsf*ele_idsf;
			weights_id[1] = mu1_idsf_up*mu2_idsf_up*ele_idsf_up;
			weights_id[2] = mu1_idsf_down*mu2_idsf_down*ele_idsf_down;
			// trigger modeling as w_trig = 1. - (1. - w_dblmu)*(1. - w_emu)
			weights_trigs_dblmu[0] = mcCorr->GetTriggerSF(electrons_loose, muons_loose, "DiMuIso_HNTopID", "");
			weights_trigs_dblmu[1] = mcCorr->GetTriggerSF(electrons_loose, muons_loose, "DiMuIso_HNTopID", "Syst_Up");
			weights_trigs_dblmu[2] = mcCorr->GetTriggerSF(electrons_loose, muons_loose, "DiMuIso_HNTopID", "Syst_Down");
			weights_trigs_emu[0] = mcCorr->GetTriggerSF(electrons_loose, muons_loose, "EMuIso_HNTopID", "");
			weights_trigs_emu[1] = mcCorr->GetTriggerSF(electrons_loose, muons_loose, "EMuIso_HNTopID", "Syst_Up");
			weights_trigs_emu[2] = mcCorr->GetTriggerSF(electrons_loose, muons_loose, "EMuIso_HNTopID", "Syst_Down");
		}
		else {
			exit(EXIT_FAILURE);
		}

		JetTagging::Parameters jtp_DeepJet_Medium
				 = JetTagging::Parameters(JetTagging::DeepJet, JetTagging::Medium, JetTagging::incl, JetTagging::mujets);
		// no b-tagging systematics yet
		weights_btag[0] = mcCorr->GetBTaggingReweight_1a(jets_lepVeto, jtp_DeepJet_Medium);
		//weights_btag[0] = 1.;
		weights_btag[1] = 0.;
		weights_btag[2] = 0.;
	}

	nMuons = muons_loose.size();
	for (unsigned int i = 0; i < nMuons; i++) {
		const Muon& mu = muons_loose.at(i);
		muons_pt[i] = mu.Pt();
		muons_eta[i] = mu.Eta();
		muons_phi[i] = mu.Phi();
		muons_mass[i] = mu.M();
		muons_miniIso[i] = mu.MiniRelIso();
		muons_charge[i] = mu.Charge();
		muons_lepType[i] = GetLeptonType(mu, gens);
		muons_isTight[i] = mu.PassID("HcToWATight");
	}
	nElectrons = electrons_loose.size();
	for (unsigned int i = 0; i < nElectrons; i++) {
		const Electron& ele = electrons_loose.at(i);
		electrons_pt[i] = ele.Pt();
		electrons_eta[i] = ele.Eta();
		electrons_phi[i] = ele.Phi();
		electrons_mass[i] = ele.M();
		electrons_miniIso[i] = ele.MiniRelIso();
		electrons_charge[i] = ele.Charge();
		electrons_lepType[i] = GetLeptonType(ele, gens);
		electrons_isTight[i] = ele.PassID("HcToWATight");
	}
	nJets = jets_lepVeto.size();
	const float bcut = mcCorr->GetJetTaggingCutValue(JetTagging::DeepJet, JetTagging::Medium);
	//const float bcut = mcCorr->GetJetTaggingCutValue(JetTagging::DeepCSV, JetTagging::Medium);
	for (unsigned int i = 0; i< nJets; i++) {
		const Jet& jet = jets_lepVeto.at(i);
		jets_pt[i] = jet.Pt();
		jets_eta[i] = jet.Eta();
		jets_phi[i] = jet.Phi();
		jets_mass[i] = jet.M();
		const float this_discr = jet.GetTaggerResult(JetTagging::DeepJet);
		//const float this_discr = jet.GetTaggerResult(JetTagging::DeepCSV);
		jets_isBtagged[i] = this_discr > bcut;
	}
	METv_pt = METv.Pt();
	METv_eta = METv.Eta();
	METv_phi= METv.Phi();
	tree->Fill();
}

