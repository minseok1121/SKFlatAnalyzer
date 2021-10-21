#include "FakeEstimator.h"

FakeEstimator::FakeEstimator() {}
FakeEstimator::~FakeEstimator() {
	outfile->cd();
	tree->Write();
}

void FakeEstimator::initializeAnalyzer(){
	// flags
	SkimSglMu = HasFlag("SkimSglMu");
	SkimSglEle = HasFlag("SkimSglEle");
	SkimDblMu = HasFlag("SkimDblMu");
	SkimDblEle = HasFlag("SkimDblEle");
	cout << "[FakeEstimator::initializeAnalyzer] SkimSglMu = " << SkimSglMu << endl;
	cout << "[FakeEstimator::initializeAnalyzer] SKimSglEle = " << SkimSglEle << endl;
	cout << "[FakeEstimator::initializeAnalyzer] SkimDblMu = " << SkimDblMu << endl;
	cout << "[FakeEstimator::initializeAnalyzer] SkimDblEle = " << SkimDblEle << endl;

	// Triggers
	if (DataYear == 2017) {
		trigs_sglmu.clear();
		trigs_sglmu.emplace_back("HLT_Mu8_TrkIsoVVL_v");	// 2.605*1.33461
		trigs_sglmu.emplace_back("HLT_Mu17_TrkIsoVVL_v"); // 70.039
		trigs_sglele.clear();
		trigs_sglele.emplace_back("HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_v");
		trigs_sglele.emplace_back("HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v");
		trigs_sglele.emplace_back("HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v");
	}
	else {
		cerr << "No trigger set for DataYear " << DataYear << endl;
		exit(EXIT_FAILURE);
	}

	// ID
	HcToWA_MuID = {"HcToWATight", "HcToWALoose"};
	HcToWA_EleID = {"HcToWATight", "HcToWALoose"};

	// B-tagging (for systematics)
	vector<JetTagging::Parameters> jtps;
	jtps.emplace_back(
			JetTagging::Parameters(JetTagging::DeepJet, JetTagging::Medium, JetTagging::incl, JetTagging::mujets));
	mcCorr->SetJetTaggingParameters(jtps);

	// Set TTree
	tree = new TTree("Events", "Events");
	tree->Branch("IsDATA", &IsDATA); tree->Branch("DataStream", &DataStream);
	tree->Branch("run", &run); tree->Branch("event", &event); tree->Branch("lumi", &lumi);
	tree->Branch("nPV", &nPV); tree->Branch("nPileUp", &nPileUp);
	tree->Branch("passMu8Path", &passMu8Path);
	tree->Branch("passMu17Path", &passMu17Path);
	tree->Branch("passEle8Path", &passEle8Path);
	tree->Branch("passEle12Path", &passEle12Path);
	tree->Branch("passEle23Path", &passEle23Path);
	tree->Branch("nMuons", &nMuons);
    tree->Branch("muons_pt", muons_pt, "muons_pt[nMuons]/F");
    tree->Branch("muons_eta", muons_eta, "muons_eta[nMuons]/F");
    tree->Branch("muons_phi", muons_phi, "muons_phi[nMuons]/F");
    tree->Branch("muons_mass", muons_mass, "muons_mass[nMuons]/F");
    tree->Branch("muons_miniIso", muons_miniIso, "muons_miniIso[nMuons]/F");
    tree->Branch("muons_scaleUp", muons_scaleUp, "muons_scaleUp[nMuons]/F");
    tree->Branch("muons_scaleDown", muons_scaleDown, "muons_scaleDown[nMuons]/F");
    tree->Branch("muons_charge", muons_charge, "muons_charge[nMuons]/I");
    tree->Branch("muons_lepType", muons_lepType, "muons_lepType[nMuons]/I");
    tree->Branch("muons_isTight", muons_isTight, "muons_isTight[nMuons]/O");
    tree->Branch("nElectrons", &nElectrons);
    tree->Branch("electrons_pt", electrons_pt, "electrons_pt[nElectrons]/F");
    tree->Branch("electrons_eta", electrons_eta, "electrons_eta[nElectrons]/F");
    tree->Branch("electrons_phi", electrons_phi, "electrons_phi[nElectrons]/F");
    tree->Branch("electrons_mass", electrons_mass, "electrons_mass[nElectrons]/F");
    tree->Branch("electrons_miniIso", electrons_miniIso, "electrons_miniIso[nElectrons]/F");
    tree->Branch("electrons_scaleUp", electrons_scaleUp, "electrons_scaleUp[nElectrons]/F");
    tree->Branch("electrons_scaleDown", electrons_scaleDown, "electrons_scaleDown[nElectrons]/F");
    tree->Branch("electrons_smearUp", electrons_smearUp, "electrons_smearUp[nElectrons]/F");
    tree->Branch("electrons_smearDown", electrons_smearDown, "electrons_smearDown[nElectrons]/F");
    tree->Branch("electrons_charge", electrons_charge, "electrons_charge[nElectrons]/I");
    tree->Branch("electrons_lepType", electrons_lepType, "electrons_lepType[nElectrons]/I");
    tree->Branch("electrons_isTight", electrons_isTight, "electrons_isTight[nElectrons]/O");
    tree->Branch("nJets", &nJets);
    tree->Branch("jets_pt", jets_pt, "jets_pt[nJets]/F");
    tree->Branch("jets_eta", jets_eta, "jets_eta[nJets]/F");
    tree->Branch("jets_phi", jets_phi, "jets_phi[nJets]/F");
    tree->Branch("jets_mass", jets_mass, "jets_mass[nJets]/F");
    tree->Branch("jets_scaleUp", jets_scaleUp, "jets_scaleUp[nJets]/F");
    tree->Branch("jets_scaleDown", jets_scaleDown, "jets_scaleDown[nJets]/F");
    tree->Branch("jets_smearUp", jets_smearUp, "jets_smearUp[nJets]/F");
    tree->Branch("jets_smearDown", jets_smearDown, "jets_smearDown[nJets]/F");
    tree->Branch("jets_isBtagged", jets_isBtagged, "jets_isBtagged[nJets]/O");
    tree->Branch("METv_pt", &METv_pt);
    tree->Branch("METv_eta", &METv_eta);
    tree->Branch("METv_phi", &METv_phi);
    tree->Branch("genWeight", &genWeight);
    tree->Branch("trigLumi", &trigLumi);
	tree->Branch("weights_L1Prefire", weights_L1Prefire, "weights_L1Prefire[3]/F");
	tree->Branch("weights_btag", weights_btag, "weights_btag[3]/F");
}

void FakeEstimator::executeEvent(){
	
	if (!PassMETFilter()) return;
	Event ev = GetEvent();
	Particle METv = ev.GetMETVector();

	// Object definition
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

	// Event selection
	if (SkimSglMu) {
		if (! ev.PassTrigger(trigs_sglmu)) return;
		if (! (muons_loose.size() == 1 && electrons_loose.size() == 0)) return;
	}
	else if (SkimDblMu) {
		if (! ev.PassTrigger(trigs_sglmu)) return;
		if (! (muons_loose.size() == 2 && electrons_loose.size() == 0)) return;
	}
	else if (SkimSglEle) {
		if (! ev.PassTrigger(trigs_sglele)) return;
		if (! (muons_loose.size() == 0 && electrons_loose.size() == 1)) return;
	}
	else if (SkimDblEle) {
		if (! ev.PassTrigger(trigs_sglele)) return;
		if (! (muons_loose.size() == 0 && electrons_loose.size() == 2)) return;
	}
	else {
		cerr << "[FakeEstimator::executeEvent] Wrong flag" << endl;
		exit(EXIT_FAILURE);
	}

	// Fill contents
	passMu8Path = ev.PassTrigger("HLT_Mu8_TrkIsoVVL_v");
	passMu17Path = ev.PassTrigger("HLT_Mu17_TrkIsoVVL_v");
	passEle8Path = ev.PassTrigger("HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_v");
	passEle12Path = ev.PassTrigger("HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v");
	passEle23Path = ev.PassTrigger("HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v");
    genWeight = 1.;
    trigLumi = 1.;
    for (unsigned int i = 0; i < 3; i++) {
        weights_L1Prefire[i] = 1.;
        weights_btag[i] = 1.;
    }

	if (!IsDATA) {
		genWeight = ev.MCweight()*weight_norm_1invpb;
        trigLumi = ev.GetTriggerLumi("Full");
        weights_L1Prefire[0] = GetPrefireWeight(0); // central
        weights_L1Prefire[1] = GetPrefireWeight(1); // up
        weights_L1Prefire[2] = GetPrefireWeight(-1); // down
        const TString ID = "HcToWATight";

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
        muons_scaleUp[i] = mu.MomentumShift(+1)/mu.Pt();
        muons_scaleDown[i] = mu.MomentumShift(-1)/mu.Pt();
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
        electrons_scaleUp[i] = ele.EnShift(+1);
        electrons_scaleDown[i] = ele.EnShift(-1);
        electrons_smearUp[i] = ele.ResShift(+1);
        electrons_smearDown[i] = ele.ResShift(-1);
        electrons_charge[i] = ele.Charge();
        electrons_lepType[i] = GetLeptonType(ele, gens);
        electrons_isTight[i] = ele.PassID("HcToWATight");
    }
    nJets = jets_lepVeto.size();
    const float bcut = mcCorr->GetJetTaggingCutValue(JetTagging::DeepJet, JetTagging::Medium);

    for (unsigned int i = 0; i< nJets; i++) {
        const Jet& jet = jets_lepVeto.at(i);
        jets_pt[i] = jet.Pt();
        jets_eta[i] = jet.Eta();
        jets_phi[i] = jet.Phi();
        jets_mass[i] = jet.M();
        jets_scaleUp[i] = jet.EnShift(+1);
        jets_scaleDown[i] = jet.EnShift(-1);
        jets_smearUp[i] = jet.ResShift(+1);
        jets_smearDown[i] = jet.ResShift(-1);
        const float this_discr = jet.GetTaggerResult(JetTagging::DeepJet);
        jets_isBtagged[i] = this_discr > bcut;
    }
    METv_pt = METv.Pt();
    METv_eta = METv.Eta();
    METv_phi= METv.Phi();
    tree->Fill();
}
