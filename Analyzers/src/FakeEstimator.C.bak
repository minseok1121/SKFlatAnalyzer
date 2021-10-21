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
	cout << "[FakeEstimator::initializeAnalyzer] RunSysts = " << RunSysts << endl;
	cout << "[FakeEstimator::initializeAnalyzer] SkimSglMu = " << SkimSglMu << endl;
	cout << "[FakeEstimator::initializeAnalyzer] SKimSglEle = " << SkimSglEle << endl;
	cout << "[FakeEstimator::initializeAnalyzer] SkimDblMu = " << SkimDblMu << endl;
	cout << "[FakeEstimator::initializeAnalyzer] SkimDblEle = " << SkimDblEle << endl;
	
	// Set TTree
	tree = new TTree("Events", "Events");
	tree->Branch("run", &run); tree->Branch("event", &event); tree->Branch("lumi", &lumi);
	tree->Branch("nPV", &nPV); tree->Branch("nPileUp", &nPileUp);
	tree->Branch("genWeight", &genWeight);
	tree->Branch("L1PrefireWeight", &L1PrefireWeight);
	tree->Branch("passMu8Path", &passMu8Path);
	tree->Branch("passMu17Path", &passMu17Path);
	tree->Branch("passEle8Path", &passEle8Path);
	tree->Branch("passEle12Path", &passEle12Path);
	tree->Branch("passEle23Path", &passEle23Path);
	tree->Branch("nMuons", &nMuons);
	tree->Branch("muons_pt", muons_pt, "muons_pt[nMuons]/D");
    tree->Branch("muons_eta", muons_eta, "muons_eta[nMuons]/D");
    tree->Branch("muons_phi", muons_phi, "muons_phi[nMuons]/D");
    tree->Branch("muons_mass", muons_mass, "muons_mass[nMuons]/D");
	tree->Branch("muons_miniIso", muons_miniIso, "muons_miniIso[nMuons]/D");
	tree->Branch("muons_relIso", muons_relIso, "muons_relIso[nMuons]/D");
    tree->Branch("muons_charge", muons_charge, "muons_charge[nMuons]/I");
    tree->Branch("muons_lepType", muons_lepType, "muons_lepType[nMuons]/I");
	tree->Branch("muons_isTight", muons_isTight, "muons_isTight[nMuons]/I");
    tree->Branch("nElectrons", &nElectrons);
    tree->Branch("electrons_pt", electrons_pt, "electrons_pt[nElectrons]/D");
    tree->Branch("electrons_eta", electrons_eta, "electrons_eta[nElectrons]/D");
    tree->Branch("electrons_phi", electrons_phi, "electrons_phi[nElectrons]/D");
    tree->Branch("electrons_mass", electrons_mass, "electrons_mass[nElectrons]/D");
	tree->Branch("electrons_miniIso", electrons_miniIso, "electrons_miniIso[nElectrons]/D");
	tree->Branch("electrons_relIso", electrons_relIso, "electrons_relIso[nElectrons]/D");
    tree->Branch("electrons_charge", electrons_charge, "electrons_charge[nElectrons]/I");
    tree->Branch("electrons_lepType", electrons_lepType, "electrons_lepType[nElectrons]/I");
    tree->Branch("electrons_isTight", electrons_isTight, "electrons_isTight[nElectrons]/I");
	tree->Branch("nJets", &nJets);
    tree->Branch("jets_pt", jets_pt, "jets_pt[nJets]/D");
    tree->Branch("jets_eta", jets_eta, "jets_eta[nJets]/D");
    tree->Branch("jets_phi", jets_phi, "jets_phi[nJets]/D");
    tree->Branch("jets_isBtagged", jets_isBtagged, "jets_isBtagged[nJets]/I");
    tree->Branch("METv_pt", &METv_pt);
    tree->Branch("METv_phi", &METv_phi);
		
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
	jtps.emplace_back(JetTagging::Parameters(JetTagging::DeepJet, JetTagging::Medium, JetTagging::incl, JetTagging::mujets));
	mcCorr->SetJetTaggingParameters(jtps);

	// Systematics
	Systematics = {"Central"};
	if (RunSysts) 
		Systematics = {
				"Central", "JetResUp", "JetResDown", "JetEnUp", "JetEnDown",
				"MuonEnUp", "MuonEnDown", "ElectronResUp", "ElectronResDown",
				"ElectronEnUp", "ElectronEnDown"};
}

void FakeEstimator::executeEvent(){
	// Get All Objects
	gens = GetGens();
	muons_all = GetAllMuons();
	electrons_all = GetAllElectrons();
	jets_all = GetAllJets();
	
	for (const auto& syst: Systematics)
		executeEventWithSystematics(syst);
}

void FakeEstimator::executeEventWithSystematics(const TString& syst) {
	if (!PassMETFilter()) return;

	Event ev = GetEvent();
	Particle METv = ev.GetMETVector();
	vector<Muon> muons = muons_all;
	vector<Electron> electrons = electrons_all;
	vector<Jet> jets = jets_all;

	// Jet Resolution, Energy Scale
	if (syst == "JetResUp")
		jets = SmearJets(jets, +1);
	if (syst == "JetResDown")
		jets = SmearJets(jets, -1);
	if (syst == "JetEnUp")
		jets = ScaleJets(jets, +1);
	if (syst == "JetEnDown")
		jets = ScaleJets(jets, -1);
	if (syst == "MuonEnUp")
		muons = ScaleMuons(muons, +1);
	if (syst == "MuonEnDown")
		muons = ScaleMuons(muons, -1);
	if (syst == "ElectronResUp")
		electrons = SmearElectrons(electrons, +1);
	if (syst == "ElectronResDown")
		electrons = SmearElectrons(electrons, -1);
	if (syst == "ElectronEnUp")
		electrons = ScaleElectrons(electrons, +1);
	if (syst == "ElectronEnDown")
		electrons = ScaleElectrons(electrons, -1);

	sort(muons.begin(), muons.end(), PtComparing);
	sort(electrons.begin(), electrons.end(), PtComparing);
	sort(jets.begin(), jets.end(), PtComparing);

	vector<Muon> muons_loose = SelectMuons(muons, "HcToWALoose", 10., 2.4);
	vector<Electron> electrons_loose = SelectElectrons(electrons, "HcToWALoose", 10., 2.5);
	vector<Jet> jets_tight = SelectJets(jets, "tight", 20., 2.4);
	vector<Jet> jets_lepVeto = JetsVetoLeptonInside(jets_tight, electrons_loose, muons_loose, 0.4);

	// Event Selection
	// Fire the trigger for the given flavor of the loose lepton
	// Exactly one loose lepton present in the event with p_T > 10 GeV
	// At least one jet of p_T > 40 GeV with |\eta| < 2.4
	// deltaR(j, l) > 1.0
	TString channel = "";
	if (! (jets_lepVeto.size() > 0)) return;
	if (! (jets_lepVeto.at(0).Pt() > 40)) return;

	if (muons_loose.size() == 1 && electrons_loose.size() == 0) channel = "QCD_SglMu";
	if (muons_loose.size() == 0 && electrons_loose.size() == 1) channel = "QCD_SglEle";
	if (muons_loose.size() == 2 && electrons_loose.size() == 0) channel = "QCD_DblMu";
	if (muons_loose.size() == 0 && electrons_loose.size() == 2) channel = "QCD_DblEle";
	if (channel == "") return;

	if (channel == "QCD_SglMu")
		if (! (muons_loose.at(0).DeltaR(jets_lepVeto.at(0)) > 1.0)) return;
	if (channel == "QCD_SglEle")
		if (! (electrons_loose.at(0).DeltaR(jets_lepVeto.at(0)) > 1.0)) return;

	if (SkimSglMu) {
		if (channel != "QCD_SglMu") return;
		passMu8Path = ev.PassTrigger("HLT_Mu8_TrkIsoVVL_v");
		passMu17Path = ev.PassTrigger("HLT_Mu17_TrkIsoVVL_v");

		if (! (passMu8Path || passMu17Path)) return;
	}
	else if (SkimSglEle) {
		if (channel != "QCD_SglEle") return;
		passEle8Path = ev.PassTrigger("HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_v");
		passEle12Path = ev.PassTrigger("HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v");
		passEle23Path = ev.PassTrigger("HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v");

		if (! (passEle8Path || passEle12Path || passEle23Path)) return;
	}
	else if (SkimDblMu) {
		if (channel != "QCD_DblMu") return;
		passMu8Path = ev.PassTrigger("HLT_Mu8_TrkIsoVVL_v");
		passMu17Path = ev.PassTrigger("HLT_Mu17_TrkIsoVVL_v");
		const Muon& mu1 = muons_loose.at(0);
		const Muon& mu2 = muons_loose.at(1);
		const Particle ZCand = mu1 + mu2;
		if (! (passMu8Path || passMu17Path)) return;
		if (fabs(ZCand.M() - 91.2) > 15.) return;
		if (mu1.Charge()+mu2.Charge() != 0) return;
	}
	else if (SkimDblEle) {
		if (channel != "QCD_DblEle") return;
		passEle8Path = ev.PassTrigger("HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_v");
		passEle12Path = ev.PassTrigger("HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v");
		passEle23Path = ev.PassTrigger("HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v");
		const Electron& ele1 = electrons_loose.at(0);
		const Electron& ele2 = electrons_loose.at(1);
		const Particle ZCand = ele1 + ele2;

		if (! (passEle8Path || passEle12Path || passEle23Path)) return;
		if (fabs(ZCand.M() - 91.2) > 15.) return;
		if (ele1.Charge()+ele2.Charge() != 0) return;
	}
	else {
		cout << "[FakeEstimator] Wrong Flag" << endl;
		exit(EXIT_FAILURE);
	}

	// Tree contents
	genWeight = ev.MCweight()*weight_norm_1invpb;
	L1PrefireWeight = GetPrefireWeight(0);
	nMuons = muons_loose.size();
	for (unsigned int i = 0; i < nMuons; i++) {
		const Muon& mu = muons_loose.at(i);
		muons_pt[i] = mu.Pt();
		muons_eta[i] = mu.Eta();
		muons_phi[i] = mu.Phi();
		muons_mass[i] = mu.M();
		muons_miniIso[i] = mu.MiniRelIso();
		muons_relIso[i] = mu.RelIso();
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
        electrons_relIso[i] = ele.RelIso();
        electrons_charge[i] = ele.Charge();
        electrons_lepType[i] = GetLeptonType(ele, gens);
        electrons_isTight[i] = ele.PassID("HcToWATight");
	}
	nJets = jets_lepVeto.size();
	for (unsigned int i = 0; i < nJets; i++) {
		const Jet& j = jets_lepVeto.at(i);
		jets_pt[i] = j.Pt();
		jets_eta[i] = j.Eta();
		jets_phi[i] = j.Phi();
		const double this_discr = 
			jets_lepVeto.at(i).GetTaggerResult(JetTagging::DeepJet);
		jets_isBtagged[i] =
			(this_discr > mcCorr->GetJetTaggingCutValue(JetTagging::DeepJet, JetTagging::Medium));
	}
	METv_pt = METv.Pt();
	METv_phi = METv.Phi();
	tree->Fill();
}


