#include "ClosureTest.h"

ClosureTest::ClosureTest() {}
ClosureTest::~ClosureTest() {
	outfile->cd();
	tree->Write();
}

void ClosureTest::initializeAnalyzer(){
	// flags
	RunTrigger = HasFlag("RunTrigger");
	RunFakeRate = HasFlag("RunFakeRate");
	cout << "[FakeEstimator::initializeAnalyzer] RunTrigger = " << RunTrigger << endl;
	cout << "[FakeEstimator::initializeAnalyzer] RunFakeRate = " << RunFakeRate << endl;

	// Set Triggers
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

	// Set Tree
	tree = new TTree("Events", "Events");
	tree->Branch("run", &run); tree->Branch("event", &event); tree->Branch("lumi", &lumi);
	tree->Branch("genWeight", &genWeight); tree->Branch("L1PrefireWeight", &L1PrefireWeight);
	tree->Branch("nPV"); tree->Branch("nPileUp");
	tree->Branch("passDblMuTrigs", &passDblMuTrigs); tree->Branch("passEMuTrigs", &passEMuTrigs);
	tree->Branch("trigDblMuEff", &trigDblMuEff);
	tree->Branch("trigDblMuEffUp", &trigDblMuEffUp);
	tree->Branch("trigDblMuEffDown", &trigDblMuEffDown);
	tree->Branch("trigEMuEff", &trigEMuEff);
	tree->Branch("trigEMuEffUp", &trigEMuEffUp);
	tree->Branch("trigEMuEffDown", &trigEMuEffDown);
	tree->Branch("nMuons", &nMuons);
    tree->Branch("muons_pt", muons_pt, "muons_pt[nMuons]/D");
    tree->Branch("muons_eta", muons_eta, "muons_eta[nMuons]/D");
    tree->Branch("muons_phi", muons_phi, "muons_phi[nMuons]/D");
    tree->Branch("muons_mass", muons_mass, "muons_mass[nMuons]/D");
	tree->Branch("nElectrons", &nElectrons);
    tree->Branch("electrons_pt", electrons_pt, "electrons_pt[nElectrons]/D");
    tree->Branch("electrons_eta", electrons_eta, "electrons_eta[nElectrons]/D");
    tree->Branch("electrons_phi", electrons_phi, "electrons_phi[nElectrons]/D");
    tree->Branch("electrons_mass", electrons_mass, "electrons_mass[nElectrons]/D");	
}

void ClosureTest::executeEvent(){

	if (!PassMETFilter()) return;
	Event ev = GetEvent();
	Particle METv = ev.GetMETVector();
	vector<Muon> muons = GetAllMuons();
	vector<Electron> electrons = GetAllElectrons();

	sort(muons.begin(), muons.end(), PtComparing);
	sort(electrons.begin(), electrons.end(), PtComparing);
	
	vector<Muon> muons_tight = SelectMuons(muons, "HcToWATight", 10., 2.4);
	vector<Muon> muons_loose = SelectMuons(muons, "HcToWALoose", 10., 2.4);
	vector<Electron> electrons_tight = SelectElectrons(electrons, "HcToWATight", 10., 2.5);
	vector<Electron> electrons_loose = SelectElectrons(electrons, "HcToWALoose", 10., 2.5);

	// just skim on nLetons base
	TString channel = "";
	if (muons_tight.size() == 3 && muons_loose.size() == 3 && electrons_tight.size() == 1 && electrons_loose.size() == 0)
		channel = "SR_3mu";
	if (muons_tight.size() == 2 && muons_loose.size() == 2 && electrons_tight.size() == 1 && electrons_loose.size() == 1)
		channel = "SR_1e2mu";
	if (channel == "") return;

	passDblMuTrigs = ev.PassTrigger(trigs_dblmu);
	passEMuTrigs = ev.PassTrigger(trigs_emu);

	// Trigger Efficiencies
	trigDblMuEff = 0.; trigDblMuEffUp = 0.; trigDblMuEffDown = 0.;
	trigEMuEff = 0.; trigEMuEffUp = 0.; trigEMuEffDown = 0.;
	trigDblMuEff 
		= mcCorr->TriggerEfficiency(electrons_tight, muons_tight, "DiMuIso_HNTopID", false, "");
	trigDblMuEffUp 
		= mcCorr->TriggerEfficiency(electrons_tight, muons_tight, "DiMuIso_HNTopID", false, "Syst_Up");
	trigDblMuEffDown 
		= mcCorr->TriggerEfficiency(electrons_tight, muons_tight, "DiMuIso_HNTopID", false, "Syst_Down");
	if (channel == "SR_1e2mu") {
		trigEMuEff
			= mcCorr->TriggerEfficiency(electrons_tight, muons_tight, "EMuIso_HNTopID", false, "");
		trigEMuEffUp
			= mcCorr->TriggerEfficiency(electrons_tight, muons_tight, "EMuIso_HNTopID", false, "Syst_Up");
		trigEMuEffDown
			= mcCorr->TriggerEfficiency(electrons_tight, muons_tight, "EMuIso_HNTopID", false, "Syst_Down");
	}
	genWeight = ev.MCweight()*weight_norm_1invpb;
	L1PrefireWeight = GetPrefireWeight(0);
	nMuons = muons_tight.size();
	for (unsigned int i = 0; i < nMuons; i++) {
		const Muon& mu = muons_tight.at(i);
		muons_pt[i] = mu.Pt();
		muons_eta[i] = mu.Eta();
		muons_phi[i] = mu.Phi();
		muons_mass[i] = mu.M();
	}
	nElectrons = electrons_tight.size();
	for (unsigned int i = 0; i < nElectrons; i++) {
		const Electron& ele = electrons_tight.at(i);
		electrons_pt[i] = ele.Pt();
		electrons_eta[i] = ele.Eta();
		electrons_phi[i] = ele.Phi();
		electrons_mass[i] = ele.M();
	}
	tree->Fill();
}
