#include "opt_pf_electron.h"

void opt_pf_electron::initializeAnalyzer(){
	
	//==== Flags for systematic sources ====
	RunSysts = HasFlag("RunSysts");
	cout << "[opt_pf_electron::initializeAnalyzer] RunSysts = " << RunSysts << endl;

	//==== Electorn IDs setting ====
	ElectronIDs = {"passMVAID_iso_WP90", "passTightID", "HcToWAT"};
	
	//==== Trigger Settings ====
	if (DataYear == 2016) {
		TrigNames = {"HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v"};
		TriggerSafePtCut = 25.;
	}
	else if (DataYear == 2017 || DataYear == 2018) {
		cout << "[opt_pf_electron::initializeAnalyzer] Trigger is not set yet" << endl;
		exit(EXIT_FAILURE);
	}
	else {
		cout << "[opt_pf_electron::initializeAnalyzer] Wrong Year" << endl;
		exit(EXIT_FAILURE);
	}

	//==== B tagging ====
	vector<JetTagging::Parameters> jtps;
	jtps.push_back( JetTagging::Parameters(JetTagging::DeepCSV, JetTagging::Medium, JetTagging::incl, JetTagging::mujets) );
	mcCorr->SetJetTaggingParameters(jtps);

	cout << "[opt_pf_electron::initializeAnalyzer] Finish Initialization" << endl;
}

void opt_pf_electron::executeEvent(){
	
	//==== Copy all objects ====
	AllMuons = GetAllMuons();
	AllElectrons = GetAllElectrons();
	AllJets = GetAllJets();

	//==== Get L1Prefire, pileup reweight ====
	w_prefire = GetPrefireWeight(0);
	w_pileup = GetPileUpWeight(nPileUp, 0);

	AnalyzerParameter param;
	
	//==== Loop over electronIDs ====
	for (const auto &ID : ElectronIDs) {
		MuonVetoID = "POGLoose";
		ElecVetoID = "passVetoID";

		ElectronID = ID;
		JetID = "tight";

		//==== No systematic setting ====
		param.Clear();
		param.syst_ = AnalyzerParameter::Central;
		param.Name = ElectronID + "_Central";

		executeEventFromParameter(param);
	}

}

void opt_pf_electron::executeEventFromParameter(AnalyzerParameter param){

	if(!PassMETFilter()) return;

	Event ev = GetEvent();
	Particle METv = ev.GetMETVector();

	//==== Trigger ====
	if (! ev.PassTrigger(TrigNames)) return;

	//==== Copy All objects ====
	vector<Muon> this_AllMuons = AllMuons;
	vector<Electron> this_AllElectrons = AllElectrons;
	vector<Jet> this_AllJets = AllJets;

	//==== ID Selection ====
	vector<Muon> muons_veto = SelectMuons(this_AllMuons, MuonVetoID, 10., 2.4);
	vector<Electron> electrons_veto = SelectElectrons(this_AllElectrons, ElecVetoID, 10, 2.5);
	vector<Electron> electrons = SelectElectrons(this_AllElectrons, ElectronID, 25, 2.5);
	vector<Jet> jets = SelectJets(this_AllJets, JetID, 25, 2.4);
	vector<Jet> jets_dR04 = JetsVetoLeptonInside(jets, electrons_veto, muons_veto, 0.4);
	vector<Jet> jets_dR10 = JetsVetoLeptonInside(jets, electrons_veto, muons_veto, 1.0);

	sort(muons_veto.begin(), muons_veto.end(), PtComparing);
	sort(electrons_veto.begin(), electrons_veto.end(), PtComparing);
	sort(jets.begin(), jets.end(), PtComparing);
	sort(jets_dR04.begin(), jets_dR04.end(), PtComparing);
	sort(jets_dR10.begin(), jets_dR10.end(), PtComparing);

	//==== preselection ====
	if (! (electrons.size() == 1 || electrons.size() == 2)) return;
	if (! (muons_veto.size() == 0)) return;
	if (! (electrons_veto.size() <= 2)) return;
	if (! (electrons.at(0).Pt() > TriggerSafePtCut)) return;

	bool isPrompt = IsPromptEvent(electrons, electrons_veto, jets_dR04, METv);
	bool isFake = IsFakeEvent(electrons, electrons_veto, jets_dR10, METv);
	bool isOnZ = IsOnZEvent(electrons, jets_dR04);
	double weight = 1.;

	TString path = param.Name + "/";

	if (isPrompt) {
		TString path_prompt = path + "PromptEvents/";
		DrawHists(path_prompt, electrons.at(0), 0, weight);
		DrawHists(path_prompt, jets_dR04.at(0), 0, weight);
		DrawHists(path_prompt, METv, electrons.at(0), weight);
	}

	if (isFake) {
		TString path_fake = path + "FakeEvents/";
		DrawHists(path_fake, electrons.at(0), 0, weight);
		DrawHists(path_fake, jets_dR10.at(0), 0, weight);
		DrawHists(path_fake, METv, electrons.at(0), weight);
	}

	if (isOnZ) {
		TString path_onZ = path + "OnZEvents/";
		Particle ZCand = electrons.at(0) + electrons.at(1);
		DrawHists(path_onZ, electrons.at(0), 0, weight);
		DrawHists(path_onZ, electrons.at(1), 1, weight);
		DrawHists(path_onZ, jets_dR04.at(0), 0, weight);
		FillHist(path_onZ + "M(ee)", ZCand.M(), weight, 60, 60., 120.);
	}

	return;
}

bool opt_pf_electron::IsPromptEvent(const vector<Electron> &electrons, const vector<Electron> &electrons_veto, const vector<Jet> &jets, const Particle &METv) {
	//==== W boson dominated region
	//==== 1. Exactly 1 electron
	//==== 2. At least 1 jet with Pt > 40 GeV
	//==== 3. Delta(e, j) > 0.4
	//==== 4. MET > 50 GeV & MT(e, METv) > 70 GeV
	if (! (electrons.size() == 1)) return false;
	if (! (electrons_veto.size() == 1)) return false;
	if (! (jets.size() > 0)) return false;
	if (! (jets.at(0).Pt() > 40)) return false;
	for (const auto &j : jets) {
		if (j.DeltaR(electrons.at(0)) < 0.4) return false;
	}
	if (! (METv.Pt() > 50)) return false;
	double Mt = MT(electrons.at(0), METv);
	if (! (Mt > 70)) return false;

	return true;
}

bool opt_pf_electron::IsFakeEvent(const vector<Electron> &electrons, const vector<Electron> &electrons_veto, const vector<Jet> &jets, const Particle &METv) {
	//==== QCD dominated region
	//==== 1. Exactly 1 electron
	//==== 2. At least 1 jet with Pt > 40 GeV
	//==== 3. Delta(e, j) > 1.0
	//==== 4. MET < 25 GeV, MT(e, METv) < 25 GeV
	if (! (electrons.size() == 1)) return false;
	if (! (electrons_veto.size() == 1)) return false;
	if (! (jets.size() > 0)) return false;
	if (! (jets.at(0).Pt() > 40)) return false;
	for (const auto &j : jets) {
		if (j.DeltaR(electrons.at(0)) < 0.4) return false;
	}
	if (! (METv.Pt() < 25)) return false;
	double Mt = MT(electrons.at(0), METv);
	if (! (Mt < 25)) return false;

	return true;
}

bool opt_pf_electron::IsOnZEvent(const vector<Electron> &electrons, const vector<Jet> &jets) {
	//==== Z boson dominated region
	//==== 1. Exactly 2 electrons
	//==== 2, At least 1 jet with Pt > 40 GeV
	//==== 3. Delta(e, j ) > 0.4
	//==== 4. |M(ee) - 91.2| < 15
	if (! (electrons.size() == 2)) return false;
	if (! (jets.size() > 0)) return false;
	for (const auto &j : jets) {
		if (j.DeltaR(electrons.at(0)) < 0.4) return false;
		if (j.DeltaR(electrons.at(1)) < 0.4) return false;
	}
	Particle ZCand = electrons.at(0) + electrons.at(1);
	if (! IsOnZ(ZCand.M(), 15)) return false;
	
	return true;
}

void opt_pf_electron::DrawHists(TString path, const Electron &e, unsigned int order, const double weight) {
	// order guard
	if (order > 2) {
		cout << "[opt_pf_electron::DrawHists] order = " << order << endl;
		cout << "[opt_pf_electron::DrawHists] Wrond order" << endl;
	}

	TString this_path = path + "e" + TString::Itoa(order+1, 10) + "/";
	FillHist(this_path + "pt", e.Pt(), weight, 48, 0., 240.);
	FillHist(this_path + "eta", e.Eta(), weight, 50, -2.5, 2.5);

	this_path += "IDVariables/";
	FillHist(this_path + "RelIso", e.RelIso(), weight, 100, 0., 1.);
    FillHist(this_path + "dXY", e.dXY(), weight, 100, -0.5, 0.5);
    FillHist(this_path + "dXYerr", e.dXYerr(), weight, 100, -0.5, 0.5);
    if (e.dXYerr() != 0) FillHist(this_path + "SIP2D", e.dXY() / e.dXYerr(), weight, 200, -10, 10);
    FillHist(this_path + "dZ", e.dZ(), weight, 100, 0., 1.);

    FillHist(this_path + "Full5x5_sigmaIetaIeta", e.Full5x5_sigmaIetaIeta(), weight, 100, 0., 1.);
    FillHist(this_path + "dEtaSeed", e.dEtaSeed(), weight, 100, 0., 1.);
    FillHist(this_path + "dPhiIn", e.dPhiIn(), weight, 100, 0., 1.);
    FillHist(this_path + "HoverE", e.HoverE(), weight, 100, 0., 1.);
    FillHist(this_path + "InvEminusInvP", e.InvEminusInvP(), weight, 100, 0., 1.);
    FillHist(this_path + "e2x5OverE5x5", e.e2x5OverE5x5(), weight, 100, 0., 1.);
    FillHist(this_path + "e1x5OverE5x5", e.e1x5OverE5x5(), weight, 100, 0., 1.);
    FillHist(this_path + "trackIso", e.TrkIso(), weight, 100, 0., 1.);
    FillHist(this_path + "dr03EcalRecHitSumEt", e.dr03EcalRecHitSumEt(), weight, 100, 0., 1.);
    FillHist(this_path + "dr03HcalDepth1TowerSumEt", e.dr03HcalDepth1TowerSumEt(), weight, 100, 0., 1.);
    FillHist(this_path + "dr03HcalTowerSumEt", e.dr03HcalTowerSumEt(), weight, 100, 0., 1.);
    FillHist(this_path + "dr03TkSumPt", e.dr03TkSumPt(), weight, 100, 0., 1.);
    FillHist(this_path + "ecalPFClusterIso", e.ecalPFClusterIso(), weight, 100, 0., 1.);
    FillHist(this_path + "hcalPFClusterIso", e.hcalPFClusterIso(), weight, 100, 0., 1.);
}
void opt_pf_electron::DrawHists(TString path, const Jet &j, unsigned int order, const double weight) {
	
	TString this_path = path + "j" + TString::Itoa(order+1, 10) + "/";
	FillHist(this_path + "pt", j.Pt(), weight, 48, 0., 240.);
	FillHist(this_path + "eta", j.Eta(), weight, 50, -2.5, 2.5);
	FillHist(this_path + "phi", j.Phi(), weight, 80, -4., 4.);
}

void opt_pf_electron::DrawHists(TString path, const Particle &METv, const Electron &e, const double weight) {

	TString this_path = path + "MET/";

	double Mt = MT(e, METv);
	FillHist(this_path + "MET", METv.Pt(), weight, 48, 0., 240.);
	FillHist(this_path + "phi", METv.Phi(), weight, 80, -4., 4.);
	FillHist(this_path + "MT", Mt, weight, 48, 0., 240.);
}
	
opt_pf_electron::opt_pf_electron() {}

opt_pf_electron::~opt_pf_electron() {}


