#include "WWW_SS.h"


// To do::
// 1. trigger SF for POGTight - 2016
// 2. need check of trigger SF
// 3. complete signal selection

void WWW_SS::initializeAnalyzer(){

	// triggers
	// only need doublemuon triggers
	if (DataYear == 2016) {
		trigs_POGTight = {
			"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v",		// BCDEF
      		"HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v",		// BCDEF
      		"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v",	// GH
      		"HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v"	// GH
		};
	}
	else if (DataYear == 2017) {
		trigs_POGTight = {
			"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v"
    	};
	}
	else if (DataYear == 2018) {
		trigs_POGTight = {
			"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v"
		};
	}
	else {
		cerr << "[WWW_SS::initializeAnalyzer] DataYear is wrong : " << DataYear << endl;
		exit(EXIT_FAILURE);
	}
	
	// B-tagging
	vector<JetTagging::Parameters> jtps;
	jtps.push_back(
			JetTagging::Parameters(JetTagging::DeepCSV, JetTagging::Medium, JetTagging::incl, JetTagging::comb));
	mcCorr->SetJetTaggingParameters(jtps);
}

void WWW_SS::executeEvent(){

	FillHist("preselection/cutflow", 0., 1., 10, 0., 10.); 
	// MET filter && trigger
	if (!PassMETFilter()) return;
	FillHist("preselection/cutflow", 1., 1., 10, 0., 10.);

	Event ev = GetEvent();
	METv = ev.GetMETVector();
	
	// objects definition
	ClearCollections();
	gens = GetGens();
	muons = GetAllMuons();
	electrons = GetAllElectrons();
	jets = GetAllJets();

	sort(muons.begin(), muons.end(), PtComparing);
	sort(electrons.begin(), electrons.end(), PtComparing);
	sort(jets.begin(), jets.end(), PtComparing);

	//muons = MuonPromptOnlyHNtype1(muons); // gen matching applied for MC
	muons_veto = SelectMuons(muons, "HNType1_POGVeto", 5., 2.4);
	muons_loose = SelectMuons(muons_veto, "HNType1_POGLoose", 10., 2.4);
	muons_tight = SelectMuons(muons_loose, "HNType1_POGTight", 10., 2.4);

	electrons_veto = SelectElectrons(electrons, "HNType1_CutBasedVeto", 10, 2.5);

	jets_tight = SelectJets(jets, "tight", 20, 2.4);
	jets_dR04 = JetsVetoLeptonInside(jets_tight, electrons_veto, muons_veto, 0.4);

	// can't use jets outside |eta| > 2.4 for b-taggings
	// no lepton veto for bjets
	bjets = SelectJets(jets_dR04, "tight", 20., 2.4);
	for (const auto& j: bjets) {
		double this_discr = j.GetTaggerResult(JetTagging::DeepCSV);
		if (this_discr >  mcCorr->GetJetTaggingCutValue(JetTagging::DeepCSV, JetTagging::Medium)) 
			bjets_dR04.push_back(j);
	}
	
	// preselection
	// 1. pass trigger
	// 2. 3 muons, include opposite-sign muon pair, no additional lepton
	// 3. pt > 20, 10, 10
	// 4. M(mumu) > 12 GeV for all os muon pair
	// 5. |M(mumu) - 91.2| > 10 GeV for all os pair
	// 6. Nj >= 2
	// 7. Nb >= 1

	// pass trigger
	bool passTrigs = ev.PassTrigger(trigs_POGTight);

	if (!passTrigs) return;
	FillHist("preselection/cutflow", 2., 1., 10, 0., 10.);

	// for fake estimation, loose id first
	if (! (muons_loose.size() == 3)) return;
	if (! (fabs(muons_loose.at(0).Charge() + muons_loose.at(1).Charge() + muons_loose.at(2).Charge()) == 1)) return;
	if (! (muons_veto.size() == 3)) return;
	if (! (electrons_veto.size() == 0)) return;
	FillHist("preselection/cutflow", 3., 1., 10, 0., 10.);

	if (! (muons_loose.at(0).Pt() > 20.)) return;
	if (! (muons_loose.at(1).Pt() > 10.)) return;
	if (! (muons_loose.at(2).Pt() > 10.)) return;
	FillHist("preselection/cutflow", 4., 1., 10, 0., 10.);

	for (unsigned int i = 0; i < muons_loose.size(); i++) {
		const Muon& first = muons_loose.at(i%3);
		const Muon& second = muons_loose.at(i%3);

		//SS, continue
		if (first.Charge()*second.Charge() > 0) continue;
		
		//OS, should above 12 GeV
		const Particle Reso = first + second;
		if (! (Reso.M() > 12)) return;
		if (! (fabs(Reso.M() - 91.2) > 10)) return;
	}
	FillHist("preselection/cutflow", 5., 1., 10, 0., 10.);
	
	if (! (jets_dR04.size() >= 2)) return;
	FillHist("preselection/cutflow", 6., 1., 10, 0., 10.);
	if (! (bjets_dR04.size() >= 1)) return;
	FillHist("preselection/cutflow", 7., 1., 10, 0., 10.);

	TString path = "fuck you";

	// objects distribution
	const double weight = 1.;

	// Categorization
	unsigned nTight = 0;
	for (auto& mu : muons_loose) {
		if (mu.PassID("HNType1_POGTight")) nTight++;
	}

	if (nTight == 0) {
		path = "preselection/3L0T/";
		DrawHists(path + "muons_loose/", muons_loose, weight);
        DrawHists(path + "muons_tight/", muons_tight, weight);
        DrawHists(path + "jets_dR04/", jets_dR04, weight);
        DrawHists(path + "bjets_dR04/", bjets_dR04, weight);
        DrawHists(path + "METv/", METv, weight);	
	}

	else if (nTight == 1) {
		path = "preselection/2L1T/";
		DrawHists(path + "muons_loose/", muons_loose, weight);
        DrawHists(path + "muons_tight/", muons_tight, weight);
        DrawHists(path + "jets_dR04/", jets_dR04, weight);
        DrawHists(path + "bjets_dR04/", bjets_dR04, weight);
        DrawHists(path + "METv/", METv, weight);
	}

	else if (nTight == 2) {
		path = "preselection/1L2T/";
		DrawHists(path + "muons_loose/", muons_loose, weight);
        DrawHists(path + "muons_tight/", muons_tight, weight);
        DrawHists(path + "jets_dR04/", jets_dR04, weight);
        DrawHists(path + "bjets_dR04/", bjets_dR04, weight);
        DrawHists(path + "METv/", METv, weight);
	}
	
	else if (nTight == 3) {
		path = "preselection/0L3T/";
		DrawHists(path + "muons_loose/", muons_loose, weight);
        DrawHists(path + "muons_tight/", muons_tight, weight);
        DrawHists(path + "jets_dR04/", jets_dR04, weight);
        DrawHists(path + "bjets_dR04/", bjets_dR04, weight);
        DrawHists(path + "METv/", METv, weight);
	}
	else {
		cerr << "[WWW_SS] Sth wrong with muonID" << endl;
		exit(EXIT_FAILURE);
	}

}

// ==== other functions ====
void WWW_SS::ClearCollections() {
	muons.clear(); muons_tight.clear(); muons_loose.clear(), muons_veto.clear();
	electrons.clear(); electrons_veto.clear();
	jets.clear(); jets_tight.clear(); jets_dR04.clear();
	bjets.clear(); bjets_dR04.clear();
	gens.clear();
}

bool WWW_SS::isSignal(const SR& sig) const {
	switch (sig) {
	case SR1:
		if (! (jets_dR04.size() >= 2)) return false;
		return true;
		break;
	case SR2:
		return true;
		break;
	default:
		cerr << "[WWW_SS::isSignal] Wrong Signal region... sig = " << sig << endl;
		exit(EXIT_FAILURE);
	}
}

bool WWW_SS::isControl(const CR& con) const {
	cout << "[WWW_SS::isControl] Control region has not been set yet" << endl;
	exit(EXIT_FAILURE);

	switch (con) {
	case CR1:
		if (! (jets_dR04.size() >= 2)) return false;
		return true;
		break;
	case CR2:
		return true;
		break;
	default:
		cerr << "[WWW_SS::isControl] Wrong control region... con = " << con << endl;
		exit(EXIT_FAILURE);
	}
}


// need to set path before executing functions
void WWW_SS::DrawHists(TString path, const vector<Muon>& muons, const double& weight) {
	TString obj_path;
	for (unsigned int i = 0; i < muons.size(); i++) {
		obj_path = path + TString::Itoa(i+1, 10) + "/";
		FillHist(obj_path + "pt", muons.at(i).Pt(), weight, 300, 0., 300.);
		FillHist(obj_path + "eta", muons.at(i).Eta(), weight, 48, -2.4, 2.4);
		FillHist(obj_path + "phi", muons.at(i).Phi(), weight, 70, -3.5, 3.5);
		double RelIso = muons.at(i).RelIso();
		if (RelIso > 0.5) RelIso = 0.49;
		FillHist(obj_path + "RelIso", RelIso, weight, 50, 0., 0.5);
		double TrkIso = muons.at(i).TrkIso() / muons.at(i).MiniAODPt();
		if (TrkIso > 0.8) TrkIso = 0.79;
		FillHist(obj_path + "TrkIso", TrkIso, weight, 80, 0., 0.8);
		FillHist(obj_path + "dXY", fabs(muons.at(i).dXY()), weight, 50, 0., 0.5);
		FillHist(obj_path + "dZ", fabs(muons.at(i).dZ()), weight, 80, 0., 0.8);

		if (muons.at(i).dXYerr() != 0.) {
			double SIP2D = fabs(muons.at(i).dXY() / muons.at(i).dXYerr());
			if (SIP2D > 10.) SIP2D = 9.99;
			FillHist(obj_path + "SIP2D", SIP2D, weight, 100, 0., 10.);
		}
		if (muons.at(i).IP3Derr() != 0.) {
			double SIP3D = fabs(muons.at(i).IP3D() / muons.at(i).IP3Derr());
			if (SIP3D > 10.) SIP3D = 9.99;
			FillHist(obj_path + "SIP3D", SIP3D, weight, 100, 0., 10.);
		}
	}
}

void WWW_SS::DrawHists(TString path, const vector<Electron>& electrons, const double& weight) {
	TString obj_path;
	for (unsigned int i = 0; i < electrons.size(); i++) {
		obj_path = path + TString::Itoa(i+1, 10) + "/";
		FillHist(obj_path + "pt", electrons.at(i).Pt(), weight, 300, 0., 300.);
		FillHist(obj_path + "eta", electrons.at(i).Eta(), weight, 50, -2.5, 2.5);
		FillHist(obj_path + "phi", electrons.at(i).Phi(), weight, 70, -3.5, 3.5);
		FillHist(obj_path + "RelIso", electrons.at(i).RelIso(), weight, 50, 0., 0.5);
		FillHist(obj_path + "dXY", fabs(electrons.at(i).dXY()), weight, 50, 0., 0.5);
		FillHist(obj_path + "dZ", fabs(electrons.at(i).dZ()), weight, 80, 0., 0.8);
	}
}
void WWW_SS::DrawHists(TString path, const vector<Jet>& jets, const double& weight) {
	TString obj_path;
	FillHist(path + "size", jets.size(), weight, 14, 0., 14.);
	for (unsigned int i = 0; i < jets.size(); i++) {
		obj_path = path + TString::Itoa(i+1, 10) + "/";
		FillHist(obj_path + "pt", jets.at(i).Pt(), weight, 300, 0., 300.);
		FillHist(obj_path + "eta", jets.at(i).Eta(), weight, 48, -2.4, 2.4);
		FillHist(obj_path + "phi", jets.at(i).Phi(), weight, 70, -3.5, 3.5);
	}	
}
void WWW_SS::DrawHists(TString path, const vector<FatJet>& jets, const double& weight) {
	TString obj_path;
	FillHist(path + "size", jets.size(), weight, 10, 0., 10.);
	for (unsigned int i = 0; i < jets.size(); i++) {
		obj_path = path + TString::Itoa(i+1, 10) + "/";
		FillHist(obj_path + "pt", jets.at(i).Pt(), weight, 500, 0., 500.);
		FillHist(obj_path + "eta", jets.at(i).Eta(), weight, 48, -2.4, 2.4);
		FillHist(obj_path + "phi", jets.at(i).Phi(), weight, 70, -3.5, 3.5);
		FillHist(obj_path + "SDMass", jets.at(i).SDMass(), weight, 300, 0., 300.);
		if (jets.at(i).PuppiTau1() != 0.) {
			double tau21 = jets.at(i).PuppiTau2() / jets.at(i).PuppiTau1();
			FillHist(obj_path + "tau21", tau21, weight, 100, 0., 1.);
		}
	}
}
void WWW_SS::DrawHists(TString path, const Particle& METv, const double& weight) {
	FillHist(path + "pt", METv.Pt(), weight, 300, 0., 300.);
	FillHist(path + "eta", METv.Eta(), weight, 48, -2.4, 2.4); // of course, 0.
	FillHist(path + "phi", METv.Phi(), weight, 70, -3.5, 3.5);
}

WWW_SS::WWW_SS(){}

WWW_SS::~WWW_SS(){}


