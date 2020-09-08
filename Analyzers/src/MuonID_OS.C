#include "MuonID_OS.h"


// To do::
// 1. trigger SF for POGTight - 2016
// 2. need check of trigger SF
// 3. complete signal selection

void MuonID_OS::initializeAnalyzer(){

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
		cerr << "[MuonID_OS::initializeAnalyzer] DataYear is wrong : " << DataYear << endl;
		exit(EXIT_FAILURE);
	}
	
	// ID settings
	MuonID_OSs = {"POGLoose", "POGTight"};

	// B-tagging
	vector<JetTagging::Parameters> jtps;
	jtps.push_back(
			JetTagging::Parameters(JetTagging::DeepCSV, JetTagging::Medium, JetTagging::incl, JetTagging::comb));
	mcCorr->SetJetTaggingParameters(jtps);
}

void MuonID_OS::executeEvent(){

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
	fatjets = GetAllFatJets();
	jets = GetAllJets();

	sort(muons.begin(), muons.end(), PtComparing);
	sort(electrons.begin(), electrons.end(), PtComparing);
	sort(fatjets.begin(), fatjets.end(), PtComparing);
	sort(jets.begin(), jets.end(), PtComparing);

	//muons = MuonPromptOnlyHNtype1(muons); // gen matching applied for MC
	muons_veto = SelectMuons(muons, "HNType1_POGVeto", 5., 2.4);
	muons_loose = SelectMuons(muons_veto, "HNType1_POGLoose", 10., 2.4);
	muons_tight = SelectMuons(muons_loose, "HNType1_POGTight", 10., 2.4);
	//muons_loose = SelectMuons(muons_veto, "HNType1_POGLoose", 10., 2.4);
	//muons_tight = SelectMuons(muons_loose, "HNType1_POGTight", 10, 2.4);

	electrons_veto = SelectElectrons(electrons, "HNType1_CutBasedVeto", 10, 2.5);

	fatjets_tight = SelectFatJets(fatjets, "tightWithSDMass", 200., 2.7);
	fatjets_dR10 = FatJetsVetoLeptonInside(fatjets_tight, electrons_veto, muons_veto, 1.0);

	jets_tight = SelectJets(jets, "tight", 20, 2.7);
	jets_dR04 = JetsVetoLeptonInside(jets_tight, electrons_veto, muons_veto, 0.4);
	jets_awayFatJet = JetsAwayFromFatJet(jets_dR04, fatjets_dR10, 0.8);

	// can't use jets outside |eta| > 2.4 for b-taggings
	// no lepton veto for bjets
	bjets = SelectJets(jets_tight, "tight", 20., 2.4);
	for (const auto& j: bjets) {
		double this_discr = j.GetTaggerResult(JetTagging::DeepCSV);
		if (this_discr >  mcCorr->GetJetTaggingCutValue(JetTagging::DeepCSV, JetTagging::Medium)) 
			bjets_tight.push_back(j);
	}
	//bjets_dR04 = JetsVetoLeptonInside(bjets_tight, electrons_veto, muons_veto, 0.4);
	
	// preselection
	// 1. pass trigger
	// 2. 2 os tight leptons. no 3rd lepton
	// 3. truth matching(?)
	// 4. muons should above pt threshold
	// 5. m(ll) > 10 GeV
	// 6. at least one jet or fatjet

	// Let's check muons_loose distributions before and after trigger
	const double weight = 1.; // no change for data
	DrawHists("beforeTrig/muons_loose/", muons_loose, weight);
	DrawHists("beforeTrig/muons_tight/", muons_tight, weight);
	// pass trigger
	bool passTrigs = ev.PassTrigger(trigs_POGTight);
	if (!passTrigs) {
		DrawHists("failTrig/muons_loose/", muons_loose, weight);
		DrawHists("failTrig/muons_tight/", muons_tight, weight);
		return;
	}
	FillHist("preselection/cutflow", 2., 1., 10, 0., 10.);
	DrawHists("afterTrig/muons_loose/", muons_loose, weight);
	DrawHists("afterTrig/muons_tight/", muons_tight, weight);

	// for fake estimation, loose id first
	if (! (muons_loose.size() == 2)) return;
	if (! (muons_loose.at(0).Charge()*muons_loose.at(1).Charge() < 0)) return;
	if (! (muons_veto.size() == 2)) return;
	if (! (electrons_veto.size() == 0)) return;
	FillHist("preselection/cutflow", 3., 1., 10, 0., 10.);

	if (! (muons_loose.at(0).Pt() > 20.)) return;
	if (! (muons_loose.at(1).Pt() > 10.)) return;
	FillHist("preselection/cutflow", 4., 1., 10, 0., 10.);

	Particle reso = muons_loose.at(0) + muons_loose.at(1);
	if (! (reso.M() > 10)) return;
	FillHist("preselection/cutflow", 5., 1., 10, 0., 10.);
	
	//if (jets_awayFatJet.size()==0 && fatjets_dR10.size()==0) return;
	//FillHist("preselection/cutflow", 6., 1., 10, 0., 10.);

	TString path;

	// objects distribution
	
	path = "preselection/loose/";
	DrawHists(path + "muons_loose/", muons_loose, weight);
	DrawHists(path + "muons_tight/", muons_tight, weight);
	DrawHists(path + "jets_awayFatJet/", jets_awayFatJet, weight);
	DrawHists(path + "bjets_tight/", bjets_tight, weight);
	DrawHists(path + "fatjets_dR10/", fatjets_dR10, weight);
	DrawHists(path + "METv/", METv, weight);

	// Categorization
	const Muon& lead = muons_loose.at(0);
	const Muon& tail = muons_loose.at(1);
	IDFlag flag;
	if (lead.PassID("HNType1_POGTight")) {
		if (tail.PassID("HNType1_POGTight")) flag = TT;
		else flag = TL;
	}
	else {
		if (tail.PassID("HNType1_POGTight")) flag = LT;
		else flag = LL;
	}

	if (flag == TT) {
		path = "preselection/TT/";
		DrawHists(path + "muons_loose/", muons_loose, weight);
		DrawHists(path + "muons_tight/", muons_tight, weight);
		DrawHists(path + "jets_awayFatJet/", jets_awayFatJet, weight);
    	DrawHists(path + "bjets_tight/", bjets_tight, weight);
    	DrawHists(path + "fatjets_dR10/", fatjets_dR10, weight);
    	DrawHists(path + "METv/", METv, weight);
	}
	else if (flag == TL) {
        path = "preselection/TL/";
        DrawHists(path + "muons_loose/", muons_loose, weight);
		DrawHists(path + "muons_tight/", muons_tight, weight);
		DrawHists(path + "jets_awayFatJet/", jets_awayFatJet, weight);
    	DrawHists(path + "bjets_tight/", bjets_tight, weight);
    	DrawHists(path + "fatjets_dR10/", fatjets_dR10, weight);
    	DrawHists(path + "METv/", METv, weight);
    }
    else if (flag == LT) {
        path = "preselection/LT/";
        DrawHists(path + "muons_loose/", muons_loose, weight);
		DrawHists(path + "muons_tight/", muons_tight, weight);
		DrawHists(path + "jets_awayFatJet/", jets_awayFatJet, weight);
    	DrawHists(path + "bjets_tight/", bjets_tight, weight);
    	DrawHists(path + "fatjets_dR10/", fatjets_dR10, weight);
    	DrawHists(path + "METv/", METv, weight);
    }
    else if (flag == LL) {
        path = "preselection/LL/";
        DrawHists(path + "muons_loose/", muons_loose, weight);
		DrawHists(path + "muons_tight/", muons_tight, weight);
		DrawHists(path + "jets_awayFatJet/", jets_awayFatJet, weight);
    	DrawHists(path + "bjets_tight/", bjets_tight, weight);
    	DrawHists(path + "fatjets_dR10/", fatjets_dR10, weight);
    	DrawHists(path + "METv/", METv, weight);
    }
    else {
        cerr << "[HNSignal::ExecuteEvent] WTF" << endl;
        exit(EXIT_FAILURE);
    }

}

// ==== other functions ====
void MuonID_OS::ClearCollections() {
	muons.clear(); muons_tight.clear(); muons_loose.clear(), muons_veto.clear();
	electrons.clear(); electrons_veto.clear();
	jets.clear(); jets_tight.clear(); jets_dR04.clear(); jets_awayFatJet.clear();
	bjets.clear(); bjets_tight.clear(); bjets_dR04.clear();
	fatjets.clear(); fatjets_tight.clear(); fatjets_dR10.clear();
	gens.clear();
}

bool MuonID_OS::isSignal(const SR& sig) const {
	switch (sig) {
	case SR1:
		if (! (jets_dR04.size() >= 2)) return false;
		if (! (fatjets_dR10.size() == 0)) return false;
		return true;
		break;
	case SR2:
		if (! (fatjets_dR10.size() > 0)) return false;
		return true;
		break;
	default:
		cerr << "[MuonID_OS::isSignal] Wrong Signal region... sig = " << sig << endl;
		exit(EXIT_FAILURE);
	}
}

bool MuonID_OS::isControl(const CR& con) const {
	cout << "[MuonID_OS::isControl] Control region has not been set yet" << endl;
	exit(EXIT_FAILURE);

	switch (con) {
	case CR1:
		if (! (jets_dR04.size() >= 2)) return false;
		if (! (fatjets_dR10.size() == 0)) return false;
		return true;
		break;
	case CR2:
		if (! (fatjets_dR10.size() > 0)) return false;
		return true;
		break;
	default:
		cerr << "[MuonID_OS::isControl] Wrong control region... con = " << con << endl;
		exit(EXIT_FAILURE);
	}
}


// need to set path before executing functions
void MuonID_OS::DrawHists(TString path, const vector<Muon>& muons, const double& weight) {
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

void MuonID_OS::DrawHists(TString path, const vector<Electron>& electrons, const double& weight) {
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
void MuonID_OS::DrawHists(TString path, const vector<Jet>& jets, const double& weight) {
	TString obj_path;
	FillHist(path + "size", jets.size(), weight, 14, 0., 14.);
	for (unsigned int i = 0; i < jets.size(); i++) {
		obj_path = path + TString::Itoa(i+1, 10) + "/";
		FillHist(obj_path + "pt", jets.at(i).Pt(), weight, 300, 0., 300.);
		FillHist(obj_path + "eta", jets.at(i).Eta(), weight, 48, -2.4, 2.4);
		FillHist(obj_path + "phi", jets.at(i).Phi(), weight, 70, -3.5, 3.5);
	}	
}
void MuonID_OS::DrawHists(TString path, const vector<FatJet>& jets, const double& weight) {
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
void MuonID_OS::DrawHists(TString path, const Particle& METv, const double& weight) {
	FillHist(path + "pt", METv.Pt(), weight, 300, 0., 300.);
	FillHist(path + "eta", METv.Eta(), weight, 48, -2.4, 2.4); // of course, 0.
	FillHist(path + "phi", METv.Phi(), weight, 70, -3.5, 3.5);
}

MuonID_OS::MuonID_OS(){}

MuonID_OS::~MuonID_OS(){}


