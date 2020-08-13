#include "PromptSelector.h"

void PromptSelector::initializeAnalyzer(){
	// flags
	RunTTLL = HasFlag("RunTTLL");
	RunTTLJ = HasFlag("RunTTLJ");
	cout << "[PromptSelector::initializeAnalyzer] RunTTLL = " << RunTTLL << endl;
	cout << "[PromptSelector::initializeAnalyzer] RunTTLJ = " << RunTTLJ << endl;

	// b tagging
	vector<JetTagging::Parameters> jtps;
	jtps.push_back(JetTagging::Parameters(JetTagging::DeepCSV, JetTagging::Medium, JetTagging::incl, JetTagging::mujets));
	mcCorr->SetJetTaggingParameters(jtps);
}

void PromptSelector::executeEvent(){
	
	// copy all objects
	AllMuons = GetAllMuons();
	AllElectrons = GetAllElectrons();
	AllJets = GetAllJets();
	gens = GetGens();
	wegiht_Prefire = GetPrefireWeight(0);

	AnalyzerParameter param;
	param.Clear();
	param.syst_ = AnalyzerParameter::Central;

	vector<Muon> muons = SelectMuons(AllMuons, "POGLoose", 5., 2.4);
	vector<Muon> muons_veto = SelectMuons(AllMuons, "POGLoose", 5., 2.4);
	vector<Electron> electrons = SelectElectrons(AllElectrons, "passLooseID", 15., 2.5);
	vector<Electron> electrons_veto = SelectElectrons(AllElectrons, "passLooseID", 15., 2.5);
	vector<Jet> jets = SelectJets(AllJets, "tight", 20, 2.4);
	vector<Jet> bjets;
    for (const auto &j : jets) {
        double this_discr = j.GetTaggerResult(JetTagging::DeepCSV);
        if (this_discr > mcCorr->GetJetTaggingCutValue(JetTagging::DeepCSV, JetTagging::Medium))
            bjets.push_back(j);
    }
	vector<Jet> jets_dR04 = JetsVetoLeptonInside(jets, electrons_veto, muons_veto, 0.4);
	vector<Jet> bjets_dR04 = JetsVetoLeptonInside(bjets, electrons_veto, muons_veto, 0.4);
	Event ev = GetEvent();
	Particle METv = ev.GetMETVector();
	Particle METv_xyCorr(pfMET_Type1_PhiCor_pt*TMath::Cos(pfMET_Type1_PhiCor_phi), pfMET_Type1_PhiCor_pt*TMath::Sin(pfMET_Type1_PhiCor_phi), 0., pfMET_Type1_PhiCor_pt);

	sort(muons.begin(), muons.end(), PtComparing);
	sort(muons_veto.begin(), muons_veto.end(), PtComparing);
	sort(electrons.begin(), electrons.end(), PtComparing);
	sort(electrons_veto.begin(), electrons_veto.end(), PtComparing);
	sort(jets.begin(), jets.end(), PtComparing);
	sort(jets_dR04.begin(), jets_dR04.end(), PtComparing);
	sort(bjets.begin(), bjets.end(), PtComparing);
	sort(bjets_dR04.begin(), bjets_dR04.end(), PtComparing);

	// preselection
	TString path = "preselection/";
	FillHist(path + "cutflow", 0., 1., 9, 0., 9.);
	if (!PassMETFilter()) return;
	FillHist(path + "cutflow", 1., 1., 9, 0., 9.);
	if (muons.size() > 2) return;
	if (electrons.size() > 2) return;
	FillHist(path + "cutflow", 2., 1., 9, 0., 9.);
	if (RunTTLL && (jets_dR04.size() < 2)) return;
	if (RunTTLJ && (jets_dR04.size() < 3)) return;
	if (bjets_dR04.size() < 1) return;
	FillHist(path + "cutflow", 3., 1., 9, 0., 9.);

	// get lepton type
	vector<Muon> muons_fake, muons_EWPrompt, muons_SigDaughter, muons_EWTauDaughter, muons_others;
	vector<Electron> electrons_fake, electrons_EWPrompt, electrons_SigDaughter, electrons_EWTauDaughter, electrons_others;
    int lep_type;
    for (const auto mu: muons) {
    	lep_type = GetLeptonType(mu, gens);
    	if (lep_type < 0) muons_fake.push_back(mu);
        else if (lep_type == 1) muons_EWPrompt.push_back(mu);
    	else if (lep_type == 2) muons_SigDaughter.push_back(mu);
        else if (lep_type == 3) muons_EWTauDaughter.push_back(mu);
        else muons_others.push_back(mu);
    }
	for (const auto e: electrons) {
        lep_type = GetLeptonType(e, gens);
        if (lep_type < 0) electrons_fake.push_back(e);
        else if (lep_type == 1) electrons_EWPrompt.push_back(e);
        else if (lep_type == 2) electrons_SigDaughter.push_back(e);
        else if (lep_type == 3) electrons_EWTauDaughter.push_back(e);
        else electrons_others.push_back(e);
    }
	
	// weight
	double weight = 1.;
	if (!IsDATA) {
		weight *= weight_norm_1invpb*ev.GetTriggerLumi("Full");
		weight *= ev.MCweight();
		weight *= wegiht_Prefire;
	}

	// TTLL
	if (RunTTLL) {
		bool DIMU = isDiMu(muons, muons_veto, electrons_veto);
		bool DIELEC = isDiElec(electrons, electrons_veto, muons_veto);
		bool EMU = isEMu(electrons, muons, electrons_veto, muons_veto);
		
		if (DIMU) {
			path = "dimu/";
			DrawHists(path, muons, gens, weight, ALL);
			DrawHists(path, muons_fake, gens, weight, FAKE);
			DrawHists(path, muons_EWPrompt, gens, weight, EWPROMPT);
			DrawHists(path, muons_SigDaughter, gens, weight, SIGDAUGHTER);
			DrawHists(path, muons_EWTauDaughter, gens, weight, EWTAUDAUGHTER);
			DrawHists(path, muons_others, gens, weight, OTHERS);
            DrawHists(path, jets_dR04, weight, JET);
            DrawHists(path, bjets_dR04, weight, BJET);
            DrawHists(path, METv_xyCorr, weight);
		}
		
		if (DIELEC) {
			path = "dielec/";
			DrawHists(path, electrons, gens, weight, ALL);
            DrawHists(path, electrons_fake, gens, weight, FAKE);
            DrawHists(path, electrons_EWPrompt, gens, weight, EWPROMPT);
            DrawHists(path, electrons_SigDaughter, gens, weight, SIGDAUGHTER);
            DrawHists(path, electrons_EWTauDaughter, gens, weight, EWTAUDAUGHTER);
            DrawHists(path, electrons_others, gens, weight, OTHERS);
			DrawHists(path, jets_dR04, weight, JET);
			DrawHists(path, bjets_dR04, weight, BJET);
			DrawHists(path, METv_xyCorr, weight);
		}

		if (EMU) {
			path = "emu/";
			DrawHists(path, muons, gens, weight, ALL);
            DrawHists(path, muons_fake, gens, weight, FAKE);
            DrawHists(path, muons_EWPrompt, gens, weight, EWPROMPT);
            DrawHists(path, muons_SigDaughter, gens, weight, SIGDAUGHTER);
            DrawHists(path, muons_EWTauDaughter, gens, weight, EWTAUDAUGHTER);
            DrawHists(path, muons_others, gens, weight, OTHERS);
			DrawHists(path, electrons, gens, weight, ALL);
            DrawHists(path, electrons_fake, gens, weight, FAKE);
            DrawHists(path, electrons_EWPrompt, gens, weight, EWPROMPT);
            DrawHists(path, electrons_SigDaughter, gens, weight, SIGDAUGHTER);
            DrawHists(path, electrons_EWTauDaughter, gens, weight, EWTAUDAUGHTER);
            DrawHists(path, electrons_others, gens, weight, OTHERS);
			DrawHists(path, jets_dR04, weight, JET);
			DrawHists(path, bjets_dR04, weight, BJET);
			DrawHists(path, METv_xyCorr, weight);
		}
	}

	// TTLJ
	if (RunTTLJ) {
		bool EJJ = isElecJJ(electrons, electrons_veto, muons_veto);
		bool MUJJ = isMuJJ(muons, muons_veto, electrons_veto);
		
		if (EJJ) {
			path = "ejj/";
			DrawHists(path, electrons, gens, weight, ALL);
            DrawHists(path, electrons_fake, gens, weight, FAKE);
            DrawHists(path, electrons_EWPrompt, gens, weight, EWPROMPT);
            DrawHists(path, electrons_SigDaughter, gens, weight, SIGDAUGHTER);
            DrawHists(path, electrons_EWTauDaughter, gens, weight, EWTAUDAUGHTER);
            DrawHists(path, electrons_others, gens, weight, OTHERS);
            DrawHists(path, jets_dR04, weight, JET);
            DrawHists(path, bjets_dR04, weight, BJET);
            DrawHists(path, METv_xyCorr, weight);
		}
		if (MUJJ) {
			path = "mujj/";
			DrawHists(path, muons, gens, weight, ALL);
            DrawHists(path, muons_fake, gens, weight, FAKE);
            DrawHists(path, muons_EWPrompt, gens, weight, EWPROMPT);
            DrawHists(path, muons_SigDaughter, gens, weight, SIGDAUGHTER);
            DrawHists(path, muons_EWTauDaughter, gens, weight, EWTAUDAUGHTER);
            DrawHists(path, muons_others, gens, weight, OTHERS);
			DrawHists(path, jets_dR04, weight, JET);
            DrawHists(path, bjets_dR04, weight, BJET);
            DrawHists(path, METv_xyCorr, weight);
		}
	}	
}

PromptSelector::PromptSelector(){}

PromptSelector::~PromptSelector(){}

// member functions
bool PromptSelector::isDiMu(const vector<Muon> &muons, const vector<Muon> &muons_veto, const vector<Electron> &electrons_veto) {
	TString path = "dimu/cutflow";
	FillHist(path, 4., 1., 9, 0., 9.);
	
	if (! (muons.size() == 2)) return false;
	FillHist(path, 5., 1., 9, 0., 9.);
	if (! (muons.at(0).Charge()*muons.at(1).Charge() < 0)) return false;
	FillHist(path, 6., 1., 9, 0., 9.);
	if (! (muons_veto.size() == 2)) return false;
	if (! (electrons_veto.size() == 0)) return false;
	FillHist(path, 7., 1., 9, 0., 9.);

	return true;
}

bool PromptSelector::isDiElec(const vector<Electron> &electrons, const vector<Electron> &electrons_veto, const vector<Muon> &muons_veto) {
	TString path = "dielec/cutflow";
	FillHist(path, 4., 1., 9, 0., 9.);

	if (! (electrons.size() == 2)) return false;
	FillHist(path, 5., 1., 9, 0., 9.);
	if (! (electrons.at(0).Charge()*electrons.at(1).Charge() < 0)) return false;
	FillHist(path, 6., 1., 9, 0., 9.);
	if (! (electrons_veto.size() == 2)) return false;
	if (! (muons_veto.size() == 0)) return false;
	FillHist(path, 7., 1., 9, 0., 9.);

	return true;
}

bool PromptSelector::isEMu(const vector<Electron> &electrons, const vector<Muon> &muons, const vector<Electron> &electrons_veto, const vector<Muon> &muons_veto) {
	TString path = "emu/cutflow";
	FillHist(path, 4., 1., 9, 0., 9.);

	if (! (muons.size()==1 && electrons.size()==1)) return false;
	FillHist(path, 5., 1., 9, 0., 9.);
	if (! (muons.at(0).Charge()*electrons.at(0).Charge() < 0)) return false;
	FillHist(path, 6., 1., 9, 0., 9.);
	if (! (muons_veto.size() == 1)) return false;
	if (! (electrons_veto.size() == 1)) return false;
	FillHist(path, 7., 1., 9, 0., 9.);

	return true;
}

bool PromptSelector::isMuJJ(const vector<Muon>& muons, const vector<Muon>& muons_veto, const vector<Electron>& electrons_veto) {
	TString path = "mujj/cutflow";
	FillHist(path, 4., 1., 9, 0., 9.);

	if (! (muons.size()==1)) return false;
	FillHist(path, 5., 1., 9, 0., 9.);
	if (! (muons_veto.size() == 1)) return false;
	if (! (electrons_veto.size() == 0)) return false;
	FillHist(path, 6., 1., 9, 0., 9.);

	return true;
}

bool PromptSelector::isElecJJ(const vector<Electron>& electrons, const vector<Electron>& electrons_veto, const vector<Muon>& muons_veto) {
	TString path = "ejj/cutflow";
	FillHist(path, 4., 1., 9, 0., 9.);

	if (! (electrons.size()==1)) return false;
	FillHist(path, 5., 1., 9, 0., 9.);
	if (! (electrons_veto.size() == 1)) return false;
	if (! (muons_veto.size() == 0)) return false;
	FillHist(path, 6., 1., 9, 0., 9.);
	
	return true;
}
// DrawHists
void PromptSelector::DrawHists(TString path, const vector<Electron>& electrons, const vector<Gen>& gens, const double& weight, PromptSelector::LEPTYPE TYPE){
	TString this_path;
	if (TYPE == ALL) this_path = path + "electrons/";
	else if (TYPE == FAKE) this_path = path + "electrons/fake/";
	else if (TYPE == EWPROMPT) this_path = path + "electrons/EWPrompt/";
	else if (TYPE == SIGDAUGHTER) this_path = path + "electrons/SignalDaughter/";
	else if (TYPE == EWTAUDAUGHTER) this_path = path + "electrons/EWTauDaughter/";
	else if (TYPE == OTHERS) this_path = path + "electrons/others/";
	
	TString obj_path;
	int lep_type = 0;
	for (unsigned int i = 0 ; i < electrons.size(); i++) {
		obj_path = this_path + TString::Itoa(i+1, 10) + "/";
		FillHist(obj_path + "pt", electrons.at(i).Pt(), weight, 300, 0., 300.);
		FillHist(obj_path + "eta", electrons.at(i).Eta(), weight, 50, -2.5, 2.5);
		FillHist(obj_path + "phi", electrons.at(i).Phi(), weight, 70, -3.5, 3.5);
		if (TYPE == ALL) {
			lep_type = GetLeptonType(electrons.at(i), gens);
			FillHist(obj_path + "lepton_type", lep_type, weight, 14, -7, 7);
		}
	}
}
void PromptSelector::DrawHists(TString path, const vector<Muon>& muons, const vector<Gen>& gens, const double& weight, PromptSelector::LEPTYPE TYPE) {
	TString this_path;
    if (TYPE == ALL) this_path = path + "muons/";
    else if (TYPE == FAKE) this_path = path + "muons/fake/";
    else if (TYPE == EWPROMPT) this_path = path + "muons/EWPrompt/";
    else if (TYPE == SIGDAUGHTER) this_path = path + "muons/SignalDaughter/";
    else if (TYPE == EWTAUDAUGHTER) this_path = path + "muons/EWTauDaughter/";
    else if (TYPE == OTHERS) this_path = path + "muons/others/";

	TString obj_path;
	int lep_type = 0;
	for (unsigned int i = 0; i < muons.size(); i++) {
		obj_path = this_path + TString::Itoa(i+1, 10) + "/";
		FillHist(obj_path + "pt", muons.at(i).Pt(), weight, 300, 0., 300.);
		FillHist(obj_path + "eta", muons.at(i).Eta(), weight, 48, -2.4, 2.4);
		FillHist(obj_path + "phi", muons.at(i).Phi(), weight, 70, -3.5, 3.5);
		if (TYPE == ALL) {
			lep_type = GetLeptonType(muons.at(i), gens);
			FillHist(obj_path + "lepton_type", lep_type, weight, 14, -7, 7);
		}
	}
}
void PromptSelector::DrawHists(TString path, const vector<Jet>& jets, const double& weight, const bool& isBJET) {
	TString this_path;
	if (isBJET) this_path = path + "bjets/";
	else this_path = path + "jets/";
	FillHist(this_path + "size", jets.size(), weight, 14, 0., 14);

	TString obj_path;
	for (unsigned int i = 0; i < jets.size(); i++) {
		obj_path = this_path + TString::Itoa(i+1, 10) + "/";
		FillHist(obj_path + "pt", jets.at(i).Pt(), weight, 300, 0., 300.);
		FillHist(obj_path + "eta", jets.at(i).Eta(), weight, 48, -2.4, 2.4);
		FillHist(obj_path + "phi", jets.at(i).Phi(), weight, 70, -3.5, 3.5);
	}
}
void PromptSelector::DrawHists(TString path, const Particle& METv, const double& weight) {
	TString this_path = path + "METv/";
	FillHist(this_path + "pt", METv.Pt(), weight, 300, 0., 300.);
	FillHist(this_path + "phi", METv.Phi(), weight, 70, -3.5, 3.5);
}
