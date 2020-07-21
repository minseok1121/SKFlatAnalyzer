#include "id_var_prompt.h"

void id_var_prompt::initializeAnalyzer(){

	RunMuon = HasFlag("RunMuon");
	RunElectron = HasFlag("RunElectron");
	cout << "[id_var_prompt::initializeAnalyzer] RunMu = " << RunMuon << endl;
	cout << "[id_var_prompt::initializeAnalyzer] RunEl = " << RunElectron << endl;

	//==== b tagging ====
	vector<JetTagging::Parameters> jtps;
    jtps.push_back( JetTagging::Parameters(JetTagging::DeepCSV, JetTagging::Medium, JetTagging::incl, JetTagging::comb));
    mcCorr->SetJetTaggingParameters(jtps);	

}

void id_var_prompt::executeEvent(){

	//==== copy all objects ====
	AllMuons = GetAllMuons();
	AllElectrons = GetAllElectrons();
	AllJets = GetAllJets();
	Gens = GetGens();

	AnalyzerParameter param;
	param.Clear();
	param.syst_ = AnalyzerParameter::Central;

	//==== No Cut ====
	FillHist("NoCut", 0., 1., 5, 0., 5.);

	//==== MET filter ====
	if (!PassMETFilter()) return;
	FillHist("NoCut", 1., 1., 5, 0., 5.);

	
	//==== copy objects and sort ====
	vector<Muon> muons = SelectMuons(AllMuons, "NOCUT", 5., 2.5);
	vector<Muon> muons_veto = SelectMuons(AllMuons, "POGLoose", 20., 2.4);
    vector<Electron> electrons = SelectElectrons(AllElectrons, "NOCUT", 5., 2.5);
    vector<Electron> electrons_veto = SelectElectrons(AllElectrons, "passLooseID", 20., 2.5);
    vector<Jet> jets = SelectJets(AllJets, "tight", 25, 2.4);
    vector<Jet> bjets;
    JetTagging::Parameters jtp_DeepCSV_Medium
        = JetTagging::Parameters(JetTagging::DeepCSV, JetTagging::Medium, JetTagging::incl, JetTagging::comb);
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

	//==== PreSelection ====
	if (muons.size() > 2) return;
	if (electrons.size() > 2) return;
	FillHist("NoCut", 2., 1., 5, 0., 5.);
	FillHist("nMuons", muons.size(), 1., 3, 0., 3.);
	FillHist("nElectrons", electrons.size(), 1., 3, 0., 3.);
	FillHist("nJets", jets.size(), 1., 7, 0., 7.);
	FillHist("nJets_dR04", jets_dR04.size(), 1., 7, 0., 7.);
	FillHist("nBJets", bjets.size(), 1., 7, 0., 7.);
	FillHist("nBJets_dR04", bjets_dR04.size(), 1., 7, 0., 7.);

	//==== event weight ====
    weight = 1.;
    w_gen = 1.; w_kfactor = 1.; w_prefire = 1.;
    w_toppt = 1.; w_lumi = 1.; w_pileup = 1.;
    sf_trig = 1.; sf_mutk = 1.; sf_muid = 1.; sf_muiso = 1.;
    sf_elreco = 1.; sf_elid = 1.; sf_btag = 1.;

    if (!IsDATA) {
        w_gen = ev.MCweight();
        w_lumi = weight_norm_1invpb*ev.GetTriggerLumi("Full");
        w_prefire = GetPrefireWeight(0);
        if (MCSample.Contains("TT") && MCSample.Contains("powheg")) {
            w_toppt = mcCorr->GetTopPtReweight(Gens);
            sf_btag = mcCorr->GetBTaggingReweight_1a(jets_dR04, jtp_DeepCSV_Medium);
        }
        w_pileup = GetPileUpWeight(nPileUp, 0);
    }

	weight *= w_gen * w_kfactor * w_lumi * w_toppt * w_pileup * w_prefire;
    weight *= sf_trig * sf_mutk * sf_muid * sf_muiso * sf_elreco * sf_elid * sf_btag;
	
	//==== what do you want to observe? ====
	if (RunMuon) {
		param.Name = "muon_nocut_central";
		TString path = param.Name + "/";

		//==== event selection ====
		bool isDY = IsDY(muons, electrons_veto);
		bool isTT = IsTT(muons, electrons_veto, jets_dR04, bjets_dR04, METv_xyCorr);
		vector<int> gen_muon_type; // 1 for EWPrompt, 3 for tau

		if (isDY) {
			gen_muon_type.clear();
			for (const auto &mu : muons) {
				int this_type = GetLeptonType(mu, Gens);
				gen_muon_type.push_back(this_type);
			}
			
			path = param.Name + "/DY/";
			FillHist(path + "METv_pt", METv.Pt(), weight, 240, 0., 240.);
			FillHist(path + "METv_xyCorr_pt", METv_xyCorr.Pt(), weight, 240, 0., 240.); 
			FillHist(path + "METv_phi", METv.Phi(), weight, 70, -3.5, 3.5);
			FillHist(path + "METv_xyCorr_phi", METv_xyCorr.Phi(), weight, 70, -3.5, 3.5);
			FillHist(path + "Nj", jets.size(), weight, 8, 0., 8.);
			FillHist(path + "Nb", bjets.size(), weight, 8, 0., 8.);
			FillHist(path + "Nj_dR04", jets_dR04.size(), weight, 8, 0., 8.);
			FillHist(path + "Nb_dR04", bjets_dR04.size(), weight, 8, 0., 8.);
			DrawHists(path, muons.at(0), 0, gen_muon_type.at(0), weight);
			DrawHists(path, muons.at(1), 1, gen_muon_type.at(1), weight);
			if (jets_dR04.size() > 0) DrawHists(path, jets_dR04.at(0), 0, weight);
			if (jets_dR04.size() > 1) DrawHists(path, jets_dR04.at(1), 1, weight);
		}

		if (isTT) {
			gen_muon_type.clear();
			for (const auto &mu : muons) {
				int this_type = GetLeptonType(mu, Gens);
				gen_muon_type.push_back(this_type);
			}

			path = param.Name + "/TT/";
			FillHist(path + "METv_pt", METv.Pt(), weight, 240, 0., 240.);
            FillHist(path + "METv_xyCorr_pt", METv_xyCorr.Pt(), weight, 240, 0., 240.);   
            FillHist(path + "METv_phi", METv.Phi(), weight, 70, -3.5, 3.5);
            FillHist(path + "METv_xyCorr_phi", METv_xyCorr.Phi(), weight, 70, -3.5, 3.5);
			FillHist(path + "Nj", jets.size(), weight, 8, 0., 8.);
            FillHist(path + "Nb", bjets.size(), weight, 8, 0., 8.);
            FillHist(path + "Nj_dR04", jets_dR04.size(), weight, 8, 0., 8.);
            FillHist(path + "Nb_dR04", bjets_dR04.size(), weight, 8, 0., 8.);
            DrawHists(path, muons.at(0), 0, gen_muon_type.at(0), weight);
            DrawHists(path, muons.at(1), 1, gen_muon_type.at(1), weight);
            if (jets_dR04.size() > 0) DrawHists(path, jets_dR04.at(0), 0, weight);
            if (jets_dR04.size() > 1) DrawHists(path, jets_dR04.at(1), 1, weight);
		}
	}

	if (RunElectron) {
		param.Name = "electron_nocut_central";
		TString path = param.Name + "/";

		//==== event selection ====
		bool isDY = IsDY(electrons, muons_veto);
		bool isTT = IsTT(electrons, muons_veto, jets_dR04, bjets_dR04, METv_xyCorr);
		vector<int> gen_elec_type; // 1 for EWPrompt, 3 for tau

		if (isDY) {
			gen_elec_type.clear();
			for (const auto &el: electrons) {
				int this_type = GetLeptonType(el, Gens);
				gen_elec_type.push_back(this_type);
			}

			path = param.Name + "/DY/";
			FillHist(path + "METv_pt", METv.Pt(), weight, 240, 0., 240.);
            FillHist(path + "METv_xyCorr_pt", METv_xyCorr.Pt(), weight, 240, 0., 240.);
            FillHist(path + "METv_phi", METv.Phi(), weight, 70, -3.5, 3.5);
            FillHist(path + "METv_xyCorr_phi", METv_xyCorr.Phi(), weight, 70, -3.5, 3.5);
			FillHist(path + "Nj", jets.size(), weight, 8, 0., 8.);
            FillHist(path + "Nb", bjets.size(), weight, 8, 0., 8.);
            FillHist(path + "Nj_dR04", jets_dR04.size(), weight, 8, 0., 8.);
            FillHist(path + "Nb_dR04", bjets_dR04.size(), weight, 8, 0., 8.);
            DrawHists(path, electrons.at(0), 0, gen_elec_type.at(0), weight);
            DrawHists(path, electrons.at(1), 1, gen_elec_type.at(1), weight);
            if (jets_dR04.size() > 0) DrawHists(path, jets_dR04.at(0), 0, weight);
            if (jets_dR04.size() > 1) DrawHists(path, jets_dR04.at(1), 1, weight);
		}

		if (isTT) {
			gen_elec_type.clear();
            for (const auto &el: electrons) {
                int this_type = GetLeptonType(el, Gens);
                gen_elec_type.push_back(this_type);
            }
			
			path = param.Name + "/TT/";
			FillHist(path + "METv_pt", METv.Pt(), weight, 240, 0., 240.);
            FillHist(path + "METv_xyCorr_pt", METv_xyCorr.Pt(), weight, 240, 0., 240.);
            FillHist(path + "METv_phi", METv.Phi(), weight, 70, -3.5, 3.5);
            FillHist(path + "METv_xyCorr_phi", METv_xyCorr.Phi(), weight, 70, -3.5, 3.5);
			FillHist(path + "Nj", jets.size(), weight, 8, 0., 8.);
            FillHist(path + "Nb", bjets.size(), weight, 8, 0., 8.);
            FillHist(path + "Nj_dR04", jets_dR04.size(), weight, 8, 0., 8.);
            FillHist(path + "Nb_dR04", bjets_dR04.size(), weight, 8, 0., 8.);
            DrawHists(path, electrons.at(0), 0, gen_elec_type.at(0), weight);
            DrawHists(path, electrons.at(1), 1, gen_elec_type.at(1), weight);
            if (jets_dR04.size() > 0) DrawHists(path, jets_dR04.at(0), 0, weight);
            if (jets_dR04.size() > 1) DrawHists(path, jets_dR04.at(1), 1, weight);
		}
	}
}

id_var_prompt::id_var_prompt(){

}

id_var_prompt::~id_var_prompt(){

}

//==== Member functions ====
bool id_var_prompt::IsDY(const vector<Muon> &muons, const vector<Electron> &electrons_veto) {
	TString path = "muon_nocut_central/DY/cutflow";
	FillHist(path, 0., 1., 5, 0., 5.);
	
	if (! (muons.size() == 2)) return false;
	if (! (muons.at(0).Charge()*muons.at(1).Charge() < 0)) return false;
	if (! (electrons_veto.size() == 0)) return false;
	FillHist(path, 1., 1., 5, 0., 5.);

	Particle ZCand = muons.at(0) + muons.at(1);
	if (! (IsOnZ(ZCand.M(), 15))) return false;
	FillHist(path, 2., 1., 5, 0., 5.);

	return true;
}

bool id_var_prompt::IsDY(const vector<Electron> &electrons, const vector<Muon> &muons_veto) {
	TString path = "electron_nocut_central/DY/cutflow";
	FillHist(path, 0., 1., 5, 0., 5.);

	if (! (electrons.size() == 2)) return false;
	if (! (electrons.at(0).Charge()*electrons.at(1).Charge() < 0)) return false;
	if (! (muons_veto.size() == 0)) return false;
	FillHist(path, 1., 1., 5, 0., 5.);

	Particle ZCand = electrons.at(0) + electrons.at(1);
	if (! (IsOnZ(ZCand.M(), 15))) return false;
	FillHist(path, 2., 1., 5, 0., 5.);

	return true;
}

bool id_var_prompt::IsTT(const vector<Muon> &muons, const vector<Electron> &electrons_veto, const vector<Jet> &jets, const vector<Jet> &bjets, const Particle &METv) {
	TString path = "muon_nocut_central/TT/cutflow";
	FillHist(path, 0., 1., 5, 0., 5.);

	if (! (muons.size() == 2)) return false;
	if (! (muons.at(0).Charge()*muons.at(1).Charge() < 0)) return false;
	if (! (muons.at(0).DeltaR(muons.at(1)) > 0.4)) return false;
	if (! (electrons_veto.size() == 0)) return false;
	FillHist(path, 1., 1., 5, 0., 5.);

	if (! (jets.size() >= 2)) return false;
	FillHist(path, 2., 1., 5, 0., 5.);
	if (! (bjets.size() >= 1)) return false;
	FillHist(path, 3., 1., 5, 0., 5.);
	if (! (METv.Pt() > 40)) return false;
	FillHist(path, 4., 1., 5, 0., 5.);

	return true;
}

bool id_var_prompt::IsTT(const vector<Electron> &electrons, const vector<Muon> &muons_veto, const vector<Jet> &jets, const vector<Jet> &bjets, const Particle &METv) {
	TString path = "electron_nocut_central/TT/cutflow";
	FillHist(path, 0., 1., 5, 0., 5.);

	if (! (electrons.size() == 2)) return false;
	if (! (electrons.at(0).Charge()*electrons.at(1).Charge() < 0)) return false;
	if (! (electrons.at(0).DeltaR(electrons.at(1)) > 0.4)) return false;
	if (! (muons_veto.size() == 0)) return false;
	FillHist(path, 1., 1., 5, 0., 5.);

	if (! (jets.size() >= 2)) return false;
    FillHist(path, 2., 1., 5, 0., 5.);
    if (! (bjets.size() >= 1)) return false;
    FillHist(path, 3., 1., 5, 0., 5.);
    if (! (METv.Pt() > 40)) return false;
    FillHist(path, 4., 1., 5, 0., 5.);

    return true;
}

void id_var_prompt::DrawHists(TString path, const Muon &mu, unsigned int order, int type, const double weight) {
	TString this_path = path + "mu" + TString::Itoa(order+1, 10) + "/";
	this_path += "type" + TString::Itoa(type, 10) + "/";
	FillHist(this_path + "pt", mu.Pt(), weight, 240, 0, 240);
    FillHist(this_path + "eta", mu.Eta(), weight, 50, -2.5, 2.5);
    FillHist(this_path + "phi", mu.Phi(), weight, 70, -3.5, 3.5);

    this_path += "IDvar/";
    FillHist(this_path + "RelIso", mu.RelIso(), weight, 100, 0, 0.6);
    FillHist(this_path + "dXY" , mu.dXY(), weight, 50, -0.04, 0.04);
    FillHist(this_path + "dXYerr", mu.dXYerr(), weight, 50, 0., 0.02);
    if (mu.dXYerr() != 0) {
        FillHist(this_path + "SIP2D", fabs(mu.dXY())/mu.dXYerr(), weight, 100, 0., 8.);
    }
    FillHist(this_path + "dZ", mu.dZ(), weight, 100, -0.06, 0.06);
    FillHist(this_path + "Chi2", mu.Chi2(), weight, 100, 0., 5.);
}

void id_var_prompt::DrawHists(TString path, const Electron &e, unsigned int order, int type, const double weight) {
	TString this_path = path + "e" + TString::Itoa(order+1, 10) + "/";
	this_path += "type" + TString::Itoa(type, 10) + "/";
	FillHist(this_path + "pt", e.Pt(), weight, 240, 0, 240);
	FillHist(this_path + "eta", e.Eta(), weight, 50, -2.5, 2.5);
	FillHist(this_path + "phi", e.Phi(), weight, 70, -3.5, 3.5);

	this_path += "IDvar/";
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
void id_var_prompt::DrawHists(TString path, const Jet &j, unsigned int order, const double weight) {

    TString this_path = path + "j" + TString::Itoa(order+1, 10) + "/";
    FillHist(this_path + "pt", j.Pt(), weight, 240, 0, 240);
    FillHist(this_path + "eta", j.Eta(), weight, 50, -2.5, 2.5);
    FillHist(this_path + "phi", j.Phi(), weight, 70, -3.5, 3.5);
}	




