#include "CRstudy.h"

CRstudy::CRstudy() {}
CRstudy::~CRstudy() {}

void CRstudy::initializeAnalyzer(){
	// Only for year 2017 now
	if (DataYear == 2017) {
		trigs_dimu.emplace_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v");
        trigs_dimu.emplace_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v");
        trigs_dimu.emplace_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v");
    
        trigs_emu.emplace_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v");
        trigs_emu.emplace_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");
    }
	else
		cerr << "DataYear " << DataYear << " is not set yet" << endl;
	
	// IDs
    HcToWA_MuID = {"HcToWATight", "HcToWALoose"};
    HcToWA_EleID = {"HcToWATight", "HcToWALoose"};

    // B-tagging
    vector<JetTagging::Parameters> jtps;
    jtps.emplace_back(
        JetTagging::Parameters(JetTagging::DeepJet, JetTagging::Medium, 
            JetTagging::incl, JetTagging::mujets));
    mcCorr->SetJetTaggingParameters(jtps);

	// initiate regions and cutflows
	regions = {"DY_DiMu", "TT_DiMu", "TT_EMu"};
	InitiateCutflow("DY_DiMu", this->getCuts("DY_DiMu"));
	InitiateCutflow("TT_DiMu", this->getCuts("TT_DiMu"));
	InitiateCutflow("TT_EMu", this->getCuts("TT_EMu"));
}

vector<TString> CRstudy::getCuts(TString region) {
	if (region == "DY_DiMu") 
		return {"nocut", "metfilter", "OSDiMu", "trigger", "passSafePtCut", "OnshellZ", "NoBjet"};
	else if (region == "TT_DiMu")
		return {"nocut", "metfilter", "OSDiMu", "trigger", "passSafePtCut", "OffshellZ", "MDiMu_ge10", "dR_ge04", "Nj_ge2", "Nb_ge1"};
	else if (region == "TT_EMu")
		return {"nocut", "metfilter", "OSEMu", "trigger", "passSafePtCut", "dR_ge04", "Nj_ge2", "Nb_ge1"};
	else {
		cerr << "[CRstudy::getCuts] Wrong Region " << region << endl;
		exit(EXIT_FAILURE);
	}
}
void CRstudy::executeEvent(){
	//FillCutflow
	for (const auto& region: regions)
		FillCutflow(region, "nocut");

	if (!PassMETFilter()) return;
	for (const auto& region: regions)
		FillCutflow(region, "metfilter");

	Event ev = GetEvent();
	vector<Gen> gens = GetGens();
	vector<Muon> muons = GetAllMuons();
	vector<Electron> electrons = GetAllElectrons();
	vector<Jet> jets = GetAllJets();
	Particle METv = ev.GetMETVector();

	// sort at the firtst time, no willing to be confused
	sort(muons.begin(), muons.end(), PtComparing);
	sort(electrons.begin(), electrons.end(), PtComparing);
	sort(jets.begin(), jets.end(), PtComparing);

	// select objects
	vector<Muon> muons_tight = SelectMuons(muons, "HcToWATight", 10., 2.4);
	vector<Muon> muons_loose = SelectMuons(muons, "HcToWALoose", 10., 2.4);
	vector<Electron> electrons_tight = SelectElectrons(electrons, "HcToWATight", 10., 2.5);
	vector<Electron> electrons_loose = SelectElectrons(electrons, "HcToWALoose", 10., 2.5);
	vector<Jet> jets_tight = SelectJets(jets, "tight", 25., 2.4);
	vector<Jet> jets_cleaned = JetsVetoLeptonInside(jets_tight, electrons_loose, muons_loose, 0.4);
	//vector<Jet> jets_puveto = SelectPUvetoJets(jets_cleaned, "medium");
	vector<Jet> bjets_cleaned;
	for (const auto& jet: jets_cleaned) {
		const double this_discr = jet.GetTaggerResult(JetTagging::DeepCSV);
		if (this_discr > mcCorr->GetJetTaggingCutValue(JetTagging::DeepJet, JetTagging::Medium))
			bjets_cleaned.emplace_back(jet);
	}

	// Note: Cutflow will automatically generated inside the SignalSelector
	const TString channel = ControlRegionSelector(
			ev, muons_tight, electrons_tight,
			muons_loose, electrons_loose,
			jets_cleaned, bjets_cleaned);
	if (channel == "")
		return;

	// set weight
	double weight = 1.;
	if (!IsDATA) {
		const double w_prefire = GetPrefireWeight(0);
		const double w_gen = ev.MCweight() * weight_norm_1invpb;
		const double w_lumi = ev.GetTriggerLumi("Full");
		const double w_pileup = GetPileUpWeight(nPileUp, 0);
		//cout << "w_prefire: " <<  w_prefire << endl;
		//cout << "w_gen: " << w_gen << endl;
		//cout << "w_lumi: " << w_lumi << endl;
		//cout << "w_pileup: " << w_pileup << endl;
		weight *= w_prefire*w_gen*w_lumi*w_pileup;

		// ID scale factors
		double w_idsf = 1.;
		const TString ID = "HcToWATight";
		double w_trigsf = 1.;
		if (channel.Contains("DiMu")) {
			const Muon& mu1 = muons_tight.at(0);
			const Muon& mu2 = muons_tight.at(1);
			const double mu1_idsf = mcCorr->MuonID_SF(ID, mu1.Eta(), mu1.MiniAODPt(), 0);
			const double mu2_idsf = mcCorr->MuonID_SF(ID, mu2.Eta(), mu2.MiniAODPt(), 0);
			w_idsf *= mu1_idsf*mu2_idsf;
			w_trigsf = mcCorr->GetTriggerSF(electrons_tight, muons_tight, "DiMuIso_HNTopID", "");
		}
		if (channel.Contains("EMu")) {
			const Muon& mu = muons_tight.at(0);
			const Electron& ele = electrons_tight.at(0);
			const double mu_idsf = mcCorr->MuonID_SF(ID, mu.Eta(), mu.MiniAODPt(), 0);
			const double ele_idsf = mcCorr->ElectronID_SF(ID, ele.scEta(), ele.Pt(), 0);
			w_idsf *= mu_idsf*ele_idsf;
			w_trigsf = mcCorr->GetTriggerSF(electrons_tight, muons_tight, "EMuIso_HNTopID", "");
		}
		//cout << "w_idsf: " << w_idsf << endl;
		//cout << "w_trigsf: " << w_trigsf << endl;
		weight *= w_idsf*w_trigsf;

		// b-tagging SF
		JetTagging::Parameters jtp_DeepJet_Medium 
			= JetTagging::Parameters(JetTagging::DeepJet, JetTagging::Medium, JetTagging::incl, JetTagging::mujets);
								
		const double w_btag = mcCorr->GetBTaggingReweight_1a(jets_cleaned, jtp_DeepJet_Medium);
		//cout << "w_btag: " << w_btag << endl;
		weight *= w_btag;
		
		// pu-veto SF
		//const double w_puveto = mcCorr->GetPUVetoSF(jets_cleaned, "medium");
		//cout << "w_puveto: " << w_puveto << endl;
		//weight *= w_puveto;
	}
	
	
	// Fill Hist
	FillObjects(channel + "/muons_tight", muons_tight, weight);
	FillObjects(channel + "/electrons_tight", electrons_tight, weight);
	FillObjects(channel + "/jets_cleaned", jets_cleaned, weight);
	//FillObjects(channel + "/jets_puveto", jets_puveto, weight);
	FillObjects(channel + "/bjets_cleaned", bjets_cleaned, weight);
	FillObject(channel + "/METv", METv, weight);
	if (channel == "DY_DiMu" || channel == "TT_DiMu") {
		const Particle ZCand = muons_tight.at(0) + muons_tight.at(1);
		FillObject(channel + "/ZCand", ZCand, weight);
	}
}

//==== Lets Select Events
TString CRstudy::ControlRegionSelector(
		Event &ev,
		vector<Muon> &muons_tight, vector<Electron> &electrons_tight,
		vector<Muon> &muons_loose, vector<Electron> &electrons_loose,
		vector<Jet> &jets, vector<Jet> &bjets) {
	
	TString region = "";
	// leptons first
	if (muons_tight.size() == 2	&& muons_loose.size() == 2 && electrons_tight.size() == 0 && electrons_loose.size() == 0)
		region = "DiMu";
	else if (muons_tight.size() == 1 && muons_loose.size() == 1 && electrons_tight.size() == 1 && electrons_loose.size() == 1)
		region = "EMu";
	else
		return "";

	if (region == "DiMu") {
		if (muons_tight.at(0).Charge() + muons_tight.at(1).Charge() != 0)
			return "";
		FillCutflow("DY_DiMu", "OSDiMu");
		FillCutflow("TT_DiMu", "OSDiMu");
	}
	if (region == "EMu") {
		if (muons_tight.at(0).Charge() + electrons_tight.at(0).Charge() != 0)
			return "";
		FillCutflow("TT_EMu", "OSEMu");
	}

	
	if (region == "DiMu") {
		if (! ev.PassTrigger(trigs_dimu)) return "";
		FillCutflow("DY_DiMu", "trigger");
		FillCutflow("TT_DiMu", "trigger");

		// safe pt cut, 20/10
		if (! (muons_tight.at(0).Pt() > 20.))
			return "";
		if (! (muons_tight.at(1).Pt() > 10.))
			return "";
		FillCutflow("DY_DiMu", "passSafePtCut");
		FillCutflow("TT_DiMu", "passSafePtCut");
		
		// Let's diverge dimuon regions
		Particle ZCand = muons_tight.at(0) + muons_tight.at(1);
		if (fabs(ZCand.M() - 91.2) < 15.) {
			region = "DY_DiMu";
			FillCutflow("DY_DiMu", "OnshellZ");
		}
		else {
			region = "TT_DiMu";
			FillCutflow("TT_DiMu", "OffshellZ");

			if (ZCand.M() < 10.)
				return "";
			FillCutflow("TT_DiMu", "MDiMu_ge10");
		}
	}
	
	// Further cuts
	if (region == "DY_DiMu") {
		if (bjets.size() > 0)
			return "";
		FillCutflow("DY_DiMu", "NoBjet");
	}
	if (region == "TT_DiMu") {
		const double dRll = muons_tight.at(0).DeltaR(muons_tight.at(1));
		if (dRll < 0.4)
			return "";
		FillCutflow("TT_DiMu", "dR_ge04");
		
		if (jets.size() < 2)
			return "";
		FillCutflow("TT_DiMu", "Nj_ge2");

		if (bjets.size() == 0)
			return "";
		FillCutflow("TT_DiMu", "Nb_ge1");
	}

	// emu
	if (region == "EMu") {
		region = "TT_EMu";
		if (! ev.PassTrigger(trigs_emu)) return "";
		FillCutflow(region, "trigger");

		bool passSafeCut = false;
		if (muons_tight.at(0).Pt() > 10. && electrons_tight.at(0).Pt() > 25.)
			passSafeCut = true;
		if (muons_tight.at(0).Pt() > 25. && electrons_tight.at(0).Pt() > 15.)
			passSafeCut = true;
		if (!passSafeCut)
			return "";
		FillCutflow(region, "passSafePtCut");

		const double dRll = muons_tight.at(0).DeltaR(electrons_tight.at(0));
		if (dRll < 0.4)
			return "";
		FillCutflow(region, "dR_ge04");

		if (jets.size() < 2)
			return "";
		FillCutflow(region, "Nj_ge2");

		if (bjets.size() == 0)
			return "";
		FillCutflow(region, "Nb_ge1");
	}
	return region;
}
