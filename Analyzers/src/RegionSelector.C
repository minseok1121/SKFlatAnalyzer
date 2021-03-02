#include "RegionSelector.h"

RegionSelector::RegionSelector() {}
RegionSelector::~RegionSelector() {}

void RegionSelector::initializeAnalyzer(){
 	// flags
    RunDeepCSV = HasFlag("RunDeepCSV");
    EMuTrigOnly = HasFlag("EMuTrigOnly");
    SkimBaseline = HasFlag("SkimBaseline");
    cout << "[RegionSelector::RegionSelector] RunDeepCSV = " << RunDeepCSV << endl;
    cout << "[RegionSelector::RegionSelector] EMuTrigOnly = " << EMuTrigOnly << endl;
    cout << "[RegionSelector::RegionSelector] SkimBaseline = " << SkimBaseline << endl;

    // TODO: prepare output TTree for SkimBaseline
    // triggers
    // only 2017 triggers for a time being
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

    // B-tagging
    vector<JetTagging::Parameters> jtps;
    if (RunDeepCSV)
        jtps.emplace_back(JetTagging::Parameters(JetTagging::DeepCSV, JetTagging::Medium, JetTagging::incl, JetTagging::mujets));
    else
         jtps.emplace_back(JetTagging::Parameters(JetTagging::DeepJet, JetTagging::Medium, JetTagging::incl, JetTagging::mujets));
    mcCorr->SetJetTaggingParameters(jtps);
	
	// initiate regions and cutflows
    regions = {"SR_3mu", "SR_1e2mu", "WZ_3mu", "WZ_1e2mu", "DY_OSdimu", "TT_OSdimu", "TT_OSemu"};
    for (const auto& region: regions)
        InitiateCutflow(region, this->getCuts(region));
}

void RegionSelector::executeEvent(){
	// Fill All Cutflows
	for (const auto& region: regions)
		FillCutflow(region, "noCut");

	if (!PassMETFilter()) return;
	for (const auto& region: regions)
		FillCutflow(region, "METFilter");

	// Object Definitions
	Event ev = GetEvent();
	vector<Gen> gens = GetGens();
	vector<Muon> muons = GetAllMuons();
	vector<Electron> electrons = GetAllElectrons();
	vector<Jet> jets = GetAllJets();
	Particle METv = ev.GetMETVector(); // phi correction applid in SKFlat

	// sort at the first time, don't want to be confused
	sort(muons.begin(), muons.end(), PtComparing);
	sort(electrons.begin(), electrons.end(), PtComparing);
	sort(jets.begin(), jets.end(), PtComparing);

	vector<Muon> muons_tight = SelectMuons(muons, "HcToWATight", 10., 2.4);
	vector<Muon> muons_loose = SelectMuons(muons, "HcToWALoose", 10., 2.4);
	vector<Electron> electrons_tight = SelectElectrons(electrons, "HcToWATgiht", 10., 2.5);
	vector<Electron> electrons_loose = SelectElectrons(electrons, "HcToWALoose", 10., 2.5);
	vector<Jet> jets_tight = SelectJets(jets, "tight", 20., 2.4);
	vector<Jet> jets_lepVeto = JetsVetoLeptonInside(jets_tight, electrons_loose, muons_loose, 0.4);
	
	vector<Jet> bjets_lepVeto, jets_nonBtag;	
	if (RunDeepCSV) {
		for (const auto& jet: jets_lepVeto) {
			const double this_discr = jet.GetTaggerResult(JetTagging::DeepCSV);
			if (this_discr > mcCorr->GetJetTaggingCutValue(JetTagging::DeepCSV, JetTagging::Medium))
				bjets_lepVeto.emplace_back(jet);
			else
				jets_nonBtag.emplace_back(jet);
		}
	}
	else {
		for (const auto& jet: jets_lepVeto) {
			const double this_discr = jet.GetTaggerResult(JetTagging::DeepJet);
			if (this_discr > mcCorr->GetJetTaggingCutValue(JetTagging::DeepJet, JetTagging::Medium))
				bjets_lepVeto.emplace_back(jet);
			else
				jets_nonBtag.emplace_back(jet);
		}
	}

	const TString channel = Selector(
			ev, muons_tight, electrons_tight,
			muons_loose, electrons_loose,
			jets_lepVeto, bjets_lepVeto, METv);
	if (channel == "")
		return;

	// initialize weight
	weight = 1.;
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
        if (channel.Contains("dimu")) {
            const Muon& mu1 = muons_tight.at(0);
            const Muon& mu2 = muons_tight.at(1);
            const double mu1_idsf = mcCorr->MuonID_SF(ID, mu1.Eta(), mu1.MiniAODPt(), 0);
            const double mu2_idsf = mcCorr->MuonID_SF(ID, mu2.Eta(), mu2.MiniAODPt(), 0);
            w_idsf *= mu1_idsf*mu2_idsf;
            w_trigsf = mcCorr->GetTriggerSF(electrons_tight, muons_tight, "DiMuIso_HNTopID", "");
        }
        if (channel.Contains("emu")) {
            const Muon& mu = muons_tight.at(0);
            const Electron& ele = electrons_tight.at(0);
            const double mu_idsf = mcCorr->MuonID_SF(ID, mu.Eta(), mu.MiniAODPt(), 0);
            const double ele_idsf = mcCorr->ElectronID_SF(ID, ele.scEta(), ele.Pt(), 0);
            //cout << "mu_idsf: " << mu_idsf << endl;
            //cout << "ele_idsf: " << ele_idsf << endl;
            w_idsf *= mu_idsf*ele_idsf;
            w_trigsf = mcCorr->GetTriggerSF(electrons_tight, muons_tight, "EMuIso_HNTopID", "");
        }
        if (channel.Contains("3mu")) {
            const Muon& mu1 = muons_tight.at(0);
            const Muon& mu2 = muons_tight.at(1);
            const Muon& mu3 = muons_tight.at(2);
            const double mu1_idsf = mcCorr->MuonID_SF(ID, mu1.Eta(), mu1.MiniAODPt(), 0);
            const double mu2_idsf = mcCorr->MuonID_SF(ID, mu2.Eta(), mu2.MiniAODPt(), 0);
            const double mu3_idsf = mcCorr->MuonID_SF(ID, mu3.Eta(), mu3.MiniAODPt(), 0);
            w_idsf *= mu1_idsf*mu2_idsf*mu3_idsf;
			w_trigsf = mcCorr->GetTriggerSF(electrons_tight, muons_tight, "DiMuIso_HNTopID", "");
        }
        if (channel.Contains("1e2mu")) {
            const Muon& mu1 = muons_tight.at(0);
            const Muon& mu2 = muons_tight.at(1);
            const Electron& ele = electrons_tight.at(0);
            const double mu1_idsf = mcCorr->MuonID_SF(ID, mu1.Eta(), mu1.MiniAODPt(), 0);
            const double mu2_idsf = mcCorr->MuonID_SF(ID, mu2.Eta(), mu2.MiniAODPt(), 0);
            const double ele_idsf = mcCorr->ElectronID_SF(ID, ele.scEta(), ele.Pt(), 0);
            w_idsf *= mu1_idsf*mu2_idsf*ele_idsf;
			if (EMuTrigOnly)
				w_trigsf = mcCorr->GetTriggerSF(electrons_tight, muons_tight, "EMuIso_HNTopID", "");
			else {
				if (ev.PassTrigger(trigs_emu))
					w_trigsf = mcCorr->GetTriggerSF(electrons_tight, muons_tight, "EMuIso_HNTopID", "");
				else if (ev.PassTrigger(trigs_dblmu))
					w_trigsf = mcCorr->GetTriggerSF(electrons_tight, muons_tight, "DiMuIso_HNTopID", "");
			}
        }
		weight *= w_idsf*w_trigsf;

		// b-tagging SF
		double w_btag = 1.;
		if (RunDeepCSV) {
			JetTagging::Parameters jtp_DeepCSV_Medium
			    = JetTagging::Parameters(JetTagging::DeepCSV, JetTagging::Medium, JetTagging::incl, JetTagging::mujets);
			w_btag = mcCorr->GetBTaggingReweight_1a(jets_lepVeto, jtp_DeepCSV_Medium);
		}
		else {
			JetTagging::Parameters jtp_DeepJet_Medium
			    = JetTagging::Parameters(JetTagging::DeepJet, JetTagging::Medium, JetTagging::incl, JetTagging::mujets);
			w_btag = mcCorr->GetBTaggingReweight_1a(jets_lepVeto, jtp_DeepJet_Medium);
		}
		//cout << "w_btag: " << w_btag << endl;
		weight *= w_btag;
	}

	// Fill Hists
	FillObjects(channel + "/muons_tight", muons_tight, weight);
	FillObjects(channel + "/electrons_tight", electrons_tight, weight);
	FillObjects(channel + "/jets_lepVeto", jets_lepVeto, weight);
	FillObjects(channel + "/bjets_lepVeto", bjets_lepVeto, weight);
	FillObject(channel + "/METv", METv, weight);
	if (channel == "DY_dimu" || channel == "TT_dimu") {
		const Particle ZCand = muons_tight.at(0) + muons_tight.at(1);
		FillObject(channel + "/ZCand", ZCand, weight);
	}
}

















