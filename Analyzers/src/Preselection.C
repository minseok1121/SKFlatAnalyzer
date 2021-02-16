#include "Preselection.h"

Preselection::Preselection() {}
Preselection::~Preselection() {
	if (SkimBaseline) {
		outfile->cd();
		tree->Write();
	}
}

void Preselection::initializeAnalyzer(){

	// flags
	RunDeepCSV = HasFlag("RunDeepCSV");
	RunPUVeto = HasFlag("RunPUVeto");
	SkimBaseline = HasFlag("SkimBaseline");
	SkipPUweight = HasFlag("SkipPUweight");
	cout << "[Preselection::initializeAnalyzer] RunDeepCSV = " << RunDeepCSV << endl;
	cout << "[Preselection::initializeAnalyzer] SkimBaseline = " << SkimBaseline << endl;

	// prepare root file for RDataFrame
	// TTree base
	if (SkimBaseline) {
		tree = new TTree("Events", "Events");
		tree->Branch("ptMu1", &ptMu1);
		tree->Branch("ptMu2", &ptMu2);
		tree->Branch("ptMu3", &ptMu3);
		tree->Branch("etaMu1", &etaMu1);
		tree->Branch("etaMu2", &etaMu2);
		tree->Branch("etaMu3", &etaMu3);
		tree->Branch("phiMu1", &phiMu1);
		tree->Branch("phiMu2", &phiMu2);
		tree->Branch("phiMu3", &phiMu3);
		tree->Branch("ptJ1", &ptJ1);
		tree->Branch("ptJ2", &ptJ2);
		tree->Branch("etaJ1", &etaJ1);
		tree->Branch("etaJ2", &etaJ2);
		tree->Branch("phiJ1", &phiJ1);
		tree->Branch("phiJ2", &phiJ2);
		tree->Branch("dRl1l2", &dRl1l2);
		tree->Branch("dRl1l3", &dRl1l3);
		tree->Branch("dRl2l3", &dRl2l3);
		tree->Branch("dRj1l1", &dRj1l1);
		tree->Branch("dRj1l2", &dRj1l2);
		tree->Branch("dRj1l3", &dRj1l3);
		tree->Branch("dRj2l1", &dRj2l3);
		tree->Branch("dRj2l2", &dRj2l2);
		tree->Branch("dRj2l3", &dRj2l3);
		tree->Branch("dRj1j2", &dRj1j2);
		tree->Branch("HT", &HT);
		tree->Branch("ST", &ST);
		tree->Branch("HToverST", &HToverST);
		tree->Branch("MET", &MET);
		//tree->Branch("mMuMu", &mMuMu);
		//tree->Branch("mMuMu2", &mMuMu2);
	}
	// triggers
	// only for 2017 currently
	if (DataYear == 2017) {
		trigs_dimu.emplace_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v");
        trigs_dimu.emplace_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v");
        trigs_dimu.emplace_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v");
    
        trigs_emu.emplace_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v");
        trigs_emu.emplace_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");
    }   
    else {
        cerr << "DataYear " << DataYear << " is not set yet" << endl;
		exit(EXIT_FAILURE);
	}

	// IDs
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
	regions = {"Pre_3mu", "Pre_1e2mu", "DY_dimu", "TT_dimu", "TT_emu"};
	for (const auto& region : regions)
		InitiateCutflow(region, this->getCuts(region));

}

vector<TString> Preselection::getCuts(TString region) {
	if (region == "Pre_3mu")
		return {"nocut", "metfilter", "3mu", "trigger", "passSafePtCut", "exist_osmu", "mass_ge12", "Nj_ge2", "Nb_ge1"};
	else if (region == "Pre_1e2mu")
		return {"nocut", "metfilter", "1e2mu", "trigger", "passSafePtCut", "exist_osmu", "mass_ge12", "Nj_ge2", "Nb_ge1"};
	else if (region == "DY_dimu")
		return {"nocut", "metfilter", "dimu", "trigger", "passSafePtCut", "os_dimu", "OnshellZ", "NoBjet"};
	else if (region == "TT_dimu")
		return {"nocut", "metfilter", "dimu", "trigger", "passSafePtCut", "os_dimu", "OffshellZ", "mass_ge12", "dR_ge04", "Nj_ge2", "Nb_ge1"};
	else if (region == "TT_emu")
		return {"nocut", "metfilter", "emu", "trigger", "passSafePtCut", "os_emu", "dR_ge04", "Nj_ge2", "Nb_ge1"};
	else {
		cerr << "[Preselection::getCuts] wrong region " << region << endl;
		exit(EXIT_FAILURE);
	}
}

void Preselection::executeEvent(){
	//FillCutflow
	for (const auto& region: regions)
		FillCutflow(region, "nocut");

	// metfilter
	if (!PassMETFilter()) return;
	for (const auto& region: regions)
		FillCutflow(region, "metfilter");

	Event ev = GetEvent();
	vector<Gen> gens = GetGens();
	vector<Muon> muons = GetAllMuons();
	vector<Electron> electrons = GetAllElectrons();
	vector<Jet> jets = GetAllJets();
	Particle METv = ev.GetMETVector();

	// sort the objects at the first time, don't want to be confused
	sort(muons.begin(), muons.end(), PtComparing);
	sort(electrons.begin(), electrons.end(), PtComparing);
	sort(jets.begin(), jets.end(), PtComparing);

	// select objects
	vector<Muon> muons_tight = SelectMuons(muons, "HcToWATight", 10., 2.4);
	vector<Muon> muons_loose = SelectMuons(muons, "HcToWALoose", 10., 2.4);
	vector<Electron> electrons_tight = SelectElectrons(electrons, "HcToWATight", 10., 2.5);
	vector<Electron> electrons_loose = SelectElectrons(electrons, "HcToWALoose", 10., 2.5);
	vector<Jet> jets_tight = SelectJets(jets, "tight", 20., 2.4);
	vector<Jet> jets_cleaned = JetsVetoLeptonInside(jets_tight, electrons_loose, muons_loose, 0.4);
	if (RunPUVeto)
		jets_cleaned = SelectPUvetoJets(jets_cleaned, "tight");
	vector<Jet> bjets_cleaned;
	if (RunDeepCSV)
		for (const auto& jet: jets_cleaned) {
			const double this_discr = jet.GetTaggerResult(JetTagging::DeepCSV);
			if (this_discr > mcCorr->GetJetTaggingCutValue(JetTagging::DeepCSV, JetTagging::Medium))
				bjets_cleaned.emplace_back(jet);
		}
	else
		for (const auto& jet: jets_cleaned) {
			const double this_discr = jet.GetTaggerResult(JetTagging::DeepJet);
			if (this_discr > mcCorr->GetJetTaggingCutValue(JetTagging::DeepJet, JetTagging::Medium))
				bjets_cleaned.emplace_back(jet);
		}

	// Cutflow will be automatically generated inside the Selector
	const TString channel = RegionSelector(
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
		double w_pileup = 1.;
		if (!SkipPUweight)
			w_pileup = GetPileUpWeight(nPileUp, 0); 
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
		}
		if (channel.Contains("1e2mu")) {
			const Muon& mu1 = muons_tight.at(0);
			const Muon& mu2 = muons_tight.at(1);
			const Electron& ele = electrons_tight.at(0);
			const double mu1_idsf = mcCorr->MuonID_SF(ID, mu1.Eta(), mu1.MiniAODPt(), 0);
			const double mu2_idsf = mcCorr->MuonID_SF(ID, mu2.Eta(), mu2.MiniAODPt(), 0);
			const double ele_idsf = mcCorr->ElectronID_SF(ID, ele.scEta(), ele.Pt(), 0);
			w_idsf *= mu1_idsf*mu2_idsf*ele_idsf;
		}
        //cout << "w_idsf: " << w_idsf << endl;
        //cout << "w_trigsf: " << w_trigsf << endl;
        weight *= w_idsf*w_trigsf;

		// b-tagging SF
		double w_btag = 1.;
		if (RunDeepCSV) {
        	JetTagging::Parameters jtp_DeepCSV_Medium
            	= JetTagging::Parameters(JetTagging::DeepCSV, JetTagging::Medium, JetTagging::incl, JetTagging::mujets);

        	w_btag = mcCorr->GetBTaggingReweight_1a(jets_cleaned, jtp_DeepCSV_Medium);
		}
		else {
			JetTagging::Parameters jtp_DeepJet_Medium
				= JetTagging::Parameters(JetTagging::DeepJet, JetTagging::Medium, JetTagging::incl, JetTagging::mujets);
			w_btag = mcCorr->GetBTaggingReweight_1a(jets_cleaned, jtp_DeepJet_Medium);
		}
        //cout << "w_btag: " << w_btag << endl;
        weight *= w_btag;

		// puveto
		if (RunPUVeto) {
			const double w_puveto = mcCorr->GetPUVetoSF(jets_cleaned, "tight");
			//cout << "w_puveto: " << w_puveto << endl;
			weight *= w_puveto;
		}
	}
	// Fill Hist
	FillObjects(channel + "/muons_tight", muons_tight, weight);
    FillObjects(channel + "/electrons_tight", electrons_tight, weight);
    FillObjects(channel + "/jets_cleaned", jets_cleaned, weight);
    //FillObjects(channel + "/jets_puveto", jets_puveto, weight);
    FillObjects(channel + "/bjets_cleaned", bjets_cleaned, weight);
    FillObject(channel + "/METv", METv, weight);
    if (channel == "DY_dimu" || channel == "TT_dimu") {
        const Particle ZCand = muons_tight.at(0) + muons_tight.at(1);
        FillObject(channel + "/ZCand", ZCand, weight);
    }
	// Fill TTree
	if (SkimBaseline) {
		if (channel == "Pre_3mu") {
			const Muon& mu1 = muons_tight.at(0);
			const Muon& mu2 = muons_tight.at(1);
			const Muon& mu3 = muons_tight.at(2);
			const Jet& j1 = jets_cleaned.at(0);
			const Jet& j2 = jets_cleaned.at(1);
			ptMu1 = mu1.Pt(); ptMu2 = mu2.Pt(); ptMu3 = mu3.Pt();
			etaMu1 = mu1.Eta(); etaMu2 = mu2.Eta(); etaMu3 = mu3.Eta();
			phiMu1 = mu1.Phi(); phiMu2 = mu2.Phi(); phiMu3 = mu3.Phi();
			ptJ1 = j1.Pt(); ptJ2 = j2.Pt();
			etaJ1 = j1.Eta(); etaJ2 = j2.Eta();
			phiJ1 = j1.Phi(); phiJ2 = j2.Phi();
			dRl1l2 = mu1.DeltaR(mu2); dRl1l3 = mu1.DeltaR(mu3); dRl2l3 = mu2.DeltaR(mu3);
			dRj1l1 = j1.DeltaR(mu1); dRj1l2 = j1.DeltaR(mu2); dRj1l3 = j1.DeltaR(mu3);
			dRj2l1 = j2.DeltaR(mu1); dRj2l2 = j2.DeltaR(mu2); dRj2l3 = j2.DeltaR(mu3);
			dRj1j2 = j1.DeltaR(j2);
			for (const auto j: jets_cleaned) {
				HT += j.Pt();
				ST += j.Pt();
			}
			for (const auto mu: muons_tight)
				ST += mu.Pt();
			HToverST = HT / ST;
			MET = METv.Pt();
			tree->Fill();
		}
	}	
}

//================================================================================//
//================================================================================//
//=======//////////===//==========//===/////////===///======//===//////////=======//
//=======//============//========//====//==========//=//====//=======//===========//
//=======//=============//======//=====//==========//==//===//=======//===========//
//=======//////////======//====//======/////////===//===//==//=======//===========//
//=======//===============//==//=======//==========//====//=//=======//===========//
//=======//================////========//==========//=====////=======//===========//
//=======//////////=========//=========/////////===//======///=======//===========//
//================================================================================//
//================================================================================//
//selection...

TString Preselection::RegionSelector(
		Event& ev,
		vector<Muon> &muons_tight, vector<Electron> &electrons_tight,
		vector<Muon> &muons_loose, vector<Electron> &electrons_loose,
		vector<Jet> &jets, vector<Jet> &bjets) {

	TString region;
	// devide by leptons first
	if (muons_tight.size() == 2 && muons_loose.size() == 2 && electrons_tight.size() == 0 && electrons_loose.size() == 0) 
		region = "dimu";
	else if (muons_tight.size() == 3 && muons_loose.size() == 3 && electrons_tight.size() == 0 && electrons_loose.size() == 0) 
		region = "3mu";
	else if (muons_tight.size() == 1 && muons_loose.size() == 1 && electrons_tight.size() == 1 && electrons_loose.size() == 1)
		region = "emu";
	else if (muons_tight.size() == 2 && muons_loose.size() == 2 && electrons_tight.size() == 1 && electrons_loose.size() == 1) 
		region = "1e2mu";
	else
		return "";

	// dimu regions
	if (region == "dimu") {
		FillCutflow("DY_dimu", region);
		FillCutflow("TT_dimu", region);
		const Muon& lead = muons_tight.at(0);
		const Muon& sub = muons_tight.at(1);
		if (!ev.PassTrigger(trigs_dimu)) 
			return "";
		FillCutflow("DY_dimu", "trigger");
		FillCutflow("TT_dimu", "trigger");

		// safe pt cut
		if (! (lead.Pt() > 20.))
			return "";
		if (! (sub.Pt() > 10.))
			return "";
		FillCutflow("DY_dimu", "passSafePtCut");
		FillCutflow("TT_dimu", "passSafePtCut");

		if (lead.Charge() + sub.Charge() != 0)
			return "";
		FillCutflow("DY_dimu", "os_dimu");
		FillCutflow("TT_dimu", "os_dimu");

		// divide DY and TT
		Particle ZCand = lead + sub;
		if (fabs(ZCand.M() - 91.2) < 15.) {
			FillCutflow("DY_dimu", "OnshellZ");

			if (bjets.size() > 0)
				return "";
			FillCutflow("DY_dimu", "NoBjet");
			return "DY_dimu";
		}
		else {
			FillCutflow("TT_dimu", "OffshellZ");

			if (ZCand.M() < 12)
				return "";
			FillCutflow("TT_dimu", "mass_ge12");
			
			if (lead.DeltaR(sub) < 0.4)
				return "";
			FillCutflow("TT_dimu", "dR_ge04");

			if (jets.size() < 2)
				return "";
			FillCutflow("TT_dimu", "Nj_ge2");

			if (bjets.size() < 1)
				return "";
			FillCutflow("TT_dimu", "Nb_ge1");
			return "TT_dimu";
		}
	}
	// emu
	else if (region == "emu") {
		FillCutflow("TT_emu", region);
		const Muon& mu = muons_tight.at(0);
		const Electron& ele = electrons_tight.at(0);
		if (!ev.PassTrigger(trigs_emu))
			return "";
		FillCutflow("TT_emu", "trigger");

		bool passSafeCut = false;
		if (mu.Pt() > 10. && ele.Pt() > 25.)
			passSafeCut = true;
		if (mu.Pt() > 25. && ele.Pt() > 15.)
			passSafeCut = true;
		if (!passSafeCut)
			return "";
		FillCutflow("TT_emu", "passSafePtCut");

		if (mu.Charge() + ele.Charge() != 0)
			return "";
		FillCutflow("TT_emu", "os_emu");

		if (mu.DeltaR(ele) < 0.4)
			return "";
		FillCutflow("TT_emu", "dR_ge04");

		if (jets.size() < 2)
			return "";
		FillCutflow("TT_emu", "Nj_ge2");

		if (bjets.size() < 1)
			return "";
		FillCutflow("TT_emu", "Nb_ge1");
		return "TT_emu";
	}
	//trimu
	else if (region == "3mu") {
		FillCutflow("Pre_3mu", region);

		if (! ev.PassTrigger(trigs_dimu)) 
			return "";
		FillCutflow("Pre_3mu", "trigger");
		
		const Muon& mu1 = muons_tight.at(0);
		const Muon& mu2 = muons_tight.at(1);
		const Muon& mu3 = muons_tight.at(2);
		if (mu1.Pt() < 20.)
			return "";
		if (mu2.Pt() < 10.)
			return "";
		if (mu3.Pt() < 10.)
			return "";
		FillCutflow("Pre_3mu", "passSafePtCut");

		const int chargeSum = mu1.Charge() + mu2.Charge() + mu3.Charge();
		if (abs(chargeSum) != 1)
			return "";
		FillCutflow("Pre_3mu", "exist_osmu");

		const Particle ZCand1 = mu1 + mu2;
		const Particle ZCand2 = mu2 + mu3;
		const Particle ZCand3 = mu1 + mu3;
		if (mu1.Charge() + mu2.Charge() == 0 && ZCand1.M() < 12.)
			return "";
		if (mu2.Charge() + mu3.Charge() == 0 && ZCand2.M() < 12.)
			return "";
		if (mu1.Charge() + mu3.Charge() == 0 && ZCand3.M() < 12.)
			return "";
		FillCutflow("Pre_3mu", "mass_ge12");

		if (jets.size() < 2)
			return "";
		FillCutflow("Pre_3mu", "Nj_ge2");

		if (bjets.size() < 1)
			return "";
		FillCutflow("Pre_3mu", "Nb_ge1");
		return "Pre_3mu";
	}
	// 1e2mu
	else if (region == "1e2mu") {
		FillCutflow("Pre_1e2mu", region);

		if (! ev.PassTrigger(trigs_emu))
			return "";
		FillCutflow("Pre_1e2mu", "trigger");

		const Muon& lead = muons_tight.at(0);
		const Muon& sub = muons_tight.at(1);
		const Electron& ele = electrons_tight.at(0);
		bool passSafeCut = false;
		if (lead.Pt() > 10. && ele.Pt() > 25.)
			passSafeCut = true;
		if (lead.Pt() > 25. && ele.Pt() > 15.)
			passSafeCut = true;
		if (!passSafeCut)
			return "";
		FillCutflow("Pre_1e2mu", "passSafePtCut");

		if (lead.Charge() + sub.Charge() != 0)
			return "";
		FillCutflow("Pre_1e2mu", "exist_osmu");

		const Particle ZCand = lead + sub;
		if (ZCand.M() < 12)
			return "";
		FillCutflow("Pre_1e2mu", "mass_ge12");

		if (jets.size() < 2)
			return "";
		FillCutflow("Pre_1e2mu", "Nj_ge2");

		if (bjets.size() < 1)
			return "";
		FillCutflow("Pre_1e2mu", "Nb_ge1");
		return "Pre_1e2mu";
	}
	else
		return "";
}
		
