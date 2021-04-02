#include "Selector.h"

Selector::Selector() {}
Selector::~Selector() {
	if (SkimSR3mu || SkimSR1e2mu || SkimWZ3mu || SkimWZ1e2mu) {
		outfile->cd();
		tree->Write();
	}
}

void Selector::initializeAnalyzer(){
 	// flags
	RunSysts = HasFlag("RunSysts");
    RunDeepCSV = HasFlag("RunDeepCSV");
    EMuTrigOnly = HasFlag("EMuTrigOnly");
    SkimSR3mu = HasFlag("SkimSR3mu");
	SkimSR1e2mu = HasFlag("SkimSR1e2mu");
	SkimWZ3mu = HasFlag("SkimWZ3mu");
	SkimWZ1e2mu = HasFlag("SkimWZ1e2mu");
	cout << "[Selector::initializeAnalyzer] RunDeepCSV = " << RunDeepCSV << endl;
    cout << "[Selector::initializeAnalyzer] EMuTrigOnly = " << EMuTrigOnly << endl;
    cout << "[Selector::initializeAnalyzer] SkimSR3mu = " << SkimSR3mu << endl;
	cout << "[Selector::initializeAnalyzer] SkimSR1e2mu = " << SkimSR1e2mu << endl;
	cout << "[Selector::initializeAnalyzer] SkimWZ3mu = " << SkimWZ3mu << endl;
	cout << "[Selector::initializeAnalyzer] SkimWZ1e2mu = " << SkimWZ1e2mu << endl;

	// Set TTree
	if (SkimSR3mu || SkimSR1e2mu || SkimWZ3mu || SkimWZ1e2mu) {
		tree = new TTree("Events", "Events");
		tree->Branch("run", &run);
		tree->Branch("event", &event);
		tree->Branch("lumi", &lumi);
		tree->Branch("weight", &weight);
		tree->Branch("nMuons", &nMuons);
		tree->Branch("muons_pt", muons_pt, "muons_pt[nMuons]/D");
		tree->Branch("muons_eta", muons_eta, "muons_eta[nMuons]/D");
		tree->Branch("muons_phi", muons_phi, "muons_phi[nMuons]/D");
		tree->Branch("muons_mass", muons_mass, "muons_mass[nMuons]/D");
		tree->Branch("muons_charge", muons_charge, "muons_charge[nMuons]/I");
		tree->Branch("muons_lepType", muons_lepType, "muons_lepType[nMuons]/I");
		tree->Branch("nElectrons", &nElectrons);
		tree->Branch("electrons_pt", electrons_pt, "electrons_pt[nElectrons]/D");
		tree->Branch("electrons_eta", electrons_eta, "electrons_eta[nElectrons]/D");
		tree->Branch("electrons_phi", electrons_phi, "electrons_phi[nElectrons]/D");
		tree->Branch("electrons_mass", electrons_mass, "electrons_mass[nElectrons]/D");
		tree->Branch("electrons_charge", electrons_charge, "electrons_charge[nElectrons]/I");
		tree->Branch("electrons_lepType", electrons_lepType, "electrons_lepType[nElectrons]/I");
		tree->Branch("nJets", &nJets);
		tree->Branch("jets_pt", jets_pt, "jets_pt[nJets]/D");
		tree->Branch("jets_eta", jets_eta, "jets_eta[nJets]/D");
		tree->Branch("jets_phi", jets_phi, "jets_phi[nJets]/D");
		tree->Branch("jets_isBtagged", jets_isBtagged, "jets_isBtagged[nJets]/I");
		tree->Branch("METv_pt", &METv_pt);
		tree->Branch("METv_phi", &METv_phi);
	}
	// Let's save CPU time, no need to dig all the regions
	AllRegions = false;
	DiLepOnly = false;
	TriLepOnly = true;

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
    for (const auto& region: regions) {
        InitiateCutflow(region, this->getCuts(region));
	}

	// Systematics
	Systematics = {"Central"};
	if (RunSysts)
		Systematics = {
			"Central", "JetResUp", "JetResDown", "JetEnUp", "JetEnDown", 
			"MuonEnUp", "MuonEnDown", "ElectronResUp", "ElectronResDown", 
			"ElectronEnUp", "ElectronEnDown", "IDSFUp", "IDSFDown", 
			"PUCorrUp", "PUCorrDown", "BtagUp", "BtagDown"
		};
}

void Selector::executeEvent() {
	// Get All objects
	gens = GetGens();
	muons_all = GetAllMuons();
	electrons_all = GetAllElectrons();
	jets_all = GetAllJets();
	
	for (const auto& syst : Systematics)
		executeEventWithSystematics(syst);
}
void Selector::executeEventWithSystematics(const TString& syst){
	
	// Fill All Cutflows
	for (const auto& region: regions)
		FillCutflow(region, "noCut", syst);

	if (!PassMETFilter()) return;
	for (const auto& region: regions)
		FillCutflow(region, "METFilter", syst);

	vector<Muon> muons = muons_all;
	vector<Electron> electrons = electrons_all;
	vector<Jet> jets = jets_all;
	// Resolution, Energy scale
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


	Event ev = GetEvent();
	Particle METv = ev.GetMETVector(); // phi correction applied
	// sort at the first time, don't want to be confused
	sort(muons.begin(), muons.end(), PtComparing);
	sort(electrons.begin(), electrons.end(), PtComparing);
	sort(jets.begin(), jets.end(), PtComparing);

	vector<Muon> muons_tight = SelectMuons(muons, "HcToWATight", 10., 2.4);
	vector<Muon> muons_loose = SelectMuons(muons, "HcToWALoose", 10., 2.4);
	vector<Electron> electrons_tight = SelectElectrons(electrons, "HcToWATight", 10., 2.5);
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

	const TString channel = RegionSelector(
			ev, muons_tight, electrons_tight,
			muons_loose, electrons_loose,
			jets_lepVeto, bjets_lepVeto, METv, syst);
	if (channel == "")
		return;

	// initialize weight
	weight = getWeight(channel, ev, muons_tight, electrons_tight, jets_lepVeto, syst);
		
	// Fill Hists
	FillObjects(channel+"/"+syst+"/muons_tight", muons_tight, weight);
	FillObjects(channel+"/"+syst+"/electrons_tight", electrons_tight, weight);
	FillObjects(channel+"/"+syst+"/jets_lepVeto", jets_lepVeto, weight);
	FillObjects(channel+"/"+syst+"/bjets_lepVeto", bjets_lepVeto, weight);
	FillObject(channel+"/"+syst+"/METv", METv, weight);
	if (channel == "DY_OSdimu" || channel == "TT_OSdimu") {
		const Particle ZCand = muons_tight.at(0) + muons_tight.at(1);
		FillObject(channel+"/"+syst+"/ZCand", ZCand, weight);
	}
	if ((channel == "SR_3mu" && SkimSR3mu) || (channel == "WZ_3mu" && SkimWZ3mu)) {
		nMuons = muons_tight.size();
		for (unsigned int i = 0; i < nMuons; i++) {
			muons_pt[i] = muons_tight.at(i).Pt();
			muons_eta[i] = muons_tight.at(i).Eta();
			muons_phi[i] = muons_tight.at(i).Phi();
			muons_mass[i] = muons_tight.at(i).M();
			muons_charge[i] = muons_tight.at(i).Charge();
			muons_lepType[i] = GetLeptonType(muons_tight.at(i), gens);
		}
		nElectrons = 0;
		nJets = jets_lepVeto.size();
		for (unsigned int i = 0; i < nJets; i++) {
			jets_pt[i] = jets_lepVeto.at(i).Pt();
			jets_eta[i] = jets_lepVeto.at(i).Eta();
			jets_phi[i] = jets_lepVeto.at(i).Phi();
			const double this_discr = 
				jets_lepVeto.at(i).GetTaggerResult(JetTagging::DeepJet);
			if (this_discr > mcCorr->GetJetTaggingCutValue(JetTagging::DeepJet, JetTagging::Medium))
				jets_isBtagged[i] = true;
			else
				jets_isBtagged[i] = false;
		}
		METv_pt = METv.Pt();
		METv_phi = METv.Phi();

		tree->Fill();
	}
	if ((channel == "SR_1e2mu" && SkimSR1e2mu) || (channel == "WZ_1e2mu" && SkimWZ1e2mu)) {
		nElectrons = electrons_tight.size();
		electrons_pt[0] = electrons_tight.at(0).Pt();
		electrons_eta[0] = electrons_tight.at(0).Eta();
		electrons_phi[0] = electrons_tight.at(0).Phi();
		electrons_charge[0] = electrons_tight.at(0).Charge();
		electrons_lepType[0] = GetLeptonType(electrons_tight.at(0), gens);

		nMuons = muons_tight.size();
        for (unsigned int i = 0; i < nMuons; i++) {
            muons_pt[i] = muons_tight.at(i).Pt();
            muons_eta[i] = muons_tight.at(i).Eta();
            muons_phi[i] = muons_tight.at(i).Phi();
            muons_mass[i] = muons_tight.at(i).M();
            muons_charge[i] = muons_tight.at(i).Charge();
			muons_lepType[i] = GetLeptonType(muons_tight.at(0), gens);
        }
        
		nJets = jets_lepVeto.size();
		for (unsigned int i = 0; i < nJets; i++) {
			jets_pt[i] = jets_lepVeto.at(i).Pt();
			jets_eta[i] = jets_lepVeto.at(i).Eta();
			jets_phi[i] = jets_lepVeto.at(i).Phi();
			const double this_discr = 
				jets_lepVeto.at(i).GetTaggerResult(JetTagging::DeepJet);
			if (this_discr > mcCorr->GetJetTaggingCutValue(JetTagging::DeepJet, JetTagging::Medium))
				jets_isBtagged[i] = true;
			else
				jets_isBtagged[i] = false;
		}
        
		METv_pt = METv.Pt();
        METv_phi = METv.Phi();
	
        tree->Fill();
	}
}
