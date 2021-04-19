#include "SignalStudy.h"

void SignalStudy::initializeAnalyzer(){
	// flags
	Skim3Mu = HasFlag("Skim3Mu");
	Skim1E2Mu = HasFlag("Skim1E2Mu");
	cout << "[SignalStudy::initializeAnalyzer] Skim3Mu = " << Skim3Mu << endl;
	cout << "[SignalStudy::initializeAnalyzer] Skim1E2Mu = " << Skim1E2Mu << endl;

	// Triggers
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

    // Jet tagger
    vector<JetTagging::Parameters> jtps;
    jtps.emplace_back(
            JetTagging::Parameters(JetTagging::DeepJet, JetTagging::Medium, JetTagging::incl, JetTagging::mujets));
    mcCorr->SetJetTaggingParameters(jtps);		
}

void SignalStudy::executeEvent(){
	
	if (!PassMETFilter()) return;
	Event ev = GetEvent();

	// Object definition
	vector<Gen> gens = GetGens();
	vector<Muon> muons = GetAllMuons();
	vector<Electron> electrons = GetAllElectrons();
	vector<Jet> jets = GetAllJets();

	// sort
	sort(muons.begin(), muons.end(), PtComparing);
	sort(electrons.begin(), electrons.end(), PtComparing);
	sort(jets.begin(), jets.end(), PtComparing);

	//vector<Muon> muons_loose = SelectMuons(muons, "POGMedium", 8., 2.4);
	//vector<Muon> muons_tight = SelectMuons(muons, "POGLoose", 8., 2.4);
	vector<Muon> muons_loose = SelectMuons(muons, "HcToWALoose", 10., 2.4);
	vector<Muon> muons_tight = SelectMuons(muons, "HcToWATight", 10., 2.4);
	//vector<Electron> electrons_loose = SelectElectrons(electrons, "passMVAID_noIso_WP90", 8., 2.5);
	//vector<Electron> electrons_tight = SelectElectrons(electrons, "passMVAID_noIso_WPLoose", 8., 2.5);
	vector<Electron> electrons_loose = SelectElectrons(electrons, "HcToWALoose", 10., 2.5);
	vector<Electron> electrons_tight = SelectElectrons(electrons, "HcToWATight", 10., 2.5);
	vector<Jet> jets_tight = SelectJets(jets, "tight", 10., 2.4);
	vector<Jet> jets_lepVeto = JetsVetoLeptonInside(jets_tight, electrons_loose, muons_loose, 0.4);
	vector<Jet> bjets_lepVeto;

	const double bcut = mcCorr->GetJetTaggingCutValue(JetTagging::DeepJet, JetTagging::Medium);
	for (const auto& jet: jets_lepVeto) {
		const double this_discr = jet.GetTaggerResult(JetTagging::DeepJet);
		if (this_discr > bcut) {
			bjets_lepVeto.emplace_back(jet);
		}
	}

	// Further classify leptons based on generator level info
	vector<Muon> muons_signal, muons_ewprompt, muons_offshellW, muons_fromtau, muons_fake;
	vector<Electron> electrons_signal, electrons_ewprompt, electrons_offshellW, 
					 electrons_fromtau, electrons_conv, electrons_fake;

	for (const auto &mu: muons_tight) {
        int LepType = GetLeptonType(mu, gens);
        if (LepType == 1)
            muons_ewprompt.emplace_back(mu);
        else if (LepType == 2)
            muons_signal.emplace_back(mu);
        else if (LepType == 6)
            muons_offshellW.emplace_back(mu);
        else if (LepType == 3)
            muons_fromtau.emplace_back(mu);
        else if (LepType < 0)
            muons_fake.emplace_back(mu);
        else
            continue;
    }
    for (const auto &ele: electrons_tight) {
        int LepType = GetLeptonType(ele, gens);
        if (LepType == 1)
            electrons_ewprompt.emplace_back(ele);
        else if (LepType == 2)
            electrons_signal.emplace_back(ele);
        else if (LepType == 6)
            electrons_offshellW.emplace_back(ele);
        else if (LepType == 3)
            electrons_fromtau.emplace_back(ele);
        else if (LepType == 4 || LepType == 5)
            electrons_conv.emplace_back(ele);
        else if (LepType < 0)
            electrons_fake.emplace_back(ele);
        else
            continue;
    }
	
	// set weight
	double weight = 1.;
    double w_gen, w_lumi;
    w_gen = ev.MCweight()*weight_norm_1invpb;
    w_lumi = ev.GetTriggerLumi("Full");
	weight = w_gen * w_lumi;

	
	// Event selection
	if (Skim1E2Mu) {
		if (! (electrons_tight.size() == 1 && muons_tight.size() == 2)) return;
		if (! (electrons_loose.size() == 1 && muons_tight.size() == 2)) return;
		if (! (jets_lepVeto.size() >= 2)) return;
		if (! (bjets_lepVeto.size() >= 1)) return;

		if (! (ev.PassTrigger(trigs_dblmu) || ev.PassTrigger(trigs_emu))) return;

		// safe pt cut
		if (! ((muons_tight.at(0).Pt() > 20. && muons_tight.at(1).Pt() > 10.) ||
			   (muons_tight.at(0).Pt() > 25. && electrons_tight.at(0).Pt() > 15.) ||
			   (muons_tight.at(0).Pt() > 10. && electrons_tight.at(0).Pt() > 25.)))
			return;

		// dimuon mass cut
		const Particle ACand = muons_tight.at(0) + muons_tight.at(1);
		if (! (ACand.M() > 12.)) return;
		
		// Fill Hists
		TString channel = "1e2mu";	
		FillObjects(channel + "/muons_tight", muons_tight, weight);
		FillObjects(channel + "/muons_signal", muons_signal, weight);
		FillObjects(channel + "/muons_ewprompt", muons_ewprompt, weight);
		FillObjects(channel + "/muons_offshellW", muons_offshellW, weight);
		FillObjects(channel + "/muons_fake", muons_fake, weight);
		FillObjects(channel + "/electrons_tight", electrons_tight, weight);
		FillObjects(channel + "/electrons_signal", electrons_signal, weight);
		FillObjects(channel + "/electrons_offshellW", electrons_offshellW, weight);
		FillObjects(channel + "/electrons_conv", electrons_conv, weight);
		FillObjects(channel + "/electrons_fake", electrons_fake, weight);
		FillObjects(channel + "/jets_lepVeto", jets_lepVeto, weight);
		FillObjects(channel + "/bjets_lepVeto", bjets_lepVeto, weight);
		FillObject(channel+"/ACand", ACand, weight);
	}
	if (Skim3Mu) {
		cerr << "[SignalStudy::ExecuteEvent] Skim3Mu not set yet" << endl;
		exit(EXIT_FAILURE);
	}


}

SignalStudy::SignalStudy(){

}

SignalStudy::~SignalStudy(){

}


