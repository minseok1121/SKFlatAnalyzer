#include "fr_muon.h"

void fr_muon::initializeAnalyzer(){
	//==== Flags for systematic sources ====
    RunSysts = HasFlag("RunSysts");
    RunPrompts = HasFlag("RunPrompts");
    RunXSecSyst = HasFlag("RunXSecSyst");
    RunNPV = HasFlag("RunNPV");

    cout << "[fr_muon::initializeAnalyzer] RunSysts = " << RunSysts << endl;
    cout << "[fr_muon::initializeAnalyzer] RunPrompts = " << RunPrompts << endl;
    cout << "[fr_muon::initializeAnalyzer] RunXSecSyst = " << RunXSecSyst << endl;
    cout << "[fr_muon::initializeAnalyzer] RunNPV = " << RunNPV << endl;

    //==== Information for histnames & systematic sources ====
    Systs = {"Central", "JetPtCut30", "JetPtCut50", "JetPtCut60", "HadFlavor"};
    Prompts = {"Central", "JetResUp", "JetResDown", "JetEnUp", "JetEnDown",
                    "ElectronResUp", "ElectronResDown", "ElectronEnUp", "ElectronEnDown", "PileUp"};
    Regions = {"QCDEnriched", "WEnriched", "ZEnriched"};

    //==== Muon ID setting ====
    MuonIDs = {"fr_muon_loose", "fr_muon_tight"};

    //==== Trigger Settings ====
    if (DataYear == 2016) {
        HLTMuonTriggerName[0] = "HLT_Mu8_TrkIsoVVL_v";
        HLTMuonTriggerName[1] = "HLT_Mu17_TrkIsoVVL_v";
        TriggerSafePtCut[0] = 10.;
        TriggerSafePtCut[1] = 20.;
    }
    else if (DataYear == 2017 || DataYear == 2018) {
        cout << "[fr_muon::initializeAnalyzer] Trigger is not set for 2017 and 2018" << endl;
        exit(EXIT_FAILURE);
    }
    else {
        cout << "[fr_muon::initializeAnalyzer] Wrong Year" << endl;
        exit(EXIT_FAILURE);
    }

	cout << "[fr_muon::initializeAnalyzer] HLTMuonTriggerName[0] = " << HLTMuonTriggerName[0] << endl;
    cout << "[fr_muon::initializeAnalyzer] HLTMuonTriggerName[1] = " << HLTMuonTriggerName[1] << endl;
    cout << "[fr_muon::initializeAnalyzer] TriggerSafePtCut[0] = " << TriggerSafePtCut[0] << endl;
    cout << "[fr_muon::initializeAnalyzer] TriggerSafePtCut[1] = " << TriggerSafePtCut[1] << endl;

    // ==== B-tagging ====
    std::vector<JetTagging::Parameters> jtps;
    jtps.push_back( JetTagging::Parameters(JetTagging::DeepCSV, JetTagging::Medium, JetTagging::incl, JetTagging::comb) );
    mcCorr->SetJetTaggingParameters(jtps);

    cout << "[fr_muon::initializeAnalyzer] finish initialization" << endl;

    // nPV setup
}

void fr_muon::executeEvent(){
	//==== Copy all objects ====
    AllMuons = GetAllMuons();
    AllElectrons = GetAllElectrons();
    AllJets = GetAllJets();

    //==== Get L1Prefire. pileup reweight ====
    weight_Prefire = GetPrefireWeight(0);
    weight_PileUp = GetPileUpWeight(nPileUp, 0);

    AnalyzerParameter param;

	//==== Loop over muon IDs, syst, prompt ====
    for (unsigned int i = 0; i < MuonIDs.size(); i++) {
        MuonID = MuonIDs.at(i);
        MuonLooseID = "fr_muon_loose";
        ElectronLooseID = "fr_elec_loose";
        JetID = "tight";
        JetPtCut = GetJetPtCut("Central");

        syst = Systs.at(0);
        prompt = Prompts.at(0);

        param.Clear();
        param.syst_ = AnalyzerParameter::Central;
        param.Name = MuonID + "_" + syst + "_" + prompt;

        executeEventFromParameter(param);

        // RunSysts, RunPrompts, RunXSecSyst,
    }

}

void fr_muon::executeEventFromParameter(AnalyzerParameter param){

  	if(!PassMETFilter()) return;

  	Event ev = GetEvent();
	Particle METv = ev.GetMETVector();

	//==== Trigger ====
    if (! (ev.PassTrigger(HLTMuonTriggerName[0])) && !(ev.PassTrigger(HLTMuonTriggerName[1]))) return;

    //==== Copy all objects ====
    vector<Muon> this_AllMuons = AllMuons;
    vector<Electron> this_AllElectrons = AllElectrons;
    vector<Jet> this_AllJets = AllJets;

    //==== Normalization Systematic Sources ====
    if (param.syst_ == AnalyzerParameter::Central) {}
    else if (param.syst_ == AnalyzerParameter::JetResUp) {
        this_AllJets = SmearJets( this_AllJets, +1);
    }
    else if (param.syst_ == AnalyzerParameter::JetResDown) {
        this_AllJets = SmearJets( this_AllJets, -1);
    }
    else if (param.syst_ == AnalyzerParameter::JetEnUp) {
        this_AllJets = ScaleJets( this_AllJets, +1);
    }
    else if (param.syst_ == AnalyzerParameter::JetEnDown) {
        this_AllJets = ScaleJets( this_AllJets, -1);
    }
    else if (param.syst_ == AnalyzerParameter::ElectronResUp) {
        this_AllElectrons = SmearElectrons( this_AllElectrons, +1);
    }
    else if (param.syst_ == AnalyzerParameter::ElectronResDown) {
        this_AllElectrons = SmearElectrons( this_AllElectrons, -1);
    }
    else if (param.syst_ == AnalyzerParameter::ElectronEnUp) {
        this_AllElectrons = ScaleElectrons( this_AllElectrons, +1);
    }
    else if (param.syst_ == AnalyzerParameter::ElectronEnDown) {
        this_AllElectrons = ScaleElectrons( this_AllElectrons, -1);
    }

	//==== baseline selection ====
    muons = SelectMuons(this_AllMuons, MuonID, 10., 2.4);
    muons_loose = SelectMuons(this_AllMuons, MuonLooseID, 10., 2.4);
    electrons_loose = SelectElectrons(this_AllElectrons, ElectronLooseID, 15, 2.5);
    jets = SelectJets(this_AllJets, JetID, 20., 2.4);

    std::sort(muons.begin(), muons.end(), PtComparing);
    std::sort(muons_loose.begin(), muons_loose.end(), PtComparing);
    std::sort(electrons_loose.begin(), electrons_loose.end(), PtComparing);
    std::sort(jets.begin(), jets.end(), PtComparing);

    int NBjets_NoSF = 0;
    int NBjets_WithSF_2a = 0;
    JetTagging::Parameters jtp_DeepCSV_Medium = JetTagging::Parameters(JetTagging::DeepCSV, JetTagging::Medium, JetTagging::incl, JetTagging::comb);
    //double btagWeight = mcCorr->GetBTaggingReweight_1a(jets, jtp_DeepCSV_Medium);

    //==== leading muon should pass trigger-safe pt cut
    if (muons.size() == 0) return;
    double corrPt = GetCorrPt(muons.at(0));

    if (corrPt < 35 && muons.at(0).Pt() <= TriggerSafePtCut[0]) return;
    else if (corrPt > 35 && muons.at(0).Pt() <= TriggerSafePtCut[1]) return;

	//==== measure nPV reweight ====
    //==== 1. no loose electron
    //==== 2. exactly 1 muon
    //==== 3. jet - ptcut 40

    if (RunNPV) {
        jets_dR04 = JetsVetoLeptonInside(jets, electrons_loose, muons, 0.4);

        if (electrons_loose.size() != 0) return;
        if (muons.size() != 1) return;
        if (jets_dR04.size() == 0) return;
        if (jets_dR04.at(0).Pt() < JetPtCut) return;

        //==== GetWeights ====
        double weight = 1.;
        if (!IsDATA) {
            if (corrPt < 35) weight *= weight_norm_1invpb*ev.GetTriggerLumi(HLTMuonTriggerName[0]);
            else if (corrPt > 35) weight *= weight_norm_1invpb*ev.GetTriggerLumi(HLTMuonTriggerName[1]);
            weight *= ev.MCweight();
            weight *= weight_Prefire;
        }

        if (nPV > 100) nPV = 100;
        FillHist("nPV_" + param.Name, nPV, weight, 100, 0, 100);

        return;
    }

	/////////////////////////////////////////////////////////////////////
    //==== Event Selection
    //==== QCD dominated region
    //==== 1. Exactly 1 electron
    //==== 2. At least 1 jet with Pt > 40 GeV
    //==== 3. Delta(e, j) > 1.0
    //==== 4. MET < 25 GeV, MT(e, METv) < 25 GeV
    /////////////////////////////////////////////////////////////////////
    //==== W boson dominated region
    //==== 1. Exactly 1 electron
    //==== 2. At least 1 jet with Pt > 40 GeV
    //==== 3. Delta(e, j) > 0.4
    //==== 4. MET > 50 GeV & MT(e, METv) > 70 GeV
    //////////////////////////////////////////////////////////////////////
    //==== Z boson dominated region
    //==== 1. Exactly 2 electrons
    //==== 2, At least 1 jet with Pt > 40 GeV
    //==== 3. Delta(e, j ) > 0.4
    //==== 4. |M(ee) - 91.2| < 15
    //////////////////////////////////////////////////////////////////////
	

}

fr_muon::fr_muon(){

}

fr_muon::~fr_muon(){

}

double fr_muon::GetCorrPt(const Muon &mu) {
	double corrPt;
	double weight;
	double rIsoTight = 0.2;

	if (mu.RelIso() > rIsoTight) {
		weight = 1. + mu.RelIso() - rIsoTight;
		corrPt = mu.Pt() * weight;
	}
	else corrPt = mu.Pt();

	return corrPt;
}

double fr_muon::GetJetPtCut(const TString syst) {
	if (syst == Systs.at(0) || syst == Systs.at(4) ) return 40;
    else if (syst == Systs.at(1)) return 30;
    else if (syst == Systs.at(2)) return 50;
    else if (syst == Systs.at(3)) return 60;
    else {
        cout << "[fr_muon::GetJetPtCUt] Wrong Syst" << endl;
        exit(EXIT_FAILURE);
    }	
}
