#include "yield_emumu.h"

void yield_emumu::initializeAnalyzer(){
	
	ElectronIDs = {"fr_elec_loose", "fr_elec_tight"};
    MuonIDs = {"fr_muon_loose", "fr_muon_tight"};

    //==== Trigger Settings ====
    if (DataYear == 2016) {
        TrigList_MuonEG_BtoG = {"HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v"};
		TrigList_MuonEG_H = {"HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v"};

		//if (IsDATA && (run > 280385)) HLTEMuTriggerNames = TrigList_MuonEG_H;
		//else HLTEMuTriggerNames = TrigList_MuonEG_BtoG;
        TriggerSafeElecPtCut = 25;
        TriggerSafeMuPtCut = 10;
    }
    else if (DataYear == 2017 || DataYear == 2018) {
        cout << "[yield_emumu::initializeAnalyzer] trigger not set" << endl;
        exit(EXIT_FAILURE);
    }
    else {
        cout << "[yield_emumu::initializeAnalyzer] Wrong Year" << endl;
        exit(EXIT_FAILURE);
    }
	cout << "[yield_emumu::initializeAnalyzer] run = " << run << endl;
	//cout << "[yield_emumu::initializeAnalyzer] HLTEMuTriggerName = " << HLTEMuTriggerNames.at(0) << endl;
    cout << "[yield_emumu::initializeAnalyzer] TriggerSafeElecPtCut = " << TriggerSafeElecPtCut << endl;
    cout << "[yield_emumu::initializeAnalyzer] TriggerSafeMuPtCut = " << TriggerSafeMuPtCut << endl;

    // ==== B-tagging ====
    std::vector<JetTagging::Parameters> jtps;
    jtps.push_back( JetTagging::Parameters(JetTagging::CSVv2, JetTagging::Medium, JetTagging::incl, JetTagging::mujets) );
    mcCorr->SetJetTaggingParameters(jtps);

    cout << "[yield_emumu::initializeAnalyzer] finish initialization" << endl;
}

void yield_emumu::executeEvent(){

	//==== Copy all objects ====
    AllMuons = GetAllMuons();
    AllElectrons = GetAllElectrons();
    AllJets = GetAllJets();

    AnalyzerParameter param;

    MuonLooseID = MuonIDs.at(0); MuonTightID = MuonIDs.at(1);
    ElecLooseID = ElectronIDs.at(0); ElecTightID = ElectronIDs.at(1);
    JetID = "tight";
    param.Clear();
    param.syst_ = AnalyzerParameter::Central;
    param.Name = "emumu_Central";


    executeEventFromParameter(param);
}

void yield_emumu::executeEventFromParameter(AnalyzerParameter param){

	bool IsPeriodH = (IsDATA && (run > 280385));

	if (IsPeriodH) HLTEMuTriggerNames = TrigList_MuonEG_H;
	else HLTEMuTriggerNames = TrigList_MuonEG_BtoG;

	//==== No cut ====
    FillHist("NoCut_" + param.Name, 0., 1., 10, 0., 10.);

    Event ev = GetEvent();
    Particle METv = ev.GetMETVector();

    //==== Trigger ====
    if (! ev.PassTrigger(HLTEMuTriggerNames)) return;
	FillHist("CutFlow_" + param.Name, 0., 1., 10, 0., 10.);

	if(!PassMETFilter()) return;
	FillHist("CutFlow_" + param.Name, 1., 1., 10, 0., 10.);

	//==== Copy all objects ====
    vector<Muon> this_AllMuons = AllMuons;
    vector<Electron> this_AllElectrons = AllElectrons;
    vector<Jet> this_AllJets = AllJets;

    //==== baseline selection ====
    muons = SelectMuons(this_AllMuons, MuonTightID, 10., 2.4);
    muons_loose = SelectMuons(this_AllMuons, MuonLooseID, 10., 2.4);
    electrons = SelectElectrons(this_AllElectrons, ElecTightID, 25., 2.5);
    electrons_loose = SelectElectrons(this_AllElectrons, ElecLooseID, 10, 2.5);
    jets = SelectJets(this_AllJets, JetID, 25., 2.4);

	int NBjets_NoSF = 0;
    int NBjets_WithSF_2a = 0;
    JetTagging::Parameters jtp_CSVv2_Medium = JetTagging::Parameters(JetTagging::CSVv2, JetTagging::Medium, JetTagging::incl, JetTagging::comb);

    std::sort(muons.begin(), muons.end(), PtComparing);
    std::sort(muons_loose.begin(), muons_loose.end(), PtComparing);
    std::sort(electrons.begin(), electrons.end(), PtComparing);
    std::sort(electrons_loose.begin(), electrons_loose.end(), PtComparing);
    std::sort(jets.begin(), jets.end(), PtComparing);

	//==== jets ====
    jets_dR04 = JetsVetoLeptonInside(jets, electrons_loose, muons_loose, 0.4);
    for (unsigned int i = 0; i < jets_dR04.size(); i++) {
        double this_discr = jets_dR04.at(i).GetTaggerResult(JetTagging::CSVv2);
        if (this_discr > mcCorr->GetJetTaggingCutValue(JetTagging::CSVv2, JetTagging::Medium)) NBjets_NoSF++;
        if (mcCorr->IsBTagged_2a(jtp_CSVv2_Medium, jets_dR04.at(i))) NBjets_WithSF_2a++;
    }

    FillHist("nMuons_" + param.Name, muons.size(), 1., 5, 0., 5.);
    FillHist("nMuons_loose_" + param.Name, muons_loose.size(), 1., 5, 0., 5.);
	FillHist("nElectrons_" + param.Name, electrons.size(), 1., 5, 0., 5.);
	FillHist("nElectrons_loose_" + param.Name, electrons_loose.size(), 1., 5, 0., 5.);
	
	//==== leading electron/muon should pass trigger safe pt cut ====
	if (electrons_loose.size() < 1) return;
	if (muons_loose.size() < 2) return;
	if (electrons_loose.at(0).Pt() < TriggerSafeElecPtCut) return;
	if (muons_loose.at(0).Pt() < TriggerSafeMuPtCut) return;
	if (muons_loose.at(1).Pt() < TriggerSafeMuPtCut) return;

	//===== no additional leptons ====
	if (electrons_loose.size() != 1) return;
	if (muons_loose.size() != 2) return;
	FillHist("CutFlow_" + param.Name, 2., 1., 10, 0., 10.);
	

	// opposite sign muon pair
	if (muons_loose.at(0).Charge() * muons_loose.at(1).Charge() > 0) return;

	bool tightFlag = false;
	if (electrons.size() == 1 && muons.size() == 2) tightFlag = true;
	if (tightFlag) FillHist("CutFlow_" + param.Name, 3., 1., 10, 0., 10.);

	FillHist("CommomLeptonVeto_" + param.Name, 0, 1., 1, 0., 1.);
	if(tightFlag) FillHist("CommonLeptonVeto_tight_" + param.Name, 0, 1., 1, 0., 1.);

	Particle Pair = muons_loose.at(0) + muons_loose.at(1);
	Particle l3 = muons_loose.at(0) + muons_loose.at(1) + electrons_loose.at(0);

	bool IsQCDLike = Pair.M() < 12;
	bool IsZLike = IsOnZ(Pair.M(), 10);
	if (tightFlag) {
		if (!IsQCDLike) FillHist("CutFlow_" + param.Name, 4., 1., 10, 0., 10.);
		if (!IsQCDLike && !IsZLike) FillHist("CutFlow_" + param.Name, 5., 1., 10, 0., 10.);
		if (!IsQCDLike && !IsZLike && NBjets_NoSF!=0) FillHist("CutFlow_" + param.Name, 6., 1., 10, 0., 10.);
		if (!IsQCDLike && !IsZLike && NBjets_NoSF!=0 && jets_dR04.size()>1) FillHist("CutFlow_" + param.Name, 7., 1., 10, 0., 10.); 
		if (!IsQCDLike && !IsZLike && NBjets_NoSF!=0 && jets_dR04.size()>1 && Pair.M()<80) FillHist("CutFlow_" + param.Name, 8., 1., 10, 0., 10.);
	}
	
	//==== CR1. ZGamma region ====
	bool IsCRZGamma = true;
	if (Pair.M() < 12 || IsOnZ(Pair.M(), 10)) IsCRZGamma = false;
	if (!IsOnZ(l3.M(), 10)) IsCRZGamma = false;
	if (METv.Pt() > 50) IsCRZGamma = false;

	if (IsCRZGamma) FillHist("PassCRZGamma_" + param.Name, 0, 1., 1, 0., 1.);
	if (IsCRZGamma && tightFlag) FillHist("PassCRZGamma_tight_" + param.Name, 0, 1., 1, 0., 1.);

	//==== CR2. OnZ trilepton region ====
	bool IsCROnZ3l = true;
	if (Pair.M() < 12 || !IsOnZ(Pair.M(), 10)) IsCROnZ3l = false;
	
	if (IsCROnZ3l) FillHist("PassCROnZ3l_" + param.Name, 0, 1., 1, 0., 1.);
	if (IsCROnZ3l && tightFlag) FillHist("PassCROnZ3l_tight_" + param.Name, 0, 1., 1, 0., 1.);

	//==== SR ====
	if (jets_dR04.size() < 2 ) return;
    if (IsData && NBjets_NoSF == 0) return;
    if (!IsData && NBjets_WithSF_2a == 0 ) return;

    //==== SRjj ====
    bool IsSRjj = true;
	if (Pair.M() < 12) IsSRjj = false;
	if (IsOnZ(Pair.M(), 10)) IsSRjj = false;

	if (IsSRjj) FillHist("PassSRjj_" + param.Name, 0, 1, 1, 0., 1.);
	if (IsSRjj && tightFlag) FillHist("PassSRjj_tight_" + param.Name, 0, 1, 1, 0., 1.);

	//==== SRjjM80 ====
	bool IsSRjjM80 = IsSRjj;
	if (Pair.M() > 80) IsSRjjM80 = false;
	if (IsSRjjM80) FillHist("PassSRjjM80_" + param.Name, 0, 1, 1, 0., 1.);
	if (IsSRjjM80 && tightFlag) FillHist("PassSRjjM80_tight_" + param.Name, 0, 1, 1, 0., 1.);

	return;
}

yield_emumu::yield_emumu(){

}

yield_emumu::~yield_emumu(){

}




