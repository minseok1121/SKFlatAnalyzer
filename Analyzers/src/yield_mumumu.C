#include "yield_mumumu.h"

void yield_mumumu::initializeAnalyzer(){
	
	ElectronIDs = {"fr_elec_loose", "fr_elec_tight"};
	MuonIDs = {"fr_muon_loose", "fr_muon_tight"};

	//==== Trigger Settings ====
	if (DataYear == 2016) {
		TrigList_DblMu_BtoG = {"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v", "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v"};
		TrigList_DblMu_H = {"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v", "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v"};
		//if (IsDATA && (run > 280385)) HLTMuonTriggerNames = TrigList_DblMu_H;
		//else HLTMuonTriggerNames = TrigList_DblMu_BtoG;
		TriggerSafePtCut1 = 20;
		TriggerSafePtCut2 = 10;
	}
	else if (DataYear == 2017 || DataYear == 2018) {
		cout << "[yield_mumumu::initializeAnalyzer] trigger not set" << endl;
		exit(EXIT_FAILURE);
	}
	else {
		cout << "[yield_mumumu::initializeAnalyzer] Wrong Year" << endl;
		exit(EXIT_FAILURE);
	}

	cout << "[yield_mumumu::initializeAnalyzer] TriggerSafePtCut1 = " << TriggerSafePtCut1 << endl;
	cout << "[yield_mumumu::initializeAnalyzer] TriggerSafePtCut2 = " << TriggerSafePtCut2 << endl;

	// ==== B-tagging ====
    std::vector<JetTagging::Parameters> jtps;
    jtps.push_back( JetTagging::Parameters(JetTagging::CSVv2, JetTagging::Medium, JetTagging::incl, JetTagging::mujets) );
    mcCorr->SetJetTaggingParameters(jtps);

	cout << "[yield_mumumu::initializeAnalyzer] finish initialization" << endl;	
}

void yield_mumumu::executeEvent(){
	
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
	param.Name = "mumumu_Central";


  	executeEventFromParameter(param);

}

void yield_mumumu::executeEventFromParameter(AnalyzerParameter param){
	
	bool IsPeriodH = (IsDATA && (run > 280385));
	if (IsPeriodH) HLTMuonTriggerNames = TrigList_DblMu_H;
	else HLTMuonTriggerNames = TrigList_DblMu_BtoG;

	//==== No cut ====
	FillHist("NoCut_" + param.Name, 0., 1., 1, 0., 1.);
	if(!PassMETFilter()) return;

	Event ev = GetEvent();
	Particle METv = ev.GetMETVector();

	//==== Trigger ====
	if (! ev.PassTrigger(HLTMuonTriggerNames)) return;

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
	
	FillHist("nMuons_" + param.Name, muons.size(), 1., 5, 0., 5.);
	FillHist("nMuons_loose_" + param.Name, muons_loose.size(), 1., 5, 0., 5.);
	FillHist("nElectrons_" + param.Name, electrons.size(), 1., 5, 0., 5.);
	FillHist("nElectrons_loose_" + param.Name, electrons_loose.size(), 1., 5, 0., 5.);

	//==== leading, subleading muons should pass trigger safe pt cut ====
	if (muons_loose.size() < 3) return;
	if (muons_loose.at(0).Pt() <= TriggerSafePtCut1) return;
	if (muons_loose.at(1).Pt() <= TriggerSafePtCut2) return;
	if (muons_loose.at(2).Pt() <= TriggerSafePtCut2) return;

	//==== mumumu events with no other leptons ====
	//==== loose criteria first =====
	if (electrons_loose.size() != 0) return;
	if (muons_loose.size() !=3) return;

	bool tightFlag = false;
	if (muons.size() == 3) tightFlag = true;

	FillHist("CommonLeptonVeto_" + param.Name, 0, 1., 1, 0., 1.);
	if(tightFlag) FillHist("CommonLeptonVeto_tight_" + param.Name, 0, 1., 1, 0., 1.);

	//==== opposite sign pair with mass constraint ====
	Pair muonPair[3];
	for (int i = 0; i < 3; i++) {
		muonPair[i].Mass = (muons_loose.at(i%3) + muons_loose.at((i+1)%3)).M();
		muonPair[i].isOS = false;
		muonPair[i].isAbove12 = false;
		muonPair[i].isOnZ = false;
		muonPair[i].isOffZ = false;

		if (muons_loose.at(i%3).Charge() * muons_loose.at((i+1)%3).Charge() < 0) muonPair[i].isOS = true;
		if (muonPair[i].Mass > 12) muonPair[i].isAbove12 = true;
		if (IsOnZ(muonPair[i].Mass, 10.)) muonPair[i].isOnZ = true;
		if (!IsOnZ(muonPair[i].Mass, 10.)) muonPair[i].isOffZ = true;
	}
	bool is3lOnZ = false;
	Particle ZCand = muons_loose.at(0) + muons_loose.at(1) + muons_loose.at(2);
	if (IsOnZ(ZCand.M(), 10.)) is3lOnZ = true;

	//==== CR1. ZGamma region ====
	//==== 1. 3 muons with 2 OS muon pair, no additional leptons
	//==== 2. offZ for all OS pairs
	//==== 3. include 12 < M(pair) < 81.2
	//==== 4. M(3l) should be onZ
	//==== 5. MET < 50
	bool IsCRZGamma = true;
	int CR1OSCnt = 0;
	int CR1OffZCnt = 0;
	int CR1AbCnt = 0;
	int CR1BelCnt = 0;
	for (int i = 0; i < 3; i++) {
		if (muonPair[i].isOS) CR1OSCnt++;
		if (muonPair[i].isOS && muonPair[i].isOffZ) CR1OffZCnt++;
		if (muonPair[i].isOS && (muonPair[i].Mass > 12)) CR1AbCnt++;
		if (muonPair[i].isOS && (muonPair[i].Mass < 81.2)) CR1BelCnt++;
	}
	if (! (CR1OSCnt == 2)) IsCRZGamma = false;
	if (! (CR1OffZCnt == 2)) IsCRZGamma = false;
	if (! (CR1AbCnt == 2)) IsCRZGamma = false; // why not all pairs???
	if (! (CR1BelCnt > 0)) IsCRZGamma = false;
	if (! is3lOnZ) IsCRZGamma = false;
	if (! (METv.Pt() < 50)) IsCRZGamma = false;

	if (IsCRZGamma) FillHist("PassCRZGamma_" + param.Name, 0, 1., 1, 0., 1.);
	if (IsCRZGamma && tightFlag) FillHist("PassCRZGamma_tight_" + param.Name, 0, 1., 1, 0., 1.);

	//==== CR2. OnZ trilepton region ====
	//==== 1. 3 muons with OS pairs, no additional leptons
	//==== 2. for all OS pair, M(mumu) > 12
	//==== 3. at least one OS pair onZ
	bool IsCR3lOnZ = true;
	int CR2OSCnt = 0;
	int CR2AbCnt = 0;
	int CR2OnZCnt = 0;
	for (int i = 0; i < 3; i++) {
		if (muonPair[i].isOS) CR2OSCnt++;
		if (muonPair[i].isOS && muonPair[i].isAbove12) CR2AbCnt++;
		if (muonPair[i].isOS && muonPair[i].isOnZ) CR2OnZCnt++;
	}
	if (! (CR2OSCnt == 2)) IsCR3lOnZ = false;
	if (! (CR2AbCnt == 2)) IsCR3lOnZ = false;
	if (! (CR2OnZCnt > 0)) IsCR3lOnZ = false;

	if (IsCR3lOnZ) FillHist("PassCR3lOnZ_" + param.Name, 0, 1., 1, 0., 1.);
	if (IsCR3lOnZ && tightFlag) FillHist("PassCROnZ3l_tight_" + param.Name, 0, 1., 1, 0., 1.);
	
	//==== jets ====
	jets_dR04 = JetsVetoLeptonInside(jets, electrons_loose, muons_loose, 0.4);
	if (jets_dR04.size() < 2 ) return;
	for (unsigned int i = 0; i < jets_dR04.size(); i++) {
		double this_discr = jets_dR04.at(i).GetTaggerResult(JetTagging::CSVv2);
		if (this_discr > mcCorr->GetJetTaggingCutValue(JetTagging::CSVv2, JetTagging::Medium)) NBjets_NoSF++;
		if (mcCorr->IsBTagged_2a(jtp_CSVv2_Medium, jets_dR04.at(i))) NBjets_WithSF_2a++;
	}

	if (IsData && NBjets_NoSF == 0) return;
	if (!IsData && NBjets_WithSF_2a == 0 ) return;

	//==== SRjj ====
	bool IsSRjj = true;
	int SRjjOSCnt = 0;
	int SRjjAbCnt = 0;
	int SRjjOffZCnt =0;
	for (int i = 0; i< 3; i++) {
		if (muonPair[i].isOS) SRjjOSCnt++;
		if (muonPair[i].isOS && muonPair[i].isAbove12) SRjjAbCnt++;
		if (muonPair[i].isOS && muonPair[i].isOffZ) SRjjOffZCnt++;
	}
	if (SRjjOSCnt == 0) IsSRjj = false;
	if (SRjjAbCnt != 2) IsSRjj = false;
	if (SRjjOffZCnt != 2) IsSRjj = false;
	if (IsSRjj) FillHist("PassSRjj_" + param.Name, 0, 1, 1, 0., 1.); 
	if (IsSRjj && tightFlag) FillHist("PassSFjj_tight_" + param.Name, 0, 1, 1, 0., 1.);

	//==== SRjjM80 ====
	bool IsSRjjM80 = IsSRjj;
	Particle ACand = ChooseACand(muons_loose, METv);
	if (! (ACand.M() < 80)) IsSRjjM80 = false;
	
	if (IsSRjjM80) FillHist("PassSRjjM80_" + param.Name, 0, 1, 1, 0., 1.);
	if (IsSRjjM80 && tightFlag) FillHist("PassSRjjM80_tight_" + param.Name, 0, 1, 1, 0., 1.);

	return;
}

yield_mumumu::yield_mumumu(){

}

yield_mumumu::~yield_mumumu(){

}

Particle yield_mumumu::ChooseACand(const vector<Muon> &muons, const Particle &METv) {
	if (muons.size() != 3) exit(EXIT_FAILURE);
	if ((fabs(muons.at(0).Charge()+muons.at(1).Charge()+muons.at(2).Charge()) > 1)) {
		return Particle(0, 0, 0, 1000);
	}
		
		
	// sort muons
	Muon Acand1, Acand2, OScand;
	for (int i = 0; i < 3; i++) {
		if (muons.at(i%3).Charge() == muons.at((i+1)%3).Charge()) {
			if (i == 2) {
				Acand1 = muons.at(0);
				Acand2 = muons.at(2);
				OScand = muons.at(1);
			}
			else {
				Acand1 = muons.at(i%3);
				Acand2 = muons.at((i+1)%3);
				OScand = muons.at((i+2)%3);
			}
			break;
		}
		else continue;
	}

	// choose muons from A
	double dPt = fabs(Acand1.Pt() - Acand2.Pt());
	double MT1 = MT(Acand1, METv), MT2 = MT(Acand2, METv);
	bool fromW1 = (50 < MT1 && MT1 < 120);
	bool fromW2 = (50 < MT2 && MT2 < 120);
	if (dPt < 25 && (!fromW1) && fromW2) return Acand1 + OScand;
	else return Acand2 + OScand;
}









