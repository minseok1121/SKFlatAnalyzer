#include "BaselineSelector.h"

BaselineSelector::BaselineSelector() {}
BaselineSelector::~BaselineSelector() {
	if (SkimBaseline) {
		outfile->cd();
		Events->Write();
	}
}
void BaselineSelector::initializeAnalyzer(){
	// flags
	RunDeepCSV = HasFlag("RunDeepCSV");
	SkimBaseline = HasFlag("SkimBaseline");
	cout << "[BaselineSelector::initializeAnalyzer] RunDeepCSV = " << RunDeepCSV << endl;
	cout << "[BaselineSelector::initializeAnalyzer] SkimBaseline = " << SkimBaseline << endl;

	// Set branch for Events tree
	if (SkimBaseline) {
		Events = new TTree("Events", "Events");
		Events->Branch("ptl1", &ptl1); Events->Branch("ptl2", &ptl2); Events->Branch("ptl3", &ptl3);
		Events->Branch("etal1", &etal1); Events->Branch("etal2", &etal2); Events->Branch("etal3", &etal3);
		Events->Branch("phil1", &phil2); Events->Branch("phil2", &phil2); Events->Branch("phil3", &phil3);
		Events->Branch("ptj1", &ptj1); Events->Branch("ptj2", &ptj2); Events->Branch("ptb1", &ptb1);
		Events->Branch("etaj1", &etaj1); Events->Branch("etaj2", &etaj2); Events->Branch("etab1", &etab1);
		Events->Branch("phij1", &phij1); Events->Branch("phij2", &phij2); Events->Branch("phib3", &phib3);
		Events->Branch("dRl1l2", &dRl1l2); Events->Branch("dRl1l3", &dRl1l3); Events->Branch("dRl2l3", &dRl2l3);
		Events->Branch("dRj1l1", &dRj1l1); Events->Branch("dRj1l2", &dRj1l2); Events->Branch("dRj1l3", &dRj1l3);
		Events->Branch("dRj2l1", &dRj2l1); Events->Branch("dRj2l2", &dRj2l2); Events->Branch("dRj2l3", &dRj2l3);
		Events->Branch("dRj1j2", &dRj1j2);
		Events->Branch("HT", &HT); Events->Branch("LT", &LT); Events->Branch("MET", &MET); Events->Branch("ST", &ST);
		Events->Branch("HToverST", &HToverST); Events->Branch("LToverST", &LToverST); Events->Branch("METoverST", &METoverST);
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
	POG_MuID = {"POGMedium", "POGLoose"};
	HcToWA_EleID = {"HcToWATight", "HcToWALoose"}; 
	POG_EleID = {"passMVAID_noIso_WP90", "passMVAID_noIso_WPLoose"};

	// B-tagging
	vector<JetTagging::Parameters> jtps;
	if (RunDeepCSV)
		jtps.emplace_back(JetTagging::Parameters(JetTagging::DeepCSV, JetTagging::Medium, JetTagging::incl, JetTagging::mujets));
	else
		jtps.emplace_back(JetTagging::Parameters(JetTagging::DeepJet, JetTagging::Medium, JetTagging::incl, JetTagging::mujets));
	mcCorr->SetJetTaggingParameters(jtps);
	
	// initiate regions and cutflows
	regions = {"SR_Base3mu", "SR_Base1e2mu", "CR_DYdimu", "CR_TTdimu", "CR_TTemu", "CR_SSdimu", "CR_SSemu", "CR_Conv3mu", "CR_Conv1e2mu", "CR_RevBjet" };
	for (const auto& region : regions)
		InitiateCutflow(region, this->getCuts(region));
}

vector<TString> BaselineSelector::getCuts(const TString& region) {
	if (region == "SR_Base3mu")
		return {"nocut", "metfilter", "3mu", "trigger", "passSafePtCut", "ExOSmu", "Mge12", "NJge2", "NBge1"};
	else if (region == "SR_Base1e2mu")
		return {"nocut", "metfilter", "1e2mu", "trigger", "passSafePtCut", "ExOSmu", "Mge12", "NJge2", "NBge1"};
	else if (region == "CR_DYdimu")
		return {"nocut", "metfilter", "OSdimu", "trigger", "passSafePtCut", "OnshellZ", "NoB"};
	else if (region == "CR_TTdimu")
		return {"nocut", "metfilter", "OSdimu", "trigger", "passSafePtCut", "OffshellZ", "Mge12", "dRge04", "NJge2", "NBge1"};
	else if (region == "CR_TTemu")
		return {"nocut", "metfilter", "OSemu", "trigger", "passSafePtCut", "dRge04", "NJge2", "NBge1"};
	else {
		cerr << "[BaselineSelector::getCuts] wrong region " << region << endl;
		exit(EXIT_FAILURE);
	}
}

void BaselineSelector::executeEvent(){
	// cutflow
	for (const auto& region: regions)
		FillCutflow(region, "nocut");
	if (!PassMETFilter())
		return;
	for (const auto& region: regions)
		FillCutflow(region, "metfilter");

	ev = GetEvent();
	gens = GetGens();
	muons = GetAllMuons();
	electrons = GetAllElectrons();
	jets = GetAllJets();
	METv = ev.GetMETVector();

	// sort the objects at the first time, don't wand to be confused
	sort(muons.begin(), muons.end(), PtComparing);
	sort(electrons.begin(), electrons.end(), PtComparing);
	sort(jets.begin(), jets.end(), PtComparing);


  AnalyzerParameter param;

  executeEventFromParameter(param);

}

void BaselineSelector::executeEventFromParameter(AnalyzerParameter param){

  if(!PassMETFilter()) return;

  Event ev = GetEvent();

}

