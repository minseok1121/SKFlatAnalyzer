#include "FakeEstimator.h"

/////////////////////////////////////////////////////////////////////////////
//==== FakeEstimator                                                   ====//
//==== This analyzer is for evaluation of the fake rate of electorns   ====//
//==== from jet activity                                               ====// 
//==== using data                                                      ====//
//==== Contact: choij@cern.ch                                          ====//
//==== Based on SKFlatAnalyzer - Run2Legacy_v4                         ====// 
/////////////////////////////////////////////////////////////////////////////

// FakeEstiamtor: Central estimation of Fake Rate

void FakeEstimator::initializeAnalyzer(){

	//==== Flags for systematic sources
	RunSysts = HasFlag("RunSysts");
	RunXsecSyst = HasFlag("RunXsecSyst");
    cout << "[FakeEstimator::initializeAnalyzer] RunSysts = " << RunSysts << endl;
    cout << "[FakeEstimator::initializeAnalyzer] RunXsecSyst = " << RunXsecSyst << endl;

	//==== informations for histnames & systematic sources
	IDs = {"passLooseID", "passTightID", "FakeLooseID", "FakeTightID"};
	Systs = {"Central", "JetPtCut30", "JetPtCut50", "JetPtCut60", "HadFlavor"};
	Prompts = {"Central", "JetResUp", "JetResDown", "JetEnUp", "JetEnDown",
		"ElectronResUp", "ElectronResDown", "ElectronEnUp", "ElectronEnDown", "PileUp"};
	Regions = {"QCDEnriched", "WEnriched", "ZEnriched"};

	//==== Electron ID setting
	ElectronIDs = IDs;

	//==== TriggerSetting
	if (DataYear == 2016) {
		HLTElecTriggerName1 = "HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v";
		HLTElecTriggerName2 = "HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v";
		TriggerSafePtCut1 = 15.;
		TriggerSafePtCut2 = 25.;
	}
	else if (DataYear == 2017) {
		HLTElecTriggerName = "HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v";
		TriggerSafePtCut = 15.;
	}
	else if (DataYear == 2018) {
		cout << "[FakeEstimator::initializeAnalyzer] Trigger is not set for 2018" << endl;
		exit(EXIT_FAILURE);
	}
	else {
		cout << "[FakeEstimator::initializeAnalyzer] Wrong Year" << endl;
		exit(EXIT_FAILURE);
	}

    cout << "[FakeEstimator::initializeAnalyzer] HLTElecTriggerName1 = " << HLTElecTriggerName1 << endl;
	cout << "[FakeEstimator::initializeAnalyzer] HLTElecTriggerName2 = " << HLTElecTriggerName2 << endl;
    cout << "[FakeEstimator::initializeAnalyzer] TriggerSafePtCut = " << TriggerSafePtCut << endl;

	//==== B-tagging
	std::vector<JetTagging::Parameters> jtps;
	jtps.push_back( JetTagging::Parameters(JetTagging::DeepCSV, JetTagging::Medium, JetTagging::incl, JetTagging::comb) );
    mcCorr->SetJetTaggingParameters(jtps);

    cout << "[FakeEstiamtor::initializeAnalyer] Finish initialization" << endl;

	if (DataYear == 2016) {
		f_nPV = new TFile("/home/choij/SKFlat/data/Run2Legacy_v4/2016/nPV/nPV_reweight.root");
	}
	else if (DataYear == 2017 || DataYear == 2018) {
		cout << "[FakeEstimator::initializeAnalyzer] nPV_reweight is not set yet" << endl;
		exit(EXIT_FAILURE);
	}
	else {
		cout << "[FakeEstiamtor::initializeAnalyzer] Wrong Year" << endl;
		exit(EXIT_FAILURE);
	}
}

void FakeEstimator::executeEvent(){

	//==== Copy all objects ====
	AllMuons = GetAllMuons();
	AllElectrons = GetAllElectrons();
	AllJets = GetAllJets();

	//==== Get L1Prefire. pileup reweight ====
	weight_Prefire = GetPrefireWeight(0);
	weight_PileUp = GetPileUpWeight(nPileUp, 0);
  
	AnalyzerParameter param;

	//==== Loop over electron IDs, Systs, prompts ====
	for (unsigned int i = 0; i < ElectronIDs.size(); i++) {
		MuonID = "POGLoose"; // to veto loose muons
		ElectronID = ElectronIDs.at(i);
		JetID = "tight";
		JetPtCut = FakeEstimator::GetJetPtCut("Central");

		if (!RunSysts) {
			syst = Systs.at(0);
			prompt = Prompts.at(0);

			param.Clear();
			param.syst_ = AnalyzerParameter::Central;
			param.Name = ElectronID + "_" + syst + "_" + prompt;
			//cout << "[FakeEstimator::executeEvent] execute " << param.Name << endl;
			param.Jet_ID = JetID;

			executeEventFromParameter(param);
		}
		else {	//run over systematic sources...turn on RunSysts
			for (unsigned int j = 0; j < Systs.size(); j++) {
				JetPtCut = GetJetPtCut(Systs.at(j));
				syst = Systs.at(j);

				for (unsigned int k = 0; k < AnalyzerParameter::NSyst; k++) {
					prompt = Prompts.at(k);
					
					param.Clear();
					param.syst_ = AnalyzerParameter::Syst(k);
					param.Name = ElectronID + "_" + syst + "_" + prompt;
					//cout << "[FakeEstimator::executeEvent] execute " << param.Name << endl;
					param.Jet_ID = JetID;
					executeEventFromParameter(param);
				}
			}
		}

		//==== Xsec syst ====
		if (RunXsecSyst) {
			cout << "[FakeEstimator::executeEvent] Xsec systematics is not set yet" << endl;
			exit(EXIT_FAILURE);
		}
	}
}

void FakeEstimator::executeEventFromParameter(AnalyzerParameter param){

	if(!PassMETFilter()) return;

	Event ev = GetEvent();
	Particle METv = ev.GetMETVector();

	//==== Trigger ====
	if (! (ev.PassTrigger(HLTElecTriggerName1)) && !(ev.PassTrigger(HLTElecTriggerName2 ))) return;

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

	//==== ID Selection ====
	muons = SelectMuons(this_AllMuons, MuonID, 10., 2.4);
	electrons = SelectElectrons(this_AllElectrons, ElectronID, 15, 2.5);
	jets = SelectJets(this_AllJets, param.Jet_ID, 20., 2.4);

	if (ElectronID.Contains("pass")) electrons_loose = SelectElectrons(this_AllElectrons, "passLooseID", 15, 2.5);
	else if (ElectronID.Contains("Fake")) electrons_loose = SelectElectrons(this_AllElectrons, "FakeLooseID", 15, 2.5);

	std::sort(muons.begin(), muons.end(), PtComparing);
	std::sort(electrons.begin(), electrons.end(), PtComparing);
	std::sort(electrons_loose.begin(), electrons_loose.end(), PtComparing);
	std::sort(jets.begin(), jets.end(), PtComparing);

	int NBjets_NoSF = 0;
	int NBjets_WithSF_2a = 0;
	JetTagging::Parameters jtp_DeepCSV_Medium = JetTagging::Parameters(JetTagging::DeepCSV, JetTagging::Medium, JetTagging::incl, JetTagging::comb);
	//double btagWeight = mcCorr->GetBTaggingReweight_1a(jets, jtp_DeepCSV_Medium);
	
	//==== leading electron should pass trigger-safe pt cut
	if (electrons.size() == 0) return;
	if (electrons.at(0).Pt() <= TriggerSafePtCut1) return;

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
	if (muons.size() != 0) return;
	double weight = 1.;

	//==== QCD dominated & W boson dominated region ====
	if (electrons.size() == 1) {
		Particle elec = Particle(electrons.at(0));
		double Mt = MT(elec, METv);
		clean_jets04 = JetsVetoLeptonInside(jets, electrons_loose, muons, 0.4);
		if (clean_jets04.size() == 0) return;
		if (clean_jets04.at(0).Pt() < JetPtCut) return;

		double corrPt = GetCorrPt(electrons.at(0));
		double elecEta = fabs(electrons.at(0).Eta());
		double ptbins[7] =  {15., 20., 30., 50., 70., 110., 200};
		double etabins[4] = {0., 0.8, 1.479, 2.5};

		if (15. < corrPt && corrPt < 30.) {
			if (!ev.PassTrigger(HLTElecTriggerName1)) return;
			if (!IsDATA) weight *= weight_norm_1invpb * ev.GetTriggerLumi(HLTElecTriggerName1);
		}
		else if (corrPt > 30.) {
			if (!ev.PassTrigger(HLTElecTriggerName2)) return;
			if (!IsDATA) weight *= weight_norm_1invpb * ev.GetTriggerLumi(HLTElecTriggerName2);
		}
		else return;

		//==== W boson dominated window ====
		if ( Mt > 70. && METv.Pt() > 50. ) {
			if (!IsDATA) {
				weight *= ev.MCweight();
				weight *= weight_Prefire;
				if (param.syst_ == AnalyzerParameter::PileUp) weight *= GetPileUpWeight(nPileUp, 0);
				else weight *= GetNPVReweight(ElectronID, "JetPtCut40");
			}
			param.Name = Regions.at(1) + "_" + param.Name;
			
			if (syst == "HadFlavor") {
				//==== B-tagging source ====
				for (unsigned int i = 0; i < clean_jets04.size(); i++) {
					double this_discr = clean_jets04.at(i).GetTaggerResult(JetTagging::DeepCSV);
					if (this_discr > mcCorr->GetJetTaggingCutValue(JetTagging::DeepCSV, JetTagging::Medium)) NBjets_NoSF++;
					if (mcCorr->IsBTagged_2a(jtp_DeepCSV_Medium, clean_jets04.at(i))) NBjets_WithSF_2a++;
				}

				if (IsData && NBjets_NoSF == 0) return;
				if (!IsData && NBjets_WithSF_2a == 0 ) return;

				FillHist("PromptEvents_" + param.Name, 1, weight, 2, 0., 2.);
				FillHist("Mt_" + param.Name, Mt, weight, 36, 60., 240.);
				FillHist("ElectronPt_" + param.Name, electrons.at(0).Pt(), weight, 48, 0., 240.);
				FillHist("corrPt_" +param.Name, corrPt, weight, 48, 0., 240.);
				FillHist("ElectronEta_" + param.Name, electrons.at(0).Eta(), weight, 25, -2.5, 2.5);
				FillHist("LeadingJetPt_" + param.Name, clean_jets04.at(0).Pt(), weight, 48, 0., 240.);
				FillHist("LeadingJetEta_" + param.Name, clean_jets04.at(0).Eta(), weight, 24, -2.4, 2.4);
				FillHist("MET_" + param.Name, METv.Pt(), weight, 48, 0., 240.);
				FillHist("MetPhi_" + param.Name, METv.Phi(), weight, 32, -4, 4);
				FillHist("passID_" + param.Name, corrPt, elecEta, weight, 6, ptbins, 3, etabins);
				FillHist("LeadingJetPhi_" + param.Name, clean_jets04.at(0).Phi(), weight, 32, -4, 4);
				FillHist("ElectronPhi_" + param.Name, electrons.at(0).Phi(), weight, 32, -4, 4);
			}

			else {
                FillHist("PromptEvents_" + param.Name, 1, weight, 2, 0., 2.);
				FillHist("Mt_" + param.Name, Mt, weight, 36, 60., 240.);
				FillHist("ElectronPt_" + param.Name, electrons.at(0).Pt(), weight, 48, 0., 240.);
				FillHist("corrPt_" +param.Name, corrPt, weight, 48, 0., 240.);
				FillHist("ElectronEta_" + param.Name, electrons.at(0).Eta(), weight, 25, -2.5, 2.5);
				FillHist("LeadingJetPt_" + param.Name, clean_jets04.at(0).Pt(), weight, 48, 0., 240.);
                FillHist("LeadingJetEta_" + param.Name, clean_jets04.at(0).Eta(), weight, 24, -2.4, 2.4);
                FillHist("MET_" + param.Name, METv.Pt(), weight, 48, 0., 240.);
                FillHist("MetPhi_" + param.Name, METv.Phi(), weight, 32, -4, 4);
                FillHist("passID_" + param.Name, corrPt, elecEta, weight, 6, ptbins, 3, etabins);
				FillHist("LeadingJetPhi_" + param.Name, clean_jets04.at(0).Phi(), weight, 32, -4, 4);
				FillHist("ElectronPhi_" + param.Name, electrons.at(0).Phi(), weight, 32, -4, 4);
			}
			return;
		}

		//==== QCD dominated window ====
		else if ( Mt < 25. && METv.Pt() < 25) {
			clean_jets10 = JetsVetoLeptonInside(jets, electrons_loose, muons, 1.0);

			if (clean_jets10.size() == 0) return;
			if (clean_jets10.at(0).Pt() < JetPtCut) return;

            if (!IsDATA) {
                weight *= ev.MCweight();
                weight *= weight_Prefire;
                if (param.syst_ == AnalyzerParameter::PileUp) weight *= GetPileUpWeight(nPileUp, 0);
                else weight *= GetNPVReweight(ElectronID, "JetPtCut40");
			}
			param.Name = Regions.at(0) + "_" + param.Name;

			if (syst == "HadFlavor") {
				//==== B-tagging source ====
				for (unsigned int i = 0; i < clean_jets10.size(); i++) {
                    double this_discr = clean_jets10.at(i).GetTaggerResult(JetTagging::DeepCSV);
                    if (this_discr > mcCorr->GetJetTaggingCutValue(JetTagging::DeepCSV, JetTagging::Medium)) NBjets_NoSF++;
					if (mcCorr->IsBTagged_2a(jtp_DeepCSV_Medium, clean_jets10.at(i))) NBjets_WithSF_2a++;
				}

				if (IsDATA && NBjets_NoSF == 0) return;
				if (!IsDATA && NBjets_WithSF_2a == 0) return;

				FillHist("corrPt_" + param.Name, corrPt, weight, 48, 0., 240.);
				FillHist("ElectronEta_" + param.Name, electrons.at(0).Eta(), weight, 25, -2.5, 2.5);
				FillHist("LeadingJetPt_" + param.Name, clean_jets10.at(0).Pt(), weight, 48, 0., 240.);
				FillHist("LeadingJetEta_" + param.Name, clean_jets10.at(0).Eta(), weight, 24, -2.4, 2.4);
				FillHist("passID_" + param.Name, corrPt, elecEta, weight, 6, ptbins, 3, etabins);
			}
			else {
                FillHist("corrPt_" + param.Name, corrPt, weight, 48, 0., 240.);
                FillHist("ElectronEta_" + param.Name, electrons.at(0).Eta(), weight, 25, -2.5, 2.5);
                FillHist("LeadingJetPt_" + param.Name, clean_jets10.at(0).Pt(), weight, 48, 0., 240.);
                FillHist("LeadingJetEta_" + param.Name, clean_jets10.at(0).Eta(), weight, 24, -2.4, 2.4);
                FillHist("passID_" + param.Name, corrPt, elecEta, weight, 6, ptbins, 3, etabins);
			}
			return;
		}
		else return;
	}
	//==== Z boson dominated window =====
	else if (electrons.size() == 2) {
		double corrPt = GetCorrPt(electrons.at(0));
		if (15. < corrPt && corrPt < 30.) {
            if (!ev.PassTrigger(HLTElecTriggerName1)) return;
			if (!IsDATA) weight *= weight_norm_1invpb*ev.GetTriggerLumi(HLTElecTriggerName1);
        }
        else if (corrPt > 30.) {
            if (!ev.PassTrigger(HLTElecTriggerName2)) return;
			if (!IsDATA) weight *= weight_norm_1invpb*ev.GetTriggerLumi(HLTElecTriggerName2);
        }
        else return;

		Particle ZCand = electrons.at(0) + electrons.at(1);
		clean_jets04 = JetsVetoLeptonInside(jets, electrons_loose, muons, 0.4);
		if (clean_jets04.size() == 0) return;
		if (clean_jets04.at(0).Pt() < JetPtCut) return;
		if (!IsOnZ(ZCand.M(), 15.)) return;

        if (!IsDATA) {
			weight *= ev.MCweight();
            weight *= weight_Prefire;
			if (param.syst_ == AnalyzerParameter::PileUp) weight *= GetPileUpWeight(nPileUp, 0);
            else weight *= GetNPVReweight(ElectronID, "JetPtCut40");
		}
		param.Name = Regions.at(2) + "_" + param.Name;

		if (syst == "HadFlavor") {
			//==== B-tagging source ====
            for (unsigned int i = 0; i < clean_jets04.size(); i++) {
                double this_discr = clean_jets04.at(i).GetTaggerResult(JetTagging::DeepCSV);
                if (this_discr > mcCorr->GetJetTaggingCutValue(JetTagging::DeepCSV, JetTagging::Medium)) NBjets_NoSF++;
                if (mcCorr->IsBTagged_2a(jtp_DeepCSV_Medium, clean_jets04.at(i))) NBjets_WithSF_2a++;
			}

			if (IsDATA && NBjets_NoSF == 0) return;
			if (!IsDATA && NBjets_WithSF_2a == 0) return;

			FillHist("M(ee)_" + param.Name, ZCand.M(), weight, 20, 70., 110.);
			FillHist("LeadingElectronPt_" + param.Name, electrons.at(0).Pt(), weight, 48, 0., 240.);
			FillHist("LeadingElectronEta_" + param.Name, electrons.at(0).Eta(), weight, 25, -2.5, 2.5);
			FillHist("SubLeadingElectronPt_" + param.Name, electrons.at(1).Pt(), weight, 48, 0., 240.);
			FillHist("SubLeadingElectronEta_" + param.Name, electrons.at(1).Eta(), weight, 25, -2.5, 2.5);
			FillHist("LeadingJetPt_" + param.Name, clean_jets04.at(0).Pt(), weight, 48, 0., 240);
			FillHist("LeadingJetEta_" + param.Name, clean_jets04.at(0).Eta(), weight, 24, -2.4, 2.4);
			FillHist("LeadingJetPhi_" + param.Name, clean_jets04.at(0).Phi(), weight, 32, -4, 4);
			FillHist("LeadingElectronPhi_" + param.Name, electrons.at(0).Phi(), weight, 32, -4, 4);
			FillHist("SubLeadingElectronPhi_" + param.Name, electrons.at(1).Phi(), weight, 32, -4, 4);
		}
		else {
            FillHist("M(ee)_" + param.Name, ZCand.M(), weight, 20, 70., 110.);
            FillHist("LeadingElectronPt_" + param.Name, electrons.at(0).Pt(), weight, 48, 0., 240.);
            FillHist("LeadingElectronEta_" + param.Name, electrons.at(0).Eta(), weight, 25, -2.5, 2.5);
            FillHist("SubLeadingElectronPt_" + param.Name, electrons.at(1).Pt(), weight, 48, 0., 240.);
            FillHist("SubLeadingElectronEta_" + param.Name, electrons.at(1).Eta(), weight, 25, -2.5, 2.5);
            FillHist("LeadingJetPt_" + param.Name, clean_jets04.at(0).Pt(), weight, 48, 0., 240);
            FillHist("LeadingJetEta_" + param.Name, clean_jets04.at(0).Eta(), weight, 24, -2.4, 2.4);
			FillHist("LeadingJetPhi_" + param.Name, clean_jets04.at(0).Phi(), weight, 32, -4, 4);
			FillHist("LeadingElectronPhi_" + param.Name, electrons.at(0).Phi(), weight, 32, -4, 4);
			FillHist("SubLeadingElectronPhi_" + param.Name, electrons.at(1).Phi(), weight, 32, -4, 4);
		}
		return;
	}

	else return;
}

FakeEstimator::FakeEstimator(){

}

FakeEstimator::~FakeEstimator(){

}


//==== member fuctions
double FakeEstimator::GetJetPtCut(TString syst) {
	double cut;
	if (syst == Systs.at(0) || syst == Systs.at(4) ) cut = 40;
	else if (syst == Systs.at(1)) cut = 30;
	else if (syst == Systs.at(2)) cut = 50;
	else if (syst == Systs.at(3)) cut = 60;
	else {
		cout << "[FakeEstimator::GetJetPtCUt] Wrong Syst" << endl;
		exit(EXIT_FAILURE);
	}

	return cut;
}

double FakeEstimator::GetNPVReweight(TString id, TString syst) {
	TDirectory* temp_dir = (TDirectory*)f_nPV->GetDirectory(id + "_Central");
	TH1D* h = (TH1D*)temp_dir->Get("nPV_reweight_" + id + "_" + syst);
	if (!h) {
		cout << "[FakeEstimator::GetNPVReweight] No Such histogram" << endl;
		exit(EXIT_FAILURE);
	}

	if (nPV > 100) nPV = 100;
	int this_bin = nPV;

	return h->GetBinContent(this_bin);
}

double FakeEstimator::GetCorrPt(Electron e) {
	double corrPt;
	double weight;
	double relIsoTight = 0.06;

	if (e.RelIso() > relIsoTight) {
		weight = 1. + e.RelIso() - relIsoTight;
		corrPt = e.Pt() * weight;
	}
	else corrPt = e.Pt();

	return corrPt;
}

TH1D* FakeEstimator::GetHist1D(TString histname) {
	TH1D* h = NULL;

	map<TString, TH1D*>::iterator mapit = maphist_TH1D.find(histname);
	
	if (mapit == maphist_TH1D.end()) return h;
	else return mapit->second;
}

TH2D* FakeEstimator::GetHist2D(TString histname) {
	TH2D* h = NULL;
	map<TString, TH2D*>::iterator mapit = maphist_TH2D.find(histname);
	
	if (mapit == maphist_TH2D.end()) return h;
	else return mapit->second;
}

void FakeEstimator::FillHist(TString histname, double value, double weight, int n_bin, double x_min, double x_max) {
	TH1D* this_hist = GetHist1D(histname);
	if (!this_hist) {
		this_hist = new TH1D(histname, "", n_bin, x_min, x_max);
		maphist_TH1D[histname] = this_hist;
	}

	this_hist->Fill(value, weight);
}

void FakeEstimator::FillHist(TString histname, double value_x, double value_y, double weight,
													double n_binx, double* xbins,
													double n_biny, double* ybins) {
	TH2D* this_hist= GetHist2D(histname);
	if (!this_hist) {
		this_hist = new TH2D(histname, "" , n_binx, xbins, n_biny, ybins);
		maphist_TH2D[histname] = this_hist;
	}

	this_hist->Fill(value_x, value_y, weight);
}

void FakeEstimator::WriteHist() {
	//==== make directories ====
	for (unsigned int i = 0; i < IDs.size(); i++) {
		outfile->cd();
		outfile->mkdir(IDs.at(i));
		auto* dirID = (TDirectory*)outfile->Get(IDs.at(i));
		for (unsigned int j = 0; j < Systs.size(); j++) {
			dirID->mkdir(Systs.at(j));
			auto* dirSyst = (TDirectory*)dirID->Get(Systs.at(j));
			for (unsigned int k = 0; k < Prompts.size(); k++) {
				dirSyst->mkdir(Prompts.at(k));
				auto* dirPrompt = (TDirectory*)dirSyst->Get(Prompts.at(k));
				for (unsigned int l = 0; l < Regions.size(); l++) {
					dirPrompt->mkdir(Regions.at(l));
					auto* dirRegion = (TDirectory*)dirPrompt->Get(Regions.at(l));

					mapDirectory[IDs.at(i)][Systs.at(j)][Prompts.at(k)][Regions.at(l)] = dirRegion;
				}
			}
		}
	}

	//==== arrange histograms to directories =====
	for (auto it = maphist_TH1D.begin(); it != maphist_TH1D.end(); it++) {
		TString id, syst, prompt, region;
		for (unsigned int i = 0; i < IDs.size(); i++) {
			if (it->first.Contains(IDs.at(i))) id = IDs.at(i);
		}
		for (unsigned int i = 0; i < Systs.size(); i++) {
			if (it->first.Contains(Systs.at(i))) syst = Systs.at(i);
		}
		for (unsigned int i = 0; i < Prompts.size(); i++) {
	        if (it->first.Contains(Prompts.at(i))) prompt = Prompts.at(i);
		}
		for (unsigned int i = 0; i < Regions.size(); i++) {
	        if (it->first.Contains(Regions.at(i))) region = Regions.at(i);
		}

		mapDirectory[id][syst][prompt][region]->cd();
		it->second->Write();
	}

	for (auto it = maphist_TH2D.begin(); it != maphist_TH2D.end(); it++) {
		TString id, syst, prompt, region;
		for (unsigned int i = 0; i < IDs.size(); i++) {
	        if (it->first.Contains(IDs.at(i))) id = IDs.at(i);
		}
		for (unsigned int i = 0; i < Systs.size(); i++) {
			if (it->first.Contains(Systs.at(i))) syst = Systs.at(i);
		}
		for (unsigned int i = 0; i < Prompts.size(); i++) {
			if (it->first.Contains(Prompts.at(i))) prompt = Prompts.at(i);
		}
		for (unsigned int i = 0; i < Regions.size(); i++) {
			if (it->first.Contains(Regions.at(i))) region = Regions.at(i);
		}

		mapDirectory[id][syst][prompt][region]->cd();
		it->second->Write();
	}
}
