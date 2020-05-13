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

	//==== Flag for systematic sources
	RunSysts = HasFlag("RunSysts");
	RunXsecSyst = HasFlag("RunXsecSyst");

	cout << "[FakeEstimator::initializeAnalyzer] RunSysts = " << RunSysts << endl;
	cout << "[FakeEstimator::initializeAnalyzer] RunXsecSyst = " << RunXsecSyst << endl;

	//==== Systematic sources
	Systs = {"BtagDep", "JetPtCut30", "JetPtCut40", "JetPtCut50", "JetPtCut60"};
	//==== ID setting for electrons
	ElectronIDs = {"passLooseID", "passTightID", "FakeLooseID", "FakeTightID"};
	
	//==== Trigger Setting
	if (DataYear == 2016) {
		HLTElecTriggerName = "HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v";
		TriggerSafePtCut = 25.;
	}
	else if (DataYear == 2017) {
		HLTElecTriggerName = "HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v";
		TriggerSafePtCut = 25.;
	}
	else if (DataYear == 2018) {
		cout << "[FakeEstimator::initializeAnalyzer] Trigger is not set for 2018" << endl;
		exit(EXIT_FAILURE);
	}
	else {
		cout << "[FakeEstimator::initializeAnalyzer] Wrong Year" << endl;
		exit(EXIT_FAILURE);
	}

	cout << "[FakeEstimator::initializeAnalyzer] HLTElecTriggerName = " << HLTElecTriggerName << endl;
	cout << "[FakeEstimator::initializeAnalyzer] TriggerSafePtCut = " << TriggerSafePtCut << endl;

	//==== B-Tagging
	std::vector<JetTagging::Parameters> jtps;
	jtps.push_back( JetTagging::Parameters(JetTagging::DeepCSV, JetTagging::Medium, JetTagging::incl, JetTagging::comb) );
	mcCorr->SetJetTaggingParameters(jtps);
	
	cout << "[FakeEstiamtor::initializeAnalyer] Finish initialization" << endl;

	f_nPV = new TFile("/home/choij/SKFlat/data/Run2Legacy_v4/2016/nPV/nPV_reweight.root");
}

void FakeEstimator::executeEvent(){

	//==== Copy all objects ====
	AllMuons = GetAllMuons();
	AllElectrons = GetAllElectrons();
	AllJets = GetAllJets();
	AllGens = GetGens();

	//==== Get L1Prefire. pileup reweight ====
	weight_Prefire = GetPrefireWeight(0);
	weight_PileUp = GetPileUpWeight(nPileUp, 0);

	AnalyzerParameter param;

	//==== Loop over electron IDs ====
	for (unsigned int i = 0; i < ElectronIDs.size(); i++) {
		MuonID = "POGLoose"; // to veto loose muons
		ElectronID = ElectronIDs.at(i);
		JetID = "tight";

		param.Clear();
		param.syst_ = AnalyzerParameter::Central;
		param.Name = ElectronID + "_Central";
		param.Jet_ID = JetID;
		

		executeEventFromParameter(param);

		//==== Systematic sources ====
		//==== need update for MET sources
		if (RunSysts) {
			for (unsigned int i = 1; i < AnalyzerParameter::NSyst; i++) {
				param.syst_ = AnalyzerParameter::Syst(i);
				param.Name = ElectronID + "_Syst_" + param.GetSystType();
				cout << "[FakeEstimator::executeEvent] execute " << param.Name << endl;
				executeEventFromParameter(param);
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
    if (! (ev.PassTrigger(HLTElecTriggerName) )) return;
	
	//==== Copy all objects ====
	vector<Muon> this_AllMuons = AllMuons;
	vector<Electron> this_AllElectrons = AllElectrons;
	vector<Jet> this_AllJets = AllJets;

	//==== Normalization Systematic Sources ====
	if (param.syst_ == AnalyzerParameter::Central) {
	}
	else if (param.syst_ == AnalyzerParameter::JetResUp) {
		this_AllJets = SmearJets( this_AllJets, +1 );
	}
	else if(param.syst_ == AnalyzerParameter::JetResDown) {
		this_AllJets = SmearJets( this_AllJets, -1 );
	}
	else if(param.syst_ == AnalyzerParameter::JetEnUp) {
		this_AllJets = ScaleJets( this_AllJets, +1 );
	}
	else if(param.syst_ == AnalyzerParameter::JetEnDown) {
		this_AllJets = ScaleJets( this_AllJets, -1 );
	}
	else if(param.syst_ == AnalyzerParameter::ElectronResUp) {
		this_AllElectrons = SmearElectrons( this_AllElectrons, +1 );
	}
	else if(param.syst_ == AnalyzerParameter::ElectronResDown) {
		this_AllElectrons = SmearElectrons( this_AllElectrons, -1 );
	}
	else if(param.syst_ == AnalyzerParameter::ElectronEnUp) {
		this_AllElectrons = ScaleElectrons( this_AllElectrons, +1 );
	}
	else if(param.syst_ == AnalyzerParameter::ElectronEnDown) {
		this_AllElectrons = ScaleElectrons( this_AllElectrons, -1 );
	}

	//==== ID Selection ====
	muons = SelectMuons(this_AllMuons, MuonID, 10., 2.4);
	electrons = SelectElectrons(this_AllElectrons, ElectronID, 25, 2.5);
	jets = SelectJets(this_AllJets, param.Jet_ID, 20., 2.4);

	if (ElectronID.Contains("pass")) electrons_loose = SelectElectrons(this_AllElectrons, "passLooseID", 25, 2.5);
	else if (ElectronID.Contains("Fake")) electrons_loose = SelectElectrons(this_AllElectrons, "FakeLooseID", 25, 2.5);
	
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
	if (electrons.at(0).Pt() <= TriggerSafePtCut) return;

	/////////////////////////////////////////////////////////////////////
	//==== Event Selection
	/////////////////////////////////////////////////////////////////////
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
	///////////////////////////////////////////////////////////////////////
	if (muons.size() != 0) return;
	
	vector<double> jetPtCuts = {0, 30, 40, 50, 60};
	double weight = 1.;

	//==== QCD dominated & W boson dominated region
	if (electrons.size() == 1) {
		Particle elec = Particle(electrons.at(0));
		double Mt = MT(elec, METv);

		//==== QCD dominated window
		if ( Mt < 25 && METv.Pt() < 25) {
			clean_jets10 = JetsVetoLeptonInside(jets, electrons_loose, muons, 1.0);
			double corrPt = GetCorrPt(electrons.at(0));

			if (clean_jets10.size() == 0) return;
			if (clean_jets10.at(0).Pt() < 30) return;

			if (!IsDATA) {
				weight *= weight_norm_1invpb*ev.GetTriggerLumi(HLTElecTriggerName);
				weight *= ev.MCweight();
				weight *= weight_Prefire;
				weight *= GetNPVReweight(ElectronID, "JetPtCut40");
			}

			if (clean_jets10.at(0).Pt() > jetPtCuts.at(1)) {
				FillHist(Systs.at(1), "corrPt_QCD_enriched_" + param.Name , corrPt, weight, 48, 0., 240.);
			}
			if (clean_jets10.at(0).Pt() > jetPtCuts.at(2)) {
				FillHist(Systs.at(2), "corrPt_QCD_enriched_" + param.Name, corrPt, weight, 48, 0., 240.);
				
				//==== B-tagging source
				for (unsigned int i = 0; i < clean_jets10.size(); i++) {
					double this_discr = clean_jets10.at(i).GetTaggerResult(JetTagging::DeepCSV);
					if (this_discr > mcCorr->GetJetTaggingCutValue(JetTagging::DeepCSV, JetTagging::Medium)) NBjets_NoSF++;
					if (mcCorr->IsBTagged_2a(jtp_DeepCSV_Medium, clean_jets10.at(i))) NBjets_WithSF_2a++;
				}

				if (IsDATA && NBjets_NoSF != 0) {
					FillHist(Systs.at(0), "corrPt_QCD_enriched_" + param.Name, corrPt, weight, 48, 0., 240.);
				}
				if (!IsDATA && NBjets_WithSF_2a != 0) {
					FillHist(Systs.at(0), "corrPt_QCD_enriched_" + param.Name, corrPt, weight, 48, 0., 240.);
				}
			}

			if (clean_jets10.at(0).Pt() > jetPtCuts.at(3)) {
				FillHist(Systs.at(3), "corrPt_QCD_enriched_" + param.Name, corrPt, weight, 48, 0., 240.);
			}
			if (clean_jets10.at(0).Pt() > jetPtCuts.at(4)) {
				FillHist(Systs.at(4), "corrPt_QCD_enriched_" + param.Name, corrPt, weight, 48, 0., 240.);
			}
			return;
		}


		//==== W boson dominated window
		else if (Mt > 70 && METv.Pt() >50) {
			clean_jets04 = JetsVetoLeptonInside(jets, electrons_loose, muons, 0.4);

			if (clean_jets04.size() == 0) return;
			if (clean_jets04.at(0).Pt() < 30) return;

			if (!IsDATA) {
				weight *= weight_norm_1invpb*ev.GetTriggerLumi(HLTElecTriggerName);
				weight *= ev.MCweight();
				weight *= weight_Prefire;
				weight *= GetNPVReweight(ElectronID, "JetPtCut40");
			}

			if (clean_jets04.at(0).Pt() > jetPtCuts.at(1)) {
				FillHist(Systs.at(1), "Mt_W_enriched_" + param.Name, Mt, weight, 36, 60., 240);
			}
			if (clean_jets04.at(0).Pt() > jetPtCuts.at(2)) {
				FillHist(Systs.at(2), "Mt_W_enriched_" + param.Name, Mt, weight, 36, 60., 240);
				
				//==== B-tagging source
				for (unsigned int i = 0; i < clean_jets04.size(); i++) {
					double this_discr = clean_jets04.at(i).GetTaggerResult(JetTagging::DeepCSV);
					if (this_discr > mcCorr->GetJetTaggingCutValue(JetTagging::DeepCSV, JetTagging::Medium)) NBjets_NoSF++;
					if (mcCorr->IsBTagged_2a(jtp_DeepCSV_Medium, clean_jets04.at(i))) NBjets_WithSF_2a++;
				}

				if (IsDATA && NBjets_NoSF != 0) {
					FillHist(Systs.at(0), "Mt_W_enriched_" + param.Name, Mt, weight, 36, 60., 240);
				}
				if (!IsDATA && NBjets_WithSF_2a != 0) {
					FillHist(Systs.at(0), "Mt_W_enriched_" + param.Name, Mt, weight, 36, 60., 240);
				}

			}
			if (clean_jets04.at(0).Pt() > jetPtCuts.at(3)) {
				FillHist(Systs.at(3), "Mt_W_enriched_" + param.Name, Mt, weight, 36, 60., 240);
			}
			if (clean_jets04.at(0).Pt() > jetPtCuts.at(4)) {
				FillHist(Systs.at(4), "Mt_W_enriched_" + param.Name, Mt, weight, 36, 60., 240);
			}

			return;
		}
		else return;
	}
	
	//==== Z boson dominated window
	else if (electrons.size() == 2) {
		Particle ZCand = electrons.at(0) + electrons.at(1);
		clean_jets04 = JetsVetoLeptonInside(jets, electrons_loose, muons, 0.4);

		if (clean_jets04.size() == 0) return;
		if (clean_jets04.at(0).Pt() < 30) return;
		if (!IsOnZ(ZCand.M(), 15.)) return;

		if (!IsDATA) {
			weight *= weight_norm_1invpb*ev.GetTriggerLumi(HLTElecTriggerName);
			weight *= ev.MCweight();
			weight *= weight_Prefire;
			weight *= GetNPVReweight(ElectronID, "JetPtCut40");
		}

		if (clean_jets04.at(0).Pt() > jetPtCuts.at(1)) {
			FillHist(Systs.at(1), "M(ee)_Z_enriched_" + param.Name, ZCand.M(), weight, 20, 70., 110.);
		}
		if (clean_jets04.at(0).Pt() > jetPtCuts.at(2)) {
			FillHist(Systs.at(2), "M(ee)_Z_enriched_" + param.Name, ZCand.M(), weight, 20, 70., 110.);
		}
		if (clean_jets04.at(0).Pt() > jetPtCuts.at(3)) {
			FillHist(Systs.at(3), "M(ee)_Z_enriched_" + param.Name, ZCand.M(), weight, 20, 70., 110.);
		}
		if (clean_jets04.at(0).Pt() > jetPtCuts.at(4)) {
			FillHist(Systs.at(4), "M(ee)_Z_enriched_" + param.Name, ZCand.M(), weight, 20, 70., 110.);
		}

		return;
	}
	else return;
}

FakeEstimator::FakeEstimator(){

}

FakeEstimator::~FakeEstimator(){

}

// Private Functions
double FakeEstimator::GetNPVReweight(TString id, TString syst) {
	TDirectory* temp_dir = (TDirectory*)f_nPV->GetDirectory(id + "_Central");
	TH1D* h = (TH1D*)temp_dir->Get("nPV_reweight_" + id + "_" + syst);

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

TH1D* FakeEstimator::GetHist1D(TString suffix, TString histname) {
	TH1D* h = NULL;

	map<TString, map<TString, TH1D*> >::iterator mapit = maphist_TH1D.find(suffix);
	if(mapit == maphist_TH1D.end()) return h;
	else {
		map<TString, TH1D*> this_maphist = mapit->second;
		map<TString, TH1D*>::iterator mapitit = this_maphist.find(histname);
		if (mapitit != this_maphist.end()) return mapitit->second;
	}

	return h;
}

TH2D* FakeEstimator::GetHist2D(TString suffix, TString histname) {
	TH2D* h = NULL;

	map<TString, map<TString, TH2D*> >::iterator mapit = maphist_TH2D.find(suffix);
	if (mapit == maphist_TH2D.end()) return h;
	else {
		map<TString, TH2D*> this_maphist = mapit->second;
		map<TString, TH2D*>::iterator mapitit = this_maphist.find(histname);
		if (mapitit != this_maphist.end()) return mapitit->second;
	}

	return h;
}

void FakeEstimator::FillHist(TString syst, TString histname, double value, double weight, int n_bin, double x_min, double x_max) {
	TString id;
	for (unsigned int i = 0; i < ElectronIDs.size(); i++) {
		if (histname.Contains(ElectronIDs.at(i))) id = ElectronIDs.at(i);
	}

	TString suffix = id + "_" + syst;
	histname = histname + "_" + syst;
	TH1D* this_hist = GetHist1D(suffix, histname);
	if (!this_hist) {
		this_hist = new TH1D(histname, "", n_bin, x_min, x_max);
		maphist_TH1D[suffix][histname] = this_hist;
	}

	this_hist->Fill(value, weight);
}

void FakeEstimator::FillHist(TString syst, TString histname, double value_x, double value_y, double weight,
			                                                double n_binx, double* xbins,
															double n_biny, double* ybins) {
	TString id;
	for (unsigned int i = 0; i < ElectronIDs.size(); i++) {
		if (histname.Contains(ElectronIDs.at(i))) id = ElectronIDs.at(i);
	}

	TString suffix = id + "_" + syst;
	histname = histname + "_" + syst;
	TH2D* this_hist = GetHist2D(suffix, histname);
	if (!this_hist) {
		this_hist = new TH2D(histname, "", n_binx, xbins, n_biny, ybins);
		maphist_TH2D[suffix][histname] = this_hist;
	}
	this_hist->Fill(value_x, value_y, weight);
}

void FakeEstimator::WriteHist() {
	outfile->cd();

	for (unsigned int i = 0; i < ElectronIDs.size(); i++) {
		outfile->cd();
		outfile->mkdir(ElectronIDs.at(i));
		auto* dirID = (TDirectory*)outfile->GetDirectory(ElectronIDs.at(i));
		for (unsigned int j = 0; j < Systs.size(); j++) {
			dirID->mkdir(Systs.at(j));
			auto* dirSyst = (TDirectory*)dirID->GetDirectory(Systs.at(j));
			mapDirectory[ElectronIDs.at(i)][Systs.at(j)] = dirSyst;
		}
	}

	for (auto it = maphist_TH1D.begin(); it != maphist_TH1D.end(); it++) {
		TString id; TString syst;
		for (unsigned int i = 0; i < ElectronIDs.size(); i++) {
			if (it->first.Contains(ElectronIDs.at(i))) id = ElectronIDs.at(i);
		}
		for (unsigned int i = 0; i < Systs.size(); i++) {
			if (it->first.Contains(Systs.at(i))) syst = Systs.at(i);
		}

		mapDirectory[id][syst]->cd();

		auto maphist = it->second;
		for (auto itit = maphist.begin(); itit != maphist.end(); itit++) {
			itit->second->Write();
		}
	}

	for (auto it = maphist_TH2D.begin(); it != maphist_TH2D.end(); it++) {
		TString id; TString syst;
		for (unsigned int i = 0; i < ElectronIDs.size(); i++) {
			if (it->first.Contains(ElectronIDs.at(i))) id = ElectronIDs.at(i);
		}
		for (unsigned int i = 0; i < Systs.size(); i++) {
			if (it->first.Contains(Systs.at(i))) syst = Systs.at(i);
		}

		mapDirectory[id][syst]->cd();

		auto maphist = it->second;
		for (auto itit = maphist.begin(); itit != maphist.end(); itit++) {
			itit->second->Write();
		}
	}
}

