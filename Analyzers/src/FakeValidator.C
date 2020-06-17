#include "FakeValidator.h"

void FakeValidator::initializeAnalyzer(){
	//==== No systemtaic setup yet =====
	RunFakeSyst = HasFlag("RunFakeSyst");
	RunPrompt = HasFlag("RunPrompt");

	cout << "[FakeValidator::initializeAnalyzer] RunFakeSyst = " << RunFakeSyst << endl;
	cout << "[FakeValidator::initializeAnalyzer] RunPrompt = " << RunPrompt << endl;
	
	//==== Electron ID setting ====
	ElectronIDs = {"passLooseID", "passTightID", "FakeLooseID", "FakeTightID"};
	idsets = {"POG", "Fake"};

	//==== Trigger setting
	if (DataYear == 2016) {
		HLTElecTriggerName = "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v";
		TriggerSafePtCut1 = 25.;
		TriggerSafePtCut2 = 15.;
	}
	else if (DataYear == 2017) {
		cout << "[FakeValidator::initializeAnalyzer] Trigger is not set for 2017" << endl;
		exit(EXIT_FAILURE);
	}
	else if (DataYear == 2018) {
		cout << "[FakeValidator::initializeAnalyzer] Trigger is not set for 2018" << endl;
		exit(EXIT_FAILURE);
	}
	else {
		cout << "[FakeValidator::initializeAnalyzer] Wrong Year" << endl;
		exit(EXIT_FAILURE);
	}
	
	cout << "[FakeValidator::initializeAnalyzer] HLTElecTriggerName  = " << HLTElecTriggerName << endl;
	cout << "[FakeValidator::initializeAnalyzer] TriggerSafePtCut1 = " << TriggerSafePtCut1 << endl;
	cout << "[FakeValidator::initializeAnalyzer] TriggerSafePtCut2 = " << TriggerSafePtCut2 << endl;

	//==== No b-tagging ====
	
}

void FakeValidator::executeEvent(){

	//==== Copy all objects ====
	AllMuons = GetAllMuons();
	AllElectrons = GetAllElectrons();
	AllJets = GetAllJets();

	//==== Get L1Prefire, pileup reweight ====
	weight_Prefire = GetPrefireWeight(0);
	weight_PileUp = GetPileUpWeight(nPileUp, 0);

	AnalyzerParameter param;

	//==== Loop over ElectornID sets =====
	for (unsigned int i = 0; i < idsets.size(); i++) {
		MuonID = "POGLoose";
		if (idsets.at(i).Contains("POG")) {
			ElectronID = ElectronIDs.at(i);
			ElectronTightID = ElectronIDs.at(i+1);
		}
		else if (idsets.at(i).Contains("Fake")) {
			ElectronID = ElectronIDs.at(2*i);
			ElectronTightID = ElectronIDs.at(2*i + 1);
		}
		JetID = "tight";
	
		param.Clear();
		param.syst_ = AnalyzerParameter::Central;
		param.Name = idsets.at(i) + "ID_Central";
		param.Jet_ID = JetID;

		executeEventFromParameter(param);

		if (RunFakeSyst) {
			param.Name = idsets.at(i) + "ID_Fake_Central";
			executeEventFromParameter(param);

			param.Name = idsets.at(i) + "ID_Fake_SystUp";
			executeEventFromParameter(param);
		
			param.Name = idsets.at(i) + "ID_Fake_SystDown";
			executeEventFromParameter(param);
		}
		
		if (RunFakeSyst && RunPrompt) {
			param.Name = idsets.at(i) + "ID_Fake_Central_WithPrompt";
			executeEventFromParameter(param);

			param.Name = idsets.at(i) + "ID_Fake_SystUp_WithPrompt";
			executeEventFromParameter(param);

			param.Name = idsets.at(i) + "ID_Fake_SystDown_WithPrompt";
			executeEventFromParameter(param);
		}
	}

}

void FakeValidator::executeEventFromParameter(AnalyzerParameter param){

	//==== No cut ====
	JSFillHist(param.Name, "NoCut_" + param.Name, 0., 1., 1, 0., 1.);

	if(!PassMETFilter()) return;

	Event ev = GetEvent();
	Particle METv = ev.GetMETVector();

	//==== Trigger ====
	if (! ev.PassTrigger(HLTElecTriggerName)) return;

	//==== Copy all objects ====
	vector<Muon> this_AllMuons = AllMuons;
	vector<Electron> this_AllElectrons = AllElectrons;
	vector<Jet> this_AllJets = AllJets;

	//==== ID selection ====
	muons = SelectMuons(this_AllMuons, MuonID, 10, 2.4);
	electrons = SelectElectrons(this_AllElectrons, ElectronID, 15, 2.5);
	jets = SelectJets(this_AllJets, param.Jet_ID, 20., 2.4);
	
	std::sort(muons.begin(), muons.end(), PtComparing);
	std::sort(electrons.begin(), electrons.end(), PtComparing);
	std::sort(jets.begin(), jets.end(), PtComparing);

	//==== Event Selectrion ====
	//==== 1. 3 loose electrons, no loose muon ====
	//==== 2. pt(1) > 25, pt(2) > 15, pt(3) > 15 ====
	//==== 3. At least 1 opposite charge electron pair with on-Z
	if (muons.size() != 0) return;
	if (electrons.size() != 3) return;
	if (electrons.at(0).Pt() < TriggerSafePtCut1) return;
	if (electrons.at(1).Pt() < TriggerSafePtCut2) return;
	if (electrons.at(2).Pt() < TriggerSafePtCut2) return;

	double charge[3];
	Particle pair[3];
	bool isOnZ[3];
	for (int i = 0; i < 3; i++) {
		isOnZ[i] = false;
		charge[i] = electrons.at(i%3).Charge() * electrons.at((i+1)%3).Charge();
		pair[i] = electrons.at(i%3) + electrons.at((i+1)%3);

		if (charge[i] < 0 && IsOnZ(pair[i].M(), 10.)) isOnZ[i] = true;
	}

	if (!isOnZ[0] && !isOnZ[1] && !isOnZ[2]) return;

	double weight = 1.;
	if (!IsDATA) {
		weight *= weight_norm_1invpb*ev.GetTriggerLumi("Full");
		weight *= ev.MCweight();
		weight *= weight_Prefire;
		weight *= weight_PileUp;
		
		//==== ID and trigger SF ====
		//==== only for POG ID ====
		if (ElectronID.Contains("pass")) {
			//==== Trigger SF ====
			double eff_DATA = 1., eff_MC = 1.;
			for (unsigned int i = 0; i < electrons.size()-1; i++) {
				if (electrons.at(i).PassID(ElectronTightID)) {
					eff_DATA *= (1. - mcCorr->ElectronTrigger_Eff(ElectronTightID, HLTElecTriggerName, 0, electrons.at(i).scEta(), electrons.at(i).Pt(), 0));
					eff_MC *= (1. - mcCorr->ElectronTrigger_Eff(ElectronTightID, HLTElecTriggerName, 1, electrons.at(i).scEta(), electrons.at(i).Pt(), 0));
				}
				else {
					eff_DATA *= (1. - mcCorr->ElectronTrigger_Eff(ElectronID, HLTElecTriggerName, 0, electrons.at(i).scEta(), electrons.at(i).Pt(), 0));
                    eff_MC *= (1. - mcCorr->ElectronTrigger_Eff(ElectronID, HLTElecTriggerName, 1, electrons.at(i).scEta(), electrons.at(i).Pt(), 0));
				}
			}
			eff_DATA = 1. - eff_DATA;
			eff_MC = 1.- eff_MC;

			weight *= eff_DATA/eff_MC;
			
			//==== ID SF ====
			for (unsigned int i = 0; i < electrons.size(); i++) {
				double this_ElecIDSF;
				if (electrons.at(i).PassID(ElectronTightID)) {
					this_ElecIDSF = mcCorr->ElectronID_SF(ElectronTightID, electrons.at(i).scEta(), electrons.at(i).Pt(), 0);
				}
				else this_ElecIDSF = mcCorr->ElectronID_SF(ElectronID, electrons.at(i).scEta(), electrons.at(i).Pt(), 0);
				weight *= this_ElecIDSF;
			}
		}
	}

	bool tightFlag = (electrons.at(0).PassID(ElectronTightID) && electrons.at(1).PassID(ElectronTightID) && electrons.at(2).PassID(ElectronTightID));

	//==== Fill hist ====
	if (param.Name.Contains("Central") && !param.Name.Contains("_Fake") && tightFlag) {
		JSFillHist(param.Name, "1st_electron_pt_" + param.Name, electrons.at(0).Pt(), weight, 20, 0., 200.);
		JSFillHist(param.Name, "2nd_electron_pt_" + param.Name, electrons.at(1).Pt(), weight, 16, 0., 160.);
		JSFillHist(param.Name, "3rd_electron_pt_" + param.Name, electrons.at(2).Pt(), weight, 16, 0., 160.);
		JSFillHist(param.Name, "1st_electron_eta_" + param.Name, electrons.at(0).Eta(), weight, 30, -3., 3.);
		JSFillHist(param.Name, "2nd_electron_eta_" + param.Name, electrons.at(1).Eta(), weight, 30, -3., 3.);
		JSFillHist(param.Name, "3rd_electron_eta_" + param.Name, electrons.at(2).Eta(), weight, 30, -3., 3.);
		JSFillHist(param.Name, "1st_electron_phi_" + param.Name, electrons.at(0).Phi(), weight, 32, -3.2, 3.2);
		JSFillHist(param.Name, "2nd_electron_phi_" + param.Name, electrons.at(1).Phi(), weight, 32, -3.2, 3.2);
		JSFillHist(param.Name, "3rd_electron_phi_" + param.Name, electrons.at(2).Phi(), weight, 32, -3.2, 3.2);
		JSFillHist(param.Name, "MET_" + param.Name, METv.Pt(), weight, 16, 0., 160.);
		JSFillHist(param.Name, "NJets_" + param.Name, jets.size(), weight, 7, -0.5, 6.5);
	}

	if (RunFakeSyst && param.Name.Contains("_Fake")) {
		double weight_fake = 1.;
		int loose_cnt = 0;
		for (unsigned int i = 0; i < electrons.size(); i++) {
			//==== Get Prompt Rate ====
			double this_prompt = -999;
			if (param.Name.Contains("Prompt") && IsDATA) 
				this_prompt = GetPromptRate(electrons.at(i), ElectronID, true);
			else if (param.Name.Contains("Prompt") && !IsData) 
				this_prompt = GetPromptRate(electrons.at(i), ElectronID, false);
			else this_prompt = 1.;
			//==== Get Fake Rate ====
			double this_fake = -999;
			if (param.Name.Contains("Fake_Central")) 
				this_fake = GetFakeRate(electrons.at(i), ElectronID, 0);
			else if (param.Name.Contains("Fake_SystUp")) 
				this_fake = GetFakeRate(electrons.at(i), ElectronID, 1);
			else if (param.Name.Contains("Fake_SystDown")) 
				this_fake = GetFakeRate(electrons.at(i), ElectronID, -1);
			else {
				cout << "[FakeValidator::executeEventFromParameter] param.Name = " << param.Name << endl;
				cout << "[FakeValidator::executeEventFromParameter] Wrong Syst" << endl;
				exit(EXIT_FAILURE);
			}
			//==== Apply as weight ====
			double this_weight = 1.;
			if (electrons.at(i).PassID(ElectronTightID)) {
				this_weight *= (((1 - this_fake) * this_prompt) / (this_prompt - this_fake));
			}
			else {
				loose_cnt++;
				this_weight *= ((-1 * this_prompt * this_fake) / (this_prompt - this_fake));
			}
			weight_fake *= this_weight;
		}
		weight_fake *= -1;
	/*double weight_fake = 1.;
	int loose_cnt = 0;
	for (unsigned int i = 0; i < electrons.size(); i++) {
		if (!electrons.at(i).PassID(ElectronTightID)) {
			double fr = -999;
			loose_cnt++;
			//cout << param.Name << endl;
			if (param.Name.Contains("Central")) fr = GetFakeRate(electrons.at(i), ElectronID, 0);
			else if (param.Name.Contains("Up")) fr = GetFakeRate(electrons.at(i), ElectronID, 1);
			else if (param.Name.Contains("Down")) fr = GetFakeRate(electrons.at(i), ElectronID, -1);
			else {
				cout << "[FakeValidator::executeEventFromParameter] Wrong Syst" << endl;
				exit(EXIT_FAILURE);
			}
			// for debug
			// fr = 0.1;
			weight_fake *= (-1 * (fr / (1 - fr)));
			
		}
		else continue;
	}
	weight_fake *= -1.;
	//cout << "loose_cnt " << loose_cnt << endl;
	//cout << "weight_fake" << weight_fake << endl;
	*/
		weight *= weight_fake;
		JSFillHist(param.Name, "nPassLooseID_" + param.Name, loose_cnt, 1, 4, -0.5, 3.5);

		//==== veto pass TTT by giving weight = 0 ====
		if (tightFlag) weight = 0;

		JSFillHist(param.Name, "1st_electron_pt_" + param.Name, electrons.at(0).Pt(), weight, 20, 0., 200.);
		JSFillHist(param.Name, "2nd_electron_pt_" + param.Name, electrons.at(1).Pt(), weight, 16, 0., 160.);
		JSFillHist(param.Name, "3rd_electron_pt_" + param.Name, electrons.at(2).Pt(), weight, 16, 0., 160.);
		JSFillHist(param.Name, "1st_electron_eta_" + param.Name, electrons.at(0).Eta(), weight, 30, -3., 3.);
		JSFillHist(param.Name, "2nd_electron_eta_" + param.Name, electrons.at(1).Eta(), weight, 30, -3., 3.);
		JSFillHist(param.Name, "3rd_electron_eta_" + param.Name, electrons.at(2).Eta(), weight, 30, -3., 3.);
		JSFillHist(param.Name, "1st_electron_phi_" + param.Name, electrons.at(0).Phi(), weight, 32, -3.2, 3.2);
		JSFillHist(param.Name, "1st_electron_phi_" + param.Name, electrons.at(0).Phi(), weight, 32, -3.2, 3.2);
		JSFillHist(param.Name, "1st_electron_phi_" + param.Name, electrons.at(0).Phi(), weight, 32, -3.2, 3.2);
		JSFillHist(param.Name, "MET_" + param.Name, METv.Pt(), weight, 16, 0., 160.);
		JSFillHist(param.Name, "Njets_" + param.Name, jets.size(), weight, 7, -0.5, 6.5);
	}
}

FakeValidator::FakeValidator(){

}

FakeValidator::~FakeValidator(){

}

//==== public functions ====
double FakeValidator::GetFakeRate(Electron &e, TString id, int sys) {
	
	double corrPt = GetCorrPt(e);
	double absEta = fabs(e.Eta());

	double value = 1.;
	double error = 0.;

	//==== corrPt and absEta should in the fiducial phase space ====
	if (corrPt <= 15.) corrPt = 16.;
	if (corrPt >= 70.) corrPt = 69.;
	if (absEta >= 2.5) absEta = 2.4;
 
	TH2D* this_hist = NULL;
	if (id.Contains("pass")) 	this_hist = (TH2D*)f->Get("Electron_fake_rate_POG");
	else if (id.Contains("Fake")) this_hist = (TH2D*)f->Get("Electron_fake_rate_FAKE");
	else {
		cout << "[FakeValidator::GetFakeRate] No fake rate histogram for " << id << endl;
		exit(EXIT_FAILURE);
	} 
	
	int this_bin = -999;
	this_bin = this_hist->FindBin(corrPt, absEta);

	value = this_hist->GetBinContent(this_bin);
	error = this_hist->GetBinError(this_bin);

	//cout << "[FakeValidator::GetFakeRate] corrPt = " << corrPt << endl;
	//cout << "[FakeValidator::GetFakeRate] absEta = " << absEta << endl;
	//cout << "[FakeValidator::GetFakeRate] value = " << value << endl;	

	return value + double(sys)*error;
}

double FakeValidator::GetPromptRate(Electron &e, TString id, bool IsData) {
	double corrPt = GetCorrPt(e);
	double absEta = fabs(e.Eta());

	double value = 1.;
	// double error = 0.;
	
	//==== corrPt and absEta should in the fiducial phase space ====
	if (corrPt <= 15.) corrPt = 16.;
	if (corrPt >= 70.) corrPt = 69.;
	if (absEta >= 2.5) absEta = 2.4;

	TH2D* this_hist = NULL;
	if (IsData) {
		if (id.Contains("pass")) this_hist = (TH2D*) f_prompt->Get("Electron_prompt_rate_POG_data");
		else if (id.Contains("Fake")) this_hist = (TH2D*) f_prompt->Get("Electron_prompt_rate_Fake_data");
		else {
			cout << "[FakeValidation::GetPromptRate] No prompt rate histogram for " << id << endl;
			exit(EXIT_FAILURE);
		}
	}
	else {
		if (id.Contains("pass")) this_hist = (TH2D*) f_prompt->Get("Electron_prompt_rate_POG_MC");
		else if (id.Contains("Fake")) this_hist = (TH2D*) f_prompt->Get("Electron_prompt_rate_Fake_MC");
		else {
			cout << "[FakeValidation::GetPromptRate] No prompt rate histogram for " << id << endl;
			exit(EXIT_FAILURE);
		}
	}

	int this_bin = -999;
	this_bin = this_hist->FindBin(corrPt, absEta);
	value = this_hist->GetBinContent(this_bin);

	return value;
}

double FakeValidator::GetCorrPt(Electron &e) {
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
