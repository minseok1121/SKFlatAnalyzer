#include "ConvEstimator.h"

void ConvEstimator::initializeAnalyzer(){
	//==== User flags ====
	RunFakeSyst = HasFlag("RunFakeSyst");
	
	cout << "[ConvEstimator::initializeAnalyzer] RunFakeSyst = " << RunFakeSyst << endl;
	
	//==== Electron ID setting ====
	ElectronIDs = {"passLooseID", "passTightID", "FakeLooseID", "FakeTightID"};
	idsets = {"POG", "Fake"};

	//==== Trigger setting ====
	if (DataYear == 2016) {
		HLTElecTriggerName = "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v";
		TriggerSafePtCut1 = 25.;
        TriggerSafePtCut2 = 15.;
    }
    else if (DataYear == 2017) {
        cout << "[ConvEstimator::initializeAnalyzer] Trigger is not set for 2017" << endl;
        exit(EXIT_FAILURE);
    }
    else if (DataYear == 2018) {
        cout << "[ConvEstimator::initializeAnalyzer] Trigger is not set for 2018" << endl;
        exit(EXIT_FAILURE);
    }
    else {
        cout << "[ConvEstimator::initializeAnalyzer] Wrong Year" << endl;
        exit(EXIT_FAILURE);
    }

	cout << "[ConvEstimator::initializeAnalyzer] HLTElecTriggerName  = " << HLTElecTriggerName << endl;
	cout << "[ConvEstimator::initializeAnalyzer] TriggerSafePtCut1 = " << TriggerSafePtCut1 << endl;
	cout << "[ConvEstimator::initializeAnalyzer] TriggerSafePtCut2 = " << TriggerSafePtCut2 << endl;
}

void ConvEstimator::executeEvent(){

	//==== Copy all objects ====
	AllMuons = GetAllMuons();
    AllElectrons = GetAllElectrons();
    AllJets = GetAllJets();

    //==== Get L1Prefire, pileup reweight ====
    weight_Prefire = GetPrefireWeight(0);
    weight_PileUp = GetPileUpWeight(nPileUp, 0);

    AnalyzerParameter param;	

	//==== Loop over ElectronID ====
	for (unsigned int i = 0; i < idsets.size(); i++) {
		MuonID = "POGLoose";
		if (idsets.at(i).Contains("POG")) {
			ElectronID = ElectronIDs.at(i);
			ElectronTightID = ElectronIDs.at(i+1);
		}
		else if (idsets.at(i).Contains("Fake")) {
			ElectronID = ElectronIDs.at(2*i);
			ElectronTightID = ElectronIDs.at(2*1 + 1);
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
	}
}

void ConvEstimator::executeEventFromParameter(AnalyzerParameter param){

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
	jets = SelectJets(this_AllJets, param.Jet_ID, 20, 2.4);

	std::sort(muons.begin(), muons.end(), PtComparing);
    std::sort(electrons.begin(), electrons.end(), PtComparing);
    std::sort(jets.begin(), jets.end(), PtComparing);

	//==== Event Selection ====
	//==== 1. 3 electrons with at least one opposite-sign electron pair ====
	//==== 2. at least one pair off from Z mass more than 10 GeV ====
	//==== 3. include OS electron pair 12 < M(e+e-) < 81.2 GeV ====
	//==== 4. |M(3e) - 91.2| < 10 GeV ====
	//==== 5. Missing ET < 50GeV ====
	if (muons.size() != 0) return;
	if (electrons.size() != 3) return;

	JSFillHist(param.Name, "leptonCut_" + param.Name, 0., 1., 1, 0., 1.);

	if (electrons.at(0).Pt() < TriggerSafePtCut1) return;
	if (electrons.at(1).Pt() < TriggerSafePtCut2) return;
	if (electrons.at(2).Pt() < TriggerSafePtCut2) return;

	Pair electronPair[3];
	int pairCnt = 0;
	bool isOnZ = false;

	for (int i = 0; i < 3; i++) {
		electronPair[i].Mass = (electrons.at(i%3) + electrons.at((i+1)%3)).M();
		electronPair[i].isOS = false;
		electronPair[i].isOffZ = false;
		electronPair[i].isBelowZ = false;
		
		// Charge of pair
		if (electrons.at(i%3).Charge() * electrons.at((i+1)%3).Charge() < 0) electronPair[i].isOS = true;
		if (!IsOnZ(electronPair[i].Mass, 10.)) electronPair[i].isOffZ = true;
		if (12 < electronPair[i].Mass && electronPair[i].Mass < 81.2) electronPair[i].isBelowZ = true;

		if (electronPair[i].isOS && electronPair[i].isOffZ && electronPair[i].isBelowZ) pairCnt++;
	}
	Particle triplet = electrons.at(0) + electrons.at(1) + electrons.at(2);
	if (IsOnZ(triplet.M(), 10)) isOnZ = true;

	if (pairCnt == 0) return;
	if (!isOnZ) return;
	if (METv.Pt() > 50) return;

	//==== Event weight ====
	double weight = 1.;
	if (!IsDATA) {
		weight *= weight_norm_1invpb * ev.GetTriggerLumi("Full");
		weight *= ev.MCweight();
		weight *= weight_Prefire;
		weight *= weight_PileUp;

		//==== ID and trigger SF ====
        //==== only for POG ID ====
		//==== clear it...====
        if (ElectronID.Contains("pass")) {
            //==== Trigger SF ====
            double eff_DATA = 1., eff_MC = 1.;
            for (unsigned int i = 0; i < electrons.size()-1; i++) {
                if (electrons.at(i).PassID(ElectronID)) {
                    eff_DATA *= (1. - mcCorr->ElectronTrigger_Eff(ElectronID, HLTElecTriggerName, 0, electrons.at(i).scEta(), electrons.at(i).Pt(), 0));
                    eff_MC *= (1. - mcCorr->ElectronTrigger_Eff(ElectronID, HLTElecTriggerName, 1, electrons.at(i).scEta(), electrons.at(i).Pt(), 0));
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
                if (electrons.at(i).PassID(ElectronID)) {
                    this_ElecIDSF = mcCorr->ElectronID_SF(ElectronID, electrons.at(i).scEta(), electrons.at(i).Pt(), 0);
                }
                else this_ElecIDSF = mcCorr->ElectronID_SF(ElectronID, electrons.at(i).scEta(), electrons.at(i).Pt(), 0);
                weight *= this_ElecIDSF;
            }
        }
	}

	bool tightFlag = (electrons.at(0).PassID(ElectronTightID) && electrons.at(1).PassID(ElectronTightID) && electrons.at(2).PassID(ElectronTightID));
	

	//==== Fill hist ====
	if (param.Name.Contains("Central") && !param.Name.Contains("_Fake") && tightFlag ) {
		JSFillHist(param.Name, "PassSelection_" + param.Name, 1., weight, 2, -0.5, 1.5);
		JSFillHist(param.Name, "pairCnt_" + param.Name, pairCnt, weight, 4, -0.5, 3.5);
		JSFillHist(param.Name, "1st_electron_pt_" + param.Name, electrons.at(0).Pt(), weight, 40, 0., 200.);
		JSFillHist(param.Name, "2nd_electron_pt_" + param.Name, electrons.at(1).Pt(), weight, 32, 0., 160.);
		JSFillHist(param.Name, "3rd_electron_pt_" + param.Name, electrons.at(2).Pt(), weight, 32, 0., 160.);
		JSFillHist(param.Name, "1st_electron_eta_" + param.Name, electrons.at(0).Eta(), weight, 30, -3., 3.);
		JSFillHist(param.Name, "2nd_electron_eta_" + param.Name, electrons.at(1).Eta(), weight, 30, -3., 3.);
		JSFillHist(param.Name, "3rd_electron_eta_" + param.Name, electrons.at(2).Eta(), weight, 30, -3., 3.);
		JSFillHist(param.Name, "1st_electron_phi_" + param.Name, electrons.at(0).Phi(), weight, 32, -3.2, 3.2);
		JSFillHist(param.Name, "2nd_electron_phi_" + param.Name, electrons.at(1).Phi(), weight, 32, -3.2, 3.2);
		JSFillHist(param.Name, "3rd_electron_phi_" + param.Name, electrons.at(2).Phi(), weight, 32, -3.2, 3.2);
		JSFillHist(param.Name, "M(3e)_" + param.Name, triplet.M(), weight, 20, 70., 110.);
		JSFillHist(param.Name, "MET_" + param.Name, METv.Pt(), weight, 32, 0., 160.);
		JSFillHist(param.Name, "NJets_" + param.Name, jets.size(), weight, 7, -0.5, 6.5);

		for (int i = 0; i < 3; i++) {
			if (electronPair[i].isOS && electronPair[i].isOffZ && electronPair[i].isBelowZ) {
				JSFillHist(param.Name, "M(e+e-)_" + param.Name, electronPair[i].Mass, weight, 50, 0, 100);
			}
		}
	}

	if (RunFakeSyst && param.Name.Contains("_Fake")) {

		double weight_fake = 1.;
		int loose_cnt = 0;
		for (unsigned int i = 0; i < electrons.size(); i++) {
			//==== Get Prompt rate ====
			//==== p = 1 approximation applied ====
			double this_prompt = 1.;

			//==== Get Fake Rate ====
			double this_fake = -999.;
			if (param.Name.Contains("Fake_Central")) 
				this_fake = GetFakeRate(electrons.at(i), ElectronID, 0);
			else if (param.Name.Contains("Fake_SystUp")) 
				this_fake = GetFakeRate(electrons.at(i), ElectronID, 1);
			else if (param.Name.Contains("Fake_SystDown")) 
				this_fake = GetFakeRate(electrons.at(i), ElectronID, -1);
			else {
				cout << "[FakeValidWithPrompt::executeEventFromParameter] Wrong Syst" << endl;
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
		if (tightFlag) weight_fake = 0;

		weight *= weight_fake;
		JSFillHist(param.Name, "PassSelection_" + param.Name, 1., weight, 2, -0.5, 1.5);
		JSFillHist(param.Name, "nPassLoose_" + param.Name, loose_cnt, 1, 4, -0.5, 3.5);
		JSFillHist(param.Name, "pairCnt_" + param.Name, pairCnt, weight, 4, -0.5, 3.5);
		JSFillHist(param.Name, "1st_electron_pt_" + param.Name, electrons.at(0).Pt(), weight, 40, 0., 200.);
		JSFillHist(param.Name, "2nd_electron_pt_" + param.Name, electrons.at(1).Pt(), weight, 32, 0., 160.);
		JSFillHist(param.Name, "3rd_electron_pt_" + param.Name, electrons.at(2).Pt(), weight, 32, 0., 160.);
		JSFillHist(param.Name, "1st_electron_eta_" + param.Name, electrons.at(0).Eta(), weight, 30, -3., 3.);
		JSFillHist(param.Name, "2nd_electron_eta_" + param.Name, electrons.at(1).Eta(), weight, 30, -3., 3.);
		JSFillHist(param.Name, "3rd_electron_eta_" + param.Name, electrons.at(2).Eta(), weight, 30, -3., 3.);
		JSFillHist(param.Name, "1st_electron_phi_" + param.Name, electrons.at(0).Phi(), weight, 32, -3.2, 3.2);
		JSFillHist(param.Name, "2nd_electron_phi_" + param.Name, electrons.at(1).Phi(), weight, 32, -3.2, 3.2);
		JSFillHist(param.Name, "3rd_electron_phi_" + param.Name, electrons.at(2).Phi(), weight, 32, -3.2, 3.2);
		JSFillHist(param.Name, "M(3e)_" + param.Name, triplet.M(), weight, 20, 70., 110.);
		JSFillHist(param.Name, "MET_" + param.Name, METv.Pt(), weight, 32, 0., 160.);
		JSFillHist(param.Name, "NJets_" + param.Name, jets.size(), weight, 7, -0.5, 6.5);
	}

}

ConvEstimator::ConvEstimator(){

}

ConvEstimator::~ConvEstimator(){

}

//==== Member functions ====
double ConvEstimator::GetFakeRate(Electron &e, TString id, int sys) {

    double corrPt = GetCorrPt(e);
    double absEta = fabs(e.Eta());

    double value = 1.;
    double error = 0.;

    //==== corrPt and absEta should in the fiducial phase space ====
    if (corrPt <= 15.) corrPt = 16.;
    if (corrPt >= 70.) corrPt = 69.;
    if (absEta >= 2.5) absEta = 2.4;

    TH2D* this_hist = NULL;
    if (id.Contains("pass"))    this_hist = (TH2D*)f->Get("Electron_fake_rate_POG");
    else if (id.Contains("Fake")) this_hist = (TH2D*)f->Get("Electron_fake_rate_FAKE");
    else {
        cout << "[ConvEstimator::GetFakeRate] No fake rate histogram for " << id << endl;
        exit(EXIT_FAILURE);
    }

    int this_bin = -999;
    this_bin = this_hist->FindBin(corrPt, absEta);

    value = this_hist->GetBinContent(this_bin);
    error = this_hist->GetBinError(this_bin);

    //cout << "[ConvEstimator::GetFakeRate] corrPt = " << corrPt << endl;
    //cout << "[ConvEstimator::GetFakeRate] absEta = " << absEta << endl;
    //cout << "[ConvEstimator::GetFakeRate] value = " << value << endl;

    return value + double(sys)*error;
}

double ConvEstimator::GetCorrPt(Electron &e) {
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
