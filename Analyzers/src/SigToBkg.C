#include "SigToBkg.h"

SigToBkg::SigToBkg() {}
SigToBkg::~SigToBkg() {}
void SigToBkg::initializeAnalyzer(){
	// TODO: Trigger is not set since I don't want to make it biased
	// Let's see how triggers affect the S/Sqrt(B) later
	
	// IDs
	HcToWA_MuID = {"HcToWATight", "HcToWALoose"};
	HcToWA_EleID = {"HcToWATight", "HcToWALoose"};

	// B-tagging
	vector<JetTagging::Parameters> jtps;
	jtps.emplace_back(
		JetTagging::Parameters(JetTagging::DeepCSV, JetTagging::Medium, 
			JetTagging::incl, JetTagging::mujets));
	mcCorr->SetJetTaggingParameters(jtps);

	// initiate regions and cutflows
	regions = {"3mu", "1e2mu"};
	cuts = {"nocut", "metfilter", "trigger", "nleptons", "passSafePtCut", "exist_osmu", "off_Zmass", "Nj_ge2", "Nb_ge1"};
	InitiateCutflow("3mu", cuts);
	InitiateCutflow("1e2mu", cuts);
}

void SigToBkg::executeEvent(){
	//FillCutflow
	for (const auto& region: regions)
		FillCutflow(region, "nocut");

	if (!PassMETFilter()) return;
	for (const auto& region: regions)
		FillCutflow(region, "metfilter");
	
	Event ev = GetEvent();
	// no triggers yet
	for (const auto& region: regions)
		FillCutflow(region, "trigger");

	// Let's define objects
	vector<Gen> gens = GetGens();
	vector<Muon> muons= GetAllMuons();
	vector<Electron> electrons = GetAllElectrons();
	vector<Jet> jets = GetAllJets();
	Particle METv = ev.GetMETVector();
	
	// sort at the first time, I don't want to get confused
	sort(muons.begin(), muons.end(), PtComparing);
	sort(electrons.begin(), electrons.end(), PtComparing);
	sort(jets.begin(), jets.end(), PtComparing);

	// select objects
	vector<Muon> muons_tight = SelectMuons(muons, "HcToWATight", 10., 2.4);
	vector<Muon> muons_loose = SelectMuons(muons, "HcToWALoose", 10., 2.4);
	vector<Electron> electrons_tight = SelectElectrons(electrons, "HcToWATight", 10., 2.5);
	vector<Electron> electrons_loose = SelectElectrons(electrons, "HcToWALoose", 10., 2.5);
	vector<Jet> jets_tight = SelectJets(jets, "tight", 20., 2.4);
	vector<Jet> jets_cleaned = JetsVetoLeptonInside(jets_tight, electrons_loose, muons_loose, 0.4);
	vector<Jet> bjets_cleaned;
	for (const auto& jet: jets_cleaned) {
		double this_discr = jet.GetTaggerResult(JetTagging::DeepCSV);
		if (this_discr > mcCorr->GetJetTaggingCutValue(JetTagging::DeepCSV, JetTagging::Medium))
			bjets_cleaned.emplace_back(jet);
	}

	// NOTE: Cutflow will automatically generated inside the SignalSelector
	TString signal = SignalRegionSelector(
			muons_tight, electrons_tight, 
			muons_loose, electrons_loose, 
			jets_cleaned, bjets_cleaned);
	
	if (signal=="") return;
	
	// Let's categorize leptons
	vector<Muon> muons_signal, muons_ewprompt, muons_fromOffW, muons_fromtau, muons_fake;
	vector<Electron> electrons_signal, electrons_ewprompt, electrons_fromOffW, electrons_fromtau, electrons_fake;

	for (const auto& mu: muons_tight) {
		const int LepType = GetLeptonType(mu, gens);
		if (LepType == 1)
			muons_ewprompt.emplace_back(mu);
		else if (LepType == 2)
			muons_signal.emplace_back(mu);
		else if (LepType == 6)
			muons_fromOffW.emplace_back(mu);
		else if (LepType == 3)
			muons_fromtau.emplace_back(mu);
		else if (LepType < 0)
			muons_fake.emplace_back(mu);
		else
			continue;
	}
	for (const auto& ele: electrons_tight) {
		const int LepType = GetLeptonType(ele, gens);
		if (LepType == 1)
			electrons_ewprompt.emplace_back(ele);
		else if (LepType == 2)
			electrons_signal.emplace_back(ele); // should be empty
		else if (LepType == 6)
			electrons_fromOffW.emplace_back(ele);
		else if (LepType == 3)
			electrons_fromtau.emplace_back(ele);
		else if (LepType < 0)
			electrons_fake.emplace_back(ele);
		else
			continue;
	}

	// set weight
	double weight = 1.;
	double w_prefire, w_gen, w_lumi;
	w_prefire = GetPrefireWeight(0);
	w_gen = ev.MCweight() * weight_norm_1invpb;
	w_lumi = ev.GetTriggerLumi("Full");
	weight *= w_prefire*w_gen*w_lumi;

	// Finally, fill histograms
	// TODO: ACnad?
	FillObjects(signal + "/muons_tight", muons_tight, weight);
	FillObjects(signal + "/muons_ewprompt", muons_ewprompt, weight);
	FillObjects(signal + "/muons_signal", muons_signal, weight);
	FillObjects(signal + "/muons_fromOffW", muons_fromOffW, weight);
	FillObjects(signal + "/muons_fromtau", muons_fromtau, weight);
	FillObjects(signal + "/muons_fake", muons_fake, weight);
	FillObjects(signal + "/electrons_tight", electrons_tight, weight);
	FillObjects(signal + "/electrons_ewprompt", electrons_ewprompt, weight);
	FillObjects(signal + "/electrons_signal", electrons_signal, weight);
	FillObjects(signal + "/electrons_fromOffW", electrons_fromOffW, weight);
	FillObjects(signal + "/electrons_fromtau", electrons_fromtau, weight);
	FillObjects(signal + "/electrons_fake", electrons_fake, weight);
	FillObjects(signal + "/jets_cleaned", jets_cleaned, weight);
	FillObjects(signal + "/bjets_cleaned", bjets_cleaned, weight);
	FillObject(signal + "/METv", METv, weight);
}

//==== This is actual selector!
//==== TODO: might divide in more functions?
TString SigToBkg::SignalRegionSelector(
		vector<Muon> &muons_tight, vector<Electron> &electrons_tight,
		vector<Muon> &muons_loose, vector<Electron> &electrons_loose,
		vector<Jet> &jets, vector<Jet> &bjets) {
	// let's count leptons first
	TString signal = "";
	if (muons_tight.size() == 3 && muons_loose.size() == 3) 
		signal = "3mu";
	if (muons_tight.size() == 2 && muons_loose.size() == 2 && electrons_tight.size() == 1 && electrons_loose.size() == 1)
		signal = "1e2mu";
	
	if (signal != "") {
		for (const auto& region: regions)
			FillCutflow(region, "nleptons");
	}

	const double ZMass = 91.2;
	if (signal == "3mu") {
		if (muons_tight.at(0).Pt() < 20.)
			return "";
		if (muons_tight.at(1).Pt() < 10.)
			return "";
		if (muons_tight.at(2).Pt() < 10.)
			return "";
		FillCutflow(signal, "passSafePtCut");

		int chargeSum
			= muons_tight.at(0).Charge() + muons_tight.at(1).Charge() + muons_tight.at(2).Charge();
		if (abs(chargeSum) != 1)
			return "";
		FillCutflow(signal, "exist_osmu");
		
		Particle ACand1 = muons_tight.at(0) + muons_tight.at(1);
		Particle ACand2 = muons_tight.at(1) + muons_tight.at(2);
		Particle ACand3 = muons_tight.at(2) + muons_tight.at(0);
		if (fabs(ACand1.M() - ZMass) < 10.) return "";
		if (fabs(ACand2.M() - ZMass) < 10.) return "";
		if (fabs(ACand3.M() - ZMass) < 10.) return "";
		FillCutflow(signal, "off_Zmass");

		if (jets.size() < 2) return "";
		FillCutflow(signal, "Nj_ge2");
		if (bjets.size() < 1) return "";
		FillCutflow(signal, "Nb_ge1");
		return signal;
	}


	else if (signal == "1e2mu") {
		if (muons_tight.at(1).Pt() < 10.)
			return "";
		if (electrons_tight.at(0).Pt() < 15.)
			return "";
		FillCutflow(signal, "passSafePtCut");

		int chargeSum = muons_tight.at(0).Charge() + muons_tight.at(1).Charge();
		if (chargeSum != 0) 
			return "";
		FillCutflow(signal, "exist_osmu");

		Particle ACand = muons_tight.at(0) + muons_tight.at(1);
		if (fabs(ACand.M() - ZMass) < 10.) 
			return "";
		FillCutflow(signal, "off_Zmass");

		if (jets.size() < 2) 
			return "";
		FillCutflow(signal, "Nj_ge2");
		if (bjets.size() < 1)
			return "";
		FillCutflow(signal, "Nb_ge1");
		return signal;
	}
	else
		return "";
}

