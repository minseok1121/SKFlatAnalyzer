#include "SigToBkg.h"

SigToBkg::SigToBkg() {} // called? let's think about the cutflow
SigToBkg::~SigToBkg() {
	// not sure deleting CutflowMakers would affect FillHist
	//map<SIGNAL, HistoMaker*>::iterator it;
	//for (it=makers.begin(); it!=makers.end(); it++)
	//	delete it->second;
}
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

	// initiate cutflows and histomakers
	cuts = {"nocut", "metfilter", "trigger", "nleptons", "passSafePtCut", "exist_osmu", "off_Zmass", "Nj_ge2", "Nb_ge1"};
	makers[SR3MU] = new HistoMaker("3mu", cuts);
	makers[SR1E2MU]= new HistoMaker("1e2mu", cuts);
}

void SigToBkg::executeEvent(){
	//FillCutflow
	map<SIGNAL, HistoMaker*>::iterator it;
	for (it=makers.begin(); it!=makers.end(); it++)
		it->second->FillCutflow("nocut");

	if (!PassMETFilter()) return;
	for (it=makers.begin(); it!=makers.end(); it++)
		it->second->FillCutflow("metfilter");
	
	Event ev = GetEvent();
	// no triggers yet
	for (it=makers.begin(); it!=makers.end(); it++)
		it->second->FillCutflow("trigger");

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
	SIGNAL signal = SignalRegionSelector(
			muons_tight, electrons_tight, 
			muons_loose, electrons_loose, 
			jets_cleaned, bjets_cleaned);
	
	if (signal == NONE) return;
	
	// set weight
	double weight = 1.;
	// Finally, fill histograms
	// auto* maker = makers[signal];	
}

//==== This is actual selector!
//==== TODO: might divide in more functions?
SigToBkg::SIGNAL SigToBkg::SignalRegionSelector(
		vector<Muon> &muons_tight, vector<Electron> &electrons_tight,
		vector<Muon> &muons_loose, vector<Electron> &electrons_loose,
		vector<Jet> &jets, vector<Jet> &bjets) {
	// let's count leptons first
	SIGNAL signal = NONE;
	if (muons_tight.size() == 3 && muons_loose.size() == 3) 
		signal = SR3MU;
	if (muons_tight.size() == 2 && muons_loose.size() == 2 && electrons_tight.size() == 1 && electrons_loose.size() == 1)
		signal = SR1E2MU;
	
	if (signal != NONE)
		makers[signal]->FillCutflow("nleptons");

	const double ZMass = 91.2;
	if (signal == SR3MU) {
		if (muons_tight.at(0).Pt() < 20.)
			return NONE;
		if (muons_tight.at(1).Pt() < 10.)
			return NONE;
		if (muons_tight.at(2).Pt() < 10.)
			return NONE;
		makers[signal]->FillCutflow("passSafePtCut");

		int chargeSum
			= muons_tight.at(0).Charge() + muons_tight.at(1).Charge() + muons_tight.at(2).Charge();
		if (abs(chargeSum) != 1)
			return NONE;
		makers[signal]->FillCutflow("exist_osmu");
		
		Particle ACand1 = muons_tight.at(0) + muons_tight.at(1);
		Particle ACand2 = muons_tight.at(1) + muons_tight.at(2);
		Particle ACand3 = muons_tight.at(2) + muons_tight.at(0);
		if (fabs(ACand1.M() - ZMass) < 10.) return NONE;
		if (fabs(ACand2.M() - ZMass) < 10.) return NONE;
		if (fabs(ACand3.M() - ZMass) < 10.) return NONE;
		makers[signal]->FillCutflow("off_Zmass");

		if (jets.size() < 2) return NONE;
		makers[signal]->FillCutflow("Nj_ge2");
		if (bjets.size() < 1) return NONE;
		makers[signal]->FillCutflow("Nb_ge1");
		return signal;
	}


	else if (signal == SR1E2MU) {
		if (muons_tight.at(1).Pt() < 10.)
			return NONE;
		if (electrons_tight.at(0).Pt() < 15.)
			return NONE;
		makers[signal]->FillCutflow("passSafePtCut");

		int chargeSum = muons_tight.at(0).Charge() + muons_tight.at(1).Charge();
		if (chargeSum != 0) 
			return NONE;
		makers[signal]->FillCutflow("exist_osmu");

		Particle ACand = muons_tight.at(0) + muons_tight.at(1);
		if (fabs(ACand.M() - ZMass) < 10.) 
			return NONE;
		makers[signal]->FillCutflow("off_Zmass");

		if (jets.size() < 2) 
			return NONE;
		makers[signal]->FillCutflow("Nj_ge2");
		if (bjets.size() < 1)
			return NONE;
		makers[signal]->FillCutflow("Nb_ge1");
		return signal;
	}
	else
		return NONE;
}

