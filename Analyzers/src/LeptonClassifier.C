#include "LeptonClassifier.h"

void LeptonClassifier::initializeAnalyzer(){
	// initialize nEvent to 0
	nEvent = 0;

	// IDs - for taggins
	MuonIDs = {"NOCUT", "POGLooseWithLooseIso", "POGTightWithTightIso"};
	ElectronIDs = {"NOCUT", "passLooseID", "passTightID"};
	JetIDs = {"NOCUT", "tight"};

	// b-tagging
	vector<JetTagging::Parameters> jtps;
	jtps.push_back(JetTagging::Parameters(JetTagging::DeepCSV, JetTagging::Medium, JetTagging::incl, JetTagging::mujets));
	mcCorr->SetJetTaggingParameters(jtps);
}

void LeptonClassifier::executeEvent(){
	nEvent++;

	TString path = "preselection/";
	FillHist(path + "cutflow", 0., 1., 9, 0., 9.);

	// clear all the containers first
	ClearCollections();
	InitializeWeights();
	
	// get objects
	Event ev = GetEvent();
	muons = GetMuons("NOCUT", 5., 2.4);
	electrons = GetElectrons("NOCUT", 5., 2.5);
	jets = GetJets("NOCUT", 5., 2.4);
	gens = GetGens();
	METv = ev.GetMETVector();

	sort(muons.begin(), muons.end(), PtComparing);
	sort(electrons.begin(), electrons.end(), PtComparing);
	sort(jets.begin(), jets.end(), PtComparing);

	if (!PassMETFilter()) return;
	FillHist(path + "cutflow", 1., 1., 9, 0., 9.);

	muons_tight = SelectMuons(muons, "POGTightWithTightIso", 30., 2.4);
	muons_veto = SelectMuons(muons, "POGLooseWithLooseIso", 25., 2.4);
	electrons_tight = SelectElectrons(electrons, "passTightID", 30., 2.5);
	electrons_veto = SelectElectrons(electrons, "passLooseID", 25., 2.5);
	jets_tight = SelectJets(jets, "tight", 25., 2.4);
	jets_dR04 = JetsVetoLeptonInside(jets_tight, electrons_veto, muons_veto, 0.4);
	bjets = jets_tight;
	for (const auto& j : bjets) {
		double this_discr = j.GetTaggerResult(JetTagging::DeepCSV);
		if (this_discr > mcCorr->GetJetTaggingCutValue(JetTagging::DeepCSV, JetTagging::Medium))
			bjets_tight.push_back(j);
	}
	bjets_dR04 = JetsVetoLeptonInside(bjets_tight, electrons_veto, muons_veto, 0.4);

	// preselection
	if (! (muons_tight.size() <= 2)) return;
	if (! (electrons_tight.size() <= 2)) return;
	FillHist(path + "cutflow", 2., 1., 9, 0., 9.);
	if (! (jets_dR04.size() >= 2)) return;
	FillHist(path + "cutflow", 3., 1., 9, 0., 9.);
	if (! (bjets_dR04.size() >= 1)) return;
	FillHist(path + "cutflow", 4., 1., 9, 0., 9.);

	// MC weights
	double weight = 1.;
	if (!IsDATA) {
		w_prefire = GetPrefireWeight(0);
		w_gen = weight_norm_1invpb*ev.MCweight();
		w_lumi = ev.GetTriggerLumi("Full");
	}
	weight *= w_prefire * w_gen * w_lumi;

	// gen level study
	// Gens = RewriteGens();
	// PrintHistory(Lep);
	//PrintHistory(Gamma);
	PrintHistory(CH);

	// To do:
	// 1. Signal / Control region selection
	// 2. DrawHists


}

void LeptonClassifier::ClearCollections() {
	muons.clear(); muons_tight.clear(); muons_veto.clear();
	electrons.clear(); electrons_tight.clear(); electrons_veto.clear();
	jets.clear(); jets_tight.clear(); jets_dR04.clear();
	bjets.clear(); bjets_tight.clear(); bjets_dR04.clear();
	gens.clear(); Gens.clear();
}

void LeptonClassifier::InitializeWeights() {
	w_prefire = 1.; w_gen = 1.; w_lumi = 1.;
}

vector<LeptonClassifier::Gen_t> LeptonClassifier::RewriteGens() {
	vector<Gen_t> out;
	for (auto& g: gens) {
		Gen_t this_gen;
		this_gen.me = g;
		this_gen.mom = gens.at(g.MotherIndex());
		for (auto& son: gens) {
			if (son.MotherIndex() == g.Index()) this_gen.sons.push_back(gens.at(son.Index()));
			else continue;
		}
		out.push_back(this_gen);
	}

	if (gens.size() != out.size()) {
		cerr << "[LeptonClassifier::RewriteGens] WARNING: missed history" << endl;
	}

	return out;
}

vector<vector<int>> LeptonClassifier::MakeHistory() {
	int my_idx = -999;
	int mom_idx = -999;

	// get final states, i.e. status == 1
	vector<Gen> endpoints;
	for (const auto& g: gens) {
		if (g.Status() == 1) endpoints.push_back(g);
		else continue;
	}

	// now trace back history
	vector<vector<int>> history;
	for (const auto& endpoint: endpoints) {
		vector<int> this_history;

		// initialize index
		my_idx = endpoint.Index();
		mom_idx = endpoint.MotherIndex();
		
		while (true) {
			this_history.push_back(my_idx);

			// go one generation up
			const Gen& mother = gens.at(mom_idx);
			my_idx = mother.Index();
			mom_idx = mother.MotherIndex();

			if (mom_idx == -1) {
				this_history.push_back(my_idx); // proton
				break;
			}
			else continue;
		}
		history.push_back(this_history);
	}

	return history;
}
	
void LeptonClassifier::PrintHistory(const FinalState_t& state) {
	vector<vector<int>> history = MakeHistory();
	cout << "[LeptonClassifier::PrintHistory] nEvent = " << nEvent << endl;

	// print
	unsigned int n_history = 0;
	for (const auto& this_history: history) {
		int endpoint_idx = this_history.at(0);
		const Gen& endpoint = gens.at(endpoint_idx);

		// check endpoint PID
		if (state == Lep) {
			if (! (endpoint.PID() == 11 || endpoint.PID() == 13)) continue;
		}
		else if (state == Gamma) {
			if (! (endpoint.PID() == 22)) continue;
		}
		else {
			cerr << "[LeptonClassifier::PrintHistory] current type is not set yet" << endl;
			exit(EXIT_FAILURE);
		}
		
		cout << "[LeptonClassiffier::PrintHistory] n_history = " << n_history << endl;
		for (const auto& idx: this_history) {
			gens.at(idx).Print();
		}
		cout << endl;
		n_history++;
	}
	cout << endl << endl;
}


LeptonClassifier::LeptonClassifier(){}

LeptonClassifier::~LeptonClassifier(){}


