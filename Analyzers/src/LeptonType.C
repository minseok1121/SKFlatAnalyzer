#include "LeptonType.h"

void LeptonType::initializeAnalyzer(){

	cout << "[LeptonType::initializeAnalyzer] No Flag" << endl;

	// nEvent initialized with 0;
	nEvent = 0;
	// IDs
	// not set for yet
	MuonIDs = {"NOCUT", "POGLooseWithTightIso", "POGTightWithTightIso"};
	ElectronIDs = {"NOCUT", "passLooseID", "passTightID"};
	JetIDs = {"NOCUT", "tight"};

	// b tagging
	vector<JetTagging::Parameters> jtps;
	jtps.push_back(JetTagging::Parameters(JetTagging::DeepCSV, JetTagging::Medium, JetTagging::incl, JetTagging::mujets));
	mcCorr->SetJetTaggingParameters(jtps);
}

void LeptonType::executeEvent(){

	nEvent++;
	// **important** clear all vector containers first
	ClearPreColl();
	muons = GetMuons("NOCUT", 5., 2.4);
	electrons = GetElectrons("NOCUT", 5., 2.5);
	jets = GetAllJets();
	gens = GetGens();
	w_prefire = GetPrefireWeight(0);

	// sort jet collections, so don't need to sort again
	sort(muons.begin(), muons.end(), PtComparing);
	sort(electrons.begin(), electrons.end(), PtComparing);
	sort(jets.begin(), jets.end(), PtComparing);

	AnalyzerParameter param;
	param.Clear();
	param.syst_ = AnalyzerParameter::Central;

	executeEventFromParameter(param);
}

void LeptonType::executeEventFromParameter(AnalyzerParameter param){

	TString path = "preselection/";
	FillHist(path + "cutflow", 0., 1., 9, 0., 9.);
	
	if(!PassMETFilter()) return;
	FillHist(path + "cutflow", 1., 1., 9, 0., 9.);

	ClearSubColl();
	Event ev = GetEvent();
	muons_tight = SelectMuons(muons, "POGTightWithTightIso", 30., 2.4);
	muons_veto = SelectMuons(muons, "POGLooseWithTightIso", 25., 2.4);
	electrons_tight = SelectElectrons(electrons, "passTightID", 30., 2.4);
	electrons_veto = SelectElectrons(electrons, "passLooseID", 25., 2.4);
	jets_tight = SelectJets(jets, "tight", 25., 2.4);
	for (const auto& j: jets_tight) {
		double this_discr = j.GetTaggerResult(JetTagging::DeepCSV);
		if (this_discr >  mcCorr->GetJetTaggingCutValue(JetTagging::DeepCSV, JetTagging::Medium)) bjets_tight.push_back(j);
	}
	jets_dR04 = JetsVetoLeptonInside(jets_tight, electrons_veto, muons_veto, 0.4);
	bjets_dR04 = JetsVetoLeptonInside(bjets_tight, electrons_veto, muons_veto, 0.4);
	METv = ev.GetMETVector();

	// preselection
	if (muons_tight.size() > 2) return;
	if (electrons_tight.size() > 2) return;
	FillHist(path + "cutflow", 2., 1., 9, 0., 9.);
	if (jets_dR04.size() < 2) return;
	if (bjets_dR04.size() < 1) return;
	FillHist(path + "cutflow", 3., 1., 9, 0., 9.);

	// weight
	double weight = 1.;
	if (!IsDATA) {
		weight *= weight_norm_1invpb*ev.GetTriggerLumi("Full");
		weight *= ev.MCweight();
		weight *= w_prefire;
	}

	// signal flag
	bool isDIMU = isSignal(DIMU);
	bool isDIELEC = isSignal(DIELEC);
	bool isEMU = isSignal(EMU);
	bool isMUJJ = isSignal(MUJJ);
	bool isEJJ = isSignal(EJJ);

	// draw histograms
	if (isDIMU) {
		path = "dimu/";
		DrawHists(path, muons, weight, "NOCUT");
		DrawHists(path, muons_veto, weight, "POGLooseWithTightIso");
		DrawHists(path, muons_tight, weight, "POGTightWithTightIso");
		DrawHists(path, electrons, weight, "NOCUT");
		DrawHists(path, electrons_veto, weight, "passLooseID");
		DrawHists(path, electrons_tight, weight, "passTightID");
		DrawHists(path, jets, weight, JET, "NOCUT");
		DrawHists(path, jets_tight, weight, JET, "tight");
		DrawHists(path, jets_dR04, weight, JET, "tight_dR04" );
		DrawHists(path, bjets_tight, weight, BJET, "tight");
		DrawHists(path, bjets_dR04, weight, BJET, "tight_dR04");
		DrawHists(path, METv, weight);
		//PrintGens();
		//vector<Gen_t> Gens = RecordHistory();
		//vector<vector<unsigned int>> history = ReturnHistory(Gens);
		PrintHistory();
	}

	if (isDIELEC) {
		path = "dielec/";
		DrawHists(path, muons, weight, "NOCUT");
        DrawHists(path, muons_veto, weight, "POGLooseWithTightIso");
        DrawHists(path, muons_tight, weight, "POGTightWithTightIso");
        DrawHists(path, electrons, weight, "NOCUT");
        DrawHists(path, electrons_veto, weight, "passLooseID");
        DrawHists(path, electrons_tight, weight, "passTightID");
        DrawHists(path, jets, weight, JET, "NOCUT");
        DrawHists(path, jets_tight, weight, JET, "tight");
        DrawHists(path, jets_dR04, weight, JET, "tight_dR04" );
        DrawHists(path, bjets_tight, weight, BJET, "tight");
        DrawHists(path, bjets_dR04, weight, BJET, "tight_dR04");
        DrawHists(path, METv, weight);
	}

	if (isEMU) {
		path = "emu/";
		DrawHists(path, muons, weight, "NOCUT");
        DrawHists(path, muons_veto, weight, "POGLooseWithTightIso");
        DrawHists(path, muons_tight, weight, "POGTightWithTightIso");
        DrawHists(path, electrons, weight, "NOCUT");
        DrawHists(path, electrons_veto, weight, "passLooseID");
        DrawHists(path, electrons_tight, weight, "passTightID");
        DrawHists(path, jets, weight, JET, "NOCUT");
        DrawHists(path, jets_tight, weight, JET, "tight");
        DrawHists(path, jets_dR04, weight, JET, "tight_dR04" );
        DrawHists(path, bjets_tight, weight, BJET, "tight");
        DrawHists(path, bjets_dR04, weight, BJET, "tight_dR04");
        DrawHists(path, METv, weight);
	}

	if (isMUJJ) {
		path = "mujj/";
		DrawHists(path, muons, weight, "NOCUT");
        DrawHists(path, muons_veto, weight, "POGLooseWithTightIso");
        DrawHists(path, muons_tight, weight, "POGTightWithTightIso");
        DrawHists(path, electrons, weight, "NOCUT");
        DrawHists(path, electrons_veto, weight, "passLooseID");
        DrawHists(path, electrons_tight, weight, "passTightID");
        DrawHists(path, jets, weight, JET, "NOCUT");
        DrawHists(path, jets_tight, weight, JET, "tight");
        DrawHists(path, jets_dR04, weight, JET, "tight_dR04" );
        DrawHists(path, bjets_tight, weight, BJET, "tight");
        DrawHists(path, bjets_dR04, weight, BJET, "tight_dR04");
        DrawHists(path, METv, weight);
	}

	if (isEJJ) {
		path = "ejj/";
		DrawHists(path, muons, weight, "NOCUT");
        DrawHists(path, muons_veto, weight, "POGLooseWithTightIso");
        DrawHists(path, muons_tight, weight, "POGTightWithTightIso");
        DrawHists(path, electrons, weight, "NOCUT");
        DrawHists(path, electrons_veto, weight, "passLooseID");
        DrawHists(path, electrons_tight, weight, "passTightID");
        DrawHists(path, jets, weight, JET, "NOCUT");
        DrawHists(path, jets_tight, weight, JET, "tight");
        DrawHists(path, jets_dR04, weight, JET, "tight_dR04" );
        DrawHists(path, bjets_tight, weight, BJET, "tight");
        DrawHists(path, bjets_dR04, weight, BJET, "tight_dR04");
        DrawHists(path, METv, weight);
	}	

}

void LeptonType::ClearPreColl() {
	muons.clear(); electrons.clear(); jets.clear(); gens.clear();
}
void LeptonType::ClearSubColl() {
	muons_tight.clear(); muons_veto.clear();
	electrons_tight.clear(); electrons_veto.clear();
	jets_tight.clear(); jets_dR04.clear(), bjets_tight.clear(); bjets_dR04.clear();
}

bool LeptonType::isSignal(const Signal& sig) {
	TString path;
	switch (sig) {
	case DIMU:
		path = "dimu/cutflow";
		FillHist(path, 4., 1., 9, 0., 9.);

		if (! (muons_tight.size() == 2)) return false;
		FillHist(path, 5., 1., 9, 0., 9.);
		if (! (muons_tight.at(0).Charge()*muons_tight.at(1).Charge() < 0)) return false;
		FillHist(path, 6., 1., 9, 0., 9.);
		
		return true;
		break;

	case DIELEC:
		path = "dielec/cutflow";
		FillHist(path, 4., 1., 9, 0., 9.);

		if (! (electrons_tight.size() == 2)) return false;
		FillHist(path, 5., 1., 9, 0., 9.);
		if (! (electrons_tight.at(0).Charge()*electrons_tight.at(1).Charge() < 0)) return false;
		FillHist(path, 6., 1., 9, 0., 9.);

		return true;
		break;

	case EMU:
		path = "emu/cutflow";
		FillHist(path, 4., 1., 9, 0., 9.);

		if (! (muons_tight.size()==1 && electrons_tight.size()==1)) return false;
		FillHist(path, 5., 1., 9, 0., 9.);
		if (! (muons_tight.at(0).Charge()*electrons_tight.at(0).Charge() < 0)) return false;
		FillHist(path, 6., 1., 9, 0., 9.);

		return true;
		break;

	case MUJJ:
		path = "mujj/cutflow";
		FillHist(path, 4., 1., 9, 0., 9.);

		if (jets_dR04.size() < 4) return false;
		if (bjets_dR04.size() < 1) return false;
		FillHist(path, 5., 1., 9, 0., 9.);
		if (! (muons_tight.size()==1)) return false;
		FillHist(path, 6., 1., 9, 0., 9.);

		return true;
		break;

	case EJJ:
		path = "ejj/cutflow";
		FillHist(path, 4., 1., 9, 0., 9.);

		if (jets_dR04.size() < 4) return false;
		if (bjets_dR04.size() < 1) return false;
		FillHist(path, 5., 1., 9, 0., 9.);
		if (! (electrons_tight.size()==1)) return false;
		FillHist(path, 6., 1., 9, 0., 9.);

		return true;
		break;

	default:
		cerr << "[LeptonType::isSignal] signal = " << sig << endl;
		cerr << "[LeptonType::isSignal] wrong signal region" << endl;
		exit(EXIT_FAILURE);
	}
}
bool LeptonType::isControl(const Control& con) {
	switch (con) {
	default:
		cerr << "[LeptonType::isControl] control = " << con << endl;
		cerr << "[LeptonType::isControl wrong control region" << endl;
		exit(EXIT_FAILURE);
	}
}

void LeptonType::DrawHists(TString path, const vector<Muon>& muons, const double& weight, const TString& muonid) {
	TString this_path = path + "muons/" +  muonid + "/";;
	FillHist(this_path + "nMuons", muons.size(), weight, 10,  0., 10.);
	
	TString obj_path;
	for (unsigned int i = 0; i < muons.size(); i++) {
		obj_path = this_path + TString::Itoa(i+1, 10) + "/";
		FillHist(obj_path + "pt", muons.at(i).Pt(), weight, 300, 0., 300.);
		FillHist(obj_path + "eta", muons.at(i).Eta(), weight, 48, -2.4, 2.4);
		FillHist(obj_path + "phi", muons.at(i).Phi(), weight, 70, -3.5, 3.5);
	}
}

void LeptonType::DrawHists(TString path, const vector<Electron>& electrons, const double& weight, const TString& elecid) {
	TString this_path = path + "electrons/" + elecid + "/";
	FillHist(this_path + "nElectrons", electrons.size(), weight, 10, 0., 10.);

	TString obj_path;
	for (unsigned int i = 0; i < electrons.size(); i++) {
		obj_path = this_path + TString::Itoa(i+1, 10) + "/";
		FillHist(obj_path + "pt", electrons.at(i).Pt(), weight, 300, 0., 300.);
		FillHist(obj_path + "eta", electrons.at(i).Eta(), weight, 50, -2.5, 2.5);
		FillHist(obj_path + "phi", electrons.at(i).Phi(), weight, 70, -3.5, 3.5);
	}
}

void LeptonType::DrawHists(TString path, const vector<Jet>& jets, const double& weight, const JetType& type, const TString& jetid) {
	TString this_path;
	if (type==JET) this_path = path + "jets/";
	else if (type==BJET) this_path = path + "bjets/";
	else {
		cerr << "[LeptonType::DrawHists] Wrong type for jets" << endl;
		exit(EXIT_FAILURE);
	}

	this_path += jetid + "/";
	FillHist(this_path + "size", jets.size(), weight, 14, 0., 14);

	TString obj_path;
	for (unsigned int i = 0; i < jets.size(); i++) {
		obj_path = this_path + TString::Itoa(i+1, 10) + "/";
		FillHist(obj_path + "pt", jets.at(i).Pt(), weight, 300, 0., 300.);
		FillHist(obj_path + "eta", jets.at(i).Eta(), weight, 48, -2.4, 2.4);
		FillHist(obj_path + "phi", jets.at(i).Phi(), weight, 70, -3.5, 3.5);
	}
}

void LeptonType::DrawHists(TString path, const Particle& METv, const double& weight) {
	TString this_path = path + "METv/";
	FillHist(this_path + "pt", METv.Pt(), weight, 300, 0., 300.);
	FillHist(this_path + "phi", METv.Phi(), weight, 70, -3.5, 3.5);
}

// test function
void LeptonType::PrintGens() {
	cout << "[LeptonType::PrintGens] nEvent = " << nEvent << endl;

	for (auto& gen: gens) {
		gen.Print();
	}
}

vector<LeptonType::Gen_t> LeptonType::RecordHistory() {
	// container for history;
	vector<Gen_t> out;
	for (const auto& gen: gens) {
		Gen_t this_gen;
		this_gen.me = &gen;
		this_gen.mother_idx = gen.MotherIndex();
		for (const auto& dau: gens) {
			if (dau.MotherIndex() == gen.Index()) this_gen.daughters_idx.push_back(dau.Index());
		}
		out.push_back(this_gen);
	}
	
	if (gens.size() != out.size()) {
		cerr << "[LeptonType::RecordHistory] WARNING:: missed history" << endl;
	}

	return out;
}

vector<vector<unsigned int>> LeptonType::ReturnHistory(const vector<LeptonType::Gen_t>& Gens) const {
	unsigned int my_idx = 9999;
	unsigned int mother_idx = 9999;
	
	// start with gens with no daughter;
	vector<LeptonType::Gen_t> status1;
	for (const auto& this_gen: Gens) {
		if (this_gen.daughters_idx.size() == 0) status1.push_back(this_gen);
		else continue;
	}
	
	// record history for every status1 particles
	unsigned int endpoints = status1.size();
	vector<vector<unsigned int>> history;

	// record history until index become proton
	unsigned int n_history = 0;
	for (const auto& endpoint: status1) {
		my_idx = endpoint.me->Index();
		mother_idx = endpoint.mother_idx;
		vector<unsigned int> this_history;
		this_history.push_back(my_idx);
		while (my_idx != 0 || my_idx != 1) {
			Gen_t mother = Gens.at(mother_idx);
			my_idx = mother.me->Index();
			mother_idx = mother.mother_idx;
			this_history.push_back(my_idx);
		}
		history.push_back(this_history);
		n_history++;

		if (n_history > endpoints) {
			cerr << "[LeptonType::ReturnHistory] WARNING::Too many histories...sth wrong" << endl;
		}
	}

	return history;
}


void LeptonType::PrintHistory(){
	// make a history first
	int my_idx = -999;
	int mother_idx = -999;

	vector<Gen> endpoints;
	for (const auto& gen: gens) {
		if (gen.Status() == 1) endpoints.push_back(gen);
		else continue;
	}
	
	// record history
	vector<vector<unsigned int>> history;
	for (const auto& endpoint: endpoints) {
		vector<unsigned int> this_history;

		// initialize index
		my_idx = endpoint.Index();
		mother_idx = endpoint.MotherIndex();

		while(mother_idx != -1) {
			this_history.push_back(my_idx);
			
			// go up one generation
			Gen mother = gens.at(mother_idx);
			my_idx = mother.Index();
			mother_idx = mother.MotherIndex();
		}
		this_history.push_back(my_idx); // which is proton
		history.push_back(this_history);
	}
	
	// now print history
	unsigned int n_history = 0;
	for (const auto& this_history: history) {
		int endpoint_idx = this_history.at(0);
		const Gen& endpoint = gens.at(endpoint_idx);

		// if endpoint is not electron or muon, continue
		if (endpoint.PID() != 11 && endpoint.PID() != 13) continue;
		
		// print history
		cout << "[LeptonType::PrintHistory] n_history = " << n_history << endl;
		for (const auto& idx: this_history) {
			gens.at(idx).Print();
		}
		n_history++;
	}
}




LeptonType::LeptonType(){}

LeptonType::~LeptonType(){}
