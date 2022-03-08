#include "Selector.h"

Selector::Selector() {}
Selector::~Selector() {
		outfile->cd();
		tree->Write();
}

void Selector::initializeAnalyzer(){
		// ArgumentGaurds
		if (DataYear != 2017) {
				cerr << "DataYear " << DataYear << " is not supported" << endl;
				exit(EXIT_FAILURE);
		}
		// flags
		Skim1E2Mu = HasFlag("Skim1E2Mu");
		Skim3Mu = HasFlag("Skim3Mu");
		// IDs
		SetMuonIDs("HcToWATight", "HcToWALoose");
		SetElectronIDs("HcToWATight", "HcToWALoose");

		// Triggers
		vector<TString> emuTrigs = {
				"HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v",
				"HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v"
		};
		vector<TString> dblmuTrigs = {
				"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v",
				"HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v",
				"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v"
		};
		SetTriggers(emuTrigs, dblmuTrigs);

		// Jet tagger
		vector<JetTagging::Parameters> jtps;
		jtps.emplace_back(
						JetTagging::Parameters(JetTagging::DeepJet, JetTagging::Medium, JetTagging::incl, JetTagging::mujets));
		mcCorr->SetJetTaggingParameters(jtps);

		// Link branches to the tree
		LinkTreeContents();
}

void Selector::executeEvent(){

		if (! PassMETFilter()) return;
		Event ev = GetEvent();
		const Particle METv = ev.GetMETVector();

		// Object Definitions
		vector<Gen>					truth_coll			= GetGens();
		vector<Muon>				muon_coll				= GetAllMuons();
		vector<Electron>		electron_coll		= GetAllElectrons();
		vector<Jet>					jet_coll				= GetAllJets();

		// to avoid any confusion, sort in the first time
		sort(muon_coll.begin(), muon_coll.end(), PtComparing);
		sort(electron_coll.begin(), electron_coll.end(), PtComparing);
		sort(jet_coll.begin(), jet_coll.end(), PtComparing);

		// muons -> loose muons, electrons -> loose electrons, jets -> lepton cleaned tight jets
		vector<Muon>				muons = SelectMuons(muon_coll, LooseMuonID, 10., 2.4);
		vector<Electron>		electrons = SelectElectrons(electron_coll, LooseElecID, 10., 2.5);
		vector<Jet>					jets = SelectJets(jet_coll, "tight", 20., 2.4);
		jets = JetsVetoLeptonInside(jets, electrons, muons, 0.4);

		// Channel definition
		if (Skim3Mu) { if (! (muons.size() == 3 && electrons.size() == 0)) return; }
		else if (Skim1E2Mu) { if (! (muons.size() == 2 && electrons.size() == 1)) return; }
		else {
				cerr << "Wrong flag" << endl;
				exit(EXIT_FAILURE);
		}

		// Initialize Contents
		// events
		passDblMuTrigs = ev.PassTrigger(trigs_dblmu);
		passEMuTrigs = ev.PassTrigger(trigs_emu);
		METv_pt = METv.Pt(); METv_phi = METv.Phi();
		genWeight = 1.; trigLumi = 1.;
		if (!IsDATA) {
				genWeight = ev.MCweight()*weight_norm_1invpb;
				trigLumi = ev.GetTriggerLumi("Full");
		}

		// muons
		nMuons = muons.size();
		for (unsigned int i = 0; i < nMuons; i++) {
				const Muon &mu = muons.at(i);
				muons_pt[i] = mu.Pt(); 
				muons_eta[i] = mu.Eta(); 
				muons_phi[i] = mu.Phi();
				muons_mass[i] = mu.M();
				muons_miniIso[i] = mu.MiniRelIso();
				muons_charge[i] = mu.Charge();
				muons_lepType[i] = GetLeptonType(mu, truth_coll);
				muons_isTight[i] = mu.PassID(TightMuonID);
		}
		// electrons
		nElectrons = electrons.size();
		for (unsigned int i = 0; i < nElectrons; i++) {
				const Electron &ele = electrons.at(i);
				electrons_pt[i] = ele.Pt();
				electrons_eta[i] = ele.Eta();
				electrons_phi[i] = ele.Phi();
				electrons_mass[i] = ele.M();
				electrons_miniIso[i] = ele.MiniRelIso();
				electrons_charge[i] = ele.Charge();
				electrons_lepType[i] = GetLeptonType(ele, truth_coll);
				electrons_isTight[i] = ele.PassID(TightElecID);
		}
		// jets
		nJets = jets.size();
		const float bcut = mcCorr->GetJetTaggingCutValue(JetTagging::DeepJet, JetTagging::Medium);
		for (unsigned int i = 0; i < nJets; i++) {
				const Jet& j = jets.at(i);
				jets_pt[i] = j.Pt();
				jets_eta[i] = j.Eta();
				jets_phi[i] = j.Phi();
				jets_mass[i] = j.M();
				jets_btagScore[i] = j.GetTaggerResult(JetTagging::DeepJet);
				jets_isBtagged[i] = j.GetTaggerResult(JetTagging::DeepJet) > bcut;
				jets_cH[i] = j.ChargedHadronEnergyFraction();
				jets_nH[i] = j.NeutralHadronEnergyFraction();
				jets_nEM[i] = j.NeutralEmEnergyFraction();
				jets_cEM[i] = j.ChargedEmEnergyFraction();
				jets_muE[i] = j.MuonEnergyFraction();
				jets_cM[i] = j.ChargedMultiplicity();
				jets_nM[i] = j.NeutralMultiplicity();
		}
		tree->Fill();
}

void Selector::SetMuonIDs(TString TightID, TString LooseID) {
		TightMuonID = TightID;
		LooseMuonID = LooseID;
}
void Selector::SetElectronIDs(TString TightID, TString LooseID) {
		TightElecID = TightID;
		LooseElecID = LooseID;
}
void Selector::SetTriggers(vector<TString> emuTrigs, vector<TString> dblmuTrigs) {
		trigs_dblmu.clear(); trigs_emu.clear();
		for (const auto trig: emuTrigs)
				trigs_emu.emplace_back(trig);
		for (const auto trig: dblmuTrigs)
				trigs_dblmu.emplace_back(trig);
}
void Selector::LinkTreeContents() {
		tree = new TTree("Events", "Events");
		// events
		tree->Branch("IsDATA", &IsDATA);
		tree->Branch("DataStream", &DataStream);
		tree->Branch("run", &run); tree->Branch("event", &event); tree->Branch("lumi", &lumi);
		tree->Branch("nPV", &nPV); tree->Branch("nPileUp", &nPileUp);
		tree->Branch("passDblMuTrigs", &passDblMuTrigs); tree->Branch("passEMuTrigs", &passEMuTrigs);
		tree->Branch("trigLumi", &trigLumi);
		tree->Branch("METv_pt", &METv_pt);
		tree->Branch("METv_phi", &METv_phi);
		tree->Branch("genWeight", &genWeight);
		// muons
		tree->Branch("nMuons", &nMuons);
		tree->Branch("muons_pt", muons_pt, "muons_pt[nMuons]/F");
		tree->Branch("muons_eta", muons_eta, "muons_eta[nMuons]/F");
		tree->Branch("muons_phi", muons_phi, "muons_phi[nMuons]/F");
		tree->Branch("muons_mass", muons_mass, "muons_mass[nMuons]/F");
		tree->Branch("muons_miniIso", muons_miniIso, "muons_miniIso[nMuons]/F");
		tree->Branch("muons_charge", muons_charge, "muons_charge[nMuons]/I");
		tree->Branch("muons_lepType", muons_lepType, "muons_lepType[nMuons]/I");
		tree->Branch("muons_isTight", muons_isTight, "muons_isTight[nMuons]/O");
		// electrons
		tree->Branch("nElectrons", &nElectrons);
		tree->Branch("electrons_pt", electrons_pt, "electrons_pt[nElectrons]/F");
		tree->Branch("electrons_eta", electrons_eta, "electrons_eta[nElectrons]/F");
		tree->Branch("electrons_phi", electrons_phi, "electrons_phi[nElectrons]/F");
		tree->Branch("electrons_mass", electrons_mass, "electrons_mass[nElectrons]/F");
		tree->Branch("electrons_miniIso", electrons_miniIso, "electrons_miniIso[nElectrons]/F");
		tree->Branch("electrons_charge", electrons_charge, "electrons_charge[nElectrons]/I");
		tree->Branch("electrons_lepType", electrons_lepType, "electrons_lepType[nElectrons]/I");
		tree->Branch("electrons_isTight", electrons_isTight, "electrons_isTight[nElectrons]/O");
		// jets
		tree->Branch("nJets", &nJets);
		tree->Branch("jets_pt", jets_pt, "jets_pt[nJets]/F");
		tree->Branch("jets_eta", jets_eta, "jets_eta[nJets]/F");
		tree->Branch("jets_phi", jets_phi, "jets_phi[nJets]/F");
		tree->Branch("jets_mass", jets_mass, "jets_mass[nJets]/F");
		tree->Branch("jets_btagScore", jets_btagScore, "jets_btagScore[nJets]/F");	
		tree->Branch("jets_isBtagged", jets_isBtagged, "jets_isBtagged[nJets]/O");
		tree->Branch("jets_cH", jets_cH, "jets_cH[nJets]/F");		 // charged hadron energy fraction
		tree->Branch("jets_nH", jets_nH, "jets_nH[nJets]/F");		 // neutral hadron energy fraction
		tree->Branch("jets_cEM", jets_cEM, "jets_cEM[nJets]/F"); // charged EM energy fraction
		tree->Branch("jets_nEM", jets_nEM, "jets_nEM[nJets]/F"); // neutral EM energy fraction
		tree->Branch("jets_muE", jets_muE, "jets_muE[nJets]/F"); // muon energy fraction
		tree->Branch("jets_cM", jets_cM, "jets_cM[nJets]/F");    // charged multiplicity
		tree->Branch("jets_nM", jets_nM, "jets_nM[nJets]/F");    // neutral multiplicity
}




























