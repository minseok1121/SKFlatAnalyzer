#include "trigEmuCuts.h"

trigEmuCuts::trigEmuCuts() {}
trigEmuCuts::~trigEmuCuts() {}

void trigEmuCuts::initializeAnalyzer(){
		// Trigger Settings
		/*
		AllEMuTrigs = {
				"HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",
				"HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v",
				"HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_v",
				"HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v",
				"HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v",
				"HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v"
		};
		*/
		if (DataEra == "2016preVFP") {
				Ele12Trigs = {"HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v",
						          "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v"
				};
				Ele23Trigs = {"HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v",
										  "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v"
				};
		}
		else if (DataEra == "2016postVFP") {
				Ele12Trigs = {"HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v",
										  "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v"
				};
				Ele23Trigs = {"HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v",
										  "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v"
				};
		}
		else if (DataEra == "2017") {
				Ele12Trigs = {"HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v"};
				Ele23Trigs = {"HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v"};
		}
		else if (DataEra == "2018") {
				Ele12Trigs = {"HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v",
										  "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v"
				};
				Ele23Trigs = {"HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v"};
		}
		else {
				cerr << "Wrong DataEra " << DataEra << "!" << endl;
				exit(EXIT_FAILURE);
		}
		AllEMuTrigs.clear();
		for (const auto &trig: Ele12Trigs) AllEMuTrigs.emplace_back(trig);
		for (const auto &trig: Ele23Trigs) AllEMuTrigs.emplace_back(trig);
		// B-tagging
		vector<JetTagging::Parameters> jtps = {
				JetTagging::Parameters(JetTagging::DeepJet, JetTagging::Medium, JetTagging::incl, JetTagging::mujets) };
		mcCorr->SetJetTaggingParameters(jtps);
}

void trigEmuCuts::executeEvent(){

		if (!PassMETFilter()) return;
		Event								ev = GetEvent();
		vector<Muon>				muons_all = GetAllMuons();
		vector<Electron>		electrons_all = GetAllElectrons();
		vector<Jet>					jets_all = GetAllJets();
		Particle						METv = ev.GetMETVector();

		// sort objects
		sort(muons_all.begin(), muons_all.end(), PtComparing);
		sort(electrons_all.begin(), electrons_all.end(), PtComparing);
		sort(jets_all.begin(), jets_all.end(), PtComparing);

		// select objects
		vector<Muon> muons_tight = SelectMuons(muons_all, "HcToWATight", 30., 2.4);
		vector<Muon> muons_loose = SelectMuons(muons_all, "HcToWALoose", 30., 2.4);
		vector<Electron> electrons = SelectElectrons(electrons_all, "NOCUT", 10., 2.5);
		vector<Jet> jets_tight = SelectJets(jets_all, "tight", 30., 2.4);
		vector<Jet> jets = JetsVetoLeptonInside(jets_tight, electrons, muons_loose, 0.4);
		vector<Jet> bjets;
		for (const auto &jet: jets) {
				const double score = jet.GetTaggerResult(JetTagging::DeepJet);
				const double bcut 
						= mcCorr->GetJetTaggingCutValue(JetTagging::DeepJet, JetTagging::Medium);
				if (score > bcut) bjets.emplace_back(jet);
		}

		// event selection
		TString channel = selectEvent(ev, muons_tight, muons_loose, electrons, jets, bjets, METv);
		if (channel == "") return;
		
		// Now fill histograms
		double weight = 1.;
		if (channel.Contains("Ele12")) {
				FillMuon("passEle12/muon", muons_tight.at(0), weight);
				FillElectron("passEle12/electron", electrons.at(0), weight);
				FillJets("passEle12/jets", jets, weight);
				FillJets("passEle12/bjets", bjets, weight);
				FillObject("passEle12/METv", METv, weight);
		}
		if (channel.Contains("Ele23")) {
				FillMuon("passEle23/muon", muons_tight.at(0), weight);
        FillElectron("passEle23/electron", electrons.at(0), weight);
        FillJets("passEle23/jets", jets, weight);
        FillJets("passEle23/bjets", bjets, weight);
        FillObject("passEle23/METv", METv, weight);
		}
		return;
}

TString trigEmuCuts::selectEvent(
				Event &ev, vector<Muon> &muons_tight, vector<Muon> &muons_loose,
				vector<Electron> &electrons, vector<Jet> &jets, vector<Jet> &bjets,
				Particle &METv) {
		if (! (muons_tight.size() == 1 && muons_loose.size() == 1 && electrons.size() == 1))
				return "";
		Muon			&mu = muons_tight.at(0);
		Electron	&ele = electrons.at(0);
		// trigger matching
		bool matched = false;
		for (const auto &trig: AllEMuTrigs) {
				// cout << trig << endl;
				if (ele.PassPath(trig)) {
						matched = true;
						break;
				}
				else continue;
		}
		if (! matched) return "";
		if (! (mu.Charge() + ele.Charge() == 0)) return "";
		if (! (mu.DeltaR(ele) > 0.4)) return "";
		if (! (jets.size() >= 3)) return "";
		if (! (bjets.size() >= 2)) return "";
		if (! (METv.Pt() > 40.)) return "";

		// now check trigger conditions
		TString out = "";
		if (ev.PassTrigger(Ele12Trigs) && mu.Pt() > 30. && ele.Pt() > 15.)
				out += "passEle12Trigs";
		if (ev.PassTrigger(Ele23Trigs) && mu.Pt() > 30. && ele.Pt() > 25.)
				out += "passEle23Trigs";

		return out;
}

