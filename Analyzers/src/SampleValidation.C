#include "SampleValidation.h"

SampleValidation::SampleValidation() {}
SampleValidation::~SampleValidation() {}

void SampleValidation::initializeAnalyzer() {

		// B-tagging
    vector<JetTagging::Parameters> jtps;
    jtps.emplace_back(
        JetTagging::Parameters(JetTagging::DeepJet, JetTagging::Medium, JetTagging::incl, JetTagging::mujets));
    mcCorr->SetJetTaggingParameters(jtps);

}

void SampleValidation::executeEvent() {

		if (!PassMETFilter()) return;
		
		Event ev = GetEvent();
		vector<Gen> truth_coll = GetGens();
		vector<Muon> muon_coll = GetAllMuons();
		vector<Electron> ele_coll = GetAllElectrons();
		vector<Jet> jet_coll = GetAllJets();
		Particle METv = ev.GetMETVector();

		// sort objects
		sort(muon_coll.begin(), muon_coll.end(), PtComparing);
		sort(ele_coll.begin(), ele_coll.end(), PtComparing);
		sort(jet_coll.begin(), jet_coll.end(), PtComparing);

		// select objects
		vector<Muon> muonA_coll, muonW_coll, muonWs_coll, muonF_coll;		// from A, from W, from W*, fake
		vector<Electron> eleW_coll, eleWs_coll, eleF_coll;							// from W, from W*, fake
		for (const auto &mu: muon_coll) {
				const int LepType = GetLeptonType(mu, truth_coll);
				if (LepType == 1) muonW_coll.emplace_back(mu);
				else if (LepType == 2) muonA_coll.emplace_back(mu);
				else if (LepType == 6) muonWs_coll.emplace_back(mu);
				else if (LepType < 0) muonF_coll.emplace_back(mu);
				else continue;
		}
		for (const auto &ele: ele_coll) {
				const int LepType = GetLeptonType(ele, truth_coll);
				if (LepType == 1) eleW_coll.emplace_back(ele);
				else if (LepType == 6) eleWs_coll.emplace_back(ele);
				else if (LepType < 0) eleF_coll.emplace_back(ele);
				else continue;
		}
		
		vector<Muon> muonP_coll;
		muonP_coll.insert(muonP_coll.end(), muonA_coll.begin(), muonA_coll.end());
		muonP_coll.insert(muonP_coll.end(), muonW_coll.begin(), muonW_coll.end());
		muonP_coll.insert(muonP_coll.end(), muonWs_coll.begin(), muonWs_coll.end());

		vector<Electron> eleP_coll;
		eleP_coll.insert(eleP_coll.end(), eleW_coll.begin(), eleW_coll.end());
		eleP_coll.insert(eleP_coll.end(), eleWs_coll.begin(), eleWs_coll.end());

		vector<Jet> jetT_coll, bjet_coll;
		for (const auto &j: jet_coll) {
				if (! (fabs(j.Eta()) < 2.4)) continue;
				if (! (j.Pt() > 10.)) continue;
				if (! j.PassID("tight")) continue;

				// remove overlap
				bool isLeptonMatched = false;
				for (const auto &mu: muonP_coll) {
						if (! (j.DeltaR(mu) > 0.4)) {
								isLeptonMatched = true;
								break;
						}
				}
				for (const auto &ele: eleP_coll) {
						if (! (j.DeltaR(ele) > 0.4)) {
								isLeptonMatched = true;
								break;
						}
				}
				if (isLeptonMatched) continue;
				jetT_coll.emplace_back(j);
		}
		// b-tagging
		for (const auto &j: jetT_coll) {
				const double score = j.GetTaggerResult(JetTagging::DeepJet);
				if (score > mcCorr->GetJetTaggingCutValue(JetTagging::DeepJet, JetTagging::Medium))
						bjet_coll.emplace_back(j);
		}

		// baseline selection
		if (! (muonA_coll.size() == 2)) return;

		TString channel = "";
		if (muonP_coll.size() == 3 && eleP_coll.size() == 0) channel = "3mu";
		else if (muonP_coll.size() == 2 && eleP_coll.size() == 1) channel = "1e2mu";
		else return;
		
		double weight = GetPrefireWeight(0)*ev.GetTriggerLumi("Full")*MCweight();
		FillMuons(channel+"/muonA", muonA_coll, weight, false);
		FillMuons(channel+"/muonW", muonW_coll, weight, false);
		FillMuons(channel+"/muonWs", muonWs_coll, weight, false);
		FillMuons(channel+"/muonF", muonF_coll, weight, false);
		FillElectrons(channel+"/eleW", eleW_coll, weight, false);
		FillElectrons(channel+"/eleWs", eleWs_coll, weight, false);
		FillJets(channel+"/jets", jetT_coll, weight);
		FillJets(channel+"/bjets", bjet_coll, weight);
		FillObject(channel+"/METv", METv, weight);
		if (muonA_coll.size() == 2) {
				Particle Acand = muonA_coll.at(0) + muonA_coll.at(1);
				FillObject(channel+"/Acand", Acand, weight);
		}

		return;
}


