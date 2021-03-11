#include "Selector.h"

double Selector::getWeight(
		const TString channel, Event& ev, 
		vector<Muon>& muons, vector<Electron>& electrons, 
		vector<Jet>& jets) {
	if (IsDATA)
		return 1.;

	double out = 1.;
	
	// genWeights
	const double w_prefire = GetPrefireWeight(0);
	const double w_gen = ev.MCweight()*weight_norm_1invpb;
	const double w_lumi = ev.GetTriggerLumi("Full");
	const double w_pileup = GetPileUpWeight(nPileUp, 0);
	//cout << "w_prefire: " <<  w_prefire << endl;
    //cout << "w_gen: " << w_gen << endl;
    //cout << "w_lumi: " << w_lumi << endl;
    //cout << "w_pileup: " << w_pileup << endl;
    out *= w_prefire*w_gen*w_lumi*w_pileup;

	double w_idsf = 1.;
	const TString ID = "HcToWATight";
	double w_trigsf = 1.;
	//cout << "channel: " << channel << endl;
	if (channel.Contains("emu")) {
        const Muon& mu = muons.at(0);
        const Electron& ele = electrons.at(0);
        const double mu_idsf = mcCorr->MuonID_SF(ID, mu.Eta(), mu.MiniAODPt(), 0);
        const double ele_idsf = mcCorr->ElectronID_SF(ID, ele.scEta(), ele.Pt(), 0);
        //cout << "mu_idsf: " << mu_idsf << endl;
        //cout << "ele_idsf: " << ele_idsf << endl;
        w_idsf *= mu_idsf*ele_idsf;
        w_trigsf = 
			mcCorr->GetTriggerSF(electrons, muons, "EMuIso_HNTopID", "");
    }
    if (channel.Contains("3mu")) {
        const Muon& mu1 = muons.at(0);
        const Muon& mu2 = muons.at(1);
        const Muon& mu3 = muons.at(2);
        const double mu1_idsf = mcCorr->MuonID_SF(ID, mu1.Eta(), mu1.MiniAODPt(), 0);
        const double mu2_idsf = mcCorr->MuonID_SF(ID, mu2.Eta(), mu2.MiniAODPt(), 0);
        const double mu3_idsf = mcCorr->MuonID_SF(ID, mu3.Eta(), mu3.MiniAODPt(), 0);
        w_idsf *= mu1_idsf*mu2_idsf*mu3_idsf;
        w_trigsf = 
			mcCorr->GetTriggerSF(electrons, muons, "DiMuIso_HNTopID", "");
        //cout << "id sf: " << w_idsf << endl;
        //cout << "trig sf: " << w_trigsf << endl;
    }
	if (channel.Contains("1e2mu")) {
        const Muon& mu1 = muons.at(0);
        const Muon& mu2 = muons.at(1);
        const Electron& ele = electrons.at(0);
        const double mu1_idsf = mcCorr->MuonID_SF(ID, mu1.Eta(), mu1.MiniAODPt(), 0);
        const double mu2_idsf = mcCorr->MuonID_SF(ID, mu2.Eta(), mu2.MiniAODPt(), 0);
        const double ele_idsf = mcCorr->ElectronID_SF(ID, ele.scEta(), ele.Pt(), 0);
        w_idsf *= mu1_idsf*mu2_idsf*ele_idsf;
        if (EMuTrigOnly)
            w_trigsf
                = mcCorr->GetTriggerSF(electrons, muons, "EMuIso_HNTopID", "");
        else {
            if (ev.PassTrigger(trigs_emu))
                w_trigsf
                    = mcCorr->GetTriggerSF(electrons, muons, "EMuIso_HNTopID", "");
            else // then it pass dblmuon trigger, see Selector::Selector
                w_trigsf
                    = mcCorr->GetTriggerSF(electrons, muons, "DiMuIso_HNTopID", "");
        }
        //cout << "id sf: " << w_idsf << endl;
        //cout << "trig sf: " << w_trigsf << endl;
    }
    out *= w_idsf*w_trigsf;

	// b-tagging SF
    double w_btag = 1.;
    if (RunDeepCSV) {
        JetTagging::Parameters jtp_DeepCSV_Medium
            = JetTagging::Parameters(JetTagging::DeepCSV, JetTagging::Medium, JetTagging::incl, JetTagging::mujets);
        w_btag = mcCorr->GetBTaggingReweight_1a(jets, jtp_DeepCSV_Medium);
    }
    else {
        JetTagging::Parameters jtp_DeepJet_Medium
            = JetTagging::Parameters(JetTagging::DeepJet, JetTagging::Medium, JetTagging::incl, JetTagging::mujets);
        w_btag = mcCorr->GetBTaggingReweight_1a(jets, jtp_DeepJet_Medium);
    }
    //cout << "w_btag: " << w_btag << endl;
    out *= w_btag;
	
	return out;
}
	
