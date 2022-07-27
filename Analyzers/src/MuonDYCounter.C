#include "MuonDYCounter.h"

void MuonDYCounter::initializeAnalyzer(){
		
}

void MuonDYCounter::executeEvent(){

		if (!PassMETFilter()) return;

		Event ev = GetEvent();
		vector<Muon> muon_coll = GetAllMuons();
		vector<Electron> ele_coll = GetAllElectrons();
		vector<Jet> jet_coll = GetAllJets();
		Particle METv = ev.GetMETVector();

		sort(muon_coll.begin(), muon_coll.end(), PtComparing);
		sort(ele_coll.begin(), ele_coll.end(), PtComparing);
		sort(jet_coll.begin(), jet_coll.end(), PtComparing);

		vector<Muon> muonT_coll = SelectMuons(muon_coll, "TrackerMuon", 10., 2.4);
		if (! (muonT_coll.size() == 2)) return;

		Muon tag, probe;
		if (muonT_coll.at(0).PassPath("HLT_IsoMu24_v")) {
				tag = muonT_coll.at(0);
				probe = muonT_coll.at(1);
		}
		else if (muonT_coll.at(1).PassPath("HLT_IsoMu24_v")) {
				tag = muonT_coll.at(1);
				probe = muonT_coll.at(0);
		}
		else return;

		if (! (tag.Pt() > 26.)) return;
		if (! (tag.RelIso() < 0.2)) return;
		if (! (tag.PassID("POGTight"))) return;
		if (! (tag.DeltaR(probe) > 0.3)) return;

		Particle pair = tag + probe;
		if (! (60 < pair.M()  && pair.M() < 140)) return;

		// now count events
		const unsigned int neta_arr = 5;
		const unsigned int npt_arr = 9;
		const unsigned int nmass_arr = (130-60)*2+1;
		double eta_bins[neta_arr] = {0., 0.9, 1.2, 2.1, 2.4};
		double pt_bins[npt_arr] = {15., 20., 25., 30., 40., 50., 60., 120., 200.};
		double mass_bins[nmass_arr] = {0};
		for (unsigned int i = 0; i < nmass_arr; i++) mass_bins[i] = 60 + i*0.5;

		FillHist("abseta_pt_zmass",
						 fabs(probe.Eta()), probe.Pt(), pair.M(), 1.,
						 neta_arr-1, eta_bins,
						 npt_arr-1, pt_bins,
						 nmass_arr, mass_bins);
						 
}

MuonDYCounter::MuonDYCounter(){

}

MuonDYCounter::~MuonDYCounter(){

}


