#include "MCCorrection.h"
double MCCorrection::GetPUVetoSF(const vector<Jet>& jets, const TString wp) {
	if (IsDATA)
		return 1.;
	
	double prob_mc = 1., prob_data = 1.;
	TH2F* hist_eff = nullptr;
	TH2F* hist_sf = nullptr;
	TString year = TString::Itoa(DataYear, 10);
	TString key_eff, key_sf;
	if (wp == "tight") {
		key_eff = "PUVeto_eff_mc" + year + "_T";
		key_sf = "PUVeto_eff_sf" + year + "_T";
	}
	else if (wp == "medium") {
		key_eff = "PUVeto_eff_mc" + year + "_M";
		key_sf = "PUVeto_eff_sf" + year + "_M";
	}
	else if (wp == "loose") {
		key_eff = "PUVeto_eff_mc" + year + "_L";
		key_sf = "PUVeto_eff_sf" + year + "_L";
	}
	else {
		std::cerr << "[MCCorrection::GetPUVetoSF] No wp " << wp << std::endl;
		exit(EXIT_FAILURE);
	}
	hist_eff = map_hist_puveto[key_eff];
	hist_sf = map_hist_puveto[key_sf];
	if (!hist_eff) {
		std::cerr << "[MCCorrection::GetPUVetoSF] No hist " << key_eff << endl;
		exit(EXIT_FAILURE);
	}
	if (!hist_sf) {
		std::cerr << "[MCCorrection::GetPUVetoSF] No hist " << key_sf << endl;
		exit(EXIT_FAILURE);
	}
	for (const auto& j : jets) {
		const double this_pt = j.Pt();
		const double this_eta = j.Eta();

		// only consider pt < 50 GeV, |eta| < 5 jets
		if (this_pt > 50. || fabs(this_eta) > 5.)
			continue;
		
		const double eff = hist_eff->GetBinContent(hist_eff->FindBin(this_pt, this_eta));
		const double sf = hist_sf->GetBinContent(hist_sf->FindBin(this_pt, this_eta));
		const bool pass_PUVeto = j.PassPileupMVA(wp);
		if (pass_PUVeto) {
			prob_mc *= eff;
			prob_data *= sf*eff;
		}
		else {
			prob_mc *= (1.-eff);
			prob_data *= (1.-sf*eff);
		}
	}
	if (prob_mc != 0)
		return prob_data/prob_mc;
	else
		return 0;
}
