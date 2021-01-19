#include "HistoMaker.h"

//==== constructor and destructor
HistoMaker::HistoMaker() {};
HistoMaker::HistoMaker(const TString &region, const vector<TString> &cuts):_region(region), _cuts(cuts) {};
HistoMaker::~HistoMaker() {};

//==== public functions
//==== TODO: update get and set functions
void HistoMaker::SetRegion(const TString &region) { _region = region; }
TString HistoMaker::GetRegion() const { return _region; }

void HistoMaker::FillCutflow(const TString &cut) {
	bool isFilled = false;
    const unsigned int ncuts = _cuts.size();
    for (unsigned int i = 0; i < ncuts; i++) {
        if (cut == _cuts.at(i)) {
            FillHist("cutflow/" + _region, static_cast<double>(i), 1., ncuts, 0., static_cast<double>(ncuts)); 
            isFilled = true;
            break;
        }
    }   
    
    if (!isFilled) {
        cerr << "[CutflowMaker::fillCutflow] no cut named " << cut << endl;
        exit(EXIT_FAILURE);
    }   
}

void HistoMaker::FillObjects(const TString path, const vector<Muon> &muons, const double &weight) {
	TString obj_path;
	FillHist(_region + "/" + path + "size", muons.size(), weight, 20, 0., 20.);
	for (unsigned int i = 0; i < muons.size(); i++) {
		obj_path = _region + "/" + path + "/" + TString::Itoa(i+1, 10) + "/";
		double this_pt = muons.at(i).Pt();
		if (this_pt > 300.)
			this_pt = 299.9;
		FillHist(obj_path + "pt", this_pt, weight, 300, 0., 300.);
		FillHist(obj_path + "eta", muons.at(i).Eta(), weight, 48, -2.4, 2.4);
		FillHist(obj_path + "phi", muons.at(i).Phi(), weight, 70, -3.5, 3.5);
		
		// Isolation
		double this_relIso = muons.at(i).RelIso();
		double this_trkIso = muons.at(i).TrkIso()/muons.at(i).Pt();
		double this_miniIso = muons.at(i).MiniRelIso();
		if (this_relIso > 0.5)
			this_relIso = 0.499;
		if (this_trkIso > 0.8)
			this_trkIso = 0.799;
		if (this_miniIso > 0.8)
			this_miniIso = 0.799;
		FillHist(obj_path + "relIso", this_relIso, weight, 50, 0., 0.5);
		FillHist(obj_path + "trkIso", this_trkIso, weight, 80, 0., 0.8);
		FillHist(obj_path + "miniRelIso", this_miniIso, weight, 80, 0., 0.8);

		// Impact Parameters
		double this_dZ = fabs(muons.at(i).dZ());
		double this_dXY = fabs(muons.at(i).dXY());
		double this_SIP3D = 9999.;
		if (muons.at(i).IP3Derr() != 0)
			this_SIP3D = fabs(muons.at(i).IP3D()/muons.at(i).IP3Derr());
		if (this_dZ > 0.8) this_dZ = 0.799;
		if (this_dXY > 0.5) this_dXY = 0.499;
		if (this_SIP3D > 10.) this_SIP3D = 9.99;
		FillHist(obj_path + "dZ", this_dZ, weight, 80, 0., 0.8);
		FillHist(obj_path + "dXY", this_dXY, weight, 50, 0., 0.5);
		FillHist(obj_path + "SIP3D", this_SIP3D, weight, 100, 0., 10.);
	}
}
