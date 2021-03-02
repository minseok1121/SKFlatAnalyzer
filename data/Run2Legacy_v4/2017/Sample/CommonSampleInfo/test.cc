void test() {
	for (int i = 0; i < 150; i++) {
		TString path = "/gv0/DATA/SKFlat/Run2Legacy_v4/2017/PrivateMC/HctoWA/FullSim_TTToHcToWA_AToMuMu_MHc70_MA15/2020_06_22_200736/SKFlatNtuple_2017_MC_" + TString::Itoa(i, 10) + ".root";
		TFile* f = new TFile(path);
		TTree* tr = (TTree*)f->Get("recoTree/SKFlat");
		int nTotal; tr->SetBranchAddress("nTotal", &nTotal);
		double gen_weight; tr->SetBranchAddress("gen_weight", &gen_weight);
		for (int event = 0; event < nTotal; event++) {
			tr->GetEvent(event);
			cout << gen_weight << endl;
		}
	}
}
