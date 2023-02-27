#include "CR_DiLepton.h"

CR_DiLepton::CR_DiLepton(){}
CR_DiLepton::~CR_DiLepton(){}

void CR_DiLepton::initializeAnalyzer(){
    // flags
    RunSyst = HasFlag("RunSyst");

    // triggers & ID settings
    if (DataEra == "2016preVFP") {
        DblMuTriggers = {
            "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v",
            "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v",
        };
        EMuTriggers = {
            "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v",
            "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",
            "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v",
            "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v"
        };
        MuonIDs = {"HcToWATight", "HcToWALoose", "HcToWAVeto"};
        ElectronIDs = {"HcToWATight16a", "HcToWALoose16a", "HcToWAVeto16a"};
    }
    else if (DataEra == "2016postVFP") {
        DblMuTriggers = {
            "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v",
            "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v",
        };
        EMuTriggers = {
            "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v",
            "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",
            "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v",
            "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v"
        };
        MuonIDs = {"HcToWATight", "HcToWALoose", "HcToWAVeto"};
        ElectronIDs = {"HcToWATight16b", "HcToWALoose16b", "HcToWAVeto16b"};
    }
    else if (DataEra == "2017") {
        DblMuTriggers = {
            "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v",
            "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v",
            "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v"
        };
        EMuTriggers = {
            "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v",
            "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v"
        };
        MuonIDs = {"HcToWATight", "HcToWALoose", "HcToWAVeto"};
        ElectronIDs = {"HcToWATight17", "HcToWALoose17", "HcToWAVeto17"};
    }
    else if (DataEra == "2018") {
        DblMuTriggers = {
            "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v",
            "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v",
            "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v"
        };
        EMuTriggers = {
            "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v",
            "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",
            "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v"
        };
        MuonIDs = {"HcToWATight", "HcToWALoose", "HcToWAVeto"};
        ElectronIDs = {"HcToWATight18", "HcToWALoose18", "HcToWAVeto18"};
    }
    else {
        cerr << "[CR_DiLepton::initializeAnalyzer] Wrong era " << DataEra << endl;
        exit(EXIT_FAILURE);
    }

    // jet tagger
    jtps = { JetTagging::Parameters(JetTagging::DeepJet, JetTagging::Medium, JetTagging::incl, JetTagging::mujets) };
    mcCorr->SetJetTaggingParameters(jtps);

    systematics = {"Central"};
    if (RunSyst) {
        systematics.emplace_back("L1PrefireUp");
        systematics.emplace_back("L1PrefireDown");
        systematics.emplace_back("PileupReweightUp");
        systematics.emplace_back("PileupReweightDown");
        systematics.emplace_back("MuonIDSFUp");
        systematics.emplace_back("MuonIDSFDown");
        systematics.emplace_back("DblMuTrigSFUp");
        systematics.emplace_back("DblMuTrigSFDown");
    }

    // load histograms for MC corrections
    TString datapath = getenv("DATA_DIR");
    
    // muonID
    TString muonIDpath = datapath + "/" + GetEra() + "/ID/Muon";
    TFile* fMuonID = new TFile(muonIDpath+"/efficiency_TopHN_IDIso.root");
    hMuonIDSF = (TH2D*)fMuonID->Get("SF_fabs(probe_eta)_probe_pt");
    hMuonIDSF->SetDirectory(0);
    fMuonID->Close();

    // doublemuon trigger
    TFile* fMu17Leg1 = new TFile(muonIDpath+"/efficiency_Mu17Leg1_DoubleMuonTriggers.root");
    hMu17Leg1_DATA = (TH2D*)fMu17Leg1->Get("muonEffi_data_fabs(probe_eta)_probe_pt");
    hMu17Leg1_MC = (TH2D*)fMu17Leg1->Get("muonEffi_mc_fabs(probe_eta)_probe_pt");
    hMu17Leg1_DATA->SetDirectory(0);
    hMu17Leg1_MC->SetDirectory(0);
    fMu17Leg1->Close();

    TFile* fMu8Leg2 = new TFile(muonIDpath+"/efficiency_Mu8Leg2_DoubleMuonTriggers.root");
    hMu8Leg2_DATA = (TH2D*)fMu8Leg2->Get("muonEffi_data_fabs(probe_eta)_probe_pt");
    hMu8Leg2_MC = (TH2D*)fMu8Leg2->Get("muonEffi_mc_fabs(probe_eta)_probe_pt");
    hMu8Leg2_DATA->SetDirectory(0);
    hMu8Leg2_MC->SetDirectory(0);
    fMu8Leg2->Close();

}

void CR_DiLepton::executeEvent(){
    

    FillHist("DYDiMu/cutflow", 0., 1., 10, 0., 10.);
    FillHist("TTDiMu/cutflow", 0., 1., 10, 0., 10.);
    FillHist("TTEMu/cutflow", 0., 1., 10, 0., 10.);
    if (! PassMETFilter()) return;
    FillHist("DYDiMu/cutflow", 1., 1., 10, 0., 10.);
    FillHist("TTDiMu/cutflow", 1., 1., 10, 0., 10.);
    FillHist("TTEMu/cutflow", 1., 1., 10, 0., 10.); 

    // Object definition
    Event ev = GetEvent();
    Particle METv = ev.GetMETVector();
    vector<Muon> rawMuons = GetAllMuons();
    vector<Electron> rawElectrons = GetAllElectrons();
    vector<Jet> rawJets = GetAllJets();

    vector<Muon> vetoMuons = SelectMuons(rawMuons, MuonIDs.at(2), 10., 2.4);
    vector<Muon> tightMuons = SelectMuons(vetoMuons, MuonIDs.at(0), 10., 2.4);
    vector<Electron> vetoElectrons = SelectElectrons(rawElectrons, ElectronIDs.at(2), 10., 2.5);
    vector<Electron> tightElectrons = SelectElectrons(vetoElectrons, ElectronIDs.at(0), 10., 2.5);
    vector<Jet> jets = SelectJets(rawJets, "tight", 20., 2.4);
    jets = JetsVetoLeptonInside(jets, vetoElectrons, vetoMuons, 0.4);
    vector<Jet> bjets;
    for (const auto &jet: jets) {
        const double btagScore = jet.GetTaggerResult(JetTagging::DeepJet);
        const double cutValue = mcCorr->GetJetTaggingCutValue(JetTagging::DeepJet, JetTagging::Medium);
        if (btagScore > cutValue) bjets.emplace_back(jet);
    }

    std::sort(tightMuons.begin(), tightMuons.end(), PtComparing);
    std::sort(tightElectrons.begin(), tightElectrons.end(), PtComparing);
    std::sort(jets.begin(), jets.end(), PtComparing);
    std::sort(bjets.begin(), bjets.end(), PtComparing);

    // Events Selection
    TString channel = "";
    if (tightMuons.size() == 2 && vetoMuons.size() == 2 && tightElectrons.size() == 0 && vetoElectrons.size() == 0)
        channel = "isDiMu";
    else if (tightMuons.size() == 1 && vetoMuons.size() == 1 && tightElectrons.size() == 1 && vetoElectrons.size() == 1)
        channel = "isEMu";
    else
        return;
    
    // Divide channels
    if (channel == "isEMu") {
        channel = "TTEMu";
        Electron &ele = tightElectrons.at(0);
        Muon &mu = tightMuons.at(0);
        if (! ev.PassTrigger(EMuTriggers)) return;         // pass trigger
        FillHist(channel+"/cutflow", 2, 1., 10, 0., 10.);
        const bool passSafeCut = ((mu.Pt() > 25. && ele.Pt() > 15.) || (mu.Pt() > 10. && ele.Pt() > 25.));
        if (! passSafeCut) return;                                 // pass safe cut
        FillHist(channel+"/cutflow", 3, 1., 10, 0., 10.); 
        if (! (ele.Charge() + mu.Charge() == 0)) return;    // OS charge condition
        FillHist(channel+"/cutflow", 4, 1., 10, 0., 10.);
        if (! (jets.size() >= 2)) return;                   // Nj >= 2
        FillHist(channel+"/cutflow", 5, 1., 10, 0., 10.);
        if (! (bjets.size() >= 1)) return;                  // Nb >= 2
        FillHist(channel+"/cutflow", 6, 1., 10, 0., 10.);

        for (const auto &syst: systematics) {
            double weight = 1.;
            if (! IsDATA) {
                weight *= MCweight();
                weight *= ev.GetTriggerLumi("Full");
                
                if (syst == "L1PrefireUp")        weight *= GetPrefireWeight(1);
                else if (syst == "L1PrefireDown") weight *= GetPrefireWeight(-1);
                else                              weight *= GetPrefireWeight(0);
                
                if (syst == "PileupReweightUp")        weight *= GetPileUpWeight(nPileUp, 1);
                else if (syst == "PileupReweightDown") weight *= GetPileUpWeight(nPileUp, -1);
                else                                   weight *= GetPileUpWeight(nPileUp, 0);
                
                weight *= mcCorr->GetBTaggingReweight_1a(jets, jtps.at(0));

                // Lepton ID efficiency correction / Trigger eff. correction
                for (auto &mu: tightMuons) {
                    if (syst == "MuonIDSFUp")        weight *= getMuonIDSF(mu, 1);
                    else if (syst == "MuonIDSFDown") weight *= getMuonIDSF(mu, -1);
                    else                             weight *= getMuonIDSF(mu, 0);
                }
                // we don't have electron ID SF & emu trigger eff. yet
            }

            FillHist(channel+"/"+syst+"/muon/pt", mu.Pt(), weight, 300, 0., 300.);
            FillHist(channel+"/"+syst+"/muon/eta", mu.Eta(), weight, 48, -2.4, 2.4);
            FillHist(channel+"/"+syst+"/muon/phi", mu.Phi(), weight, 64, -3.2, 3.2);
            FillHist(channel+"/"+syst+"/electron/pt", ele.Pt(), weight, 300, 0., 300.);
            FillHist(channel+"/"+syst+"/electron/eta", ele.Eta(), weight, 50, -2.5, 2.5);
            FillHist(channel+"/"+syst+"/electron/phi", ele.Phi(), weight, 64, -3.2, 3.2);
            FillHist(channel+"/"+syst+"/jets/size", jets.size(), weight, 20, 0., 20.);
            FillHist(channel+"/"+syst+"/bjets/size", bjets.size(), weight, 20, 0., 20.);
            for (unsigned int i = 0; i < jets.size(); i++) {
                TString histkey = channel+"/"+syst+"/jets/"+TString::Itoa(i+1, 10);
                FillHist(histkey+"/pt", jets.at(i).Pt(), weight, 300, 0., 300.);
                FillHist(histkey+"/eta", jets.at(i).Eta(), weight, 48, -2.4, 2.4);
                FillHist(histkey+"/phi", jets.at(i).Phi(), weight, 64, -3.2, 3.2);
            }
            for (unsigned int i = 0; i < bjets.size(); i++) {
                TString histkey = channel+"/"+syst+"/bjets/"+TString::Itoa(i+1, 10);
                FillHist(histkey+"/pt", bjets.at(i).Pt(), weight, 300, 0., 300.);
                FillHist(histkey+"/eta", bjets.at(i).Eta(), weight, 48, -2.4, 2.4);
                FillHist(histkey+"/phi", bjets.at(i).Phi(), weight, 64, -3.2, 3.2);
            } 
            FillHist(channel+"/"+syst+"/METv/pt", METv.Pt(), weight, 300, 0., 300.);
            FillHist(channel+"/"+syst+"/METv/phi", METv.Phi(), weight, 64, -3.2, 3.2);
        }
    }
    else if (channel == "isDiMu") {
        Muon &mu1 = tightMuons.at(0);
        Muon &mu2 = tightMuons.at(1);
        Particle pair = mu1+mu2;
        if (! ev.PassTrigger(DblMuTriggers)) return;
        FillHist("DYDiMu/cutflow", 2., 1., 10, 0., 10.);
        FillHist("TTDiMu/cutflow", 2., 1., 10, 0., 10.);
        const bool passSafeCut = (mu1.Pt() > 20. && mu2.Pt() > 10.);
        if (! passSafeCut) return;
        FillHist("DYDiMu/cutflow", 3., 1., 10, 0., 10.);
        FillHist("TTDiMu/cutflow", 3., 1., 10, 0., 10.);
        if (! (mu1.Charge() + mu2.Charge() == 0)) return;
        FillHist("DYDiMu/cutflow", 4., 1., 10, 0., 10.);
        FillHist("TTDiMu/cutflow", 4., 1., 10, 0., 10.);
        
        // DYDiMu
        if ( fabs(pair.M() - 91.2) < 15.) {
            channel = "DYDiMu";
            FillHist(channel+"/cutflow", 5., 1., 10, 0., 10.);
            if (! (bjets.size() == 0)) return;
            FillHist(channel+"/cutflow", 6., 1., 10, 0., 10.);
        }
        else {
            channel = "TTDiMu";
            if (! (pair.M() > 12.)) return;
            FillHist(channel+"/cutflow", 5., 1., 10, 0., 10.);
            if (! (jets.size() >= 2)) return;
            FillHist(channel+"/cutflow", 6., 1., 10, 0., 10.);
            if (! (bjets.size() >= 1)) return;
            FillHist(channel+"/cutflow", 7., 1., 10, 0., 10.);
            if (! (METv.Pt() > 40.)) return;
            FillHist(channel+"/cutflow", 8., 1., 10, 0., 10.);
        }

        for (const auto &syst: systematics) {
            double weight = 1.;
            if (! IsDATA) {
                weight *= MCweight();
                weight *= ev.GetTriggerLumi("Full");

                if (syst == "L1PrefireUp")        weight *= GetPrefireWeight(1);
                else if (syst == "L1PrefireDown") weight *= GetPrefireWeight(-1);
                else                              weight *= GetPrefireWeight(0);

                if (syst == "PileupReweightUp")        weight *= GetPileUpWeight(nPileUp, 1);
                else if (syst == "PileupReweightDown") weight *= GetPileUpWeight(nPileUp, -1);
                else                                   weight *= GetPileUpWeight(nPileUp, 0);

                weight *= mcCorr->GetBTaggingReweight_1a(jets, jtps.at(0));

                // Lepton ID efficiency correction / Trigger eff. correction
                for (const auto &mu: tightMuons) {
                    if (syst == "MuonIDSFUp")        weight *= getMuonIDSF(mu, 1);
                    else if (syst == "MuonIDSFDown") weight *= getMuonIDSF(mu, -1);
                    else                             weight *= getMuonIDSF(mu, 0);
                }
                if (syst == "DblMuTrigSFUp")        weight *= getDblMuTriggerSF(tightMuons, 1);
                else if (syst == "DblMuTrigSFDown") weight *= getDblMuTriggerSF(tightMuons, -1);
                else                                weight *= getDblMuTriggerSF(tightMuons, 0);
            }

            FillHist(channel+"/"+syst+"/muons/1/pt", mu1.Pt(), weight, 300, 0., 300.);
            FillHist(channel+"/"+syst+"/muons/1/eta", mu1.Eta(), weight, 48, -2.4, 2.4);
            FillHist(channel+"/"+syst+"/muons/1/phi", mu1.Phi(), weight, 64, -3.2, 3.2);
            FillHist(channel+"/"+syst+"/muons/2/pt", mu2.Pt(), weight, 300, 0., 300.);
            FillHist(channel+"/"+syst+"/muons/2/eta", mu2.Eta(), weight, 50, -2.5, 2.5);
            FillHist(channel+"/"+syst+"/muons/2/phi", mu2.Phi(), weight, 64, -3.2, 3.2);
            FillHist(channel+"/"+syst+"/pair/pt", pair.Pt(), weight, 300, 0., 300.);
            FillHist(channel+"/"+syst+"/pair/eta", pair.Eta(), weight, 100, -5., 5.);
            FillHist(channel+"/"+syst+"/pair/phi", pair.Phi(), weight, 64, -3.2, 3.2);
            FillHist(channel+"/"+syst+"/pair/mass", pair.M(), weight, 300, 0., 300.);
            FillHist(channel+"/"+syst+"/jets/size", jets.size(), weight, 20., 0., 20.);
            FillHist(channel+"/"+syst+"/bjets/size", bjets.size(), weight, 20., 0., 20.);
            for (unsigned int i = 0; i < jets.size(); i++) {
                TString histkey = channel+"/"+syst+"/jets/"+TString::Itoa(i+1, 10);
                FillHist(histkey+"/pt", jets.at(i).Pt(), weight, 300, 0., 300.);
                FillHist(histkey+"/eta", jets.at(i).Eta(), weight, 48, -2.4, 2.4);
                FillHist(histkey+"/phi", jets.at(i).Phi(), weight, 64, -3.2, 3.2);
            }
            for (unsigned int i = 0; i < bjets.size(); i++) {
                TString histkey = channel+"/"+syst+"/bjets/"+TString::Itoa(i+1, 10);
                FillHist(histkey+"/pt", bjets.at(i).Pt(), weight, 300, 0., 300.);
                FillHist(histkey+"/eta", bjets.at(i).Eta(), weight, 48, -2.4, 2.4);
                FillHist(histkey+"/phi", bjets.at(i).Phi(), weight, 64, -3.2, 3.2);
            }
            FillHist(channel+"/"+syst+"/METv/pt", METv.Pt(), weight, 300, 0., 300.);
            FillHist(channel+"/"+syst+"/METv/phi", METv.Phi(), weight, 64, -3.2, 3.2);       
        }
    }
    else { 
        return; 
    }
}


double CR_DiLepton::getMuonIDSF(const Muon &mu, int sys) {
    double pt = mu.Pt();
    double eta = fabs(mu.Eta());

    if (pt < 10.) pt = 10;
    if (pt > 200.) pt = 199.;
    if (eta > 2.4) eta = 2.39;
    int thisBin = hMuonIDSF->FindBin(eta, pt);
    double value = hMuonIDSF->GetBinContent(thisBin);
    double error = hMuonIDSF->GetBinError(thisBin);

    return value + int(sys)*error;
}

double CR_DiLepton::getTriggerEff(const Muon &mu, TString histkey, bool isDataEff, int sys) {
    TH2D *h = nullptr;
    double pt = mu.Pt();
    double eta = fabs(mu.Eta());
    if (histkey == "Mu17Leg1" && isDataEff) {
        h = hMu17Leg1_DATA;
        if (pt < 16.) pt = 16.;
        if (pt > 200.) pt = 199.;
        if (eta > 2.4) eta = 2.39;
    }
    else if (histkey == "Mu17Leg1" && (!isDataEff)) {
        h = hMu17Leg1_MC;
        if (pt < 16.) pt = 16.;
        if (pt > 200.) pt = 199.;
        if (eta > 2.4) eta = 2.39;
    }
    else if (histkey == "Mu8Leg2" && isDataEff) {
        h = hMu8Leg2_DATA;
        if (pt < 10.) pt = 10.;
        if (pt > 200.) pt = 199.;
        if (eta > 2.4) eta = 2.39;
    }
    else if (histkey == "Mu8Leg2" && (!isDataEff)) {
        h = hMu8Leg2_MC;
        if (pt < 10.) pt = 10.;
        if (pt > 200.) pt = 199.;
        if (eta > 2.4) eta = 2.39;
    }
    else {
        cerr << "[TriLeptonBase::getTriggerEff] Wrong combination of histkey and isDataEff" << endl;
        cerr << "[TriLeptonBase::getTriggerEff] histkey = " << histkey << endl;
        cerr << "[TriLeptonBase::getTriggerEff] isDataEff = " << isDataEff << endl;
    }

    int thisBin = h->FindBin(eta, pt);
    double value = h->GetBinContent(thisBin);
    double error = h->GetBinError(thisBin);

    return value + int(sys)*error;
} 

double CR_DiLepton::getDblMuTriggerEff(vector<Muon> &muons, bool isDATA, int sys) {
    // check no. of muons
    if (! (muons.size() == 2)) {
        cerr << "[CR_DiLepton::getDblMuTriggerEff] Wrong no. of muons " << muons.size() << endl;
        exit(EXIT_FAILURE);
    }

    const double leadingMuEff = getTriggerEff(muons.at(0), "Mu17Leg1", isDATA, sys);
    const double subleadingMuEff = getTriggerEff(muons.at(1), "Mu8Leg2", isDATA, sys);

    return leadingMuEff * subleadingMuEff;
}

double CR_DiLepton::getDblMuTriggerSF(vector<Muon> &muons, int sys) {
    const double effData = getDblMuTriggerEff(muons, true, sys);
    const double effMC   = getDblMuTriggerEff(muons, false, sys);
    if (effMC == 0 || effData == 0)
        return 1.;

    return effData / effMC;
}

