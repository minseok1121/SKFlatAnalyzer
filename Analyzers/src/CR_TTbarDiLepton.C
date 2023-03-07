#include "CR_TTbarDiLepton.h"

CR_TTbarDiLepton::CR_TTbarDiLepton(){}
CR_TTbarDiLepton::~CR_TTbarDiLepton(){}

void CR_TTbarDiLepton::initializeAnalyzer(){
    // flags
    TTDiMu = HasFlag("TTDiMu");
    TTEMu  = HasFlag("TTEMu");
    DYDiMu = HasFlag("DYDiMu");

    // triggers & ID settings
    if (DataEra == "2016preVFP") {
        DblMuTriggers = {
            "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v",
            "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v",
        };
        EMuTriggers = {
            //"HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v",
            //"HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",
            "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v",
            "HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v",
            //"HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v"
        };
        MuonIDs = {"HcToWATight", "HcToWALoose", "HcToWAVeto"};
        ElectronIDs = {"HcToWATight16a", "HcToWALoose16a", "HcToWAVeto16a"};
    }
    else if (DataEra == "2016postVFP") {
        DblMuTriggers = {
            "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v",
            "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v",
            "HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v"
        };
        EMuTriggers = {
            //"HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v",
            "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",
            //"HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v",
            "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v"
        };
        MuonIDs = {"HcToWATight", "HcToWALoose", "HcToWAVeto"};
        ElectronIDs = {"HcToWATight16b", "HcToWALoose16b", "HcToWAVeto16b"};
    }
    else if (DataEra == "2017") {
        DblMuTriggers = {
            "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v",
            //"HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v",
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
            //"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v",
            //"HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v",
            "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v"
        };
        EMuTriggers = {
            //"HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v",
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
}

void CR_TTbarDiLepton::executeEvent(){
    

    FillHist("cutflow", 0., 1., 10, 0., 10.);
    if (! PassMETFilter()) return;
    FillHist("cutflow", 1., 1., 10, 0., 10.);

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
        channel = "TTDiMu";
    else if (tightMuons.size() == 1 && vetoMuons.size() == 1 && tightElectrons.size() == 1 && vetoElectrons.size() == 1)
        channel = "TTEMu";
    else
        return;


    
    if (TTEMu && (channel == "TTEMu")) {
        Electron &ele = tightElectrons.at(0);
        Muon &mu = tightMuons.at(0);
        if (! ev.PassTrigger(EMuTriggers)) return;         // pass trigger
        FillHist("cutflow", 2, 1., 10, 0., 10.);
        const bool passSafeCut = ((mu.Pt() > 25. && ele.Pt() > 15.) || (mu.Pt() > 10. && ele.Pt() > 25.));
        if (! passSafeCut) return;                                 // pass safe cut
        FillHist("cutflow", 3, 1., 10, 0., 10.); 
        if (! (ele.Charge() + mu.Charge() == 0)) return;    // OS charge condition
        FillHist("cutflow", 4, 1., 10, 0., 10.);
        if (! (jets.size() >= 2)) return;                   // Nj >= 2
        FillHist("cutflow", 5, 1., 10, 0., 10.);
        if (! (bjets.size() >= 1)) return;                  // Nb >= 2
        FillHist("cutflow", 6, 1., 10, 0., 10.);

        double weight = 1.;
        if (! IsDATA) {
            weight *= MCweight();
            weight *= ev.GetTriggerLumi("Full");
            weight *= GetPrefireWeight(0);
            weight *= GetPileUpWeight(nPileUp, 0);
            weight *= mcCorr->GetBTaggingReweight_1a(jets, jtps.at(0));

            // Lepton ID efficiency correction / Trigger eff. correction
        }

        FillHist(channel+"/muon/pt", mu.Pt(), weight, 300, 0., 300.);
        FillHist(channel+"/muon/eta", mu.Eta(), weight, 48, -2.4, 2.4);
        FillHist(channel+"/muon/phi", mu.Phi(), weight, 64, -3.2, 3.2);
        FillHist(channel+"/electron/pt", ele.Pt(), weight, 300, 0., 300.);
        FillHist(channel+"/electron/eta", ele.Eta(), weight, 50, -2.5, 2.5);
        FillHist(channel+"/electron/phi", ele.Phi(), weight, 64, -3.2, 3.2);
        FillHist(channel+"/jets/size", jets.size(), weight, 20, 0., 20.);
        for (unsigned int i = 0; i < jets.size(); i++) {
            TString histkey = channel+"/jets/"+TString::Itoa(i+1, 10);
            FillHist(histkey+"/pt", jets.at(i).Pt(), weight, 300, 0., 300.);
            FillHist(histkey+"/eta", jets.at(i).Eta(), weight, 48, -2.4, 2.4);
            FillHist(histkey="/phi", jets.at(i).Phi(), weight, 64, -3.2, 3.2);
        }
        for (unsigned int i = 0; i < bjets.size(); i++) {
            TString histkey = channel+"/bjets/"+TString::Itoa(i+1, 10);
            FillHist(histkey+"/pt", bjets.at(i).Pt(), weight, 300, 0., 300.);
            FillHist(histkey+"/eta", bjets.at(i).Eta(), weight, 48, -2.4, 2.4);
            FillHist(histkey="/phi", bjets.at(i).Phi(), weight, 64, -3.2, 3.2);
        } 
        FillHist(channel+"/METv/pt", METv.Pt(), weight, 300, 0., 300.);
        FillHist(channel+"/METv/phi", METv.Phi(), weight, 64, -3.2, 3.2);
    }
    else if (TTDiMu && (channel == "TTDiMu")) {
        Muon &mu1 = tightMuons.at(0);
        Muon &mu2 = tightMuons.at(1);
        Particle pair = mu1+mu2;
        if (! ev.PassTrigger(DblMuTriggers)) return;
        FillHist("cutflow", 2, 1., 10, 0., 10.);
        const bool passSafeCut = (mu1.Pt() > 20. && mu2.Pt() > 10.);
        if (! passSafeCut) return;
        FillHist("cutflow", 3, 1., 10, 0., 10.);
        if (! (mu1.Charge() + mu2.Charge() == 0)) return;
        FillHist("cutflow", 4, 1., 10, 0., 10);
        if (! (fabs(pair.M() - 91.2) > 15.)) return;
        FillHist("cutflow", 5, 1., 10, 0., 10.);
        if (! (jets.size() >= 2)) return;
        FillHist("cutflow", 6, 1., 10, 0., 10.);
        if (! (bjets.size() >= 1)) return;
        FillHist("cutflow", 7, 1., 10, 0., 10.);
        if (! (pair.M() > 12.)) return;
        FillHist("cutflow", 8, 1., 10, 0., 10.);
        if (! (METv.Pt() > 40.)) return;
        FillHist("cutflow", 9, 1., 10, 0., 10.);

        /*
        AnalyzerParameter param;
        TString MuonIDSFKey = "NUM_MediumID_DEN_TrackerMuons";
        param.Clear();
        //param.syst_ = AnalyzerParameter::Central;
        //param.Name = MuonIDs+"_"+"Central";
        //param.Muon_Tight_ID = MuonID;
        param.Muon_ID_SF_Key = MuonIDSFKey;
        //param.Jet_ID = "tight";
        //executeEventFromParameter(param);
        */


        double weight = 1.;
        if (! IsDATA) {
            weight *= MCweight();
            weight *= ev.GetTriggerLumi("Full");
            weight *= GetPrefireWeight(0);
            weight *= GetPileUpWeight(nPileUp, 0);
            weight *= mcCorr->GetBTaggingReweight_1a(jets, jtps.at(0));
            /*
            for(unsigned int i=0; i<tightMuons.size(); i++){
                //double this_idsf = 1.;
                double this_idsf = mcCorr->MuonID_SF_MS ("TightID", tightMuons.at(i).Eta(), tightMuons.at(i).MiniAODPt());
                weight *= this_idsf;
            }
            */
        }

        FillHist(channel+"/muon/1/pt", mu1.Pt(), weight, 300, 0., 300.);
        FillHist(channel+"/muon/1/eta", mu1.Eta(), weight, 48, -2.4, 2.4);
        FillHist(channel+"/muon/1/phi", mu1.Phi(), weight, 64, -3.2, 3.2);
        FillHist(channel+"/muon/2/pt", mu2.Pt(), weight, 300, 0., 300.);
        FillHist(channel+"/muon/2/eta", mu2.Eta(), weight, 50, -2.5, 2.5);
        FillHist(channel+"/muon/2/phi", mu2.Phi(), weight, 64, -3.2, 3.2);
        FillHist(channel+"/pair/pt", pair.Pt(), weight, 300, 0., 300.);
        FillHist(channel+"/pair/eta", pair.Eta(), weight, 100, -5., 5.);
        FillHist(channel+"/pair/phi", pair.Phi(), weight, 64, -3.2, 3.2);
        FillHist(channel+"/pair/mass", pair.M(), weight, 300, 0., 300.);
        FillHist(channel+"/jets/size", jets.size(), weight, 30, 0., 30.);
        FillHist(channel+"/bjets/size", bjets.size(), weight, 30, 0., 30.);
        for (unsigned int i = 0; i < jets.size(); i++) {
            TString histkey = channel+"/jets/"+TString::Itoa(i+1, 10);
            FillHist(histkey+"/pt", jets.at(i).Pt(), weight, 300, 0., 300.);
            FillHist(histkey+"/eta", jets.at(i).Eta(), weight, 48, -2.4, 2.4);
            FillHist(histkey+"/phi", jets.at(i).Phi(), weight, 64, -3.2, 3.2);
        }
        for (unsigned int i = 0; i < bjets.size(); i++) {
            TString histkey = channel+"/bjets/"+TString::Itoa(i+1, 10);
            FillHist(histkey+"/pt", bjets.at(i).Pt(), weight, 300, 0., 300.);
            FillHist(histkey+"/eta", bjets.at(i).Eta(), weight, 48, -2.4, 2.4);
            FillHist(histkey+"/phi", bjets.at(i).Phi(), weight, 64, -3.2, 3.2);
        }
        FillHist(channel+"/METv/pt", METv.Pt(), weight, 300, 0., 300.);
        FillHist(channel+"/METv/phi", METv.Phi(), weight, 64, -3.2, 3.2);       
    }
    else if (DYDiMu && (channel == "TTDiMu")) {
        Muon &mu1 = tightMuons.at(0);
        Muon &mu2 = tightMuons.at(1);
        Particle pair = mu1+mu2;
        if (! ev.PassTrigger(DblMuTriggers)) return;
        FillHist("cutflow", 2, 1., 10, 0., 10.);
        const bool passSafeCut = (mu1.Pt() > 20. && mu2.Pt() > 10.);
        if (! passSafeCut) return;
        FillHist("cutflow", 3, 1., 10, 0., 10.);
        if (! (mu1.Charge() + mu2.Charge() == 0)) return;
        FillHist("cutflow", 4, 1., 10, 0., 10);
        if (! (fabs(pair.M() - 91.2) < 15.)) return;
        FillHist("cutflow", 5, 1., 10, 0., 10.);
        if ( (bjets.size() >= 1)) return;
        FillHist("cutflow", 6, 1., 10, 0., 10.);
        
        /*
        AnalyzerParameter param;
        TString MuonIDSFKey = "NUM_MediumID_DEN_TrackerMuons";
        param.Clear();
        //param.syst_ = AnalyzerParameter::Central;
        //param.Name = MuonIDs+"_"+"Central";
        //param.Muon_Tight_ID = MuonID;
        param.Muon_ID_SF_Key = MuonIDSFKey;
        //param.Jet_ID = "tight";
        //executeEventFromParameter(param);
        */


        double weight = 1.;
        if (! IsDATA) {
            weight *= MCweight();
            weight *= ev.GetTriggerLumi("Full");
            weight *= GetPrefireWeight(0);
            weight *= GetPileUpWeight(nPileUp, 0);
            weight *= mcCorr->GetBTaggingReweight_1a(jets, jtps.at(0));
            /*
            for(unsigned int i=0; i<tightMuons.size(); i++){
                //double this_idsf = 1.;
                double this_idsf = mcCorr->MuonID_SF_MS ("TightID", tightMuons.at(i).Eta(), tightMuons.at(i).MiniAODPt());
                weight *= this_idsf;
            }
            */
        }

        FillHist(channel+"/muon/1/pt", mu1.Pt(), weight, 300, 0., 300.);
        FillHist(channel+"/muon/1/eta", mu1.Eta(), weight, 48, -2.4, 2.4);
        FillHist(channel+"/muon/1/phi", mu1.Phi(), weight, 64, -3.2, 3.2);
        FillHist(channel+"/muon/2/pt", mu2.Pt(), weight, 300, 0., 300.);
        FillHist(channel+"/muon/2/eta", mu2.Eta(), weight, 50, -2.5, 2.5);
        FillHist(channel+"/muon/2/phi", mu2.Phi(), weight, 64, -3.2, 3.2);
        FillHist(channel+"/pair/pt", pair.Pt(), weight, 300, 0., 300.);
        FillHist(channel+"/pair/eta", pair.Eta(), weight, 100, -5., 5.);
        FillHist(channel+"/pair/phi", pair.Phi(), weight, 64, -3.2, 3.2);
        FillHist(channel+"/pair/mass", pair.M(), weight, 300, 0., 300.);
        FillHist(channel+"/jets/size", jets.size(), weight, 30, 0., 30.);
        FillHist(channel+"/bjets/size", bjets.size(), weight, 30, 0., 30.);
        for (unsigned int i = 0; i < jets.size(); i++) {
            TString histkey = channel+"/jets/"+TString::Itoa(i+1, 10);
            FillHist(histkey+"/pt", jets.at(i).Pt(), weight, 300, 0., 300.);
            FillHist(histkey+"/eta", jets.at(i).Eta(), weight, 48, -2.4, 2.4);
            FillHist(histkey+"/phi", jets.at(i).Phi(), weight, 64, -3.2, 3.2);
        }
        for (unsigned int i = 0; i < bjets.size(); i++) {
            TString histkey = channel+"/bjets/"+TString::Itoa(i+1, 10);
            FillHist(histkey+"/pt", bjets.at(i).Pt(), weight, 300, 0., 300.);
            FillHist(histkey+"/eta", bjets.at(i).Eta(), weight, 48, -2.4, 2.4);
            FillHist(histkey+"/phi", bjets.at(i).Phi(), weight, 64, -3.2, 3.2);
        }
        FillHist(channel+"/METv/pt", METv.Pt(), weight, 300, 0., 300.);
        FillHist(channel+"/METv/phi", METv.Phi(), weight, 64, -3.2, 3.2);       
    }
    else { 
        return; 
    }
}
