#include "TTAnalyzer.h"

TTAnalyzer::TTAnalyzer(){
    //==== constructor of this analyzer ====
}

void TTAnalyzer::initializeAnalyzer(){
    
    //=============================================
    //==== TTbar -> W+W-bbbar -> llvvbb events ====
    //=============================================
    
    //===== if you use "--userflags RunSyst" with SKFlat.py, HasFlag("RunSyst") will return "true"
    RunSyst = HasFlag("RunSyst");
    cout << "[TTAnalyzer::initializeAnalyzer] Runsyst = " << RunSyst << endl;
    
    //==== I defined "vector<TString> MuonIDs;" in Analyzers/include/TTAnalyzer.h
    MuonIDs = { "POGTightWithTightIso"};
    //==== corresponding Muon ID & Iso SF Keys for mcCorr->MuonID_SF()
    MuonIDSFKeys = {"NUM_TightID_DEN_genTracks"};
    MuonIsoSFKeys = {"NUM_TightRelIso_DEN_TightIDandIPCut"};
    
    //=== Also for electrons
    ElectronIDs = {"passTightID"};
    ElectronIDSFKeys = {"passTightID"};
    
    //==== At this point, sample informations (e.g., IsDATA, DataStream, MCSample, or DataYear) are all set
    //==== You can define sample-dependent or year-dependent variables here
    //==== (Example) Year-dependent variables
    //==== I defined "TString IsoMuTriggerName;" and "double TriggerSafePtCut;" in Analyzers/include/TTAnalyzer.h
    //==== IsoMuTriggerName is a year-dependent variable, and you don't want to do "if(DatYaer==~~)" for every event (let's save cpu time).
    //==== Then, do it here, which only ran once for each macro
    if(DataYear == 2016) {
        IsoMuTriggerName = "HLT_IsoMu24_v";
        TriggerSafePtCut = 26.;
    }
    else if(DataYear == 2017) {
        IsoMuTriggerName = "HLT_IsoMu27_v";
        TriggerSafePtCut = 29.;
    }
    else if(DataYear == 2018) {
      IsoMuTriggerName = "HLT_IsoMu24_v";
        TriggerSafePtCut = 26.;
    }
    
    cout << "[TTAnalyzer::initializeAnalyzer] IsoMuTriggerName = " << IsoMuTriggerName << endl;
    cout << "[TTAnalyzer::initializeAnalyzer TriggerSafePtCut = " << TriggerSafePtCut << endl;
    
    //==== Test btagging code
    //==== add taggers and WP that you want to use in analysis
    std::vector<JetTagging::Parameters> jtps;
    jtps.push_back( JetTagging::Parameters(JetTagging::DeepCSV, JetTagging::Medium, JetTagging::incl, JetTagging::comb));
    mcCorr->SetJetTaggingParameters(jtps);
 
    //================================
    //==== Using new PDF
    //==== It consumes so much time, so only being actiavted with --userflags RunNewPDF
    //================================

    RunNewPDF = HasFlag("RunNewPDF");
    cout << "[TTAnalyzer::initializeAnalyzer] RunNewPDF = " << RunNewPDF << endl;
    if(RunNewPDF && !IsDATA){

        LHAPDFHandler LHAPDFHandler_Prod;
        LHAPDFHandler_Prod.CentralPDFName = "NNPDF31_nnlo_hessian_pdfas";
        LHAPDFHandler_Prod.init();

        LHAPDFHandler LHAPDFHandler_New;
        LHAPDFHandler_New.CentralPDFName = "NNPDF31_nlo_hessian_pdfas";
        LHAPDFHandler_New.ErrorSetMember_Start = 1;
        LHAPDFHandler_New.ErrorSetMember_End = 100;
        LHAPDFHandler_New.AlphaSMember_Down = 101;
        LHAPDFHandler_New.AlphaSMember_Up = 102;
        LHAPDFHandler_New.init();

        pdfReweight->SetProdPDF( LHAPDFHandler_Prod.PDFCentral );
        pdfReweight->SetNewPDF( LHAPDFHandler_New.PDFCentral );
        pdfReweight->SetNewPDFErrorSet( LHAPDFHandler_New.PDFErrorSet );
        pdfReweight->SetNewPDFAlphaS( LHAPDFHandler_New.PDFAlphaSDown, LHAPDFHandler_New.PDFAlphaSUp );
    }

    //================================================
    //==== How to estimate xsec errors (PDF & Scale)
    //==== For example, MET
    //================================================

    RunXSecSyst = HasFlag("RunXSecSyst");
    cout << "[TTAnalyzer::initializeAnalyzer] RunXSecSyst = " << RunXSecSyst << endl;
}

TTAnalyzer::~TTAnalyzer() {
    //==== Destructor of this Analyzer ====
}


void TTAnalyzer::executeEvent() {
    //==== *IMPORTANT TO SAVE CPU TIME*
    //==== Every GetMuon() funtion first collect ALL MINIAOD muons with GetAllMuons(),
    //==== and then check ID booleans.
    //==== GetAllMuons not only loops over all MINIAOD muons, but also actually CONSTRUCT muon objects for each muons.
    //==== We are now running systematics, and you don't want to do this for every systematic sources
    //==== So, I defined "vector<Muon> AllMuons;" in Analyzers/include/TTAnalyzer.h,
    //==== and save muons objects at the very beginning of executeEvent().
    //==== Later, do "SelectMuons(AllMuons, ID, pt, eta)" to get muons with ID cuts
    AllMuons = GetAllMuons();
    AllJets = GetAllJets();
    AllElectrons = GetAllElectrons();
    AllGens = GetGens();

    //==== Get L1Prefire reweight
    //==== If data, 1.;
    //==== If MC && DataYear > 2017, 1.;
    //==== If MC && DataYear <= 2017, we have to reweight the event with this value
    //==== I defined "double weight_Prefire;" in Analyzers/include/TTAnalyzer.h
    weight_Prefire = GetPrefireWeight(0);
    
    //==== PileUp reweight
    weight_PileUp = GetPileUpWeight(nPileUp, 0);
    
    //==== Top Pt reweight
    weight_TopPt = mcCorr->GetTopPtReweight(AllGens);
    
    //==== Declare AnalyzerParameter

    AnalyzerParameter param;

    //==== Loop over muon IDs
    
    for(unsigned int i = 0; i < MuonIDs.size(); i++){

        TString MuonID = MuonIDs.at(i);
        TString MuonIDSFKey = MuonIDSFKeys.at(i);
        TString MuonIsoSFKey = MuonIsoSFKeys.at(i);
        TString ElectronID = ElectronIDs.at(i);
        TString ElectronIDSFKey = ElectronIDSFKeys.at(i);
        TString WP = "Tight";
        
        //==== 1) First, let's run Central values of the systematics

        //==== clear parameter set
        param.Clear();
        //==== set which systematic sources you want to run this time
        //==== default syst_ is AnalyzerParameter::Central
        param.syst_ = AnalyzerParameter::Central;
        //==== set name of the parameter set
        //==== this will be used for the directory name of histograms
        param.Name = WP+"_"+"Central";
        //==== You can define lepton ID string here
        param.Muon_Tight_ID = MuonID;
        param.Muon_ID_SF_Key = MuonIDSFKey;
        param.Muon_ISO_SF_Key = MuonIsoSFKey;
        //==== And, Jet ID, Electron ID
        param.Jet_ID = "tight";
        param.Electron_Tight_ID = ElectronID;
        param.Electron_ID_SF_Key = ElectronIDSFKey;
        //==== Now, all parameters are set. Run executeEventFromParameter() with this parameter set
        executeEventFromParameter(param);
        //==== 2) Now, loop over systematic sources
        //==== without --userflag RunSyst, this will not be ran

        if(RunSyst){

          for(int it_syst=1; it_syst<AnalyzerParameter::NSyst; it_syst++){

            //==== Everything else remains same, but only change syst_ and parameter name

            param.syst_ = AnalyzerParameter::Syst(it_syst);
            param.Name = MuonID+"_"+"Syst_"+param.GetSystType();
            executeEventFromParameter(param);
          }
        }
    }
    
    //================================
    //==== Using new PDF
    //================================

    if(RunNewPDF && !IsDATA){
        FillHist("NewPDF_PDFReweight", GetPDFReweight(), 1., 2000, 0.90, 1.10);
        for(int i=0; i<pdfReweight->NErrorSet; i++){
            //cout << "[TTAnalyzer::executeEvent]   " << GetPDFReweight(i) << endl;
            JSFillHist("NewPDF_PDFErrorSet", "PDFReweight_Member_"+TString::Itoa(i,10), GetPDFReweight(i), 1., 2000, 0.90, 1.10);
        }
    }

    //================================================
    //==== How to estimate xsec errors (PDF & Scale)
    //==== For example, MET
    //================================================

    if(RunXSecSyst && !IsDATA){

        Event ev = GetEvent();
        double MET = ev.GetMETVector().Pt();

        //==== 1) PDF Error
        //==== Obtain RMS of the distribution later
        for(unsigned int i=0; i<PDFWeights_Error->size(); i++){
            JSFillHist("XSecError", "MET_PDFError_"+TString::Itoa(i,10), MET, PDFWeights_Error->at(i), 200, 0., 200.);
        }

        //==== 2) PDF AlphaS
        //==== Look for PDF4LHC paper..
        //==== https://arxiv.org/abs/1510.03865
        if(PDFWeights_AlphaS->size()==2){
            JSFillHist("XSecError", "MET_PDFAlphaS_Down", MET, PDFWeights_AlphaS->at(0), 200, 0., 200.);
            JSFillHist("XSecError", "MET_PDFAlphaS_Up", MET, PDFWeights_AlphaS->at(1), 200, 0., 200.);
        }

        //==== 3) Scale
        //==== Obtain the envelop of the distribution later
        for(unsigned int i=0; i<PDFWeights_Scale->size(); i++){
            //==== i=5 and 7 are unphysical
            if(i==5) continue;
            if(i==7) continue;
            JSFillHist("XSecError", "MET_Scale_"+TString::Itoa(i,10), MET, PDFWeights_Scale->at(i), 200, 0., 200.);
        }
    }
}

void TTAnalyzer::executeEventFromParameter(AnalyzerParameter param) {
    
    //=============
    //==== No Cut
    //=============
    
    JSFillHist(param.Name, "NoCut_"+param.Name, 0., 1., 1, 0., 1.);

    //========================
    //==== MET Filter
    //========================

    if(!PassMETFilter()) return;

    Event ev = GetEvent();
    Particle METv = ev.GetMETVector();
    //==============
    //==== Trigger
    //==============
    if(! (ev.PassTrigger(IsoMuTriggerName) )) return;
    // if(! (ev.PassTrigger(IsoEleTriggerName))) return;
    
    //======================
    //==== Copy AllObjects
    //======================

    vector<Muon> this_AllMuons = AllMuons;
    vector<Jet> this_AllJets = AllJets;
    vector<Electron> this_AllElectrons = AllElectrons;
    

    //==== Then, for each systematic sources
    //==== 1) Smear or scale them
    //==== 2) Then apply ID selections
    //==== This order should be explicitly followed
    //==== Below are all variables for available systematic sources
    
    if(param.syst_ == AnalyzerParameter::Central){

    }
    else if(param.syst_ == AnalyzerParameter::JetResUp){
      this_AllJets = SmearJets( this_AllJets, +1 );
      //this_AllFatJets = SmearFatJets( this_AllFatJets, +1 );
    }
    else if(param.syst_ == AnalyzerParameter::JetResDown){
      this_AllJets = SmearJets( this_AllJets, -1 );
      //this_AllFatJets = SmearFatJets( this_AllFatJets, -1 );
    }
    else if(param.syst_ == AnalyzerParameter::JetEnUp){
      this_AllJets = ScaleJets( this_AllJets, +1 );
      //this_AllFatJets = ScaleFatJets( this_AllFatJets, +1 );
    }
    else if(param.syst_ == AnalyzerParameter::JetEnDown){
      this_AllJets = ScaleJets( this_AllJets, -1 );
      //this_AllFatJets = ScaleFatJets( this_AllFatJets, -1 );
    }
	/*
    else if(param.syst_ == AnalyzerParameter::MuonEnUp){
      this_AllMuons = ScaleMuons( this_AllMuons, +1 );
    }
    else if(param.syst_ == AnalyzerParameter::MuonEnDown){
      this_AllMuons = ScaleMuons( this_AllMuons, -1 );
    }
    */
    else if(param.syst_ == AnalyzerParameter::ElectronResUp){
      //this_AllElectrons = SmearElectrons( this_AllElectrons, +1 );
    }
    else if(param.syst_ == AnalyzerParameter::ElectronResDown){
      //this_AllElectrons = SmearElectrons( this_AllElectrons, -1 );
    }
    else if(param.syst_ == AnalyzerParameter::ElectronEnUp){
      //this_AllElectrons = ScaleElectrons( this_AllElectrons, +1 );
    }
    else if(param.syst_ == AnalyzerParameter::ElectronEnDown){
      //this_AllElectrons = ScaleElectrons( this_AllElectrons, -1 );
    }
    else{
      cout << "[TTAnalyzer::executeEventFromParameter] Wrong syst" << endl;
      exit(EXIT_FAILURE);
    }
    
    //==================================================
    //==== Then, apply ID selections using this_AllXXX
    //==================================================

    vector<Muon> loose_muons = SelectMuons(this_AllMuons, "POGLoose", 20., 2.4);
    vector<Muon> muons = SelectMuons(this_AllMuons, param.Muon_Tight_ID, 20., 2.4);
    vector<Jet> jets = SelectJets(this_AllJets, param.Jet_ID, 30., 2.4);
    vector<Electron> loose_electrons = SelectElectrons(this_AllElectrons, "passLooseID", 20., 2.4);
    vector<Electron> electrons = SelectElectrons(this_AllElectrons, param.Electron_Tight_ID, 20., 2.4);
      
    //=======================
    //==== Sort in pt-order
    //=======================

    //==== 1) leptons : after scaling/smearing, pt ordring can differ from MINIAOD
    std::sort(muons.begin(), muons.end(), PtComparing);
    std::sort(electrons.begin(), electrons.end(), PtComparing);
    std::sort(loose_muons.begin(), loose_muons.end(), PtComparing);
    std::sort(loose_electrons.begin(), loose_electrons.end(), PtComparing);
    //==== 2) jets : similar, but also when applying new JEC, ordering is changes. This is important if you use leading jets
    std::sort(jets.begin(), jets.end(), PtComparing);
    
    //jet cleaning
    vector<Jet> clean_jets = JetsVetoLeptonInside(jets, loose_electrons, loose_muons, 0.4);

    int NBjets_NoSF(0), NBjets_WithSF_2a(0);
    JetTagging::Parameters jtp_DeepCSV_Medium = JetTagging::Parameters(JetTagging::DeepCSV, JetTagging::Medium, JetTagging::incl, JetTagging::comb);
  
    for (unsigned int ij = 0; ij < jets.size(); ij++) {
      double this_discr = jets.at(ij).GetTaggerResult(JetTagging::DeepCSV);
      if (this_discr > mcCorr->GetJetTaggingCutValue(JetTagging::DeepCSV, JetTagging::Medium) ) NBjets_NoSF++;
      if (mcCorr->IsBTagged_2a(jtp_DeepCSV_Medium, jets.at(ij)) ) NBjets_WithSF_2a++;
    }

     //==== leading muon has trigger-safe pt
     if( muons.size() == 0) return;
     if( muons.at(0).Pt() <= TriggerSafePtCut ) return;
     
     //==========================
     //==== Event selections ====
     //==========================

     if((muons.size() != 1) || (electrons.size() != 1)) return;
     if(muons.at(0).Charge() * electrons.at(0).Charge() > 0) return;
     if(muons.at(0).DeltaR(electrons.at(0)) < 0.4) return;

     if(clean_jets.size() == 0 || clean_jets.size() == 1) return;
     if (!IsData) {
       if(NBjets_WithSF_2a == 0) return; // for SF applied to MC
       //if(n_bjet_deepcsv_m_noSF == 0) return; // for no SF applied to MC
     }
     else {
       if (NBjets_WithSF_2a == 0) return;
     }
    
     //==== mass of electron + muon
     Particle DiLepton = muons.at(0) + electrons.at(0);

     //===================
     //==== Event weight
     //===================

     double weight = 1.;
     //==== If MC
    
    if(!IsDATA){

      //==== weight_norm_1invpb is set to be event weight normalized to 1 pb-1
      //==== So, you have to multiply trigger luminosity
      //==== you can pass trigger names to ev.GetTriggerLumi(), but if you are using unprescaled trigger, simply pass "Full"

      weight *= weight_norm_1invpb*ev.GetTriggerLumi("Full");

      //==== MCweight is +1 or -1. Should be multiplied if you are using e.g., aMC@NLO NLO samples
      weight *= ev.MCweight();

      //==== L1Prefire reweight
      weight *= weight_Prefire;

      //==== Muon ID SF ====
      double this_MuonIDSF  = mcCorr->MuonID_SF (param.Muon_ID_SF_Key,  muons.at(0).Eta(), muons.at(0).MiniAODPt());
        
      weight *= this_MuonIDSF;

      //==== Electron ID SF ====
      double this_ElecIDSF = mcCorr->ElectronID_SF (param.Electron_ID_SF_Key, electrons.at(0).Eta(), electrons.at(0).Pt());
        
      weight *= this_ElecIDSF;
        
      //==== Muon ISO SF ====
      double this_isosf = mcCorr->MuonISO_SF(param.Muon_ISO_SF_Key, muons.at(0).Eta(), muons.at(0).MiniAODPt());
        
      weight *= this_isosf;
        
 
      //==== Muon Trigger SF ====
      double this_MuonTrigSF = 0.;
      if (DataYear == 2017)
        this_MuonTrigSF = mcCorr->MuonTrigger_SF ("POGTight", "IsoMu27", muons, 0);
      else if (DataYear == 2016 || 2018)
        this_MuonTrigSF = mcCorr->MuonTrigger_SF ("POGTight", "IsoMu24", muons, 0);
      else {
        cout << "[TTAnlayzer::executeEventFromParameter] Wrong Trigger Year" << endl;
        exit(EXIT_FAILURE);
      }
      
      weight *= this_MuonTrigSF;
      
      //==== Pileup reweight
      weight *= weight_PileUp;
      
      //==== Top pt reweight
      weight *= weight_TopPt;
    }
    
    //==========================
    //==== Now fill histograms
    //==========================
      
    // Leptons
    JSFillHist(param.Name, "Electron_pt_"+param.Name, electrons.at(0).Pt(), weight, 38, 20, 400);
    JSFillHist(param.Name, "Muon_pt_"+param.Name, muons.at(0).Pt(), weight, 38, 20, 400);
    JSFillHist(param.Name, "Electron_eta_"+param.Name, electrons.at(0).Eta(), weight, 24, -2.4, 2.4);
    JSFillHist(param.Name, "Muon_eta_"+param.Name, muons.at(0).Eta(), weight, 24, -2.4, 2.4);
    JSFillHist(param.Name, "Electorn_phi_" + param.Name, electrons.at(0).Phi(), weight, 16, -4, 4);
    JSFillHist(param.Name, "Muon_phi_" + param.Name, muons.at(0).Phi(), weight, 16, -4, 4);

    // Jets
    JSFillHist(param.Name, "N_jets_"+param.Name, clean_jets.size(), weight, 8, -0.5, 7.5);
    JSFillHist(param.Name, "N_bjets_"+param.Name, NBjets_WithSF_2a, weight, 5, -0.5, 4.5);
    if (clean_jets.size() > 0) {
      JSFillHist(param.Name, "1st_jet_pt_" + param.Name, clean_jets.at(0).Pt(), weight, 38, 20, 400);
      JSFillHist(param.Name, "1st_jet_eta_" + param.Name, clean_jets.at(0).Eta(), weight, 24, -2.4, 2.4);
      JSFillHist(param.Name, "1st_jet_phi_" + param.Name, clean_jets.at(0).Phi(), weight, 16, -4, 4);
    }
    if (clean_jets.size() > 1) {
    JSFillHist(param.Name, "2nd_jet_pt_" + param.Name, clean_jets.at(1).Pt(), weight, 38, 20, 400);
    JSFillHist(param.Name, "2nd_jet_eta_" + param.Name, clean_jets.at(1).Eta(), weight, 24, -2.4, 2.4);
    JSFillHist(param.Name, "2nd_jet_phi_" + param.Name, clean_jets.at(1).Phi(), weight, 16, -4, 4);
    }
    if (clean_jets.size() > 2) JSFillHist(param.Name, "3rd_jet_pt_" + param.Name, clean_jets.at(2).Pt(), weight, 38, 20, 400);
    if (clean_jets.size() > 3) JSFillHist(param.Name, "4th_jet_pr_" + param.Name, clean_jets.at(3).Pt(), weight, 38, 20, 400);

    // Other objects...
    JSFillHist(param.Name, "MET_pt_"+param.Name, METv.Pt(), weight, 38, 20, 400);
    JSFillHist(param.Name, "MET_phi_"+param.Name, METv.Phi(), weight, 18, -4, 4);
    JSFillHist(param.Name, "Delta_R_"+param.Name, muons.at(0).DeltaR(electrons.at(0)), weight, 20, 0, 5);
    JSFillHist(param.Name, "Delta_phi_" + param.Name, muons.at(0).DeltaPhi(electrons.at(0)), weight, 16, -4, 4);
    JSFillHist(param.Name, "No_of_PV_" + param.Name, nPV, weight, 100, 0, 100);
}
