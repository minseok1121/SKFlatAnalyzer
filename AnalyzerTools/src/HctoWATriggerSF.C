#include "MCCorrection.h"

float MCCorrection::GetTriggerSF(vector<Electron>& EleColl, vector<Muon>& MuColl, TString SFKey, TString Option){

  if(IsDATA) return 1.;

  TString NominalOpt=Option; NominalOpt.ReplaceAll("Syst","");
  bool SystRun=Option.Contains("Syst");
  float SystDir=0., RelSystData=0., RelSystMC=0., TotRelSyst=0.;

  float TriggerEff_Data = TriggerEfficiency(EleColl, MuColl, SFKey, true,  NominalOpt);
  float TriggerEff_MC   = TriggerEfficiency(EleColl, MuColl, SFKey, false, NominalOpt);
  if(SystRun){
    float TriggerEff_Data_syst = TriggerEfficiency(EleColl, MuColl, SFKey, true,  Option);
    float TriggerEff_MC_syst   = TriggerEfficiency(EleColl, MuColl, SFKey, false, Option);
    RelSystData = TriggerEff_Data!=0.? (TriggerEff_Data_syst-TriggerEff_Data)/TriggerEff_Data:0.;
    RelSystMC   = TriggerEff_MC  !=0.? (TriggerEff_MC_syst  -TriggerEff_MC  )/TriggerEff_MC  :0.;
    TotRelSyst  = sqrt(pow(RelSystData,2.)+pow(RelSystMC,2.));
    if     (Option.Contains("Up"))   SystDir =  1.;
    else if(Option.Contains("Down")) SystDir = -1.;
  }
  
  float TriggerScaleFactor = TriggerEff_MC!=0.? TriggerEff_Data/TriggerEff_MC:0.;
   if(TriggerScaleFactor<0) TriggerScaleFactor=0.;
   if(SystRun) TriggerScaleFactor *= (1.+SystDir*TotRelSyst);

  return TriggerScaleFactor;

}


float MCCorrection::TriggerEfficiency(vector<Electron>& EleColl, vector<Muon>& MuColl, TString SFKey, bool ReturnDataEff, TString Option){
  //DataorMC : T: Return DataEff, F: Return MCEff

  if(IsDATA) return 1.;

  TString StrMCorData = ReturnDataEff? "DATA":"MC";
  int SystDir=0;
  if(Option.Contains("Syst")){
    if     (Option.Contains("Up"))   SystDir= 1.;
    else if(Option.Contains("Down")) SystDir=-1.;
  }

  bool SiglMuTrig=false, SiglElTrig=false, DiMuTrig=false, DiElTrig=false, EMuTrig=false;
  float MinPt1=-1, MaxPt1=-1, MinPt2=-1, MaxPt2=-1, MinPt3=-1, MaxPt3=-1., MinPt4=-1, MaxPt4=-1, MaxfEta1=-1, MaxfEta2=-1;
  TH2F* HistEff1=NULL; TH2F* HistEff2=NULL; TH2F* HistEff3=NULL; TH2F* HistEff4=NULL;
  if(DataYear==2016 && SFKey.Contains("IsoORTkIsoMu24_POGTight")){
    SiglMuTrig=true, MinPt1=26., MaxPt1=500., MaxfEta1=2.4; 
    HistEff1 = map_hist_Muon["Trigger_Eff_"+StrMCorData+"_IsoMu24_POGTight"];
  }
  else if(DataYear==2017 && SFKey.Contains("IsoMu27_POGTight")){
    SiglMuTrig=true, MinPt1=29., MaxPt1=1200., MaxfEta1=2.4; 
    HistEff1 = map_hist_Muon["Trigger_Eff_"+StrMCorData+"_IsoMu27_POGTight"];
  }
  else if(DataYear==2018 && SFKey.Contains("IsoMu24_POGTight")){
    SiglMuTrig=true, MinPt1=26., MaxPt1=1200., MaxfEta1=2.4; 
    HistEff1 = map_hist_Muon["Trigger_Eff_"+StrMCorData+"_IsoMu24_POGTight"];
  }
  else if(DataYear==2016 && SFKey.Contains("Ele27WPTight_POGMVAIsoWP90")){
    SiglElTrig=true, MinPt1=35., MaxPt1=500., MaxfEta1=2.5;
    HistEff1 = map_hist_Electron["Trigger_Eff_"+StrMCorData+"_Ele27WPTight_POGMVAIsoWP90"];
  }
  else if(DataYear==2017 && SFKey.Contains("Ele32WPTight1OR2_POGMVAIsoWP90")){
    SiglElTrig=true, MinPt1=35., MaxPt1=500., MaxfEta1=2.5;
    HistEff1 = map_hist_Electron["Trigger_Eff_"+StrMCorData+"_Ele32WPTight1OR2_POGMVAIsoWP90"];
  }
  else if(DataYear==2018 && SFKey.Contains("Ele32WPTight_POGMVAIsoWP90")){
    SiglElTrig=true, MinPt1=35., MaxPt1=500., MaxfEta1=2.5;
    HistEff1 = map_hist_Electron["Trigger_Eff_"+StrMCorData+"_Ele32WPTight_POGMVAIsoWP90"];
  }
  else if(SFKey.Contains("DiMuIso_HNTopID")){
    DiMuTrig=true; MinPt1=20., MinPt2=10., MaxPt1=200., MaxPt2=200., MaxfEta1=2.4;
    HistEff1 = map_hist_Muon["Trigger_Eff_"+StrMCorData+"_DiMuIsoMu17_HNTopID"];
    HistEff2 = map_hist_Muon["Trigger_Eff_"+StrMCorData+"_DiMuIsoMu8_HNTopID"];
  }
  else if(SFKey=="DiElIso_HNTopID" or SFKey=="DiElIso_HNTopIDSS"){
    DiElTrig=true; MinPt1=25., MinPt2=15., MaxPt1=200., MaxPt2=200., MaxfEta1=2.5;
    TString SSLabel = SFKey.Contains("IDSS")? "SS":"";
    HistEff1 = map_hist_Electron["Trigger_Eff_"+StrMCorData+"_DiElIsoEl23_HNTopID"+SSLabel];
    HistEff2 = map_hist_Electron["Trigger_Eff_"+StrMCorData+"_DiElIsoEl12_HNTopID"+SSLabel];
  }
  else if(SFKey=="EMuIso_HNTopID" or SFKey=="EMuIso_HNTopIDSS"){
    EMuTrig=true;
    TString SSLabel = SFKey.Contains("IDSS")? "SS":"";
    MinPt1=25., MinPt2=15., MaxPt1=200., MaxPt2=200., MaxfEta1=2.5;
    MinPt3=25., MinPt4=10., MaxPt3=200., MaxPt4=200., MaxfEta2=2.4;
    HistEff1 = map_hist_Electron["Trigger_Eff_"+StrMCorData+"_EMuIsoEl23_HNTopID"+SSLabel];
    HistEff2 = map_hist_Electron["Trigger_Eff_"+StrMCorData+"_EMuIsoEl12_HNTopID"+SSLabel];
    HistEff3 = map_hist_Muon["Trigger_Eff_"+StrMCorData+"_EMuIsoMu23_HNTopID"];
    HistEff4 = map_hist_Muon["Trigger_Eff_"+StrMCorData+"_EMuIsoMu8_HNTopID"];
  }
  MinPt1+=1E-5, MinPt2+=1E-5, MaxPt1-=1E-5, MaxPt2-=1E-5, MaxfEta1-=1E-5;
  MinPt3+=1E-5, MinPt4+=1E-5, MaxPt3-=1E-5, MaxPt4-=1E-5, MaxfEta2-=1E-5;

  if(!HistEff1){cerr<<"[MCCorrection::MuonTrigger_Eff] No eff file for "<<SFKey<<endl; exit(EXIT_FAILURE);}
  if( (DiMuTrig or DiElTrig) && (!HistEff2) ){cerr<<"[MCCorrection::MuonTrigger_Eff] No eff file for leg2 "<<SFKey<<endl; exit(EXIT_FAILURE);}
  if( EMuTrig && !(HistEff2 && HistEff3 && HistEff4)){ cerr<<"[MCCorrection::MuonTrigger_Eff] No eff file for "<<SFKey<<endl; exit(EXIT_FAILURE); }


  float TriggerEff=0.;
  int NMu=MuColl.size(), NEl=EleColl.size();
  if(SiglMuTrig){
    for(unsigned int it_m=0; it_m<MuColl.size(); it_m++){
      float pt   = MuColl.at(it_m).Pt();
      float feta = fabs(MuColl.at(it_m).Eta());
      if     (pt<MinPt1) return 1.;
      else if(pt>MaxPt1) return 1.;
      if     (feta>MaxfEta1) return 1.;

      int BinIdx = HistEff1->FindBin(pt, feta);
      TriggerEff = HistEff1->GetBinContent(BinIdx);
      if(SystDir!=0){ TriggerEff += float(SystDir)*HistEff1->GetBinError(BinIdx); }
    }
  }
  else if(SiglElTrig){
    for(unsigned int it_e=0; it_e<EleColl.size(); it_e++){
      float pt   = EleColl.at(it_e).Pt();
      float feta = fabs(EleColl.at(it_e).Eta());
      if     (pt<MinPt1) return 1.;
      else if(pt>MaxPt1) return 1.;
      if     (feta>MaxfEta1) return 1.;

      int BinIdx = HistEff1->FindBin(pt, feta);
      TriggerEff = HistEff1->GetBinContent(BinIdx);
      if(SystDir!=0){ TriggerEff += float(SystDir)*HistEff1->GetBinError(BinIdx); }
    }
  }
  else if(DiMuTrig){
    if(NMu==2){
      float pt1  = MuColl.at(0).Pt(), pt2 = MuColl.at(1).Pt();
      float eta1 = MuColl.at(0).Eta(), eta2 = MuColl.at(1).Eta();
      pt1 = min(max(pt1,MinPt1),MaxPt1), pt2 = min(max(pt2,MinPt2),MaxPt2);
      eta1 = min(max(eta1,((float)-1.)*MaxfEta1),MaxfEta1), eta2 = min(max(eta2,((float)-1.)*MaxfEta1),MaxfEta1);

      float EffLeg1_Mu1 = HistEff1->GetBinContent(HistEff1->FindBin(eta1, pt1));
      float EffLeg2_Mu2 = HistEff2->GetBinContent(HistEff2->FindBin(eta2, pt2));
      float ErrLeg1_Mu1 = HistEff1->GetBinError(HistEff1->FindBin(eta1, pt1));
      float ErrLeg2_Mu2 = HistEff2->GetBinError(HistEff2->FindBin(eta2, pt2));
      if(SystDir!=0){ EffLeg1_Mu1+=float(SystDir)*ErrLeg1_Mu1; EffLeg2_Mu2+=float(SystDir)*ErrLeg2_Mu2; }
      TriggerEff = EffLeg1_Mu1*EffLeg2_Mu2;
    }
  }
  else if(DiElTrig){
    if(NEl==2){
      float pt1  = EleColl.at(0).Pt(), pt2 = EleColl.at(1).Pt();
      float eta1 = EleColl.at(0).Eta(), eta2 = EleColl.at(1).Eta();
      pt1 = min(max(pt1,MinPt1),MaxPt1), pt2 = min(max(pt2,MinPt2),MaxPt2);
      eta1 = min(max(eta1,((float)-1.)*MaxfEta1),MaxfEta1), eta2 = min(max(eta2,((float)-1.)*MaxfEta1),MaxfEta1);

      float EffLeg1_El1 = HistEff1->GetBinContent(HistEff1->FindBin(eta1, pt1));
      float EffLeg2_El2 = HistEff2->GetBinContent(HistEff2->FindBin(eta2, pt2));
      float ErrLeg1_El1 = HistEff1->GetBinError(HistEff1->FindBin(eta1, pt1));
      float ErrLeg2_El2 = HistEff2->GetBinError(HistEff2->FindBin(eta2, pt2));
      if(SystDir!=0){ EffLeg1_El1+=float(SystDir)*ErrLeg1_El1; EffLeg2_El2+=float(SystDir)*ErrLeg2_El2; }
      TriggerEff = EffLeg1_El1*EffLeg2_El2;
    }
  }
  else if(EMuTrig){
    if(NEl==1 && NMu==1){
      float pt_m = MuColl.at(0).Pt(), pt_e = EleColl.at(0).Pt();
      float eta_m = fabs(MuColl.at(0).Eta()), eta_e = fabs(EleColl.at(0).Eta());//efficiency folded.
      eta_m = min(eta_m,MaxfEta2), eta_e = min(eta_e,MaxfEta1);
      //eta_m = min(max(eta_m,((float)-1.)*MaxfEta2),MaxfEta2), eta_e = min(max(eta_e,((float)-1.)*MaxfEta1),MaxfEta1);

      float Eff_Mu = 0., Eff_El=0., Err_Mu=0., Err_El=0.;;
      if(pt_e>MinPt1){
        pt_m   = min(max(pt_m,MinPt4),MaxPt4);
        Eff_Mu = HistEff4->GetBinContent(HistEff4->FindBin(pt_m, eta_m));
        Err_Mu = HistEff4->GetBinError(HistEff4->FindBin(pt_m, eta_m));
      }
      else{
        pt_m   = min(max(pt_m,MinPt3),MaxPt3);
        Eff_Mu = HistEff3->GetBinContent(HistEff3->FindBin(pt_m, eta_m));
        Err_Mu = HistEff3->GetBinError(HistEff3->FindBin(pt_m, eta_m));
      }
      if(pt_m>MinPt3){  
        pt_e   = min(max(pt_e,MinPt2),MaxPt2);
        Eff_El = HistEff2->GetBinContent(HistEff2->FindBin(pt_e, eta_e));
        Err_El = HistEff2->GetBinError(HistEff2->FindBin(pt_e, eta_e));
      }
      else{
        pt_e   = min(max(pt_e,MinPt1),MaxPt1);
        Eff_El = HistEff1->GetBinContent(HistEff1->FindBin(pt_e, eta_e));
        Err_El = HistEff1->GetBinError(HistEff1->FindBin(pt_e, eta_e));
      }
      if(SystDir!=0){ Eff_Mu+=float(SystDir)*Err_Mu; Eff_El+=float(SystDir)*Err_El; }
      TriggerEff = Eff_Mu*Eff_El;
    }
  }

  return TriggerEff;
}
