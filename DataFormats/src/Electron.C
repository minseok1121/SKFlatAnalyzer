#include "Electron.h"

ClassImp(Electron)

Electron::Electron(){

  j_En_up=1.;
  j_En_down=1.;;
  j_Res_up = 1.;
  j_Res_down = 1.;

  j_scEta = -999.;
  j_scPhi = -999.;
  j_scE = -999.;
  j_mvaiso = -999.;
  j_mvanoiso = -999.;
  j_EnergyUnCorr = -999.;
  j_passConversionVeto = false;
  j_NMissingHits = 0;
  j_Full5x5_sigmaIetaIeta = -999.;
  j_dEtaSeed = -999.;
  j_dPhiIn = -999.;
  j_HoverE  = -999.;
  j_InvEminusInvP = -999.;
  j_e2x5OverE5x5 = -999.;
  j_e1x5OverE5x5 = -999.;
  j_trkiso = -999.;
  j_dr03EcalRecHitSumEt = -999.;
  j_dr03HcalDepth1TowerSumEt = -999.;
  j_dr03HcalTowerSumEt = -999.;
  j_dr03TkSumPt = -999.;
  j_ecalPFClusterIso = -999.;
  j_hcalPFClusterIso = -999.;
  j_isEcalDriven = false;
  j_IDBit = 0;
  j_IDCutBit.clear();
  j_Rho = -999.;
  j_isGsfCtfScPixChargeConsistent = false;
  this->SetLeptonFlavour(ELECTRON);
}

Electron::~Electron(){

}

void Electron::SetEnShift(double en_up, double en_down){
  j_En_up = en_up;
  j_En_down = en_down;
}

void Electron::SetResShift(double res_up, double res_down){
  j_Res_up = res_up;
  j_Res_down = res_down;
}

void Electron::SetSC(double sceta, double scphi, double sce){
  j_scEta = sceta;
  j_scPhi = scphi;
  j_scE = sce;
}

void Electron::SetMVA(double mvaiso, double mvanoiso){
  j_mvaiso = mvaiso;
  j_mvanoiso = mvanoiso;
}

void Electron::SetUncorrE(double une){
  j_EnergyUnCorr = une;
}

void Electron::SetPassConversionVeto(bool b){
  j_passConversionVeto = b;
}

void Electron::SetNMissingHits(int n){
  j_NMissingHits = n;
}

void Electron::SetCutBasedIDVariables(
    double Full5x5_sigmaIetaIeta,
    double dEtaSeed,
    double dPhiIn,
    double HoverE,
    double InvEminusInvP,
    double e2x5OverE5x5,
    double e1x5OverE5x5,
    double trackIso,
    double dr03EcalRecHitSumEt,
    double dr03HcalDepth1TowerSumEt,
    double dr03HcalTowerSumEt,
    double dr03TkSumPt,
    double ecalPFClusterIso,
    double hcalPFClusterIso,
    int ecalDriven
  ){
  j_Full5x5_sigmaIetaIeta = Full5x5_sigmaIetaIeta;
  j_dEtaSeed = dEtaSeed;
  j_dPhiIn = dPhiIn;
  j_HoverE = HoverE;
  j_InvEminusInvP = InvEminusInvP;
  j_e2x5OverE5x5 = e2x5OverE5x5;
  j_e1x5OverE5x5 = e1x5OverE5x5;
  j_trkiso = trackIso;
  j_dr03EcalRecHitSumEt = dr03EcalRecHitSumEt;
  j_dr03HcalDepth1TowerSumEt = dr03HcalDepth1TowerSumEt;
  j_dr03HcalTowerSumEt = dr03HcalTowerSumEt;
  j_dr03TkSumPt = dr03TkSumPt;
  j_ecalPFClusterIso = ecalPFClusterIso;
  j_hcalPFClusterIso = hcalPFClusterIso;

  if(ecalDriven==0) j_isEcalDriven = false;
  else j_isEcalDriven = true;
}

void Electron::SetIDBit(unsigned int idbit){
  j_IDBit = idbit;
}

void Electron::SetIDCutBit(vector<int> idcutbit){
  j_IDCutBit = idcutbit;
}

void Electron::SetRelPFIso_Rho(double r){
  j_RelPFIso_Rho = r;
  this->SetRelIso(r);
}

double Electron::EA(){

  //==== RecoEgamma/ElectronIdentification/data/Fall17/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_94X.txt
  
  double eta = fabs(this->scEta());

  if     (eta<1.000) return 0.1440;
  else if(eta<1.479) return 0.1562;
  else if(eta<2.000) return 0.1032;
  else if(eta<2.200) return 0.0859;
  else if(eta<2.300) return 0.1116;
  else if(eta<2.400) return 0.1321;
  else if(eta<2.500) return 0.1654;
  else return 0.1654;

}

bool Electron::PassID(TString ID) const{

  //==== XXX veto Gap Always
  if(etaRegion()==GAP) return false;

  //==== POG
  if(ID=="passVetoID") return passVetoID();
  if(ID=="passLooseID") return passLooseID();
  if(ID=="passMediumID") return passMediumID();
  if(ID=="passTightID") return passTightID();
  if(ID=="passHEEPID") return passHEEPID();
  if(ID=="passMVAID_noIso_WP80") return passMVAID_noIso_WP80();
  if(ID=="passMVAID_noIso_WP90") return passMVAID_noIso_WP90();
  if(ID=="passMVAID_iso_WP80") return passMVAID_iso_WP80();
  if(ID=="passMVAID_iso_WP90") return passMVAID_iso_WP90();
  //==== Customized
  if(ID=="FakeTightID") return Pass_FakeTight();
  if(ID=="FakeLooseID") return Pass_FakeLoose();
  if(ID=="SUSYTight") return Pass_SUSYTight();
  if(ID=="SUSYLoose") return Pass_SUSYLoose();
  if(ID=="NOCUT") return true;
  if(ID=="TEST") return Pass_TESTID();
  if(ID=="fr_elec_tight") return Pass_fr_elec_tight();
  if(ID=="fr_elec_loose") return Pass_fr_elec_loose();

  cout << "[Electron::PassID] No id : " << ID << endl;
  exit(EXIT_FAILURE);

  return false;
}

bool Electron::Pass_SUSYMVAWP(TString wp) const{

  double sceta = fabs(scEta());

    double cutA = 0.77;
    double cutB = 0.52;

  if(wp=="Tight"){
    if     (sceta<0.8)  { cutA = 0.77; cutB = 0.52; }
    else if(sceta<1.479){ cutA = 0.56; cutB = 0.11; } 
    else                { cutA = 0.48; cutB =-0.01; }
  }
  else if(wp=="Loose"){
    if     (sceta<0.8)  { cutA =-0.48; cutB =-0.85; }
    else if(sceta<1.479){ cutA =-0.67; cutB =-0.91; }
    else                { cutA =-0.49; cutB =-0.83; }
  }
  else{}

  double cutSlope = (cutB-cutA) / 10.;
  double cutFinal = std::min( cutA, std::max(cutB , cutA + cutSlope*(this->Pt()-15.) ) );

  //==== Using NoIso MVA, because we apply MiniIso later
  if(MVANoIso()>cutFinal) return true;
  else return false;

}

bool Electron::Pass_SUSYTight() const{
  if(! Pass_SUSYMVAWP("Tight") ) return false;
  if(! (MiniRelIso()<0.1) ) return false;	
  if(! (fabs(dXY())<0.05 && fabs(dZ())<0.1 && fabs(IP3D()/IP3Derr())<8.) ) return false;
  if(! PassConversionVeto() ) return false;
  if(! (NMissingHits()==0) ) return false;

  return true;
}

bool Electron::Pass_SUSYLoose() const{
  if(! Pass_SUSYMVAWP("Loose") ) return false;
  if(! (MiniRelIso()<0.4) ) return false;
  if(! (fabs(dXY())<0.05 && fabs(dZ())<0.1 && fabs(IP3D()/IP3Derr())<8.) ) return false;
  if(! PassConversionVeto() ) return false;
  if(! (NMissingHits()==0) ) return false;

  return true;
}

//==== TEST ID

bool Electron::Pass_TESTID() const{
  return true;
}



bool Electron::Pass_CutBasedLooseNoIso() const{

  if( fabs(scEta()) <= 1.479 ){

    if(! (Full5x5_sigmaIetaIeta() < 0.0112) ) return false;
    if(! (fabs(dEtaSeed()) < 0.00377) ) return false;
    if(! (fabs(dPhiIn()) < 0.0884) ) return false;
    if(! (HoverE() < 0.05 + 1.16/scE() + 0.0324*Rho()/scE()) ) return false;
    if(! (fabs(InvEminusInvP()) < 0.193) ) return false;
    if(! (NMissingHits() <= 1) ) return false;
    if(! (PassConversionVeto()) ) return false;

    return true;

  }
  else{

    if(! (Full5x5_sigmaIetaIeta() < 0.0425) ) return false;
    if(! (fabs(dEtaSeed()) < 0.00674) ) return false;
    if(! (fabs(dPhiIn()) <  0.169 ) ) return false;
    if(! (HoverE() < 0.0441 + 2.54/scE() + 0.183*Rho()/scE()) ) return false;
    if(! (fabs(InvEminusInvP()) < 0.111) ) return false;
    if(! (NMissingHits() <= 1) ) return false;
    if(! (PassConversionVeto()) ) return false;

    return true;

  }

}

bool Electron::Pass_CutBasedVetoNoIso() const{
  
  if( fabs(scEta()) <= 1.479 ){
    
    if(! (Full5x5_sigmaIetaIeta() < 0.0126 ) ) return false;
    if(! (fabs(dEtaSeed()) < 0.00463 ) ) return false;
    if(! (fabs(dPhiIn()) < 0.148) ) return false;
    if(! (HoverE() < 0.05 + 1.16/scE() + 0.0324*Rho()/scE()) ) return false;
    if(! (fabs(InvEminusInvP()) < 0.209) ) return false;
    if(! (NMissingHits() <= 2) ) return false;
    if(! (PassConversionVeto()) ) return false;
    
    return true;
  
  }
  else{
    
    if(! (Full5x5_sigmaIetaIeta() < 0.0457) ) return false;
    if(! (fabs(dEtaSeed()) < 0.00814) ) return false;
    if(! (fabs(dPhiIn()) < 0.19) ) return false;
    if(! (HoverE() < 0.05 + 2.54/scE() + 0.183*Rho()/scE()) ) return false;
    if(! (fabs(InvEminusInvP()) < 0.132) ) return false;
    if(! (NMissingHits() <= 3) ) return false;
    if(! (PassConversionVeto()) ) return false;
    
    return true;
  
  }

}

bool Electron::Pass_CutBasedLoose() const{

  if( fabs(scEta()) <= 1.479 ){

    if(! (Full5x5_sigmaIetaIeta() < 0.0112) ) return false;
    if(! (fabs(dEtaSeed()) < 0.00377) ) return false;
    if(! (fabs(dPhiIn()) < 0.0884) ) return false;
    if(! (HoverE() < 0.05 + 1.16/scE() + 0.0324*Rho()/scE()) ) return false;
    if(! (RelIso() < 0.112+0.506/UncorrPt()) ) return false;
    if(! (fabs(InvEminusInvP()) < 0.193) ) return false;
    if(! (NMissingHits() <= 1) ) return false;
    if(! (PassConversionVeto()) ) return false;

    return true;

  }
  else{

    if(! (Full5x5_sigmaIetaIeta() < 0.0425) ) return false;
    if(! (fabs(dEtaSeed()) < 0.00674) ) return false;
    if(! (fabs(dPhiIn()) <  0.169 ) ) return false;
    if(! (HoverE() < 0.0441 + 2.54/scE() + 0.183*Rho()/scE()) ) return false;
    if(! (RelIso() < 0.108+0.963/UncorrPt()) ) return false;
    if(! (fabs(InvEminusInvP()) < 0.111) ) return false;
    if(! (NMissingHits() <= 1) ) return false;
    if(! (PassConversionVeto()) ) return false;

    return true;

  }

}

bool Electron::Pass_CutBasedVeto() const{

  if( fabs(scEta()) <= 1.479 ){

    if(! (Full5x5_sigmaIetaIeta() < 0.0126 ) ) return false;
    if(! (fabs(dEtaSeed()) < 0.00463 ) ) return false;
    if(! (fabs(dPhiIn()) < 0.148) ) return false;
    if(! (HoverE() < 0.05 + 1.16/scE() + 0.0324*Rho()/scE()) ) return false;
    if(! (RelIso() < 0.198+0.506/UncorrPt()) ) return false;
    if(! (fabs(InvEminusInvP()) < 0.209) ) return false;
    if(! (NMissingHits() <= 2) ) return false;
    if(! (PassConversionVeto()) ) return false;

    return true;

  }
  else{

    if(! (Full5x5_sigmaIetaIeta() < 0.0457) ) return false;
    if(! (fabs(dEtaSeed()) < 0.00814) ) return false;
    if(! (fabs(dPhiIn()) < 0.19) ) return false;
    if(! (HoverE() < 0.05 + 2.54/scE() + 0.183*Rho()/scE()) ) return false;
    if(! (RelIso() < 0.203+0.963/UncorrPt()) ) return false;
    if(! (fabs(InvEminusInvP()) < 0.132) ) return false;
    if(! (NMissingHits() <= 3) ) return false;
    if(! (PassConversionVeto()) ) return false;

    return true;

  }

}

//==== Customized IDs ====
bool Electron::Pass_fr_elec_tight() const {
	const double D2SIP = fabs(dXY() / dXYerr());
	if (!passMVAID_iso_WP90()) return false;
	if (! (RelIso() < 0.06)) return false;
	if (! (fabs(dXY()) < 0.025 && fabs(dZ()) < 0.1)) return false;
	if (! ((D2SIP < 4) && (dXYerr() != 0.))) return false;
	if (! PassConversionVeto()) return false;
	if (! (fabs(Eta()) < 2.5)) return false;

	// trigger emulation cuts
	if (Pt() > 15) {
		if (fabs(scEta()) < 1.479) {
			if (!(Full5x5_sigmaIetaIeta() < 0.012)) return false;
			if (!(fabs(dEtaSeed()) < 0.0095)) return false;
			if (!(fabs(dPhiIn()) < 0.065)) return false;
			if (!(HoverE() < 0.09)) return false;
			if (!(ecalPFClusterIso() < 0.37*Pt())) return false;
			if (!(hcalPFClusterIso() < 0.25*Pt())) return false;
			if (!(TrkIso() < 0.18*Pt())) return false;
		}

		else if (fabs(scEta()) < 2.5) {
			if (!(Full5x5_sigmaIetaIeta() < 0.033)) return false;
			if (!(HoverE() < 0.09)) return false;
			if (!(ecalPFClusterIso() < 0.45*Pt())) return false;
			if (!(hcalPFClusterIso() < 0.28*Pt())) return false;
			if (!(TrkIso() < 0.18*Pt())) return false;
		}
	}

	return true;
}

bool Electron::Pass_fr_elec_loose() const {
	const double D2SIP = fabs((dXY() / dXYerr()));
	//if (!Pass_FakeMVAWP("Loose")) return false;
	double fEta = fabs(Eta());
	if     (fEta<0.8  ){ if(MVAIso()<-0.92) return false; }
	else if(fEta<1.479){ if(MVAIso()<-0.88) return false; }
	else if(fEta<2.5  ){ if(MVAIso()<-0.78) return false; }

	if (! (RelIso() < 0.4)) return false;
    if (! ((fabs(dXY()) < 0.025) && (fabs(dZ()) < 0.1))) return false;
    if (! ((D2SIP < 4) && (dXYerr() != 0.))) return false;
    if (! PassConversionVeto()) return false;
	if (! (fabs(Eta()) < 2.5)) return false;

    // trigger emulation cuts
	// only for pt > 15 GeV for loose ID
    if (Pt() > 15) {
		if (fabs(scEta()) < 1.479) {
			if (!(Full5x5_sigmaIetaIeta() < 0.012)) return false;
			if (!(fabs(dEtaSeed()) < 0.0095)) return false;
			if (!(fabs(dPhiIn()) < 0.065)) return false;
			if (!(HoverE() < 0.09)) return false;
			if (!(ecalPFClusterIso() < 0.37*Pt())) return false;
			if (!(hcalPFClusterIso() < 0.25*Pt())) return false;
			if (!(TrkIso() < 0.18*Pt())) return false;
		}
    
		else if (fabs(scEta()) < 2.5) {
			if (!(Full5x5_sigmaIetaIeta() < 0.033)) return false;
			if (!(HoverE() < 0.09)) return false;
			if (!(ecalPFClusterIso() < 0.45*Pt())) return false;
			if (!(hcalPFClusterIso() < 0.28*Pt())) return false;
			if (!(TrkIso() < 0.18*Pt())) return false;   
		}
	}
	return true;
}

bool Electron::Pass_FakeMVAWP(TString wp) const {
  
  //cout << "[Electron::Pass_FakeMAVWP] need reoptimization for Fall17v2 dataset" << endl;

  double eta = fabs(Eta());
  double cutA = -999., cutB = -999., cutC = -999.;

  if (wp == "Tight") {
	cutA = 0.837; cutB = 0.715; cutC = 0.357;
	if (eta < 0.8) {
	  if (MVANoIso() > cutA) return true;
	  else return false;
	}
	else if (eta < 1.479) {
	  if (MVANoIso() > cutB) return true;
	  else return false;
	}
	else if (eta < 2.5) {
	  if (MVANoIso() > cutC) return true;
	  else return false;
    }
    else return false; 
  }

  else if (wp == "Loose") {
    cutA = -0.92; cutB = -0.88; cutC = -0.78;
    if (eta < 0.8) {
	  if (MVANoIso() > cutA) return true;
	  else return false;
    }
	else if (eta < 1.479) {
	  if (MVANoIso() > cutB) return true;
	  else return false;
	}
	else if (eta < 2.5) {
	  if (MVANoIso() > cutC) return true;
	  else return false;
    }
    else return false;
  }

  else {
	cout << "[Electron::Pass_FakeMVAWP] Wrong WP" << endl;
	exit(EXIT_FAILURE);
  }
}

bool Electron::Pass_FakeTight() const {
  // Electron_ID = "FakeTightID"
  // 2DSIP has not been set yet.... not important
  if (! Pass_FakeMVAWP("Tight")) return false;
  if (! (RelIso() < 0.06)) return false;
  if (! (fabs(dXY()) < 0.025 && fabs(dZ()) < 0.1)) return false;
  if (! PassConversionVeto()) return false;
  
  return true;
}

bool Electron::Pass_FakeLoose() const {
  // Electron_ID = "PassLooseID"
  // 2DSIP has not been set yet... not important
  if(! Pass_FakeMVAWP("Loose")) return false;
  if (! (RelIso() < 0.4)) return false;
  if (! (fabs(dXY()) < 0.025 && fabs(dZ()) < 0.1)) return false;
  if (! PassConversionVeto()) return false;

  return true;
}



void Electron::SetRho(double r){
  j_Rho = r;
}

void Electron::SetIsGsfCtfScPixChargeConsistent(bool b){
  j_isGsfCtfScPixChargeConsistent = b;
}
