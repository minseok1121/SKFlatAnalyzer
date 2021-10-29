#include "TutorialRun.h"

//==== Constructor and Destructor
TutorialRun::TutorialRun() {}
TutorialRun::~TutorialRun() {}

//==== Initialize variables
// NOTE: Every local varaible declared in executeEvent() will be initialized for every event, 
// but some variables do not have to be reinitialized for every event.
// for example, we will use the same trigger and trigger safe pt cut throughout all events.
// For these varibles(which are called global variables) are declared in the TutorialRun.h
// and initialized in initializeAnalyzer step
void TutorialRun::initializeAnalyzer(){
		//==== Dimuon Z-peak events with two muons IDs
		// One can define customized Muon IDs in DataFrommat/src/Muon.C
		// For this tutorial, let's use POG based IDs
		// https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonSelection
		MuonIDs = {"POGMedium", "POGTight"};
		MuonVetoID = "POGLoose";
		ElectronVetoID = "passLooseID"; // to ensure there is no electrons

		//==== Trigger settings
		// In this tutorial, we will use HLT_IsoMu27_v (HighLevelTrigger_IsolatedMuon ptcut 27)
		// which is the recommended singlemuon tirgger for year 2017
		// https://twiki.cern.ch/twiki/bin/view/CMS/MuonHLT2017
		// NOTE: for each dataset, you should use corresponding triggers
		// for example, SingleMuon -> HLT_IsoMu27_v
		//              DoubleMuon -> HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v
		if (DataYear == 2017) {
				IsoMuTriggerName = "HLT_IsoMu27_v";
				TriggerSafePtCut = 29.;
		}
		else {
				cerr << "[TutorialRun::InitializeAnalyzer] This tutorail is for year 2017 only" << endl;
				exit(EXIT_FAILURE);
		}

		cout << "[TutorialRun::initailizeAnalyzer] IsoMuTriggerName = " << IsoMuTriggerName << endl;
		cout << "[TutorailRun::initializeAnalyzer] TriggerSafePtCut = " << TriggerSafePtCut << endl;

}

void TutorialRun::executeEvent(){
		//==== *IMPORTANT TO SAVE CPU TIME*
		//==== Every GetMuon() function first collect ALL MINIAOD muons with GetAllMuons() 
		//==== and then check ID booleans
		AllMuons = GetAllMuons();
		AllElectrons = GetAllElectrons();
		AllJets = GetAllJets();

		for (const auto &MuonID: MuonIDs) {
				//==== No cut
				//==== FillHist function is defined in Analyzers/src/AnalyzerCore.C
				//==== void FillHist(histkey, value, weight, nbins, left_edge, right_edge);
				FillHist(MuonID+"/Cutflow", 0., 1., 10, 0., 10.);

				//==== MET (Missing Transverse Energy) Filter
				if (! PassMETFilter()) return;
				FillHist(MuonID+"/Cutflow", 1., 1., 10, 0., 10.);

				//==== Trigger
				//==== check DataFormat/src/Event.C
				Event ev = GetEvent();
				if (! (ev.PassTrigger(IsoMuTriggerName))) return;
				FillHist(MuonID+"/Cutflow", 2., 1., 10, 0., 10.);

				//==== Object definition
				//==== muon/electron/jet candidates are in AllMuons, AllElectrons, AllJets
				//==== select good quality particles with certain IDs
				vector<Muon> muons = SelectMuons(AllMuons, MuonID, 30., 2.4);
				vector<Muon> muons_veto = SelectMuons(AllMuons, MuonVetoID, 25., 2.4);
				vector<Electron> electrons_veto = SelectElectrons(AllElectrons, ElectronVetoID, 25., 2.5);
				vector<Jet> jets = SelectJets(AllJets, "tight", 30., 2.4);
				// lepton cleaning, to avoid double counting between leptons and jets
				jets = JetsVetoLeptonInside(jets, electrons_veto, muons_veto, 0.4);
				const Particle METv = ev.GetMETVector();

				//==== sort in Pt order
				std::sort(muons.begin(), muons.end(), PtComparing);
				std::sort(jets.begin(), jets.end(), PtComparing);

				//==== b-tagging
				//==== in this tutorial, we will use DeepJet algorithm
				//==== see https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation
				vector<Jet> bjets, non_bjets;
				for (const auto &jet: jets) {
						const double this_discr = jet.GetTaggerResult(JetTagging::DeepJet);
						if (this_discr > mcCorr->GetJetTaggingCutValue(JetTagging::DeepJet, JetTagging::Medium))
								bjets.emplace_back(jet);
						else
								non_bjets.emplace_back(jet);
				}

				//==== Event Selection
				//==== Problem 1. Complete the cutflow!
				//==== OS dimuon, no additional veto leptons
				if (! (muons.size() == 2 && muons_veto.size() == 2)) return;
				if (! (muons.at(0).Charge() + muons.at(1).Charge() == 0)) return;

				//==== leading muon should pass trigger-safe pt cut
				if (! (muons.at(0).Pt() > TriggerSafePtCut)) return;

				//==== On-Z condition
				//==== Reconstruct Z candidates using two muons
				//==== the mass of reconstructed Z candidates do not have to be
				//==== exact the same with Z mass, since it's not a real particle
				//==== so let's give a Z-mass window within 15 GeV
				// Question: what is the definition of adding two partilces?
				const Particle ZCand = muons.at(0) + muons.at(1);
				const double mZ = 91.2; // GeV
				if (! (fabs(mZ - ZCand.M()) < 15.)) return;

				//==== For MC, the events should be normalized to the xsec * lumi
				double weight = 1.;
				if (!IsDATA) {
						const double gen_weight = ev.MCweight()*weight_norm_1invpb;
						const double lumi = ev.GetTriggerLumi("Full");
						weight *= gen_weight*lumi;
				}

				//==== Now fill histograms
				FillHist(MuonID+"/ZCand/mass", ZCand.M(), weight, 40, 70., 110.);
				FillHist(MuonID+"/ZCand/pt", ZCand.Pt(), weight, 300, 0., 300.);
				FillHist(MuonID+"/ZCand/eta", ZCand.Eta(), weight, 100, -5., 5.);
				FillHist(MuonID+"/ZCand/phi", ZCand.Phi(), weight, 64, -3.2, 3.2);
				FillHist(MuonID+"/MET", METv.Pt(), weight, 300, 0., 300.);
				FillHist(MuonID+"/bjets/size", bjets.size(), weight, 20, 0., 20.);
				//==== Problem 2. Let's see the kinematic distributions of decayed products. 
				//==== Fill histograms for the leading & subleading muons.
				
				//==== Problem 3. After you compared the data and MCs, 
				//==== you can check many kinematic distributions of DY(Drell-Yan) and other sources.
				//==== Some variables can be used to suppress non-DY backgrounds (such as ttbar process)
				//==== implement additional cuts to suppress backgrounds further.	
		}
}
