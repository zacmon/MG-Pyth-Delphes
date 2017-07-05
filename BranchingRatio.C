//////////////////////////////////////////////////////////////////////////////
///                  Delphes Tutorial w/ Branching Ratio                   ///
///                  Zachary Montague 27.06.17                             ///
///                                                                        ///
///  This takes a Madgraph, Pythia, Delphes simulated pp to ttbar          ///
///  and analyzes the branching ratios of the fully hadronic,              ///
///  semi-leptonic, and fully leptonic decays.                             ///
///  PYTHIA Status codes:                                                  ///
///  http://home.thep.lu.se/~torbjorn/pythia81html/ParticleProperties.html ///
///  PDG PID:                                                              ///
///  http://pdg.lbl.gov/2007/reviews/montecarlorpp.pdf                     ///
/////////////////////////////////////////////////////////////////////////////

#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes)
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#endif

#include <iostream>
#include <cmath>
#include <cstdio>

#include <TSystem.h>
#include <TH1D.h>
#include <TChain.h>
#include <TString.h>

//  Delcare histograms.
TH1D* muonPtHist = new TH1D("muonPt", "Muon P_{T}; P_{T}; Frequency", 50, -10, 500);
TH1D* muonEtaHist = new TH1D("muonEta", "Muon #eta; #eta; Frequency", 50, -2.5, 2.5);

TH1D* electronPtHist = new TH1D("electronPt", "Electron P_{T}; P_{T}; Frequency", 50, 0, 400);
TH1D* electronEtaHist = new TH1D("electronEta", "Electron #eta; #eta; Frequency", 50, -2.5, 2.5);

TH1D* jetPtHist = new TH1D("jetMass", "Jet P_{T}; Jet P_{T}; Frequency", 50, -1, 450);
TH1D* jetEtaHist = new TH1D("jetEta", "Jet Eta; Jet #eta; Frequency", 50, -5, 5);

TH1D* metMETHist = new TH1D("metMET", "Missing E_{T}; Missing E_{T}; Frequency", 50, -10, 300);
TH1D* metEtaHist = new TH1D("metEta", "Missing E_{T} #eta; #eta; Frequency", 50, -8, 8);

TH1D* daughterPIDHist = new TH1D("daughterPID", "Distribution of Daughter Particles; PID; Frequency", 80, -20, 20);
TH1D* hadronPtHist = new TH1D("hadronPt", "Daughter Hadron P_{T}; P_{T}; Frequency", 50, -10, 300);
TH1D* hadronEtaHist = new TH1D("hadronEta", "Daughter Hadron #eta; #eta; Frequency", 50, -5, 5);
TH1D* leptonPtHist = new TH1D("leptonPt", "Daughter Lepton P_{T}; P_{T}; Frequency", 50, -5, 300);
TH1D* leptonEtaHist = new TH1D("leptonEta", "Daughter Lepton #eta; #eta; Frequency", 50, -5, 5);

//  Declare histogram array.
TH1D* histogramArray[] = {muonPtHist, muonEtaHist, electronPtHist, jetPtHist, jetEtaHist, metMETHist, metEtaHist,
			  daughterPIDHist, hadronPtHist, hadronEtaHist, leptonPtHist, leptonEtaHist};

//  Delcare save directory.
TString saveDirectory = "~/Work/MG-Pyth-Delphes/";

void printParticleInfo(GenParticle* particle, Int_t index) {
  printf("\n %3d %6d %4d %4d %4d %5d %6d %6d %3d",
	 index, particle -> PID, particle -> Status, particle -> IsPU, particle -> M1, particle -> M2,
	 particle -> D1, particle -> D2, particle -> Charge);
  printf("%9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %7.2f %7.2f %7.2f %7.2f %7.2f",
	 particle -> Mass,
	 particle -> E,particle -> Px, particle -> Py, particle -> Pz,
	 particle -> P, particle -> PT, particle -> Eta, particle -> Phi,
	 particle -> Rapidity, particle -> CtgTheta,
	 particle -> D0, particle -> DZ,
	 particle -> T, particle -> X, particle -> Y, particle -> Z);
}

void fillMuonHistograms(TClonesArray* branchMuon) {
  if (branchMuon -> GetEntries() > 0) {
    for (Int_t i = 0; i < branchMuon -> GetEntries(); ++i) {
      Muon* muon = (Muon*) branchMuon -> At(i);
      muonPtHist -> Fill(muon -> PT);
      muonEtaHist -> Fill(muon -> Eta);
    }
  }
}

void fillElectronHistograms(TClonesArray* branchElectron) {
  if (branchElectron -> GetEntries() > 0) {
    for (Int_t i = 0; i < branchElectron -> GetEntries(); ++i) {
      Electron* electron = (Electron*) branchElectron -> At(i);
      electronPtHist -> Fill(electron -> PT);
      electronEtaHist -> Fill(electron -> Eta);
    }
  }
}

void fillJetHistograms(TClonesArray* branchJet) {
  if (branchJet -> GetEntries() > 0) {
    for (Int_t i = 0; i < branchJet -> GetEntries(); ++i) {
      Jet* jet = (Jet*) branchJet -> At(i);
      jetPtHist -> Fill(jet -> PT);
      jetEtaHist -> Fill(jet -> Eta);
    }
  }
}

void fillMETHistograms(TClonesArray* branchMET) {
  if (branchMET -> GetEntries() > 0) {
    for (Int_t i = 0; i < branchMET -> GetEntries(); ++i) {
      MissingET* met = (MissingET*) branchMET -> At(i);
      metMETHist -> Fill(met -> MET);
      metEtaHist -> Fill(met -> Eta);
    }
  }
}

bool hasMother(GenParticle* particle) {
  return particle -> M1 != -1;
}  

bool isPhoton(GenParticle* particle) {
  return particle -> PID == 22;
}

bool isWBoson(GenParticle* particle) {
  return abs(particle -> PID) == 24;
}

bool hasDaughters(TClonesArray* branchParticle, GenParticle* particle) {
  GenParticle* daughter1 = (GenParticle*) branchParticle -> At(particle -> D1);
  return abs(daughter1 -> PID) != 24;
}

int getDaughterPID(TClonesArray* branchParticle, int index) {
  GenParticle* daughter = (GenParticle*) branchParticle -> At(index);
  return daughter -> PID;
}

bool isQuark(int PID) {
  return abs(PID) < 6;
}

bool isLepton(int PID) {
  return abs(PID) > 10 && abs(PID) < 17;
}

void saveHistograms() {
  TCanvas* c = new TCanvas("c", "c");
  for (auto histogram : histogramArray) {
    histogram -> Draw();
    c -> SaveAs(saveDirectory + histogram -> GetName() + ".png");
  }
}

//  Print observed branching ratios.
//  int type corresponds to 0 for measurement,
//  1 for lower bound measurement, 2 for upper
//  bound measurement.
void printObservation(double hadron, double lepton, int type) {
  double full_hadron = hadron * hadron;
  double full_lepton = lepton * lepton;
  double semi_lepton = 2 * hadron * lepton;
  
  TString string_fh = "\nObserved fully hadronic";
  TString string_fl = "Observed fully leptonic";
  TString string_sl = "Observed semi-leptonic";
  TString endstring;

  if (type == 0) {
    endstring = ": ";
  }
  else if (type == 1) {
    endstring = " lower bound: ";
  }
  else if (type == 2) {
    endstring = " upper bound: ";
  }
    
  std::cout << string_fh << endstring << full_hadron << std::endl;
  std::cout << string_fl << endstring << full_lepton << std::endl;
  std::cout << string_sl << endstring << semi_lepton << std::endl;
}

//  Main function.
void BranchingRatio(const char *inputFile) {
  gSystem -> Load("libDelphes");
  
  //  Create chain of ROOT trees.
  TChain chain("Delphes");
  chain.Add(inputFile);

  //  Create object of class ExRootTreeReader.
  ExRootTreeReader* treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader -> GetEntries();
  
  //  Get pointers to branches used in this analysis.
  TClonesArray* branchParticle = treeReader -> UseBranch("Particle");
  TClonesArray* branchMuon = treeReader -> UseBranch("Muon");
  TClonesArray* branchElectron = treeReader -> UseBranch("Electron");
  TClonesArray* branchJet = treeReader -> UseBranch("Jet");
  TClonesArray* branchMET = treeReader -> UseBranch("MissingET");

  //  Bool to enable printing of particle info.
  bool printParticle(0);
  
  if (printParticle) {
    printf("\n \n %6s %3s %4s %1s %3s %4s %6s %6s %6s",
	   "Index", "PID", "Status", "IsPU", "M1", "M2", "D1", "D2", "Charge");
    printf("%6s %7s %9s %10s %8s %10s %9s %8s %10s %9s %12s %6s %7s %6s %7s %7s %7s",
	   "Mass", "E", "Px", "Py", "Pz", "P", "PT", "Eta", "Phi", "Rapid", "CtgTheta",
	   "D0", "DZ", "T", "X", "Y", "Z");
  std:cout << "\n--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------";
  }

  //  Ints to track if decay is fully hardonic, semi-leptonic, or fully leptonic.
  int numFullyHadronic(0);
  int numSemiLeptonic(0);
  int numFullyLeptonic(0);
  
  for(Int_t entry = 0; entry < numberOfEntries; ++entry) {
    treeReader -> ReadEntry(entry);

    fillMuonHistograms(branchMuon);
    fillElectronHistograms(branchElectron);
    fillJetHistograms(branchJet);
    fillMETHistograms(branchMET);

    //  Look in Particle branch.
    if (branchParticle -> GetEntries() > 0) {
      //  Leptons are given a value of -1; hadrons are given a value of +1.
      //  If fully hadronic, type_sum = 4; fully leptonic, type_sum = 4.
      //  If semi-leptonic, -4 < type_sum < 4.
      int decayTypeSum(0);

      for (Int_t i = 0; i < branchParticle -> GetEntries(); ++i) {
	GenParticle *particle = (GenParticle*) branchParticle -> At(i);
       
	if (printParticle) {
	  printParticleInfo(particle, i);
	}
	//  Get W bosons that are not produced by final state showers.
	//  Neglecting the latter bool would yield a W boson as a daughter.
	if (isWBoson(particle) && hasDaughters(branchParticle, particle)) {
	  GenParticle* daughterParticle1 = (GenParticle*) branchParticle -> At(particle -> D1);
	  GenParticle* daughterParticle2 = (GenParticle*) branchParticle -> At(particle -> D2);
	  daughterPIDHist -> Fill(daughterParticle1 -> PID);
	  daughterPIDHist -> Fill(daughterParticle2 -> PID);
	  
	  if (isQuark(daughterParticle1 -> PID)) {
	    decayTypeSum++;
	    hadronPtHist -> Fill(particle -> PT);
	    hadronEtaHist -> Fill(particle -> Eta);
	  }
	  
	  else if (isLepton(daughterParticle1 -> PID)) {
	    decayTypeSum = decayTypeSum + 3;
	    leptonPtHist -> Fill(particle -> PT);
	    leptonEtaHist -> Fill(particle -> Eta);
	  }
	}
	
	//  If all the daughter particles have been recovered,
	//  then stop searching through this event to suppress bottlenecking.
	if (decayTypeSum == 4) { 
	  numFullyHadronic++;
	  break;
	}
	else if (decayTypeSum == 12) {
	  numFullyLeptonic++;
	  break;
	}
	else if (decayTypeSum == 8) {
	  numSemiLeptonic++;
	  break;
	}
      }
    }
  }

  //  Save histograms.
  saveHistograms();

  //  Print fractions.
  std::cout << "Simulated fraction fully hadronic: " << numFullyHadronic / (double)  numberOfEntries <<std::endl;
  std::cout << "Simulated fraction fraction full leptonic: " << numFullyLeptonic / (double) numberOfEntries << std::endl;
  std::cout << "Simulated fraction semi-leptonic: " << numSemiLeptonic / (double) numberOfEntries << std::endl;

  //  Experimental observations of W decay from PDG 2016.
  double decay_hadron_obs = 0.6741;
  double dho_error = 0.0027;
  double decay_lept_obs = 0.1071 + 0.1063 + 0.1138;
  double dlo_error = 0.0016 + 0.0015 + 0.0021;

  printObservation(decay_hadron_obs, decay_lept_obs, 0);
  printObservation(decay_hadron_obs - dho_error, decay_lept_obs - dlo_error, 1);
  printObservation(decay_hadron_obs + dho_error, decay_lept_obs + dlo_error, 2);
}
