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

TH1D* muonPt = new TH1D("muonPt", "Muon P_{T}; P_{T}; Frequency", 50, -10, 500);
TH1D* muonEta = new TH1D("muonEta", "Muon #eta; #eta; Frequency", 50, -2.5, 2.5);

TH1D* electronPt = new TH1D("electronPt", "Electron P_{T}; P_{T}; Frequency", 50, 0, 400);
TH1D* electronEta = new TH1D("electronEta", "Electron #eta; #eta; Frequency", 50, -2.5, 2.5);

TH1D* jetPt = new TH1D("jetMass", "Jet P_{T}; Jet P_{T}; Frequency", 50, -1, 450);
TH1D* jetEta = new TH1D("jetEta", "Jet Eta; Jet #eta; Frequency", 50, -5, 5);

TH1D* metMET = new TH1D("metMET", "Missing E_{T}; Missing E_{T}; Frequency", 50, -10, 300);
TH1D* metEta = new TH1D("metEta", "Missing E_{T} #eta; #eta; Frequency", 50, -8, 8);

TH1D* daughterType = new TH1D("daughterType", "Distribution of Daughter Particles; PID; Frequency", 80, -20, 20);
TH1D* hadronPt = new TH1D("hadronPt", "Daughter Hadron P_{T}; P_{T}; Frequency", 50, -10, 300);
TH1D* hadronEta = new TH1D("hadronEta", "Daughter Hadron #eta; #eta; Frequency", 50, -5, 5);
TH1D* leptonPt = new TH1D("leptonPt", "Daughter Lepton P_{T}; P_{T}; Frequency", 50, -5, 300);
TH1D* leptonEta = new TH1D("leptonEta", "Daughter Lepton #eta; #eta; Frequency", 50, -5, 5);

//  Prints particle information in a table.
void printParticle(GenParticle* particle, Int_t i) {
  printf("\n %3d %6d %4d %4d %4d %5d %6d %6d %3d",
	 i, particle -> PID, particle -> Status, particle -> IsPU, particle -> M1, particle -> M2,
	 particle -> D1, particle -> D2, particle -> Charge);
  printf("%9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %7.2f %7.2f %7.2f %7.2f %7.2f",
	 particle -> Mass,
	 particle -> E,particle -> Px, particle -> Py, particle -> Pz,
	 particle -> P, particle -> PT, particle -> Eta, particle -> Phi,
	 particle -> Rapidity, particle -> CtgTheta,
	 particle -> D0, particle -> DZ,
	 particle -> T, particle -> X, particle -> Y, particle -> Z);
}

//  Fills muons histograms.
void muonHist(TClonesArray* branchMuon) {
  if (branchMuon -> GetEntries() > 0) {
    for (Int_t i = 0; i < branchMuon -> GetEntries(); ++i) {
      Muon* muon = (Muon*) branchMuon -> At(i);
      muonPt -> Fill(muon -> PT);
      muonEta -> Fill(muon -> Eta);
    }
  }
}

//  Fills electron histograms.
void electronHist(TClonesArray* branchElectron) {
  if (branchElectron -> GetEntries() > 0) {
    for (Int_t i = 0; i < branchElectron -> GetEntries(); ++i) {
      Electron* electron = (Electron*) branchElectron -> At(i);
      electronPt -> Fill(electron -> PT);
      electronEta -> Fill(electron -> Eta);
    }
  }
}

//  Fills jet histograms.
void jetHist(TClonesArray* branchJet) {
  if (branchJet -> GetEntries() > 0) {
    for (Int_t i = 0; i < branchJet -> GetEntries(); ++i) {
      Jet* jet = (Jet*) branchJet -> At(i);
      jetPt -> Fill(jet -> PT);
      jetEta -> Fill(jet -> Eta);
    }
  }
}

//  Fills MET histograms.
void metHist(TClonesArray* branchMET) {
  if (branchMET -> GetEntries() > 0) {
    for (Int_t i = 0; i < branchMET -> GetEntries(); ++i) {
      MissingET* met = (MissingET*) branchMET -> At(i);
      metMET -> Fill(met -> MET);
      metEta -> Fill(met -> Eta);
    }
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

  //  Bool to enable printing for particles.
  bool print(0);
  
  if (print) {
    printf("\n \n %6s %3s %4s %1s %3s %4s %6s %6s %6s",
	   "Index", "PID", "Status", "IsPU", "M1", "M2", "D1", "D2", "Charge");
    printf("%6s %7s %9s %10s %8s %10s %9s %8s %10s %9s %12s %6s %7s %6s %7s %7s %7s",
	   "Mass", "E", "Px", "Py", "Pz", "P", "PT", "Eta", "Phi", "Rapid", "CtgTheta",
	   "D0", "DZ", "T", "X", "Y", "Z");
  std:cout << "\n--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------";
  }
  
  //  Vector to store PIDs of all W+- daughter particles of each event.
  std::vector<std::vector<Int_t>> final_PID(numberOfEntries);

  //  Ints to track if decay is fully hardonic, semi-leptonic, or fully leptonic.
  int full_hadron(0);
  int semi_lept(0);
  int full_lept(0);
  
  //  Loop over all events.
  for(Int_t entry = 0; entry < numberOfEntries; ++entry) {
    //  Load selected branches with data from specified event.
    treeReader->ReadEntry(entry);

    //  Fill muon pT and eta histograms.
    muonHist(branchMuon);

    //  Fill electron pT and eta histograms.
    electronHist(branchElectron);

    //  Fill jet pT and eta.
    jetHist(branchJet);

    //  Fill MET pT and eta.
    metHist(branchMET);

    //  Leptons are given a value of -1; hadrons are given a value of +1.
    //  If fully hadronic, type_sum = 4; fully leptonic, type_sum = 4.
    //  If semi-leptonic, -4 < type_sum < 4.
    int type_sum(0);
    
    //  Look in Particle branch.
    if (branchParticle -> GetEntries() > 0) {
      for (Int_t i = 0; i < branchParticle -> GetEntries(); ++i) {
	GenParticle *particle = (GenParticle*) branchParticle -> At(i);
       
	if (print) {
	  printParticle(particle, i);
	}
	
	//  Use only particles with mothers.
	if (particle -> M1 != -1) {
	  GenParticle *mom = (GenParticle*) branchParticle -> At(particle -> M1);

	  //  Get particles with W boson mom and that are not produced by final state showers.
	  //  Neglecting the latter bool would yield a W boson as a daughter.
	  if ((abs(mom -> PID) == 24 && particle -> Status < 50)) {
	    final_PID[entry].push_back(particle -> PID);
	    daughterType -> Fill(particle -> PID);

	    //  Hadron daughters.
	    if ((abs(particle -> PID)) < 9) {
	      type_sum++;
	      hadronPt -> Fill(particle -> PT);
	      hadronEta -> Fill(particle -> Eta);
	    }
	    //  Lepton daughters.
	    else if (abs(particle -> PID) > 9) {
	      type_sum--;
	      leptonPt -> Fill(particle -> PT);
	      leptonEta -> Fill(particle -> Eta);
	    }
	  }
	}

	//  If all the daughter particles have been recovered,
	//  then stop searching through this event to suppress bottlenecking.
	if (final_PID[entry].size() == 4) {
	  break;
	}
      }
    }

    //  Classify event.
    if (type_sum == 4) {
      full_hadron++;
    }
    else if (type_sum == -4) {
      full_lept++;
    }
    else {
      semi_lept++;
    }

    //  Reset type_sum for next event.
    type_sum = 0;
  }

  //  Save directory
  TString directory = "~/Work/MG-Pyth-Delphes/";

  //  Save histograms.
  TCanvas* c = new TCanvas("c");
  muonPt -> Draw();
  c -> SaveAs(directory + "muonPt.png");
  muonEta -> Draw();
  c -> SaveAs(directory + "muonEta.png");
  electronPt -> Draw();
  c -> SaveAs(directory + "electronPt.png");
  electronEta -> Draw();
  c -> SaveAs(directory + "electronEta.png");
  jetPt -> Draw();
  c -> SaveAs(directory + "jetPt.png");
  jetEta -> Draw();
  c -> SaveAs(directory + "jetEta.png");
  metMET -> Draw();
  c -> SaveAs(directory + "metMET.png");
  metEta -> Draw();
  c -> SaveAs(directory + "metEta.png");
  daughterType -> Draw();
  c -> SaveAs(directory + "DaughterType.png");
  hadronPt -> Draw();
  c -> SaveAs(directory + "hadronPt.png");
  hadronEta -> Draw();
  c -> SaveAs(directory + "hadronEta.png");
  leptonPt -> Draw();
  c -> SaveAs(directory + "leptonPt.png");
  leptonEta -> Draw();
  c -> SaveAs(directory + "leptonEta.png");

  //  Print fractions.
  std::cout << numberOfEntries << std::endl;
  std::cout << "Simulated fraction fully hadronic: " << full_hadron / (double)  numberOfEntries <<std::endl;
  std::cout << "Simulated fraction fraction full leptonic: " << full_lept / (double) numberOfEntries << std::endl;
  std::cout << "Simulated fraction semi-leptonic: " << semi_lept / (double) numberOfEntries << std::endl;

  //  Experimental observations of W decay from PDG 2016.
  double decay_hadron_obs = 0.6741;
  double dho_error = 0.0027;
  double decay_lept_obs = 0.1071 + 0.1063 + 0.1138;
  double dlo_error = 0.0016 + 0.0015 + 0.0021;

  printObservation(decay_hadron_obs, decay_lept_obs, 0);
  printObservation(decay_hadron_obs - dho_error, decay_lept_obs - dlo_error, 1);
  printObservation(decay_hadron_obs + dho_error, decay_lept_obs + dlo_error, 2);
}
