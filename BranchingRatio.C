//  Intro

#include "TSystem.h"

#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes)
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#endif

#include <iostream>
#include <cmath>
#include <cstdio>

#include <TH1.h>
#include <TChain.h>

//  Ints to keep track of if decay is fully hardonic, semi-leptonic, or fully leptonic.
int full_hadron(0);
int semi_lept(0);
int full_lept(0);

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

void BranchingRatio(const char *inputFile) {
  gSystem -> Load("libDelphes");
  
  //  Create chain of ROOT trees.
  TChain chain("Delphes");
  chain.Add(inputFile);

  //  Create object of class ExRootTreeReader.
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader -> GetEntries();
  
  //  Get pointers to branches used in this analysis.
  TClonesArray *branchMuon = treeReader -> UseBranch("Muon");
  TClonesArray *branchJet = treeReader -> UseBranch("Jet");
  TClonesArray *branchElectron = treeReader -> UseBranch("Electron");
  TClonesArray *branchMET = treeReader -> UseBranch("MissingET");
  TClonesArray *branchParticle = treeReader -> UseBranch("Particle");

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
  
  //  Loop over all events.
  for(Int_t entry = 0; entry < numberOfEntries; ++entry) {
    //  Load selected branches with data from specified event.
    treeReader->ReadEntry(entry);

    //  Leptons are given a value of -1; hadrons are given a value of +1.
    //  If fully hadronic, type_sum = 4; fully leptonic, type_sum = 4.
    //  If semi-leptonic, -4 < type_sum < 4.
    int type_sum(0);

    //  Check to see if event is nonempty.
    if(branchParticle->GetEntries() > 0) {
      for (Int_t i = 0; i < branchParticle -> GetEntries(); ++i) {
	GenParticle *particle = (GenParticle*) branchParticle -> At(i);
       
	if (print) {
	  printParticle(particle, i);
	}
	
	//  Use only particles with mothers.
	if (particle -> M1 != -1) {
	  GenParticle *mom = (GenParticle*) branchParticle -> At(particle -> M1);
	  if ((abs(mom -> PID) == 24 && particle -> Status < 50)) {
	    final_PID[entry].push_back(particle -> PID);
	    if ((abs(particle -> PID)) < 9) {
	      type_sum++;
	    }
	    else if (abs(particle -> PID) > 9) {
	      type_sum--;
	    }
	  }
	}
	//  If all the daughter particles have been recovered,
	//  then stop searching through this event.
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

    //  Resent type_sum.
    type_sum = 0;
  }


  int total = full_hadron + full_lept + semi_lept;
  std::cout << "Fraction fully hadronic: " << full_hadron / (double) total <<std::endl;
  std::cout << "Fraction full leptonic: " << full_lept / (double) total << std::endl;
  std::cout << "Fraction semi-leptonic: " << semi_lept / (double) total << std::endl;

  double decay_hadron_obs = 0.6741;
  double decay_lept_obs = 0.1071 + 0.1063 + 0.1138;

  std::cout << "Observed fraction fully hadronic: " << decay_hadron_obs * decay_hadron_obs << std::endl;
  std::cout << "Observed fraction fully leptonic: " << decay_lept_obs * decay_lept_obs << std::endl;
  std::cout << "Observed fraction semi-leptonic: " << 2 * decay_hadron_obs * decay_lept_obs << std::endl;
}
