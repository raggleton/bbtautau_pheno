#include <iostream>
#include <fstream>
#include <string>
#include <thread>
#include "Pythia8/Pythia.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
// #include "Pythia8Plugins/HepMC2.h"
#include "Pythia8/Pythia8ToHepMC.h"
#include "HepMC/GenEvent.h"
#include "HepMC/IO_GenEvent.h"
#include "MyJet.cc"
#include "MyEvent.cc"
#include "MyProcess.cc"
#include "bbmuHist.cc"
#include "bbtauHist.cc"
#include "Simul.cc"
using namespace std;
using namespace Pythia8;

void runSim(int ma1) {
  Simul *simul;
  stringstream Sa1;
  Sa1 << ma1;
  simul = new Simul(Sa1.str());
  simul->run(50);
  delete simul;
}

int main(int argc, char *argv[]) {
  vector < thread > tradar;
  tradar.resize(0);
  // Run over diff masses of a1
  // for (int a1 = 15; a1 < 64; a1++){
  //   if (tradar.size() < 30) {
  //     tradar.push_back(thread(runSim, a1));
  //   }else{
  //     for (int i = 0; i < tradar.size(); i++){tradar[i].join();}
  //     tradar.resize(0);
  //     tradar.push_back(thread(runSim, a1));
  //   }
  // }
  // for (int i = 0; i < tradar.size(); i++){tradar[i].join();}
  runSim(15);
}
