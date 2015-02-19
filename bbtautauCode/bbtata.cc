#include <iostream>
#include <fstream>
#include <string>
#include <thread>
#include "Pythia8/Pythia.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
/* for Pythia 8.2X */
#include "Pythia8Plugins/HepMC2.h"
/* These 3 lines for Pythia 8.1X */
// #include "Pythia8/Pythia8ToHepMC.h"
// #include "HepMC/GenEvent.h"
// #include "HepMC/IO_GenEvent.h"
// #include <boost/lexical_cast.hpp>

#include "Simul.cc"
#include "PythiaProgramOpts.h"

using namespace std;
using namespace Pythia8;

int main(int argc, char *argv[]) {

  PythiaProgramOpts opts(argc, argv);
  opts.printProgramOptions();

  Simul simul(opts);
  simul.run(opts.nEvents());

}
