#include <iostream>
#include <fstream>
#include <string>
#include <thread>
#include "Pythia8/Pythia.h"
// #include "fastjet/PseudoJet.hh"
// #include "fastjet/ClusterSequence.hh"
#include "Pythia8Plugins/HepMC2.h"
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

// ROOT headers
#include "TH1.h"
#include "TFile.h"

// Own headers
#include "Simul.cc"
#include "PythiaProgramOpts.h"

using namespace std;
using namespace Pythia8;


/**
 * @brief Main function for generating MC events
 *
 * @param argc [description]
 * @param argv [description]
 *
 * @return [description]
 */
int main(int argc, char *argv[]) {

  PythiaProgramOpts opts(argc, argv);
  opts.printProgramOptions();

  ///////////////////////
  // SETUP PYTHIA HERE //
  ///////////////////////

  Pythia pythia;
  pythia.readString("Random:setSeed = on");
  pythia.readString("Random:seed = " + boost::lexical_cast<std::string>(opts.seed()));
  pythia.readString("Next:numberShowEvent = 00");

  /**
  * Setup pythia to produce gg->h(125)->AA,
  * and A->bb or A->tautau
  */
  pythia.readString("HiggsSM:gg2H = on");
  pythia.readString("25:m0 = 450.");
  for (int i = 0; i < 76; i++) {
    // turn off all h(125) decay modes
    stringstream SSt;
    SSt << i;
    pythia.readString("25:" + SSt.str() + ":onMode = off");
  }
  pythia.readString("25:addChannel = 1 1 100 23 36"); // h -> Z a
  // pythia.readString("25:addChannel = 1 1 100 36 36"); // h -> a a
  pythia.readString("36:m0 = " + boost::lexical_cast<std::string>(opts.mass()));
  pythia.readString("36:mMin = 3.5");
  pythia.readString("36:mWidth = 0.1");

  // Note that pythia auto rescales the total BR to 1.
  // So even if we add tautau wiht BR = 0.5, if this is the only channel it
  // will auto rescale to 1.
  if (opts.bbDecay())
    pythia.readString("36:addChannel = 1 0.5 100 5 -5");
  if  (opts.tautauDecay())
    pythia.readString("36:addChannel = 1 0.5 100 15 -15");
  pythia.readString("Beams:eCM = 13000.");

  // CUETP8M1 tune a la CMS
  pythia.readString("Tune:pp = 18");

  // pythia.readString("PartonLevel:all = off");
  // pythia.readString("PartonLevel:ISR = on");
  // pythia.readString("PartonLevel:FSR = on");

  pythia.init();

  // Any PYTHIA histograms here
  Hist hTransverseMomentum("h transverse momentum", 40, 0, 200);
  Hist a1Separation("a1 separation deltaR", 100, 0, 5);
  Hist a1PhiSeparation("a1 separation deltaPhi", 31, 0, 3.1);
  Hist a1Momentum("a1 momentum", 80, 0, 400);
  Hist a1TransverseMomentum("a1 transverse momentum", 80, 0, 400);
  Hist bbDR("a1 product separation deltaR", 100, 0, 5);
  Hist bbDPhi("a1 product separation deltaPhi", 31, 0, 3.1);

  // Interface for conversion from Pythia8::Event to HepMC event.
  HepMC::Pythia8ToHepMC ToHepMC;
  HepMC::IO_GenEvent ascii_io(opts.filenameHEPMC(), std::ios::out);
  if (opts.writeToHEPMC()) {
    cout << "Writing HepMC to " << opts.filenameHEPMC() << endl;
  }

  // Create an LHAup object that can access relevant information in pythia for writing to LHE
  LHAupFromPYTHIA8 * myLHA;
  if (opts.writeToLHE()) {
    cout << "Writing LHE to " << opts.filenameLHE() << endl;
    // Create an LHAup object that can access relevant information in pythia.
    myLHA = new LHAupFromPYTHIA8(&pythia.process, &pythia.info);
    // Open a file on which LHEF events should be stored, and write header.
    myLHA->openLHEF(opts.filenameLHE());
    // Store initialization info in the LHAup object.
    myLHA->setInit();
    // Write out this initialization info on the file.
    myLHA->initLHEF();
  }

  // Text file to write progress - handy for monitoring during jobs
  ofstream progressFile;
  std::string ext = ".hepmc";
  std::string stem = opts.filenameHEPMC().substr(0, opts.filenameHEPMC().size() - ext.size());
  progressFile.open(stem + "_progress.txt");

  //////////////////////
  // SETUP ROOT STUFF //
  //////////////////////
  TFile* outFile = new TFile("hist.root", "RECREATE");

  TH1F * h_bbDr = new TH1F("bbDr","bb DeltaR", 100, 0, 5);


  //////////////////////////
  // GENERATE EVENTS HERE //
  //////////////////////////

  int outputEvery = 50;
  for (int iEvent = 0; iEvent < opts.nEvents(); ++iEvent) {
    // output progress info
    if (iEvent % outputEvery == 0) {
      cout << "iEvent: " << iEvent << " - " << getCurrentTime() << endl;
      progressFile << "iEvent: " << iEvent << " - " << getCurrentTime() << endl;
    }

    // Generate event safely
    if (!pythia.next()) {
      break;
    }

    // Output to screen if wanted
    if (iEvent < 2 && opts.printEvent()) {
      pythia.info.list();
      pythia.event.list();
      pythia.process.list();
    }

    // find h decay products and look at separation
    for (int i = 0; i < pythia.event.size(); ++i) {

      if (abs(pythia.event[i].id()) == 25 && pythia.event[i].status() == -62) {

        // plot h1 daughter variables
        int d1 = pythia.event[i].daughter1();
        int d2 = pythia.event[i].daughter2();
        a1Separation.fill(REtaPhi(pythia.event[d1].p(), pythia.event[d2].p()));
        a1PhiSeparation.fill(phi(pythia.event[d1].p(), pythia.event[d2].p()));

        // figure out which is the a1
        int a1Ind = d1;
        if (pythia.event[d1].id() != 36) {
          a1Ind = d2;
        }
        // plot h1 variables
        hTransverseMomentum.fill(pythia.event[i].pT());

        // now plot a1 variables, and for its decay products
        a1Momentum.fill(pythia.event[a1Ind].pAbs());
        a1TransverseMomentum.fill(pythia.event[a1Ind].pT());

        // look at a1 daughter particles
        Vec4 daughter1Mom = pythia.event[pythia.event[a1Ind].daughter1()].p();
        Vec4 daughter2Mom = pythia.event[pythia.event[a1Ind].daughter2()].p();
        bbDR.fill(REtaPhi(daughter1Mom, daughter2Mom));
        h_bbDr->Fill(REtaPhi(daughter1Mom, daughter2Mom));
        bbDPhi.fill(phi(daughter1Mom, daughter2Mom));

      }
    }

    // Construct new empty HepMC event and fill it.
    // Write the HepMC event to file. Done with it.
    if (opts.writeToHEPMC()) {
      HepMC::GenEvent* hepmcevt = new HepMC::GenEvent(HepMC::Units::GEV, HepMC::Units::MM);
      ToHepMC.fill_next_event(pythia, hepmcevt);
      ascii_io << hepmcevt;
      delete hepmcevt;
    }

    if (opts.writeToLHE()) {
      // Store event info in the LHAup object.
      myLHA->setEvent();
      // Write out this event info on the file.
      // With optional argument (verbose =) false the file is smaller.
      myLHA->eventLHEF();
    }

  }

  progressFile.close();
  pythia.stat();

  cout << hTransverseMomentum << endl;
  cout << a1PhiSeparation << endl;
  cout << a1Separation << endl;
  cout << a1Momentum << endl;
  cout << a1TransverseMomentum << endl;
  cout << bbDR << endl;
  cout << bbDPhi << endl;

  if (opts.writeToLHE()) {
    // Update the cross section info based on Monte Carlo integration during run.
    myLHA->updateSigma();
    // Write endtag. Overwrite initialization info with new cross sections.
    myLHA->closeLHEF(true);
  }

  h_bbDr->Write();
  outFile->Close();
  delete outFile;

}
