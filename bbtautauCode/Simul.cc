// #include <iostream>
// #include <fstream>
// #include <string>
// #include "Pythia.h"
// #include "MyEvent.cc"
using namespace std;
using namespace Pythia8;

class Simul {
private:
  // Interface for conversion from Pythia8::Event to HepMC event.
  HepMC::Pythia8ToHepMC ToHepMC;
  // Specify file where HepMC events will be stored.
  HepMC::IO_GenEvent ascii_io;
  string file = "test.txt";
  Pythia pythia;
  MyEvent *event;
  MyProcess *process;
  bbmuHist *bbmudata;
  bbtauHist *bbtaudata;
  void bbmuanalyse();
  void bbtauanalyse();
public:
  Simul(string model);
  Simul(string model, int back);
  Simul(string model, string ptmin, string ptmax);
  ~Simul();
  void run(int nEvnt);
};

Simul::Simul(string model): ascii_io("hepmcout41.dat", std::ios::out) {
  pythia.readString("HiggsSM:gg2H = on");
  pythia.readString("25:m0 = 125.");
  for (int i = 0; i < 76; i++) {
    stringstream SSt;
    SSt << i;
    pythia.readString("25:" + SSt.str() + ":onMode = off");
  }
  pythia.readString("25:addChannel = 1 .05 100 36 36");
  pythia.readString("36:m0 = " + model + ".0");
  pythia.readString("36:mMin = 9.5");
  pythia.readString("36:mWidth = 0.1");
  pythia.readString("36:oneChannel = 1 0.5 100 5 -5");
  pythia.readString("36:addChannel = 1 0.5 100 15 -15");
  pythia.readString("36:1:bRatio = 0.5");
  pythia.readString("36:2:bRatio = 0.5");
  //  pythia.readString("PartonLevel:all = off");    
  //pythia.readString("PartonLevel:ISR = on");    
  //pythia.readString("PartonLevel:FSR = on");    
  pythia.readString("Beams:eCM = 14000.");
  //  pythia.readString("PartonLevel:all = off");
  pythia.init();
  bbmudata = new bbmuHist(model, 0.25 * 0.25 * 2);
  bbtaudata = new bbtauHist(model, 0.25 * 0.25 * 2);
}

Simul::~Simul() {
  delete bbmudata;
  delete bbtaudata;
}

void Simul::run(int nEvnt) {
  for (int iEvent = 0; iEvent < nEvnt; ++iEvent) {
    cout << "Event start" << endl;
    if (!pythia.next()) {
      break;
    }
    if (iEvent < 2) {
      pythia.info.list();
      pythia.event.list();
      pythia.process.list();
    }
    cout << "Event end" << endl;
    // process = new MyProcess(&pythia);
    cout << "Event end2" << endl;

    // event = new MyEvent(&pythia, 0.4);
    cout << "Event end3" << endl;
    // bbmuanalyse();
    // bbtauanalyse();
    cout << "Event end4" << endl;
    // delete event;
    // delete process;

    // Construct new empty HepMC event and fill it.
    // Units will be as chosen for HepMC build; but can be changed
    // by arguments, e.g. GenEvt( HepMC::Units::GEV, HepMC::Units::MM)
    HepMC::GenEvent* hepmcevt = new HepMC::GenEvent( HepMC::Units::GEV, HepMC::Units::MM);
    cout << "Event end5" << endl;
    ToHepMC.fill_next_event(pythia, hepmcevt );
    cout << "Event end6" << endl;

    // Write the HepMC event to file. Done with it.
    ascii_io << hepmcevt;
    delete hepmcevt;
  }
  // ofstream OutFile;
  // OutFile.open(file.c_str());
  // OutFile << bbmudata->GetN() << endl;
  // OutFile << bbtaudata->GetN() << endl;
  // OutFile.close();

  // cout << bbmudata->GetN() << endl;
  // cout << bbtaudata->GetN() << endl;
  // bbmudata->Save();
  // bbtaudata->Save();
  pythia.stat();
}

void Simul::bbmuanalyse() {
  //  bbmudata->AddEvent();
  bbmudata->Fill(process, event);
}

void Simul::bbtauanalyse() {
  //  Zbbdata->AddEvent();
  bbtaudata->Fill(process, event);
}

