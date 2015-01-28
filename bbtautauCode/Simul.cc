// #include <iostream>
// #include <fstream>
// #include <string>
// #include "Pythia.h"
// #include "MyEvent.cc"
using namespace std;
using namespace Pythia8; 

class Simul{
private:
  string file;
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

Simul::Simul(string model){
  pythia.readString("HiggsSM:gg2H = on");    
  pythia.readString("25:m0 = 125.");
  for (int i = 0; i < 76; i++){
    stringstream SSt;
    SSt << i;
    pythia.readString("25:"+SSt.str()+":onMode = off");
  }
  pythia.readString("25:addChannel = 1 .05 100 36 36");
  pythia.readString("36:m0 = "+model+".0");
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
  bbmudata = new bbmuHist(model,0.25*0.25*2);
  bbtaudata = new bbtauHist(model,0.25*0.25*2);
}


Simul::~Simul(){
  delete bbmudata;
  delete bbtaudata;
}

void Simul::run(int nEvnt){
  for (int iEvent = 0; iEvent < nEvnt; ++iEvent) {
    if (!pythia.next()){ break;}
    if (iEvent < 2) {pythia.info.list(); pythia.event.list(); pythia.process.list();} 
    process = new MyProcess(&pythia);
    
    event = new MyEvent(&pythia,0.4);
    bbmuanalyse(); 
    bbtauanalyse(); 
    delete event;
    delete process;
  }
  ofstream OutFile;
  OutFile.open(file.c_str());
  OutFile << bbmudata->GetN() << endl;
  OutFile << bbtaudata->GetN() << endl;
  OutFile.close();
  
  cout << bbmudata->GetN() << endl;
  cout << bbtaudata->GetN() << endl;
  bbmudata->Save();
  bbtaudata->Save();
  pythia.statistics();
}



void Simul::bbmuanalyse(){
  //  bbmudata->AddEvent();
  bbmudata->Fill(process,event);
}

void Simul::bbtauanalyse(){
  //  Zbbdata->AddEvent();
  bbtaudata->Fill(process,event);
}

