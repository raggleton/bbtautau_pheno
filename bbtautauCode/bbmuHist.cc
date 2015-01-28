using namespace std;
using namespace Pythia8;

class bbmuHist {
private:
  Hist bSp;
  Hist BmSp;
  Hist bJSp;
  Hist FjSp;
  Hist Fjm;
  Hist Nb;
  Hist NBm;
  Hist NbJ;
  Hist NFj;
  Hist bbtauHad;
  Hist bbHad;
  Hist tataHad;
  Hist a1a1Had;
  Hist h1a1a1Had;
  Hist a1h1Had;
  Hist bbtauSub;
  Hist a1a1Sub;
  Hist h1a1a1Sub;
  Hist a1h1Sub;

  Hist Fat2bpT;
  Hist Fat2bm;
  Hist FatMDpT;
  Hist FatMDm;
  Hist Fat1tpT;
  Hist Fat1tm;
  Hist FatinpT;
  Hist Fatinm;
  string name;
  double ScaleFactor;
  long N;
  void SaveHist(Hist *h, string file, double delta);

public:
  bbmuHist(string model, double fact);
  ~bbmuHist();
  void Fill(MyProcess *Mp, MyEvent *Me);
  void AddEvent();
  long GetN();

  void Save();

};

bbmuHist::bbmuHist(string model, double fact) {

  bSp.book("", 100, 0., 100.);
  bSp.null();
  BmSp.book("", 100, 0., 100.);
  BmSp.null();
  bJSp.book("", 100, 0., 100.);
  bJSp.null();
  FjSp.book("", 100, 0., 100.);
  FjSp.null();

  Fjm.book("", 100, 0., 100.);
  Fjm.null();

  Nb.book("", 15, 0.5, 15.5);
  Nb.null();
  NBm.book("", 15, 0.5, 15.5);
  NBm.null();
  NbJ.book("", 15, 0.5, 15.5);
  NbJ.null();
  NFj.book("", 15, 0.5, 15.5);
  NFj.null();

  bbtauHad.book("", 100, 0., 150.);
  bbtauHad.null();
  bbHad.book("", 100, 0., 100.);
  bbHad.null();
  tataHad.book("", 100, 0., 150.);
  tataHad.null();
  a1a1Had.book("", 100, 0., 150.);
  a1a1Had.null();
  h1a1a1Had.book("", 100, 0., 100.);
  h1a1a1Had.null();
  a1h1Had.book("", 100, 0., 150.);
  a1h1Had.null();

  bbtauSub.book("", 100, 0., 150.);
  bbtauSub.null();
  a1a1Sub.book("", 100, 0., 150.);
  a1a1Sub.null();
  h1a1a1Sub.book("", 100, 0., 100.);
  h1a1a1Sub.null();
  a1h1Sub.book("", 100, 0., 150.);
  a1h1Sub.null();

  Fat2bpT.book("", 100, 0., 200.);
  Fat2bpT.null();
  Fat2bm.book("", 100, 0., 200.);
  Fat2bm.null();
  FatMDpT.book("", 100, 0., 200.);
  FatMDpT.null();
  FatMDm.book("", 100, 0., 200.);
  FatMDm.null();
  Fat1tpT.book("", 100, 0., 200.);
  Fat1tpT.null();
  Fat1tm.book("", 100, 0., 200.);
  Fat1tm.null();
  FatinpT.book("", 100, 0., 200.);
  FatinpT.null();
  Fatinm.book("", 100, 0., 200.);
  Fatinm.null();

  name = model;
  N = 0;
  ScaleFactor = fact;
}

bbmuHist::~bbmuHist() {
}

void bbmuHist::Fill(MyProcess *Mp, MyEvent *Me) {
  Mp->Spec(&bSp, 5, 1);
  Me->Spec(&bJSp, 5, 1);
  Me->Spec(&BmSp, 500, 1);
  Me->Spec(&FjSp, 55, 1);

  Me->InvM(&Fjm, 55, 1);

  Nb.fill(Mp->MaxI(5) + Mp->MaxI(-5), 1);
  NbJ.fill(Me->MaxI(5) + Me->MaxI(-5), 1);
  NBm.fill(Me->MaxI(500), 1);
  NFj.fill(Me->MaxI(55), 1);

  N++;

  Me->FatEff(&Fat2bpT, &Fat2bm, &FatMDpT, &FatMDm, &Fat1tpT, &Fat1tm, &FatinpT,
      &Fatinm);

  Me->bbmuSub(&bbtauSub, &a1a1Sub, &h1a1a1Sub, &a1h1Sub, 10);
  if (Me->bbmu(&bbtauHad, &a1a1Had, &h1a1a1Had, &a1h1Had, 10)) {
    Me->InvM(&bbHad, 5, -5, 1);
    Me->InvM(&bbHad, 5, 5, 0.5);
    Me->InvM(&bbHad, -5, -5, 0.5);
    Me->InvM(&tataHad, 13, -13, 1);
    Me->InvM(&tataHad, 13, 13, 0.5);
    Me->InvM(&tataHad, -13, -13, 0.5);
  }
}

void bbmuHist::AddEvent() {
  N++;
}

long bbmuHist::GetN() {
  return N;
}

void bbmuHist::Save() {

  SaveHist(&bSp, "bbmudata/bSp_" + name, 1.);
  SaveHist(&BmSp, "bbmudata/BmSp_" + name, 1.);
  SaveHist(&bJSp, "bbmudata/bJSp_" + name, 1.);
  SaveHist(&FjSp, "bbmudata/FjSp_" + name, 1.);

  SaveHist(&Fjm, "bbmudata/Fjm_" + name, 1.);

  SaveHist(&Nb, "bbmudata/Nb_" + name, 1.);
  SaveHist(&NBm, "bbmudata/NBm_" + name, 1.);
  SaveHist(&NbJ, "bbmudata/NbJ_" + name, 1.);
  SaveHist(&NFj, "bbmudata/NFj_" + name, 1.);

  SaveHist(&bbtauHad, "bbmudata/bbmuH_" + name, 1.5);
  SaveHist(&bbHad, "bbmudata/bbH_" + name, 1.);
  SaveHist(&tataHad, "bbmudata/mumuH_" + name, 1.5);
  SaveHist(&a1a1Had, "bbmudata/a1a1H_" + name, 1.5);
  SaveHist(&h1a1a1Had, "bbmudata/h1a1a1H_" + name, 1.);
  SaveHist(&a1h1Had, "bbmudata/a1h1H_" + name, 1.5);

  SaveHist(&bbtauSub, "bbmudata/bbmuS_" + name, 1.5);
  SaveHist(&a1a1Sub, "bbmudata/a1a1S_" + name, 1.5);
  SaveHist(&h1a1a1Sub, "bbmudata/h1a1a1S_" + name, 1.);
  SaveHist(&a1h1Sub, "bbmudata/a1h1S_" + name, 1.5);

  SaveHist(&Fat2bpT, "bbmudata/Fat2bpT_" + name, 2.);
  SaveHist(&Fat2bm, "bbmudata/Fat2bm_" + name, 2.);
  SaveHist(&FatMDpT, "bbmudata/FatMDpT_" + name, 2.);
  SaveHist(&FatMDm, "bbmudata/FatMDm_" + name, 2.);
  SaveHist(&Fat1tpT, "bbmudata/Fat1tpT_" + name, 2.);
  SaveHist(&Fat1tm, "bbmudata/Fat1tm_" + name, 2.);
  SaveHist(&FatinpT, "bbmudata/FatinpT_" + name, 2.);
  SaveHist(&Fatinm, "bbmudata/Fatinm_" + name, 2.);

}

void bbmuHist::SaveHist(Hist *h, string file, double delta) {
  *h *= ScaleFactor / (delta * (double) N);
  ofstream OutFile;
  OutFile.open(file.c_str());
  h->table(OutFile);
  OutFile.close();
}
