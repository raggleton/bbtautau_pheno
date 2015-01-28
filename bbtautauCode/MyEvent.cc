// #include <iostream>
// #include <fstream>
// #include <string>
// #include "Pythia.h"
// #include "Simul.cc"
using namespace std;
using namespace Pythia8;
//const double pi = 3.14159265358979323846;

class MyEvent {
private:
  vector<int> em;
  vector<int> ep;
  vector<int> mum;
  vector<int> mup;
  vector<int> gamma;
  vector<int> tauJet;
  vector<int> taupJet;
  vector<int> Fatbb;
  vector<int> Bm;
  vector<int> bJet;
  vector<int> bpJet;
  vector<int> Jet20;
  vector<int> lJJ;
  Sphericity *sph;
  Sphericity *linsph;
  Thrust *thr;
  MyJet *Jet;
  Pythia *pythia;
  Vec4 tau1p;
  Vec4 tau2p;
  double CollinearM(int tau1, int tau2);
  int I(int sp, int i);
  bool isJet(int sp);
  bool AcceptLepton(int lep);
  double deltaR(int i1, int i2);
  double deltaRlJ(int i1, int i2);
  double deltaRJJ(int i1, int i2);
  double sign(double d);
  double phi2(double phi);
  bool Same(int temp1, int temp2, int sp1, int sp2);
  double m(const Vec4& v1);
  double m(const Vec4& v1, const Vec4& v2);
  double m(const Vec4& v1, const Vec4& v2, const Vec4& v3);
  double m(const Vec4& v1, const Vec4& v2, const Vec4& v3, const Vec4& v4);
  int Bmes(int i);
public:
  MyEvent(Pythia *p, double Rparam);
  ~MyEvent();
  bool Fourtau(Hist *ta4, Hist *a1a1, Hist *h1a1a1, Hist *a1h1, double tares);
  bool bbtau(Hist *bbta, Hist *a1a1, Hist *h1a1a1, Hist *a1h1, double tares);
  bool bbtauSub(Hist *bbta, Hist *a1a1, Hist *h1a1a1, Hist *a1h1, double tares);
  bool bbtauhlSub(Hist *bbta, Hist *a1a1, Hist *h1a1a1, Hist *a1h1,
      double tares);
  bool bbmu(Hist *bbta, Hist *a1a1, Hist *h1a1a1, Hist *a1h1, double tares);
  bool bbmuSub(Hist *bbta, Hist *a1a1, Hist *h1a1a1, Hist *a1h1, double tares);
  bool bbgamma(Hist *bbta, Hist *a1a1, Hist *h1a1a1, Hist *a1h1, double tares);
  bool bbgammaSub(Hist *bbta, Hist *a1a1, Hist *h1a1a1, Hist *a1h1,
      double tares);
  bool Fourb(Hist *b4, Hist *a1a1, Hist *h1a1a1, Hist *a1h1, double bres);
  bool FourbJetSub(Hist *b4, Hist *a1a1, Hist *h1a1a1, Hist *a1h1, double bres);
  bool FourbJSemi(Hist *b4, Hist *a1a1, Hist *h1a1a1, Hist *a1h1, double bres);
  bool Ztata(Hist *tatah1, Hist *Za1, double tares);
  bool Zbb(Hist *bbh1, Hist *Za1, double bres, Hist *dRbb, Hist *dRlb,
      Hist *dRll, Hist *ZH);
  bool ZbbJetSub(Hist *bbh1, Hist *Za1, double bres, Hist *dRbb, Hist *dRlb,
      Hist *dRll, Hist *ZH);
  void FatEff(Hist *F2bpT, Hist *F2bm, Hist *FMDpT, Hist *FMDm, Hist *F1tpT,
      Hist *F1tm, Hist *FinpT, Hist *Finm);

  void Spec(Hist *h, int sp, double w);
  void eta(Hist *h, int sp, double w);
  void InvMlJJ(Hist *h);
  void InvM(Hist *h, int sp1, double w);
  void InvM(Hist *h, int sp1, int sp2, double w);
  void InvM(Hist *h, int sp1, int sp2, int sp3, double w);
  void InvM(Hist *h, int sp1, int sp2, int sp3, int sp4, double w);
  void InvMpTMiss(Hist *h, double w);
  double mTC();
  int NlepJet(double cutoff);
  int NJet(double cutoff);
  double MJet();
  double Spher();
  double LinSpher();
  double Thru();
  double TranSph();
  int MaxI(int sp);
  int pT20lep;
  int Ntau;
  int pT15lep;
  int SSisolep;
  int NlJJ;
  double pTmaxlep;
  double pTMiss;
  double pxMiss;
  double pyMiss;
};

MyEvent::MyEvent(Pythia *p, double Rparam) {
  pythia = p;
  for (int i = 0; i < pythia->event.size(); i++) {
    if (abs(pythia->event[i].eta()) > 2.5) {
      pythia->event[i].statusNeg();
    }
  }
  //  for (int i = 0; i < pythia->event.size(); i++){ if (abs(Bmes(i)) == 5) {pythia->event[i].statusNeg();}}
  //  sph->analyze(pythia->event);
  //linsph->analyze(pythia->event);
  // thr->analyze(pythia->event);
  Jet = new MyJet(pythia, Rparam);
  tauJet = Jet->GettauJet();
  taupJet = Jet->GettaupJet();
  bJet = Jet->GetbJet();
  bpJet = Jet->GetbpJet();
  Jet20 = Jet->GetJet20();
  Fatbb = Jet->GetFatbb();
  Ntau = 0;
  for (unsigned int i = 0; i < tauJet.size(); i++) {
    if (Jet->pT(tauJet[i]) > 20) {
      Ntau++;
    }
  }
  for (unsigned int i = 0; i < taupJet.size(); i++) {
    if (Jet->pT(taupJet[i]) > 20) {
      Ntau++;
    }
  }
  lJJ.resize(0);
  double px = 0;
  double py = 0;
  pT15lep = 0;
  for (int i = 0; i < pythia->event.size(); i++) {
    if (abs(Bmes(i)) == 5 && abs(abs(pythia->event[i].status()) - 83) <= 3) {
      Bm.push_back(i);
    }
    switch (pythia->event[i].id()) {
    case 11:
      if (AcceptLepton(i)) {
        em.push_back(i);
      }
      if (pythia->event[i].pT() > 15. && pythia->event[i].status() > 0) {
        pT15lep++;
      }
      break;
    case -11:
      if (AcceptLepton(i)) {
        ep.push_back(i);
      }
      if (pythia->event[i].pT() > 15. && pythia->event[i].status() > 0) {
        pT15lep++;
      }
      break;
    case 13:
      if (AcceptLepton(i)) {
        mum.push_back(i);
      }
      if (pythia->event[i].pT() > 15. && pythia->event[i].status() > 0) {
        pT15lep++;
      }
      break;
    case -13:
      if (AcceptLepton(i)) {
        mup.push_back(i);
      }
      if (pythia->event[i].pT() > 15. && pythia->event[i].status() > 0) {
        pT15lep++;
      }
      break;
    case 22:
      if (AcceptLepton(i)) {
        gamma.push_back(i);
      }
      break;
    case 12:
    case -12:
    case 14:
    case -14:
    case 16:
    case -16:
    case 1000022:
    case -1000022:
      if (pythia->event[i].status() > 0) {
        px += pythia->event[i].px();
        py += pythia->event[i].py();
      }
      break;
    }
  }
  pTMiss = sqrt(px * px + py * py);
  pxMiss = px;
  pyMiss = py;
  pT20lep = 0;
  pTmaxlep = 0;
  int poslep = 0;
  int neglep = 0;
  for (unsigned int i = 0; i < em.size(); i++) {
    if (pythia->event[em[i]].pT() > 20) {
      int Jco = 0;
      int temp1 = 0;
      int temp2 = 0;
      for (unsigned int j = 0; j < Jet20.size(); j++) {
        if (deltaRlJ(em[i], Jet20[j]) < 1) {
          Jco++;
          (Jco == 1) ? temp1 = j : temp2 = j;
        }
      }
      if (Jco == 2) {
        lJJ.push_back(em[i]);
        lJJ.push_back(Jet20[temp1]);
        lJJ.push_back(Jet20[temp2]);
      }
      pT20lep++;
      neglep++;
    }
    if (pythia->event[em[i]].pT() > pTmaxlep) {
      pTmaxlep = pythia->event[em[i]].pT();
    }
  }
  for (unsigned int i = 0; i < ep.size(); i++) {
    if (pythia->event[ep[i]].pT() > 20) {
      int Jco = 0;
      int temp1 = 0;
      int temp2 = 0;
      for (unsigned int j = 0; j < Jet20.size(); j++) {
        if (deltaRlJ(ep[i], Jet20[j]) < 1) {
          Jco++;
          (Jco == 1) ? temp1 = j : temp2 = j;
        }
      }
      if (Jco == 2) {
        lJJ.push_back(ep[i]);
        lJJ.push_back(Jet20[temp1]);
        lJJ.push_back(Jet20[temp2]);
      }
      pT20lep++;
      poslep++;
    }
    if (pythia->event[ep[i]].pT() > pTmaxlep) {
      pTmaxlep = pythia->event[ep[i]].pT();
    }
  }
  for (unsigned int i = 0; i < mum.size(); i++) {
    if (pythia->event[mum[i]].pT() > 20) {
      int Jco = 0;
      int temp1 = 0;
      int temp2 = 0;
      for (unsigned int j = 0; j < Jet20.size(); j++) {
        if (deltaRlJ(mum[i], Jet20[j]) < 1) {
          Jco++;
          (Jco == 1) ? temp1 = j : temp2 = j;
        }
      }
      if (Jco == 2) {
        lJJ.push_back(mum[i]);
        lJJ.push_back(Jet20[temp1]);
        lJJ.push_back(Jet20[temp2]);
      }
      pT20lep++;
      neglep++;
    }
    if (pythia->event[mum[i]].pT() > pTmaxlep) {
      pTmaxlep = pythia->event[mum[i]].pT();
    }
  }
  for (unsigned int i = 0; i < mup.size(); i++) {
    if (pythia->event[mup[i]].pT() > 20) {
      int Jco = 0;
      int temp1 = 0;
      int temp2 = 0;
      for (unsigned int j = 0; j < Jet20.size(); j++) {
        if (deltaRlJ(mup[i], Jet20[j]) < 1) {
          Jco++;
          (Jco == 1) ? temp1 = j : temp2 = j;
        }
      }
      if (Jco == 2) {
        lJJ.push_back(mup[i]);
        lJJ.push_back(Jet20[temp1]);
        lJJ.push_back(Jet20[temp2]);
      }
      pT20lep++;
      poslep++;
    }
    if (pythia->event[mup[i]].pT() > pTmaxlep) {
      pTmaxlep = pythia->event[mup[i]].pT();
    }
  }
  if (poslep > neglep) {
    SSisolep = poslep;
  } else {
    SSisolep = neglep;
  }
  NlJJ = lJJ.size() / 3;
}

MyEvent::~MyEvent() {
  delete Jet;
}

bool MyEvent::FourbJSemi(Hist *b4, Hist *a1a1, Hist *h1a1a1, Hist *a1h1,
    double bres) {
  int bArray[2];
  if (bJet.size() + bpJet.size() != 2 || Fatbb.size() != 1) {
    return false;
  } else {
    for (int i = 0; i < bJet.size(); i++) {
      bArray[i] = bJet[i];
    }
    for (int i = 0; i < bpJet.size(); i++) {
      bArray[bJet.size() + i] = bpJet[i];
    }
  }
  double Inv4b = m(Jet->p(bArray[0]), Jet->p(bArray[1]), Jet->p(Fatbb[0]));
  b4->fill(Inv4b);
  double Inv1, Inv2;
  Inv1 = m(Jet->p(Fatbb[0]));
  Inv2 = m(Jet->p(bArray[0]), Jet->p(bArray[1]));
  if (abs(Inv1 - 125) < bres) {
    a1h1->fill(Inv2);
  }
  if (abs(Inv2 - 125) < bres) {
    a1h1->fill(Inv1);
  }
  if (abs(Inv1 - Inv2) < bres) {
    a1a1->fill(Inv1);
    a1a1->fill(Inv2);
    if (abs(Inv4b - 125) < bres) {
      h1a1a1->fill(Inv1);
      h1a1a1->fill(Inv2);
    }
  }
  return true;
}
bool MyEvent::FourbJetSub(Hist *b4, Hist *a1a1, Hist *h1a1a1, Hist *a1h1,
    double bres) {
  int bArray[2];
  if (Fatbb.size() != 2) {
    return false;
  } else {
    for (int i = 0; i < Fatbb.size(); i++) {
      bArray[i] = Fatbb[i];
    }
  }
  double Inv4b = m(Jet->p(bArray[0]), Jet->p(bArray[1]));
  b4->fill(Inv4b);
  double Inv1, Inv2;
  Inv1 = m(Jet->p(bArray[0]));
  Inv2 = m(Jet->p(bArray[1]));
  if (abs(Inv1 - 125) < bres) {
    a1h1->fill(Inv2);
  }
  if (abs(Inv2 - 125) < bres) {
    a1h1->fill(Inv1);
  }
  if (abs(Inv1 - Inv2) < bres) {
    a1a1->fill(Inv1);
    a1a1->fill(Inv2);
    if (abs(Inv4b - 125) < 2 * bres) {
      h1a1a1->fill(Inv1);
      h1a1a1->fill(Inv2);
    }
  }
  return true;
}

bool MyEvent::Fourb(Hist *b4, Hist *a1a1, Hist *h1a1a1, Hist *a1h1,
    double bres) {
  int bArray[4];
  if (bJet.size() + bpJet.size() != 4) {
    return false;
  } else {
    for (int i = 0; i < bJet.size(); i++) {
      bArray[i] = bJet[i];
    }
    for (int i = 0; i < bpJet.size(); i++) {
      bArray[bJet.size() + i] = bpJet[i];
    }
  }
  double Inv4b = m(Jet->p(bArray[0]), Jet->p(bArray[1]), Jet->p(bArray[2]),
      Jet->p(bArray[3]));
  b4->fill(Inv4b);
  int tempA[2][3] = { { 2, 1, 1 }, { 3, 3, 2 } };
  double Inv1, Inv2;
  for (int i = 1; i < 4; i++) {
    Inv1 = m(Jet->p(bArray[0]), Jet->p(bArray[i]));
    Inv2 = m(Jet->p(bArray[tempA[0][i - 1]]), Jet->p(bArray[tempA[1][i - 1]]));
    if (abs(Inv1 - 125) < bres) {
      a1h1->fill(Inv2);
    }
    if (abs(Inv2 - 125) < bres) {
      a1h1->fill(Inv1);
    }
    if (abs(Inv1 - Inv2) < bres) {
      a1a1->fill(Inv1);
      a1a1->fill(Inv2);
      if (abs(Inv4b - 125) < 2 * bres) {
        h1a1a1->fill(Inv1);
        h1a1a1->fill(Inv2);
      }
    }
  }
  return true;
}

bool MyEvent::bbtau(Hist *bbta, Hist *a1a1, Hist *h1a1a1, Hist *a1h1,
    double tares) {
  int taArray[2];
  if (tauJet.size() + taupJet.size() != 2) {
    return false;
  } else {
    for (int i = 0; i < tauJet.size(); i++) {
      taArray[i] = tauJet[i];
    }
    for (int i = 0; i < taupJet.size(); i++) {
      taArray[tauJet.size() + i] = taupJet[i];
    }
  }
  int bArray[2];
  if (bJet.size() + bpJet.size() != 2) {
    return false;
  } else {
    for (int i = 0; i < bJet.size(); i++) {
      bArray[i] = bJet[i];
    }
    for (int i = 0; i < bpJet.size(); i++) {
      bArray[bJet.size() + i] = bpJet[i];
    }
  }
  double Invtau = CollinearM(taArray[0], taArray[1]);
  double Invbb = m(Jet->p(bArray[0]), Jet->p(bArray[1]));
  double Inv4tau = m(Jet->p(bArray[0]), Jet->p(bArray[1]), tau1p, tau2p);
  bbta->fill(Inv4tau);
  if (abs(Invtau - 125) < tares) {
    a1h1->fill(Invbb);
  }
  if (abs(Invbb - 125) < tares) {
    a1h1->fill(Invtau);
  }
  if (abs(Invtau - Invbb) < tares) {
    a1a1->fill(Invtau);
    a1a1->fill(Invbb);
    if (abs(Inv4tau - 125) < 2 * tares) {
      h1a1a1->fill(Invtau);
      h1a1a1->fill(Invbb);
    }
  }
  return true;
}
bool MyEvent::bbtauSub(Hist *bbta, Hist *a1a1, Hist *h1a1a1, Hist *a1h1,
    double tares) {
  int taArray[2];
  if (tauJet.size() + taupJet.size() != 2 || Fatbb.size() != 1) {
    return false;
  } else {
    for (int i = 0; i < tauJet.size(); i++) {
      taArray[i] = tauJet[i];
    }
    for (int i = 0; i < taupJet.size(); i++) {
      taArray[tauJet.size() + i] = taupJet[i];
    }
  }
  double Invtau = CollinearM(taArray[0], taArray[1]);
  double Invbb = m(Jet->p(Fatbb[0]));
  double Inv4tau = m(Jet->p(Fatbb[0]), tau1p, tau2p);
  bbta->fill(Inv4tau);
  if (abs(Invtau - 125) < tares) {
    a1h1->fill(Invbb);
  }
  if (abs(Invbb - 125) < tares) {
    a1h1->fill(Invtau);
  }
  if (abs(Invtau - Invbb) < tares) {
    a1a1->fill(Invtau);
    a1a1->fill(Invbb);
    if (abs(Inv4tau - 125) < 2 * tares) {
      h1a1a1->fill(Invtau);
      h1a1a1->fill(Invbb);
    }
  }
  return true;
}

bool MyEvent::bbtauhlSub(Hist *bbta, Hist *a1a1, Hist *h1a1a1, Hist *a1h1,
    double tares) {
  int taArray[2];
  if (tauJet.size() + taupJet.size() != 1
      || em.size() + ep.size() + mum.size() + mup.size() != 1
      || Fatbb.size() != 1) {
    return false;
  } else {
    if (tauJet.size() == 1) {
      taArray[0] = tauJet[0];
    } else if (taupJet.size() == 1) {
      taArray[0] = taupJet[0];
    }
    if (em.size() == 1) {
      taArray[1] = em[0];
    } else if (ep.size() == 1) {
      taArray[1] = ep[0];
    } else if (mum.size() == 1) {
      taArray[1] = mum[0];
    } else if (mup.size() == 1) {
      taArray[1] = mup[0];
    }
  }
  double Invtau = m(Jet->p(taArray[0]), pythia->event[taArray[1]].p());
  double Invbb = m(Jet->p(Fatbb[0]));
  double Inv4tau = m(Jet->p(Fatbb[0]), Jet->p(taArray[0]),
      pythia->event[taArray[1]].p());
  bbta->fill(Inv4tau);
  if (abs(Invtau - 125) < tares) {
    a1h1->fill(Invbb);
  }
  if (abs(Invbb - 125) < tares) {
    a1h1->fill(Invtau);
  }
  if (abs(Invtau - Invbb) < tares) {
    a1a1->fill(Invtau);
    a1a1->fill(Invbb);
    if (abs(Inv4tau - 125) < 2 * tares) {
      h1a1a1->fill(Invtau);
      h1a1a1->fill(Invbb);
    }
  }
  return true;
}

bool MyEvent::bbmu(Hist *bbta, Hist *a1a1, Hist *h1a1a1, Hist *a1h1,
    double tares) {
  int taArray[2];
  if (mum.size() + mup.size() != 2) {
    return false;
  } else {
    for (int i = 0; i < mum.size(); i++) {
      taArray[i] = mum[i];
    }
    for (int i = 0; i < mup.size(); i++) {
      taArray[mum.size() + i] = mup[i];
    }
  }
  int bArray[2];
  if (bJet.size() + bpJet.size() != 2) {
    return false;
  } else {
    for (int i = 0; i < bJet.size(); i++) {
      bArray[i] = bJet[i];
    }
    for (int i = 0; i < bpJet.size(); i++) {
      bArray[bJet.size() + i] = bpJet[i];
    }
  }
  double Invtau = m(pythia->event[taArray[0]].p(),
      pythia->event[taArray[1]].p());
  double Invbb = m(Jet->p(bArray[0]), Jet->p(bArray[1]));
  double Inv4tau = m(Jet->p(bArray[0]), Jet->p(bArray[1]),
      pythia->event[taArray[0]].p(), pythia->event[taArray[1]].p());
  bbta->fill(Inv4tau);
  if (abs(Invtau - 125) < tares) {
    a1h1->fill(Invbb);
  }
  if (abs(Invbb - 125) < tares) {
    a1h1->fill(Invtau);
  }
  if (abs(Invtau - Invbb) < tares) {
    a1a1->fill(Invtau);
    a1a1->fill(Invbb);
    if (abs(Inv4tau - 125) < 2 * tares) {
      h1a1a1->fill(Invtau);
      h1a1a1->fill(Invbb);
    }
  }
  return true;
}
bool MyEvent::bbmuSub(Hist *bbta, Hist *a1a1, Hist *h1a1a1, Hist *a1h1,
    double tares) {
  int taArray[2];
  if (mum.size() + mup.size() != 2 || Fatbb.size() != 1) {
    return false;
  } else {
    for (int i = 0; i < mum.size(); i++) {
      taArray[i] = mum[i];
    }
    for (int i = 0; i < mup.size(); i++) {
      taArray[mum.size() + i] = mup[i];
    }
  }
  double Invtau = m(pythia->event[taArray[0]].p(),
      pythia->event[taArray[1]].p());
  double Invbb = m(Jet->p(Fatbb[0]));
  double Inv4tau = m(Jet->p(Fatbb[0]), pythia->event[taArray[0]].p(),
      pythia->event[taArray[1]].p());
  bbta->fill(Inv4tau);
  if (abs(Invtau - 125) < tares) {
    a1h1->fill(Invbb);
  }
  if (abs(Invbb - 125) < tares) {
    a1h1->fill(Invtau);
  }
  if (abs(Invtau - Invbb) < tares) {
    a1a1->fill(Invtau);
    a1a1->fill(Invbb);
    if (abs(Inv4tau - 125) < 2 * tares) {
      h1a1a1->fill(Invtau);
      h1a1a1->fill(Invbb);
    }
  }
  return true;
}

bool MyEvent::bbgamma(Hist *bbta, Hist *a1a1, Hist *h1a1a1, Hist *a1h1,
    double tares) {
  if (gamma.size() != 2) {
    return false;
  }
  int bArray[2];
  if (bJet.size() + bpJet.size() != 2) {
    return false;
  } else {
    for (int i = 0; i < bJet.size(); i++) {
      bArray[i] = bJet[i];
    }
    for (int i = 0; i < bpJet.size(); i++) {
      bArray[bJet.size() + i] = bpJet[i];
    }
  }
  double Invtau = m(pythia->event[gamma[0]].p(), pythia->event[gamma[1]].p());
  double Invbb = m(Jet->p(bArray[0]), Jet->p(bArray[1]));
  double Inv4tau = m(Jet->p(bArray[0]), Jet->p(bArray[1]),
      pythia->event[gamma[0]].p(), pythia->event[gamma[1]].p());
  bbta->fill(Inv4tau);
  if (abs(Invtau - 125) < 0.25 * tares) {
    a1h1->fill(Invbb);
  }
  if (abs(Invbb - 125) < tares) {
    a1h1->fill(Invtau);
  }
  if (abs(Invtau - Invbb) < tares) {
    a1a1->fill(Invtau);
    a1a1->fill(Invbb);
    if (abs(Inv4tau - 125) < tares) {
      h1a1a1->fill(Invtau);
      h1a1a1->fill(Invbb);
    }
  }
  return true;
}
bool MyEvent::bbgammaSub(Hist *bbta, Hist *a1a1, Hist *h1a1a1, Hist *a1h1,
    double tares) {
  if (gamma.size() != 2 || Fatbb.size() != 1) {
    return false;
  }
  double Invtau = m(pythia->event[gamma[0]].p(), pythia->event[gamma[1]].p());
  double Invbb = m(Jet->p(Fatbb[0]));
  double Inv4tau = m(Jet->p(Fatbb[0]), pythia->event[gamma[0]].p(),
      pythia->event[gamma[1]].p());
  bbta->fill(Inv4tau);
  if (abs(Invtau - 125) < 0.25 * tares) {
    a1h1->fill(Invbb);
  }
  if (abs(Invbb - 125) < tares) {
    a1h1->fill(Invtau);
  }
  if (abs(Invtau - Invbb) < tares) {
    a1a1->fill(Invtau);
    a1a1->fill(Invbb);
    if (abs(Inv4tau - 125) < tares) {
      h1a1a1->fill(Invtau);
      h1a1a1->fill(Invbb);
    }
  }
  return true;
}

bool MyEvent::Fourtau(Hist *ta4, Hist *a1a1, Hist *h1a1a1, Hist *a1h1,
    double tares) {
  int taArray[4];
  if (tauJet.size() + taupJet.size() != 4) {
    return false;
  } else {
    for (int i = 0; i < tauJet.size(); i++) {
      taArray[i] = tauJet[i];
    }
    for (int i = 0; i < taupJet.size(); i++) {
      taArray[tauJet.size() + i] = taupJet[i];
    }
  }
  double Inv4tau = m(Jet->p(taArray[0]), Jet->p(taArray[1]), Jet->p(taArray[2]),
      Jet->p(taArray[3]));
  ta4->fill(Inv4tau);
  int tempA[2][3] = { { 2, 1, 1 }, { 3, 3, 2 } };
  double Inv1, Inv2;
  for (int i = 1; i < 4; i++) {
    Inv1 = m(Jet->p(taArray[0]), Jet->p(taArray[i]));
    Inv2 = m(Jet->p(taArray[tempA[0][i - 1]]),
        Jet->p(taArray[tempA[1][i - 1]]));
    if (abs(Inv1 - 125) < tares) {
      a1h1->fill(Inv2);
    }
    if (abs(Inv2 - 125) < tares) {
      a1h1->fill(Inv1);
    }
    if (abs(Inv1 - Inv2) < tares) {
      a1a1->fill(Inv1);
      a1a1->fill(Inv2);
      if (abs(Inv4tau - 125) < 2 * tares) {
        h1a1a1->fill(Inv1);
        h1a1a1->fill(Inv2);
      }
    }
  }
  return true;
}

bool MyEvent::Ztata(Hist *tatah1, Hist *Za1, double tares) {
  int taArray[2];
  if (tauJet.size() + taupJet.size() != 2) {
    return false;
  } else {
    for (int i = 0; i < tauJet.size(); i++) {
      taArray[i] = tauJet[i];
    }
    for (int i = 0; i < taupJet.size(); i++) {
      taArray[tauJet.size() + i] = taupJet[i];
    }
  }
  int Zl = -1;
  int Zlp = -1;
  if (em.size() == 1 && ep.size() == 1) {
    Zl = em[0];
    Zlp = ep[0];
  } else if (mum.size() == 1 && mup.size() == 1) {
    Zl = mum[0];
    Zlp = mup[0];
  } else {
    return false;
  }
  if (abs(m(pythia->event[Zl].p(), pythia->event[Zlp].p()) - 90) > 10) {
    return false;
  }
  double Inv2ta = CollinearM(taArray[0], taArray[1]);
  double InvZa1 = m(pythia->event[Zl].p(), pythia->event[Zlp].p(), tau1p,
      tau2p);
  Za1->fill(InvZa1);
  if (abs(InvZa1 - 125) < tares) {
    tatah1->fill(Inv2ta);
  }
  return true;
}

bool MyEvent::Zbb(Hist *bbh1, Hist *Za1, double bres, Hist *dRbb, Hist *dRlb,
    Hist *dRll, Hist *ZH) {
  int bArray[2];
  if (bJet.size() + bpJet.size() != 2) {
    return false;
  } else {
    for (int i = 0; i < bJet.size(); i++) {
      bArray[i] = bJet[i];
    }
    for (int i = 0; i < bpJet.size(); i++) {
      bArray[bJet.size() + i] = bpJet[i];
    }
  }
  int Zl = -1;
  int Zlp = -1;
  if (em.size() == 1 && ep.size() == 1) {
    Zl = em[0];
    Zlp = ep[0];
  } else if (mum.size() == 1 && mup.size() == 1) {
    Zl = mum[0];
    Zlp = mup[0];
  } else {
    return false;
  }
  if (abs(m(pythia->event[Zl].p(), pythia->event[Zlp].p()) - 90) > 10) {
    return false;
  }
  dRbb->fill(deltaRJJ(bArray[0], bArray[1]));
  dRlb->fill(deltaRlJ(Zl, bArray[0]));
  dRlb->fill(deltaRlJ(Zl, bArray[1]));
  dRlb->fill(deltaRlJ(Zlp, bArray[0]));
  dRlb->fill(deltaRlJ(Zlp, bArray[1]));
  dRll->fill(deltaR(Zl, Zlp));
  double Inv2b = m(Jet->p(bArray[0]), Jet->p(bArray[1]));
  double InvZa1 = m(pythia->event[Zl].p(), pythia->event[Zlp].p(),
      Jet->p(bArray[0]), Jet->p(bArray[1]));
  Za1->fill(InvZa1);
  if (abs(Inv2b - 125) < bres) {
    ZH->fill(InvZa1);
  }
  if (abs(InvZa1 - 125) < bres) {
    bbh1->fill(Inv2b);
  }
  return true;
}

bool MyEvent::ZbbJetSub(Hist *bbh1, Hist *Za1, double bres, Hist *dRbb,
    Hist *dRlb, Hist *dRll, Hist *ZH) {
  if (Fatbb.size() != 1) {
    return false;
  }
  int Zl = -1;
  int Zlp = -1;
  if (em.size() == 1 && ep.size() == 1) {
    Zl = em[0];
    Zlp = ep[0];
  } else if (mum.size() == 1 && mup.size() == 1) {
    Zl = mum[0];
    Zlp = mup[0];
  } else {
    return false;
  }
  if (abs(m(pythia->event[Zl].p(), pythia->event[Zlp].p()) - 90) > 10) {
    return false;
  }
  double Inv2b = m(Jet->p(Fatbb[0]));
  double InvZa1 = m(pythia->event[Zl].p(), pythia->event[Zlp].p(),
      Jet->p(Fatbb[0]));
  dRbb->fill(Jet->DeltaRfj[0]);
  dRlb->fill(deltaRlJ(Zl, Fatbb[0]));
  dRlb->fill(deltaRlJ(Zlp, Fatbb[0]));
  dRll->fill(deltaR(Zl, Zlp));
  Za1->fill(InvZa1);
  if (abs(Inv2b - 125) < bres) {
    ZH->fill(InvZa1);
  }
  if (abs(InvZa1 - 125) < bres) {
    bbh1->fill(Inv2b);
  }
  return true;
}

double MyEvent::CollinearM(int tau1, int tau2) {
  double q2 = Jet->px(tau2) / Jet->px(tau1) - Jet->py(tau2) / Jet->py(tau1);
  double q = pxMiss / Jet->px(tau1) - pyMiss / Jet->py(tau1);
  double a2 = 0;  //q2/q;
  q2 = Jet->px(tau1) / Jet->px(tau2) - Jet->py(tau1) / Jet->py(tau2);
  q = pxMiss / Jet->px(tau2) - pyMiss / Jet->py(tau2);
  double a1 = 0;  //q2/q;
  tau1p = (1 + a1) * Jet->p(tau1);
  tau2p = (1 + a2) * Jet->p(tau2);
  return m(tau1p, tau2p);
}

void MyEvent::FatEff(Hist *F2bpT, Hist *F2bm, Hist *FMDpT, Hist *FMDm,
    Hist *F1tpT, Hist *F1tm, Hist *FinpT, Hist *Finm) {
  Jet->GetFatSp(F2bpT, F2bm, 1);
  Jet->GetFatSp(FMDpT, FMDm, 2);
  Jet->GetFatSp(F1tpT, F1tm, 3);
  Jet->GetFatSp(FinpT, Finm, 4);

}

void MyEvent::Spec(Hist *h, int sp, double w) {
  int max = this->MaxI(sp);
  int temp = 0;
  for (int i = 0; i < max; i++) {
    h->fill(isJet(sp) ? Jet->pT(I(sp, i)) : pythia->event[I(sp, i)].pT(), w);
  }
  max = this->MaxI(-sp);
  temp = 0;
  for (int i = 0; i < max; i++) {
    h->fill(isJet(-sp) ? Jet->pT(I(-sp, i)) : pythia->event[I(-sp, i)].pT(), w);
  }
}

void MyEvent::eta(Hist *h, int sp, double w) {
  int max = this->MaxI(sp);
  int temp = 0;
  for (int i = 0; i < max; i++) {
    h->fill(isJet(sp) ? Jet->eta(I(sp, i)) : pythia->event[I(sp, i)].eta(), w);
  }
  max = this->MaxI(-sp);
  temp = 0;
  for (int i = 0; i < max; i++) {
    h->fill(isJet(-sp) ? Jet->eta(I(-sp, i)) : pythia->event[I(-sp, i)].eta(),
        w);
  }
}

void MyEvent::InvMlJJ(Hist *h) {
  for (unsigned int i = 0; i < lJJ.size(); i += 3) {
    h->fill(
        this->m(pythia->event[lJJ[i]].p(), Jet->p(lJJ[i + 1]),
            Jet->p(lJJ[i + 2])), 1);
  }
}

void MyEvent::InvM(Hist *h, int sp1, double w) {
  int max1 = this->MaxI(sp1);
  int temp1 = 0;
  for (int i1 = 0; i1 < max1; i1++) {
    temp1 = this->I(sp1, i1);
    h->fill(this->m(isJet(sp1) ? Jet->p(temp1) : pythia->event[temp1].p()), w);
  }
}

void MyEvent::InvM(Hist *h, int sp1, int sp2, double w) {
  int max1 = this->MaxI(sp1);
  int temp1 = 0;
  int max2 = this->MaxI(sp2);
  int temp2 = 0;
  for (int i1 = 0; i1 < max1; i1++) {
    for (int i2 = 0; i2 < max2; i2++) {
      temp1 = this->I(sp1, i1);
      temp2 = this->I(sp2, i2);
      if (!Same(temp1, temp2, sp1, sp2)) {   //temp1 != temp2  || sp1 != sp2
        h->fill(
            this->m((isJet(sp1) ? Jet->p(temp1) : pythia->event[temp1].p()),
                (isJet(sp2) ? Jet->p(temp2) : pythia->event[temp2].p())), w);
      }
    }
  }
}

void MyEvent::InvM(Hist *h, int sp1, int sp2, int sp3, double w) {
  int max1 = this->MaxI(sp1);
  int temp1 = 0;
  int max2 = this->MaxI(sp2);
  int temp2 = 0;
  int max3 = this->MaxI(sp3);
  int temp3 = 0;
  for (int i1 = 0; i1 < max1; i1++) {
    for (int i2 = 0; i2 < max2; i2++) {
      for (int i3 = 0; i3 < max3; i3++) {
        temp1 = this->I(sp1, i1);
        temp2 = this->I(sp2, i2);
        temp3 = this->I(sp3, i3);
        if (!Same(temp1, temp2, sp1, sp2) && !Same(temp1, temp3, sp1, sp3)
            && !Same(temp2, temp3, sp2, sp3)) {
          h->fill(
              this->m((isJet(sp1) ? Jet->p(temp1) : pythia->event[temp1].p()),
                  (isJet(sp2) ? Jet->p(temp2) : pythia->event[temp2].p()),
                  (isJet(sp3) ? Jet->p(temp3) : pythia->event[temp3].p())), w);
        }
      }
    }
  }
}

void MyEvent::InvM(Hist *h, int sp1, int sp2, int sp3, int sp4, double w) {
  int max1 = this->MaxI(sp1);
  int temp1 = 0;
  int max2 = this->MaxI(sp2);
  int temp2 = 0;
  int max3 = this->MaxI(sp3);
  int temp3 = 0;
  int max4 = this->MaxI(sp4);
  int temp4 = 0;
  for (int i1 = 0; i1 < max1; i1++) {
    for (int i2 = 0; i2 < max2; i2++) {
      for (int i3 = 0; i3 < max3; i3++) {
        for (int i4 = 0; i4 < max4; i4++) {
          temp1 = this->I(sp1, i1);
          temp2 = this->I(sp2, i2);
          temp3 = this->I(sp3, i3);
          temp4 = this->I(sp4, i4);
          if (!Same(temp1, temp2, sp1, sp2) && !Same(temp1, temp3, sp1, sp3)
              && !Same(temp2, temp3, sp2, sp3) && !Same(temp3, temp4, sp3, sp4)
              && !Same(temp1, temp4, sp1, sp4)
              && !Same(temp2, temp4, sp2, sp4)) {
            h->fill(
                this->m((isJet(sp1) ? Jet->p(temp1) : pythia->event[temp1].p()),
                    (isJet(sp2) ? Jet->p(temp2) : pythia->event[temp2].p()),
                    (isJet(sp3) ? Jet->p(temp3) : pythia->event[temp3].p()),
                    (isJet(sp4) ? Jet->p(temp4) : pythia->event[temp4].p())),
                w);
          }
        }
      }
    }
  }
}

void MyEvent::InvMpTMiss(Hist *h, double w) {
  vector<int> jet100 = Jet->GetJet100();
  Vec4 ptemp;
  double ax = 0;
  double ay = 0;
  double bx = 0;
  double by = 0;
  double temp = 0;
  double A = 0;
  double B = 0;
  if (jet100.size() > 1) {
    for (unsigned int i1 = 0; i1 < jet100.size() - 1; i1++) {
      for (unsigned int i2 = i1 + 1; i2 < jet100.size(); i2++) {
        ax = Jet->px(i1);
        ay = Jet->py(i1);
        bx = Jet->px(i2);
        by = Jet->py(i2);
        temp = sqrt(ax * ax + ay * ay);
        ax = ax / temp;
        ay = ay / temp;
        temp = sqrt(bx * bx + by * by);
        bx = bx / temp;
        by = by / temp;
        A = (ax * pxMiss + ay * pyMiss
            - (bx * pxMiss + by * pyMiss) * (ax * bx + ay * by))
            / (1 - (ax * bx + ay * by) * (ax * bx + ay * by));
        B = (bx * pxMiss + by * pyMiss
            - (ax * pxMiss + ay * pyMiss) * (ax * bx + ay * by))
            / (1 - (ax * bx + ay * by) * (ax * bx + ay * by));
        ptemp = Jet->p(i1);
        temp = m(ptemp) * m(ptemp) + 2 * A * (ptemp.e() - ptemp.pAbs());
        h->fill(sqrt(abs(temp)), w);
        ptemp = Jet->p(i2);
        temp = m(ptemp) * m(ptemp) + 2 * B * (ptemp.e() - ptemp.pAbs());
        h->fill(sqrt(abs(temp)), w);
      }
    }
  }
}

double MyEvent::mTC() {
  double temp = 0;
  for (unsigned int i = 0; i < Jet20.size(); i++) {
    if (abs(Jet->eta(Jet20[i])) < 2) {
      temp += Jet->pT(Jet20[i]);
    }
  }
  for (unsigned int i = 0; i < em.size(); i++) {
    if (abs(pythia->event[em[i]].eta()) < 2) {
      temp += pythia->event[em[i]].pT();
    }
  }
  for (unsigned int i = 0; i < ep.size(); i++) {
    if (abs(pythia->event[ep[i]].eta()) < 2) {
      temp += pythia->event[ep[i]].pT();
    }
  }
  for (unsigned int i = 0; i < mum.size(); i++) {
    if (abs(pythia->event[mum[i]].eta()) < 2) {
      temp += pythia->event[mum[i]].pT();
    }
  }
  for (unsigned int i = 0; i < mup.size(); i++) {
    if (abs(pythia->event[mup[i]].eta()) < 2) {
      temp += pythia->event[mup[i]].pT();
    }
  }
  return temp;
}

int MyEvent::NlepJet(double cutoff) {
  return Jet->Nlep(cutoff);
}
int MyEvent::NJet(double cutoff) {
  return Jet->NJet(cutoff);
}
double MyEvent::MJet() {
  return Jet->MJet();
}
double MyEvent::Spher() {
  return sph->sphericity();
}
double MyEvent::LinSpher() {
  return linsph->sphericity();
}
double MyEvent::Thru() {
  return thr->thrust();
}
double MyEvent::TranSph() {
  double mat[4] = { 0, 0, 0, 0 };
  for (int i = 0; i < pythia->event.size(); i++) {
    if (pythia->event[i].isFinal() && pythia->event[i].isVisible()) {
      mat[0] += pythia->event[i].px() * pythia->event[i].px();
      mat[1] += pythia->event[i].px() * pythia->event[i].py();
      mat[2] += pythia->event[i].py() * pythia->event[i].py();
      mat[3] += pythia->event[i].py() * pythia->event[i].py()
          + pythia->event[i].px() * pythia->event[i].px();
    }
  }
  if (mat[3] == 0) {
    return 0;
  } else {
    mat[0] = mat[0] / mat[3];
    mat[1] = mat[1] / mat[3];
    mat[2] = mat[2] / mat[3];
  }
  return 1 - sqrt(1 + 4. * mat[1] * mat[1] - 4. * mat[0] * mat[2]);
}

bool MyEvent::isJet(int sp) {
  switch (sp) {
  case 11:
  case -11:
  case 13:
  case -13:
  case 22:
  case 500:
    return false;
    break;
  case 15:
  case -15:
  case 5:
  case -5:
  case 55:
  case 0:
    return true;
    break;
  default:
    return false;

  }
}

int MyEvent::I(int sp, int i) {
  switch (sp) {
  case 11:
    return em[i];
    break;
  case -11:
    return ep[i];
    break;
  case 13:
    return mum[i];
    break;
  case -13:
    return mup[i];
    break;
  case 22:
    return gamma[i];
    break;
  case 15:
    return tauJet[i];
    break;
  case -15:
    return taupJet[i];
    break;
  case 5:
    return bJet[i];
    break;
  case -5:
    return bpJet[i];
    break;
  case 500:
    return Bm[i];
    break;
  case 55:
    return Fatbb[i];
    break;
  case 0:
    return Jet20[i];
    break;
  default:
    return 0;
  }
}

int MyEvent::MaxI(int sp) {
  switch (sp) {
  case 11:
    return em.size();
    break;
  case -11:
    return ep.size();
    break;
  case 13:
    return mum.size();
    break;
  case -13:
    return mup.size();
    break;
  case 22:
    return gamma.size();
    break;
  case 15:
    return tauJet.size();
    break;
  case -15:
    return taupJet.size();
    break;
  case 5:
    return bJet.size();
    break;
  case -5:
    return bpJet.size();
    break;
  case 500:
    return Bm.size();
    break;
  case 55:
    return Fatbb.size();
    break;
  case 0:
    return Jet20.size();
    break;
  default:
    return 0;
  }
}

bool MyEvent::AcceptLepton(int lep) {
  if (pythia->event[lep].pT() < 5. || pythia->event[lep].status() < 0) {
    return false;
  }
  double Etemp = 0;
  for (int i = 0; i < pythia->event.size(); i++) {
    if (pythia->event[i].status() > 0 && pythia->event[i].isVisible()
        && i != lep) {
      if (deltaR(lep, i) < 0.3) {
        Etemp += pythia->event[i].pT();
      }
    }
  }
  int CloseJet = 0;
  //  for(int i = 0; i < Jet->NJet(); i++){if (deltaRlJ(lep,i) < 0.4 && Jet->pT(i) > 10 && !Jet->lep(i)) {CloseJet++;}} 
  if (Etemp < 0.1 * pythia->event[lep].pT() && CloseJet == 0) {
    return true;
  } else {
    return false;
  }
}

double MyEvent::deltaR(int i1, int i2) {
  return sqrt(
      pow(pythia->event[i1].eta() - pythia->event[i2].eta(), 2)
          + phi2(pythia->event[i1].phi() - pythia->event[i2].phi()));
}

double MyEvent::deltaRlJ(int i1, int i2) {
  return sqrt(
      pow(pythia->event[i1].eta() - Jet->eta(i2), 2)
          + phi2(pythia->event[i1].phi() - Jet->phi(i2)));
}

double MyEvent::deltaRJJ(int i1, int i2) {
  return sqrt(
      pow(Jet->eta(i1) - Jet->eta(i2), 2) + phi2(Jet->phi(i1) - Jet->phi(i2)));
}

double MyEvent::sign(double d) {
  if (d >= 0) {
    return 1;
  } else {
    return -1;
  }
}

double MyEvent::phi2(double phi) {
  phi = abs(phi);
  if (phi > pi) {
    phi -= 2 * pi;
  }
  return pow(phi, 2);
}

bool MyEvent::Same(int temp1, int temp2, int sp1, int sp2) {
  if (temp1 != temp2) {
    return false;
  } else {
    if (sp1 == sp2) {
      return true;
    } else {
      if ((sp1 == 0 && isJet(sp2)) || (sp2 == 0 && isJet(sp1))) {
        return true;
      } else {
        return false;
      }
    }
  }
}

double MyEvent::m(const Vec4& v1) {
  double m2 = pow2(v1.e()) - pow2(v1.px()) - pow2(v1.py()) - pow2(v1.pz());
  return (m2 > 0.) ? sqrt(m2) : 0.;
}

double MyEvent::m(const Vec4& v1, const Vec4& v2) {
  double m2 = pow2(v1.e() + v2.e()) - pow2(v1.px() + v2.px())
      - pow2(v1.py() + v2.py()) - pow2(v1.pz() + v2.pz());
  return (m2 > 0.) ? sqrt(m2) : 0.;
}

double MyEvent::m(const Vec4& v1, const Vec4& v2, const Vec4& v3) {
  double m2 = pow2(v1.e() + v2.e() + v3.e()) - pow2(v1.px() + v2.px() + v3.px())
      - pow2(v1.py() + v2.py() + v3.py()) - pow2(v1.pz() + v2.pz() + v3.pz());
  return (m2 > 0.) ? sqrt(m2) : 0.;
}

double MyEvent::m(const Vec4& v1, const Vec4& v2, const Vec4& v3,
    const Vec4& v4) {
  double m2 = pow2(v1.e() + v2.e() + v3.e() + v4.e())
      - pow2(v1.px() + v2.px() + v3.px() + v4.px())
      - pow2(v1.py() + v2.py() + v3.py() + v4.py())
      - pow2(v1.pz() + v2.pz() + v3.pz() + v4.pz());
  return (m2 > 0.) ? sqrt(m2) : 0.;
}

int MyEvent::Bmes(int i) {
  int id = pythia->event[i].id();
  if (abs(id) < 500) {
    return id;
  }
  int signid = id / abs(id);
  id = abs(id);
  id = ((int) id / 100) % 10;
  if (id == 5) {
    return signid * 5;
  }
  return 0;
}
