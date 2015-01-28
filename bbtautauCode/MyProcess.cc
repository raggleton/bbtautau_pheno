// #include <iostream>
// #include <fstream>
// #include <string>
// #include "Pythia.h"
// #include "Simul.cc"
using namespace std;
using namespace Pythia8; 
//const double pi = 3.14159265358979323846;

class MyProcess{
private:
  vector <int> em;
  vector <int> ep;
  vector <int> mum;
  vector <int> mup;
  vector <int> tau;
  vector <int> taup;
  vector <int> b;
  vector <int> bp;
  Pythia *pythia;
  int I(int sp, int i);
  bool AcceptLepton(int lep);
  double deltaR(int i1, int i2);
  double sign(double d);
  double phi2(double phi);
  bool Same(int temp1, int temp2, int sp1, int sp2);
  double m(const Vec4& v1);
  double m(const Vec4& v1, const Vec4& v2);
  double m(const Vec4& v1, const Vec4& v2, const Vec4& v3);
  double m(const Vec4& v1, const Vec4& v2, const Vec4& v3, const Vec4& v4);
public:
  MyProcess(Pythia *p);
  ~MyProcess();
  bool Fourb(Hist *b4, Hist *a1a1, Hist *h1a1a1, Hist *a1h1,  double bres);
  double DRaa();
  bool Zbb(Hist *bbh1, Hist *Za1, double bres);
  void Zsp(Hist *pT, Hist *eta, Hist *boost, Hist *a1bst);
  void Spec(Hist *h, int sp, double w);
  void InvM(Hist *h, int sp1, double w);
  void InvM(Hist *h, int sp1, int sp2, double w);
  void InvM(Hist *h, int sp1, int sp2, int sp3, double w);
  void InvM(Hist *h, int sp1, int sp2, int sp3, int sp4, double w);
  int MaxI(int sp);
  int pT20lep;
  int pT15lep;
  double pTmaxlep;
  double pTMiss;
  double pxMiss;
  double pyMiss;
};


MyProcess::MyProcess(Pythia *p){
  pythia = p;
  for (int i = 0; i < pythia->process.size(); i ++){ if (abs(pythia->process[i].eta()) > 2.5) {pythia->process[i].statusNeg();}}
  double px = 0;
  double py = 0;
  pT15lep = 0;
  for(int i = 0; i < pythia->process.size(); i++){
    switch (pythia->process[i].id()){
    case 5:  if (AcceptLepton(i)) {b.push_back(i);}
      if (pythia->process[i].pT() > 15. && pythia->process[i].status() > 0) {pT15lep++;}
      break;
    case -5:  if (AcceptLepton(i)) {bp.push_back(i);}
      if (pythia->process[i].pT() > 15. && pythia->process[i].status() > 0) {pT15lep++;}
      break;
    case 11:   if (AcceptLepton(i)) {em.push_back(i);}
      if (pythia->process[i].pT() > 15. && pythia->process[i].status() > 0) {pT15lep++;}
      break;
    case -11:  if (AcceptLepton(i)) {ep.push_back(i);}
      if (pythia->process[i].pT() > 15. && pythia->process[i].status() > 0) {pT15lep++;}
      break;
    case 13:   if (AcceptLepton(i)) {mum.push_back(i);}
      if (pythia->process[i].pT() > 15. && pythia->process[i].status() > 0) {pT15lep++;}
      break;
    case -13:  if (AcceptLepton(i)) {mup.push_back(i);}
      if (pythia->process[i].pT() > 15. && pythia->process[i].status() > 0) {pT15lep++;}
      break;
    case 15:  if (AcceptLepton(i)) {tau.push_back(i);}
      if (pythia->process[i].pT() > 15. && pythia->process[i].status() > 0) {pT15lep++;}
      break;
    case -15:  if (AcceptLepton(i)) {taup.push_back(i);}
      if (pythia->process[i].pT() > 15. && pythia->process[i].status() > 0) {pT15lep++;}
      break;
    case 12: case -12: case 14: case -14: case 16: case -16: case 1000022: case -1000022: 
      if (pythia->process[i].status() > 0) {px += pythia->process[i].px();  py += pythia->process[i].py();}
      break;
    }
  }
  pTMiss = sqrt(px*px+py*py);
  pxMiss = px;     pyMiss = py;
}

MyProcess::~MyProcess(){
}

bool MyProcess::Fourb(Hist *b4, Hist *a1a1, Hist *h1a1a1, Hist *a1h1,  double bres){
  int bArray[4];
  if (b.size() + bp.size() != 4) {return false;} else {
    for (int i = 0; i < b.size(); i++){ bArray[i] = b[i];}
    for (int i = 0; i < bp.size(); i++){ bArray[b.size()+i] = bp[i];}
  }
  double Inv4b = m(pythia->process[bArray[0]].p(),pythia->process[bArray[1]].p(),pythia->process[bArray[2]].p(),pythia->process[bArray[3]].p());
  b4->fill(Inv4b);
  int tempA[2][3] = {{2,1,1},{3,3,2}};
  double Inv1,Inv2;
  for (int i = 1; i < 4; i++){
    Inv1 = m(pythia->process[bArray[0]].p(),pythia->process[bArray[i]].p());
    Inv2 = m(pythia->process[bArray[tempA[0][i-1]]].p(),pythia->process[bArray[tempA[1][i-1]]].p());
    if (abs(Inv1 - 125) < bres) {a1h1->fill(Inv2);}
    if (abs(Inv2 - 125) < bres) {a1h1->fill(Inv1);}
    if (abs(Inv1 - Inv2) < bres) {
      a1a1->fill(Inv1);     a1a1->fill(Inv2);
      if (abs(Inv4b - 125) < bres) {h1a1a1->fill(Inv1);   h1a1a1->fill(Inv2);}
    }
  }
  return true;
}

double MyProcess::DRaa(){
  vector <int> a1;
  a1.resize(0);
  for (int i = 0; i < pythia->process.size(); i++){if (pythia->process[i].id() == 36) {a1.push_back(i);}}
  if (a1.size() == 2) {return deltaR(a1[0],a1[1]);} else {return 0;}  
}

bool MyProcess::Zbb(Hist *bbh1, Hist *Za1, double bres){
  int bArray[2];
  if (b.size() + bp.size() != 2) {return false;} else {
    for (int i = 0; i < b.size(); i++){ bArray[i] = b[i];}
    for (int i = 0; i < bp.size(); i++){ bArray[b.size()+i] = bp[i];}
  }
  int Zint = -1;
  for (int i = 0; i < pythia->process.size(); i++){if (pythia->process[i].id() == 23) {Zint = i;}}
  if (Zint == -1) {return false;}
  double Inv2b = m(pythia->process[bArray[0]].p(),pythia->process[bArray[1]].p());
  double InvZa1 = m(pythia->process[Zint].p(),pythia->process[bArray[0]].p(),pythia->process[bArray[1]].p());
  Za1->fill(InvZa1);
  if (abs(InvZa1 - 125) < bres) {bbh1->fill(Inv2b);}
  return true;
}

void MyProcess::Zsp(Hist *pT, Hist *eta, Hist *boost, Hist *a1bst){
  int Zint = -1;
  for (int i = 0; i < pythia->process.size(); i++){if (pythia->process[i].id() == 23) {Zint = i;}}
  if (Zint >= 0) {
    pT->fill(pythia->process[Zint].pT());
    eta->fill(pythia->process[Zint].eta());
    boost->fill(pythia->process[Zint].e()/pythia->process[Zint].m());
  }
  Zint = -1;
  for (int i = 0; i < pythia->process.size(); i++){if (pythia->process[i].id() == 36) {Zint = i;}}
  if (Zint >= 0) {
    a1bst->fill(pythia->process[Zint].e()/pythia->process[Zint].m());
  }
}



void MyProcess::Spec(Hist *h, int sp, double w){
  int max = this->MaxI(sp);  int temp = 0;
  for (int i = 0; i < max; i++){h->fill(pythia->process[I(sp,i)].pT(),w);}
  max = this->MaxI(-sp);  temp = 0;
  for (int i = 0; i < max; i++){h->fill(pythia->process[I(-sp,i)].pT(),w);}
}


void MyProcess::InvM(Hist *h, int sp1, double w){
  int max1 = this->MaxI(sp1);  int temp1 = 0;
  for (int i1 = 0; i1 < max1; i1++){
    temp1 = this->I(sp1,i1);
    h->fill(this->m(pythia->process[temp1].p()),w);
  }
}  


void MyProcess::InvM(Hist *h, int sp1, int sp2, double w){
  int max1 = this->MaxI(sp1);  int temp1 = 0;
  int max2 = this->MaxI(sp2);  int temp2 = 0;
  for (int i1 = 0; i1 < max1; i1++){
    for (int i2 = 0; i2 < max2; i2++){
      temp1 = this->I(sp1,i1);   temp2 = this->I(sp2,i2);
      if (!Same(temp1,temp2,sp1,sp2)) {   //temp1 != temp2  || sp1 != sp2
	h->fill(this->m(pythia->process[temp1].p(),pythia->process[temp2].p()),w);
      }
    }
  }
}  



void MyProcess::InvM(Hist *h, int sp1, int sp2, int sp3, double w){
  int max1 = this->MaxI(sp1);  int temp1 = 0;
  int max2 = this->MaxI(sp2);  int temp2 = 0;
  int max3 = this->MaxI(sp3);  int temp3 = 0;
  for (int i1 = 0; i1 < max1; i1++){
    for (int i2 = 0; i2 < max2; i2++){
      for (int i3 = 0; i3 < max3; i3++){
	temp1 = this->I(sp1,i1);   temp2 = this->I(sp2,i2);    temp3 = this->I(sp3,i3); 
	if (!Same(temp1,temp2,sp1,sp2) && !Same(temp1,temp3,sp1,sp3) && !Same(temp2,temp3,sp2,sp3)) {
	  h->fill(this->m( pythia->process[temp1].p(),pythia->process[temp2].p(),pythia->process[temp3].p()),w);
	}
      }
    }
  }
}  

void MyProcess::InvM(Hist *h, int sp1, int sp2, int sp3, int sp4, double w){
  int max1 = this->MaxI(sp1);  int temp1 = 0;
  int max2 = this->MaxI(sp2);  int temp2 = 0;
  int max3 = this->MaxI(sp3);  int temp3 = 0;
  int max4 = this->MaxI(sp4);  int temp4 = 0;
  for (int i1 = 0; i1 < max1; i1++){
    for (int i2 = 0; i2 < max2; i2++){
      for (int i3 = 0; i3 < max3; i3++){
	for (int i4 = 0; i4 < max4; i4++){
	  temp1 = this->I(sp1,i1);   temp2 = this->I(sp2,i2);    temp3 = this->I(sp3,i3);    temp4 = this->I(sp4,i4); 
	  if (!Same(temp1,temp2,sp1,sp2) && !Same(temp1,temp3,sp1,sp3) && !Same(temp3,temp4,sp3,sp4) && !Same(temp3,temp2,sp3,sp2) && !Same(temp2,temp4,sp2,sp4) && !Same(temp3,temp4,sp3,sp4)) {
	    h->fill(this->m(pythia->process[temp1].p(),pythia->process[temp2].p(),pythia->process[temp3].p(),pythia->process[temp4].p()),w);
	  }
	}
      }
    }
  }
}  





int MyProcess::I(int sp, int i){
  switch (sp){
  case 11:   return em[i];      break;
  case -11:  return ep[i];      break;
  case 13:   return mum[i];     break;
  case -13:  return mup[i];     break;
  case 15:   return tau[i];  break;
  case -15:  return taup[i]; break;
  case 5:    return b[i];    break;
  case -5:   return bp[i];   break;
  default:   return 0;
  }
}
    

int MyProcess::MaxI(int sp){
  switch (sp){
  case 11:   return em.size();      break;
  case -11:  return ep.size();      break;
  case 13:   return mum.size();     break;
  case -13:  return mup.size();     break;
  case 15:   return tau.size();  break;
  case -15:  return taup.size(); break;
  case 5:    return b.size();    break;
  case -5:   return bp.size();   break;
  default:   return 0;
  }
}


bool MyProcess::AcceptLepton(int lep){
  if (pythia->process[lep].pT() < 15. || abs(pythia->process[lep].eta()) > 2.5) {return false;}
  return true;
}

double MyProcess::deltaR(int i1, int i2){
  return sqrt(pow(pythia->process[i1].eta() - pythia->process[i2].eta(),2) + phi2(pythia->process[i1].phi() - pythia->process[i2].phi()));
}


double MyProcess::sign(double d){
  if (d >= 0) {return 1;}
  else {return -1;}
}

double MyProcess::phi2(double phi){
  phi = abs(phi);
  if (phi > pi) {phi -= 2*pi;}
  return pow(phi,2);
} 

bool MyProcess::Same(int temp1, int temp2, int sp1, int sp2){
  if (temp1 != temp2){ return false;}else{
    if (sp1 == sp2){return true;}
  }
}

double MyProcess::m(const Vec4& v1) {
  double m2 = pow2(v1.e()) - pow2(v1.px()) - pow2(v1.py()) - pow2(v1.pz());
  return (m2 > 0.) ? sqrt(m2) : 0.; 
}

double MyProcess::m(const Vec4& v1, const Vec4& v2) {
  double m2 = pow2(v1.e() + v2.e()) - pow2(v1.px() + v2.px())
     - pow2(v1.py() + v2.py()) - pow2(v1.pz() + v2.pz());
  return (m2 > 0.) ? sqrt(m2) : 0.; 
}

double MyProcess::m(const Vec4& v1, const Vec4& v2, const Vec4& v3){
  double m2 = pow2(v1.e() + v2.e() + v3.e()) - pow2(v1.px() + v2.px() + v3.px())
    - pow2(v1.py() + v2.py() + v3.py()) - pow2(v1.pz() + v2.pz() + v3.pz());
  return (m2 > 0.) ? sqrt(m2) : 0.; 
}

double MyProcess::m(const Vec4& v1, const Vec4& v2, const Vec4& v3, const Vec4& v4){
  double m2 = pow2(v1.e() + v2.e() + v3.e() + v4.e()) - pow2(v1.px() + v2.px() + v3.px() + v4.px())
    - pow2(v1.py() + v2.py() + v3.py() + v4.py()) - pow2(v1.pz() + v2.pz() + v3.pz() + v4.pz());
  return (m2 > 0.) ? sqrt(m2) : 0.; 
}
