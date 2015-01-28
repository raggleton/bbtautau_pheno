using namespace std;
using namespace Pythia8; 
using namespace fastjet;

class MyJet{
private:
  double R;
  Pythia *pythia;
  JetDefinition *jetDef;
  vector <PseudoJet> Input;
  vector <PseudoJet> Jets;
  vector <PseudoJet> FatbJ;
  vector <PseudoJet> Fat1tag;
  vector <PseudoJet> Fat2b;
  vector <PseudoJet> FatMD;
  vector <int> taub;
  vector <int> tauJet;
  vector <int> taupJet;
  vector <int> bJet;
  vector <int> bpJet;
  vector <int> Jet20;
  vector <int> Jet100;
  vector <int> FatbbJet;
  bool IsHadron(int id);
  double deltaR(int i1, int i2);
  double deltaR(PseudoJet J1, int i2);
  double deltaR(PseudoJet J1, PseudoJet J2);
  double phi2(double phi);
  int Bmes(int i);
  int Jettype(PseudoJet *J, ClusterSequence *CS);
  void SubJets(ClusterSequence *CS);
  bool SameJet(PseudoJet J1, PseudoJet J2);
  bool Between(PseudoJet J, PseudoJet J1, PseudoJet J2);


public:
  MyJet(Pythia *p, double Rparam);
  ~MyJet();
  double MJet();
  int Nlep(double cutoff);
  int NJet(double cutoff);
  int NJet();
  int part(int i);
  bool lep(int i);
  double eta(int i);
  double phi(int i);
  double pT(int i);
  double e(int i);
  double px(int i);
  double py(int i);
  Vec4 p(int i);
  vector <double> DeltaRfj;
  vector <int> GettauJet();
  vector <int> GettaupJet();
  vector <int> GetbJet();
  vector <int> GetbpJet();
  vector <int> GetJet20();
  vector <int> GetJet100();
  vector <int> GetFatbb();
  void GetFatSp(Hist *pT, Hist *m, int dist);
};


MyJet::MyJet(Pythia *p, double Rparam){
  R = Rparam;
  pythia = p;         Input.resize(0);     Jets.resize(0);   FatbJ.resize(0);   taub.resize(0);     
  tauJet.resize(0);   taupJet.resize(0);   bJet.resize(0);   bpJet.resize(0);   Jet20.resize(0);   Jet100.resize(0);
  FatbbJet.resize(0); Fat1tag.resize(0);   Fat2b.resize(0);  FatMD.resize(0);   DeltaRfj.resize(0);
  for(int i = 0; i < pythia->event.size(); i++){
    if (abs(pythia->event[i].id()) == 15 && abs(pythia->event[i].status()) == 23) {taub.push_back(i);}
    if (Bmes(i) == 5) {taub.push_back(i);}
  }
  for(int i = 0; i < pythia->event.size(); i++){
    if (pythia->event[i].isFinal() && IsHadron(pythia->event[i].id()) && abs(pythia->event[i].eta()) < 2.5){ 
      PseudoJet temp(pythia->event[i].px(),pythia->event[i].py(),pythia->event[i].pz(),pythia->event[i].e());
      for(unsigned int itau = 0; itau < taub.size(); itau++){
	if (pythia->event.isAncestor(i,taub[itau]) && deltaR(i,taub[itau]) < R) {temp.set_user_index(Bmes(taub[itau]));}
      }
      if (abs(abs(pythia->event[i].id())-12) == 1) {temp.set_user_index(pythia->event[i].id());}
      Input.push_back(temp);
    }
  }
  //  double Rparam = 0.2;
  jetDef = new fastjet::JetDefinition(cambridge_algorithm, 1.2, E_scheme, Best);
  cout << "1" << flush;
  ClusterSequence clust(Input, *jetDef);
  Jets = clust.inclusive_jets();  

  SubJets(&clust);
  cout << "2" << flush;
  jetDef = new fastjet::JetDefinition(kt_algorithm, R, E_scheme, Best);
  Jets.resize(0);
  ClusterSequence clustSeq(Input, *jetDef);
  Jets = clustSeq.inclusive_jets();  


  cout << "3" << flush;
  for(unsigned int i = 0; i < Jets.size(); i++){
    Jets[i].set_user_index(Jettype(&Jets[i],&clustSeq));
    if (Jets[i].perp() >= 20 && Jets[i].user_index() == 0) {Jet20.push_back(i);}
    if (Jets[i].perp() >= 100 && abs(Jets[i].user_index()) <  10) {Jet100.push_back(i);}
  }
  for(unsigned int i = 0; i < Jets.size(); i++){if (Jets[i].perp() > 15){
    switch ((int)Jets[i].user_index()){
    case 15:    tauJet.push_back(i);     break;
    case -15:   taupJet.push_back(i);    break;
    case 5:     bJet.push_back(i);       break;
    case -5:    bpJet.push_back(i);      break;
    }
  }}
  cout << "4" << flush;
  int iJs = Jets.size();
  FatbbJet.resize(0);
  for (int i = 0; i < FatbJ.size(); i++){
    FatbbJet.push_back(iJs+i);
    Jets.push_back(FatbJ[i]);
  }
  cout << "5" << flush;
}


MyJet::~MyJet(){
  delete jetDef;
}

double MyJet::MJet(){
  double Max = 0;
  for (unsigned int i = 0; i < Jets.size(); i++) if (Jets[i].perp() > Max && !lep(i)) Max = Jets[i].perp();
  return Max;
}

int    MyJet::Nlep(double cutoff){
  int Nl = 0;
  for (unsigned int i = 0; i < Jets.size(); i++) if (Jets[i].perp() > cutoff && lep(i)) Nl++;
  return Nl;
}
int    MyJet::NJet(double cutoff){
  int Nj = 0;
  for (unsigned int i = 0; i < Jets.size(); i++) if (Jets[i].perp() > cutoff && !lep(i)) Nj++;
  return Nj;
}
int    MyJet::NJet(){      return Jets.size();}
int    MyJet::part(int i){ return Jets[i].user_index();}
bool   MyJet::lep(int i){  if (abs(abs(Jets[i].user_index())-12) == 1){return true;} else {return false;}}
//double MyJet::y1(int i){   return dmerge()*R*R/Jets[i].m2();
double MyJet::eta(int i){  return Jets[i].eta();}
double MyJet::phi(int i){  return Jets[i].phi_std();}
double MyJet::pT(int i){   return Jets[i].perp();}
double MyJet::e(int i){    return Jets[i].e();}
double MyJet::px(int i){   return Jets[i].px();}
double MyJet::py(int i){   return Jets[i].py();}
Vec4   MyJet::p(int i){    return Vec4(Jets[i].px(),Jets[i].py(),Jets[i].pz(),Jets[i].e());}
vector <int> MyJet::GettauJet(){   return tauJet;}
vector <int> MyJet::GettaupJet(){  return taupJet;}
vector <int> MyJet::GetbJet(){     return bJet;}
vector <int> MyJet::GetbpJet(){    return bpJet;}
vector <int> MyJet::GetJet20(){    return Jet20;}
vector <int> MyJet::GetJet100(){   return Jet100;}
vector <int> MyJet::GetFatbb(){    return FatbbJet;}

void MyJet::GetFatSp(Hist *pT, Hist *m, int dist){
  vector <PseudoJet> temp;
  temp.resize(0);
  switch (dist){
  case 1: temp = Fat2b;   break;
  case 2: temp = FatMD;   break;
  case 3: temp = Fat1tag; break;
  case 4: temp = FatbJ;
  }
  for (int i = 0; i < temp.size(); i++){
    pT->fill(temp[i].perp());
    m->fill(temp[i].m());
  }
}


bool MyJet::IsHadron(int id){
  switch (id){
  case 12:      case -12:      case 14:      case -14:     case 16:    case -16:  return false; break;
  case 1000022: case -1000022: case 1000039: case -1000039:                       return false; break;
    // case 11:      case -11:      case 13:      case -13:                            return false; break;   
  default: return true;
  }
}

double MyJet::deltaR(int i1, int i2){
  return sqrt(pow(pythia->event[i1].eta() - pythia->event[i2].eta(),2) + phi2(pythia->event[i1].phi() - pythia->event[i2].phi()));
}

double MyJet::deltaR(PseudoJet J1, int i2){
  return sqrt(pow(J1.eta() - pythia->process[i2].eta(),2) + phi2(J1.phi() - pythia->process[i2].phi()));
}

double MyJet::deltaR(PseudoJet J1, PseudoJet J2){
  return sqrt(pow(J1.eta() - J2.eta(),2) + phi2(J1.phi() - J2.phi()));
}


double MyJet::phi2(double phi){
  phi = abs(phi);
  if (phi > pi) {phi -= 2*pi;}
  return pow(phi,2);
} 

int MyJet::Bmes(int i){
  int id = pythia->event[i].id();
  if (abs(id) < 500) {return id;}
  int signid = id/abs(id);
  id = abs(id);
  id = ((int)id/100) % 10;
  if (id == 5) {return signid*5;}
  return 0;
}

int MyJet::Jettype(PseudoJet *J, ClusterSequence *CS){
  int ret = 0;
  vector <PseudoJet> tautemp = CS->constituents(*J);
  double taumx = 0;   double taupx = 0;   double bmx = 0;   double bpx = 0;
  double taumy = 0;   double taupy = 0;   double bmy = 0;   double bpy = 0;
  double emx = 0;     double epx = 0;     double mumx = 0;  double mupx = 0;
  double emy = 0;     double epy = 0;     double mumy = 0;  double mupy = 0;
  for(unsigned int iCo = 0; iCo < tautemp.size(); iCo++){ 
    switch (tautemp[iCo].user_index()){
    case 15:   taumx += tautemp[iCo].px();   taumy += tautemp[iCo].py();    break;
    case -15:  taupx += tautemp[iCo].px();   taupy += tautemp[iCo].py();    break;
    case 5:    bmx   += tautemp[iCo].px();   bmy   += tautemp[iCo].py();    break;
    case -5:   bpx   += tautemp[iCo].px();   bpy   += tautemp[iCo].py();    break;
    case 11:   emx   += tautemp[iCo].px();   emy   += tautemp[iCo].py();    break;
    case -11:  epx   += tautemp[iCo].px();   epy   += tautemp[iCo].py();    break;
    case 13:   mumx  += tautemp[iCo].px();   mumy  += tautemp[iCo].py();    break;
    case -13:  mupx  += tautemp[iCo].px();   mupy  += tautemp[iCo].py();    break;
    }
  }
  if (taumx*taumx+taumy*taumy > 0.36*J->kt2()) {ret = 15;}
  if (taupx*taupx+taupy*taupy > 0.36*J->kt2()) {ret = -15;}
  if (bmx*bmx+bmy*bmy >         0.36*J->kt2()) {ret = 5;}
  if (bpx*bpx+bpy*bpy >         0.36*J->kt2()) {ret = -5;}
  if (emx*emx+emy*emy >         0.49*J->kt2()) {ret = 11;}
  if (epx*epx+epy*epy >         0.49*J->kt2()) {ret = -11;}
  if (mumx*mumx+mumy*mumy >     0.49*J->kt2()) {ret = 13;}
  if (mupx*mupx+mupy*mupy >     0.49*J->kt2()) {ret = -13;}

  return ret;
}

void MyJet::SubJets(ClusterSequence *CS){
  PseudoJet p;
  PseudoJet p1;
  PseudoJet p2;
  double muP = 0.67;
  double mu  = 1.0;
  double yP  = 0.09;
  double y   = 0.0;
  double mj,mj1,mj2;
  for (int i = 0; i < Jets.size(); i++){ if (Jets[i].pt() > 30. && Jets[i].m() > 12.){
      p1 = Jets[i];
      int bNum = 0;
      for (int i1 = 0; i1 < pythia->process.size(); i1++){
	double DRtemp = deltaR(p1,i1);
	bool closest = true;
	if (abs(pythia->process[i1].id()) == 5 && DRtemp < 1.2 && abs(pythia->process[i1].eta()) < 2.5) {
	  for (int i2 = 0; i2 < Jets.size(); i2++){ if (DRtemp > deltaR(Jets[i2],i1) && i != i2) { closest = false;}}
	  if (closest) bNum++;
	}
      }
      if (bNum >= 2) {Fat2b.push_back(p1); }
      do{
	p = p1;
	if (p.has_parents(p1,p2)){
	  mj  = p.m();
	  mj1 = p1.m();
	  mj2 = p2.m();
	  if (mj > 0.0){
	    mu = min(mj1,mj2)/mj;
	    y  = min(p1.pt2(),p2.pt2())*deltaR(p1,p2)/pow(mj,2);
	  }else{
	    mu = 0.0;
	    y  = 1.0;
	  }
	}else {
	  mu = 0.0;
	  y  = 1.0;
	}
      }while (mu > muP || y < yP);
      if (mu != 0.0 || y != 1.0){
	bNum = 0;
	for (int i1 = 0; i1 < pythia->process.size(); i1++){
	  double DRtemp = deltaR(p,i1);
	  bool closest = true;
	  if (abs(pythia->process[i1].id()) == 5 && DRtemp < 1.2 && abs(pythia->process[i1].eta()) < 2.5) {
	    for (int i2 = 0; i2 < Jets.size(); i2++){ if (DRtemp > deltaR(Jets[i2],i1) && i != i2) { closest = false;}}
	    if (closest) bNum++;
	  }
	}
	if (bNum >= 2) {FatMD.push_back(p); }
	double Rp = max(min(0.3,deltaR(p1,p2)/2.0),0.2);
	cout << "1.3" << flush;
	JetDefinition *JD = new fastjet::JetDefinition(cambridge_algorithm, Rp, E_scheme, Best);
	vector <PseudoJet> temp;
	ClusterSequence clustSeq(CS->constituents(p), *JD);
	temp = clustSeq.inclusive_jets();
	temp = sorted_by_pt(temp);
	cout << "1.4" << flush;
	if (temp.size() >= 3){
	  p = join(temp[0],temp[1],temp[2]);
	  if ((abs(Jettype(&temp[0],&clustSeq)) == 5 || abs(Jettype(&temp[1],&clustSeq)) == 5) && p.m() > 12.) {
	    Fat1tag.push_back(p);
	  }
	  if (abs(Jettype(&temp[0],&clustSeq)) == 5 && abs(Jettype(&temp[1],&clustSeq)) == 5 && p.m() > 12.) {
	    DeltaRfj.push_back(deltaR(temp[0],temp[1]));
	    FatbbJet.push_back(i);
	    FatbJ.push_back(p);
	    cout << "1.5" << flush;
	    vector <PseudoJet> con; 
	    for (int ico = 0; ico < 3; ico++){
	      con.resize(0);
	      con = clustSeq.constituents(temp[ico]);
	      cout << "consize" << con.size()<< "|"<< (clustSeq.constituents(temp[0])).size() << flush;
	      for (int i1 = 0; i1 < con.size(); i1++){
		for (int inp = 0; inp < Input.size(); inp++){
		  if (SameJet(con[i1],Input[inp])){
		    Input.erase(Input.begin()+inp);
		    cout << "1.6" << flush;
		    break;
		  }
		}
	      }
	    }
	  }
	  cout << "1.7" << flush;
	}
      }
    }
  }
}

bool MyJet::SameJet(PseudoJet J1, PseudoJet J2){
  if (J1.pz() == J2.pz() && J1.px() == J2.px() && J1.py() == J2.py()){return true;} else {return false;}
}


bool MyJet::Between(PseudoJet J, PseudoJet J1, PseudoJet J2){
  double dR1 = deltaR(J,J1);
  double dR2 = deltaR(J,J2);
  double dR  = deltaR(J1,J2);
  if (dR1 < dR && dR2 < dR){return true;} else {return false;}
}
