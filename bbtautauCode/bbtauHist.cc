using namespace std;
using namespace Pythia8; 

class bbtauHist{
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

  Hist bbtauhlSub;
  Hist a1a1hlSub;
  Hist h1a1a1hlSub;
  Hist a1h1hlSub; 

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
  bbtauHist(string model, double fact);
  ~bbtauHist();
  void Fill(MyProcess *Mp, MyEvent *Me);
  void AddEvent();
  long GetN();
  
  void Save();


};

bbtauHist::bbtauHist(string model, double fact){

  bSp.book("",100,0.,100.);    bSp.null();
  BmSp.book("",100,0.,100.);   BmSp.null();
  bJSp.book("",100,0.,100.);   bJSp.null();
  FjSp.book("",100,0.,100.);   FjSp.null();

  Fjm.book("",100,0.,100.);    Fjm.null();

  Nb.book("",15,0.5,15.5);     Nb.null();
  NBm.book("",15,0.5,15.5);    NBm.null();
  NbJ.book("",15,0.5,15.5);    NbJ.null();
  NFj.book("",15,0.5,15.5);    NFj.null();

  bbtauHad.book("",100,0.,150.);   bbtauHad.null();
  bbHad.book("",100,0.,100.);      bbHad.null();
  tataHad.book("",100,0.,150.);    tataHad.null();
  a1a1Had.book("",100,0.,150.);    a1a1Had.null();
  h1a1a1Had.book("",100,0.,100.);  h1a1a1Had.null();
  a1h1Had.book("",100,0.,150.);    a1h1Had.null();

  bbtauSub.book("",100,0.,150.);   bbtauSub.null();
  a1a1Sub.book("",100,0.,150.);    a1a1Sub.null();
  h1a1a1Sub.book("",100,0.,100.);  h1a1a1Sub.null();
  a1h1Sub.book("",100,0.,150.);    a1h1Sub.null();

  bbtauhlSub.book("",100,0.,150.); bbtauhlSub.null();
  a1a1hlSub.book("",100,0.,150.);  a1a1hlSub.null();
  h1a1a1hlSub.book("",100,0.,100.);h1a1a1hlSub.null();
  a1h1hlSub.book("",100,0.,150.);  a1h1hlSub.null();


  Fat2bpT.book("",100,0.,200.);   Fat2bpT.null();
  Fat2bm.book("",100,0.,200.);    Fat2bm.null();
  FatMDpT.book("",100,0.,200.);   FatMDpT.null();
  FatMDm.book("",100,0.,200.);    FatMDm.null();
  Fat1tpT.book("",100,0.,200.);   Fat1tpT.null();
  Fat1tm.book("",100,0.,200.);    Fat1tm.null();
  FatinpT.book("",100,0.,200.);   FatinpT.null();
  Fatinm.book("",100,0.,200.);    Fatinm.null();

  name = model;
  N = 0;
  ScaleFactor = fact;
}


bbtauHist::~bbtauHist(){
}

void bbtauHist::Fill(MyProcess *Mp, MyEvent *Me){
  Mp->Spec(&bSp,5,1);   
  Me->Spec(&bJSp,5,1);  
  Me->Spec(&BmSp,500,1);
  Me->Spec(&FjSp,55,1);

  Me->InvM(&Fjm,55,1);

  Nb.fill(Mp->MaxI(5)+Mp->MaxI(-5),1);
  NbJ.fill(Me->MaxI(5)+Me->MaxI(-5),1);
  NBm.fill(Me->MaxI(500),1);
  NFj.fill(Me->MaxI(55),1);
  
  N++;

  Me->FatEff(&Fat2bpT,&Fat2bm,&FatMDpT,&FatMDm,&Fat1tpT,&Fat1tm,&FatinpT,&Fatinm);
  
  Me->bbtauSub(&bbtauSub,&a1a1Sub,&h1a1a1Sub,&a1h1Sub,10);
  Me->bbtauhlSub(&bbtauhlSub,&a1a1hlSub,&h1a1a1hlSub,&a1h1hlSub,10);
  if (Me->bbtau(&bbtauHad,&a1a1Had,&h1a1a1Had,&a1h1Had,10)){ 
    Me->InvM(&bbHad,5,-5,1);    Me->InvM(&bbHad,5,5,0.5);    Me->InvM(&bbHad,-5,-5,0.5);
    Me->InvM(&tataHad,15,-15,1);    Me->InvM(&tataHad,15,15,0.5);    Me->InvM(&tataHad,-15,-15,0.5);
  }
}

void bbtauHist::AddEvent(){N++;}

long bbtauHist::GetN(){return N;}

void bbtauHist::Save(){

  SaveHist(&bSp,"bbtaudata/bSp_"+name,1.);
  SaveHist(&BmSp,"bbtaudata/BmSp_"+name,1.);
  SaveHist(&bJSp,"bbtaudata/bJSp_"+name,1.);
  SaveHist(&FjSp,"bbtaudata/FjSp_"+name,1.);

  SaveHist(&Fjm,"bbtaudata/Fjm_"+name,1.);

  SaveHist(&Nb,"bbtaudata/Nb_"+name,1.);
  SaveHist(&NBm,"bbtaudata/NBm_"+name,1.);
  SaveHist(&NbJ,"bbtaudata/NbJ_"+name,1.);
  SaveHist(&NFj,"bbtaudata/NFj_"+name,1.);

  SaveHist(&bbtauHad,"bbtaudata/bbtauH_"+name,1.5);
  SaveHist(&bbHad,"bbtaudata/bbH_"+name,1.);
  SaveHist(&tataHad,"bbtaudata/tataH_"+name,1.5);
  SaveHist(&a1a1Had,"bbtaudata/a1a1H_"+name,1.5);
  SaveHist(&h1a1a1Had,"bbtaudata/h1a1a1H_"+name,1.);
  SaveHist(&a1h1Had,"bbtaudata/a1h1H_"+name,1.5);

  SaveHist(&bbtauSub,"bbtaudata/bbtauS_"+name,1.5);
  SaveHist(&a1a1Sub,"bbtaudata/a1a1S_"+name,1.5);
  SaveHist(&h1a1a1Sub,"bbtaudata/h1a1a1S_"+name,1.);
  SaveHist(&a1h1Sub,"bbtaudata/a1h1S_"+name,1.5);

  SaveHist(&bbtauhlSub,"bbtaudata/bbtauShl_"+name,1.5);
  SaveHist(&a1a1hlSub,"bbtaudata/a1a1Shl_"+name,1.5);
  SaveHist(&h1a1a1hlSub,"bbtaudata/h1a1a1Shl_"+name,1.);
  SaveHist(&a1h1hlSub,"bbtaudata/a1h1Shl_"+name,1.5);

  SaveHist(&Fat2bpT,"bbtaudata/Fat2bpT_"+name,2.);
  SaveHist(&Fat2bm,"bbtaudata/Fat2bm_"+name,2.);
  SaveHist(&FatMDpT,"bbtaudata/FatMDpT_"+name,2.);
  SaveHist(&FatMDm,"bbtaudata/FatMDm_"+name,2.);
  SaveHist(&Fat1tpT,"bbtaudata/Fat1tpT_"+name,2.);
  SaveHist(&Fat1tm,"bbtaudata/Fat1tm_"+name,2.);
  SaveHist(&FatinpT,"bbtaudata/FatinpT_"+name,2.);
  SaveHist(&Fatinm,"bbtaudata/Fatinm_"+name,2.);

}

void bbtauHist::SaveHist(Hist *h, string file, double delta){
  *h *= ScaleFactor/(delta*(double)N);
  ofstream OutFile;
  OutFile.open(file.c_str());
  h->table(OutFile);
  OutFile.close();
}
