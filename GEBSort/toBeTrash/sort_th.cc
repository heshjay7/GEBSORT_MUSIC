#include <iostream>
#include "TFile.h"
#include "stdio.h"
#include <fstream>
#include "string.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TCutG.h"
#include "TChain.h"
//#include "functions_dfma.h"

using namespace std;

#define newtrees 1

#define DSSDTRLEN 256
#define MAXNUMGE 15
#define MAXNUMSI 10

void sort_th()

{
  ofstream myfile;
  myfile.open ("new.txt");

  //   TFile *f = new TFile("sorted_th.root","recreate"); // all runs
    TFile *f = new TFile("./TREE_FILES/Sum_1_95.root","recreate"); 
  //  TFile *f = new TFile("./TREE_FILES/Sum37.root","recreate");
  TTree tr("tr","254Rf");

  int k1=10;
  int k2=10; 
  float Eaf[6],Ear[6],Ea[6],Ta[6],Ei,Etot[6];
  unsigned short int Na,T[256];
  Short_t Traf[256]={0},Trab[256]={0},Traba[256]={0},Traba2[256]={0},Tradb[256]={0},Tradba[256]={0};
  int N; 
  long long int Eg[15],Ebox[6][10],Egc[15];
  unsigned long long int dTg[15],Tg[15],Tgc[15];
  unsigned long long int Ts,Tds;
  int Geid[15],Nsi,Nge,Nbox[6],Nsi2,Geidc[15],Ngec;
  long long int Ptofl,Pdetof,Ptofr;
  float Pde;

  tr.Branch("N", &N,  "N/I");
  tr.Branch("Na", &Na,  "Na/s");
  tr.Branch("Nsi", &Nsi,  "Nsi/I");
  tr.Branch("Nsi2", &Nsi2,  "Nsi2/I");
  tr.Branch("Nbox", Nbox[6],  "Nbox[6]/L");
  tr.Branch("Nge", &Nge,  "Nge/I");
  //  tr.Branch("Ngec", &Ngec,  "Ngec/i");
  tr.Branch("Ei", &Ei,  "Ei/F");
  tr.Branch("Ts", &Ts,  "Ts/l");
  tr.Branch("Tds", &Tds,  "Tds/l");
  tr.Branch("Tg", Tg,  "Tg[15]/l");
  tr.Branch("Tgc", Tgc,  "Tgc[15]/l");
  //  tr.Branch("dTg", dTg,  "dTg[15]/l");
  tr.Branch("Eaf", Eaf,  "Eaf[6]/F");
  tr.Branch("Ear", Ear,  "Ear[6]/F");
  tr.Branch("Ea", Ea,  "Ea[6]/F");
  tr.Branch("Etot", Etot,  "Etot[6]/F");
  tr.Branch("Ta", Ta,  "Ta[6]/F");
  tr.Branch("Ebox", Ebox,  "Ebox[6][10]/F");
  tr.Branch("Eg", Eg, "Eg[15]/L");
  tr.Branch("Geid", Geid, "Geid[15]/I");
  tr.Branch("Egc", Egc, "Egc[15]/L");
  tr.Branch("Geidc", Geidc, "Geidc[15]/I");
  tr.Branch("Pde", &Pde, "Pde/F");
  tr.Branch("Ptofl", &Ptofl, "Ptofl/L");
  tr.Branch("Ptofr", &Ptofr, "Ptofr/L");
  tr.Branch("Pdetof", &Pdetof, "Pdetof/L");
  // tr.Branch("Traf", Traf,  "Traf[256]/S");
  //  tr.Branch("Trab", Trab,  "Trab[256]/S");
  //  tr.Branch("Traba", Traba,  "Traba[256]/S");
  //  tr.Branch("Traba2", Traba2,  "Traba2[256]/S");
  //  tr.Branch("Tradb", Tradb,  "Tradb[256]/S");
  //  tr.Branch("Tradba", Tradba,  "Tradba[256]/S");
  //tr.Branch("T", T,  "T[256]/s");

 TChain *chain = new TChain("tree");

  char fth[100];
  int ith;
  /*   
  for(ith=19;ith<=26;ith++)
    {
      //      if(ith==36||ith==37||(ith>=53&&ith<=60)||ith==63||ith==67||ith==68||(ith>=88&&ith<=90)||ith==93||ith==95||ith==99||(ith>=46&&ith<=49)||ith==118){continue;}
     if(ith==21||ith==22||ith==24){continue;}
    sprintf(fth,"./TREE_FILES/new_presort.run%d.root",ith);
    chain->Add(fth);
    }

    chain->Add("./TREE_FILES/new_presort.run30.root"); 
    chain->Add("./TREE_FILES/new_presort.run31.root");
   chain->Add("./TREE_FILES/new_presort.run37.root");

  for(ith=39;ith<=48;ith++)
    {
    sprintf(fth,"./TREE_FILES/new_presort.run%d.root",ith);
    chain->Add(fth);
    }
  

   chain->Add("./TREE_FILES/new_presort.run54.root");
   chain->Add("./TREE_FILES/new_presort.run61.root");
   chain->Add("./TREE_FILES/new_presort.run66.root");*/
  //   chain->Add("./TREE_FILES/new_presort.run51.root");
  /*
 for(ith=55;ith<=65;ith++)
    {
    sprintf(fth,"./TREE_FILES/new_presort.run%d.root",ith);
    chain->Add(fth);
    }   
  */
  //   chain->Add("./TREE_FILES/new_presort.run37.root");
        
for(ith=1;ith<=62;ith++)
    {
    sprintf(fth,"./TREE_FILES/new_presort.run%d.root",ith);
    chain->Add(fth);
    }

for(ith=68;ith<=95;ith++)
    {
    sprintf(fth,"./TREE_FILES/new_presort.run%d.root",ith);
    chain->Add(fth);
    }
   
struct recoil_type {

   unsigned long long int ts;
   int en;
   int enb;
   int pu;
   float left;
   float right;
   signed long long int ptofl;
   signed long long int ptofr;
   signed long long int pdetof;
   float pde;
   double x;
   int d2t0;
   int d2t1;
   int d2t2;
   int nge;
   long long int geehi[MAXNUMGE];   
   unsigned long long int tge[MAXNUMGE];
   int getid[MAXNUMGE];
   int traceLen;
   short int trace_fr[DSSDTRLEN];
   short int trace_ba[DSSDTRLEN];

};

struct decay_type {

   unsigned long long int ts;
   int en;
   int enb;
   int pu_fr;
   int pu_ba;
   unsigned long long int time;
   int traceLen;
   short int trace_fr[DSSDTRLEN];
   short int trace_ba[DSSDTRLEN];
   int nge;
   long long int geehi[MAXNUMGE];   
   unsigned long long int tge[MAXNUMGE];
   int getid[MAXNUMGE];
   int nsi;
   long long int siehi[MAXNUMSI];  
   unsigned long long int tsi[MAXNUMSI];
   int sitid[MAXNUMSI];
   int d2t0;
   int d2t1;
   int d2t2;
};


// correlated chain


struct chain_type {
  int s_fr;
  int s_ba;
  recoil_type recoil;
  int ndec;
  decay_type decay[6];
  int corr_type;
};


chain_type chain_event;


chain->SetBranchAddress("s_fr",&chain_event.s_fr);
chain->SetBranchAddress("s_ba",&chain_event.s_ba);
chain->SetBranchAddress("recoil",&chain_event.recoil);
chain->SetBranchAddress("ndec",&chain_event.ndec);
chain->SetBranchAddress("dec1",&chain_event.decay[0]);
chain->SetBranchAddress("dec2",&chain_event.decay[1]);
chain->SetBranchAddress("dec3",&chain_event.decay[2]);
//chain->SetBranchAddress("dec3",&chain_event.decay[2]);
chain->SetBranchAddress("dec4",&chain_event.decay[3]);
chain->SetBranchAddress("dec5",&chain_event.decay[4]);
chain->SetBranchAddress("dec6",&chain_event.decay[5]);
//chain->SetBranchAddress("dec7",&chain_event.decay[6]);
//chain->SetBranchAddress("dec8",&chain_event.decay[7]);
//chain->SetBranchAddress("dec9",&chain_event.decay[8]);
//chain->SetBranchAddress("dec10",&chain_event.decay[9]);

// variables

// TH1F *Ea = new TH1F("Ea","Ea",15000,0,15000);
 TH2D *h2_traces = new TH2D("traces","traces",300,0,300,2000,0,2000);
 TH2D *h2_traces2 = new TH2D("traces2","traces2",300,0,300,100,0,100);
 TH1F *hs = new TH1F("hs","test",1000,0,1000);

 float logt, logt2;

 int dec_pu_traces_c=0;
 int dec_pud2_traces_c=0;


 int i,j,k,t,ii;
short int  s,d;
/*
 FILE *file;
 char str[255];
 float a[161],b[161];
 file =fopen("F.cal","r");
 fgets(str,255,file);
 fgets(str,255,file);
 for(i=1;i<161;i++)
   {
     fscanf(file,"%d",&j);
     if(j==i)fscanf(file,"%f %f",&a[j],&b[j]);
     else {return -3;}
    fgets(str,255,file); 
   }
 fclose(file);
*/
 // EVENT LOOP
 
 for (int i=0; i<chain->GetEntries();i++) {
//for(i=2150964;i<chain->GetEntries();i++){
//for (int i=0; i<1000;i++) {
   
   chain->GetEntry(i);
   //printf ("Event number %d\n",i);
   /*     
   
   if(chain_event.decay[0].en>0&&chain_event.decay[0].en<1000&&0.4343*log(10.0*chain_event.decay[0].time)>0&&0.4343*log(10.0*chain_event.decay[0].time)<4.5&&i>70826)
     {
       for(j=0;j<256-k1;j++)
	 {
	   //           s = chain_event.recoil.trace_fr[j] & 0x3FFF;
	   //	   s = chain_event.recoil.trace_fr[j];
	   s=chain_event.decay[0].trace_fr[j];
	   d=chain_event.decay[0].trace_fr[j+k1]-s;	  
	   h2_traces->Fill(j, s );
	   h2_traces2->Fill(j, d );
	     }
   cout<<i<<endl;
   h2_traces->Draw();
   h2_traces2->Draw();
   return;
     }
   */   
   s=0;  
   N=i;
   Na=chain_event.ndec;
   Nsi=chain_event.s_fr;
   Nsi2=chain_event.s_ba;
   Nge=chain_event.recoil.nge;
   Ngec=chain_event.decay[0].nge;
   if(chain_event.recoil.en>chain_event.recoil.enb)
   {Ei=chain_event.recoil.en;}
      else {Ei=chain_event.recoil.enb;}
    Ts=chain_event.recoil.ts;
    Tds=chain_event.decay[0].ts;
    Pde= chain_event.recoil.pde;
    Ptofl= chain_event.recoil.ptofl;
    Ptofr= chain_event.recoil.ptofr;
    Pdetof= chain_event.recoil.pdetof;
        for(j=0;j<6;j++)
	  {
	    if(j<Na)
	      {
	      Eaf[j]=chain_event.decay[j].en;
	      Ear[j]=chain_event.decay[j].enb;
	      if(Eaf[j]>Ear[j]){Ea[j]=Eaf[j];}
	      else {Ea[j]=Ear[j];}
	      Ta[j]=log10(10.0*chain_event.decay[j].time);
	      Nbox[j]=chain_event.decay[j].nsi;
	      Etot[j]=0;
	      for(k=0;k<10;k++)
		{
		  if(k<Nbox[j])
		    { Ebox[j][k]=chain_event.decay[j].siehi[k];
		      if(Etot[j]==0){Etot[j]=Ebox[j][k];}
		      else if(Etot[j]<Ebox[j][k]){Etot[j]=Ebox[j][k];}
		    }
		  else
		    { Ebox[j][k]=0;}
		}
	      Etot[j]=Etot[j]+Ea[j];
	      }
	    else
	      {
	      Eaf[j]=0;
	      Ear[j]=0;
	      Ea[j]=0;
	      Ta[j]=0;
	      Nbox[j]=0;
	      for(k=0;k<10;k++)
		{
		   Ebox[j][k]=0;
		}	
	      }
	  }
	for(j=0;j<15;j++)
	  {
	    if(j<Nge)
	      {
	    Eg[j]=chain_event.recoil.geehi[j];
	    Geid[j]=chain_event.recoil.getid[j];
	    Tg[j]=chain_event.recoil.tge[j];
	    //	    dTg[j]=chain_event.recoil.tge[j]-chain_event.recoil.ts;
	    //	    if(Tg[j]-Ts>-7&&Tg[j]-Ts<6&&Eg[j]>0&&Eg[j]<1000)
	    //	      {hs->Fill(Eg[j]);}
	      }
	    else
	      {
	    Eg[j]=0;
	    Geid[j]=0;
	    Tg[j]=0;	
	      }
	  }

	for(j=0;j<15;j++)
	  {
	    if(j<Ngec)
	      {
	    Egc[j]=chain_event.decay[0].geehi[j];
	    Geidc[j]=chain_event.decay[0].getid[j];
	    Tgc[j]=chain_event.decay[0].tge[j];
	    //	    dTg[j]=chain_event.recoil.tge[j]-chain_event.recoil.ts;
	      }
	    else
	      {
	    Egc[j]=0;
	    Geidc[j]=0;
	    Tgc[j]=0;	
	      }
	  }

	

	/*
	if(Ea[0]>0&&Ea[0]<400&&Ta[0]<9.5&&Ta[0]>6.5)
	  {
	    s=1;
	  }
	else if(Ea[0]>0&&Ea[0]<1000&&Ta[0]<=6.5&&Ta[0]>0)
	  {
	    s=1;
	  }
	else
	  {
	    for(j=0;j<6;j++)
	      {
		if((Ea[j]>8000&&Ea[j]<8150)||(Ea[j]>7300&&Ea[j]<7500))
		  {
		    s=1;
		    break;
		  }
	      }
	  }
	if(s==1)
	  {
        for(j=0;j<223-k1-k2;j++)
	  {
	   //	     Traf[j]=chain_event.recoil.trace_fr[j];
	   Trab[j]=chain_event.recoil.trace_ba[j];
	   Traba[j]=chain_event.recoil.trace_ba[j+k1]-Trab[j];
	   Traba2[j]=(chain_event.recoil.trace_ba[j+k1+k2]-chain_event.recoil.trace_ba[j+k2])-Traba[j];
	   T[j]=j;
	  }
	*/
	//	if(N==297570){cout<<chain_event.decay[0].time<<endl;}
	tr.Fill(); 
	//	  }
	/*	
	if(Ea[0]>0&&Ea[0]<1000&&Ta[0]>0&&Ta[0]<4.5)
	  {
	   for(j=0;j<223-k1-k2;j++)
	     {
	       Tradb[j]=chain_event.decay[0].trace_ba[j];
	       Tradba[j]=chain_event.decay[0].trace_ba[j+k1]-Trab[j];
	     }	     
	  }
	*/
   

   /*  
   if(chain_event.decay[0].en>8000&&chain_event.decay[0].en<8150&&chain_event.ndec>1)
     {
       Ea->Fill(chain_event.decay[1].en);
     }
   */
   

   /*
   for(j=0;j<=5;j++)
     {
       if(chain_event.decay[j].en>8000&&chain_event.decay[j].en<8150&&chain_event.ndec>j+1)
	 {
	   for(k=j+1;k<chain_event.ndec;k++)
	     {
	      Ea->Fill(chain_event.decay[k].en);
	     }
	 }
     }
   */
   
  } // event loop

 //printf("fastpileup1: %i, fastpileup2: %i, fastpileup3: %i\n",fastpileup1, fastpileup2, fastpileup3);

 printf("total events =  %d\n",i);


 // hs->Draw();
 f->cd();
 tr.Write();
 // f.Close();
   // f->Write();
 // f->Close();

} // macro end

// *********************************************************************


