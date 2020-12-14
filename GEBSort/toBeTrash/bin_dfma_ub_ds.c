#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include "Rtypes.h"
#include "TROOT.h"
#include "TFile.h"
#include "TRandom.h"
#include "TH1.h"
#include "TH2.h"
#include "TObjArray.h"
#include "TObject.h"
#include "TKey.h"
#include "TSystem.h"
#include "TCutG.h"
#include "TTree.h"
#include "gdecomp.h"
#include "limits.h"

#include "GEBSort.h"
#include "GTMerge.h"

#define DSSDTRLEN 256
#define MCPTRLEN 10

//#define SICH 6

#define RCH  4
#define LCH  3
#define UPCH  5
#define DOWNCH  6
#define PPACDE 7
#define TOFCH 17
#define EMONRCH 8
#define WHEEL2CH 1  
#define EMONLCH 9
#define SICH 10
#define WHEELCH 2

#define EnergyFromTrace 0
#define SmoothOrNot 1
#define WITHDGS 1

// Change DSSD tid info here-> 

#define FRTIDLOW 1	//1
#define FRTIDHIGH 160	//160
#define BATIDLOW 161	//161
#define BATIDHIGH 320	//320
#define NUMFR 160	//160
#define NUMBA 160	//160
#define NUMSTRIPS 320	//320
#define NUMSTRIPSBOX 320	//320

/* Doppler correction and Ge calibrations */

#define NGSGE 110
#define MAXNUMGE 15
#define MAXNUMDEC 4
#define MAXNUMSI 10



struct strip_type {
   int phystrip;
   int thr;
   float off;
   float gain;
   int baseline;
};

struct strip_type map_fr[NUMSTRIPS+1];
struct strip_type map_ba[NUMSTRIPS+1];

struct strip_type map_box[NUMSTRIPSBOX+1];

struct clover_map {

  float gain;
  float off;

};

struct clover_map clmap[21];

struct recoil_pu_type {

  //short int trace[MAXTRACELEN];
   short int trace[DSSDTRLEN];
   unsigned short int traceLen;
   int s_fr;
   int s_ba;
   unsigned long long int ts;

};



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

struct decay_type_new {

   unsigned long long int ts;
   int en;
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
   int fdecEv;
   int d2t0;
   int d2t1;
   int s_ba;
   int s_fr;
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

struct focal_plane {

   int left;
   int right;
   int icde;
   int ppacde;
   int ppacde2;
   unsigned long long int left_ts;
   unsigned long long int right_ts;
   unsigned long long int icde_ts;
   unsigned long long int ppacde_ts;
   unsigned long long int ppacde2_ts;
 

};

struct chain_type {
   int s_fr;
   int s_ba;
   recoil_type recoil;
   int ndec;
  decay_type decay[6];//*4thGen
   int corr_type;//*26th Apr

};

struct pixel_type {

   int status;
   chain_type chain;

};

struct clover_type {

  int nge;
  int tid[MAXNUMGE];
  long long int ehi[MAXNUMGE];
  unsigned long long int tge[MAXNUMGE];

};

struct clover_type clover;

struct chain_type chain;
struct recoil_pu_type recoil_pu[100];
struct recoil_type recoil;
struct decay_type decay;
struct decay_type decay_fr;
struct decay_type decay_ba;
struct focal_plane fplane;
struct pixel_type dssd_corr[NUMFR+1][NUMBA+1];
struct pixel_type dssd_front[NUMFR+1];
struct pixel_type dssd_back[NUMBA+1];

int ncl=0;
int d_fr_basesample_av[161];
int d_fr_basesample_first[161];

unsigned long long int ts, t_first, t_last, t_firstge;
int CurSubEvNo;
bool first=true;
bool firstge = true;

int RealCount = 0;
int trn1 = 0;
int trn2 = 0;
int trn3 = 0;
int TestTr1 = 0;
int TestTr2 = 0;
int TestTr3 = 0;
int TestTr4 = 0;
int leTr = 0;
int numDFMA = 0;
int numDGS = 0;
// external parameters 

extern TFile *treef;//TREES...
extern TTree *tree;//TREES...

extern PARS Pars;
extern int ng;
extern DGSEVENT DGSEvent[MAXCOINEV];
extern int tlkup[NCHANNELS];
extern int tid[NCHANNELS];
extern int DFMAEvDecompose_v3 (unsigned int *ev, int len, DFMAEVENT * DFMAEvent); 
extern DFMAEVENT DFMAEvent[MAXCOINEV];

extern DGSEVENT XAEvent[MAXCOINEV];
extern int XAng;

// DEFINE HISTOGRAMS HERE!!!!

TH1D *h1_cl, *h1_cr, *h1_cup, *h1_cdown, *h1_esi, *h1_emonl, *h1_emonr, *h1_wheel_total, *h1_wheel_mon, *h1_wheel_left,  *h1_wheel_up, *h1_wheel_ge, *h1_wheel, *h1_wheel2, *h1_ppacde, *h1_ppacdeg, *h1_cx, *h1_cy, *h1_cxn, *h1_csum, *h1_csumy;
TH1D *h1_emonl_rate, *h1_emonl_rate2;
TH1D *h1_emonr_rate, *h1_emonr_rate2;
TH1D *h1_cl_pu, *h1_cr_pu, *h1_ppacde_pu;

TH2F *h2_clr, *h2_cud, *h2_cxy;
TH1D *h1_clg, *h1_crg, *h1_cxg, *h1_cxng;
TH1D *h1_cx_en_g1, *h1_cxn_en_g1;
TH1D *h1_cx_en_g2, *h1_cxn_en_g2;

TH2F *h2_xde, *h2_xdeg; 

TH2F *h2_clrg;
TH2F *h2_clrg_int;

// DSSD stuff

TH2F *h2_dssd_en, *h2_dssd_en_raw;
TH2F *h2_dssd_fb;
TH2F *h2_dssd_fr_emax, *h2_dssd_fr_emax_phy;
TH2F *h2_dssd_ba_emax, *h2_dssd_ba_emax_phy;

TH2F *h2_dssd_rate, *h2_FP_rate;
TH1D *h1_ndssd;

TH1D *h1_cl_int;
TH1D *h1_cr_int;
TH1D *h1_cxg_int;
TH1D *h1_cxng_int;
TH2F *h2_clr_int;
TH2F *h2_clrg_en;
TH2F *h2_clrg_en_g1, *h2_clrg_en_g2;
TH2F *h2_r_hitxy, *h2_d_hitxy, *h2_dssd_hitxy_phys, *h2_dssd_hitxy_physg, *h2_dssd_hitxy_physg251Md, *h2_dssd_hitxy_phys_corr, *h2_dssd_hitxy_phys_corrf;
TH2F *h2_dssd_fb_corr;
TH2F *h2_dssd_fr_p2;
TH2F *h2_dssd_ba_p2;
TH2F *h2_clr_p2;
TH2F *h2_dssd_fr_p1;
TH2F *h2_dssd_ba_p1;
TH2F *h2_clr_p1;

TH2F *h2_dssd_traces_fr, *h2_dssd_traces_ba;

TH1D *h1_trlen;

TH2F *h2_xehi, *h2_xehig;
TH2F *h2_x_d_fr_e, *h2_x_d_fr_e_short;


extern FILE *trf;

/*-----------------------------------------------------*/

int exit_dfma()
{

}

/*-----------------------------------------------------*/

int
sup_dfma ()
{

#if(1)
  // declarations \\

  char str1[STRLEN], str2[STRLEN];
  float pi;
  int i,j;

  TH1D *mkTH1D (char *, char *, int, double, double);
  TH2F *mkTH2F (char *, char *, int, double, double, int, double, double);

 char presortFile[100];
 char compstr[10];
 char onechar[10];

 t_first = 0;
 t_firstge = 0;
 t_last = 0;
 CurSubEvNo = 0;

//**************************\\
// Histogram initialisation \\
//**************************\\

// LIKE USERINIT.H!!!

//Focal plane

h1_ppacdeg = mkTH1D((char *)"ppacdeg",(char *)"ppacdeg",4000,0,40000);
h1_cl = mkTH1D((char *)"cl",(char *)"cl",4000,0,80000);
h1_cr = mkTH1D((char *)"cr",(char *)"cr",4000,0,80000);
h1_cx = mkTH1D((char *)"cx",(char *)"cx",4000,-100000,100000);
h1_cup = mkTH1D((char *)"cup",(char *)"cup",4000,0,80000);
h1_cdown = mkTH1D((char *)"cdown",(char *)"cdown",4000,0,80000);
h1_cy = mkTH1D((char *)"cy",(char *)"cy",4000,-100000,100000);

h1_cl_pu = mkTH1D((char *)"cl_pu",(char *)"cl_pu",10,0,10);
h1_cr_pu = mkTH1D((char *)"cr_pu",(char *)"cr_pu",10,0,10);

h1_csum = mkTH1D((char *)"csum",(char *)"csum",4000, 0,200000);
h1_csumy = mkTH1D((char *)"csumy",(char *)"csumy",4000, 0,200000);

h1_cxn = mkTH1D((char *)"cxn",(char *)"cxn",4000,0,4000);
h2_clr = mkTH2F((char *)"clr",(char *)"clr",2000,0,80000,2000,0,80000);
h2_cud = mkTH2F((char *)"cud",(char *)"cud",2000,0,80000,2000,0,80000);
h2_cxy = mkTH2F((char *)"cxy",(char *)"cxy",2000,-100000,100000,2000,-100000,100000);

h2_xde = mkTH2F((char *)"xde",(char *)"xde",1000,-100000,100000,1000,0,400000);
h2_xdeg = mkTH2F((char *)"xdeg",(char *)"xdeg",1000,-100000,100000,1000,0,400000);

// DSSD stuff

h2_dssd_fr_emax= mkTH2F((char *)"dssd_fr_emax",(char *)"dssd_fr_emax",5000,0,50000,400,0,400);
h2_dssd_ba_emax= mkTH2F((char *)"dssd_ba_emax",(char *)"dssd_ba_emax",5000,0,50000,400,0,400);

h2_dssd_fr_emax_phy= mkTH2F((char *)"dssd_fr_emax_phy",(char *)"dssd_fr_emax_phy",5000,0,500000,400,0,400);
h2_dssd_ba_emax_phy= mkTH2F((char *)"dssd_ba_emax_phy",(char *)"dssd_ba_emax_phy",5000,0,500000,400,0,400);


corr_time_short = 10000000000;
corr_time = 1000000000000;


//<><><><><><><><><><><><>\\
//        MAPFILES        \\ 
//<><><><><><><><><><><><>\\


FILE *fmap1;
char fmapname1[32];

int strip, phystrip, thr, baseline;
float off, gain;


fmap1=fopen("MAP_FILES/dssd_fr_calib_th.map","read");
//fmap1=fopen("MAP_FILES/dssd_fr_test.map","read");

if(fmap1==0) printf("Failed to read mapfile...");
if(fmap1!=0) printf("Read mapfile ok... ");

for(i=0;i<161;i++){

   fscanf(fmap1,"%d %d %d %f %f %i", &strip, &phystrip, &thr, &off, &gain, &baseline);

   if (phystrip!=0) {  

   map_fr[strip].thr = thr;
   map_fr[strip].phystrip = phystrip;
   map_fr[strip].off = off;
   map_fr[strip].gain = gain;  
   map_fr[strip].baseline = baseline;
   
   } else {

   map_fr[strip].thr = 100;
   map_fr[strip].phystrip = phystrip;
   map_fr[strip].off = 0;
   map_fr[strip].gain = 1.0;  
   map_fr[strip].baseline = 0;

  }

   printf("map front strip %d (phys strip %d) has thr %d, offset %f and gain %f and baseline %i\n",strip,phystrip,thr,off,gain,baseline);

}

fclose(fmap1);


fmap1=fopen("MAP_FILES/dssd_ba_calib_th.map","read");
//fmap1=fopen("MAP_FILES/dssd_ba_test.map","read");

 for(i=0;i<161;i++){

   fscanf(fmap1,"%d %d %d %f %f %i", &strip, &phystrip, &thr, &off, &gain, &baseline);
 

   if (phystrip!=0) {  

   map_ba[strip].thr = thr;
   map_ba[strip].phystrip = phystrip;
   map_ba[strip].off = off; 
   map_ba[strip].gain = gain; 
   map_ba[strip].baseline = baseline;

  } else {

   map_ba[strip].thr = 100;
   map_ba[strip].phystrip = phystrip;
   map_ba[strip].off = 0;
   map_ba[strip].gain = 1.0;  
   map_ba[strip].baseline = 0;

  }


  printf("map back strip %d (phys strip %d) has thr %d, offset %f and gain %f and baseline %i\n",strip,phystrip,thr,off,gain,baseline);

}

fclose(fmap1);


// window

// TFile *g = new TFile("eftof","read");
 TFile *g = new TFile("eftofG.root","read");

 // eftof = (TCutG*) g->Get("eftof");
 eftof = (TCutG*) g->Get("eftofGate"); 
 printf("Got the gate HLC drew/n");

 g->Close();
 
#endif

//******************\\
// MAP FILE map.dat \\
//******************\\

  char str[STRLEN];
  int i1, i2, i7, i8;
  FILE *fp;
  for (i = 0; i < NCHANNELS; i++)
    {
      tlkup[i] = NOTHING;
      tid[i] = NOTHING;
    };


  printf("MAPPING H\n");
  fp = fopen ("./map_z115.dat", "r");
  //fp = fopen ("./map_gsfma333.dat", "r");
 // fp = fopen ("./map_LBNL_final4.dat", "r");
  if (fp == NULL)
    {
      printf ("need a map file to run\n");
      printf ("do - g++ mkMap_LBNL.c -o mkMap_LBNL; ./mkMap_LBNL > map_LBNL.dat\n");
      exit (1);
    };

  printf ("\nmapping H\n");
  i2 = fscanf (fp, "\n%i %i %i %s", &i1, &i7, &i8, str);
  printf ("%i %i %i %s\n", i1, i7, i8, str);
  while (i2 == 4)
    {
      tlkup[i1] = i7;
      tid[i1] = i8;
      i2 = fscanf (fp, "\n%i %i %i %s", &i1, &i7, &i8, str);
      printf ("%i %i %i %s\n", i1, i7, i8, str);
    };
  fclose (fp);

  Pars.wlist = gDirectory->GetList ();
  Pars.wlist->Print ();



#endif

  //trf = fopen("pucm_traces.txt","w"); 
    //if (trf==NULL) (printf("TRACE FILE DOES NOT OPEN\n")); 

  load_ranges();

  for (i=0;i<161;i++) d_fr_basesample_av[i]=0;
  for (i=0;i<161;i++) d_fr_basesample_first[i]=1;

};

/* ----------------------------------------------------------------- */

int
DFMAEvDecompose_v3 (unsigned int *ev, int len, DFMAEVENT * DFMAEvent)
{

 /* firmware circa Sept 2014 */

  /* declarations */

  int i, k, i1;
  unsigned int ui0 = 0, ui1 = 0, ui2 = 0;
  unsigned int PRE_RISE_SUM = 0, POST_RISE_SUM = 0;
  int rawE;
  unsigned int t1 = 0, t2 = 0, t3 = 0, t4 = 0;
  unsigned long long int ulli1;


  if (Pars.CurEvNo <= Pars.NumToPrint){
  printf ("entered DFMAEvDecompose_v3:\n");
  printf ("\nmapping\n");


  for (i = 2010; i < 2020; i++)
    {
      printf("lkup %d tid %d\n",tlkup[i],tid[i]);

    };
  }




  /* swap the bytes */

  i = 0;
  while (i < len)
    {

      /* before 4 3 2 1 */
      /*        | | | | */
      /* after  1 2 3 4 */

      t1 = (*(ev + i) & 0x000000ff) << 24;
      t2 = (*(ev + i) & 0x0000ff00) << 8;
      t3 = (*(ev + i) & 0x00ff0000) >> 8;
      t4 = (*(ev + i) & 0xff000000) >> 24;
      *(ev + i) = t1 + t2 + t3 + t4;

      i++;
    }


  /* debug print */

  if (Pars.CurEvNo <= Pars.NumToPrint)
    {
      printf ("event len=%i (%i bytes) >\n", len, len * sizeof (unsigned int));
      for (i = 0; i < len; i++)
        {
          printf ("%3i[doc: %3i]: %12u, 0x%8.8x; ", i, i + 1, *(ev + i), *(ev + i));
          ui0 = 0x80000000;
          printf ("|");
          for (k = 0; k <= 31; k++)
            {
              if ((*(ev + i) & ui0) == ui0)
                printf ("1");
              else
                printf ("0");
              ui0 = ui0 / 2;
              if ((k + 1) % 4 == 0)
                printf ("|");
            };
          printf ("\n");
        };
    };

  /* extract IDs */

  DFMAEvent->chan_id = (*(ev + 0) & 0x0000000f);
  DFMAEvent->board_id = ((*(ev + 0) >> 4) & 0xfff);
  DFMAEvent->id = DFMAEvent->board_id * 10 + DFMAEvent->chan_id;

  /* store the type and type ID */

  DFMAEvent->tpe = tlkup[DFMAEvent->id];
  DFMAEvent->tid = tid[DFMAEvent->id];

  if (Pars.CurEvNo <= Pars.NumToPrint)
    {
      printf ("chan_id = %i, board_id=%i, id=%i\n", DFMAEvent->chan_id, DFMAEvent->board_id, DFMAEvent->id);
    }

  /* extract wheel */

  DFMAEvent->wheel=  ( *(ev + 6) & 0xffff0000 ) >> 16;

  /* extract the energy */

  PRE_RISE_SUM = *(ev + 7) & 0x00ffffff;

  if (Pars.CurEvNo <= Pars.NumToPrint)
    {
      printf ("PRE_RISE_SUM =  0x%8.8x %12i  ", PRE_RISE_SUM, PRE_RISE_SUM);
      ui0 = 0x80000000;
      printf ("|");
      for (k = 0; k <= 31; k++)
        {
          if ((PRE_RISE_SUM & ui0) == ui0)
            printf ("1");
          else
            printf ("0");
          ui0 = ui0 / 2;
          if ((k + 1) % 4 == 0)
            printf ("|");
        };
      printf ("\n");
    };

  ui1 = *(ev + 7) & 0xff000000;
  ui1 >>= 24;
  ui2 = *(ev + 8) & 0x0000ffff;
  ui2 <<= 8;
  POST_RISE_SUM = ui1 + ui2;

  if (Pars.CurEvNo <= Pars.NumToPrint)
    {
      printf ("POST_RISE_SUM = 0x%8.8x %12i  ", POST_RISE_SUM, POST_RISE_SUM);
      ui0 = 0x80000000;
      printf ("|");
      for (k = 0; k <= 31; k++)
        {
          if ((POST_RISE_SUM & ui0) == ui0)
            printf ("1");
          else
            printf ("0");
          ui0 = ui0 / 2;
          if ((k + 1) % 4 == 0)
            printf ("|");
        };
      printf ("\n");
    };

  /* note: POS/PRE_RISE_SUM only have 24 bits */
  /* so making ints of them is not a problem */

  rawE = (int) POST_RISE_SUM - (int) PRE_RISE_SUM;
  DFMAEvent->ehi = rawE/10;
  if (DFMAEvent->ehi<0)  DFMAEvent->ehi=-DFMAEvent->ehi;

  if (Pars.CurEvNo <= Pars.NumToPrint)
    printf ("rawE = 0x%8.8x %i, DFMAEvent->ehi= %i\n", rawE, rawE, DFMAEvent->ehi);

  /* extract the LED time stamp */

  DFMAEvent->LEDts = 0;
  DFMAEvent->LEDts = (unsigned long long int) *(ev + 1);
  ulli1 = (unsigned long long int) (*(ev + 2) & 0x0000ffff);
  ulli1 = (ulli1 << 32);
  DFMAEvent->LEDts += ulli1;

  /* extract pulse shap parameters 7/15/2016 assuming Event Type = 4 */

  /*
  m2_begin_sample = *(ev + 10) & 0x00003fff;
  m2last_end_sample = *(ev + 9) & 0x00003fff;

  m2last_begin_sample = ( *(ev + 10) & 0x3fff0000 ) >> 16;

  m1_begin_sample = *(ev + 11) & 0x00003fff;
  m1_end_sample = ( *(ev + 11) & 0x3fff0000 ) >> 16;
  */

  DFMAEvent->header_type = ((*(ev + 2) & 0x000f0000)  >> 16); 


  /* extract trace */

  
  DFMAEvent->traceLen=0;
  for (i=0; i<len-HDRLENINTS; i++) {

   DFMAEvent->trace[DFMAEvent->traceLen] = (unsigned long long int) (*(ev + 13 + i) & 0x0000ffff);
   //DFMAEvent->trace[DFMAEvent->traceLen] = (short int) (*(ev + 13 + i) & 0x0000ffff);
   DFMAEvent->traceLen++;
   DFMAEvent->trace[DFMAEvent->traceLen] = (unsigned long long int) (*(ev + 13 + i) >> 16 & 0x0000ffff);
   //DFMAEvent->trace[DFMAEvent->traceLen] = (short int) (*(ev + 13 + i) >> 16 & 0x0000ffff);
   DFMAEvent->traceLen++;

   

  }

  /* DS changes 3/1/2014 */

  /*PARSE*/  

  unsigned long long int prevTS;
  int baseline, basesample;
  int postrisebeg, postriseend, prerisebeg, preriseend; 
  int peaksample;
  int prerisesum;
  int postrisesum;
  int m2_last_beg, m2_last_end;

  // previous TS

  //prevTS = 0;
  // correctoin DS 7/15/2016
  //prevTS = (unsigned long long int) ( *(ev + 3) & 0xffff);
  prevTS = (unsigned long long int) (( *(ev + 3) & 0xffff0000) >> 16);
  ulli1 = (unsigned long long int) (*(ev + 4) );
  ulli1 = (ulli1 << 16);
  prevTS += ulli1;

  // baseline

  baseline = (int) ( *(ev + 5) & 0xffffff);


  // event type 3
  // m2_last_end

  m2_last_end = (int) ( *(ev + 9) & 0x3fff);

  // post rise begin sample
  //postrisebeg = (int) ( *(ev + 10) & 0x3fff);

 // event type 3
  // m2_last_beg

  postrisebeg = (int) ( *(ev + 10) & 0x3fff);

  // pre rise begin sample
  prerisebeg = (int) ( *(ev + 11) & 0x3fff);

  // post rise end sample
  m2_last_beg = (int) ( *(ev + 10) & 0x3fff0000);
 

  m2_last_beg = m2_last_beg >> 16; 

  // pre rise end sample

  preriseend = (int) ( *(ev + 11) & 0x3fff0000);
  
  preriseend = preriseend >> 16; 

  // peak sample

   peaksample = (int) ( *(ev + 12) & 0x3fff); 

  // base sample

  basesample = (int) ( *(ev + 12) & 0x3fff0000);
  
  basesample = basesample >> 16; 

  //prerisesum = PRE_RISE_SUM/400;
  
  prerisesum = PRE_RISE_SUM;
  postrisesum = POST_RISE_SUM;

  DFMAEvent->baseline = baseline;
  DFMAEvent->m2_last_beg = m2_last_beg;
  DFMAEvent->m2_last_end = m2_last_end;

  DFMAEvent->postrisebeg = postrisebeg;
  DFMAEvent->prerisebeg = prerisebeg;
  DFMAEvent->postriseend = postriseend;
  DFMAEvent->preriseend = preriseend;
  DFMAEvent->peaksample = peaksample;
  DFMAEvent->basesample = basesample;
  DFMAEvent->prevTS = prevTS;
  DFMAEvent->prerisesum = prerisesum;
  DFMAEvent->postrisesum = postrisesum;

  /* done */

  if (Pars.CurEvNo <= Pars.NumToPrint)
    printf ("exit DFMAEvDecompose_v3:\n");

  return (0);

}

//*********************************************************\\


int
bin_dfma (GEB_EVENT * GEB_event)
{

  if (Pars.CurEvNo <= 100){//Pars.NumToPrint){
    printf ("entered bin_dfma:\n");
  }


double dssdfrsi_t, dssdbasi_t;
double dssdfrsi_t2, dssdbasi_t2;

double ratecomp=1e8;

#if(1)

  //if(Pars.CurEvNo > 0){
  //if(((DFMAEvent[0].LEDts - t_first)/1E8)>30){
  if(Pars.CurEvNo > 0){


  // printf("\nTest1, event number %i\n",Pars.CurEvNo);

  float cal_e;

  char strg[128];
  int i, j, ii, jj;
  int ndssd;
  int ndfma;
  int nfp;
  int nsubev;
  int trn[100];

  for(i=0;i<100;i++){
    trn[i] = 0;
  }

  int GebTypeStr (int type, char strg[]);


  if (Pars.CurEvNo <= Pars.NumToPrint){
    printf ("entered bin_dfma:\n");
  }

  //********************************************************************\\
  // loop through the coincidence event and fish out GEB_TYPE_DFMA data \\
  //********************************************************************\\

  ndfma = 0;
  ndssd = 0;
  nsubev = 0;
  nfp = 0;

  for (i = 0; i < GEB_event->mult; i++){
    
    GebTypeStr (GEB_event->ptgd[i]->type, strg);

    if(GEB_event->ptgd[i]->type == 16){

      if (Pars.CurEvNo <= Pars.NumToPrint){
            printf ("bin_dfma, %2i> %2i, %s, TS=%lli\n", i, GEB_event->ptgd[i]->type, strg,GEB_event->ptgd[i]->timestamp);
      }

      DFMAEvDecompose_v3 ((unsigned int *) GEB_event->ptinp[i], GEB_event->ptgd[i]->length / sizeof (unsigned int),
                             &DFMAEvent[nsubev]);

        
        h1_wheel_total->Fill(DFMAEvent[nsubev].wheel);

      if(DFMAEvent[nsubev].tpe == DSSD){
             ndssd++;
             ndfma++;
      }
      if(DFMAEvent[nsubev].tpe == FP){
             nfp++;
             ndfma++;
      }
           nsubev++;
     }

           //nsubev++;

   }



h1_ndssd->Fill(ndssd);

if(DFMAEvent[0].LEDts != 0) numDFMA++;
if(DGSEvent[0].event_timestamp != 0) numDGS++;


if (first && (numDFMA > 100)) { t_first=DFMAEvent[0].LEDts; first=false; }
if(Pars.CurEvNo%10000 == 0) printf("Processing event number %i with timestamp %llu. Time elapsed: %i seconds\n",Pars.CurEvNo,DFMAEvent[0].LEDts,int(float(DFMAEvent[0].LEDts - t_first)/1E8));if(firstge) { t_firstge = DGSEvent[0].event_timestamp; firstge = false; }

long long int tmp1, tmp2;

// Declarations for variables

float cr_int;
cr_int = 0;

float csum, csumy;
csum=0; csumy=0;

float crat_int, cratn_int;
crat_int = 0.0;
cratn_int = 0.0;

//>>>>>>>>>>>>>>>>>>>>>>>>>

 float cl, cl_pu, cr, cr_pu, esi, emonl, emonr, wheel, wheel2, wheel_ge, wheel_left, wheel_up, ppacde, ppacde_pu, crat, cratn;
 float cup, cdown, craty;

float dssd_fr_emax;
float dssd_ba_emax;
dssd_fr_emax = 0.0;
dssd_ba_emax = 0.0;
int dssd_fr_subev;
int dssd_ba_subev;
dssd_fr_subev = 100000;
dssd_ba_subev = 100000;

cl=0;
cr=0;
cl_pu=0;
cr_pu=0;
esi=0;
emonl=0; emonr=0;
ppacde=0;
ppacde_pu=0;
crat=0;
cratn=0;

int gate = 0;

//<><><><><><><><><><><><><><><><><>\\
// Dig out DSSD event and calibrate \\
//<><><><><><><><><><><><><><><><><>\\

for (i=0;i<nsubev;i++) {

   if(DFMAEvent[i].tpe == DSSD){
    
        cal_e = .0;
        //cal_e = DFMAEvent[i].ehi;

  	h2_dssd_en_raw->Fill(DFMAEvent[i].ehi,DFMAEvent[i].tid);

   	if (DFMAEvent[i].ehi>50) h2_dssd_rate->Fill((DFMAEvent[i].LEDts - t_first)/ratecomp, DFMAEvent[i].tid);

        // FRONT - selecting maximum energy strip and calibrating

	if(DFMAEvent[i].tid < BATIDLOW){

	  // CORRECT FOR TE TAIL FROM THE PREVIOUS EVENT

          // GE correction

	  //bl_corr=(float)(DGSEvent[i].post_rise_energy*DGSEvent[i].m1_begin_sample-DGSEvent[i].pre_rise_energy*DGSEvent[i].m2_begin_sample)/
	  //(DGSEvent[i].post_rise_energy-DGSEvent[i].pre_rise_energy+M*(DGSEvent[i].m1_begin_sample-DGSEvent[i].m2_begin_sample));

          // M=250?
          // ehiPZ?

          //Energy = ((float)(DGSEvent[i].post_rise_energy)-(float)(DGSEvent[i].pre_rise_energy)*ehiPZ[gsid])/M*ehigain[gsid];
          //Energy = Energy - bl_corr_av2[gsid]*(1.-ehiPZ[gsid])*ehigain[gsid] + ehioffset[gsid]; 

           cal_e = double(map_fr[DFMAEvent[i].tid].gain)*(DFMAEvent[i].ehi + double(rand())/RAND_MAX-0.5) + double(map_fr[DFMAEvent[i].tid].off);  
	   
           if(cal_e > dssd_fr_emax) {
              dssd_fr_emax = cal_e;
              dssd_fr_subev = i;
           }
        }        

        // BACK - selecting maximum energy strip and calibrating

	if(DFMAEvent[i].tid > FRTIDHIGH){
           cal_e = double(map_ba[DFMAEvent[i].tid].gain)*(DFMAEvent[i].ehi + double(rand())/RAND_MAX-0.5) + double(map_ba[DFMAEvent[i].tid].off);  

           if(cal_e > dssd_ba_emax){
              dssd_ba_emax = cal_e;
              dssd_ba_subev = i;
           }
        }
        // UNCOMMENT THIS IF YOU WANT CALIBRATED ENERGIES PROPAGATED THROUGH REMAINDER OF CODE
        DFMAEvent[i].ehi = cal_e;
  	h2_dssd_en->Fill(DFMAEvent[i].ehi,DFMAEvent[i].tid);
   }
}

 h2_dssd_fb->Fill(dssd_fr_emax, dssd_ba_emax);

//<><><><><><><><><><><><><><><><><><><><><><><\\
// Time difference between DSSD front and back \\
//<><><><><><><><><><><><><><><><><><><><><><><\\

signed long long int tdssd_fr;
signed long long int tdssd_ba;

tdssd_fr = 0;
tdssd_ba = 0;

double frba_t;
frba_t = 0.0;


if(dssd_fr_emax != 0 && dssd_ba_emax != 0){

   frba_t = double(DFMAEvent[dssd_fr_subev].LEDts) - double(DFMAEvent[dssd_ba_subev].LEDts);
   h2_frba_t->Fill(frba_t,DFMAEvent[dssd_fr_subev].tid);
   h2_frba_t->Fill(frba_t,DFMAEvent[dssd_ba_subev].tid);

   tdssd_fr = DFMAEvent[dssd_fr_subev].LEDts;
   tdssd_ba = DFMAEvent[dssd_ba_subev].LEDts;
}
 
//<><><><><><><><><><><><><>\\
// Dig out focal plane data \\
//<><><><><><><><><><><><><>\\

signed long long int tdssdmcp_fr;
signed long long int tdssdmcp_ba;
signed long long int tdssdmcp_fr_r;
signed long long int tdssdmcp_ba_r;
signed long long int tdssdppacde_fr;

tdssdmcp_fr = 0;
tdssdmcp_ba = 0;
tdssdmcp_fr_r = 0;
tdssdmcp_ba_r = 0;

tdssdppacde_fr=0;

int left_subev = 0;
int right_subev = 0;
int up_subev = 0;
int down_subev = 0;

int letrace[MCPTRLEN+1];
int let0 = 0;
int let1 = 0;

#if(1)

 int n_sib;
 float sib_e[56];
 int sib_tid[56];
 int sib_wn[56];
 int sib_detn[56];
 int sib_stripn[56];
 
 int sib_wn_table[56]={4,4,4,4,4,4,4,4,4,4,4,4,4,4,
                 3,3,3,3,3,3,3,3,3,3,3,3,3,3,
                 2,2,2,2,2,2,2,2,2,2,2,2,2,2,
                 1,1,1,1,1,1,1,1,1,1,1,1,1,1};
 int sib_detn_table[56]={1,1,1,1,1,1,1,2,2,2,2,2,2,2,
                   1,1,1,1,1,1,1,2,2,2,2,2,2,2,
                   2,2,2,2,2,2,2,1,1,1,1,1,1,1,
		   2,2,2,2,2,2,2,1,1,1,1,1,1,1};
 int sib_stripn_table[56]={1,2,3,4,5,6,7,1,2,3,4,5,6,7,
                     1,2,3,4,5,6,7,1,2,3,4,5,6,7,
                     1,2,3,4,5,6,7,1,2,3,4,5,6,7,
                     1,2,3,4,5,6,7,1,2,3,4,5,6,7};

 float sib_emax;
 signed long long int sib_ts[56];
 int idnotpresent;

 n_sib=0;
 sib_emax=0;

for (i=0;i<nsubev;i++) {
  
  switch (DFMAEvent[i].tpe) {
    
         //-------------\\
         // Focal Plane \\
         //-------------\\


	 case FP:

        if(DFMAEvent[i].tid == WHEELCH && DFMAEvent[i].ehi > 0) 
	  { wheel = DFMAEvent[i].ehi; h1_wheel->Fill(wheel); ts_wheel=DFMAEvent[i].LEDts; }

        if(DFMAEvent[i].tid == WHEEL2CH && DFMAEvent[i].ehi > 0) 
	  { wheel2 = DFMAEvent[i].ehi; h1_wheel2->Fill(wheel2); ts_wheel2=DFMAEvent[i].LEDts; }

	 if(DFMAEvent[i].tid == SICH && DFMAEvent[i].ehi > 0) 
           { esi = DFMAEvent[i].ehi; h1_esi->Fill(esi); ts_esi=DFMAEvent[i].LEDts; 
             h1_esidt->Fill(ts_esi/10000-ts_wheel/10000);
             h1_esidt2->Fill(ts_esi/10000-ts_wheel2/10000);
             h2_esidtr->Fill(DFMAEvent[i].LEDts/10000000000-t_first/10000000000,DFMAEvent[i].LEDts/10000-ts_wheel/10000);
             if ((esi>3000)&&(esi<6000)) h2_esidtrg->Fill(DFMAEvent[i].LEDts/10000000000-t_first/10000000000,DFMAEvent[i].LEDts/10000-ts_wheel/10000);
           }

	 if(DFMAEvent[i].tid == EMONLCH && DFMAEvent[i].ehi > 0) 
           { emonl = DFMAEvent[i].ehi; h1_emonl->Fill(emonl); ts_emonl=DFMAEvent[i].LEDts; 
             if (emonl>20000 && emonl<50000) {
              h1_emonl_rate->Fill((ts_emonl-t_first)/1e8);
              h1_emonl_rate2->Fill((ts_emonl-t_first)/1e7);
              h1_wheel_mon->Fill(DFMAEvent[i].wheel);
             }
           }

	 if(DFMAEvent[i].tid == EMONRCH && DFMAEvent[i].ehi > 0) 
           { emonr = DFMAEvent[i].ehi; h1_emonr->Fill(emonr); ts_emonr=DFMAEvent[i].LEDts; 
             h1_emonr_rate->Fill((ts_emonr-t_first)/1e8);
              h1_emonr_rate2->Fill((ts_emonr-t_first)/1e7);
              h1_wheel_mon->Fill(DFMAEvent[i].wheel);
           }

	 // MWPC LEFT
  	 if (DFMAEvent[i].tid==LCH && DFMAEvent[i].ehi > 0) {

                left_subev = i;
               
                 cl_pu = DFMAEvent[i].pu;

         	h2_FP_rate->Fill((DFMAEvent[i].LEDts - t_first)/ratecomp, DFMAEvent[i].tid);
                if((DFMAEvent[i].LEDts - t_first)/ratecomp < 100000000){
         	    h1_FP_rate_left->Fill((DFMAEvent[i].LEDts - t_first)/ratecomp);
		    }
                tdssdmcp_fr = (DFMAEvent[i].LEDts) -  (tdssd_fr);
                tdssdmcp_ba = (DFMAEvent[i].LEDts) -  (tdssd_ba);

                if(tdssd_fr != 0) h2_tdssdmcp_fr->Fill(tdssdmcp_fr,DFMAEvent[dssd_fr_subev].tid);
                if(tdssd_ba != 0) h2_tdssdmcp_ba->Fill(tdssdmcp_ba,DFMAEvent[dssd_ba_subev].tid-160);


    		cl=DFMAEvent[i].ehi;
 
                h1_leftdt->Fill(DFMAEvent[i].LEDts/10000-ts_wheel/10000);
                h1_leftdt2->Fill(DFMAEvent[i].LEDts/10000-ts_wheel2/10000);
                h2_leftdtr->Fill(DFMAEvent[i].LEDts/10000000000-t_first/10000000000,DFMAEvent[i].LEDts/10000-ts_wheel/10000);

                ts_left=DFMAEvent[i].LEDts;

                wheel_left=DFMAEvent[i].wheel;
                h1_wheel_left->Fill(DFMAEvent[i].wheel);

	 }

	 // MWPC RIGHT
	 if (DFMAEvent[i].tid==RCH && DFMAEvent[i].ehi > 0){

                right_subev = i;
                cr_pu = DFMAEvent[i].pu;
             	h2_FP_rate->Fill((DFMAEvent[i].LEDts - t_first)/ratecomp, DFMAEvent[i].tid);
                if((DFMAEvent[i].LEDts - t_first)/ratecomp < 100000){
         	    h1_FP_rate_right->Fill((DFMAEvent[i].LEDts - t_first)/ratecomp);
                }
                tdssdmcp_fr_r = (DFMAEvent[i].LEDts) -  (tdssd_fr);
                tdssdmcp_ba_r = (DFMAEvent[i].LEDts) -  (tdssd_ba);
                if(tdssd_fr != 0) h2_tdssdmcp_fr_r->Fill(tdssdmcp_fr_r,DFMAEvent[dssd_fr_subev].tid);
                if(tdssd_ba != 0) h2_tdssdmcp_ba_r->Fill(tdssdmcp_ba_r,DFMAEvent[dssd_ba_subev].tid-160);
  	    	cr=DFMAEvent[i].ehi;

	}

	 // MWPC UP
	 if (DFMAEvent[i].tid==UPCH && DFMAEvent[i].ehi > 0){

                up_subev = i;
  	    	cup=DFMAEvent[i].ehi;
		h1_cup->Fill(cup);
               if((DFMAEvent[i].LEDts - t_first)/ratecomp < 100000){
                  if ((cup>100)&&(cup<100000)) h1_FP_rate_up->Fill((DFMAEvent[i].LEDts - t_first)/ratecomp);
                }
	}

	 // MWPC UP
	 if (DFMAEvent[i].tid==DOWNCH && DFMAEvent[i].ehi > 0){

                down_subev = i;
  	    	cdown=DFMAEvent[i].ehi;
		h1_cdown->Fill(cdown);

	}

	 // MWPC DE
	 if (DFMAEvent[i].tid==PPACDE && DFMAEvent[i].ehi > 0) {

           ppacde =  DFMAEvent[i].ehi;
           tdssdppacde_fr = (DFMAEvent[i].LEDts) -  (tdssd_fr);     
           h2_tdssdppacde_fr->Fill(tdssdppacde_fr,DFMAEvent[dssd_fr_subev].tid);
           h2_tdssdppacde->Fill(tdssdppacde_fr,ppacde);

	 }

        break;

        default:
	break;
     
  } 
  
}

#endif


//<><><>\\
// MWPC \\
//<><><>\\

  h1_esi->Fill(esi);

  h1_ppacde->Fill(ppacde);

  h1_cl_pu->Fill(cl_pu);
  h1_cr_pu->Fill(cr_pu);

  h2_clr->Fill(cl,cr);

  if(cl > 0 && cr > 0) { 

    csum = cl + cr;
    crat = cl - cr;
 
    h1_cx->Fill(crat);
    h2_xde->Fill(crat,ppacde);

    h1_cl->Fill(cl);
    h1_cr->Fill(cr);

    h1_csum->Fill(csum);

    if(cup > 0 || cdown > 0) { 
     craty=cup-cdown;
     csumy=cup+cdown;
     h1_cy->Fill(craty);
     h2_cud->Fill(cup,cdown);
     h1_csumy->Fill(csumy);

  }

  if ((cup > 0) && (cdown > 0) && (cl > 0) && (cr > 0)) { 
      h2_cxy->Fill(crat,craty);
  }
 
 } // cl and  cr not zero


//<><><><><><><><><>\\
// Print statements \\
//<><><><><><><><><>\\

if (Pars.CurEvNo <= Pars.NumToPrint){

printf("Print statements at end of bin_dfma\n");

 for(i=0;i<1;i++){

   
      printf ("\n\n\nwe have %i gamma rays\n", ng);
      printf ("%2i> ", i);
      printf ("chan_id=%i; ", DGSEvent[i].chan_id);
      printf ("board_id=%i; ", DGSEvent[i].board_id);
      printf ("id =%i; ", DGSEvent[i].id);
      printf ("tpe=%i; ", DGSEvent[i].tpe);
      printf ("tid=%i; ", DGSEvent[i].tid);
      printf ("EventTS=%llu; ", DGSEvent[i].event_timestamp);
      printf ("ehi=%i ", DGSEvent[i].ehi);
      printf ("\n\n\n");

   
 }
}

  // debug list the dssd events we found 

  if (Pars.CurEvNo <= Pars.NumToPrint)
    for (i = 0; i < ndssd; i++)
      {
        printf ("we have %i DSSD event(s)\n", ndssd);
        printf ("%2i> ", i);
        printf ("chan_id=%i; ", DFMAEvent[i].chan_id);
        printf ("board_id=%i; ", DFMAEvent[i].board_id);
        printf ("id =%i; ", DFMAEvent[i].id);
        printf ("tpe=%i; ", DFMAEvent[i].tpe);
        printf ("tid=%i; ", DFMAEvent[i].tid);
        printf ("LEDTS=%llu; ", DFMAEvent[i].LEDts);
        printf ("ehi=%8i ", DFMAEvent[i].ehi);
        printf ("\n\n\n\n");
        fflush (stdout);
      };


  // done


  if (Pars.CurEvNo <= Pars.NumToPrint)
    printf ("exit bin_dfma\n");

  return (0);
  }

 

#endif
}

//*******************************************************************************************************

#if(0)
void
SetBeta ()
{

  /* declarations */

  int i;
  double d1;

  /*-------------------------------------*/
  /* find Doppler and aberration factors */
  /*-------------------------------------*/

  for (i = 0; i < NGSGE; i++)
    {
      //printf("det %3.3i-> ", i);
      d1 = angle[i] / 57.29577951;
      DopCorFac[i] = (1 - Pars.beta * cos (d1)) / sqrt (1 - Pars.beta * Pars.beta);
      //printf("dop cor fac: %6.4f; ", DopCorFac[i]);
      ACFac[i] = DopCorFac[i] * DopCorFac[i];
      //printf("aberration cor fac: %6.4f\n", ACFac[i]);

    };
  fflush (stdout);



}

#endif
