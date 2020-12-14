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
#include "functions_dfma.h"
#include "get_dead_layer_corrections2.cpp"

//#define LEFT 1
//#define RIGHT 2
//#define ICDE 4
//#define PPACDE2 5

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
   //unsigned short int traceLen;
   //int dssdbase;
   //short int trace[MAXTRACELEN];
   int traceLen;
   //short int trace[MAXTRACELEN];
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



//>>>>>>>>>>>>>>>>>>\\
// RING BUFFER CODE \\

#define RBUFFSIZE 5

#if(0)


struct dssd_type {

   decay_type_new east;
   decay_type_new west;
   decay_type_new top;

};


struct dssd_buffer_type {

   dssd_type event[RBUFFSIZE];
   int head;
   int tail;
   int count;

};

struct dssd_type dssd_event;
struct dssd_buffer_type dssd_buffer;

//struct global_type global;

struct dssd_buffer_type dssd_tree;



struct decay_type_new east;



struct decay_type_new empty;


struct recoil_buffer_type {

  recoil_type recoil[RBUFFSIZE];
  int head;
  int tail;
  int count;
  int s_fr;
  int s_ba;

}; 

struct recoil_buffer_type recbuffer;

struct decay_buffer_type {

  decay_type decay[RBUFFSIZE];
  int head;
  int tail;
  int count;
  int s_fr;
  int s_ba;
}; 

struct decay_buffer_type decbuffer;


struct pixel_buffer_type {
   recoil_buffer_type recbuffer;
   decay_buffer_type decbuffer;
   int status_flag;
};

struct pixel_buffer_type pixel_buff[101][101];

struct pixel_buffer_type front_buff[101];
struct pixel_buffer_type back_buff[101];

//*********************************************************\\


void DSSDBuffer_put (dssd_buffer_type *_this, dssd_type &c)
{

    dssd_type newrecoil = c;

    _this->head = ((_this->head)%(RBUFFSIZE)+1);
    _this->event[(_this->head)-1] = newrecoil;	// WAS _this->recoil[(_this->head)-1] = newrecoil;

    ++_this->count;
   
}

//*********************************************************\\

//*********************************************************\\

void RBuffer_init (recoil_buffer_type *_this)
{
         //The following clears:
        //-> buf
       // -> head
      //  -> tail
     //  -> count
    //  and sets head = tail
    
    memset (_this, 0, sizeof (*_this));

}
//*********************************************************\\

void DBuffer_init (decay_buffer_type *_this)
{
         //The following clears:
        //-> buf
       // -> head
      //  -> tail
     //  -> count
    //  and sets head = tail
    
    memset (_this, 0, sizeof (*_this));

}

//*********************************************************\\

int RBuffer_empty (recoil_buffer_type *_this)
{
    // If the buffer is empty, returns 0. If not, returns 1.
    return (0==_this->count);
}


//*********************************************************\\

int RBuffer_full (recoil_buffer_type *_this)
{    
    // If the buffer is NOT full, returns 0. If FULL, returns 1.
    return (_this->count>=RBUFFSIZE);
}

//*********************************************************\\

void RBuffer_put (recoil_buffer_type *_this, recoil_type &c)
{

    recoil_type newrecoil = c;

    _this->head = ((_this->head)%(RBUFFSIZE)+1);
    _this->recoil[(_this->head)-1] = newrecoil;	// WAS _this->recoil[(_this->head)-1] = newrecoil;

    ++_this->count;
   
}

//*********************************************************\\

void DBuffer_put (decay_buffer_type *_this, decay_type &c)
{

    decay_type newrecoil = c;

    _this->head = ((_this->head)%(RBUFFSIZE)+1);
    _this->decay[(_this->head)-1] = newrecoil; // WAS_this->decay[(_this->head)-1] = newrecoil;

    ++_this->count;
   
}

//*********************************************************\\

int lookback (int P, int M){

   int f = 0;
   
   f = ((P+RBUFFSIZE-(M%RBUFFSIZE))%RBUFFSIZE);
   

   return f;

}

#endif

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

// DEFINE HISTOGRAMS HERE!!!! (LIKE USERDECLARE.H)

TH1D *h1_short_dec_en;
TH1D *h1_short_rec_en;
TH1D *h1_short_rec_bkg_en;
TH2F *h2_le_output_left;
TH2F *h2_rec_gammas;
TH2F *h2_all_gammas;
TH2F *h2_dssd_gammas;
TH2F *h2_mcp_gammas;

TH1D *h1_nge;

//TH2F *h2_te107_gammas;

TH2F *h2_corr_gammas, *h2_corr_gammas_box, *h2_corr_gammasf, *h2_corr_gammasff, *h2_corr_gammasfff, *h2_corr_gammasffff; 

TH2F *h2_corr_cl, *h2_corr_cla, *h2_corr_alpha_cl, *h2_corr_alpha_cla, *h2_corr_proton_cl;
TH2F *h2_corr_gammas_id;
TH2F *h2_corr_gammas_short;
TH2F *h2_gamma_gamma_shortdec;
TH2F *h2_dTgdssd, *h2_dTgdssd_corr, *h2_dTgdssd_corr_ehi, *h2_dTgdssd_ch40, *h2_dTcldssd_corr, *h2_dTcldssda_corr, *h2_dTcldssdall, *h2_dTcldssdp_corr, *h2_dTgdssdr_clenergy, *h2_corr_proton_cl_times, *h2_dTXadssdE, *h2_dTXaclE, *h2_dssdEclE;
TH2F *h2_dTgdssd_allge;
TH1D *h1_dTgl, *h1_dTgr;
TH2F *h2_stripfr_cxng, *h2_stripba_cxng, *h2_stripba_cxng_frontgated;

TH1D *h1_cl, *h1_cr, *h1_cup, *h1_cdown, *h1_esi, *h1_emonl, *h1_emonr, *h1_wheel_total, *h1_wheel_mon, *h1_wheel_left,  *h1_wheel_up, *h1_wheel_ge, *h1_wheel, *h1_wheel2, *h1_ppacde, *h1_ppacdeg, *h1_cx, *h1_cy, *h1_cxn, *h1_csum, *h1_csumy;
TH1D *h1_emonl_rate, *h1_emonl_rate2;
TH1D *h1_emonr_rate, *h1_emonr_rate2;
TH1D *h1_cl_pu, *h1_cr_pu, *h1_ppacde_pu;

//TH2F *h2_ppacde_E, *h2_ppacde_E_corr;

TH2F *h2_clr, *h2_cud, *h2_cxy;
TH1D *h1_clg, *h1_crg, *h1_cxg, *h1_cxng;
TH1D *h1_cx_en_g1, *h1_cxn_en_g1;
TH1D *h1_cx_en_g2, *h1_cxn_en_g2;

TH2F *h2_xde, *h2_xdeg; 


TH2F *h2_clrg;
TH2F *h2_clrg_int; //new
TH2F *h2_traceg;
TH2F *h2_tracegb14;
TH2F *h2_tracegb15;
TH2F *h2_traceg2;
TH2F *h2_traceg3;
TH1D *h1_trgn, *h1_trg2n, *h1_trg3n;

TH1D *h1_ch0, *h1_ch1, *h1_ch2, *h1_ch3, *h1_ch4, *h1_ch5;

TH2F *h2_frba_t;


// DSSD stuff

TH2F *h2_dssd_en, *h2_dssd_en_raw;
TH2F *h2_dssd_fb;
TH2F *h2_dssd_fr_emax, *h2_dssd_fr_emax_phy;
TH2F *h2_dssd_ba_emax, *h2_dssd_ba_emax_phy;
TH2F *h2_d_fr_e, *h2_d_ba_e;
TH2F *h2_r_fr_e, *h2_r_ba_e;

/*POINT*/

TH1D *h1_basesample;
TH2F *h2_blav;
TH2F *h2_basesample;
TH2F *h2_baseline;
TH2F *h2_prerisebeg;
TH2F *h2_preriseend;
TH2F *h2_postrisebeg;
TH2F *h2_postriseend;

TH2F *h2_pre2;
TH2F *h2_post2;
TH2F *h2_last2;

TH2F *h2_e_bs;
TH2F *h2_e_dt;

TH1D *h1_header_type;

TH2F *h2_dssd_rate, *h2_FP_rate;
//TH1D *h1_dssd_rate_1us;
TH1D *h1_dssd_rate_1D;
TH1D *h1_decay_rate, *h1_recoil_rate, *h1_corr_rate;
TH1D *h1_GE_rate;
TH1D *h1_FP_rate_left, *h1_FP_rate_right, *h1_FP_rate_up;
TH2F *h2_tdssdmcp_fr, *h2_tdssdmcp_ba;
TH2F *h2_tdssdmcp_fr_r, *h2_tdssdmcp_ba_r, *h2_tdssdppacde_fr ;
TH2F  *h2_tdssdppacde ;
TH1D *h1_ndssd;

TH2F *h2_lefttraces;
TH2F *h2_righttraces;

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

// Correlation histograms

TH2F *h2_e1t1log, *h2_e1t1log_box, *h2_e1t1log_box1, *h2_e2t2log, *h2_e3t3log,*h2_e4t4log,*h2_e5t5log, *h2_e6t6log;
TH2F *h2_e1t1logf, *h2_e2t2logf, *h2_e3t3logf,*h2_e4t4logf,*h2_e5t5logf, *h2_e6t6logf;

TH2F *h2_e1t1c1, *h2_e2t2c1, *h2_e3t3c1;
TH2F *h2_e1t1c2, *h2_e2t2c2, *h2_e3t3c2;
TH2F *h2_e1t1c3, *h2_e2t2c3, *h2_e3t3c3;
TH2F *h2_e1t1c4, *h2_e2t2c4, *h2_e3t3c4;


float corr_time, corr_time_short;

TH2F *h2_e1e2, *h2_e1e2_short, *h2_e2e3, *h2_e3e4, *h2_e4e5, *h2_e5e6,  *h2_e2e3f, *h2_e3e4f, *h2_e4e5f, *h2_e5e6f;

TH2F *h2_eftof, *h2_eftofg, *h2_eftof_corr;


// clovers

TH2F *h2_cle_id, *h2_cle_id_raw, *h2_IDT;

/*
TH2F *h2_cle_id_betag, *h2_cle_id_keV;
TH2F *h2_cl_ehi_basesample[21];
TH2F *h2_cl_ehi_slope[21];
*/
TH1D *h1_clover_en;
/*
TH1D *h1_dec_clover;
*/



TH2F *h2_cl_E4_ehi_baseline;
TH2F *h2_cl_E4_ehi_prerisebeg;
TH2F *h2_cl_E4_ehi_preriseend;
TH2F *h2_cl_E4_preriseend_prerisebeg;
TH2F *h2_cl_E4_preriseend_diff1;
TH2F *h2_cl_E4_preriseend_diff2;
TH2F *h2_cl_E4_ehi_sampledbase;
TH2F *h2_cl_E4_ehi_diff1;
TH2F *h2_cl_E4_ehi_diff2;
TH2F *h2_cl_E4_ehi_diff3;
TH2F *h2_cl_E4_ehi_diff4;
TH2F *h2_cl_E4_ehi_trueBL;

TH2F *h2_electron_traces_ch100_fr;
TH2F *h2_alpha_traces_ch100_fr;
TH2F *h2_electron_traces_ch100_ba;

// Si Box

TH2F *h2_sibox_raw;
TH2F *h2_sibox;

TH2F *h2_sisi, *h2_sidssd, *h2_dtdssdfrsi,  *h2_dtdssdbasi;
TH2F *h2_sisi_corr, *h2_sidssd_corr,  *h2_dtdssdfrsi_corr,  *h2_dtdssdbasi_corr;

TH2F *h2_sidssdcor, *h2_sidssdcor2, *h2_sisicor2esc1, *h2_sisicor2esc2, *h2_aa1esc, *h2_aa2esc;  

TH2F *h2_ada1,  *h2_ada2, *h2_ada3,  *h2_ada4, *h2_ada5,  *h2_ada6, *h2_ada7,  *h2_ada8 ;

TH1D *h1_atot, *h1_atotcor;

TH1D *h1_fdecEv;

// wheel

long long int ts_wheel=0;
long long int ts_wheel2=0;
long long int ts_esi=0;
long long int ts_emonl=0;
long long int ts_emonr=0;

long long int ts_left;

TH1D* h1_esidt, *h1_esidt2;
TH1D* h1_leftdt, *h1_leftdt2;
TH2F* h2_leftdtr, *h2_leftdtam, *h2_esidtr, *h2_esidtrg;

TH2F *h2_ts2, *h2_ts2c;

 TCutG *eftof;

extern FILE *trf;

/*-----------------------------------------------------*/


int exit_dfma()
{

int i,j;

for (i=1;i<161;i++) {
for (j=1;j<161;j++) {

  if (dssd_corr[i][j].status>1) {
       chain=dssd_corr[i][j].chain;     
        tree->Fill();
  }

 }
}

#if(1)
treef->cd();

tree->Write("tree");
Pars.f1->cd();
treef->Close();
#endif
//if (trf==NULL) (printf("TRACE FILE DOES NOT CLOSE\n")); 
//fclose(trf); 
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

// //TREES...
strcpy(compstr,"2");
strcpy(presortFile,"TREE_FILES/new_presort.");
//strcpy(presortFile,"/home/dgs/data2/user/254No/TREE_FILES/new_presort.");

   i=strlen(Pars.ROOTFile);

  while (Pars.ROOTFile[i]!='/') {
    i--;
  }

for (j=i+1;j<strlen(Pars.ROOTFile);j++) {
  onechar[0]=Pars.ROOTFile[j];
  onechar[1]='\0';
  strcat(presortFile,onechar);
 }

printf("ROOT file %s\n",Pars.ROOTFile);
printf("presort file %s\n",presortFile);



#if(1)
treef = new TFile(presortFile,"recreate");
tree = new TTree("tree","presorted data");

Pars.f1->cd();

tree->Branch("s_fr",&chain.s_fr,"s_fr/I");
tree->Branch("s_ba",&chain.s_ba,"s_ba/I");

//tree->Branch("recoil",&chain.recoil,"ts/l:en/i:pu/i:left/F:right/F:x/D:d2t0/i:d2t1/i:d2t2/i:nge/i:geehi[15]/L:tge[15]/l:getid[15]/i:traceLen/i:trace_fr[256]/S:trace_ba[256]/S");

 tree->Branch("recoil",&chain.recoil,"ts/l:en/I:enb/I:pu/I:left/F:right/F:ptofl/L:ptofr/L:pdetof/L:pde/F:x/D:d2t0/I:d2t1/I:d2t2/I:nge/I:geehi[15]/L:tge[15]/l:getid[15]/I:traceLen/I:trace_fr[256]/S:trace_ba[256]/S");   //Adding the back side energy, PPACtof, PPACdetof and PPACde by Tianheng.

//tree->Branch("recoil",&chain.recoil,"ts/l:en/i:toten/i:pu:right:left:ppacde:ppacde2:icde:nge:geehi[15]/L:tge[15]/l:getid[15]/i:tof/l:d2t0/i:d2t1/i:fdecEv/i:traceLen/s:trace_fr[DSSDTRLEN]/S:trace_ba[DSSDTRLEN]/S");
tree->Branch("ndec",&chain.ndec,"ndec/I");

tree->Branch("dec1",&chain.decay[0],"ts/l:en/I:enb/I:pu_fr/I:pu_ba/I:time/l:traceLen/I:trace_fr[256]/S:trace_ba[256]/S:nge/I:geehi[15]/L:tge[15]/l:getid[15]/I:nsi/I:siehi[10]/L:tsi[10]/l:sitid[10]/I:d2t0/I:d2t1/I:d2t2/I");
tree->Branch("dec2",&chain.decay[1],"ts/l:en/I:enb/I:pu_fr/I:pu_ba/I:time/l:traceLen/I:trace_fr[256]/S:trace_ba[256]/S:nge/I:geehi[15]/L:tge[15]/l:getid[15]/I:nsi/I:siehi[10]/L:tsi[10]/l:sitid[10]/I:d2t0/I:d2t1/I:d2t2/I");
tree->Branch("dec3",&chain.decay[2],"ts/l:en/I:enb/I:pu_fr/I:pu_ba/I:time/l:traceLen/I:trace_fr[256]/S:trace_ba[256]/S:nge/I:geehi[15]/L:tge[15]/l:getid[15]/I:nsi/I:siehi[10]/L:tsi[10]/l:sitid[10]/I:d2t0/I:d2t1/I:d2t2/I");
tree->Branch("dec4",&chain.decay[3],"ts/l:en/I:enb/I:pu_fr/I:pu_ba/I:time/l:traceLen/I:trace_fr[256]/S:trace_ba[256]/S:nge/I:geehi[15]/L:tge[15]/l:getid[15]/I:nsi/I:siehi[10]/L:tsi[10]/l:sitid[10]/I:d2t0/I:d2t1/I:d2t2/I");
tree->Branch("dec5",&chain.decay[4],"ts/l:en/I:enb/I:pu_fr/I:pu_ba/I:time/l:traceLen/I:trace_fr[256]/S:trace_ba[256]/S:nge/I:geehi[15]/L:tge[15]/l:getid[15]/I:nsi/I:siehi[10]/L:tsi[10]/l:sitid[10]/I:d2t0/I:d2t1/I:d2t2/I");
tree->Branch("dec6",&chain.decay[5],"ts/l:en/I:enb/I:pu_fr/I:pu_ba/I:time/l:traceLen/I:trace_fr[256]/S:trace_ba[256]/S:nge/I:geehi[15]/L:tge[15]/l:getid[15]/I:nsi/I:siehi[10]/L:tsi[10]/l:sitid[10]/I:d2t0/I:d2t1/I:d2t2/I");



//tree->Branch("dec1",&chain.decay[0],"ts/l:en/i:pu_fr/i:pu_ba/i:time/l:traceLen/i:trace_fr[1000]/S:trace_ba[1000]/S:nge/i:geehi[15]/L:tge[15]/l:getid[15]/i");
//tree->Branch("dec2",&chain.decay[1],"ts/l:en/i:pu_fr/i:pu_ba/i:time/l:traceLen/i:trace_fr[1000]/S:trace_ba[1000]/S:nge/i:geehi[15]/L:tge[15]/l:getid[15]/i");
//tree->Branch("dec3",&chain.decay[2],"ts/l:en/i:pu_fr/i:pu_ba/i:time/l:traceLen/i:trace_fr[1000]/S:trace_ba[1000]/S:nge/i:geehi[15]/L:tge[15]/l:getid[15]/i");


#endif


////TREES...

//**************************\\
// Histogram initialisation \\
//**************************\\

// LIKE USERINIT.H!!!


h1_short_dec_en = mkTH1D((char *)"short_dec_en",(char *)"short_dec_en",5000,0,50000);
h1_short_rec_en = mkTH1D((char *)"short_rec_en",(char *)"short_rec_en",5000,0,50000);
h1_short_rec_bkg_en = mkTH1D((char *)"short_rec_bkg_en",(char *)"short_rec_bkg_en",5000,0,50000);

h2_le_output_left = mkTH2F((char *)"le_output_left",(char *)"le_output_left",100,0,100,1000,0,1000);

//h2_rec_gammas = mkTH2F((char *)"rec_gammas",(char *)"rec_gammas",10000,0,10000,120,0,120);
//h2_all_gammas = mkTH2F((char *)"all_gammas",(char *)"all_gammas",10000,0,10000,120,0,120);
//h2_dssd_gammas = mkTH2F((char *)"dssd_gammas",(char *)"dssd_gammas",10000,0,10000,120,0,120);

//h2_corr_gammas_id = mkTH2F((char *)"corr_gammas_id",(char *)"corr_gammas_id",10000,0,10000,120,0,120);
//h2_mcp_gammas = mkTH2F((char *)"mcp_gammas",(char *)"mcp_gammas",10000,0,10000,120,0,120);
//h2_corr_gammas = mkTH2F((char *)"corr_gammas",(char *)"corr_gammas",3000,0,30000,2000,0,2000);
//h2_corr_gammas_box = mkTH2F((char *)"corr_gammas_box",(char *)"corr_gammas_box",3000,0,300000,2000,0,2000);
//h2_corr_gammasf = mkTH2F((char *)"corr_gammasf",(char *)"corr_gammasf",3000,0,300000,2000,0,2000);
//h2_corr_gammasff = mkTH2F((char *)"corr_gammasff",(char *)"corr_gammasff",3000,0,300000,2000,0,2000);
//h2_corr_gammasfff = mkTH2F((char *)"corr_gammasfff",(char *)"corr_gammasfff",3000,0,300000,2000,0,2000);
//h2_corr_gammasffff = mkTH2F((char *)"corr_gammasffff",(char *)"corr_gammasffff",3000,0,300000,2000,0,2000);
//h2_corr_gammas_short = mkTH2F((char *)"corr_gammas_short",(char *)"corr_gammas_short",3000,0,300000,2000,0,2000);
h2_corr_cl = mkTH2F((char *)"corr_cl",(char *)"corr_cl",3000,0,30000,2000,0,2000);
h2_corr_cla = mkTH2F((char *)"corr_cla",(char *)"corr_cla",3000,0,30000,2000,0,2000);

h2_corr_alpha_cl = mkTH2F((char *)"corr_alpha_cl",(char *)"corr_alpha_cl",3000,0,30000,2000,0,2000);
h2_corr_alpha_cla = mkTH2F((char *)"corr_alpha_cla",(char *)"corr_alpha_cla",3000,0,30000,2000,0,2000);
h2_corr_proton_cl = mkTH2F((char *)"corr_proton_cl",(char *)"corr_proton_cl",1000,0,10000,2000,0,2000);
h2_corr_proton_cl_times = mkTH2F((char *)"corr_proton_cl_times",(char *)"corr_proton_cl_times",2000,0,2000,4000,-2000,2000);

h2_IDT = mkTH2F((char *)"IDT",(char *)"IDT",2000,0,2000,2000,0,2000);

//h2_gamma_gamma_shortdec = mkTH2F((char *)"gamma_gamma_shortdec",(char *)"gamma_gamma_shortdec",2000,0,2000,2000,0,2000);
h1_nge = mkTH1D((char *)"nge",(char *)"nge",500,0,500);

h2_dTgdssd = mkTH2F((char *)"dTgdssd",(char *)"dTgdssd",4000,-2000,2000,400,0,400);
h2_dTgdssd_corr = mkTH2F((char *)"dTgdssd_corr",(char *)"dTgdssd_corr",4000,-2000,2000,400,0,400);
h2_dTcldssd_corr = mkTH2F((char *)"dTcldssd_corr",(char *)"dTcldssd_corr",4000,-2000,2000,400,0,400);
h2_dTcldssda_corr = mkTH2F((char *)"dTcldssda_corr",(char *)"dTcldssda_corr",4000,-2000,2000,400,0,400);
h2_dTcldssdp_corr = mkTH2F((char *)"dTcldssdp_corr",(char *)"dTcldssdp_corr",4000,-2000,2000,400,0,400);
h2_dTgdssdr_clenergy= mkTH2F((char *)"dTgdssdr_clenergy",(char *)"dTgdssdr_clenergy",8000,-2000,6000,2000,0,2000);
h2_dTgdssd_corr_ehi = mkTH2F((char *)"dTgdssd_corr_ehi",(char *)"dTgdssd_corr_ehi",2000,-1000,1000,1000,0,1000);
h2_dTgdssd_ch40 = mkTH2F((char *)"dTgdssd_ch40",(char *)"dTgdssd_ch40",4000,-2000,2000,400,0,400);
h2_dTgdssd_allge = mkTH2F((char *)"dTgdssd_allge",(char *)"dTgdssd_allge",4000,-2000,2000,400,0,400);
h2_dTXadssdE= mkTH2F((char *)"dTXadssdE",(char *)"dTXadssdE",2000,0,10000,4000,-2000,2000);
h2_dTXaclE= mkTH2F((char *)"dTXaclE",(char *)"dTXaclE",2000,0,20000,4000,-2000,2000);
h2_dssdEclE= mkTH2F((char *)"dssdEclE",(char *)"dssdEclE",2000,0,10000,2000,0,2000);

h2_dTcldssdall = mkTH2F((char *)"dTcldssdall",(char *)"dTcldssdall",4000,-2000,2000,400,0,400);

h1_dTgl = mkTH1D((char *)"dTgl",(char *)"dTgl",4000,-2000,2000);
h1_dTgr = mkTH1D((char *)"dTgr",(char *)"dTgr",4000,-2000,2000);

h2_stripfr_cxng = mkTH2F((char *)"stripfr_cxng",(char *)"stripfr_cxng",200,0,200,4000,0,4000);
h2_stripba_cxng = mkTH2F((char *)"stripba_cxng",(char *)"stripba_cxng",200,0,200,4000,0,4000);
h2_stripba_cxng_frontgated = mkTH2F((char *)"stripba_cxng_frontgated",(char *)"stripba_cxng_frontgated",200,0,200,4000,0,4000);
h1_esi = mkTH1D((char *)"esi",(char *)"esi",4000,0,40000);
h1_emonl = mkTH1D((char *)"emonl",(char *)"emonl",40000,0,400000);
h1_emonr = mkTH1D((char *)"emonr",(char *)"emonr",40000,0,400000);
h1_wheel = mkTH1D((char *)"wheel",(char *)"wheel",4000,0,40000);
h1_wheel2 = mkTH1D((char *)"wheel2",(char *)"wheel2",4000,0,40000);
h1_wheel_total = mkTH1D((char *)"wheel_total",(char *)"wheel_total",4000,0,4000);
h1_wheel_mon = mkTH1D((char *)"wheel_mon",(char *)"wheel_mon",4000,0,4000);
h1_wheel_left = mkTH1D((char *)"wheel_left",(char *)"wheel_left",4000,0,4000);
h1_wheel_up = mkTH1D((char *)"wheel_up",(char *)"wheel_up",4000,0,4000);
h1_wheel_ge = mkTH1D((char *)"wheel_ge",(char *)"wheel_ge",4000,0,4000);
h1_ppacde = mkTH1D((char *)"ppacde",(char *)"ppacde",4000,0,40000);
h1_emonl_rate = mkTH1D((char *)"emonl_rate",(char *)"emonl_rate",100000,0,100000);
h1_emonl_rate2 = mkTH1D((char *)"emonl_rate2",(char *)"emonl_rate2",1000000,0,1000000);
h1_emonr_rate = mkTH1D((char *)"emonr_rate",(char *)"emonr_rate",100000,0,100000);
h1_emonr_rate2 = mkTH1D((char *)"emonr_rate2",(char *)"emonr_rate2",1000000,0,1000000);


//h2_ppacde_E= mkTH2F((char *)"ppacde_E",(char *)"ppacde_E",5000,0,50000,4000,0,40000);
//h2_ppacde_E_corr= mkTH2F((char *)"ppacde_E_corr",(char *)"ppacde_E_corr",5000,0,50000,4000,0,40000);

h1_ppacdeg = mkTH1D((char *)"ppacdeg",(char *)"ppacdeg",4000,0,40000);
h1_cl = mkTH1D((char *)"cl",(char *)"cl",4000,0,80000);
h1_cr = mkTH1D((char *)"cr",(char *)"cr",4000,0,80000);
h1_cx = mkTH1D((char *)"cx",(char *)"cx",4000,-100000,100000);
h1_cup = mkTH1D((char *)"cup",(char *)"cup",4000,0,80000);
h1_cdown = mkTH1D((char *)"cdown",(char *)"cdown",4000,0,80000);
h1_cy = mkTH1D((char *)"cy",(char *)"cy",4000,-100000,100000);

h1_cl_pu = mkTH1D((char *)"cl_pu",(char *)"cl_pu",10,0,10);
h1_cr_pu = mkTH1D((char *)"cr_pu",(char *)"cr_pu",10,0,10);
h1_ppacde_pu = mkTH1D((char *)"ppacde_pu",(char *)"ppacde_pu",10,0,10);

h1_csum = mkTH1D((char *)"csum",(char *)"csum",4000, 0,200000);
h1_csumy = mkTH1D((char *)"csumy",(char *)"csumy",4000, 0,200000);

h1_cxn = mkTH1D((char *)"cxn",(char *)"cxn",4000,0,4000);
h2_clr = mkTH2F((char *)"clr",(char *)"clr",2000,0,80000,2000,0,80000);
h2_cud = mkTH2F((char *)"cud",(char *)"cud",2000,0,80000,2000,0,80000);
h2_cxy = mkTH2F((char *)"cxy",(char *)"cxy",2000,-100000,100000,2000,-100000,100000);

h2_xde = mkTH2F((char *)"xde",(char *)"xde",1000,-100000,100000,1000,0,400000);
h2_xdeg = mkTH2F((char *)"xdeg",(char *)"xdeg",1000,-100000,100000,1000,0,400000);


h1_cx_en_g1 = mkTH1D((char *)"cx_en_g1",(char *)"cx_en_g1",4000,0,4000);
h1_cxn_en_g1 = mkTH1D((char *)"cxn_en_g1",(char *)"cxn_en_g1",4000,0,4000);
h1_cx_en_g2 = mkTH1D((char *)"cx_en_g2",(char *)"cx_en_g2",4000,0,4000);
h1_cxn_en_g2 = mkTH1D((char *)"cxn_en_g2",(char *)"cxn_en_g2",4000,0,4000);

h2_x_d_fr_e = mkTH2F((char *)"x_d_fr_e",(char *)"x_d_fr_e",2000,-100000,100000,5000,0,50000);
h2_x_d_fr_e_short = mkTH2F((char *)"x_d_fr_e_short",(char *)"x_d_fr_e_short",2000,-100000,100000,5000,0,50000);
h2_xehi = mkTH2F((char *)"xehi",(char *)"xehi",1000,-100000,100000,10000,0,10000);
h2_xehig = mkTH2F((char *)"xehig",(char *)"xehig",1000,-100000,100000,10000,0,10000);
h1_clg = mkTH1D((char *)"clg",(char *)"clg",4000,0,40000);
h1_crg = mkTH1D((char *)"crg",(char *)"crg",4000,0,40000);
//h2_clrg = mkTH2F((char *)"clrg",(char *)"clrg",4000,0,40000,4000,0,40000);
//h2_clrg_int = mkTH2F((char *)"clrg_int",(char *)"clrg_int",4000,0,40000,4000,0,40000); //new
h1_cxg = mkTH1D((char *)"cxg",(char *)"cxg",4000,0,4000);
h1_cxng = mkTH1D((char *)"cxng",(char *)"cxng",4000,0,4000);
//h2_clrg_en = mkTH2F((char *)"clrg_en",(char *)"clrg_en",4000,0,40000,4000,0,40000);
//h2_clrg_en_g1 = mkTH2F((char *)"clrg_en_g1",(char *)"clrg_en_g1",4000,0,40000,4000,0,40000);
//h2_clrg_en_g2 = mkTH2F((char *)"clrg_en_g2",(char *)"clrg_en_g2",4000,0,40000,4000,0,40000);
h2_traceg = mkTH2F((char *)"traceg",(char *)"traceg",1000,0,1000,1000,0,1000);
h2_tracegb14 = mkTH2F((char *)"tracegb14",(char *)"tracegb14",1000,0,1000,1000,0,1000);
h2_tracegb15 = mkTH2F((char *)"tracegb15",(char *)"tracegb15",1000,0,1000,1000,0,1000);

h2_traceg2 = mkTH2F((char *)"traceg2",(char *)"traceg2",1000,0,1000,1000,0,1000);
h2_traceg3 = mkTH2F((char *)"traceg3",(char *)"traceg3",1000,0,1000,1000,0,1000);
h1_trgn=mkTH1D((char *)"trgn",(char *)"trgn",2000,0,2000);
h1_trg2n=mkTH1D((char *)"trg2n",(char *)"trg2n",2000,0,2000);
h1_trg3n=mkTH1D((char *)"trg3n",(char *)"trg3n",2000,0,2000);

// DSSD stuff

h2_dssd_fr_emax= mkTH2F((char *)"dssd_fr_emax",(char *)"dssd_fr_emax",5000,0,50000,400,0,400);
h2_dssd_ba_emax= mkTH2F((char *)"dssd_ba_emax",(char *)"dssd_ba_emax",5000,0,50000,400,0,400);

h2_dssd_fr_emax_phy= mkTH2F((char *)"dssd_fr_emax_phy",(char *)"dssd_fr_emax_phy",5000,0,500000,400,0,400);
h2_dssd_ba_emax_phy= mkTH2F((char *)"dssd_ba_emax_phy",(char *)"dssd_ba_emax_phy",5000,0,500000,400,0,400);

h2_d_fr_e= mkTH2F((char *)"d_fr_e",(char *)"d_fr_e",7500,0,30000,400,0,400);
h2_d_ba_e= mkTH2F((char *)"d_ba_e",(char *)"d_ba_e",7500,0,30000,400,0,400);
h2_r_fr_e= mkTH2F((char *)"r_fr_e",(char *)"r_fr_e",7500,0,30000,400,0,400);
h2_r_ba_e= mkTH2F((char *)"r_ba_e",(char *)"r_ba_e",7500,0,30000,400,0,400);

/*INIT*/

h1_basesample=mkTH1D((char *)"bs",(char *)"bs",2000,0,20000);

h1_header_type=mkTH1D((char *)"header",(char *)"header",100,0,100);

 h2_basesample= mkTH2F((char *)"basesample",(char *)"basesample",2000,0,20000,200,0,200);
 h2_baseline= mkTH2F((char *)"baseline",(char *)"baseline",2000,0,20000,200,0,200);

 h2_blav= mkTH2F((char *)"blav",(char *)"blav",2000,0,20000,200,0,200);

 h2_postrisebeg= mkTH2F((char *)"postrisebeg",(char *)"postrisebeg",2000,0,20000,200,0,200);
 h2_postriseend= mkTH2F((char *)"postriseend",(char *)"postriseend",2000,0,20000,200,0,200);
 h2_prerisebeg= mkTH2F((char *)"prerisebeg",(char *)"prerisebeg",2000,0,20000,200,0,200);
 h2_preriseend= mkTH2F((char *)"preriseend",(char *)"preriseend",2000,0,20000,200,0,200);

 h2_pre2= mkTH2F((char *)"pre2",(char *)"pre2",2000,0,20000,2000,0,20000);
 h2_post2= mkTH2F((char *)"post2",(char *)"post2",2000,0,20000,2000,0,20000);
h2_last2= mkTH2F((char *)"last2",(char *)"last2",2000,0,20000,2000,0,20000);

h2_e_bs= mkTH2F((char *)"e_bs",(char *)"e_bs",200,-1000,1000,5000,0,60000);
h2_e_dt= mkTH2F((char *)"e_dt",(char *)"e_dt",2000,0,20000,5000,0,50000);

h2_dssd_en_raw= mkTH2F((char *)"dssd_en_raw",(char *)"dssd_en_raw",30000,0,300000,400,0,400);
//h2_dssd_en = mkTH2F((char *)"dssd_en",(char *)"dssd_en",20000,0,20000,400,0,400);
h2_dssd_en = mkTH2F((char *)"dssd_en",(char *)"dssd_en",2500,0,10000,400,0,400);
h2_dssd_fb = mkTH2F((char *)"dssd_fb",(char *)"dssd_fb",2000,0,40000,2000,0,40000);

h2_dssd_rate = mkTH2F((char *)"dssd_rate",(char *)"dssd_rate",5000,0,5000,500,0,500);
//h1_dssd_rate_1us = mkTH1D((char *)"dssd_rate_1us",(char *)"dssd_rate_1us",1000000,0,1000000);
h1_decay_rate = mkTH1D((char *)"decay_rate",(char *)"decay_rate",100000,0,100000);
h1_recoil_rate = mkTH1D((char *)"recoil_rate",(char *)"recoil_rate",100000,0,100000);
h1_corr_rate = mkTH1D((char *)"corr_rate",(char *)"corr_rate",100000,0,100000);

h1_dssd_rate_1D = mkTH1D((char *)"dssd_rate_1D",(char *)"dssd_rate_1D",100000,0,100000);
h2_FP_rate = mkTH2F((char *)"FP_rate",(char *)"FP_rate",5000,0,5000,10,0,10);
h1_FP_rate_left = mkTH1D((char *)"FP_rate_left",(char *)"FP_rate_left",100000,0,100000);
h1_FP_rate_right = mkTH1D((char *)"FP_rate_right",(char *)"FP_rate_right",100000,0,100000);
h1_FP_rate_up = mkTH1D((char *)"FP_rate_up",(char *)"FP_rate_up",100000,0,100000);
h1_GE_rate = mkTH1D((char *)"GE_rate",(char *)"GE_rate",100000,0,100000);
h1_cl_int = mkTH1D((char *)"cl_int",(char *)"cl_int",4000,0,40000);
h1_cr_int = mkTH1D((char *)"cr_int",(char *)"cr_int",4000,0,40000);
h1_cxg_int = mkTH1D((char *)"cxg_int",(char *)"cxg_int",4000,0,4000);
h1_cxng_int = mkTH1D((char *)"cxng_int",(char *)"cxng_int",4000,0,4000);
h2_clr_int = mkTH2F((char *)"clr_int",(char *)"clr_int",4000,0,40000,4000,0,40000);

h2_tdssdmcp_fr = mkTH2F((char *)"tdssdmcp_fr",(char *)"tdssdmcp_fr",1000,-5000,5000,200,0,200);
h2_tdssdmcp_ba = mkTH2F((char *)"tdssdmcp_ba",(char *)"tdssdmcp_ba",1000,-5000,5000,200,0,200);
h2_tdssdmcp_fr_r = mkTH2F((char *)"tdssdmcp_fr_r",(char *)"tdssdmcp_fr_r",1000,-5000,5000,200,0,200);
h2_tdssdmcp_ba_r = mkTH2F((char *)"tdssdmcp_ba_r",(char *)"tdssdmcp_ba_r",1000,-5000,5000,200,0,200);
h2_tdssdppacde_fr = mkTH2F((char *)"tdssdppacde_fr",(char *)"tdssdppacde_fr",1000,-5000,5000,200,0,200);

h2_tdssdppacde = mkTH2F((char *)"tdssdppacde",(char *)"tdssdppacde",1000,-5000,5000,5000,0,50000);

h1_ndssd = mkTH1D((char *)"ndssd",(char *)"ndssd",350,0,350);

h2_lefttraces = mkTH2F((char *)"lefttraces",(char *)"lefttraces",100,0,100,1000,0,1000);
h2_righttraces = mkTH2F((char *)"righttraces",(char *)"righttraces",100,0,100,1000,0,1000);

h2_dssd_traces_fr = mkTH2F((char *)"dssd_traces_fr",(char *)"dssd_traces_fr",100,0,100,1000,0,1000);
h2_dssd_traces_ba = mkTH2F((char *)"dssd_traces_ba",(char *)"dssd_traces_ba",100,0,100,1000,0,1000);


h1_trlen = mkTH1D((char *)"trlen",(char *)"trlen",4000,0,4000);

h2_frba_t = mkTH2F((char *)"frba_t",(char *)"frba_t",1000,-500,500,350,0,350);


h1_ch0 = mkTH1D((char *)"ch0",(char *)"ch0",4000,0,40000);
h1_ch1 = mkTH1D((char *)"ch1",(char *)"ch1",4000,0,40000);
h1_ch2 = mkTH1D((char *)"ch2",(char *)"ch2",4000,0,40000);
h1_ch3 = mkTH1D((char *)"ch3",(char *)"ch3",4000,0,40000);
h1_ch4 = mkTH1D((char *)"ch4",(char *)"ch4",4000,0,40000);
h1_ch5 = mkTH1D((char *)"ch5",(char *)"ch5",4000,0,40000);

h2_r_hitxy = mkTH2F((char *)"r_hitxy",(char *)"r_hitxy",200,0,200,200,0,200);
h2_d_hitxy = mkTH2F((char *)"d_hitxy",(char *)"d_hitxy",200,0,200,200,0,200);

h2_dssd_hitxy_phys = mkTH2F((char *)"dssd_hitxy_phys",(char *)"dssd_hitxy_phys",200,0,200,200,0,200);
h2_dssd_hitxy_physg = mkTH2F((char *)"dssd_hitxy_physg",(char *)"dssd_hitxy_physg",200,0,200,200,0,200);
//h2_dssd_hitxy_physg251Md = mkTH2F((char *)"dssd_hitxy_physg251Md",(char *)"dssd_hitxy_physg251Md",200,0,200,200,0,200);
h2_dssd_hitxy_phys_corr= mkTH2F((char *)"dssd_hitxy_phys_corr",(char *)"dssd_hitxy_phys_corr",200,0,200,200,0,200);
h2_dssd_hitxy_phys_corrf= mkTH2F((char *)"dssd_hitxy_phys_corrf",(char *)"dssd_hitxy_phys_corrf",200,0,200,200,0,200);

h2_dssd_fb_corr= mkTH2F((char *)"dssd_fb_corr",(char *)"dssd_fb_corr",1000,0,10000,1000,0,10000);

h2_dssd_fr_p2 = mkTH2F((char *)"dssd_fr_p2",(char *)"dssd_fr_p2",5000,0,50000,400,0,400);
h2_dssd_ba_p2 = mkTH2F((char *)"dssd_ba_p2",(char *)"dssd_ba_p2",5000,0,50000,400,0,400);
h2_clr_p2 = mkTH2F((char *)"clr_p2",(char *)"clr_p2",4000,0,40000,4000,0,40000);

h2_dssd_fr_p1 = mkTH2F((char *)"dssd_fr_p1",(char *)"dssd_fr_p1",5000,0,50000,400,0,400);
h2_dssd_ba_p1 = mkTH2F((char *)"dssd_ba_p1",(char *)"dssd_ba_p1",5000,0,50000,400,0,400);
//h2_clr_p1 = mkTH2F((char *)"clr_p1",(char *)"clr_p1",4000,0,40000,4000,0,40000);


//Correlations


h2_e1t1log = mkTH2F((char *)"e1t1log",(char *)"e1t1log",5000,0,15000,200,0,20);
h2_e1t1log_box = mkTH2F((char *)"e1t1log_box",(char *)"e1t1log_box",6000,-10000,50000,200,0,20);
h2_e1t1log_box1 = mkTH2F((char *)"e1t1log_box1",(char *)"e1t1log_box1",5000,0,15000,200,0,20);
//h2_e1t1c1 = mkTH2F((char *)"e1t1c1",(char *)"e1t1c1",600,-30000,30000,1000,0,1000);
//h2_e1t1c2 = mkTH2F((char *)"e1t1c2",(char *)"e1t1c2",600,-30000,30000,1000,0,1000);
//h2_e1t1c3 = mkTH2F((char *)"e1t1c3",(char *)"e1t1c3",600,-30000,30000,1000,0,1000);
//h2_e1t1c4 = mkTH2F((char *)"e1t1c4",(char *)"e1t1c4",600,-30000,30000,1000,0,1000);


h2_e2t2log = mkTH2F((char *)"e2t2log",(char *)"e2t2log",6000,-10000,50000,200,0,20);
//h2_e2t2c1 = mkTH2F((char *)"e2t2c1",(char *)"e2t2c1",600,-30000,30000,1000,0,1000);
//h2_e2t2c2 = mkTH2F((char *)"e2t2c2",(char *)"e2t2c2",600,-30000,30000,1000,0,1000);
//h2_e2t2c3 = mkTH2F((char *)"e2t2c3",(char *)"e2t2c3",600,-30000,30000,1000,0,1000);
//h2_e2t2c4 = mkTH2F((char *)"e2t2c4",(char *)"e2t2c4",600,-30000,30000,1000,0,1000);


h2_e3t3log = mkTH2F((char *)"e3t3log",(char *)"e3t3log",600,-30000,30000,200,0,20);
//h2_e3t3c1 = mkTH2F((char *)"e3t3c1",(char *)"e3t3c1",600,-30000,30000,1000,0,1000);
//h2_e3t3c2 = mkTH2F((char *)"e3t3c2",(char *)"e3t3c2",600,-30000,30000,1000,0,1000);
//h2_e3t3c3 = mkTH2F((char *)"e3t3c3",(char *)"e3t3c3",600,-30000,30000,1000,0,1000);
//h2_e3t3c4 = mkTH2F((char *)"e3t3c4",(char *)"e3t3c4",600,-30000,30000,1000,0,1000);

h2_e4t4log = mkTH2F((char *)"e4t4log",(char *)"e4t4log",6000,-10000,50000,200,0,20);
h2_e5t5log = mkTH2F((char *)"e5t5log",(char *)"e5t5log",6000,-10000,50000,200,0,20);
h2_e6t6log = mkTH2F((char *)"e6t6log",(char *)"e6t6log",6000,-10000,50000,200,0,20);


//h2_e1t1logf = mkTH2F((char *)"e1t1logf",(char *)"e1t1logf",6000,-10000,50000,200,0,20);
//h2_e2t2logf = mkTH2F((char *)"e2t2logf",(char *)"e2t2logf",6000,-10000,50000,200,0,20);
//h2_e3t3logf = mkTH2F((char *)"e3t3logf",(char *)"e3t3logf",600,-30000,30000,200,0,20);
//h2_e4t4logf = mkTH2F((char *)"e4t4logf",(char *)"e4t4logf",60,-10000,50000,200,0,20);//was 6000 x 200, now 60 x 200 bins, kalle
//h2_e5t5logf = mkTH2F((char *)"e5t5logf",(char *)"e5t5logf",60,-10000,50000,200,0,20);//was 6000 x 200, now 60 x 200 bins, kalle
//h2_e6t6logf = mkTH2F((char *)"e6t6logf",(char *)"e6t6logf",60,-10000,50000,200,0,20);//was 6000 x 200, now 60 x 200 bins, kalle


h2_eftof = mkTH2F((char *)"eftof",(char *)"eftof",5000,0,200000,2000,-1000,1000);
h2_eftofg = mkTH2F((char *)"eftofg",(char *)"eftofg",5000,0,200000,2000,-1000,1000);
h2_eftof_corr = mkTH2F((char *)"eftof_corr",(char *)"eftof_corr",5000,0,200000,2000,-1000,1000);


h2_e1e2 = mkTH2F((char *)"e1e2",(char *)"e1e2",2000, 0, 20000,2000, 0, 20000);
h2_e1e2_short = mkTH2F((char *)"e1e2_short",(char *)"e1e2_short", 2000, 0, 20000, 2000, 0, 20000);

h2_e2e3 = mkTH2F((char *)"e2e3",(char *)"e2e3",2000, 0, 20000,2000, 0, 20000);
h2_e3e4 = mkTH2F((char *)"e3e4",(char *)"e3e4",2000, 0, 20000,2000, 0, 20000);
h2_e4e5 = mkTH2F((char *)"e4e5",(char *)"e4e5",2000, 0, 20000,2000, 0, 20000);
h2_e5e6 = mkTH2F((char *)"e5e6",(char *)"e5e6",2000, 0, 20000,2000, 0, 20000);
h2_e2e3f = mkTH2F((char *)"e2e3f",(char *)"e2e3f",20, 0, 20000,20, 0, 20000);
h2_e3e4f = mkTH2F((char *)"e3e4f",(char *)"e3e4f",20, 0, 20000,20, 0, 20000);//was 2000 x 2000, now 20 x 20 bins, kalle
h2_e4e5f = mkTH2F((char *)"e4e5f",(char *)"e4e5f",20, 0, 20000,20, 0, 20000);//was 2000 x 2000, now 20 x 20 bins, kalle
h2_e5e6f = mkTH2F((char *)"e5e6f",(char *)"e5e6f",20, 0, 20000,20, 0, 20000);//was 2000 x 2000, now 20 x 20 bins, kalle

corr_time_short = 10000000000;
corr_time = 1000000000000;

//Clovers
char str_1[50];

h1_clover_en = mkTH1D((char *)"clover_en",(char *)"clover_en",400000,0,400000);



//h1_dec_clover = mkTH1D((char *)"dec_clover",(char *)"dec_clover",40000,0,40000);
//h2_cle_id_betag = mkTH2F((char *)"cle_id_betag",(char *)"cle_id_betag",40000,0,40000,31,0,31);

h2_cle_id = mkTH2F((char *)"cle_id",(char *)"cle_id",4000,0,4000,30,1,31);

//h2_cle_id_keV = mkTH2F((char *)"cle_id_keV",(char *)"cle_id_keV",10000,0,10000,31,0,31);

 h2_cle_id_raw = mkTH2F((char *)"cle_id_raw",(char *)"cle_id_raw",15000,0,15000,31,0,31);

/*
for (i = 1; i < 21; i++ ){
    sprintf(str_1, "cl_ehi_basesample%i", i) ;
    h2_cl_ehi_basesample[i] = mkTH2F(str_1,str_1,5000,0,5000,1000,7500,8500);
}
for (i = 1; i < 21; i++ ){
    sprintf(str_1, "cl_ehi_slope%i", i) ;
    h2_cl_ehi_slope[i] = mkTH2F(str_1,str_1,5000,0,5000,1000,-500,500);
}
*/

// BOX histograms

h2_sibox_raw = mkTH2F((char *)"sibox_raw",(char *)"sibox_raw",5000,0,50000,60,0,60);
h2_sibox = mkTH2F((char *)"sibox",(char *)"sibox",5000,0,50000,60,0,60);

h2_dtdssdfrsi = mkTH2F((char *)"dtdssdfrsi",(char *)"dtdssdfrsi",1000,-500,500,200,0,200);
h2_dtdssdbasi = mkTH2F((char *)"dtdssdbasi",(char *)"dtdssdbasi",1000,-500,500,200,0,200);
h2_dtdssdfrsi_corr = mkTH2F((char *)"dtdssdfrsi_corr",(char *)"dtdssdfrsi_corr",1000,-500,500,200,0,200);
h2_dtdssdbasi_corr = mkTH2F((char *)"dtdssdbasi_corr",(char *)"dtdssdbasi_corr",1000,-500,500,200,0,200);


h2_sisi = mkTH2F((char *)"sisi",(char *)"sisi",1000,0,20000,1000,0,20000);
h2_sidssd = mkTH2F((char *)"sidssd",(char *)"sidssd",1000,0,20000,1000,0,20000);

h2_sidssdcor = mkTH2F((char *)"sidssdcor",(char *)"sidssdcor",1000,0,20000,1000,0,20000);
h2_sidssdcor2 = mkTH2F((char *)"sidssdcor2",(char *)"sidssdcor2",1000,0,20000,1000,0,20000);

h2_ada1 = mkTH2F((char *)"ada1",(char *)"ada1",1000,0,20000,1000,0,20000);
h2_ada2 = mkTH2F((char *)"ada2",(char *)"ada2",1000,0,20000,1000,0,20000);
h2_ada3 = mkTH2F((char *)"ada3",(char *)"ada3",1000,0,20000,1000,0,20000);
h2_ada4 = mkTH2F((char *)"ada4",(char *)"ada4",1000,0,20000,1000,0,20000);
h2_ada5 = mkTH2F((char *)"ada5",(char *)"ada5",1000,0,20000,1000,0,20000);
h2_ada6 = mkTH2F((char *)"ada6",(char *)"ada6",1000,0,20000,1000,0,20000);
h2_ada7 = mkTH2F((char *)"ada7",(char *)"ada7",1000,0,20000,1000,0,20000);
h2_ada8 = mkTH2F((char *)"ada8",(char *)"ada8",1000,0,20000,1000,0,20000);


h2_sisicor2esc1 = mkTH2F((char *)"sisicor2esc1",(char *)"sisicor2esc1",1000,0,20000,1000,0,20000);
h2_sisicor2esc2 = mkTH2F((char *)"sisicor2esc2",(char *)"sisicor2esc2",1000,0,20000,1000,0,20000);

h2_aa1esc = mkTH2F((char *)"aa1esc",(char *)"aa1esc",1000,0,20000,1000,0,20000);
h2_aa2esc = mkTH2F((char *)"aa2esc",(char *)"aa2esc",1000,0,20000,1000,0,20000);

h1_atot = mkTH1D((char *)"atot",(char *)"atot",40000,0,40000);
h1_atotcor = mkTH1D((char *)"atotcor",(char *)"atotcor",40000,0,40000);

// correlated histograms

h2_sisi_corr = mkTH2F((char *)"sisi_corr",(char *)"sisi_corr",1000,0,20000,1000,0,20000);
h2_sidssd_corr = mkTH2F((char *)"sidssd_corr",(char *)"sidssd_corr",1000,0,20000,1000,0,20000);


// WHEEL

h1_esidt = mkTH1D((char *)"esidt",(char *)"esidt",10000,0,10000);
h1_esidt2 = mkTH1D((char *)"esidt2",(char *)"esidt2",10000,0,10000);
h1_leftdt = mkTH1D((char *)"leftdt",(char *)"leftdt",10000,0,10000);
h1_leftdt2 = mkTH1D((char *)"leftdt2",(char *)"leftdt2",10000,0,10000);
h2_leftdtr = mkTH2F((char *)"leftdtr",(char *)"leftdtr",1000,0,1000,2000,0,2000);
h2_esidtr = mkTH2F((char *)"esidtr",(char *)"esidtr",1000,0,1000,2000,0,2000);
h2_esidtrg = mkTH2F((char *)"esidtrg",(char *)"esidtrg",1000,0,1000,2000,0,2000);
h2_leftdtam = mkTH2F((char *)"leftdtam",(char *)"leftdtam",1000,0,1000,2000,0,2000);


h2_ts2 = mkTH2F((char *)"ts2",(char *)"ts2",10000,0,10000,1000,0,1000);
h2_ts2c = mkTH2F((char *)"ts2c",(char *)"ts2c",10000,0,10000,1000,0,1000);

h1_fdecEv =  mkTH1D((char *)"fdecEv",(char *)"fdecEv",100,0,100);

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

/*

  //recalibration for the DSSD by Tianheng
 FILE *fileth;
 char strth[255];
 float ath[161],bth[161]; 
 fileth =fopen("F.cal","r");
 fgets(strth,255,fileth);
 fgets(strth,255,fileth);
 for(i=1;i<161;i++)
   {
     fscanf(fileth,"%d",&strip);
     if(strip==i)fscanf(fileth,"%f %f",&ath[i],&bth[i]);
     else {return -3;}
    fgets(strth,255,fileth); 
   }
 fclose(fileth);

  for(i=1;i<161;i++)
   {
     map_fr[i].gain = map_fr[i].gain * ath[i];
     map_fr[i].off = map_fr[i].off *ath[i] + bth[i];
   }

 fileth =fopen("B.cal","r");
 fgets(strth,255,fileth);
 fgets(strth,255,fileth);
 for(i=1;i<161;i++)
   {
     fscanf(fileth,"%d",&strip);
     if(strip==i)fscanf(fileth,"%f %f",&ath[i],&bth[i]);
     else {return -3;}
    fgets(strth,255,fileth); 
   }
 fclose(fileth);

  for(i=1;i<161;i++)
   {
     map_ba[i].gain = map_ba[i].gain * ath[i];
     map_ba[i].off = map_ba[i].off *ath[i] + bth[i];
   }

*/

// CLOVER MAP FILE

#if(1)
 fmap1=fopen("MAP_FILES/clovers.map","read");

 for(i=0;i<20;i++){

   fscanf(fmap1,"%d %f %f", &strip, &gain, &off);
   clmap[strip].gain = gain;
   clmap[strip].off = off;

  printf("Clover map channel %d has offset %f and gain %f \n",strip,off,gain);

 }

fclose(fmap1);

#endif

// SIBOX MAP FILE

#if(1)

//fmap1=fopen("MAP_FILES/sibox_calib_z115.map","read");
fmap1=fopen("MAP_FILES/sibox_calib_th.map","read");
//fmap1=fopen("MAP_FILES/sibox_test.map","read");

 for(i=0;i<NUMSTRIPSBOX;i++){

   fscanf(fmap1,"%d %d %d %f %f %i", &strip, &phystrip, &thr, &off, &gain, &baseline);
 

   if (phystrip!=0) {  

   map_box[strip].thr = thr;
   map_box[strip].phystrip = phystrip;
   map_box[strip].off = off; 
   map_box[strip].gain = gain; 
   map_box[strip].baseline = baseline;

  } else {

   map_box[strip].thr = 100;
   map_box[strip].phystrip = phystrip;
   map_box[strip].off = 0;
   map_box[strip].gain = 1.0;  
   map_box[strip].baseline = 0;

  }


  printf("map back strip %d (phys strip %d) has thr %d, offset %f and gain %f and baseline %i\n",strip,phystrip,thr,off,gain,baseline);

}

fclose(fmap1);

// window



// = new TCUT("eftof");

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



//if (first && (numDFMA > 1000) && (numDGS > 1000)) { t_first=DFMAEvent[0].LEDts; first=false; }
if (first && (numDFMA > 100)) { t_first=DFMAEvent[0].LEDts; first=false; }
if(Pars.CurEvNo%10000 == 0) printf("Processing event number %i with timestamp %llu. Time elapsed: %i seconds\n",Pars.CurEvNo,DFMAEvent[0].LEDts,int(float(DFMAEvent[0].LEDts - t_first)/1E8));


 if ((numDFMA > 100) ) {

 h2_ts2->Fill((DFMAEvent[0].LEDts-t_first)/100, 990);

 h2_ts2->Fill((DGSEvent[0].event_timestamp-t_first)/100, 970);

for (i=0;i<nsubev;i++) {

  if(DFMAEvent[i].tpe == DSSD){ h2_ts2->Fill((DFMAEvent[i].LEDts-t_first)/100, DFMAEvent[i].tid+400); }

  if(DFMAEvent[i].tpe == FP ){ h2_ts2->Fill((DFMAEvent[i].LEDts-t_first)/100, DFMAEvent[i].tid+800); }

  if(DFMAEvent[i].tpe == SIBOX ){ h2_ts2->Fill((DFMAEvent[i].LEDts-t_first)/100, DFMAEvent[i].tid+900); }



  if(DFMAEvent[i].tpe == DSSD){ h2_ts2c->Fill((DFMAEvent[i].LEDts-t_first)/100000, DFMAEvent[i].tid+400); }

  if(DFMAEvent[i].tpe == FP ){ h2_ts2c->Fill((DFMAEvent[i].LEDts-t_first)/100000, DFMAEvent[i].tid+800); }

  if(DFMAEvent[i].tpe == SIBOX ){ h2_ts2c->Fill((DFMAEvent[i].LEDts-t_first)/100000, DFMAEvent[i].tid+900); }

}

for (i=0;i<ng;i++) {

  if(DGSEvent[i].tpe == GE) { h2_ts2->Fill((DGSEvent[i].event_timestamp-t_first)/100, DGSEvent[i].tid); }

  if(DGSEvent[i].tpe == BGO) { h2_ts2->Fill((DGSEvent[i].event_timestamp-t_first)/100, DGSEvent[i].tid+200); }

  if(DGSEvent[i].tpe == GE) { h2_ts2c->Fill((DGSEvent[i].event_timestamp-t_first)/100000, DGSEvent[i].tid); }

  if(DGSEvent[i].tpe == BGO) { h2_ts2c->Fill((DGSEvent[i].event_timestamp-t_first)/100000, DGSEvent[i].tid+200); }


}

 }


if(firstge) { t_firstge = DGSEvent[0].event_timestamp; firstge = false; }

long long int tmp1, tmp2;


// Declarations for variables


   // NEW STUFF HERE>>>>>>>>>>
int l_num_bkgd, l_num_sig;
l_num_bkgd = 95;
l_num_sig = 400;
float l_bkgdsum, l_signal, l_avgbkgd;
l_bkgdsum = 0.0;
l_signal = 0.0;
l_avgbkgd = 0.0;
float cl_int;
cl_int = 0;

int r_num_bkgd, r_num_sig;
r_num_bkgd = 95;
r_num_sig = 350;
float r_bkgdsum, r_signal, r_avgbkgd;
r_bkgdsum = 0.0;
r_signal = 0.0;
r_avgbkgd = 0.0;
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

	   if (DFMAEvent[i].tid==1) {
             if (trn1<1000) {

	       for (jj=0;jj<DSSDTRLEN;jj++) {

	       short int f;
               float ff;               

               f= (short int) (DFMAEvent[i].trace[jj] & 0x3fff);
               ff=f;
               //printf("nob %d %d %f\n", jj, f, ff);
               //h2_traceg->Fill(trn1,jj,ff);
               h2_traceg->Fill(trn1,jj,(DFMAEvent[i].trace[jj] & 0x3fff));
               f= (short int) (DFMAEvent[i].trace[jj]);
               ff=f;
               //printf("withb %d %d %f\n", jj, f, ff);
	       //h2_traceg2->Fill(trn1,jj,ff);
               h2_traceg2->Fill(trn1,jj,DFMAEvent[i].trace[jj]);
              //h2_traceg->SetBinContent(trn1,999,s_fd_fr);
	      if ((DFMAEvent[i].trace[jj] & 0x4000) != 0) h2_tracegb14->Fill(trn1,jj,1);
	      if ((DFMAEvent[i].trace[jj] & 0x8000) != 0) h2_tracegb15->Fill(trn1,jj,1);
	       }
               trn1=trn1+1;
             }
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

   //if (trf==NULL) (printf("TRACE FILE DOES NOT WRITE\n")); 

   /*
   fprintf(trf,"%4d %4d ",DFMAEvent[dssd_fr_subev].tid, DFMAEvent[dssd_fr_subev].traceLen);

   for (i=0;i<DFMAEvent[dssd_fr_subev].traceLen;i++) {
     fprintf(trf,"%6d",DFMAEvent[dssd_fr_subev].trace[i]);
   }

   fprintf(trf,"\n");

   fprintf(trf,"%4d %4d",DFMAEvent[dssd_ba_subev].tid, DFMAEvent[dssd_fr_subev].traceLen);

   for (i=0;i<DFMAEvent[dssd_ba_subev].traceLen;i++) {
     fprintf(trf,"%6d",DFMAEvent[dssd_ba_subev].trace[i]);
   }

   fprintf(trf,"\n");
   */

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
 //idnotpresent=1;

 

for (i=0;i<nsubev;i++) {
  
  switch (DFMAEvent[i].tpe) {
    
         //-------------\\
         // Focal Plane \\
         //-------------\\

         case SIBOX:
 
           // check if this id was already present
           idnotpresent=1;
	   for (j=0;j<n_sib;j++) {
	     if (sib_tid[j]==DFMAEvent[i].tid) idnotpresent=0; 
           }         

	   if (idnotpresent==1) {
            cal_e = double(map_box[DFMAEvent[i].tid].gain)*(DFMAEvent[i].ehi + double(rand())/RAND_MAX-0.5) + double(map_box[DFMAEvent[i].tid].off); 
            h2_sibox->Fill(cal_e,DFMAEvent[i].tid);
            h2_sibox_raw->Fill(DFMAEvent[i].ehi,DFMAEvent[i].tid);
            //h2_dssd_rate->Fill((DFMAEvent[i].LEDts - t_first)/1E5, DFMAEvent[i].tid);
            sib_ts[n_sib]=DFMAEvent[i].LEDts;
            sib_e[n_sib]=cal_e;
            sib_tid[n_sib]=DFMAEvent[i].tid;
            if (cal_e>sib_emax) sib_emax=cal_e;
            sib_wn[n_sib]=sib_wn_table[DFMAEvent[i].tid];
            sib_stripn[n_sib]=sib_stripn_table[DFMAEvent[i].tid];
            sib_detn[n_sib]=sib_detn_table[DFMAEvent[i].tid];

            n_sib++;
	    printf("ind box %d %d %d %f\n",n_sib,DFMAEvent[i].tid,DFMAEvent[i].ehi,cal_e);
	   }


         break; 

	 case FP:

           //h2_dssd_rate->Fill((DFMAEvent[i].LEDts - t_first)/1E5, DFMAEvent[i].tid);


  	 if(DFMAEvent[i].tid == 1) h1_ch0->Fill(DFMAEvent[i].ehi);
   	 if(DFMAEvent[i].tid == 2) h1_ch1->Fill(DFMAEvent[i].ehi);
   	 if(DFMAEvent[i].tid == 3) h1_ch2->Fill(DFMAEvent[i].ehi);
   	 if(DFMAEvent[i].tid == 4) h1_ch3->Fill(DFMAEvent[i].ehi);
   	 if(DFMAEvent[i].tid == 5) h1_ch4->Fill(DFMAEvent[i].ehi);
   	 if(DFMAEvent[i].tid == 6) h1_ch5->Fill(DFMAEvent[i].ehi);
   	 if(DFMAEvent[i].tid == SICH && DFMAEvent[i].ehi > 0) gate = 1;


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
		/*		
		if(tdssdmcp_fr<1 && tdssdmcp_fr>-1)
		  {
	        printf("Testing the TOF! %f %f %f\n",tdssdmcp_fr,DFMAEvent[i].LEDts,tdssd_fr);
		  }
		*/
                //left_ts=DFMAEvent[i].LEDts;

                if(tdssd_fr != 0) h2_tdssdmcp_fr->Fill(tdssdmcp_fr,DFMAEvent[dssd_fr_subev].tid);
                if(tdssd_ba != 0) h2_tdssdmcp_ba->Fill(tdssdmcp_ba,DFMAEvent[dssd_ba_subev].tid-160);


    		cl=DFMAEvent[i].ehi;
 
                h1_leftdt->Fill(DFMAEvent[i].LEDts/10000-ts_wheel/10000);
                h1_leftdt2->Fill(DFMAEvent[i].LEDts/10000-ts_wheel2/10000);
                h2_leftdtr->Fill(DFMAEvent[i].LEDts/10000000000-t_first/10000000000,DFMAEvent[i].LEDts/10000-ts_wheel/10000);

                ts_left=DFMAEvent[i].LEDts;

                wheel_left=DFMAEvent[i].wheel;
                h1_wheel_left->Fill(DFMAEvent[i].wheel);

                /*


 		if (trn1<1000) {
    		   for (jj=0;jj<1000;jj++) h2_traceg->Fill(trn1,jj,DFMAEvent[i].trace[jj]);
       		   trn1=trn1+1;
  		}


		// Trace integration >>>>>>>>>>

                // Firstly, find trigger point...

                

                for(j=0;j<1000;j++) letrace[j] = 0;
                let0 = 0;
                let1 = 0;

                leFilter(DFMAEvent[i].trace, 1000, 20, letrace, let0, let1, 200); 


                if(leTr<100){
                   for(j=0;j<1000;j++) h2_le_output_left->Fill(leTr,j,letrace[j]);
                   leTr++;
                }


		for(j=5;j<(5+l_num_bkgd);j++){
		   l_bkgdsum = l_bkgdsum + DFMAEvent[i].trace[j];
		}
            
		l_avgbkgd = (l_bkgdsum/l_num_bkgd);

                for(j=let0; j<(let0+l_num_sig); j++){
		   l_signal = l_signal + DFMAEvent[i].trace[j];
		}

		cl_int = l_signal - (l_avgbkgd*l_num_sig);

                */

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

                /*


 	     	if (trn2<1000) {
    		   for (jj=0;jj<1000;jj++) h2_traceg2->Fill(trn2,jj,DFMAEvent[i].trace[jj]);
    		   trn2=trn2+1;
  	        }


	        // Trace integration>>>>>>>>>>
                // Firstly, find trigger point...



                for(j=0;j<1000;j++) letrace[j] = 0;
                let0 = 0;
                let1 = 0;

                leFilter(DFMAEvent[i].trace, 1000, 20, letrace, let0, let1, 200); 

		for(j=5;j<(5+r_num_bkgd);j++){
		   r_bkgdsum = r_bkgdsum + DFMAEvent[i].trace[j];
		}
            
		r_avgbkgd = (r_bkgdsum/r_num_bkgd);

                for(j=let0;j<(let0+r_num_sig);j++){
		   r_signal = r_signal + DFMAEvent[i].trace[j];
		}

		cr_int = r_signal - (r_avgbkgd*r_num_sig);

                */

		//>>>>>>>>>>>>>>>>>>>>>>>>>>

	}

	 // MWPC UP
	 if (DFMAEvent[i].tid==UPCH && DFMAEvent[i].ehi > 0){

                up_subev = i;

                //cr_pu = DFMAEvent[i].pu;

             	//h2_FP_rate->Fill((DFMAEvent[i].LEDts - t_first)/1E8, DFMAEvent[i].tid);
                //if((DFMAEvent[i].LEDts - t_first)/1E5 < 100000){
         	//    h1_FP_rate_right->Fill((DFMAEvent[i].LEDts - t_first)/1E5);
                //}

                //tdssdmcp_fr_r = (DFMAEvent[i].LEDts) -  (tdssd_fr);
                //tdssdmcp_ba_r = (DFMAEvent[i].LEDts) -  (tdssd_ba);

                //f(tdssd_fr != 0) h2_tdssdmcp_fr_r->Fill(tdssdmcp_fr_r,DFMAEvent[dssd_fr_subev].tid);
		  //if(tdssd_ba != 0) h2_tdssdmcp_ba_r->Fill(tdssdmcp_ba_r,DFMAEvent[dssd_ba_subev].tid-160);

  	    	cup=DFMAEvent[i].ehi;
		h1_cup->Fill(cup);

                wheel_up=DFMAEvent[i].wheel;
		if (cup>15000) h1_wheel_up->Fill(DFMAEvent[i].wheel);
 
               if((DFMAEvent[i].LEDts - t_first)/ratecomp < 100000){
		 //if ((cup>25000)&&(cup<30000)) 
                  if ((cup>100)&&(cup<100000)) h1_FP_rate_up->Fill((DFMAEvent[i].LEDts - t_first)/ratecomp);
                }


	}

	 // MWPC UP
	 if (DFMAEvent[i].tid==DOWNCH && DFMAEvent[i].ehi > 0){

                down_subev = i;

                //cr_pu = DFMAEvent[i].pu;

             	//h2_FP_rate->Fill((DFMAEvent[i].LEDts - t_first)/1E8, DFMAEvent[i].tid);
                //if((DFMAEvent[i].LEDts - t_first)/1E5 < 100000){
         	//    h1_FP_rate_right->Fill((DFMAEvent[i].LEDts - t_first)/1E5);
                //}

                //tdssdmcp_fr_r = (DFMAEvent[i].LEDts) -  (tdssd_fr);
                //tdssdmcp_ba_r = (DFMAEvent[i].LEDts) -  (tdssd_ba);

                //f(tdssd_fr != 0) h2_tdssdmcp_fr_r->Fill(tdssdmcp_fr_r,DFMAEvent[dssd_fr_subev].tid);
		  //if(tdssd_ba != 0) h2_tdssdmcp_ba_r->Fill(tdssdmcp_ba_r,DFMAEvent[dssd_ba_subev].tid-160);

  	    	cdown=DFMAEvent[i].ehi;
		h1_cdown->Fill(cdown);
 

	}



	 // MWPC DE
	 if (DFMAEvent[i].tid==PPACDE && DFMAEvent[i].ehi > 0) {

           ppacde =  DFMAEvent[i].ehi;

           ppacde_pu = DFMAEvent[i].pu;

           tdssdppacde_fr = (DFMAEvent[i].LEDts) -  (tdssd_fr);     
      
           h2_tdssdppacde_fr->Fill(tdssdppacde_fr,DFMAEvent[dssd_fr_subev].tid);

           //h2_ppacde_E->Fill(ppacde,dssd_fr_emax);
           
           h2_tdssdppacde->Fill(tdssdppacde_fr,ppacde);


	 }





        break;

       

        default:
	break;
      
 

  } 
  
}





#endif


float d;

for (i=0;i<nsubev;i++) {
  
  switch (DFMAEvent[i].tpe) {
    
         //---------\\
         // X array \\
         //---------\\

	 case XARRAY:


      #if(0)


      //Energy = ((float)(DFMAEvent[i].postrise_energy)-(float)(DFMAEvent[i].prerise_energy))/20;//*ehigain[gsid];

      //h2_cle_id_raw->Fill(Energy, DFMAEvent[i].tid);
	   //    h2_cle_id_raw->Fill(DFMAEvent[i].ehi, DFMAEvent[i].tid);
      //h1_cl_tracelen->Fill(DFMAEvent[i].traceLen);



      d = (float) (RAND_MAX) + 1.0;  
      //Energy = Energy*clmap[DFMAEvent[i].tid].gain+clmap[DFMAEvent[i].tid].off+(float)rand()/d-0.5;

        DFMAEvent[i].ehi = DFMAEvent[i].ehi*clmap[DFMAEvent[i].tid].gain+clmap[DFMAEvent[i].tid].off+(float)rand()/d-0.5;

      //DFMAEvent[i].ehi = Energy; // SAVING RAW CLOVERS TO TREE

	//     h2_cle_id->Fill(DFMAEvent[i].ehi, DFMAEvent[i].tid);

 


/*
         if(DFMAEvent[i].tid < 21){
            h2_cl_ehi_basesample[DFMAEvent[i].tid]->Fill(DFMAEvent[i].ehi,DFMAEvent[i].basesample);
            //h2_cl_ehi_baseline[DFMAEvent[i].tid]->Fill(DFMAEvent[i].ehi,(DFMAEvent[i].prerisebeg + DFMAEvent[i].preriseend)/2);
  	    h2_cl_ehi_slope[DFMAEvent[i].tid]->Fill(DFMAEvent[i].ehi,float(DFMAEvent[i].preriseend) - float(DFMAEvent[i].prerisebeg));
         }

             // Traces for clovers...

            // printf("\nDebug #1: DFMAEvent[i].tid: %i, TrIn1: %i, DFMAEvent[i].traceLen: %i, DFMAEvent[i].trace[100]: %i, Energy: %f\n",DFMAEvent[i].tid,TrIn1,DFMAEvent[i].traceLen,DFMAEvent[i].trace[100],Energy);

             if(DFMAEvent[i].tid == 17 && TrIn1 < 1000 && Energy < 21774 && Energy > 21659){
                for(kk=0;kk<DFMAEvent[i].traceLen;kk++){
                   h2_cl_E1_traces->Fill(TrIn1,kk,DFMAEvent[i].trace[kk]-DFMAEvent[i].preriseend);
                }
                TrIn1++;
             } 

             if(DFMAEvent[i].tid == 5 && TrIn2 < 1000 && Energy < 21600 && Energy > 21500){
                for(kk=0;kk<DFMAEvent[i].traceLen;kk++){
                   h2_cl_B1_traces_middle->Fill(TrIn2,kk,DFMAEvent[i].trace[kk]-DFMAEvent[i].preriseend);
                }
                TrIn2++;
             } 

             if(DFMAEvent[i].tid == 5 && TrIn3 < 1000 && Energy < 21700 && Energy > 21600){
                for(kk=0;kk<DFMAEvent[i].traceLen;kk++){
                   h2_cl_B1_traces_end->Fill(TrIn3,kk,DFMAEvent[i].trace[kk]-DFMAEvent[i].preriseend);
                }
                TrIn3++;
             } 

             if(DFMAEvent[i].tid == 2 && TrIn4 < 1000 && Energy < 20847 && Energy > 20766){
                for(kk=0;kk<DFMAEvent[i].traceLen;kk++){
                   h2_cl_A2_traces->Fill(TrIn4,kk,DFMAEvent[i].trace[kk]-DFMAEvent[i].preriseend);
                }
                TrIn4++;
             } 



         if(DFMAEvent[i].tid == 20){

              for(j=0;j<DFMAEvent[i].traceLen;j++){
                   traceAverage = traceAverage + DFMAEvent[i].trace[j];
              }

              traceAverage = traceAverage/(DFMAEvent[i].traceLen);


              //printf("\nPoint one: %i. Point two: %f. Point three: %f. Point four: %i\n",DFMAEvent[i].prerisebeg,float(DFMAEvent[i].prerisesum/300),float(DFMAEvent[i].trace[0]+DFMAEvent[i].trace[1])/2, DFMAEvent[i].preriseend);             


              //printf("Prerisebeg is : %i, preriseend is : %i\n",DFMAEvent[i].prerisebeg,DFMAEvent[i].preriseend);
              //printf("Postrisebeg is : %i, postriseend is : %i\n",DFMAEvent[i].postrisebeg,DFMAEvent[i].postriseend);
              //printf("prerisesum is %i\n",DFMAEvent[i].prerisesum);

             h2_cl_E4_ehi_baseline->Fill(DFMAEvent[i].ehi,(DFMAEvent[i].baseline/1000));
             h2_cl_E4_ehi_prerisebeg->Fill(DFMAEvent[i].ehi,DFMAEvent[i].prerisebeg);
             h2_cl_E4_ehi_preriseend->Fill(DFMAEvent[i].ehi,DFMAEvent[i].preriseend);
             h2_cl_E4_preriseend_prerisebeg->Fill(DFMAEvent[i].preriseend,DFMAEvent[i].prerisebeg);

             h2_cl_E4_preriseend_diff1->Fill(DFMAEvent[i].preriseend,float(DFMAEvent[i].preriseend - float(DFMAEvent[i].prerisesum)/300));
             h2_cl_E4_preriseend_diff2->Fill(DFMAEvent[i].preriseend,float(DFMAEvent[i].preriseend - DFMAEvent[i].prerisebeg));

           if(DFMAEvent[i].basesample > 3400 && DFMAEvent[i].basesample < 3575){

             h2_cl_E4_ehi_diff1->Fill(DFMAEvent[i].ehi,float(DFMAEvent[i].preriseend - float(DFMAEvent[i].prerisesum)/300));
             h2_cl_E4_ehi_diff2->Fill(DFMAEvent[i].ehi,float(DFMAEvent[i].preriseend - DFMAEvent[i].prerisebeg));
             h2_cl_E4_ehi_diff3->Fill(DFMAEvent[i].ehi,(traceAverage - float(DFMAEvent[i].prerisesum)/300));
             postriseend_corrected = DFMAEvent[i].postriseend - (DFMAEvent[i].ehi*(-0.0145)-14);
             diffsum = (postriseend_corrected - DFMAEvent[i].postrisebeg) + (DFMAEvent[i].preriseend - DFMAEvent[i].prerisebeg) + (traceAverage - float(DFMAEvent[i].prerisesum)/300);

	     h2_cl_E4_ehi_diff4->Fill(DFMAEvent[i].ehi,diffsum);
           }



             //h2_cl_E4_ehi_diff4->Fill(DFMAEvent[i].ehi,(DFMAEvent[i].postriseend - DFMAEvent[i].postrisebeg));

             h2_cl_E4_ehi_sampledbase->Fill(DFMAEvent[i].ehi,float(DFMAEvent[i].sampledbase)/1000);

             h2_cl_E4_ehi_trueBL->Fill(float(DFMAEvent[i].trace[0] + DFMAEvent[i].trace[1])/2);




         }


*/

         //printf("\nclover %i base sample: %i, baseline: %i\n",DFMAEvent[i].tid,DFMAEvent[i].basesample,DFMAEvent[i].baseline);



         //h2_cle_id_keV->Fill(float(DFMAEvent[i].ehi)/3.0, DFMAEvent[i].tid);

         //*****************************\\
         // Allocating clover structure \\
         //*****************************\\
    
         if(ncl < MAXNUMGE){

         clover.tid[ncl] = DFMAEvent[i].tid;
         clover.ehi[ncl] = DFMAEvent[i].ehi;
         clover.tge[ncl] = DFMAEvent[i].LEDts;
         //clover.slope[ncl] = float(DFMAEvent[i].preriseend) - float(DFMAEvent[i].prerisebeg); 



         ncl++;
         }

#endif



      break;

       

        default:
	break;
      
 

  } 
  
}

// reading XArray datas from bin_xa

clover.nge = ncl;

 double dtcldssd;
 double clts;
 double dssdts;

 clover.nge = XAng;

for(j=0;j<clover.nge;j++){

   clover.ehi[j]=XAEvent[j].ehi;
   clover.tge[j]=XAEvent[j].event_timestamp;
   clover.tid[j]=XAEvent[j].tid;
   h1_clover_en->Fill(clover.ehi[j]);
   h2_cle_id->Fill(XAEvent[j].ehi, XAEvent[j].tid);

 clts=double(XAEvent[j].event_timestamp);
   dssdts=double(tdssd_fr);

   dtcldssd=double(XAEvent[j].event_timestamp)-double(tdssd_fr);
   //  dtcldssd=clts-dssdts;

   h2_dTcldssdall->Fill(dtcldssd, XAEvent[j].tid);

}






for(j=0;j<ncl;j++){
   h1_clover_en->Fill(clover.ehi[j]);
}


//<><><><><><><><><><><><><><><><><><>\\
//  Make deltaT between DGS and DFMA  \\
//<><><><><><><><><><><><><><><><><><>\\

double dTgdssd;
dTgdssd = 0.0;

double dTgdssda;
dTgdssda = 0.0;

double dTgdssdp;
dTgdssdp = 0.0;

double dTgdssdr;
dTgdssdr = 0.0;

int firstge = 1000;

for(i=0;i<ng;i++){

   if(DGSEvent[i].tpe == GE && firstge == 1000) firstge = i;

}


int ch40subev = 1000;


for(j=0;j<ng;j++){
    if(DGSEvent[j].tid == 40 && DGSEvent[j].tpe == GE && DGSEvent[j].flag == 0){
             ch40subev = j;
    }
}



if(ch40subev < 1000){
   for(i=0;i<nsubev;i++){
     if(DFMAEvent[i].tpe == DSSD){
        dTgdssd = double(DGSEvent[ch40subev].event_timestamp) - double(DFMAEvent[i].LEDts);
        h2_dTgdssd_ch40->Fill(dTgdssd,DFMAEvent[i].tid);
     }
   }
}

dTgdssd = 0.0;

for(i=0;i<nsubev;i++){
   if((DFMAEvent[i].LEDts > 0) && (DFMAEvent[i].tpe == DSSD) && (ng > 0)){

        dTgdssd = double(DGSEvent[firstge].event_timestamp) - double(DFMAEvent[i].LEDts);
        h2_dTgdssd->Fill(dTgdssd,DFMAEvent[i].tid);

        if(dTgdssd > 1000 && Pars.CurEvNo < 10000){
 
           //printf("\n\nDebugging!");
           //printf("\nEvent number %i",Pars.CurEvNo);
           //printf("\n DFMAEvent[%i].LEDts: %llu",i,DFMAEvent[i].LEDts);
           //printf("\n DGSEvent[%i].LEDts: %llu",firstge,DGSEvent[firstge].event_timestamp);

        } 


     }
}


dTgdssd = 0.0;


if(dssd_fr_emax >12000){

   for(i=0;i<ng;i++){
     if(DGSEvent[i].tpe == GE && DGSEvent[i].flag == 0){
        dTgdssd = double(DGSEvent[i].event_timestamp) - double(DFMAEvent[dssd_fr_subev].LEDts);
        h2_dTgdssd_allge->Fill(dTgdssd,DFMAEvent[dssd_fr_subev].tid);
     }
   }

}





double dTgl;
dTgl = 0.0;
double dTgr;
dTgr = 0.0;

if(DFMAEvent[left_subev].LEDts > 0){
   dTgl = double(DGSEvent[0].event_timestamp) - double(DFMAEvent[left_subev].LEDts);
   h1_dTgl->Fill(dTgl);
} 
if(DFMAEvent[right_subev].LEDts > 0){
   dTgr = double(DGSEvent[0].event_timestamp) - double(DFMAEvent[right_subev].LEDts);
   h1_dTgr->Fill(dTgr);
} 


//<><><>\\
// MWPC \\
//<><><>\\

h1_esi->Fill(esi);
//h1_ppacde->Fill(ppacde);

/*
if(cl != 0){
       if(TestTr1 < 100){
          for(i=0;i<1000;i++){
             h2_lefttraces->Fill(TestTr1,i,DFMAEvent[left_subev].trace[i]);
            
          }
          TestTr1++;
       }
}

if(cr != 0){
       if(TestTr4 < 100){
          for(i=0;i<1000;i++){
             h2_righttraces->Fill(TestTr4,i,DFMAEvent[right_subev].trace[i]);
            
          }
          TestTr4++;
       }
}
*/


//if(cl != 0 & cr != 0){
  
/*
   if (cr>0) crat=cl/cr*1000.0;
   if ((cl>0)||(cr>0)) cratn=cl/(cl+cr)*1000.0;

   if (cr_int>0) crat_int=cl_int/cr_int*1000.0;	//new
   if ((cl_int>0)||(cr_int>0)) cratn_int=cl_int/(cl_int+cr_int)*1000.0;	//new

   h1_cl_int->Fill(cl_int);// new
   h1_cr_int->Fill(cr_int);// new

   h1_cl->Fill(cl);
   h1_cr->Fill(cr);
   h2_clr->Fill(cl/10,cr/10);
   h2_clr_int->Fill(cl_int/50,cr_int/50); //new

   h1_cxn->Fill(cratn);

*/

  h1_ppacde->Fill(ppacde);


  h1_cl_pu->Fill(cl_pu);
  h1_cr_pu->Fill(cr_pu);
  h1_ppacde_pu->Fill(ppacde_pu);

  h2_clr->Fill(cl,cr);

  if(cl > 0 && cr > 0) { 

    csum = cl + cr;
    crat = cl - cr;
 
    

    h1_cx->Fill(crat);
    h2_xde->Fill(crat,ppacde);

   h1_cl->Fill(cl);
   h1_cr->Fill(cr);

 

   h1_csum->Fill(csum);






   //  } // cl and cr greater than 0 zero

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
 



  
   if(ndssd >0){

      h1_ppacdeg->Fill(ppacde);


      h1_clg->Fill(cl);
      h1_crg->Fill(cr);
      //h2_clrg->Fill(cl,cr);

      if(dssd_fr_emax > 4500 && dssd_ba_emax > 4500){// Gating on Alpha energies

	//h2_clrg_en->Fill(cl,cr);

         if(map_fr[DFMAEvent[dssd_fr_subev].tid].phystrip >= 20 && map_fr[DFMAEvent[dssd_fr_subev].tid].phystrip < 40){
	   //h2_clrg_en_g1->Fill(cl,cr);
            h1_cx_en_g1->Fill(cl/cr*1000.0);
            h1_cxn_en_g1->Fill(cl/(cl+cr)*1000.0);
         }

         if(map_fr[DFMAEvent[dssd_fr_subev].tid].phystrip >= 100 && map_fr[DFMAEvent[dssd_fr_subev].tid].phystrip < 120){
	   //h2_clrg_en_g2->Fill(cl,cr);
            h1_cx_en_g2->Fill(cl/cr*1000.0);
            h1_cxn_en_g2->Fill(cl/(cl+cr)*1000.0);
         }

         h2_stripfr_cxng->Fill(map_fr[DFMAEvent[dssd_fr_subev].tid].phystrip,cratn);
         h2_stripba_cxng->Fill(map_ba[DFMAEvent[dssd_ba_subev].tid].phystrip,cratn);

         if(map_fr[DFMAEvent[dssd_fr_subev].tid].phystrip > 74 && map_fr[DFMAEvent[dssd_fr_subev].tid].phystrip < 86){

            h2_stripba_cxng_frontgated->Fill(map_ba[DFMAEvent[dssd_ba_subev].tid].phystrip,cratn);

         }



      }

      //h2_clrg_int->Fill(cl_int/50,cr_int/50);//new

      h1_cxg->Fill(crat);
      h1_cxng->Fill(cratn);
 
      h1_cxg_int->Fill(crat_int);//new
      h1_cxng_int->Fill(cratn_int);//new

   } // n dssd more than zero

 } // cl and  cr not zero

  // printf("\nTest2, event number %i\n",Pars.CurEvNo);


//********************************\\
// Now decide if recoil or decay! \\
//********************************\\

int r_fr_emax, r_ba_emax;
int s_r_fr, s_r_ba;
int r_fr_subev, r_ba_subev;
int r_fr_PU, r_ba_PU;
int r_fr_trace_len, r_ba_trace_len;
int r_fr_d2t0, r_ba_d2t0;
int r_fr_d2t1, r_ba_d2t1;
int r_fr_d2t2, r_ba_d2t2;
int r_fr_trace[1000], r_ba_trace[1000];
unsigned long long int r_fr_ts, r_ba_ts;

int d_fr_emax, d_ba_emax;
int s_d_fr, s_d_ba;
int s_d_fr_phys, s_d_ba_phys;
int d_fr_subev, d_ba_subev;
int d_fr_PU, d_ba_PU;
int d_fr_trace_len, d_ba_trace_len;
int d_fr_d2t0, d_ba_d2t0;
int d_fr_d2t1, d_ba_d2t1;
int d_fr_d2t2, d_ba_d2t2;
int d_fr_trace[1000], d_ba_trace[1000];
unsigned long long int d_fr_ts, d_ba_ts;

int rec_pu;
int fdecEv;
int s_fd_fr;
int s_fd_ba;

int fd_fr_emax, fd_ba_emax;
unsigned long long int fd_fr_ts, fd_fr_ts2;
unsigned long long int fd_ba_ts;
int fd_fr_emax2;



fd_fr_ts = 0;
fd_fr_ts2;
fd_fr_emax2 = 0;
fd_fr_emax = 0;
fd_ba_emax = 0;

fd_ba_ts = 0;

r_fr_emax = 0;
r_ba_emax = 0;
s_r_fr = 0;
s_r_ba = 0;
r_fr_subev = 0;
r_ba_subev = 0;
r_fr_PU = 0;
r_ba_PU = 0;
r_fr_trace_len = 0;
r_ba_trace_len = 0;
r_fr_d2t0 = 0;
r_ba_d2t0 = 0;
r_fr_d2t1 = 0;
r_ba_d2t1 = 0;
r_fr_d2t2 = 0;
r_ba_d2t2 = 0;
for(i=0;i<DSSDTRLEN;i++){
  r_fr_trace[i] = 0;
  r_ba_trace[i] = 0;
  d_fr_trace[i] = 0;
  d_ba_trace[i] = 0;
} 
r_fr_ts = 0;
r_ba_ts = 0;

d_fr_emax = 0;
d_ba_emax = 0;
s_d_fr = 0;
s_d_ba = 0;
s_d_fr_phys = 0;
s_d_ba_phys = 0;
d_fr_subev = 0;
d_ba_subev = 0;
d_fr_PU = 0;
d_ba_PU = 0;
d_fr_trace_len = 0;
d_ba_trace_len = 0;
d_fr_d2t0 = 0;
d_ba_d2t0 = 0;
d_fr_d2t1 = 0;
d_ba_d2t1 = 0;
d_fr_d2t2 = 0;
d_ba_d2t2 = 0;
d_fr_ts = 0;
d_ba_ts = 0;

rec_pu = 0;
fdecEv = 0;
s_fd_fr = 0;
s_fd_ba = 0;



#if(1)

   short int pud2[1000];
     
     int dec_PU, rec_PU;

     int t0, t1;
     int d2t0thr, d2t1thr;


     //THRESHOLDS

     d2t0thr=500;
     d2t1thr=500;

     int mm,kk,nn;
     float trise;

     mm=100;
     kk=15;
     nn=0;
     trise=5;


     float e0,e1,ef0,ef1,eb0,eb1;

      int kdd = 10; 
      int d2dd= 5;

int dtmin = 190; // time gate for coincidence beteween GS-dssd
int dtmax = 230;  

if((dssd_fr_emax != 0) && (dssd_ba_emax != 0)){

   //**********************\\
   // Condition for recoil \\
   //**********************\\

       h2_eftof->Fill(dssd_fr_emax,tdssdmcp_fr_r);


      //if(((cl != 0) || (cr != 0)) && (dssd_fr_emax > 7000) && (dssd_fr_emax < 30000) && (tdssdmcp_fr>365) && (tdssdmcp_fr<374) ){

      // RECOIL CONDITION FOR THE HOSTOGRAM BELOW

      if(((cl != 0) || (cr != 0)) && (eftof->IsInside(dssd_fr_emax,tdssdmcp_fr)) ){

	//if(((cl != 0) || (cr != 0))) {

        h2_leftdtam->Fill(ts_left/10000000000-t_first/10000000000,ts_left/10000-ts_wheel/10000);
      }

  //if(((cl != 0) || (cr != 0)) && (dssd_fr_emax > 6000) && (dssd_fr_emax < 30000) && (tdssdmcp_fr>365) && (tdssdmcp_fr<376) ){

      // ACTUAL RECOIL CONDITION

      int tr,tl;

      tr = int(tdssdmcp_fr_r);
      tl=int(tdssdmcp_fr);

      if(  
	 //	 ((cl != 0) || (cr != 0) || (ppacde!= 0)) && (dssd_fr_emax > 10000) && (dssd_ba_emax > 10000) && (dssd_fr_emax < 50000) && (dssd_ba_emax < 50000)
	 ((cl != 0) || (cr != 0) || (ppacde!= 0)) 
	 //&& (dssd_fr_emax > 1000) && (dssd_ba_emax > 1000) && (dssd_fr_emax < 50000) && (dssd_ba_emax < 50000) 
         &&(ppacde>16000)
	 &&(cl-cr<14000)
	 // && (((tr<-16) && (tr>-30))||((tr<34) && (tr>25)))
	 //&& (tr<34) && (tr>25)
	 //((cl != 0) || (cr != 0) || (ppacde!= 0)) && (dssd_ba_emax > 10000) && (dssd_ba_emax < 25000) &&
	 //(eftof->IsInside(dssd_fr_emax,tdssdmcp_fr)) 
	 ){

	 // && (eftof->IsInside(dssd_fr_emax,tdssdmcp_fr))
	 //	 (dssd_fr_emax > 10000) && (dssd_ba_emax > 10000) && (dssd_fr_emax < 50000) && (dssd_ba_emax < 50000)//217Pa
	 //	 		 (dssd_fr_emax > 17000) && (dssd_ba_emax > 17000) && (dssd_fr_emax < 33000) && (dssd_ba_emax < 33000) //254Rf
		 //	 (((dssd_fr_emax > 16000) || (dssd_ba_emax > 16000))&&((dssd_fr_emax < 24000) && (dssd_ba_emax < 24000)))
         // && (dssd_fr_emax < 25000)

	//if(((cl != 0) || (cr != 0))) {
	//h2_eftofg->Fill(dssd_fr_emax,tdssdmcp_fr);

  //if(((cl != 0) || (cr != 0)) && (ppacde > 8000)) {

     //&&(dssd_fr_emax >= 15000)&&(dssd_ba_emax >= 15000)&&(ppacde<75000) ) {
     //if((cl == 0) || (cr == 0) || (ppacde == 0)) {

    //if ((dssd_fr_emax >= 20000) && (dssd_ba_emax >= 20000) &&
    //    (dssd_fr_emax <= 40000) && (dssd_ba_emax <= 40000) )
    //{

     h2_xdeg->Fill(crat,ppacde);
     h2_dssd_hitxy_physg->Fill(map_ba[DFMAEvent[dssd_ba_subev].tid].phystrip,161-map_fr[DFMAEvent[dssd_fr_subev].tid].phystrip);

      //FRONT
      r_fr_emax = dssd_fr_emax;
      s_r_fr = DFMAEvent[dssd_fr_subev].tid;
      r_fr_ts = DFMAEvent[dssd_fr_subev].LEDts;	// Take recoil timestamp from FRONT...
      r_fr_subev = dssd_fr_subev;               
      r_fr_PU = DFMAEvent[dssd_fr_subev].pu;
      r_fr_trace_len = DFMAEvent[dssd_fr_subev].traceLen;
      r_fr_d2t0 = DFMAEvent[dssd_fr_subev].d2t0;
      r_fr_d2t1 = DFMAEvent[dssd_fr_subev].d2t1;
      r_fr_d2t2 = DFMAEvent[dssd_fr_subev].d2t2;
      for (j=0;j<DFMAEvent[dssd_fr_subev].traceLen;j++) r_fr_trace[j]= DFMAEvent[dssd_fr_subev].trace[j];
      int r_fr_basesample; 
      r_fr_basesample=DFMAEvent[dssd_fr_subev].basesample; 
      
      h2_eftofg->Fill(r_fr_emax,tdssdmcp_fr_r);

      //BACK
      r_ba_emax = dssd_ba_emax;
      s_r_ba = DFMAEvent[dssd_ba_subev].tid - NUMFR;
      r_ba_ts = DFMAEvent[dssd_ba_subev].LEDts;	// Take recoil timestamp from BACK...
      r_ba_subev = dssd_ba_subev;
      r_ba_PU = DFMAEvent[dssd_ba_subev].pu;
      r_ba_trace_len = DFMAEvent[dssd_ba_subev].traceLen;
      r_ba_d2t0 = DFMAEvent[dssd_ba_subev].d2t0;
      r_ba_d2t1 = DFMAEvent[dssd_ba_subev].d2t1;
      r_ba_d2t2 = DFMAEvent[dssd_ba_subev].d2t2;
      for (j=0;j<DFMAEvent[dssd_ba_subev].traceLen;j++) r_ba_trace[j]= DFMAEvent[dssd_ba_subev].trace[j];
      int r_ba_basesample; 
      r_ba_basesample=DFMAEvent[dssd_ba_subev].basesample;

      // gammas gated by recoils



       /*      if(ng > 0){
	      for(j=0;j<ng;j++) {
                dTgdssd = DGSEvent[i].event_timestamp - r_fr_ts;
                //h2_dTgdssd_corr->Fill(dTgdssd,dssd_corr[s_d_fr][s_d_ba].chain.recoil.getid[j]);
                if(dTgdssd > dtmin && dTgdssd < dtmax) {
                  h2_rec_gammas->Fill(DGSEvent[j].ehi,DGSEvent[j].tid);
		}
	      }
            }	  */ 
 
               //****************************\\
               // DEALING WITH RECOIL PILEUP \\
               //****************************\\

               //****************************\\
               // TWO-SIDED CORRELATIONS...
               //****************************\\

     if ((s_r_fr!=0)&&(s_r_ba!=0)) {


     // front side

      // flip negative traces 

      for (i=0;i<r_fr_trace_len;i++) { 
	r_fr_trace[i] = r_fr_trace[i] - r_fr_basesample+100;
      }

      // fix the spikes
      
      for (i=50;i<r_fr_trace_len;i++) { 
	//if (r_fr_trace[i] < r_fr_trace[i-1]/2 ) r_fr_trace[i] = r_fr_trace[i-1];
      }   


     //printf("FRONT DDFILTER\n");



     //ddFilterNEW(r_fr_trace,r_fr_trace_len,pud2);
      ddFilterNEW2(r_fr_trace,r_fr_trace_len,kdd,d2dd,pud2);


     t0=0;
     t1=0;
          
     for(i=30;i<(r_fr_trace_len-5);i++){
   
       if ((pud2[i]>=d2t1thr)&&(pud2[i-1]<d2t1thr)&&(t0!=0)&&(t1==0)) {
           t1=i;
       }

       if ((pud2[i]>=d2t0thr)&&(pud2[i-1]<d2t0thr)&&(t0==0)) {
           t0=i;
       }

     }


     if (t1!=0) { 

     //printf("FRONT T1 NOT ZERO\n");

     e0=0;
     e1=0;

     GetPUEn(t1,t0,r_fr_trace,r_fr_trace_len,e0,e1,mm,kk,nn,trise);

     //printf("e0=%f e1=%f\n",e0,e1);


 


     //e0 = e0/10;
     //e1 = e1/10;
     e0 = e0/5;
     e1 = e1/5;
     e0 = map_fr[s_r_fr].gain*(e0 + float(rand())/RAND_MAX-0.5) + map_fr[s_r_fr].off;     
     e1 = map_fr[s_r_fr].gain*(e1 + float(rand())/RAND_MAX-0.5) + map_fr[s_r_fr].off;   
     ef0=e0; ef1=e1;

     //printf("ef0=%f ef1=%f t1=%i t0=%i\n",ef0,ef1,t1,t0);

     //d_fr_emax=ef0;
     fd_fr_emax=ef1;

     fd_fr_ts=r_fr_ts+t1-t0;
      
     r_fr_PU=2;

     }

     t0=0;
     t1=0;

     // back side

      for (i=0;i<r_ba_trace_len;i++) { 
            r_ba_trace[i] = -r_ba_trace[i];
            r_ba_trace[i] = r_ba_trace[i] + r_ba_basesample+100;
      }

     // fix the spikes
      
      for (i=50;i<r_ba_trace_len;i++) { 
	if (r_ba_trace[i] < r_ba_trace[i-1]/2 ) r_ba_trace[i] = r_ba_trace[i-1];
      }   
      
     //printf("BACK DDFILTER\n");


     //ddFilterNEW(r_ba_trace,r_ba_trace_len,pud2);
      ddFilterNEW2(r_ba_trace,r_ba_trace_len,kdd,d2dd,pud2);
/*
          if (trn1<1000) {

              for (jj=0;jj<1000;jj++) h2_traceg->Fill(trn1,jj,d_fr_trace[jj]);
              trn1=trn1+1;
          }

        if (trn2<1000) {

              for (jj=0;jj<1000;jj++) h2_traceg2->Fill(trn1,jj,pud2[jj]);
              trn2=trn2+1;
          }
*/

     for(i=30;i<(r_ba_trace_len-5);i++){

       if ((pud2[i]>=d2t1thr)&&(pud2[i-1]<d2t1thr)&&(t0!=0)&&(t1==0)) {
           t1=i;
       }
       if ((pud2[i]>=d2t0thr)&&(pud2[i-1]<d2t0thr)&&(t0==0)) {
           t0=i;
       }
     }


     if (t1!=0) { 

     //printf("BACK T1 NOT ZERO\n");



     e0=0;
     e1=0;

     GetPUEn(t1,t0,r_ba_trace,r_ba_trace_len,e0,e1,mm,kk,nn,trise);

     //printf("e0=%f e1=%f\n",e0,e1);

     //e0 = e0/10;
     //e1 = e1/10;
     e0 = e0/5;
     e1 = e1/5;
     //printf("e0=%f e1=%f\n",e0,e1);

     e0 = map_ba[s_r_ba+NUMFR].gain*(e0 + float(rand())/RAND_MAX-0.5) + map_ba[s_r_ba+NUMFR].off;     
     e1 = map_ba[s_r_ba+NUMFR].gain*(e1 + float(rand())/RAND_MAX-0.5) + map_ba[s_r_ba+NUMFR].off;  

     //printf("e0=%f e1=%f\n",e0,e1);

     eb0=e0; eb1=e1;

     //printf("eb0=%f eb1=%f t1=%i t0=%i\n",eb0,eb1,t1,t0);

     //d_ba_emax=eb0;
     fd_ba_emax=eb1;

     fd_ba_ts=r_ba_ts+t1-t0;
      
     r_ba_PU=2;

     }


       int diff;

       if ((r_fr_PU > 1)&&(r_ba_PU > 1)) diff= fd_ba_ts-fd_fr_ts;


       if ((r_fr_PU > 1)&&(r_ba_PU > 1)&&(diff<5)&&(diff>-5)) {
       //if ((r_fr_PU > 1)&&(r_ba_PU > 1)) {

 
       //printf("DIFF =%d\n", diff);
      
       rec_PU=1;
       fdecEv=1;

       s_fd_fr=s_r_fr;
       s_fd_ba=s_r_ba;

       /*
    if (trn1<1000) {

      //for (jj=0;jj<1000;jj++) h2_traceg->Fill(trn1,jj,r_fr_trace[jj]);
      //h2_traceg->SetBinContent(trn1,999,s_fd_fr);
      //trn1=trn1+1;
          }

        if (trn2<1000) {

	  //for (jj=0;jj<1000;jj++) h2_traceg2->Fill(trn2,jj,pud2[jj]);
	  //for (jj=0;jj<1000;jj++) h2_traceg2->Fill(trn2,jj,r_ba_trace[jj]);
	  //h2_traceg2->SetBinContent(trn2,999,s_fd_ba);
	  //trn2=trn2+1;
          }
       */
       //printf("FDECEV\n");

     }

   } // s_r_fr and s_r_ba


               /*

	       if ((r_fr_PU > 1)&&(r_ba_PU > 1)) rec_PU = 1;
               if ((rec_PU!=0)&&(s_r_fr!=0)&&(s_r_ba!=0)) {


               // true PU condition    
               //if (((curt1s[s_r_fr]-curt1s[s_r_ba+NUMFR])<5)&&((curt1s[s_r_fr]-curt1s[s_r_ba+NUMFR])>-5)&&(curt1s[s_r_fr]>0)&&(curt1s[s_r_ba+NUMFR]>0)) {

		  fdecEv = 1;
                  if((r_fr_PU == 3) && (r_ba_PU == 3)) fdecEv = 3;
                  

                  s_fd_fr=s_r_fr;
                  s_fd_ba=s_r_ba;
      
                  t0 = curt0s[s_r_fr];
                  t1 = curt1s[s_r_fr];
                  t2 = curt2s[s_r_fr];


                  GetPUEn(t1,t0,r_fr_trace,r_fr_trace_len,e0,e1,m,kk,nn,trise);
                  

                  e0 = e0/100;
                  e1 = e1/100;

                  e0 = map_fr[s_r_fr].gain*(e0 + float(rand())/RAND_MAX-0.5) + map_fr[s_r_fr].off;     
                  e1 = map_fr[s_r_fr].gain*(e1 + float(rand())/RAND_MAX-0.5) + map_fr[s_r_fr].off;  

                  ef0=e0; ef1=e1;

                  r_fr_emax=ef0;
                  fd_fr_emax=ef1;

                  fd_fr_ts=r_fr_ts+curt1s[s_r_fr]-curt0s[s_r_fr];

                  if(fdecEv == 3){

                     GetPUEn(t2,t1,r_fr_trace,r_fr_trace_len,e1_b,e2,m,kk,nn,trise);
                     
                     e2 = e2/100;

                     e2 = map_fr[s_r_fr].gain*(e2 + float(rand())/RAND_MAX-0.5) + map_fr[s_r_fr].off;  
                     ef2 = e2;
                     
                     fd_fr_emax2 = ef2; 
                     fd_fr_ts2 = r_fr_ts + curt2s[s_r_fr] - curt0s[s_r_fr];

                    
                  }


                  // back energies
                  t0 = curt0s[s_r_ba+NUMFR];
                  t1 = curt1s[s_r_ba+NUMFR];
                  t2 = curt2s[s_r_ba+NUMFR];

                  GetPUEn(t1,t0,r_ba_trace,r_ba_trace_len,e0,e1,m,kk,nn,trise);

   
                  e0 = e0/10;
                  e1 = e1/10;
                  e0 = map_ba[s_r_ba+NUMFR].gain*(e0 + float(rand())/RAND_MAX-0.5) + map_ba[s_r_ba+NUMFR].off;     
                  e1 = map_ba[s_r_ba+NUMFR].gain*(e1 + float(rand())/RAND_MAX-0.5) + map_ba[s_r_ba+NUMFR].off; 
                  eb0=e0; eb1=e1;
    
                  r_ba_emax=eb0;
                  fd_ba_emax=eb1;
                  fd_ba_ts=r_ba_ts+curt1s[s_r_ba+NUMFR]-curt0s[s_r_ba+NUMFR];


                  if(fdecEv == 3){

                     GetPUEn(t2,t1,r_ba_trace,r_ba_trace_len,e1_b,e2,m,kk,nn,trise);
                     
                     e2 = e2/100;

                     e2 = map_ba[s_r_ba+NUMFR].gain*(e2 + float(rand())/RAND_MAX-0.5) + map_ba[s_r_ba+NUMFR].off;  
                     ef2 = e2;
                     
                     fd_ba_emax2 = ef2; 
                     fd_ba_ts2 = r_ba_ts + curt2s[s_r_ba+NUMFR] - curt0s[s_r_ba+NUMFR];

                    
                  }

     
		} */
    
        //}
}
  /*
if(n_sib>0){

  for (i=0;i<n_sib;i++) {
     

    dssdfrsi_t = double(DFMAEvent[dssd_fr_subev].LEDts) - double(sib_ts[i]);

    dssdbasi_t = double(DFMAEvent[dssd_ba_subev].LEDts) - double(sib_ts[i]);

  h2_dtdssdfrsi->Fill(dssdfrsi_t, DFMAEvent[dssd_fr_subev].tid);
  h2_dtdssdbasi->Fill(dssdbasi_t, DFMAEvent[dssd_ba_subev].tid-160);

  h2_sidssd->Fill(sib_e[i],d_fr_emax);

  printf("ind box mat %d %f\n",sib_tid[i],sib_e[i]);

  if (n_sib>1) {
  for (j=i+1;j<n_sib;j++) {
    h2_sisi->Fill(sib_e[i],sib_e[j]);
    h2_sisi->Fill(sib_e[j],sib_e[i]);
    }
  }

  }

}
    
  */  

   //*************************\\
   // Condition for decays..  \\
   //*************************\\

      if(
	 (cl == 0) && (cr == 0) && (ppacde ==0) 
	 //	 && (dssd_fr_emax < 220000) && (dssd_ba_emax < 220000)
	 //		 		 (dssd_fr_emax > 5000) && (dssd_ba_emax > 5000)&&
	 //	  (dssd_fr_emax < 10000) && (dssd_ba_emax < 10000)
		 //	      	 && ((dssd_fr_emax > 110000)||(dssd_ba_emax>110000))
         ) {

   //if( (dssd_fr_emax < 12000) ){

   //if( (dssd_fr_emax < 10000) && (dssd_fr_emax > 100) ){

      //FRONT
      d_fr_emax = dssd_fr_emax;
      s_d_fr = DFMAEvent[dssd_fr_subev].tid;
      s_d_fr_phys = map_fr[DFMAEvent[dssd_fr_subev].tid].phystrip;
      d_fr_subev = dssd_fr_subev;
      d_fr_ts = DFMAEvent[dssd_fr_subev].LEDts;
      d_fr_PU = DFMAEvent[dssd_fr_subev].pu;
      d_fr_trace_len = DFMAEvent[dssd_fr_subev].traceLen;
      h1_trlen->Fill(d_fr_trace_len);
      for(jj=0;jj<DFMAEvent[dssd_fr_subev].traceLen;jj++) d_fr_trace[jj]=DFMAEvent[dssd_fr_subev].trace[jj];
      int d_fr_basesample; 
      d_fr_basesample=DFMAEvent[dssd_fr_subev].basesample;
 
     	      
      /*WORK*/
 
        long long int dts;

	//if (d_fr_emax>2000) {
        
        int basesample;


        //basesample=DFMAEvent[dssd_fr_subev].basesample;
      
        h1_header_type->Fill(DFMAEvent->header_type);

        //basesample=DFMAEvent[dssd_fr_subev].prerisesum/200;


        //if (DFMAEvent[dssd_fr_subev].tid==40) {

        h1_basesample->Fill(DFMAEvent[dssd_fr_subev].basesample);
        h2_basesample->Fill(DFMAEvent[dssd_fr_subev].basesample,s_d_fr);
 
        h2_prerisebeg->Fill(DFMAEvent[dssd_fr_subev].prerisebeg,s_d_fr);
        h2_postrisebeg->Fill(DFMAEvent[dssd_fr_subev].postrisebeg,s_d_fr);
        h2_preriseend->Fill(DFMAEvent[dssd_fr_subev].preriseend,s_d_fr);
        h2_postriseend->Fill(DFMAEvent[dssd_fr_subev].postriseend,s_d_fr);
        h2_baseline->Fill(DFMAEvent[dssd_fr_subev].baseline,s_d_fr);
 
        h2_pre2->Fill(DFMAEvent[dssd_fr_subev].prerisebeg,DFMAEvent[dssd_fr_subev].preriseend);
        //h2_post2->Fill(DFMAEvent[dssd_fr_subev].postrisebeg,DFMAEvent[dssd_fr_subev].postriseend);
        h2_post2->Fill(DFMAEvent[dssd_fr_subev].postrisebeg,DFMAEvent[dssd_fr_subev].m2_last_beg);
        //h2_last2->Fill(DFMAEvent[dssd_fr_subev].m2_last_beg,DFMAEvent[dssd_fr_subev].m2_last_end);
        h2_last2->Fill(DFMAEvent[dssd_fr_subev].prerisebeg,DFMAEvent[dssd_fr_subev].m2_last_end);

        //}

        dts=DFMAEvent[dssd_fr_subev].LEDts-DFMAEvent[dssd_fr_subev].prevTS;
        if (d_fr_basesample_first[s_d_fr]==1)  { 
         d_fr_basesample_av[s_d_fr]=DFMAEvent[dssd_fr_subev].basesample; 
         d_fr_basesample_first[s_d_fr]++;
        } 
        else {
	  if (dts>20000) d_fr_basesample_av[s_d_fr]=0.99*d_fr_basesample_av[s_d_fr]+0.01*DFMAEvent[dssd_fr_subev].basesample;
         d_fr_basesample_first[s_d_fr]++;
	}
	if ((d_fr_emax>3000)&&(d_fr_emax<3300))  h2_e_bs->Fill(basesample+200-d_fr_basesample_av[s_d_fr],d_fr_emax-3000+s_d_fr*300);
         h2_e_dt->Fill(dts,d_fr_emax);
         if (d_fr_basesample_first[s_d_fr]>100) h2_blav->Fill(d_fr_basesample_av[s_d_fr],s_d_fr);	

	 //} // d_fr_max


      //BACK
      d_ba_emax = dssd_ba_emax;
      s_d_ba = DFMAEvent[dssd_ba_subev].tid-NUMFR;    
      s_d_ba_phys = map_ba[DFMAEvent[dssd_ba_subev].tid].phystrip; 
      d_ba_ts = DFMAEvent[dssd_ba_subev].LEDts;		// Take decay timestamp from BACK...
      d_ba_subev = dssd_ba_subev;
      d_ba_PU = DFMAEvent[dssd_ba_subev].pu;
      d_ba_trace_len = DFMAEvent[dssd_ba_subev].traceLen;
      for(jj=0;jj<DFMAEvent[dssd_ba_subev].traceLen;jj++) d_ba_trace[jj]=DFMAEvent[dssd_ba_subev].trace[jj]; 
      int d_ba_basesample; 
      d_ba_basesample=DFMAEvent[dssd_ba_subev].basesample;
   
         // Si box

      //double dssdfrsi_t, dssdbasi_t;


if(n_sib>0){

printf("n_sib %d\n",n_sib);

//std::array<double, 4> get_dead_layer_corrections(double dssd_E, int dssd_strip_x, int dssd_strip_y, double box_E, int box_wall, int box_strip, int box_detector);
 
  float sibox_emax1, sibox_emax2;
  int s_box1, s_box2;

  // find 2 box strips with the highest energy

  /*

  s_box1=-1; s_box2=-1; sibox_emax1=0; sibox_emax2=0;
  
  for (i=0;i<n_sib;i++) {
    if (sib_e[i]>sibox_emax1) {  
      sibox_emax2=sibox_emax1; 
      s_box2=s_box1; 
      sibox_emax1=sib_e[i]; 
      s_box1=i; } 
    else { 
      if (sib_e[i]>sibox_emax2)  { 
         sibox_emax2=sib_e[i]; 
         s_box2=i; 
      };
    }
  }


  float ae1, ae2, aedssd;
  float aebox1, aebox2;
  float dedssda1, dedssddla1, deboxdla1; 
  float dedssda2, dedssddla2, deboxdla2; 
  float atot, atotcor;
 

  // energy compression

  aedssd=d_fr_emax/1000.0;
  sibox_emax1=sibox_emax1/1000.0;
  sibox_emax2=sibox_emax2/1000.0;

  int dssdx, dssdy;

  dssdx=161-map_fr[DFMAEvent[dssd_fr_subev].tid].phystrip;
  dssdy=161-map_ba[DFMAEvent[dssd_ba_subev].tid].phystrip;

  // single escape condition

  if ((aedssd>0.5)&&(aedssd<5.0)&&(sibox_emax1>0.5)&&(sibox_emax2<0.1)) {
   
   printf("single_escape dssde=%f x=%d y=%d sieb=%f wn=%d sn=%d detn=%d \n", aedssd, dssdx, dssdy, sibox_emax1, sib_wn[s_box1], sib_stripn[s_box1], sib_detn[s_box1]); 

   int siw1, sis1, sid1;

   siw1=sib_wn[s_box1]; sis1= sib_stripn[s_box1]; sid1=sib_detn[s_box1];
    
   std::array<double, 4> dead_layer_corrections =  get_dead_layer_corrections(aedssd, dssdx, dssdy, sibox_emax1, siw1, sis1, sid1);

    ae1=dead_layer_corrections[0]*1000.0;
    dedssddla1=dead_layer_corrections[1]*1000.0;
    deboxdla1=dead_layer_corrections[2]*1000.0;

    //printf("ae1=%f dedssddla1=%f deboxdla1=%f \n", ae1, dedssddla1, deboxdla1); 
    printf("ae1=%f dedssddla1=%f deboxdla1=%f \n", dead_layer_corrections[0], dead_layer_corrections[1], dead_layer_corrections[2]); 

    float  dela1tot; 
    dela1tot=dedssddla1+deboxdla1;
    aebox1= dedssddla1+deboxdla1+sibox_emax1*1000.0;
    atotcor=aebox1+aedssd*1000.0;

    h2_sidssdcor->Fill(aebox1,aedssd*1000.0);

    h1_atotcor->Fill(atotcor);
    
    if ((sid1+siw1*2)==3) h2_ada1->Fill(aedssd*1000+sibox_emax1*1000,dela1tot);
    if ((sid1+siw1*2)==4) h2_ada2->Fill(aedssd*1000+sibox_emax1*1000,dela1tot);
    if ((sid1+siw1*2)==5) h2_ada3->Fill(aedssd*1000+sibox_emax1*1000,dela1tot);
    if ((sid1+siw1*2)==6) h2_ada4->Fill(aedssd*1000+sibox_emax1*1000,dela1tot);
    if ((sid1+siw1*2)==7) h2_ada5->Fill(aedssd*1000+sibox_emax1*1000,dela1tot);
    if ((sid1+siw1*2)==8) h2_ada6->Fill(aedssd*1000+sibox_emax1*1000,dela1tot);
    if ((sid1+siw1*2)==9) h2_ada7->Fill(aedssd*1000+sibox_emax1*1000,dela1tot);
    if ((sid1+siw1*2)==10) h2_ada8->Fill(aedssd*1000+sibox_emax1*1000,dela1tot);

  }

  // single escape condition with one alpha going in

  if ((aedssd>5.0)&&(aedssd<10.0)&&(sibox_emax1>0.5)&&(sibox_emax2<0.1)) {

    printf("single_escape_one_nonescape dssde=%f x=%d y=%d sieb=%f wn=%d sn=%d detn=%d \n", aedssd, dssdx, dssdy, sibox_emax1, sib_wn[s_box1], sib_stripn[s_box1], sib_detn[s_box1]); 

    std::array<double, 4> se =  single_escape_one_nonescape(aedssd, dssdx, dssdy, sibox_emax1, sib_wn[s_box1], sib_stripn[s_box1], sib_detn[s_box1]);

    dedssda1=se[0]*1000.0;
    dedssda2=se[1]*1000.0;
    dedssddla2=se[2]*1000.0;
    deboxdla2=se[3]*1000.0;

    printf("dedssda1=%f dedssda2=%f  dedssddla2=%f deboxdla2=%f \n", dedssda1, dedssda2, dedssddla2, deboxdla2); 

    ae1=dedssda1;
    ae2= dedssda2+dedssddla2+deboxdla2+sibox_emax1*1000.0;
    aebox2= dedssddla2+deboxdla2+sibox_emax1*1000.0;

    h2_sidssdcor2->Fill(dedssddla2,ae2);
    h2_aa1esc->Fill(ae1, ae2);

  }

  // double escape condition

  if ((aedssd>0.5)&&(aedssd<10.0)&&(sibox_emax1>0.5)&&(sibox_emax2>0.5)) {

   printf("double_escape dssde=%f x=%d y=%d sieb1=%f wn1=%d sn1=%d detn1=%d sieb2=%f wn2=%d sn2=%d detn2=%d \n", aedssd, dssdx, dssdy, sibox_emax1,  
	  sib_wn[s_box1], sib_stripn[s_box1], sib_detn[s_box1], sibox_emax2, sib_wn[s_box2], sib_stripn[s_box2], sib_detn[s_box2]); 

   std::array<double, 6> de = double_escape(aedssd, dssdx, dssdy, sibox_emax1, sibox_emax2, sib_stripn[s_box1], sib_detn[s_box1],  
                                                         sib_e[s_box2], sib_wn[s_box2], sib_stripn[s_box2], sib_detn[s_box2]);
    
  dedssda1=de[0]*1000.0;
  dedssda2=de[1]*1000.0;
  dedssddla1=de[2]*1000.0;
  deboxdla1=de[3]*1000.0;
  dedssddla2=de[4]*1000.0;
  deboxdla2=de[5]*1000.0;  

  printf("dedssda1=%f dedssda2=%f  dedssddla1=%f deboxdla1=%f dedssddla2=%f deboxdla2=%f \n", dedssda1, dedssda2, dedssddla1, deboxdla1, dedssddla2, deboxdla2); 

  ae1=dedssda1+dedssddla1+deboxdla1+sibox_emax1*1000.0;
  ae2=dedssda2+dedssddla2+deboxdla2+sibox_emax2*1000.0;

  aebox1=dedssddla1+deboxdla1+sibox_emax1*1000.0;
  aebox2=dedssddla2+deboxdla2+sibox_emax2*1000.0;

  h2_sisicor2esc1->Fill(dedssda1, aebox1);
  h2_sisicor2esc2->Fill(dedssda2, aebox2);
  h2_aa2esc->Fill(ae1, ae2);

  }

  */

 // si-dssd and si-si histograms

  float atot;

  for (i=0;i<n_sib;i++) {
     
    dssdfrsi_t = double(DFMAEvent[dssd_fr_subev].LEDts) - double(sib_ts[i]);
    dssdbasi_t = double(DFMAEvent[dssd_ba_subev].LEDts) - double(sib_ts[i]);

    if (sib_e[i]>1000) h2_dtdssdfrsi->Fill(dssdfrsi_t, DFMAEvent[dssd_fr_subev].tid);
    if (sib_e[i]>1000) h2_dtdssdbasi->Fill(dssdbasi_t, DFMAEvent[dssd_ba_subev].tid-160);

    atot=sib_e[i]+d_fr_emax;

    h1_atot->Fill(atot);
    h2_sidssd->Fill(sib_e[i],d_fr_emax);

    //printf("ind box mat %d %f\n",sib_tid[i],sib_e[i]);

    if (n_sib>1) {
    for (j=i+1;j<n_sib;j++) {
      h2_sisi->Fill(sib_e[i],sib_e[j]);
      h2_sisi->Fill(sib_e[j],sib_e[i]);
    }
    }

  } // si box


}


            
               //***************************\\
               // DEALING WITH DECAY PILEUP \\
 
               //***************************\\
 
  
     if ((s_d_fr!=0)&&(s_d_ba!=0)) {


     // front side

      // flip negative traces 

      for (i=0;i<d_fr_trace_len;i++) { 
            d_fr_trace[i] = d_fr_trace[i] - d_fr_basesample+100;
      }
      


     //printf("FRONT DDFILTER\n");



      //ddFilterNEW2(d_fr_trace,d_fr_trace_len, pud2);
      ddFilterNEW2(d_fr_trace,d_fr_trace_len,kdd,d2dd,pud2);



     t0=0;
     t1=0;
          
     for(i=30;i<(d_fr_trace_len-5);i++){
   
       if ((pud2[i]>=d2t1thr)&&(pud2[i-1]<d2t1thr)&&(t0!=0)&&(t1==0)) {
           t1=i;
       }

       if ((pud2[i]>=d2t0thr)&&(pud2[i-1]<d2t0thr)&&(t0==0)) {
           t0=i;
       }

     }



     if (t1!=0) { 

     //printf("FRONT T1 NOT ZERO\n");

     e0=0;
     e1=0;

     GetPUEn(t1,t0,d_fr_trace,d_fr_trace_len,e0,e1,mm,kk,nn,trise);

     //printf("e0=%f e1=%f\n",e0,e1);

     //e0 = e0/10;
     //e1 = e1/10;
     e0 = e0/5;
     e1 = e1/5;
     e0 = map_fr[s_d_fr].gain*(e0 + float(rand())/RAND_MAX-0.5) + map_fr[s_d_fr].off;     
     e1 = map_fr[s_d_fr].gain*(e1 + float(rand())/RAND_MAX-0.5) + map_fr[s_d_fr].off;   
     ef0=e0; ef1=e1;

     //printf("ef0=%f ef1=%f t1=%i t0=%i\n",ef0,ef1,t1,t0);

     //d_fr_emax=ef0;
     fd_fr_emax=ef1;

     fd_fr_ts=d_fr_ts+t1-t0;
      


     d_fr_PU=2;

     }

     t0=0;
     t1=0;

     // back side

      for (i=0;i<d_ba_trace_len;i++) { 
            d_ba_trace[i] = -d_ba_trace[i];
            d_ba_trace[i] = d_ba_trace[i] + d_ba_basesample+100;
      }
      
     //printf("BACK DDFILTER\n");

     //ddFilterNEW(d_ba_trace,d_ba_trace_len,pud2);
     ddFilterNEW2(d_ba_trace,d_ba_trace_len,kdd,d2dd,pud2);

     for(i=30;i<(d_ba_trace_len-5);i++){

       if ((pud2[i]>=d2t1thr)&&(pud2[i-1]<d2t1thr)&&(t0!=0)&&(t1==0)) {
           t1=i;
       }
       if ((pud2[i]>=d2t0thr)&&(pud2[i-1]<d2t0thr)&&(t0==0)) {
           t0=i;
       }
     }


     if (t1!=0) { 

     //printf("BACK T1 NOT ZERO\n");



     //e0=0;
     //e1=0;
     e0 = e0/5;
     e1 = e1/5;

     GetPUEn(t1,t0,d_ba_trace,d_ba_trace_len,e0,e1,mm,kk,nn,trise);

     //printf("e0=%f e1=%f\n",e0,e1);

     e0 = e0/10;
     e1 = e1/10;

     //printf("e0=%f e1=%f\n",e0,e1);

     e0 = map_ba[s_d_ba+NUMFR].gain*(e0 + float(rand())/RAND_MAX-0.5) + map_ba[s_d_ba+NUMFR].off;     
     e1 = map_ba[s_d_ba+NUMFR].gain*(e1 + float(rand())/RAND_MAX-0.5) + map_ba[s_d_ba+NUMFR].off;  

     //printf("e0=%f e1=%f\n",e0,e1);

     eb0=e0; eb1=e1;

     //printf("eb0=%f eb1=%f t1=%i t0=%i\n",eb0,eb1,t1,t0);

     //d_ba_emax=eb0;
     fd_ba_emax=eb1;

     fd_ba_ts=d_ba_ts+t1-t0;
      
     d_ba_PU=2;

     }


     if ((d_fr_PU > 1)&&(d_ba_PU > 1)&&(fd_ba_ts-fd_fr_ts<5)&&(fd_ba_ts-fd_fr_ts>-5)) {
       dec_PU=1;
       fdecEv=2;

       s_fd_fr=s_d_fr;
       s_fd_ba=s_d_ba;

       //printf("FDECEV\n");

     }

   } // s_d_fr and s_d_ba

/*

               int curt0s[1000];
               int curt1s[1000];
               
               if ((d_fr_PU > 1)&&(d_ba_PU > 1))  dec_PU = 1; 

               if ((dec_PU!=0)&&(s_d_fr!=0)&&(s_d_ba!=0)) {

               // true PU condition 
               if (((curt1s[s_d_fr]-curt1s[s_d_ba+NUMFR])<5)&&((curt1s[s_d_fr]-curt1s[s_d_ba+NUMFR])>-5)&&(curt1s[s_d_fr]>0)&&(curt1s[s_d_ba+NUMFR]>0)) {

                  fdecEv=2;

                  s_fd_fr=s_d_fr;
                  s_fd_ba=s_d_ba;


                  t0 = curt0s[s_d_fr];
                  t1 = curt1s[s_d_fr];

                  GetPUEn(t1,t0,d_fr_trace,d_fr_trace_len,e0,e1,m,kk,n,trise);


                  e0 = e0/100;
                  e1 = e1/100;
                  e0 = map_fr[s_d_fr].gain*(e0 + float(rand())/RAND_MAX-0.5) + map_fr[s_d_fr].off;     
                  e1 = map_fr[s_d_fr].gain*(e1 + float(rand())/RAND_MAX-0.5) + map_fr[s_d_fr].off;   
                  ef0=e0; ef1=e1;

                  d_fr_emax=ef0;
                  fd_fr_emax=ef1;

                  fd_fr_ts=d_fr_ts+curt1s[s_d_fr]-curt0s[s_d_fr];


                 //******************
                 // back energies


                 t0 = curt0s[s_d_ba+NUMFR];
                 t1 = curt1s[s_d_ba+NUMFR];

                 GetPUEn(t1,t0,d_ba_trace,d_ba_trace_len,e0,e1,m,kk,n,trise);

                 e0 = e0/10;
                 e1 = e1/10;
                 e0 = map_ba[s_d_ba+NUMFR].gain*(e0 + float(rand())/RAND_MAX-0.5) + map_ba[s_d_ba+NUMFR].off;     
                 e1 = map_ba[s_d_ba+NUMFR].gain*(e1 + float(rand())/RAND_MAX-0.5) + map_ba[s_d_ba+NUMFR].off;   
                 eb0=e0; eb1=e1;

                 d_ba_emax=eb0;
                 fd_ba_emax=eb1;
                 fd_ba_ts=d_ba_ts+curt1s[s_d_ba+NUMFR]-curt0s[s_d_ba+NUMFR];



		 }}
  */     
	//}

 
   } // DECAY condition

#endif

//<><><><><><><><><><><><><><><><>\\
// DSSD hit patterns and energies \\
//<><><><><><><><><><><><><><><><>\\

if(s_r_fr != 0 && s_r_ba != 0) h2_r_hitxy->Fill(s_r_fr,s_r_ba);
if(s_d_fr != 0 && s_d_ba != 0) h2_d_hitxy->Fill(s_d_fr,s_d_ba);

if(dssd_fr_emax > 0 && dssd_ba_emax > 0){
  //h2_dssd_hitxy_phys->Fill(map_fr[DFMAEvent[dssd_fr_subev].tid].phystrip,map_ba[DFMAEvent[dssd_ba_subev].tid].phystrip);

h2_dssd_hitxy_phys->Fill(map_ba[DFMAEvent[dssd_ba_subev].tid].phystrip,161-map_fr[DFMAEvent[dssd_fr_subev].tid].phystrip);


 h2_dssd_fr_emax_phy->Fill(dssd_fr_emax, 161-map_fr[DFMAEvent[dssd_fr_subev].tid].phystrip);
 h2_dssd_ba_emax_phy->Fill(dssd_ba_emax, map_ba[DFMAEvent[dssd_ba_subev].tid].phystrip);

   if(TestTr3 < 100){
      for(i=0;i<DSSDTRLEN;i++){
         h2_dssd_traces_fr->Fill(TestTr3,i,DFMAEvent[dssd_fr_subev].trace[i]);
         h2_dssd_traces_ba->Fill(TestTr3,i,DFMAEvent[dssd_ba_subev].trace[i]);
      }
      TestTr3++;
   }
}

if(dssd_fr_emax > 0){

   h2_dssd_fr_emax->Fill(dssd_fr_emax,DFMAEvent[dssd_fr_subev].tid);
   h2_d_fr_e->Fill(d_fr_emax,s_d_fr);
   h2_r_fr_e->Fill(r_fr_emax,s_r_fr);

   //   h2_ppacde_E->Fill(d_fr_emax,ppacde);

   //if((DFMAEvent[dssd_fr_subev].LEDts - t_first)/1E5 < 100000){
         h1_dssd_rate_1D->Fill((DFMAEvent[dssd_fr_subev].LEDts - t_first)/ratecomp);
         //h1_dssd_rate_1us->Fill((DFMAEvent[dssd_fr_subev].LEDts - t_first)/100);
        if(s_d_fr != 0) h1_decay_rate->Fill((DFMAEvent[dssd_fr_subev].LEDts - t_first)/1E8);

        if(s_r_fr != 0) h1_recoil_rate->Fill((DFMAEvent[dssd_fr_subev].LEDts - t_first)/1E8);
   
	//}

}
}

if(dssd_ba_emax > 0){

   h2_dssd_ba_emax->Fill(dssd_ba_emax,DFMAEvent[dssd_ba_subev].tid-NUMFR);
   h2_d_ba_e->Fill(d_ba_emax,s_d_ba);
   h2_r_ba_e->Fill(r_ba_emax,s_r_ba);

}






/*
if(Pars.CurEvNo > 1925000){
printf("\n\nEvent number: %i, r_fr_subev: %i, s_r_fr: %i, s_r_ba : %i, ng: %i\n",Pars.CurEvNo,r_fr_subev,s_r_fr,s_r_ba,ng);
printf("status: %i\n",dssd_corr[s_r_fr][s_r_ba].status);
printf("s_d_fr: %i, s_d_ba: %i\n",s_d_fr,s_d_ba);

}
*/


#if(1) // comment out correlations code

/*
 if (s_r_fr==0) {
   s_r_fr=161;
   r_fr_emax=r_ba_emax;
 };

 if (s_r_ba==0) {
   s_r_ba=161;
   r_ba_emax=r_fr_emax;
 };
*/

//  s_r_fr=1;
//  r_fr_emax=r_ba_emax;
//  r_fr_ts=r_ba_ts;

//   s_r_ba=1;
//   r_ba_emax=r_fr_emax;
//   r_ba_ts=r_fr_ts;


//<><><><><><><><><><><><><><><><><><><>\\
//				        \\
// ******** Correlations Code ********* \\
//				        \\
//<><><><><><><><><><><><><><><><><><><>\\

//<><><><><><><><><><>\\
// Recoil correlation \\
//<><><><><><><><><><>\\

int decaycond;

decaycond=0;

if((s_r_fr != 0) && (s_r_ba != 0)){

   //if(Pars.CurEvNo > 1925000)printf("EvNum %i Are we in this loop?\n",Pars.CurEvNo);

   switch (dssd_corr[s_r_fr][s_r_ba].status) {

   case 0: 

   // no recoils so far

 recoil.pu = DFMAEvent[r_fr_subev].pu;
   recoil.en = r_fr_emax;
   recoil.enb = r_ba_emax;
   recoil.ts = r_fr_ts;
   recoil.d2t2 = r_fr_d2t2;
   //   recoil.d2t2 = tdssdmcp_fr_r;
   recoil.d2t0 = r_fr_d2t0;
   recoil.d2t1 = r_fr_d2t1;
   //recoil.fdecEv = fdecEv;
   recoil.x = crat;
   //   recoil.right=ppacde;
   recoil.ptofl= tdssdmcp_fr;
   recoil.ptofr= tdssdmcp_fr_r;
   recoil.pdetof= tdssdppacde_fr;
   recoil.pde=ppacde;
   //recoil.nge = ng;
   //for(i=0;i<ng;i++){
     //recoil.geehi[i] = clover.ehi[i];
     //recoil.tge[i] = clover.tge[i];
     //recoil.getid[i] = clover.tid[i];
   //}

   for(i=0;i<DSSDTRLEN;i++){
      recoil.trace_ba[i] = r_ba_trace[i];
      recoil.trace_fr[i] = r_fr_trace[i];
   }




   #if(WITHDGS)
   // Storing DGS info for this recoil...
   //recoil.nge = ng;
   jj = 0;
   for(i=0;i<ng;i++){
      if(DGSEvent[i].tpe == GE &&  DGSEvent[i].flag == 0 && jj < MAXNUMGE){
          recoil.geehi[jj] = DGSEvent[i].ehi;
          recoil.tge[jj] = DGSEvent[i].event_timestamp;
          recoil.getid[jj] = DGSEvent[i].tid;
          jj++;
          h1_wheel_ge->Fill(wheel_left);
      }
   }
  recoil.nge = jj;

  h1_nge->Fill(recoil.nge);

   #endif


   dssd_corr[s_r_fr][s_r_ba].status = 1;
   dssd_corr[s_r_fr][s_r_ba].chain.recoil = recoil;
   //   dssd_corr[s_r_fr][s_r_ba].chain.s_fr=s_r_fr;
   //   dssd_corr[s_r_fr][s_r_ba].chain.s_ba=s_r_ba;
 dssd_corr[s_r_fr][s_r_ba].chain.s_fr=161-map_fr[DFMAEvent[dssd_fr_subev].tid].phystrip;
 dssd_corr[s_r_fr][s_r_ba].chain.s_ba=map_ba[DFMAEvent[dssd_ba_subev].tid].phystrip;
 

   // if(Pars.CurEvNo > 1925800)printf("\nTest middle event number %i\n",Pars.CurEvNo);


   break;

   case 1:

   // recoil last

   //if ((r_fr_ts-dssd_corr[s_r_fr][s_r_ba].chain.recoil.ts)<10000000) { // decay time condition, 10 ms
   //decaycond=1;
   //}
   //else
   //{
   recoil.pu = DFMAEvent[r_fr_subev].pu;
   recoil.en = r_fr_emax;
   recoil.enb = r_ba_emax;
   recoil.ts = r_fr_ts;
   recoil.d2t2 = r_fr_d2t2;
   //   recoil.d2t2 = tdssdmcp_fr_r;
   recoil.d2t0 = r_fr_d2t0;
   recoil.d2t1 = r_fr_d2t1;
   //recoil.fdecEv = fdecEv;
   recoil.x = crat;
   //   recoil.right=ppacde;
   recoil.ptofl= tdssdmcp_fr;
   recoil.ptofr= tdssdmcp_fr_r;
   recoil.pdetof= tdssdppacde_fr;
   recoil.pde=ppacde;
   //recoil.nge = ng;
   //for(i=0;i<ng;i++){
     //recoil.geehi[i] = clover.ehi[i];
     //recoil.tge[i] = clover.tge[i];
     //recoil.getid[i] = clover.tid[i];
   //}

   for(i=0;i<DSSDTRLEN;i++){
      recoil.trace_ba[i] = r_ba_trace[i];
      recoil.trace_fr[i] = r_fr_trace[i];

      //printf("N %i INT %i SHORT %hi\n",i, r_fr_trace[i], recoil.trace_fr[i]);

   }




   #if(WITHDGS)
   // Storing DGS info for this recoil...

   jj = 0;
   for(i=0;i<ng;i++){
      if(DGSEvent[i].tpe == GE &&  DGSEvent[i].flag == 0 && jj < MAXNUMGE){
          recoil.geehi[jj] = DGSEvent[i].ehi;
          recoil.tge[jj] = DGSEvent[i].event_timestamp;
          recoil.getid[jj] = DGSEvent[i].tid;
          jj++;
      }
   }
  recoil.nge = jj;  

   #endif


   dssd_corr[s_r_fr][s_r_ba].status = 1;
   dssd_corr[s_r_fr][s_r_ba].chain.recoil = recoil;
   //   dssd_corr[s_r_fr][s_r_ba].chain.s_fr=s_r_fr;
   //   dssd_corr[s_r_fr][s_r_ba].chain.s_ba=s_r_ba;
 dssd_corr[s_r_fr][s_r_ba].chain.s_fr=161-map_fr[DFMAEvent[dssd_fr_subev].tid].phystrip;
 dssd_corr[s_r_fr][s_r_ba].chain.s_ba=map_ba[DFMAEvent[dssd_ba_subev].tid].phystrip;

   //} // decay time condition

   // if(Pars.CurEvNo > 1925800)printf("\nTest middle event number %i\n",Pars.CurEvNo);

   break;


   case 2:
  

     //if ((r_fr_ts-dssd_corr[s_r_fr][s_r_ba].chain.decay[0].ts)<10000000) { // decay time condition, 10 ms
     //decaycond=1;
     //}
     //else
     //{

        // store previous chain in the tree

        chain=dssd_corr[s_r_fr][s_r_ba].chain;     
        chain.corr_type = 1;
        tree->Fill();

	// reset the pixel

        recoil.pu = DFMAEvent[r_fr_subev].pu;
        recoil.en = r_fr_emax;
	recoil.enb = r_ba_emax;
        recoil.ts = r_fr_ts;
        recoil.d2t2 = r_fr_d2t2;
	//        recoil.d2t2 = tdssdmcp_fr_r;
        recoil.d2t0 = r_fr_d2t0;
        recoil.d2t1 = r_fr_d2t1;
        recoil.x = crat;
//recoil.ptof;
	//        recoil.right=ppacde;
	recoil.ptofl= tdssdmcp_fr;
	recoil.ptofr= tdssdmcp_fr_r;
        recoil.pdetof= tdssdppacde_fr;
	recoil.pde=ppacde;
        //recoil.fdecEv = fdecEv;

  	#if(WITHDGS)
   	// Storing DGS info for this recoil...

        jj = 0;
        for(i=0;i<ng;i++){
          if(DGSEvent[i].tpe == GE &&  DGSEvent[i].flag == 0 && jj < MAXNUMGE){
            recoil.geehi[jj] = DGSEvent[i].ehi;
            recoil.tge[jj] = DGSEvent[i].event_timestamp;
            recoil.getid[jj] = DGSEvent[i].tid;
            jj++;
        }
        }

        recoil.nge = jj;

   	#endif

        for(i=0;i<DSSDTRLEN;i++){
          recoil.trace_ba[i] = r_ba_trace[i];
          recoil.trace_fr[i] = r_fr_trace[i];
        }


  	dssd_corr[s_r_fr][s_r_ba].chain.recoil=recoil;


        // reinitialising decays...
        for (i=0; i<dssd_corr[s_r_fr][s_r_ba].chain.ndec; i++) {
         if(i<MAXNUMDEC){
          dssd_corr[s_r_fr][s_r_ba].chain.decay[i].ts=0;
          dssd_corr[s_r_fr][s_r_ba].chain.decay[i].en=0;
          dssd_corr[s_r_fr][s_r_ba].chain.decay[i].enb=0;
          dssd_corr[s_r_fr][s_r_ba].chain.decay[i].pu_fr=0;
          dssd_corr[s_r_fr][s_r_ba].chain.decay[i].pu_ba=0;
          dssd_corr[s_r_fr][s_r_ba].chain.decay[i].time=0;
          for(j=0;j<dssd_corr[s_r_fr][s_r_ba].chain.decay[i].traceLen;j++){
             dssd_corr[s_r_fr][s_r_ba].chain.decay[i].trace_fr[j]=0;
             dssd_corr[s_r_fr][s_r_ba].chain.decay[i].trace_ba[j]=0;
           }
          dssd_corr[s_r_fr][s_r_ba].chain.decay[i].traceLen=0;
	
	
        dssd_corr[s_r_fr][s_r_ba].chain.decay[i].nge = 0;
  	for(j=0;j<MAXNUMGE;j++){
         
   	  dssd_corr[s_r_fr][s_r_ba].chain.decay[i].geehi[j] = 0;
    	  dssd_corr[s_r_fr][s_r_ba].chain.decay[i].tge[j] = 0;
    	  dssd_corr[s_r_fr][s_r_ba].chain.decay[i].getid[j] = 0;
         
   	}

        dssd_corr[s_r_fr][s_r_ba].chain.decay[i].nsi = 0;
  	for(j=0;j<MAXNUMSI;j++){
         
   	  dssd_corr[s_r_fr][s_r_ba].chain.decay[i].siehi[j] = 0;
    	  dssd_corr[s_r_fr][s_r_ba].chain.decay[i].tsi[j] = 0;
    	  dssd_corr[s_r_fr][s_r_ba].chain.decay[i].sitid[j] = 0;
         
   	}

        }}

        dssd_corr[s_r_fr][s_r_ba].chain.ndec=0;
        dssd_corr[s_r_fr][s_r_ba].status=1;
	//       dssd_corr[s_r_fr][s_r_ba].chain.s_fr=s_r_fr;
	//     dssd_corr[s_r_fr][s_r_ba].chain.s_ba=s_r_ba;
 dssd_corr[s_r_fr][s_r_ba].chain.s_fr=161-map_fr[DFMAEvent[dssd_fr_subev].tid].phystrip;
 dssd_corr[s_r_fr][s_r_ba].chain.s_ba=map_ba[DFMAEvent[dssd_ba_subev].tid].phystrip;
	//} // decay time condition

    //} 

    break;

   default: // there is chain already

       // store previous chain in the tree

        chain=dssd_corr[s_r_fr][s_r_ba].chain;     
        chain.corr_type = 1;
        tree->Fill();

	// reset the pixel

        recoil.pu = DFMAEvent[r_fr_subev].pu;
        recoil.en = r_fr_emax;
	recoil.enb = r_ba_emax;
        recoil.ts = r_fr_ts;
        recoil.d2t2 = r_fr_d2t2;
	//      recoil.d2t2 = tdssdmcp_fr_r;
        recoil.d2t0 = r_fr_d2t0;
        recoil.d2t1 = r_fr_d2t1;
        recoil.x = crat;
	//        recoil.right=ppacde;
	recoil.ptofl= tdssdmcp_fr;
	recoil.ptofr= tdssdmcp_fr_r;
        recoil.pdetof= tdssdppacde_fr;
	recoil.pde=ppacde;
        //recoil.fdecEv = fdecEv;

  	#if(WITHDGS)
   	// Storing DGS info for this recoil...
   	recoil.nge = ng;
   	jj = 0;
   	for(i=0;i<ng;i++){
      	   if(DGSEvent[i].tpe == GE && jj < MAXNUMGE){
              recoil.geehi[jj] = DGSEvent[i].ehi;
              recoil.tge[jj] = DGSEvent[i].event_timestamp;
              recoil.getid[jj] = DGSEvent[i].tid;
              jj++;
           }
   	}
        recoil.nge = jj;
   	#endif

        for(i=0;i<DSSDTRLEN;i++){
          recoil.trace_ba[i] = r_ba_trace[i];
          recoil.trace_fr[i] = r_fr_trace[i];
        }


  	dssd_corr[s_r_fr][s_r_ba].chain.recoil=recoil;


        // reinitialising decays...
        for (i=0; i<dssd_corr[s_r_fr][s_r_ba].chain.ndec; i++) {
         if(i<MAXNUMDEC){
          dssd_corr[s_r_fr][s_r_ba].chain.decay[i].ts=0;
          dssd_corr[s_r_fr][s_r_ba].chain.decay[i].en=0;
          dssd_corr[s_r_fr][s_r_ba].chain.decay[i].enb=0;
          dssd_corr[s_r_fr][s_r_ba].chain.decay[i].pu_fr=0;
          dssd_corr[s_r_fr][s_r_ba].chain.decay[i].pu_ba=0;
          dssd_corr[s_r_fr][s_r_ba].chain.decay[i].time=0;
          for(j=0;j<dssd_corr[s_r_fr][s_r_ba].chain.decay[i].traceLen;j++){
             dssd_corr[s_r_fr][s_r_ba].chain.decay[i].trace_fr[j]=0;
             dssd_corr[s_r_fr][s_r_ba].chain.decay[i].trace_ba[j]=0;
           }
          dssd_corr[s_r_fr][s_r_ba].chain.decay[i].traceLen=0;
	
	
        dssd_corr[s_r_fr][s_r_ba].chain.decay[i].nge = 0;
  	for(j=0;j<MAXNUMGE;j++){
         
   	  dssd_corr[s_r_fr][s_r_ba].chain.decay[i].geehi[j] = 0;
    	  dssd_corr[s_r_fr][s_r_ba].chain.decay[i].tge[j] = 0;
    	  dssd_corr[s_r_fr][s_r_ba].chain.decay[i].getid[j] = 0;
         
   	}

        dssd_corr[s_r_fr][s_r_ba].chain.decay[i].nsi = 0;
  	for(j=0;j<MAXNUMSI;j++){
         
   	  dssd_corr[s_r_fr][s_r_ba].chain.decay[i].siehi[j] = 0;
    	  dssd_corr[s_r_fr][s_r_ba].chain.decay[i].tsi[j] = 0;
    	  dssd_corr[s_r_fr][s_r_ba].chain.decay[i].sitid[j] = 0;
         
   	}

        }}

        dssd_corr[s_r_fr][s_r_ba].chain.ndec=0;
        dssd_corr[s_r_fr][s_r_ba].status=1;
	//     dssd_corr[s_r_fr][s_r_ba].chain.s_fr=s_r_fr;
	//     dssd_corr[s_r_fr][s_r_ba].chain.s_ba=s_r_ba;
 dssd_corr[s_r_fr][s_r_ba].chain.s_fr=161-map_fr[DFMAEvent[dssd_fr_subev].tid].phystrip;
 dssd_corr[s_r_fr][s_r_ba].chain.s_ba=map_ba[DFMAEvent[dssd_ba_subev].tid].phystrip;

   break;


   } // switch

} // f&b




//<><><><><><><><><><>\\
// Decay correlations \\
//<><><><><><><><><><>\\

// int dtmin = -15; // 3.5 us 
// int dtmax = 10;  // 3.5 us

//int dtmin = -160; // 2.0 us
//int dtmax = -135;  // 2.0 us

//int dtmin = 0; // 2.0 us
//int dtmax = 30;  // 2.0 us 


#define T1COMP1 1.0
#define T1COMP2 1000.0
#define T1COMP3 1000000.0
#define T1COMP4 100000000.0

#define T2COMP1 1.0
#define T2COMP2 1000.0
#define T2COMP3 1000000.0

#define T3COMP1 1.0
#define T3COMP2 1000.0
#define T3COMP3 1000000.0

#define T1GATE 10000.0
#define T2GATE 100000.0
#define T3GATE 10000.0

float logtime;
logtime = .0;
int decEv = 0;
unsigned long long int decay0_time;
decay0_time = 0;

 h1_fdecEv->Fill(fdecEv);


/*
 if (s_d_fr==0) {
   s_d_fr=161;
   d_fr_emax=d_ba_emax;
 };

 if (s_d_ba==0) {
   s_d_ba=161;
   d_ba_emax=d_fr_emax;
 };
*/

//  s_d_fr=1;
//  d_fr_emax=d_ba_emax;
//  d_fr_ts=d_ba_ts;   

//  s_d_ba=1;
//  d_ba_emax=d_fr_emax;
//  d_ba_ts=d_fr_ts;


if((s_d_fr != 0) && (s_d_ba != 0) ) decEv=1;

while (((decEv==1)||(fdecEv!=0))) {

  if (decEv==1) {
 
    decay.en = d_fr_emax;
    decay.enb = d_ba_emax;
    //    decay.en = d_ba_emax;


    decay.ts = d_fr_ts;
 
    //decay.decpos = 0;
    decay.d2t0 = d_fr_d2t0;
    decay.d2t1 = d_fr_d2t1;
    decay.d2t2 = d_fr_d2t2;

    // Decay-Decay PU
    decay.pu_fr=d_fr_PU;
    decay.pu_ba=d_ba_PU;
    for (i=0; i<DSSDTRLEN; i++) {
      decay.trace_fr[i]=d_fr_trace[i];
      decay.trace_ba[i]=d_ba_trace[i];
    }
    decay.nge = XAng;
    for(j=0;j<XAng;j++){
       decay.geehi[j] = clover.ehi[j];
       decay.tge[j] = clover.tge[j];
       decay.getid[j] = clover.tid[j];
    }
   
    if (n_sib>0) {
 
    decay.nsi = n_sib; 
    if (decay.nsi>MAXNUMSI) decay.nsi=MAXNUMSI;
    for(j=0;j<decay.nsi;j++){
      //printf("reload 1 %d  %f \n",sib_tid[j],sib_e[j]);
       decay.siehi[j] = sib_e[j];
       decay.tsi[j] = sib_ts[j];
       decay.sitid[j] = sib_tid[j];
       //printf("reload 11 %d  %d \n",decay.sitid[j],decay.siehi[j]);
    }
    }


    decEv=0;


  } else {

    if (fdecEv!=0) {


     s_d_fr=s_fd_fr;
     s_d_ba=s_fd_ba;

     decay.en = fd_fr_emax;
     decay.enb = fd_ba_emax;
     decay.ts = fd_fr_ts;


     if ((fdecEv==1) || (fdecEv == 3)) { // recoil PU

       if(fdecEv == 1){


          decay.pu_fr=r_fr_PU;
          decay.pu_ba=r_ba_PU;   
          decay.d2t0 = r_fr_d2t0;
          decay.d2t1 = r_fr_d2t1;
          decay.d2t2 = r_fr_d2t2;

          for (i=0; i<DSSDTRLEN; i++) {
             decay.trace_fr[i]=r_fr_trace[i];
             decay.trace_ba[i]=r_ba_trace[i];
          }
        decay.nge = XAng;
  	for(j=0;j<XAng;j++){
   	  decay.geehi[j] = clover.ehi[j];
    	  decay.tge[j] = clover.tge[j];
    	  decay.getid[j] = clover.tid[j];
   	}
       decay.nsi = n_sib;
    if (n_sib>0) {
      //printf("reload 2 %d  %f \n",sib_tid[j],sib_e[j]);
       if (decay.nsi>MAXNUMSI) decay.nsi=MAXNUMSI;
       for(j=0;j<decay.nsi;j++){
         //printf("reload 2 %d  %f \n",sib_tid[j],sib_e[j]);
         decay.siehi[j] = sib_e[j];
         decay.tsi[j] = sib_ts[j];
         decay.sitid[j] = sib_tid[j];
	 //printf("reload 22 %d  %d \n",decay.sitid[j],decay.siehi[j]);
       }
       
    }

          fdecEv = 0;
       }
       else if(fdecEv == 3){

          decay.pu_fr=r_fr_PU;
          decay.pu_ba=r_ba_PU;   
          decay.d2t0 = r_fr_d2t0;
          decay.d2t1 = r_fr_d2t1;
          decay.d2t2 = r_fr_d2t2;

          for (i=0; i<DSSDTRLEN; i++) {
             decay.trace_fr[i]=r_fr_trace[i];
             decay.trace_ba[i]=r_ba_trace[i];
          }
        decay.nge = XAng;
  	for(j=0;j<XAng;j++){
   	  decay.geehi[j] = clover.ehi[j];
    	  decay.tge[j] = clover.tge[j];
    	  decay.getid[j] = clover.tid[j];
   	}
        decay.nsi = n_sib;
       if (decay.nsi>MAXNUMSI) decay.nsi=MAXNUMSI;
        if (n_sib>0) {
	  //printf("reload 3 %d  %f \n",sib_tid[j],sib_e[j]);
 
        for(j=0;j<decay.nsi;j++){
         decay.siehi[j] = sib_e[j];
         decay.tsi[j] = sib_ts[j];
         decay.sitid[j] = sib_tid[j];
        }
    }

          fd_fr_emax = fd_fr_emax2;
          fd_fr_ts = fd_fr_ts2;
          fdecEv = 1;

       }
     }

      if (fdecEv==2) { // decay PU
  
      //decay.decpos = 1;

      decay.pu_fr=d_fr_PU;
      decay.pu_ba=d_ba_PU;
      decay.d2t0 = d_fr_d2t0;
      decay.d2t1 = d_fr_d2t1;
      decay.d2t2 = d_fr_d2t2;

      for (i=0; i<DSSDTRLEN; i++) {
       decay.trace_fr[i]=d_fr_trace[i];
       decay.trace_ba[i]=d_ba_trace[i];
       }
        decay.nge = XAng;
  	for(j=0;j<XAng;j++){
   	  decay.geehi[j] = clover.ehi[j];
    	  decay.tge[j] = clover.tge[j];
    	  decay.getid[j] = clover.tid[j];
   	}
	decay.nsi = n_sib;
       if (decay.nsi>MAXNUMSI) decay.nsi=MAXNUMSI;
    if (n_sib>0) {
      //printf("reload 4 %d  %f \n",sib_tid[j],sib_e[j]);
       for(j=0;j<decay.nsi;j++){
	 decay.siehi[j] = (long long int) sib_e[j];
        decay.tsi[j] = sib_ts[j];
        decay.sitid[j] = sib_tid[j];
       }
  
    }

       fdecEv=0;

     }
 
    
    }
  }



  float dece_tot;

   switch(dssd_corr[s_d_fr][s_d_ba].status){

      case 0:	// No recoil before - do nothing.
         break;
      case 1:   // First decay generation

	if (decay.en>0) {

         decay.time = decay.ts - dssd_corr[s_d_fr][s_d_ba].chain.recoil.ts;
         dssd_corr[s_d_fr][s_d_ba].chain.decay[0] = decay;
         dssd_corr[s_d_fr][s_d_ba].status++;
         dssd_corr[s_d_fr][s_d_ba].chain.ndec=1;

         if (decay.time>0){
            logtime=log10(10.0*decay.time);
			h2_e1t1log->Fill(decay.en,logtime) ;
			//-464 => -458
			// 10k => 20k

            dece_tot=decay.en+sib_emax;
            if (sib_emax>1000) h2_e1t1log_box->Fill(dece_tot,logtime);
			if (sib_emax>1000) h2_e1t1log_box1->Fill(dece_tot,logtime);

		//if(logtime > 8 && logtime < 11.4 && decay.en > 8020 && decay.en < 8140){  //previous 254No cuts
			//	if(logtime > 7 && logtime < 11.4 && decay.en > 8000 && decay.en < 8600){    //255Lr
			if(logtime > 9 && logtime < 10.5 && decay.en > 6310 && decay.en < 6360){                  //251Md
			//h2_dssd_hitxy_physg251Md->Fill(map_ba[DFMAEvent[dssd_ba_subev].tid].phystrip,161-map_fr[DFMAEvent[dssd_fr_subev].tid].phystrip);
            h2_eftof_corr->Fill(dssd_corr[s_d_fr][s_d_ba].chain.recoil.en,dssd_corr[s_d_fr][s_d_ba].chain.recoil.d2t2);
		}
         }

	 /*
         if (dssd_corr[s_d_fr][s_d_ba].chain.decay[0].time<1000000000) { 
11
           h2_e1t1logf->Fill(decay.en,logtime);

	 }
	 */

        //if ((logtime<12)&&(logtime>9)&&(decay.en>8000)) { 
	if ((logtime<11.4)&&(logtime>7)&&(decay.en>8000)&&(decay.en<8600)
) { 

             h2_dssd_hitxy_phys_corr->Fill(map_ba[DFMAEvent[dssd_ba_subev].tid].phystrip,161-map_fr[DFMAEvent[dssd_fr_subev].tid].phystrip);
             h2_dssd_fb_corr->Fill(dssd_fr_emax,dssd_ba_emax);
	     //h2_ppacde_E_corr->Fill(dssd_corr[s_d_fr][s_d_ba].chain.recoil.right,dssd_corr[s_d_fr][s_d_ba].chain.recoil.en);
         }
	/*
	if ((logtime<9.0)&&(logtime>3)&&(decay.en>50000)&&(decay.en<250000)) { 

             h2_dssd_hitxy_phys_corrf->Fill(map_ba[DFMAEvent[dssd_ba_subev].tid].phystrip,161-map_fr[DFMAEvent[dssd_fr_subev].tid].phystrip);

         }
	*/

  /*      // 188Bi 6.99 MeV

	if ((logtime>6.5)&&(logtime<8.5)&&(decay.en>6950)&&(decay.en<7050)) {
             h2_dssd_hitxy_phys_corrf->Fill(map_ba[DFMAEvent[dssd_ba_subev].tid].phystrip,161-map_fr[DFMAEvent[dssd_fr_subev].tid].phystrip);
             if(dssd_corr[s_d_fr][s_d_ba].chain.recoil.nge > 0){
	      for(j=0;j<dssd_corr[s_d_fr][s_d_ba].chain.recoil.nge;j++) {
                dTgdssd = double(dssd_corr[s_d_fr][s_d_ba].chain.recoil.tge[j]) - double(dssd_corr[s_d_fr][s_d_ba].chain.recoil.ts);
                h2_dTgdssd_corr->Fill(dTgdssd,dssd_corr[s_d_fr][s_d_ba].chain.recoil.getid[j]);
                if(dTgdssd > dtmin && dTgdssd < dtmax) {
                  h2_corr_gammasf->Fill(decay.en,dssd_corr[s_d_fr][s_d_ba].chain.recoil.geehi[j]);
		}
	      }
            } 
	}

       // 188Bi 6.81 MeV

	if ((logtime>6.5)&&(logtime<10.0)&&(decay.en>6750)&&(decay.en<6900)) {
             h2_dssd_hitxy_phys_corrf->Fill(map_ba[DFMAEvent[dssd_ba_subev].tid].phystrip,161-map_fr[DFMAEvent[dssd_fr_subev].tid].phystrip);
             if(dssd_corr[s_d_fr][s_d_ba].chain.recoil.nge > 0){
	      for(j=0;j<dssd_corr[s_d_fr][s_d_ba].chain.recoil.nge;j++) {
                dTgdssd = double(dssd_corr[s_d_fr][s_d_ba].chain.recoil.tge[j]) - double(dssd_corr[s_d_fr][s_d_ba].chain.recoil.ts);
                //h2_dTgdssd_corr->Fill(dTgdssd,dssd_corr[s_d_fr][s_d_ba].chain.recoil.getid[j]);
                if(dTgdssd > dtmin && dTgdssd < dtmax) {
                  h2_corr_gammasff->Fill(decay.en,dssd_corr[s_d_fr][s_d_ba].chain.recoil.geehi[j]);
		}
	      }
            }	   
	}

      // 189Bi 6.7 MeV

	if ((logtime>7)&&(logtime<10.0)&&(decay.en>6650)&&(decay.en<6750)) {
             h2_dssd_hitxy_phys_corrf->Fill(map_ba[DFMAEvent[dssd_ba_subev].tid].phystrip,161-map_fr[DFMAEvent[dssd_fr_subev].tid].phystrip);
             if(dssd_corr[s_d_fr][s_d_ba].chain.recoil.nge > 0){
	      for(j=0;j<dssd_corr[s_d_fr][s_d_ba].chain.recoil.nge;j++) {
                dTgdssd = double(dssd_corr[s_d_fr][s_d_ba].chain.recoil.tge[j]) - double(dssd_corr[s_d_fr][s_d_ba].chain.recoil.ts);
                //h2_dTgdssd_corr->Fill(dTgdssd,dssd_corr[s_d_fr][s_d_ba].chain.recoil.getid[j]);
                if(dTgdssd > dtmin && dTgdssd < dtmax) {
                  h2_corr_gammasfff->Fill(decay.en,dssd_corr[s_d_fr][s_d_ba].chain.recoil.geehi[j]);
		}
	      }
            }	   
	}


      // 188Po 8 MeV (fast)

	if ((logtime>4)&&(logtime<6.5)&&(decay.en>7800)&&(decay.en<8050)) {
             h2_dssd_hitxy_phys_corrf->Fill(map_ba[DFMAEvent[dssd_ba_subev].tid].phystrip,161-map_fr[DFMAEvent[dssd_fr_subev].tid].phystrip);
             if(dssd_corr[s_d_fr][s_d_ba].chain.recoil.nge > 0){
	      for(j=0;j<dssd_corr[s_d_fr][s_d_ba].chain.recoil.nge;j++) {
                dTgdssd = double(dssd_corr[s_d_fr][s_d_ba].chain.recoil.tge[j]) - double(dssd_corr[s_d_fr][s_d_ba].chain.recoil.ts);
                //h2_dTgdssd_corr->Fill(dTgdssd,dssd_corr[s_d_fr][s_d_ba].chain.recoil.getid[j]);
                if(dTgdssd > dtmin && dTgdssd < dtmax) {
                  h2_corr_gammasffff->Fill(decay.en,dssd_corr[s_d_fr][s_d_ba].chain.recoil.geehi[j]);
		}
	      }
            }	   
	} 



/*	// recoil - X array (test)
	if (logtime > 4 && logtime < 13) {
	   if(XAng > 0){
	      for(j=0;j<XAng;j++) {
		dTgdssdr = double(clover.tge[j]) - double(dssd_corr[s_d_fr][s_d_ba].chain.recoil.ts);
                h2_dTgdssdr_clenergy->Fill(dTgdssdr,clover.ehi[j]);
		}
		}
	} */

	if (logtime > 0 && logtime < 10.5 ) {
             if(XAng > 0){
	      for(j=0;j<XAng;j++) {
		dTgdssd = double(clover.tge[j]) - double(decay.ts);
                  h2_dTXadssdE->Fill(decay.en,dTgdssd);
		  // if(clover.ehi[j]>10)
		    if(clover.ehi[j]>10&&dTgdssd>427&&dTgdssd<460) h2_dssdEclE->Fill(decay.en,clover.ehi[j]);
		  if(decay.en>6300&&decay.en<6500&&clover.ehi[j]>10)
		  h2_dTXaclE->Fill(clover.ehi[j],dTgdssd);

	      }
            }	   
	}

	// electron gamma at FP

	if (logtime > 4 && logtime < 10.5 ) {
             if(XAng > 0){
	      for(j=0;j<XAng;j++) {
                //dTgdssd = double(dssd_corr[s_d_fr][s_d_ba].chain.recoil.tge[j]) - double(dssd_corr[s_d_fr][s_d_ba].chain.recoil.ts);
		dTgdssd = double(clover.tge[j]) - double(decay.ts);
		if (decay.en>6000) h2_dTcldssd_corr->Fill(dTgdssd,clover.tid[j]);
		   //if(dTgdssd > 330 && dTgdssd < 360) {
                  if(dTgdssd > 310 && dTgdssd < 360) {
                  h2_corr_cl->Fill(decay.en,clover.ehi[j]);
		  }
		if(dTgdssd <310) {
                  h2_corr_cla->Fill(decay.en,clover.ehi[j]);
		  }

	      }
            }	   
	}


	// alpha gamma at FP

	if (logtime > 4 && logtime < 10.5 ) {
             if(XAng > 0){
	      for(j=0;j<XAng;j++) {
                //dTgdssd = double(dssd_corr[s_d_fr][s_d_ba].chain.recoil.tge[j]) - double(dssd_corr[s_d_fr][s_d_ba].chain.recoil.ts);
		dTgdssda = double(clover.tge[j]) - double(decay.ts);
                   h2_dTcldssda_corr->Fill(dTgdssda,clover.tid[j]);
		   //if(dTgdssd > 330 && dTgdssd < 360) {
                  if(dTgdssda > 310 && dTgdssda < 360) {
                  h2_corr_alpha_cl->Fill(decay.en,clover.ehi[j]);
		  }
		if(dTgdssda < 310) {
                  h2_corr_alpha_cla->Fill(decay.en,clover.ehi[j]);
		  }

	      }
            }	   
	}

	// proton - gamma events
	if (logtime > 3.5 && logtime < 5.5 && decay.en < 2000 && decay.en>250) {
             if(XAng > 0){
	      for(j=0;j<XAng;j++) {
		dTgdssdp = double(clover.tge[j]) - double(decay.ts);
                   h2_dTcldssdp_corr->Fill(dTgdssdp,clover.tid[j]);

		   if(dTgdssdp > -1000 && dTgdssdp < 310) {
		     h2_corr_proton_cl->Fill(decay.en,clover.ehi[j]);
		   }
		   h2_corr_proton_cl_times->Fill(clover.ehi[j], dTgdssdp);
		}
            }	   
	}



	/*

	if (logtime > 4 && logtime < 10) {
             if(dssd_corr[s_d_fr][s_d_ba].chain.recoil.nge > 0){
	      for(j=0;j<dssd_corr[s_d_fr][s_d_ba].chain.recoil.nge;j++) {
                dTgdssd = double(dssd_corr[s_d_fr][s_d_ba].chain.recoil.tge[j]) - double(dssd_corr[s_d_fr][s_d_ba].chain.recoil.ts);
                //h2_dTgdssd_corr->Fill(dTgdssd,dssd_corr[s_d_fr][s_d_ba].chain.recoil.getid[j]);
                if(dTgdssd > dtmin && dTgdssd < dtmax) {
                  h2_corr_gammasff->Fill(decay.en,dssd_corr[s_d_fr][s_d_ba].chain.recoil.geehi[j]);
		}
	      }
            }	   
	}
	*/

         h2_x_d_fr_e->Fill(dssd_corr[s_d_fr][s_d_ba].chain.recoil.x,decay.en);

         if(logtime < 5.5 && logtime > 3.5 && decay.en >250 && decay.en<2000){

            h2_x_d_fr_e_short->Fill(dssd_corr[s_d_fr][s_d_ba].chain.recoil.x,decay.en);

      
            h1_short_dec_en->Fill(decay.en);

            //if(decay.en > 0 && decay.en < 2000){
              h1_short_rec_en->Fill(dssd_corr[s_d_fr][s_d_ba].chain.recoil.en);
            //}

	      //if(decay.en > 0 ){
           for (i=0;i<decay.nsi;i++) {
  
	     //printf("ind box corr %d %d\n",decay.sitid[i], decay.siehi[i]);

            dssdfrsi_t = double(DFMAEvent[dssd_fr_subev].LEDts) - double(decay.tsi[i]);
            dssdbasi_t = double(DFMAEvent[dssd_ba_subev].LEDts) - double(decay.tsi[i]);

            h2_dtdssdfrsi_corr->Fill(dssdfrsi_t, DFMAEvent[dssd_fr_subev].tid);
            h2_dtdssdbasi_corr->Fill(dssdbasi_t, DFMAEvent[dssd_ba_subev].tid-160);

            if ((dssdfrsi_t>-20)&&(dssdfrsi_t<10)&&(dssdbasi_t>-20)&&(dssdbasi_t<10)) h2_sidssd_corr->Fill(decay.siehi[i],decay.en);

            if (decay.en<10000) {
            if (decay.nsi>1) {
	     if   ((dssdfrsi_t>-20)&&(dssdfrsi_t<10)&&(dssdbasi_t>-20)&&(dssdbasi_t<10)) {
             for (j=i+1;j<decay.nsi;j++) {
              dssdfrsi_t2 = double(DFMAEvent[dssd_fr_subev].LEDts) - double(decay.tsi[j]);
              dssdbasi_t2 = double(DFMAEvent[dssd_ba_subev].LEDts) - double(decay.tsi[j]);           
              if ((dssdfrsi_t2>-20)&&(dssdfrsi_t2<10)&&(dssdbasi_t2>-20)&&(dssdbasi_t2<10)) {   
              h2_sisi_corr->Fill(decay.siehi[i],decay.siehi[j]);
              h2_sisi_corr->Fill(decay.siehi[j],decay.siehi[i]);
	      }
             }
	    }
            }
	    }



          }
	   //}

         } 

           
         // gs alpha decay

         if(logtime < 10 && logtime > 4.0){

          //if(decay.en > 5500 && decay.en < 6500){
              h1_short_rec_bkg_en->Fill(dssd_corr[s_d_fr][s_d_ba].chain.recoil.en);
            //}

	      
          if (trn3<1000) {
	    //if(decay.en > 8400 && decay.en < 8550){
              for (jj=0;jj<DSSDTRLEN;jj++) h2_traceg3->Fill(trn3,jj,decay.trace_fr[jj]);
              trn3=trn3+1;
	      //}
          }
	 
          // ALPHA GAMMA CORRELATION MATRIX

            if(dssd_corr[s_d_fr][s_d_ba].chain.recoil.nge > 0){
	      for(j=0;j<dssd_corr[s_d_fr][s_d_ba].chain.recoil.nge;j++) {
                dTgdssd = double(dssd_corr[s_d_fr][s_d_ba].chain.recoil.tge[j]) - double(dssd_corr[s_d_fr][s_d_ba].chain.recoil.ts);
                //h2_dTgdssd_corr->Fill(dTgdssd,dssd_corr[s_d_fr][s_d_ba].chain.recoil.getid[j]);
                h1_corr_rate->Fill((double(dssd_corr[s_d_fr][s_d_ba].chain.recoil.ts)-t_first)/ratecomp);
			    h2_dTgdssd_corr_ehi->Fill(dTgdssd,dssd_corr[s_d_fr][s_d_ba].chain.recoil.geehi[j]);
                if(dTgdssd > dtmin && dTgdssd < dtmax) {
               //   h2_corr_gammas->Fill(decay.en,dssd_corr[s_d_fr][s_d_ba].chain.recoil.geehi[j]);
                  dece_tot=sib_emax+decay.en;
                  if (sib_emax>1000)  h2_corr_gammas_box->Fill(dece_tot,dssd_corr[s_d_fr][s_d_ba].chain.recoil.geehi[j]);
		}
	      }
           } 


            if(dssd_corr[s_d_fr][s_d_ba].chain.recoil.nge > 0){

	      if(decay.en > 8200 && decay.en < 8700 ) {
              // if(decay.en > 100000  ) {


        /*        for(j=0;j<dssd_corr[s_d_fr][s_d_ba].chain.recoil.nge;j++){
                     dTgdssd = double(dssd_corr[s_d_fr][s_d_ba].chain.recoil.tge[j]) - double(dssd_corr[s_d_fr][s_d_ba].chain.recoil.ts);
		     if(dssd_corr[s_d_fr][s_d_ba].chain.recoil.geehi[j]>0)
		       //h2_dTgdssd_corr->Fill(dTgdssd,dssd_corr[s_d_fr][s_d_ba].chain.recoil.getid[j]);
                     if(dTgdssd > dtmin && dTgdssd < dtmax) {
                     //  h2_corr_gammas_short->Fill(decay.en,dssd_corr[s_d_fr][s_d_ba].chain.recoil.geehi[j]);
                     //  h2_corr_gammas_id->Fill(dssd_corr[s_d_fr][s_d_ba].chain.recoil.geehi[j],dssd_corr[s_d_fr][s_d_ba].chain.recoil.getid[j]);
		     }
                 } */

                   for(j=0;j<dssd_corr[s_d_fr][s_d_ba].chain.recoil.nge-1;j++){

                     dTgdssd = double(dssd_corr[s_d_fr][s_d_ba].chain.recoil.tge[j]) - double(dssd_corr[s_d_fr][s_d_ba].chain.recoil.ts);

                     if(dTgdssd > dtmin && dTgdssd < dtmax){ 

                        for(int k = j+1;k<dssd_corr[s_d_fr][s_d_ba].chain.recoil.nge;k++){
			   dTgdssd = double(dssd_corr[s_d_fr][s_d_ba].chain.recoil.tge[k]) - double(dssd_corr[s_d_fr][s_d_ba].chain.recoil.ts);
                               if(dTgdssd > dtmin && dTgdssd < dtmax){ 

				 //h2_gamma_gamma_shortdec->Fill(dssd_corr[s_d_fr][s_d_ba].chain.recoil.geehi[j],dssd_corr[s_d_fr][s_d_ba].chain.recoil.geehi[k]);
				 //h2_gamma_gamma_shortdec->Fill(dssd_corr[s_d_fr][s_d_ba].chain.recoil.geehi[k],dssd_corr[s_d_fr][s_d_ba].chain.recoil.geehi[j]);
 
                               }
                        }


                     }
                   }
                 } // alpha energy gate 
	    } // Ge loop

         } // decay time gate

       } //decay.en>0

         //h2_e1t1c1->Fill(decay.en,decay.time/T1COMP1);
         //h2_e1t1c2->Fill(decay.en,decay.time/T1COMP2);
         //h2_e1t1c3->Fill(decay.en,decay.time/T1COMP3);
         //h2_e1t1c4->Fill(decay.en,decay.time/T1COMP4);
    
  

         break;

      case 2:   // Second decay generation

         decay.time = decay.ts - dssd_corr[s_d_fr][s_d_ba].chain.decay[0].ts;
         decay0_time = dssd_corr[s_d_fr][s_d_ba].chain.decay[0].ts - dssd_corr[s_d_fr][s_d_ba].chain.recoil.ts;
         dssd_corr[s_d_fr][s_d_ba].chain.decay[1] = decay;
         dssd_corr[s_d_fr][s_d_ba].status++;
         dssd_corr[s_d_fr][s_d_ba].chain.ndec=2;
	 
         if (decay.time>0){
            logtime=log10(10.0*decay.time);
            h2_e2t2log->Fill(decay.en,logtime);
         }
	 

         //if (dssd_corr[s_d_fr][s_d_ba].chain.decay[0].time<corr_time_short)  h2_e2t2logf->Fill(decay.en,logtime);

         /*
          if (trn2<1000) {
          if ((decay.en > 8300)&&(dssd_corr[s_d_fr][s_d_ba].chain.decay[0].en > 8400)&&(dssd_corr[s_d_fr][s_d_ba].chain.decay[0].en < 8550)) {
              for (jj=0;jj<1000;jj++) h2_traceg2->Fill(trn2,jj,decay.trace_fr[jj]);
              trn2=trn2+1;
           }
          }
	 */


         //h2_e2t2c1->Fill(decay.en,decay.time/T1COMP1);
         //h2_e2t2c2->Fill(decay.en,decay.time/T1COMP2);
         //h2_e2t2c3->Fill(decay.en,decay.time/T1COMP3);
         //h2_e2t2c4->Fill(decay.en,decay.time/T1COMP4);

/*
         if(decay.time < corr_time_short && dssd_corr[s_d_fr][s_d_ba].chain.decay[0].time  < corr_time_short)  
            h2_e1e2_short->Fill(decay.en, dssd_corr[s_d_fr][s_d_ba].chain.decay[0].en);
         if(decay.time < corr_time &&  dssd_corr[s_d_fr][s_d_ba].chain.decay[0].time < corr_time)	            
            h2_e1e2->Fill(decay.en, dssd_corr[s_d_fr][s_d_ba].chain.decay[0].en);
*/

        //  if(decay.time < corr_time_short && dssd_corr[s_d_fr][s_d_ba].chain.decay[0].time  < corr_time_short) 
	  if(decay.time < corr_time_short && dssd_corr[s_d_fr][s_d_ba].chain.decay[0].time  < corr_time_short) 
	  { 
            h2_e1e2_short->Fill(dssd_corr[s_d_fr][s_d_ba].chain.decay[0].en, decay.en);
	  }
	  if(decay.time < corr_time &&  dssd_corr[s_d_fr][s_d_ba].chain.decay[0].time < corr_time)	
	  {            
            h2_e1e2->Fill(dssd_corr[s_d_fr][s_d_ba].chain.decay[0].en, decay.en);
          }

         break;

      case 3:   // Third decay generation

         decay.time = decay.ts - dssd_corr[s_d_fr][s_d_ba].chain.decay[1].ts;
         dssd_corr[s_d_fr][s_d_ba].chain.decay[2] = decay;
         dssd_corr[s_d_fr][s_d_ba].status++;
         dssd_corr[s_d_fr][s_d_ba].chain.ndec=3;

         
         if (decay.time>0){
            logtime=log10(10.0*decay.time);
            h2_e3t3log->Fill(decay.en,logtime);
         }

         //if (dssd_corr[s_d_fr][s_d_ba].chain.decay[0].time<corr_time_short)  h2_e3t3logf->Fill(decay.en,logtime);

         //h2_e3t3c1->Fill(decay.en,decay.time/T1COMP1);
         //h2_e3t3c2->Fill(decay.en,decay.time/T1COMP2);
         //h2_e3t3c3->Fill(decay.en,decay.time/T1COMP3);
         //h2_e3t3c4->Fill(decay.en,decay.time/T1COMP4);

       //if(decay.time < corr_time_short && dssd_corr[s_d_fr][s_d_ba].chain.decay[1].time < corr_time_short)  
         // h2_e2e3_short->Fill(decay.en, dssd_corr[s_d_fr][s_d_ba].chain.decay[1].en);
       //if(decay.time < corr_time &&  dssd_corr[s_d_fr][s_d_ba].chain.decay[1].time < corr_time)	            
            h2_e2e3->Fill(decay.en, dssd_corr[s_d_fr][s_d_ba].chain.decay[1].en);
          if(  dssd_corr[s_d_fr][s_d_ba].chain.decay[1].time < corr_time_short && dssd_corr[s_d_fr][s_d_ba].chain.decay[0].time  < corr_time_short)  
            h2_e2e3f->Fill(decay.en, dssd_corr[s_d_fr][s_d_ba].chain.decay[1].en);

         break;

      case 4:   // Fourth decay generation


         dssd_corr[s_d_fr][s_d_ba].status++;
         dssd_corr[s_d_fr][s_d_ba].chain.ndec=4;
         
         decay.time = decay.ts - dssd_corr[s_d_fr][s_d_ba].chain.decay[2].ts;
         dssd_corr[s_d_fr][s_d_ba].chain.decay[3] = decay; // was [4]
         if (decay.time>0){
            logtime=log10(10.0*decay.time);
            h2_e4t4log->Fill(decay.en,logtime);
         }

         //if (dssd_corr[s_d_fr][s_d_ba].chain.decay[0].time<corr_time_short)  h2_e4t4logf->Fill(decay.en,logtime);

          h2_e3e4->Fill(decay.en, dssd_corr[s_d_fr][s_d_ba].chain.decay[2].en);
          if(  dssd_corr[s_d_fr][s_d_ba].chain.decay[1].time < corr_time_short && dssd_corr[s_d_fr][s_d_ba].chain.decay[0].time  < corr_time_short)  
            h2_e3e4f->Fill(decay.en, dssd_corr[s_d_fr][s_d_ba].chain.decay[2].en);

         /*
         h2_e4t4c1->Fill(decay.en,decay.time/T1COMP1);
         h2_e4t4c2->Fill(decay.en,decay.time/T1COMP2);
         h2_e4t4c3->Fill(decay.en,decay.time/T1COMP3);
         h2_e4t4c4->Fill(decay.en,decay.time/T1COMP4);
         */
     

         break;
      case 5:   // Fifth decay generation

         dssd_corr[s_d_fr][s_d_ba].status++;
         dssd_corr[s_d_fr][s_d_ba].chain.ndec=5;

         decay.time = decay.ts - dssd_corr[s_d_fr][s_d_ba].chain.decay[3].ts;
         dssd_corr[s_d_fr][s_d_ba].chain.decay[4] = decay; // was [5]
         if (decay.time>0){
            logtime=log10(10.0*decay.time);
            h2_e5t5log->Fill(decay.en,logtime);
         }


         //if (dssd_corr[s_d_fr][s_d_ba].chain.decay[0].time<corr_time_short)  h2_e5t5logf->Fill(decay.en,logtime);
         h2_e4e5->Fill(decay.en, dssd_corr[s_d_fr][s_d_ba].chain.decay[3].en);
          if(  dssd_corr[s_d_fr][s_d_ba].chain.decay[1].time < corr_time_short && dssd_corr[s_d_fr][s_d_ba].chain.decay[0].time  < corr_time_short)  
            h2_e4e5f->Fill(decay.en, dssd_corr[s_d_fr][s_d_ba].chain.decay[3].en);

         break;
      case 6:   // Sixth decay generation

         dssd_corr[s_d_fr][s_d_ba].status++;
         dssd_corr[s_d_fr][s_d_ba].chain.ndec=6;

         decay.time = decay.ts - dssd_corr[s_d_fr][s_d_ba].chain.decay[4].ts;
         dssd_corr[s_d_fr][s_d_ba].chain.decay[5] = decay; // was [6]
         if (decay.time>0){
            logtime=log10(10.0*decay.time);
            h2_e6t6log->Fill(decay.en,logtime);
         }

         //if (dssd_corr[s_d_fr][s_d_ba].chain.decay[0].time<corr_time_short)  h2_e6t6logf->Fill(decay.en,logtime);

         h2_e5e6->Fill(decay.en, dssd_corr[s_d_fr][s_d_ba].chain.decay[4].en);
          if(  dssd_corr[s_d_fr][s_d_ba].chain.decay[1].time < corr_time_short && dssd_corr[s_d_fr][s_d_ba].chain.decay[0].time  < corr_time_short)  
            h2_e5e6f->Fill(decay.en, dssd_corr[s_d_fr][s_d_ba].chain.decay[4].en);

         break;
      default:
         dssd_corr[s_d_fr][s_d_ba].status++;
         dssd_corr[s_d_fr][s_d_ba].chain.ndec=7;
         break;

   }

   
 } //end while


#endif



//<><><><><><><><>\\
// Singles gammas \\
//<><><><><><><><>\\

  // printf("\nTest4, event number %i\n",Pars.CurEvNo);

#if(1)

for(i=0;i<ng;i++){
  if((DGSEvent[i].tpe == GE) && (DGSEvent[i].flag == 0)){ // Check for clean GEs
//   h2_all_gammas->Fill(DGSEvent[i].ehi,DGSEvent[i].tid);


  if((DGSEvent[i].event_timestamp - t_first)/1E5 < 100000){
     h1_GE_rate->Fill((DGSEvent[i].event_timestamp - t_first)/1E5);
  }


  }
}

#endif


//<><><><><><><><><><><><><>\\
// Recoil-correlated gammas \\
//<><><><><><><><><><><><><>\\



#if(1)

//printf("\n\n1) Event number %i has s_r_fr %i, s_r_ba %i, ng %i",Pars.CurEvNo,s_r_fr,s_r_ba,ng);
//printf("\nndssd: %i, cl: %f, cr: %f\n",ndssd,cl,cr);

if(s_r_fr != 0 && s_r_ba != 0){

 for(i=0;i<ng;i++) {
  if((DGSEvent[i].tpe == GE) && (DGSEvent[i].flag == 0)){ // Check for clean GEs
   h2_xehi->Fill(crat,DGSEvent[i].ehi);
   dTgdssd = double(DGSEvent[i].event_timestamp) - double(DFMAEvent[dssd_fr_subev].LEDts);
   if(dTgdssd > dtmin && dTgdssd < dtmax){
      h2_xehig->Fill(crat,DGSEvent[i].ehi);
   }
  }
 }

}


if(ndssd > 0){



  for(i=0;i<ng;i++) {
     if((DGSEvent[i].tpe == GE) && (DGSEvent[i].flag == 0)){ // Check for clean GEs
//        h2_dssd_gammas->Fill(DGSEvent[i].ehi,DGSEvent[i].tid);
     }
  }
}


if(cl != 0 || cr != 0){

  for(i=0;i<ng;i++) {
     if((DGSEvent[i].tpe == GE) && (DGSEvent[i].flag == 0)){ // Check for clean GEs
   //     h2_mcp_gammas->Fill(DGSEvent[i].ehi,DGSEvent[i].tid);
     }
  }
}





#endif

	// recoil - X array (test)
	   if(XAng > 0){
	      for(j=0;j<XAng;j++) {
		dTgdssdr = double(clover.tge[j]) - double(dssd_corr[s_d_fr][s_d_ba].chain.recoil.ts);
                h2_dTgdssdr_clenergy->Fill(dTgdssdr,clover.ehi[j]);
		
		// IDT (test)
		if( dTgdssdr > 235 && dTgdssdr < 535){  // 0 to 3 us
		for(i=0;i<ng;i++) {
     		  if((DGSEvent[i].tpe == GE) && (DGSEvent[i].flag == 0)){
		dTgdssd = double(DGSEvent[i].event_timestamp) - double(dssd_corr[s_d_fr][s_d_ba].chain.recoil.ts);
  		 if(dTgdssd > dtmin && dTgdssd < dtmax){
		 h2_IDT->Fill(clover.ehi[j],DGSEvent[i].ehi);
		}
		}
		}

		}
		}
		}
XAng=0.0;
ng = 0.0;


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
