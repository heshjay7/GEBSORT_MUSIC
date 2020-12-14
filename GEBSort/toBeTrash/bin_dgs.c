#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <limits.h>

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

#include "GEBSort.h"
#include "GTMerge.h"
#define NGSGE NGE
#include "gsang.h"

//useful bit mask values (JTA)

const unsigned short usCLR_BIT0 = 0xFFFE;
const unsigned short usCLR_BIT1 = 0xFFFD;
const unsigned short usCLR_BIT2 = 0xFFFB;
const unsigned short usCLR_BIT3 = 0xFFF7;
const unsigned short usCLR_BIT4 = 0xFFEF;
const unsigned short usCLR_BIT5 = 0xFFDF;
const unsigned short usCLR_BIT6 = 0xFFBF;
const unsigned short usCLR_BIT7 = 0xFF7F;
const unsigned short usCLR_BIT8 = 0xFEFF;
const unsigned short usCLR_BIT9 = 0xFDFF;
const unsigned short usCLR_BIT10 = 0xFBFF;
const unsigned short usCLR_BIT11 = 0xF7FF;
const unsigned short usCLR_BIT12 = 0xEFFF;
const unsigned short usCLR_BIT13 = 0xDFFF;
const unsigned short usCLR_BIT14 = 0xBFFF;
const unsigned short usCLR_BIT15 = 0x7FFF;

const unsigned short usSET_BIT0 = 0x0001;
const unsigned short usSET_BIT1 = 0x0002;
const unsigned short usSET_BIT2 = 0x0004;
const unsigned short usSET_BIT3 = 0x0008;
const unsigned short usSET_BIT4 = 0x0010;
const unsigned short usSET_BIT5 = 0x0020;
const unsigned short usSET_BIT6 = 0x0040;
const unsigned short usSET_BIT7 = 0x0080;
const unsigned short usSET_BIT8 = 0x0100;
const unsigned short usSET_BIT9 = 0x0200;
const unsigned short usSET_BIT10 = 0x0400;
const unsigned short usSET_BIT11 = 0x0800;
const unsigned short usSET_BIT12 = 0x1000;
const unsigned short usSET_BIT13 = 0x2000;
const unsigned short usSET_BIT14 = 0x4000;
const unsigned short usSET_BIT15 = 0x8000;

const unsigned int CLR_BIT0 = 0xFFFFFFFE;
const unsigned int CLR_BIT1 = 0xFFFFFFFD;
const unsigned int CLR_BIT2 = 0xFFFFFFFB;
const unsigned int CLR_BIT3 = 0xFFFFFFF7;
const unsigned int CLR_BIT4 = 0xFFFFFFEF;
const unsigned int CLR_BIT5 = 0xFFFFFFDF;
const unsigned int CLR_BIT6 = 0xFFFFFFBF;
const unsigned int CLR_BIT7 = 0xFFFFFF7F;
const unsigned int CLR_BIT8 = 0xFFFFFEFF;
const unsigned int CLR_BIT9 = 0xFFFFFDFF;
const unsigned int CLR_BIT10 = 0xFFFFFBFF;
const unsigned int CLR_BIT11 = 0xFFFFF7FF;
const unsigned int CLR_BIT12 = 0xFFFFEFFF;
const unsigned int CLR_BIT13 = 0xFFFFDFFF;
const unsigned int CLR_BIT14 = 0xFFFFBFFF;
const unsigned int CLR_BIT15 = 0xFFFF7FFF;
const unsigned int CLR_BIT16 = 0xFFFEFFFF;
const unsigned int CLR_BIT17 = 0xFFFDFFFF;
const unsigned int CLR_BIT18 = 0xFFFBFFFF;
const unsigned int CLR_BIT19 = 0xFFF7FFFF;
const unsigned int CLR_BIT20 = 0xFFEFFFFF;
const unsigned int CLR_BIT21 = 0xFFDFFFFF;
const unsigned int CLR_BIT22 = 0xFFBFFFFF;
const unsigned int CLR_BIT23 = 0xFF7FFFFF;
const unsigned int CLR_BIT24 = 0xFEFFFFFF;
const unsigned int CLR_BIT25 = 0xFDFFFFFF;
const unsigned int CLR_BIT26 = 0xFBFFFFFF;
const unsigned int CLR_BIT27 = 0xF7FFFFFF;
const unsigned int CLR_BIT28 = 0xEFFFFFFF;
const unsigned int CLR_BIT29 = 0xDFFFFFFF;
const unsigned int CLR_BIT30 = 0xBFFFFFFF;
const unsigned int CLR_BIT31 = 0x7FFFFFFF;

const unsigned int SET_BIT0 = 0x00000001;
const unsigned int SET_BIT1 = 0x00000002;
const unsigned int SET_BIT2 = 0x00000004;
const unsigned int SET_BIT3 = 0x00000008;
const unsigned int SET_BIT4 = 0x00000010;
const unsigned int SET_BIT5 = 0x00000020;
const unsigned int SET_BIT6 = 0x00000040;
const unsigned int SET_BIT7 = 0x00000080;
const unsigned int SET_BIT8 = 0x00000100;
const unsigned int SET_BIT9 = 0x00000200;
const unsigned int SET_BIT10 = 0x00000400;
const unsigned int SET_BIT11 = 0x00000800;
const unsigned int SET_BIT12 = 0x00001000;
const unsigned int SET_BIT13 = 0x00002000;
const unsigned int SET_BIT14 = 0x00004000;
const unsigned int SET_BIT15 = 0x00008000;
const unsigned int SET_BIT16 = 0x00010000;
const unsigned int SET_BIT17 = 0x00020000;
const unsigned int SET_BIT18 = 0x00040000;
const unsigned int SET_BIT19 = 0x00080000;
const unsigned int SET_BIT20 = 0x00100000;
const unsigned int SET_BIT21 = 0x00200000;
const unsigned int SET_BIT22 = 0x00400000;
const unsigned int SET_BIT23 = 0x00800000;
const unsigned int SET_BIT24 = 0x01000000;
const unsigned int SET_BIT25 = 0x02000000;
const unsigned int SET_BIT26 = 0x04000000;
const unsigned int SET_BIT27 = 0x08000000;
const unsigned int SET_BIT28 = 0x10000000;
const unsigned int SET_BIT29 = 0x20000000;
const unsigned int SET_BIT30 = 0x40000000;
const unsigned int SET_BIT31 = 0x80000000;

const unsigned int BIT0_MASK = 0x00000001;
const unsigned int BIT1_MASK = 0x00000002;
const unsigned int BIT2_MASK = 0x00000004;
const unsigned int BIT3_MASK = 0x00000008;
const unsigned int BIT4_MASK = 0x00000010;
const unsigned int BIT5_MASK = 0x00000020;
const unsigned int BIT6_MASK = 0x00000040;
const unsigned int BIT7_MASK = 0x00000080;
const unsigned int BIT8_MASK = 0x00000100;
const unsigned int BIT9_MASK = 0x00000200;
const unsigned int BIT10_MASK = 0x00000400;
const unsigned int BIT11_MASK = 0x00000800;
const unsigned int BIT12_MASK = 0x00001000;
const unsigned int BIT13_MASK = 0x00002000;
const unsigned int BIT14_MASK = 0x00004000;
const unsigned int BIT15_MASK = 0x00008000;
const unsigned int BIT16_MASK = 0x00010000;
const unsigned int BIT17_MASK = 0x00020000;
const unsigned int BIT18_MASK = 0x00040000;
const unsigned int BIT19_MASK = 0x00080000;
const unsigned int BIT20_MASK = 0x00100000;
const unsigned int BIT21_MASK = 0x00200000;
const unsigned int BIT22_MASK = 0x00400000;
const unsigned int BIT23_MASK = 0x00800000;
const unsigned int BIT24_MASK = 0x01000000;
const unsigned int BIT25_MASK = 0x02000000;
const unsigned int BIT26_MASK = 0x04000000;
const unsigned int BIT27_MASK = 0x08000000;
const unsigned int BIT28_MASK = 0x10000000;
const unsigned int BIT29_MASK = 0x20000000;
const unsigned int BIT30_MASK = 0x40000000;
const unsigned int BIT31_MASK = 0x80000000;

//===========================================================
//  pseudo-functions (JTA)
//===========================================================

#define EXTRACT_BIT0(x) (x & BIT0_MASK)
#define EXTRACT_BIT1(x) ((x & BIT1_MASK) >> 1)
#define EXTRACT_BIT2(x) ((x & BIT2_MASK) >> 2)
#define EXTRACT_BIT3(x) ((x & BIT3_MASK) >> 3)
#define EXTRACT_BIT4(x) ((x & BIT4_MASK) >> 4)
#define EXTRACT_BIT5(x) ((x & BIT5_MASK) >> 5)
#define EXTRACT_BIT6(x) ((x & BIT6_MASK) >> 6)
#define EXTRACT_BIT7(x) ((x & BIT7_MASK) >> 7)
#define EXTRACT_BIT8(x) ((x & BIT8_MASK) >> 8)
#define EXTRACT_BIT9(x) ((x & BIT9_MASK) >> 9)
#define EXTRACT_BIT10(x) ((x & BIT10_MASK) >> 10)
#define EXTRACT_BIT11(x) ((x & BIT11_MASK) >> 11)
#define EXTRACT_BIT12(x) ((x & BIT12_MASK) >> 12)
#define EXTRACT_BIT13(x) ((x & BIT13_MASK) >> 13)
#define EXTRACT_BIT14(x) ((x & BIT14_MASK) >> 14)
#define EXTRACT_BIT15(x) ((x & BIT15_MASK) >> 15)
#define EXTRACT_BIT16(x) ((x & BIT16_MASK) >> 16)
#define EXTRACT_BIT17(x) ((x & BIT17_MASK) >> 17)
#define EXTRACT_BIT18(x) ((x & BIT18_MASK) >> 18)
#define EXTRACT_BIT19(x) ((x & BIT19_MASK) >> 19)
#define EXTRACT_BIT20(x) ((x & BIT20_MASK) >> 20)
#define EXTRACT_BIT21(x) ((x & BIT21_MASK) >> 21)
#define EXTRACT_BIT22(x) ((x & BIT22_MASK) >> 22)
#define EXTRACT_BIT23(x) ((x & BIT23_MASK) >> 23)
#define EXTRACT_BIT24(x) ((x & BIT24_MASK) >> 24)
#define EXTRACT_BIT25(x) ((x & BIT25_MASK) >> 25)
#define EXTRACT_BIT26(x) ((x & BIT26_MASK) >> 26)
#define EXTRACT_BIT27(x) ((x & BIT27_MASK) >> 27)
#define EXTRACT_BIT28(x) ((x & BIT28_MASK) >> 28)
#define EXTRACT_BIT29(x) ((x & BIT29_MASK) >> 22)
#define EXTRACT_BIT30(x) ((x & BIT30_MASK) >> 30)
#define EXTRACT_BIT31(x) ((x & BIT31_MASK) >> 31)

#define READ_MOD_WRITE_BIT0(x)  CLR_BIT0, ((x & 0x1) << 0)
#define READ_MOD_WRITE_BIT1(x)  CLR_BIT1, ((x & 0x1) << 1)
#define READ_MOD_WRITE_BIT2(x)  CLR_BIT2, ((x & 0x1) << 2)
#define READ_MOD_WRITE_BIT3(x)  CLR_BIT3, ((x & 0x1) << 3)
#define READ_MOD_WRITE_BIT4(x)  CLR_BIT4, ((x & 0x1) << 4)
#define READ_MOD_WRITE_BIT5(x)  CLR_BIT5, ((x & 0x1) << 5)
#define READ_MOD_WRITE_BIT6(x)  CLR_BIT6, ((x & 0x1) << 6)
#define READ_MOD_WRITE_BIT7(x)  CLR_BIT7, ((x & 0x1) << 7)
#define READ_MOD_WRITE_BIT8(x)  CLR_BIT8, ((x & 0x1) << 8)
#define READ_MOD_WRITE_BIT9(x)  CLR_BIT9, ((x & 0x1) << 9)
#define READ_MOD_WRITE_BIT10(x) CLR_BIT10,((x & 0x1) << 10)
#define READ_MOD_WRITE_BIT11(x) CLR_BIT11,((x & 0x1) << 11)
#define READ_MOD_WRITE_BIT12(x) CLR_BIT12,((x & 0x1) << 12)
#define READ_MOD_WRITE_BIT13(x) CLR_BIT13,((x & 0x1) << 13)
#define READ_MOD_WRITE_BIT14(x) CLR_BIT14,((x & 0x1) << 14)
#define READ_MOD_WRITE_BIT15(x) CLR_BIT15,((x & 0x1) << 15)
#define READ_MOD_WRITE_BIT16(x) CLR_BIT16,((x & 0x1) << 16)
#define READ_MOD_WRITE_BIT17(x) CLR_BIT17,((x & 0x1) << 17)
#define READ_MOD_WRITE_BIT18(x) CLR_BIT18,((x & 0x1) << 18)
#define READ_MOD_WRITE_BIT19(x) CLR_BIT19,((x & 0x1) << 19)
#define READ_MOD_WRITE_BIT20(x) CLR_BIT20,((x & 0x1) << 20)
#define READ_MOD_WRITE_BIT21(x) CLR_BIT21,((x & 0x1) << 21)
#define READ_MOD_WRITE_BIT22(x) CLR_BIT22,((x & 0x1) << 22)
#define READ_MOD_WRITE_BIT23(x) CLR_BIT23,((x & 0x1) << 23)
#define READ_MOD_WRITE_BIT24(x) CLR_BIT24,((x & 0x1) << 24)
#define READ_MOD_WRITE_BIT25(x) CLR_BIT25,((x & 0x1) << 25)
#define READ_MOD_WRITE_BIT26(x) CLR_BIT26,((x & 0x1) << 26)
#define READ_MOD_WRITE_BIT27(x) CLR_BIT27,((x & 0x1) << 27)
#define READ_MOD_WRITE_BIT28(x) CLR_BIT28,((x & 0x1) << 28)
#define READ_MOD_WRITE_BIT29(x) CLR_BIT29,((x & 0x1) << 29)
#define READ_MOD_WRITE_BIT30(x) CLR_BIT30,((x & 0x1) << 30)
#define READ_MOD_WRITE_BIT31(x) CLR_BIT31,((x & 0x1) << 31)

/* Gain Calibration */

float M = 200.0;                // changed from 350.0
//float M = 350.0;
//float M = 344.0;            //changed by Tianheng for 3.5us
float ehigain[NGE + 1];
float ehioffset[NGE + 1];
float ehibase[NGE + 1];
float PZ[NGE + 1];
void SetBeta ();

float ath[NGE + 1], bth[NGE + 1], slth[NGE + 1],asth[NGE + 1], bsth[NGE + 1],xth[NGE + 1], wth[NGE + 1];
float Pth[NGE + 1], BLPth[NGE + 1], Pdtth[NGE + 1][9]; 
float Kth=146;
float dBLth=500;  //70;
int Navth=3000;
int ith[NGE + 1];
int isth[NGE + 1];
float BLTth[NGE + 1];
float BLCth[NGE + 1][3000];   //3000 = Navth;
int BLMth[NGE+1][1201];       //1201=dBLth*2+1
int BLMaxth[NGE+1];
int BLMCth[NGE+1];
float BLSth[NGE + 1];
float BLth,BLavfth,BLavth[NGE + 1]; 
float BLth1,BLth2;

float Lth,Lbth,Leth;
float Eth,S1th,S2th,PLth,P1th,P2th;
unsigned short int Idth,Bsth,S1bth,S1eth,S2bth,S2eth;
unsigned long long int Tth,TPth,dTth;
int S0th,Tsth[1025];
int Nth=0;

unsigned short int Trath[1025];

#define th 0
#if(th)
/*Root File for the calibration_Tianheng*/
//  TFile f1("/home/dgs/data2/user/254No/ROOT_FILES/Bi503_th.root","RECREATE");
//  TFile f1("/home/dgs/data2/user/data/Eu0_th.root","RECREATE");
 TFile f1("./Calibration/Ge_thbiAG2.root","RECREATE");
  TTree tr("tr","gamsET_cal");
#endif

/* Other variables */

unsigned long long int EvTimeStam0 = 0;

/* pointers to ROOT spectra */

TH1D *hEventCounter;
TH2F *hGeCounter, *hBGOCounter;
TH2F *hEhiRaw, *hEhiCln, *hEhiDrty, *hEhiCln_nodop;
TH2F *hBLRaw;
TH2F *hGeBGO_DT;
TH2F *hTrB, *hFwB, *hE_TrB[NGE + 1], *hE_TrBcorr40, *hEhi_TrBcorr40;
TH2F *hCorrB1, *hCorrB2, *hCorrB3, *hts_bl, *hts_blc, *hts_bl2, *hts_blc2, *hts_sbl, *hts_sblc, *hts_E, *hts_Ecorr;
TH2F *hm1, *hm2, *hm2last, *hm12;
TH2F *pzraw;
TH2F *base1_diff;
TH2F *base2_diff;
TH2F *h2_ebl, *h2_blts, *h2_sblts;
TH2F *h2_ebl1, *h2_ebl2, *h2_ebl3, *h2_ebl4, *h2_ebl5, *h2_ebl6;
TH1D *hrr;
TH2F *e2_e1vse1[NGE + 1];
TH2F *SZe2_e1vse1[NGE + 1];
TH2F *hbase1,*hbase2;


TH1D *hdT;
TH2F *hdT_E;

extern TH1D *ehi[MAXDETPOS + 1];

/* parameters */

extern DGSEVENT DGSEvent[MAXCOINEV];
extern int ng;
extern PARS Pars;
int tlkup[NCHANNELS];
int tid[NCHANNELS];

extern double angle[NGSGE], anglephi[NGSGE];

extern double DopCorFac[NGSGE], ACFac[NGSGE];

int en = 0;
float bl_array[5000];

float base1, base2, base1_av = 0;
float ave_base[NGE + 1];

long long int dgsHeaderID[20];


int firstbl = 1;
int firstbl2[NGE + 1];

/*-----------------------------------------------------*/

int
exit_dgs ()
{
  /* declarations */

  int i;

  /* done */

  printf ("\n");
  printf ("last running baseline values:\n");
  for (i = 0; i < NGE; i++)
    if (ave_base[i] > 0)
      printf ("ave_base[%3i]=%8.1f\n", i, ave_base[i]);
  printf ("\n");

#if(th)
  f1.cd();
  tr.Write();
  f1.Close();
#endif

  return (0);
};

/*-----------------------------------------------------*/

int
sup_dgs ()
{
  /* declarations */

  char str[STRLEN];
  int i, i1, i2, i7, i8;
  FILE *fp;

  void getcal (char *);

  char file_name[] = "./dgscal.dat";    // place this is sort directory


// functions for making root histograms 

  TH2F *make2D (const char *, int, int, int, int, int, int);
  TH1D *make1D (const char *, int, int, int);

/* base diff spectra */

//  base1_diff = make2D ("base1_diff", NGE, 1, NGE, 2048, -1024, 1024);
//  base2_diff = make2D ("base2_diff", NGE, 1, NGE, 2048, -1024, 1024);
//  hbase1      = make2D ("hbase1", NGE, 1, NGE, 4096, 0, 16384);
//  hbase2      = make2D ("hbase2", NGE, 1, NGE, 4096, 0, 16384);
#if(th)
  tr.Branch("Geid",  &Idth,   "Idth/s");
  tr.Branch("N",  &Nth,   "Nth/I");
  //  tr.Branch("Bs",  &Bsth,   "Bsth/s");
   tr.Branch("S1b",  &S1bth,   "S1bth/s");
   tr.Branch("S1e",  &S1eth,   "S1eth/s");
   tr.Branch("S2b",  &S2bth,   "S2bth/s");
   tr.Branch("S2e",  &S2eth,   "S2eth/s");
   tr.Branch("T",  &Tth,   "Tth/l");
//  tr.Branch("TP",  &TPth,   "TPth/l");
  tr.Branch("dT",  &dTth,   "dTth/l");
  tr.Branch("E",  &Eth,  "Eth/F" );
  tr.Branch("S1", &S1th,  "S1th/F");
  tr.Branch("S2", &S2th,  "S2th/F");
  tr.Branch("S0", &S0th,  "S0th/I");
  tr.Branch("BLav", &BLavfth,  "BLavfth/F");
  tr.Branch("P", &PLth,  "PLth/F");
  //  tr.Branch("L", &Lth,  "Lth/F");
  //   tr.Branch("Trace", Trath,  "Trath[1025]/s");
  //   tr.Branch("T", Tsth,  "Tsth[1025]/I");
  //  tr.Branch("Lb", &Lbth,  "Lbth/F");
  //  tr.Branch("Le", &Leth,  "Leth/F");
#endif

  h2_ebl = make2D ("ebl", 2000, 0, 2000, 4000, 0, 8000);
  h2_sblts = make2D ("sblts", 2000, 0, 2000, 2000, 0, 2000);
  h2_blts = make2D ("blts", 2000, 0, 2000, 2000, 0, 2000);

  hdT=make1D("dT",1000,0,20000);
  hdT_E=make2D("dT_E",1000,0,20000,2000,0,2000);


// 2-D's for Rate

  hEventCounter = make1D ("EvntCounter", 14400, 0, 14400);      // Good for 4 hours if Counts/sec
  hrr = make1D ("rr", 1012, 0, 2);
  hGeCounter = make2D ("GeCounter", 14400, 0, 14400, NGE, 1, NGE);
  hBGOCounter = make2D ("BGOCounter", 14400, 0, 14400, NGE, 1, NGE);
  for (i = 1; i <= NGE; i++)
    {
      sprintf (str, "e2_e1vse1_%3.3i", i);
//      e2_e1vse1[i] = make2D (str, 1024, 1, 16384, 1024, 1, 16384);
      sprintf (str, "SZe2_e1vse1_%3.3i", i);
//      SZe2_e1vse1[i] = make2D (str, 2048, 1, 2048, 1024, 1, 12000);
    };

// 2-D's for Energy

  hEhiRaw = make2D ("EhiRaw", LENSP, 0, LENSP + 1, NGE+1, 0, NGE+1);
  hBLRaw = make2D ("BLRaw", 20000, -10000, 10000, NGE+1, 0, NGE+1);
  hEhiCln = make2D ("EhiCln", LENSP, 0, LENSP + 1, NGE+1, 0, NGE+1);
  hEhiCln_nodop = make2D ("EhiCln_nodop", LENSP, 0, LENSP + 1, NGE+1, 0, NGE+1);
  hEhiDrty = make2D ("EhiDrty", LENSP, 0, LENSP + 1, NGE+1, 0, NGE+1);

// 2-D's for Tacs

  hGeBGO_DT = make2D ("GeBGO_DT", 400, -200, 200, NGE, 1, NGE);

// 2-D's for PZ and Baseline

  pzraw = make2D ("pzraw", NGE, 1, NGE, 2000, 0, 2.0);

/* list what we have */

  //printf (" we have define the following spectra:\n");

  Pars.wlist = gDirectory->GetList ();
//  Pars.wlist->Print ();

  /* -------------------- */
  /* read in the map file */
  /* -------------------- */

  for (i = 0; i < NCHANNELS; i++)
    {
      tlkup[i] = NOTHING;
      tid[i] = NOTHING;
    };

  fp = fopen ("map.dat", "r");
  if (fp == NULL)
    {
      printf ("need a \"map.dat\" file to run\n");
      exit (1);
      //system("./mkMap > map.dat");
      //printf("just made you one...\n");
      //fp = fopen ("map.dat", "r");
      //assert(fp != NULL);
    };

  printf ("\nmapping - started\n");

  i2 = fscanf (fp, "\n%i %i %i %s", &i1, &i7, &i8, str);
  printf ("%i %i %i %s\n", i1, i7, i8, str);
  while (i2 == 4)
    {
      tlkup[i1] = i7;
      tid[i1] = i8;
      i2 = fscanf (fp, "\n%i %i %i %s", &i1, &i7, &i8, str);
      //printf ("%i %i %i %s\n", i1, i7, i8, str);
    };
  fclose (fp);

  printf ("\nmapping - complete!!\n");

  /* Set Default Calibration and PZ parameters */

  for (i = 0; i <= NGE + 1; i++)
    {
//      printf ("\nsup_dgs %i \n", i);
      ehigain[i] = 1.0;
      ehioffset[i] = 0.0;
      PZ[i] = 1.0;
      ehibase[i] = 0.0;
    };

  for (i = 0; i < 5000; i++)
    {
      bl_array[i] = 0;
    }

  for (i = 0; i < NGE + 1; i++)
    {
      ave_base[i] = 0;
      firstbl2[i] = 1;
    }


  /* get the DGS calibration file */

  getcal (file_name);

  /* done */

  printf ("\nsup_dgs done!!\n");

  return (0);

};

/*------------------------------------------------------------------*/

int
DGSEvDecompose_v3_ds (unsigned int *ev, int len, DGSEVENT * DGSEvent)
{

  /* declarations */

  int i, k, i1;
  unsigned int ui0 = 0, ui1 = 0, ui2 = 0;
  int rawE;
  unsigned int t1 = 0, t2 = 0, t3 = 0, t4 = 0;
  unsigned long long int ulli1;


  if (Pars.CurEvNo <= Pars.NumToPrint)
    printf ("entered DGSEvDecompose_v3:\n");

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

  // Decode the generic part of the header.

  DGSEvent->chan_id = (*(ev + 0) & 0x0000000f);
  DGSEvent->board_id = ((*(ev + 0) & 0x0000fff0) >> 4); // USER_DEF
  DGSEvent->id = DGSEvent->board_id * 10 + DGSEvent->chan_id;
  DGSEvent->packet_length = ((*(ev + 0) & 0x07ff0000) >> 16);
  DGSEvent->geo_addr = ((*(ev + 0) & 0xf8000000) >> 27);

  DGSEvent->header_type = ((*(ev + 2) & 0x000f0000) >> 16);
  DGSEvent->event_type = ((*(ev + 2) & 0x03800000) >> 23);      // hope this is right.
  DGSEvent->header_length = ((*(ev + 2) & 0xFC000000) >> 26);   // hope this is right.

  /* extract the LED time stamp */

#if(0)
  DGSEvent->event_timestamp = 0;
  DGSEvent->event_timestamp = (unsigned long long int)(*(ev + 1) & 0x3FFFFFFF);
#endif
#if(1)
  DGSEvent->event_timestamp = (unsigned long long int) *(ev + 1);
  ulli1 = (unsigned long long int) (*(ev + 2) & 0x0000ffff);
  ulli1 = (ulli1 << 32);
  DGSEvent->event_timestamp += ulli1;
 #endif

  /* store the type and type ID */

  DGSEvent->tpe = tlkup[DGSEvent->id];
  DGSEvent->tid = tid[DGSEvent->id];
  DGSEvent->flag = 0;

  /* extract the header values based on header type */

  if (Pars.CurEvNo <= Pars.NumToPrint)
    printf ("DGSEvent->header_type=%i\n", DGSEvent->header_type);
  switch (DGSEvent->header_type)
    {

    case 0:                    // Original LED header (Compatibility mode)
      DGSEvent->timestamp_match_flag = 0;
      DGSEvent->external_disc_flag = 0;
      DGSEvent->cfd_valid_flag = 0;
      DGSEvent->pileup_only_flag = 0;
      DGSEvent->cfd_sample_0 = 0;
      DGSEvent->cfd_sample_1 = 0;
      DGSEvent->cfd_sample_2 = 0;

      DGSEvent->peak_valid_flag = ((*(ev + 3) & 0x00000200) >> 9);      // Word 3: 9
      DGSEvent->offset_flag = ((*(ev + 3) & 0x00000400) >> 10); // Word 3: 10
      DGSEvent->sync_error_flag = ((*(ev + 3) & 0x00001000) >> 12);     // Word 3: 12
      DGSEvent->general_error_flag = ((*(ev + 3) & 0x00002000) >> 13);  // Word 3: 13
      DGSEvent->pileup_flag = ((*(ev + 3) & 0x00008000) >> 15); // Word 3: 15
      DGSEvent->last_disc_timestamp = (((unsigned long long int) (*(ev + 3) & 0xFFFF0000)) >> 16) |     // Word 3: 31..16 & 
        (((unsigned long long int) (*(ev + 4) & 0xFFFFFFFF)) << 16);    // Word 4 :31..0  
      DGSEvent->sampled_baseline = ((*(ev + 5) & 0x00FFFFFF) >> 0);     // Word 5: 23..0
      DGSEvent->sum1 = ((*(ev + 7) & 0x00FFFFFF) >> 0); // Word 7: 23..0
      DGSEvent->sum2 = ((*(ev + 7) & 0xFF000000) >> 28) |       // Word 7: 31..24 & 
        ((*(ev + 8) & 0x0000FFFF) << 8);        // Word 8: 15..0 
      DGSEvent->peak_timestamp = (((unsigned long long int) (*(ev + 8) & 0xFFFF0000)) >> 16) |  // Word 8: 31..16 & 
        (((unsigned long long int) (*(ev + 9) & 0xFFFFFFFF)) << 16);    // Word 9 :31..0  
      DGSEvent->m2begin = ((*(ev + 10) & 0x00003FFF) >> 0);     // Word 10:13..0
      DGSEvent->m2end = ((*(ev + 10) & 0x3FFF0000) >> 16);      // Word 10:29..16
      DGSEvent->m1begin = ((*(ev + 11) & 0x00003FFF) >> 0);     // Word 11:13..0
      DGSEvent->m1end = ((*(ev + 11) & 0x3FFF0000) >> 16);      // Word 11:29..16
      DGSEvent->peak_sample = ((*(ev + 12) & 0x00003FFF) >> 0); // Word 12:13..0
      DGSEvent->base_sample = ((*(ev + 12) & 0x3FFF0000) >> 16);        // Word 12:29..16
      break;

    case 1:                    // New LED Header
      DGSEvent->timestamp_match_flag = 0;
      DGSEvent->cfd_valid_flag = 0;
      DGSEvent->cfd_sample_0 = 0;
      DGSEvent->cfd_sample_1 = 0;
      DGSEvent->cfd_sample_2 = 0;

      DGSEvent->external_disc_flag = ((*(ev + 3) & 0x00000100) >> 8);   // Word 3: 8
      DGSEvent->peak_valid_flag = ((*(ev + 3) & 0x00000200) >> 9);      // Word 3: 9
      DGSEvent->offset_flag = ((*(ev + 3) & 0x00000400) >> 10); // Word 3: 10
      DGSEvent->sync_error_flag = ((*(ev + 3) & 0x00001000) >> 12);     // Word 3: 12
      DGSEvent->general_error_flag = ((*(ev + 3) & 0x00002000) >> 13);  // Word 3: 13
      DGSEvent->pileup_only_flag = ((*(ev + 3) & 0x00004000) >> 14);    // Word 3: 14
      DGSEvent->pileup_flag = ((*(ev + 3) & 0x00008000) >> 15); // Word 3: 15  
      DGSEvent->last_disc_timestamp = (((unsigned long long int) (*(ev + 3) & 0xFFFF0000)) >> 16) |     // Word 3: 31..16 & 
        (((unsigned long long int) (*(ev + 4) & 0xFFFFFFFF)) << 16);    // Word 4 :31..0  
      DGSEvent->sampled_baseline = ((*(ev + 5) & 0x00FFFFFF) >> 0);     // Word 5: 23..0
      DGSEvent->sum1 = ((*(ev + 7) & 0x00FFFFFF) >> 0); // Word 7: 23..0
      DGSEvent->sum2 = ((*(ev + 7) & 0xFF000000) >> 28) |       // Word 7: 31..24 & 
        ((*(ev + 8) & 0x0000FFFF) << 8);        // Word 8: 15..0 
      DGSEvent->peak_timestamp = (((unsigned long long int) (*(ev + 8) & 0xFFFF0000)) >> 16) |  // Word 8: 31..16 & 
        (((unsigned long long int) (*(ev + 9) & 0xFFFFFFFF)) << 16);    // Word 9: 31..0  
      DGSEvent->m1begin = ((*(ev + 10) & 0x00003FFF) >> 0);     // Word 10:13..0
      //DGSEvent->m1end          = ((*(ev + 10) & 0x3FFF0000) >> 16);                             // Word 10:29..16

      DGSEvent->m2last_end_sample = ((*(ev + 9) & 0x00003FFF) >> 0);    // Word 9:13..0  
      DGSEvent->m2last_begin_sample = ((*(ev + 10) & 0x3FFF0000) >> 16);        // Word 10:13..0

      DGSEvent->m2end = ((*(ev + 11) & 0x00003FFF) >> 0);       // Word 11:13..0
      DGSEvent->m2begin = ((*(ev + 11) & 0x3FFF0000) >> 16);    // Word 11:29..16
      DGSEvent->peak_sample = ((*(ev + 12) & 0x00003FFF) >> 0); // Word 12:13..0
      DGSEvent->base_sample = ((*(ev + 12) & 0x3FFF0000) >> 16);        // Word 12:29..16
      break;

    case 4:                    // was 2              // CFD Header
      DGSEvent->timestamp_match_flag = ((*(ev + 3) & 0x00000080) >> 7); // Word 3: 7
      DGSEvent->external_disc_flag = ((*(ev + 3) & 0x00000100) >> 8);   // Word 3: 8
      DGSEvent->peak_valid_flag = ((*(ev + 3) & 0x00000200) >> 9);      // Word 3: 9
      DGSEvent->offset_flag = ((*(ev + 3) & 0x00000400) >> 10); // Word 3: 10
      DGSEvent->cfd_valid_flag = ((*(ev + 3) & 0x00000800) >> 11);      // Word 3: 11
      DGSEvent->sync_error_flag = ((*(ev + 3) & 0x00001000) >> 12);     // Word 3: 12
      DGSEvent->general_error_flag = ((*(ev + 3) & 0x00002000) >> 13);  // Word 3: 13
      DGSEvent->pileup_only_flag = ((*(ev + 3) & 0x00004000) >> 14);    // Word 3: 14
      DGSEvent->pileup_flag = ((*(ev + 3) & 0x00008000) >> 15); // Word 3: 15
      //      DGSEvent->last_disc_timestamp = (((unsigned long long int) (*(ev + 3) & 0xFFFF0000)) >> 16) |     // Word 3: 31..16 &
      //        (((unsigned long long int) (*(ev + 4) & 0xFFFFFFFF)) << 16);    // Word 4 :31..0
      DGSEvent->last_disc_timestamp = (((unsigned long long int) (*(ev + 3) & 0xFFFF0000)) >> 16) |     // Word 3: 31..16 &
        (((unsigned long long int) (*(ev + 4) & 0x00003FFF)) << 16);    // Word 4 :31..0
      DGSEvent->cfd_sample_0 = ((*(ev + 4) & 0x3FFF0000) >> 16);        // Word 4: 29..16
      //     DGSEvent->sampled_baseline = ((*(ev + 5) & 0x00FFFFFF) >> 0);     // Word 5: 23..0
     DGSEvent->sampled_baseline = ((*(ev + 5) & 0x00FFFFFF) >> 6);     // Word 5: 23..0
     DGSEvent->cfd_sample_1 = ((*(ev + 6) & 0x00003FFF) >> 0); // Word 6: 13..0
      DGSEvent->cfd_sample_2 = ((*(ev + 6) & 0x3FFF0000) >> 16);        // Word 6: 29..16
      DGSEvent->sum1 = ((*(ev + 7) & 0x00FFFFFF) >> 0); // Word 7: 23..0
      DGSEvent->sum2 = ((*(ev + 7) & 0xFF000000) >> 24) |       //was 28                      // Word 7: 31..24 & 
        ((*(ev + 8) & 0x0000FFFF) << 8);        // Word 8: 15..0 
      DGSEvent->peak_timestamp = (((unsigned long long int) (*(ev + 8) & 0xFFFF0000)) >> 16) |  // Word 8: 31..16 & 
        (((unsigned long long int) (*(ev + 9) & 0xFFFFFFFF)) << 16);    // Word 9: 31..0  

      DGSEvent->m2end = ((*(ev + 10) & 0x00003FFF) >> 0);       // Word 10:13..0   POST_RISE_SAMPLE
      DGSEvent->m2begin = ((*(ev + 10) & 0x3FFF0000) >> 16);    // Word 10:29..16  LAST_POST_RISE_SAMPLE
      DGSEvent->m1begin = ((*(ev + 11) & 0x00003FFF) >> 0);     // Word 11:13..0   PRE_RISE_ENTER_SAMPLE
      DGSEvent->m1end = ((*(ev + 11) & 0x3FFF0000) >> 16);      // Word 11:29..16  PRE_RISE_LEAVE_SAMPLE

      DGSEvent->peak_sample = ((*(ev + 12) & 0x00003FFF) >> 0); // Word 12:13..0
      DGSEvent->base_sample = ((*(ev + 12) & 0x3FFF0000) >> 16);        // Word 12:29..16
      break;

    case 5:                    // New LED Header //TL: copy of 1 for now
      DGSEvent->timestamp_match_flag = 0;
      DGSEvent->cfd_valid_flag = 0;
      DGSEvent->cfd_sample_0 = 0;
      DGSEvent->cfd_sample_1 = 0;
      DGSEvent->cfd_sample_2 = 0;

      DGSEvent->external_disc_flag = ((*(ev + 3) & 0x00000100) >> 8);   // Word 3: 8
      DGSEvent->peak_valid_flag = ((*(ev + 3) & 0x00000200) >> 9);      // Word 3: 9
      DGSEvent->offset_flag = ((*(ev + 3) & 0x00000400) >> 10); // Word 3: 10
      DGSEvent->sync_error_flag = ((*(ev + 3) & 0x00001000) >> 12);     // Word 3: 12
      DGSEvent->general_error_flag = ((*(ev + 3) & 0x00002000) >> 13);  // Word 3: 13
      DGSEvent->pileup_only_flag = ((*(ev + 3) & 0x00004000) >> 14);    // Word 3: 14
      DGSEvent->pileup_flag = ((*(ev + 3) & 0x00008000) >> 15); // Word 3: 15  
      DGSEvent->last_disc_timestamp = (((unsigned long long int) (*(ev + 3) & 0xFFFF0000)) >> 16) |     // Word 3: 31..16 & 
        (((unsigned long long int) (*(ev + 4) & 0xFFFFFFFF)) << 16);    // Word 4 :31..0  
      DGSEvent->sampled_baseline = ((*(ev + 5) & 0x00FFFFFF) >> 0);     // Word 5: 23..0
      DGSEvent->sum1 = ((*(ev + 7) & 0x00FFFFFF) >> 0); // Word 7: 23..0
      DGSEvent->sum2 = ((*(ev + 7) & 0xFF000000) >> 28) |       // Word 7: 31..24 & 
        ((*(ev + 8) & 0x0000FFFF) << 8);        // Word 8: 15..0 
      DGSEvent->peak_timestamp = (((unsigned long long int) (*(ev + 8) & 0xFFFF0000)) >> 16) |  // Word 8: 31..16 & 
        (((unsigned long long int) (*(ev + 9) & 0xFFFFFFFF)) << 16);    // Word 9: 31..0  
      DGSEvent->m1begin = ((*(ev + 10) & 0x00003FFF) >> 0);     // Word 10:13..0
      //DGSEvent->m1end          = ((*(ev + 10) & 0x3FFF0000) >> 16);                             // Word 10:29..16

      DGSEvent->m2last_end_sample = ((*(ev + 9) & 0x00003FFF) >> 0);    // Word 9:13..0  
      DGSEvent->m2last_begin_sample = ((*(ev + 10) & 0x3FFF0000) >> 16);        // Word 10:13..0

      DGSEvent->m2end = ((*(ev + 11) & 0x00003FFF) >> 0);       // Word 11:13..0
      DGSEvent->m2begin = ((*(ev + 11) & 0x3FFF0000) >> 16);    // Word 11:29..16
      DGSEvent->peak_sample = ((*(ev + 12) & 0x00003FFF) >> 0); // Word 12:13..0
      DGSEvent->base_sample = ((*(ev + 12) & 0x3FFF0000) >> 16);        // Word 12:29..16
      break;

    case 6:                    // was 2              // CFD Header   //TL: copy of 2 for now
      DGSEvent->timestamp_match_flag = ((*(ev + 3) & 0x00000080) >> 7); // Word 3: 7
      DGSEvent->external_disc_flag = ((*(ev + 3) & 0x00000100) >> 8);   // Word 3: 8
      DGSEvent->peak_valid_flag = ((*(ev + 3) & 0x00000200) >> 9);      // Word 3: 9
      DGSEvent->offset_flag = ((*(ev + 3) & 0x00000400) >> 10); // Word 3: 10
      DGSEvent->cfd_valid_flag = ((*(ev + 3) & 0x00000800) >> 11);      // Word 3: 11
      DGSEvent->sync_error_flag = ((*(ev + 3) & 0x00001000) >> 12);     // Word 3: 12
      DGSEvent->general_error_flag = ((*(ev + 3) & 0x00002000) >> 13);  // Word 3: 13
      DGSEvent->pileup_only_flag = ((*(ev + 3) & 0x00004000) >> 14);    // Word 3: 14
      DGSEvent->pileup_flag = ((*(ev + 3) & 0x00008000) >> 15); // Word 3: 15
      DGSEvent->last_disc_timestamp = (((unsigned long long int) (*(ev + 3) & 0xFFFF0000)) >> 16) |     // Word 3: 31..16 &
        (((unsigned long long int) (*(ev + 4) & 0xFFFFFFFF)) << 16);    // Word 4 :31..0
      DGSEvent->cfd_sample_0 = ((*(ev + 4) & 0x3FFF0000) >> 16);        // Word 4: 29..16
      DGSEvent->sampled_baseline = ((*(ev + 5) & 0x00FFFFFF) >> 0);     // Word 5: 23..0
      DGSEvent->cfd_sample_1 = ((*(ev + 6) & 0x00003FFF) >> 0); // Word 6: 13..0
      DGSEvent->cfd_sample_2 = ((*(ev + 6) & 0x3FFF0000) >> 16);        // Word 6: 29..16
      DGSEvent->sum1 = ((*(ev + 7) & 0x00FFFFFF) >> 0); // Word 7: 23..0
      DGSEvent->sum2 = ((*(ev + 7) & 0xFF000000) >> 24) |       //was 28                      // Word 7: 31..24 & 
        ((*(ev + 8) & 0x0000FFFF) << 8);        // Word 8: 15..0 
      DGSEvent->peak_timestamp = (((unsigned long long int) (*(ev + 8) & 0xFFFF0000)) >> 16) |  // Word 8: 31..16 & 
        (((unsigned long long int) (*(ev + 9) & 0xFFFFFFFF)) << 16);    // Word 9: 31..0  

      DGSEvent->m2end = ((*(ev + 10) & 0x00003FFF) >> 0);       // Word 10:13..0   POST_RISE_SAMPLE
      DGSEvent->m2begin = ((*(ev + 10) & 0x3FFF0000) >> 16);    // Word 10:29..16  LAST_POST_RISE_SAMPLE
      DGSEvent->m1begin = ((*(ev + 11) & 0x00003FFF) >> 0);     // Word 11:13..0   PRE_RISE_ENTER_SAMPLE
      DGSEvent->m1end = ((*(ev + 11) & 0x3FFF0000) >> 16);      // Word 11:29..16  PRE_RISE_LEAVE_SAMPLE

      DGSEvent->peak_sample = ((*(ev + 12) & 0x00003FFF) >> 0); // Word 12:13..0
      DGSEvent->base_sample = ((*(ev + 12) & 0x3FFF0000) >> 16);        // Word 12:29..16
      break;

    default:
      printf ("DGSEvent->header_type= %i is unknown to this version of bin_dgs.c;\n", DGSEvent->header_type);
      return (1);
      break;

    }

  DGSEvent->baseline = (DGSEvent->m2end + DGSEvent->m2begin) / 2;


  /* Now load Trace into DGSEvent Structure */

  DGSEvent->traceLen = 0;
  for (i = 13; i < len - 1; i++)
    {
      if (i < 1037)
        {
          DGSEvent->trace[2 * (i - 13)] = (unsigned short int) (*(ev + i) & 0x3fff);
          DGSEvent->trace[2 * (i - 13) + 1] = (unsigned short int) ((*(ev + i) >> 16) & 0x3fff);
          DGSEvent->traceLen += 2;
        }
    }

  //rawE = (int) POST_RISE_SUM - (int) PRE_RISE_SUM;
  //DGSEvent->ehi = rawE / 800;

  if (Pars.CurEvNo <= Pars.NumToPrint)
    printf ("rawE = 0x%8.8x %i, DGSEvent->ehi= %i\n", rawE, rawE, DGSEvent->ehi);

  /* done */

  if (Pars.CurEvNo <= Pars.NumToPrint)
    printf ("exit DGSEvDecompose_v3:\n");

  return (0);

}

/* ----------------------------------------------------------------- */

int
bin_dgs (GEB_EVENT * GEB_event)
{

  /* declarations */
  char str[128];
  int i, j, i1, gsid, e, st;
  unsigned int ui1;
  double d1, d2;

  int RelEvT = 0, DEvT = 0, DTrT = 0, tdiff = 0;
  float Energy, Energy_nodop, top, bot, r1, rr;

  /* prototypes */

  int GebTypeStr (int type, char str[]);

  /* Print debug */

  if (Pars.CurEvNo <= Pars.NumToPrint)
    printf ("entereD bin_dgs:\n");

  /* loop through the coincidence event and fish out DGS data */
  /* (gamma rays) count in ng */

  ng = 0;
  for (i = 0; i < GEB_event->mult; i++)
    {
      if (GEB_event->ptgd[i]->type == GEB_TYPE_DGS)
        {
          if (Pars.CurEvNo <= Pars.NumToPrint)
            {
              GebTypeStr (GEB_event->ptgd[i]->type, str);
              printf ("bin_mode1, %2i> %2i, %s, TS=%lli\n", i, GEB_event->ptgd[i]->type, str,
                      GEB_event->ptgd[i]->timestamp);
            }

          st=DGSEvDecompose_v3_ds ((unsigned int *) GEB_event->ptinp[i], GEB_event->ptgd[i]->length / sizeof (unsigned int),
                             &DGSEvent[ng]);
          if(st!=0) return(0);
          ng++;
        }
    }

  // Initialize EvTimeStam0

  //printf("Stat 1 \n");
  if (Pars.CurEvNo <= Pars.NumToPrint)
    {
      printf ("\n\nCurEvNo: %i\n", Pars.CurEvNo);
      for (i = 0; i < ng; i++)
        printf ("\nDGSEvent[%i].event_timestamp: %llu", i, DGSEvent[i].event_timestamp);
    }

  if (EvTimeStam0 == 0)
    EvTimeStam0 = DGSEvent[0].event_timestamp;
  RelEvT = (int) ((DGSEvent[0].event_timestamp - EvTimeStam0) / 100000000);     // overflow?
  hEventCounter->Fill (RelEvT);

  /* Loop */
  //Zero clearing before the Tianheng's calibration!
  
  ith[NGE + 1]={0};
  isth[NGE + 1]={0};
  BLSth[NGE + 1]={0};
  BLTth[NGE + 1]={0};
  
  
  // BLMth[NGE + 1][1201]={0};
  //  BLTth[NGE + 1]={0};
  // BLMCth[NGE + 1]={0};
   
  for(i=1;i<=110;i++)
    {
      BLavth[i]=BLPth[i];
      //     BLSth[i]=BLPth[i]-dBLth;
      // BLMaxth[i]=BLPth[i];
      }
  /*                //for the test
  BLavth[81]=1174;
  BLavth[83]=1014;
  BLavth[85]=1121;
  BLavth[87]=1016;
  BLavth[89]=2965;
  */
  for (i = 0; i < ng; i++)
    {
      if (DGSEvent[i].tpe == GE)
        {
          gsid = DGSEvent[i].tid;
          hGeCounter->Fill ((int) ((DGSEvent[0].event_timestamp - EvTimeStam0) / 100000000), gsid);

          /* current baseline; SZ's Method1 implementation */

          if (Pars.CurEvNo <= Pars.NumToPrint)
            {
              printf ("\n");
              printf ("m1begin= DGSEvent[%i].m1begin  = %8i w11,lobit\n", i, DGSEvent[i].m1begin);
              printf ("m1end  = DGSEvent[%i].m1end    = %8i w11,upbit ", i, DGSEvent[i].m1end);
              printf ("diff %i\n", DGSEvent[i].m1begin - DGSEvent[i].m1end);
              printf ("m2begin= DGSEvent[%i].m2begin  = %8i w10,upbit\n", i, DGSEvent[i].m2begin);
              printf ("m2end  = DGSEvent[%i].m2end    = %8i w10,lobit ", i, DGSEvent[i].m2end);
              printf ("diff %i\n", DGSEvent[i].m2begin - DGSEvent[i].m2end);
            };
//          base1_diff->Fill ((double) gsid, DGSEvent[i].m1begin - DGSEvent[i].m1end);
//          base2_diff->Fill ((double) gsid, DGSEvent[i].m2begin - DGSEvent[i].m2end);


          top = (float) (DGSEvent[i].sum2 * DGSEvent[i].m1begin - DGSEvent[i].sum1 * DGSEvent[i].m2begin);
          bot = DGSEvent[i].sum2 - DGSEvent[i].sum1;
          bot -= M * (DGSEvent[i].m2begin - DGSEvent[i].m1begin);
          base1 = top / bot;
//          hbase1->Fill((double) gsid, (double)base1);

          top = (float) (DGSEvent[i].m1end * DGSEvent[i].m2begin - DGSEvent[i].m1begin * DGSEvent[i].m2end);
          bot = (DGSEvent[i].m2begin - DGSEvent[i].m2end) - (DGSEvent[i].m1begin - DGSEvent[i].m1end);
          base2 = top / bot;
//          hbase2->Fill((double) gsid, (double)base2);

          /* keep a running average of the baseline */
          /* avoid extreme values (per Darek) */
          /* add constraints from SZ as well */

          rr = base1 / base2;
          if (rr > 0 && rr < 2.0)
            hrr->Fill (double (rr));
          if ((base1 < 4000) && (base1 > 0))
            if ((DGSEvent[i].sum2 - DGSEvent[i].sum1) > 200)
              if (rr > 0.95 && rr < 1.05)
                {
                  if (firstbl2[gsid] == 1)
                    {
                      ave_base[gsid] = base2;
                      firstbl2[gsid] = 0;
                    }
                  else
                    {
                      ave_base[gsid] = ave_base[gsid] * 0.99 + base1 * 0.01;
                    }
                }

          /* this is SZ's baseline restore correction, see */
          /* soffice /home/tl/d6/keep/2018_zhu_baseline.ppt */
          /* ave_base: the running base. note: M is a float */

#if(1)
          Energy = DGSEvent[i].sum2 / M - DGSEvent[i].sum1 * PZ[gsid] / M;
          Energy -= (1. - PZ[gsid]) * ave_base[gsid];
#endif
#if(0)
          Energy = DGSEvent[i].sum2 / M - DGSEvent[i].sum1 * PZ[gsid] / M;
          Energy -= (1. - PZ[gsid]) * base2;
#endif

          if (Pars.CurEvNo <= Pars.NumToPrint)
            {
              printf ("\n");
              printf ("DGSEvent[i].sum2= %i\n", DGSEvent[i].sum2);
              printf ("DGSEvent[i].sum1= %i\n", DGSEvent[i].sum1);
              printf ("M= %f\n", M);
              printf ("PZ[gsid]= %f\n", PZ[gsid]);
              printf ("base1= %8.1f base2= %8.1f, diff= %8.1f rat= %5.2f\n", base1, base2, base1 - base2,
                      base1 / base2);
              printf ("ave_base[gsid]= %f\n", ave_base[gsid]);
              printf ("Energy = %f\n", Energy);
            };


          /* gain match */

          Energy = Energy * ehigain[gsid] + ehioffset[gsid];
          if (Pars.CurEvNo <= Pars.NumToPrint)
            {
              printf ("Energy = %f (calibrated)\n", Energy);
            };

#if(1)
	  /*Overwrite the energy with Tianheng's calibration*/
	  /*	  	  
	  if(slth[gsid]>0)
	    {
	      Energy = (DGSEvent[i].sum2 - DGSEvent[i].sum1)/M;
	      Energy = Energy + DGSEvent[i].sum1/M/slth[gsid];
	      Energy = Energy*ath[gsid] + bth[gsid];
	    }
	  	  
	  */	  
	  /*For dT>1024 steps*/
	  Bsth= DGSEvent[i].base_sample;
//	  Tth= DGSEvent[i].event_timestamp; 
          Tth =  DGSEvent[i].event_timestamp & 0x3FFFFFFF; 
	  TPth=  DGSEvent[i].last_disc_timestamp;
	  if(Tth<TPth){Tth=Tth+1073741823;}  //1073741823=3FFFFFFF
	  //          Eth = (DGSEvent[i].sum2- DGSEvent[i].sum1)/M;
	  //          Blth = DGSEvent[i].sum1/M;
          S1th = DGSEvent[i].sum1;
          S2th = DGSEvent[i].sum2;
	  dTth= Tth-TPth;
	  S0th=  DGSEvent[i].sampled_baseline;
	  S0th=S0th*64;
	  S1bth=DGSEvent[i].m1begin;
          S1eth=DGSEvent[i].m1end;
	  S2bth=DGSEvent[i].m2begin;
          S2eth=DGSEvent[i].m2end;
	  /*
	  if(S1bth+S2eth-S1eth-S2bth!=0)
	    {
	      BLth=float(S1bth*S2eth-S1eth*S2bth)/float(S1bth+S2eth-S1eth-S2bth);
	      if(BLth<BLPth[gsid]-dBLth || BLth>BLPth[gsid]+dBLth)
		{
		  if(BLTth[gsid]!=0){BLth=BLTth[gsid];}
		    else {BLth=BLMaxth[gsid];}
		}
	    }
	  else if(BLTth[gsid]!=0){BLth=BLTth[gsid];}
	  else {BLth=BLMaxth[gsid];}
	  
	   BLTth[gsid]=BLth;

	   BLMth[gsid][int(BLth)-int(BLSth[gsid])]++;

	   if(BLMCth[gsid]==0||BLMth[gsid][int(BLth)-int(BLSth[gsid])]>=BLMCth[gsid])
	     {
	       BLMCth[gsid]=BLMth[gsid][int(BLth)-int(BLSth[gsid])];
	       BLMaxth[gsid]=BLth;
	     }

	   Energy =  S2th/M - (BLMaxth[gsid]+0.5) - (S1th/M - (BLMaxth[gsid]+0.5)) * pow(Pth[gsid], (M+Kth[gsid])/M);
	   Energy = Energy*ath[gsid] + bth[gsid];
*/

	  	  
	  //	  	  if(S1bth+S2eth-S1eth-S2bth!=0)
	  //	             if(2*S1th/M-S1bth-S1eth!=0)

	  
	if(dTth<=500) Pth[gsid] =  Pdtth[gsid][0];
	else if(dTth>500&&dTth<=1000) Pth[gsid] =  Pdtth[gsid][1];
	else if(dTth>1000&&dTth<=1500) Pth[gsid] =  Pdtth[gsid][2];
	else if(dTth>1500&&dTth<=2000) Pth[gsid] =  Pdtth[gsid][3];
	else if(dTth>2000&&dTth<=2500) Pth[gsid] =  Pdtth[gsid][4];
	else if(dTth>2500&&dTth<=3000) Pth[gsid] =  Pdtth[gsid][5];
	else if(dTth>3000&&dTth<=3500) Pth[gsid] =  Pdtth[gsid][6];
	else if(dTth>3500&&dTth<=4000) Pth[gsid] =  Pdtth[gsid][7];
	else  Pth[gsid] =  Pdtth[gsid][8];
	   
      	hBLRaw->Fill ((S1eth-S1bth*Pth[gsid])/(1-Pth[gsid]), gsid);
	  
		  if((S1eth-S1bth*Pth[gsid])/(1-Pth[gsid])>-1000)
	//  	if(2*S1th/M-S1bth-S1eth!=0)
	    {
	      //    	      	      BLth=float(S1bth*S2eth-S1eth*S2bth)/float(S1bth+S2eth-S1eth-S2bth);
	      //        BLth=(S1th/M*S1th/M-S1bth*S1eth)/(2*S1th/M-S1bth-S1eth);
	      	      BLth=(S1eth-S1bth*Pth[gsid])/(1-Pth[gsid]);
		      //	      if(BLth<BLavth[gsid]-dBLth || BLth>BLavth[gsid]+dBLth)
		       if(BLth<BLPth[gsid]-dBLth || BLth>BLPth[gsid]+dBLth)
		{
		  if(BLTth[gsid]!=0){BLth=BLTth[gsid];}
		    else {BLth=BLavth[gsid];}
		}
	    }
	  else if(BLTth[gsid]!=0){BLth=BLTth[gsid];}
	  else {BLth=BLavth[gsid];}
	
	  /*	  
	  if(2*S1th/M-S1bth-S1eth!=0){BLth1=(S1th/M*S1th/M-S1bth*S1eth)/(2*S1th/M-S1bth-S1eth);}
	  if(2*S2th/M-S2bth-S2eth!=0){BLth2=(S2th/M*S2th/M-S2bth*S2eth)/(2*S2th/M-S2bth-S2eth);}

	  if(2*S1th/M-S1bth-S1eth!=0&&2*S2th/M-S2bth-S2eth!=0)
	    {
	      BLth=(BLth1+BLth2)/2;
	      if(BLth<BLavth[gsid]-dBLth || BLth>BLavth[gsid]+dBLth)
		{
		  if(BLTth[gsid]!=0){BLth=BLTth[gsid];}
		    else {BLth=BLavth[gsid];}
		}
	    }
	  else if(2*S1th/M-S1bth-S1eth==0&&2*S2th/M-S2bth-S2eth!=0)
	    {
	      BLth=BLth2;
	     	      if(BLth<BLavth[gsid]-dBLth || BLth>BLavth[gsid]+dBLth)
		{
		  if(BLTth[gsid]!=0){BLth=BLTth[gsid];}
		    else {BLth=BLavth[gsid];}
		} 
	    }
	  else if(2*S2th/M-S2bth-S2eth==0&&2*S1th/M-S1bth-S1eth!=0)
	    {
	      BLth=BLth1;
	      	      if(BLth<BLavth[gsid]-dBLth || BLth>BLavth[gsid]+dBLth)
		{
		  if(BLTth[gsid]!=0){BLth=BLTth[gsid];}
		    else {BLth=BLavth[gsid];}
		}
	    }
	  else if(BLTth[gsid]!=0){BLth=BLTth[gsid];}
	  else {BLth=BLavth[gsid];}
*/

	  BLTth[gsid]=BLth;

	  if(ith[gsid]==Navth)
	    {
	      ith[gsid]=0;
	      if(isth[gsid]==0){isth[gsid]=1;}
	    }

	 if(isth[gsid]==0)
	   {
	     BLCth[gsid][ith[gsid]]=BLth;
	     BLSth[gsid]=BLSth[gsid]+BLCth[gsid][ith[gsid]];
	     BLavth[gsid]=BLSth[gsid]/(ith[gsid]+1);
	   }
	 else
	   {
	     BLSth[gsid]=BLSth[gsid]-BLCth[gsid][ith[gsid]];
	     BLCth[gsid][ith[gsid]]=BLth;
	     BLSth[gsid]=BLSth[gsid]+BLCth[gsid][ith[gsid]];
	     BLavth[gsid]=BLSth[gsid]/Navth;
	   }

	BLavfth=BLavth[gsid];
	ith[gsid]++;

		Energy =  S2th/M - BLavfth - (S1th/M - BLavfth) * pow(Pth[gsid], (M+Kth)/M);
	/*
	if(dTth<=500) Energy =  S2th/M - BLavfth - (S1th/M - BLavfth) * pow(Pdtth[gsid][0], (M+Kth)/M);
	else if(dTth>500&&dTth<=1000) Energy =  S2th/M - BLavfth - (S1th/M - BLavfth) * pow(Pdtth[gsid][1], (M+Kth)/M);
	else if(dTth>1000&&dTth<=1500) Energy =  S2th/M - BLavfth - (S1th/M - BLavfth) * pow(Pdtth[gsid][2], (M+Kth)/M);
	else if(dTth>1500&&dTth<=2000) Energy =  S2th/M - BLavfth - (S1th/M - BLavfth) * pow(Pdtth[gsid][3], (M+Kth)/M);
	else if(dTth>2000&&dTth<=2500) Energy =  S2th/M - BLavfth - (S1th/M - BLavfth) * pow(Pdtth[gsid][4], (M+Kth)/M);
	else if(dTth>2500&&dTth<=3000) Energy =  S2th/M - BLavfth - (S1th/M - BLavfth) * pow(Pdtth[gsid][5], (M+Kth)/M);
	else if(dTth>3000&&dTth<=3500) Energy =  S2th/M - BLavfth - (S1th/M - BLavfth) * pow(Pdtth[gsid][6], (M+Kth)/M);
	else if(dTth>3500&&dTth<=4000) Energy =  S2th/M - BLavfth - (S1th/M - BLavfth) * pow(Pdtth[gsid][7], (M+Kth)/M);
	else  Energy =  S2th/M - BLavfth - (S1th/M - BLavfth) * pow(Pdtth[gsid][8], (M+Kth)/M);
*/
        Energy = Energy*ath[gsid] + bth[gsid];
	 
	  /*	  	      
	  //	  	     if(dTth>=1280&&ath[gsid]<1)  //for 2us
	    	  if(dTth>=1434&&ath[gsid]<1)               //for 3.5us
	    {         
	      Lth = (S0th*M*M+S1th*1024*1024)/(2*(S0th*M-S1th*1024));
	      Energy=(S2th-S1th+S1th*(2*(M+80))/(2*Lth-M))/M;
	      Energy=Energy*ath[gsid]+bth[gsid];
	    }
	  
		  //      else if(dTth<1280&&asth[gsid]<1)    //for 2us
	      else if(dTth<1434&&asth[gsid]<1)             //for 3.5us
	    {
	      Lth=(S0th*M*M+S1th*(dTth-M-60)*(dTth-M-60)-S1th*2*xth[gsid]*wth[gsid]*(1024-dTth+M+60))/(2*(S0th*M-S1th*(dTth-M-60)-S1th*xth[gsid]*(1024-dTth+M+60)));     //wth[i]>0;
	      Energy=(S2th-S1th+S1th*(2*(M+80))/(2*Lth-M))/M;
	      Energy=Energy*asth[gsid]+bsth[gsid];
	    }
	  
         else {Energy=0;}
*/
	 
#endif


#if(0)
          /* overwrite the energy with just the simple energy */
          /* so Mike's calibration coefficients can be used */
          /* NO PZ CORRECTION */

          Energy = ((float) (DGSEvent[i].sum2) - (float) (DGSEvent[i].sum1)) / 350.0 * ehigain[gsid] + ehioffset[gsid];
#endif

          /* dopler correct the energy */

          Energy_nodop = Energy;
          if (Pars.beta != 0)
            {
              d1 = angtheta[gsid - 1] / 57.29577951;
              Energy = Energy * (1 - Pars.beta * cos (d1)) / sqrt (1 - Pars.beta * Pars.beta);
            }
          DGSEvent[i].ehi = Energy;
          DGSEvent[i].ehi_nodop = Energy_nodop;


          d1 = (DGSEvent[i].sum2 - DGSEvent[i].sum1) / M;       /* uncorrected energy */
          d2 = DGSEvent[i].sum1 / M;    /* ~ baseline */
          if (gsid <= NGE)
            if (d1 < (double) 8192)
              if (d1 < (double) 8192)
                {
//                 e2_e1vse1[gsid]->Fill (d1, d2, 1);
                  d1 = DGSEvent[i].ehi; /* SZ corrected */
//                SZe2_e1vse1[gsid]->Fill (d1, d2, 1);
                };


          /* fill SZ's pol zero (raw) spectrum */

          top = DGSEvent[i].m2end - DGSEvent[i].m1end;
          bot = DGSEvent[i].m2begin - DGSEvent[i].m1begin;
          r1 = top / bot;
          if (r1 > 0 && r1 < 2.0);
          pzraw->Fill (gsid, r1);

#if(th)        
          Idth = gsid;
	  Bsth= DGSEvent[i].base_sample;
	  //Tth= DGSEvent[i].event_timestamp;
          Tth =  DGSEvent[i].event_timestamp & 0x3FFFFFFF; 
	  TPth=  DGSEvent[i].last_disc_timestamp;
	  if(Tth<TPth){Tth=Tth+1073741823;}  //1073741823=3FFFFFFF
	  //          Eth = (DGSEvent[i].sum2- DGSEvent[i].sum1)/M;
	  //	  Tth=(Tth-893713441670)/100000000;
	  Eth = Energy;
	  //          Blth = DGSEvent[i].sum1/M;
          S1th = DGSEvent[i].sum1;
          S2th = DGSEvent[i].sum2;
      	  dTth= Tth-TPth;
	  S0th=  DGSEvent[i].sampled_baseline;
	  S0th=S0th*64;
	  S1bth=DGSEvent[i].m1begin;
          S1eth=DGSEvent[i].m1end;
	  S2bth=DGSEvent[i].m2begin;
          S2eth=DGSEvent[i].m2end;
	  //         Lth = (S0th*M*M+S1th*1024*1024)/(2*(S0th*M-S1th*1024));  //Lth uses the value calculated before.
	  //	  Lbth=(M*M)/(2*(M-S1th/S1bth));
	  //  Leth=(2*M*(S1th/S1eth)-M*M)/(2*(S1th/S1eth-M));
	  /*  
	    for(j=0;j<1025;j++)
	    {
	      Trath[j]= DGSEvent[i].trace[j];
	      Tsth[j]=j;
	    }
	  */
	  //	  if(dTth<1434)     //for 3.5us
	  //         if(dTth<1280)     //for 2us
	  //	P1th=(S1eth-BLavfth)/(S1bth-BLavfth);
	  //	P2th=(S2eth-BLavfth)/(S2bth-BLavfth);
	  //     	P1th=(S1th/M-BLavfth)/(S1bth-BLavfth);
	  //	P2th=(S1eth-BLavfth)/(S1th/M-BLavfth);
	  //	PLth=(P1th+P2th)/2;
	  	PLth=(S1eth-BLavfth)/(S1bth-BLavfth);
	      tr.Fill();
	      Nth++;         
#endif


	      //	  	  hdT->Fill(dTth);
	      //  hdT_E->Fill(dTth,Energy);
	  //         h2_ebl->Fill (Energy, DGSEvent[i].sum1/M);
		  //        if (gsid == 105)
		  if (gsid == 110)	  
            {
	      //              h2_ebl->Fill (Energy, ave_base[gsid]);
	      //             Energy =  (DGSEvent[i].sum2- DGSEvent[i].sum1)/M;
	      //              h2_ebl->Fill (Energy, DGSEvent[i].sum1/M);
              h2_blts->Fill ((DGSEvent[i].event_timestamp - EvTimeStam0) / 1E8, ave_base[gsid]);
              h2_sblts->Fill ((DGSEvent[i].event_timestamp - EvTimeStam0) / 1E8, DGSEvent[i].sampled_baseline);
            }


          /* do the Compton suppression by setting a flag */

          for (j = 0; j < ng; j++)
            {
              if (DGSEvent[j].tpe == BGO && DGSEvent[j].tid == gsid)
                {               // BGO & GE in coincidence
                  tdiff = (int) (DGSEvent[i].event_timestamp - DGSEvent[j].event_timestamp);
                  hGeBGO_DT->Fill (tdiff, gsid);
                  if (abs (tdiff) <= 50)
                    DGSEvent[i].flag = 1;       // Mark as Dirty Ge
                }
            }
        }
      if (DGSEvent[i].tpe == BGO)
        {
          hBGOCounter->Fill ((int) ((DGSEvent[0].event_timestamp - EvTimeStam0) / 100000000), DGSEvent[i].tid);
          DGSEvent[i].ehi = (float) (DGSEvent[i].sum2) - (float) (DGSEvent[i].sum1);
        }
    }

  /* Energy Histogram loop */

  for (i = 0; i < ng; i++)
    {
      if (DGSEvent[i].tpe == GE)
        {
          e = (int) DGSEvent[i].ehi;
          if (e > 0 && e < LENSP)
            {
              gsid = DGSEvent[i].tid;
              hEhiRaw->Fill (e, gsid);
	      //     	      hBLRaw->Fill ((DGSEvent[i].sum1/M*DGSEvent[i].sum1/M-DGSEvent[i].m1begin*DGSEvent[i].m1end)/(2*DGSEvent[i].sum1/M-DGSEvent[i].m1begin-DGSEvent[i].m1end), gsid);
	      //	      hBLRaw->Fill (((DGSEvent[i].sum1/M*DGSEvent[i].sum1/M-DGSEvent[i].m1begin*DGSEvent[i].m1end)/(2*DGSEvent[i].sum1/M-DGSEvent[i].m1begin-DGSEvent[i].m1end)+(DGSEvent[i].sum2/M*DGSEvent[i].sum2/M-DGSEvent[i].m2begin*DGSEvent[i].m2end)/(2*DGSEvent[i].sum2/M-DGSEvent[i].m2begin-DGSEvent[i].m2end))/2, gsid);
              if (DGSEvent[i].flag == 0)
                {
                  hEhiCln->Fill (e, gsid);
                  hEhiCln_nodop->Fill (DGSEvent[i].ehi_nodop, gsid);
                };
              if (DGSEvent[i].flag == 1)
                hEhiDrty->Fill (e, gsid);
            };
        }
    }

/* debug list the gamma rays we found */

  if (Pars.CurEvNo <= Pars.NumToPrint)
    for (i = 0; i < ng; i++)
      {
        printf ("we have %i gamma rays\n", ng);
        printf ("%2i> ", i);
        printf ("chan_id=%i; ", DGSEvent[i].chan_id);
        printf ("board_id=%i; ", DGSEvent[i].board_id);
        printf ("id =%i; ", DGSEvent[i].id);
        printf ("tpe=%i; ", DGSEvent[i].tpe);
        printf ("tid=%i; ", DGSEvent[i].tid);
        printf ("EventTS=%llu; ", DGSEvent[i].event_timestamp);
        printf ("ehi=%8i ", DGSEvent[i].ehi);
        printf ("\n");
        fflush (stdout);
      };

  /* done */

  if (Pars.CurEvNo <= Pars.NumToPrint)
    printf ("exit bin_dgs\n");

  return (0);

}

//////////////////////////////////////////////////////////////////////////////
TH2F *
make2D (const char *txt, int xln, int xlo, int xhi, int yln, int ylo, int yhi)
{
  char str[STRLEN];
  char fn[STRLEN];

  TH2F *mkTH2F (char *, char *, int, double, double, int, double, double);


  TH2F *h2D;

  sprintf (str, txt);
  h2D = mkTH2F (str, str, xln, xlo, xhi, yln, ylo, yhi);

  return h2D;
}

//////////////////////////////////////////////////////////////////////////////
TH1D *
make1D (const char *txt, int xln, int xlo, int xhi)
{
  char str[STRLEN];
  char fn[STRLEN];
  double xlod, xhid;
  TH1D *mkTH1D (char *, char *, int, double, double);
  TH1D *h1D;

  xlod = xlo;
  xhid = xhi;

  sprintf (str, txt);
  h1D = mkTH1D (str, str, xln, xlod, xhid);
  return h1D;
}

//////////////////////////////////////////////////////////////////////////////
void
getcal (char *file)
{
  int i, ii, j, k, l, ret = 0;
  float a, b, c, d, e,tp[9];
  char mystring[1000];
  FILE *fp;

  /* get pol zero */


  fp = fopen ("dgs_pz_6.cal", "r");       // read mode
  if (fp == NULL)
    {
      printf ("could not open the DGS cal file: dgs_pz.cal; use default == 1\n", file);
    }
  else
    {

      // read file and parse

      while (fgets (mystring, 100, fp) != NULL)
        {
          ret = sscanf (mystring, "%i %f ", &i, &b);
          PZ[i] = b;
//      printf ("ge %3i has pz of %8.4f\n", i, PZ[i]);
        }
      fclose (fp);
    };


  /* get energy cal file */

  fp = fopen ("dgs_ehi_6.cal", "r");      // read mode
  if (fp == NULL)
    {
      printf ("could not open the DGS cal file: dgs_ehi.cal; set 0 offset and 1.0 gain\n");
      for (i = 1; i <= 110; i++)
        {
          ehigain[i] = 1;
          ehioffset[i] = 0;
        };
    }
  else
    {

      // read file and parse

      while (fgets (mystring, 100, fp) != NULL)
        {
          ret = sscanf (mystring, "%i %f %f ", &i, &c, &d);
          ehigain[i] = d;
          ehioffset[i] = c;
          printf ("ge %3i has offset and gain of: %8.2f %8.4f and a PZ of %8.4f\n", i, ehioffset[i], ehigain[i], PZ[i]);
        }
      fclose (fp);
    };

  /*  
 fp = fopen ("Para.cal", "r");      // read mode
  if (fp == NULL)
    {
      printf ("could not open the Para.cal file\n");
      return;
    }
  else
    {

      // read file and parse
      fgets (mystring, 255, fp);
      fgets (mystring, 255, fp);

      while (fgets (mystring, 100, fp) != NULL)
        {
          ret = sscanf (mystring, "%i %i %f %f %f", &i, &c, &a, &b,&d);
          ath[i] = a;
          bth[i] = b;
	  slth[i]=d;
          printf ("ge %3i has offset and gain of: %8.2f %8.4f and a SL of %8.4f\n", i, bth[i], ath[i], slth[i]);
        }
      fclose (fp);
    }
  */  
  
  /*  
  //       fp = fopen ("Lt_t.cal", "r");      // read mode  //for 2us and dT>1280
             fp = fopen ("Lt3_t.cal", "r");     //for 3.5us and dT>1434
  if (fp == NULL)
    {
      printf ("could not open the Para.cal file\n");
      return;
    }
  else
    {

      // read file and parse
      fgets (mystring, 255, fp);
      fgets (mystring, 255, fp);

      while (fgets (mystring, 100, fp) != NULL)
        {
          ret = sscanf (mystring, "%i %f %f", &i, &a, &b);
          ath[i] = a;
          bth[i] = b;
	  //	  slth[i]=d;
          printf ("ge %3i has offset and gain of: %8.2f %8.4f\n", i, bth[i], ath[i]);
        }
      fclose (fp);
    }

  //       fp = fopen ("P_2SL.cal", "r");              //for 2us and dT<1280
         fp = fopen ("P_3SL.cal", "r");                //for 3.5us and dT<1434
  if (fp == NULL)
    {
      printf ("could not open the Para.cal file\n");
      return;
    }
  else
    {

      // read file and parse
      fgets (mystring, 255, fp);

      while (fgets (mystring, 100, fp) != NULL)
        {
          ret = sscanf (mystring, "%i %f %f %f %f", &i, &a, &b, &c, &d);
          asth[i] = a;
          bsth[i] = b;
	  wth[i]=c;
          xth[i]=d;
          printf ("ge %3i has offset and gain of: %8.2f %8.4f\n", i, bsth[i], asth[i]);
        }
      fclose (fp);
    }   
*/
  
  //       fp = fopen ("Ge3L_002.cal", "r");
  fp = fopen ("Ge_Md.cal", "r");
  //  	        fp = fopen ("Ge158.cal", "r");              //for 2us 
  //    fp = fopen ("Ge185.cal", "r");
  //          fp = fopen ("Ge4504.cal", "r");                   // for 2us
  //        fp = fopen ("Ge503.cal", "r");                //for 3.5us 
  if (fp == NULL)
    {
      printf ("could not open the Para.cal file\n");
      return;
    }
  else
    {

      // read file and parse
      fgets (mystring, 255, fp);
      //      fgets (mystring, 255, fp);
      while (fgets (mystring, 100, fp) != NULL)
        {
          ret = sscanf (mystring, "%i %f %f %f %f", &i, &a, &b, &c, &d);
          Pth[i] = a;
          BLPth[i] = b;
          ath[i] = c;
	  bth[i] = d;
          printf ("ge %3i has offset and gain of: %8.2f %8.4f\n", i, bth[i], ath[i]);
        }
      fclose (fp);
    }

  
    fp = fopen ("dT_P2.cal", "r");
  if (fp == NULL)
    {
      printf ("could not open the Para.cal file\n");
      return;
    }
  else
    {

      // read file and parse
      fgets (mystring, 255, fp);
      //      fgets (mystring, 255, fp);
      for(ii=1;ii<=110;ii++)
        {
          fscanf (fp,"%i %f %f %f %f %f %f %f %f %f", &i, &tp[0], &tp[1], &tp[2], &tp[3], &tp[4], &tp[5], &tp[6], &tp[7], &tp[8]);
	  for(j=0;j<=8;j++) {Pdtth[i][j] = tp[j];}
          printf ("Id %3i tp[0] %8.6f  tp[1] %8.6f  tp[2] %8.6f  tp[3] %8.6f  tp[4] %8.6f  tp[5] %8.6f tp[6] %8.6f tp[7] %8.6f tp[8] %8.6f\n", i, tp[0], tp[1], tp[2], tp[3], tp[4], tp[5], tp[6], tp[7], tp[8]);
        }
      fclose (fp);
    }


}




#if(0)
/*--------------------------------------------------------*/
void
SetBeta ()
{

  /* delarations */

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

  /* done */

}
#endif

/*-----------------------------------------------------------*/
