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
#include "TH3.h"
#include "TObjArray.h"
#include "TObject.h"
#include "TKey.h"
#include "TSystem.h"
#include "TCutG.h"
#include "TTree.h"
#include "gdecomp.h"
#include "veto_pos.h"
#include "GEBSort.h"
#include "GTMerge.h"


#define MAXPAD 1024


float egamBinWidth;

/* pointers to ROOT spectra */
/* __spectra affected by TS+T0 time gating */

TH1D *CCsum_nodop;

TH1D *CCsum;
TH1D *CCadd;
TH2F *SEGadd_FPGA;
TH2F *SEGadd_PSA;
TH2F *SEGsum_FPGA;
TH2F *SEGsum_PSA;

TH1D *radius_all;
TH2F *CCe;
TH2F *ggCC;
TH2F *SMAP_allhits;
TH2F *SMAP_allhits_sq;
TH1D *hitpat;
TH1D *numHitArray2_sp;
TH1D *CCmult;
TH1D *rate_mode2_min;
TH1D *rate_mode2_sec;
TH1D *rate_mode2_intp;
TH2F *evsr_all;
TH2F *z_plot;
TH2F *seghitpat;
TH1D *T0;
TH3F *crystal_hit;
TH1D *hit2;
TH1D *hit3;
TH1D *hit4;
TH2F *cdtbtev;


TH2F *crXcr;
TH2F *SegE;

/* parameters */

extern PARS Pars;
extern TH1D *ehi[MAXDETPOS + 1];

/* for mean polar and azimuth angles */

double pol[MAXDETPOS + 1];
double azi[MAXDETPOS + 1];
long long int ndethits[MAXDETPOS + 1];
long long int ClastTS[MAXDETNO + 1];

unsigned int modhit[MAXMODNO + 1][MAXCRYSTALNO + 1];

extern unsigned int *veto_cube;

double meanCCenergy_sum = 0;
long long int nCCenergy_sum = 0;
double meanCCenergy_add = 0;
long long int nCCenergy_add = 0;

/* counters */

long long int pad_sp[MAXPAD];


/* ----------------------------------------------------------------- */

int
sup_mode2 ()
{
  /* declarations */

  char str1[STRLEN], str2[STRLEN];
  float pi;
  int i, j, i1;

  TH1D *mkTH1D (char *, char *, int, double, double);
  TH2F *mkTH2F (char *, char *, int, double, double, int, double, double);
  TH3F *mkTH3F (char *, char *, int, double, double, int, double, double, int, double, double);

  /* initialize */

  for (i = 0; i <= MAXDETPOS; i++)
    {
      pol[i] = 0;
      azi[i] = 0;
      ndethits[i] = 0;
    };

  for (i = 0; i <= MAXMODNO; i++);
  for (j = 0; j <= MAXCRYSTALNO; j++)
    modhit[i][j] = 0;

  for (i = 0; i <= MAXDETNO; i++)
    ClastTS[i] = 0;

  for (i = 0; i < MAXPAD; i++)
    pad_sp[i] = 0;

  /* define spectra */

  sprintf (str1, "T0");
  sprintf (str2, "T0");
  T0 = mkTH1D (str1, str2, 2048, 0, 150);
  sprintf (str1, "");
  T0->SetXTitle ("10 nsec");

  sprintf (str1, "hitpat");
  sprintf (str2, "hitpat");
  hitpat = mkTH1D (str1, str2, 201, 0, 200);
  sprintf (str1, "det number");
  hitpat->SetXTitle (str1);

  sprintf (str1, "numHitArray2");
  sprintf (str2, "numHitArray2");
  numHitArray2_sp = mkTH1D (str1, str2, 200, 1, 200);
  sprintf (str1, "# array interactions");
  numHitArray2_sp->SetXTitle (str1);

  sprintf (str1, "CCmult");
  sprintf (str2, "CC Multiplicity");
  CCmult = mkTH1D (str1, str2, 21, 0, 20);
  sprintf (str1, "CCmult");
  CCmult->SetXTitle (str1);

  sprintf (str1, "CCsum");
  sprintf (str2, "CCsum");
  CCsum = mkTH1D (str1, str2, LONGLEN, 1, LONGLEN);
  sprintf (str1, "(keV)");
  CCsum->SetXTitle (str1);

  sprintf (str1, "hit2");
  sprintf (str2, "hit2");
  hit2 = mkTH1D (str1, str2, LONGLEN, 1, LONGLEN);
  sprintf (str1, "(keV)");
  CCsum->SetXTitle (str1);

  sprintf (str1, "hit3");
  sprintf (str2, "hit3");
  hit3 = mkTH1D (str1, str2, LONGLEN, 1, LONGLEN);
  sprintf (str1, "(keV)");
  CCsum->SetXTitle (str1);

  sprintf (str1, "hit4");
  sprintf (str2, "hit4");
  hit4 = mkTH1D (str1, str2, LONGLEN, 1, LONGLEN);
  sprintf (str1, "(keV)");
  CCsum->SetXTitle (str1);


  sprintf (str1, "CCsum_nodop");
  sprintf (str2, "CCsum_nodop");
  CCsum_nodop = mkTH1D (str1, str2, LONGLEN, 1, LONGLEN);
  sprintf (str1, "(keV)");
  CCsum_nodop->SetXTitle (str1);

  sprintf (str1, "CCadd");
  sprintf (str2, "CCadd");
  CCadd = mkTH1D (str1, str2, LONGLEN, 1, LONGLEN);
  sprintf (str1, "(keV)");
  CCadd->SetXTitle (str1);

  sprintf (str1, "SEGadd_FPGA");
  sprintf (str2, "SEGadd_FPGA");
  SEGadd_FPGA = mkTH2F (str1, str2, MAXDETNO + 1, 0, MAXDETNO, LONGLEN, 1, LONGLEN);
  sprintf (str1, "crystal #");
  SEGadd_FPGA->SetXTitle (str1);
  sprintf (str1, "SEGadd_FPGA");
  SEGadd_FPGA->SetYTitle (str1);

  sprintf (str1, "SEGadd_PSA");
  sprintf (str2, "SEGadd_PSA");
  SEGadd_PSA = mkTH2F (str1, str2, MAXDETNO + 1, 0, MAXDETNO, LONGLEN, 1, LONGLEN);
  sprintf (str1, "crystal #");
  SEGadd_PSA->SetXTitle (str1);
  sprintf (str1, "SEGadd_PSA");
  SEGadd_PSA->SetYTitle (str1);

  sprintf (str1, "SEGsum_FPGA");
  sprintf (str2, "SEGsum_FPGA");
  SEGsum_FPGA = mkTH2F (str1, str2, MAXDETNO + 1, 0, MAXDETNO, LONGLEN, 1, LONGLEN);
  sprintf (str1, "crystal #");
  SEGsum_FPGA->SetXTitle (str1);
  sprintf (str1, "SEGsum_FPGA");
  SEGsum_FPGA->SetYTitle (str1);

  sprintf (str1, "SEGsum_PSA");
  sprintf (str2, "SEGsum_PSA");
  SEGsum_PSA = mkTH2F (str1, str2, MAXDETNO + 1, 0, MAXDETNO, LONGLEN, 1, LONGLEN);
  sprintf (str1, "crystal #");
  SEGsum_PSA->SetXTitle (str1);
  sprintf (str1, "SEGsum_PSA");
  SEGsum_PSA->SetYTitle (str1);

  sprintf (str1, "radius_all");
  sprintf (str2, "radius (all points)");
  radius_all = mkTH1D (str1, str2, 4096, RMIN, RMAX);
  radius_all->SetXTitle (str1);

  sprintf (str1, "rate_mode2_min");
  sprintf (str2, "rate_mode2_min headers");
  rate_mode2_min = mkTH1D (str1, str2, RATELEN, 0, RATELEN);
  rate_mode2_min->SetXTitle (str1);
  rate_mode2_min->SetXTitle ("minutes");
  rate_mode2_min->SetYTitle ("Hz");

  sprintf (str1, "rate_mode2_sec");
  sprintf (str2, "rate_mode2_sec headers");
  rate_mode2_sec = mkTH1D (str1, str2, RATELEN, 0, RATELEN);
  rate_mode2_sec->SetXTitle (str1);
  rate_mode2_sec->SetXTitle ("seconds");
  rate_mode2_sec->SetYTitle ("Hz");


  sprintf (str1, "rate_mode2_intp");
  sprintf (str2, "rate_mode2 interaction points");
  rate_mode2_intp = mkTH1D (str1, str2, RATELEN, 0, RATELEN);
  rate_mode2_intp->SetXTitle (str1);
  rate_mode2_intp->SetXTitle ("minutes");
  rate_mode2_intp->SetYTitle ("Hz");

  sprintf (str1, "cdtbtev");
  sprintf (str2, "cdtbtev");
  cdtbtev = mkTH2F (str1, str2, DTBTEVLEN, 0, DTBTEVLEN, MAXDETNO, 0, MAXDETNO - 1);
  sprintf (str1, "dt [usec]");
  cdtbtev->SetXTitle (str1);
  sprintf (str1, "detector number");
  cdtbtev->SetYTitle (str1);


  /* star map of GRETA */

  sprintf (str1, "evsr_all");
  sprintf (str2, "evsr_all");
  evsr_all = mkTH2F (str1, str2, MEDIUMLEN, RMIN, RMAX, MEDIUMLEN, 1, MEDIUMLEN);
  sprintf (str1, "energy");
  evsr_all->SetXTitle (str1);
  sprintf (str1, "radius");
  evsr_all->SetYTitle (str1);

  sprintf (str1, "SMAP_allhits");
  sprintf (str2, "SMAP_allhits");
  SMAP_allhits = mkTH2F (str1, str2, 720, -180, 180, 360, 0, 180);
  sprintf (str1, "Azimuth");
  SMAP_allhits->SetXTitle (str1);
  sprintf (str1, "Polar");
  SMAP_allhits->SetYTitle (str1);

  sprintf (str1, "SMAP_allhits_sq");
  sprintf (str2, "SMAP_allhits_sq");
  SMAP_allhits_sq = mkTH2F (str1, str2, 720, -180, 180, 360, 0, 180);
  sprintf (str1, "Azimuth");
  SMAP_allhits_sq->SetXTitle (str1);
  sprintf (str1, "Polar");
  SMAP_allhits_sq->SetYTitle (str1);

  sprintf (str1, "seghitpat");
  sprintf (str2, "seghitpat");
  seghitpat = mkTH2F (str1, str2, MAXDETNO, 1, MAXDETNO, 38, 0, 37);
  sprintf (str1, "detector");
  seghitpat->SetXTitle (str1);
  sprintf (str1, "segment");
  seghitpat->SetYTitle (str1);

  sprintf (str1, "CCe");
  sprintf (str2, "CCe");
  CCe = mkTH2F (str1, str2, MAXDETPOS + 1, 0, MAXDETPOS, LONGLEN, 1, LONGLEN);
  sprintf (str1, "crystal #");
  CCe->SetXTitle (str1);
  sprintf (str1, "cc energy");
  CCe->SetYTitle (str1);

  sprintf (str1, "ggCC");
  sprintf (str2, "CC gg matrix");
  ggCC = mkTH2F (str1, str2, Pars.GGMAX, 1, Pars.GGMAX, Pars.GGMAX, 1, Pars.GGMAX);
  sprintf (str1, "g1");
  ggCC->SetXTitle (str1);
  sprintf (str1, "g2");
  ggCC->SetYTitle (str1);

  sprintf (str1, "z_plot");
  sprintf (str2, "z_plot");
  z_plot = mkTH2F (str1, str2, NUMAGATAPOS + 1, 0, NUMAGATAPOS, 1000, 0, 100);
  sprintf (str1, "crystal");
  z_plot->SetXTitle (str1);
  sprintf (str1, "z values");
  z_plot->SetYTitle (str1);

  sprintf (str1, "crystal_hit");
  crystal_hit = mkTH3F (str1, str1, 182, -45, 45, 182, -45, 45, 192, 0, 95);

  sprintf (str1, "SegE");
  sprintf (str2, "SegE");
  SegE = mkTH2F (str1, str2, MAXSEGNO + 1, 0, MAXSEGNO, SHORTLEN, 1, LONGLEN);
  SegE->SetXTitle ("Crystal#*36+seg#");
  SegE->SetYTitle ("segment energy");



  sprintf (str1, "crXcr");
  sprintf (str2, "crXcr");
  crXcr = mkTH2F (str1, str2, MAXDETNO + 1, 0, MAXDETNO, MAXDETNO + 1, 0, MAXDETNO);
  sprintf (str1, "crystal 1");
  crXcr->SetXTitle (str1);
  sprintf (str1, "crystal 2");
  crXcr->SetYTitle (str1);

  /* list what we have */

//  Pars.wlist = gDirectory->GetList ();
//  Pars.wlist->Print ();

  for (i = 1; i <= MAXDETPOS; i++)
    Pars.detpolang[i] = 0;
  printf ("MAXDETPOS=%i\n", MAXDETPOS);
  fflush (stdout);

}

/* ----------------------------------------------------------------- */

int
exit_bin_mode2 (void)
{

/* declarations*/

  int i, j, mod, crystal, ncrystals = 0;
  double csum, mcsum, d1, d2;
  int k, l, m;
  unsigned int ibad, icsum, nbadcubes = 0, ntotcubes = 0, lli;
  FILE *fp;
  unsigned int ui1;
  float xxi, yyi, zzi, xxj, yyj, zzj, rpol, dotProduct, minpol, maxang;
  int crystalNumber, holeNum, side = 0, noside = 0, nn = 0;

  /* normalize and report */

  printf ("\n");
  printf ("mean observed polar and azimuth angles\n");
  printf ("for the detectors we saw (module,crystal)\n");
  printf ("\n");
  printf ("------------------------------------------\n");

  for (i = 0; i <= MAXDETPOS; i++)
    if (ndethits[i] > 0)
      {
        ncrystals++;
        pol[i] /= ndethits[i];
        azi[i] /= ndethits[i];
        printf ("mean pol/azi angle for detector %3i are %7.2f/%7.2f\n", i,
                pol[i] * 180.0 / M_PI, azi[i] * 180.0 / M_PI);
      };

  printf ("we found %i crystals\n", ncrystals);
  printf ("------------------------------------------\n");
  printf ("\n");

  minpol = 1000000;
  for (i = 0; i <= MAXDETPOS; i++)
    for (j = i + 1; j <= MAXDETPOS; j++)
      if (ndethits[i] > 10 && ndethits[j] > 10)
        {
          if (Pars.AGATA_data)
            {
              xxi = Pars.TrX[i];
              yyi = Pars.TrY[i];
              zzi = Pars.TrZ[i];
//          printf("i %i:: %f,%f,%f\n",i,xxi,yyi,zzi);
              xxj = Pars.TrX[j];
              yyj = Pars.TrY[j];
              zzj = Pars.TrZ[j];
//          printf("j %i:: %f,%f,%f\n",j,xxj,yyj,zzj);
              maxang = 19 / 57.2958;
            }

          else
            {

              crystalNumber = i % 4;
              holeNum = holeNum = i / 4;
              xxi = Pars.crmat[holeNum][crystalNumber][0][3];
              yyi = Pars.crmat[holeNum][crystalNumber][1][3];
              zzi = Pars.crmat[holeNum][crystalNumber][2][3];
//          printf("%i::%i,%i: %f,%f,%f    ",i,bin_mode2holeNum,crystalNumber,xxi,yyi,zzi);

              crystalNumber = j % 4;
              holeNum = j / 4;
              xxj = Pars.crmat[holeNum][crystalNumber][0][3];
              yyj = Pars.crmat[holeNum][crystalNumber][1][3];
              zzj = Pars.crmat[holeNum][crystalNumber][2][3];
//          printf("%i::%i,%i: %f,%f,%f\n",j,holeNum,crystalNumber,xxj,yyj,zzj);

              maxang = 25 / 57.2958;
            }
          d1 = xxi * xxi + yyi * yyi + zzi * zzi;
          d1 = sqrt (d1);
          d2 = xxj * xxj + yyj * yyj + zzj * zzj;
          d2 = sqrt (d2);
//          printf("d1=%f,d2=%f\n", d1,d2);
          dotProduct = xxi * xxj + yyi * yyj + zzi * zzj;
          dotProduct /= d1;
          dotProduct /= d2;
          if (dotProduct < -1.0)
            dotProduct = -1.0;
          if (dotProduct > 1.0)
            dotProduct = 1.0;
          rpol = acosf (dotProduct);

          nn++;
//          printf ("%3i:: %3i %3i %6.2f ; ", nn, i, j, rpol * 57.2958);
          if (rpol < maxang)
            {
              side++;
//              printf (" *\n");
            }
          else
            {
              noside++;
//              printf ("\n");
            }

          if (rpol < minpol)
            minpol = rpol;
        }
  printf ("minimum pol angle between crystals is %6.2f\n", minpol * 57.2958);
  d1 = 100 * 2 * (float) side / (6.0 * ncrystals);
  printf ("compactness = %6.2f%%, %i of %i\n", d1, side, 6 * ncrystals);

  if (Pars.vetoSpots)
    {

      printf ("we may have to evaluate as a fuction of Z\n");
      printf ("rather than base on mean from entire crystal...\n");
      printf ("\n");


      fp = fopen (Pars.vetoSpotsfn, "w");
      if (fp == NULL)
        {
          printf ("could not open %s\n", Pars.vetoSpotsfn);
          exit (1);
        }
      else
        printf ("open %s for bad cube listings\n", Pars.vetoSpotsfn);

      for (i = 0; i <= MAXGTMODNO; i++)
        for (j = 0; j <= MAXCRYSTALNO; j++)
          {

#if(0)
            csum = 0;
            for (k = 0; k <= VETO_NX; k++)
              for (l = 0; l <= VETO_NY; l++)
                for (m = 0; m <= VETO_NZ; m++)
                  csum += *(veto_cube + VETO_INDX (i, j, k, l, m));
#else
            csum = modhit[i][j];
#endif

            if (csum > 0)
              {
                mcsum = csum / VETO_NX / VETO_NY / VETO_NZ;
                printf ("mod %2i, crystal %1i, has %6.2f counts, mean %f, ", i, j, csum, (float) mcsum);
                fflush (stdout);
                ntotcubes += VETO_NX * VETO_NY * VETO_NZ;
              };

            if (csum > 0)
              {
                icsum = (int) (Pars.vetocutfac * mcsum + 0.5);
                ibad = 0;

                /* now write out the bad cubes only */

                for (k = 0; k <= VETO_NX; k++)
                  for (l = 0; l <= VETO_NY; l++)
                    for (m = 0; m <= VETO_NZ; m++)
                      if ((int) *(veto_cube + VETO_INDX (i, j, k, l, m)) > icsum)
                        {
                          ibad++;
                          nbadcubes++;
                          fprintf (fp, "%2i %2i %3i %3i %3i ; %3i\n", i, j, k, l, m,
                                   *(veto_cube + VETO_INDX (i, j, k, l, m)));
//                          fprintf (fp, "%5.1f %5.1f %5.1f\n", VETO_X_XEDNI (k), VETO_Y_XEDNI (l), VETO_Z_XEDNI (m));
                        };
                d1 = (double) ibad / VETO_NX / VETO_NY / VETO_NZ;
                printf ("there are %6i bad cubes; %6.2f%%\n", ibad, 100 * d1);
              };

          };
      fclose (fp);
      printf ("\n");
      printf ("bin_mode2: closed %s\n", Pars.vetoSpotsfn);
      printf ("\n");

      printf ("VETO_INDX (7, 3, 10, 10, 10) =%i\n", VETO_INDX (7, 3, 10, 10, 10));
      printf ("VETO_INDX (8, 0, 10, 10, 10) =%i\n", VETO_INDX (8, 0, 10, 10, 10));
      printf ("VETO_INDX (8, 1, 10, 10, 10) =%i\n", VETO_INDX (8, 1, 10, 10, 10));
      printf ("VETO_INDX (8, 2, 10, 10, 10) =%i\n", VETO_INDX (8, 2, 10, 10, 10));
      printf ("VETO_INDX (8, 3, 10, 10, 10) =%i\n", VETO_INDX (8, 3, 10, 10, 10));

      ui1 = MAXGTMODNO * MAXCRYSTALNO * VETO_NX * VETO_NY * VETO_NZ;
      printf ("max index...: %10u\n", ui1);
      printf ("ntotcubes...= %10i\n", ntotcubes);
      d1 = 100 * (double) nbadcubes / ntotcubes;
      printf ("nbadcubes...= %10i or %6.2f%%\n", nbadcubes, d1);
      printf ("VETO_X_D = %f\n", VETO_X_D);
      printf ("VETO_Y_D = %f\n", VETO_Y_D);
      printf ("VETO_Z_MAX = %f\n", VETO_Z_MAX);
      printf ("VETO_BINWIDTH = %f mm\n", VETO_BINWIDTH);
      printf ("VETO_NX = %i\n", VETO_NX);
      printf ("VETO_NY = %i\n", VETO_NY);
      printf ("VETO_NZ = %i\n", VETO_NZ);
      printf ("\n");
    };


  printf ("module hit statistics\n");
  for (i = 0; i <= MAXMODNO; i++)
    for (j = 0; j <= MAXCRYSTALNO; j++)
      if (modhit[i][j] > 0)
        printf ("module %3i, crystal %i was hit %10i times\n", i, j, modhit[i][j]);
  printf ("\n");

  /* report pad distribution */

  printf ("bin_mode2: pad distribution\n");
  lli = 0;
  for (i = 0; i < MAXPAD; i++)
    lli += pad_sp[i];
  for (i = 0; i < MAXPAD; i++)
    if (pad_sp[i] > 0)
      printf ("__id %3i, %12lli, %5.2f%%\n", i, pad_sp[i], 100 * (double) pad_sp[i] / lli);
  printf ("\n");

  /* report mean CC energy */

  meanCCenergy_sum /= nCCenergy_sum;
  printf ("bin_mode2 calculated mean CCsum energy to be %8.3f keV\n", (float) meanCCenergy_sum);
  meanCCenergy_add /= nCCenergy_add;
  printf ("bin_mode2 calculated mean CCadd energy to be %8.3f keV\n", (float) meanCCenergy_add);


  return (0);
}

/* ----------------------------------------------------------------- */

int
bin_mode2 (GEB_EVENT * GEB_event)
{

  /* declarations */

  int i, j, k, crystalno, moduleno, detno, nn = 0;
  float sX, sY, polAng, aziAng, rr, xx, yy, zz, r1, rmax, rmin;
  double d1, d2;
  static int nev = 0, xnev = 0;
  static float sumcc = 0, addcc = 0, mcc;
  char str[128];
  int GebTypeStr (int type, char str[]);
  float detDopFac, dp, addedEnergy = 0, addedSEGEnergy = 0, r2;
  float addedSEGEnergy_FPGA[MAXDETNO], addedSEGEnergy_PSA[MAXDETNO];
  int segProcessed[MAXDETNO][NUMSEGS];
  float RAD2DEG = 0.0174532925;
  float CCenergies[MAX_GAMMA_RAYS];
  static int firsttime = 1;
  static long long int t0;
  float polang[MAX_GAMMA_RAYS];
  float doppler_factor[MAX_GAMMA_RAYS];
  float xar[MAXCOINEV], yar[MAXCOINEV], zar[MAXCOINEV];
  int detectorPosition, crystalNumber;
  int i1, i2, i3, ok;
  int nCCenergies, crystaltype;
  static int nperrors = 0;
  unsigned int ui1;
  int numHitArray2 = 0, ntracked = 0;
  float fom = 0;
  int crhit[MAXDETNO], ncr;
  long long int lli;
  int segHit[MAXDETNO + 1][TOT_SEGS + 1];

  TRACKED_GAMMA_HIT *grh;
  CRYS_INTPTS *ptinp;
  GEBDATA *ptgd;

  /* prototypes */

  float findAzimuthFromCartesian (float, float, float);
  float findPolarFromCartesian (float, float, float, float *);

  if (Pars.CurEvNo <= Pars.NumToPrint)
    printf ("entered bin_mode2: %i/%i\n", Pars.CurEvNo, Pars.NumToPrint);

  /* scan for tracked data */

  ntracked = 0;
  for (i = 0; i < GEB_event->mult; i++)
    {
      if (GEB_event->ptgd[i]->type == GEB_TYPE_TRACK)
        {
          grh = (TRACKED_GAMMA_HIT *) GEB_event->ptinp[i];
          ntracked += grh->ngam;
          for (j = 0; j < grh->ngam; j++)
            fom = grh->gr[j].fom;
        };
    };
//  printf("# tracked gamma rays: %i, FOM= %f\n", ntracked, fom);

  if (Pars.requiretracked)
    {

      /* require tracked data before we bin mode2 data */

      if (ntracked == 0)
        return (0);

    }


  /* bail out if the multiplicity is too high */

  assert (GEB_event->mult < MAX_GAMMA_RAYS);


  /*--------------------------------*/
  /* apply calibration coefficients */
  /*--------------------------------*/

//TBD. should be done before next step

#if(0)

  /*---------------------------------*/
  /* recalibrate the PSA energies to */
  /* the better segment energies     */
  /* THIS IS UNFINISHED!!            */
  /*---------------------------------*/

  for (i = 0; i < MAXDETNO; i++)
    for (j = 0; j < TOT_SEGS; j++)
      segHit[i][j] = 0;

  for (i = 0; i < GEB_event->mult; i++)
    if (GEB_event->ptgd[i]->type == GEB_TYPE_DECOMP)
      {
        /* cast */

        ptinp = (CRYS_INTPTS *) GEB_event->ptinp[i];

        /* extract IDs */

        crystalno = (ptinp->crystal_id & 0x0003);
        moduleno = ((ptinp->crystal_id & 0xfffc) >> 2);

        if (Pars.AGATA_data == 0)
          detno = moduleno * 4 + crystalno;
        else if (Pars.AGATA_data == 1)
          detno = moduleno * 3 + crystalno;

        /* check */

        if (detno < 0 || detno > MAXDETNO)
          {
            printf ("bad detector number: %i, set to 1\n", detno);
            detno = 1;
          };

        /* count hits in segments */

        for (k = 0; k < ptinp->num; k++)
          {
            if (ptinp->intpts[k].seg < TOT_SEGS && ptinp->intpts[k].seg >= 0)
              segHit[detno][ptinp->intpts[k].seg]++;
            else
              printf ("bad segment number: %i, k=%i\n", ptinp->intpts[k].seg, k);
          };

      };

  if (Pars.CurEvNo <= Pars.NumToPrint)
    for (i = 0; i < MAXDETNO; i++)
      for (j = 0; j < TOT_SEGS; j++)
        if (segHit[i][j] > 0)
          printf ("det/seg %3i/%2i has hits %i\n", i, j, segHit[i][j]);

#endif

  /*-----------*/
  /* main loop */
  /*-----------*/

  addedEnergy = 0;
  for (i = 0; i < MAXDETNO; i++)
    {
      addedSEGEnergy_FPGA[i] = 0;
      addedSEGEnergy_PSA[i] = 0;
      for (j = 0; j < NUMSEGS; j++)
        segProcessed[i][j] = 0;
    };
  for (i = 0; i < GEB_event->mult; i++)
    CCenergies[i] = 0;
  nCCenergies = 0;
  numHitArray2 = 0;
  ncr = 0;


  for (i = 0; i < GEB_event->mult; i++)
    {

      if (GEB_event->ptgd[i]->type == GEB_TYPE_DECOMP)
        {

          /* cast */

          ptinp = (CRYS_INTPTS *) GEB_event->ptinp[i];

          /* manipulate pad value and bin */

//          if (ptinp->pad>=128) ptinp->pad-=128;
          if ((ptinp->pad & 0x000000ff) < MAXPAD)
            pad_sp[(ptinp->pad & 0x000000ff)]++;

          if (Pars.CurEvNo <= Pars.NumToPrint)
            printf ("ptinp->pad= %i, (ptinp->pad&0x000000ff)= %i\n", ptinp->pad, (ptinp->pad & 0x000000ff));

          /* check for good pad value */

          if ((ptinp->pad & 0x000000ff) == 0)
            {

              if (Pars.CurEvNo <= Pars.NumToPrint)
                printf ("good pad value\n");

              /* number of hits in the array */

              numHitArray2 += ptinp->num;

              /* debug print */

              if (Pars.CurEvNo <= Pars.NumToPrint)
                {
                  GebTypeStr (GEB_event->ptgd[i]->type, str);
                  printf ("mode2 data: %2i> %2i, %s, TS=%lli\n", i, GEB_event->ptgd[i]->type, str,
                          GEB_event->ptgd[i]->timestamp);
                  printf ("number of hits: %i\n", ptinp->num);
//              printf("total number of hits in array: %i\n", numHitArray2);
                }

              /* extract IDs */

              crystalno = (ptinp->crystal_id & 0x0003);
              moduleno = ((ptinp->crystal_id & 0xfffc) >> 2);
              modhit[moduleno][crystalno] += ptinp->num;

              if (Pars.AGATA_data == 0)
                {
                  detno = moduleno * 4 + crystalno;
                  if (crystalno == 0 || crystalno == 2)
                    crystaltype = 0;
                  else
                    crystaltype = 1;
                }
              else if (Pars.AGATA_data == 1)
                {
                  detno = moduleno * 3 + crystalno;
                  crystaltype = 0;
                };

              /* update cdtbtev array (channel deadtime) */

              d1 = (double) GEB_event->ptgd[i]->timestamp - (double) ClastTS[detno];
              d1 /= 100;
              if (d1 >= 0 && d1 < (double) DTBTEVLEN)
                cdtbtev->Fill (d1, detno);
              ClastTS[detno] = GEB_event->ptgd[i]->timestamp;
//          printf("new ClastTS[detno]=%lli\n",ClastTS[detno]);

              /* store away hit chrystals */

              crhit[ncr] = detno;
              ncr++;

              /* 2D segment hit matrix */

              for (j = 0; j < ptinp->num; j++)
                {
                  seghitpat->Fill (detno, ptinp->intpts[j].seg, 1);
                };

              /* fill 3D crystal cube */

              if (Pars.crystalID3D_fomlo > 0 || Pars.crystalID3D_fomhi < 2.0)
                {
                  if (ntracked == 1)
                    ok = ((fom >= Pars.crystalID3D_fomlo) && (fom <= Pars.crystalID3D_fomhi));
                  else
                    ok = 0;
                }
              else
                ok = 1;

              if ((detno == Pars.crystalID3D) && ok)
                for (j = 0; j < ptinp->num; j++)
                  if (ptinp->intpts[j].e > Pars.crystalID3D_elo && ptinp->intpts[j].e < Pars.crystalID3D_ehi)
                    crystal_hit->Fill (ptinp->intpts[j].x, ptinp->intpts[j].y, ptinp->intpts[j].z);

              /* bin if enabled */

              if (Pars.enabled[detno])
                {
                  /* mode 2 rate spectrum, x=minute, y=Hz */

                  if (firsttime)
                    {
                      firsttime = 0;
                      t0 = GEB_event->ptgd[i]->timestamp;
                    };

                  /* update the rate spectra */

                  d1 = (double) (GEB_event->ptgd[i]->timestamp - t0);
                  d1 /= 100000000;

                  if (d1 > 0 && d1 < (double) RATELEN)
                    rate_mode2_sec->Fill (d1, 1);

                  d1 /= 60;
                  if (d1 > 0 && d1 < (double) RATELEN)
                    {
                      rate_mode2_min->Fill (d1, 1 / 60.0);
                      rate_mode2_intp->Fill (d1, ptinp->num / 60.0);
                    };


                  /* update veto_cube */

                  if (Pars.vetoSpots)
                    for (j = 0; j < ptinp->num; j++)
                      if (ptinp->intpts[j].e <= Pars.vetoecut)
                        {
                          i1 = VETO_X_INDEX (ptinp->intpts[j].x);
                          i2 = VETO_Y_INDEX (ptinp->intpts[j].y);
                          i3 = VETO_Z_INDEX (ptinp->intpts[j].z);

                          if (i1 >= 0 && i2 >= 0 && i3 >= 0)
                            if (i1 < VETO_NX && i2 < VETO_NY && i3 < VETO_NZ)
                              {
                                (*(veto_cube + VETO_INDX (moduleno, crystalno, i1, i2, i3)))++;
                              };
                        };


                  if (Pars.CurEvNo <= Pars.NumToPrint)
                    {
                      printf ("* %i/%i, is GEB_TYPE_DECOMP: num=%i\n", i, GEB_event->mult, ptinp->num);
                      printf ("__detno: %i, module: %i, crystalNumber: %i\n", detno, moduleno, crystalno);
                    }

                  /* calibrate mode 2 CC data */

                  ptinp->tot_e = ptinp->tot_e * Pars.CCcal_gain[detno] + Pars.CCcal_offset[detno];
                  addedEnergy += ptinp->tot_e;
                  CCenergies[nCCenergies] = ptinp->tot_e / Pars.modCCdopfac[detno];
                  if (Pars.CurEvNo <= Pars.NumToPrint)
                    printf ("CCenergies[%i]=%7.2f\n", nCCenergies, CCenergies[nCCenergies]);

#if(0)
                  xnev++;
                  printf ("\n");
                  printf ("gamma ray # %4i (loop), ", xnev);
                  printf (" %8.3f modCCdopfac[%i]=%8.3f", ptinp->tot_e, detno, Pars.modCCdopfac[detno]);
                  printf ("\n");
#endif
                  nCCenergies++;

                  /* calibrate mode2 segment energies */

                  assert (ptinp->num < MAX_INTPTS);
                  for (j = 0; j < ptinp->num; j++)
                    {
                      ptinp->intpts[j].e = ptinp->intpts[j].e * Pars.SEGcal_gain[detno][ptinp->intpts[j].seg]
                        + Pars.SEGcal_offset[detno][ptinp->intpts[j].seg];
                    }

                  /* update SEGadd_FPGA and SEGadd_PSA.  */
                  /* as ell as sum versions */


                  for (j = 0; j < ptinp->num; j++)
                    {
                      if (segProcessed[detno][ptinp->intpts[j].seg] == 0)
                        {
                          addedSEGEnergy_FPGA[detno] += ptinp->intpts[j].seg_ener;
                          SEGsum_FPGA->Fill (detno, ptinp->intpts[j].seg_ener);
                          segProcessed[detno][ptinp->intpts[j].seg] = 1;
                        };
                      addedSEGEnergy_PSA[detno] += ptinp->intpts[j].e;
                      SEGsum_PSA->Fill (detno, ptinp->intpts[j].e);
                    };

                  /* quietly rescale all interaction energies to the CC energy */

                  r1 = 0;
                  for (j = 0; j < ptinp->num; j++)
                    r1 += ptinp->intpts[j].e;
                  for (j = 0; j < ptinp->num; j++)
                    ptinp->intpts[j].e = ptinp->intpts[j].e / r1 * ptinp->tot_e;


                  ptinp->tot_e = 0;
                  for (j = 0; j < ptinp->num; j++)
                    ptinp->tot_e += (ptinp->intpts[j].e);

                  /* hit pattern */

                  hitpat->Fill ((double) detno, 1);

                  /* unsuppressed segment energy spectra */

                  for (j = 0; j < ptinp->num; j++)
                    {
                      d1 = detno * 36 + ptinp->intpts[j].seg;
                      d2 = ptinp->intpts[j].seg_ener;
                      if (d1 >= 0 && d1 < MAXSEGNO)
                        if (d1 > 0 && d2 < LONGLEN)
                          SegE->Fill (d1, d2);
                    };

                  /* z_plot */

                  for (j = 0; j < ptinp->num; j++)
                    z_plot->Fill (detno, ptinp->intpts[j].z, 1);

                  /* Now rotate into world coordinates */
                  /* from now on x,y,z will be world coordinates */

                  for (j = 0; j < ptinp->num; j++)
                    {

                      if (Pars.nocrystaltoworldrot == 0)
                        {

                          if (Pars.AGATA_data == 0)
                            {

                              /* rotate into world coordinates first */
                              /* and make it cm rather than mm because */
                              /* crmat needs it in cm */

                              if (Pars.CurEvNo <= Pars.NumToPrint)
                                {
                                  printf ("* %i: ", j);
                                  printf ("%7.2f,%7.2f,%7.2f --> ", ptinp->intpts[j].x, ptinp->intpts[j].y,
                                          ptinp->intpts[j].z);
                                };

                              xx = ptinp->intpts[j].x / 10.0;
                              yy = ptinp->intpts[j].y / 10.0;
                              zz = ptinp->intpts[j].z / 10.0;


                              ptinp->intpts[j].x = Pars.crmat[moduleno][crystalno][0][0] * xx
                                + Pars.crmat[moduleno][crystalno][0][1] * yy
                                + Pars.crmat[moduleno][crystalno][0][2] * zz + Pars.crmat[moduleno][crystalno][0][3];

                              ptinp->intpts[j].y = Pars.crmat[moduleno][crystalno][1][0] * xx
                                + Pars.crmat[moduleno][crystalno][1][1] * yy
                                + Pars.crmat[moduleno][crystalno][1][2] * zz + Pars.crmat[moduleno][crystalno][1][3];

                              ptinp->intpts[j].z = Pars.crmat[moduleno][crystalno][2][0] * xx
                                + Pars.crmat[moduleno][crystalno][2][1] * yy
                                + Pars.crmat[moduleno][crystalno][2][2] * zz + Pars.crmat[moduleno][crystalno][2][3];

                              /* apply y-offsets [are in cm] */

                              ptinp->intpts[j].x += Pars.xyzoffset[moduleno * 4 + crystalno][0];
                              ptinp->intpts[j].y += Pars.xyzoffset[moduleno * 4 + crystalno][1];
                              ptinp->intpts[j].z += Pars.xyzoffset[moduleno * 4 + crystalno][2];

                              if (Pars.CurEvNo <= Pars.NumToPrint)
                                {
                                  printf ("GT: %7.2f,%7.2f,%7.2f\n", ptinp->intpts[j].x, ptinp->intpts[j].y,
                                          ptinp->intpts[j].z);
                                  printf("applied xyzoffsets: %7.4f %7.4f %7.4f\n", Pars.xyzoffset[detno][0], Pars.xyzoffset[detno][1], Pars.xyzoffset[detno][2]);
                               }

                            }
                          else if (Pars.AGATA_data == 1)
                            {
                              detno = moduleno * 3 + crystalno;


                              xx = ptinp->intpts[j].x;
                              yy = ptinp->intpts[j].y;
                              zz = ptinp->intpts[j].z;

                              ptinp->intpts[j].x =
                                Pars.rotxx[detno] * xx + Pars.rotxy[detno] * yy + Pars.rotxz[detno] * zz +
                                Pars.TrX[detno];
                              ptinp->intpts[j].y =
                                Pars.rotyx[detno] * xx + Pars.rotyy[detno] * yy + Pars.rotyz[detno] * zz +
                                Pars.TrY[detno];;
                              ptinp->intpts[j].z =
                                Pars.rotzx[detno] * xx + Pars.rotzy[detno] * yy + Pars.rotzz[detno] * zz +
                                Pars.TrZ[detno];;


                              if (Pars.CurEvNo <= Pars.NumToPrint)
                                {
                                  printf ("AG::x: %9.2f --> %9.2f\n", xx, ptinp->intpts[j].x);
                                  printf ("AG::y: %9.2f --> %9.2f\n", yy, ptinp->intpts[j].y);
                                  printf ("AG::z: %9.2f --> %9.2f\n", zz, ptinp->intpts[j].z);
                                  r1 = xx * xx + yy * yy + zz * zz;
                                  r1 = sqrtf (r1);
                                  r2 = ptinp->intpts[j].x * ptinp->intpts[j].x
                                    + ptinp->intpts[j].y * ptinp->intpts[j].y + ptinp->intpts[j].z * ptinp->intpts[j].z;
                                  r2 = sqrtf (r2);
                                  printf ("AG::radius %f --> %f\n", r1, r2);
                                }

                              ptinp->intpts[j].x /= 10;
                              ptinp->intpts[j].y /= 10;
                              ptinp->intpts[j].z /= 10;

                              /* apply xyz-offsets [are in cm] */

                              ptinp->intpts[j].x += Pars.xyzoffset[detno][0];
                              ptinp->intpts[j].y += Pars.xyzoffset[detno][1];
                              ptinp->intpts[j].z += Pars.xyzoffset[detno][2];

                              if (Pars.CurEvNo <= Pars.NumToPrint)
                                 printf("applied xyzoffsets: %7.4f %7.4f %7.4f\n", Pars.xyzoffset[detno][0], Pars.xyzoffset[detno][1], Pars.xyzoffset[detno][2]);

                            }

                          else
                            {
                              /* no rotation case, just make it cm */

                              xx = ptinp->intpts[j].x / 10.0;
                              yy = ptinp->intpts[j].y / 10.0;
                              zz = ptinp->intpts[j].z / 10.0;

                            };

                        };

                    };




                  /* doppler correction (this is not the only way to do this!) */

                  for (j = 0; j < ptinp->num; j++)
                    {

                      rr =
                        (ptinp->intpts[j].x - Pars.target_x / 10) * (ptinp->intpts[j].x - Pars.target_x / 10) +
                        (ptinp->intpts[j].y - Pars.target_y / 10) * (ptinp->intpts[j].y - Pars.target_y / 10) +
                        (ptinp->intpts[j].z - Pars.target_z / 10) * (ptinp->intpts[j].z - Pars.target_z / 10);
                      rr = sqrtf (rr);

//printf("ptinp->intpts[j].x=%f, ptinp->intpts[j].x - Pars.target_x=%f\n",ptinp->intpts[j].x,ptinp->intpts[j].x - Pars.target_x);
//printf("ptinp->intpts[j].y=%f, ptinp->intpts[j].y - Pars.target_y=%f\n",ptinp->intpts[j].y,ptinp->intpts[j].y - Pars.target_y);
//printf("ptinp->intpts[j].z=%f, ptinp->intpts[j].z - Pars.target_z=%f\n",ptinp->intpts[j].z,ptinp->intpts[j].z - Pars.target_z);


                      if (rr > RMIN && rr < RMAX)
                        {
                          radius_all->Fill ((double) rr, 1);
                          if (ptinp->intpts[j].e > 0 && ptinp->intpts[j].e < MEDIUMLEN);
                          evsr_all->Fill ((double) rr, ptinp->intpts[j].e);
                        };

                      dp = ((ptinp->intpts[j].x - Pars.target_x / 10) * Pars.beamdir[0] +
                            (ptinp->intpts[j].y - Pars.target_y / 10) * Pars.beamdir[1] +
                            (ptinp->intpts[j].z - Pars.target_z / 10) * Pars.beamdir[2]) / rr;

                      if (dp < -1.0)
                        dp = -1.0;
                      if (dp > 1.0)
                        dp = 1.0;
                      polang[j] = acosf (dp);

                      rr = 1.0 - Pars.beta * Pars.beta;
                      doppler_factor[j] = sqrt (rr) / (1.0 - Pars.beta * cos (polang[j]));

                      /* central contact energy matrix and total energy */
                      /* the 1.0/ptinp->num preserves the counts */

                      if (detno >= 0 && detno < MAXDETPOS)
                        if (ptinp->tot_e > 0 && ptinp->tot_e < (LONGLEN + 30))
                          {
                            CCsum_nodop->Fill (ptinp->tot_e, 1.0 / ptinp->num);
                            CCsum->Fill ((double) ptinp->tot_e / doppler_factor[j], 1.0 / ptinp->num);
                            CCe->Fill ((double) detno, (double) ptinp->tot_e / doppler_factor[j], 1.0 / ptinp->num);
                          };

                    };

                  /* worldmap all hits */
                  /* target position does not get into this */

                  for (j = 0; j < ptinp->num; j++)
                    {

                      polAng = findPolarFromCartesian (ptinp->intpts[j].x, ptinp->intpts[j].y, ptinp->intpts[j].z, &rr);
                      aziAng = findAzimuthFromCartesian (ptinp->intpts[j].x, ptinp->intpts[j].y, ptinp->intpts[j].z);

                      ndethits[detno]++;
                      pol[detno] += polAng;
                      azi[detno] += aziAng;


                      /* SMAP coordinates */

                      SMAP_allhits_sq->Fill ((double) aziAng / RAD2DEG, (double) polAng / RAD2DEG, 1);

                      sX = aziAng * sinf (polAng) / RAD2DEG;
                      sY = polAng / RAD2DEG;    /* + 1.5; */

                      if (Pars.CurEvNo <= Pars.NumToPrint && 0)
                        {
                          printf ("%i [type %i] ", j, GEB_event->ptgd[i]->type);
                          printf ("e: %9.2f/%9.2f ", ptinp->intpts[j].e, ptinp->tot_e);
                          printf ("(%6.2f,%6.2f,%6.2f)cry --> ", xx, yy, zz);
                          printf ("(%6.2f,%6.2f,%6.2f)world(cm); ", ptinp->intpts[j].x, ptinp->intpts[j].y,
                                  ptinp->intpts[j].z);
//                    printf (" sX,sY=%6.2f,%6.2f ", sX, sY);
                          printf ("\n");
                        };

                      /* update */

                      if (sX >= -180 && sX <= 180 && sY >= 0 && sY <= 180)
                        SMAP_allhits->Fill (sX, sY, 1);
                      else
                        {
                          if (nperrors < 10)
                            {
                              nperrors++;
                              printf ("error: sX,sY= ( %11.6f , %11.6f )\n", sX, sY);
//                          exit (1);
                            };
                        };


                    };          /* for (j = 0; j < ptinp->num; j++) */


                };              /* if (Pars.enabled[detno]) */

            };                  /* pad==0 */

        };                      /* if (GEB_event->ptgd[i]->type == GEB_TYPE_DECOMP) */

    };                          /* for (i = 0; i < GEB_event->mult; i++) */

// if (Pars.CurEvNo >2000) exit(0);

  /* crystal crystal hit matrix, who sees who */

  if (ncr >= 2)
    for (i = 0; i < ncr; i++)
      for (j = i + 1; j < ncr; j++)
        crXcr->Fill (crhit[i], crhit[j], 1);

  /* update added energy spectrum */
  addedEnergy = 0;
  for (i = 0; i < nCCenergies; i++)
    addedEnergy += CCenergies[i];

  CCadd->Fill ((double) addedEnergy, 1);
  if (Pars.CurEvNo <= Pars.NumToPrint)
    printf ("updated CCadd with %7.2f\n", addedEnergy);
  for (i = 0; i < MAXDETNO; i++)
    {
      SEGadd_FPGA->Fill ((double) i, (double) addedSEGEnergy_FPGA[i], 1);
      SEGadd_PSA->Fill ((double) i, (double) addedSEGEnergy_PSA[i], 1);
    }


  /* update mean energy */

  for (i = 0; i < GEB_event->mult; i++)
    {
      if (GEB_event->ptgd[i]->type == GEB_TYPE_DECOMP)
        {

          /* cast */

          ptinp = (CRYS_INTPTS *) GEB_event->ptinp[i];
          if (ptinp->pad == 0)
            {
              meanCCenergy_sum += ptinp->tot_e;
              nCCenergy_sum++;
            };
        };
    };
  meanCCenergy_add += addedEnergy;
  nCCenergy_add++;

  /* update the special spectra used to */
  /* find the scattering probability */

  switch (ncr)
    {
    case 2:
      hit2->Fill ((double) addedEnergy, 1);
      break;
    case 3:
      hit3->Fill ((double) addedEnergy, 1);
      break;
    case 4:
      hit4->Fill ((double) addedEnergy, 1);
      break;
    }

  /* fill the ggCC martrix */


  CCmult->Fill (nCCenergies, 1);
  if (nCCenergies >= Pars.multlo && nCCenergies <= Pars.multhi)
    for (i = 0; i < nCCenergies; i++)
      for (j = i + 1; j < nCCenergies; j++)
        {
          ggCC->Fill (CCenergies[i], CCenergies[j], 1.0);
          ggCC->Fill (CCenergies[j], CCenergies[i], 1.0);
        };

  /* numbers of hits in the array */

  if (Pars.CurEvNo <= Pars.NumToPrint)
    printf ("total number of hits in array: %i\n", numHitArray2);
  numHitArray2_sp->Fill (numHitArray2, 1);

  /* done */

  if (Pars.CurEvNo <= Pars.NumToPrint)
    printf ("exit bin_mode2\n");

  return (0);

}
