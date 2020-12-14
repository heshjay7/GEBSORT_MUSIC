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

#include "veto_pos.h"
#include "GEBSort.h"
#include "GTMerge.h"

#define MAXNGAM 20

#define HKLEN 127
#define MAXH HKLEN
#define MAXK HKLEN

#define MINDOPFAC 0.5
#define MAXDOPFAC 1.5

/* pointers to ROOT spectra */


TH2F *SMAP_firsthits;
TH1D *fmsp;
TH2F *ecos_raw;
TH2F *ecos_dopcor;
TH1D *dopfac;
TH1D *polangle;
TH1D *numHitArray1_sp;
TH1D *gmult;
TH1D *sumTrackE;
TH1D *TrackE_1gate;
TH1D *TrackE_2gates;
TH1D *TrackE_3gates;
TH2F *fomXe;
TH2F *ndetXfom;
TH2F *HK;
TH2F *HK1;
TH2F *HK2;
TH2F *HK3;
TH1D *rate_mode1;
TH2F *gg;
TH2F *ndet_e;
TH2F *rad_e;
TH2F *rad_e_raw;
TH1D *radius_first;
TH2F *evsr_first;
TH2F *ngamXsumTrackE;
TH2F *cr_trackedE;
TH2F *nCC_e;
TH2F *firstScatAng;
TH2F *firstScatRange;

typedef struct PAYLOAD
{
  char p[MAXDATASIZE];
} PAYLOAD;

typedef struct TRACK_STRUCT
{
  int n;
  GEBDATA *gd;
  PAYLOAD *payload;
} TRACK_STRUCT;

/* parameters */

extern PARS Pars;
extern EXCHANGE exchange;

/* gate lookup spectrum */

int gate_spe[4000];

/*-----------------------------------------------------*/

int
gggate (int mult_g, float e_good[40], int *n1gated, int *n2gated, int *n3gated)
{

  /* declarations */

  int ifocus, i, ngates;

  *n1gated = 0;
  *n2gated = 0;
  *n3gated = 0;

  /* trivial return? */

  if (mult_g == 0 || mult_g >= 40)
    return (0);

  /* loop over gamma rays */

  for (ifocus = 0; ifocus < mult_g; ifocus++)
    {

      /* how many gates has this gamma ray */

      ngates = 0;
      for (i = 0; i < mult_g; i++)
        if (i != ifocus)
          {
            ngates += gate_spe[(int) e_good[i]];
          };

      /* exchange */

      exchange.ngates = ngates;

      /* update */

      if (ngates >= 1)
        {
          (*n1gated)++;
          TrackE_1gate->Fill ((double) e_good[ifocus], 1);
        };

      if (ngates >= 2)
        {
          (*n2gated)++;
          TrackE_2gates->Fill ((double) e_good[ifocus], 1);
        };

      if (ngates >= 3)
        {
          (*n3gated)++;
          TrackE_3gates->Fill ((double) e_good[ifocus], 1);
        };


    }

  /* done */

  return (0);
};


/*-----------------------------------------------------*/

int
print_tracked_gamma_rays (FILE * fp, TRACKED_GAMMA_HIT * grh)
{

/* declarations */

  int j;
  float r1;

  fprintf (fp, "number of gamma rays: %i\n", grh->ngam);
  for (j = 0; j < grh->ngam; j++)
    {
      fprintf (fp, "esum=%6.3f, ", grh->gr[j].esum);
      fprintf (fp, "ndet (#interactions)=%2i, ", grh->gr[j].ndet);
      fprintf (fp, "fom=%6.3f, ", grh->gr[j].fom);
      fprintf (fp, "tracked=%i\n", grh->gr[j].tracked);
      fprintf (fp, "1'th hit: (%6.3f, %6.3f, %6.3f), e= %6.1f keV\n", grh->gr[j].x0, grh->gr[j].y0, grh->gr[j].z0,
               grh->gr[j].e0);
      fprintf (fp, "2'nd hit: (%6.3f, %6.3f, %6.3f), e= %6.1f keV\n", grh->gr[j].x1, grh->gr[j].y1, grh->gr[j].z1,
               grh->gr[j].e1);
      r1 = grh->gr[j].x0 * grh->gr[j].x0 + grh->gr[j].y0 * grh->gr[j].y0 + grh->gr[j].z0 * grh->gr[j].z0;
      r1 = sqrtf (r1);
      fprintf (fp, "1'th interaction radius: %9.1f mm\n", r1);
    };

  /* done */

  return (0);

};

/*-----------------------------------------------------*/

int
sup_mode1 ()
{
  /* declarations */

  char str1[STRLEN], str2[STRLEN];
  int jj, ii;
  unsigned int seed;
  double nn, lo, hi;
  int ener_l[40], ener_h[40];

  TH1D *mkTH1D (char *, char *, int, double, double);
  TH2F *mkTH2F (char *, char *, int, double, double, int, double, double);
  int get_a_seed (unsigned int *);

  /* define spectra */

  sprintf (str1, "ecos_raw");
  sprintf (str2, "ecos_raw");
  ecos_raw = mkTH2F (str1, str2, 2048, 0, 2047, 1000, -1, 1);
  sprintf (str1, "uncorrected Energy");
  ecos_raw->SetXTitle (str1);
  sprintf (str1, "cos()");
  ecos_raw->SetYTitle (str1);


  sprintf (str1, "ecos_dopcor");
  sprintf (str2, "ecos_dopcor");
  ecos_dopcor = mkTH2F (str1, str2, 2048, 0, 2047, 1000, -1, 1);
  sprintf (str1, "doppler corrected Energy");
  ecos_dopcor->SetXTitle (str1);
  sprintf (str1, "cos()");
  ecos_dopcor->SetYTitle (str1);


  sprintf (str1, "SMAP_firsthits");
  sprintf (str2, "SMAP_firsthits");
  SMAP_firsthits = mkTH2F (str1, str2, 256, -180, 180, 256, 0, 180);
  sprintf (str1, "horizontal");
  SMAP_firsthits->SetXTitle (str1);
  sprintf (str1, "vertical");
  SMAP_firsthits->SetYTitle (str1);


  sprintf (str1, "fm");
  sprintf (str2, "figure of merit");
  fmsp = mkTH1D (str1, str2, 1024, 0, 2.1);
  fmsp->SetXTitle (str1);


  sprintf (str1, "polangle");
  sprintf (str2, "polangle");
  polangle = mkTH1D (str1, str2, 1024, 0, 180);
  polangle->SetXTitle (str1);


  sprintf (str1, "dopfac");
  sprintf (str2, "dopfac");
  dopfac = mkTH1D (str1, str2, 2048, MINDOPFAC, MAXDOPFAC);
  dopfac->SetXTitle (str1);


  sprintf (str1, "gmult");
  sprintf (str2, "tracked gamma ray mutiplicity");
  gmult = mkTH1D (str1, str2, 21, 0, 20);
  gmult->SetXTitle (str1);


  sprintf (str1, "sumTrackE");
  sumTrackE = mkTH1D (str1, str1, LONGLEN, 1, LONGLEN);
  sumTrackE->SetXTitle (str1);


  sprintf (str1, "TrackE_1gate");
  TrackE_1gate = mkTH1D (str1, str1, LONGLEN, 1, LONGLEN);
  TrackE_1gate->SetXTitle (str1);


  sprintf (str1, "TrackE_2gates");
  TrackE_2gates = mkTH1D (str1, str1, LONGLEN, 1, LONGLEN);
  TrackE_2gates->SetXTitle (str1);

  sprintf (str1, "TrackE_3gates");
  TrackE_3gates = mkTH1D (str1, str1, LONGLEN, 1, LONGLEN);
  TrackE_3gates->SetXTitle (str1);


  sprintf (str1, "fomXe");
  sprintf (str2, "fomXe");
  fomXe = mkTH2F (str1, str2, LONGLEN, 1, LONGLEN, 211, 0, 2.1);
  fomXe->SetXTitle ("e");
  fomXe->SetYTitle ("fom");


  sprintf (str1, "cr_trackedE");
  sprintf (str2, "cr_trackedE");
  cr_trackedE = mkTH2F (str1, str2, MAXDETNO + 1, 0, MAXDETNO, Pars.GGMAX, 1, Pars.GGMAX);
  cr_trackedE->SetXTitle ("crystal no");
  cr_trackedE->SetYTitle ("tracked spectrum");

  sprintf (str1, "firstScatAng");
  sprintf (str2, "firstScatAng");
  firstScatAng = mkTH2F (str1, str2, Pars.GGMAX, 1, Pars.GGMAX, 181, 0, 180);
  firstScatAng->SetXTitle ("energy [keV]");
  firstScatAng->SetYTitle ("scatter angle [deg]");

  sprintf (str1, "firstScatRange");
  sprintf (str2, "firstScatRange");
  firstScatRange = mkTH2F (str1, str2, Pars.GGMAX, 1, Pars.GGMAX, 200, 1, 200);
  firstScatRange->SetXTitle ("energy [keV]");
  firstScatRange->SetYTitle ("scatter range [mm]");

  sprintf (str1, "ndetXfom");
  ndetXfom = mkTH2F (str1, str1, 8, 1, 8, 211, 0, 2.1);
  ndetXfom->SetXTitle ("ndet (# interaction points)");
  ndetXfom->SetYTitle ("fom");


  sprintf (str1, "ngamXsumTrackE");
  sprintf (str2, "# gammarays vs gamma ray energy");
  ngamXsumTrackE = mkTH2F (str1, str2, MAXNGAM, 1, MAXNGAM, LONGLEN, 1, LONGLEN);
  sprintf (str1, "ngam");
  ngamXsumTrackE->SetXTitle (str1);
  sprintf (str1, "sumTrackE");
  ngamXsumTrackE->SetYTitle (str1);

  /* these values make the mean work */
  /* properly for the HK matrices */

  nn = 127;
  lo = -0.5;
  hi = nn + lo;

  sprintf (str1, "HK");
  sprintf (str2, "ungated tracked HK with FOM cut");
  HK = mkTH2F (str1, str2, nn, lo, hi, nn, lo, hi);
  sprintf (str1, "K, # gamma rays +(rn-0.5)");
  HK->SetXTitle (str1);
  sprintf (str1, "H, summed energy");
  HK->SetYTitle (str1);

  sprintf (str1, "HK1");
  sprintf (str2, "1 gated tracked HK with FOM cut");
  HK1 = mkTH2F (str1, str2, nn, lo, hi, nn, lo, hi);
  sprintf (str1, "K, # gamma rays +(rn-0.5)");
  HK1->SetXTitle (str1);
  sprintf (str1, "H, summed energy");
  HK1->SetYTitle (str1);

  sprintf (str1, "HK2");
  sprintf (str2, "2 gated tracked HK with FOM cut");
  HK2 = mkTH2F (str1, str2, nn, lo, hi, nn, lo, hi);
  sprintf (str1, "K, # gamma rays +(rn-0.5)");
  HK2->SetXTitle (str1);
  sprintf (str1, "H, summed energy");
  HK2->SetYTitle (str1);

  sprintf (str1, "HK3");
  sprintf (str2, "3 gated tracked HK with FOM cut");
  HK3 = mkTH2F (str1, str2, nn, lo, hi, nn, lo, hi);
  sprintf (str1, "K, # gamma rays +(rn-0.5)");
  HK3->SetXTitle (str1);
  sprintf (str1, "H, summed energy");
  HK3->SetYTitle (str1);






  sprintf (str1, "rate_mode1");
  sprintf (str2, "rate_mode1 gamma rays");
  rate_mode1 = mkTH1D (str1, str2, RATELEN, 0, RATELEN);
  rate_mode1->SetXTitle ("minutes");
  rate_mode1->SetYTitle ("Hz");


  sprintf (str1, "gg");
  sprintf (str2, "tracked gg matrix");
  gg = mkTH2F (str1, str2, Pars.GGMAX, 1, Pars.GGMAX, Pars.GGMAX, 1, Pars.GGMAX);
  sprintf (str1, "g1");
  gg->SetXTitle (str1);
  sprintf (str1, "g2");
  gg->SetYTitle (str1);


  sprintf (str1, "numHitArray1");
  sprintf (str2, "numHitArray1");
  numHitArray1_sp = mkTH1D (str1, str2, 201, 0, 200);
  sprintf (str1, "# array interactions");
  numHitArray1_sp->SetXTitle (str1);


  sprintf (str1, "ndet_e");
  sprintf (str2, "interaction points vs gamma energy");
  ndet_e = mkTH2F (str1, str2, 10, 1, 10, Pars.GGMAX, 1, (double) Pars.GGMAX);
  sprintf (str1, "ndet");
  ndet_e->SetXTitle (str1);
  sprintf (str1, "tracked e");
  ndet_e->SetYTitle (str1);

  sprintf (str1, "nCC_e");
  sprintf (str2, "tracked gamma rays vs number of CC");
  nCC_e = mkTH2F (str1, str2, 20, 1, 20, Pars.GGMAX, 1, (double) Pars.GGMAX);
  sprintf (str1, "nCC");
  nCC_e->SetXTitle (str1);
  sprintf (str1, "tracked e");
  nCC_e->SetYTitle (str1);

  sprintf (str1, "rad_e");
  sprintf (str2, "rad_e with FOM cut");
  rad_e = mkTH2F (str1, str2, Pars.GGMAX, 1, Pars.GGMAX, 2048, 0, 340);
  sprintf (str1, "Energy (keV)");
  rad_e->SetXTitle (str1);
  sprintf (str1, "Radius (mm)");
  rad_e->SetYTitle (str1);

  sprintf (str1, "rad_e_raw");
  sprintf (str2, "rad_e_raw with FOM cut");
  rad_e_raw = mkTH2F (str1, str2, Pars.GGMAX, 1, Pars.GGMAX, 2048, 0, 340);
  sprintf (str1, "Energy (keV)");
  rad_e_raw->SetXTitle (str1);
  sprintf (str1, "Radius (mm)");
  rad_e_raw->SetYTitle (str1);

  sprintf (str1, "radius_first");
  sprintf (str2, "radius (first points)");
  radius_first = mkTH1D (str1, str2, 4096, RMIN, RMAX);
  radius_first->SetXTitle (str1);

  sprintf (str1, "evsr_first");
  sprintf (str2, "evsr_first");
  evsr_first = mkTH2F (str1, str2, MEDIUMLEN, RMIN, RMAX, MEDIUMLEN, 1, MEDIUMLEN);
  sprintf (str1, "radius (cm)");
  evsr_first->SetXTitle (str1);
  sprintf (str1, "energy");
  evsr_first->SetYTitle (str1);


  /* list what we have */

//  printf (" we have define the following spectra:\n");

//  Pars.wlist = gDirectory->GetList ();
//  Pars.wlist->Print ();

  get_a_seed (&seed);
  srand (seed);

  for (jj = 0; jj < 4000; jj++)
    {
      gate_spe[jj] = 0;
    };

  for (ii = 0; ii < Pars.numgggates; ii++)
    {
      ener_l[ii] = Pars.gg_gate_pos[ii] - Pars.gg_gate_width[ii];
      ener_h[ii] = Pars.gg_gate_pos[ii] + Pars.gg_gate_width[ii];

      for (jj = ener_l[ii]; jj < ener_h[ii]; jj++)
        {

          gate_spe[jj] = 1;
        };

    };



  return (0);

};

/* ----------------------------------------------------------------- */

int
exit_bin_mode1 ()
{

  int dim, i;
  float rr[LONGLEN];
  char str[512];
//  int wr_spe (char *, int *, float *);

  /* write the tracked spectrum in spe format as well */

  dim = LONGLEN;
  for (i = 0; i < dim; i++)
    {
      rr[i] = (float) sumTrackE->GetBinContent (i);
    };
  sprintf (str, "sumTrackE.spe");
//  wr_spe (str, &dim, rr);

  /* write the gated spectra out */

  dim = LONGLEN;
  for (i = 0; i < dim; i++)
    {
      rr[i] = (float) TrackE_1gate->GetBinContent (i);
    };
  sprintf (str, "TrackE_1gate.spe");
//  wr_spe (str, &dim, rr);

  dim = LONGLEN;
  for (i = 0; i < dim; i++)
    {
      rr[i] = (float) TrackE_2gates->GetBinContent (i);
    };
  sprintf (str, "TrackE_2gates.spe");
//  wr_spe (str, &dim, rr);

  dim = LONGLEN;
  for (i = 0; i < dim; i++)
    {
      rr[i] = (float) TrackE_3gates->GetBinContent (i);
    };
  sprintf (str, "TrackE_3gates.spe");
//  wr_spe (str, &dim, rr);

  return (0);

};

/* ----------------------------------------------------------------- */

int
bin_mode1 (GEB_EVENT * GEB_event)
{

  /* declarations */

  TRACKED_GAMMA_HIT *grh;
  int k, l, i, j, nMode1 = 0;
  float sX, sY, polAng, aziAng, rr, dp;
  double d1, d2;
  float RAD2DEG = 0.01745329;
  char str[128];
  static long long int t0;
  static int firsttime = 1;
  float polang[MAX_GAMMA_RAYS];
  float doppler_factor[MAX_GAMMA_RAYS];
  int numHitArray1 = 0;
  int nTrackedGammas = 0;
  int n1gated, n2gated, n3gated;
  int KK, nCC;
  float HH;
  float xx, yy, zz, x1, y1, z1, dd, x2, y2, z2;

/* define the gates */

  int numgate, mult_g, num_sp1;
  float e_good[40];
  num_sp1 = 1;
  numgate = Pars.numgggates;


  /* prototypes */

  int GebTypeStr (int type, char str[]);
  float findAzimuthFromCartesian (float, float, float);
  float findPolarFromCartesian (float, float, float, float *);



  if (Pars.CurEvNo <= Pars.NumToPrint)
    printf ("--------------------\nentered bin_mode1:\n");


#if(0)
  /* find the number of tracked headers in event */


  nMode1 = 0;
  for (i = 0; i < GEB_event->mult; i++)
    {
      if (GEB_event->ptgd[i]->type == GEB_TYPE_TRACK)
        nMode1++;
      if (GEB_event->ptgd[i]->type == GEB_TYPE_DECOMP)
        nMode2++;
    };
  printf ("nMode2=%i, nMode1=%i\n", nMode2, nMode1);
  if (nMode1 != 1)
    printf ("nMode1=%i\n", nMode1);
  assert (nMode1 <= 1);
#endif

  /* process the total number of interactions in the array */
  /* find the number of tracked headers in event */
  /* and gate on them */

  nTrackedGammas = 0;
  numHitArray1 = 0;
  for (i = 0; i < GEB_event->mult; i++)
    if (GEB_event->ptgd[i]->type == GEB_TYPE_TRACK)
      {
        grh = (TRACKED_GAMMA_HIT *) GEB_event->ptinp[i];
        nTrackedGammas += grh->ngam;
        for (j = 0; j < grh->ngam; j++)
          numHitArray1 += grh->gr[j].ndet;
      };

  /* num interactions and gamma ray gating */

  if (numHitArray1 < Pars.minnumHitArray || numHitArray1 > Pars.maxnumHitArray)
    return (0);

  if (nTrackedGammas < Pars.minNumGammas || nTrackedGammas > Pars.maxNumGammas)
    return (0);


  if (numHitArray1 > 0 && numHitArray1 < 200)
    numHitArray1_sp->Fill (numHitArray1, 1);


  /* find number of central contact in the event */
  /* this reflects the trigger of events */

  nCC = 0;
  for (i = 0; i < GEB_event->mult; i++)
    if (GEB_event->ptgd[i]->type == GEB_TYPE_DECOMP)
      nCC++;

  /* process mode 1 data */

  nMode1 = 0;
  for (i = 0; i < GEB_event->mult; i++)
    {

      if (GEB_event->ptgd[i]->type == GEB_TYPE_TRACK)
        {

          if (Pars.CurEvNo <= Pars.NumToPrint)
            {
              GebTypeStr (GEB_event->ptgd[i]->type, str);
              printf ("bin_mode1, %2i> %2i, %s, TS=%lli;\n", i, GEB_event->ptgd[i]->type, str,
                      GEB_event->ptgd[i]->timestamp);
            }

          /* mode 1 rate spectrum, x=minute, y=Hz */

          if (firsttime)
            {
              firsttime = 0;
              t0 = GEB_event->ptgd[i]->timestamp;
            };
          d1 = (double) (GEB_event->ptgd[i]->timestamp - t0);
          d1 /= 100000000;
          d1 /= 60;
          if (d1 > 0 && d1 < (double) RATELEN)
            rate_mode1->Fill (d1, 1 / 60.0);
          nMode1++;

          grh = (TRACKED_GAMMA_HIT *) GEB_event->ptinp[i];

          gmult->Fill (grh->ngam, 1);

          /* truncate fom in overflow channel */

          for (j = 0; j < grh->ngam; j++)
            if (grh->gr[j].tracked)
              if (grh->gr[j].fom > 2.0)
                grh->gr[j].fom = 2.0;

          /* check for multiplicity requirements */

          if (grh->ngam < Pars.multlo || grh->ngam > Pars.multhi)
            return (0);

          if (Pars.CurEvNo <= Pars.NumToPrint)
            {
              printf ("ngam=%i\n", grh->ngam);
              for (j = 0; j < grh->ngam; j++)
                {
                  printf ("%2i: esum=%8.0f, fom=%4.2f, tracked %i, gate %i\n", j, grh->gr[j].esum, grh->gr[j].fom,
                          grh->gr[j].tracked, gate_spe[(int) grh->gr[j].esum]);
                };
            };

          /* modify FOM for single inteaction */

          for (j = 0; j < grh->ngam; j++)
            if (grh->gr[j].tracked)
              {
                if (grh->gr[j].ndet == 1)
                  if (grh->gr[j].esum > Pars.maxsnglintrE)
                    grh->gr[j].fom = Pars.maxsnglintrEFOM;
              };

          /* find doppler correction factors */

          if (Pars.beta >= 0.000)
            {
              for (j = 0; j < grh->ngam; j++)
                if (grh->gr[j].tracked)
                  {

                    /* note: x,y,z and rr are in mm */
                    /* so are Pars.target_x|y|z */

                    /* find normalized vector from target position */

                    rr = (grh->gr[j].x0 - Pars.target_x) * (grh->gr[j].x0 - Pars.target_x)
                      + (grh->gr[j].y0 - Pars.target_y) * (grh->gr[j].y0 - Pars.target_y)
                      + (grh->gr[j].z0 - Pars.target_z) * (grh->gr[j].z0 - Pars.target_z);
                    rr = sqrtf (rr);

                    /* bin rr but in cm */

                    if (j == 0)
                      if (rr / 10 > RMIN && rr / 10 < RMAX)
                        {
                          radius_first->Fill ((double) rr / 10, 1);
                          if (grh->gr[j].esum > 0 && grh->gr[j].esum < MEDIUMLEN)
                            evsr_first->Fill ((double) rr / 10, grh->gr[j].esum);
                        };

                    /* dot with beam direction */

                    dp = ((grh->gr[j].x0 - Pars.target_x) * Pars.beamdir[0] +
                          (grh->gr[j].y0 - Pars.target_y) * Pars.beamdir[1] +
                          (grh->gr[j].z0 - Pars.target_z) * Pars.beamdir[2]) / rr;

                    /* find polar angle and bin it */

                    if (dp < -1.0)
                      dp = -1.0;
                    if (dp > 1.0)
                      dp = 1.0;
                    polang[j] = acosf (dp);
                    polangle->Fill (polang[j] / M_PI * 180, 1);

                    d1 = (double) (grh->gr[j].esum);
                    if (d1 < 0)
                      d1 = 0;
                    if (d1 > Pars.GGMAX)
                      d1 = Pars.GGMAX - 1;

                    rad_e_raw->Fill (d1, rr, 1.0);

                    if (grh->gr[j].ndet >= Pars.ndetlimlo && grh->gr[j].ndet <= Pars.ndetlimhi)
                      if (grh->gr[j].fom >= Pars.fomlo[grh->gr[j].ndet]
                          && grh->gr[j].fom <= Pars.fomhi[grh->gr[j].ndet])
                        {
                          d1 = (double) (grh->gr[j].esum);
                          if (d1 < 0)
                            d1 = 0;
                          if (d1 > Pars.GGMAX)
                            d1 = Pars.GGMAX - 1;

                          rad_e->Fill (d1, rr, 1.0);
                        };

                    /* find a bin doppler factor */

                    rr = 1.0 - Pars.beta * Pars.beta;
                    doppler_factor[j] = sqrt (rr) / (1.0 - Pars.beta * cos (polang[j]));
                    if (doppler_factor[j] > MINDOPFAC && doppler_factor[j] < MAXDOPFAC)
                      dopfac->Fill (doppler_factor[j], 1);

                    if (Pars.CurEvNo <= Pars.NumToPrint)
                      {
                        printf ("doppler_factor[%i]=%f\n", j, doppler_factor[j]);
                      }
                  };

              /* apply doppler corrections */

              for (j = 0; j < grh->ngam; j++)
                if (grh->gr[j].tracked)
                  {

                    ecos_raw->Fill (grh->gr[j].esum, cos (polang[j]));

                    if (Pars.CurEvNo <= Pars.NumToPrint)
                      printf ("dopler correction: esum %f e0 %f e1 %f -> ", grh->gr[j].esum, grh->gr[j].e0,
                              grh->gr[j].e1);

                    grh->gr[j].esum /= doppler_factor[j];
                    ecos_dopcor->Fill (grh->gr[j].esum, cos (polang[j]));
                    grh->gr[j].e0 /= doppler_factor[j];
                    grh->gr[j].e1 /= doppler_factor[j];
                    if (Pars.CurEvNo <= Pars.NumToPrint)
                      printf ("%f %f %f\n", grh->gr[j].esum, grh->gr[j].e0, grh->gr[j].e1);

                  };
            };

          if (Pars.CurEvNo <= Pars.NumToPrint)
            print_tracked_gamma_rays (stdout, grh);

          /* tracked spectra vs CC's that fired */

          for (j = 0; j < grh->ngam; j++)
            if (grh->gr[j].tracked)
              if (grh->gr[j].ndet >= Pars.ndetlimlo && grh->gr[j].ndet <= Pars.ndetlimhi)
                if (grh->gr[j].esum < Pars.GGMAX)
                  if (nCC < 20)
                    nCC_e->Fill ((double) nCC, (double) grh->gr[j].esum, 1);

          /* energy spectra vs ndet */

          for (j = 0; j < grh->ngam; j++)
            if (grh->gr[j].tracked)
              if (grh->gr[j].ndet >= Pars.ndetlimlo && grh->gr[j].ndet <= Pars.ndetlimhi)
                {
                  if (grh->gr[j].ndet < 10)
                    if (grh->gr[j].esum < Pars.GGMAX)
                      ndet_e->Fill ((double) grh->gr[j].ndet, (double) grh->gr[j].esum, 1);
                };


          for (j = 0; j < grh->ngam; j++)
            if (grh->gr[j].tracked)
              {
                polAng = findPolarFromCartesian (grh->gr[j].x0, grh->gr[j].y0, grh->gr[j].z0, &rr);
                aziAng = findAzimuthFromCartesian (grh->gr[j].x0, grh->gr[j].y0, grh->gr[j].z0);


                /* SMAP coordinates */

                sX = aziAng * sinf (polAng) / RAD2DEG;
                sY = polAng / RAD2DEG;  /* + 1.5; */

                /* update */

                if (Pars.CurEvNo <= Pars.NumToPrint)
                  {
                    printf ("aziAng=%6.2f\n", aziAng / M_PI * 180);
                    printf ("polAng=%6.2f\n", polAng / M_PI * 180);
                    printf (" sX,sY=%6.2f,%6.2f\n", sX, sY);
                    fflush (stdout);
                  };

                if (sX > -180 && sX < 180)
                  if (sY > 0 && sY < 180)
                    SMAP_firsthits->Fill (sX, sY, 1);
              };


          /* esum vs fom matrix */

          for (j = 0; j < grh->ngam; j++)
            if (grh->gr[j].tracked)
              {
                d1 = (double) (grh->gr[j].esum);
                if (grh->gr[j].fom < 2.1 && grh->gr[j].fom >= 0)
                  if (d1 > 0 && d1 < LONGLEN)
                    if (grh->gr[j].ndet >= Pars.ndetlimlo && grh->gr[j].ndet <= Pars.ndetlimhi)
                      {

                        fomXe->Fill (d1, grh->gr[j].fom, 1);
//                        printf("*** %f %f\n", (float)d1, (float)grh->gr[j].fom);

                        if (grh->gr[j].fom >= Pars.fomlo[grh->gr[j].ndet])
                          if (grh->gr[j].fom <= Pars.fomhi[grh->gr[j].ndet])
                            {
                              sumTrackE->Fill (d1, 1);
                              if (grh->ngam < MAXNGAM)
                                ngamXsumTrackE->Fill (grh->ngam, d1, 1);
#if (TRACK2==1)
                              if (grh->gr[j].fhcrID < MAXDETNO && d1 < Pars.GGMAX)
                                cr_trackedE->Fill ((double) grh->gr[j].fhcrID, d1, 1);
#endif

                            };
                      };
              };
//  fomXe = mkTH2F (str1, str2,MAXDETNO+1 , 0, MAXDETNO, Pars.GGMAX, 1, Pars.GGMAX);

          /* fom spectrum */

          for (j = 0; j < grh->ngam; j++)
            if (grh->gr[j].tracked)
              if (grh->gr[j].fom >= 0 && grh->gr[j].fom < 2.1)
                {
                  fmsp->Fill (grh->gr[j].fom, 1);
                  ndetXfom->Fill (grh->gr[j].ndet, grh->gr[j].fom);
                }

          /* make list of FOM++ gated gamma rays */
          /* find H and K at the same time */

          KK = 0;
          HH = 0;
          mult_g = 0;
          for (k = 0; k < grh->ngam; k++)
            if (grh->gr[k].tracked)
              if (grh->gr[k].ndet >= Pars.ndetlimlo && grh->gr[k].ndet <= Pars.ndetlimhi)
                if (grh->gr[k].fom >= Pars.fomlo[grh->gr[k].ndet] && grh->gr[k].fom <= Pars.fomhi[grh->gr[k].ndet])
                  if (grh->gr[k].esum > 0 && grh->gr[k].esum <= 4000)
                    {
                      e_good[mult_g] = grh->gr[k].esum;
                      mult_g++;
                      KK++;
                      HH += grh->gr[k].esum;
                    };
          HH /= Pars.Hresolution;

          /* process the gated spectra */

          gggate (mult_g, e_good, &n1gated, &n2gated, &n3gated);
          if (Pars.CurEvNo <= Pars.NumToPrint)
            {
              printf ("n1gated=%i ,n2gated=%i, n3gated=%i\n", n1gated, n2gated, n3gated);
            };

          /* HK matrices, d2: summed energy, d1: number of gamma rays */
          /* HK is ungated, HK1 and HK2 are one and two gated */
          /* we only sum up tracked and FOM OK events for the HK!! */
          /* note: d1 is not necessarely grh->ngam because of FOM cut */

          if (KK >= 0 && KK < 100)
            if (HH >= 0 && HH < 100)
              {

                if (Pars.CurEvNo <= Pars.NumToPrint)
                  {
                    printf ("HK: KK= %i HH= %f\n", KK, HH);
                    fflush (stdout);
                  }

                /* fill ungated HK matrix */

                HK->Fill ((double) KK, (double) HH, 1);

                /* fill gated HK matrices */

                if (n1gated > 0)
                  {
                    HK1->Fill ((double) KK, (double) HH, 1);
                  };
                if (n2gated > 0)
                  {
                    HK2->Fill ((double) KK, (double) HH, 1);
                  };
                if (n3gated > 0)
                  {
                    HK3->Fill ((double) KK, (double) HH, 1);
                  };

              };


          /* gamma gamma matrix */

          if (grh->ngam >= 2)
            for (k = 0; k < grh->ngam; k++)
              if (grh->gr[k].tracked)
                if (grh->gr[k].ndet >= Pars.ndetlimlo && grh->gr[k].ndet <= Pars.ndetlimhi)
                  if (grh->gr[k].fom >= Pars.fomlo[grh->gr[k].ndet] && grh->gr[k].fom <= Pars.fomhi[grh->gr[k].ndet])
                    for (l = k + 1; l < grh->ngam; l++)
                      if (grh->gr[l].tracked)
                        if (grh->gr[l].ndet >= Pars.ndetlimlo && grh->gr[l].ndet <= Pars.ndetlimhi)
                          if (grh->gr[l].fom >= Pars.fomlo[grh->gr[l].ndet]
                              && grh->gr[l].fom <= Pars.fomhi[grh->gr[l].ndet])
                            {

                              d1 = (double) grh->gr[k].esum;
                              if (d1 < 0)
                                d1 = 0;
                              if (d1 > Pars.GGMAX)
                                d1 = Pars.GGMAX;

                              d2 = (double) grh->gr[l].esum;
                              if (d2 < 0)
                                d2 = 0;
                              if (d2 > Pars.GGMAX)
                                d2 = Pars.GGMAX;

                              gg->Fill (d1, d2, 1.0);
                              gg->Fill (d2, d1, 1.0);


                            };


          /* firstScatAng */

          for (k = 0; k < grh->ngam; k++)
            if (grh->gr[k].tracked)
              if (grh->gr[k].ndet >= Pars.ndetlimlo && grh->gr[k].ndet <= Pars.ndetlimhi)
                if (grh->gr[k].fom >= Pars.fomlo[grh->gr[k].ndet] && grh->gr[k].fom <= Pars.fomhi[grh->gr[k].ndet])
                  if (grh->gr[k].ndet > 0)
                    {

                      /* find first scattering angle */

                      xx = grh->gr[k].x0 - Pars.target_x;
                      yy = grh->gr[k].y0 - Pars.target_y;
                      zz = grh->gr[k].z0 - Pars.target_z;
                      if (Pars.CurEvNo <= Pars.NumToPrint)printf("first interaction (.x0): %4.2f  %4.2f %4.2f\n", xx,yy,zz);
                      dd = xx * xx + yy * yy + zz * zz;
                      dd = sqrtf (dd);
                      x1 = xx / dd;
                      y1 = yy / dd;
                      z1 = zz / dd;

                      xx = grh->gr[k].x1;
                      yy = grh->gr[k].y1 ;
                      zz = grh->gr[k].z1 ;
                      if (Pars.CurEvNo <= Pars.NumToPrint)printf("second interaction(.x1): %4.2f  %4.2f %4.2f\n", xx,yy,zz);
                      xx -= grh->gr[k].x0;
                      yy -= grh->gr[k].y0;
                      zz -= grh->gr[k].z0;
                      if (Pars.CurEvNo <= Pars.NumToPrint)printf("second vector: %4.2f  %4.2f %4.2f\n", xx,yy,zz);
                      dd = xx * xx + yy * yy + zz * zz;
                      dd = sqrtf (dd);
                      if (Pars.CurEvNo <= Pars.NumToPrint)printf("second vector length: %4.2f \n", dd);
                      x2 = xx/dd;
                      y2 = yy/dd;
                      z2 = zz/dd;
                      if (Pars.CurEvNo <= Pars.NumToPrint)printf("second vector: %4.2f  %4.2f %4.2f (normalized)\n", x2,y2,z2);

                      /* scatter range here */

                      if (dd<200)
                      firstScatRange->Fill (grh->gr[k].esum-grh->gr[k].e0, (double) dd);

                      /* first scatter angle */

                      dd = x1 * x2 + y1 * y2 + z1 * z2;
                      if (dd >= -1.0 && dd <= 1.0)
                        {
                          dd = acosf (dd);
                          dd *= 180.0/M_PI;
                          if (Pars.CurEvNo <= Pars.NumToPrint)printf("ang =%4.2f\n", dd);
                          firstScatAng->Fill (grh->gr[k].esum, dd);
                        };

                    };

        };                      /* if (GEB_event->ptgd[i]->type == GEB_TYPE_TRACK) */

    };


/* done */

        if (Pars.CurEvNo <= Pars.NumToPrint)
          printf ("exit bin_mode1\n");

  return (0);

}
