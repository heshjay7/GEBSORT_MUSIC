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

#define PPLO 1
#define PPHI 2
#define BBLO 3
#define BBHI 4

/* pointers to ROOT spectra */


TH2F *mpolang;
TH2F *mpolang1;
TH2F *mpolang2;
TH2F *mpolang3;

/* parameters */

extern PARS Pars;
extern EXCHANGE exchange;

#define NBINS 180

/*-----------------------------------------------------*/

int
exit_angdis ()
{


  return (0);

};

/*-----------------------------------------------------*/

int
sup_angdis ()
{
  /* declarations */

  char str1[STRLEN], str2[STRLEN];
  TH1D *mkTH1D (char *, char *, int, double, double);
  TH2F *mkTH2F (char *, char *, int, double, double, int, double, double);


  sprintf (str1, "polAng");
  sprintf (str2, "polAng");
  mpolang = mkTH2F (str1, str2, NBINS+1,0,NBINS, Pars.GGMAX, 1, Pars.GGMAX);
  sprintf (str1, "polar angle [deg]");
  mpolang->SetXTitle (str1);
  sprintf (str1, "energy [kev]");
  mpolang->SetYTitle (str1);

  sprintf (str1, "polAng1");
  sprintf (str2, "polAng1");
  mpolang1 = mkTH2F (str1, str2, NBINS+1,0,NBINS, Pars.GGMAX, 1, Pars.GGMAX);
  sprintf (str1, "polar angle [deg]");
  mpolang1->SetXTitle (str1);
  sprintf (str1, "energy [kev]");
  mpolang1->SetYTitle (str1);

  sprintf (str1, "polAng2");
  sprintf (str2, "polAng2");
  mpolang2 = mkTH2F (str1, str2, NBINS+1,0,NBINS, Pars.GGMAX, 1, Pars.GGMAX);
  sprintf (str1, "polar angle [deg]");
  mpolang2->SetXTitle (str1);
  sprintf (str1, "energy [kev]");
  mpolang2->SetYTitle (str1);

  sprintf (str1, "polAng3");
  sprintf (str2, "polAng3");
  mpolang3 = mkTH2F (str1, str2, NBINS+1,0,NBINS, Pars.GGMAX, 1, Pars.GGMAX);
  sprintf (str1, "polar angle [deg]");
  mpolang3->SetXTitle (str1);
  sprintf (str1, "energy [kev]");
  mpolang3->SetYTitle (str1);


  /* list what we have */

  printf (" we have defined the following spectra:\n");


  Pars.wlist = gDirectory->GetList ();
  Pars.wlist->Print ();

  /* make sure bin_mode1 is on so the */
  /* Doppler correction is taken care of */

  if (Pars.do_bin_mode1 != 1)
    {
      printf ("For the doppler correction, you must have bin_mode1 enabled.\n");
      printf ("bin_linpol does not Doppler correct the mode1 data.\n");
      printf ("quit!\n");
      exit (1);
    };

  return(0);

};

/* ----------------------------------------------------------------- */

int
bin_angdis (GEB_EVENT * GEB_event)
{
  int i, j;
  char str[128];
  TRACKED_GAMMA_HIT *grh;
  double rr, dp, pang;

  int GebTypeStr (int type, char str[]);
  int print_tracked_gamma_rays (FILE *, TRACKED_GAMMA_HIT *);

  for (i = 0; i < GEB_event->mult; i++)
    {

      if (Pars.CurEvNo <= Pars.NumToPrint)
        {
          GebTypeStr (GEB_event->ptgd[i]->type, str);
          printf ("bin_angdis, %2i> %2i, %s, TS=%lli\n", i, GEB_event->ptgd[i]->type, str,
                  GEB_event->ptgd[i]->timestamp);
        }

      if (GEB_event->ptgd[i]->type == GEB_TYPE_TRACK)
        {

          grh = (TRACKED_GAMMA_HIT *) GEB_event->ptinp[i];
          if (Pars.CurEvNo <= Pars.NumToPrint)
            print_tracked_gamma_rays (stdout, grh);

          /*------------------------------------*/
          /* update angular distribution matrix */
          /*------------------------------------*/

          for (j = 0; j < grh->ngam; j++)
            {

              /* find normalized vector from target position */

              rr = (grh->gr[j].x0 - Pars.target_x) * (grh->gr[j].x0 - Pars.target_x)
                + (grh->gr[j].y0 - Pars.target_y) * (grh->gr[j].y0 - Pars.target_y)
                + (grh->gr[j].z0 - Pars.target_z) * (grh->gr[j].z0 - Pars.target_z);
              rr = sqrtf (rr);

              /* dot with beam direction (which is normalized already) */

              dp = ((grh->gr[j].x0 - Pars.target_x) * Pars.beamdir[0] +
                    (grh->gr[j].y0 - Pars.target_y) * Pars.beamdir[1] +
                    (grh->gr[j].z0 - Pars.target_z) * Pars.beamdir[2]) / rr;

              /* find polar angle and bin it */

              if (dp < -1.0)
                dp = -1.0;
              if (dp > 1.0)
                dp = 1.0;
              pang = acosf (dp);
              pang *= (180 / M_PI);
              if (grh->gr[j].ndet >= Pars.ndetlimlo && grh->gr[j].ndet <= Pars.ndetlimhi)
                if (grh->gr[j].fom >= Pars.fomlo[grh->gr[j].ndet] && grh->gr[j].fom <= Pars.fomhi[grh->gr[j].ndet])
                  if (grh->gr[j].esum < Pars.GGMAX && grh->gr[j].esum >= 1)
                    if (pang > 0 && pang < NBINS)
                      {
                        mpolang->Fill ((double) pang, (double) grh->gr[j].esum, 1);
                        if (exchange.ngates >= 1)
                          mpolang1->Fill ((double) pang, (double) grh->gr[j].esum, 1);
                        if (exchange.ngates >= 2)
                          mpolang2->Fill ((double) pang, (double) grh->gr[j].esum, 1);
                        if (exchange.ngates >= 3)
                          mpolang3->Fill ((double) pang, (double) grh->gr[j].esum, 1);
                      };

            };

        };

    };


  /* done */

  return (0);
};
