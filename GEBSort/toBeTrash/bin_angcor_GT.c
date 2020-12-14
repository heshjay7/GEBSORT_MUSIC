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

#define LOQ 25

/* pointers to ROOT spectra */

TH3F *angcor_cube;
TH3F *angcor_cube_oo;

/* parameters */

extern PARS Pars;
extern EXCHANGE exchange;


int nn_sp[100];

/*-----------------------------------------------------*/

int
exit_angcor_GT ()
{

  int i;
  double mm = 0;

  for (i = 0; i < 100; i++)
    mm += nn_sp[i];
  mm /= 100;

  printf ("\nexit_angcor: type statistics\n");
  for (i = 0; i < 100; i++)
    if (nn_sp[i] > 0)
      printf ("nn_sp[%2i]=%10i, %5.1f%%\n", i, nn_sp[i], (float) nn_sp[i] / mm);

  return (0);

};

/*-----------------------------------------------------*/

int
sup_angcor_GT ()
{
  /* declarations */

  int i;
  double d1;
  TH1D *mkTH1D (char *, char *, int, double, double);
  TH2F *mkTH2F (char *, char *, int, double, double, int, double, double);
  TH3F *mkTH3F (char *, char *, int, double, double, int, double, double, int, double, double);

  /* the catch all 3D arrays */

  angcor_cube = mkTH3F ((char *) "angcor_cube", (char *) "angcor_cube",
                        Pars.GGMAX, 1, Pars.GGMAX, Pars.GGMAX, 1, Pars.GGMAX, 36, 0, 180);
  angcor_cube_oo = mkTH3F ((char *) "angcor_cube_oo", (char *) "angcor_cube_oo",
                           Pars.GGMAX, 1, Pars.GGMAX, Pars.GGMAX, 1, Pars.GGMAX, 36, 0, 180);

  d1 = Pars.GGMAX * Pars.GGMAX * 36 * sizeof (float);
  printf ("angcur cubes are %10.0f bytes long, that is %10.3f GBytes\n", d1, d1 / 1024 / 1024 / 1024);
  if (d1 > 1073741822)
    {
      printf ("ROOT does not like object larger than 1073741822, you have %10.0f, QUIT\n", d1);
      fflush (stdout);
      exit (1);
    };

  /* list what we have */

  printf (" we have defined the following spectra:\n");


  Pars.wlist = gDirectory->GetList ();
  Pars.wlist->Print ();


  printf ("Note: always run through bin_mode1 before \n");
  printf ("using this function so that dopper \n");
  printf ("corrections can be handled properly\n");

  for (i = 0; i < 100; i++)
    nn_sp[i] = 0;

  return (0);

};

/* ----------------------------------------------------------------- */

int
bin_angcor_GT (GEB_EVENT * GEB_event)
{

  /* declarations */

  char str[128];
  int i, j, k = 0, l, nn;
  double xx[MAX_GAMMA_RAYS], yy[MAX_GAMMA_RAYS], zz[MAX_GAMMA_RAYS];
  TRACKED_GAMMA_HIT *grh;
  double d1, d2;
  double x1, y1, z1;
  double x2, y2, z2;

  static double xx_oo[LOQ][MAX_GAMMA_RAYS], yy_oo[LOQ][MAX_GAMMA_RAYS], zz_oo[LOQ][MAX_GAMMA_RAYS],
    ee_oo[LOQ][MAX_GAMMA_RAYS];
  static int nn_oo[LOQ];


  /* prototypes */

  int GebTypeStr (int type, char str[]);
  int print_tracked_gamma_rays (FILE *, TRACKED_GAMMA_HIT *);
  double dotproductangle (double, double, double, double, double, double);
  double rad2deg (double);
  double unitvector (double *, double *, double *);
  int crossprod (double, double, double, double, double, double, double *, double *, double *);



  if (Pars.CurEvNo <= Pars.NumToPrint)
    printf ("entered bin_angcor_GT: mult=%i\n", GEB_event->mult);

  /* list event fundamentals */

  if (Pars.CurEvNo <= Pars.NumToPrint)
    {
      printf ("GEB_event->mult=%i\n", GEB_event->mult);
      for (i = 0; i < GEB_event->mult; i++)
        {
          GebTypeStr (GEB_event->ptgd[i]->type, str);
          printf ("%2i> %s, TS=%lli; ", i, str, GEB_event->ptgd[i]->timestamp);
          if (GEB_event->ptgd[i]->type == GEB_TYPE_TRACK)
            {
              grh = (TRACKED_GAMMA_HIT *) GEB_event->ptinp[i];
              printf ("gammas: ");
              for (j = 0; j < grh->ngam; j++)
                printf ("[%i] %5.1f ", j, grh->gr[j].esum);
            };
          printf ("\n");
        };
    };

  /*----------------------*/
  /* fill the angcor_cube */
  /*----------------------*/

  /* count number of tracked headers */

  nn = 0;
  for (i = 0; i < GEB_event->mult; i++)
    if (GEB_event->ptgd[i]->type == GEB_TYPE_TRACK)
      nn++;

  /* we don't want to deal with double loops */
  /* we can, but it is more messy... */
  /* there should only be one tracked head/pay */

  if (nn != 1)
    {
      if (Pars.CurEvNo <= Pars.NumToPrint)
        printf ("# tracked headers=%i, quit\n", nn);
      return (0);
    };
  if (Pars.CurEvNo <= Pars.NumToPrint)
    printf ("# tracked headers=%i, proceed\n", nn);

  for (i = 0; i < GEB_event->mult; i++)
    if (GEB_event->ptgd[i]->type == GEB_TYPE_TRACK)
      {

        grh = (TRACKED_GAMMA_HIT *) GEB_event->ptinp[i];

        /* if only one gamma ray, no need to continue */

        if (grh->ngam < 2)
          return (0);

        /* debug list */

        if (Pars.CurEvNo <= Pars.NumToPrint)
          for (j = 0; j < grh->ngam; j++)
            {
              printf ("{%i} %5.1f ", j, grh->gr[j].esum);
              if (Pars.CurEvNo <= Pars.NumToPrint)
                printf ("\n");
            };

        /* find unit vectors */

        for (j = 0; j < grh->ngam; j++)
          {
            xx[j] = grh->gr[j].x0;
            yy[j] = grh->gr[j].y0;
            zz[j] = grh->gr[j].z0;
            unitvector (&xx[j], &yy[j], &zz[j]);
          };

        /* update angcor_cube */

        for (j = 0; j < grh->ngam; j++)
          if (grh->gr[j].esum > 0 && grh->gr[j].esum < Pars.GGMAX)
            for (k = j + 1; k < grh->ngam; k++)
              if (grh->gr[k].esum > 0 && grh->gr[k].esum < Pars.GGMAX)

                {
                  if (Pars.angcor_useplaneang)
                    {
                      /* use angle between gamma-zaxis planes (for in-beam data, DCO) */

                      /* find the two normal vectors to the planes */

                      crossprod ((double) 0, (double) 0, (double) 1, xx[j], yy[j], zz[j], &x1, &y1, &z1);
                      d1 = unitvector (&x1, &y1, &z1);
                      crossprod ((double) 0, (double) 0, (double) 1, xx[k], yy[k], zz[k], &x2, &y2, &z2);
                      d1 = unitvector (&x2, &y2, &z2);

                      /* find angle bt the planes */

                      d1 = dotproductangle (x1, y1, z1, x2, y2, z2);

                      /* for the log file */

                      if (Pars.CurEvNo <= Pars.NumToPrint)
                        {
                          printf ("using angle bt gam/z-axis planes: %6.1f deg\n", rad2deg (d1));
                          d2 = dotproductangle (xx[j], yy[j], zz[j], xx[k], yy[k], zz[k]);
                          printf ("angle bt gamma rays: %6.1f deg\n", rad2deg (d2));
                          d2 = dotproductangle ((double) 0, (double) 0, (double) 1, xx[j], yy[j], zz[j]);
                          printf ("gamma ray 1 polar angle: %6.1f deg\n", rad2deg (d2));
                          d2 = dotproductangle ((double) 0, (double) 0, (double) 1, xx[k], yy[k], zz[k]);
                          printf ("gamma ray 2 polar angle: %6.1f deg\n", rad2deg (d2));
                        };

                    }
                  else
                    {
                      /* use angle between the gamma rays (for sources data) */

                      d1 = dotproductangle (xx[j], yy[j], zz[j], xx[k], yy[k], zz[k]);
                      if (Pars.CurEvNo <= Pars.NumToPrint)
                        printf ("using angle bt gamma rays: %6.1f deg\n", rad2deg (d1));
                    };

                  d1 = rad2deg (d1);
                  angcor_cube->Fill (grh->gr[j].esum, grh->gr[k].esum, d1, 1);
                  angcor_cube->Fill (grh->gr[k].esum, grh->gr[j].esum, d1, 1);
                  if (Pars.CurEvNo <= Pars.NumToPrint)
                    {
                      printf ("angcor_cube j: %f %f %f %f\n", xx[j], yy[j], zz[j], grh->gr[j].esum);
                      printf ("angcor_cube k: %f %f %f %f\n", xx[k], yy[k], zz[k], grh->gr[k].esum);
                      printf ("angcor_cube: update %f %f %f\n", grh->gr[j].esum, grh->gr[k].esum, d1);
                    };
                };

        /* update angcor_cube_oo */
        /* note: this event against everything in the old queue */
        /* think already symmetrized!? */

        for (l = 0; l < LOQ; l++)
          for (j = 0; j < grh->ngam; j++)
            if (grh->gr[j].esum > 0 && grh->gr[j].esum < Pars.GGMAX)
              for (k = 0; k < nn_oo[l]; k++)
                if (ee_oo[l][k] > 0.0 && ee_oo[l][k] < Pars.GGMAX)
                  {
                    if (Pars.angcor_useplaneang)
                      {
                        /* use angle between gamma-zaxis planes (for in-beam data, DCO) */

                        /* find the two normal vectors to the planes */

                        crossprod ((double) 0, (double) 0, (double) 1, xx[j], yy[j], zz[j], &x1, &y1, &z1);
                        d1 = unitvector (&x1, &y1, &z1);
                        crossprod ((double) 0, (double) 0, (double) 1, xx_oo[l][k], yy_oo[l][k], zz_oo[l][k], &x2, &y2,
                                   &z2);
                        d1 = unitvector (&x2, &y2, &z2);

                        /* find angle bt the planes */

                        d1 = dotproductangle (x1, y1, z1, x2, y2, z2);
                        if (Pars.CurEvNo <= Pars.NumToPrint)
                          printf ("using angle bt gam/z axis planes: %6.1f deg\n", rad2deg (d1));

                      }
                    else
                      {
                        /* use angle between the gamma rays (for sources data) */

                        d1 = dotproductangle (xx[j], yy[j], zz[j], xx_oo[l][k], yy_oo[l][k], zz_oo[l][k]);
                        if (Pars.CurEvNo <= Pars.NumToPrint)
                          printf ("using angle bt gamma rays: %6.1f deg\n", rad2deg (d1));
                      };
                    d1 = rad2deg (d1);
                    angcor_cube_oo->Fill (grh->gr[j].esum, ee_oo[l][k], d1, 1);
                    if (Pars.CurEvNo <= Pars.NumToPrint)
                      {
                        printf ("angcor_cube_oo j: %f %f %f %f\n", xx[j], yy[j], zz[j], grh->gr[j].esum);
                        printf ("angcor_cube_oo k: %f %f %f %f\n", xx_oo[l][k], yy_oo[l][k], zz_oo[l][k], ee_oo[l][k]);
                        printf ("angcor_cube_oo: update %f %f %f\n", grh->gr[j].esum, ee_oo[l][k], d1);
                      };
                  };

        /* store current to be next old */
        /* already unit vectors */

        /* move everyone one down */

        l = LOQ - 1;
        while (l > 0)
          {
            nn_oo[l] = nn_oo[l - 1];
            for (k = 0; k < nn_oo[l - 1]; k++)
              {
                xx_oo[l][k] = xx_oo[l - 1][k];
                yy_oo[l][k] = yy_oo[l - 1][k];
                zz_oo[l][k] = zz_oo[l - 1][k];
                ee_oo[l][k] = ee_oo[l - 1][k];
              };
            l--;
          };

        /* store this event */
        /* at the top */

        nn_oo[0] = grh->ngam;
        for (k = 0; k < nn_oo[0]; k++)
          {
            xx_oo[0][k] = xx[k];
            yy_oo[0][k] = yy[k];
            zz_oo[0][k] = zz[k];
            ee_oo[0][k] = grh->gr[k].esum;
          };

      };

  /* done */

  if (Pars.CurEvNo <= Pars.NumToPrint)
    printf ("exit bin_angcor_GT\n");

  return (0);

}
