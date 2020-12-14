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
#include "gsang.h"

#define NBINS 180

/* pointers to ROOT spectra */

TH1D *Gsangdiff;
TH3F *Angcor_cube;
TH3F *Angcor_cube_oo;
TH2F *SMAP_DGS;

/* parameters */

extern DGSEVENT DGSEvent[MAXCOINEV];
extern int ng;

extern PARS Pars;
extern EXCHANGE exchange;

#define LOQ 15

int Nn_sp[100];
double angdif[MAX_GES + 1][MAX_GES + 1];
long long int dgshit[MAX_GES + 1][MAX_GES + 1];
double pol_angdis[MAX_GES + 1], azi_angdis[MAX_GES + 1];
long long int angdis_hitp[MAX_GES + 1];

/*-----------------------------------------------------*/

int
exit_angcor_DGS ()
{

  int i, nn = 0;
  double mm = 0, d1;

  for (i = 0; i < 100; i++)
    mm += Nn_sp[i];
  mm /= 100;

  printf ("\nexit_angcor: type statistics\n");
  for (i = 0; i < 100; i++)
    if (Nn_sp[i] > 0)
      printf ("Nn_sp[%2i]=%10i, %5.1f%%\n", i, Nn_sp[i], (float) Nn_sp[i] / mm);

  /* special angcor hit pattern */

  mm = 0;
  nn = 0;
  for (i = 1; i < MAX_GES; i++)
    if (angdis_hitp[i] > 0)
      {
        nn++;
        mm += angdis_hitp[i];
      };
  mm /= nn;

  printf ("\n");
  printf ("special angcor counting:\n");
  printf ("\n");
  for (i = 1; i < MAX_GES; i++)
    if (angdis_hitp[i] > 0)
      {
        printf ("det %3i: ", i);
        printf (" counts %10lli ", angdis_hitp[i]);
        d1 = 100.0 * angdis_hitp[i] / mm;
        printf (" %7.1f %% ", d1);
        if (d1 < 85.0)
          printf (" too low");
        if (d1 > 115.0)
          printf (" too high");
        printf ("\n");
      }
  printf ("\n");

  return (0);

};

/*-----------------------------------------------------*/

int
sup_angcor_DGS ()
{
  /* declarations */

  int i, k, l;
  char str1[STRLEN], str2[STRLEN];
  double d1, d2;
  unsigned int seed;

  double xx[MAX_GES + 1], yy[MAX_GES + 1], zz[MAX_GES + 1];

  FILE *fp;

  /* prototypes */

  TH1D *mkTH1D (char *, char *, int, double, double);
  TH2F *mkTH2F (char *, char *, int, double, double, int, double, double);
  TH3F *mkTH3F (char *, char *, int, double, double, int, double, double, int, double, double);
  double dotproductangle (double, double, double, double, double, double);
  double rad2deg (double);
  int check_unitvector (double, double, double);
  int get_a_seed (unsigned int *);

  /* make sure bin_dgs is on for doppler correction etc */

  if (Pars.do_bin_angcor_DGS != 1)
    {
      printf ("you must have bin_dgs on so the data\n");
      printf ("seen in bin_angcor_DGS is energy \n");
      printf ("calibrated and Doppler corrected, quit\n");
      exit (1);
    };

  Gsangdiff = mkTH1D ((char *) "Gsangdiff", (char *) "GS detector angcor angles", NBINS, 0, 180);
  Gsangdiff->SetXTitle ((char *) "detno");
  Gsangdiff->SetYTitle ((char *) "yield");

  SMAP_DGS = mkTH2F ((char *) "SMAP_DGS", (char *) "SMAP_DGS (schematic det sizes)", 361, 0, 360, 181, 0, 180);
  SMAP_DGS->SetXTitle ((char *) "aximuth");
  SMAP_DGS->SetYTitle ((char *) "polar");

  /* list what we have */

  printf (" we have defined the following spectra:\n");

  Pars.wlist = gDirectory->GetList ();
  Pars.wlist->Print ();

  for (i = 0; i < 100; i++)
    Nn_sp[i] = 0;

  /* find normal vectors for all detectors */

  /* The h file has the 1-110 detectors mapped into  */
  /* angle lookup tables that go from 0-109. Here we  */
  /* quitly take care of that when we create vectors  */
  /* so that we can ignore the offset later on and refer */
  /* to the detectors as 1,2,....109,110 */


#if(0)
  for (i = 1; i <= MAX_GES; i++)
    {
      printf ("ge%3i: ", i);
      printf ("pol= %7.3f azi= %8.3f \n", angtheta[i - 1], angphi[i - 1]);
      printf ("cos(azi)= %7.4f sin (azi)= %7.4f\n", cos (angphi[i - 1]), sin (angphi[i - 1]));
      if (angphi[i - 1] > 180.0)
        angphi[i - 1] = angphi[i - 1] - 360;
      printf ("cos(azi)= %7.4f sin (azi)= %7.4f\n", cos (angphi[i - 1]), sin (angphi[i - 1]));
      /* ^^^^ clearly shows this does not hold true */
    };
  exit (0);
#endif

  for (i = 1; i <= MAX_GES; i++)
    {

      printf ("ge%3i: pol= %7.3f azi= %8.3f ", i, angtheta[i - 1], angphi[i - 1]);
      d1 = 180.0 + (angphi[i - 1] - 180.0) * sin ((double) angtheta[i - 1] * M_PI / 180.0);
      printf (" (SMAP %8.3f sin: %f); ", d1, sin ((double) angtheta[i - 1] * M_PI / 180.0));
      xx[i] = sin ((double) angtheta[i - 1] * M_PI / 180.0) * cos ((double) angphi[i - 1] * M_PI / 180.0);
      yy[i] = sin ((double) angtheta[i - 1] * M_PI / 180.0) * sin ((double) angphi[i - 1] * M_PI / 180.0);;
      zz[i] = cos ((double) angtheta[i - 1] * M_PI / 180.0);
      check_unitvector (xx[i], yy[i], zz[i]);
      printf (" ( %7.4f  %7.4f  %7.4f ) \n", xx[i], yy[i], zz[i]);
      pol_angdis[i] = (double) angtheta[i - 1] * M_PI / 180.0;
      azi_angdis[i] = (double) angphi[i - 1] * M_PI / 180.0;
    };

  /* find and plot angle differences */
  /* make symmetric and include zero diagonal */

  for (k = 1; k <= MAX_GES; k++)
    for (l = 1; l <= MAX_GES; l++)
      {
        d1 = dotproductangle (xx[k], yy[k], zz[k], xx[l], yy[l], zz[l]);
        d1 = rad2deg (d1);
        angdif[k][l] = d1;
//        printf ("%3i %3i - %7.2f\n", k, l, angdif[k][l]);
        Gsangdiff->Fill (angdif[k][l], 1);
      };

#if(1)

  /* one could also use this formula */
  /* but I like to see the vectors */
  /* https://www.physicsforums.com/threads/the-angle-between-two-vectors.856065/ */

  for (k = 1; k <= MAX_GES; k++)
    for (l = k + 1; l <= MAX_GES; l++)
      {
        d1 = dotproductangle (xx[k], yy[k], zz[k], xx[l], yy[l], zz[l]);
        d1 *= 180 / M_PI;
        d2 = cos (pol_angdis[k]) * cos (pol_angdis[l]);
        d2 += sin (pol_angdis[k]) * sin (pol_angdis[l]) * cos (azi_angdis[k] - azi_angdis[l]);
        d2 = acos (d2);
        d2 *= 180 / M_PI;
//        printf ("%3i %3i - %7.2f %7.2f (alt) diff= %f\n", k, l, d1, d2, d1-d2);
        assert (fabs (d1 - d2) < 0.0001);
      };

#endif

  /* list in file (only one way) */

  fp = fopen ("GS_ancor_angles.txt", "w");
  for (i = 1; i <= MAX_GES; i++)
    {
      fprintf (fp, "ge%3i: pol= %7.3f azi= %7.3f ; ", i, angtheta[i - 1], angphi[i - 1]);
      fprintf (fp, " ( %7.4f  %7.4f  %7.4f ) \n", xx[i], yy[i], zz[i]);
    }
  fprintf (fp, "----\n");
  for (k = 1; k <= MAX_GES; k++)
    for (l = k + 1; l <= MAX_GES; l++)
      {
        fprintf (fp, "%3.3i %3.3i - %7.2f\n", k, l, angdif[k][l]);
        if (fabs (angdif[k][l] - angdif[l][k]) > 0.0001)
          {
            printf ("ooops, %3i %3i - %7.2f != %7.2f \n", k, l, angdif[k][l], angdif[l][k]);
            fflush (stdout);
            exit (1);
          }
      };
  fclose (fp);

  /* checked some polar angles by hand, were OK */

  printf ("sort -n -k 4 GS_ancor_angles.txt | more \n");
  printf ("wc GS_ancor_angles.txt; better be 110*109/2 = 5995 \n");
  printf ("sort -n -k 4 GS_ancor_angles.txt | awk '{print $4}' | uniq | wc ; 275 unique angles\n");
  printf ("\n");
  printf ("we take the angle information from gsang.h\n");
  printf ("we do not know where that file came from... I.Y. Lee?\n");
  printf ("\n");

  /* the catch all 3D arrays */

  Angcor_cube = mkTH3F ((char *) "Angcor_cube", (char *) "Angcor_cube",
                        Pars.GGMAX, 1, Pars.GGMAX, Pars.GGMAX, 1, Pars.GGMAX, 36, 0, 180);
  Angcor_cube_oo = mkTH3F ((char *) "Angcor_cube_oo", (char *) "Angcor_cube_oo",
                           Pars.GGMAX, 1, Pars.GGMAX, Pars.GGMAX, 1, Pars.GGMAX, 36, 0, 180);

  d1 = Pars.GGMAX * Pars.GGMAX * 36 * sizeof (float);
  printf ("angcor cubes are %10.0f bytes long, that is %10.3f GBytes\n", d1, d1 / 1024 / 1024 / 1024);
  if (d1 > 1073741822)
    {
      printf ("ROOT does not like object larger than 1073741822, you have %10.0f, reduce gglen, QUIT\n", d1);
      fflush (stdout);
      exit (1);
    };

  get_a_seed (&seed);
  srand (seed);

  for (i = 1; i < MAX_GES; i++)
    angdis_hitp[i] = 0;

  /* done */

  return (0);

};

/* ----------------------------------------------------------------- */

int
bin_angcor_DGS (GEB_EVENT * GEB_event)
{

  /* declarations */

  char str[128];
  int i, j, k = 0, l, m, nn;
  double ee[MAX_GAMMA_RAYS];
  int id[MAX_GAMMA_RAYS];
  float sX, sY, RAD2DEG = 0.017453;
  float dx, dy, rr;

  static double ee_oo[LOQ][MAX_GAMMA_RAYS];
  static int id_oo[LOQ][MAX_GAMMA_RAYS];
  static int nn_oo[LOQ];

  /* prototypes */

  int GebTypeStr (int type, char str[]);
  int print_tracked_gamma_rays (FILE *, TRACKED_GAMMA_HIT *);
  double dotproductangle (double, double, double, double, double, double);
  double rad2deg (double);
  double unitvector (double *, double *, double *);

  if (Pars.CurEvNo <= Pars.NumToPrint)
    printf ("entered bin_angcor_DGS: mult=%i\n", GEB_event->mult);

  /* list event fundamentals */

  if (Pars.CurEvNo <= Pars.NumToPrint)
    {
      printf ("GEB_event->mult=%i\n", GEB_event->mult);
      for (i = 0; i < GEB_event->mult; i++)
        {
          GebTypeStr (GEB_event->ptgd[i]->type, str);
          printf ("%2i> %s, TS=%lli; \n", i, str, GEB_event->ptgd[i]->timestamp);
        };
      printf ("\n");
      for (i = 0; i < ng; i++)
        printf ("DGS: tid %i energy %f\n", DGSEvent[i].tid, DGSEvent[i].ehi);
    };

  /* check for bad IDs */

  for (i = 0; i < ng; i++)
    if (DGSEvent[i].tpe == GE)
    if (DGSEvent[i].tid < 1 || DGSEvent[i].tid > MAX_GES)
      {
//        printf ("bin_dgs: %i> tid %i energy %f\n", i, DGSEvent[i].tid, DGSEvent[i].ehi);
        return (0);
      };


  /* assemble good event so we */
  /* do not have to check later */

  nn = 0;
  for (j = 0; j < ng; j++)
    if (DGSEvent[j].tpe == GE)
      if (DGSEvent[j].flag == 0)
        if (Pars.enabled[DGSEvent[j].tid])
          if (DGSEvent[j].ehi < Pars.GGMAX)
            if (DGSEvent[j].ehi > 0)
              if (DGSEvent[j].tid >= 1)
                if (DGSEvent[j].tid < MAX_GES)
                  {
                    ee[nn] = DGSEvent[j].ehi;
                    id[nn] = DGSEvent[j].tid;
                    nn++;
                  };

  /* we need at least two to continue */

  if (nn < 2)
    return (0);

  if (Pars.CurEvNo <= Pars.NumToPrint)
    {
      printf ("ev# %5i: ", Pars.CurEvNo);
      for (j = 0; j < nn; j++)
        printf ("%6.1f [%3i] ", ee[j], id[j]);
      printf ("\n");
    };

  /*----------------------*/
  /* fill the angcor_cube */
  /*----------------------*/

  for (k = 0; k < nn; k++)
    for (l = k + 1; l < nn; l++)
      {
        Angcor_cube->Fill (ee[k], ee[l], angdif[id[k]][id[l]]);
        Angcor_cube->Fill (ee[l], ee[k], angdif[id[l]][id[k]]);
        if (Pars.CurEvNo <= Pars.NumToPrint)
          printf ("cube update: %6.1f [%3i] %6.1f [%3i] ang= %5.1f\n", ee[k], id[k], ee[l], id[l],
                  angdif[id[k]][id[l]]);
      };

  /* SMAP coordinates (limit updates) */

  if (Pars.CurEvNo <= 1000000)
    for (k = 0; k < nn; k++)
      {
        sX = M_PI + (azi_angdis[id[k]] - M_PI) * sin (pol_angdis[id[k]]);
        sX /= RAD2DEG;
        sY = pol_angdis[id[k]] / RAD2DEG;
        if (Pars.CurEvNo <= Pars.NumToPrint)
          {
            printf ("tpe=%i, tid=%i\n", DGSEvent[k].tpe, id[k]);
            printf ("azi/pol: %5.3f %5.3f [rad] ", azi_angdis[id[k]], pol_angdis[id[k]]);
            printf ("%5.1f %5.1f [deg]\n", azi_angdis[id[k]] / RAD2DEG, pol_angdis[id[k]] / RAD2DEG);
            printf ("sX=%f, sY=%f\n", sX, sY);
          };
        rr = 15;
        while (rr > 5)
          {
            dx = 10.0 * ((float) rand () / RAND_MAX - 0.5);
            dy = 10.0 * ((float) rand () / RAND_MAX - 0.5);
            rr = sqrtf (dx * dx + dy * dy);
          };
        sX += dx;
        sY += dy;
        if (Pars.CurEvNo <= Pars.NumToPrint)
          printf ("dx = %f dy = %f\n", dx, dy);

        if (sX > 0 && sX < 360)
          if (sY > 0 && sY < 180)
            SMAP_DGS->Fill (sX, sY, 1);

        /* count hits in detectors that get through */

        angdis_hitp[id[k]]++;

      };


  /*-------------------------*/
  /* fill the angcor_cube_oo */
  /*-------------------------*/

  for (m = 0; m < LOQ; m++)
    for (k = 0; k < nn; k++)
      for (l = 0; l < nn_oo[m]; l++)
        if (nn_oo[m]>0)
        {
          Angcor_cube_oo->Fill (ee[k], ee_oo[m][l], angdif[id[k]][id_oo[m][l]]);
          if (Pars.CurEvNo <= Pars.NumToPrint)
            printf ("cube_oo update: %6.1f [%3i] %6.1f [%3i] ang= %5.1f\n", ee[k], id[k], ee_oo[m][l], id_oo[m][l],
                    angdif[id[k]][id_oo[m][l]]);
        };

  if (Pars.CurEvNo <= Pars.NumToPrint)
    {
      printf ("current mix queue:\n");
      for (m = 0; m < LOQ; m++)
        {
          printf ("%2i [l=%2i]: ", m, nn_oo[m]);
          for (k = 0; k < nn_oo[m]; k++)
            printf ("%6.1f [%3i], ", ee_oo[m][k], id_oo[m][k]);
          printf ("\n");
        };
    };

  /* move everyone one down in the back mix list. */
  /* This could be made more efficient with a */
  /* circular queue and a pointer... */

  l = LOQ - 1;
  while (l > 0)
    {
      nn_oo[l] = nn_oo[l - 1];
      for (k = 0; k < nn_oo[l - 1]; k++)
        {
          ee_oo[l][k] = ee_oo[l - 1][k];
          id_oo[l][k] = id_oo[l - 1][k];
        };
      l--;
    };

  /* store current as first in back--mix list */

  nn_oo[0] = nn;
  for (i = 0; i < nn; i++)
    {
      id_oo[0][i] = id[i];
      ee_oo[0][i] = ee[i];
    };


  /* done */

  if (Pars.CurEvNo <= Pars.NumToPrint)
    printf ("exit bin_angcor_DGS\n");

  return (0);

}
