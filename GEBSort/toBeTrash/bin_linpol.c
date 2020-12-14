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
#include "GEBSort.h"
#include "GTMerge.h"


unsigned int seed;

#define NBINS 180
#define LOQ 1000

/* pointers to ROOT spectra */

TH2F *lp_azin;
TH2F *lp_azin_ref;

/* parameters */

extern PARS Pars;
extern EXCHANGE exchange;

unsigned int linpol_good_update;
unsigned int linpol_good_update_ref;
unsigned int linpol_good_emission_ang;
unsigned int linpol_total;
unsigned int linpol_good_basic;
unsigned int linpol_good_rr_scat;
unsigned int linpol_good_scatang;
int qlen;

/*-----------------------------------------------------*/

int
exit_linpol ()
{
  printf ("\n");
  printf ("bin_linpol statistics\n");
  printf ("\n");
  printf ("LOQ= %i (mixing depth)\n", LOQ);
  printf ("\n");
  printf ("of linpol_total\n");
  printf ("linpol_total             %10i, %6.2f%%\n", linpol_total, 100.0 * linpol_total / linpol_total);
  printf ("linpol_good_basic        %10i, %6.2f%%\n", linpol_good_basic, 100.0 * linpol_good_basic / linpol_total);
  printf ("__means: two gamma rays\n");
  printf ("  tracked, fom in range, ndet limits OK, ndet>=2\n");
  printf ("  and the 'beam axis' gamma ray was found\n");
  printf ("\n");

  printf ("of linpol_good_basic count\n");
  printf ("linpol_good_emission_ang %10i, %6.2f%%\n", linpol_good_emission_ang,
          100.0 * linpol_good_emission_ang / linpol_good_basic);
  printf ("linpol_good_rr_scat      %10i, %6.2f%%\n", linpol_good_rr_scat,
          100.0 * linpol_good_rr_scat / linpol_good_basic);
  printf ("linpol_good_scatang      %10i, %6.2f%%\n", linpol_good_scatang,
          100.0 * linpol_good_scatang / linpol_good_basic);
  printf ("\n");

  printf ("of linpol_good_basic count\n");
  printf ("linpol_good_update       %10i, %6.2f%%\n", linpol_good_update,
          100.0 * linpol_good_update / linpol_good_basic);
  printf ("of linpol_total\n");
  printf ("linpol_good_update       %10i, %6.2f%%\n", linpol_good_update, 100.0 * linpol_good_update / linpol_total);
  printf ("\n");

  printf ("reference contribution to error: %6.2f%%\n",
          100.0 * sqrt ((double) linpol_good_update / (double) linpol_good_update_ref));
  printf ("\n");

  /* done */

  return (0);

}

/*-----------------------------------------------------*/

int
sup_linpol ()
{
  /* declarations */

  TH1D *mkTH1D (char *, char *, int, double, double);
  TH2F *mkTH2F (char *, char *, int, double, double, int, double, double);

  lp_azin = mkTH2F ((char *) "lp_azin", (char *) "lp_azin", Pars.GGMAX, 1, Pars.GGMAX, 2 * NBINS + 1, -NBINS, NBINS);
  lp_azin->SetXTitle ((char *) "energy");
  lp_azin->SetYTitle ((char *) "linear polarization: azimuth angle of normal");

  lp_azin_ref =
    mkTH2F ((char *) "lp_azin_ref", (char *) "lp_azin_ref", Pars.GGMAX, 1, Pars.GGMAX, 2 * NBINS + 1, -NBINS, NBINS);
  lp_azin_ref->SetXTitle ((char *) "energy");
  lp_azin_ref->SetYTitle ((char *) "linear polarization: azimuth angle of normal REFERENCE");

  linpol_good_update = 0;
  linpol_good_update_ref = 0;
  linpol_good_emission_ang = 0;
  linpol_total = 0;
  linpol_good_rr_scat = 0;
  linpol_good_scatang = 0;
  linpol_good_basic = 0;

  qlen = 1;

  /* list what we have */

  printf (" we have define the following spectra:\n");

  Pars.wlist = gDirectory->GetList ();
  Pars.wlist->Print ();


  return (0);

};

/* ----------------------------------------------------------------- */

int
find_linpol (double n0_x, double n0_y, double n0_z,
             double n1_x, double n1_y, double n1_z, double n2_x, double n2_y, double n2_z, float ee, TH2F * lp_azin)
{

  /* declarations */

  static int nn = 0;
  double d1, d2;

  double xx, yy, zz;
  double pol_new, azi_new;

  /* normal vector to reaction and scatter plane */

  double x1, y1, z1;
  double x2, y2, z2;

  /* new coordinate system */

  double newXaxis_x, newXaxis_y, newXaxis_z;
  double newYaxis_x, newYaxis_y, newYaxis_z;
  double newZaxis_x, newZaxis_y, newZaxis_z;

  /* prototypes */

  int crossprod (double, double, double, double, double, double, double *, double *, double *);
  double unitvector (double *, double *, double *);
  int check_coord_sys (double, double, double, double, double, double, double, double, double);
  int coord_in_new (double, double, double,
                    double, double, double,
                    double, double, double, double, double, double, double *, double *, double *);
  double AzimuthFromCartesian (double, double, double);
  double rad2deg (double);
  int check_unitvector (double, double, double);

  if (Pars.CurEvNo <= Pars.NumToPrint)
    {
      printf ("find_linpol @ event # %i:\n", Pars.CurEvNo);
      printf ("  n0  (%10.7f %10.7f %10.7f) \n", n0_x, n0_y, n0_z);
      printf ("  n1  (%10.7f %10.7f %10.7f) \n", n1_x, n1_y, n1_z);
      printf ("  n2  (%10.7f %10.7f %10.7f) \n", n2_x, n2_y, n2_z);
    };


  /* Normal vectors for reaction/scatter planes */
  /* Remember to make unit vector; because a cross */
  /* does not make a unit vector from unit vectors */

  crossprod (n0_x, n0_y, n0_z, n1_x, n1_y, n1_z, &x1, &y1, &z1);
  d1 = unitvector (&x1, &y1, &z1);
  crossprod (n1_x, n1_y, n1_z, n2_x, n2_y, n2_z, &x2, &y2, &z2);
  d1 = unitvector (&x2, &y2, &z2);

  /* new cordsys x-axis vector is a vector in the reaction plane */
  /* perpendicular to the gamma emission direction. */
  /* So we make a cross product (remember to make unit vector) */

  crossprod (x1, y1, z1, n1_x, n1_y, n1_z, &newXaxis_x, &newXaxis_y, &newXaxis_z);
  d1 = unitvector (&newXaxis_x, &newXaxis_y, &newXaxis_z);

  /* new cordsys z-axis gamma emission direction */

  newZaxis_x = n1_x;
  newZaxis_y = n1_y;
  newZaxis_z = n1_z;

  /* new cordsys y-axis is reaction plane normal vector */

  newYaxis_x = x1;
  newYaxis_y = y1;
  newYaxis_z = z1;

  if (Pars.CurEvNo <= Pars.NumToPrint)
    {
      printf ("new coord system:\n");
      printf ("  x-axis: (%10.7f %10.7f %10.7f)\n", newXaxis_x, newXaxis_y, newXaxis_z);
      printf ("  y-axis: (%10.7f %10.7f %10.7f)\n", newYaxis_x, newYaxis_y, newYaxis_z);
      printf ("  z-axis: (%10.7f %10.7f %10.7f)\n", newZaxis_x, newZaxis_y, newZaxis_z);
      check_coord_sys (newXaxis_x, newXaxis_y, newXaxis_z, newYaxis_x, newYaxis_y,
                       newYaxis_z, newZaxis_x, newZaxis_y, newZaxis_z);
    };

  coord_in_new (newXaxis_x, newXaxis_y, newXaxis_z,
                newYaxis_x, newYaxis_y, newYaxis_z, newZaxis_x, newZaxis_y, newZaxis_z, x2, y2, z2, &xx, &yy, &zz);

  if (Pars.CurEvNo <= Pars.NumToPrint)
    {
      printf ("normal vector to scattering plane in new system:\n");
      printf ("  (%10.7f %10.7f %10.7f)\n", xx, yy, zz);
      check_unitvector (xx, yy, zz);
    };

  /* polar angle in new coordinate system is */
  /* [trivial since we found zz above] */
  /* [also trivial because it is zero */

  pol_new = acos (zz);
//                                    lp_scatn->Fill (grh->gr[j].esum, rad2deg (pol_new), 1);

  /* azimuth angle, in new coordinate system, is determined by xx and yy */
  /* CRC page 385 + trigfun.c */

  azi_new = AzimuthFromCartesian (xx, yy, zz);
  lp_azin->Fill (ee, rad2deg (azi_new), 1);
  if (Pars.CurEvNo <= Pars.NumToPrint)
    {
      printf ("for scatter plane normal vector:\n");
      printf ("pol_new = %5.2f (deg)\n", rad2deg (pol_new));
      printf ("azi_new = %5.2f (deg)\n", rad2deg (azi_new));
    };

#if(0)
  nn++;
  if (nn > 100)
    exit (0);
#endif

  /* done */

  return (0);

};

/* ----------------------------------------------------------------- */

int
bin_linpol (GEB_EVENT * GEB_event)
{

  /* declarations */

  int i, j, k, m, ok[MAX_GAMMA_RAYS];
  TRACKED_GAMMA_HIT *grh;
  double rr, rr_scat[MAX_GAMMA_RAYS];
  int nv, beamDir, nitgr;
  double xx_0[MAX_GAMMA_RAYS], yy_0[MAX_GAMMA_RAYS], zz_0[MAX_GAMMA_RAYS];
  double xx_1[MAX_GAMMA_RAYS], yy_1[MAX_GAMMA_RAYS], zz_1[MAX_GAMMA_RAYS];
  float ee[MAX_GAMMA_RAYS], see[MAX_GAMMA_RAYS], emissionAng[MAX_GAMMA_RAYS];
  float semissionAng[MAX_GAMMA_RAYS], secondScatterAngg[MAX_GAMMA_RAYS];
  double sn0_x[MAX_GAMMA_RAYS], sn0_y[MAX_GAMMA_RAYS], sn0_z[MAX_GAMMA_RAYS];
  double sn1_x[MAX_GAMMA_RAYS], sn1_y[MAX_GAMMA_RAYS], sn1_z[MAX_GAMMA_RAYS];
  double sn2_x[MAX_GAMMA_RAYS], sn2_y[MAX_GAMMA_RAYS], sn2_z[MAX_GAMMA_RAYS];

  /* the list of old good events we will mix with for the reference */
  /* index [0] is the current */

  static double n0_x[MAX_GAMMA_RAYS], n0_y[MAX_GAMMA_RAYS], n0_z[MAX_GAMMA_RAYS];
  static double n1_x[MAX_GAMMA_RAYS], n1_y[MAX_GAMMA_RAYS], n1_z[MAX_GAMMA_RAYS];
  static double n2_x[MAX_GAMMA_RAYS], n2_y[MAX_GAMMA_RAYS], n2_z[MAX_GAMMA_RAYS];

  /* prototypes */

  double unitvector (double *, double *, double *);
  double dotproductangle (double, double, double, double, double, double);
  double rad2deg (double);
  int check_unitvector (double, double, double);
  int ranVectorOnSphere (double *, double *, double *);

  if (Pars.CurEvNo <= Pars.NumToPrint)
    printf ("entered bin_linpol\n");

  /* run through event and collect */
  /* events with two or more interactions */
  /* with good FOM */

  nitgr = 0;
  for (i = 0; i < GEB_event->mult; i++)
    if (GEB_event->ptgd[i]->type == GEB_TYPE_TRACK)
      {
        grh = (TRACKED_GAMMA_HIT *) GEB_event->ptinp[i];
        linpol_total += grh->ngam;
        for (j = 0; j < grh->ngam; j++)
          if (grh->gr[j].ndet >= Pars.linpol_nlo && grh->gr[j].ndet <= Pars.linpol_nhi)
            if (grh->gr[j].fom >= Pars.fomlo[grh->gr[j].ndet])
              if (grh->gr[j].fom <= Pars.fomhi[grh->gr[j].ndet])
                if (grh->gr[j].tracked)
                  {
                    if (Pars.CurEvNo <= Pars.NumToPrint)
                      {
                        printf ("[nitgr=%i] interesting event\n", nitgr);
                        printf ("gam %i, esum= %6.1f keV ndet= %i::\n", grh->ngam, grh->gr[j].esum, grh->gr[j].ndet);
                        printf ("__e0=%6.1f  (%10.7f %10.7f %10.7f) \n", grh->gr[j].e0, grh->gr[j].x0, grh->gr[j].y0,
                                grh->gr[j].z0);
                        printf ("__e1=%6.1f  (%10.7f %10.7f %10.7f) \n", grh->gr[j].e1, grh->gr[j].x1, grh->gr[j].y1,
                                grh->gr[j].z1);
                      };

                    /* if ndet>=2 it is a potentially good event. */
                    /* But also if the gamma ray that polarizes is */
                    /* fulfilled, then that too is a good event even if ndet=1 */

                    if (grh->gr[j].ndet > 1 || (fabs (ee[j] - Pars.linpol_source_ee) < Pars.linpol_source_de))
                      {
                        /* store event */

                        ee[nitgr] = grh->gr[j].esum;
                        xx_0[nitgr] = grh->gr[j].x0;
                        yy_0[nitgr] = grh->gr[j].y0;
                        zz_0[nitgr] = grh->gr[j].z0;
                        xx_1[nitgr] = grh->gr[j].x1;
                        yy_1[nitgr] = grh->gr[j].y1;
                        zz_1[nitgr] = grh->gr[j].z1;

                        nitgr++;

                      };
                  };
      };

  if (Pars.CurEvNo <= Pars.NumToPrint)
    printf ("nitgr=%i\n", nitgr);

  /* drop out if there is no hope */

  if (Pars.linpol_usesource == 1 && nitgr < 2)
    return (0);
  if (Pars.linpol_usesource == 0 && nitgr < 1)
    return (0);

  if (Pars.CurEvNo <= Pars.NumToPrint)
    printf ("passed nitgr test\n");


  /* list event */

  if (Pars.CurEvNo <= Pars.NumToPrint)
    {
      for (j = 0; j < nitgr; j++)
        {
          printf ("  gam %i, esum= %6.1f keV::\n", j, ee[j]);
          printf ("    (%10.7f %10.7f %10.7f) \n", xx_0[j], yy_0[j], zz_0[j]);
          printf ("    (%10.7f %10.7f %10.7f) \n", xx_1[j], yy_1[j], zz_1[j]);
        };
    };

  /* search for polarization gamma ray in event */

  nv = 0;
  if (Pars.linpol_usesource == 1)
    {
      for (j = 0; j < nitgr; j++)
        {
          if (fabs (ee[j] - Pars.linpol_source_ee) < Pars.linpol_source_de)
            {
              nv++;
              beamDir = j;
            };
        };
    }
  else
    {
      /* fool the rest of the code */

      nv = 1;
      beamDir = -1;
    };

  /* drop out if no 'first' gamma ray or beam axis */

  if (Pars.CurEvNo <= Pars.NumToPrint)
    printf ("nv=%i, beamDir=%i\n", nv, beamDir);
  if (nv != 1)
    return (0);
  linpol_good_basic += (nitgr - 1);

  if (Pars.CurEvNo <= Pars.NumToPrint)
    {
      for (j = 0; j < nitgr; j++)
        {
          printf ("  gam %i, esum= %6.1f keV::\n", j, ee[j]);
          printf ("    (%10.7f %10.7f %10.7f) \n", xx_0[j], yy_0[j], zz_0[j]);
          printf ("    (%10.7f %10.7f %10.7f) \n", xx_1[j], yy_1[j], zz_1[j]);
        };
    };

  /* make unit vectors */

  for (j = 0; j < nitgr; j++)
    if (j == beamDir)
      {
        /* this is the 'beam direction' */
        /* there is only one */

        n0_x[0] = xx_0[j] - Pars.target_x;
        n0_y[0] = yy_0[j] - Pars.target_y;
        n0_z[0] = zz_0[j] - Pars.target_z;
        rr = unitvector (&n0_x[0], &n0_y[0], &n0_z[0]);
        xx_1[j] = yy_1[j] = zz_1[j] = 0;
      }
    else
      {
        /* first interaction point unit vector */

        n1_x[j] = xx_0[j] - Pars.target_x;
        n1_y[j] = yy_0[j] - Pars.target_y;
        n1_z[j] = zz_0[j] - Pars.target_z;
        rr = unitvector (&n1_x[j], &n1_y[j], &n1_z[j]);

        /* second interaction point unit vector */

        n2_x[j] = xx_1[j] - xx_0[j];
        n2_y[j] = yy_1[j] - yy_0[j];
        n2_z[j] = zz_1[j] - zz_0[j];
        rr_scat[j] = unitvector (&n2_x[j], &n2_y[j], &n2_z[j]);

      };

  /* if we use beam axis (already normalized) */

  if (Pars.linpol_usesource == 0)
    {
      n0_x[0] = Pars.beamdir[0];
      n0_y[0] = Pars.beamdir[1];
      n0_z[0] = Pars.beamdir[2];
      if (Pars.CurEvNo <= Pars.NumToPrint)
        printf ("use beam dir:  (%10.7f %10.7f %10.7f) \n", n0_x[0], n0_y[0], n0_z[0]);
    };

  /* find emission direction - */
  /* We have to do after loop above  */
  /* since the beam axis might  */
  /* not come first */

  for (j = 0; j < nitgr; j++)
    if (j != beamDir)
      {
        emissionAng[j] = dotproductangle (n0_x[0], n0_y[0], n0_z[0], n1_x[j], n1_y[j], n1_z[j]);
        secondScatterAngg[j] = dotproductangle (n1_x[j], n1_y[j], n1_z[j], n2_x[j], n2_y[j], n2_z[j]);
      };
  emissionAng[beamDir] = -1;
  secondScatterAngg[beamDir] = -1;

  /* filter event on rr_scat and secondScatterAng */
  /* to purify the event for further use */
  /* we do not filter on emissionAng because of the */
  /* mixing with old events later on */

  m = 0;
  for (j = 0; j < nitgr; j++)
    {
      ok[j] = 0;
      if (j != beamDir)
        {

          if ((float) rr_scat[j] >= Pars.linpol_rrlo)
            if ((float) rr_scat[j] <= Pars.linpol_rrhi)
              {
                linpol_good_rr_scat++;
                ok[j] = 1;
              };

          if ((float) secondScatterAngg[j] >= Pars.linpol_scatmin)
            if ((float) secondScatterAngg[j] <= Pars.linpol_scatmax)
              {
                linpol_good_scatang++;
                if (ok[j])
                  ok[j] = 1;
              };

          /* just for counting, no veto on emissionAng */

          if (emissionAng[j] >= Pars.linpol_polmin)
            if (emissionAng[j] <= Pars.linpol_polmax)
              linpol_good_emission_ang++;

        }
      else
        ok[j] = 1;
      if (ok[j])
        m++;
    };

  if (Pars.CurEvNo <= Pars.NumToPrint)
    printf ("#m=%i\n", m);

  /* m test */

  if (Pars.linpol_usesource == 1 && m < 2)
    return (0);
  if (Pars.linpol_usesource == 0 && m < 1)
    return (0);
  if (Pars.CurEvNo <= Pars.NumToPrint)
    printf ("passed m test, #m=%i\n", m);

  if (Pars.CurEvNo <= Pars.NumToPrint)
    {
      printf ("[4]event [nitgr=%i]\n", nitgr);
      for (j = 0; j < nitgr; j++)
        if (j == beamDir)
          {
            printf ("  n0  (%10.7f %10.7f %10.7f) beam dir e=%6.1f\n", n0_x[0], n0_y[0], n0_z[0], ee[j]);
            check_unitvector (n0_x[0], n0_y[0], n0_z[0]);
          }
        else
          {
            printf ("  gam %i, esum= %6.1f keV OK=%i\n", j, ee[j], rr_scat[j], ok[j]);
            printf ("  rr_scat %5.1f mm ", rr_scat[j]);
            printf ("emissionAng %5.1f deg ", rad2deg (emissionAng[j]));
            printf ("secondScatterAngg %5.1f deg\n", rad2deg (secondScatterAngg[j]));
            printf ("  n1  (%10.7f %10.7f %10.7f) \n", n1_x[j], n1_y[j], n1_z[j]);
            check_unitvector (n1_x[j], n1_y[j], n1_z[j]);
            printf ("  n2  (%10.7f %10.7f %10.7f) \n", n2_x[j], n2_y[j], n2_z[j]);
            check_unitvector (n2_x[j], n2_y[j], n2_z[j]);
          };
    };

  /* condence the event. */
  /* mostly for the storage for mixing. */
  /* We do this to make the stored events  */
  /* for the mixing as relevant as possible.  */
  /* The only thing we don't store on is the  */
  /* emission angle (which here means the ang  */
  /* wrt to the n0 axis) whis has to be  */
  /* calculated for each mixed event */

  for (j = 0; j < nitgr; j++)
    {
      semissionAng[j] = emissionAng[j];
      see[j] = ee[j];
      sn0_x[j] = n0_x[j];
      sn0_y[j] = n0_y[j];
      sn0_z[j] = n0_z[j];
      sn1_x[j] = n1_x[j];
      sn1_y[j] = n1_y[j];
      sn1_z[j] = n1_z[j];
      sn2_x[j] = n2_x[j];
      sn2_y[j] = n2_y[j];
      sn2_z[j] = n2_z[j];
    };

  m = 0;
  for (j = 0; j < nitgr; j++)
    if (ok[j])
      {
        emissionAng[j] = semissionAng[j];
        ee[m] = see[j];
        n0_x[m] = sn0_x[j];
        n0_y[m] = sn0_y[j];
        n0_z[m] = sn0_z[j];
        n1_x[m] = sn1_x[j];
        n1_y[m] = sn1_y[j];
        n1_z[m] = sn1_z[j];
        n2_x[m] = sn2_x[j];
        n2_y[m] = sn2_y[j];
        n2_z[m] = sn2_z[j];
        if (j == beamDir)
          beamDir = m;
        m++;
      };
  nitgr = m;

  if (Pars.CurEvNo <= Pars.NumToPrint)
    {
      printf ("[4]event [nitgr=%i]\n", nitgr);
      for (j = 0; j < nitgr; j++)
        if (j == beamDir)
          {
            printf ("  n0  (%10.7f %10.7f %10.7f) beam dir e=%6.1f\n", n0_x[0], n0_y[0], n0_z[0], ee[j]);
            check_unitvector (n0_x[0], n0_y[0], n0_z[0]);
          }
        else
          {
            printf ("  gam %i, esum= %6.1f keV\n", j, ee[j]);
            printf ("  n1  (%10.7f %10.7f %10.7f) \n", n1_x[j], n1_y[j], n1_z[j]);
            check_unitvector (n1_x[j], n1_y[j], n1_z[j]);
            printf ("  n2  (%10.7f %10.7f %10.7f) \n", n2_x[j], n2_y[j], n2_z[j]);
            check_unitvector (n2_x[j], n2_y[j], n2_z[j]);
          };
    };

  /* now do the linpol for the signal lp_azin */
  /* if conditions are OK */

  for (j = 0; j < nitgr; j++)
    if (j != beamDir)
      if (emissionAng[j] >= Pars.linpol_polmin)
        if (emissionAng[j] <= Pars.linpol_polmax)
          {
            find_linpol (n0_x[0], n0_y[0], n0_z[0], n1_x[j], n1_y[j], n1_z[j], n2_x[j], n2_y[j], n2_z[j], ee[j],
                         lp_azin);
            linpol_good_update++;
          };


  /* now do the linpol for the reference lp_azin_ref */
  /* if conditions are OK */
  /* everythin is the same, but we use */
  /* the 'beam' direction of the stored */
  /* events so the reaction plane becomes uncorrelated */

  /* make beam axis a random direction */

  if (Pars.linpol_random)
    ranVectorOnSphere (&n0_x[0], &n0_y[0], &n0_z[0]);


//  if (linpol_good_update > LOQ)
  for (m = 1; m < qlen; m++)
    for (j = 0; j < nitgr; j++)
      if (j != beamDir)
        {
          emissionAng[j] = dotproductangle (n0_x[m], n0_y[m], n0_z[m], n1_x[j], n1_y[j], n1_z[j]);

          if (emissionAng[j] >= Pars.linpol_polmin)
            if (emissionAng[j] <= Pars.linpol_polmax)
              {
                find_linpol (n0_x[m], n0_y[m], n0_z[m], n1_x[j], n1_y[j], n1_z[j], n2_x[j], n2_y[j], n2_z[j],
                             ee[j], lp_azin_ref);
                linpol_good_update_ref++;
              };
        };


  /* store LOQ back of beam directions */
  /* [0] entry will be used next time */

  if (qlen < LOQ)
    qlen++;

  j = qlen - 1;
  while (j > 0)
    {
      n0_x[j] = n0_x[j - 1];
      n0_y[j] = n0_y[j - 1];
      n0_z[j] = n0_z[j - 1];
      j--;
    };


/* done */

  if (Pars.CurEvNo <= Pars.NumToPrint)
    printf ("exit bin_linpol\n");

  return (0);
}
