#include <stdlib.h>
#include <fcntl.h>
#include <stddef.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/param.h>
#include <sys/resource.h>
#include <sys/times.h>
#include <assert.h>
#include "GTMerge.h"

#define DIM 16384

#define BI207LO  569.702
#define BI207HI 1063.662
#define Y88LO    898.0450
#define Y88HI   1836.0630
#define CO60LO  1173.228
#define CO60HI  1332.490

/*-----------------------------------------------------------------------------*/

int
main (int argc, char **argv)
{
  /* declarations */

  int i, j, dim, mm, kk, nn = 0;
  char str[128], *p;
  float sp[DIM], sum, mean, spz = 0;
  FILE *fp, *fp1;
  float d1, pz, r1, r2, dplo, dphi, kevch = 1.0;
  float max, thr, aphi, aplo;
  float meanrawhi = 0, meanrawlo = 0, ww, w1 = 5;
  int maxch, i1, i2, ihi, ilo, lowch;
  float off, gain, desgain;
  int n, st;
  float area1[NGE], area2[NGE], ptval[NGE], peak1[NGE], peak2[NGE], sumarea[NGE];
  int deton[NGE], det[NGE];
  float sigma1[NGE], sigma2[NGE], skew1[NGE], skew2[NGE];

  /* prototypes */

  int
    pt_1sp (int, float sp[], int, float *, float *, float *, float *, float *, float *, float *, float *, float *,
            float *);
  int rd_spe (char *fn, int *dim, float *sp);

  /* help ? */

  if (argc != 5)
    {
      printf ("use: dgs_ecal   cal_file_name   source lowch  desgain\n");
      //               0           1            2       3       4
      printf ("will work on ehinnn.spe files, nnn=001...110\n");
      printf ("use \".x get_ecln.cc\" to extract those spe files \nn");
      exit (0);
    }

  /* get pars */

  lowch = atoi (argv[3]);
  desgain = atof (argv[4]);
  printf ("desired gain: %f\n", desgain);
  fflush(stdout);

  /* find desired peak positions for selected source */

  if ((p = strstr (argv[2], "207Bi")) != NULL)
    {
      r1 = BI207LO;
      r2 = BI207HI;
      printf ("assuming 207Bi source calibration\n");
      printf ("__with peaks at %7.3f and %7.3f keV!!!!!\n", r1, r2);
      fflush (stdout);

      printf ("will calibrate for %9.6f kev/channel\n", (float) kevch);
      dplo = BI207LO;
      dphi = BI207HI;
      printf ("==>desired peak positions: %9.3f and %9.3f\n", (float) dplo, (float) dphi);
      fflush (stdout);
    }
  else if ((p = strstr (argv[2], "88Y")) != NULL)
    {
      r1 = Y88LO;
      r2 = Y88HI;
      printf ("assuming 88Y source calibration\n");
      printf ("__with peaks at %7.3f and %7.3f keV!!!!!\n", r1, r2);
      fflush (stdout);

      printf ("will calibrate for %9.6f kev/channel\n", (float) kevch);
      dplo = Y88LO;
      dphi = Y88HI;
      printf ("==>desired peak positions: %9.3f and %9.3f\n", (float) dplo, (float) dphi);
      fflush (stdout);
    }
  else if ((p = strstr (argv[2], "60Co")) != NULL)
    {
      r1 = CO60LO;
      r2 = CO60HI;
      printf ("assuming 60Co source calibration\n");
      printf ("__with peaks at %7.3f and %7.3f keV!!!!!\n", r1, r2);
      fflush (stdout);

      printf ("will calibrate for %9.6f kev/channel\n", (float) kevch);
      dplo = CO60LO;
      dphi = CO60HI;
      printf ("==>desired peak positions: %9.3f and %9.3f\n", (float) dplo, (float) dphi);
      fflush (stdout);
    }
  else
    {
      printf ("***error, source \"%s\" not recognized\n", argv[2]);
      printf ("   must be one of: 207Bi, 88Y or 60Co\n");
      printf ("\n");
      return (1);
    };


  fp = fopen ((char *) argv[1], "w");
  if (fp != NULL)
    printf ("%s is open for write\n", argv[1]);
  else
    {
      printf ("could not open %s for write\n", argv[1]);
      exit (0);
    };

  fp1 = fopen ("d_ecal.cmd", "w");
  if (fp1 != NULL)
    printf ("%s is open for write\n", "d_ecal.cmd");
  else
    {
      printf ("could not open %s for write\n", "d_ecal.cmd");
      exit (0);
    };

  n = 0;
  for (i = 1; i <= 110; i++)
    {

      /* get spectrum */

      sprintf (str, "ehi%3.3i.spe", i);
      printf ("-- %s ", str);
  fflush(stdout);
      dim = DIM;
      rd_spe (str, &dim, sp);
//      printf ("dim= %i ", dim);

      /* find ROIs */

      /* pt analyze */

      det[n] = i;
      st = pt_1sp (dim, sp, lowch,
                   &peak1[n], &peak2[n],
                   &area1[n], &area2[n], &ptval[n], &sigma1[n], &sigma2[n], &skew1[n], &skew2[n], &sumarea[n]);

      if (st == 0)
        {

          /* process */

//                ave += ptval[n];

//                printf ("[%3i] -> ", i);
              printf ("peaks: %6.1f %6.1f  ", peak1[n] / kevch, peak2[n] / kevch);
//                printf ("FWHM: %5.2f  %5.2f  ", 2.355 * sigma1[n] / kevch, 2.355 * sigma2[n] / kevch);
//                printf ("areas: %8.0f %8.0f  ", area1[n], area2[n]);
//                printf ("skew: %5.2f, %5.2f", skew1[n], skew2[n]);

//                fwhm_ave += sigma2[n];
//                printf (" ok\n");

              /* find new calibration */

              gain = (dphi - dplo) / (peak2[n] - peak1[n]);
              off = dphi - gain * peak2[n];
              printf ("offset= %8.2f gain=%8.4f", off / desgain, gain / desgain);

              fprintf (fp, "%3i %8.2f %8.4f\n", i, off / desgain, gain / desgain);
              fprintf (fp1, "sp %s\nds\ncr\n", str);     
              printf ("\n");

              n++;
        };

    };

  /* done */

  printf ("found %i cal values\n", n);
  printf ("display: d_ecal.cmd\n");

  exit (0);

}

