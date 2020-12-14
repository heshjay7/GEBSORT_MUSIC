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

  lowch=atoi(argv[3]);
  desgain=atof(argv[4]);
  printf("desired gain: %f\n",desgain );

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
      return(1);
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

  for (i = 1; i <= 110; i++)
    {

      /* get spectrum */

      sprintf (str, "ehi%3.3i.spe", i);
      printf ("-- %s ", str);
      dim = DIM;
      rd_spe (str, &dim, sp);
      printf ("dim= %i ", dim);

      /* find counts */

      sum = 0;
      for (j = 10; j < (dim - 10); j++)
        sum += sp[j];
      printf ("sum= %8.0f ", sum);

      if (sum < 1000)
        {
          /* write dummy */

          off = 0;
          gain = 1.0;
          fprintf (fp, "%3i %8.2f %8.4f\n", i, off, gain);
        }
      else
        {
          nn++;

          /* find max channel */

          max = 0;
          for (j = lowch; j < (dim - 10); j++)
            if (sp[j] > max)
              max = sp[j];
          printf ("max %8.0f ", max);

          if (max > 100.0)
            {

              /* set threshold */

              thr = max / 4.0;
              printf ("thr %8.0f ", thr);

              /* find approx hi  peak pos */

              ihi = dim - 100;
              while (sp[ihi] < thr && ihi > 0)
                ihi--;
              printf("ihi=%i  ", ihi);

              /* find max hi peak positions */

              max = 0;
              maxch = 0;
              for (j = ihi - 10; j < ihi + 10; j++)
                if (sp[j] > max)
                  {
                    max = sp[j];
                    maxch = j;
                  };
              ihi = maxch;

              /* find mean hi peak position */

              mean = 0;
              ww = 0;
              for (j = ihi - w1; j < ihi + w1; j++)
                {
                  mean += sp[j] * j;
                  ww += sp[j];
                };
              mean /= ww;
              aphi = mean;

              /* find approx lo  peak pos */
              /* (continue down...)       */

              ilo = ihi - 50;
              while (sp[ilo] < thr && ilo > 0)
                ilo--;

              /* find max hlo peak positions */

              max = 0;
              maxch = 0;
              for (j = ilo - 10; j < ilo + 10; j++)
                if (sp[j] > max)
                  {
                    max = sp[j];
                    maxch = j;
                  };
              ilo = maxch;

              /* find mean lo peak position */

              mean = 0;
              ww = 0;
              for (j = ilo - w1; j < ilo + w1; j++)
                {
                  mean += sp[j] * j;
                  ww += sp[j];
                };

              mean /= ww;
              aplo = mean;

              printf ("p1,p2: %7.2f, %7.2f ", aplo, aphi);

              /* find new calibration */

              gain = (dphi - dplo) / (aphi - aplo);
              off = dphi - gain * aphi;
              printf ("offset= %8.2f gain=%8.4f", off/desgain, gain/desgain);

              fprintf (fp, "%3i %8.2f %8.4f\n", i, off/desgain, gain/desgain);
              
              /* display cmd file */
          
              fprintf (fp1, "sp %s\nds\ncr\n", str);

            }
          else
            printf ("not enough counts, will not change cal... ");

          /* report and set offset/gain to 0,1 */

//          fprintf (fp, "%3i     0.0   %8.4f   0.00   1.000000\n", i, pz);

        };

      printf ("\n");

    };

  /* done */

  printf ("found %i cal values\n", nn);
  printf ("display: d_ecal.cmd\n");

  exit (0);

}
