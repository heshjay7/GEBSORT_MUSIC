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

#define DIM 4096
#define WIDTH 50

/*-----------------------------------------------------------------------------*/

int
main (int argc, char **argv)
{
  /* declarations */

  int i, j, dim, mm, kk, nn = 0;
  int i1, i2, jmax;
  char str[128];
  float sp[DIM], spm[DIM], sum, mean, spz = 0, max;
  FILE *fp, *fp1;
  float power_factor, pz, SZfactor;

  /* help ? */

  if (argc != 5)
    {
      printf ("use: dgs_pz M K cal_file_name factor\n");
      printf ("will work on pznnn.spe files, nnn=001...110\n");
      printf ("use \".x get_pz.cc\" to extract those spe files \nn");
      printf ("\n");
      printf ("with factor=1.0, you get the PZ values that Shoufei Zhu\n");
      printf ("determines. The factor is multiplied to the PZ values\n");
      printf ("before written out to the cal_file_name\n");
      printf ("\n");
      exit (0);
    }

  /* get pars */

  mm = atoi (argv[1]);
  kk = atoi (argv[2]);
  power_factor = ((float) mm + (float) kk) / (float) mm;
  SZfactor= atof(argv[4]);

  fp = fopen ((char *) argv[3], "w");
  if (fp != NULL)
    printf ("%s is open for write\n", argv[3]);
  else
    {
      printf ("could not open %s for write\n", argv[3]);
      exit (0);
    };

  fp1 = fopen ("d_pz.cmd", "w");
  if (fp1 != NULL)
    printf ("%s is open for write\n", "d_pz.cmd");
  else
    {
      printf ("could not open %s for write\n", "d_pz.cmd");
      exit (0);
    };

  for (i = 1; i <= 110; i++)
    {

      /* get spectrum */

      sprintf (str, "pz%3.3i.spe", i);
      printf ("-- %s ", str);
      dim = DIM;
      rd_spe (str, &dim, sp);
      printf ("dim= %i ", dim);

      /* make a smoothed version of the spectrum */

      
      for (j = 1; j < dim; j++)
        spm[j]=0;
      for (j = 2; j < (dim-1); j++)
        spm[j]=0.5*sp[j]+0.25*sp[j-1]+0.25*sp[j+1];

      /* inspect smoothed spectrum */

      sum = 0;
      jmax=-1;
      max=-1;
      for (j = 850; j < 950; j++)
        {
        sum += spm[j];
        if (sp[j]>max)
          {
          max=spm[j];
          jmax=j;
          };
        };
      printf ("sum= %8.0f, jmax= %3i, ", sum, jmax);

      if (sum < 1000)
        {
          /* write dummy */

          fprintf (fp, "%3i     1.0 \n", i);
        }
      else
        {
          nn++;

          /* find mean value == pol zero */

          i1=jmax-WIDTH;
          i2=jmax+WIDTH;
          mean = 0;
          sum = 0;
          for (j = i1; j < i2; j++)
            {
            mean += (sp[j] * j);
            sum  += sp[j];
           };
          mean /= sum;
          mean *= (2.0 / 2000);
          if (mean > 1.0)
            mean = 1.0;
          printf ("mean= %8.4f ", mean);
          pz = powf (mean, power_factor);
          pz *= SZfactor;
          printf ("pz= %8.4f ", pz);
          spz += pz;

          /* report  */

          fprintf (fp, "%3i   %8.4f \n", i, pz);
          
          /* display cmd file */
          
          fprintf (fp1, "sp %s\nds\ncr\n", str);

        };


      printf ("\n");

    };

  /* done */
  
  fclose(fp);
  fclose(fp1);

  printf ("M=%i and K=%i, power_factor=%f\n", mm, kk, power_factor);
  printf ("found %i pz values, mean pz is %f\n", nn, spz / nn);
  printf ("display: d_pz.cmd\n");

  exit (0);

}
