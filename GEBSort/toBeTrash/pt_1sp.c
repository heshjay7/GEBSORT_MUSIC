/* peak_total_1sp.c */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define DEBUG 0


/* ------------------------------------------------- */
int
f_peak (float sp[], float sps[], int cut, int lo, int hi, float *peak, float *area, float *sig, float *skew)
{
  /* declarations */

  int j, maxj, ww, jump;
  float max;
  float b1, b2, b;
  int i1, i2;
  double u2, u3, f, d1;

  /* initialize */

  *peak=-1;
  *area=-1;
  *sig=-1;
  *skew=-1;

  /* find new trigger level in local area */

  if (DEBUG)
    printf ("search in area: %i to %i\n", lo, hi);
  max = 0;
  for (j = lo; j <= hi; j++)
    {
      if (sps[j] > max)
        {
          max = sps[j];
          maxj = j;
          if (DEBUG && 0)
            printf ("%i : %f %f%i\n", j, sp[j], max);
        };
    };
  if (DEBUG)
    printf ("counts in strongest is %f at ch %i \n", max, maxj);
  max = max / 4;
  if (DEBUG)
    printf (" --> peak trigger level: %f\n", max);
  fflush (stdout);
  if (max < 30)
    {
      printf ("no peaks... max=%f\n", max);
      return (-1);
    };

  /* redetermine lo and hi with local max */

  jump=(hi-lo)/4;
  j = hi + 2;
  while (sps[j] < max && j > 0)
    j--;
  hi = j;
  j-=jump;
  while (sps[j] > max && j > 0)
    j--;
  lo = j;
  ww = hi - lo;
  if (DEBUG)
    printf ("redefined lo, hi: from %4i to %4i (width=%i)\n", lo, hi, ww);

  /* set ROI */

  hi += ww / 2;
  lo -= ww / 2;
  if (DEBUG)
    printf ("high peak search area: from %4i to %4i (width=%i)\n", lo, hi, ww);
  if (lo < 200 || lo > 15000)
    {
      printf ("high peak search area: from %4i to %4i; ", lo, hi);
      printf ("unable to find high search region\n");
      return (-1);
    };
  if (hi < 200 || hi > 15000)
    {
      printf ("high peak search area: from %4i to %4i; ", lo, hi);
      printf ("unable to find high search region\n");
      return (-1);
    };

  /* find background parameters */

  b1 = 0;
  for (j = (lo - 3); j <= (lo); j++)
    b1 += (float) sp[j];
  b1 /= 3.0;
  i1 = lo - 2;

  b2 = 0;
  for (j = (hi); j <= (hi + 3); j++)
    b2 += (float) sp[j];
  b2 /= 3.0;
  i2 = hi + 2;

  if (DEBUG)
    printf ("b1=%f b2=%f i1=%i,i2=%i\n", b1, b2, i1, i2);

  /* find peak position */

  *peak = 0;
  *area = 0;
  for (j = lo; j <= hi; j++)
    {
      b = (b2 - b1) * (j - i1) / (i2 - i1) + b1;
      if (DEBUG &&0)
        printf ("j=%i, sp[j]=%f, b=%f sp[j]-b=%f\n", j, sp[j], b, sp[j] - b);
      *peak += ((float) sp[j] - b) * (float) j;
      *area += ((float) sp[j] - b);
    };
  *peak /= *area;
  if (DEBUG)
    printf ("peak pos: %f\n", *peak);

  /* find the second and third moments */

  u2 = 0;
  u3 = 0;
  for (j = lo; j <= hi; j++)
    {
      b = (b2 - b1) * (j - i1) / (i2 - i1) + b1;
      f = sp[j] - b;
      d1 = (double) j - *peak;
      u2 += f * d1 * d1;
      u3 += f * d1 * d1 * d1;
    }
  u2 /= *area;
  u3 /= *area;
  *sig = sqrt (u2);
  *skew = u3 / (*sig) / (*sig) / (*sig);

  /* done */

  return (0);

};

/* ------------------------------------------------- */

int
pt_1sp (int dim, float sp[], int cut, float *peak1, float *peak2, float *area1, float *area2, float *ptval,
        float *sigma1, float *sigma2, float *skew1, float *skew2, float *sumarea)
{


  /* declarations */

  float max;
  int j, ipos, maxj, ww1, ww2, ww, st1, st2;
  int h1, l1, h2, l2;
  int i1, i2, w, lo, hi;
  float b1, b2, b, area;
  float sps[20000];

  if (DEBUG)
    printf ("\npt_1sp called with spectrum dim= %i\n", dim);

  /* make a smoothed spectrum for limits */

  for (j=0;j<2000; j++)
    sps[j]=0;
  for (j=2;j<(dim-2); j++)
    {
    sps[j]=0.5*sp[j];
    sps[j]+=0.25*sp[j-1];
    sps[j]+=0.25*sp[j+1];
//    printf("%f --> %f\n",sp[j], sps[j]);
    }
  

  /***********************************************/
  /* find the trigger levels and search channels */
  /***********************************************/

  /* find overall trigger level */

  max=-2000000;
  for (j = cut; j < dim; j++)
    {
      if (sps[j] > max)
        {
          max = sps[j];
          maxj = j;
        };
    };
  if (DEBUG)
    printf ("counts in strongest peak: %f at ch %i", max, maxj);
  max = max / 4;
  if (DEBUG)
    printf (" --> peak trigger level: %f\n", max);
  fflush (stdout);
//if(1)exit(0);
  if (max < 50)
    {
      printf ("no peaks... max=%f\n", max);
      return (-1);
    };

  /* find ROIs */

  while (sps[j] < max && j > 0)
    j--;
  h2 = j;
  while (sps[j] > max && j > 0)
    j--;
  l2 = j;
  ww2 = h2 - l2;
  h2 += ww2;
  l2 -= ww2;

  st2 = f_peak (sp, sps, cut, l2, h2, peak2, area2, sigma2, skew2);
//if(1)exit(0);
  j = l2 - 3 * ww2;

  while (sps[j] < max && j > 0)
    j--;
  h1 = j;
  while (sps[j] > max && j > 0)
    j--;
  l1 = j;
  ww1 = h1 - l1;
  h1 += ww1;
  l1 -= ww1;

  st1 = f_peak (sp, sps, cut, l1, h1, peak1, area1, sigma1, skew1);
  if (DEBUG)
    printf ("peak1 peak2 pos: %9.2, sigma2, skew2f %9.2f\n", *peak1, *peak2);

//if(1)exit(0);


  area = 0;
  for (j = cut; j < dim; j++)
    area += sp[j];
  *sumarea = area;

  /* find peak to total */

  *ptval = (*area1 + *area2) / area;
  /* printf("p/t: %9.4f\n",*ptval); */

  /* done */

  return (st1 + st2);

}
