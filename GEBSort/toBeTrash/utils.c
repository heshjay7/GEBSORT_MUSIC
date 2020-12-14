#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <math.h>
#include <string.h>
#include <assert.h>

/* ----------------------------------------------------------------- */

void
pprint_32 (char *str, unsigned int vv)
{
  int i, k;
  unsigned int val;

  printf ("%s", str);
  printf ("%10u %10i 0x%8.8x ", vv, (unsigned int) vv, vv);
  printf ("|");
  val = 0x80000000;
  for (k = 32; k > 0; k--)
    {
      if ((vv & val) == val)
        printf ("1");
      else
        printf ("0");
      if ((k % 4 - 1) == 0)
        printf ("|");
      val /= 2;
    };
//  printf ("|");
  printf ("\n");

}


/* ----------------------------------------------------------------- */

unsigned int
c32bit24bit (int vv)
{

  /* declarations */

  unsigned int val;
  int sign, base = 0x00800000;

  /* process the sign */

  if (vv < 0)
    {
      val = base;
      val += ((unsigned int) vv) & 0x007fffff;
    }
  else
    val = (unsigned int) vv;

  /* convert */

  if (sign == 0)
    val = vv;
  else
    {
      vv -= 0x00800000;
      val = -(base - vv);
    };

  /* done */

  return (val);

}


/* ----------------------------------------------------------------- */

int
c24bit32bit (unsigned int vv)
{

  /* declarations */

  int val, sign, base = 8388608;

  /* find the sign */

  sign = (vv & 0x00800000) >> 23;

  /* convert */

  if (sign == 0)
    val = vv;
  else
    {
      vv -= 0x00800000;
      val = -(base - vv);
    };

  /* done */

  return (val);

}

/* ----------------------------------------------------------------- */

int
twoscomp_to_int_24 (unsigned int tempE)
{

/* this function/method comes from Shofei Zhu */

  /*declarations */

  int rawE;
  int i, k;
  unsigned int val;


  /* convert */

  tempE = (tempE << 8);
  rawE = (int) tempE;
  rawE = rawE >> 8;



  /* done */

  return (rawE);

}

/* ----------------------------------------------------------------- */

int
test_convert ()
{

  /* declarations */

  int ii1, ii2, i;
  unsigned int ui1, ui2;
  int base = 8388608;

  printf ("\n");
  printf ("\n");
  printf ("\n");
  printf ("convert to 24 bit signed two's complement\n");
  printf ("\n");

  for (ii1 = 0; ii1 < 5; ii1++)
    {
      printf ("%15i converts to \n", ii1);
      pprint_32 ("__ ", c32bit24bit (ii1));
    };

  printf ("\n");

  for (ii1 = (base - 5); ii1 < base; ii1++)
    {
      printf ("%15i converts to \n", ii1);
      pprint_32 ("__ ", c32bit24bit (ii1));
    };
  printf ("\n");

  for (ii1 = 0; ii1 > -5; ii1--)
    {
      printf ("%15i converts to \n", ii1);
      pprint_32 ("__ ", c32bit24bit (ii1));
    };

  printf ("\n");

  for (ii1 = (-base + 5); ii1 >= -base; ii1--)
    {
      printf ("%15i converts to \n", ii1);
      pprint_32 ("__ ", c32bit24bit (ii1));
    };

  printf ("\n");
  /* done */

  printf ("\n");
  printf ("\n");
  exit (0);

};
