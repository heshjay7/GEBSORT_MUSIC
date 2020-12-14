
#include <stdlib.h>
#include <fcntl.h>
#include <stddef.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <limits.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/param.h>
#include <sys/resource.h>
#include <sys/times.h>
#include <assert.h>

#include "GTMerge.h"

int
main (int argc, char **argv)
{

  /* declarations */

  int nge, nbgo, nside, naux;
  int i, j, i1;


  /* DGS  */

  nge = 1;
  nbgo = 1;
  nside = 1;
  naux = 1;
  j = 0;
  for (i = 1010; i < (41 * 11 + 1010); i++)
    {

      printf ("%4i ", i);
      if (j < 5)
        {
          printf ("%4i %4i BGO\n", BGO, nbgo);
          nbgo++;
        }
      else if (j < 10)
        {
          printf ("%4i %4i GE\n", GE, nge);
          nge++;
        }
      else if (j < 15)
        {
          printf ("%4i %4i SIDE\n", SIDE, nside);
          nside++;
        }
      else if (j < 20)
        {
          printf ("%4i %4i AUX\n", AUX, naux);
          naux++;
        };
      j++;
      if (j >= 20)
        j = 0;

    };

  /* DSSD IDs */

  i1 = 1;
  for (i = 2000; i <= 2319; i++)
    {
      printf ("%4i %4i %4i DSSD\n", i, DSSD, i1);
      i1++;
    };

  /* Focal Plane IDs */

  i1 = 1;
  for (i = 2320; i <= 2339; i++)
    {
      printf ("%4i %4i %4i FP\n", i, FP, i1);
      i1++;
    };

  /* Focal Plane IDs */

  i1 = 1;
  for (i = 2340; i <= 2359; i++)
    {
      printf ("%4i %4i %4i XARRAY\n", i, XARRAY, i1);
      i1++;
    };

  /* Chico2 */

  i1 = 1;
  for (i = 3000; i <= 3010; i++)
    {
      printf ("%4i %4i %4i CHICO2\n", i, CHICO2, i1);
      i1++;
    };


  /* done */

  return (0);

}
