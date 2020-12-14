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

#include "ctk.h"
#include "veto_pos.h"
#include "GEBFilter.h"

#define MAXGTMODNO 30

extern unsigned int *veto_cube;
extern GEBMERGEPAR Pars;

int
setup_veto_cube ()
{
  unsigned int ui1, i1;
  int st, i, j, k, l, m, nr = 0;
  float xx, yy, zz;
  FILE *fp;

  printf ("setup_veto_cube: attempt to setup veto cube\n");
  fflush (stdout);


/* open veto_cube file */

  fp = fopen (Pars.vetocubefn, "r");
  if (fp == NULL)
    {
      printf ("could not open %s, set Pars.vetocubes=0 and go on\n", Pars.vetocubefn);
      Pars.vetocubes=0;
      return(0);
    }
  else
    printf ("\"%s\" is open\n", Pars.vetocubefn);
  fflush (stdout);

  /* allocate and zero out the veto cube */

  ui1 = (MAXGTMODNO + 1) * (MAXCRYSTALNO + 1) * (VETO_NX + 1) * (VETO_NY + 1) * (VETO_NZ + 1) * sizeof (unsigned int);
  veto_cube = (unsigned int *) calloc (ui1, 1);
  printf ("allocated %11u bytes  for veto_cube\n", ui1);
  printf ("       or %11u Mbytes for veto_cube\n", ui1 / 1024 / 1024);
  ui1 = MAXGTMODNO * MAXCRYSTALNO * VETO_NX * VETO_NY * VETO_NZ;
  printf ("max index: %u\n", ui1);
  printf ("VETO_X_D = %f\n", VETO_X_D);
  printf ("VETO_Y_D = %f\n", VETO_Y_D);
  printf ("VETO_Z_MAX = %f\n", VETO_Z_MAX);
  printf ("VETO_BINWIDTH = %f\n", VETO_BINWIDTH);
  printf ("VETO_NX = %i\n", VETO_NX);
  printf ("VETO_NY = %i\n", VETO_NY);
  printf ("VETO_NZ = %i\n", VETO_NZ);


  st = fscanf (fp, "%i %i %i %i %i ; %i ", &i, &j, &k, &l, &m, &i1);
  while (st == 6)
    {

      /* read and fill veto_cube with markers (not counts) */

      if (i1 > 0)
        *(veto_cube + VETO_INDX (i, j, k, l, m)) = 1;
      nr++;

      if (nr < 100)
        {
          printf (" %2i %2i %3i %3i %3i = %3i INDX[%i]=", i, j, k, l, m, i1, VETO_INDX (i, j, k, l, m));
          printf ("%i\n", *(veto_cube + VETO_INDX (i, j, k, l, m)));
        };

      /* next */

      st = fscanf (fp, "%i %i %i %i %i ; %i ", &i, &j, &k, &l, &m, &i1);

    };

  printf ("read %i values from \"%s\"\n", nr, Pars.vetocubefn);

  return(0);

};


