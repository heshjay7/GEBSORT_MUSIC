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

#include "gdecomp.h"

#define MAX_SEGS   8            /* max. number of segments to take in events */
#define MAX_PARS   (8*MAX_SEGS+1)       /* max. number of parameters to fit */
#define MAX_INTPTS (2*MAX_SEGS) /* max. number of interaction points */
#define PMODE 0644
#define GEB_TYPE_DECOMP         1
#define MAXDETPOS 30
#define MAXCRYSTALNO 3
#define MAXDETNO 180+2
#define NOMOREDATA 99

#define MAXL 1000
#define SPLEN 2000
#define STRLEN 1024
#define NBIN 10
#define MAX_SIM_GAMMAS 10
#define MAXPAYLOADSIZE 30000

/* #define NEWFORMAT 2 */

float minr, emin, eSmear_seg, eSmear_CC, pSmear, maxdist, vc;
int maxevents, minNoInteractions, maxNoInteractions;
int minNoCrystals, maxNoCrystals;
int nThresholds, evmod, noworldtocrystalrot = 0;
float e_th[MAXDETNO], e_0[MAXDETNO];
int nResolution;
float enRes_0[120], enRes_1[1000];
int NPRINT = 20;
off_t inData;
int enabled[MAXDETNO + 1];

double xx_mean[1000];
double yy_mean[1000];
long long int nn_mean[1000];

int ngammas;
double etot;



/*-----------------------------------------------------------------------------*/

int
get_event_GT_ascii3 (FILE * fp, int *nhits, float *egamma, int *evno, int id[MAXL], float eg[MAXL], float xx[MAXL],
                     float yy[MAXL], float zz[MAXL])
{

  /* for UCGretina, Jul 2014 ascii format lr3.3 version for source data */

  /* declarations */

  char str[STRLEN], str1[STRLEN], str2[STRLEN], str3[STRLEN], str4[STRLEN], *pc;
  int i, j, ndp = 0, nip = 0, idummy, crystid;
  static int nn = 0, nlines = 0;

  bzero (str, STRLEN);

  if (nn < 200)
    printf ("entered get_event_GT_ascii3\n");

  /* search for a line starting with D */

  i = 0;
  pc = fgets (str, STRLEN, fp);
  nlines++;
  while (str[0] != 68)
    {
      i++;
      pc = fgets (str, STRLEN, fp);
      nlines++;
      if (pc == 0)
        return (NOMOREDATA);
//      printf ("__%s", str);
      if (str[0] == 69)
        {
        sscanf (str, "E %i %i %i", &idummy, &idummy, evno);
//        printf("== evno= %i\n", *evno);
        };
    };
  if (nn < 200 && 0)
    printf ("skipped %i lines to get to a D line \n", i);
  sscanf (str, "S %s %s %s %s %i", str1, str2, str3, str4, &idummy);

  /* interpret the D line */

  sscanf (str, "D %i %i", &ndp, evno);

  *nhits = 0;
  *egamma = 0;
  for (i = 0; i < ndp; i++)
    {
      pc = fgets (str, STRLEN, fp);
      nlines++;
//    printf("%s", str);
      sscanf (str, "C  %i %i", &crystid, &nip);
      for (j = 0; j < nip; j++)
        {
          id[*nhits] = crystid;
//          printf("crystid=%i, nip=%i, id[*nhits]=%i\n",crystid, nip,id[*nhits]);
          pc = fgets (str, STRLEN, fp);
          nlines++;
//      printf("%s", str);
          sscanf (str, "%i %f %f %f %f %f", &idummy, &eg[*nhits], &xx[*nhits], &yy[*nhits], &zz[*nhits]);
//      printf("## %i\n", *nhits);
//      printf("== id %i e %f x %f y %f z %f \n",id[*nhits],eg[*nhits],xx[*nhits],yy[*nhits],zz[*nhits]);
          etot += eg[*nhits];
          (*nhits)++;
        }

    };

  if (nn < 200)
    {
      printf ("++get_event_GT_ascii3: hit # %i, nlines= %i\n", nn, nlines);
      for (i = 0; i < *nhits; i++)
        {
          printf ("++%2i: [%3i] e %8.2f x %8.2f y %8.2f z %8.2f\n", i, id[i], eg[i], xx[i], yy[i], zz[i]);
        };
    };
  nn++;
//  if (*evno > 10)exit (0);
//  if (nn > 3)exit (0);
//   printf("%c, %i, *egamma=%f\n", pc,*evno,*egamma);

  /* keep record of the mean xy positions */
  /* for ID check */

  for (i = 0; i < *nhits; i++)
    {
      nn_mean[id[i]]++;
      xx_mean[id[i]] += xx[i];
      yy_mean[id[i]] += yy[i];
    };

  for (i = 0; i < *nhits; i++)
    if (id[i] <= 0 || id[i] > 123)
      {
        printf ("bad id= %i at lines read %i\n", id[i], nlines);
        assert (id[i] > 0 && id[i] <= 123);
      };

  for (i = 0; i < *nhits; i++)
    assert (id[i] > 0 && id[i] <= 123);

  return (0);
}

/*-----------------------------------------------------------------------------*/

int
get_event_GT_ascii2 (FILE * fp, int *nhits, float *egamma, int *evno, int id[MAXL], float eg[MAXL], float xx[MAXL],
                     float yy[MAXL], float zz[MAXL])
{

  /* for UCGretina, Jul 2014 ascii format lr3.3 version for in beam data */

  /* declarations */

  char str[STRLEN], str1[STRLEN], str2[STRLEN], str3[STRLEN], str4[STRLEN], *pc;
  int i, j, ndp = 0, nip = 0, idummy, crystid;
  static int nn = 0, nlines = 0;

  bzero (str, STRLEN);

  if (nn < 200)
    printf ("entered get_event_GT_ascii2\n");

  /* search for a line starting with S */

  i = 0;
  pc = fgets (str, STRLEN, fp);
  nlines++;
  while (str[0] != 83)
    {
      i++;
      pc = fgets (str, STRLEN, fp);
      nlines++;
      if (pc == 0)
        return (NOMOREDATA);
//      printf ("__%s", str);
      if (str[0] == 69)
        {
        sscanf (str, "E %i %i %i", &idummy, &idummy, evno);
//        printf("== evno= %i\n", *evno);
        };
    };
  if (nn < 200 && 0)
    printf ("skipped %i lines to get to a S line \n", i);
  sscanf (str, "S %s %s %s %s %i", str1, str2, str3, str4, &idummy);

  /* next line is  D */

  pc = fgets (str, STRLEN, fp);
  ndp = 0;
  sscanf (str, "D %i %i", &ndp, evno);
//  printf ("*ndp=%i, evno=%i\n", ndp, *evno);


  *nhits = 0;
  *egamma = 0;
  for (i = 0; i < ndp; i++)
    {
      pc = fgets (str, STRLEN, fp);
      nlines++;
//    printf("%s", str);
      sscanf (str, "C  %i %i", &crystid, &nip);
      for (j = 0; j < nip; j++)
        {
          id[*nhits] = crystid;
//          printf("crystid=%i, nip=%i, id[*nhits]=%i\n",crystid, nip,id[*nhits]);
          pc = fgets (str, STRLEN, fp);
          nlines++;
//      printf("%s", str);
          sscanf (str, "%i %f %f %f %f %f", &idummy, &eg[*nhits], &xx[*nhits], &yy[*nhits], &zz[*nhits]);
//      printf("## %i\n", *nhits);
//      printf("== id %i e %f x %f y %f z %f \n",id[*nhits],eg[*nhits],xx[*nhits],yy[*nhits],zz[*nhits]);
          etot += eg[*nhits];
          (*nhits)++;
        }

    };

  if (nn < 200)
    {
      printf ("++get_event_GT_ascii2: hit # %i, nlines= %i\n", nn, nlines);
      for (i = 0; i < *nhits; i++)
        {
          printf ("++%2i: [%3i] e %8.2f x %8.2f y %8.2f z %8.2f\n", i, id[i], eg[i], xx[i], yy[i], zz[i]);
        };
    };
  nn++;
//  if (*evno > 10)exit (0);
//  if (nn > 3)exit (0);
//   printf("%c, %i, *egamma=%f\n", pc,*evno,*egamma);

  /* keep record of the mean xy positions */
  /* for ID check */

  for (i = 0; i < *nhits; i++)
    {
      nn_mean[id[i]]++;
      xx_mean[id[i]] += xx[i];
      yy_mean[id[i]] += yy[i];
    };

  for (i = 0; i < *nhits; i++)
    if (id[i] <= 0 || id[i] > 123)
      {
        printf ("bad id= %i at lines read %i\n", id[i], nlines);
        assert (id[i] > 0 && id[i] <= 123);
      };

  for (i = 0; i < *nhits; i++)
    assert (id[i] > 0 && id[i] <= 123);

  return (0);
}

/*-----------------------------------------------------------------------------*/

int
get_event_AGT (FILE * fp, int *nhits, float *egamma, int *evno, int id[MAXL], float eg[MAXL], float xx[MAXL],
               float yy[MAXL], float zz[MAXL])
{

  int i, i1, nret, type;
  static int firsttime = 1, nn = 0;;
  static char str[STRLEN], header[STRLEN];
  char *p, *pc;
  float r1, r2, r3;

  /* skip the long header the first time we are called */

  if (firsttime)
    {
      printf ("get_event_AGT, firsttime activity\n");
      firsttime = 0;
      pc = fgets (str, STRLEN, fp);
      while ((p = strstr (str, "$")) == NULL)
        {
          printf ("%s", str);
          pc = fgets (str, STRLEN, fp);
          nn++;
        };
      printf ("get_event_AGT: skipped %i lines\n", nn);

      /* now read the header */

      pc = fgets (header, STRLEN, fp);
      printf ("get_event_AGT: first header: %s\n", header);

    };

  /* interpret header */

  nret = sscanf (header, "%i %f %f %f %f %i ", &type, egamma, &r1, &r2, &r3, evno);
//  *evno=i1;
//  printf ("%i %f %f %f %f %i \n", type, *egamma, r1, r2, r3, i1); fflush(stdout);
  if (nret != 6)
    {
      printf ("get_event_AGT: failed to read header...\n");
      fflush (stdout);
      return (NOMOREDATA);
    };
  assert (nret == 6);

  /* read this event (until we see a -1) */

  i = 0;
  while (1)
    {

      /*read next line */

//      printf ("get_event_AGT: attempting to read next line, fp=%p, str=%p\n", fp, str);
//      fflush (stdout);
      pc = fgets (str, STRLEN, fp);
      nn++;
//      printf ("get_event_AGT: %i, %p: %s", nn, pc, str);
//      fflush (stdout);

      /* done */

      if (pc == NULL)
        {
          printf ("get_event_AGT: failed to read more data...\n");
          fflush (stdout);
          return (NOMOREDATA);
        };


      /* interpret */

//      printf ("get_event_AGT: i=%i\n", i);
//      fflush (stdout);
      nret = sscanf (str, "%i %f %f %f %f %i", &id[i], &eg[i], &xx[i], &yy[i], &zz[i], &type);
//      printf ("get_event_AGT: ->  %i %f %f %f %f %i\n", id[i], eg[i], xx[i], yy[i], zz[i], type);
//      fflush (stdout);
      if (id[i] == -1)
        {

          /* we hit the next header */

//          printf ("get_event_AGT: done reading interaction points\n");
//          fflush (stdout);
          strcpy (header, str);
          *nhits = i;
//          printf ("get_event_AGT: return 0\n");
//          TS.listfflush (stdout);
          return (0);

        }
      else
        {

          /* this was just an interaction point */

          i++;
//          printf ("get_event_AGT: i now %i\n", i);
//          fflush (stdout);

          id[i] += 4;
        };

    };


  return (0);


}


/*-----------------------------------------------------------------------------*/

int
get_event_GT_ascii (FILE * fp, int *nhits, float *egamma, int *evno, int id[MAXL], float eg[MAXL], float xx[MAXL],
                    float yy[MAXL], float zz[MAXL])
{

  /* this function reads the ascii output */
  /* from LRs UCGRETINA. Notice the ascii  */
  /* format changed at soem point, this will  */
  /* read the new format */


  /* declarations */

  int nret, i, dummy;
  static unsigned int nlines = 0;

  /* get header */

#if(NEWFORMAT==1)
  dummy = 0;
  nret = fscanf (fp, "$ %i %f %i %i\n", nhits, egamma, evno, &dummy);
//  printf("%i %f %i %i\n", *nhits, *egamma, *evno, dummy);
  if (nret != 4)
    return (NOMOREDATA);
#else
  nret = fscanf (fp, "$ %i %f %i\n", nhits, egamma, evno);
//  printf("%i %f %i\n", *nhits, *egamma, *evno);
  if (nret != 3)
    return (NOMOREDATA);
#endif
//  printf ("nret=%i, nhits=%3i, e=%9.2f, evno= %i\n", nret, *nhits, *egamma, *evno);
  nlines++;

  /* read event */

  for (i = 0; i < *nhits; i++)
    {
#if(NEWFORMAT==1)
      /* latest has segment id */

      nret = fscanf (fp, "%i %i %f %f %f %f\n", &id[i], &dummy, &eg[i], &xx[i], &yy[i], &zz[i]);
      if (nret != 6)
        {
          printf ("data read error, need nret=6, but got %i, on line %i\n", nret, nlines);
          exit (1);
        };
#else
      nret = fscanf (fp, "%i %f %f %f %f\n", &id[i], &eg[i], &xx[i], &yy[i], &zz[i]);
//      printf("id[i]=%i\n",id[i]);

      if (nret != 5)
        {
          printf ("data read error, need nret=5, but got %i, on line %i\n", nret, nlines);
          exit (1);
        };

#endif
      nlines++;
//     printf ("id=%i e=%9.2f x=%8.2f y=%8.2f z=%8.2f\n", id[i], eg[i], xx[i], yy[i], zz[i]);


    };

  return (0);
}

/*-----------------------------------------------------------------------------*/

int
get_event_GT_bin (FILE * fp, int *nhits, float *egamma, int *evno, int id[MAXL], float eg[MAXL], float xx[MAXL],
                  float yy[MAXL], float zz[MAXL], long long int *TS)
{

  /* This reader handles the binary format  */
  /* ID==11 from LR's UCGRETINA code */

  /* declarations */

  int nret, i, dummy, siz, j;
  static unsigned int nlines = 0, nn = 0;
  GEBDATA ptgd;
  int *pi;
  char dustbin[100000];
  CRYS_INTPTS mode2;

  /* LRs G4 data stuctures */

  typedef struct g4sim_emitted_gamma
  {
    float e;
    float x, y, z;
    float phi, theta;
  } EG;

  typedef struct g4sim_abcd1234
  {
    int type;                   /* defined as abcd1234 */
    int num;                    /* # of emitted gammas */
    EG gammas[MAX_SIM_GAMMAS];
  } G4SIM_EGS;

  G4SIM_EGS ptpl;

  nn++;
  if (nn < 20)
    {
      printf ("entered get_event_GT_bin\n");
    }

  *nhits = 0;

  /* read the geb header */

  siz = read (inData, (char *) &ptgd, sizeof (GEBDATA));
  if (siz != sizeof (GEBDATA))
    {
      printf ("failed to read %i bytes for header, got %i, quit\n", sizeof (GEBDATA), siz);
      return (NOMOREDATA);
    };

  if (nn < 20)
    {
      printf ("read GEB header of size=%i; ", siz);
      printf ("TS=%lli, ", ptgd.timestamp);
      printf ("length= %i, type=%i\n", ptgd.length, ptgd.type);
    }

  /* pass back the timestamp */

  *TS = ptgd.timestamp;

  /* read the payload */

  if (ptgd.type == GEB_TYPE_G4SIM)
    {

      if (nn < 20)
        printf ("got GEB_TYPE_G4SIM data\n");
//      assert (ptgd.length == sizeof (G4SIM_EGS));
      siz = read (inData, (char *) &ptpl, ptgd.length);
      if (nn < 20)
        printf ("read GEB payload of size=%i\n", siz);

      if (nn < 20)
        {
          printf ("type=0x%x\n", ptpl.type);
          printf ("num=%i\n", ptpl.num);
          for (i = 0; i < ptpl.num; i++)
            {
              printf ("%i: ", i);
              printf ("%5.1f; ", ptpl.gammas[i].e);
              printf ("%5.1f, ", ptpl.gammas[i].x);
              printf ("%5.1f, ", ptpl.gammas[i].y);
              printf ("%5.1f ... ", ptpl.gammas[i].z);
              printf ("%5.1f, ", ptpl.gammas[i].phi);
              printf ("%5.1f ", ptpl.gammas[i].theta);
              printf ("\n");
            };
        }
      return (1);
    }
  else if (ptgd.type == GEB_TYPE_DECOMP)
    {
      if (nn < 20)
        printf ("got GEB_TYPE_DECOMP data\n");
      siz = read (inData, (char *) &mode2, ptgd.length);
      if (nn < 20)
        printf ("read GEB payload of size=%i into mode2 structure\n", siz);

      assert (mode2.type == 0xabcd5678);

      /* return info in old format for processing */

      *egamma = 0;
      *evno = 0;
      *nhits = mode2.num;
      assert (mode2.num > 0);
      for (i = 0; i < mode2.num; i++)
        {
          id[i] = mode2.crystal_id;
          xx[i] = mode2.intpts[i].x;
          yy[i] = mode2.intpts[i].y;
          zz[i] = mode2.intpts[i].z;
          eg[i] = mode2.intpts[i].e;
          *egamma += mode2.intpts[i].e;
//          printf ("%f %f %f %f\n", eg[i], xx[i], yy[i], zz[i]);
        };

      return (0);

    }
  else if (ptgd.type == GEB_TYPE_S800PHYSDATA)
    {
      if (nn < 20)
        printf ("got GEB_TYPE_S800PHYSDATA data, drop it\n");
      siz = read (inData, (char *) dustbin, ptgd.length);
      if (nn < 20)
        printf ("read GEB payload of size=%i\n", siz);
      return (1);
    }
  else
    {
      printf ("dont know what to do with data of type %i, drop it\n", ptgd.type);
      siz = read (inData, (char *) dustbin, ptgd.length);
      printf ("read GEB payload of size=%i\n", siz);
      return (1);
    }


  return (1);
}

/*-----------------------------------------------------------------------------*/

float
ranGauss (void)
{

  double x1, x2, w, val;
  static int nn = 0;
  static float y1, y2;

  if (nn == 0)
    {

      do
        {
//          x1 = 2.0 * (double) rand () / RAND_MAX - 1.0;
//          x2 = 2.0 * (double) rand () / RAND_MAX - 1.0;
          x1 = 2.0 * drand48 () - 1.0;
          x2 = 2.0 * drand48 () - 1.0;
          w = x1 * x1 + x2 * x2;
        }
      while (w >= (double) 1.0);

      w = sqrt ((-2.0 * log (w)) / w);
      y1 = x1 * w;
      y2 = x2 * w;

//  printf("%f %f, %f %f; %f\n", x1, x2, y1, y2, w);

      nn = 1;
      return ((float) y1);
    }
  else
    {
      nn = 0;
      return ((float) y2);
    }


}

/*--------------------------------------------------------*/

void
CheckNoArgs (int required, int actual, char *str)
{

  if (required < actual)
    {
      printf ("argument problem with chat option\n");
      printf ("--> %s\n", str);
      printf ("required # arguments: %i\n", required - 1);
      printf ("actual   # arguments: %i\n", actual - 1);
      printf ("Please fix and try again, quitting...\n");
      exit (1);
    };

}

/*--------------------------------------------------------*/

int
readChatFile (char *name)
{

  /* declarations */

  FILE *fp, *fp1;
  char *p, *pc, str[STRLEN], str1[STRLEN], str2[STRLEN];
  int echo = 0, nret = 0, nni = 0, nn = 0;
  int i1, i;
  float r2, r3;
  int detID;
  float et, e0;
  float res0, res1;

  /* prototypes */

  int str_decomp (char *, int, int *, int);

  /* open chat file */

  if ((fp = fopen (name, "r")) == NULL)
    {
      printf ("error: could not open chat file: <%s>\n", name);
      exit (0);
    };
  printf ("chat file: <%s> open\n", name);
  printf ("\n");
  fflush (stdout);

  /* read content and act */

  pc = fgets (str, STRLEN, fp);

  while (pc != NULL)
    {
      if (echo)
        printf ("chat->%s", str);
      fflush (stdout);

      /* attemp to interpret the line */

      if ((p = strstr (str, "echo")) != NULL)
        {
          echo = 1;
          if (echo)
            printf ("will echo command lines\n");
          fflush (stdout);

        }
      else if (str[0] == 35)
        {
          /* '#' comment line, do nothing */
          nni--;                /* don't count as instruction */

        }
      else if (str[0] == 59)
        {
          /* ';' comment line, do nothing */
          nni--;                /* don't count as instruction */

        }
      else if (str[0] == 10)
        {
          /* empty line, do nothing */
          nni--;                /* don't count as instruction */

        }
      else if ((p = strstr (str, "minr")) != NULL)
        {
          nret = sscanf (str, "%s %f", str1, &minr);
          CheckNoArgs (nret, 2, str);
        }
      else if ((p = strstr (str, "vc")) != NULL)
        {
          nret = sscanf (str, "%s %f", str1, &vc);
          CheckNoArgs (nret, 2, str);
        }
      else if ((p = strstr (str, "emin")) != NULL)
        {
          nret = sscanf (str, "%s %f", str1, &emin);
          CheckNoArgs (nret, 2, str);
        }
      else if ((p = strstr (str, "enabled")) != NULL)
        {
          nret = sscanf (str, "%s %s", str1, str2);
          CheckNoArgs (nret, 2, str);
          str_decomp (str2, MAXDETNO + 1, enabled, 1);
          printf ("enabled crystals\n");
          for (i = 0; i <= MAXDETNO; i++)
            if (enabled[i])
              printf ("[%i]", i);
          printf ("\n");
        }
      else if ((p = strstr (str, "eSmear_seg")) != NULL)
        {
          nret = sscanf (str, "%s %f", str1, &eSmear_seg);
          CheckNoArgs (nret, 2, str);
        }
      else if ((p = strstr (str, "eSmear_CC")) != NULL)
        {
          nret = sscanf (str, "%s %f", str1, &eSmear_CC);
          CheckNoArgs (nret, 2, str);
        }
      else if ((p = strstr (str, "pSmear")) != NULL)
        {
          nret = sscanf (str, "%s %f", str1, &pSmear);
          CheckNoArgs (nret, 2, str);
        }
      else if ((p = strstr (str, "maxdist")) != NULL)
        {
          nret = sscanf (str, "%s %f", str1, &maxdist);
          CheckNoArgs (nret, 2, str);
        }
      else if ((p = strstr (str, "evmod")) != NULL)
        {
          nret = sscanf (str, "%s %i", str1, &evmod);
          CheckNoArgs (nret, 2, str);
        }
      else if ((p = strstr (str, "maxevents")) != NULL)
        {
          nret = sscanf (str, "%s %i", str1, &maxevents);
          CheckNoArgs (nret, 2, str);
        }
      else if ((p = strstr (str, "noworldtocrystalrot")) != NULL)
        {
          nret = sscanf (str, "%s", str1);
          noworldtocrystalrot = 1;
          CheckNoArgs (nret, 1, str);
        }
      else if ((p = strstr (str, "nprint")) != NULL)
        {
          nret = sscanf (str, "%s %i", str1, &NPRINT);
          CheckNoArgs (nret, 2, str);
        }
      else if ((p = strstr (str, "minNoInteractions")) != NULL)
        {
          nret = sscanf (str, "%s %i", str1, &minNoInteractions);
          CheckNoArgs (nret, 2, str);
        }
      else if ((p = strstr (str, "maxNoInteractions")) != NULL)
        {
          nret = sscanf (str, "%s %i", str1, &maxNoInteractions);
          CheckNoArgs (nret, 2, str);
        }
      else if ((p = strstr (str, "minNoCrystals")) != NULL)
        {
          nret = sscanf (str, "%s %i", str1, &minNoCrystals);
          CheckNoArgs (nret, 2, str);
        }
      else if ((p = strstr (str, "maxNoCrystals")) != NULL)
        {
          nret = sscanf (str, "%s %i", str1, &maxNoCrystals);
          CheckNoArgs (nret, 2, str);
          printf ("maxNoCrystals set to %i\n", maxNoCrystals);
        }
      else if ((p = strstr (str, "thresholds")) != NULL)
        {
          nret = sscanf (str, "%s %i", str1, &nThresholds);
          CheckNoArgs (nret, 2, str);
          for (i = 0; i < nThresholds; i++)
            {
              nn++;             /* line counter */
              pc = fgets (str, STRLEN, fp);     /* read thresholds */
              nret = sscanf (str, "%i %f %f", &detID, &et, &e0);
              CheckNoArgs (nret, 3, str);
              e_th[detID] = et;
              e_0[detID] = e0;
            }
        }
      else if ((p = strstr (str, "crystalResPar")) != NULL)
        {
          nret = sscanf (str, "%s %i", str1, &nResolution);
          CheckNoArgs (nret, 2, str);
          for (i = 0; i < nResolution; i++)
            {
              nn++;             /* line counter */
              pc = fgets (str, STRLEN, fp);     /* read thresholds */
              nret = sscanf (str, "%i %f %f", &detID, &res0, &res1);
              CheckNoArgs (nret, 3, str);
              enRes_0[detID] = res0;
              enRes_1[detID] = res1;
            }
        }
      else
        {

          /*-----------------------------*/
          /* chatscript read error point */
          /*-----------------------------*/

          printf ("line %2.2i in chat script, option :%s \n__not understood\n", nn, str);
          printf ("%i\n", str[0]);
          printf ("aborting\n");
          fflush (stdout);
          exit (0);
        };

      /* read next line in chat script */

      nn++;                     /* line counter */
      nni++;                    /* instruction counter */
      pc = fgets (str, STRLEN, fp);

    };



  /* done */

  fclose (fp);
  printf ("\n");
  printf ("chat file: <%s> closed\n", name);
  printf ("__processed %i sort instructions and %i lines\n", nni, nn);
  printf ("\n");
  fflush (stdout);
  return (0);

}


/*-----------------------------------------------------------------------------*/

int
main (int argc, char **argv)
{
  FILE *fp;
  int nhits, evno, nret, evcount = 0, i1, i2, ndet, tooFar, ndet_sp[20], in;
  int i, j, k, l, m, ng = 0, readcount = 0, nPoints = 0, nit;
  float egamma, sp_obs[SPLEN], sp_emit[SPLEN], spg[SPLEN], scatter_l[SPLEN], r1, r2, r3;
  float sp_obs_nosmear[SPLEN], sp_emit_nosmear[SPLEN];
  float esum, esum_smeared, egnow;
  float eg[MAXL], xx[MAXL], yy[MAXL], zz[MAXL], maxPackDist[MAXL];
  int id[MAXL], pos[MAXL], hitnum[MAXL], valid[MAXL];
  float realx[MAXL], realy[MAXL], realz[MAXL], dr;
  int eminok = 0, eminnotok = 0;
  GEBDATA *gd;
  CRYS_INTPTS *payload;
  long long int ts = 1000000, TS;
  off_t outData;
  float lx, ly, lz;
  unsigned int seed;
  int siz = 0, ilok = 0, ilnotok = 0, producedNumIntPoints = 0, readNumIntPoints = 0;
  int trueNomGammas = 0;
  double nphoto, ncompton;
  float polang, doppler_factor, dp, radius, emax;
  FILE *fp1;
  char str[512];
  float thresh;
  int numclusters, nn, j1, j2;
  float x1, y1, z1;
  float crmat[MAXDETPOS + 1][MAXCRYSTALNO + 1][4][4];
  float mat[3][3];
  int holeNum, crystalNumber, eOK = 1;
  float sppd[SPLEN];
  int CChitrate[MAXDETNO], numG4Events, zero = 0;
  char buffer[512];
  FILE *fp0;
  long long int id_hit[MAXDETNO];


  int wr_spe (char *, int *, float *);
  int printCRYS_INTPTS (FILE *, CRYS_INTPTS *, GEBDATA *);
  int get_a_seed (unsigned int *);

  if (argc != 5)
    {
      printf ("use: %s input output chatfile #events\n", argv[0]);
      exit (0);
    };

  numG4Events = atoi (argv[4]);
  printf ("G4 data has %i events [INPUT, NOT COUNTED]\n", numG4Events);

  /* init */

  etot=0;
  for (i = 0; i < MAXDETNO; i++)
    {
      xx_mean[i] = 0;
      yy_mean[i] = 0;
      nn_mean[i] = 0;
      id_hit[i] = 0;
    };

  for (i = 0; i < 20; i++)
    ndet_sp[i] = 0;

  for (i = 0; i < MAXDETNO; i++);
  CChitrate[i] = 0;

  for (i = 0; i < SPLEN; i++)
    {
      sppd[i] = 0;
      sp_obs[i] = 0;
      sp_emit[i] = 0;
      sp_obs_nosmear[i] = 0;
      sp_emit_nosmear[i] = 0;
      spg[i] = 0;
      scatter_l[i] = 0;
    };

  printf ("doing (GEBDATA *) calloc (1, sizeof (GEBDATA)); ... ");
  fflush (stdout);
  gd = (GEBDATA *) calloc ((size_t) 1, sizeof (GEBDATA));
  printf (" done\n");
  fflush (stdout);
  i1 = sizeof (CRYS_INTPTS);
  i1 = 10000;                   /* make MAXPAYLOADSIZEsure there is enough room */
  //  payload = (CRYS_INTPTS *) calloc (1, 1);  //LR
  printf ("doing (CRYS_INTPTS *) calloc (1, i1); ... ");
  fflush (stdout);
  payload = (CRYS_INTPTS *) calloc ((size_t) 1, (size_t) i1);
  printf (" done\n");
  fflush (stdout);
  printf ("allocated %i bytes for payload\n", i1);
  printf ("zero payload...");
  fflush (stdout);
  bzero ((void *) payload, MAXPAYLOADSIZE);
  printf ("done\n");
  fflush (stdout);

  /* initialize random number generator */

  get_a_seed (&seed);
  srand (seed);

  ngammas = 0;
  vc = 0;
  minr = 5.0;
  emin = 5;
  eSmear_seg = 1.0;
  eSmear_CC = 1.0;
  pSmear = 1.0;
  maxdist = 300;
  maxevents = 2000000000;
  minNoInteractions = 1;
  maxNoInteractions = 8;
  minNoCrystals = 1;
  maxNoCrystals = 10;
  nThresholds = 0;
  for (i = 0; i < 120; i++)
    {
      e_th[i] = 0;
      e_0[i] = 0;
      enabled[i] = 1;
    }
  evmod = 1;

  /* open chat file and read */

  printf ("reading chatfile...");
  fflush (stdout);
  readChatFile (argv[3]);
  printf ("done\n");
  fflush (stdout);

  /* to rotate or to rotate not */

#if(0)
  printf ("[1]noworldtocrystalrot=%i\n", noworldtocrystalrot);
#if(USGEANTFORMAT==1 && ASCIIDATA==0)
  noworldtocrystalrot = 1;
#endif
#if(USGEANTFORMAT==0 && ASCIIDATA==1)
  noworldtocrystalrot = 0;
#endif
#endif

  printf ("USGEANTFORMAT=%i\n", USGEANTFORMAT);
  printf ("ASCIIDATA=%i\n", ASCIIDATA);
  printf ("[2]noworldtocrystalrot=%i\n", noworldtocrystalrot);

  i1 = 0;
  for (i = 121; i < MAXDETNO; i++)
    if (enabled[i])
      i1++;
  if (i1 > 0 && noworldtocrystalrot == 0)
    {
      printf ("%s:\n", argv[0]);
      printf ("cannot handle crystal numbers > 120\n");
      printf ("without the noworldtocrystalrot option on\n");
      printf ("you have %i crystal enabled above 120\n", i1);
      printf ("quit\n");
      exit (1);
    }



  /* log parameters */

  printf ("vc               =%f\n", vc);
  printf ("minr             =%f\n", minr);
  printf ("emin             =%f\n", emin);
  printf ("eSmear_seg       =%f\n", eSmear_seg);
  printf ("eSmear_CC        =%f\n", eSmear_CC);
  printf ("pSmear           =%f\n", pSmear);
  printf ("maxdist          =%f\n", maxdist);
  printf ("maxevents        =%i\n", maxevents);
  printf ("minNoInteractions=%i\n", minNoInteractions);
  printf ("maxNoInteractions=%i\n", maxNoInteractions);
  printf ("minNoCrystals    =%i\n", minNoCrystals);
  printf ("maxNoCrystals    =%i\n", maxNoCrystals);
  printf ("nThresholds      =%i\n", nThresholds);
  if (nThresholds > 0)
    {
      for (i = 0; i < 120; i++)
        {
          if (e_th[i] > 0)
            printf ("                  %i  %f  %f\n", i, e_th[i], e_0[i]);
        }
    }
  printf ("crystalResPar    =%i\n", nResolution);
  if (nResolution > 0)
    {
      for (i = 0; i < 120; i++)
        {
          if (enRes_0[i] > 0)
            printf ("                  %i  %f  %f\n", i, enRes_0[i], enRes_1[i]);
        }
    }

#if(0)

  /* test gaussian random number generator */

  for (i = 0; i < 10000000; i++)
    {
      emin r1 = 300.0 * ranGauss ();
      i1 = (int) (r1 + 0.5);
//      if (i1 )
//       printf ("%i, %f\n", i1,spg[i1]);
      i1 += SPLEN / 2;
      if (i1 > 0 && i1 < SPLEN)
        spg[i1]++;
    };
  i1 = SPLEN;
  wr_spe ("gauss.spe", &i1, spg);
  printf ("wrote \"gauss.spe\"\n");
  if (1)
    exit (0);
#endif

  /* open input file from G4 */

#if(ASCIIDATA==1)

  /* ascii input file */

  printf ("will attempt to open file \"%s\" \n", argv[1]);
  fflush (stdout);
  if ((fp = fopen (argv[1], "r")) == NULL)
    {
      printf ("could not open \"%s\" as input file\n", argv[1]);
      exit (1);
    }
  else
    printf ("ascii \"%s\" is open as input file\n", argv[1]);
#endif

#if(ASCIIDATA==0)

  /* binary ID=11 LR UCGRETINA file */

  inData = open (argv[1], O_RDONLY, 0);
  if (inData == -1)
    {
      printf ("could not open\"%s\"; quit\n", argv[1]);
      exit (1);
    }
  else
    printf ("binary input file \"%s\" is open, inData=%i\n", argv[1], (int) inData);

#endif

  /* open output file */

  outData = 0;
  outData = open ((char *) argv[2], O_WRONLY | O_CREAT | O_TRUNC, PMODE);
  if (outData == 0)
    {
      printf ("could not open output data file %s, quit!\n", argv[2]);
      exit (1);
    }
  else
    printf ("\"%s\" is open as out file\n", argv[2]);

#if(0)
  nret = 0;
  nn = 0;
  while (nret == 0)
    {
      nret = get_event_AGT (fp, &nhits, &egamma, &evno, id, eg, xx, yy, zz);
      nn++;
      printf ("nhits=%3i, e=%9.2f, evno= %i\n", nhits, egamma, evno);
    };
  if (1)
    exit (0);
#endif

  /* get the rotation matrices */

#if(USGEANTFORMAT==1)
  sprintf (str, "crmat.LINUX");
  in = open (str, O_RDONLY, 0);
  if (in > 0)
    printf ("%s is open (input) binary format\n", str);
  else
    {
      printf ("could not open %s\n", str);
      exit (1);
    };
  siz = read (in, (char *) crmat, sizeof (crmat));
  printf ("read %i bytes into crmat\n", siz);
  close (in);
#endif

#if(USGEANTFORMAT==0)
  printf ("this should really be AGATA crmats...\n");
  sprintf (str, "crmat.LINUX");
  in = open (str, O_RDONLY, 0);
  if (in > 0)
    printf ("%s is open (input) binary format\n", str);
  else
    {
      printf ("could not open %s\n", str);
      exit (1);
    };
  siz = read (in, (char *) crmat, sizeof (crmat));
  printf ("read %i bytes into crmat\n", siz);
  close (in);
#endif

  /* invert the rotation part of the matrix */
  /* seems inverted matices are the same as original matrix */
  /* why is that? Thus, it is unnecessary to invert the */
  /* matrices in the code below */



  /* read until we drop */

  while (1)
    {

    head:
      if (evcount < NPRINT)
        printf ("\n-------------\n");

      /* get the next event, read until we see good data */

      nret = 1;
      while (nret != 0)
        {

          if (nret == NOMOREDATA)
            break;

#if(USGEANTFORMAT==1)

#if(ASCIIDATA==1)

#if(NEWFORMAT==1 || NEWFORMAT==0)
          nret = get_event_GT_ascii (fp, &nhits, &egamma, &evno, id, eg, xx, yy, zz);
#endif

#if(NEWFORMAT==2)
          /* in-beam data */
          nret = get_event_GT_ascii2 (fp, &nhits, &egamma, &evno, id, eg, xx, yy, zz);
#endif

#if(NEWFORMAT==3)
          /* source data */
          nret = get_event_GT_ascii3 (fp, &nhits, &egamma, &evno, id, eg, xx, yy, zz);
#endif

#endif
#if(ASCIIDATA==0)
          nret = get_event_GT_bin (fp, &nhits, &egamma, &evno, id, eg, xx, yy, zz, &TS);
#endif


#else

          nret = get_event_AGT (fp, &nhits, &egamma, &evno, id, eg, xx, yy, zz);
//      printf ("main: nretemin=%i\n", nret);
//      fflush (stdout);
#endif

        };

//for (i=0;i<nhits;i++)
//  assert(id[i]>0);

#if(0)
#if(USGEANTFORMAT==1)
      /* modify id from LR G4 code */
      /* this realy needs to be understood!!!! */

      for (i = 0; i < nhits; i++)
        id[i] -= 4;
#endif
#endif

      if (evcount < NPRINT)
        printf ("===hit, egamma=%f, event number= %i\n", egamma, evno);

      readcount++;
      readcount=evno; /* <-- this is the way it should be?? tl/06/11/2018 */

      /* next timestamp */

      if (readcount % evmod == 0)
        {
          ts += 500;
          if (evcount < NPRINT)
            printf ("incremented timestamp to %lli\n", ts);
        }
      else if (evcount < NPRINT)
        printf ("did not increment timestamp\n");

      if (evcount < NPRINT)
        printf ("TS increment: %i (== readcount %% evmod)\n", readcount % evmod);

//      if (evcount >= 2)
//        exit (0);

      /* are we done? */

      if (nret == NOMOREDATA || evcount > maxevents)
        {
          printf ("header read error at event no %i/%i, nret=%i\n", evcount, maxevents, nret);

          printf ("...in total\n");
          printf ("made %i attempts to read events from the data file [readcount]\n", readcount);
          printf ("produced  %5i interaction points [producedNumIntPoints]\n", producedNumIntPoints);
          printf ("read      %5i interaction points [readNumIntPoints]\n", readNumIntPoints);
          r1 = (float) producedNumIntPoints / (float) readNumIntPoints;
          printf ("compacting reduction %9.2f%%\n", 100 * r1);
          printf ("processed %5i gamma rays [trueNomGammas]\n", trueNomGammas);
          printf ("nphoto  =%10.0f\n", (float) nphoto);
          printf ("ncompton=%10.0f\n", (float) ncompton);
          printf ("total   =%10.0f\n", (float) (ncompton + nphoto));
          r1 = nphoto / (nphoto + ncompton);
          printf ("true P/T= %4.2f\n", r1);
          printf ("but not all gamma rays emitted are full energy, so...\n");
          fflush (stdout);

          /* write out spectra */

          i1 = SPLEN;
//          wr_spe ("G4toMode2_obs.spe", &i1, sp_obs);
//         printf ("wrote \"G4toMode2_obs.spe\"\n");
          wr_spe ("G4toMode2_emit.spe", &i1, sp_emit);
          printf ("wrote \"G4toMode2_emit.spe\"\n");
//          wr_spe ("G4toMode2_obs_nosmear.spe", &i1, sp_obs_nosmear);
//          printf ("wrote \"G4toMode2_obs_nosmear.spe\"\n");
          wr_spe ("G4toMode2_emit_nosmear.spe", &i1, sp_emit_nosmear);
          printf ("wrote \"G4toMode2_emit_nosmear.spe\"\n");

          i1 = SPLEN;
          wr_spe ("scatter_l.spe", &i1, scatter_l);
          printf ("wrote \"scatter_l.spe\"\n");

          r1 = (float) ilok / ((float) ilok + (float) ilnotok);
          printf ("fraction within %9.2f is %9.3f%%\n", maxdist, 100 * r1);


          i1 = 0;
          for (i = 0; i < 20; i++)
            i1 += ndet_sp[i];

          r2 = 0;
          for (i = 0; i < 20; i++)
            if (ndet_sp[i] > 0)
              {
                r1 = 100.0 * (float) ndet_sp[i] / (float) i1;
                printf ("ndet == %2i; %10i, %9.2f %%\n", i, ndet_sp[i], r1);
                r2 += i * r1 / 100;
              };
          printf ("mean number of interactions %9.2f MAX_INTPTS=%i\n", r2, MAX_INTPTS);

#if(0)
          for (i = 0; i < 100; i++)
            printf ("%f mm, %f\n", (float) i / NBIN, sppd[i]);
          i1 = SPLEN;
          wr_spe ("intdist.spe", &i1, sppd);
#endif
          printf ("\naccounting of compacted points rejected on emin\n");
          i1 = eminok + eminnotok;
          r1 = 100.0 * eminok / i1;
          printf ("eminok   =%12i, %7.3f%%\n", eminok, r1);
          r1 = 100.0 * eminnotok / i1;
          printf ("eminnotok=%12i, %7.3f%%\n", eminnotok, r1);

          /* keep track of hitrate in CC */

          i1 = 0;
          for (i = 0; i <= MAXDETNO; i++)
            {
//            printf("-- %3i, %f\n",i,CChitrate[i]);
              if (CChitrate[i] > 0)
                i1++;
            };
          printf ("\n");
          printf ("we hit %i detectors with %i events [INPUT, NOT COUNTED]\n", i1, numG4Events);
          r2 = 0;
          for (k = 0; k <= MAXDETNO; k++)
            if (CChitrate[k] > 0)
              {
                r1 = (float) CChitrate[k] / (float) numG4Events;
                printf ("%3i: %10i, CChitrate=%f, ", k, CChitrate[k], r1 * MAXDETNO / i1);
                printf ("apparent detector success: %f\n", r1 * MAXDETNO / i1);
                r2 += r1;
              };
          r2 /= i1;

          printf ("CC average total efficiency per detector....: %f\n", r2);
          if (i1 > 0)
            printf ("apparent geometric efficiency per detector..: %f\n", r2 * MAXDETNO / i1);
          else
            printf ("apparent geometric efficiency per detector..: %i\n", 0);
          printf ("total CC apparent effienciency..............: %f\n", r2 * i1);


          printf ("\n");
          printf ("mean detector positions\n");
          printf ("\n");

          for (i = 0; i < 200; i++)
            if (nn_mean[i] > 0)
              {
                xx_mean[i] /= nn_mean[i];
                yy_mean[i] /= nn_mean[i];
                printf ("det %3i has mean x and y of %9.2f, %9.2f\n", i, xx_mean[i], yy_mean[i]);
              };

          for (i = 0; i < 200; i++)
            if (nn_mean[i] > 0)
              {
                i1 = 0;
                r1 = 100000;
                for (j = 0; j < 120; j++)
                  {
                    r2 = xx_mean[i] / 10 - crmat[j / 4][j % 4][0][3] * xx_mean[i] / 10 - crmat[j / 4][j % 4][0][3]
                      + yy_mean[i] / 10 - crmat[j / 4][j % 4][1][3] * yy_mean[i] / 10 - crmat[j / 4][j % 4][1][3];
                    r2 = sqrtf (r2);
                    if (r2 < r1)
                      {
                        i1 = j;
                        r1 = r2;
                      };
                  };
                printf ("det %3i nearest GT %3i (%3i), %9.2f\n", i, i1, i - i1, r1);
              };

          /* id hit pattern */

          printf ("\n");
          for (i = 0; i < MAXDETNO; i++)
            if (id_hit[i] > 0)
              printf ("crystal %3i: %10i hits\n", i, id_hit[i]);
          printf ("\n");


          /* full stop */

          printf ("ngammas=%i, evno=%i\n", ngammas, evno);

printf ("etot= %f, mean absorted energy = %f\n", etot, etot/(evno+1));

          printf ("\n%s is done\n\n", argv[0]);
          exit (0);

        };


      if (nhits == 0)
        {
          if (evcount < NPRINT)
            printf ("skip as there were no hits\n");
        }
      else
        {

          /* open detailed listing file */

          if (evcount < NPRINT)
            {

              /* open event file */

              sprintf (str, "G4toMode2_%3.3i.txt", evcount + 1);
              fp1 = fopen (str, "w");
              if (fp1 != NULL)
                printf ("%s is open\n", str);
              else
                {
                  printf ("could not open %s\n", str);
                  exit (1);
                };
              fprintf (fp1, "----\n----file: %s\n\n", str);
#if(0)
              for (i = 0; i < nhits; i++)
                printf ("--++ %i\n", id[i]);
              if (1)
                exit (0);
#endif

            };

          if (fp1 != NULL)
            {
              printf ("nhits=%3i, e=%9.2f, evno= %i\n", nhits, egamma, evno);
              fprintf (fp1, "[1]nhits=%3i, e=%9.2f, evno= %i\n", nhits, egamma, evno);
              esum = 0;
              for (i = 0; i < nhits; i++)
                {
                  if (i == 0)
                    fprintf (fp1, "*");
                  else
                    fprintf (fp1, " ");
                  fprintf (fp1, "read[%2i]: id=%i e=%9.2f (%8.1f %8.1f %8.1f)\n", i, id[i], eg[i], xx[i], yy[i], zz[i]);
                };
            }


          if (nhits > 0)
            {

              /* find the dopler correction factor base on */
              /* first interaction point reported */

#if(1)
              j = 0;
#else
              emax = 0;
              j = 0;
              for (k = 0; k < nhits; k++)
                if (eg[k] > emax)
                  {
                    emax = eg[k];
                    j = k;
                  };
#endif
              radius = xx[j] * xx[j] + yy[j] * yy[j] + zz[j] * zz[j];
              radius = sqrtf (radius);

              dp = (zz[j]) / radius;

              if (dp < -1.0)
                dp = -1.0;
              if (dp > 1.0)
                dp = 1.0;
              polang = acosf (dp);

              radius = 1.0 - vc * vc;
              doppler_factor = sqrtf (radius) / (1.0 - vc * cosf (polang));

              eOK = 1;
              if (fp1 != NULL)
                {
                  fprintf (fp1, "true doppler factor= %9.6f, %9.2f-->%9.2f, eOK=%i\n", doppler_factor, egamma,
                           egamma / doppler_factor, eOK);
                };

              /* keep track of how many interaction points we have encountered */
              /* and how many gamma rays we have seen */

              readNumIntPoints += nhits;
              trueNomGammas += 1;

              /* resort according to detector number, */
              /* which is what I need to simulate decomposed data */


              for (i = 0; i < nhits; i++)
                for (j = i + 1; j < nhits; j++)
                  if (id[i] > id[j])
                    {

                      /*swap them */

                      i1 = id[i];
                      id[i] = id[j];
                      id[j] = i1;
                      r1 = eg[i];
                      eg[i] = eg[j];
                      eg[j] = r1;
                      r1 = xx[i];
                      xx[i] = xx[j];
                      xx[j] = r1;
                      r1 = yy[i];
                      yy[i] = yy[j];
                      yy[j] = r1;
                      r1 = zz[i];
                      zz[i] = zz[j];
                      zz[j] = r1;

                    };

#if(1)
              /* check radius */

              for (i = 0; i < nhits; i++)
                {
                  r1 = xx[i] * xx[i];
                  r1 += yy[i] * yy[i];
                  r1 += zz[i] * zz[i];
                  r1 = sqrtf (r1);
//                assert (r1 < 290.0);
                  if (r1 > 350.0)
                    {
                      printf ("ooops, radius= %9.2f", r1);
                      printf ("read[%2i]: id=%i e=%9.2f ", i, id[i], eg[i]);
                      printf ("x,y,z=(%8.1f %8.1f %8.1f) ", xx[i], yy[i], zz[i]);
                      printf ("\n");
                    };

                };
#endif

              /* debug print */

              if (fp1 != NULL && r1 > 350.0)
                {
                  esum = 0;
                  fprintf (fp1, "after reordering, r1=%f\n", r1);
                  for (i = 0; i < nhits; i++)
                    {
                      fprintf (fp1, "read[%2i]: id=%i e=%9.2f x,y,z=(%8.1f %8.1f %8.1f), ", i, id[i], eg[i], xx[i],
                               yy[i], zz[i]);
                      esum += eg[i];
                      fprintf (fp1, "esum=%9.2f; ", esum);
                      r1 = xx[i] * xx[i];
                      r1 += yy[i] * yy[i];
                      r1 += zz[i] * zz[i];
                      r1 = sqrtf (r1);
                      if (i > 0)
                        {
                          r1 = (xx[i] - xx[i - 1]);
                          r2 = (yy[i] - yy[i - 1]);
                          r3 = (zz[i] - zz[i - 1]);
                          r1 = sqrtf (r1 * r1 + r2 * r2 + r3 * r3);
                          fprintf (fp1, "raw scatter len= %9.2f", r1);
                        };
                      fprintf (fp1, "\n");
                    };

                };

              /* keep record of photo to compton */

              esum = 0;
              for (i = 0; i < nhits; i++)
                esum += eg[i];
              if ((egamma - esum) > 3.0)
                ncompton++;
              else
                nphoto++;

              /* bin the observed energy */

              i1 = (int) (esum / doppler_factor);
              if (i1 >= 0 && i1 < SPLEN)
                {
                  sp_obs_nosmear[i1]++;
                };


              /* bin the emitted energy */

              i1 = (int) (egamma / doppler_factor);
              if (i1 >= 0 && i1 < SPLEN)
                {
                  sp_emit_nosmear[i1]++;
                };

              /* adapted AGATA energy smearing */

              r1 = sqrtf (1.0 + 3.7 * (egamma / 1000.0));
              r1 /= 2.3548;
              if (fp1 != NULL)
                fprintf (fp1, "AGATA e smearing val=%f\n", r1);
              r1 *= ranGauss ();
              if (fp1 != NULL)
                fprintf (fp1, "smearing egamma=%f with %f\n", egamma, (eSmear_seg * r1));
              egamma += (eSmear_seg * r1);
              i1 = (int) (egamma / doppler_factor);
              if (i1 >= 0 && i1 < SPLEN)
                {
                  sp_emit[i1]++;
//                  printf("%i %f\n", i1,sp_emit[i1]);
                };


              /* find new detector points in the list */

              ndet = 0;
              i1 = -1;
              for (i = 0; i < nhits; i++)
                if (id[i] != i1)
                  {

                    /* new detector in play */

                    pos[ndet] = i;
                    ndet++;
                    i1 = id[i];

                  };
              pos[ndet] = nhits;

              /* only process if within requested limits */

              if (ndet >= minNoCrystals && ndet <= maxNoCrystals)
                {



                  if (fp1 != NULL)
                    {
                      fprintf (fp1, "\nwe have %i detectors in play\n", ndet);

                      k = 0;
                      for (i = 0; i < ndet; i++)
                        {
                          fprintf (fp1, "------");
                          fprintf (fp1, "range %i to %i\n", pos[i], pos[i + 1] - 1);
                          esum = 0;
                          for (j = pos[i]; j < pos[i + 1]; j++)
                            {
                              fprintf (fp1, "%2i: id=%i e=%9.2f (%8.1f,%8.1f,%8.1f) ", k, id[j], eg[j], xx[j], yy[j],
                                       zz[j]);
                              k++;
                              esum += eg[j];
                              fprintf (fp1, "esum=%9.2f; ", esum);
                              fprintf (fp1, "\n");
                            }
                        };
                      fprintf (fp1, "\n");
                    };


                  for (i = 0; i < ndet; i++)
                    {

                      /* limits for hist in this detector */

                      i1 = pos[i];
                      i2 = pos[i + 1];

                      if (fp1 != NULL)
                        {
                          printf ("det # %i; ", i);
                          printf ("range %i to %i\n", pos[i], pos[i + 1] - 1);
                        };

                      /* reset assignments */

                      for (j = i1; j < i2; j++)
                        hitnum[j] = -1;

                     /*----------------------------------------------*/
                      /* repack is a packing pameter > 0 is specified */
                     /*----------------------------------------------*/

                      /* packing using the philosophy of AGATA */

                      /* prime */

                      for (j = i1; j < i2; j++)
                        {
                          valid[j] = 1;
                          realx[j] = xx[j];
                          realy[j] = yy[j];
                          realz[j] = zz[j];
                        };


                      for (nit = 0; nit < 3; nit++)
                        for (j = i1; j < i2; j++)
                          for (k = i1; k < i2; k++)
                            if (j != k)
                              if (valid[j])
                                if (valid[k])
                                  {

                                    /* find relative distance */

                                    if (fp1 != NULL)
                                      {
                                        fprintf (fp1, "(%2i) %5.2f %5.2f %5.2f ... ", j, realx[j], realy[j], realz[j]);
                                        fprintf (fp1, "(%2i) %5.2f %5.2f %5.2f ", k, realx[k], realy[k], realz[k]);
                                      };
                                    dr = (realx[j] - realx[k]) * (realx[j] - realx[k])
                                      + (realy[j] - realy[k]) * (realy[j] - realy[k])
                                      + (realz[j] - realz[k]) * (realz[j] - realz[k]);
                                    dr = sqrtf (dr);
                                    if (fp1 != NULL)
                                      fprintf (fp1, "dr= %5.2f\n", dr);

                                    /*--------------------------------------------------*/
                                    /* two points are too close and need to be combined */
                                    /*--------------------------------------------------*/

                                    if (dr < minr && minr > 0)
                                      {

                                        if (fp1 != NULL)
                                          fprintf (fp1, "pack2: combine %i and %i dr=%f < minr=%f\n", j, k, dr, minr);

                                        /* invalidate j */

                                        valid[j] = 0;

                                        /* let k be the energy weighted position */

                                        realx[k] = (realx[k] * eg[k] + realx[j] * eg[j]) / (eg[k] + eg[j]);
                                        realy[k] = (realy[k] * eg[k] + realy[j] * eg[j]) / (eg[k] + eg[j]);
                                        realz[k] = (realz[k] * eg[k] + realz[j] * eg[j]) / (eg[k] + eg[j]);

                                        /* k now has the summed energy */

                                        eg[k] = eg[k] + eg[j];

                                      }

                                  };

                      numclusters = 0;
                      for (j = i1; j < i2; j++)
                        if (valid[j])
                          {
                            if (eg[j] > emin)
                              {
                                numclusters++;
                                xx[j] = realx[j];
                                yy[j] = realy[j];
                                zz[j] = realz[j];
                              }
                            else
                              {
                                valid[j] = 0;
                                if (fp1 != NULL)
                                  fprintf (fp1, "%2i: invalidated, e=%9.1f < %9.1f\n", j, eg[j], emin);
                              };
                          };
//                     printf("mmm_numclusters=%i\n",numclusters);

                      if (fp1 != NULL)
                        {
                          fprintf (fp1, "result of packing:\n");
                          for (j = i1; j < i2; j++)
                            {
                              if (valid[j])
                                fprintf (fp1, "%2i: (%5.2f %5.2f %5.2f) e=%9.3f\n", j, xx[j], yy[j], zz[j], eg[j]);
                              else
                                fprintf (fp1, "%2i: invalid (was packed)\n", j);
                            };
                        };



                      if (fp1 != NULL)
                        {
                          /* debug list what we have */

                          fprintf (fp1, "there are %2i clusters for detector %2i\n", numclusters, i);
                        };

                      /* write them out */

//                      assert (numclusters > 0);

                      if (numclusters > 0)
                        {
                          if (fp1 != NULL)
                            {
                              fprintf (fp1, "writing these %i interactions out\n", numclusters);
                              for (j = i1; j < i2; j++)
                                if (valid[j])
                                  fprintf (fp1, "%2i: (%5.2f %5.2f %5.2f) e=%9.3f\n", j, xx[j], yy[j], zz[j], eg[j]);
                            }

                          /* prepare geb header */

                          bzero ((void *) payload, MAXPAYLOADSIZE);
                          gd->length = numclusters * sizeof (CRYS_INTPTS);
                          if (gd->length >= MAXPAYLOADSIZE)
                            {
                              printf ("gd->length=%i is greater than MAXPAYLOADSIZE=%i; QUIT\n", gd->length,
                                      MAXPAYLOADSIZE);
                              exit (1);
                            };
                          gd->type = GEB_TYPE_DECOMP;
                          payload->type = (int) 0xabcd5678;
#if (USGEANTFORMAT==1 &&  ASCIIDATA==0)
                          gd->timestamp = TS;
                          payload->timestamp = TS;
#else
                          gd->timestamp = ts;
                          payload->timestamp = ts;
#endif

                          /* keep track of hitrate in CC */

                          CChitrate[id[i1]]++;
//                          printf("CChitrate[%i]=%i\n",id[i1],CChitrate[id[i1]]);

                          /* hide G4 event number in pad upper bits */

                          if (fp1 != NULL)
                            fprintf (fp1, "[2]event number %i, pad=%i, %i\n", evno, payload->pad);
                          payload->pad |= (evno << 8);
                          if (fp1 != NULL)
                            fprintf (fp1, "new pad=%i\n", evno, payload->pad);


                          payload->num = 0;
                          for (l = i1; l < i2; l++)
                            if (valid[l])
                              {
//                                id[l]+=1;
                                holeNum = id[l] / 4;

                                id_hit[id[l]]++;

                                crystalNumber = id[l] % 4;
//                                printf ("aa: l=%i,id=%i holeNum=%i crystalNumber=%i\n", l, id[l],
//                                        holeNum, crystalNumber);

                                payload->crystal_id = crystalNumber;
                                payload->crystal_id |= (holeNum << 2);



                                if (fp1 != NULL)
                                  {
                                    fprintf (fp1, "id[l]=%i, ", id[l]);
                                    fprintf (fp1, "holeNum=%i, ", holeNum);
                                    fprintf (fp1, "crystalNumber=%i, ", crystalNumber);
                                    fprintf (fp1, "payload->crystal_id=0x%x,%i\n", payload->crystal_id,
                                             payload->crystal_id);
                                  }

                                payload->intpts[payload->num].x = xx[l];
                                payload->intpts[payload->num].y = yy[l];
                                payload->intpts[payload->num].z = zz[l];
                                payload->intpts[payload->num].e = eg[l];
                                payload->tot_e += payload->intpts[payload->num].e;
                                payload->num++;
//                                assert (payload->num < MAX_INTPTS);


                              };

                          /* add some gaussian width to energies and postitions */
                          /* CC */

                          r1 = sqrtf (1.0 + 3.7 * (payload->tot_e / 1000.0));
                          r1 /= 2.3548;
                          if (fp1 != NULL)
                            fprintf (fp1, "== tot_e %f (%f %f)--> ", payload->tot_e, r1, (eSmear_CC * r1));
                          r1 *= ranGauss ();
                          payload->tot_e += (eSmear_CC * r1);
                          if (fp1 != NULL)
                            fprintf (fp1, "%f\n", payload->tot_e);

                          /* add some gaussian width to segment energies and postitions */
                          /* segments */

                          for (j = 0; j < payload->num; j++)
                            {
                              if (payload->intpts[j].e > 10)
                                {

                                  /* Simulate energy resolution */

                                  if (enRes_0[payload->crystal_id] > 0)
                                    {
                                      /* Use crystal resolution parameters, if defined in the chat file ... */
                                      r1 =
                                        enRes_0[payload->crystal_id] * sqrtf (1.0 +
                                                                              payload->intpts[j].e *
                                                                              enRes_1[payload->crystal_id]);
                                      if (fp1 != NULL)
                                        fprintf (fp1, "Crystal %i resolution=%f\n", payload->crystal_id, r1);
                                    }
                                  else
                                    {

                                      /* ... otherwise: adapted AGATA energy smearing */

                                      r1 = sqrtf (1.0 + 3.7 * (payload->intpts[j].e / 1000.0));
                                      r1 /= 2.3548;
                                      if (fp1 != NULL)
                                        fprintf (fp1, "AGATA e smearing val=%f\n", r1);
                                    }

                                  /* scale by eSmear_seg factor from the chat file */

                                  if (fp1 != NULL)
                                    fprintf (fp1, "== seg_e %f (%f %f)--> ", payload->intpts[j].e, r1,
                                             (eSmear_seg * r1));
                                  r1 *= ranGauss ();
                                  payload->intpts[j].e += (eSmear_seg * r1);
                                  if (fp1 != NULL)
                                    fprintf (fp1, "%f\n", payload->intpts[j].e);
                                  esum_smeared += payload->intpts[j].e;

                                  /* adapted AGATA position smearing */

                                  r1 = 0.5 * sqrtf (0.1 / (payload->intpts[j].e / 1000.0));
                                  r1 /= 2.3548;
                                  r1 *= pSmear;
                                  payload->intpts[j].x += (r1 * ranGauss ());
                                  payload->intpts[j].y += (r1 * ranGauss ());
                                  payload->intpts[j].z += (r1 * ranGauss ());

                                };
                            };

                          if (noworldtocrystalrot == 0)
                            for (j = 0; j < payload->num; j++)
                              {

                                if (fp1 != NULL)
                                  {
                                    fprintf (fp1, "* holeNum=%i, crystalNumber=%i\n", holeNum, crystalNumber);
                                    fprintf (fp1, "to crystal * %i: ", j);
                                    fprintf (fp1, "%7.2f,%7.2f,%7.2f --> \n", payload->intpts[j].x,
                                             payload->intpts[j].y, payload->intpts[j].z);
                                  }


                                /* positions need to be in cm for this version of crmat */

                                x1 = payload->intpts[j].x / 10;
                                y1 = payload->intpts[j].y / 10;
                                z1 = payload->intpts[j].z / 10;
                                if (fp1 != NULL)
                                  fprintf (fp1, "1- %f %f %f\n", x1, y1, z1);

                                /* first subtract the translation */

                                x1 -= crmat[holeNum][crystalNumber][0][3];
                                y1 -= crmat[holeNum][crystalNumber][1][3];
                                z1 -= crmat[holeNum][crystalNumber][2][3];
                                if (fp1 != NULL)
                                  fprintf (fp1, "2- %f %f %f\n",
                                           crmat[holeNum][crystalNumber][0][3],
                                           crmat[holeNum][crystalNumber][1][3], crmat[holeNum][crystalNumber][2][3]);
                                if (fp1 != NULL)
                                  fprintf (fp1, "3- %f %f %f\n", x1, y1, z1);

                                /* then invert (transpose) the rotation */

                                payload->intpts[j].x = crmat[holeNum][crystalNumber][0][0] * x1
                                  + crmat[holeNum][crystalNumber][1][0] * y1 + crmat[holeNum][crystalNumber][2][0] * z1;

                                payload->intpts[j].y = crmat[holeNum][crystalNumber][0][1] * x1
                                  + crmat[holeNum][crystalNumber][1][1] * y1 + crmat[holeNum][crystalNumber][2][1] * z1;

                                payload->intpts[j].z = crmat[holeNum][crystalNumber][0][2] * x1
                                  + crmat[holeNum][crystalNumber][1][2] * y1 + crmat[holeNum][crystalNumber][2][2] * z1;
                                if (fp1 != NULL)
                                  fprintf (fp1, "4- %f %f %f\n", x1, y1, z1);


                                /* make it mm again */

                                payload->intpts[j].x *= 10;
                                payload->intpts[j].y *= 10;
                                payload->intpts[j].z *= 10;
//assert(payload->intpts[j].x<45);

                                if (fp1 != NULL)
                                  fprintf (fp1, "%7.2f,%7.2f,%7.2f\n", payload->intpts[j].x, payload->intpts[j].y,
                                           payload->intpts[j].z);

                              };


                          /* write out if enabled */

                          if (enabled[payload->crystal_id])
                            {
                              siz = write (outData, (char *) gd, sizeof (GEBDATA));
                              assert (siz == sizeof (GEBDATA));
                              siz = write (outData, (char *) payload, gd->length);
                              assert (siz == gd->length);
                              nPoints += payload->num;
                            };

                          /* keep record of how many packed */
                          /* interaction points we write out */

                          producedNumIntPoints += payload->num;


                          if (fp1 != NULL)
                            {
                              fprintf (fp1, " OK: payload->tot_e=%f emin=%f, pad=%i\n", payload->tot_e, emin,
                                       payload->pad);
                              fprintf (fp1, "^^^^^wrote %i interaction points out for this detector, ts=%lli\n",
                                       payload->num, gd->timestamp);
                              fprintf (fp1, "in crystal coordinates\n");
                            };
                          if (fp1 != NULL)
                            {
                              fprintf (fp1, "evcount=%3i, wrote:\n", evcount);
                              printCRYS_INTPTS (fp1, payload, gd);
                              fprintf (fp1, "...so far we have\n");
                              fprintf (fp1, "made %i attempts to read events from the data file [readcount]\n",
                                       readcount);
                              fprintf (fp1, "read      %7i interaction points [readNumIntPoints]\n", readNumIntPoints);
                              fprintf (fp1, "produced  %7i interaction points [producedNumIntPoints]\n",
                                       producedNumIntPoints);
                              fprintf (fp1, "processed %7i gamma rays [trueNomGammas]\n", trueNomGammas);
                            };


                        };
                    };

                };



//              if (evcount >= 7) exit (0);

              /* output the clusters */


              evcount++;


#if(0)
              if (evcount >= 7)
                if (1)
                  exit (0);
#endif

              if (fp1 != NULL)
                {
                  fprintf (fp1, "\nend of interactions\n");
                  fclose (fp1);
                  fp1 = NULL;
                };

#if(0)
              /* stop for debugging */

              if ((int) ts > 33500)
                if (1)
                  exit (0);
#endif
            };
        };

    };



  /* done */

  return (0);
};
