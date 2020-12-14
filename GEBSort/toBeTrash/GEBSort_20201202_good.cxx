
/* sorter for data from the GEB tab (or a file) */
/* totally universal in its scope */
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <netinet/in.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <signal.h>
#include <time.h>
#include <stddef.h>
#include <zlib.h>
//#include <values.h>

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TFile.h"
#include "TCutG.h"
#include "TTree.h"
#include "TMapFile.h"
#include "TSystem.h"

#include "gdecomp.h"
#include "veto_pos.h"
#include "GEBSort.h"
#include "GTMerge.h"

#if(HAVE_VXWORKS)
#include <GEBClient.h>
#include <gretTapClient.h>
#endif

/* global variables */

FILE *inData;

/* buffer of data for sorting */

char *rbuf;

extern TH1D *sumTrackE;

int nn2 = 0;
int nn3 = 0;

//NprintEvNo = 0;
int bufPos = 0;
int nBadTestPat = 0;
int egemin = 2;
double oldLEDTs = 0;
double oldCFDTs = 0;
time_t tdmp = 0, tdmplast;
float ehiGeOffset[NCHANNELS];
float ehiGeGain[NCHANNELS];
int ehiDoGainCor = 0;

//char CommandFileName[STRLEN] = "GEBSort.command";

int time_stamp (FILE *);
TH1D *mkTH1D (char *, char *, int, double, double);
TH2F *mkTH2F (char *, char *, int, double, double, int, double, double);

TH1D *ehi[MAXDETPOS + 1];

/* user */

#define DEBUG1 0
#define DEBUG2 0

/* Doppler correction and Ge calibrations */

#define NGSGE 110

float cal_off[NGSGE], cal_gain[NGSGE];
double angle[NGSGE], anglephi[NGSGE];

double DopCorFac[NGSGE], ACFac[NGSGE];

void SetBeta ();

/* parameters */

//TFile *treef;
//TTree *tree;

PARS Pars;
EXCHANGE exchange;

/* common storage of bin_xxx analysis */

int ng;
DGSEVENT DGSEvent[MAXCOINEV];
DFMAEVENT DFMAEvent[MAXCOINEV];
int XAng;
DGSEVENT XAEvent[MAXCOINEV];

static const int ConnectionRetryWait = 2 * 100; /* usec between con. retries. */
static const int ConnectionRetryCount = 150;    /* Times we retry connect. */

struct gretTap *pTap;
struct GEBData *pData;

int nn1 = 0;

unsigned int *veto_cube;

/*----------------------------------------------------*/
int
fleft (FILE * fp)
{
  int prev = ftell (fp);
  fseek (fp, 0L, SEEK_END);
  int sz = ftell (fp);
  fseek (fp, prev, SEEK_SET);   //go back to where we were
  return (sz - prev);
}

/*----------------------------------------------------*/

TH1D *
mkTH1D (char *str1, char *str2, int nn, double lo, double hi)
{
  TH1D *tmppt;

  if (!Pars.UpdateRootFile)
    {
      tmppt = new TH1D (str1, str2, nn, lo, hi);
      printf ("Created Object \"%s\", %p\n, \"%s\"", str1, tmppt, str2);
    }
  else
    {
      tmppt = (TH1D *) gROOT->FindObject (str1);
      printf ("Found Object \"%s\", %p\n", str1, tmppt);
    }

  return (tmppt);

}

/*----------------------------------------------------*/

TH2F *
mkTH2F (char *str1, char *str2, int nn1, double lo1, double hi1, int nn2, double lo2, double hi2)
{
  TH2F *tmppt;

  if (!Pars.UpdateRootFile)
    {
      tmppt = new TH2F (str1, str2, nn1, lo1, hi1, nn2, lo2, hi2);
      printf ("Created Object \"%s\", %p\n", str1, tmppt);
    }
  else
    {
      tmppt = (TH2F *) gROOT->FindObject (str1);
      printf ("Found Object \"%s\", %p\n", str1, tmppt);
    };

  return (tmppt);

}

/*----------------------------------------------------*/

TH3F *
mkTH3F (char *str1, char *str2,
        int nn1, double lo1, double hi1, int nn2, double lo2, double hi2, int nn3, double lo3, double hi3)
{
  TH3F *tmppt;

  if (!Pars.UpdateRootFile)
    {
      tmppt = new TH3F (str1, str2, nn1, lo1, hi1, nn2, lo2, hi2, nn3, lo3, hi3);
      printf ("Created Object \"%s\", %p\n", str1, tmppt);
    }
  else
    {
      tmppt = (TH3F *) gROOT->FindObject (str1);
      printf ("Found Object \"%s\", %p\n", str1, tmppt);
    };

  return (tmppt);

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
void
SetBeta ()
{

  /* delarations */

  int i;
  double d1;

  /*-------------------------------------*/
  /* find Doppler and aberration factors */
  /*-------------------------------------*/

  for (i = 0; i < NGSGE; i++)
    {
      //printf("det %3.3i-> ", i);
      d1 = angle[i] / 57.29577951;
      DopCorFac[i] = (1 - Pars.beta * cos (d1)) / sqrt (1 - Pars.beta * Pars.beta);
      //printf("dop cor fac: %6.4f; ", DopCorFac[i]);
      ACFac[i] = DopCorFac[i] * DopCorFac[i];
      //printf("aberration cor fac: %6.4f\n", ACFac[i]);

    };
  fflush (stdout);

  /* done */

}

/*-----------------------------------------------------------*/

float
findPolarFromCartesian (float xx, float yy, float zz, float *rr)
{
  float d1;

  *rr = sqrtf (xx * xx + yy * yy + zz * zz);
  d1 = acosf (zz / *rr);

  return (d1);
}

/*-----------------------------------------------------------*/

float
findAzimuthFromCartesian (float xx, float yy, float zz)
{

  float d1;

#if(0)
  if (xx > 0 && yy >= 0)
    d1 = atanf (yy / xx);
  else if (xx > 0 && yy < 0)
    d1 = atanf (yy / xx) + 2. * M_PI;
  else if (xx < 0)
    d1 = atanf (yy / xx) + M_PI;
  else if (xx == 0 && yy > 0)
    d1 = M_PI / 2.;
  else if (xx == 0 && yy < 0)
    d1 = 3. * M_PI / 2.;
  else
    d1 = 0.0;
#endif
#if(1)
  d1 = atan2f (yy, xx);
//  if (d1<0) d1+=M_PI;
//  if (d1>M_PI) d1-=M_PI;
#endif


  return (d1);

}



/*----------------------------------------------------------------------------*/

int GebTypeStr (int type, char str[]){
//   printf("got type %i, %i\n", type,GEB_TYPE_DECOMP);
  if (type == GEB_TYPE_DECOMP)
    sprintf (str, "GEB_TYPE_DECOMP      ");
  else if (type == GEB_TYPE_BGS)
    sprintf (str, "GEB_TYPE_BGS         ");
  else if (type == GEB_TYPE_RAW)
    sprintf (str, "GEB_TYPE_RAW         ");
  else if (type == GEB_TYPE_TRACK)
    sprintf (str, "GEB_TYPE_TRACK       ");
  else if (type == GEB_TYPE_S800_RAW)
    sprintf (str, "GEB_TYPE_S800_RAW    ");
  else if (type == GEB_TYPE_NSCLnonevent)
    sprintf (str, "GEB_TYPE_NSCLnonevent");
  else if (type == GEB_TYPE_GT_SCALER)
    sprintf (str, "GEB_TYPE_GT_SCALER   ");
  else if (type == GEB_TYPE_GT_MOD29)
    sprintf (str, "GEB_TYPE_GT_MOD29    ");
  else if (type == GEB_TYPE_S800PHYSDATA)
    sprintf (str, "GEB_TYPE_S800PHYSDATA");
  else if (type == GEB_TYPE_G4SIM)
    sprintf (str, "GEB_TYPE_G4SIM       ");
  else if (type == GEB_TYPE_CHICO)
    sprintf (str, "GEB_TYPE_CHICO       ");
  else if (type == GEB_TYPE_NSCLNONEVTS)
    sprintf (str, "GEB_TYPE_NSCLNONEVTS ");
  else if (type == GEB_TYPE_DGS)
    sprintf (str, "GEB_TYPE_DGS         ");
  else if (type == GEB_TYPE_DGSTRIG)
    sprintf (str, "GEB_TYPE_DGSTRIG     ");
  else if (type == GEB_TYPE_DFMA)
    sprintf (str, "GEB_TYPE_DFMA        ");
  else if (type == GEB_TYPE_PHOSWICH)
    sprintf (str, "GEB_TYPE_PHOSWICH    ");
  else if (type == GEB_TYPE_PHOSWICHAUX)
    sprintf (str, "GEB_TYPE_PHOSWICHAUX ");
  else if (type == GEB_TYPE_GODDESS)
    sprintf (str, "GEB_TYPE_GODDESS ");
  else if (type == GEB_TYPE_LABR)
    sprintf (str, "GEB_TYPE_LABR ");
  else if (type == GEB_TYPE_GODDESSAUX)
    sprintf (str, "GEB_TYPE_GODDESSAUX ");
  else if (type == GEB_TYPE_XA)
    sprintf (str, "GEB_TYPE_XA ");
  else
    {
      sprintf (str, "%i is unknown, maybe update 'GebTypeStr' function?", type);
      return (1);
    };
//      printf("type: %s\n",str);

  return (0);

};


/*----------------------------------------------------------------------------*/
int
fbuf_read (FILE * inData, char *buf, int nbytes)
{

  /* declarations */

  int st, i;
  long long int li1;

//    printf("fbuf_read, asked for %i bytes\n", nbytes);

  /* check if there is enough data to read */
  /* wait for data? */

  if (Pars.waitfordataseconds > 0)
    {
      li1 = fleft (inData);
      if (li1 < (long long int) nbytes)
        {
          /* enter wait phase */

          for (i = 0; i < Pars.waitfordataseconds; i++)
            {
              printf (".");
              fflush (stdout);
              sleep (1);
              li1 = fleft (inData);
              if (li1 > 2 * (long long int) nbytes)
                {
                  printf ("\n");
                  fflush (stdout);
                  break;
                };
            };
        };
    };


  /* get the data */

  st = fread (buf, 1, nbytes, inData);

  /* done */

  return (st);

};

/*----------------------------------------------------------------------------*/
int
buf_read (FILE * inData, char *buf, int nbytes)
{

/* buffered read of data from disk or network */
/* 'rbuf' is the buffer we maintain, it is global */
/* 'buf' is the buffer of data we return */

/* declarations */

  ssize_t st;
  int partread, l, ntries;
  struct GEBData *tmpData, *tmpStore;
  static int bpos = RBUFSIZE, ActualBufSize;
  static int bleft = 0;
  char str[64];
  unsigned int *pi4;
  int i, bstart;
  struct stat fst;
  long long int li1;

  if (Pars.CurEvNo <= Pars.NumToPrint)
    printf ("buf_read called, asked for %i byte, have %i bytes left in buffer, bpos=%i\n", nbytes, bleft, bpos);

  if (buf == NULL)
    {
      printf ("buf_read: pointer to buf is NULL, cannot continue\n");
      return (0);
    };

#if(0)

  /*---------------------------------*/
  /* get the data with a simple read */
  /*---------------------------------*/
  /* we never want to do that, too slow */

//  st = read (inData, buf, nbytes);
  st = fread (buf, 1, nbytes, inData);

#else

  /*-----------------------*/
  /* get data via a buffer */
  /*-----------------------*/

//  printf("entered buffer read, asking for %i bytes\n",nbytes);

  /* do we need to read a new buffer in ? */

  if (nbytes > bleft)
    {
      if (Pars.CurEvNo <= Pars.NumToPrint)
        printf ("__forced to read data because nbytes=%i > bleft=%i\n", nbytes, bleft);

      /* first read what is left in the buffer */

      partread = bleft;
      if (bleft > 0)
        {
          memcpy (buf, rbuf + bpos, partread);

          bleft = 0;
          bpos = ActualBufSize;
          if (Pars.CurEvNo <= Pars.NumToPrint)
            printf ("partial transfer of %i bytes, set bleft=%i, set bpos=%i\n", partread, bleft, bpos);
        }
      else
        partread = 0;

      /* then read in a new buffer from disk or GEB */

      if (Pars.InputSrc == DISK)
        {

          /* check if there is enogh left in the file */
          /* to read a buffer. If not, we wait and see */
          /* if it gets filled before we go on */

          if (Pars.waitfordataseconds > 0)
            {

              /* how many bytes left? */

//              fstat (inData, &fst);
//              li1 = fst.st_size - Pars.nbytes;
              li1 = fleft (inData);

              /* wait? */

              if (li1 < (long long int) RBUFSIZE)
                {
                  /* enter wait phase */

                  for (i = 0; i < Pars.waitfordataseconds; i++)
                    {
                      printf (".");
                      fflush (stdout);
                      sleep (1);
//                      fstat (inData, &fst);
//                      li1 = fst.st_size - Pars.nbytes;
                      li1 = fleft (inData);
                      if (li1 > 2 * (long long int) RBUFSIZE)
                        {
                          printf ("\n");
                          fflush (stdout);
                          break;
                        };
                    };
                };
            };


//          ActualBufSize = read (inData, rbuf, RBUFSIZE);
          ActualBufSize = fread (rbuf, 1, RBUFSIZE, inData);
          nn1 += ActualBufSize;
          bpos = 0;
          bleft = ActualBufSize;
          if (Pars.CurEvNo <= Pars.NumToPrint)
            printf ("read new buffer of size %i bytes, bpos=%i, bleft=%i\n", ActualBufSize, bpos, bleft);
          if (ActualBufSize != RBUFSIZE)
            {
              printf ("could not read a full buffer of size %i, only got %i, bleft=%i\n", RBUFSIZE, ActualBufSize,
                      bleft);
              printf ("read new buffer of size %i bytes, bpos=%i, bleft=%i\n", ActualBufSize, bpos, bleft);
              if (ActualBufSize <= 0)
                return (0);
            };

        }
      else if (Pars.InputSrc == GEB)
        {
#if (HAVE_VXWORKS)
          bzero (rbuf, RBUFSIZE + 1);

          if (Pars.CurEvNo <= Pars.NumToPrint)
            {
              printf ("attempting to get %i gebheader/payloads from GEB\n", Pars.grouping);
              fflush (stdout);
            }

          ntries = 0;
          pData = gretTapData (pTap, Pars.grouping, Pars.timeout);

//          printf("pData = %p\n", pData);fflush(stdout);

          /* failed to get data */

          if (pData == NULL)
            {
              printf ("failed to get more data from GEBTap, return 0 bytes\n");
              return (0);
            }



          if (Pars.CurEvNo <= Pars.NumToPrint)
            printf ("... got them, we think\n");

          /* what comes back is a linked list of events */
          /* which we quietly read and stuff into the  */
          /* buffer. It would probably be better to use  */
          /* the list directly, but the buffer method  */
          /* exists and we are not as clever as Carl, by far */

          /* FYI: struct GEBData *pData; GEBData defined in GEBLink.h */
          /*      dont know what to do with refCount and refIndex yet */


          bpos = 0;
          bleft = st;
          tmpData = pData;
          bstart = 0;
          for (l = 0; l < Pars.grouping; l++)
            {

//              printf ("***type=%i\n", tmpData->type);
//              fflush (stdout);

              /* we do not take RAW data nor bloated payloads */

//              printf("l=%i, %i, %i\n", l, tmpData->type,tmpData->length);
//              fflush(stdout);
//              if (tmpData->type != GEB_TYPE_RAW && tmpData->length<MAXPAYLOADSIZE)

              /* crash and burn if the payload is too big */
              /* better to know than to proceed */

              if (tmpData->length >= MAXPAYLOADSIZE)
                {
                  printf ("\n");
                  printf ("PROBLEM in function buf_read:\n");
                  printf ("tmpData->length is %i\n", tmpData->length);
                  printf ("that is larger than MAXPAYLOADSIZE=%i\n", MAXPAYLOADSIZE);
                  printf ("quitting\n");
                  printf ("\n");
                  exit (1);
                };

//              if (tmpData->length < MAXPAYLOADSIZE)
              if (tmpData->payload != NULL)
                {

                  bstart = bpos;
                  if (Pars.CurEvNo <= Pars.NumToPrint)
                    {
                      printf ("processing link # %i\n", l);
                      fflush (stdout);
                    }

                  /* copy GEB header */

                  if (Pars.CurEvNo <= Pars.NumToPrint)
                    {
                      printf ("will copy header; len %i, ", GEB_HEADER_BYTES);
                      printf ("TS=%lli;", tmpData->timestamp);
                      printf ("type=%i", tmpData->type);
                      GebTypeStr (tmpData->type, str);
                      printf ("==%s", str);
                      printf ("to bpos=%i(Byte) ", bpos);
                      printf ("... ");
                    };
                  memcpy (rbuf + bpos, (char *) tmpData, GEB_HEADER_BYTES);
                  bpos += GEB_HEADER_BYTES;
                  assert (bpos < RBUFSIZE);
                  if (Pars.CurEvNo <= Pars.NumToPrint)
                    {
                      printf ("done\n");
                    };

                  /* copy payload */

#if(0)
                  printf ("tmpData->length=%i, ", tmpData->length);
                  printf ("bpos=%i, ", bpos);
                  printf ("RBUFSIZE=%i, ", RBUFSIZE);
                  printf ("tmpData->payload=0x%p\n", tmpData->payload);
                  fflush (stdout);
#endif
                  memcpy (rbuf + bpos, (char *) tmpData->payload, tmpData->length);
                  bpos += tmpData->length;
                  assert (bpos < RBUFSIZE);

                  if (Pars.CurEvNo <= Pars.NumToPrint)
                    printf ("copied payload of length %i\n", tmpData->length);

                }
              else
                {
                  if (Pars.CurEvNo <= Pars.NumToPrint && 0)
                    {
                      printf ("skipped data of type ");
                      GebTypeStr (tmpData->type, str);
                      printf ("%s", str);
                      printf ("with` ll=ength of %i\n", tmpData->length);
                      fflush (stdout);

                    }
                };

              /* next link and free the previous one */

              tmpStore = tmpData;
              tmpData = tmpData->next;
              free (tmpStore);

            };

          /* done transferring, set counters */
          /* for reading again */

          bleft = bpos;
          bpos = 0;
#endif
        };

      /* then read the rest */

      if (Pars.CurEvNo <= Pars.NumToPrint)
        {
          printf ("buf=0x%p ", buf);
          printf ("rbuf=0x%p ", rbuf);
          printf ("partread:%i ", partread);
          printf ("bpos:%i ", bpos);
          printf ("nbytes:%i ", nbytes);
          printf ("nbytes - partread:%i ", nbytes - partread);
          printf ("bleft:%i\n", bleft);
          fflush (stdout);
        }

      /* trap */

      if (bleft < (nbytes - partread))
        {
          printf ("GEBSort error: the buffer must be define big enough to hold the biggest event we have\n");
          printf ("#define RBUFSIZE %i (at least) in GEBSort.h\n", nbytes);
          printf ("buf=0x%p ", buf);
          printf ("rbuf=0x%p ", rbuf);
          printf ("partread:%i ", partread);
          printf ("bpos:%i ", bpos);
          printf ("nbytes:%i ", nbytes);
          printf ("nbytes - partread:%i ", nbytes - partread);
          printf ("bleft:%i\n", bleft);
          fflush (stdout);
          exit (1);
        };

      memcpy (buf + partread, rbuf + bpos, nbytes - partread);
      bpos += nbytes - partread;
      bleft -= (nbytes - partread);

    }
  else
    {

      /* just transfer what was aked for */

//    printf("%p, %p - %p, %i\n", buf, rbuf, rbuf+bpos,nbytes);
      if (Pars.CurEvNo <= Pars.NumToPrint)
        printf ("simple transfer: now: bleft=%i, bpos=%i\n", bleft, bpos);

      memcpy (buf, rbuf + bpos, nbytes);
      bpos += nbytes;
      bleft -= nbytes;
      if (Pars.CurEvNo <= Pars.NumToPrint)
        printf ("transfer of %i bytes, set bleft=%i, set bpos=%i\n", nbytes, bleft, bpos);

    };
  if (Pars.CurEvNo <= Pars.NumToPrint)
    printf ("buf_read done:bpos=%i,bleft=%i\n", bpos, bleft);

  assert (bleft >= 0);


#endif

  /* done */

  return (nbytes);

};

/*----------------------------------------------------------------------------*/

int
GEBGetEv (GEB_EVENT * GEV_event, int curEvNo)
{

  /* declarations */

  static int nn = 0, ii = 0, nx = 0, nbadTS = 0, firsttime = 1;
  int siz, val, i1, i, stType;
  static int newbuf = 1, *pos;
  long long int TS, dTS;
  char str[256];

  /* prototypes */

  int bread (int, int *, int *, int *, int *);

#if(DEBUG2)
  printf ("GEBGetEv: called, nx=%i\n", nx);
  fflush (stdout);
#endif

  if (firsttime)
    {
      firsttime = 0;

      /* get the initial header, into position 1, not 0 */

      ii = 1;
      siz = fbuf_read (inData, (char *) GEV_event->ptgd[ii], sizeof (GEBDATA));
      if (siz != sizeof (GEBDATA))
        {
          printf ("failed to read %i bytes for header, got %i\n", sizeof (GEBDATA), siz);
          return (1);
        };
      Pars.nbytes += siz;
      nn2++;
      printf ("got initial header, TS=%lli ", GEV_event->ptgd[ii]->timestamp);
      printf ("ID: %i ", GEV_event->ptgd[ii]->type);
      printf ("length: %i\n", GEV_event->ptgd[ii]->length);

      /* get the initial payload */

      i1 = GEV_event->ptgd[ii]->length;
      siz = fbuf_read (inData, (char *) GEV_event->ptinp[ii], i1);
      if (siz != i1)
        {
          printf ("failed to read %i bytes for payload, got %i\n", i1, siz);
          return (2);
        };
      nn3++;
      Pars.nbytes += siz;
      printf ("read initial payload of siz=%i into event position %i\n", siz, ii);
      fflush (stdout);

      printf ("__ptgd[0]->type=%2i; ", GEV_event->ptgd[0]->type);
      printf ("ptgd[0]->length=%4i; ", GEV_event->ptgd[0]->length);
      printf ("ptgd[0]->timestamp=%lli\n", GEV_event->ptgd[0]->timestamp);
      fflush (stdout);

      ii = 1;
      nn = 1;
      printf ("initial ii=%i, nn=%i\n", ii, nn);

    };

  Pars.curTS = GEV_event->ptgd[0]->timestamp;
  TS = Pars.curTS;

  /* process leftovers from the last read */

  if (nn > 0)
    {

#if(DEBUG2)
      printf ("we have old geb/payload left over at position %i: \n", nn);
      fflush (stdout);
      printf ("__ptgd[nn]->type=%2i; ", GEV_event->ptgd[nn]->type);
      printf ("ptgd[nn]->length=%4i; ", GEV_event->ptgd[nn]->length);
      printf ("ptgd[nn]->timestamp=%lli\n", GEV_event->ptgd[nn]->timestamp);
      fflush (stdout);
#endif
      /* move the last (at pos nn) to the first position */

      memcpy ((char *) GEV_event->ptgd[0], (char *) GEV_event->ptgd[nn], sizeof (GEBDATA));
      memcpy ((char *) GEV_event->ptinp[0], (char *) GEV_event->ptinp[nn], GEV_event->ptgd[0]->length);

    }

  /* reset, ii=0 alwasy taken for the first leftover from last time */

  ii = 1;
  nn = 1;
  Pars.curTS = GEV_event->ptgd[0]->timestamp;

  GEV_event->mult = 0;
#if(DEBUG2)
  printf ("ii=%i, nn=%i\n", ii, nn);
  printf ("Pars.curTS=%lli, TS=%lli\n", Pars.curTS, TS);
#endif

  while ((TS - Pars.curTS) < Pars.dTS)
    {
      /*read geb header */

//      assert (i < MAXGEBS);
//      assert (GEV_event->ptgd[i] != NULL);

#if (0)
      printf ("for ii=%i:: trying to get a geb header of size %i\n", ii, sizeof (GEBDATA));
      printf ("GEV_event->ptgd[%i]=0x%p\n", ii, GEV_event->ptgd[ii]);
      fflush (stdout);
#endif

      /* trap for too long events */

      if (ii >= MAXGEBS)
        {
          printf ("error: this event is too long > %i\n", ii);
          return (1);
        };

      siz = fbuf_read (inData, (char *) GEV_event->ptgd[ii], sizeof (GEBDATA));
      if (siz != sizeof (GEBDATA))
        {
          printf ("failed to read %i bytes for header, got %i\n", sizeof (GEBDATA), siz);
          return (1);
        };
      nn2++;
      Pars.nbytes += siz;
      GEV_event->mult++;


#if (0)
      printf ("ii=%i, found header with TS=%lli, payload length=%i\n", ii, GEV_event->ptgd[ii]->timestamp,
              GEV_event->ptgd[ii]->length);
      fflush (stdout);
#endif
      TS = GEV_event->ptgd[ii]->timestamp;
#if(DEBUG2)
      printf ("ii=%i, nn=%i\n", ii, nn);
      printf ("Pars.curTS=%lli, TS=%lli\n", Pars.curTS, TS);
#endif

      /* periodic write out */
      /* shows the next event, but we don't care */

      if ((Pars.CurEvNo % Pars.modwrite == 0) )
        {
          printf ("GEBGetEv: ev# %5i ", Pars.CurEvNo);
          stType = GebTypeStr (GEV_event->ptgd[ii]->type, str);
          printf ("%s ", str);
          printf ("%4iBytes ", GEV_event->ptgd[ii]->length);
          printf ("TS=%lli ", GEV_event->ptgd[ii]->timestamp);
          printf ("curTS=%lli ", Pars.curTS);
          dTS = TS - Pars.curTS;
          printf ("dTS= %lli\n", dTS);
          fflush (stdout);
        };

      /* trap for bad timestamps */

      if (TS < Pars.curTS)
        {
          if (nbadTS < 100)
            {
              printf ("batflag:: TS<Pars.curTS, reset it at event # %i\n", Pars.CurEvNo);
              printf ("TS=%lli, Pars.curTS=%lli, DT=%lli\n", TS, Pars.curTS, TS - Pars.curTS);
              fflush (stdout);
            };
          Pars.curTS = TS;
          if (nbadTS < 100)
            {
              printf ("new Pars.curTS=%lli\n", Pars.curTS);
              printf ("we have read %lli bytes so far\n", Pars.nbytes);
            };
          nbadTS++;

#if(0)
          if (nbadTS > 1000)
            {
              printf ("too many bad TS, quit withe error code 3\n");
              fflush (stdout);
//              return (3);
            };
#endif
        };


      /* read payload */

      i1 = GEV_event->ptgd[ii]->length;
      siz = fbuf_read (inData, (char *) GEV_event->ptinp[ii], i1);
      if (siz != i1)
        {
          printf ("failed to read %i bytes for payload, got %i\n", i1, siz);
          return (2);
        };
      nn3++;
      Pars.nbytes += siz;
#if (DEBUG2)
      printf ("__read payload of siz=%i into event position %i\n", siz, ii);
      fflush (stdout);
#endif

      ii++;
      nn++;
#if (DEBUG2)
      printf ("1: ii=%i, nn=%i\n", ii, nn);
      fflush (stdout);
#endif

     }

  ii--;
  nn--;
#if (DEBUG2)
  printf ("2: ii=%i, nn=%i\n", ii, nn);
  fflush (stdout);
#endif

  /* return the mutiplicity */

  GEV_event->mult = ii;

#if(DEBUG2)
  printf ("complete event, next TS is %lli or  %lli out\n", TS, TS - Pars.curTS);
  printf ("we found %i events in coincidence, timestamps are\n", GEV_event->mult);
  for (i = 0; i < GEV_event->mult; i++)
    printf ("[%i] TS=%lli\n", i, GEV_event->ptgd[i]->timestamp);
  printf ("we have read %lli bytes so far\n", Pars.nbytes);
  printf ("GEV_event->mult=%i\n", GEV_event->mult);
  fflush (stdout);
  nx++;
  if (nx > 3)
    exit (0);
#endif

#if(DEBUG2)
  printf ("GEBGetEv: done %i\n");
  fflush (stdout);
#endif

  return (0);
}


/*----------------------------------------------------------------------------*/
int main (int argc, char **argv) {

  /*--------------*/
  /* declarations */
  /*--------------*/

  int j, i, HaveChatFile = 0;
  char *p;
  char ChatFileName[STRLEN];
  int GEBacq (char *);
  int time_stamp (FILE *);
  char str2[STRLEN], str3[STRLEN], str4[STRLEN];
  struct stat st;


  /*------*/
  /* help */
  /*------*/

  if (argc < 2)
    {
      printf ("\n");
      printf ("use: %s -chat file [-help] [-version] .... TBD\n", argv[0]);
      printf ("\n");
      return (0);
    };

  /* initialize */

  printf ("started GEBSort at ");
  time_stamp (stderr);

  Pars.InputSrc = NOTDEF;
  Pars.HaveRootFileName = 0;
  sprintf (Pars.ROOTFileOption, "UNDEFINED");
  Pars.GGMAX = 2000;
  Pars.ndetlimlo = 1;
  Pars.ndetlimhi = 8;
  for (i = 0; i < MAXNOSEG; i++)
    {
      Pars.fomlo[i] = 0;
      Pars.fomhi[i] = 2.0;
    }
  Pars.UpdateRootFile = 0;
//  Pars.UseShareMemFile = FALSE;
  Pars.StartMapAddress = 0;
  sprintf (Pars.ShareMemFile, "GTSort.map");
  Pars.maxsnglintrE = 2000.0;
  Pars.maxsnglintrEFOM = 1.81;
  Pars.minnumHitArray = 0;
  Pars.maxnumHitArray = 200;
  Pars.minNumGammas = 0;
  Pars.maxNumGammas = 200;
  Pars.minNumCC = 0;
  Pars.maxNumCC = 200;
  Pars.maxTS = MAXLONG;

/*
  for (i = 0; i < MAXDETNO; i++)
    {
      Pars.CCcal_offset[i] = 0;
      Pars.CCcal_gain[i] = 1.0;
      Pars.enabled[i] = 1;
      for (j = 0; j <= MAXCRYSTALNO; j++)
        {
          Pars.SEGcal_gain[i][j] = 1.0;
          Pars.SEGcal_offset[i][j] = 0.0;
        };
    }
*/
  /*--------------------*/
  /* Parse command line */
  /* and call GEBacq     */
  /*--------------------*/

  j = 1;                        /* current command line arg position */

  printf ("we have %i arguments\n", argc);
  fflush (stdout);

  if (argc == 1)
    {
      printf ("help: see go file for examples of use of GEBSort");
      exit (0);
    };

  if (argc > 1)
    while (j < argc)
      {
        printf ("%i... %s\n", j, argv[j]);
        fflush (stdout);

        if ((p = strstr (argv[j], "-version")) != NULL)
          {
            printf ("try: svn info\n");
            exit (0);
          }
        else if ((p = strstr (argv[j], "-help")) != NULL)
          {
            printf ("\n");
            printf ("GEBSort is documented on the WWW, URL:\n");
            printf ("\n");
            printf ("http://wiki.anl.gov/gretina_at_anl \n");
            printf ("\n");
            exit (0);

          }
        else if ((p = strstr (argv[j], "-input")) != NULL)
          {
            /* check that user specified enough arguments */
            j++;

            /* determine input source */
            strcpy (str2, argv[j++]);
            if (strcmp ("disk", str2) == 0)
              {
                strcpy (str3, argv[j++]);
//                strcpy (Pars.ROOTFileOption, argv[j++]);
                printf ("will take input from disk\n");
                strcpy (Pars.GTSortInputFile, str3);
                Pars.InputSrc = DISK;
                printf("file : %s \n", str3);
                stat (Pars.GTSortInputFile, &st);
                printf ("file size is %lli bytes or %lli MBytes\n", st.st_size, st.st_size / 1024 / 1000);
                fflush (stdout);
              }
            else if (strcmp ("geb", str2) == 0)
              {
                strcpy (Pars.pHost, argv[j++]);
                Pars.grouping = atol (argv[j++]);
                Pars.type = atol (argv[j++]);
                Pars.timeout = (float) atof (argv[j++]);

                printf ("Pars.pHost=%s\n", Pars.pHost);
                printf ("Pars.grouping=%i\n", Pars.grouping);
                printf ("Pars.type=%i\n", Pars.type);
                printf ("Pars.timeout=%f\n", Pars.timeout);
                Pars.InputSrc = GEB;
//                strcpy (Pars.ROOTFileOption, argv[j++]);
                printf ("root file option: %s\n", Pars.ROOTFileOption);
#if(HAVE_VXWORKS==0)
                printf ("oppsie... you cannot specify this option unless\n");
                printf ("you have #define HAVE_VXWORKS 1 in GEBSort.h\n");
                printf ("and have a VxWorks license, quit\n");
                exit (0);
#endif
              }
            else
              {
                printf ("unknown input option: %s\n", str2);
                printf ("aborting\n");
                fflush (stdout);
                exit (0);
              };

            printf ("\n");

          }
        else if ((p = strstr (argv[j], "-chat")) != NULL)
          {
            j++;
            strcpy (ChatFileName, argv[j++]);
            printf ("will read sorting instructions from chatfile: %s\n", ChatFileName);
            system ("ulimit -a");
            HaveChatFile = 1;
            printf ("\n");
          }
        else if ((p = strstr (argv[j], "-rootfile")) != NULL)
          {
            j++;
            strcpy (Pars.ROOTFile, argv[j++]);
            printf ("rootfile name specified on command line\n");
            printf ("--> will store spectra in rootfile: %s\n", Pars.ROOTFile);
            
            Pars.HaveRootFileName = 1;
            Pars.UseRootFile = 1;
            if ((p = strstr (argv[j], "RECREATE")) != NULL)
              {
                Pars.UpdateRootFile = FALSE;
                printf ("will recreate root file\n");
                sprintf (Pars.ROOTFileOption, "RECREATE");
                j++;
              }
            else if ((p = strstr (argv[j], "UPDATE")) != NULL)
              {
                Pars.UpdateRootFile = TRUE;
                printf ("will update root file\n");
                sprintf (Pars.ROOTFileOption, "UPDATE");
                j++;
              }
            else
              {
                printf (" you must specify RECREATE or UPDATE after -rootfile\n");
                exit (0);
              }
            printf("\n");  
          }
        else
          {
            printf ("command line argument not understood!\n");

            printf ("%s: I was called as: \n--->[", argv[0]);
            for (i = 0; i < argc; i++)
              {
                printf ("%s ", argv[i]);
                fflush (stdout);
              }
            printf ("]\non ");
            time_stamp (stderr);
            exit (0);

          }
      };

  /* now start the sorter */
  if (HaveChatFile){
    printf( "------- run GEBacq(ChatFile) in GEBSort.cxx\n");
    GEBacq (ChatFileName);
  }else
    {
      printf ("you must specify a chat script\n");
      exit (0);
    }

}

/*--------------------------------------------------------------------------*/

void
signal_catcher (int sigval)
{
  int time_stamp (FILE *);
  printf ("\nGSSort/GEBacq received signal <%i> on ", sigval);
  time_stamp (stderr);
  Pars.WeWereSignalled = TRUE;
  fflush (stdout);
}

/*---------------------------------------------------------------------------*/

#if(HAVE_VXWORKS)
void
sdummyload (Long_t size)
{

  /* dummy load a shared memory to find out what */
  /* start address it chooses................... */

  TMapFile *m;
  m = TMapFile::Create ("dummy.map", "recreate", size);
#if(MAC==1)
  Pars.StartMapAddress = (unsigned long long int) m->GetMmallocDesc ();
#else
  Pars.StartMapAddress = (unsigned int) m->GetMmallocDesc ();
#endif
  m->Print ();

  /* close and remove dummy map file */

  m->Close ();
  gSystem->Exec ("\\rm dummy.map");
}
#endif

/*---------------------------------------------------------------------------*/
int GEBSort_read_chat (char *name) {

  /* declarations */

  FILE *fp, *fp1;
  char *pc, *pc1, str[STRLEN] = { '0' }, str1[STRLEN] = {'0'}, str2[STRLEN] = {'0'};
  char str3[STRLEN], str4[STRLEN], str5[STRLEN], str6[STRLEN];
  int nn = 0, nni = 0, st, PType;
  char *p;
  int i, j, k, l, m, n, i1, i2, i3, i4, i5, i6;
  int j1, j2, j3, j4, j5, j6, j7;
  float f1, f2, f3, f4, pi, r1, r2, r3, r4, rr;
  int echo = 0, nret;
  double d1;
  unsigned int ui1;
  double pol, azi;

  /* prototypes */
  int str_decomp (char *, int, int *, int);

/*
  TCutG *rd2dwin (Char_t *);
  int FindPEvMatNo (char *, int *);
  void FindCondNo (char *, int *);
  int SetFERAVSN (char *, int);
  void InitializeFERAvsntable ();
  void ZeroFERAtypeVSN ();
  void PrintFERATypes ();
  void SetNPosWarn (int);
  void SetRecordVer_tape (int);
  void SetRecordVer_disk (int);
  void SetRecordVer_net (int);

  void FindTimeMaskNo (char *, int *);
  int RdOffFile (char *, int *);
  int RdGeCalFile (char *, float *, float *);
  void CheckNoArgs (int, int, char *);
//  void SetSpecial (char *str);
//  void SetExportModNWords (char *, int);
//  void SetlongRaNCHANNELSTDCNWords (char *, int);
 // void setIsomerIDs (int);
  void SetClockPar (int, float, float);
  void SetFERAvsntable (int, int);
  void SetnFeraDebug (int);
*/

  /* open chat file */

  if ((fp = fopen (name, "r")) == NULL){
      printf ("error: could not open chat file: <%s>\n", name);
      exit (0);
  };
  printf ("chat file: <%s> open\n", name);
  printf ("\n");
  fflush (stdout);

  /* read content and act */

  pc = fgets (str, STRLEN, fp);

  /* rmEndComment(str, STRLEN); */

  while (pc != NULL)
    {
      if (echo)
        printf ("chat->%s", str);
      fflush (stdout);

      /* attemp to interpret the line */

      if ((p = strstr (str, "nevents")) != NULL)
        {
          nret = sscanf (str, "%s %i", str1, &Pars.nEvents);
          CheckNoArgs (nret, 2, str);
          printf ("will sort a max of %i events\n", Pars.nEvents);
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
      else if ((p = strstr (str, "enabled")) != NULL)
        {
          nret = sscanf (str, "%s %s", str1, str2);
          CheckNoArgs (nret, 2, str);
          printf ("%s\n", str2);
          for (j = 0; j < MAXDETNO; j++)
            Pars.enabled[j] = 0;
          str_decomp (str2, MAXDETNO + 1, Pars.enabled, 1);
          for (j = 0; j < MAXDETNO; j++)
            if (Pars.enabled[j])
              printf ("detector %3i is ENABLED\n", j);
        }
      else if ((p = strstr (str, "firstEvent")) != NULL)
        {
          nret = sscanf (str, "%s %i", str1, &Pars.firstEvent);
          CheckNoArgs (nret, 2, str);
          printf ("will start sorting at event %d\n", Pars.firstEvent);
          fflush (stdout);

        }
      else if ((p = strstr (str, "echo_data")) != NULL)
        {
          nret = sscanf (str, "%s %s", str1, &Pars.echo_data_fn);
          CheckNoArgs (nret, 2, str);
          Pars.echo_data = 1;
          printf ("will echo data to %s\n", Pars.echo_data_fn);
          fflush (stdout);
          Pars.echo_data_pipe = open (Pars.echo_data_fn, O_WRONLY | O_CREAT | O_TRUNC, PMODE);
          if (Pars.echo_data_pipe == 0)
            {
              printf ("could not open input data file \"%s\", quit!\n", Pars.echo_data_fn);
              exit (1);
            };
          printf ("echo data file \"%s\" is open\n", Pars.echo_data_fn);
        }
      else if ((p = strstr (str, "bin_none")) != NULL)
        {
          Pars.do_bin_mode3 = 0;
          Pars.do_bin_mode2 = 0;
          Pars.do_bin_mode1 = 0;
          Pars.do_bin_dgs = 0;
          Pars.do_bin_XA = 0;
        }
      else if ((p = strstr (str, "bin_mode3")) != NULL)
        {
          Pars.do_bin_mode3 = 1;
          printf ("Pars.do_bin_mode3=%i\n", Pars.do_bin_mode3);
        }
      else if ((p = strstr (str, "bin_mode2")) != NULL)
        {
          Pars.do_bin_mode2 = 1;
          printf ("Pars.do_bin_mode2=%i\n", Pars.do_bin_mode2);
        }
      else if ((p = strstr (str, "bin_mode1")) != NULL)
        {
          Pars.do_bin_mode1 = 1;
          printf ("Pars.do_bin_mode1=%i\n", Pars.do_bin_mode1);
        }
      else if ((p = strstr (str, "bin_dgs")) != NULL)
        {
          Pars.do_bin_dgs = 1;
          printf ("Pars.do_bin_dgs=%i\n", Pars.do_bin_dgs);
        }
      else if ((p = strstr (str, "bin_XA")) != NULL)
        {
          Pars.do_bin_XA = 1;
          printf ("Pars.do_bin_XA=%i\n", Pars.do_bin_XA);
        }
      else if ((p = strstr (str, "waitfordata")) != NULL)
        {
          nret = sscanf (str, "%s %i", str1, &Pars.waitfordataseconds);
          CheckNoArgs (nret, 2, str);
        }
     else if ((p = strstr (str, "xa_MM")) != NULL)
        {
          nret = sscanf (str, "%s %f", str1, &Pars.xa_MM);
          CheckNoArgs (nret, 2, str);
          printf ("will use M= %f in bin_dgs\n", Pars.xa_MM);
        }
      else if ((p = strstr (str, "xa_PZ")) != NULL)
        {
          nret = sscanf (str, "%s %s", str1, Pars.xa_PZfn);
          CheckNoArgs (nret, 2, str);
          printf ("will use PZ file \"%s\" in bin_dgs\n", Pars.xa_PZfn);
        }
      else if ((p = strstr (str, "xa_ecal")) != NULL)
        {
          nret = sscanf (str, "%s %s", str1, Pars.xa_ecalfn);
          CheckNoArgs (nret, 2, str);
          printf ("will use ecal file \"%s\" in bin_dgs\n", Pars.xa_ecalfn);
        }
     else if ((p = strstr (str, "maxDataTime")) != NULL)
        {
          nret = sscanf (str, "%s %i", str1, &i1);
          CheckNoArgs (nret, 2, str);
          Pars.tmpmaxTS = (long long int) i1;
          printf ("will sort %lli minutes of data\n", Pars.tmpmaxTS);
          // NOTE: will be re-calculated when we have the first timestamp 
        }
      else if ((p = strstr (str, "timewin")) != NULL)
        {
          nret = sscanf (str, "%s %lli", str1, &Pars.dTS);
          CheckNoArgs (nret, 2, str);
          printf ("will sort using a time window of %lli\n", Pars.dTS);
          fflush (stdout);
        }
      else if ((p = strstr (str, "ndetlimits")) != NULL)
        {
          nret = sscanf (str, "%s %i %i", str1, &Pars.ndetlimlo, &Pars.ndetlimhi);
          CheckNoArgs (nret, 3, str);
          printf ("Pars.ndetlimlo=%i,", Pars.ndetlimlo);
          printf ("Pars.ndetlimhi=%i\n", Pars.ndetlimhi);
        }
      else if ((p = strstr (str, "modwrite")) != NULL)
        {
          nret = sscanf (str, "%s %i", str1, &Pars.modwrite);
          CheckNoArgs (nret, 2, str);
        }
      else if ((p = strstr (str, "tsnumwrites")) != NULL)
        {
          nret = sscanf (str, "%s %i", str1, &Pars.tsnumwrites);
          CheckNoArgs (nret, 2, str);

        }
      else if ((p = strstr (str, "printevents")) != NULL)
        {
          nret = sscanf (str, "%s %i", str1, &Pars.NumToPrint);
          CheckNoArgs (nret, 2, str);
          printf ("will print details of the first %i events\n", Pars.NumToPrint);
          fflush (stdout);

        }
      else if ((p = strstr (str, "multlims")) != NULL)
        {
          nret = sscanf (str, "%s %i %i", str1, &Pars.multlo, &Pars.multhi);
          CheckNoArgs (nret, 2, str);
          printf ("will require mult bt: %i %i\n", Pars.multlo, Pars.multhi);
          fflush (stdout);
        }
      else if ((p = strstr (str, "gglen")) != NULL)
        {
          nret = sscanf (str, "%s %i", str1, &Pars.GGMAX);
          CheckNoArgs (nret, 2, str);
        }
      else if ((p = strstr (str, "exit")) != NULL)
        {

          printf ("will skip rest of chat file\n");
          fclose (fp);
          printf ("\n");
          printf ("chat file: <%s> closed\n", name);
          printf ("processed %i sort instructions and %i lines\n", nni, nn);
          printf ("\n");
          fflush (stdout);
          return (0);

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
  printf ("processed %i sort instructions and %i lines\n", nni, nn);
  printf ("\n");
  fflush (stdout);
  return (0);

};


/*----------------------------------------------------------------------------*/

int showStatus(){
  printf ("read %i events; ", (Pars.CurEvNo - Pars.firstEvent));
  printf ("Pars.beta=%6.4f; ", (float) Pars.beta);
  printf ("time since last update %i minutes\n", (int) tdmp);
//  printf ("CommandFileName=\"%s\"\n", CommandFileName);

  fflush (stdout);

  return (0);

};

/* ----------------------------------------------------------------- */

int
echo_data (GEB_EVENT * GEB_event)
{

  /* Repeat the current event to a file after validation */
  /* per tradition, we do not write GEB_event->mult out */

  int i, siz;

  for (i = 0; i < GEB_event->mult; i++)
    {

      /* write geb header */

      siz = write (Pars.echo_data_pipe, (char *) GEB_event->ptgd[i], sizeof (GEBDATA));
      assert (siz == sizeof (GEBDATA));

      /* write payload */

      siz = write (Pars.echo_data_pipe, (char *) GEB_event->ptinp[i], GEB_event->ptgd[i]->length);
      assert (siz == GEB_event->ptgd[i]->length);


    };

  /* done */

  return (0);

};

/*----------------------------------------------------------------------------*/
int GEBacq (char *ChatFileName) {

  /* declarations */

  int NprintEvNo = 0, in, zero = 0;
  GEB_EVENT GEB_event;
  int st = 0, eov = 0, nWords = 0, i1, i2, i, j, nret, siz;
  int ii, jj;
  char str[256], str1[256], str2[246];
  FILE *fp;
  time_t t1, t2;
  Int_t ComPressLevel = NOTDEF;
  char *p, buffer[512];
  FILE *fp0;
  int ir[NUMAGATAPOS], dummy_i[NUMAGATAPOS];;
  double rn, dtmpehi, d1, nsec;
  DGSHEADER dgsHeader;
  long long int tac, tFP;
  int nFP;
  float r1, r2, r3;
  static int firsttime = 1;
  static long long int t0, TSprev = 0;
  long long int tcur;
  unsigned int typehit[MAX_GEB_TYPE];
  FILE *TSfile;
  TH2F *dtbtev;
  long long int firtsTSinEvent, dTS;
  int dim;
  float rr[LONGLEN + 1];
  int ok_validate;
  long long int nok_validate = 0, ntot_validate = 0;
  float azi, pol;
  struct stat fst;
  long long int li1;

  /* root spectra pointers */

  TH1 *hhtemp;

  TList *zlist;
  TIterator *hiterator;

  TMapFile *mfile;

  TH1D *tmpTH1D = NULL;
  TH2F *tmpTH2F = NULL;

  /* prototypes */

  int GEBSort_read_chat (char *);
  int wr_spe (char *, int *, float *);

  /* data type binners sup==setup, bin_==binner */

  int exit_bin_mode1 ();
  int exit_bin_mode2 ();
  int exit_bin_mode3 ();
  int exit_dgs ();
  int exit_XA ();
  int sup_mode1 ();
  int sup_mode2 ();
  int sup_mode3 ();
  int sup_dgs ();
  int sup_XA ();
  int bin_mode1 (GEB_EVENT *);
  int bin_mode2 (GEB_EVENT *);
  int bin_mode3 (GEB_EVENT *);
  int bin_dgs (GEB_EVENT *);
  int bin_XA (GEB_EVENT *);
  int validate (GEB_EVENT *);

  /*-------*/
  /* Hello */
  /*-------*/

  printf ("\n");
  printf ("GEBsort running on: ");
  time_stamp (stdout);
  printf ("\n");

  /*------------*/
  /* initialize */
  /*------------*/

  /* Note: do not zero out Pars as there are */
  /* command line inputs here already */
  /*  bzero((char *) &Pars, sizeof (PARS)); */

  printf ("\n");
  printf ("initializing Pars, defined at GEBSort.h\n");
  printf ("\n");
  fflush (stdout);

  // MAX_GEB_TYPE = 25, is defined in gdecomp.h 
  for (i = 0; i < MAX_GEB_TYPE; i++) typehit[i] = 0;

  for (i = 0; i <= MAXDETNO; i++){  // MAXDETNO = 200 defined in ctk.h
    Pars.xyzoffset[i][0] = 0;
    Pars.xyzoffset[i][1] = 0;
    Pars.xyzoffset[i][2] = 0;
  };

  Pars.nEvents = 2000000000;
  Pars.WeWereSignalled = FALSE; /* signal  */
  Pars.UseRootFile = FALSE;
  Pars.SizeShareMemFile = FALSE;
  Pars.spname[STRLEN];
  Pars.firstEvent = 0;
  Pars.GSudpPort = 1101;
  Pars.NumToPrint = 25;
  Pars.DumpEvery = 10;
  Pars.dTS = 500;
  Pars.nbytes = 0;
  Pars.beta = 0;
  Pars.modwrite = 1000;
  Pars.tsnumwrites = 100;
  Pars.nocrystaltoworldrot = 0;
  Pars.multlo = 1;
  Pars.multhi = 20;
  Pars.requiretracked = 0;
  Pars.vetoSpots = 0;
  Pars.target_x = 0;
  Pars.target_y = 0;
  Pars.target_z = 0;
  Pars.crystalID3D = -1;
  Pars.Hresolution = 1000;
  Pars.havexyzoffset = 0;

  Pars.do_bin_mode3 = 0;
  Pars.do_bin_mode2 = 0;
  Pars.do_bin_mode1 = 0;
  Pars.do_bin_dgs = 0;

  Pars.echo_data = 0;

  for (i = 0; i < MAXGEBS; i++) {// MAXGEBS = 1000, defined at GEBSort.h
    GEB_event.ptgd[i] = (GEBDATA *) calloc (2 * sizeof (GEBDATA), 1);
    GEB_event.ptinp[i] = (CRYS_INTPTS *) calloc (MAXPAYLOADSIZE, 1);
  };

  rbuf = (char *) calloc (RBUFSIZE + 1, 1); //RBUFSIZE = 500000, defined at GEBSort.h

  /*------------------*/
  /* read chat script */
  /*------------------*/

  GEBSort_read_chat (ChatFileName);

  printf ("\nsorting this:\n");
  printf ("Pars.do_bin_mode3 = %i\n", Pars.do_bin_mode3);
  printf ("Pars.do_bin_mode2 = %i\n", Pars.do_bin_mode2);
  printf ("Pars.do_bin_mode1 = %i\n", Pars.do_bin_mode1);
  printf ("Pars.do_bin_dgs = %i\n", Pars.do_bin_dgs);
  printf ("Pars.do_bin_XA = %i\n", Pars.do_bin_XA);
  printf ( "\n");
  

  printf ("checking proper input of chat file...\n");
  if (Pars.InputSrc == NOTDEF){
    printf ("you must specify an input source\n");
    exit (1);
  }else if (Pars.InputSrc == DISK){
    /* attempt to open input file or input stream */
    if ((p = strstr (Pars.GTSortInputFile, "STDIN")) != NULL){
      inData = stdin;
    }else{
#if(ISMAC)
        inData = fopen (Pars.GTSortInputFile, "r");
#else
        inData = fopen64 (Pars.GTSortInputFile, "r");
#endif
    }
      
    if (inData == NULL){
      printf ("could not open\"%s\"; quit\n", Pars.GTSortInputFile);
      exit (1);
    }else{
      printf ("input file \"%s\" is open, inData=%p\n", Pars.GTSortInputFile, inData);
    }

    /* get the file size */

    //  fstat (inData, &fst);
    li1 = fleft (inData);
    printf ("file size is %lli bytes or %lli MBs \n", li1, li1 / 1024 / 1000);
    fflush (stdout);

    /* check for badstate */
    if ((p = strstr (Pars.GTSortInputFile, "STDIN")) != NULL && Pars.waitfordataseconds > 0){
      printf ("you cannot specify waifordata and STDIN at the same time, quit\n");
      exit (0);
    };

  }else if (Pars.InputSrc == GEB){

#if(HAVE_VXWORKS)
    printf ("will take input from GEB with these parameters:\n");
    printf ("Pars.pHost=%s; ", Pars.pHost);
    printf ("Pars.grouping=%i; ", Pars.grouping);
    printf ("Pars.type=%i; ", Pars.type);
    printf ("Pars.timeout=%f\n", Pars.timeout);
    printf ("connecting to tap\n");

    for (i = 0; i < ConnectionRetryCount; i++){
      printf ("connecting, %i/%i, %i\n", i, ConnectionRetryCount, ConnectionRetryWait / 1000);
      fflush (stdout);

      pTap = gretTapConnect (Pars.pHost, GRETTAP_GEB, Pars.type);
      if (pTap || (gretTapClientError != GTC_TAPCONN)){
        /* Either success or a failure other than connection failure */
        printf ("got here, that is good, break out and go on\n");
        break;
      }
      usleep (ConnectionRetryWait);

#if(DEBUG)
      fprintf (stderr, "Retry number %d\n", i);
      fflush (stderr);
#endif
    }
    
    fprintf (stderr, "Retries: %d\n", i);
    fflush (stderr);
    if (!pTap){
      fprintf (stderr, "Unable to connect to tap server at %s : %s\n", Pars.pHost, gretTapClientErrorStrings[gretTapClientError]);
      exit (1);
    };
#endif
  }else{
    printf ("input source not recognized, quit\n");
    exit (1);
  };

  printf ("input source is set up, supposedly\n");

  /*---------------------*/
  /* other sanity checks */
  /*---------------------*/

  if (Pars.InputSrc == NOTDEF){
    printf ("you must specify an input source!\n");
    printf ("quitting...\n");
    exit (1);
  };

  NprintEvNo = 0;
  Pars.CurEvNo = 0;

  if (!Pars.ShareMemFile){
    Pars.DumpEvery = 2000000000;
    printf ("\n");
    printf ("_since rootfile: setting `Pars.DumpEvery` to infinity..!\n");
    printf ("\n");
  };

  /*------------------------------------------*/
  /* if we are using root file, then either   */
  /* read in old rootfile or create a nev one */
  /*------------------------------------------*/

  if (Pars.UpdateRootFile){  //who want to update a root file? just create a new one please...
    /* check here whether the old root file exists */
    fp = fopen (Pars.ROOTFile, "r");
    if (fp == NULL){
      printf ("could not open old rootfile: %s\n", Pars.ROOTFile);
      printf ("the old rootfile must exist if you  \n");
      printf ("want to use the UPDATE option\n");
      printf ("aborting...\n");
      exit (0);
    };
    fclose (fp);

    // read in old root file 

    Pars.f1 = NULL;
    Pars.f1 = new TFile (Pars.ROOTFile, "UPDATE");
    printf ("read old root file <%s>\n", Pars.ROOTFile);
    if (!Pars.f1->IsOpen ()){
      printf ("could not open file....\n\n");
      exit (-1);
    };
    printf ("base=<%s>\n", Pars.f1->GetPath ());
    Pars.f1->Print ();

  }else{ // create the rootfile 

    Pars.f1 = NULL;
    Pars.f1 = new TFile (Pars.ROOTFile, "RECREATE");
    printf ("root file <%s>\n", Pars.ROOTFile);
    if (!Pars.f1->IsOpen ()){
      printf ("could not open file....\n\n");
      exit (-1);
    };
    printf ("base=<%s>\n", Pars.f1->GetPath ());
    Pars.f1->Print ();
  };


  TSfile = fopen ("TS.list", "w");

  /*--------------------------------*/
  /* setup the root spectra we need */
  /*--------------------------------*/

  /* spectra that are always there */

  sprintf (str1, "dtbtev");
  sprintf (str2, "dtbtev");
  dtbtev = mkTH2F (str1, str2, DTBTEVLEN / 2, 0, DTBTEVLEN, MAX_GEB_TYPE, 1, MAX_GEB_TYPE);
  sprintf (str1, "delta t");
  dtbtev->SetXTitle (str1);
  sprintf (str1, "type");
  dtbtev->SetYTitle (str1);

  /* spectra for different types of data */
  
  printf("\n---------- Run sup_XXX () \n");
  if (Pars.do_bin_mode2 == 1)  sup_mode2 ();
  if (Pars.do_bin_mode1 == 1)  sup_mode1 ();
  if (Pars.do_bin_mode3 == 1) sup_mode3 ();
  if (Pars.do_bin_dgs == 1) sup_dgs ();
  if (Pars.do_bin_XA == 1) sup_XA ();

  printf ("we have define the following ROOT spectra:\n");

  Pars.wlist = gDirectory->GetList ();
  Pars.wlist->Print ();

  /*----------------------*/
  /* setup signal catcher */
  /*----------------------*/

  printf ("\nsetup signal catcher\n");

#ifdef LINUX
  signal (SIGHUP, signal_catcher);
#endif
#ifdef SOLARIS
  sigset (SIGHUP, signal_catcher);
#endif


  /*---------------*/
  /* start sorting */
  /*---------------*/

  printf ("started sorting... ");
  if (Pars.InputSrc == DISK){
    printf ("from disk...\n");
  }else if (Pars.InputSrc == GEB){
    printf ("from GEB...\n");
  }else if (Pars.InputSrc == NET){
    printf ("from net... SHOULD NOT HAPPEN\n");
    exit (1);
  };
  printf ("\n");
  fflush (stdout);

  /* make sure we make it to the first event */
  /* before Pars.maxTS can cut the sort off */

  printf ("Pars.maxTS=%lli\n", Pars.maxTS);

  tdmplast = time (NULL);

  while (st >= 0 && ((Pars.CurEvNo - Pars.firstEvent) < Pars.nEvents) && (eov == 0) && (Pars.curTS <= Pars.maxTS)){

    //memset ((char *) &CEvent, 0, sizeof (COINEV));

    /*----------------*/
    /* get next event */
    /*----------------*/

    Pars.CurEvNo++;
    NprintEvNo++;

#if(DEBUG2)
    printf ("calling GEBGetEv, Pars.CurEvNo=%i\n", Pars.CurEvNo);
#endif

    st = GEBGetEv (&GEB_event, Pars.CurEvNo);

    /* debug print some events */

    if (Pars.CurEvNo <= Pars.NumToPrint){
      printf ("\n+++++++++++++++++++++++++++++++\n");
      printf ("*start event # %i with multiplicity %i looks like this:\n", Pars.CurEvNo, GEB_event.mult);
      for (i = 0; i < GEB_event.mult; i++){
        GebTypeStr (GEB_event.ptgd[i]->type, str);
        printf ("%2i> %2i, %s, TS=%lli\n", i, GEB_event.ptgd[i]->type, str, GEB_event.ptgd[i]->timestamp);
      };
    };

    if (st == 0 && Pars.CurEvNo < Pars.tsnumwrites){
      for (i = 0; i < GEB_event.mult; i++){
        if (i == 0) fprintf (TSfile, "\n");
        GebTypeStr (GEB_event.ptgd[i]->type, str);
        fprintf (TSfile, "%4i/%2i: (%2i,%s) TS=%20lli; ", Pars.CurEvNo, i, GEB_event.ptgd[i]->type, str, GEB_event.ptgd[i]->timestamp);
        fprintf (TSfile, "dT=%lli\n", GEB_event.ptgd[i]->timestamp - TSprev);
        TSprev = GEB_event.ptgd[i]->timestamp;
      }
    };

#if(DEBUG2)
        printf ("st=%i\n", st);
        printf ("GEB_event.mult=%i\n", GEB_event.mult);
        fflush (stdout);
        if (1) exit (0);
#endif

    if (st == 0){
      if (firsttime){
        firsttime = 0;
        t0 = GEB_event.ptgd[0]->timestamp;
        printf ("t0=%lli\n", t0);
        printf ("first event: GEB_event.mult=%i\n", GEB_event.mult);
        for (i = 0; i < GEB_event.mult; i++){
          GebTypeStr (GEB_event.ptgd[i]->type, str);
          printf ("%4i/%2i: (%2i,%s) TS=%20lli; ", Pars.CurEvNo, i, GEB_event.ptgd[i]->type, str, GEB_event.ptgd[i]->timestamp);
          printf ("dT=%lli\n", GEB_event.ptgd[i]->timestamp - TSprev);
        };

        /* calculate the termination timestamp */
        /* now that we know the first timestamp */

        if (Pars.tmpmaxTS > 0)
            Pars.maxTS = t0 + Pars.tmpmaxTS * 100000000;
        else
            Pars.maxTS = MAXLONG;

        printf ("maxTS=%lli\n", Pars.maxTS);

        d1 = (double) (Pars.maxTS - t0);
        d1 /= 100000000;
        printf ("that is %.1f seconds or ", (float) d1);
        nsec = d1;
        i1 = (unsigned int) d1 / 3600;
        d1 -= i1 * 3600;
        i2 = (unsigned int) d1 / 60;
        d1 -= i2 * 60;
        printf ("%ih%im%is minutes of sorting\n", i1, i2, (int) d1);

      }

      tcur = GEB_event.ptgd[0]->timestamp;


      /* count data types */
      for (i = 0; i < GEB_event.mult; i++){
        if (GEB_event.ptgd[i]->type > 0 && GEB_event.ptgd[i]->type < MAX_GEB_TYPE) typehit[GEB_event.ptgd[i]->type]++;
      };

      /* fill dtbtev spectrum */
      firtsTSinEvent = LLONG_MAX;
      for (i = 0; i < GEB_event.mult; i++){
          if (GEB_event.ptgd[i]->timestamp < firtsTSinEvent) firtsTSinEvent = GEB_event.ptgd[i]->timestamp;
      }
      
      for (i = 0; i < GEB_event.mult; i++){
        dTS = GEB_event.ptgd[i]->timestamp - firtsTSinEvent;
        d1 = (double) dTS;
        if (d1 >= (double) 0 && d1 < RATELEN) dtbtev->Fill (d1, GEB_event.ptgd[i]->type, 1);
      };
    };

    /*----------------------------------------*/
    /* allow user to manipulate raw data here */
    /*----------------------------------------*/

    if (st != 0){
      printf ("GEBSort: GEBGetEv returned %i, FAIL\n", st);
      printf ("we have read %lli bytes; ", Pars.nbytes);
      printf ("CurEvNo=%i\n", Pars.CurEvNo);
      fflush (stdout);

      /* terminate sort */

      eov = 1;

      /* note: */
      /* we might want to wait and try GEBGetEv */
      /* later to give the impresssion of interactivity */
      /* here in some future version... */
    }
    
    
//      if (st == 0 && (GEB_event.mult <30))
    if (st == 0 ) {

      /*----------------------------*/
      /* good event, now process it */
      /*----------------------------*/

      /* statistics */

      if (Pars.CurEvNo <= Pars.NumToPrint && 0){
        printf ("GEBSort: GEBGetEv returned st=%i, OK\n", st);
        printf ("we have read %lli bytes; ", Pars.nbytes);
        printf ("CurEvNo=%i\n", Pars.CurEvNo);
        fflush (stdout);
      };

      /* zap out the exchange structure so we */
      /* do not catch anything from last coincidence */

      bzero ((void *) &exchange, sizeof (EXCHANGE));

      /* search for external */
      /* if not needed, make it return 1 */

      ok_validate = validate (&GEB_event);
      ntot_validate++;

      if (ok_validate){
        nok_validate++;

        /* echo data to another file */
        /* after validation, do this before */
        /* the bin functions who can change the data */
        /* e.g., rotate into world coordinates */

        if (Pars.echo_data)  echo_data (&GEB_event);

        /* bin GT mode 3 data  (== raw data with traces) */
        if (Pars.do_bin_mode3 == 1) bin_mode3 (&GEB_event);

        /* bin GT mode 2 data  (== decomposed data) */
        if (Pars.do_bin_mode2 == 1) bin_mode2 (&GEB_event);

        /* bin mode 1 data (==tracked data) */
        if (Pars.do_bin_mode1 == 1) bin_mode1 (&GEB_event);

        /* bin DGS data */
        if (Pars.do_bin_dgs == 1) bin_dgs (&GEB_event);

        /* bin XA data */
        if (Pars.do_bin_XA == 1){
          bin_XA (&GEB_event);
          extern TTree *tree;
          tree->Fill();
        }
      };
      if (Pars.CurEvNo <= Pars.NumToPrint){
        printf ("*end of event # %i\n", Pars.CurEvNo);
        printf ("+++++++++++++++++++++++++++++++\n");
      };
    };


    /*---------------------*/
    /* house keeping...... */
    /* done every so often */
    /*---------------------*/

    if (Pars.CurEvNo % 100 == 0){
      /* calc time since last dump */
      tdmp = time (NULL);
      tdmp -= tdmplast;
      tdmp /= 60;           /* now minutes */
    };

    /*-----------------------------------------------------------*/
    /* dump all spectra on signal or dump every Pars.DumpEvery events */
    /* or respond to 'interactive' command...................... */
    /*-----------------------------------------------------------*/
    /*
    if (Pars.WeWereSignalled || (int) tdmp >= Pars.DumpEvery){

      // disarm signal 
      Pars.WeWereSignalled = FALSE;

      // check for command file 

      fp = fopen (CommandFileName, "r");
      if (fp != NULL){

        printf ("found command file: %s\n", CommandFileName);
        fgets (str, STRLEN, fp);
        printf ("with command: %s\n", str);

        if ((p = strstr (str, "DumpEvery")) != NULL){
          sscanf (str, "%s %i", str1, &Pars.DumpEvery);
          printf ("will dump to output file every %i minutes\n", Pars.DumpEvery);
          fflush (stdout);
        }else if ((p = strstr (str, "printevents")) != NULL){
          // reset print event counter 

          nret = sscanf (str, "%s %i", str1, &i1);
          if (nret == 2) Pars.NumToPrint = i1;
          printf ("will print %i events\n", Pars.NumToPrint);
          NprintEvNo = 0;

        }else if ((p = strstr (str, "status")) != NULL){

          showStatus ();
          
        }else if ((p = strstr (str, "stopsort")) != NULL){
          // simulate end of event to stop sort 

          eov = 1;

        }else if ((p = strstr (str, "zap")) != NULL){
          // extract spectrum name 

          sscanf (str, "%s %s", str1, Pars.spname);
          hhtemp = (TH1D *) gROOT->FindObject (Pars.spname);
        }else{
          printf ("command not understood\n");
        }
        // delete command file 

        fclose (fp);
        sprintf (str, "\\rm %s", CommandFileName);
        system (str);
        printf ("%s\n", str);

      }else{
        printf ("\"%s\" was not found\n", CommandFileName);

        // update sh mem or writeout root file 

        printf ("time since last dump: %i minute(s)\n", (int) tdmp);
        tdmp = 0;
        
        //if (!Pars.UseShareMemFile)
        {
          printf ("*---------------------------------\n");
          printf ("* you cannot update a disk file.  \n");
          printf ("  you must wait for sort to finish\n");
          printf ("  or stop the sort! Ignoring you...\n");
          printf ("*---------------------------------\n");
        }

        tdmplast = time (NULL);

      };

      printf ("continuing the sort...\n");
      fflush (stdout);

    };
  
    */
  };

  /*-----------------------*/
  /* we are done sorting!! */
  /* save all ROOT spectra */
  /*-----------------------*/

  printf ("\n");
  printf ("Sorting is done\n");
  printf ("attempting to save root or map file\n");
  printf ("\n");
  fflush (stdout);

  if (Pars.InputSrc == DISK){
      fclose (inData);
  }else if (Pars.InputSrc == GEB){
#if(HAVE_VXWORKS)
    gretTapClose (pTap);
#endif
  };
  printf ("\n");
  fflush (stdout);

  /* execute the exit scripts */

  if (Pars.do_bin_mode1) exit_bin_mode1 ();
  if (Pars.do_bin_mode2)  exit_bin_mode2 ();
  if (Pars.do_bin_mode3)  exit_bin_mode3 ();
  if (Pars.do_bin_dgs)   exit_dgs ();
  if (Pars.do_bin_XA)  exit_XA ();

  
  printf ("\nattempting to close root file...\n"); fflush (stdout);
  printf ("Pars.f1->Write();\n"); fflush (stdout);
  Pars.f1->Write ();
  printf ("Pars.f1->Print();\n"); fflush (stdout);
  Pars.f1->Print ();
  printf ("Pars.f1->Close();\n"); fflush (stdout);
  Pars.f1->Close ();
  printf ("done saving rootfile \"%s\n\n", Pars.ROOTFile); fflush (stdout);
  printf ("\n");

  /*-------------------------*/
  /* print simple statistics */
  /*-------------------------*/

  showStatus ();

  /* done */

  printf ("\n");
  printf ("sorted timestamp range %lli-->%lli: %lli\n", t0, tcur, tcur - t0);
  d1 = (double) (tcur - t0);
  d1 /= 100000000;
  printf ("that is %.1f seconds or ", (float) d1);
  nsec = d1;
  i1 = (unsigned int) d1 / 3600;
  d1 -= i1 * 3600;
  i2 = (unsigned int) d1 / 60;
  d1 -= i2 * 60;
  printf ("%ih%im%is\n", i1, i2, (int) d1);
  printf ("^^^^^ any timestamp jumps will upset this accounting\n");
  printf ("\n");

  printf ("hit statistics per type\n");
  i1 = 0;
  for (i = 1; i < MAX_GEB_TYPE; i++){
    if (typehit[i] > 0){
      GebTypeStr (i, str);
      printf ("%2i %s %10i ;", i, str, typehit[i]);
      i1 += typehit[i];
      d1 = (double) typehit[i] / nsec;
      printf (" %9.2f Hz ", (float) d1);
      printf ("\n");
    };
  }
  printf ("read a total of              %i ; header/payloads\n", i1);
  printf ("\n");

/*
  printf ("\n");
  printf ("nn1=%i\n", nn1);
  printf ("nn2=%i\n", nn2);
  printf ("nn3=%i\n", nn3);
  printf ("\n");

*/
  printf ("nok_validate= %lli, ntot_validate=%lli; ", nok_validate, ntot_validate);
  printf ("validate fraction %5.2f%%\n", 100 * (double) nok_validate / (double) ntot_validate);
  printf ("bytes read=%lli, %f MB\n", Pars.nbytes, (double) Pars.nbytes / 1024.0 / 1024.0);
  printf ("\n");

  printf ("boniva sancta! ");
//  printf ("...GEBSort (unexpectedly) did not crash!\n");
  printf ("\n ** GEBSort is done at ");
  time_stamp (stdout);
  printf ("\n");
  return (0);

}

/*----------------------------------------------------*/
