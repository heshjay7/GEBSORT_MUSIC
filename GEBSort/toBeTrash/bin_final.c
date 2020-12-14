#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <limits.h>

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

#include "veto_pos.h"
#include "GEBSort.h"
#include "GTMerge.h"

#define SIMPLE 0


typedef struct PAYLOAD
{
  char p[MAXDATASIZE];
} PAYLOAD;

typedef struct TRACK_STRUCT
{
  int n;
  GEBDATA *gd;
  PAYLOAD *payload;
} TRACK_STRUCT;


/* pointers to ROOT spectra */

TH2F *firstType;
TH2F *dgsdfma_hit;

/* parameters */

extern PARS Pars;
extern EXCHANGE exchange;

/* common storage of bin_xxx analysis */

extern int ng;
extern DGSEVENT DGSEvent[MAXCOINEV];
extern DFMAEVENT DFMAEvent[MAXCOINEV];

/*-----------------------------------------------------*/

int
sup_final ()
{
  /* declarations */

  char str1[128], str2[128];

  /* prototypes */

  TH1D *mkTH1D (char *, char *, int, double, double);
  TH2F *mkTH2F (char *, char *, int, double, double, int, double, double);

  /* root spectra */

  sprintf (str1, "firstType");
  sprintf (str2, "first type vs multiplicity");
  firstType = mkTH2F (str1, str2, MAX_GEB_TYPE + 1, 1, MAX_GEB_TYPE + 1, 100, 1, 100);
  sprintf (str1, "first type");
  firstType->SetXTitle (str1);
  sprintf (str1, "event mult");
  firstType->SetYTitle (str1);

  sprintf (str1, "dgsdfma_hit");
  sprintf (str2, "dfma vs dgs header numbers");
  dgsdfma_hit = mkTH2F (str1, str2, 100, 0, 99, 100, 0, 99);
  dgsdfma_hit->SetXTitle ("ndgs");
  dgsdfma_hit->SetYTitle ("event mult");


  printf ("spectra define in bin_final:\n");

  Pars.wlist = gDirectory->GetList ();
  Pars.wlist->Print ();

  return (0);
};

/*-----------------------------------------------------*/

int
exit_final ()
{
  /* declarations */

};

/* ----------------------------------------------------------------- */

/* The idea of bin_final is that it takes results from other  */
/* bin_functions.c where these functions have deposited their  */
/* individual information in the exchange structure and makes  */
/* the final cross coincidence evaluations. It shold not itself  */
/* really look at the datastream. */

int
bin_final (GEB_EVENT * GEB_event)
{

  /* declarations */

  int i, j, mintype, t1, t2, ndgs, ndfma;
  long long int minTS;
  int dts[100], tp[100];
  char str[20];
  int GebTypeStr (int type, char str[]);

  /* prototypes */


  if (Pars.CurEvNo <= Pars.NumToPrint)
    printf ("entered bin_final:\n");

#if(1)

  /* some trigger statistics */

  minTS = LLONG_MAX;
  mintype = 0;
  ndgs=0;
  ndfma=0;
  for (i = 0; i < GEB_event->mult; i++)
    {
      if (GEB_event->ptgd[i]->type == GEB_TYPE_DGS) ndgs++;
      if (GEB_event->ptgd[i]->type == GEB_TYPE_DFMA) ndfma++;
      if (GEB_event->ptgd[i]->timestamp < minTS)
        {
          minTS = GEB_event->ptgd[i]->timestamp;
          mintype = GEB_event->ptgd[i]->type;
//          printf ("xx%2i: TS= %lli, mintype=%i\n", i, minTS, mintype);
        };
    };
  assert (mintype > 0);
  if (GEB_event->mult < 100 && mintype > 0)
    firstType->Fill (mintype, GEB_event->mult);
  dgsdfma_hit->Fill(ndgs,ndfma);

  if (Pars.CurEvNo <= Pars.NumToPrint)
    {
      printf ("mintype,GEB_event->mult=%i %i\n", mintype, GEB_event->mult);

      /* store */

      for (i = 0; i < GEB_event->mult; i++)
        {
          printf ("yy%2i: TS= %lli, type=%i\n", i,GEB_event->ptgd[i]->timestamp , GEB_event->ptgd[i]->type);
          dts[i] = (int) (GEB_event->ptgd[i]->timestamp - minTS);
          printf("dts: %i\n", dts[i]);
          tp[i] = GEB_event->ptgd[i]->type;
        }

      /* sort */

      for (i = 0; i < GEB_event->mult; i++)
        for (j = i + 1; j < GEB_event->mult; j++)
          if (dts[i] > dts[j])
            {
              t1 = dts[i];
              t2 = tp[i];
              dts[i] = dts[j];
              tp[i] = tp[j];
              dts[j] = t1;
              tp[j] = t2;
            }

      /* list */

      for (i = 0; i < GEB_event->mult; i++)
        {
//          printf ("zz%2i: TS= %lli, type=%i\n", i,GEB_event->ptgd[i]->timestamp , GEB_event->ptgd[i]->type);
        printf(" @ dt=%4i, tp=%2i ", dts[i], tp[i]);
        if (tp[i]==14)
          printf("DGS  ");
        else if
          (tp[i]==16)
          printf("DFMA ");
        };
      printf("\n");

    };

#endif


  /* done */

  if (Pars.CurEvNo <= Pars.NumToPrint)
    printf ("exit bin_final\n");

  return (0);

}
