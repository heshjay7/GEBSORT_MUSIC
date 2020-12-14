#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <math.h>
#include <string.h>
#include <assert.h>

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

#define SIMPLE 1


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


#if(SIMPLE)
TH1D *simple_CCsum;
TH1D *simple_CCadd;
#endif

/* parameters */

extern PARS Pars;
extern EXCHANGE exchange;


#define NBINS 180

#if(0)
/*-----------------------------------------------------------*/

TCutG *
rd2dwin (Char_t * winname)
{

  /* declarations */

  char str[STRLEN];
  TCutG *mycutg;

  TFile *f = new TFile (winname, "read");
  mycutg = (TCutG *) f->Get (winname);

  if (mycutg != NULL)
    {
      printf ("2dwin read from file:\n");
      sprintf (str, "ls -l %s", winname);
      gSystem->Exec (str);
      /* mycutg->Print(); */
      fflush (stdout);
    }
  else
    {
      mycutg = NULL;
      printf ("could not read 2dwin file %s\n", winname);
      printf ("TFile error number: %i\n", f->GetErrno ());
      printf ("abort\n\n");
      exit (-1);
    };

  /* done */

  f->Delete ();
  f->Close ();
  return (mycutg);

}
#endif


/*-----------------------------------------------------*/

int
sup_template ()
{
  /* declarations */

  char str1[STRLEN], str2[STRLEN];
  float pi;
  int i;
  TCutG *mywin;
  TH1D *mkTH1D (char *, char *, int, double, double);
  TH2F *mkTH2F (char *, char *, int, double, double, int, double, double);


#if(SIMPLE)
  sprintf (str1, "simple_CCadd");
  sprintf (str2, "simple_CCadd");
  simple_CCadd = mkTH1D (str1, str2, LONGLEN, 1, LONGLEN);
  sprintf (str1, "(keV)");
  simple_CCadd->SetXTitle (str1);

  sprintf (str1, "simple_CCsum");
  sprintf (str2, "simple_CCsum");
  simple_CCsum = mkTH1D (str1, str2, LONGLEN, 1, LONGLEN);
  sprintf (str1, "(keV)");
  simple_CCsum->SetXTitle (str1);
#endif

  /* list what we have */

  printf (" we have define the following spectra:\n");

#if(0)
  mywin = rd2dwin("test.win");
  mywin->Print();
if(1)exit(0);
#endif

  Pars.wlist = gDirectory->GetList ();
  Pars.wlist->Print ();



  return(0);
};

/* ----------------------------------------------------------------- */

int
bin_template (GEB_EVENT * GEB_event)
{

  /* declarations */

  char str[128];
  int i, j, k, nn, isGGLOW, isGGHI;
  int ee[1];
  float xx[10], yy[10], zz[10];
  static float xx_o[10], yy_o[10], zz_o[10];
  TRACKED_GAMMA_HIT *grh;
  float dd, dotProduct, polang;
  CRYS_INTPTS *ptinp;
  float simple_add=0;

  /* prototypes */

  int GebTypeStr (int type, char str[]);

//if(1)return(0);

  if (Pars.CurEvNo <= Pars.NumToPrint)
    printf ("entered bin_template:\n");
 
  /*------------------------------------------*/
  /* generate very simple add and sum spectra */
  /* for debugging really, works on mode2 CC  */
  /*------------------------------------------*/

#if(SIMPLE)

  simple_add=0;
  for (i = 0; i < GEB_event->mult; i++)
    {
      if (GEB_event->ptgd[i]->type == GEB_TYPE_DECOMP)
        {

          /* cast */

          ptinp = (CRYS_INTPTS *) GEB_event->ptinp[i];

         if (ptinp->pad == 0)
            {
            simple_CCsum->Fill ((double) ptinp->tot_e , 1.0);
            simple_add+=ptinp->tot_e;
            };

         };
      };
   simple_CCadd->Fill (simple_add,1.0);

#endif



  /* done */

  if (Pars.CurEvNo <= Pars.NumToPrint)
    printf ("exit bin_template\n");

  return (0);

}
