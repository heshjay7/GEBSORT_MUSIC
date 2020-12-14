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


/* parameters */

extern PARS Pars;


/* ----------------------------------------------------------------- */

int
validate (GEB_EVENT * GEB_event)
{

  /* declarations */

  char str[128];
  int i, j, k, rval = 0, nGod, nBad;
  CRYS_INTPTS *ptinp;
  GEBDATA *ptgd;
  int crystalno, moduleno, detno;
  static int nn = 0;
  TRACKED_GAMMA_HIT *grh;
  float r1, esum;
  int nmode1, nmode2, emode2;

  /* prototypes */

  int GebTypeStr (int type, char str[]);


  /*-------------------------------------*/
  /* Validate the trigger conditions     */
  /* we allways do that in this function */
  /*-------------------------------------*/

  emode2 = 0;
  for (i = 0; i < GEB_event->mult; i++)
    {

//      if (GEB_event->ptgd[i]->type == GEB_TYPE_DECOMP)
//        {
//          /* count how many CC energies are greater than threshold */
//          /* soo we can minic the GT trigger system */
//          /* we assume one crystal can only have one header */
//
//          ptinp = (CRYS_INTPTS *) GEB_event->ptinp[i];
//          if (ptinp->tot_e > Pars.minCCe)
//            emode2++;
//
//        };

    };

//  if (emode2<Pars.minNumCC) return(0);
//  if (emode2>Pars.maxNumCC) return(0);






  /* trivial return. Ignores anything below */

  return (1);




  if (Pars.CurEvNo <= Pars.NumToPrint)
    printf ("entered validate:\n");

  /* put a check on the event multiplicity */
  /* to avoid problems in the bin_ functions */

  if (GEB_event->mult >= MAX_GAMMA_RAYS)
    return (0);












#if(0)

  /* kill events where there was more than one interaction */
  /* in a segment in the coincidence event */

  for (i = 0; i < GEB_event->mult; i++)
    {

      if (GEB_event->ptgd[i]->type == GEB_TYPE_DECOMP)
        {

          if (Pars.CurEvNo <= Pars.NumToPrint)
            {
              GebTypeStr (GEB_evsp CCsum ent->ptgd[i]->type, str);
              printf ("validate, %2i> %2i, %s, TS=%lli\n", i, GEB_event->ptgd[i]->type, str,
                      GEB_event->ptgd[i]->timestamp);
            }

          /* cast */

          ptinp = (CRYS_INTPTS *) GEB_event->ptinp[i];



          /* loop over interaction points */

          for (j = 0; j < ptinp->num; j++)
            for (k = j + 1; k < ptinp->num; k++)
              if (ptinp->intpts[j].seg == ptinp->intpts[k].seg)
                {
                  if (nn < 100)
                    {
                      printf ("ptinp->intpts[%i].seg=%i\n", j, ptinp->intpts[j].seg);
                      printf ("ptinp->intpts[%i].seg=%i\n", k, ptinp->intpts[k].seg);
                      printf ("validate: found two hits in same segment\n");
                    }
                  nn++;
                  return (0);
                };

        };

    };

  /* if we did not kickout above then return good */

  if (nn < 100)
    printf ("validate: no double hits in same segment\n");
  return (1);

#endif





















#if(0)

  /* external detector triggering */

  rval = 0;
  for (i = 0; i < GEB_event->mult; i++)
    {

      if (GEB_event->ptgd[i]->type == GEB_TYPE_DECOMP)
        {

          if (Pars.CurEvNo <= Pars.NumToPrint)
            {
              GebTypeStr (GEB_event->ptgd[i]->type, str);
              printf ("validate, %2i> %2i, %s, TS=%lli\n", i, GEB_event->ptgd[i]->type, str,
                      GEB_event->ptgd[i]->timestamp);
            }

          /* extract ID and energy */

          /* cast */

          ptinp = (CRYS_INTPTS *) GEB_event->ptinp[i];

          /* extract IDs */

          crystalno = (ptinp->crystal_id & 0x0003);
          moduleno = ((ptinp->crystal_id & 0xfffc) >> 2);

          if (Pars.AGATA_data == 0)
            {
              detno = moduleno * 4 + crystalno;
            }
          else if (Pars.AGATvalidate.cA_data == 1)
            {
              detno = moduleno * 3 + crystalno;
            };

          if (Pars.CurEvNo <= Pars.NumToPrint)
            {
              GebTypeStr (GEB_event->ptgd[i]->type, str);
              printf ("validate, %2i> %2i, %s, TS=%lli\n", i, GEB_event->ptgd[i]->type, str,
                      GEB_event->ptgd[i]->timestamp);
            }

          if (detno == 33)
            {
              if (Pars.CurEvNo <= Pars.NumToPrint)
                printf ("found external det with energy %f \n", ptinp->tot_e);
//             if ( (ptinp->tot_e > 1170.2 ) && (ptinp->tot_e < 1176.6 ) )
              if ((ptinp->tot_e > (1173 - 2)) && (ptinp->tot_e < (1173 + 2)))
//             if ( (ptinp->tot_e > 1180.0 ) && (ptinp->tot_e < 1186.4 ) )
                {
                  rval = 1;
                  ptinp->totvalidate.c_e = 0;
                };
            };

        };

    };

#endif





#if(1)

  /* require to see a specific summed energy in the array */

  rval = 0;
  esum = 0;
  for (i = 0; i < GEB_event->mult; i++)
    {

      if (GEB_event->ptgd[i]->type == GEB_TYPE_DECOMP)
        {

          if (Pars.CurEvNo <= Pars.NumToPrint)
            {
              GebTypeStr (GEB_event->ptgd[i]->type, str);
              printf ("validate, %2i> %2i, %s, TS=%lli\n", i, GEB_event->ptgd[i]->type, str,
                      GEB_event->ptgd[i]->timestamp);
            }

          /* extract ID and energy */

          /* cast */

          ptinp = (CRYS_INTPTS *) GEB_event->ptinp[i];

          /* sum up to calorimetric energy */

          esum += ptinp->tot_e;

        };

    };

  /* evaluate */

//  printf("esum=%7.1f\n",esum);
  if (abs (esum - 1706.3) < 2.5)
    rval = 1;

#endif



















#if(0)

  /* only valiate if a particular tracked energy is seen */
  /* that way we can filter on summed peaks */
  /* or good photopeaks */

  nGod = 0;
  nBad = 0;
  for (i = 0; i < GEB_event->mult; i++)
    {

      if (GEB_event->ptgd[i]->type == GEB_TYPE_TRACK)
        {
          grh = (TRACKED_GAMMA_HIT *) GEB_event->ptinp[i];

          if (Pars.CurEvNo <= Pars.NumToPrint || nn < 200)
            printf ("ngam=%i\n", grh->ngam);
          for (j = 0; j < grh->ngam; j++)
            {
              if (Pars.CurEvNo <= Pars.NumToPrint || nn < 200)
                printf ("%2i: esum=%8.0f, fom=%4.2f, tracked %i\n", j, grh->gr[j].esum, grh->gr[j].fom,
                        grh->gr[j].tracked);

              /* inspect energies */

              /* photopeaks */

              r1 = abs (183.8 - grh->gr[j].esum);
              if (r1 < 3.0)
                nGod++;
              r1 = abs (711.7 - grh->gr[j].esum);
              if (r1 < 3.0)
                nGod++;
              r1 = abs (810.2 - grh->gr[j].esum);
              if (r1 < 3.0)
                nGod++;

              /* summed peaks */

              r1 = abs (895.5 - grh->gr[j].esum);
              if (r1 < 3.0)
                nBad++;
              r1 = abs (994.3 - grh->gr[j].esum);
              if (r1 < 3.0)
                nBad++;
              r1 = abs (1705.7 - grh->gr[j].esum);
              if (r1 < 3.0)
                nBad++;
              r1 = abs (1522.3 - grh->gr[j].esum);
              if (r1 < 3.0)
                nBad++;


            };

          if (Pars.CurEvNo <= Pars.NumToPrint || nn < 200)
            {
              printf ("nGod=%i, nBad=%i\n", nGod, nBad);
              if (nGod == grh->ngam)
                printf ("return (1);\n");
              else
                printf ("return (0);\n");
            }
          nn++;

          if (nGod == grh->ngam)
//          if (nBad == grh->ngam)
            return (1);
          else
            return (0);


        };

    };


  /* return */



#endif

  /* done */


  if (Pars.CurEvNo <= Pars.NumToPrint)
    printf ("exit validate\n");

  return (rval);

}
