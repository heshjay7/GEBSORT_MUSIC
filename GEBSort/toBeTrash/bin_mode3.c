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
#include "GTmode3.h"

#define MAXINT   2147483647

#define DEBUGTWOCOMPLEMENT 1

/* pointers to ROOT spectra */

TH2F *eCC_a;
TH2F *eCC_b;
TH2F *eCC_c;
TH2F *eCC_d;
TH2F *eSeg;

TH1D *mode3_hitpat;
TH1D *mode3_hitpat_chan;

/* parameters */

extern PARS Pars;
extern EXCHANGE exchange;

/* ----------------------------------------------------------------- */

int
sup_mode3 ()
{
  /* declarations */

  char str1[STRLEN], str2[STRLEN];
  int i;

  TH1D *mkTH1D (char *, char *, int, double, double);
  TH2F *mkTH2F (char *, char *, int, double, double, int, double, double);
  int test_convert();

#if(0)
  for (i = 0; i < 20; i++)
    {
      printf ("%15u 0x%8.8x --> ", i, i);
      //printf ("%15i 0x%8.8x \n", twoscomp_to_int_24 ((unsigned int) i), twoscomp_to_int_24 ((unsigned int) i));
    }
  if (1)
    exit (0);
#endif

  /* initialize */

  /*----------------*/
  /* define spectra */

  sprintf (str1, "eCC_a");
  sprintf (str2, "eCC_a");
  eCC_a = mkTH2F (str1, str2, NCRYSTALS, 1, NCRYSTALS, ELENGTH, 1, ELENGTH);
  eCC_a->SetXTitle (str1);

  sprintf (str1, "eCC_b");
  sprintf (str2, "eCC_b");
  eCC_b = mkTH2F (str1, str2, NCRYSTALS, 1, NCRYSTALS, ELENGTH, 1, ELENGTH);
  eCC_b->SetXTitle (str1);

  sprintf (str1, "eCC_c");
  sprintf (str2, "eCC_c");
  eCC_c = mkTH2F (str1, str2, NCRYSTALS, 1, NCRYSTALS, ELENGTH, 1, ELENGTH);
  eCC_c->SetXTitle (str1);

  sprintf (str1, "eCC_d");
  sprintf (str2, "eCC_d");
  eCC_d = mkTH2F (str1, str2, NCRYSTALS, 1, NCRYSTALS, ELENGTH, 1, ELENGTH);
  eCC_d->SetXTitle (str1);

  sprintf (str1, "eSeg");
  sprintf (str2, "eSeg");
  eSeg = mkTH2F (str1, str2, NSEG, 1, NSEG, ELENGTH, 1, ELENGTH);
  eSeg->SetXTitle (str1);

  sprintf (str1, "mode3_hitpat");
  sprintf (str2, "mode3_hitpat");
  mode3_hitpat = mkTH1D (str1, str2, NCRYSTALS, 1, NCRYSTALS);
  mode3_hitpat->SetXTitle ("chystal ID");
  mode3_hitpat->SetYTitle ("counts");

  /* this one should be trivially flat  */
  /* because we read all segments out; but  */
  /* we have the spectrum just to verify this */

  sprintf (str1, "mode3_hitpat_chan");
  sprintf (str2, "mode3_hitpat_chan");
  mode3_hitpat_chan = mkTH1D (str1, str2, 45, 0, 44);
  mode3_hitpat_chan->SetXTitle ("channel ID");
  mode3_hitpat_chan->SetYTitle ("counts");

  /* list what we have */

//  Pars.wlist = gDirectory->GetList ();
//  Pars.wlist->Print ();

#if(0)
  /* test 32 bit to 24 bit converter */

  test_convert();
#endif

  return (0);

}

/* ----------------------------------------------------------------- */

int
exit_bin_mode3 ()
{

  return (0);

};


/* ----------------------------------------------------------------- */

int
bin_mode3 (GEB_EVENT * GEB_event)
{

  /* for some documentation, see: */
  /* evince ~/d6/keep/2017_GT_mode3_format.odp  */

  /* declarations */

  GTEVENT Event;
  unsigned int t1, t2, t3, t4;
  int ii, i, nBadTestPat = 0;
  volatile unsigned int *bit32Pointer;
  unsigned short int *bit16Pointer;
  FILE *fp1 = NULL, *fp2 = NULL;
  char str[128];
  int pos = 0, ncrystal = 0, sppos = 0, ccpos = 0;
  unsigned int tempE;
  int rawE;
  int gtheaders=0, newcrystal=0;

  /* prototypes */

  int GebTypeStr (int type, char str[]);
  int c24bit32bit (unsigned int);
//  void pprint_32 (char *, unsigned int);
//  int twoscomp_to_int_24 (unsigned int );

  if (Pars.CurEvNo <= Pars.NumToPrint)
    printf ("entered bin_mode3 @ event number: %i \n", Pars.CurEvNo);

  for (ii = 0; ii < GEB_event->mult; ii++)
    {

      /* pos keeps record of how far we have */
      /* proceeded in the payload */
      /* and count how many crystals in payload */

      pos = 0;
      ncrystal = 0;

      if (GEB_event->ptgd[ii]->type == GEB_TYPE_RAW || GEB_event->ptgd[ii]->type == GEB_TYPE_GT_MOD29)
        {

          if (Pars.CurEvNo <= Pars.NumToPrint)
            {
              GebTypeStr (GEB_event->ptgd[ii]->type, str);
              printf ("\nbin_mode3: %2i> %2i, %s, TS=%lli, 0x%llx; ", ii, GEB_event->ptgd[ii]->type, str,
                      GEB_event->ptgd[ii]->timestamp, GEB_event->ptgd[ii]->timestamp);
              printf ("payload length: %i bytes", GEB_event->ptgd[ii]->length);
              printf ("\n");
            }

          /* byteswap the entire payload */

          bit32Pointer = (unsigned int *) GEB_event->ptinp[ii];
          for (i = 0; i < GEB_event->ptgd[ii]->length / 4; i++)
            {
              t1 = (*(bit32Pointer + i) & 0x000000ff) << 24;
              t2 = (*(bit32Pointer + i) & 0x0000ff00) << 8;
              t3 = (*(bit32Pointer + i) & 0x00ff0000) >> 8;
              t4 = (*(bit32Pointer + i) & 0xff000000) >> 24;
              *(bit32Pointer + i) = t1 + t2 + t3 + t4;
            };

          /* inside the payload we have the a number */
          /* of header/trace data sets. These are the  */
          /* crystals that were in coincidence. Here  */
          /* we loop over these header/trace data sets. */
          /* We may have data from more than one crystal  */
          /* in this payload. Thus, there will be at  */
          /* least 40 traces, but there can also be  */
          /* 80, 120... etc. Be sure MAXPAYLOADSIZE  */
          /* is big enough to handle GT mode3 payloads. */

          /* loop over the header/traces of the mode 3 data */

          while (pos < GEB_event->ptgd[ii]->length)
            {

              /* start of event (Event.len known from last event) */

              if (pos == 0)
                bit32Pointer = (unsigned int *) GEB_event->ptinp[ii];
              else
                bit32Pointer += (Event.len / 4);

              /* check the EOE situation */

              if (*bit32Pointer != EOE)
                {
                  nBadTestPat++;

                  if (nBadTestPat == 10)
                    printf ("** suspending warnings about bad EOE markers...\n");

                  if (nBadTestPat < 10)
                    {
                      printf ("ooops: bit32Pointer=%8.8x after event # %i trace # %i\n", *bit32Pointer, Pars.CurEvNo,
                              ncrystal);
                      fflush (stdout);
                      exit (1);
                    };

                };

            if ( (gtheaders%40) == 0)
               {
               ncrystal++;
               newcrystal=1;
               }
             gtheaders++;
//             printf("%i %i %i\n", gtheaders, ncrystal, newcrystal);
                
              /* fill event header, skip EOE */

              bit16Pointer = (unsigned short int *) (bit32Pointer + 1);
              for (i = 0; i < HDRLENWORDS; i++)
                {
                  Event.hdr[i] = *(bit16Pointer + i);
                };

              /* debug list */

              if (Pars.CurEvNo <= Pars.NumToPrint)
                  {
                  printf("\n---------\nnew GTheader at pos=%i or %i words (%i):\n",pos,pos/2, pos/2+8);
                for (i = 0; i < HDRLENWORDS; i++)
                  printf ("Event.hdr[%2i]=0x%4.4x, %6i\n", i, Event.hdr[i], Event.hdr[i]);
                  };

              /* eventlength/tracelength in bytes */
              /* the +4 comes from the EOE 4 bit word */

              Event.len = 4 * (Event.hdr[1] & 0x7ff) + 4;
              Event.traceLen = Event.len - HDRLENBYTES - sizeof (unsigned int);

              if (Pars.CurEvNo <= Pars.NumToPrint)
                {
                  printf ("Event.len=%4i, Event.traceLen=%4i (in Bytes)\n", Event.len, Event.traceLen);
                  printf ("Event.len=%4i, Event.traceLen=%4i (in 16 bit words)\n", Event.len / sizeof (short int),
                          Event.traceLen / sizeof (short int));
                  printf ("Event.len=%4i, Event.traceLen=%4i (in 32 bit words)\n", Event.len / sizeof (int),
                          Event.traceLen / sizeof (int));
                };


              /* check the next EOE is there */
              /* before we go on (unless last set) */

              pos += Event.len;

              if (Pars.CurEvNo <= Pars.NumToPrint && pos < GEB_event->ptgd[ii]->length)
                printf ("next start: 0x%8.8x, pos= %i\n", *(bit32Pointer + (Event.len / 4)), pos);

              /* potential debug info */

              if (*(bit32Pointer + (Event.len / 4)) != EOE && pos < GEB_event->ptgd[ii]->length)
                {
                  printf ("next EOE not found at pointer %p, seek manually\n", bit32Pointer + (Event.len / 4));
                  printf ("Event.len=%i, Event.traceLen=%i\n", Event.len, Event.traceLen);
                  i = 1;
                  while (i < (Event.len + 10))
                    {
                      printf ("%3i has 0x%8.8x\n", 4 * i, *(bit32Pointer + i));
                      fflush (stdout);
                      if (*(bit32Pointer + i) == EOE)
                        {
                          printf ("found it\n");
                          break;
                        };
                      i++;
                    };
                  i--;
                  printf ("%3i has 0x%8.8x\n", i, *(bit32Pointer + i));
                };

              /* extract the Board IDs etc */
              /* chan_id is 0-9 on digitizer */

              Event.chan_id = (Event.hdr[0] & 0x000f);
              Event.digitizer_id = (Event.hdr[0] >> 4) & 0x0003;
              Event.crystal_id = (Event.hdr[0] >> 6) & 0x0003;
              Event.module_id = (Event.hdr[0] >> 8) & 0x001f;   /* zhu use 0x00ff */

              /* construct detector/channel/segment ID */

              Event.ge_id = Event.module_id * 4 + Event.crystal_id;
              mode3_hitpat->Fill (Event.ge_id, 1);

              Event.channel = Event.digitizer_id * 10 + Event.chan_id;
              mode3_hitpat_chan->Fill (Event.channel, 1);

              Event.seg_id = Event.ge_id * 40 + Event.channel;

              if (Pars.CurEvNo <= Pars.NumToPrint)
                {
                  printf ("0x%4.4x ", Event.hdr[0]);
                  printf ("module_id = %i ", Event.module_id);
                  printf ("crystal_id = %i ", Event.crystal_id);
                  printf ("digitizer_id = %i ", Event.digitizer_id);
                  printf ("chan_id = %i ", Event.chan_id);
                  printf ("ge_id = %i ", Event.ge_id);
                  printf ("channel = %i ", Event.channel);
                  if (((int) (Event.channel) % 10) == 9)
                    printf (" -> is CC");
                  else
                    printf (" -> is segment");
                  printf ("\n");
                };

              /* count the crystals we have processed */



                  if (Pars.CurEvNo <= Pars.NumToPrint)
                  if (newcrystal==1)
                    {
                      newcrystal=0;

                  /* write first few superpulse traces out */
                  /* open files here, fill later */

                      if (fp1 != NULL)
                        fclose (fp1);
                      if (fp2 != NULL)
                        fclose (fp2);

                      sprintf (str, "superpulse_%3.3i_%2.2i_%2.2i.xy", Pars.CurEvNo, ii, ncrystal);
                      printf ("open %s at Event.chan_id = %i\n", str, Event.chan_id);
                      fp1 = fopen (str, "w");
                      sppos = 0;

                      sprintf (str, "CCpulse_%3.3i_%2.2i_%2.2i.xy", Pars.CurEvNo, ii, ncrystal);
                      printf ("open %s at Event.chan_id = %i\n", str, Event.chan_id);
                      fp2 = fopen (str, "w");
                      ccpos = 0;

                    };

              if (Pars.CurEvNo <= Pars.NumToPrint)
                printf ("ncrystal=%i\n", ncrystal);

              /* copy trace into 16 bit words */
              /* bit16Pointer points first word after first EOE */

              for (i = 0; i < (Event.traceLen / 2); i++)
                if (((int) (Event.chan_id) % 10) == 9)
                  Event.trace[i] = *(bit16Pointer + i + HDRLENWORDS);
                else
                  {
                    Event.trace[i] = *(bit16Pointer + i + HDRLENWORDS);
                    Event.trace[i] = -(short int) Event.trace[i];
                  };


              if ((Pars.CurEvNo <= Pars.NumToPrint))
                if (fp1 != NULL && fp2 != NULL)
                  {
                    assert (fp1 != NULL);
                    if (((int) (Event.chan_id) % 10) != 9)
                      {
                        for (i = 0; i < (Event.traceLen / 2); i++)
                          fprintf (fp1, "%i %i\n", sppos++, Event.trace[i]);
                      }
                    else
                      {
                        for (i = 0; i < (Event.traceLen / 2); i++)
                          fprintf (fp2, "%i %i\n", ccpos++, Event.trace[i]);
                      }

                    printf ("wrote trace to %s\n", str);
                    fflush (fp1);
                    fflush (fp2);
                  };

              /* extract the energies like Shoufei does in muxTest.cc */
              /* needs some translating.... */

              tempE = (((unsigned int) Event.hdr[6]) & 0x00ff) << 16;
              tempE += (unsigned int) Event.hdr[5] & 0xffff;
              if (Pars.CurEvNo <= Pars.NumToPrint)
//                pprint_32 ("two's comp tempE  :: ", tempE);

//              rawE = twoscomp_to_int_24 (tempE);
//              rawE =c24bit32bit(tempE);
              if (Pars.CurEvNo <= Pars.NumToPrint)
 //               pprint_32 ("rawE              :: ", rawE);

              /* change the sign as the the Digs are neg */

              if (rawE <= 0)
                rawE = -rawE;
              if (Pars.CurEvNo <= Pars.NumToPrint)
 //                pprint_32 ("rawE (sign change):: ", rawE);

              /* downscale */

              Event.ehi = (int) ((float) rawE / Mvalue);
              if (Pars.CurEvNo <= Pars.NumToPrint)
                printf ("Event.ehi= %15i 0x%8.8x\n", Event.ehi, Event.ehi);

              /* for the energy extraction method from: */
              /* http://gretina.lbl.gov/tools-etc/GEBHeaderTypes.pdf */
              /* see svn version 733 */


              /* bin central contact (four gains) and segments  energies */

              Event.ehi /= 10;
              if (Pars.CurEvNo <= Pars.NumToPrint)
                printf ("Event.ehi=%6i (gain changed) for channel %2i\n", Event.ehi, Event.channel);

              if (Event.ehi > 10 && Event.ehi < ELENGTH)
                {

                  if ((Event.chan_id % 10) != 9)
                    if (Event.seg_id < NSEG)
                      eSeg->Fill ((double) Event.seg_id, (double) Event.ehi, 1);

                  if (Event.channel == 9)
                    eCC_a->Fill ((double) Event.ge_id, (double) Event.ehi, 1);
                  if (Event.channel == 19)
                    eCC_b->Fill ((double) Event.ge_id, (double) Event.ehi, 1);
                  if (Event.channel == 29)
                    eCC_c->Fill ((double) Event.ge_id, (double) Event.ehi, 1);
                  if (Event.channel == 39)
                    eCC_d->Fill ((double) Event.ge_id, (double) Event.ehi, 1);

                };

              /* extract LED external time, per documentation, works */

              Event.LEDts = (unsigned long long int) Event.hdr[2] +
                ((unsigned long long int) Event.hdr[3] << 16) + ((unsigned long long int) Event.hdr[4] << 32);
              if (Pars.CurEvNo <= Pars.NumToPrint)
                printf ("Event.LEDts = %20lli\n", Event.LEDts);

              /* extract CFD external time, per documentation, seems odd...not working? */

              Event.CFDts = (unsigned long long int) Event.hdr[7] +
                ((unsigned long long int) Event.hdr[8] << 16) + ((unsigned long long int) Event.hdr[9] << 32);
              if (Pars.CurEvNo <= Pars.NumToPrint)
                printf ("Event.CFDts = %20lli\n", Event.CFDts);


            };                  /* while (pos<=GEB_event->ptgd[ii]->length) */


        };                      /* if(GEB_event->ptgd[ii]->type == GEB_TYPE_RAW) */


    };


  /* done */

  if (Pars.CurEvNo <= Pars.NumToPrint)
    printf ("exit bin_mode3\n");

  return (0);

}
