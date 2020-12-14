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
#include "GEBFilter.h"
#include "veto_pos.h"

#define MAXGTMODNO 30

off_t inData, outData;
unsigned int *veto_cube;
GEBMERGEPAR Pars;
FILE *afp;

float crmat[MAXDETPOS + 1][MAXCRYSTALNO + 1][4][4];
double TrX[NUMAGATAPOS], TrY[NUMAGATAPOS], TrZ[NUMAGATAPOS];
double rotxx[NUMAGATAPOS], rotxy[NUMAGATAPOS], rotxz[NUMAGATAPOS];
double rotyx[NUMAGATAPOS], rotyy[NUMAGATAPOS], rotyz[NUMAGATAPOS];
double rotzx[NUMAGATAPOS], rotzy[NUMAGATAPOS], rotzz[NUMAGATAPOS];

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

/*---------------------------------------------------------------------------------*/

int
readChatFile (char *name)
{

  /* declarations */

  FILE *fp, *fp1, *fp0;
  char *p, *pc, str[STRLEN], str1[STRLEN], buffer[512];;
  int echo, nn = 0, nni = 0, nret, in, siz;
  int i, j, zero = 0;
  int ir[NUMAGATAPOS], dummy_i[NUMAGATAPOS];

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
      else if ((p = strstr (str, "nevents")) != NULL)
        {
          nret = sscanf (str, "%s %i", str1, &Pars.nevents);
          CheckNoArgs (nret, 2, str);
          printf ("will process %i events\n", Pars.nevents);
        }
      else if ((p = strstr (str, "waitusec")) != NULL)
        {
          nret = sscanf (str, "%s %i", str1, &Pars.waitusec);
          CheckNoArgs (nret, 2, str);
          printf ("will sleep %i usec between events\n", Pars.waitusec);
        }
      else if ((p = strstr (str, "vetocube")) != NULL)
        {
          nret = sscanf (str, "%s %s", str1, Pars.vetocubefn);
          CheckNoArgs (nret, 2, str);
          Pars.vetocubes = 1;
          printf ("will read %s\n", Pars.vetocubefn);
          printf ("  and process veto_cube\n");
          printf ("\n");
        }
      else if ((p = strstr (str, "addT0")) != NULL)
        {
          printf ("will do TS --> TS*10+(long long int)(T0+0.5)\n");
          printf ("for mode2 data and just\n");
          printf ("        TS --> TS*10\n");
          printf ("for any other GEBHEader type\n");
          Pars.addT0 = 1;
        }
      else if ((p = strstr (str, "GT2AGG4")) != NULL)
        {
          nret = sscanf (str, "%s %s", str1, Pars.GT2AGG4_fn);
          Pars.GT2AGG4 = 1;
          printf ("will write GT data in AGATA G4 ascii format\n");
          afp = fopen (Pars.GT2AGG4_fn, "w");
          if (afp != NULL)
            printf ("__\"%s\" is open for write\n", Pars.GT2AGG4_fn);
          fprintf(afp,"$\n");

          /* read the GT rotational matrices in */

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

          printf ("read in AGATA rotation and translation data\n");

          fp0 = fopen ("GANIL_AGATA_crmat.dat", "r");
          if (fp0 != NULL)
            {
              printf ("GANIL_AGATA_crmat.dat is open for reading\n");

              j = 0;
              for (i = 0; i < 180; i++)
                {

                  memset (buffer, zero, sizeof (buffer));
                  fgets (buffer, 150, fp0);
                  sscanf (buffer, "%d %d %lf %lf %lf ", &ir, &dummy_i, &TrX[j], &TrY[j], &TrZ[j]);

                  memset (buffer, zero, sizeof (buffer));
                  fgets (buffer, 150, fp0);
                  sscanf (buffer, "%d %lf %lf %lf  ", &dummy_i, &rotxx[j], &rotxy[j], &rotxz[j]);

                  memset (buffer, zero, sizeof (buffer));
                  fgets (buffer, 150, fp0);
                  sscanf (buffer, "%d %lf %lf %lf  ", &dummy_i, &rotyx[j], &rotyy[j], &rotyz[j]);

                  memset (buffer, zero, sizeof (buffer));
                  fgets (buffer, 150, fp0);
                  sscanf (buffer, "%d %lf %lf %lf  ", &dummy_i, &rotzx[j], &rotzy[j], &rotzz[j]);

                  j++;
                };
              printf ("read %i AGATA rotational/translational matrices\n", j);
              fclose (fp0);
            };

        }
      else if ((p = strstr (str, "xyz_smear")) != NULL)
        {
          nret = sscanf (str, "%s %f %f %f", str1, &Pars.smear_x, &Pars.smear_y, &Pars.smear_z);
          CheckNoArgs (nret, 4, str);
          printf ("will smear x,y,z with +/- (%5.2f,%5.2f,%5.2f) mm\n", Pars.smear_x, Pars.smear_y, Pars.smear_z);
          Pars.xyz_smear = 1;
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

  return (0);

};

/*---------------------------------------------------------------------*/

int
main (int argc, char **argv)
{

  /* declarations */

  GEBDATA *ptgd;
  PAYLOAD *ptinp;
  int i, siz, ip, nb, j, holeNum, crystalNumber;
  long long int ngood = 0, nbad = 0, last_timestamp = 0;
  int haveTrack = 0;
  int nread = 0, nwrite = 0;
  int writeOK;
  int curev = 0;
  CRYS_INTPTS *cip, tmpcip;
  int bip[MAX_INTPTS];
  int i1, i2, i3, moduleno, crystalno;
  double d1, rr, min, r1, r2, r3, ee;
  float ebad, esum, egood, ff;
  float rn, xx, yy, zz;
  unsigned int seed;
  int AG_crystal_no[100], detectornumber, evno=0;

  /* prototypes */

  int setup_veto_cube ();

  /* help */

  if (argc != 4)
    {
      printf ("use: %s chatfile infile outfile\n", argv[0]);
      exit (0);
    };

  /* initialize */

  bzero ((void *) &Pars, sizeof (GEBMERGEPAR));
  Pars.xyz_smear = 0;

  /* initialize random number generator */

  get_a_seed (&seed);
  srand (seed);

  /* read chat file */

  readChatFile (argv[1]);

  /* allocate space for track structure */

  ptgd = (GEBDATA *) calloc (1, sizeof (GEBDATA));
  ptinp = (PAYLOAD *) calloc (1, sizeof (PAYLOAD));

  /* open data file */

  inData = 0;
  /*                        + name of data file */
  /*                        |                   */
  inData = open ((char *) argv[2], O_RDONLY);
  if (inData == 0)
    {
      printf ("could not open input data file %s, quit!\n", argv[2]);
      exit (1);
    };
  printf ("input  data file %s is open\n", argv[2]);

  /* open output file */

  outData = 0;
  /*                        + name of data file */
  /*                        |                   */
  outData = open ((char *) argv[3], O_WRONLY | O_CREAT | O_TRUNC, PMODE);
  if (outData == 0)
    {
      printf ("could not open output data file %s, quit!\n", argv[3]);
      exit (1);
    };
  printf ("output data file %s is open\n", argv[3]);

  /* setup */

  if (Pars.vetocubes)
    setup_veto_cube ();

  if (Pars.GT2AGG4)
    {

      /* find AGATA norm vectors to center of crystals */

      for (i = 0; i <= 180; i++)
        {
          rr = TrX[i] * TrX[i] + TrY[i] * TrY[i] + TrZ[i] * TrZ[i];
          rr = sqrt (rr);
          TrX[i] /= rr;
          TrY[i] /= rr;
          TrZ[i] /= rr;
        };

    };



  /* loop through the data and filter */

  while (1)
    {

      /* read GEB header and payload */

      siz = read (inData, (char *) ptgd, sizeof (GEBDATA));
      if (siz != sizeof (GEBDATA))
        goto errordone;

      siz = read (inData, (char *) ptinp, ptgd->length);
      if (siz != ptgd->length)
        goto errordone;

      nread++;
      if (nread < 20)
        printf ("evno %4i, type %i TS= %lli\n", nread, ptgd->type,ptgd->timestamp);

      writeOK = 1;


      /* write GT data in AGATA G4 format */

      if (Pars.GT2AGG4)
        {

          if (ptgd->type == GEB_TYPE_DECOMP)
            {

              /* cast */

              cip = (CRYS_INTPTS *) ptinp;

              holeNum = ((cip->crystal_id & 0xfffc) >> 2);
              crystalNumber = (cip->crystal_id & 0x0003);
              detectornumber = 4 * holeNum + crystalNumber;

              if (nread < 20)
                {
                  printf ("-----------\n");
                  printf ("holeNum=%i, crystalNumber=%i, detectornumber=%i\n", holeNum, crystalNumber, detectornumber);
                  for (ip = 0; ip < cip->num; ip++)
                    {
                      printf ("[%3i]mode2: (%5.1f ", nread, cip->intpts[ip].x);
                      printf ("%5.1f ", cip->intpts[ip].y);
                      printf ("%5.1f) ", cip->intpts[ip].z);
                      printf ("%5.1f ", cip->intpts[ip].e);
                      printf (" ID(raw) %3i ", cip->crystal_id);
                      printf (" TS %3lli ", cip->timestamp);
                      printf ("\n");
                    };
                };

              /* is this a new event? */

              if ((cip->timestamp - last_timestamp) > (long long int) 30)
                {

                  /* write G4 header out since it must be a new event */

                  ee = 1000;
                  r1 = 1.0;
                  r2 = 1.0;
                  r3 = 1.0;
                  fprintf (afp, "-1 %f %f %f %f %i\n", ee, r1, r2, r3, evno);
                  evno++;

                };
              last_timestamp = cip->timestamp;

              /* rotate into world coordinates */

              for (j = 0; j < cip->num; j++)
                {


                  xx = cip->intpts[j].x / 10;
                  yy = cip->intpts[j].y / 10;
                  zz = cip->intpts[j].z / 10;


                  cip->intpts[j].x = crmat[holeNum][crystalNumber][0][0] * xx
                    + crmat[holeNum][crystalNumber][0][1] * yy
                    + crmat[holeNum][crystalNumber][0][2] * zz + crmat[holeNum][crystalNumber][0][3];

                  cip->intpts[j].y = crmat[holeNum][crystalNumber][1][0] * xx
                    + crmat[holeNum][crystalNumber][1][1] * yy
                    + crmat[holeNum][crystalNumber][1][2] * zz + crmat[holeNum][crystalNumber][1][3];

                  cip->intpts[j].z = crmat[holeNum][crystalNumber][2][0] * xx
                    + crmat[holeNum][crystalNumber][2][1] * yy
                    + crmat[holeNum][crystalNumber][2][2] * zz + crmat[holeNum][crystalNumber][2][3];

                  cip->intpts[j].x *= 10;
                  cip->intpts[j].y *= 10;
                  cip->intpts[j].z *= 10;

                  if (nread < 20)
                    {
                      printf ("crmat: hit no (j): %i of %i\n", j, cip->num);
                      printf ("crystal coords: %9.4f %9.4f %9.4f \n", xx, yy, zz);
                      printf ("world coords: %9.4f %9.4f %9.4f (mm)\n", cip->intpts[j].x, cip->intpts[j].y,
                              cip->intpts[j].z);
                    };

                  /* find interaction point AG crystal number */

                  r1 =
                    cip->intpts[j].x * cip->intpts[j].x + cip->intpts[j].y * cip->intpts[j].y +
                    cip->intpts[j].z * cip->intpts[j].z;
                  r1 = sqrt (r1);
                  if (nread < 20)
                    {
                      printf ("finding nearest AG crystal for interaction point\n");
                      printf ("__ %f %f %f\n", cip->intpts[j].x / r1, cip->intpts[j].y / r1, cip->intpts[j].z / r1);
                    };
                  min = 1000000000.0;
                  for (i = 0; i < 180; i++)
                    {
                      rr =
                        cip->intpts[j].x / r1 * TrX[i] + cip->intpts[j].y / r1 * TrY[i] +
                        cip->intpts[j].z / r1 * TrZ[i];
                      if (nread < 20)
                        {
//                        printf("AG # %3i norm vec %f %f %f\n",i, TrX[i],TrY[i],TrZ[i]);
                        };
                      rr = acos (rr);
                      if (rr < min)
                        {
                          min = rr;
                          AG_crystal_no[j] = i;
                        };
                    };

                  if (nread < 20)
                    {
                      printf ("GT # %i, AG # %i\n", detectornumber, AG_crystal_no[j]);
                    }


                };

              /* rotate back into AGATA crystal coordinates */
              /* NOT NECESSARY: OFT expects world cooordinates */

              /* write out in AGATA G4 ascii format */

              for (j = 0; j < cip->num; j++)
                {
                  fprintf (afp, "%i %f %f %f %f %i\n", AG_crystal_no[j], cip->intpts[j].e, cip->intpts[j].x,
                           cip->intpts[j].y, cip->intpts[j].z, 1);
                  if (nread < 20)
                    printf ("%i %f %f %f %f %i\n", AG_crystal_no[j], cip->intpts[j].e, cip->intpts[j].x,
                            cip->intpts[j].y, cip->intpts[j].z, 1);
                };



            };


          /* do not filter to the regular output file */

          writeOK = 0;

        };

      /* T0 adjustment ? */

      if (Pars.addT0)
        {

          /* make the TS ns rather than 10 ns */
          /* add allow for mode2 T0 adjustment */

          if (curev < 100)
            printf ("TS= %lli, ", ptgd->timestamp);
          ptgd->timestamp *= 10;

          /* if mode2 data, add T0 to nearest ns */

          if (ptgd->type == GEB_TYPE_DECOMP)
            {

              /* cast */

              cip = (CRYS_INTPTS *) ptinp;
              if (curev < 100)
                printf ("cip->t0=%f --> ", cip->t0);

              /* add T0 to the nearest 1 ns */

              ptgd->timestamp += (long long int) (10.0 * cip->t0 + 0.5);

            };
          if (curev < 100)
            printf ("%lli\n", ptgd->timestamp);

        };

      /* smear the xyz data */

      if (Pars.xyz_smear)
        {
          if (ptgd->type == GEB_TYPE_DECOMP)
            {

              /* cast */

              cip = (CRYS_INTPTS *) ptinp;
              assert (cip->num < MAX_INTPTS);

              for (ip = 0; ip < cip->num; ip++)
                {
                  rn = Pars.smear_x * (drand48 () - 0.5);
                  cip->intpts[ip].x += rn;
                  rn = Pars.smear_y * (drand48 () - 0.5);
                  cip->intpts[ip].y += rn;
                  rn = Pars.smear_z * (drand48 () - 0.5);
                  cip->intpts[ip].z += rn;
                };

            };

        };

      /* remove bad interaction */

      if (Pars.vetocubes)
        if (ptgd->type == GEB_TYPE_DECOMP)
          {

            /* cast */

            cip = (CRYS_INTPTS *) ptinp;
            assert (cip->num < MAX_INTPTS);

            /* init marker array */

            for (ip = 0; ip < cip->num; ip++)
              bip[ip] = 0;


            /* loop over interactions */

            ebad = 0;
            esum = 0;
            for (ip = 0; ip < cip->num; ip++)
              {
//                printf ("[]%i %i\n", cip->num, tmpcip.num);
                fflush (stdout);

                nb = 0;

                /* find vetocube index */

                crystalno = (cip->crystal_id & 0x0003);
                moduleno = ((cip->crystal_id & 0xfffc) >> 2);
                i1 = VETO_X_INDEX (cip->intpts[ip].x);
                i2 = VETO_Y_INDEX (cip->intpts[ip].y);
                i3 = VETO_Z_INDEX (cip->intpts[ip].z);
                bip[ip] = *(veto_cube + VETO_INDX (moduleno, crystalno, i1, i2, i3));
                assert (cip->num < MAX_INTPTS);

                if (curev < 10)
                  {
                    printf ("%i: %6.2f %6.2f %6.2f --> bip=", ip, cip->intpts[ip].x, cip->intpts[ip].y,
                            cip->intpts[ip].z);
                    printf ("%i\n", bip[ip]);
                    fflush (stdout);
                  };

                /* accounting */

                esum += cip->intpts[ip].e;
                if (bip[ip] == 0)
                  ngood++;
                else
                  {
                    nb++;
                    nbad++;
                    ebad += cip->intpts[ip].e;
                  };
              };
            egood = esum - ebad;
            ff = ebad / egood;
            assert (cip->num < MAX_INTPTS);


            /* zap the bad interaction points */

            if (nb > 0)
              {
                bcopy ((char *) cip, (char *) &tmpcip, sizeof (CRYS_INTPTS));
//                printf ("[]%i %i\n", cip->num, tmpcip.num);
//                fflush (stdout);
                assert (tmpcip.num == cip->num);
                assert (tmpcip.num < MAX_INTPTS);

                /* transfer only valid points back */
                /* redistribute the bad energy */

                j = 0;
                assert (cip->num < MAX_INTPTS);
                for (ip = 0; ip < cip->num; ip++)
                  {
                    if (bip[ip] == 0)
                      {
                        if (curev < 10)
                          {
                            printf ("ip=%i\n", ip);
                            fflush (stdout);
                            printf ("tmpcip.intpts[ip].x=%f\n", tmpcip.intpts[ip].x);
                            fflush (stdout);
                            printf ("j=%i\n", j);
                            fflush (stdout);
                            printf ("cip->intpts[j].x=%f\n", cip->intpts[j].x);
                            fflush (stdout);
                          };

                        cip->intpts[j].x = tmpcip.intpts[ip].x;
                        cip->intpts[j].y = tmpcip.intpts[ip].y;
                        cip->intpts[j].z = tmpcip.intpts[ip].z;
                        cip->intpts[j].e = tmpcip.intpts[ip].e * (1 + ff);
                        cip->intpts[j].seg = tmpcip.intpts[ip].seg;
                        cip->intpts[j].seg_ener = tmpcip.intpts[ip].seg_ener;
                        j++;
                      };
                  };


                /* update number of interaction points */

                cip->num = j;

                /* prevent writeout if none survived */

                if (cip->num == 0)
                  writeOK = 0;

              };

          };


      /* filter */

      /* slow down writeout for on-line simulation? */

      if (Pars.waitusec>0);
        usleep(Pars.waitusec);

      /* only repeat what filter OKed */

      if (writeOK)
        {
          siz = write (outData, (char *) ptgd, sizeof (GEBDATA));
//          assert (siz == sizeof (GEBDATA));
          siz = write (outData, (char *) ptinp, ptgd->length);
//          assert (siz == ptgd->length);
          nwrite++;
        };

      /* are we done? */

      curev++;
      if (curev > Pars.nevents)
        goto done;

    };

  /* done */

errordone:
  printf ("GEBFilter: could not read payload\n");
  printf ("could not read more data\n");
done:

  if (Pars.vetocubes)
    {
      printf ("veto_cube:\n");
      d1 = 100.0 * (double) ngood / ((double) ngood + (double) nbad);
      printf ("ngood= %12lli %6.2f%%\n", ngood, d1);
      d1 = 100.0 * (double) nbad / ((double) ngood + (double) nbad);
      printf ("nbad= %12lli %6.2f%%\n", nbad, d1);
    };

  printf ("read  header/payloads: %10i\n", nread);
  printf ("wrote header/payloads: %10i\n", nwrite);

  exit (0);

}
