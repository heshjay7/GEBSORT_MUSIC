
#include <time.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#include <cadef.h>
#include <stdlib.h>

#define STRLEN 128

#define TEST 0

/*-----------------------------------------------------------------*/

float
CAgetval_float (char *ss)
{

  chid mychid;
  char str[128];
  float rval;

  /* find the channel identifier */

  SEVCHK (ca_search (ss, &mychid), "ca_search failure");
  SEVCHK (ca_pend_io (5.0), "ca_pend_io failure");

  /* get the value */

  SEVCHK (ca_get (DBR_FLOAT, mychid, (void *) &rval), "ca_get failure");
  SEVCHK (ca_pend_io (10.0), "ca_pend_io failure");

  return (rval);

};


/*-----------------------------------------------------------------*/

int
CAgetval_int (char *ss)
{

  int ival = -1;
  chid mychid;
  char str[128];

  /* find the channel identifier */

  SEVCHK (ca_search (ss, &mychid), "ca_search failure");
  SEVCHK (ca_pend_io (5.0), "ca_pend_io failure");

  /* get the value */

  SEVCHK (ca_get (DBR_LONG, mychid, (void *) &ival), "ca_get failure");
  SEVCHK (ca_pend_io (10.0), "ca_pend_io failure");

  return (ival);

};

/*-----------------------------------------------------------------*/

int
CAgetval_string (char *ss, char *rstr)
{

  int ival;
  chid mychid;
  char str[128];

  /* find the channel identifier */

  SEVCHK (ca_search (ss, &mychid), "ca_search failure");
  SEVCHK (ca_pend_io (5.0), "ca_pend_io failure");

  /* get the value */

  SEVCHK (ca_get (DBR_STRING, mychid, (void *) rstr), "ca_get failure");
  SEVCHK (ca_pend_io (10.0), "ca_pend_io failure");

  return (0);

};

/*-----------------------------------------------------------------*/

int
send_message (char *message)
{

/* declarations */

  FILE *fp;
  char address[128], command[512];
  int st;

  fp = fopen ("/home/dgs/.LNALARM", "r");
  st = fscanf (fp, "%s", address);
  while (st == 1)
    {
      sprintf (command, "echo \"%s\" | mailx -s LN_ALARM -r gs@anl.gov %s ", message, address);
      if (TEST == 0)
        {
          printf ("*");
          system (command);
          sleep (2);
        }
      else
        {
          printf ("*");
//                          printf ("</br>\n");
//                          printf ("\ntest: %s", command);
        };

      /* next in alarm list */

      st = fscanf (fp, "%s", address);

    };


/* done */

  return (0);

}

/*-----------------------------------------------------------------*/

main (int argc, char **argv)
{
  /* declarations */

  int ival, n = 0, i, stat, nret, i1, nalarm = 0;
  chid mychid;
  char str[STRLEN], message[128], command[512], address[128], *pc, str2[STRLEN];
  float rval, temperature, hi, lo;
  struct tm *local;
  time_t t;
  char hostnam[8], timestamp[32];
  FILE *fp;
  int st;
  char hose[1000][STRLEN];

  /* init */

  if (argc != 1)
    {
      printf ("ooops, shoud have no arguments\n");
      exit (1);
    }

  bzero (&hose[0][0], 1000 * STRLEN);

  /* get time for time stamp */

  t = time (NULL);
  local = localtime (&t);
  sprintf (timestamp, "%s", asctime (local));

  /* get host name */

  gethostname (hostnam, 8);

  /* header */

  printf ("<tt>\n");
  printf ("</br>\n");
  printf ("Gammasphere LN status/alarm system</br>\n");
  printf ("</br>\n");
  printf ("%s: @ %s</br>\n", hostnam, timestamp);

  /* read in the hose list */

  fp = fopen ("det.list", "r");
  if(fp==NULL)
    {
    printf("failed to open \"det.list\"\n");
    exit(1);
    }
  else
    {
  pc = fgets (str, STRLEN, fp);
  while (pc != NULL)
    {
      nret = sscanf (str, "%i %s", &i1, str2);
      strcpy (&hose[i1][0], str2);
      pc = fgets (str, STRLEN, fp);
    };
  fclose (fp);
    };

#if(0)
  for (i = 1; i <= 110; i++)
    printf ("det %3i has hose %s\n", i, &hose[i][0]);
#endif

  /* fire up channel access */

  SEVCHK (ca_task_initialize (), "ca_task_initialize");

  /* first check crates */

  printf ("</br>\n");

  for (i = 1; i <= 6; i++)
    {
      sprintf (str, "HEARTBEAT%1.1i", i);
      ival = CAgetval_int (str);
      if (ival == -1)
        {
          printf ("VXI crate%1.1i is down\n", i);
          sprintf (message, "vxi crate%1.1i is down", i);
          send_message (message);
          printf ("</br>\n");
        }
      else
        printf ("VXI crate%1.1i is up</br>\n", i);
    };

  printf ("</br>\n");

  ival = CAgetval_int ("detector_display");
  if (ival == 0)
    printf ("LN: detector filling: idle\n</br>");
  else
    printf ("LN: detectors are being filled*\n</br>");

  ival = CAgetval_int ("tank_display");
  if (ival == 0)
    printf ("LN: tank filling....: idle\n</br>");
  else
    printf ("LN: tanks are being filled*\n</br>");

  printf ("</br>\n");

  ival = CAgetval_string ("LN_ATLF:XC", str2);
  printf ("last fill: %s </br>\n", str2);
  ival = CAgetval_string ("LN_ATNF:XC", str2);
  printf ("next fill: %s </br>\n", str2);
  ival = CAgetval_string ("LN_TSLF:XC", str2);
  printf ("time since last fill: %s </br>\n", str2);
  ival = CAgetval_string ("LN_TSFS:XC", str2);
  printf ("last fill took: %s </br>\n", str2);

  ival = CAgetval_string ("LN_ALM:XC", str2);
  printf ("number of alarms: %s </br>\n", str2);

  ival = CAgetval_string ("LN_ATLTF:XC", str2);
  printf ("last tank fill: %s </br>\n", str2);

  ival = CAgetval_string ("LN_TSTFS:XC", str2);
  printf ("last tank fill took: %s </br>\n", str2);

  printf ("</br>\n");
  printf ("--list of monitored detectors:\n");
  printf ("</br>\n");

  /* loop through detectors */

  for (i = 1; i <= 110; i++)
    {

      sprintf (str, "MOD%3.3i_DV_EN", i);
      stat = CAgetval_int (str);

      /* process the detectors that are on */

      if (stat == 1)
        {
          n++;
          printf ("ge %3.3i ", i);
          if (stat == 1)
            printf ("_on ");
          else
            printf ("off ");

          sprintf (str, "MOD%3.3i_DV_GEHV", i);
          ival = CAgetval_int (str);
          printf (" %4.4iV ", ival);

          sprintf (str, "MOD%3.3i_DV_TEMP", i);
          temperature = CAgetval_float (str);

          if (i == TEST)
            temperature = 120;

          sprintf (str, "MOD%3.3i_DV_TEMP.HIGH", i);
          hi = CAgetval_float (str);
          sprintf (str, "MOD%3.3i_DV_TEMP.LOW", i);
          lo = CAgetval_float (str);

          printf ("T[K]: [%5.1f - %6.1f]  %5.1f ", lo, hi, temperature);
          printf (" {%5.1f} ", temperature - hi);

          if (temperature >= hi)
            {
              nalarm++;
              printf (" WARM, hose: %s ", &hose[i][0]);

              /* send Email alarm */

              sprintf (message, "GS det %i is warm: %5.1f {%5.1f} hose: %s", i, temperature, temperature - hi,
                       &hose[i][0]);
              send_message (message);

            }
          else if ( (hi-temperature) < 3.0)
            printf (" <--- watch, hose: %s ", &hose[i][0]);

          printf ("\n");
          printf ("</br>\n");
        };
    }

  /* shutdown CA */

  SEVCHK (ca_task_exit (), "ca_task_exit");

  /* done */

  printf ("</br>\n");
  printf ("a total of %i detectors are monitored\n", n);
  printf ("</br>\n");
  printf ("</br>\n");
  printf ("</tt>\n");
  exit (nalarm);

};
