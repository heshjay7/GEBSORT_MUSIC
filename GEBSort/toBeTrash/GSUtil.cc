
/* root utilities for GammaSphere - tl */

#include <stdio.h>
#include <math.h>

#include <TFile.h>
#include <TMapFile.h>
#include <TH1.h>
#include <TF1.h>
#include <TH2.h>
#include <TH3.h>
#include <TCanvas.h>
#include <TSystem.h>
#include <TCutG.h>
#include <TKey.h>
#include <TPolyLine3D.h>
#include <TLegend.h>





#if(0)

/* not needed we think!? */

#define MAXFLOAT 3.40282347e+38
#define MAXINT   2147483647
#endif


/* definitions */

#define MAXGE 400
#define STRLEN 132
#define MINCOUNTS 100
#define PMODE 0644
#define MAXSPELEN 20000
#define FAIL -1
#define UNDETERMINED 0
#define SOLARIS 1
#define LINUX 2
#define CMDTRIES 10
#define MAXDETPOS 120
#define NUMSEGS 36
#define MAXSEGINDEX (MAXDETPOS+1)*(NUMSEGS+1)
#define MAXCRYSTALNO 3

#define BI207LO  569.702
#define BI207HI 1063.662
#define Y88LO    898.0450
#define Y88HI   1836.0630
#define CO60LO  1173.228
#define CO60HI  1332.490

/* global pointers etc */

TMapFile *mfile = 0;            /* shared memory file pointer */
TFile *dfile = 0;               /* disk file pointer */
TCutG *mycutg;                  /* graphical 2dwin window  */
int SystemType = UNDETERMINED;

TCanvas *c1 = 0;
Int_t LineColor = kRed;
int HaveyLimits = 0;
Double_t Ylowlim = 0, Yhilim = 0;
int HavexLimits = 0;
Double_t Xlowlim = 0, Xhilim = 0;

/* prototypes */

int wr_spe (char *, int *, float *);
int rd_spe (char *, int *, float *);
TH1D *FindSpectrum (char *);
void ealign (Char_t *, Char_t *, Char_t *, Double_t, Double_t, Int_t, Char_t *, Int_t);
void talign (Char_t *, Char_t *, Char_t *, float, int, int, int);
void dod1 (Char_t *, Char_t *, Float_t, Float_t, Int_t);
void d1 (char *, Float_t xmin = -999999., Float_t xmax = 999999.);
void d1o (char *, Float_t xmin = -999999., Float_t xmax = 999999.);
void str_decomp (Char_t *, Int_t, Int_t *);
TCutG *rd2dwin (Char_t *);
int read_dm (char *, int *, int *, char *, char *);
int wr_dm (char *, int, int, char *, char *);
int rd_mat (char *, int, int, float *, char *);
int wr_mat (char *, int, int, float *, char *);
TH2D *FindMatrix (Char_t *);
int wr_Radmat (char *, unsigned short int *);
int mkGSSortCmdFile (Char_t *);

/*--------------------------------------*/

int
sm3 (int dim, float *y)
{
  /* three point smooth of an array, from Bevington */

  /* declarations */

  int i;
  float z[30000];

  /* handle start point */

  z[0] = (2. * y[0] + y[1] + y[2]) / 4.;

  /* handle intermediate points */

  for (i = 1; i < dim - 1; i++)
    z[i] = (y[i - 1] + 2. * y[i] + y[i + 1]) / 4.;

  /* handle end point */

  z[dim - 1] = (2. * y[dim - 1] + y[dim - 2] + y[dim - 3]) / 4.;

  /* transfer smoothed array back to input array */

  for (i = 0; i < dim; i++)
    y[i] = z[i];

  /* done */

  return (0);

}

/*-------------------------------------------------------------------*/

void
setylimits (Double_t d1, Double_t d2)
{
  /* specify ylimits */

  HaveyLimits = 1;
  Ylowlim = d1;
  Yhilim = d2;

}

/*-------------------------------------------------------------------*/

void
unsetylimits ()
{
  HaveyLimits = 0;
}

/*-------------------------------------------------------------------*/

void
setxlimits (Double_t d1, Double_t d2)
{
  /* specify ylimits */

  HavexLimits = 1;
  Xlowlim = d1;
  Xhilim = d2;

}

/*-------------------------------------------------------------------*/

void
unsetxlimits ()
{
  HavexLimits = 0;
}

/*-------------------------------------------------------------------*/

void
showlimits ()
{
  if (HavexLimits)
    printf ("x limits from %f to %f\n", Xlowlim, Xhilim);
  else
    printf ("auto x limits\n");
  if (HaveyLimits)
    printf ("y limits from %f to %f\n", Ylowlim, Yhilim);
  else
    printf ("auto y limits\n");
}

/*-------------------------------------------------------------------*/

int
GetSystemType (void)
{
  /* determine what system we are on   */
  /* is there a better way to do this? */


  /* declarations */

  FILE *fp;
  char *p;
  char *pc, str[STRLEN];
  int type;

  /* create a file with the info */

  gSystem->Exec ("uname > /tmp/tmpGSSortID");

  /* open/read/delete */

  if ((fp = fopen ("/tmp/tmpGSSortID", "r")) == NULL)
    {
      printf ("GetSystemType/error: could not open /tmp/tmpGSSortID\n");
      return (FAIL);
    };
  pc = fgets (str, STRLEN, fp);
  fclose (fp);
  gSystem->Exec ("rm /tmp/tmpGSSortID");

  /* evaluate */

  if ((p = strstr (str, "SunOS")) != NULL)
    type = SOLARIS;
  else if ((p = strstr (str, "Linux")) != NULL)
    type = LINUX;
  else
    type = UNDETERMINED;

  /* done */

  return (type);

}


/*-------------------------------------------------------------------*/

void
rdmat (Char_t * fn)
  /** Read a two-dimensional matrix from a file in to Root. 
      Example: rdmat("matrix_example")
        will read the following files:
          matrix_example.dm  (contains dimensions, machinetype & datatype)
          matrix_example     (contains data)
          
  **/
{
  /* declarations */

  int nx, ny, k, i, j;
  char machinetype[32], datatype[32];
  float *FloatMat1;
  TH2F *hist2 = NULL;

  if (1)
    {
      printf ("please use rdmat_ascii instead\n");
      return;
    };

  /* read info file */

  /* read_dm is in 2d_fun.cc */
  if (read_dm (fn, &nx, &ny, machinetype, datatype) == -1)
    {
      printf ("rdmat: aborted because read_dm failed.\n");
      fflush (stdout);
      return;
    }

  k = nx * ny * sizeof (float);
  FloatMat1 = (float *) malloc (k);
  k = nx * ny;
  for (i = 0; i < k; i++)
    *(FloatMat1 + i) = 0;

  /* read data */

  rd_mat (fn, nx, ny, FloatMat1, datatype);

  /* create root 2D histogram */

  printf ("creating %s ...", fn);
  fflush (stdout);
  hist2 = new TH2F (fn, fn, nx, 0, nx, ny, 0, ny);
  hist2->GetXaxis ()->SetLimits (0, nx);
  hist2->GetYaxis ()->SetLimits (0, ny);
  printf ("done\n");
  fflush (stdout);

  /* xfer data */

  {
    TAxis *x = hist2->GetXaxis ();
    TAxis *y = hist2->GetYaxis ();

    for (i = 0; i < nx; i++)
      for (j = 0; j < ny; j++)
        hist2->SetBinContent (x->FindBin (i), y->FindBin (j), *(FloatMat1 + nx * i + j));
  }

  /* done */

  free (FloatMat1);

};


/*-------------------------------------------------------------------*/
int
rdmat_ascii (Char_t * name)
{
  /* read a matrix in simple ascii format */

  TH2F *hist2;
  FILE *fp;
  Int_t nx, ny;                 /* width and height of matrix */
  Int_t i, j;                   /* iterators */
  float val;
  int st;
  char str[512], machinetype[100], datatype[100];
  double dnn, lo, hi, numEntries;

  /* read the dm file */

  sprintf (str, "%s.dm", name);
  fp = fopen (str, "r");
  if (fp == NULL)
    printf ("failed to open \"%s\"\n", str);
  else
    printf ("openened \"%s\"\n", str);
  fflush(stdout);

  st = fscanf (fp, "\n%i", &nx);
  st = fscanf (fp, "\n%i", &ny);
  st = fscanf (fp, "\n%s", machinetype);
  st = fscanf (fp, "\n%s", datatype);
  st = fclose (fp);

  printf ("x dim.......: %4i\n", nx);
  printf ("y dim.......: %4i\n", ny);
  printf ("machine type: %s\n", machinetype);
  printf ("data type...: %s\n", datatype);


  /* delete any old version of matrix */

  hist2 = (TH2F *) gROOT->FindObject (name);

  if (hist2 != NULL)
    {
      printf ("\"%s\" exist, deleting it\n", name);
      hist2->Delete ();
      fflush (stdout);
    }

  /* create root 2D histogram */

  printf ("creating %s ...", name);
  fflush (stdout);

  if (nx != ny)
    {
      printf ("can only handle square matrices\n");
      return (0);
    };

  /* create the same way as we made the HK matrices */

  dnn = nx;
  lo = -0.5;
  hi = dnn + lo;

  sprintf (str, "%s [rdmat_ascii input]", name);
  hist2 = new TH2F (name, str, dnn, lo, hi, dnn, lo, hi);
  printf ("done\n");
  fflush (stdout);

  /* open the ascii data file */

  fp = fopen (name, "r");
  if (fp == NULL)
    printf ("failed to open \"%s\"\n", name);
  else
    printf ("openened \"%s\"\n", name);

  /* read and fill the name matrix */

  numEntries = 0;
  st = fscanf (fp, "\n%i %i %f", &i, &j, &val);
  while (st > 0)
    {

      printf ("%3i %3i %f\n", i, j, val);
      numEntries += val;

      /* fill matrix */

      hist2->SetBinContent ((double) i, (double) j, (double) val);

      /* read next */

      st = fscanf (fp, "\n%i %i %f", &i, &j, &val);

    }
  hist2->SetEntries (numEntries);


  fclose (fp);
  return (0);

};

/*-------------------------------------------------------------------*/
int
wrmat_ascii (Char_t * name)
{

  /* write a matrix in simple ascii format */

  TH2F *hist2;
  Int_t nx, ny;                 /* width and height of matrix */
  Int_t i, j;                   /* iterators */
  float val;
  FILE *fp;
  char str[512];

  hist2 = (TH2F *) gROOT->FindObject (name);

  if (hist2 == NULL)
    {
      printf ("wrmat: Matrix %s does not exist.\n", name);
      fflush (stdout);
      return (-2);
    }

  nx = hist2->GetNbinsX ();
  ny = hist2->GetNbinsY ();
  printf ("GetNbinsX ()= %i,GetNbinsY () = %i\n", nx, ny);

  fp = fopen (name, "w");
  if (fp == NULL)
    printf ("failed to open \"%s\"\n", name);
  else
    printf ("openened \"%s\"\n", name);

  {
    TAxis *x = hist2->GetXaxis ();
    TAxis *y = hist2->GetYaxis ();
    for (i = 0; i < nx; i++)
      for (j = 0; j < ny; j++)
        {
          val = hist2->GetBinContent (x->FindBin (i), y->FindBin (j));
          if (val > 0 && i < 127 && j < 127)
            {
              printf ("%3i ", x->FindBin (i));
              fprintf (fp, "%3i ", x->FindBin (i));
              printf ("[%6.2f-%6.2f] ", x->GetBinLowEdge (i), x->GetBinLowEdge (i) + x->GetBinWidth (i));
              printf ("%3i ", y->FindBin (j));
              fprintf (fp, "%3i ", y->FindBin (j));
              printf ("[%6.2f-%6.2f] ", y->GetBinLowEdge (j), y->GetBinLowEdge (j) + y->GetBinWidth (j));
              printf ("%f ", hist2->GetBinContent (x->FindBin (i), y->FindBin (j)));
              fprintf (fp, "%f ", hist2->GetBinContent (x->FindBin (i), y->FindBin (j)));
              printf ("\n");
              fprintf (fp, "\n");
            };
        };

    fclose (fp);
    printf ("closed \"%s\"\n", name);

    sprintf (str, "%s.dm", name);
    fp = fopen (str, "w");
    if (fp == NULL)
      printf ("failed to open \"%s\"\n", str);
    else
      printf ("openened \"%s\"\n", str);

    fprintf (fp, "%i\n", nx);
    fprintf (fp, "%i\n", ny);
    fprintf (fp, "linux\n");
    fprintf (fp, "ascii\n");

    fclose (fp);
    printf ("closed \"%s\"\n", str);

    printf ("x BinWidth = %f\n", x->GetBinWidth (0));
    printf ("y BinWidth = %f\n", y->GetBinWidth (0));
  }

  return (0);

};

/*-------------------------------------------------------------------*/

int
wrmat (Char_t * name)
{
  /* declarations */

  TH2F *hist2;
  float *FloatMat1;
  Int_t nx, ny;                 /* width and height of matrix */
  Int_t i, j;                   /* iterators */
  Int_t need_bytes;

  if (1)
    {
      printf ("please use wrmat_ascii instead\n");
      return (0);
    };

  hist2 = (TH2F *) gROOT->FindObject (name);

  if (hist2 == NULL)
    {
      printf ("wrmat: Matrix %s does not exist.\n", name);
      fflush (stdout);
      return (-2);
    }

  nx = hist2->GetNbinsX ();
  ny = hist2->GetNbinsY ();
  printf ("GetNbinsX ()= %i,GetNbinsY () = %i\n", nx, ny);


  /** allocate memory **/

  need_bytes = (nx) * (ny) * sizeof (float);
//  FloatMat1 = (float *) malloc (need_bytes);
  FloatMat1 = (float *) calloc (1, need_bytes);

  /** Extract the matrix in to a 2-dimensional array
      (inverse of rdmat) **/

  /* This part is very slow; for a 2048 by 2048 matrix,
     it has to call GetBinContent 4.2 million times!   */

  {
    TAxis *x = hist2->GetXaxis ();
    TAxis *y = hist2->GetYaxis ();
    for (i = 0; i < nx; i++)
      for (j = 0; j < ny; j++)
        {
          *(FloatMat1 + ny * i + j) = hist2->GetBinContent (x->FindBin (i), y->FindBin (j));
          if (*(FloatMat1 + ny * i + j) > 0 && i < 127 && j < 127)
            {
              printf ("%i ", x->FindBin (i));
              printf ("[%6.2f-%6.2f] ", x->GetBinLowEdge (i), x->GetBinLowEdge (i) + x->GetBinWidth (i));
              printf ("%i ", y->FindBin (j));
              printf ("[%6.2f-%6.2f] ", y->GetBinLowEdge (j), y->GetBinLowEdge (j) + y->GetBinWidth (j));
              printf ("%f ", hist2->GetBinContent (x->FindBin (i), y->FindBin (j)));
              printf ("\n");
            };

        };
    printf ("x BinWidth = %f\n", x->GetBinWidth (0));
    printf ("y BinWidth = %f\n", y->GetBinWidth (0));
  }

  {
    char datatype[32];
    char machinetype[32];
    sprintf (datatype, "ascii");
    sprintf (machinetype, "unknown");

    /* wr_mat and wr_dm are defined in 2d_fun.c */
    printf ("wrmat: Writing file %s in %s format...", name, datatype);
    fflush (stdout);

    wr_mat (name, nx, ny, FloatMat1, datatype);

    printf ("done.\nwrmat: Writing file %s.dm...", name);

    wr_dm (name, nx, ny, machinetype, datatype);

//    free(FloatMat1);


    printf ("done.\n");
    fflush (stdout);
  }

  free (FloatMat1);

  return (0);                   /* Success */
}


/*-------------------------------------------------------------------*/

void
confirm ()
{
  /* attempt to allow user to pause in script files */
  /* unfortunately it doesn't work for mk2dwin      */
  /* which was originally why I invented confirm()  */
  /* but it might still be useful                   */

  char ans;

  printf ("press return to continue ");
  ans = getchar ();
  printf ("\n");

}

/*-------------------------------------------------------------------*/

void
mkcanvas ()
{
  /* make a default simple canvas with */
  /* some useful properties........... */

  c1 = (TCanvas *) gROOT->FindObject ("c1");
  if (c1 != NULL)
    c1->Delete ();
  c1 = new TCanvas ("c1");
  c1->SetCrosshair (1);
  c1->SetFillColor (0);
  c1->ToggleEventStatus ();
  c1->SetTitle ("c1 GT canvas");
  printf ("__pointer to canvas is c1\n");
}

/*-------------------------------------------------------------------*/

void
mkcanvas_square ()
{
  /* make a default simple canvas with */
  /* some useful properties........... */

  c1 = (TCanvas *) gROOT->FindObject ("c1");
  if (c1 != NULL)
    c1->Delete ();
  c1 = new TCanvas ("c1", "c1", 700, 700);
  c1->SetCrosshair (1);
  c1->SetFillColor (0);
  c1->ToggleEventStatus ();
  c1->SetTitle ("c1 GT canvas");
  printf ("__pointer to canvas is c1\n");
}


/*-------------------------------------------------------------------*/

void
version ()
{
  printf ("GammaSphere root display utility\n");
  printf ("extract from svn\n");
  printf ("this version has Float 2D matrices\n");
}

/*-------------------------------------------------------------------*/

void
md1 (Char_t * List, Char_t * GenNam, Float_t xmin, Float_t xmax, Int_t wait = 1)
{
  /* declarations */

  Char_t str[STRLEN];
  Int_t i, firsttime = 1;
  Int_t DisplayIt[MAXGE + 1];

  if (c1 != NULL)
    c1->SetFillColor (0);       /* white background */
  else
    printf ("no canvas\n");

  /* find list of detectors to display */

  str_decomp (List, MAXGE + 1, DisplayIt);

  /* display them */

  firsttime = 1;
  for (i = 1; i <= MAXGE; i++)
    if (DisplayIt[i])
      {
        if (firsttime || wait == 2 || wait == 3)
          {
            firsttime = 0;
            sprintf (str, "%s%3.3i", GenNam, i);
            LineColor = 1;
            d1 (str, xmin, xmax);
            if (c1 != NULL && wait)
              c1->Update ();
            if (wait == 3)
              {
                printf ("hit return for next spectrum...");
                getc (stdin);
              };
          }
        else
          {
            sprintf (str, "%s%3.3i", GenNam, i);
            if (wait)
              {
                printf ("hit return for next spectrum...");
                getc (stdin);
              };
            d1o (str, xmin, xmax);
            if (c1 != NULL && wait)
              c1->Update ();
          };
      };
  c1->Update ();

  /* reset line color */

  LineColor = kRed;
  if (c1 != NULL)
    c1->SetFillColor (0);

};

/*-------------------------------------------------------------------*/

void
d1 (char *histname, Float_t xmin, Float_t xmax)
{
  /* declarations */

  Char_t str[STRLEN];

  /* simple display, clear old spectrum */

//  sprintf (str, "");
  str[0] = 0;
  LineColor = 2;
  dod1 (str, histname, xmin, xmax, 0);
}

/*-------------------------------------------------------------------*/

void
d1o (char *histname, Float_t xmin, Float_t xmax)
{
  /* declarations */

  Char_t str[STRLEN];

  /* cycle line colors */

  LineColor++;                  /* new color */
  if (LineColor == 10)
    LineColor = 1;

  /* simple display, clear old spectrum */

  sprintf (str, "SAME");
  dod1 (str, histname, xmin, xmax, 1);
}

/*-------------------------------------------------------------------*/

void
tgealign (Char_t * oldf, Char_t * newf, float dpos = 4000, int di = 4, int lo = 10, int bp = 0)
{
  /* declarations */

  Char_t str[STRLEN];

  /* calibrate germanium time signal */

  sprintf (str, "tge");
  talign (str, oldf, newf, dpos, di, lo, bp);

}

/*-------------------------------------------------------------------*/

void
tbgoalign (Char_t * oldf, Char_t * newf, float dpos = 4000, int di = 8, int lo = 10, int bp = 5)
{
  /* declarations */

  Char_t str[STRLEN];

  /* calibrate germanium time signal */

  sprintf (str, "tbgo");
  talign (str, oldf, newf, dpos, di, lo, bp);

}

/*-------------------------------------------------------------------*/

void
ehialign (Char_t * oldf, Char_t * newf, Double_t kevch, Double_t sgain, Int_t w1, Char_t * source, Int_t minchan)
{
  /* declarations */

  Char_t str[STRLEN];

  /* calibrate hi res germanium signal */

  sprintf (str, "EhiCln");
  ealign (str, oldf, newf, kevch, sgain, w1, source, minchan);

}

/*-------------------------------------------------------------------*/

void
eloalign (Char_t * oldf, Char_t * newf, Double_t kevch, Int_t w1, Char_t * source, Int_t minchan)
{
  /* declarations */

  Char_t str[STRLEN];
  Double_t sgain = 1.0;

  /* calibrate lo res germanium signal */

  sprintf (str, "elo");
  ealign (str, oldf, newf, kevch, sgain, w1, source, minchan);

}

/*-------------------------------------------------------------------*/

void
esidealign (Char_t * oldf, Char_t * newf, Double_t kevch, Int_t w1, Char_t * source, Int_t minchan)
{
  /* declarations */

  Char_t str[STRLEN];
  Double_t sgain = 1.0;

  /* calibrate side channel germanium signal */

  sprintf (str, "eside");
  ealign (str, oldf, newf, kevch, sgain, w1, source, minchan);

}

/*-----------------------------------------------------------------*/

void
add1_old (Char_t * h1, Char_t * h2, Double_t af = 1.0)
{

  /* declarations */

  Int_t ok1, ok2, lo1, hi1, nn1, lo2, hi2, nn2;
  Axis_t xlo, xhi;
  TH1D *hist1 = NULL, *hist2 = NULL;

  /* get parameters for h2 */

  printf ("checking for spectrum %s:\n", h2);
  fflush (stdout);
  ok2 = 1;
  if (mfile != 0)
    {
      hist2 = (TH1D *) mfile->Get (h2);
      if (hist2 == NULL)
        {
          printf ("spectrum not in shared memory\n");
          hist2 = (TH1D *) gROOT->FindObject (h2);
          if (hist2 == NULL)
            {
              printf ("spectrum not local either...\n");
              ok2 = 0;
            }
          else
            printf ("spectrum found locally\n");

        }
      else if (dfile != 0)
        hist2 = (TH1D *) gROOT->FindObject (h2);
    }
  else
    {
      printf ("no shared mem or disk file loaded!\n");
      ok2 = 0;
    }

  if (ok2)
    {
      printf ("found it ");
      hist2->Print ();
      lo2 = hist2->GetXaxis ()->GetFirst ();
      hi2 = hist2->GetXaxis ()->GetLast ();
      nn2 = hist2->GetXaxis ()->GetNbins ();
      xlo = hist2->GetXaxis ()->GetXmin ();
      xhi = hist2->GetXaxis ()->GetXmax ();
      printf ("__%i bins, from %i to %i\n", nn2, lo2, hi2);
      printf ("__x limits from %i to %i\n", (int) xlo, (int) xhi);
      fflush (stdout);

      /* see if h1 exists */

      printf ("checking for spectrum %s:\n", h1);
      fflush (stdout);
      ok1 = 1;
      if (mfile != 0)
        {
          hist1 = (TH1D *) mfile->Get (h1);
          if (hist1 == NULL)
            {
              printf ("spectrum not in shared memory\n");
              hist1 = (TH1D *) gROOT->FindObject (h1);
              if (hist1 == NULL)
                {
                  printf ("spectrum not local either...\n");
                  ok1 = 0;
                }
              else
                printf ("spectrum found locally\n");

            }
          else if (dfile != 0)
            hist1 = (TH1D *) gROOT->FindObject (h1);
        }
      else
        {
          printf ("no shared mem or disk file loaded!\n");
          ok1 = 0;
        };
      fflush (stdout);

      if (ok1 == 0)
        {
          /* create h1 */

          printf ("creating %s ...", h1);
          fflush (stdout);
          hist1 = new TH1D (h1, "sum", nn2, lo2, hi2);
          hist1->GetXaxis ()->SetLimits (xlo, xhi);
          printf ("done\n");
          fflush (stdout);

        }
      else
        printf ("%s already exists:\n", h1);

      /* print info on h1 */

      hist1->Print ();
      lo1 = hist1->GetXaxis ()->GetFirst ();
      hi1 = hist1->GetXaxis ()->GetLast ();
      nn1 = hist1->GetXaxis ()->GetNbins ();
      printf ("__%i bins, from %i to %i\n", nn1, lo1, hi1);
      fflush (stdout);

      /* add */

      hist1->Add (hist2, af);
      printf ("added %s to %s with factor: %f\n", h2, h1, (float) af);
      fflush (stdout);

    };

};

/*-----------------------------------------------------------------------*/

void
add1 (Char_t * h1, Char_t * h2, Double_t af = 1.0)
{

  /* declarations */

  Int_t lo1, hi1, nn1, lo2, hi2, nn2;
  Axis_t xlo, xhi;
  TH1D *hist1 = NULL, *hist2 = NULL;

  /* find spectra */

  hist1 = FindSpectrum (h1);
  if (hist1 != NULL)
    printf ("spectrum %s found\n", h1);
  else
    printf ("spectrum %s not found\n", h1);
  hist2 = FindSpectrum (h2);
  if (hist2 != NULL)
    printf ("spectrum %s found\n", h2);
  else
    printf ("spectrum %s not found\n", h2);


  if (hist2 != NULL)
    {
      printf ("found %s ", h2);
      hist2->Print ();
      lo2 = hist2->GetXaxis ()->GetFirst ();
      hi2 = hist2->GetXaxis ()->GetLast ();
      nn2 = hist2->GetXaxis ()->GetNbins ();
      xlo = hist2->GetXaxis ()->GetXmin ();
      xhi = hist2->GetXaxis ()->GetXmax ();
      printf ("__%i bins, from %i to %i\n", nn2, lo2, hi2);
      printf ("__x limits from %i to %i\n", (int) xlo, (int) xhi);
      fflush (stdout);


      if (hist1 == NULL)
        {
          /* create h1 */

          printf ("creating %s ...", h1);
          fflush (stdout);
          hist1 = new TH1D (h1, "sum", nn2, lo2, hi2);
          hist1->GetXaxis ()->SetLimits (xlo, xhi);
          printf ("done\n");
          fflush (stdout);

        }
      else
        printf ("%s already exists:\n", h1);

      /* print info on h1 */

      hist1->Print ();
      lo1 = hist1->GetXaxis ()->GetFirst ();
      hi1 = hist1->GetXaxis ()->GetLast ();
      nn1 = hist1->GetXaxis ()->GetNbins ();
      printf ("__%i bins, from %i to %i\n", nn1, lo1, hi1);
      fflush (stdout);

      /* add */

      hist1->Add (hist2, af);
      printf ("added %s to %s with factor: %f\n", h2, h1, (float) af);
      fflush (stdout);

    };

};


/*-----------------------------------------------------------------*/

void
del (Char_t * h1)
{

  /* declarations */

  TH1D *hist1;

  /* check where spectrum is */

  hist1 = FindSpectrum (h1);

  /* delete */

  if (hist1 != NULL)
    {
      hist1->Print ();
      hist1->Delete ();
      printf ("%s deleted!\n", h1);
    }
  else
    printf ("%s not found\n", h1);

};

/*-----------------------------------------------------------------*/
/* experimental */
/* doesn 't really work for shared, why? */

void
find ()
{

  TList *zlist;
  TIterator *hiterator;
  TH1 *htemp;


#if(0)
  /* shared mem spectra */

  zlist = mfile->GetList ();
  hiterator = zlist->MakeIterator ();
  while (htemp = (TH1 *) hiterator->Next ())
    htemp->Print ();
#endif


  /* local spectra */

  zlist = gDirectory->GetList ();
  hiterator = zlist->MakeIterator ();
  while ((htemp = (TH1 *) hiterator->Next ()))
    htemp->Print ();

};

/*-----------------------------------------------------------------*/

void
update ()
{

  /* signal the GSSort sorter to update shared memory */

  /* declarations */

  char str[STRLEN];

  /* execute command */

  sprintf (str, "pkill -HUP GEBSort");
  printf ("%s\n", str);
  fflush (stdout);
  gSystem->Exec (str);
  gSystem->Exec ("date");

  /* done */

};

/*-----------------------------------------------------------------*/

void
zapcounters ()
{


  /* declarations */

  Char_t str[STRLEN];

  /* signal the GSSort sorter to zap its counters */

  printf ("will stop GSSort gracefully\n");

  /* write command file */

  sprintf (str, "zapcounters");
  mkGSSortCmdFile (str);

  /* execute command */

  sprintf (str, "pkill -HUP GSSort; pkill -HUP GTSort");
  printf ("%s\n", str);
  fflush (stdout);
  gSystem->Exec (str);
  gSystem->Exec ("date");

  /* done */

}

/*-----------------------------------------------------------------*/

void
stopsort ()
{


  /* declarations */

  Char_t str[STRLEN];

  /* signal the GSSort sorter to stop gracefully */

  printf ("will stop GSSort gracefully\n");

  /* write command file */

  sprintf (str, "stopsort");
  mkGSSortCmdFile (str);

  /* execute command */

  sprintf (str, "pkill -HUP GSSort; pkill -HUP GTSort");
  printf ("%s\n", str);
  fflush (stdout);
  gSystem->Exec (str);
  gSystem->Exec ("date");

  /* done */

}

/*-----------------------------------------------------------------*/

void
killsort ()
{

  /* stop the GSSort sorter cold! */

  /* declarations */

  char str[STRLEN];

  /* execute command */

  sprintf (str, "pkill -KILL GSSort; pkill -KILL GTSort;");
  printf ("%s\n", str);
  fflush (stdout);
  gSystem->Exec (str);
  gSystem->Exec ("date");

  /* done */

}

/*-----------------------------------------------------------------*/

void
submitsort (Char_t * SortEngine, Char_t * chatfile, Char_t * logfile)
{
  /* declarations */

  char str[STRLEN];

  /* delete old logfile */

  sprintf (str, "\\rm %s", logfile);
  printf ("%s\n", str);
  gSystem->Exec (str);

  /* start the sort */

  sprintf (str, "%s -chat %s > %s &", SortEngine, chatfile, logfile);
  printf ("%s\n", str);
  fflush (stdout);
  gSystem->Exec (str);
  if (SystemType == SOLARIS)
    gSystem->Exec ("ps -ef | grep GSSort | grep -v grep");
  else if (SystemType == SOLARIS)
    gSystem->Exec ("ps -axe | grep GSSort | grep -v grep");
  gSystem->Exec ("date");

}

/*-----------------------------------------------------------------*/

void
dttermlog (Char_t * logfile)
{

  /* start sort logging in dtterm window */

  /* declarations */

  char str[STRLEN];

  sprintf (str, "dtterm -sl 2000 -e tail -f %s &", logfile);
  printf ("%s\n", str);
  fflush (stdout);
  gSystem->Exec (str);
  gSystem->Exec ("date");

};

/*-----------------------------------------------------------------*/

void
xtermlog (Char_t * logfile)
{

  /* start sort logging in dtterm window */

  /* declarations */

  char str[STRLEN];

  sprintf (str, "xterm -sl 2000 -T %s -e tail -f %s &", logfile, logfile);
  printf ("%s\n", str);
  fflush (stdout);
  gSystem->Exec (str);
  gSystem->Exec ("date");

};

/*-----------------------------------------------------------------*/

void
startsort (Char_t * SortEngine, Char_t * chatfile, Char_t * logfile)
{
  SystemType = GetSystemType ();
  submitsort (SortEngine, chatfile, logfile);
  if (SystemType == SOLARIS)
    dttermlog (logfile);
  else if (SystemType == LINUX)
    xtermlog (logfile);
  else
    printf ("don't know how to start a logging window\n");

  /* remind user to reload */

  printf ("Remember to (re)load shared memory mapfile or root file\n");
  printf ("sload(\"mapfile\") or dload(\"rootfile\")\n");

};

/*-----------------------------------------------------------------*/

void
te (Char_t * fn)
{
  /* edit a file */

  /* declarations */

  char str[STRLEN];

  /* execute command */

  sprintf (str, "/usr/dt/bin/dtpad -standAlone -statusLine %s &", fn);
  printf ("%s\n", str);
  fflush (stdout);
  gSystem->Exec (str);
  gSystem->Exec ("date");

  /* done */

}

/*-----------------------------------------------------------------*/

void
setdumpinterval (Int_t dt)
{

  /* set dump time in GSSort */

  /* declarations */

  Char_t str[STRLEN];

  printf ("will set GSSort program to update shared\n");
  printf ("memory every %i minutes\n", (int) dt);

  /* write command file */

  sprintf (str, "dumpevery %i", (int) dt);
  mkGSSortCmdFile (str);

  /* execute command */

  sprintf (str, "pkill -HUP GSSort; pkill -HUP GTSort");
  printf ("%s\n", str);
  fflush (stdout);
  gSystem->Exec (str);
  gSystem->Exec ("date");

  /* done */

}

/*-----------------------------------------------------------------*/

void
printevents (int i)
{

  /* tell GSSort to print events */

  /* declarations */

  Char_t str[STRLEN];

  /* write command file */

  sprintf (str, "printevents %i", i);
  mkGSSortCmdFile (str);

  /* execute command */

  sprintf (str, "pkill -HUP GSSort; pkill -HUP GTSort");
  printf ("%s\n", str);
  fflush (stdout);
  gSystem->Exec (str);
  gSystem->Exec ("date");

  /* done */

}

/*-----------------------------------------------------------------*/

void
zapall ()
{

  /* zero all spectra in shared memory....... */
  /* by signaling GSSort through command file */
  /* (we cannot zero spectra in shared memory */
  /* from rootn.exe since the shared memory. */
  /* is open in read mode).................. */

  /* declarations */

  Char_t str[STRLEN];

  /* write command file */

  sprintf (str, "zapall");
  mkGSSortCmdFile (str);

  /* execute command */

  sprintf (str, "pkill -HUP GSSort; pkill -HUP GTSort");
  printf ("%s\n", str);
  fflush (stdout);
  gSystem->Exec (str);
  gSystem->Exec ("date");

  /* done */

}

/*-----------------------------------------------------------------*/

void
zap (Char_t * spname)
{

  /* zero spname in shared memory............ */
  /* by signaling GSSort through command file */
  /* (we cannot zero spectra in shared memory */
  /* from rootn.exe since the shared memory. */
  /* is open in read mode).................. */

  /* declarations */

  Char_t str[STRLEN];

  /* write command file */

  sprintf (str, "zap %s", spname);
  mkGSSortCmdFile (str);

  /* execute command */

  sprintf (str, "pkill -HUP GSSort; pkill -HUP GTSort");
  printf ("%s\n", str);
  fflush (stdout);
  gSystem->Exec (str);
  gSystem->Exec ("date");

  /* done */

}

/*-----------------------------------------------------------------*/

void
setbeta (float beta)
{
  /* declarations */

  char str[STRLEN];

  /* write command file */

  sprintf (str, "beta %f", beta);
  mkGSSortCmdFile (str);

  /* execute command */

  sprintf (str, "pkill -HUP GSSort; pkill -HUP GTSort");
  printf ("%s\n", str);
  fflush (stdout);
  gSystem->Exec (str);
  gSystem->Exec ("date");

  /* done */

}

/*-----------------------------------------------------------------*/

void
sett0 (float t0)
{
  /* declarations */

  char str[STRLEN];

  /* write command file */

  sprintf (str, "evtimet0 %f", t0);
  mkGSSortCmdFile (str);

  /* execute command */

  sprintf (str, "pkill -HUP GSSort; pkill -HUP GTSort");
  printf ("%s\n", str);
  fflush (stdout);
  gSystem->Exec (str);
  gSystem->Exec ("date");

  /* done */

}

/*-----------------------------------------------------------------*/

void
resett0 (void)
{
  /* declarations */

  char str[STRLEN];

  /* write command file */

  sprintf (str, "resett0");
  mkGSSortCmdFile (str);

  /* execute command */

  sprintf (str, "pkill -HUP GSSort; pkill -HUP GTSort");
  printf ("%s\n", str);
  fflush (stdout);
  gSystem->Exec (str);
  gSystem->Exec ("date");

  /* done */

}

/*-----------------------------------------------------------------*/

void
status (void)
{
  /* declarations */

  char str[STRLEN];

  /* write command file */

  sprintf (str, "status");
  mkGSSortCmdFile (str);

  /* execute command */

  sprintf (str, "pkill -HUP GSSort; pkill -HUP GTSort");
  printf ("%s\n", str);
  fflush (stdout);
  gSystem->Exec (str);
  gSystem->Exec ("date");

  /* done */

}

/*-----------------------------------------------------------------*/

void
segmentHitPattern (void)
{
  /* declarations */

  char str[STRLEN];

  /* write command file */

  sprintf (str, "segmentHitPattern");
  mkGSSortCmdFile (str);

  /* execute command */

  sprintf (str, "pkill -HUP GSSort; pkill -HUP GTSort");
  printf ("%s\n", str);
  fflush (stdout);
  gSystem->Exec (str);
  gSystem->Exec ("date");

  /* done */

}

/*-----------------------------------------------------------------*/


int
rdspe (Char_t * dfn, Char_t * rfn)
{

  /* declarations */

  int dim = 20000, st, i;
  float sp[20000];
  TH1D *tmppt;
  char str1[132];

  /* read the spectrum */

  st = rd_spe (dfn, &dim, sp);
  printf ("read spectrum \"%s\" of length %i\n", dfn, dim);
  fflush (stdout);

  /* make spectrum */

  sprintf (str1, "%s", rfn);
  tmppt = new TH1D (str1, str1, dim, 1, dim);
  printf ("Created Object \"%s\", %p\n", str1, tmppt);

  /* transfer */

  for (i = 1; i < dim; i++)
    tmppt->Fill (i, (double) sp[i]);

  /* done */

  return (0);

}


/*-----------------------------------------------------------------*/

void
wrspe (Char_t * histname, Char_t * fn)
{

  /* write a spectrum out in Radford spe format */

  TH1D *hist1;
  Int_t i, nn;
  Int_t i1, i2;
  float sp[MAXSPELEN];
  double_t d1,d2;

  /* get spectrum */

  hist1 = FindSpectrum (histname);

  if (hist1 != NULL)
    {

      /* get x dimension information */

      i1 = (int) hist1->GetXaxis ()->GetXmin ();
      printf ("x from %4i ", i1);
      i2 = (int) hist1->GetXaxis ()->GetXmax ();
      printf ("to %4i; ", i2);
      fflush (stdout);
      nn = hist1->GetXaxis ()->GetNbins ();
      printf (" %i channels\n", nn);

      d1=hist1->GetXaxis ()->GetBinCenter(1);
      d2=hist1->GetXaxis ()->GetBinCenter(nn);
      printf("%f to %f (for center) in %i channels\n", d1, d2,nn);


//      if (ul1> i1) i1=ul1;
//      if (ul2< i2) i2=ul2;
//      nn=i2-i1;

      /* fill spe float spectrum */

//      for (i = 0; i < i1; i++)
//        sp[i] = 0;
      for (i = 1; i <= nn; i++)
        sp[i] = (float) hist1->GetBinContent (i);

//      sp[nn + i1] = 0;          /* just to be safe */

      /* write spe file */

      if (nn >= 16384)
        nn = 16383;             /* so gf3 can handle spectrum */
      wr_spe (fn, &nn, sp);
      printf ("wrote %i channels to \"%s\"\n", nn, fn);

    }
  else
    printf ("spectrum %s not found\n", histname);

  /* done */

}

/*-----------------------------------------------------------------*/

void
pjx (Char_t * histname, Char_t * pname, Float_t ymin = -999999., Float_t ymax = 999999.)
{
  TH1D *xproj;
  TH2F *hist2 = 0;
  Int_t ok;
  Int_t ybin1, ybin2;

  if (c1 == NULL)
    mkcanvas ();

  /* get spectrum */

  ok = 1;
  if (mfile != 0)
    hist2 = (TH2F *) mfile->Get (histname);
  else if (dfile != 0)
    hist2 = (TH2F *) gROOT->FindObject (histname);
  else
    {
      printf ("no shared mem or disk file loaded!\n");
      ok = 0;
    };

  /* project */

  if (ok)
    {
      xproj = (TH1D *) gROOT->FindObject (pname);
      printf ("%p; ", xproj);
      if (xproj != NULL)
        xproj->Delete ();

      if (ymin == -999999. && ymax == 999999.)
        {
          ybin2 = hist2->GetNbinsY ();
          ybin1 = 0;
        }
      else if (ymax == 999999. && ymin != -999999.)
        {
          ybin1 = hist2->GetYaxis ()->FindBin (ymin);
          ybin2 = ybin1;
        }
      else
        {
          ybin1 = hist2->GetYaxis ()->FindBin (ymin);
          ybin2 = hist2->GetYaxis ()->FindBin (ymax);
        }
      printf ("projection from bin %i to bin %i\n", ybin1, ybin2);
      fflush (stdout);
      hist2->ProjectionX (pname, ybin1, ybin2)->Draw ();

      /* clean up */

      if (hist2 != NULL)
        hist2->Delete ();

    }
}

/*-----------------------------------------------------------------*/

void
pj3z (Char_t * histname, Char_t * pname, Float_t xmin = -999999., Float_t xmax = 999999., Float_t ymin = -999999., Float_t ymax = 999999.)
{
  TH1D *xproj;
  TH3F *hist2 = 0;
  Int_t ok;
  Int_t xbin1=0, xbin2=0;
  Int_t ybin1=0, ybin2=0;

  if (c1 == NULL)
    mkcanvas ();

  /* get spectrum */

  ok = 1;
  if (mfile != 0)
    hist2 = (TH3F *) mfile->Get (histname);
  else if (dfile != 0)
    hist2 = (TH3F *) gROOT->FindObject (histname);
  else
    {
      printf ("no shared mem or disk file loaded!\n");
      ok = 0;
    };

  /* project */

  if (ok)
    {
      xproj = (TH1D *) gROOT->FindObject (pname);
      printf ("%p; ", xproj);
      if (xproj != NULL)
        xproj->Delete ();

      if (ymin == -999999. && ymax == 999999.)
        {
          ybin2 = hist2->GetNbinsY ();
          ybin1 = 0;
        }
      else if (ymax == 999999. && ymin != -999999.)
        {
          ybin1 = hist2->GetYaxis ()->FindBin (ymin);
          ybin2 = ybin1;
        }
      else
        {
          xbin1 = hist2->GetXaxis ()->FindBin (xmin);
          xbin2 = hist2->GetXaxis ()->FindBin (xmax);
          ybin1 = hist2->GetYaxis ()->FindBin (ymin);
          ybin2 = hist2->GetYaxis ()->FindBin (ymax);
        }
      printf ("projection from bin %i to bin %i\n", xbin1, xbin2);
      printf ("projection from bin %i to bin %i\n", ybin1, ybin2);
      fflush (stdout);
      hist2->ProjectionZ (pname, xbin1, xbin2, ybin1, ybin2)->Draw ();

      /* clean up */

      if (hist2 != NULL)
        hist2->Delete ();

    }
}

/*-----------------------------------------------------------------*/

#if(0)
DOES NOT WORK ON LINUX void
checkmax (Char_t * histname)
{
  TH2F *hist2 = 0;
  Int_t ok, i, j;
  Int_t ybin1, ybin2;
  Float_t chcont;
  Double_t maxcount, maxmax, d1;

  printf ("FSIGNIF: %i\n", FSIGNIF);
  maxmax = pow ((double) 2, (double) FSIGNIF);
  printf ("max counts: %d\n", maxmax);

  /* get spectrum */

  ok = 1;
  if (mfile != 0)
    hist2 = (TH2F *) mfile->Get (histname);
  else if (dfile != 0)
    hist2 = (TH2F *) gROOT->FindObject (histname);
  else
    {
      printf ("no shared mem or disk file loaded!\n");
      ok = 0;
    };

  maxcount = 0;
  if (ok)
    {

      ybin2 = hist2->GetNbinsY ();
      ybin1 = 0;
      for (i = 0; i < ybin2; i++)
        for (j = 0; j < ybin2; j++)
          {
            chcont = hist2->GetBinContent (i, j);
            if (chcont > maxcount)
              maxcount = chcont;
          };

      printf ("\n");
      printf ("max count is %f\n", (float) maxcount);
      printf ("\n");
      d1 = 100 * maxcount / maxmax;
      printf ("used %8.4f%% of bitspace\n", float (d1));
      printf ("(unless you have already wrapped...)\n");
      printf ("\n");

    };

}

#endif

/*-----------------------------------------------------------------*/

void
pjy (char *histname, char *pname, Float_t xmin = -999999., Float_t xmax = 999999.)
{

  TH1D *yproj;
  TH2F *hist2 = 0;
  Int_t ok;
  Int_t xbin1, xbin2;

  if (c1 == NULL)
    mkcanvas ();

  /* get spectrum */

  ok = 1;
  if (mfile != 0)
    hist2 = (TH2F *) mfile->Get (histname);
  else if (dfile != 0)
    hist2 = (TH2F *) gROOT->FindObject (histname);
  else
    {
      printf ("no shared mem or disk file loaded!\n");
      ok = 0;
    }

  /* project */

  if (ok)
    {

      yproj = (TH1D *) gROOT->FindObject (pname);
//      printf ("%p; ", yproj);
      if (yproj != NULL)
        yproj->Delete ();

      if (xmin == -999999. && xmax == 999999.)
        {
          xbin2 = hist2->GetNbinsX ();
          xbin1 = 0;
        }
      else if (xmax == 999999. && xmin != -999999.)
        {
          xbin1 = hist2->GetXaxis ()->FindBin (xmin);
          xbin2 = xbin1;
        }
      else
        {
          xbin1 = hist2->GetXaxis ()->FindBin (xmin);
          xbin2 = hist2->GetXaxis ()->FindBin (xmax);
        }
//      printf ("projection from bin %i to bin %i\n", xbin1, xbin2);
      fflush (stdout);
      hist2->ProjectionY (pname, xbin1, xbin2)->Draw ();
    }
  /* clean up */

  if (hist2 != NULL)
    hist2->Delete ();

}

/*-----------------------------------------------------------------*/

void
gate (Char_t * histname, Char_t * pname, Float_t p1, Float_t w1, Double_t bgf = 1.0)
{

  /* extract a gate and subtract a fraction */
  /* of the total projection as background. */

  /* declarations */

  Float_t ylo, yhi;
  TH1D *gpro, *tpro;
  Int_t ok;
  Int_t yblo, ybhi;
  Double_t gc, tc, bgff;
  TH2F *whist = NULL;

  /* convert to lo hi channels */

  ylo = p1 - w1;
  yhi = p1 + w1;
  printf ("gate: %f to %f\n", ylo, yhi);

  /* get spectrum */

  ok = 1;
  if (mfile != 0)
    whist = (TH2F *) mfile->Get (histname);
  else if (dfile != 0)
    whist = (TH2F *) gROOT->FindObject (histname);
  else
    {
      printf ("no shared mem or disk file loaded!\n");
      ok = 0;
    }

  /* do the projection */

  if (ok)
    {

      gpro = (TH1D *) gROOT->FindObject (pname);
      printf ("%p; ", gpro);
      if (gpro != NULL)
        gpro->Delete ();
      tpro = (TH1D *) gROOT->FindObject ("tpro");
      printf ("%p; \n", tpro);
      if (tpro != NULL)
        tpro->Delete ();

      /* find bin limits */

      yblo = whist->GetYaxis ()->FindBin (ylo);
      ybhi = whist->GetYaxis ()->FindBin (yhi);

      /* project gate */

      printf ("projecting from %4i to %4i, ", yblo, ybhi);
      whist->ProjectionX (pname, yblo, ybhi);
      TH1 *hgpro = (TH1 *) gDirectory->Get (pname);
      hgpro->Draw ();
      //hgpro->Sumw2();
      gc = hgpro->Integral ();
      printf ("count: %10i\n", (int) gc);

      /* total projection */

      yblo = 0;
      ybhi = whist->GetNbinsY ();
      printf ("projecting from %4i to %4i, ", yblo, ybhi);
      fflush (stdout);
      whist->ProjectionX ("tpro", yblo, ybhi);
      TH1 *htpro = (TH1 *) gDirectory->Get ("tpro");
      htpro->Draw ();
      //htpro->Sumw2();
      tc = htpro->Integral ();
      printf ("count: %10i\n", (int) tc);
      fflush (stdout);

      /* subtract background from gate */

      bgff = -bgf * gc / tc;
      printf ("background factor: %f\n", (float) bgff);
      fflush (stdout);
      hgpro->Add (htpro, bgff);

      /* draw background subtracted gate */

      hgpro->Draw ();

      /* clean up */

      if (whist != NULL)
        whist->Delete ();
#if(0)
      if (gpro != NULL)
        gpro->Delete ();
      if (tpro != NULL)
        tpro->Delete ();
#endif

    };

}

/*-----------------------------------------------------------------*/

void
close ()
{

  /* close chared memory file or disk file */

  if (mfile != 0)
    {
      mfile->Close ();
      mfile = 0;
      printf ("closed/detached shared memory\n");
    }
  else if (dfile != 0)
    {
      dfile->Close ();
      dfile = 0;
      printf ("closed disk file\n");
    }
};

/*-----------------------------------------------------------------*/

void
sload (Char_t * mapname, Long_t startadr = 0)
{
  close ();

  /* attach shared memory */
  /* mfile is global...... */

  /* set start address */

  if (startadr > 0)
    TMapFile::SetMapAddress (startadr);

  /* attach mapfile, use default read mode */

  mfile = TMapFile::Create (mapname);

  /* print standard info */

  mfile->Print ();
  printf ("pointer name: mfile\n");

}

/*-----------------------------------------------------------------*/

void
dload (Char_t * fname)
{
  close ();

  /* load disk file. */
  /* dfile is global */

  dfile = TFile::Open (fname);

  /* print standard info */

  dfile->Print ();
  printf ("pointer name: dfile\n");

}

/*-----------------------------------------------------------------*/

void
sdummyload (Long_t size)
{

  /* dummy load a shared memory to find out what */
  /* start address it chooses................... */

  TMapFile *m;
  m = TMapFile::Create ("dummy.map", "recreate", size);
  m->Print ();

  /* close and remove dummy map file */

  m->Close ();
  gSystem->Exec ("\\rm dummy.map");

}

/*-----------------------------------------------------------------*/

void
ls ()
{
  if (mfile != 0)
    {
      mfile->Print ();
      mfile->ls ();
    }
  else if (dfile != 0)
    {
      dfile->Print ();
      dfile->ls ();
    }
  else
    printf ("no shared mem or disk file loaded!\n");
}

/*-----------------------------------------------------------------*/

void
d2 (Char_t * histname, Float_t xmin = -999999., Float_t xmax = 999999., Float_t ymin = -999999., Float_t ymax = 999999.)
{
  /* display 2D spectrum */

  /* declarations */

  Int_t minxbin, maxxbin;
  Int_t lowxbin, highxbin;
  Int_t minybin, maxybin;
  Int_t lowybin, highybin;
  static TH2F *hist2 = NULL;
  Int_t ok;

  /* clean up */

//  if (hist2 != NULL)
//    hist2->Delete ();
//  ^^^^^ found this crashes program when
//        we switch bt d1 and d2

  /* get spectrum */

  ok = 1;
  if (mfile != 0)
    {
      hist2 = (TH2F *) mfile->Get (histname);
      if (hist2 == NULL)
        {
          printf ("d2: could not find the histogram <%s> in mem file!\n", histname);
          fflush (stdout);
          ok = 0;
        }
    }
  else if (dfile != 0)
    {
      hist2 = (TH2F *) gROOT->FindObject (histname);
      if (hist2 == NULL)
        {
          printf ("d2: could not find the histogram <%s> in disk file!\n", histname);
          fflush (stdout);
          ok = 0;
        }
    }
  else
    {
      printf ("no shared mem or disk file loaded!\n");
      ok = 0;
    }

  /* display spectrum */

  if (ok)
    {

      /* quietly make canvas if it is not there */

      if (c1 == NULL)
        mkcanvas ();

      lowxbin = hist2->GetXaxis ()->GetFirst ();
      highxbin = hist2->GetXaxis ()->GetLast ();
      lowybin = hist2->GetYaxis ()->GetFirst ();
      highybin = hist2->GetYaxis ()->GetLast ();

      if (xmin == -999999. && xmax == 999999.)
        {
          hist2->GetXaxis ()->SetRange (lowxbin, highxbin);
        }
      else if (xmax == 999999. && xmin != -999999.)
        {
          minxbin = hist2->GetXaxis ()->FindBin (xmin);
          hist2->GetXaxis ()->SetRange (minxbin, highxbin);
        }
      else
        {
          minxbin = hist2->GetXaxis ()->FindBin (xmin);
          maxxbin = hist2->GetXaxis ()->FindBin (xmax);
          hist2->GetXaxis ()->SetRange (minxbin, maxxbin);
        }
      if (ymin == -999999. && ymax == 999999.)
        {
          hist2->GetYaxis ()->SetRange (lowybin, highybin);
        }
      else if (ymax == 999999. && ymin != -999999.)
        {
          minybin = hist2->GetYaxis ()->FindBin (ymin);
          hist2->GetYaxis ()->SetRange (minybin, highybin);
        }
      else
        {
          minybin = hist2->GetYaxis ()->FindBin (ymin);
          maxybin = hist2->GetYaxis ()->FindBin (ymax);
          hist2->GetYaxis ()->SetRange (minybin, maxybin);
        };

      hist2->Draw ("col2");
    };
}

/*-----------------------------------------------------------------*/

void
d3 (Char_t * histname)
{
  /* display 2D spectrum */

  /* declarations */

  Int_t lowxbin, highxbin;

  Int_t lowybin, highybin;

  Int_t lowzbin, highzbin;

  static TH3F *hist3 = NULL;
  Int_t ok;

  /* clean up */

  /* get spectrum */

  ok = 1;
  if (mfile != 0)
    {
      hist3 = (TH3F *) mfile->Get (histname);
      if (hist3 == NULL)
        {
          printf ("d3: could not find the histogram <%s> in mem file!\n", histname);
          fflush (stdout);
          ok = 0;
        }
    }
  else if (dfile != 0)
    {
      hist3 = (TH3F *) gROOT->FindObject (histname);
      if (hist3 == NULL)
        {
          printf ("d3: could not find the histogram <%s> in disk file!\n", histname);
          fflush (stdout);
          ok = 0;
        }
    }
  else
    {
      printf ("no shared mem or disk file loaded!\n");
      ok = 0;
    }

  /* display spectrum */

  if (ok)
    {

      /* quietly make canvas if it is not there */

      if (c1 == NULL)
        mkcanvas ();

      lowxbin = hist3->GetXaxis ()->GetFirst ();
      highxbin = hist3->GetXaxis ()->GetLast ();
      lowybin = hist3->GetYaxis ()->GetFirst ();
      highybin = hist3->GetYaxis ()->GetLast ();
      lowzbin = hist3->GetZaxis ()->GetFirst ();
      highzbin = hist3->GetZaxis ()->GetLast ();

#if(0)
      if (xmin == -999999. && xmax == 999999.)
        {
          hist2->GetXaxis ()->SetRange (lowxbin, highxbin);
        }
      else if (xmax == 999999. && xmin != -999999.)
        {
          minxbin = hist2->GetXaxis ()->FindBin (xmin);
          hist2->GetXaxis ()->SetRange (minxbin, highxbin);
        }
      else
        {
          minxbin = hist2->GetXaxis ()->FindBin (xmin);
          maxxbin = hist2->GetXaxis ()->FindBin (xmax);
          hist2->GetXaxis ()->SetRange (minxbin, maxxbin);
        }
      if (ymin == -999999. && ymax == 999999.)
        {
          hist2->GetYaxis ()->SetRange (lowybin, highybin);
        }
      else if (ymax == 999999. && ymin != -999999.)
        {
          minybin = hist2->GetYaxis ()->FindBin (ymin);
          hist2->GetYaxis ()->SetRange (minybin, highybin);
        }
      else
        {
          minybin = hist2->GetYaxis ()->FindBin (ymin);
          maxybin = hist2->GetYaxis ()->FindBin (ymax);
          hist2->GetYaxis ()->SetRange (minybin, maxybin);
        };
#endif

      hist3->Draw ("col2");
    };
}

/*----------------------------------------------------------------*/

void
mk2dwin (Char_t * histname, Float_t xmin = -999999.,
         Float_t xmax = 999999., Float_t ymin = -999999., Float_t ymax = 999999.)
{

  /* declarations */

//  TCutG *mycutg;

  /* delete old 2dwin */

  mycutg = (TCutG *) gROOT->GetListOfSpecials ()->FindObject ("CUTG");
  if (mycutg != NULL)
    mycutg->Delete ();

  /* quietly set up canvas */

  mkcanvas ();

  /* display 2D spectrum */

  d2 (histname, xmin, xmax, ymin, ymax);

  /* ask user to make 2D window */

  gROOT->SetEditorMode ("CutG");
  c1->Update ();

  printf ("-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+\n");
  printf ("o use left hand button to mark the 2D window\n");
  printf ("o double click on the last point to complete 2D window\n");
  printf ("o use save2dwin(winname) to save this window to file\n");
  printf ("-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+\n");

}

/*-----------------------------------------------------------------*/

void
save2dwin (Char_t * winname)
{
  /* rename TCutG window and save to file */

  /* declarations */

  char str[STRLEN];

  /* find CUTG window */

  mycutg = (TCutG *) gROOT->GetListOfSpecials ()->FindObject ("CUTG");
  if (mycutg == NULL)
    {
      printf ("no 2D window present!\n");
      printf ("use mk2dwin to make a window you can save to a file\n");
      return;
    };
  mycutg->SetName (winname);
  mycutg->Print ();

  /* write it to a file */

  TFile *f = new TFile (winname, "recreate");
  mycutg->Write (winname);
  f->Print ();
  f->ls ();
  f->Close ();
  printf ("2dwin written to file:\n");
  sprintf (str, "ls -l %s", winname);
  gSystem->Exec (str);

  /* get rid of the 2dwin */

  mycutg->Delete ();
  mycutg = 0;
  printf ("CUTG has been deleted\n");

  /* done */

}

/*----------------------------------------------------------------*/

void
d2dwin (Char_t * winname, Char_t * histname, Float_t xmin = -999999.,
        Float_t xmax = 999999., Float_t ymin = -999999., Float_t ymax = 999999.)
{
  /* read 2dwin from file and display */
  /* it on top of a 2d histogram      */

  /* declarations */

//  TCutG *mycutg = NULL;

  /* warn if nothing is loaded */

  if (mfile == 0 && dfile == 0)
    printf ("warning: no disk file or shared memory loaded\n");

  /* quietly set up canvas if it is not there */

  if (c1 == NULL)
    mkcanvas ();

  /* display 2D spectrum */

  d2 (histname, xmin, xmax, ymin, ymax);
  c1->Update ();

  /* read 2dwin */

  mycutg = rd2dwin (winname);

  /* display if we can */

  if (mycutg != NULL)
    mycutg->Draw ();

  /* done */

  return;

}

/*-----------------------------------------------------------------*/

void
lsrootfile (Char_t * histname)
{

  /* list the content of a root file */

  /* declarations */

  TFile *f = NULL;
  TKey *key;
  TObject *obj;

  /* open rootfile */

  f = new TFile (histname, "read");

  /* list (this is really equiv to ls() ) */

  TIter nextkey (f->GetListOfKeys ());
  while ((key = (TKey *) nextkey ()))
    {
      obj = key->ReadObj ();
      obj->Print ();
    };

  /* done */

  f->Close ();

}

/*--------------------------------------------------------------*/
//does not work since spectra not read in...

#if(0)
void
allTH1tospe (Char_t * histname)
{

  /* declarations */

  TFile *f = NULL;
  TH1 *htemp;
  TKey *key;
  TObject *obj;
  char str[30];

  /* open rootfile */

  f = new TFile (histname, "read");

  /* list (this is really equiv to ls() ) */

  TIter nextkey (f->GetListOfKeys ());
  while ((key = (TKey *) nextkey ()))
    {
      obj = key->ReadObj ();
      if (obj->IsA ()->InheritsFrom ("TH1"))
        {
          obj->Print ();
          htemp = (TH1 *) obj;
          sprintf (str, "%s\n", htemp->GetName ());
          wrspe (str, str);
        };
    };

  /* done */

  f->Close ();

}
#endif

/*-----------------------------------------------------------------*/

#if(0)
Double_t
FITft1 (Double_t * x, Double_t * par)
{
  Double_t arg = 0;
  if (par[2] != 0)
    arg = (x[0] - par[1]) / par[2];
  Double_t fitval = par[0] * TMath::Exp (-0.5 * arg * arg) + par[3] + par[4] * x[0] + par[5] * x[0] * x[0];
  return fitval;
}

void
ft1 (Char_t * HName, Float_t xLO, Float_t xHI, Float_t xPK = 0)
{
  /* fit single peak with a quadratic background */

  /* declarations */

  TH1D *hist1 = 0;

  /* get spectrum */

  hist1 = 0;
  hist1 = (TH1D *) gROOT->FindObject (HName);

  /* fit */

  if (hist1 != 0)
    {
      TF1 *func = new TF1 ("FITft1", FITft1, xLO, xHI, 6);
      func->SetParameters (hist1->GetMaximum (), xPK, 1, hist1->GetMinimum (), 0, 0);
      func->SetParNames ("Constant", "Mean_value", "Sigma", "a0", "b0", "c0");
      hist1->Fit (func, "r");
    };

}

#endif

void
wrRadmat (Char_t * histname)
{

  /* write a spectrum out in Radford format */
  /* Darek contribution 8/03                */

  TH2D *hist2;
  Int_t i, j, nx, ny, size;
  unsigned short int *sp;
  char str[STRLEN];

  /* get spectrum */

  hist2 = FindMatrix (histname);

  if (hist2 != NULL)
    {

      /* get x and y dimension information */

      nx = hist2->GetXaxis ()->GetNbins ();
      printf (" %i channels in x\n", nx);
      fflush (stdout);

      ny = hist2->GetYaxis ()->GetNbins ();
      printf (" %i channels in y\n", ny);
      fflush (stdout);

      size = 4096 * 4096 * sizeof (unsigned short int);
      sp = (unsigned short int *) malloc (size);

      /* fill short int matrix */

      for (i = 0; i < 4096; i++)
        for (j = 0; j < 4096; j++)
          if (i < nx && j < ny)
            {
              *(sp + i * 4096 + j) = (unsigned short int) hist2->GetBinContent (i, j);
            }
          else
            *(sp + i * 4096 + j) = 0;

      /* write mat file */

      sprintf (str, "%s.mat", histname);
      wr_Radmat (str, sp);
      printf ("wrote %i X %i channels to %s\n", nx, ny, str);

      free (sp);

    }
  else
    printf ("spectrum %s not found\n", histname);

  /* done */

}







/* ================================================================ */
/* FUNCTIONS not called directly */

/*--------------------------------------------------------*/

TH1D *
FindSpectrum (Char_t * hnam)
{

  /* this function will search both any shared memory */
  /* as well as local memory for TH1D spectrum hnam */

  /* declarations */

  Int_t FoundSpectrum = 0;
  TH1D *hpointer = NULL;

  /* warning */

  if (mfile == 0 && dfile == 0)
    printf ("warning: no disk file or shared memory available\n");

  /* try finding the spectrum in shared memory */

  if (mfile != 0)
    {

      hpointer = (TH1D *) mfile->Get (hnam);

      /* if it wasn't there, it could be local */

      if (hpointer == NULL)
        FoundSpectrum = 0;
      else
        FoundSpectrum = 1;

    };

  /* try finding spectrum in disk file or locally */

  if (FoundSpectrum == 0 || dfile != 0)
    {

      hpointer = (TH1D *) gROOT->FindObject (hnam);

      if (hpointer == NULL)
        FoundSpectrum = 0;
      else
        FoundSpectrum = 1;

    };

  /* return pointer or error message */

  if (FoundSpectrum)
    return (hpointer);
  else
    return (NULL);

}

/*-------------------------------------------------------------------*/

void
findGT_cal_SEG (Char_t * source, Char_t * newf)
{
  /* declarations */

  FILE *fp;
  int i,j;
  float off[MAXSEGINDEX], gain[MAXSEGINDEX];
  char *p, str[128], str1[128];
  TH1D *hist1 = 0;
  Stat_t counts;

  printf("will calibrate the GT segment energies for a \"%s\" source\n",source);
  printf("will write the calibration coefficients to the file \"%s\"\n",newf);

  if(1)
    {
    printf("this function is under development\n");
    return;
    }

  /* open cal file */

  fp = fopen (newf, "r");
  if (fp != NULL)
    {
      printf ("\"%s\" already exist! please delete or move and try again\n", newf);
//      return;
    }
  fp = fopen (newf, "w");
  if (fp == NULL)
    {
      printf ("could not open: \"%s\"\n", newf);
      return;
    }
  else
    printf ("\"%s\" open\n", newf);

/*HEREWEARE*/

  /* loop through the SEGe 2D matrix */

//  for (i=0;i<=MAXDETPOS;i++)
    for (i=20;i<21;i++)
    for (j=0;j<= NUMSEGS;j++)
      {
      printf("det %3i seg %2i, index %4i\n", i,j,i*NUMSEGS+j);

      /* make a projected spectrum */

      sprintf (str, "SegE");
      sprintf (str1, "x");
      pjy (str, str1,i*NUMSEGS+j ,i*NUMSEGS+j );

      /* look at projection */

      sprintf (str, "x");
      hist1 = FindSpectrum (str1);
      counts = hist1->Integral ();

      if ((int) counts > 10)
        {

          printf ("crystal %3i has %8i counts\n", i, (int) counts);


      /* find the two peaks */

      /* find the calibration */

         };

      };

  /*done */

  fclose(fp);
  return;

};

/*-------------------------------------------------------------------*/

void
findGT_cal_CC (Char_t * source, Char_t * newf)
{

  /* declarations */

  int i, i1, i2, j, nn = 0;
  TH1D *hist1 = 0;
  Stat_t counts;
  float r1, r2, max = 0, thr;
  Double_t dplo, dphi, kevch = 1, aphi, aplo, ww, w1 = 3, mean;
  char *p, str[128], str1[128];
  int minchan = 100;
  Int_t maxch, ihi, ilo;
  float off[MAXDETPOS + 1], gain[MAXDETPOS + 1];
  FILE *fp;

  /* hello */

  printf ("entered findGT_cal_CC function\n");
  if (c1 == NULL)
    mkcanvas ();

  /* open cal file */

  fp = fopen (newf, "r");
  if (fp != NULL)
    {
      printf ("\"%s\" already exist! please delete or move and try again\n", newf);
      return;
    }
  fp = fopen (newf, "w");
  if (fp == NULL)
    {
      printf ("could not open: \"%s\"\n", newf);
      return;
    }
  else
    printf ("\"%s\" open\n", newf);


  /* find desired peak positions for selected source */

  if ((p = strstr (source, "207Bi")) != NULL)
    {
      r1 = BI207LO;
      r2 = BI207HI;
      printf ("assuming 207Bi source calibration\n");
      printf ("__with peaks at %7.3f and %7.3f keV!!!!!\n", r1, r2);
      fflush (stdout);

      printf ("will calibrate for %9.6f kev/channel\n", (float) kevch);
      dplo = BI207LO / kevch;
      dphi = BI207HI / kevch;
      printf ("==>desired peak positions: %9.3f and %9.3f\n", (float) dplo, (float) dphi);
      fflush (stdout);
    }
  else if ((p = strstr (source, "88Y")) != NULL)
    {
      r1 = Y88LO;
      r2 = Y88HI;
      printf ("assuming 88Y source calibration\n");
      printf ("__with peaks at %7.3f and %7.3f keV!!!!!\n", r1, r2);
      fflush (stdout);

      printf ("will calibrate for %9.6f kev/channel\n", (float) kevch);
      dplo = Y88LO / kevch;
      dphi = Y88HI / kevch;
      printf ("==>desired peak positions: %9.3f and %9.3f\n", (float) dplo, (float) dphi);
      fflush (stdout);
    }
  else if ((p = strstr (source, "60Co")) != NULL)
    {
      r1 = CO60LO;
      r2 = CO60HI;
      printf ("assuming 60Co source calibration\n");
      printf ("__with peaks at %7.3f and %7.3f keV!!!!!\n", r1, r2);
      fflush (stdout);

      printf ("will calibrate for %9.6f kev/channel\n", (float) kevch);
      dplo = CO60LO / kevch;
      dphi = CO60HI / kevch;
      printf ("==>desired peak positions: %9.3f and %9.3f\n", (float) dplo, (float) dphi);
      fflush (stdout);
    }
  else
    {
      printf ("***error, source \"%s\" not recognized\n", source);
      printf ("   must be one of: 207Bi, 88Y or 60Co\n");
      printf ("\n");
      return;
    };

  for (i = 1; i <= MAXDETPOS; i++)
    {

      /* make projection */

//     printf("i=%i\n", i);
//     fflush(stdout);
      sprintf (str, "CCe");
      sprintf (str1, "x");
      pjy (str, str1, i, i);

      /* look at projection */

      sprintf (str, "x");
      hist1 = FindSpectrum (str);
      counts = hist1->Integral ();
      if ((int) counts > 10)
        {

          printf ("crystal %3i has %8i counts\n", i, (int) counts);


          /* find max channel (avoid lowest and highest) */

          i1 = hist1->GetXaxis ()->GetFirst ();
          i1 += 20;
          if (i1 < minchan)
            i1 = minchan;
          i2 = hist1->GetXaxis ()->GetLast ();
          i2 -= 20;
          max = 0;
//          printf("i1=%i, i2=%i\n", i1,i2);
          for (j = i1; j < i2; j++)
            if (hist1->GetBinContent (j) > max)
              {
                max = hist1->GetBinContent (j);
//              printf("%i -- %f\n", j, max);
              };
          printf ("crystal %3i has a max count of %9.0f\n", i, max);

          /* set threshold */

          thr = max / 3.0;

          /* find approx hi  peak pos */

          ihi = i2;
          while (hist1->GetBinContent (ihi) < thr && ihi > i1)
            ihi--;

          /* find max hi peak positions */

          max = 0;
          maxch = 0;
          for (j = ihi - 10; j < ihi + 10; j++)
            if (hist1->GetBinContent (j) > max)
              {
                max = hist1->GetBinContent (j);
                maxch = j;
              };
          ihi = maxch;

          /* find mean hi peak position */

          mean = 0;
          ww = 0;
          for (j = ihi - w1; j < ihi + w1; j++)
            {
              mean += hist1->GetBinContent (j) * j;
              ww += hist1->GetBinContent (j);
            };
          mean /= ww;
          aphi = mean;

          /* find approx lo  peak pos */
          /* (continue down...)       */

          ilo = ihi - 30;
          while (hist1->GetBinContent (ilo) < thr && ilo > i1)
            ilo--;

          /* find max hlo peak positions */

          max = 0;
          maxch = 0;
          for (j = ilo - 10; j < ilo + 10; j++)
            if (hist1->GetBinContent (j) > max)
              {
                max = hist1->GetBinContent (j);
                maxch = j;
              };
          ilo = maxch;

          /* find mean lo peak position */

          mean = 0;
          ww = 0;
          for (j = ilo - w1; j < ilo + w1; j++)
            {
              mean += hist1->GetBinContent (j) * j;
              ww += hist1->GetBinContent (j);
            };

          mean /= ww;
          aplo = mean;

          printf ("(%7.2f,%7.2f); ", aplo, aphi);

          gain[i] = (dphi - dplo) / (aphi - aplo);
          off[i] = dphi - gain[i] * aphi;

          /* report and store */

          printf ("off= %6.2f, gain= %7.5f\n", off[i], gain[i]);
          fprintf (fp, "%3i %f %f\n", i, off[i], gain[i]);
          nn++;

          printf ("check: aplo*gain[i]+off[i]=%f\n", aplo * gain[i] + off[i]);
          printf ("check: aphi*gain[i]+off[i]=%f\n", aphi * gain[i] + off[i]);

          /* plot */

          sprintf (str, "x");
          d1 (str, 500, 1100);

        };

    };


  /* done */

  fclose (fp);
  printf (" determined %i calibration coefficients\n", nn);
  printf ("exited findGT_cal_CC function\n");
  return;

}

/*-------------------------------------------------------------------*/

void
SegE_to_spe (Int_t ilo, Int_t ihi)
{

  /* declarations */

  int i, k, i1;
  TH1D *hist1 = 0;
  Stat_t counts;
  char str[128], str1[128];

  i1 = MAXDETPOS * (MAXCRYSTALNO + 1) * 6 * 6;
  for (i = ilo; i < ihi; i++)
    for (k = 0; k < 36; k++)
      {
        i1 = i * 36 + k;
        sprintf (str, "SegE");
        sprintf (str1, "x");
        pjy (str, str1, i1, i1);
        printf ("projecting %3i_%2i\n", i, k);
        fflush (stdout);

        /* look at projection */

        sprintf (str1, "x");
        hist1 = FindSpectrum (str1);
        counts = hist1->Integral ();
        if ((int) counts > 10)
          {
            sprintf (str, "seg_%3.3i_%2.2i.spe", i, k);
            printf ("%3i_%2i has %8i counts --> %s\n", i, k, (int) counts, str);
            sprintf (str1, "x");
            wrspe (str1, str);
            d1 (str1);
          };
      };

}

/*-------------------------------------------------------------------*/

void
ealign (Char_t * GenSpNam, Char_t * oldf, Char_t * newf, Double_t kevch, Double_t sgain, Int_t w1, Char_t * source,
        Int_t minchan)
{

  /* find (new) germanium energy calibration */


  /* declarations */

  Int_t i, j;
  char str[STRLEN], *pc;
  TH1D *hist1;
  Double_t mean, ww;
  Int_t maxch, nn, i1, i2, ihi, ilo;
  FILE *fp;
  float off[MAXGE + 1], gain[MAXGE + 1], r1, r2;
  float oldoff, oldgain;
  Double_t dplo, dphi, max, thr, aphi, aplo;
  Double_t meanrawhi = 0, meanrawlo = 0, d1;
  Int_t nmeanraw = 0;
  char *p;
  int have_spectrum[MAXGE];

  /* check we have shared memory loaded */

#if(0)
  if (mfile == 0)
    {
      printf ("\nyou must load shared memory with sload first\n\n");
      return;
    };
#endif

  /* find desired peak positions for selected source */

  if ((p = strstr (source, "207Bi")) != NULL)
    {
      r1 = BI207LO;
      r2 = BI207HI;
      printf ("assuming 207Bi source calibration\n");
      printf ("__with peaks at %7.3f and %7.3f keV!!!!!\n", r1, r2);
      fflush (stdout);

      printf ("will calibrate for %9.6f kev/channel\n", (float) kevch);
      dplo = BI207LO / kevch;
      dphi = BI207HI / kevch;
      printf ("==>desired peak positions: %9.3f and %9.3f\n", (float) dplo, (float) dphi);
      fflush (stdout);
    }
  else if ((p = strstr (source, "88Y")) != NULL)
    {
      r1 = Y88LO;
      r2 = Y88HI;
      printf ("assuming 88Y source calibration\n");
      printf ("__with peaks at %7.3f and %7.3f keV!!!!!\n", r1, r2);
      fflush (stdout);

      printf ("will calibrate for %9.6f kev/channel\n", (float) kevch);
      dplo = Y88LO / kevch;
      dphi = Y88HI / kevch;
      printf ("==>desired peak positions: %9.3f and %9.3f\n", (float) dplo, (float) dphi);
      fflush (stdout);
    }
  else if ((p = strstr (source, "60Co")) != NULL)
    {
      r1 = CO60LO;
      r2 = CO60HI;
      printf ("assuming 60Co source calibration\n");
      printf ("__with peaks at %7.3f and %7.3f keV!!!!!\n", r1, r2);
      fflush (stdout);

      printf ("will calibrate for %9.6f kev/channel\n", (float) kevch);
      dplo = CO60LO / kevch;
      dphi = CO60HI / kevch;
      printf ("==>desired peak positions: %9.3f and %9.3f\n", (float) dplo, (float) dphi);
      fflush (stdout);
    }
  else
    {
      printf ("***error, source \"%s\" not recognized\n", source);
      printf ("   must be one of: 207Bi, 88Y or 60Co\n");
      printf ("\n");
      return;
    };

  printf ("assuming a data multipication factor of %9.3f in GSSort\n", (float) sgain);
  fflush (stdout);

  printf ("will read old calibration parameters from: %s\n", oldf);
  fflush (stdout);
  printf ("and write new calibration parameters to  : %s\n", newf);
  fflush (stdout);

  printf ("will use +/- %i channels around peak to find mean\n", w1);
  fflush (stdout);

  printf ("will look from channel %i (minchan) and up\n", minchan);
  fflush (stdout);

  /*-------------------------------*/
  /* read in old calibrations */
  /*-------------------------------*/

  if ((p = strstr (oldf, "DUMMY")) != NULL)
    {
      for (i = 0; i <= MAXGE; i++)
        {
          off[i] = 0;
          gain[i] = 1;
        };
    }
  else
    {
      fp = fopen (oldf, "r");
      if (fp == NULL)
        {
          printf ("could not open: %s\n", oldf);
          return;
        }
      else
        printf ("%s open\n", oldf);

      /* read first line (ID like) */

      pc = fgets (str, STRLEN, fp);
      printf ("%s\n", str);

      /* read old offsets */

      nn = 0;
      pc = fgets (str, STRLEN, fp);
      while (pc != NULL)
        {
          nn++;
          sscanf (str, "%i %f %f", &i1, &r1, &r2);
          off[i1] = r1;
          gain[i1] = r2;
          pc = fgets (str, STRLEN, fp);
        };
      printf ("read %i calibration parameters\n", nn);
      fflush (stdout);

      fclose (fp);

    };

  /*---------------------------------------------*/
  /* loop through detectors and find new offsets */
  /*---------------------------------------------*/

  for (i = 1; i <= MAXGE; i++)
    have_spectrum[i] = 0;

  for (i = 1; i <= MAXGE; i++)
    {

      /* fetch spectrum  */

      sprintf (str, "%s%3.3i", GenSpNam, i);
      if (mfile != 0)
        hist1 = (TH1D *) mfile->Get (str);
      else
        hist1 = (TH1D *) gROOT->FindObject (str);
      if (hist1 != NULL)
        printf ("[%s]; ", str);
      else
        {
          printf ("could not find spectrum %s \n", str);
          /*return; */
        };

      if (hist1 != NULL)
        {

          /* find max channel (avoid lowest and highest) */

          i1 = hist1->GetXaxis ()->GetFirst ();
          i1 += 20;
          if (i1 < minchan)
            i1 = minchan;
          i2 = hist1->GetXaxis ()->GetLast ();
          i2 -= 20;
          max = 0;
          for (j = i1; j < i2; j++)
            if (hist1->GetBinContent (j) > max)
              max = hist1->GetBinContent (j);
          /*printf("max channel bt %i and %i hand counts %f\n", i1,i2,max); */

          if (max > 100.0)
            {

              have_spectrum[i] = 1;

              /* set threshold */

              thr = max / 3.0;

              /* find approx hi  peak pos */

              ihi = i2;
              while (hist1->GetBinContent (ihi) < thr && ihi > i1)
                ihi--;

              /* find max hi peak positions */

              max = 0;
              maxch = 0;
              for (j = ihi - 10; j < ihi + 10; j++)
                if (hist1->GetBinContent (j) > max)
                  {
                    max = hist1->GetBinContent (j);
                    maxch = j;
                  };
              ihi = maxch;

              /* find mean hi peak position */

              mean = 0;
              ww = 0;
              for (j = ihi - w1; j < ihi + w1; j++)
                {
                  mean += hist1->GetBinContent (j) * j;
                  ww += hist1->GetBinContent (j);
                };
              mean /= ww;
              aphi = mean;

              /* find approx lo  peak pos */
              /* (continue down...)       */

              ilo = ihi - 30;
              while (hist1->GetBinContent (ilo) < thr && ilo > i1)
                ilo--;

              /* find max hlo peak positions */

              max = 0;
              maxch = 0;
              for (j = ilo - 10; j < ilo + 10; j++)
                if (hist1->GetBinContent (j) > max)
                  {
                    max = hist1->GetBinContent (j);
                    maxch = j;
                  };
              ilo = maxch;

              /* find mean lo peak position */

              mean = 0;
              ww = 0;
              for (j = ilo - w1; j < ilo + w1; j++)
                {
                  mean += hist1->GetBinContent (j) * j;
                  ww += hist1->GetBinContent (j);
                };

              mean /= ww;
              aplo = mean;

              printf ("(%7.2f,%7.2f)sort;", aplo, aphi);

              /* undo sort compression */

              aplo /= sgain;
              aphi /= sgain;
              printf ("(%7.2f,%7.2f)eff;", aplo, aphi);

              /* undo current calibration */

              aplo = (aplo - off[i]) / gain[i];
              aphi = (aphi - off[i]) / gain[i];
              printf ("(%7.2f,%7.2f)raw\n", aplo, aphi);

              if (!isnan (aplo) && !isnan (aphi))
                {
                  meanrawlo += aplo;
                  meanrawhi += aphi;
                  nmeanraw++;
                };

              /* find new calibration */

              oldoff = off[i];
              oldgain = gain[i];

              gain[i] = (dphi - dplo) / (aphi - aplo);
              off[i] = dphi - gain[i] * aphi;

              /* check for reasonable values */

              if (isnan (gain[i]) || isnan (off[i]))
                {
                  printf ("unreasonable new parameters, set to 0,1 \n");
                  gain[i] = 1.0;
                  off[i] = 0;
                };

              /* report change */

              printf ("________offset change: ");
              printf ("       %6.1f --> %6.1f; ", oldoff, off[i]);
              printf ("%6.1f\n", off[i] - oldoff);

              printf ("________gain change: ");
              printf ("%7.5f --> %7.5f; ", oldgain, gain[i]);
              printf ("%7.5f %% \n", 100 * (gain[i] - oldgain) / oldgain);

            }
          else
            printf ("not enough counts, will not change cal... ");

          printf ("\n");

        };
    };

  /*-----------------------------*/
  /* write new alignment offsets */
  /*-----------------------------*/

  fp = fopen (newf, "w");
  if (fp == NULL)
    {
      printf ("could not open: %s\n", newf);
      return;
    }
  else
    printf ("%s open\n", newf);

  /* write ID line */

//  fprintf (fp, "%s aligment parameters root/GSSort\n", GenSpNam);

  /* write offsets */

  for (i = 1; i <= MAXGE; i++)
    if (have_spectrum[i] == 1)
      fprintf (fp, "%3i %10.4f %10.8f\n", i, (float) off[i], (float) gain[i]);
  fclose (fp);
  printf ("%s energy calibration file %s written\n", GenSpNam, newf);

  /* find 'natural' gain */

  meanrawhi /= nmeanraw;
  meanrawlo /= nmeanraw;
  d1 = (r2 - r1) / (meanrawhi - meanrawlo);
  printf ("natural gain seems to be: %9.5f keV/ch\n", (float) d1);

  /* done */

  fflush (stdout);

}

/*----------------------------------------------------------------*/

void
talign (Char_t * GenSpNam, Char_t * oldf, Char_t * newf, float dpos, int di, int lo, int bp)
{

  /* find (new) germanium/BGO time alignment offsets */

  /* di:   number of channels on earch side of max ch to use */
  /* ch to use for calculating mean peak position      */
  /* dpos: desired channel position in displayed spectrum    */

  /* declarations */

  int i, j, Ok[MAXGE + 1];
  char str[STRLEN], *pc;
  TH1D *hist1;
  double maxval, val, mean, ww, bg;
  int maxch, nn, i1, i2;
  FILE *fp;
  int toff[MAXGE + 1];
  char *p;
  float sp[20000];

  /* check we have shared memory loaded */

  if (!(mfile == 0 || dfile == 0))
    {
      printf ("\nyou must load shared memory or rootfile first\n\n");
      return;
    };

  /* default all offsets to zero */

  for (i = 0; i <= MAXGE; i++)
    toff[i] = 0;

  printf ("generic spectra names: [%s]\n", GenSpNam);
  printf ("old cal file name: [%s]\n", oldf);
  printf ("new cal file name: [%s]\n", newf);
  printf ("will attempt to align times at %9.2f\n", dpos);
  printf ("__+/- %i used for mean\n", di);
  printf ("will not look below %i\n", lo);
  fflush (stdout);

  /* is old file DUMMY? */

  if ((p = strstr (oldf, "DUMMY")) != NULL)
    printf ("will assume old alignment offsets are zero!\n");
  else
    {
    /*-------------------------------*/
      /* read in old alignment offsets */
    /*-------------------------------*/

      fp = fopen (oldf, "r");
      if (fp == NULL)
        {
          printf ("could not open: %s\n", oldf);
          return;
        }
      else
        printf ("%s open\n", oldf);

      /* read first line (ID like) */

      pc = fgets (str, STRLEN, fp);
      printf ("%s\n", str);

      /* read old offsets */

      nn = 0;
      pc = fgets (str, STRLEN, fp);
      while (pc != NULL)
        {
          nn++;
          sscanf (str, "%i %i", &i1, &i2);
          toff[i1] = i2;
          pc = fgets (str, STRLEN, fp);
        };
      printf ("read %i old offsets\n", nn);

      fclose (fp);

    };

  /*---------------------------------------------*/
  /* loop through detectors and find new offsets */
  /*---------------------------------------------*/

  for (i = 1; i <= MAXGE; i++)
    {

      /* fetch spectrum */

      sprintf (str, "%s%3.3i", GenSpNam, i);
      hist1 = FindSpectrum (str);
      if (hist1 != NULL)
        printf ("[%s]%9i; ", str, (int) hist1->GetEntries ());
      else
        printf ("spectrum %s not found\n", str);
      fflush (stdout);

      /* zap local spectrum */

      for (j = 0; j < 20000; j++)
        sp[j] = 0;

      /* transfer time spectrum to local spectrum */

      i1 = hist1->GetXaxis ()->GetFirst ();
      i2 = hist1->GetXaxis ()->GetLast ();
#if(0)
      printf ("fetching from %i to %i\n", i1, i2);
#endif
      for (j = i1; j < i2; j++)
        sp[j] = (float) hist1->GetBinContent (j);

      /* zap out bottom */

      for (j = 0; j < lo; j++)
        sp[j] = 0;

#if(1)
      /* write raw spectrum out */

      sprintf (str, "t_raw_%3.3i.spe", i);
      wr_spe (str, &i2, sp);
#endif

      /* remove spikes */

      if (bp)
        {
          printf ("despiking; ");
          for (j = i1 + 3; j < (i2 - 3); j++)
            {
              bg = sp[j - 3];
              bg += sp[j - 2];
              bg += sp[j + 2];
              bg += sp[j + 3];
              bg = bg / 4;
              if (sp[j] > bp * bg)
                {
                  sp[j] = bg;
                  /*printf("[%i]", j); */
                };
            };

          /* smooth the spectrum a few time */

          for (j = 0; j < 3; j++)
            sm3 (i2, sp);

        };

#if(1)
      /* write processed spectrum out */

      sprintf (str, "t_mod_%3.3i.spe", i);
      wr_spe (str, &i2, sp);
#endif

      /* find max channel */

      maxch = 0;
      maxval = 0;
      if (hist1->GetEntries () > MINCOUNTS)
        for (j = i1 + 5; j < i2 - 5; j++)
          {
            val = sp[j];
            if (val > maxval)
              {
                maxval = sp[j];
                maxch = j;
              };

          };
      if (maxch > 10)
        {
          Ok[i] = 1;
        }
      else
        Ok[i] = 0;

      /* find a mean position of the peak */

      if (Ok[i])
        {

          mean = 0;
          ww = 0;
          if (maxch > 10)
            for (j = maxch - di; j < maxch + di; j++)
              {
                mean += j * hist1->GetBinContent (j);
                ww += hist1->GetBinContent (j);
              };
          mean /= ww;
          printf ("pos %7.2f; ", mean);

          /* find new offset */

          printf ("toff %4i --> ", toff[i]);
          i1 = toff[i];
          toff[i] -= (int) (dpos - mean);
          printf ("%4i ", toff[i]);
          i1 = toff[i] - i1;
          printf (" _change: %4i ", i1);

        };

      printf ("\n");

    };

  /*-----------------------------*/
  /* write new alignment offsets */
  /*-----------------------------*/

  fp = fopen (newf, "w");
  if (fp == NULL)
    {
      printf ("could not open: %s\n", newf);
      return;
    }
  else
    printf ("%s open\n", newf);

  /* write ID line */

  fprintf (fp, "#%s alignment root/GSSort\n", GenSpNam);

  /* write offsets */

  for (i = 1; i <= MAXGE; i++)
    fprintf (fp, "%3i %4i\n", i, toff[i]);

  fclose (fp);

  /* done */

  printf ("%s time alignment file \"%s\" written\n", GenSpNam, newf);
  fflush (stdout);
  return;

}

/*-----------------------------------------------------------------*/

void
dod1 (Char_t * opt, Char_t * histname, Float_t xmin, Float_t xmax, Int_t rescale)
{
  /* display 1D spectrum with options */

  /* rescale=1: overlay plot */

  /* delarations */

  Int_t minbin, maxbin;
  Int_t lowbin, highbin;
  TH1D *hist1 = 0;
  Double_t d1, scale = 1;
  Stat_t counts;

  if (c1 == NULL)
    mkcanvas ();

  /* get spectrum */

  hist1 = FindSpectrum (histname);

  /* display spectrum */

  if (hist1 != NULL)
    {

      hist1->Print ();


      lowbin = hist1->GetXaxis ()->GetFirst ();
      highbin = hist1->GetXaxis ()->GetLast ();
      if (xmin == -999999. && xmax == 999999.)
        {
          hist1->GetXaxis ()->SetRange (lowbin, highbin);
        }
      else if (xmax == 999999. && xmin != -999999.)
        {
          minbin = hist1->FindBin (xmin);
          hist1->GetXaxis ()->SetRange (minbin, highbin);
        }
      else
        {
          minbin = hist1->FindBin (xmin);
          maxbin = hist1->FindBin (xmax);
          hist1->GetXaxis ()->SetRange (minbin, maxbin);
        }


      /* display only if meaning full # entries */

      counts = hist1->Integral ();
      if ((float) counts > 0)
        {
          hist1->SetLineColor (LineColor);

          if (rescale == 1)
            {
              d1 = 1.1 * hist1->GetMaximum ();
              scale = gPad->GetUymax () / d1;
              hist1->Scale (scale);
              printf ("WARNING: rescaled data by %f\n", scale);

            };

          /* overwrite limits ? */

          if (HaveyLimits)
            {
              hist1->SetMaximum (Yhilim);
              hist1->SetMinimum (Ylowlim);
            };
          if (HavexLimits)
            {
              minbin = hist1->FindBin (Xlowlim);
              maxbin = hist1->FindBin (Xhilim);
              hist1->GetXaxis ()->SetRange (minbin, maxbin);
            };

          /* draw with options */

          hist1->Draw (opt);
          c1->Update ();

          /* set back to original scale */

          if (rescale == 1)
            {
              scale = 1.0 / scale;
              hist1->Scale (scale);
            };

        }
      else
        printf ("no counts, not displayed...\n");

    }
  else
    printf ("spectrum %s not found\n", histname);
}

/*----------------------------------------------- */

void
str_decomp (Char_t * str, Int_t dim, Int_t yy[])
{
  /* decompose a detector list into a logical array */

  /* declarations */

  int i, j, pos, nn, ok, lo, hi;
  char lstr[100];

  /* zero array */

  for (i = 0; i < dim; i++)
    yy[i] = 0;

  pos = 0;
  ok = 1;
  while (ok == 1)
    {

      /* search for the next sub range */

      nn = 0;
      while ((str[pos + nn] != ',') && (ok == 1))
        {
          nn++;
          if (str[pos + nn] == 0)
            ok = 0;
        };

      /* create the local string */

      j = 0;
      for (i = pos; i < (pos + nn); i++)
        {
          lstr[j] = str[i];
          j++;
        };
      lstr[j] = 0;

      /* remove non numeric characters in this string */

      for (i = 0; i < nn; i++)
        if ((lstr[i] < 48) || (lstr[i] > 57))
          lstr[i] = ' ';

      /* extract limits */

      i = sscanf (lstr, "%i%i", &lo, &hi);

      /* fill array */

      if (i == 1)
        hi = lo;
      for (j = lo; j <= hi; j++)
        yy[j] = 1;

      /* move pos pointer forward */

      pos += nn;
      pos++;

    };

  /* done */

  return;

}

/*--------------------------------------------------------------*/

TCutG *
rd2dwin (Char_t * winname)
{

  /* declarations */

  char str[STRLEN];
//  TCutG *mycutg;

  TFile *f = new TFile (winname, "read");
  mycutg = (TCutG *) f->Get (winname);
  f->Close ();

  if (mycutg != NULL)
    {
      printf ("2dwin read from file:\n");
      sprintf (str, "ls -l %s", winname);
      gSystem->Exec (str);
      mycutg->Print ();
      fflush (stdout);
    }
  else
    {
      printf ("could not read 2dwin file %s\n", winname);
      mycutg = NULL;
    };

  /* done */

  return (mycutg);

}


/*-----------------------------------------------------------------*/

int
wr_Radmat (char *fn, unsigned short int *sp)
{

  /* declarations */

  int st, i1;
  int exa;
  int siz;

  /* open file */

  exa = creat (fn, PMODE);
  if (exa <= 0)
    {
      printf ("wr_Radmat: could not open spectrum file %s\n", fn);
      return (0);
    };

  /* write spectrum */

  //i1 = *dimx * *dimy * 2;
  i1 = 4096 * 4096 * 2;
  siz = write (exa, (char *) sp, i1);
  if (i1 != siz)
    {
      printf ("tried to write %i\n", i1);
      printf ("actually wrote %i\n", siz);
      close (exa);
      return (-1);
    };

  /* done */

  st = close (exa);
  if (st != 0)
    printf ("could not close [%s]\n", fn);
  return (0);

}

/*--------------------------------------------------------*/

TH2D *
FindMatrix (Char_t * hnam)
{

  /* this function will search both any shared memory */
  /* as well as local memory for TH2D spectrum hnam */

  /* declarations */

  Int_t FoundMatrix = 0;
  TH2D *hpointer = NULL;

  /* warning */

  if (mfile == 0 && dfile == 0)
    printf ("warning: no disk file or shared memory available\n");

  /* try finding the spectrum in shared memory */

  if (mfile != 0)
    {

      hpointer = (TH2D *) mfile->Get (hnam);

      /* if it wasn't there, it could be local */

      if (hpointer == NULL)
        FoundMatrix = 0;
      else
        FoundMatrix = 1;

    };

  /* try finding spectrum in disk file or locally */

  if (FoundMatrix == 0 || dfile != 0)
    {

      hpointer = (TH2D *) gROOT->FindObject (hnam);

      if (hpointer == NULL)
        FoundMatrix = 0;
      else
        FoundMatrix = 1;

    };

  /* return pointer or error message */

  if (FoundMatrix)
    return (hpointer);
  else
    return (NULL);

}

/*------------------------------------------------------*/

int
mkGSSortCmdFile (Char_t * cmd)
{

  /* safe version of command file for GSSort */

  /* declarations */

  FILE *fp;
  int ntry;

  /* first check if cmd file already exists */

  fp = fopen ("GSSort.command", "r");
  if (fp != NULL)
    {

      /*fclose(fp); linux machines get upset if I do this */
      printf ("will wait for previous command to clear");
      fflush (stdout);

      ntry = 0;
      while (fp != NULL && ntry <= CMDTRIES)
        {
          /*fclose(fp); linux machines get upset if I do this */
          ntry++;
          sleep (1);
          printf (".");
          fflush (stdout);
          fp = fopen ("GSSort.command", "r");
        };

      /* check status */

      if (ntry >= CMDTRIES)
        {
          printf ("\n\n**GSSort seems to be stuck or not running\n\n");
          fflush (stdout);
          return (1);
        };

    };

  /* write command */

  fp = fopen ("GSSort.command", "w");
  if (fp == NULL)
    {
      printf ("could not open: \"GSSort.command\"\n");
      return (1);
    };
  fprintf (fp, "%s", cmd);
  fclose (fp);
  printf ("[%s] written to \"GSSort.command\"\n", cmd);

  return (0);

}

/*-----------------------------------------------------------------------------*/

void
write_ascii (char *hist1_name, char *name)
{
  int i, st;
  int nbins;
  char fname[500];
  char buffer[500];
  int file;
  TH1D *hist1;

  hist1 = (TH1D *) gROOT->FindObject (hist1_name);
  if (hist1 == NULL)
    {
      printf ("write_ascii: Error: Could not find 1D histogram %s\n", hist1_name);
      fflush (stdout);
      return;
    }

//  sprintf (fname, "%s.ascii", hist1_name);
  file = open (name, O_WRONLY | O_CREAT | O_TRUNC, 0644);

  nbins = hist1->GetNbinsX ();

  for (i = 1; i <= nbins; i++)
    {
      sprintf (buffer, "%f %f\n", (double) hist1->GetBinLowEdge (i), (double) hist1->GetBinContent (i));
      st = write (file, buffer, strlen (buffer));
    }

  close (file);

}

/*-----------------------------------------------------------------------------*/

void
trace3d (Char_t * fn, Int_t Igam = 0)
{
  Int_t nn;
  int i, j, igam=0;
  static int np = 0;

  Double_t px[10];
  Double_t py[10];
  Double_t pz[10];
  int order[10], tmpi;
  double u1, u2, u3, rr, tmp;
  double pol, azi, xylim;
  double s1, s2, s3, v1, v2, v3;
  double dx, dy, dz;
  FILE *fp;
  char *pc, *p;
  char str[STRLEN], str1[STRLEN];
  float xx, yy, zz;

  /* open file */

  printf ("try to open file: \"%s\"\n", fn);
  fp = fopen (fn, "r");
  if (fp == NULL)
    {
      printf ("\"%s\" not found\n", fn);
      return;
    };

  /* read file */

  nn = 0;
#if(1)
  px[nn] = 0;
  py[nn] = 0;
  pz[nn] = 0;
  order[nn] = 0;
  nn++;
#endif
  pc = fgets (str, STRLEN, fp);
  while (pc != NULL)
    {
      printf ("%s", str);

      if ((p = strstr (str, "valid")) != NULL)
        {
          sscanf (str, "%i: ", &igam);
//         printf("igam=%i, Igam=%i\n", igam,Igam);
        };
      if (str[0] == 35 && igam == Igam)
        {
//        printf("%s", str);
          sscanf (str, "%s %s ( %f %f %f ) %s %i", str1, str1, &xx, &yy, &zz, str1, &order[nn]);
          px[nn] = (double) xx;
          py[nn] = (double) yy;
          pz[nn] = (double) zz;
          order[nn]++;
//          printf ("%f %f %f order:%i\n", px[nn], py[nn], pz[nn], order[nn]);
          nn++;
        }

      pc = fgets (str, STRLEN, fp);
    };
  fclose (fp);

  /* check */

  if (nn == 1)
    {
      printf ("single interaction or error, cannot plot, nn=%i\n", nn);
      return;
    };

  /* order them using bubble sort */

  for (i = 0; i < nn; i++)
    for (j = 0; j < nn; j++)
      if (order[i] < order[j])
        {
          tmp = px[i];
          px[i] = px[j];
          px[j] = tmp;

          tmp = py[i];
          py[i] = py[j];
          py[j] = tmp;

          tmp = pz[i];
          pz[i] = pz[j];
          pz[j] = tmp;

          tmpi = order[i];
          order[i] = order[j];
          order[j] = tmpi;

        };
#if(0)
  for (i = 0; i < nn; i++)
    {
      printf ("%7.2f %7.2f %7.2f order: %i ; ", px[i], py[i], pz[i], order[i]);
      rr = px[i] * px[i] + py[i] * py[i] + pz[i] * pz[i];
      rr = sqrt (rr);
      printf ("radius= %f; ", rr);
      if (i > 0)
        {
          rr =
            (px[i] - px[i - 1]) * (px[i] - px[i - 1]) + (py[i] - py[i - 1]) * (py[i] - py[i - 1]) + (pz[i] -
                                                                                                     pz[i -
                                                                                                        1]) * (pz[i] -
                                                                                                               pz[i -
                                                                                                                  1]);
          printf ("distance: %f\n", sqrt (rr));
        }
      else
        printf ("first interaction point\n");
    }
#endif

  /* find mean of cluster */

  u1 = 0;
  u2 = 0;
  u3 = 0;
  for (i = 0; i < nn; i++)
    {
      u1 += px[i] / nn;
      u2 += py[i] / nn;
      u3 += pz[i] / nn;
    };
//  printf ("mean of cluster %f %f %f\n", u1, u2, u3);

  /* new z axis */

  rr = u1 * u1 + u2 * u2 + u3 * u3;
  rr = sqrt (rr);
//  printf ("u1 radius= %f\n", rr);
  u1 /= rr;
  u2 /= rr;
  u3 /= rr;

  pol = acos (u3);
  azi = acos (u1 / sin (pol));
  if (u2 < 0)
    azi = -azi;
//  printf ("new z axis: pol = %f, azi= %f\n", pol * 57.2958, azi * 57.2958);

  /* new y axis (arbitrary) */

  v1 = -u1;
  v2 = -u2;
  v3 = (u1 * u1 + u2 * u2) / u3;

  rr = v1 * v1 + v2 * v2 + v3 * v3;
  rr = sqrt (rr);
//  printf ("v1 radius= %f\n", rr);
  v1 /= rr;
  v2 /= rr;
  v3 /= rr;

//  printf("check: %f should be zero\n", v1*u1+v2*u2+v3*u3);

  /* new x axis */

  s1 = u2 * v3 - u3 * v2;
  s2 = u3 * v1 - u1 * v3;
  s3 = u1 * v2 - u2 * v1;

  rr = s1 * s1 + s2 * s2 + s3 * s3;
  rr = sqrt (rr);
//  printf ("s1 radius= %f\n", rr);
  s1 /= rr;
  s2 /= rr;
  s3 /= rr;

//  printf("check: %f should be zero\n", v1*s1+v2*s2+v3*s3);
//  printf("check: %f should be zero\n", u1*s1+u2*s2+u3*s3);

  /* find new coordinates */

  printf ("displaying gamma ray # %i\n", Igam);
  for (i = 0; i < nn; i++)
    {
      dx = px[i];
      dy = py[i];
      dz = pz[i];
      px[i] = dx * v1 + dy * v2 + dz * v3;
      py[i] = dx * s1 + dy * s2 + dz * s3;
      pz[i] = dx * u1 + dy * u2 + dz * u3 - 18.5;
      printf ("new: %7.2f %7.2f %7.2f order: %i \n", px[i], py[i], pz[i], order[i]);

    }


  TPolyLine3D *line3D_1 = new TPolyLine3D (nn, px, py, pz);
  line3D_1->SetLineColor (kRed);

  xylim = 14.0;

  Double_t xmin = -xylim, xmax = xylim;
  Double_t ymin = -xylim, ymax = xylim;
  Double_t zmin = .0, zmax = 10.;
  Double_t resolution = 0.01;
  Int_t nBins = Int_t ((xmax - xmin) / resolution);

  /// The binning should be adjusted to the interactive zoom resolution desired

//if (histo != NULL)        histo->Delete ();      

//      histo = (TH1D *) gROOT->FindObject ("histo");
//      printf ("%p; ", );
//      if (histo != NULL)
//       histo->Delete ();

  sprintf (str, "track%i", np);
  np++;
  TH2F *histo = new TH2F (str, "", nBins, xmin, xmax, 10, ymin, ymax);
  histo->SetStats (kFALSE);
  histo->SetMinimum (zmin);
  histo->SetMaximum (zmax);
  histo->SetXTitle ("x[cm]");
  histo->GetXaxis ()->CenterTitle ();
  histo->SetYTitle ("y[cm]");
  histo->GetYaxis ()->CenterTitle ();
  histo->SetZTitle ("z[cm]");
  histo->GetZaxis ()->CenterTitle ();

  /* make canvas and draw lines */

  TCanvas *canvas = new TCanvas ("canvas", "canvas", 0, 0, 1000, 1000);
  canvas->SetTheta (23.);
  canvas->SetPhi (-23.);
  gPad->SetLeftMargin (0.10);
  gPad->SetRightMargin (0.10);
  gPad->SetTopMargin (0.10);
  gPad->SetBottomMargin (0.10);
  histo->Draw ("lego0,fb");
  line3D_1->Draw ();

  /* legend */

  TLegend *leg = new TLegend (0.81, 0.86, 0.99, 0.94);
  leg->SetTextSize (0.06);
  leg->SetBorderSize (0.);
  leg->AddEntry (line3D_1, "Track", "L");
  leg->Draw ();
  //canvas->SaveAs("test.png") ;

  return;
}

/*---------------------------------------------------------------*/

#include "spe_fun.c"

#include "2d_fun.c"
