
#ifndef _GEBSort_h
#define _GEBSort_h


#include "TFile.h"
#include "TTree.h"
#include "gdecomp.h"

#define TAPE 0
#define NET  1
#define DISK 2
#define GEB  3

#define MEDIUMLEN 2048
#define LONGLEN 8000
#define SHORTLEN 1024
#define RBUFSIZE 500000
#define MAXPAYLOADSIZE 500000
#define STRLEN 256

//#define MAXDETPOS 180
//#define MAXCRYSTALNO 3
//#define MAXDETNO 200
//#define MAXMODNO 60
//#define NUMAGATAPOS 180

/* veto cube */
//#define MAXGTMODNO 30
//
//#define RATELEN 60*24
//#define DTBTEVLEN 1000
//#define MAXNOSEG 9
//#define NUMSEGS 36
//#define MAXSEGNO MAXDETNO*NUMSEGS
//#define RMIN 10
//#define RMAX 35
//
#define MAXGEBS 1000
#define MAXLONG (long long)1<<62

//
//#define MAX_GAMMA_RAYS 1000
//#define GEB_BITS GEB_HEADER_BYTES*8

/* max values for # of bits */

#define M14BITS 0x3fff 
#define M13BITS 0x1fff
#define M12BITS 0x0fff
#define M11BITS 0x07ff
#define M10BITS 0x03ff

/* basic spectra lengths */

#define L14BITS  M14BITS+1
#define L13BITS  M13BITS+1
#define L12BITS  M12BITS+1
#define L11BITS  M11BITS+1
#define L10BITS  M10BITS+1

typedef struct EXCHANGE
{

  int ngates;

} EXCHANGE;


typedef struct GEB_event
{
  int mult;
  CRYS_INTPTS *ptinp[MAXGEBS];
  GEBDATA *ptgd[MAXGEBS];
} GEB_EVENT;

typedef struct PARS
{
  char ROOTFile[STRLEN];
  int nEvents;
  char GTSortInputFile[STRLEN];
  unsigned int StartMapAddress;
  int InputSrc;
  int WeWereSignalled;
//  char spname[STRLEN];
  int firstEvent;
  int GSudpPort;
  int NumToPrint;
  int DumpEvery;
  TFile *f1;
  TList *wlist;
  long long int curTS;
  long long int maxTS, tmpmaxTS;
  long long int dTS;
  long long int nbytes;
  int CurEvNo;
  char pHost[16];
  int grouping;
  int type;
  float timeout;
  int modwrite;
  int tsnumwrites;
  int do_bin_XA;
  int echo_data;
  char echo_data_fn[512];
  off_t echo_data_pipe;

} PARS;



/* macros */

/*
#define WRITEALLHISTS  \
  gDirectory->cd("root:/"); \
  wlist = gDirectory->GetList(); \
  if (ComPressLevel>NOTDEF) f1->SetCompressionLevel(ComPressLevel); \
  printf("writing all spectra to [%s]\n", Pars.RootFile); \
  printf("be patient... "); \
  fflush(stdout); \
  t1 = time(NULL); \
  wlist->Write(0,TObject::kOverwrite); \
  t2 = time(NULL); \
  printf("DONE! on "); \
  time_stamp(stderr); \
  printf("file size: %i, ",f1->GetSize()); \
  printf("compression level: %i, ",f1->GetCompressionLevel()); \
  printf("and factor: %f\n",f1->GetCompressionFactor()); \
  printf("uncompressed root file size: %f\n",f1->GetSize()*f1->GetCompressionFactor()); \
  printf("writeout time: %i seconds\n", t2 - t1); \
  printf("at %7.2f Mbytes/sec\n", (float) f1->GetSize() / (t2 - t1) / 1000000.0); \
  printf("on "); \
  time_stamp(stderr); \
  fflush(stdout);

#define UPDSSHMEM \
  t1 = time(NULL); \
  mfile->Update(); \
  t2 = time(NULL); \
  printf("done! "); \
  printf("shared memory size: %i\n", mfile->GetSize()); \
  printf("update time: %i seconds ", t2 - t1); \
  printf("at %7.2f Mbytes/sec\n", (float) mfile->GetSize() / (t2 - t1) / 1000000.0); \
  printf("to mapfile [%s] on ",Pars.ShareMemFile); \
  time_stamp(stderr); \
  fflush(stdout);
*/

#endif	/* _GEBSort_h */
