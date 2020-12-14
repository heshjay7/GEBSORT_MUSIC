
#define PMODE 0644
#define LENSP 16384

#define USEBREAD 0
#define BREAD_BUFSIZE 16384
#define USEZLIB 0

#define MAXFILES 600
#define MAXCOINEV 500
#define DEBUG1 0
#define NCHANNELS 1000*11
#define MAXVALIDHIE 32000
#define NCRYSTALS 120+1
#define MAXBOARDID 100000
#define NSEG 31*4*40

#define TRUE 1
#define FALSE 0
#define NOTDEF -1

/* GRETINA has 7 32 bit ints header */
/* HDRLENWORDS counts 16 bit words */

#define HDRLENINTS    7
#define HDRLENBYTES  4*HDRLENINTS
#define HDRLENWORDS  2*HDRLENINTS

#define EOE          0xaaaaaaaa
#define MAXLENINTS  519

#define LENEOVWORDS  2

#define MAXTRACELEN 8192
#define EOE          0xaaaaaaaa

#define CC_ID1 10      /* instead of 29 as default */
#define CC_ID2 11      /* instead of 39 as default */

#define Mvalue 32

#define ELENGTH 4096

typedef struct GTEVENT
{
  /* raw */

  unsigned short int      len;
  short int               ehi;
  short int               ge_id;
  short int               module_id;
  unsigned short int      tpe, tid;
  unsigned short int      digitizer_id;
  unsigned short int      crystal_id;
  unsigned short int      chan_id;
  unsigned short int      seg_id;
  short int               channel;
  unsigned long long int  LEDts;
  unsigned long long int  CFDts;
  unsigned long long int  PEAKts;
  char                    flag;  
  short int               baseline;
  int                     rawE;
  unsigned short int      hdr[HDRLENINTS];
  unsigned short int      traceLen;
  short int               trace[MAXTRACELEN];
}  GTEVENT;







