
#define PMODE 0644
#define LENSP 16384

#define USEBREAD 0
#define BREAD_BUFSIZE 16384
#define USEZLIB 0

#define MAXFILES 1000
#define MAXCOINEV 1000
#define DEBUG1 0
#define NCHANNELS 1000*11
#define MAXVALIDHIE 32000
#define NGE 110
#define MAXBOARDID 100000

#define TRUE 1
#define FALSE 0
#define NOTDEF -1

#define EOE          0xaaaaaaaa
#define HDRLENBYTES  28
//#define HDRLENINTS    7
#define HDRLENINTS  13  
#define HDRLENWORDS  2*HDRLENINTS
#define MAXLENINTS  519

#define LENEOVWORDS  2

#define MAXTRACELEN 8192
#define EOE          0xaaaaaaaa

#define CC_ID1 10      /* instead of 29 as default */
#define CC_ID2 11      /* instead of 39 as default */

/* instrument types */

#define NOTHING 0
#define GE 1
#define BGO 2
#define SIDE 3
#define AUX 4
#define DSSD 5
#define FP 6
#define XARRAY 7
#define CHICO2 8
#define SSD 9
#define CLOVER 10
#define SPARE 11
#define SIBOX 12

#define MAXTPE 12
#define MAXTID 1000

/*---------------*/
/* single events */
/*---------------*/


typedef struct DGSEVENT
{
  float                   ehi;
  float                   ehi_nodop;
  short int               id;
  //short int               id_pw[MAXPW];
  unsigned short int      tpe, tid;
  unsigned short int      flag;
  unsigned short int      board_id;
  unsigned short int      chan_id;
  unsigned short int      geo_addr;
  unsigned short int      packet_length;

  unsigned short int      header_type;
  unsigned short int      event_type;
  unsigned short int      header_length;

  unsigned long long int  event_timestamp;
  unsigned long long int  last_disc_timestamp;
  unsigned long long int  peak_timestamp;

  unsigned short int      timestamp_match_flag;
  unsigned short int      external_disc_flag;
  unsigned short int      cfd_valid_flag;
  unsigned short int      pileup_only_flag;
  unsigned short int      offset_flag;
  unsigned short int      sync_error_flag;
  unsigned short int      general_error_flag;

  unsigned short int      peak_valid_flag;
  unsigned short int      pileup_flag ;

  int                     sampled_baseline;
  int                     cfd_sample_0;
  int                     cfd_sample_1;
  int                     cfd_sample_2;
  int                     sum1;
  int                     sum2;
  
  unsigned short int      m2end;
  float                   cfd_interpolate; // should be added to TS for additional accuracy - in 10ns units.
  unsigned short int      m2begin;
  unsigned short int      m2last_begin_sample;
  unsigned short int      m2last_end_sample;
  unsigned short int      m1begin;
  unsigned short int      m1end;
  unsigned short int      peak_sample;
  unsigned short int      base_sample;
  
  int                     baseline;
  
  unsigned short int      traceLen;
  short int               trace[MAXTRACELEN];  


  //unsigned short int      A[MAXPW];
  //unsigned short int      B[MAXPW];
  //unsigned short int      C[MAXPW];
  //unsigned short int      T[MAXPW];


  unsigned long long int  LEDts;
  unsigned long long int  CFDts;
  //unsigned long long int  PEAKts;
  //char                    flag;  
  //short int               baseline;
  //unsigned short int      traceLen;
  //short int               trace[MAXTRACELEN];
}  DGSEVENT;

#define DGSEVENT_BASELEN sizeof(DGSEVENT)-MAXTRACELEN*sizeof(short int)

typedef struct DFMAEVENT
{
  int               ehi;    // WAS SHORT INT
  short int               id;
  unsigned short int      tpe, tid;
  unsigned short int      board_id;
  unsigned short int      chan_id;
  unsigned long long int  LEDts;
  unsigned long long int  CFDts;
  unsigned long long int  PEAKts;
  char                    flag;  
  char 			  pu;
  int			  d2t0;
  int			  d2t1;
  int                     d2t2;
  int wheel;
  unsigned long long int  prevTS;
  int               baseline;
  int               postrisebeg;
  int               prerisebeg;
  int               postriseend;
  int               preriseend;
  int               peaksample;
  int               basesample;
  int		    postrisesum;
  int		    prerisesum;
  int               header_type;
  int               m2_last_beg;
  int               m2_last_end;
  unsigned short int      traceLen;
  short int               trace[MAXTRACELEN];
}  DFMAEVENT;

#define DFMAEVENT_BASELEN sizeof(DFMAEVENT)-MAXTRACELEN*sizeof(short int)

typedef struct GTEVENT
{
  /* raw */

  unsigned short int      len;
  short int               ehi;
  short int               id;
  short int               module;
  unsigned short int      tpe, tid;
  unsigned short int      board_id;
  unsigned short int      chan_id;
  unsigned long long int  LEDts;
  unsigned long long int  CFDts;
  unsigned long long int  PEAKts;
  char                    flag;  
  short int               baseline;
  int                     rawE;
  unsigned int            hdr[HDRLENINTS];
  unsigned short int      traceLen;
  short int               trace[MAXTRACELEN];
}  GTEVENT;


/*--------------------*/
/* coincidence events */
/*--------------------*/

typedef struct COINEV_struct
  {
  unsigned short int      len;
  unsigned char      lenclean;
  unsigned char      lendirty;
  unsigned char      lenbgo;
  unsigned char      lenaux;
  GTEVENT  GTEvent[MAXCOINEV];
  DGSEVENT DGSEvent[MAXCOINEV];
  DFMAEVENT DFMAEvent[MAXCOINEV];
  } COINEV;


typedef struct CONTROL_struct
  {
  int nOpenFiles;
  int fileActive[MAXFILES];
  int fileEventsRead[MAXFILES];
  int filesiz[MAXFILES];
  int nwritten;
  int nwritten_chunk;
  int chunkno;
  long long int chunksiz;
  int nread;
  int nzeroehi;
  int minGE;
  int minFP;
  int minDSSD;
  int minmult;
  int TSlistelen;
  int noverflowehi;
  int CurEvNo;
  long long int dts_min;
  long long int dts_max;
  int suppressBadAtRead;
  int zzipout;
  int dtsfabort;
  int dtsbabort;
  long long int startTS_lo;
  long long int startTS_hi;
  int startTS;
  int TSlist_lo;
  int TSlist_hi;
  unsigned int waitfordataseconds;
  } CONTROL;


typedef struct stat_struct
  {
  long long int inbytes;
  long long int badid;
  long long int outbytes;
  unsigned int ge_hit[NGE+1];
  unsigned int ge_cs[NGE+1];
  int nTSjumprecover_f[MAXTID];
  int nTSjumprecover_b[MAXTID];
  unsigned int in_hit[MAXTPE][MAXTID];
  unsigned int out_hit[MAXTPE][MAXTID];
  int nbigbufreads;
  int nswaps;
  unsigned int id_hit[NCHANNELS];
  unsigned int GEBIds[30];
  long long int GEBlen[30];
  } MSTAT;


/*-----------------*/
/* header for file */
/*-----------------*/

typedef struct DGSHEADER_struct
  {
  unsigned int      id;
  char              date[19];
  unsigned int      nfiles;
  char              dummy;
  } DGSHEADER;





