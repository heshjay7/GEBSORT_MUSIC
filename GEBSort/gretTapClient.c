
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <netdb.h>
#include <errno.h>

#include <epicsMutex.h>

#include "GEBLink.h"
#include "gretTapClient.h"

enum gretTapClientErrorValues gretTapClientError;

char *gretTapClientErrorStrings[] = { "GTC_NOERROR", "GTC_WRITE", "GTC_READ",
     "GTC_M_CREATE", "GTC_HOSTNAME", "GTC_INADDR", "GTC_INSOCK", "GTC_TAPCONN",
     "GTC_INITWRITE", "GTC_INITREAD", "GTC_NOT_FOUND", "GTC_NOT_RUNNING",
     "GTC_UNKNOWN", "GTC_CLOSED", "GTC_HEADER_READ", "GTC_HEADER_WRITE", 
     "GTC_TIMEOUT",
     "GTC_GEB_HEADER", "GTC_LEN_0", "GTC_LEN_HUGE", "GTC_DATA_READ" };


#define TAP_HEADER_LEN (2 * sizeof(int))

/* Do a write with proper handling of partial writes 
 * return byteswrit (0 or -1) on failure, amount written on success 
 */
int fdWrite(int where, void *dat, int len) {
   int towrite, numwrit = 0, byteswrit = 0;
   unsigned char *outp;

   towrite = len;
   outp = (unsigned char *)dat;
   while (numwrit < towrite) {
      byteswrit = write(where, outp + numwrit, len - numwrit);
      if (byteswrit <= 0) {
         return byteswrit;
      }
      numwrit += byteswrit;
   }
   return numwrit;
}

/* Do a read with proper handling of partial reads
 * return bytesread (0 or -1) on failure, amount read on success 
 */
int fdRead(int where, void *dat, int len) {
   int toread, numread = 0, bytesread = 0;
   unsigned char *inp;

   toread = len;
   inp = (unsigned char *)dat;
   while (numread < toread) {
      bytesread = read(where, inp + numread, len - numread);
      if (bytesread <= 0) {
         return bytesread;
      }
      numread += bytesread;
   }
   return numread;
}

struct gretTap *gretTapConnect(char *addr, int position, unsigned int typeMask)
{
   struct sockaddr_in adr_srvr;  /* AF_INET */
   struct hostent *hp;
   struct gretTap *retval;
   int numrec, firstret[2];
   unsigned int outbuf[2];
   int numsent;

   gretTapClientError = GTC_NOERROR;
   retval = (struct gretTap *) calloc(1, sizeof(struct gretTap));
   retval->connectionMutex = epicsMutexCreate();
   if (!retval->connectionMutex) {
      gretTapClientError = GTC_M_CREATE;
      free(retval);
      return 0;
   }
   memset(&adr_srvr,0,sizeof adr_srvr);
   adr_srvr.sin_family = AF_INET;
   adr_srvr.sin_port = htons(GRETTAP_PORT);

   hp = gethostbyname(addr);
   if (!hp) {
      gretTapClientError = GTC_HOSTNAME;
      free(retval);
      return 0;
   }
   adr_srvr.sin_addr.s_addr = inet_addr(inet_ntoa(*((struct in_addr *)
                                (hp->h_addr_list[0]))));

   if ( adr_srvr.sin_addr.s_addr == INADDR_NONE ) {
      gretTapClientError = GTC_INADDR;
      free(retval);
      return 0;
   }

   retval->inSock = socket(AF_INET, SOCK_STREAM, 0);
   if (retval->inSock == -1) {
      gretTapClientError = GTC_INSOCK;
      free(retval);
      return 0;
   }
   if (connect(retval->inSock,
        (struct sockaddr *)&adr_srvr, sizeof(adr_srvr)) < 0) {
      gretTapClientError = GTC_TAPCONN;
      close(retval->inSock);
      free(retval);
      return 0;
   }
   /* send two ints: type and position */
   outbuf[0] = typeMask;
   outbuf[1] = position;
   
   numsent = fdWrite(retval->inSock, outbuf, 2 * sizeof(unsigned int));
   if (numsent != 2 * sizeof(unsigned int)){
      gretTapClientError = GTC_INITWRITE;
      close(retval->inSock);
      free(retval);
      return 0;
   }
   numrec = fdRead(retval->inSock, firstret, TAP_HEADER_LEN);
   if (numrec != TAP_HEADER_LEN)  {
      gretTapClientError = GTC_INITREAD;
      close(retval->inSock);
      free(retval);
      return 0;
   }
   
   if (firstret[0] != TAP_ACK) {
      switch (firstret[0]) {
         case TAP_NOT_FOUND:
            gretTapClientError = GTC_NOT_FOUND;
            break;
         case TAP_NOT_RUNNING:
            gretTapClientError = GTC_NOT_RUNNING;
            break;
         default:
            gretTapClientError = GTC_UNKNOWN;
      }
      close(retval->inSock);
      free(retval);
      return 0;
   }

   return (retval);
}

void gretTapDataFree(struct GEBData *in) {
   struct GEBData *temp;

   while ( in ) {
      temp = in->next;
      free(in->payload);
      free(in);
      in = temp;
   }
}
     
struct GEBData *gretTapData(struct gretTap *tap, int nreq, float timeout) {
    int numrec, i;
    int outbuf[2];
    int firstret[2] = {0,0};
    struct GEBData *retval=0, *inbuf = 0;

    gretTapClientError = GTC_NOERROR;

    if (tap->inSock == 0) {
       gretTapClientError = GTC_CLOSED;
       return 0;
    }

    outbuf[0] = nreq;

    if (timeout < 0) { 
       outbuf[1] = 0;
    } else {
       outbuf[1] = (int)(timeout * 1000.0);
    }

   if (fdWrite(tap->inSock, outbuf, TAP_HEADER_LEN) != TAP_HEADER_LEN) {
      gretTapClientError = GTC_HEADER_WRITE;
      close(tap->inSock);
      tap->inSock = 0;
      return 0;
   }
    
   numrec = fdRead(tap->inSock, firstret, TAP_HEADER_LEN);
   if (numrec != TAP_HEADER_LEN)  {
      if (numrec) {
         gretTapClientError = GTC_HEADER_READ;
      } else {
         gretTapClientError = GTC_CLOSED;
      }
      close(tap->inSock);
      tap->inSock = 0;
      return 0;
   }
   if (firstret[0] == TAP_TIMEOUT) {
      gretTapClientError = GTC_TIMEOUT;
   }
   if (firstret[0] == TAP_NOT_RUNNING) {
      gretTapClientError = GTC_NOT_RUNNING;
   }
   if (firstret[1] == 0) {
      return 0;
   }


   for (i = 0; i < firstret[1]; i++) {
      if (!retval) {
         inbuf = (struct GEBData *) calloc(1, sizeof(struct GEBData));
         retval = inbuf;
      } else {
         inbuf->next = (struct GEBData *) calloc(1, sizeof(struct GEBData));
         inbuf = inbuf->next;
      }
      numrec = fdRead(tap->inSock, inbuf, GEB_HEADER_BYTES);
      if (numrec != GEB_HEADER_BYTES) {
          gretTapClientError = GTC_GEB_HEADER;
          gretTapDataFree(retval);
          close(tap->inSock);
          tap->inSock = 0;
      }
#if(1)
      if (inbuf == NULL) {
          printf("inbuf=0x%p, return GTC_LEN_0\n",inbuf);
          fflush(stdout);
          gretTapClientError = GTC_LEN_0;
          inbuf->next = (struct GEBData *) calloc(1, sizeof(struct GEBData));
          inbuf = inbuf->next;
          continue;
      }
#endif
      if (inbuf->length == 0) {
          gretTapClientError = GTC_LEN_0;
          inbuf->next = (struct GEBData *) calloc(1, sizeof(struct GEBData));
          inbuf = inbuf->next;
          continue;
      }
      inbuf->payload = calloc(1, inbuf->length);
      if (!inbuf->payload) {
          /* if length rediculous this could happen ? */
          gretTapClientError = GTC_LEN_HUGE;
          inbuf->next = (struct GEBData *) calloc(1, sizeof(struct GEBData));
          inbuf = inbuf->next;
          continue;
      }
      numrec = fdRead(tap->inSock, inbuf->payload, inbuf->length);
      if (numrec != inbuf->length) {
          gretTapClientError = GTC_DATA_READ;
          gretTapDataFree(retval);
          close(tap->inSock);
          tap->inSock = 0;
      }
   }
   return retval;
}

void gretTapClose(struct gretTap *tap) {
   
    if (!tap) return;

    if (!tap->inSock) {
       free(tap);
       return;
    }

    close(tap->inSock);
    free(tap);
}
       
int gretTapCheck(struct gretTap *tap) {
   if (!tap->inSock) {
      return 1;
   }
   return 0;
}
   
