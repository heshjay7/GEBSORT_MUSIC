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

#include "GEBClient.h"


struct gebClient *GEBClientInit() {
   struct gebClient *retval;
   int testvar = 1;

   retval = (struct gebClient *)calloc(1, sizeof(struct gebClient));
   if (!retval) return 0;

   retval->connectionMutex = epicsMutexCreate();
   if (!retval->connectionMutex) {
      free(retval);
      return 0;
   }

   if (htonl(testvar) == 1) retval->bigendian = 1;

   return retval;
}


/* geb_addr is a network name, like node0005-cluster */
int setGEBClient(struct gebClient *i, char * geb_addr, int port) {

   struct sockaddr_in adr_srvr;  /* AF_INET */
   struct hostent *hp;
   int numrec, reject;

   epicsMutexLock(i->connectionMutex);

   if (i->outSock != 0) {
      close(i->outSock);
   }
   memset(&adr_srvr,0,sizeof adr_srvr);
   adr_srvr.sin_family = AF_INET;
   adr_srvr.sin_port = htons(port);

   hp = gethostbyname(geb_addr);
   if (!hp) {
      printf("hostname %s not resolved\n", geb_addr);
      epicsMutexUnlock(i->connectionMutex);
      return 4;
   }
   adr_srvr.sin_addr.s_addr = inet_addr(inet_ntoa(*((struct in_addr *) 
                                (hp->h_addr_list[0]))));

   if ( adr_srvr.sin_addr.s_addr == INADDR_NONE ) {
      printf("bad geb address %s port %d\n", geb_addr, port);
      epicsMutexUnlock(i->connectionMutex);
      return 1;
   }

   i->outSock = socket(AF_INET, SOCK_STREAM, 0);
   if (i->outSock == -1) {
      printf("Unable to open output socket.\n");
      epicsMutexUnlock(i->connectionMutex);
      return 2;
   }

   if (connect(i->outSock, 
        (struct sockaddr *)&adr_srvr, sizeof(adr_srvr)) < 0) {
      printf("Connect to GEB failed\n");
      close(i->outSock);
      i->outSock = 0;
      epicsMutexUnlock(i->connectionMutex);
      return 3;
   }

   numrec = recv(i->outSock, &reject, sizeof(int),0);
   if (numrec != 4)  {
      printf("initial read returned %d\n", numrec);
      close(i->outSock);
      i->outSock = 0;
      epicsMutexUnlock(i->connectionMutex);
      return 5;
   }

   if (reject) {
      printf("GEB not receiving at this time\n");
      close(i->outSock);
      i->outSock = 0;
      epicsMutexUnlock(i->connectionMutex);
      return 6;
   }

   epicsMutexUnlock(i->connectionMutex);

   return 0;

}

void closeGEBClient(struct gebClient *i) {

   if (i->outSock) {
      epicsMutexLock(i->connectionMutex);
      close(i->outSock);
      i->outSock=0;
      epicsMutexUnlock(i->connectionMutex);
   }
}

/* CheckGEBClient returns basic information about the state of the GEB 
   connection.

   return value of 0 means all is well; 1 means the connection is fine but 
   the GEB is not accepting data because of memory depletion; 2 means that 
   the connection is down and needs to be restarted with setGEBClient(). A
   return value of 3 indicates that GEBClientInit() must be used. 
 */
int checkGEBClient(struct gebClient *i) {

   if (!i || !i->connectionMutex) return 3;
   return sendGEBData(i, 0);
}

/* sendGEBData(i, outmsg) sends data in the outmsg to client i. It byteswaps
   the outmsg elements to littleendian but not the payload before sending.  It 
   sends the auxData payload as a string of bytes.  Freeing of outmsg and the
   contained payload are the responsibility of the caller.

   A return value of soft error indicates that the GEB data buffers are full 
   and of hard error means the GEB has disconnected.  A null value for outmsg
   sends a keepalive message and has the same return value meanings.

   return values: 2 = hard error (need to reconnect GEB, and resend data from 
                      call)
                  1 = soft error (need to keep resending same data until a
                      successful return)
                  0 = successful return
 */
int sendGEBData(struct gebClient *i, struct GEBData *outmsg) {

   unsigned int *outdata;
   int numrec, reject, outsize, sendsize;

   epicsMutexLock(i->connectionMutex);
   if (!i->outSock) {
      epicsMutexUnlock(i->connectionMutex);
      return 2;
   }
   numrec = recv(i->outSock, &reject, sizeof(int),0);
   if (numrec != 4)  {
      printf("ack from GEB failed %s\n", strerror(errno));
      close(i->outSock);
      i->outSock = 0;
      epicsMutexUnlock(i->connectionMutex);
      return 2;
   }
   /* If we get reject, it means that the buffers were full in the GEB as of
      the last send, which was graciously accepted anyway. So this time, just 
      send a keepalive to keep the protocol cranking and return soft error.
    */
   if (reject) reject = 1;
   if (0 == outmsg || reject) {
     outmsg = &i->zeromsg;
     outdata = (unsigned int *)calloc(GEB_HEADER_BYTES/4, 4);
     if (!(i->bigendian)) {
        bcopy(outmsg, &outdata[0], GEB_HEADER_BYTES);
     } else {
        swab(outmsg, &outdata[0], GEB_HEADER_BYTES);
     }
     sendsize = GEB_HEADER_BYTES;
   } else {
     if (!outmsg->payload || outmsg->length < 1) {
       printf("Garbage into GEB Client. dat 0x%p len %d \n", outmsg->payload,
                                                              outmsg->length);
       close(i->outSock);
       i->outSock = 0;
       epicsMutexUnlock(i->connectionMutex);
       return 2;
     }
     sendsize = GEB_HEADER_BYTES + outmsg->length;
     outsize = sendsize/4 + 1;
     outdata = (unsigned int *)calloc(outsize, 4);

     if (!(i->bigendian)) {
        bcopy(outmsg, &outdata[0], GEB_HEADER_BYTES);
        bcopy(outmsg->payload, &outdata[4], outmsg->length);
     } else {
        swab(outmsg, &outdata[0], GEB_HEADER_BYTES);
        swab(outmsg->payload, &outdata[4], outmsg->length);
     }
   }

  if ( write(i->outSock, outdata, sendsize) < 0) {
     printf("output to GEB failed %s\n", strerror(errno));
     close(i->outSock);
     i->outSock = 0;
     printf("GEB disconnected\n");
     epicsMutexUnlock(i->connectionMutex);
     free(outdata);
     return 2;
   }
   epicsMutexUnlock(i->connectionMutex);
   free(outdata);
   return reject;
}


