#include <epicsMutex.h>
#include "GEBLink.h"

struct gebClient {
   epicsMutexId connectionMutex;
   int outSock;
   int bigendian;
   struct GEBData zeromsg;
};

struct gebClient *GEBClientInit();

int sendGEBData(struct gebClient *, struct GEBData *);
int setGEBClient(struct gebClient *, char *addr, int port);
void closeGEBClient(struct gebClient *);
int checkGEBClient(struct gebClient * );

