/*
 * Positions. Types are in GEBLink.h */
#define GRETTAP_RAW		1
#define GRETTAP_CRYSTAL_EVENT	2
#define GRETTAP_POSITION	3
#define GRETTAP_GEB		4
#define GRETTAP_TRACK		5

/* Tap network status (in reply word 1) in GEBLink.h */

/* 
 * Error values
 */
enum gretTapClientErrorValues { GTC_NOERROR, GTC_WRITE, GTC_READ, GTC_M_CREATE,
     GTC_HOSTNAME, GTC_INADDR, GTC_INSOCK, GTC_TAPCONN, GTC_INITWRITE,
     GTC_INITREAD, GTC_NOT_FOUND, GTC_NOT_RUNNING, GTC_UNKNOWN, GTC_CLOSED,
     GTC_HEADER_READ, GTC_HEADER_WRITE, GTC_TIMEOUT, GTC_GEB_HEADER, GTC_LEN_0, GTC_LEN_HUGE,
     GTC_DATA_READ };
extern char *gretTapClientErrorStrings[];

extern enum gretTapClientErrorValues gretTapClientError;

#define GRETTAP_PORT	9305
#include <epicsMutex.h>

struct gretTap {
   int inSock;
   epicsMutexId connectionMutex;
}; 
/* 
 * Connect to a tap.  
 * addr is the name of the node to connect to. 
 * position indicates which of the available taps are 
 * to be used, from the list above.  GRETTAP_RAW can be found in digitizer 
 * crates, CRYSTAL_EVENT and POSITION are in signal decomposition nodes, 
 * GEB in the global event builder, and TRACK in tracking nodes.
 * gebType indicates what type of data is desired; this is only necessary in 
 * the case of GRETTAP_GEB where various auxiliary detector data may be 
 * available.
 *
 * The struct returned is used only for calling the other functions below and
 * is not for user manipulation.
 */
struct gretTap *gretTapConnect(char *addr, int position, unsigned int typeMask);
/* 
 * get tap data
 * Returns a gebData struct containing some data.  The gebData struct itself
 * gives a type, timestamp and length of the attached buffer. Note that a 
 * type of 0 returns all
 * 
 * nreq specifies amount of data to get in events.  This call actually goes 
 * and gets the info.
 * 
 */
struct GEBData *gretTapData(struct gretTap *, int nreq, float timeout);
/* 
 * close tap
 */
void gretTapClose(struct gretTap *);
/* 
 * check tap.  return of 0 indicates things are fine, 1 that there is some
 * fatal problem.  In that case, gretTapClose() should be called and a new
 * connection established.
 */
int gretTapCheck(struct gretTap *);

/*
 * a utility function to free a list of the sort returned by gretTapData 
 */
void gretTapDataFree(struct GEBData *);
