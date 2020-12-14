struct GEBData {
int type;
int length;     /* of payload in bytes */
long long timestamp;
void *payload;
short refCount;		/* for data tap upgrade, all programs */
short refIndex;
struct GEBData *next;	/* used in TrackIF.c */
};
#define GEB_TYPE_KEEPALIVE 	0
#define GEB_TYPE_DECOMP 	1
#define GEB_TYPE_RAW		2
#define GEB_TYPE_TRACK		3
#define GEB_TYPE_BGS		4
#define GEB_TYPE_S800		5

#define GEB_PORT 9005
#define GEB_HEADER_BYTES 16

/*
 * These are a few things used in gretTap
 */

#define TAP_DATA                0
#define TAP_ACK                 1
#define TAP_TIMEOUT             2
#define TAP_NOT_FOUND           4
#define TAP_NOT_RUNNING         8
#define TAP_ERROR		16

