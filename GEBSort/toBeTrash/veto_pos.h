/* position veto */

/* veto cube parameters */

#define VETO_X_D (2*42.0)
#define VETO_Y_D (2*42.0)
#define VETO_Z_MAX 92.0
#define VETO_BINWIDTH 1.0

/* grid size evaluations */

#define VETO_NX (int)(VETO_X_D/VETO_BINWIDTH+1)
#define VETO_NY VETO_NX
#define VETO_NZ (int)(VETO_Z_MAX/VETO_BINWIDTH+1)

//unsigned int veto_cube[MAXGTMODNO][MAXCRYSTALNO][VETO_NX][VETO_NY][VETO_NZ];
//calloc in program instead

/* index functions */

#define VETO_X_INDEX(x)  (x*VETO_NX/VETO_X_D+VETO_NX/2.0)
#define VETO_Y_INDEX(x)  (x*VETO_NY/VETO_X_D+VETO_NY/2.0)
#define VETO_Z_INDEX(x)  (x*VETO_NZ/VETO_Z_MAX)
 
/* reverse index functions */

#define VETO_X_XEDNI(i)   ( ((float)i+0.5) -VETO_NX/2.0) *VETO_X_D/VETO_NX 
#define VETO_Y_XEDNI(i)   ( ((float)i+0.5) -VETO_NY/2.0) *VETO_Y_D/VETO_NY 
#define VETO_Z_XEDNI(i)     ((float)i+0.5) *VETO_Z_MAX/VETO_NZ

/* careful with this index function... since */
/* we go from i0, <=i, we need the +1 everywhere */

#define VETO_INDX(M,C,X,Y,Z) \
                             (M  *(MAXCRYSTALNO+1)*(VETO_NX+1)*(VETO_NY+1)*(VETO_NZ+1) \
                             +C                   *(VETO_NX+1)*(VETO_NY+1)*(VETO_NZ+1) \
                             +X                               *(VETO_NY+1)*(VETO_NZ+1) \
                             +Y                                           *(VETO_NZ+1) \
                             +Z)


