#include <stdlib.h>
#include <fcntl.h>
#include <stddef.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <limits.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/param.h>
#include <sys/resource.h>
#include <sys/times.h>
#include <assert.h>
#include <limits.h>

/* collection of functions for vector manipulation */

#define MAXERR 0.0001

#if(0)
  /* prototypes [for copy and paste] */

  double unitvector (double *, double *, double *);
  double dotproductangle (double, double, double, double, double, double);
  int crossprod (double, double, double, double, double, double, double *, double *, double *);
  int crossprod1 (double, double, double, double, double, double, double *, double *, double *);
  int crossprod2 (double, double, double, double, double, double, double *, double *, double *);
  double rad2deg (double);
  double find_azimuth (double, double);
  double AzimuthFromCartesian (double, double, double);
  double PolarFromCartesian (double, double, double, double *);
  double vectorlen (double, double, double);
  int coord_in_new (double, double, double,
                    double, double, double,
                    double, double, double, double, double, double, double *, double *, double *);
  int check_coord_sys (double, double, double, double, double, double, double, double, double);
  int check_unitvector (double, double, double);
  int ranVectorOnSphere (double *, double *, double *);
#endif

/*-----------------------------------------------------------*/

double
PolarFromCartesian (double xx, double yy, double zz, double *rr)
{
  double d1;

  *rr = sqrt (xx * xx + yy * yy + zz * zz);
  d1 = acos (zz / *rr);

  return (d1);
}

/*-----------------------------------------------------------*/

double
AzimuthFromCartesian (double xx, double yy, double zz)
{

  double d1;

#if(0)
  if (xx > 0 && yy >= 0)
    d1 = atan (yy / xx);
  else if (xx > 0 && yy < 0)
    d1 = atan (yy / xx) + 2. * M_PI;
  else if (xx < 0)
    d1 = atan (yy / xx) + M_PI;
  else if (xx == 0 && yy > 0)
    d1 = M_PI / 2.;
  else if (xx == 0 && yy < 0)
    d1 = 3. * M_PI / 2.;
  else
    d1 = 0.0;
#endif
#if(1)
  d1 = atan2 (yy, xx);
//  if (d1<0) d1+=M_PI;
//  if (d1>M_PI) d1-=M_PI;
#endif


  return (d1);

}

/*-----------------------------------------------------*/

int
ranInsideCirle (double *u0, double *u1)
{
  /* declarations */

  double rr;

  rr = 2;
  while (rr > 1)
    {
      *u0 = 2 * (drand48 () - 0.5);
      *u1 = 2 * (drand48 () - 0.5);
      rr = (*u0) * (*u0) + (*u1) * (*u1);
    };

  /* done */

  return (0);


};

/*-----------------------------------------------------*/

int
ranVectorOnSphere (double *x0, double *x1, double *x2)
{

  /* see section 21.5.1 in Numerical Recipes book, page 1129 */
  /* generating vector on sphere is non trivial */

  /* declarations */

  double u0, u1;

  /* get a point within a cirle */

  ranInsideCirle (&u0, &u1);

  /* generate random points on sphere */

  *x0 = 2 * u0 * sqrt (1 - u0 * u0 - u1 * u1);
  *x1 = 2 * u1 * sqrt (1 - u0 * u0 - u1 * u1);
  *x2 = 1 - 2 * (u0 * u0 + u1 * u1);

  /* done */

  return (0);

};

/*-----------------------------------------------------*/

int
crossprod (double u1, double u2, double u3, double v1, double v2, double v3, double *s1, double *s2, double *s3)
{


//  printf ("vector  1: (%10.7f %10.7f %10.7f)\n", u1,u2,u3);
//  printf ("vector  2: (%10.7f %10.7f %10.7f)\n", v1,v2,v3);

  *s1 = u2 * v3 - u3 * v2;
  *s2 = u3 * v1 - u1 * v3;
  *s3 = u1 * v2 - u2 * v1;
//  printf ("1 cross 2: (%10.7f %10.7f %10.7f)\n",*s1, *s2, *s3);

  return (0);

};

/*-----------------------------------------------------*/

int
crossprod1 (double a1, double b1, double c1, double a2, double b2, double c2, double *s1, double *s2, double *s3)
{

  /* CRC 22 ed, page 556 */

  *s1 = b1 * c2 - b2 * c1;
  *s2 = c1 * a2 - c2 * a1;
  *s3 = a1 * b2 - a2 * b1;

  return (0);

}

/*-----------------------------------------------------*/

int
crossprod2 (double a1, double a2, double a3, double b1, double b2, double b3, double *s1, double *s2, double *s3)
{

  *s1 = a2 * b2 - a3 * b2;
  *s2 = a3 * b1 - a1 * b3;
  *s3 = a1 * b2 - a2 * b1;

  return(0);

}

/*-----------------------------------------------------*/

double
dotproductangle (double u1, double u2, double u3, double v1, double v2, double v3)
{

  double polarAngle, dp;

  dp = u1 * v1 + u2 * v2 + u3 * v3;
  if (dp < -1.0)
    dp = -1.0;
  if (dp > 1.0)
    dp = 1.0;
  polarAngle = acos (dp);

//  assert (polarAngle >= 0 && polarAngle <= (2 * M_PI));

  return (polarAngle);

}

/*-----------------------------------------------------*/

double
vectorlen (double u1, double u2, double u3)
{

  double rr;

  rr = u1 * u1 + u2 * u2 + u3 * u3;
  rr = sqrt (rr);

  return (rr);

}


/*-----------------------------------------------------*/

double
unitvector (double *u1, double *u2, double *u3)
{

  double rr;

  rr = (*u1) * (*u1) + (*u2) * (*u2) + (*u3) * (*u3);
  rr = sqrt (rr);
  *u1 /= rr;
  *u2 /= rr;
  *u3 /= rr;

  return (rr);

}

/*-----------------------------------------------------*/


int
test_fun ()
{

  double a1, b1, c1;
  double a2, b2, c2;
  double a3, b3, c3;

  a1 = 1, b1 = 0, c1 = 0;
  a2 = 0, b2 = 1, c2 = 0;

  crossprod (a1, b1, c1, a2, b2, c2, &a3, &b3, &c3);
  printf ("1: %10.6f %10.6f %10.6f\n", a1, b1, c1);
  printf ("2: %10.6f %10.6f %10.6f\n", a2, b2, c2);
  printf ("3: %10.6f %10.6f %10.6f\n", a3, b3, c3);

/* done */

  exit (0);

}

/* ----------------------------------------------------------------- */

double
rad2deg (double ang)
{
  return (ang * 180 / M_PI);
};

/* ----------------------------------------------------------------- */

int
check_coord_sys (double x1, double x2, double x3, double y1, double y2, double y3, double z1, double z2, double z3)
{

  /* declarations */

  double d1;
  double zz1, zz2, zz3;

  /* length */

  d1 = vectorlen (x1, x2, x3);
  if (fabs (1.0 - d1) > MAXERR)
    {
      printf ("x-axis does not have unit length, is %f\n", d1);
      exit (1);
    };
  d1 = vectorlen (y1, y2, y3);
  if (fabs (1.0 - d1) > MAXERR)
    {
      printf ("y-axis does not have unit length, is %f\n", d1);
      exit (1);
    };
  d1 = vectorlen (z1, z2, z3);
  if (fabs (1.0 - d1) > MAXERR)
    {
      printf ("z-axis does not have unit length, is %f\n", d1);
      exit (1);
    };

  /* dot products, angles between the three axis */

  d1 = dotproductangle (x1, x2, x3, y1, y2, y3);
  if (fabs (M_PI / 2.0 - d1) > MAXERR)
    {
      printf ("x-axis and y-axis not perpendicular, angle is %f (deg)\n", rad2deg (d1));
      exit (1);
    };

  d1 = dotproductangle (x1, x2, x3, z1, z2, z3);
  if (fabs (M_PI / 2.0 - d1) > MAXERR)
    {
      printf ("x-axis and z-axis not perpendicular, angle is %f (deg)\n", rad2deg (d1));
      exit (1);
    };

  d1 = dotproductangle (y1, y2, y3, z1, z2, z3);
  if (fabs (M_PI / 2.0 - d1) > MAXERR)
    {
      printf ("y-axis and z-axis not perpendicular, angle is %f (deg)\n", rad2deg (d1));
      exit (1);
    };

  /* check right handed */

  crossprod (x1, x2, x3, y1, y2, y3, &zz1, &zz2, &zz3);
  d1 = dotproductangle (z1, z2, z3, zz1, zz2, zz3);
  if (fabs (d1) > MAXERR)
    {
      printf ("coordinate system is not right handed, z-axis off by %f deg\n", rad2deg (d1));
     fflush(stdout);
     exit (1);
    };

  return(0);

};

/* ----------------------------------------------------------------- */

int
check_unitvector (double x1, double x2, double x3)
{
  double d1;

  d1 = vectorlen (x1, x2, x3);
  if (fabs (1.0 - d1) > MAXERR)
    {
      printf ("vector does not have unit length, is %f\n", d1);
      fflush(stdout);
      exit (1);
    };

  return(0);

};

/* ----------------------------------------------------------------- */

int
coord_in_new (double x1, double x2, double x3,
              double y1, double y2, double y3,
              double z1, double z2, double z3, double xx, double yy, double zz, double *xx1, double *yy1, double *zz1)
{
  *xx1 = xx * x1 + yy * x2 + zz * x3;
  *yy1 = xx * y1 + yy * y2 + zz * y3;
  *zz1 = xx * z1 + yy * z2 + zz * z3;
  return (0);
};
