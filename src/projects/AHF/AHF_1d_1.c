/* AHF_1d_1.c */
/* Copyright (C) 1998 Bernd Bruegmann, 27.5.2009 */
/* Produced with Mathematica */

#include "bam.h"
#include "AHF.h"

#define Power(x,y) (pow((double) (x), (double) (y)))
#define Log(x)     log((double) (x)))
#define pow2(x)    ((x)*(x))
#define pow2inv(x) (1.0/((x)*(x)))

#define Cal(x,y,z) ((x)?(y):(z))




void AHF_1d_1(tVarList *vars, tVarList *ahfvars)
{

tL *level = ahfvars->level;
double *xp1  = Ptr(level, "x");
double *xp2  = Ptr(level, "y");
double *xp3  = Ptr(level, "z");

double *g11 = vldataptr(vars, 0);
double *g12 = vldataptr(vars, 1);
double *g13 = vldataptr(vars, 2);
double *g22 = vldataptr(vars, 3);
double *g23 = vldataptr(vars, 4);
double *g33 = vldataptr(vars, 5);
double *A11 = vldataptr(vars, 6);
double *A12 = vldataptr(vars, 7);
double *A13 = vldataptr(vars, 8);
double *A22 = vldataptr(vars, 9);
double *A23 = vldataptr(vars, 10);
double *A33 = vldataptr(vars, 11);
double *K = vldataptr(vars, 12);
double *phi = vldataptr(vars, 13);
double *s1 = vldataptr(ahfvars, 0);
double *s2 = vldataptr(ahfvars, 1);
double *s3 = vldataptr(ahfvars, 2);
double *H = vldataptr(ahfvars, 3);


double dx = level->dx;
double dy = level->dy;
double dz = level->dz;
double oo2dx = 1/(2*dx), oo2dy = 1/(2*dy), oo2dz = 1/(2*dz);
double oodx2 = 1/(dx*dx), oody2 = 1/(dy*dy), oodz2 = 1/(dz*dz);
double oo4dxdy = 1/(4*dx*dy);
double oo4dxdz = 1/(4*dx*dz);
double oo4dydz = 1/(4*dy*dz);
long ccc;

double norm;
double S1;
double S2;
double S3;



forinner25(level,oo2dx,oo2dy,oo2dz) {


S1
=
xp1[ccc]
;

S2
=
xp2[ccc]
;

S3
=
xp3[ccc]
;

norm
=
Power(2.7182818284590452354,4.*phi[ccc])*
  (2.*(S1*(S2*g12[ccc] + S3*g13[ccc]) + S2*S3*g23[ccc]) + 
    g11[ccc]*pow2(S1) + g22[ccc]*pow2(S2) + g33[ccc]*pow2(S3))
;

s1[ccc]
=
xp1[ccc]/sqrt(norm)
;

s2[ccc]
=
xp2[ccc]/sqrt(norm)
;

s3[ccc]
=
xp3[ccc]/sqrt(norm)
;



} endforinner; /* loop i, j, k */

}  /* function */

/* AHF_1d_1.c */
/* nvars = 21, nauxs = 4, n* = 59,  n/ = 23,  n+ = 9, n = 91, O = 0 */
