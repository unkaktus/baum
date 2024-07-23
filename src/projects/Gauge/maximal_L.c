/* maximal_L.c */
/* Copyright (C) 1998 Bernd Bruegmann, 12.3.2012 */
/* Produced with Mathematica */

#include "bam.h"
#include "Gauge.h"

#define Power(x,y) (pow((double) (x), (double) (y)))
#define Log(x)     log((double) (x)))
#define pow2(x)    ((x)*(x))
#define pow2inv(x) (1.0/((x)*(x)))

#define Cal(x,y,z) ((x)?(y):(z))




void maximal_L(tL *level, tVarList *vlLu, tVarList *vlu, tVarList *vlgi, tVarList *vlG, tVarList *vlKK)
{

double *Lu = vldataptr(vlLu, 0);
double *u = vldataptr(vlu, 0);
double *gi11 = vldataptr(vlgi, 0);
double *gi12 = vldataptr(vlgi, 1);
double *gi13 = vldataptr(vlgi, 2);
double *gi22 = vldataptr(vlgi, 3);
double *gi23 = vldataptr(vlgi, 4);
double *gi33 = vldataptr(vlgi, 5);
double *G1 = vldataptr(vlG, 0);
double *G2 = vldataptr(vlG, 1);
double *G3 = vldataptr(vlG, 2);
double *KK = vldataptr(vlKK, 0);


double dx = level->dx;
double dy = level->dy;
double dz = level->dz;
double oo2dx = 1/(2*dx), oo2dy = 1/(2*dy), oo2dz = 1/(2*dz);
double oodx2 = 1/(dx*dx), oody2 = 1/(dy*dy), oodz2 = 1/(dz*dz);
double oo4dxdy = 1/(4*dx*dy);
double oo4dxdz = 1/(4*dx*dz);
double oo4dydz = 1/(4*dy*dz);

double ddu11 = 0.;
double ddu12 = 0.;
double ddu13 = 0.;
double ddu22 = 0.;
double ddu23 = 0.;
double ddu33 = 0.;
double du1 = 0.;
double du2 = 0.;
double du3 = 0.;



forinner19(level) {


du1
=
oo2dx*(-u[-di + ijk] + u[di + ijk])
;

du2
=
oo2dy*(-u[-dj + ijk] + u[dj + ijk])
;

du3
=
oo2dz*(-u[-dk + ijk] + u[dk + ijk])
;

ddu11
=
oodx2*(-2.*u[ijk] + u[-di + ijk] + u[di + ijk])
;

ddu12
=
oo4dxdy*(u[-di - dj + ijk] - u[di - dj + ijk] - u[-di + dj + ijk] + 
    u[di + dj + ijk])
;

ddu13
=
oo4dxdz*(u[-di - dk + ijk] - u[di - dk + ijk] - u[-di + dk + ijk] + 
    u[di + dk + ijk])
;

ddu22
=
oody2*(-2.*u[ijk] + u[-dj + ijk] + u[dj + ijk])
;

ddu23
=
oo4dydz*(u[-dj - dk + ijk] - u[dj - dk + ijk] - u[-dj + dk + ijk] + 
    u[dj + dk + ijk])
;

ddu33
=
oodz2*(-2.*u[ijk] + u[-dk + ijk] + u[dk + ijk])
;

Lu[ijk]
=
-(du1*G1[ijk]) - du2*G2[ijk] - du3*G3[ijk] + ddu11*gi11[ijk] + 
  ddu22*gi22[ijk] + 2.*(ddu12*gi12[ijk] + ddu13*gi13[ijk] + 
     ddu23*gi23[ijk]) + ddu33*gi33[ijk] - KK[ijk]*(1. + u[ijk])
;



} endfor;  /* loop i, j, k */

}  /* function */

/* maximal_L.c */
/* nvars = 12, nauxs = 9, n* = 65,  n/ = 20,  n+ = 86, n = 171, O = 0 */
