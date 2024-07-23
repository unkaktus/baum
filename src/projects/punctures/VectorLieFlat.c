/* VectorLieFlat.c */
/* Copyright (C) 1998 Bernd Bruegmann, 25.6.2012 */
/* Produced with Mathematica */

#include "bam.h"

#define Power(x,y) (pow((double) (x), (double) (y)))
#define Log(x)     log((double) (x)))
#define pow2(x)    ((x)*(x))
#define pow2inv(x) (1.0/((x)*(x)))

#define Cal(x,y,z) ((x)?(y):(z))




void VectorLieFlat(tL *level, tVarList *vllu, tVarList *vlu)
{

double *lv11 = vldataptr(vllu, 0);
double *lv12 = vldataptr(vllu, 1);
double *lv13 = vldataptr(vllu, 2);
double *lv22 = vldataptr(vllu, 3);
double *lv23 = vldataptr(vllu, 4);
double *lv33 = vldataptr(vllu, 5);
double *v1 = vldataptr(vlu, 0);
double *v2 = vldataptr(vlu, 1);
double *v3 = vldataptr(vlu, 2);


double dx = level->dx;
double dy = level->dy;
double dz = level->dz;
double oo2dx = 1/(2*dx), oo2dy = 1/(2*dy), oo2dz = 1/(2*dz);
double oodx2 = 1/(dx*dx), oody2 = 1/(dy*dy), oodz2 = 1/(dz*dz);
double oo4dxdy = 1/(4*dx*dy);
double oo4dxdz = 1/(4*dx*dz);
double oo4dydz = 1/(4*dy*dz);

double dv11 = 0.;
double dv12 = 0.;
double dv13 = 0.;
double dv21 = 0.;
double dv22 = 0.;
double dv23 = 0.;
double dv31 = 0.;
double dv32 = 0.;
double dv33 = 0.;
double trdv = 0.;



forinner19(level) {


dv11
=
oo2dx*(-v1[-di + ijk] + v1[di + ijk])
;

dv12
=
oo2dx*(-v2[-di + ijk] + v2[di + ijk])
;

dv13
=
oo2dx*(-v3[-di + ijk] + v3[di + ijk])
;

dv21
=
oo2dy*(-v1[-dj + ijk] + v1[dj + ijk])
;

dv22
=
oo2dy*(-v2[-dj + ijk] + v2[dj + ijk])
;

dv23
=
oo2dy*(-v3[-dj + ijk] + v3[dj + ijk])
;

dv31
=
oo2dz*(-v1[-dk + ijk] + v1[dk + ijk])
;

dv32
=
oo2dz*(-v2[-dk + ijk] + v2[dk + ijk])
;

dv33
=
oo2dz*(-v3[-dk + ijk] + v3[dk + ijk])
;

trdv
=
0.66666666666666666667*(dv11 + dv22 + dv33)
;

lv11[ijk]
=
2.*dv11 - trdv
;

lv12[ijk]
=
dv12 + dv21
;

lv13[ijk]
=
dv13 + dv31
;

lv22[ijk]
=
2.*dv22 - trdv
;

lv23[ijk]
=
dv23 + dv32
;

lv33[ijk]
=
2.*dv33 - trdv
;



} endfor;  /* loop i, j, k */

}  /* function */

/* VectorLieFlat.c */
/* nvars = 9, nauxs = 10, n* = 49,  n/ = 20,  n+ = 56, n = 125, O = 0 */
