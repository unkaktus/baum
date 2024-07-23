/* VectorLaplaceFlatGS.c */
/* Copyright (C) 1998 Bernd Bruegmann, 5.3.2012 */
/* Produced with Mathematica */

#include "bam.h"

#define Power(x,y) (pow((double) (x), (double) (y)))
#define Log(x)     log((double) (x)))
#define pow2(x)    ((x)*(x))
#define pow2inv(x) (1.0/((x)*(x)))

#define Cal(x,y,z) ((x)?(y):(z))




void VectorLaplaceFlatGS(tL *level, tVarList *vlv, tVarList *vlu, tVarList *vlf)
{

double *v1 = vldataptr(vlv, 0);
double *v2 = vldataptr(vlv, 1);
double *v3 = vldataptr(vlv, 2);
double *u1 = vldataptr(vlu, 0);
double *u2 = vldataptr(vlu, 1);
double *u3 = vldataptr(vlu, 2);
double *f1 = vldataptr(vlf, 0);
double *f2 = vldataptr(vlf, 1);
double *f3 = vldataptr(vlf, 2);


double dx = level->dx;
double dy = level->dy;
double dz = level->dz;
double oo2dx = 1/(2*dx), oo2dy = 1/(2*dy), oo2dz = 1/(2*dz);
double oodx2 = 1/(dx*dx), oody2 = 1/(dy*dy), oodz2 = 1/(dz*dz);
double oo4dxdy = 1/(4*dx*dy);
double oo4dxdz = 1/(4*dx*dz);
double oo4dydz = 1/(4*dy*dz);

double ddu111;
double ddu112;
double ddu113;
double ddu121;
double ddu122;
double ddu123;
double ddu131;
double ddu132;
double ddu133;
double ddu221;
double ddu222;
double ddu223;
double ddu231;
double ddu232;
double ddu233;
double ddu331;
double ddu332;
double ddu333;
double Lii;
double Lu1;
double Lu2;
double Lu3;
double oot;



forinner19(level) {


oot
=
0.33333333333333333333
;

ddu111
=
oodx2*(-2.*u1[ijk] + u1[-di + ijk] + u1[di + ijk])
;

ddu112
=
oodx2*(-2.*u2[ijk] + u2[-di + ijk] + u2[di + ijk])
;

ddu113
=
oodx2*(-2.*u3[ijk] + u3[-di + ijk] + u3[di + ijk])
;

ddu121
=
oo4dxdy*(u1[-di - dj + ijk] - u1[di - dj + ijk] - u1[-di + dj + ijk] + 
    u1[di + dj + ijk])
;

ddu122
=
oo4dxdy*(u2[-di - dj + ijk] - u2[di - dj + ijk] - u2[-di + dj + ijk] + 
    u2[di + dj + ijk])
;

ddu123
=
oo4dxdy*(u3[-di - dj + ijk] - u3[di - dj + ijk] - u3[-di + dj + ijk] + 
    u3[di + dj + ijk])
;

ddu131
=
oo4dxdz*(u1[-di - dk + ijk] - u1[di - dk + ijk] - u1[-di + dk + ijk] + 
    u1[di + dk + ijk])
;

ddu132
=
oo4dxdz*(u2[-di - dk + ijk] - u2[di - dk + ijk] - u2[-di + dk + ijk] + 
    u2[di + dk + ijk])
;

ddu133
=
oo4dxdz*(u3[-di - dk + ijk] - u3[di - dk + ijk] - u3[-di + dk + ijk] + 
    u3[di + dk + ijk])
;

ddu221
=
oody2*(-2.*u1[ijk] + u1[-dj + ijk] + u1[dj + ijk])
;

ddu222
=
oody2*(-2.*u2[ijk] + u2[-dj + ijk] + u2[dj + ijk])
;

ddu223
=
oody2*(-2.*u3[ijk] + u3[-dj + ijk] + u3[dj + ijk])
;

ddu231
=
oo4dydz*(u1[-dj - dk + ijk] - u1[dj - dk + ijk] - u1[-dj + dk + ijk] + 
    u1[dj + dk + ijk])
;

ddu232
=
oo4dydz*(u2[-dj - dk + ijk] - u2[dj - dk + ijk] - u2[-dj + dk + ijk] + 
    u2[dj + dk + ijk])
;

ddu233
=
oo4dydz*(u3[-dj - dk + ijk] - u3[dj - dk + ijk] - u3[-dj + dk + ijk] + 
    u3[dj + dk + ijk])
;

ddu331
=
oodz2*(-2.*u1[ijk] + u1[-dk + ijk] + u1[dk + ijk])
;

ddu332
=
oodz2*(-2.*u2[ijk] + u2[-dk + ijk] + u2[dk + ijk])
;

ddu333
=
oodz2*(-2.*u3[ijk] + u3[-dk + ijk] + u3[dk + ijk])
;

Lu1
=
ddu111 + ddu221 + ddu331 + (ddu111 + ddu122 + ddu133)*oot
;

Lu2
=
ddu112 + ddu222 + ddu332 + (ddu121 + ddu222 + ddu233)*oot
;

Lu3
=
ddu113 + ddu223 + ddu333 + (ddu131 + ddu232 + ddu333)*oot
;

Lii
=
-20.*oo2dx*oot
;

v1[ijk]
=
(-Lu1 + f1[ijk])/Lii + u1[ijk]
;

v2[ijk]
=
(-Lu2 + f2[ijk])/Lii + u2[ijk]
;

v3[ijk]
=
(-Lu3 + f3[ijk])/Lii + u3[ijk]
;



} endfor;  /* loop i, j, k */

}  /* function */

/* VectorLaplaceFlatGS.c */
/* nvars = 9, nauxs = 23, n* = 69,  n/ = 23,  n+ = 199, n = 291, O = 0 */
