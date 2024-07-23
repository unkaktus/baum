/* VectorLaplaceFlat.c */
/* Copyright (C) 1998 Bernd Bruegmann, 5.3.2012 */
/* Produced with Mathematica */

#include "bam.h"

#define Power(x,y) (pow((double) (x), (double) (y)))
#define Log(x)     log((double) (x)))
#define pow2(x)    ((x)*(x))
#define pow2inv(x) (1.0/((x)*(x)))

#define Cal(x,y,z) ((x)?(y):(z))




void VectorLaplaceFlat(tL *level, tVarList *vllu, tVarList *vlu)
{

double *Lv1 = vldataptr(vllu, 0);
double *Lv2 = vldataptr(vllu, 1);
double *Lv3 = vldataptr(vllu, 2);
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

double ddv111;
double ddv112;
double ddv113;
double ddv121;
double ddv122;
double ddv123;
double ddv131;
double ddv132;
double ddv133;
double ddv221;
double ddv222;
double ddv223;
double ddv231;
double ddv232;
double ddv233;
double ddv331;
double ddv332;
double ddv333;
double oot;



forinner19(level) {


oot
=
0.33333333333333333333
;

ddv111
=
oodx2*(-2.*v1[ijk] + v1[-di + ijk] + v1[di + ijk])
;

ddv112
=
oodx2*(-2.*v2[ijk] + v2[-di + ijk] + v2[di + ijk])
;

ddv113
=
oodx2*(-2.*v3[ijk] + v3[-di + ijk] + v3[di + ijk])
;

ddv121
=
oo4dxdy*(v1[-di - dj + ijk] - v1[di - dj + ijk] - v1[-di + dj + ijk] + 
    v1[di + dj + ijk])
;

ddv122
=
oo4dxdy*(v2[-di - dj + ijk] - v2[di - dj + ijk] - v2[-di + dj + ijk] + 
    v2[di + dj + ijk])
;

ddv123
=
oo4dxdy*(v3[-di - dj + ijk] - v3[di - dj + ijk] - v3[-di + dj + ijk] + 
    v3[di + dj + ijk])
;

ddv131
=
oo4dxdz*(v1[-di - dk + ijk] - v1[di - dk + ijk] - v1[-di + dk + ijk] + 
    v1[di + dk + ijk])
;

ddv132
=
oo4dxdz*(v2[-di - dk + ijk] - v2[di - dk + ijk] - v2[-di + dk + ijk] + 
    v2[di + dk + ijk])
;

ddv133
=
oo4dxdz*(v3[-di - dk + ijk] - v3[di - dk + ijk] - v3[-di + dk + ijk] + 
    v3[di + dk + ijk])
;

ddv221
=
oody2*(-2.*v1[ijk] + v1[-dj + ijk] + v1[dj + ijk])
;

ddv222
=
oody2*(-2.*v2[ijk] + v2[-dj + ijk] + v2[dj + ijk])
;

ddv223
=
oody2*(-2.*v3[ijk] + v3[-dj + ijk] + v3[dj + ijk])
;

ddv231
=
oo4dydz*(v1[-dj - dk + ijk] - v1[dj - dk + ijk] - v1[-dj + dk + ijk] + 
    v1[dj + dk + ijk])
;

ddv232
=
oo4dydz*(v2[-dj - dk + ijk] - v2[dj - dk + ijk] - v2[-dj + dk + ijk] + 
    v2[dj + dk + ijk])
;

ddv233
=
oo4dydz*(v3[-dj - dk + ijk] - v3[dj - dk + ijk] - v3[-dj + dk + ijk] + 
    v3[dj + dk + ijk])
;

ddv331
=
oodz2*(-2.*v1[ijk] + v1[-dk + ijk] + v1[dk + ijk])
;

ddv332
=
oodz2*(-2.*v2[ijk] + v2[-dk + ijk] + v2[dk + ijk])
;

ddv333
=
oodz2*(-2.*v3[ijk] + v3[-dk + ijk] + v3[dk + ijk])
;

Lv1[ijk]
=
ddv111 + ddv221 + ddv331 + (ddv111 + ddv122 + ddv133)*oot
;

Lv2[ijk]
=
ddv112 + ddv222 + ddv332 + (ddv121 + ddv222 + ddv233)*oot
;

Lv3[ijk]
=
ddv113 + ddv223 + ddv333 + (ddv131 + ddv232 + ddv333)*oot
;



} endfor;  /* loop i, j, k */

}  /* function */

/* VectorLaplaceFlat.c */
/* nvars = 6, nauxs = 19, n* = 63,  n/ = 20,  n+ = 189, n = 272, O = 0 */
