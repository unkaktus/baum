/* z4_algcon.c */
/* Copyright (C) 1998 Bernd Bruegmann, 24.1.2024 */
/* Produced with Mathematica */

#include "bam.h"
#include "z4.h"

#define Power(x,y) (pow((double) (x), (double) (y)))
#define Log(x)     log((double) (x)))
#define pow2(x)    ((x)*(x))
#define pow2inv(x) (1.0/((x)*(x)))

#define Cal(x,y,z) ((x)?(y):(z))

#define PROBLEM    printf("  %d pts away from boundary\n",boundaryaway(6)); \
                   printf("  %e %e %e \n  detg=%e,trA=%e \n", \
                   Ptr(level,"x")[ijk],Ptr(level,"y")[ijk],Ptr(level,"z")[ijk], \
                   detg, trA); stop++ 



void z4_algcon(tVarList *u)
{

double *g11 = vldataptr(u, 0);
double *g12 = vldataptr(u, 1);
double *g13 = vldataptr(u, 2);
double *g22 = vldataptr(u, 3);
double *g23 = vldataptr(u, 4);
double *g33 = vldataptr(u, 5);
double *A11 = vldataptr(u, 6);
double *A12 = vldataptr(u, 7);
double *A13 = vldataptr(u, 8);
double *A22 = vldataptr(u, 9);
double *A23 = vldataptr(u, 10);
double *A33 = vldataptr(u, 11);
double *trAstore = vldataptr(u, 12);
double *detgstore = vldataptr(u, 13);
tL *level = u->level;

const int subtractA      = Getv("z4_subtractA", "yes");
const int normalizedetg  = Getv("z4_normalizedetg", "yes");
const double oot = 1.0/3.0;
const int store          = Getv("z4_register_algcon", "store");



bampi_openmp_start


double dx = level->dx;
double dy = level->dy;
double dz = level->dz;
double oo2dx = 1/(2*dx), oo2dy = 1/(2*dy), oo2dz = 1/(2*dz);
double oodx2 = 1/(dx*dx), oody2 = 1/(dy*dy), oodz2 = 1/(dz*dz);
double oo4dxdy = 1/(4*dx*dy);
double oo4dxdz = 1/(4*dx*dz);
double oo4dydz = 1/(4*dy*dz);

double aux = 0.;
double detg = 0.;
double detginv = 0.;
double epsilon = 0.;
double trA = 0.;



int stop=0; forallpoints_ijk_openmp(level) {


detg
=
2.*g12[ijk]*g13[ijk]*g23[ijk] + g11[ijk]*g22[ijk]*g33[ijk] - 
  g33[ijk]*pow2(g12[ijk]) - g22[ijk]*pow2(g13[ijk]) - g11[ijk]*pow2(g23[ijk])
;


if (detg<=0.) { PROBLEM; /*errorexit(" detg <= 0.");*/ detg=1.;}  

detginv
=
1/detg
;



/* conditional */
if (normalizedetg) {

epsilon
=
-1. + detg
;


if (fabs(epsilon)<1e-12) aux = 1.-oot*epsilon; else  

aux
=
Power(detg,-oot)
;

g11[ijk]
=
aux*g11[ijk]
;

g12[ijk]
=
aux*g12[ijk]
;

g13[ijk]
=
aux*g13[ijk]
;

g22[ijk]
=
aux*g22[ijk]
;

g23[ijk]
=
aux*g23[ijk]
;

g33[ijk]
=
aux*g33[ijk]
;

}
/* if (normalizedetg) */


detg
=
2.*g12[ijk]*g13[ijk]*g23[ijk] + g11[ijk]*g22[ijk]*g33[ijk] - 
  g33[ijk]*pow2(g12[ijk]) - g22[ijk]*pow2(g13[ijk]) - g11[ijk]*pow2(g23[ijk])
;


if (detg<=0.) { PROBLEM; /*errorexit(" detg <= 0.");*/ detg=1.;}  

detginv
=
1/detg
;

trA
=
detginv*(-2.*A23[ijk]*g11[ijk]*g23[ijk] + A11[ijk]*g22[ijk]*g33[ijk] + 
    g11[ijk]*(A33[ijk]*g22[ijk] + A22[ijk]*g33[ijk]) + 
    2.*(g13[ijk]*(A23[ijk]*g12[ijk] - A13[ijk]*g22[ijk] + 
          A12[ijk]*g23[ijk]) + 
       g12[ijk]*(A13[ijk]*g23[ijk] - A12[ijk]*g33[ijk])) - 
    A33[ijk]*pow2(g12[ijk]) - A22[ijk]*pow2(g13[ijk]) - 
    A11[ijk]*pow2(g23[ijk]))
;



/* conditional */
if (subtractA) {

aux
=
-(oot*trA)
;

A11[ijk]
=
A11[ijk] + aux*g11[ijk]
;

A12[ijk]
=
A12[ijk] + aux*g12[ijk]
;

A13[ijk]
=
A13[ijk] + aux*g13[ijk]
;

A22[ijk]
=
A22[ijk] + aux*g22[ijk]
;

A23[ijk]
=
A23[ijk] + aux*g23[ijk]
;

A33[ijk]
=
A33[ijk] + aux*g33[ijk]
;

}
/* if (subtractA) */


trA
=
detginv*(-2.*A23[ijk]*g11[ijk]*g23[ijk] + A11[ijk]*g22[ijk]*g33[ijk] + 
    g11[ijk]*(A33[ijk]*g22[ijk] + A22[ijk]*g33[ijk]) + 
    2.*(g13[ijk]*(A23[ijk]*g12[ijk] - A13[ijk]*g22[ijk] + 
          A12[ijk]*g23[ijk]) + 
       g12[ijk]*(A13[ijk]*g23[ijk] - A12[ijk]*g33[ijk])) - 
    A33[ijk]*pow2(g12[ijk]) - A22[ijk]*pow2(g13[ijk]) - 
    A11[ijk]*pow2(g23[ijk]))
;


if (store) {  

detgstore[ijk]
=
detg
;

trAstore[ijk]
=
trA
;


}  



} endfor_ijk_openmp;/* if (stop) errorexit(""); */ /* loop ccc */



bampi_openmp_stop


}  /* function */

/* z4_algcon.c */
/* nvars = 14, nauxs = 5, n* = 124,  n/ = 37,  n+ = 50, n = 211, O = 0 */
