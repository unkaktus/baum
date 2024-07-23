/* compute_hydrovars.c */
/* Copyright (C) 1998 Bernd Bruegmann, 24.8.2023 */
/* Produced with Mathematica */

#include "bam.h"
#include "hydroanalysis.h"

#define Power(x,y) pow((double) (x), (double) (y))
#define Sqrt(x)    sqrt((double) (x))
#define Log(x)     log((double) (x))
#define pow2(x)    ((x)*(x))
#define pow4(x)    ((x)*(x)*(x)*(x))
#define pow2inv(x) (1.0/((x)*(x)))

#define Tanh(x)    tanh(x)
#define Sech(x)    (1/cosh(x))
#define Cal(x,y,z) ((x)?(y):(z))


#define ALPHAMIN 1e-6



void compute_hydrovars(tVarList *u)
{

double *conD = vldataptr(u, 0);
double *pres = vldataptr(u, 1);
double *rho = vldataptr(u, 2);
double *epsl = vldataptr(u, 3);
double *v1 = vldataptr(u, 4);
double *v2 = vldataptr(u, 5);
double *v3 = vldataptr(u, 6);
double *g11 = vldataptr(u, 7);
double *g12 = vldataptr(u, 8);
double *g13 = vldataptr(u, 9);
double *g22 = vldataptr(u, 10);
double *g23 = vldataptr(u, 11);
double *g33 = vldataptr(u, 12);
double *alpha = vldataptr(u, 13);
double *beta1 = vldataptr(u, 14);
double *beta2 = vldataptr(u, 15);
double *beta3 = vldataptr(u, 16);
double *xp = vldataptr(u, 17);
double *yp = vldataptr(u, 18);
double *zp = vldataptr(u, 19);
double *Dh = vldataptr(u, 20);
double *Db = vldataptr(u, 21);
double *Du = vldataptr(u, 22);
double *etot = vldataptr(u, 23);
double *uesc = vldataptr(u, 24);
double *vor1 = vldataptr(u, 25);
double *vor2 = vldataptr(u, 26);
double *vor3 = vldataptr(u, 27);
double *ivor1 = vldataptr(u, 28);
double *ivor2 = vldataptr(u, 29);
double *ivor3 = vldataptr(u, 30);
double *P1 = vldataptr(u, 31);
double *P2 = vldataptr(u, 32);
double *P3 = vldataptr(u, 33);
double *Pu1 = vldataptr(u, 34);
double *Pu2 = vldataptr(u, 35);
double *Pu3 = vldataptr(u, 36);
double *ut = vldataptr(u, 37);
tL *level = u->level;
int index_matter_mask = Ind("matter_mask");
double *mask = level->v[index_matter_mask + 0];

const double Wmax = Getd("hydroanalysis_wMAX");



bampi_openmp_start


double dx = level->dx;
double dy = level->dy;
double dz = level->dz;
double oo2dx = 1/(2*dx), oo2dy = 1/(2*dy), oo2dz = 1/(2*dz);
double oodx2 = 1/(dx*dx), oody2 = 1/(dy*dy), oodz2 = 1/(dz*dz);
double oo4dxdy = 1/(4*dx*dy);
double oo4dxdz = 1/(4*dx*dz);
double oo4dydz = 1/(4*dy*dz);

double auupt = 0.;
double betadown1 = 0.;
double betadown2 = 0.;
double betadown3 = 0.;
double betaduu = 0.;
double betasq = 0.;
double h = 0.;
double ooalpha = 0.;
double udown1 = 0.;
double udown2 = 0.;
double udown3 = 0.;
double udownt = 0.;
double ur = 0.;
double uup1 = 0.;
double uup2 = 0.;
double uup3 = 0.;
double uupt = 0.;
double vsq = 0.;
double vup1 = 0.;
double vup2 = 0.;
double vup3 = 0.;
double W = 0.;



forinnerpoints_ijk_openmp(level) {



 if ((MATTER.USEMASK) && (mask[ijk]>.99)) { 

Db[ijk]
=
0
;

Du[ijk]
=
0
;

Dh[ijk]
=
0
;

ivor1[ijk]
=
0
;

ivor2[ijk]
=
0
;

ivor3[ijk]
=
0
;


 continue; 


 } 

vup1
=
v1[ijk]
;

vup2
=
v2[ijk]
;

vup3
=
v3[ijk]
;

vsq
=
2.*(vup1*(vup2*g12[ijk] + vup3*g13[ijk]) + vup2*vup3*g23[ijk]) + 
  g11[ijk]*pow2(vup1) + g22[ijk]*pow2(vup2) + g33[ijk]*pow2(vup3)
;

ooalpha
=
1/alpha[ijk]
;


if (alpha[ijk]<=ALPHAMIN){ 


 printf(" alpha<%e in hydroanalysis, reset\n",ALPHAMIN); ooalpha=1./ALPHAMIN; 


} 

W
=
1/Sqrt(1. - vsq)
;


if (vsq>=1.) W = Wmax; 

uup1
=
W*(vup1 - ooalpha*beta1[ijk])
;

uup2
=
W*(vup2 - ooalpha*beta2[ijk])
;

uup3
=
W*(vup3 - ooalpha*beta3[ijk])
;

uupt
=
ooalpha*W
;

betadown1
=
beta1[ijk]*g11[ijk] + beta2[ijk]*g12[ijk] + beta3[ijk]*g13[ijk]
;

betadown2
=
beta1[ijk]*g12[ijk] + beta2[ijk]*g22[ijk] + beta3[ijk]*g23[ijk]
;

betadown3
=
beta1[ijk]*g13[ijk] + beta2[ijk]*g23[ijk] + beta3[ijk]*g33[ijk]
;

udown1
=
betadown1*uupt + uup1*g11[ijk] + uup2*g12[ijk] + uup3*g13[ijk]
;

udown2
=
betadown2*uupt + uup1*g12[ijk] + uup2*g22[ijk] + uup3*g23[ijk]
;

udown3
=
betadown3*uupt + uup1*g13[ijk] + uup2*g23[ijk] + uup3*g33[ijk]
;

h
=
1. + epsl[ijk] + pres[ijk]/rho[ijk]
;

betasq
=
betadown1*beta1[ijk] + betadown2*beta2[ijk] + betadown3*beta3[ijk]
;

betaduu
=
betadown1*uup1 + betadown2*uup2 + betadown3*uup3
;

udownt
=
betaduu + uupt*(betasq - pow2(alpha[ijk]))
;

ut[ijk]
=
udownt
;

ur
=
vup1*xp[ijk] + vup2*yp[ijk] + vup3*zp[ijk]
;

P1[ijk]
=
h*udown1*conD[ijk]
;

P2[ijk]
=
h*udown2*conD[ijk]
;

P3[ijk]
=
h*udown3*conD[ijk]
;

ivor1[ijk]
=
h*udown1
;

ivor2[ijk]
=
h*udown2
;

ivor3[ijk]
=
h*udown3
;

Db[ijk]
=
conD[ijk]
;

Du[ijk]
=
conD[ijk]
;

Dh[ijk]
=
conD[ijk]
;


 if ( (- udownt > 1)  && (ur > 0) ){ 

Db[ijk]
=
0
;

auupt
=
uupt*alpha[ijk]
;

etot[ijk]
=
conD[ijk]*(auupt*h - pres[ijk]/(auupt*rho[ijk]))
;

uesc[ijk]
=
conD[ijk]*epsl[ijk]
;

Pu1[ijk]
=
P1[ijk]
;

Pu2[ijk]
=
P2[ijk]
;

Pu3[ijk]
=
P3[ijk]
;


 } else { 

Du[ijk]
=
0
;

uesc[ijk]
=
0
;

etot[ijk]
=
0
;

Pu1[ijk]
=
0
;

Pu2[ijk]
=
0
;

Pu3[ijk]
=
0
;


 }  



} endfor_ijk_openmp; /* loop i, j, k */



bampi_openmp_stop


}  /* function */

/* compute_hydrovars.c */
/* nvars = 38, nauxs = 22, n* = 129,  n/ = 26,  n+ = 43, n = 198, O = 0 */
