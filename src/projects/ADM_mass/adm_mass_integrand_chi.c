/* adm_mass_integrand_chi.c */
/* Copyright (C) 1998 Bernd Bruegmann, 12.3.2012 */
/* Produced with Mathematica */

#include "bam.h"
#include "ADM_mass.h"

#define Power(x,y) (pow((double) (x), (double) (y)))
#define pow2(x)    ((x)*(x))
#define pow3(x)    ((x)*(x)*(x))
#define Sqrt(x)    sqrt(x)
#define Tanh(x)    tanh(x)
#define Sech(x)    (1/cosh(x))
#define PI  3.1415926535897932




void adm_mass_integrand_chi(tL *level, int i_x, int i_chi, int i_integrand)
{

double *x1 = level->v[i_x+0];
double *x2 = level->v[i_x+1];
double *x3 = level->v[i_x+2];
double *chi = level->v[i_chi+0];
double *integrand = level->v[i_integrand+0];
const double *xp = level->v[Ind("x")];
const double *yp = level->v[Ind("y")];
const double *zp = level->v[Ind("z")];
const int useShellsTransfo = level->shells;
const double *rp = level->v[IndLax("shells_R")];
const double *rr = level->v[IndLax("shells_r")];
const double shellsS = GetdLax("amr_shells_stretch");
const double shellsR = GetdLax("amr_shells_r0");
const double shellsE = GetdLax("amr_shells_eps");

double dx = level->dx;
double dy = level->dy;
double dz = level->dz;
double oo2dx = 1/(2*dx), oo2dy = 1/(2*dy), oo2dz = 1/(2*dz);
double oodx2 = 1/(dx*dx), oody2 = 1/(dy*dy), oodz2 = 1/(dz*dz);
double oo4dxdy = 1/(4*dx*dy);
double oo4dxdz = 1/(4*dx*dz);
double oo4dydz = 1/(4*dy*dz);

double ddRdr = 0.;
double delchi1 = 0.;
double delchi2 = 0.;
double delchi3 = 0.;
double delchiSST1 = 0.;
double delchiSST2 = 0.;
double delchiSST3 = 0.;
double dRdr = 0.;
double Jac11 = 0.;
double Jac12 = 0.;
double Jac13 = 0.;
double Jac21 = 0.;
double Jac22 = 0.;
double Jac23 = 0.;
double Jac31 = 0.;
double Jac32 = 0.;
double Jac33 = 0.;
double NN1 = 0.;
double NN2 = 0.;
double NN3 = 0.;
double R = 0.;



forinnerpoints_ijk(level) {


delchi1
=
oo2dx*(-chi[-di + ijk] + chi[di + ijk])
;

delchi2
=
oo2dy*(-chi[-dj + ijk] + chi[dj + ijk])
;

delchi3
=
oo2dz*(-chi[-dk + ijk] + chi[dk + ijk])
;


if (useShellsTransfo) { 

dRdr
=
shellsS + 0.5*(1. + Power(2.7182818284590452354,(-2.*shellsR)/shellsE))*
   (1. - shellsS + (-1 + shellsS)*Tanh((-shellsR + rr[ijk])/shellsE))
;

ddRdr
=
(0.5*(1. + Power(2.7182818284590452354,(-2.*shellsR)/shellsE))*
    (-1 + shellsS)*pow2(Sech((-shellsR + rr[ijk])/shellsE)))/shellsE
;


if (nbox<=1 ) { 

Jac11
=
xp[ijk]/(dRdr*rp[ijk])
;

Jac12
=
yp[ijk]/(dRdr*rp[ijk])
;

Jac13
=
zp[ijk]/(dRdr*rp[ijk])
;

Jac21
=
-(yp[ijk]/(pow2(xp[ijk]) + pow2(yp[ijk])))
;

Jac22
=
xp[ijk]/(pow2(xp[ijk]) + pow2(yp[ijk]))
;

Jac23
=
0
;

Jac31
=
-(zp[ijk]/(pow2(xp[ijk]) + pow2(zp[ijk])))
;

Jac32
=
0
;

Jac33
=
xp[ijk]/(pow2(xp[ijk]) + pow2(zp[ijk]))
;


} else if (nbox<=3) { 

Jac11
=
xp[ijk]/(dRdr*rp[ijk])
;

Jac12
=
yp[ijk]/(dRdr*rp[ijk])
;

Jac13
=
zp[ijk]/(dRdr*rp[ijk])
;

Jac21
=
yp[ijk]/(pow2(xp[ijk]) + pow2(yp[ijk]))
;

Jac22
=
-(xp[ijk]/(pow2(xp[ijk]) + pow2(yp[ijk])))
;

Jac23
=
0
;

Jac31
=
0
;

Jac32
=
-(zp[ijk]/(pow2(yp[ijk]) + pow2(zp[ijk])))
;

Jac33
=
yp[ijk]/(pow2(yp[ijk]) + pow2(zp[ijk]))
;


} else if (nbox<=5) { 

Jac11
=
xp[ijk]/(dRdr*rp[ijk])
;

Jac12
=
yp[ijk]/(dRdr*rp[ijk])
;

Jac13
=
zp[ijk]/(dRdr*rp[ijk])
;

Jac21
=
zp[ijk]/(pow2(xp[ijk]) + pow2(zp[ijk]))
;

Jac22
=
0
;

Jac23
=
-(xp[ijk]/(pow2(xp[ijk]) + pow2(zp[ijk])))
;

Jac31
=
0
;

Jac32
=
zp[ijk]/(pow2(yp[ijk]) + pow2(zp[ijk]))
;

Jac33
=
-(yp[ijk]/(pow2(yp[ijk]) + pow2(zp[ijk])))
;


} 

delchiSST1
=
delchi1*Jac11 + delchi2*Jac21 + delchi3*Jac31
;

delchiSST2
=
delchi1*Jac12 + delchi2*Jac22 + delchi3*Jac32
;

delchiSST3
=
delchi1*Jac13 + delchi2*Jac23 + delchi3*Jac33
;

delchi1
=
delchiSST1
;

delchi2
=
delchiSST2
;

delchi3
=
delchiSST3
;


} 

R
=
sqrt(pow2(x1[ijk]) + pow2(x2[ijk]) + pow2(x3[ijk]))
;

NN1
=
x1[ijk]/R
;

NN2
=
x2[ijk]/R
;

NN3
=
x3[ijk]/R
;

integrand[ijk]
=
(0.125*(delchi1*NN1 + delchi2*NN2 + delchi3*NN3)*pow2(R))/
  (PI*Power(chi[ijk],1.25))
;



} endfor_ijk; /* loop i, j, k */

}  /* function */

/* adm_mass_integrand_chi.c */
/* nvars = 10, nauxs = 21, n* = 71,  n/ = 50,  n+ = 77, n = 198, O = 0 */
