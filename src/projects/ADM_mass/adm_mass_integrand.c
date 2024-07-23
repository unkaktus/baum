/* adm_mass_integrand.c */
/* Copyright (C) 1998 Bernd Bruegmann, 12.3.2012 */
/* Produced with Mathematica */

#include "bam.h"
#include "ADM_mass.h"

#define Power(x,y) (pow((double) (x), (double) (y)))
#define Sqrt(x)    sqrt(x)
#define Log(x)     log((double) (x)))
#define pow2(x)    ((x)*(x))
#define pow3(x)    ((x)*(x)*(x))
#define pow4(x)    ((x)*(x)*(x)*(x))
#define pow2inv(x) (1.0/((x)*(x)))

#define Tanh(x)    tanh(x)
#define Sech(x)    (1/cosh(x))
#define Cal(x,y,z) ((x)?(y):(z))

#define PI  3.1415926535897932




void adm_mass_integrand(tL *level, int i_x, int i_g, int i_integrand)
{

double *x1 = level->v[i_x+0];
double *x2 = level->v[i_x+1];
double *x3 = level->v[i_x+2];
double *g11 = level->v[i_g+0];
double *g12 = level->v[i_g+1];
double *g13 = level->v[i_g+2];
double *g22 = level->v[i_g+3];
double *g23 = level->v[i_g+4];
double *g33 = level->v[i_g+5];
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
double delg111 = 0.;
double delg112 = 0.;
double delg113 = 0.;
double delg122 = 0.;
double delg123 = 0.;
double delg133 = 0.;
double delg211 = 0.;
double delg212 = 0.;
double delg213 = 0.;
double delg222 = 0.;
double delg223 = 0.;
double delg233 = 0.;
double delg311 = 0.;
double delg312 = 0.;
double delg313 = 0.;
double delg322 = 0.;
double delg323 = 0.;
double delg333 = 0.;
double delgSST111 = 0.;
double delgSST112 = 0.;
double delgSST113 = 0.;
double delgSST122 = 0.;
double delgSST123 = 0.;
double delgSST133 = 0.;
double delgSST211 = 0.;
double delgSST212 = 0.;
double delgSST213 = 0.;
double delgSST222 = 0.;
double delgSST223 = 0.;
double delgSST233 = 0.;
double delgSST311 = 0.;
double delgSST312 = 0.;
double delgSST313 = 0.;
double delgSST322 = 0.;
double delgSST323 = 0.;
double delgSST333 = 0.;
double det = 0.;
double detginvf = 0.;
double dRdr = 0.;
double ginv11 = 0.;
double ginv12 = 0.;
double ginv13 = 0.;
double ginv22 = 0.;
double ginv23 = 0.;
double ginv33 = 0.;
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


delg111
=
oo2dx*(-g11[-di + ijk] + g11[di + ijk])
;

delg112
=
oo2dx*(-g12[-di + ijk] + g12[di + ijk])
;

delg113
=
oo2dx*(-g13[-di + ijk] + g13[di + ijk])
;

delg122
=
oo2dx*(-g22[-di + ijk] + g22[di + ijk])
;

delg123
=
oo2dx*(-g23[-di + ijk] + g23[di + ijk])
;

delg133
=
oo2dx*(-g33[-di + ijk] + g33[di + ijk])
;

delg211
=
oo2dy*(-g11[-dj + ijk] + g11[dj + ijk])
;

delg212
=
oo2dy*(-g12[-dj + ijk] + g12[dj + ijk])
;

delg213
=
oo2dy*(-g13[-dj + ijk] + g13[dj + ijk])
;

delg222
=
oo2dy*(-g22[-dj + ijk] + g22[dj + ijk])
;

delg223
=
oo2dy*(-g23[-dj + ijk] + g23[dj + ijk])
;

delg233
=
oo2dy*(-g33[-dj + ijk] + g33[dj + ijk])
;

delg311
=
oo2dz*(-g11[-dk + ijk] + g11[dk + ijk])
;

delg312
=
oo2dz*(-g12[-dk + ijk] + g12[dk + ijk])
;

delg313
=
oo2dz*(-g13[-dk + ijk] + g13[dk + ijk])
;

delg322
=
oo2dz*(-g22[-dk + ijk] + g22[dk + ijk])
;

delg323
=
oo2dz*(-g23[-dk + ijk] + g23[dk + ijk])
;

delg333
=
oo2dz*(-g33[-dk + ijk] + g33[dk + ijk])
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

delgSST111
=
delg111*Jac11 + delg211*Jac21 + delg311*Jac31
;

delgSST112
=
delg112*Jac11 + delg212*Jac21 + delg312*Jac31
;

delgSST113
=
delg113*Jac11 + delg213*Jac21 + delg313*Jac31
;

delgSST122
=
delg122*Jac11 + delg222*Jac21 + delg322*Jac31
;

delgSST123
=
delg123*Jac11 + delg223*Jac21 + delg323*Jac31
;

delgSST133
=
delg133*Jac11 + delg233*Jac21 + delg333*Jac31
;

delgSST211
=
delg111*Jac12 + delg211*Jac22 + delg311*Jac32
;

delgSST212
=
delg112*Jac12 + delg212*Jac22 + delg312*Jac32
;

delgSST213
=
delg113*Jac12 + delg213*Jac22 + delg313*Jac32
;

delgSST222
=
delg122*Jac12 + delg222*Jac22 + delg322*Jac32
;

delgSST223
=
delg123*Jac12 + delg223*Jac22 + delg323*Jac32
;

delgSST233
=
delg133*Jac12 + delg233*Jac22 + delg333*Jac32
;

delgSST311
=
delg111*Jac13 + delg211*Jac23 + delg311*Jac33
;

delgSST312
=
delg112*Jac13 + delg212*Jac23 + delg312*Jac33
;

delgSST313
=
delg113*Jac13 + delg213*Jac23 + delg313*Jac33
;

delgSST322
=
delg122*Jac13 + delg222*Jac23 + delg322*Jac33
;

delgSST323
=
delg123*Jac13 + delg223*Jac23 + delg323*Jac33
;

delgSST333
=
delg133*Jac13 + delg233*Jac23 + delg333*Jac33
;

delg111
=
delgSST111
;

delg112
=
delgSST112
;

delg113
=
delgSST113
;

delg122
=
delgSST122
;

delg123
=
delgSST123
;

delg133
=
delgSST133
;

delg211
=
delgSST211
;

delg212
=
delgSST212
;

delg213
=
delgSST213
;

delg222
=
delgSST222
;

delg223
=
delgSST223
;

delg233
=
delgSST233
;

delg311
=
delgSST311
;

delg312
=
delgSST312
;

delg313
=
delgSST313
;

delg322
=
delgSST322
;

delg323
=
delgSST323
;

delg333
=
delgSST333
;


} 

det
=
2.*g12[ijk]*g13[ijk]*g23[ijk] + g11[ijk]*g22[ijk]*g33[ijk] - 
  g33[ijk]*pow2(g12[ijk]) - g22[ijk]*pow2(g13[ijk]) - g11[ijk]*pow2(g23[ijk])
;

detginvf
=
1/(2.*g12[ijk]*g13[ijk]*g23[ijk] + g11[ijk]*g22[ijk]*g33[ijk] - 
    g33[ijk]*pow2(g12[ijk]) - g22[ijk]*pow2(g13[ijk]) - 
    g11[ijk]*pow2(g23[ijk]))
;

ginv11
=
detginvf*(g22[ijk]*g33[ijk] - pow2(g23[ijk]))
;

ginv12
=
detginvf*(g13[ijk]*g23[ijk] - g12[ijk]*g33[ijk])
;

ginv13
=
detginvf*(-(g13[ijk]*g22[ijk]) + g12[ijk]*g23[ijk])
;

ginv22
=
detginvf*(g11[ijk]*g33[ijk] - pow2(g13[ijk]))
;

ginv23
=
detginvf*(g12[ijk]*g13[ijk] - g11[ijk]*g23[ijk])
;

ginv33
=
detginvf*(g11[ijk]*g22[ijk] - pow2(g12[ijk]))
;

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
((ginv12*(0.0625*((-delg112 + delg211)*NN1 + (delg122 - delg212)*NN2) + 
         (0.0625*(delg123 + delg213) - 0.125*delg312)*NN3) + 
      ginv33*(0.0625*((-delg133 + delg313)*NN1 + 
            (-delg233 + delg323)*NN2) + 0.*delg333*NN3) + 
      ginv11*(0.*delg111*NN1 + 
         0.0625*((delg112 - delg211)*NN2 + (delg113 - delg311)*NN3)) + 
      ginv13*((0.0625*delg123 - 0.125*delg213)*NN2 + 
         0.0625*((-delg113 + delg311)*NN1 + delg312*NN2 + 
            (delg133 - delg313)*NN3)) + 
      ginv22*(0.*delg222*NN2 + 
         0.0625*((-delg122 + delg212)*NN1 + (delg223 - delg322)*NN3)) + 
      ginv23*((-0.125*delg123 + 0.0625*(delg213 + delg312))*NN1 + 
         0.0625*(-(delg223*NN2) + delg322*NN2 + delg233*NN3 - delg323*NN3))\
)*pow2(R)*sqrt(det))/PI
;



} endfor_ijk; /* loop i, j, k */

}  /* function */

/* adm_mass_integrand.c */
/* nvars = 19, nauxs = 59, n* = 209,  n/ = 52,  n+ = 247, n = 508, O = 0 */
