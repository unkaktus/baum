/* dphi_int2_psi4_N.c */
/* Copyright (C) 1998 Bernd Bruegmann, 12.11.2012 */
/* Produced with Mathematica */

#include "bam.h"
#include "Invariants.h"

#define Power(x,y) pow((double) (x), (double) (y))
#define Sqrt(x)    sqrt((double) (x))
#define Log(x)     log((double) (x))
#define pow2(x)    ((x)*(x))
#define pow4(x)    ((x)*(x)*(x)*(x))
#define pow2inv(x) (1.0/((x)*(x)))

#define Tanh(x)    tanh(x)
#define Sech(x)    (1/cosh(x))
#define Cal(x,y,z) ((x)?(y):(z))




void dphi_int2_psi4_N(tVarList *u)
{

double *I2rpsi4 = vldataptr(u, 0);
double *I2ipsi4 = vldataptr(u, 1);
double *dphiI2rpsi4 = vldataptr(u, 2);
double *dphiI2ipsi4 = vldataptr(u, 3);
tL *level = u->level;

double oosqrt2 = 1.0/sqrt(2);
double sqrt2 = sqrt(2);

const int order_centered = Geti("order_centered");
const double x0 = Getd("sphere_x0");
const double y0 = Getd("sphere_y0");
const double z0 = Getd("sphere_z0");
if (!(x0==y0==z0==0.)) errorexit("need sphere center at 0");

const double *xp = level->v[Ind("x")];
const double *yp = level->v[Ind("y")];
const double *zp = level->v[Ind("z")];
const int useShellsTransfo = level->shells;
const double *rp = level->v[IndLax("shells_R")];
const double *rr = level->v[IndLax("shells_r")];
const double shellsS = GetdLax("amr_shells_stretch");
const double shellsR = GetdLax("amr_shells_r0");
const double shellsE = GetdLax("amr_shells_eps");




bampi_openmp_start


double dx = level->dx;
double dy = level->dy;
double dz = level->dz;
double oo2dx = 1/(2*dx), oo2dy = 1/(2*dy), oo2dz = 1/(2*dz);
double oodx2 = 1/(dx*dx), oody2 = 1/(dy*dy), oodz2 = 1/(dz*dz);
double oo4dxdy = 1/(4*dx*dy);
double oo4dxdz = 1/(4*dx*dz);
double oo4dydz = 1/(4*dy*dz);

double cosphi = 0.;
double costheta = 0.;
double ddRdr = 0.;
double dI2ipsi41 = 0.;
double dI2ipsi42 = 0.;
double dI2ipsi43 = 0.;
double dI2ipsi4SST1 = 0.;
double dI2ipsi4SST2 = 0.;
double dI2ipsi4SST3 = 0.;
double dI2rpsi41 = 0.;
double dI2rpsi42 = 0.;
double dI2rpsi43 = 0.;
double dI2rpsi4SST1 = 0.;
double dI2rpsi4SST2 = 0.;
double dI2rpsi4SST3 = 0.;
double DJac111 = 0.;
double DJac112 = 0.;
double DJac113 = 0.;
double DJac122 = 0.;
double DJac123 = 0.;
double DJac133 = 0.;
double DJac211 = 0.;
double DJac212 = 0.;
double DJac213 = 0.;
double DJac222 = 0.;
double DJac223 = 0.;
double DJac233 = 0.;
double DJac311 = 0.;
double DJac312 = 0.;
double DJac313 = 0.;
double DJac322 = 0.;
double DJac323 = 0.;
double DJac333 = 0.;
double dRdr = 0.;
double dxyzdphi1 = 0.;
double dxyzdphi2 = 0.;
double dxyzdphi3 = 0.;
double Jac11 = 0.;
double Jac12 = 0.;
double Jac13 = 0.;
double Jac21 = 0.;
double Jac22 = 0.;
double Jac23 = 0.;
double Jac31 = 0.;
double Jac32 = 0.;
double Jac33 = 0.;
double r = 0.;
double sinphi = 0.;
double sintheta = 0.;
double x = 0.;
double y = 0.;
double z = 0.;



forinnerpoints_ijk_openmp(level) {


x
=
-x0 + xp[ijk]
;

y
=
-y0 + yp[ijk]
;

z
=
-z0 + zp[ijk]
;

r
=
Sqrt(pow2(x) + pow2(y) + pow2(z))
;


if (order_centered == 2 || boundary1away) { 

dI2rpsi41
=
oo2dx*(-I2rpsi4[-di + ijk] + I2rpsi4[di + ijk])
;

dI2rpsi42
=
oo2dy*(-I2rpsi4[-dj + ijk] + I2rpsi4[dj + ijk])
;

dI2rpsi43
=
oo2dz*(-I2rpsi4[-dk + ijk] + I2rpsi4[dk + ijk])
;

dI2ipsi41
=
oo2dx*(-I2ipsi4[-di + ijk] + I2ipsi4[di + ijk])
;

dI2ipsi42
=
oo2dy*(-I2ipsi4[-dj + ijk] + I2ipsi4[dj + ijk])
;

dI2ipsi43
=
oo2dz*(-I2ipsi4[-dk + ijk] + I2ipsi4[dk + ijk])
;


} else if (order_centered == 4 || boundaryNaway(2)) { 

dI2rpsi41
=
0.16666666666666666667*oo2dx*(I2rpsi4[-2*di + ijk] + 
    8.*(-I2rpsi4[-di + ijk] + I2rpsi4[di + ijk]) - I2rpsi4[2*di + ijk])
;

dI2rpsi42
=
0.16666666666666666667*oo2dy*(I2rpsi4[-2*dj + ijk] + 
    8.*(-I2rpsi4[-dj + ijk] + I2rpsi4[dj + ijk]) - I2rpsi4[2*dj + ijk])
;

dI2rpsi43
=
0.16666666666666666667*oo2dz*(I2rpsi4[-2*dk + ijk] + 
    8.*(-I2rpsi4[-dk + ijk] + I2rpsi4[dk + ijk]) - I2rpsi4[2*dk + ijk])
;

dI2ipsi41
=
0.16666666666666666667*oo2dx*(I2ipsi4[-2*di + ijk] + 
    8.*(-I2ipsi4[-di + ijk] + I2ipsi4[di + ijk]) - I2ipsi4[2*di + ijk])
;

dI2ipsi42
=
0.16666666666666666667*oo2dy*(I2ipsi4[-2*dj + ijk] + 
    8.*(-I2ipsi4[-dj + ijk] + I2ipsi4[dj + ijk]) - I2ipsi4[2*dj + ijk])
;

dI2ipsi43
=
0.16666666666666666667*oo2dz*(I2ipsi4[-2*dk + ijk] + 
    8.*(-I2ipsi4[-dk + ijk] + I2ipsi4[dk + ijk]) - I2ipsi4[2*dk + ijk])
;


} else if (order_centered == 6 || boundaryNaway(3)) { 

dI2rpsi41
=
0.033333333333333333333*oo2dx*(-I2rpsi4[-3*di + ijk] + 
    45.*(-I2rpsi4[-di + ijk] + I2rpsi4[di + ijk]) + 
    9.*(I2rpsi4[-2*di + ijk] - I2rpsi4[2*di + ijk]) + I2rpsi4[3*di + ijk])
;

dI2rpsi42
=
0.033333333333333333333*oo2dy*(-I2rpsi4[-3*dj + ijk] + 
    45.*(-I2rpsi4[-dj + ijk] + I2rpsi4[dj + ijk]) + 
    9.*(I2rpsi4[-2*dj + ijk] - I2rpsi4[2*dj + ijk]) + I2rpsi4[3*dj + ijk])
;

dI2rpsi43
=
0.033333333333333333333*oo2dz*(-I2rpsi4[-3*dk + ijk] + 
    45.*(-I2rpsi4[-dk + ijk] + I2rpsi4[dk + ijk]) + 
    9.*(I2rpsi4[-2*dk + ijk] - I2rpsi4[2*dk + ijk]) + I2rpsi4[3*dk + ijk])
;

dI2ipsi41
=
0.033333333333333333333*oo2dx*(-I2ipsi4[-3*di + ijk] + 
    45.*(-I2ipsi4[-di + ijk] + I2ipsi4[di + ijk]) + 
    9.*(I2ipsi4[-2*di + ijk] - I2ipsi4[2*di + ijk]) + I2ipsi4[3*di + ijk])
;

dI2ipsi42
=
0.033333333333333333333*oo2dy*(-I2ipsi4[-3*dj + ijk] + 
    45.*(-I2ipsi4[-dj + ijk] + I2ipsi4[dj + ijk]) + 
    9.*(I2ipsi4[-2*dj + ijk] - I2ipsi4[2*dj + ijk]) + I2ipsi4[3*dj + ijk])
;

dI2ipsi43
=
0.033333333333333333333*oo2dz*(-I2ipsi4[-3*dk + ijk] + 
    45.*(-I2ipsi4[-dk + ijk] + I2ipsi4[dk + ijk]) + 
    9.*(I2ipsi4[-2*dk + ijk] - I2ipsi4[2*dk + ijk]) + I2ipsi4[3*dk + ijk])
;


} else errorexit("order is not implemented yet"); 


if (useShellsTransfo) { 

r
=
rp[ijk]
;

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

DJac111
=
(pow2(rp[ijk]) - pow2(xp[ijk]) - ddRdr*rp[ijk]*pow2(xp[ijk])*pow2inv(dRdr))/
  (dRdr*Power(rp[ijk],3))
;

DJac112
=
-((xp[ijk]*yp[ijk]*(ddRdr*rp[ijk] + pow2(dRdr)))/
    (Power(dRdr,3)*Power(rp[ijk],3)))
;

DJac113
=
-((xp[ijk]*zp[ijk]*(ddRdr*rp[ijk] + pow2(dRdr)))/
    (Power(dRdr,3)*Power(rp[ijk],3)))
;

DJac122
=
(pow2(rp[ijk]) - pow2(yp[ijk]) - ddRdr*rp[ijk]*pow2(yp[ijk])*pow2inv(dRdr))/
  (dRdr*Power(rp[ijk],3))
;

DJac123
=
-((yp[ijk]*zp[ijk]*(ddRdr*rp[ijk] + pow2(dRdr)))/
    (Power(dRdr,3)*Power(rp[ijk],3)))
;

DJac133
=
(pow2(rp[ijk]) - pow2(zp[ijk]) - ddRdr*rp[ijk]*pow2(zp[ijk])*pow2inv(dRdr))/
  (dRdr*Power(rp[ijk],3))
;

DJac211
=
2.*xp[ijk]*yp[ijk]*pow2inv(pow2(xp[ijk]) + pow2(yp[ijk]))
;

DJac212
=
(-pow2(xp[ijk]) + pow2(yp[ijk]))*pow2inv(pow2(xp[ijk]) + pow2(yp[ijk]))
;

DJac213
=
0
;

DJac222
=
-2.*xp[ijk]*yp[ijk]*pow2inv(pow2(xp[ijk]) + pow2(yp[ijk]))
;

DJac223
=
0
;

DJac233
=
0
;

DJac311
=
2.*xp[ijk]*zp[ijk]*pow2inv(pow2(xp[ijk]) + pow2(zp[ijk]))
;

DJac312
=
0
;

DJac313
=
(-pow2(xp[ijk]) + pow2(zp[ijk]))*pow2inv(pow2(xp[ijk]) + pow2(zp[ijk]))
;

DJac322
=
0
;

DJac323
=
0
;

DJac333
=
-2.*xp[ijk]*zp[ijk]*pow2inv(pow2(xp[ijk]) + pow2(zp[ijk]))
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

DJac111
=
(pow2(rp[ijk]) - pow2(xp[ijk]) - ddRdr*rp[ijk]*pow2(xp[ijk])*pow2inv(dRdr))/
  (dRdr*Power(rp[ijk],3))
;

DJac112
=
-((xp[ijk]*yp[ijk]*(ddRdr*rp[ijk] + pow2(dRdr)))/
    (Power(dRdr,3)*Power(rp[ijk],3)))
;

DJac113
=
-((xp[ijk]*zp[ijk]*(ddRdr*rp[ijk] + pow2(dRdr)))/
    (Power(dRdr,3)*Power(rp[ijk],3)))
;

DJac122
=
(pow2(rp[ijk]) - pow2(yp[ijk]) - ddRdr*rp[ijk]*pow2(yp[ijk])*pow2inv(dRdr))/
  (dRdr*Power(rp[ijk],3))
;

DJac123
=
-((yp[ijk]*zp[ijk]*(ddRdr*rp[ijk] + pow2(dRdr)))/
    (Power(dRdr,3)*Power(rp[ijk],3)))
;

DJac133
=
(pow2(rp[ijk]) - pow2(zp[ijk]) - ddRdr*rp[ijk]*pow2(zp[ijk])*pow2inv(dRdr))/
  (dRdr*Power(rp[ijk],3))
;

DJac211
=
-2.*xp[ijk]*yp[ijk]*pow2inv(pow2(xp[ijk]) + pow2(yp[ijk]))
;

DJac212
=
(pow2(xp[ijk]) - pow2(yp[ijk]))*pow2inv(pow2(xp[ijk]) + pow2(yp[ijk]))
;

DJac213
=
0
;

DJac222
=
2.*xp[ijk]*yp[ijk]*pow2inv(pow2(xp[ijk]) + pow2(yp[ijk]))
;

DJac223
=
0
;

DJac233
=
0
;

DJac311
=
0
;

DJac312
=
0
;

DJac313
=
0
;

DJac322
=
2.*yp[ijk]*zp[ijk]*pow2inv(pow2(yp[ijk]) + pow2(zp[ijk]))
;

DJac323
=
(-pow2(yp[ijk]) + pow2(zp[ijk]))*pow2inv(pow2(yp[ijk]) + pow2(zp[ijk]))
;

DJac333
=
-2.*yp[ijk]*zp[ijk]*pow2inv(pow2(yp[ijk]) + pow2(zp[ijk]))
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

DJac111
=
(pow2(rp[ijk]) - pow2(xp[ijk]) - ddRdr*rp[ijk]*pow2(xp[ijk])*pow2inv(dRdr))/
  (dRdr*Power(rp[ijk],3))
;

DJac112
=
-((xp[ijk]*yp[ijk]*(ddRdr*rp[ijk] + pow2(dRdr)))/
    (Power(dRdr,3)*Power(rp[ijk],3)))
;

DJac113
=
-((xp[ijk]*zp[ijk]*(ddRdr*rp[ijk] + pow2(dRdr)))/
    (Power(dRdr,3)*Power(rp[ijk],3)))
;

DJac122
=
(pow2(rp[ijk]) - pow2(yp[ijk]) - ddRdr*rp[ijk]*pow2(yp[ijk])*pow2inv(dRdr))/
  (dRdr*Power(rp[ijk],3))
;

DJac123
=
-((yp[ijk]*zp[ijk]*(ddRdr*rp[ijk] + pow2(dRdr)))/
    (Power(dRdr,3)*Power(rp[ijk],3)))
;

DJac133
=
(pow2(rp[ijk]) - pow2(zp[ijk]) - ddRdr*rp[ijk]*pow2(zp[ijk])*pow2inv(dRdr))/
  (dRdr*Power(rp[ijk],3))
;

DJac211
=
-2.*xp[ijk]*zp[ijk]*pow2inv(pow2(xp[ijk]) + pow2(zp[ijk]))
;

DJac212
=
0
;

DJac213
=
(pow2(xp[ijk]) - pow2(zp[ijk]))*pow2inv(pow2(xp[ijk]) + pow2(zp[ijk]))
;

DJac222
=
0
;

DJac223
=
0
;

DJac233
=
2.*xp[ijk]*zp[ijk]*pow2inv(pow2(xp[ijk]) + pow2(zp[ijk]))
;

DJac311
=
0
;

DJac312
=
0
;

DJac313
=
0
;

DJac322
=
-2.*yp[ijk]*zp[ijk]*pow2inv(pow2(yp[ijk]) + pow2(zp[ijk]))
;

DJac323
=
(pow2(yp[ijk]) - pow2(zp[ijk]))*pow2inv(pow2(yp[ijk]) + pow2(zp[ijk]))
;

DJac333
=
2.*yp[ijk]*zp[ijk]*pow2inv(pow2(yp[ijk]) + pow2(zp[ijk]))
;


} 

dI2rpsi4SST1
=
dI2rpsi41*Jac11 + dI2rpsi42*Jac21 + dI2rpsi43*Jac31
;

dI2rpsi4SST2
=
dI2rpsi41*Jac12 + dI2rpsi42*Jac22 + dI2rpsi43*Jac32
;

dI2rpsi4SST3
=
dI2rpsi41*Jac13 + dI2rpsi42*Jac23 + dI2rpsi43*Jac33
;

dI2ipsi4SST1
=
dI2ipsi41*Jac11 + dI2ipsi42*Jac21 + dI2ipsi43*Jac31
;

dI2ipsi4SST2
=
dI2ipsi41*Jac12 + dI2ipsi42*Jac22 + dI2ipsi43*Jac32
;

dI2ipsi4SST3
=
dI2ipsi41*Jac13 + dI2ipsi42*Jac23 + dI2ipsi43*Jac33
;

dI2rpsi41
=
dI2rpsi4SST1
;

dI2rpsi42
=
dI2rpsi4SST2
;

dI2rpsi43
=
dI2rpsi4SST3
;

dI2ipsi41
=
dI2ipsi4SST1
;

dI2ipsi42
=
dI2ipsi4SST2
;

dI2ipsi43
=
dI2ipsi4SST3
;


} 

costheta
=
z/r
;

sintheta
=
Sqrt(1. - pow2(z)/Power(r,2))
;

cosphi
=
x/(r*sintheta)
;

sinphi
=
y/(r*sintheta)
;

dxyzdphi1
=
0.
;

dxyzdphi2
=
0.
;

dxyzdphi3
=
0.
;


dxyzdphi1 = -r*sinphi*sintheta; 


dxyzdphi2 =  r*cosphi*sintheta; 

dphiI2rpsi4[ijk]
=
dI2rpsi41*dxyzdphi1 + dI2rpsi42*dxyzdphi2 + dI2rpsi43*dxyzdphi3
;

dphiI2ipsi4[ijk]
=
dI2ipsi41*dxyzdphi1 + dI2ipsi42*dxyzdphi2 + dI2ipsi43*dxyzdphi3
;



} endfor_ijk_openmp; /* loop i, j, k */



bampi_openmp_stop


}  /* function */

/* dphi_int2_psi4_N.c */
/* nvars = 9, nauxs = 52, n* = 283,  n/ = 70,  n+ = 324, n = 677, O = 0 */
