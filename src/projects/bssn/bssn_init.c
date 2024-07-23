/* bssn_init.c */
/* Copyright (C) 1998 Bernd Bruegmann, 24.1.2024 */
/* Produced with Mathematica */

#include "bam.h"
#include "bssn.h"

#define Power(x,y) (pow((double) (x), (double) (y)))
#define Sqrt(x)    sqrt((double) (x))
#define Log(x)     log((double) (x))
#define pow2(x)    ((x)*(x))
#define pow2inv(x) (1.0/((x)*(x)))

#define Tanh(x)    tanh(x)
#define Sech(x)    (1/cosh(x))
#define Cal(x,y,z) ((x)?(y):(z))




void bssn_init(tVarList *ucur)
{

tL *level = ucur->level;
double *gb11 = vldataptr(ucur, 0);
double *gb12 = vldataptr(ucur, 1);
double *gb13 = vldataptr(ucur, 2);
double *gb22 = vldataptr(ucur, 3);
double *gb23 = vldataptr(ucur, 4);
double *gb33 = vldataptr(ucur, 5);
double *Kb11 = vldataptr(ucur, 6);
double *Kb12 = vldataptr(ucur, 7);
double *Kb13 = vldataptr(ucur, 8);
double *Kb22 = vldataptr(ucur, 9);
double *Kb23 = vldataptr(ucur, 10);
double *Kb33 = vldataptr(ucur, 11);
double *gt11 = vldataptr(ucur, 12);
double *gt12 = vldataptr(ucur, 13);
double *gt13 = vldataptr(ucur, 14);
double *gt22 = vldataptr(ucur, 15);
double *gt23 = vldataptr(ucur, 16);
double *gt33 = vldataptr(ucur, 17);
double *chi = vldataptr(ucur, 18);
double *At11 = vldataptr(ucur, 19);
double *At12 = vldataptr(ucur, 20);
double *At13 = vldataptr(ucur, 21);
double *At22 = vldataptr(ucur, 22);
double *At23 = vldataptr(ucur, 23);
double *At33 = vldataptr(ucur, 24);
double *K = vldataptr(ucur, 25);
double *G1 = vldataptr(ucur, 26);
double *G2 = vldataptr(ucur, 27);
double *G3 = vldataptr(ucur, 28);
double *gtinv11 = vldataptr(ucur, 29);
double *gtinv12 = vldataptr(ucur, 30);
double *gtinv13 = vldataptr(ucur, 31);
double *gtinv22 = vldataptr(ucur, 32);
double *gtinv23 = vldataptr(ucur, 33);
double *gtinv33 = vldataptr(ucur, 34);
double *shiftDr = vldataptr(ucur, 35);

double chipsipower  = Getd("bssn_chi_psipower");

const double *xp = level->v[Ind("x")];
const double *yp = level->v[Ind("y")];
const double *zp = level->v[Ind("z")];
const int useShellsTransfo = level->shells;
const double *rp = level->v[IndLax("shells_R")];
const double *rr = level->v[IndLax("shells_r")];
const double shellsS = GetdLax("amr_shells_stretch");
const double shellsR = GetdLax("amr_shells_r0");
const double shellsE = GetdLax("amr_shells_eps");
const int use_eta   = Getv("bssn_use_eta", "yes");
bssn_eta_init(level->grid);



bampi_openmp_start


double dx = level->dx;
double dy = level->dy;
double dz = level->dz;
double oo2dx = 1/(2*dx), oo2dy = 1/(2*dy), oo2dz = 1/(2*dz);
double oodx2 = 1/(dx*dx), oody2 = 1/(dy*dy), oodz2 = 1/(dz*dz);
double oo4dxdy = 1/(4*dx*dy);
double oo4dxdz = 1/(4*dx*dz);
double oo4dydz = 1/(4*dy*dz);

double A11 = 0.;
double A12 = 0.;
double A13 = 0.;
double A21 = 0.;
double A22 = 0.;
double A23 = 0.;
double A31 = 0.;
double A32 = 0.;
double A33 = 0.;
double absdpsim2 = 0.;
double dchi1 = 0.;
double dchi2 = 0.;
double dchi3 = 0.;
double dchiSST1 = 0.;
double dchiSST2 = 0.;
double dchiSST3 = 0.;
double ddRdr = 0.;
double detgb = 0.;
double detgt = 0.;
double detgtinv = 0.;
double dgtinv111 = 0.;
double dgtinv112 = 0.;
double dgtinv113 = 0.;
double dgtinv122 = 0.;
double dgtinv123 = 0.;
double dgtinv133 = 0.;
double dgtinv211 = 0.;
double dgtinv212 = 0.;
double dgtinv213 = 0.;
double dgtinv222 = 0.;
double dgtinv223 = 0.;
double dgtinv233 = 0.;
double dgtinv311 = 0.;
double dgtinv312 = 0.;
double dgtinv313 = 0.;
double dgtinv322 = 0.;
double dgtinv323 = 0.;
double dgtinv333 = 0.;
double dgtinvSST111 = 0.;
double dgtinvSST112 = 0.;
double dgtinvSST113 = 0.;
double dgtinvSST121 = 0.;
double dgtinvSST122 = 0.;
double dgtinvSST123 = 0.;
double dgtinvSST131 = 0.;
double dgtinvSST132 = 0.;
double dgtinvSST133 = 0.;
double dgtinvSST211 = 0.;
double dgtinvSST212 = 0.;
double dgtinvSST213 = 0.;
double dgtinvSST221 = 0.;
double dgtinvSST222 = 0.;
double dgtinvSST223 = 0.;
double dgtinvSST231 = 0.;
double dgtinvSST232 = 0.;
double dgtinvSST233 = 0.;
double dgtinvSST311 = 0.;
double dgtinvSST312 = 0.;
double dgtinvSST313 = 0.;
double dgtinvSST321 = 0.;
double dgtinvSST322 = 0.;
double dgtinvSST323 = 0.;
double dgtinvSST331 = 0.;
double dgtinvSST332 = 0.;
double dgtinvSST333 = 0.;
double dpsim21 = 0.;
double dpsim22 = 0.;
double dpsim23 = 0.;
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
double Kt11 = 0.;
double Kt12 = 0.;
double Kt13 = 0.;
double Kt22 = 0.;
double Kt23 = 0.;
double Kt33 = 0.;
double psi = 0.;
double psim2 = 0.;
double psim4 = 0.;



forallpoints_ijk_openmp(level) {


detgb
=
2.*gb12[ijk]*gb13[ijk]*gb23[ijk] + gb11[ijk]*gb22[ijk]*gb33[ijk] - 
  gb33[ijk]*pow2(gb12[ijk]) - gb22[ijk]*pow2(gb13[ijk]) - 
  gb11[ijk]*pow2(gb23[ijk])
;

psim4
=
Power(detgb,-0.3333333333333333)
;

psi
=
Power(detgb,0.08333333333333333)
;

gt11[ijk]
=
psim4*gb11[ijk]
;

gt12[ijk]
=
psim4*gb12[ijk]
;

gt13[ijk]
=
psim4*gb13[ijk]
;

gt22[ijk]
=
psim4*gb22[ijk]
;

gt23[ijk]
=
psim4*gb23[ijk]
;

gt33[ijk]
=
psim4*gb33[ijk]
;

Kt11
=
psim4*Kb11[ijk]
;

Kt12
=
psim4*Kb12[ijk]
;

Kt13
=
psim4*Kb13[ijk]
;

Kt22
=
psim4*Kb22[ijk]
;

Kt23
=
psim4*Kb23[ijk]
;

Kt33
=
psim4*Kb33[ijk]
;

detgt
=
2.*gt12[ijk]*gt13[ijk]*gt23[ijk] + gt11[ijk]*gt22[ijk]*gt33[ijk] - 
  gt33[ijk]*pow2(gt12[ijk]) - gt22[ijk]*pow2(gt13[ijk]) - 
  gt11[ijk]*pow2(gt23[ijk])
;

detgtinv
=
1/detgt
;

gtinv11[ijk]
=
detgtinv*(gt22[ijk]*gt33[ijk] - pow2(gt23[ijk]))
;

gtinv12[ijk]
=
detgtinv*(gt13[ijk]*gt23[ijk] - gt12[ijk]*gt33[ijk])
;

gtinv13[ijk]
=
detgtinv*(-(gt13[ijk]*gt22[ijk]) + gt12[ijk]*gt23[ijk])
;

gtinv22[ijk]
=
detgtinv*(gt11[ijk]*gt33[ijk] - pow2(gt13[ijk]))
;

gtinv23[ijk]
=
detgtinv*(gt12[ijk]*gt13[ijk] - gt11[ijk]*gt23[ijk])
;

gtinv33[ijk]
=
detgtinv*(gt11[ijk]*gt22[ijk] - pow2(gt12[ijk]))
;

K[ijk]
=
Kt11*gtinv11[ijk] + Kt22*gtinv22[ijk] + 
  2.*(Kt12*gtinv12[ijk] + Kt13*gtinv13[ijk] + Kt23*gtinv23[ijk]) + 
  Kt33*gtinv33[ijk]
;

A11
=
-0.33333333333333333333*gb11[ijk]*K[ijk] + Kb11[ijk]
;

A12
=
-0.33333333333333333333*gb12[ijk]*K[ijk] + Kb12[ijk]
;

A13
=
-0.33333333333333333333*gb13[ijk]*K[ijk] + Kb13[ijk]
;

A21
=
-0.33333333333333333333*gb12[ijk]*K[ijk] + Kb12[ijk]
;

A22
=
-0.33333333333333333333*gb22[ijk]*K[ijk] + Kb22[ijk]
;

A23
=
-0.33333333333333333333*gb23[ijk]*K[ijk] + Kb23[ijk]
;

A31
=
-0.33333333333333333333*gb13[ijk]*K[ijk] + Kb13[ijk]
;

A32
=
-0.33333333333333333333*gb23[ijk]*K[ijk] + Kb23[ijk]
;

A33
=
-0.33333333333333333333*gb33[ijk]*K[ijk] + Kb33[ijk]
;

At11[ijk]
=
A11*psim4
;

At12[ijk]
=
A21*psim4
;

At13[ijk]
=
A31*psim4
;

At22[ijk]
=
A22*psim4
;

At23[ijk]
=
A32*psim4
;

At33[ijk]
=
A33*psim4
;

chi[ijk]
=
Power(psi,chipsipower)
;


} endfor_ijk_openmp; 


forinnerpoints_ijk_openmp(level) { 

dgtinv111
=
oo2dx*(-gtinv11[-di + ijk] + gtinv11[di + ijk])
;

dgtinv112
=
oo2dx*(-gtinv12[-di + ijk] + gtinv12[di + ijk])
;

dgtinv113
=
oo2dx*(-gtinv13[-di + ijk] + gtinv13[di + ijk])
;

dgtinv122
=
oo2dx*(-gtinv22[-di + ijk] + gtinv22[di + ijk])
;

dgtinv123
=
oo2dx*(-gtinv23[-di + ijk] + gtinv23[di + ijk])
;

dgtinv133
=
oo2dx*(-gtinv33[-di + ijk] + gtinv33[di + ijk])
;

dgtinv211
=
oo2dy*(-gtinv11[-dj + ijk] + gtinv11[dj + ijk])
;

dgtinv212
=
oo2dy*(-gtinv12[-dj + ijk] + gtinv12[dj + ijk])
;

dgtinv213
=
oo2dy*(-gtinv13[-dj + ijk] + gtinv13[dj + ijk])
;

dgtinv222
=
oo2dy*(-gtinv22[-dj + ijk] + gtinv22[dj + ijk])
;

dgtinv223
=
oo2dy*(-gtinv23[-dj + ijk] + gtinv23[dj + ijk])
;

dgtinv233
=
oo2dy*(-gtinv33[-dj + ijk] + gtinv33[dj + ijk])
;

dgtinv311
=
oo2dz*(-gtinv11[-dk + ijk] + gtinv11[dk + ijk])
;

dgtinv312
=
oo2dz*(-gtinv12[-dk + ijk] + gtinv12[dk + ijk])
;

dgtinv313
=
oo2dz*(-gtinv13[-dk + ijk] + gtinv13[dk + ijk])
;

dgtinv322
=
oo2dz*(-gtinv22[-dk + ijk] + gtinv22[dk + ijk])
;

dgtinv323
=
oo2dz*(-gtinv23[-dk + ijk] + gtinv23[dk + ijk])
;

dgtinv333
=
oo2dz*(-gtinv33[-dk + ijk] + gtinv33[dk + ijk])
;

dchi1
=
oo2dx*(-chi[-di + ijk] + chi[di + ijk])
;

dchi2
=
oo2dy*(-chi[-dj + ijk] + chi[dj + ijk])
;

dchi3
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

dgtinvSST111
=
dgtinv111*Jac11 + dgtinv211*Jac21 + dgtinv311*Jac31
;

dgtinvSST112
=
dgtinv112*Jac11 + dgtinv212*Jac21 + dgtinv312*Jac31
;

dgtinvSST113
=
dgtinv113*Jac11 + dgtinv213*Jac21 + dgtinv313*Jac31
;

dgtinvSST121
=
dgtinv112*Jac11 + dgtinv212*Jac21 + dgtinv312*Jac31
;

dgtinvSST122
=
dgtinv122*Jac11 + dgtinv222*Jac21 + dgtinv322*Jac31
;

dgtinvSST123
=
dgtinv123*Jac11 + dgtinv223*Jac21 + dgtinv323*Jac31
;

dgtinvSST131
=
dgtinv113*Jac11 + dgtinv213*Jac21 + dgtinv313*Jac31
;

dgtinvSST132
=
dgtinv123*Jac11 + dgtinv223*Jac21 + dgtinv323*Jac31
;

dgtinvSST133
=
dgtinv133*Jac11 + dgtinv233*Jac21 + dgtinv333*Jac31
;

dgtinvSST211
=
dgtinv111*Jac12 + dgtinv211*Jac22 + dgtinv311*Jac32
;

dgtinvSST212
=
dgtinv112*Jac12 + dgtinv212*Jac22 + dgtinv312*Jac32
;

dgtinvSST213
=
dgtinv113*Jac12 + dgtinv213*Jac22 + dgtinv313*Jac32
;

dgtinvSST221
=
dgtinv112*Jac12 + dgtinv212*Jac22 + dgtinv312*Jac32
;

dgtinvSST222
=
dgtinv122*Jac12 + dgtinv222*Jac22 + dgtinv322*Jac32
;

dgtinvSST223
=
dgtinv123*Jac12 + dgtinv223*Jac22 + dgtinv323*Jac32
;

dgtinvSST231
=
dgtinv113*Jac12 + dgtinv213*Jac22 + dgtinv313*Jac32
;

dgtinvSST232
=
dgtinv123*Jac12 + dgtinv223*Jac22 + dgtinv323*Jac32
;

dgtinvSST233
=
dgtinv133*Jac12 + dgtinv233*Jac22 + dgtinv333*Jac32
;

dgtinvSST311
=
dgtinv111*Jac13 + dgtinv211*Jac23 + dgtinv311*Jac33
;

dgtinvSST312
=
dgtinv112*Jac13 + dgtinv212*Jac23 + dgtinv312*Jac33
;

dgtinvSST313
=
dgtinv113*Jac13 + dgtinv213*Jac23 + dgtinv313*Jac33
;

dgtinvSST321
=
dgtinv112*Jac13 + dgtinv212*Jac23 + dgtinv312*Jac33
;

dgtinvSST322
=
dgtinv122*Jac13 + dgtinv222*Jac23 + dgtinv322*Jac33
;

dgtinvSST323
=
dgtinv123*Jac13 + dgtinv223*Jac23 + dgtinv323*Jac33
;

dgtinvSST331
=
dgtinv113*Jac13 + dgtinv213*Jac23 + dgtinv313*Jac33
;

dgtinvSST332
=
dgtinv123*Jac13 + dgtinv223*Jac23 + dgtinv323*Jac33
;

dgtinvSST333
=
dgtinv133*Jac13 + dgtinv233*Jac23 + dgtinv333*Jac33
;

dchiSST1
=
dchi1*Jac11 + dchi2*Jac21 + dchi3*Jac31
;

dchiSST2
=
dchi1*Jac12 + dchi2*Jac22 + dchi3*Jac32
;

dchiSST3
=
dchi1*Jac13 + dchi2*Jac23 + dchi3*Jac33
;

dgtinv111
=
dgtinvSST111
;

dgtinv112
=
dgtinvSST121
;

dgtinv113
=
dgtinvSST131
;

dgtinv122
=
dgtinvSST122
;

dgtinv123
=
dgtinvSST132
;

dgtinv133
=
dgtinvSST133
;

dgtinv211
=
dgtinvSST211
;

dgtinv212
=
dgtinvSST221
;

dgtinv213
=
dgtinvSST231
;

dgtinv222
=
dgtinvSST222
;

dgtinv223
=
dgtinvSST232
;

dgtinv233
=
dgtinvSST233
;

dgtinv311
=
dgtinvSST311
;

dgtinv312
=
dgtinvSST321
;

dgtinv313
=
dgtinvSST331
;

dgtinv322
=
dgtinvSST322
;

dgtinv323
=
dgtinvSST332
;

dgtinv333
=
dgtinvSST333
;

dchi1
=
dchiSST1
;

dchi2
=
dchiSST2
;

dchi3
=
dchiSST3
;


 } 

G1[ijk]
=
-dgtinv111 - dgtinv212 - dgtinv313
;

G2[ijk]
=
-dgtinv112 - dgtinv222 - dgtinv323
;

G3[ijk]
=
-dgtinv113 - dgtinv223 - dgtinv333
;


if (use_eta) { 

psim2
=
pow2inv(psi)
;

dpsim21
=
(-2.*chipsipower*dchi1*psim2)/chi[ijk]
;

dpsim22
=
(-2.*chipsipower*dchi2*psim2)/chi[ijk]
;

dpsim23
=
(-2.*chipsipower*dchi3*psim2)/chi[ijk]
;

absdpsim2
=
Sqrt(2.*(dpsim21*(dpsim22*gtinv12[ijk] + dpsim23*gtinv13[ijk]) + 
      dpsim22*dpsim23*gtinv23[ijk]) + gtinv11[ijk]*pow2(dpsim21) + 
   gtinv22[ijk]*pow2(dpsim22) + gtinv33[ijk]*pow2(dpsim23))
;


shiftDr[ijk] = bssn_eta_set(xp[ijk],yp[ijk],zp[ijk], absdpsim2,psim2); 


} 



} endfor_ijk_openmp; /* loop i, j, k */



bampi_openmp_stop


}  /* function */

/* bssn_init.c */
/* nvars = 41, nauxs = 87, n* = 287,  n/ = 51,  n+ = 265, n = 603, O = 0 */
