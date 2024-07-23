/* adm_pj_integrand.c */
/* Copyright (C) 1998 Bernd Bruegmann, 16.9.2011 */
/* Produced with Mathematica */

#include "bam.h"
#include "ADM_mass.h"

#define Power(x,y) (pow((double) (x), (double) (y)))
#define pow2(x)    ((x)*(x))
#define Sqrt(x)    sqrt(x)
#define PI  3.1415926535897932




void adm_pj_integrand(tL *level, int i_x, int i_g, int i_K, int i_TrK)
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
double *K11 = level->v[i_K+0];
double *K12 = level->v[i_K+1];
double *K13 = level->v[i_K+2];
double *K22 = level->v[i_K+3];
double *K23 = level->v[i_K+4];
double *K33 = level->v[i_K+5];
double *TrK = level->v[i_TrK+0];
double *IP1 = level->v[Ind("ADM_mass_Pxint")];
double *IP2 = level->v[Ind("ADM_mass_Pyint")];
double *IP3 = level->v[Ind("ADM_mass_Pzint")];
double *IJ1 = level->v[Ind("ADM_mass_Jxint")];
double *IJ2 = level->v[Ind("ADM_mass_Jyint")];
double *IJ3 = level->v[Ind("ADM_mass_Jzint")];

double dx = level->dx;
double dy = level->dy;
double dz = level->dz;
double oo2dx = 1/(2*dx), oo2dy = 1/(2*dy), oo2dz = 1/(2*dz);
double oodx2 = 1/(dx*dx), oody2 = 1/(dy*dy), oodz2 = 1/(dz*dz);
double oo4dxdy = 1/(4*dx*dy);
double oo4dxdz = 1/(4*dx*dz);
double oo4dydz = 1/(4*dy*dz);
long ccc;

double det;
double detginvf;
double ginv11;
double ginv12;
double ginv13;
double ginv22;
double ginv23;
double ginv33;
double NN1;
double NN2;
double NN3;
double R;



forinner19(level) {


det
=
2.*g12[ccc]*g13[ccc]*g23[ccc] + g11[ccc]*g22[ccc]*g33[ccc] - 
  g33[ccc]*pow2(g12[ccc]) - g22[ccc]*pow2(g13[ccc]) - g11[ccc]*pow2(g23[ccc])
;

detginvf
=
1/(2.*g12[ccc]*g13[ccc]*g23[ccc] + g11[ccc]*g22[ccc]*g33[ccc] - 
    g33[ccc]*pow2(g12[ccc]) - g22[ccc]*pow2(g13[ccc]) - 
    g11[ccc]*pow2(g23[ccc]))
;

ginv11
=
detginvf*(g22[ccc]*g33[ccc] - pow2(g23[ccc]))
;

ginv12
=
detginvf*(g13[ccc]*g23[ccc] - g12[ccc]*g33[ccc])
;

ginv13
=
detginvf*(-(g13[ccc]*g22[ccc]) + g12[ccc]*g23[ccc])
;

ginv22
=
detginvf*(g11[ccc]*g33[ccc] - pow2(g13[ccc]))
;

ginv23
=
detginvf*(g12[ccc]*g13[ccc] - g11[ccc]*g23[ccc])
;

ginv33
=
detginvf*(g11[ccc]*g22[ccc] - pow2(g12[ccc]))
;

R
=
sqrt(pow2(x1[ccc]) + pow2(x2[ccc]) + pow2(x3[ccc]))
;

NN1
=
x1[ccc]/R
;

NN2
=
x2[ccc]/R
;

NN3
=
x3[ccc]/R
;

IP1[ccc]
=
(0.125*(ginv22*NN2*K12[ccc] + ginv12*(NN2*K11[ccc] + NN1*K12[ccc]) + 
      ginv33*NN3*K13[ccc] + ginv13*(NN3*K11[ccc] + NN1*K13[ccc]) + 
      ginv23*(NN3*K12[ccc] + NN2*K13[ccc]) + 
      NN1*(ginv11*K11[ccc] - TrK[ccc]))*pow2(R)*sqrt(det))/PI
;

IP2[ccc]
=
(0.125*(ginv11*NN1*K12[ccc] + ginv12*(NN2*K12[ccc] + NN1*K22[ccc]) + 
      ginv33*NN3*K23[ccc] + ginv13*(NN3*K12[ccc] + NN1*K23[ccc]) + 
      ginv23*(NN3*K22[ccc] + NN2*K23[ccc]) + 
      NN2*(ginv22*K22[ccc] - TrK[ccc]))*pow2(R)*sqrt(det))/PI
;

IP3[ccc]
=
(0.125*(ginv11*NN1*K13[ccc] + ginv22*NN2*K23[ccc] + 
      ginv12*(NN2*K13[ccc] + NN1*K23[ccc]) + 
      ginv13*(NN3*K13[ccc] + NN1*K33[ccc]) + 
      ginv23*(NN3*K23[ccc] + NN2*K33[ccc]) + 
      NN3*(ginv33*K33[ccc] - TrK[ccc]))*pow2(R)*sqrt(det))/PI
;

IJ1[ccc]
=
(0.125*(TrK[ccc]*(-(NN3*x2[ccc]) + NN2*x3[ccc]) + 
      ginv11*NN1*(K13[ccc]*x2[ccc] - K12[ccc]*x3[ccc]) + 
      ginv22*NN2*(K23[ccc]*x2[ccc] - K22[ccc]*x3[ccc]) + 
      ginv12*((NN2*K13[ccc] + NN1*K23[ccc])*x2[ccc] - 
         (NN2*K12[ccc] + NN1*K22[ccc])*x3[ccc]) + 
      ginv33*NN3*(K33[ccc]*x2[ccc] - K23[ccc]*x3[ccc]) + 
      ginv13*((NN3*K13[ccc] + NN1*K33[ccc])*x2[ccc] - 
         (NN3*K12[ccc] + NN1*K23[ccc])*x3[ccc]) + 
      ginv23*((NN3*K23[ccc] + NN2*K33[ccc])*x2[ccc] - 
         (NN3*K22[ccc] + NN2*K23[ccc])*x3[ccc]))*pow2(R)*sqrt(det))/PI
;

IJ2[ccc]
=
(0.125*(TrK[ccc]*(NN3*x1[ccc] - NN1*x3[ccc]) + 
      ginv11*NN1*(-(K13[ccc]*x1[ccc]) + K11[ccc]*x3[ccc]) + 
      ginv22*NN2*(-(K23[ccc]*x1[ccc]) + K12[ccc]*x3[ccc]) + 
      ginv12*(-((NN2*K13[ccc] + NN1*K23[ccc])*x1[ccc]) + 
         (NN2*K11[ccc] + NN1*K12[ccc])*x3[ccc]) + 
      ginv33*NN3*(-(K33[ccc]*x1[ccc]) + K13[ccc]*x3[ccc]) + 
      ginv13*(-((NN3*K13[ccc] + NN1*K33[ccc])*x1[ccc]) + 
         (NN3*K11[ccc] + NN1*K13[ccc])*x3[ccc]) + 
      ginv23*(-((NN3*K23[ccc] + NN2*K33[ccc])*x1[ccc]) + 
         (NN3*K12[ccc] + NN2*K13[ccc])*x3[ccc]))*pow2(R)*sqrt(det))/PI
;

IJ3[ccc]
=
(0.125*(TrK[ccc]*(-(NN2*x1[ccc]) + NN1*x2[ccc]) + 
      ginv11*NN1*(K12[ccc]*x1[ccc] - K11[ccc]*x2[ccc]) + 
      ginv22*NN2*(K22[ccc]*x1[ccc] - K12[ccc]*x2[ccc]) + 
      ginv12*((NN2*K12[ccc] + NN1*K22[ccc])*x1[ccc] - 
         (NN2*K11[ccc] + NN1*K12[ccc])*x2[ccc]) + 
      ginv33*NN3*(K23[ccc]*x1[ccc] - K13[ccc]*x2[ccc]) + 
      ginv13*((NN3*K12[ccc] + NN1*K23[ccc])*x1[ccc] - 
         (NN3*K11[ccc] + NN1*K13[ccc])*x2[ccc]) + 
      ginv23*((NN3*K22[ccc] + NN2*K23[ccc])*x1[ccc] - 
         (NN3*K12[ccc] + NN2*K13[ccc])*x2[ccc]))*pow2(R)*sqrt(det))/PI
;



} endforinner; /* loop i, j, k */

}  /* function */

/* adm_pj_integrand.c */
/* nvars = 22, nauxs = 12, n* = 248,  n/ = 29,  n+ = 150, n = 427, O = 0 */
