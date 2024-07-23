/* puncture_properties_integrand.c */
/* Copyright (C) 1998 Bernd Bruegmann, 7.8.2020 */
/* Produced with Mathematica */

#include "bam.h"
#include "Gauge.h"

#define Power(x,y) (pow((double) (x), (double) (y)))
#define pow2(x)    ((x)*(x))
#define pow3(x)    ((x)*(x)*(x))
#define pow2inv(x) (1.0/((x)*(x)))

#define Sqrt(x)    sqrt(x)
#define PI  3.1415926535897932




void puncture_properties_integrand(tL *level, int i_x, int i_g, int i_K, int i_TrK, double x0,double y0, double z0)
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
double *IM  = level->v[Ind("puncture_properties_Mint")];
double *IP1 = level->v[Ind("puncture_properties_Pxint")];
double *IP2 = level->v[Ind("puncture_properties_Pyint")];
double *IP3 = level->v[Ind("puncture_properties_Pzint")];
double *IS1 = level->v[Ind("puncture_properties_Sxint")];
double *IS2 = level->v[Ind("puncture_properties_Syint")];
double *IS3 = level->v[Ind("puncture_properties_Szint")];

double dx = level->dx;
double dy = level->dy;
double dz = level->dz;
double oo2dx = 1/(2*dx), oo2dy = 1/(2*dy), oo2dz = 1/(2*dz);
double oodx2 = 1/(dx*dx), oody2 = 1/(dy*dy), oodz2 = 1/(dz*dz);
double oo4dxdy = 1/(4*dx*dy);
double oo4dxdz = 1/(4*dx*dz);
double oo4dydz = 1/(4*dy*dz);

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
double det = 0.;
double detginvf = 0.;
double ginv11 = 0.;
double ginv12 = 0.;
double ginv13 = 0.;
double ginv22 = 0.;
double ginv23 = 0.;
double ginv33 = 0.;
double NN1 = 0.;
double NN2 = 0.;
double NN3 = 0.;
double ones1 = 0.;
double ones2 = 0.;
double ones3 = 0.;
double phi11 = 0.;
double phi12 = 0.;
double phi13 = 0.;
double phi21 = 0.;
double phi22 = 0.;
double phi23 = 0.;
double phi31 = 0.;
double phi32 = 0.;
double phi33 = 0.;
double R = 0.;
double rho = 0.;
double xp1 = 0.;
double xp2 = 0.;
double xp3 = 0.;



forinner19(level) {


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
0
;

rho
=
0
;


R =sqrt((x1[ijk]-x0)*(x1[ijk]-x0)+(x2[ijk]-y0)*(x2[ijk]-y0)+(x3[ijk]-z0)*(x3[ijk]-z0)); 


rho =sqrt((x1[ijk]-x0)*(x1[ijk]-x0)+(x2[ijk]-y0)*(x2[ijk]-y0)); 

xp1
=
0
;

xp2
=
0
;

xp3
=
0
;

ones1
=
1.
;

ones2
=
1.
;

ones3
=
1.
;


xp1 = (x1[ijk]-x0); 


xp2 = (x2[ijk]-y0); 


xp3 = (x3[ijk]-z0); 

NN1
=
xp1/R
;

NN2
=
xp2/R
;

NN3
=
xp3/R
;

phi11
=
0
;

phi12
=
xp3
;

phi13
=
-xp2
;

phi21
=
-xp3
;

phi22
=
0
;

phi23
=
xp1
;

phi31
=
xp2
;

phi32
=
-xp1
;

phi33
=
0
;

IM[ijk]
=
0. + (R*(ginv12*(0.0625*((-delg112 + delg211)*xp1 + 
             (delg122 - delg212)*xp2) + 
          (0.0625*(delg123 + delg213) - 0.125*delg312)*xp3) + 
       ginv13*((0.0625*delg123 - 0.125*delg213)*xp2 + 
          0.0625*((-delg113 + delg311)*xp1 + delg312*xp2 + 
             (delg133 - delg313)*xp3)) + 
       0.0625*(ginv33*((-delg133 + delg313)*xp1 + 
             (-delg233 + delg323)*xp2) + 
          ginv11*((delg112 - delg211)*xp2 + (delg113 - delg311)*xp3) + 
          ginv22*((-delg122 + delg212)*xp1 + (delg223 - delg322)*xp3)) + 
       ginv23*((-0.125*delg123 + 0.0625*(delg213 + delg312))*xp1 + 
          0.0625*(-(delg223*xp2) + delg322*xp2 + delg233*xp3 - delg323*xp3)\
))*sqrt(det))/PI
;

IP1[ijk]
=
0. + (0.125*R*(ginv22*xp2*K12[ijk] + 
       ginv12*(xp2*K11[ijk] + xp1*K12[ijk]) + ginv33*xp3*K13[ijk] + 
       ginv13*(xp3*K11[ijk] + xp1*K13[ijk]) + 
       ginv23*(xp3*K12[ijk] + xp2*K13[ijk]) + 
       xp1*(ginv11*K11[ijk] - TrK[ijk]))*sqrt(det))/PI
;

IP2[ijk]
=
0. + (0.125*R*(ginv11*xp1*K12[ijk] + 
       ginv12*(xp2*K12[ijk] + xp1*K22[ijk]) + ginv33*xp3*K23[ijk] + 
       ginv13*(xp3*K12[ijk] + xp1*K23[ijk]) + 
       ginv23*(xp3*K22[ijk] + xp2*K23[ijk]) + 
       xp2*(ginv22*K22[ijk] - TrK[ijk]))*sqrt(det))/PI
;

IP3[ijk]
=
0. + (0.125*R*(ginv11*xp1*K13[ijk] + ginv22*xp2*K23[ijk] + 
       ginv12*(xp2*K13[ijk] + xp1*K23[ijk]) + 
       ginv13*(xp3*K13[ijk] + xp1*K33[ijk]) + 
       ginv23*(xp3*K23[ijk] + xp2*K33[ijk]) + 
       xp3*(ginv33*K33[ijk] - TrK[ijk]))*sqrt(det))/PI
;

IS1[ijk]
=
(0.125*R*(ginv11*xp1*(phi11*K11[ijk] + phi21*K12[ijk] + phi31*K13[ijk]) + 
      ginv22*xp2*(phi11*K12[ijk] + phi21*K22[ijk] + phi31*K23[ijk]) + 
      ginv12*(phi11*(xp2*K11[ijk] + xp1*K12[ijk]) + 
         xp2*(phi21*K12[ijk] + phi31*K13[ijk]) + 
         xp1*(phi21*K22[ijk] + phi31*K23[ijk])) + 
      ginv33*xp3*(phi11*K13[ijk] + phi21*K23[ijk] + phi31*K33[ijk]) + 
      ginv13*(xp3*(phi11*K11[ijk] + phi21*K12[ijk]) + 
         (phi11*xp1 + phi31*xp3)*K13[ijk] + 
         xp1*(phi21*K23[ijk] + phi31*K33[ijk])) + 
      ginv23*(phi11*(xp3*K12[ijk] + xp2*K13[ijk]) + 
         phi21*(xp3*K22[ijk] + xp2*K23[ijk]) + 
         phi31*(xp3*K23[ijk] + xp2*K33[ijk])) - 
      (phi11*xp1 + phi21*xp2 + phi31*xp3)*TrK[ijk])*sqrt(det))/PI
;

IS2[ijk]
=
(0.125*R*(ginv11*xp1*(phi12*K11[ijk] + phi22*K12[ijk] + phi32*K13[ijk]) + 
      ginv22*xp2*(phi12*K12[ijk] + phi22*K22[ijk] + phi32*K23[ijk]) + 
      ginv12*(phi12*(xp2*K11[ijk] + xp1*K12[ijk]) + 
         xp2*(phi22*K12[ijk] + phi32*K13[ijk]) + 
         xp1*(phi22*K22[ijk] + phi32*K23[ijk])) + 
      ginv33*xp3*(phi12*K13[ijk] + phi22*K23[ijk] + phi32*K33[ijk]) + 
      ginv13*(xp3*(phi12*K11[ijk] + phi22*K12[ijk]) + 
         (phi12*xp1 + phi32*xp3)*K13[ijk] + 
         xp1*(phi22*K23[ijk] + phi32*K33[ijk])) + 
      ginv23*(phi12*(xp3*K12[ijk] + xp2*K13[ijk]) + 
         phi22*(xp3*K22[ijk] + xp2*K23[ijk]) + 
         phi32*(xp3*K23[ijk] + xp2*K33[ijk])) - 
      (phi12*xp1 + phi22*xp2 + phi32*xp3)*TrK[ijk])*sqrt(det))/PI
;

IS3[ijk]
=
(0.125*R*(ginv11*xp1*(phi13*K11[ijk] + phi23*K12[ijk] + phi33*K13[ijk]) + 
      ginv22*xp2*(phi13*K12[ijk] + phi23*K22[ijk] + phi33*K23[ijk]) + 
      ginv12*(phi13*(xp2*K11[ijk] + xp1*K12[ijk]) + 
         xp2*(phi23*K12[ijk] + phi33*K13[ijk]) + 
         xp1*(phi23*K22[ijk] + phi33*K23[ijk])) + 
      ginv33*xp3*(phi13*K13[ijk] + phi23*K23[ijk] + phi33*K33[ijk]) + 
      ginv13*(xp3*(phi13*K11[ijk] + phi23*K12[ijk]) + 
         (phi13*xp1 + phi33*xp3)*K13[ijk] + 
         xp1*(phi23*K23[ijk] + phi33*K33[ijk])) + 
      ginv23*(phi13*(xp3*K12[ijk] + xp2*K13[ijk]) + 
         phi23*(xp3*K22[ijk] + xp2*K23[ijk]) + 
         phi33*(xp3*K23[ijk] + xp2*K33[ijk])) - 
      (phi13*xp1 + phi23*xp2 + phi33*xp3)*TrK[ijk])*sqrt(det))/PI
;



} endforinner; /* loop i, j, k */

}  /* function */

/* puncture_properties_integrand.c */
/* nvars = 23, nauxs = 46, n* = 350,  n/ = 31,  n+ = 323, n = 704, O = 0 */
