/* z4_boundary_box_rad.c */
/* Copyright (C) 1998 Bernd Bruegmann, 7.6.2022 */
/* Produced with Mathematica */

#include "bam.h"
#include "z4.h"

#define Power(x,y) (pow((double) (x), (double) (y)))
#define Sqrt(x)    sqrt(x)
#define Log(x)     log((double) (x))
#define pow2(x)    ((x)*(x))
#define pow3(x)    ((x)*(x)*(x))
#define pow4(x)    ((x)*(x)*(x)*(x))
#define pow2inv(x) (1.0/((x)*(x)))

#define Cal(x,y,z) ((x)?(y):(z))

#define Tan(x)     tan(x)
#define ArcTan(x)  atan(x)
#define Sin(x)     sin(x)
#define Cos(x)     cos(x)
#define Csc(x)     (1./sin(x))
#define Abs(x)     (fabs(x))
#define sqrt2      (sqrt(2))
#define Tanh(x)    tanh(x)
#define Sech(x)    (1/cosh(x))



void z4_boundary_box_rad(tVarList *unew, tVarList *upre, double c, tVarList *ucur)
{

tL *level = ucur->level;
int addlinear = (c != 0.0l);
double time = level->time;

double *g11 = vldataptr(ucur, METRIC_z4_INDX_VAR);
double *g12 = vldataptr(ucur, 1 + METRIC_z4_INDX_VAR);
double *g13 = vldataptr(ucur, 2 + METRIC_z4_INDX_VAR);
double *g22 = vldataptr(ucur, 3 + METRIC_z4_INDX_VAR);
double *g23 = vldataptr(ucur, 4 + METRIC_z4_INDX_VAR);
double *g33 = vldataptr(ucur, 5 + METRIC_z4_INDX_VAR);
double *chi = vldataptr(ucur, 6 + METRIC_z4_INDX_VAR);
double *A11 = vldataptr(ucur, 7 + METRIC_z4_INDX_VAR);
double *A12 = vldataptr(ucur, 8 + METRIC_z4_INDX_VAR);
double *A13 = vldataptr(ucur, 9 + METRIC_z4_INDX_VAR);
double *A22 = vldataptr(ucur, 10 + METRIC_z4_INDX_VAR);
double *A23 = vldataptr(ucur, 11 + METRIC_z4_INDX_VAR);
double *A33 = vldataptr(ucur, 12 + METRIC_z4_INDX_VAR);
double *Khat = vldataptr(ucur, 13 + METRIC_z4_INDX_VAR);
double *G1 = vldataptr(ucur, 14 + METRIC_z4_INDX_VAR);
double *G2 = vldataptr(ucur, 15 + METRIC_z4_INDX_VAR);
double *G3 = vldataptr(ucur, 16 + METRIC_z4_INDX_VAR);
double *Theta = vldataptr(ucur, 17 + METRIC_z4_INDX_VAR);
double *alpha = vldataptr(ucur, 18 + METRIC_z4_INDX_VAR);
double *beta1 = vldataptr(ucur, 19 + METRIC_z4_INDX_VAR);
double *beta2 = vldataptr(ucur, 20 + METRIC_z4_INDX_VAR);
double *beta3 = vldataptr(ucur, 21 + METRIC_z4_INDX_VAR);
double *B1 = vldataptr(ucur, 22 + METRIC_z4_INDX_VAR);
double *B2 = vldataptr(ucur, 23 + METRIC_z4_INDX_VAR);
double *B3 = vldataptr(ucur, 24 + METRIC_z4_INDX_VAR);
double *ng11 = vldataptr(unew, METRIC_z4_INDX_VAR);
double *ng12 = vldataptr(unew, 1 + METRIC_z4_INDX_VAR);
double *ng13 = vldataptr(unew, 2 + METRIC_z4_INDX_VAR);
double *ng22 = vldataptr(unew, 3 + METRIC_z4_INDX_VAR);
double *ng23 = vldataptr(unew, 4 + METRIC_z4_INDX_VAR);
double *ng33 = vldataptr(unew, 5 + METRIC_z4_INDX_VAR);
double *nchi = vldataptr(unew, 6 + METRIC_z4_INDX_VAR);
double *nA11 = vldataptr(unew, 7 + METRIC_z4_INDX_VAR);
double *nA12 = vldataptr(unew, 8 + METRIC_z4_INDX_VAR);
double *nA13 = vldataptr(unew, 9 + METRIC_z4_INDX_VAR);
double *nA22 = vldataptr(unew, 10 + METRIC_z4_INDX_VAR);
double *nA23 = vldataptr(unew, 11 + METRIC_z4_INDX_VAR);
double *nA33 = vldataptr(unew, 12 + METRIC_z4_INDX_VAR);
double *nKhat = vldataptr(unew, 13 + METRIC_z4_INDX_VAR);
double *nG1 = vldataptr(unew, 14 + METRIC_z4_INDX_VAR);
double *nG2 = vldataptr(unew, 15 + METRIC_z4_INDX_VAR);
double *nG3 = vldataptr(unew, 16 + METRIC_z4_INDX_VAR);
double *nTheta = vldataptr(unew, 17 + METRIC_z4_INDX_VAR);
double *nalpha = vldataptr(unew, 18 + METRIC_z4_INDX_VAR);
double *nbeta1 = vldataptr(unew, 19 + METRIC_z4_INDX_VAR);
double *nbeta2 = vldataptr(unew, 20 + METRIC_z4_INDX_VAR);
double *nbeta3 = vldataptr(unew, 21 + METRIC_z4_INDX_VAR);
double *nB1 = vldataptr(unew, 22 + METRIC_z4_INDX_VAR);
double *nB2 = vldataptr(unew, 23 + METRIC_z4_INDX_VAR);
double *nB3 = vldataptr(unew, 24 + METRIC_z4_INDX_VAR);
double *pg11 = vldataptr(upre, METRIC_z4_INDX_VAR);
double *pg12 = vldataptr(upre, 1 + METRIC_z4_INDX_VAR);
double *pg13 = vldataptr(upre, 2 + METRIC_z4_INDX_VAR);
double *pg22 = vldataptr(upre, 3 + METRIC_z4_INDX_VAR);
double *pg23 = vldataptr(upre, 4 + METRIC_z4_INDX_VAR);
double *pg33 = vldataptr(upre, 5 + METRIC_z4_INDX_VAR);
double *pchi = vldataptr(upre, 6 + METRIC_z4_INDX_VAR);
double *pA11 = vldataptr(upre, 7 + METRIC_z4_INDX_VAR);
double *pA12 = vldataptr(upre, 8 + METRIC_z4_INDX_VAR);
double *pA13 = vldataptr(upre, 9 + METRIC_z4_INDX_VAR);
double *pA22 = vldataptr(upre, 10 + METRIC_z4_INDX_VAR);
double *pA23 = vldataptr(upre, 11 + METRIC_z4_INDX_VAR);
double *pA33 = vldataptr(upre, 12 + METRIC_z4_INDX_VAR);
double *pKhat = vldataptr(upre, 13 + METRIC_z4_INDX_VAR);
double *pG1 = vldataptr(upre, 14 + METRIC_z4_INDX_VAR);
double *pG2 = vldataptr(upre, 15 + METRIC_z4_INDX_VAR);
double *pG3 = vldataptr(upre, 16 + METRIC_z4_INDX_VAR);
double *pTheta = vldataptr(upre, 17 + METRIC_z4_INDX_VAR);
double *palpha = vldataptr(upre, 18 + METRIC_z4_INDX_VAR);
double *pbeta1 = vldataptr(upre, 19 + METRIC_z4_INDX_VAR);
double *pbeta2 = vldataptr(upre, 20 + METRIC_z4_INDX_VAR);
double *pbeta3 = vldataptr(upre, 21 + METRIC_z4_INDX_VAR);
double *pB1 = vldataptr(upre, 22 + METRIC_z4_INDX_VAR);
double *pB2 = vldataptr(upre, 23 + METRIC_z4_INDX_VAR);
double *pB3 = vldataptr(upre, 24 + METRIC_z4_INDX_VAR);

const int N_points_away = 1;
const int order_dissipation   = Geti("order_dissipation");
const double dissfactor       = get_dissipation_factor(level);
const double kappa1      = Getd("z4_kappa1");
const double kappa2      = Getd("z4_kappa2");
const double chiDivFloor = Getd("z4_chi_div_floor");
const double chipsipower = Getd("z4_chi_psipower");
const double shiftdriver = Getd("z4_shiftdriver") * Getv("z4_bc_use_eta","yes");
const double *xp = level->v[Ind("x")];
const double *yp = level->v[Ind("y")];
const double *zp = level->v[Ind("z")];


double dx = level->dx;
double dy = level->dy;
double dz = level->dz;
double oo2dx = 1/(2*dx), oo2dy = 1/(2*dy), oo2dz = 1/(2*dz);
double oodx2 = 1/(dx*dx), oody2 = 1/(dy*dy), oodz2 = 1/(dz*dz);
double oo4dxdy = 1/(4*dx*dy);
double oo4dxdz = 1/(4*dx*dz);
double oo4dydz = 1/(4*dy*dz);

double dA111 = 0.;
double dA112 = 0.;
double dA113 = 0.;
double dA122 = 0.;
double dA123 = 0.;
double dA133 = 0.;
double dA211 = 0.;
double dA212 = 0.;
double dA213 = 0.;
double dA222 = 0.;
double dA223 = 0.;
double dA233 = 0.;
double dA311 = 0.;
double dA312 = 0.;
double dA313 = 0.;
double dA322 = 0.;
double dA323 = 0.;
double dA333 = 0.;
double dG11 = 0.;
double dG12 = 0.;
double dG13 = 0.;
double dG21 = 0.;
double dG22 = 0.;
double dG23 = 0.;
double dG31 = 0.;
double dG32 = 0.;
double dG33 = 0.;
double dKhat1 = 0.;
double dKhat2 = 0.;
double dKhat3 = 0.;
double dTheta1 = 0.;
double dTheta2 = 0.;
double dTheta3 = 0.;
double r = 0.;
double rA11 = 0.;
double rA12 = 0.;
double rA13 = 0.;
double rA22 = 0.;
double rA23 = 0.;
double rA33 = 0.;
double rG1 = 0.;
double rG2 = 0.;
double rG3 = 0.;
double rKhat = 0.;
double rTheta = 0.;
double sup1 = 0.;
double sup2 = 0.;
double sup3 = 0.;
double x = 0.;
double y = 0.;
double z = 0.;



forinnerpoints_ijk_openmp(level) {



if (!(boundaryNaway(N_points_away))) continue; 


if (CheckForNANandINF(25,                                                   
    g11[ijk],g12[ijk],g13[ijk],g22[ijk],g23[ijk],g33[ijk],
    A11[ijk],A12[ijk],A13[ijk],A22[ijk],A23[ijk],A33[ijk], 
    G1[ijk],G2[ijk],G3[ijk], Khat[ijk],chi[ijk], Theta[ijk],
    alpha[ijk],beta1[ijk],beta2[ijk],beta3[ijk],B1[ijk],B2[ijk],B3[ijk])) {
    printf("problem with vars in Z4d_boundary.m\n");
    printf("x=%2.5e, y=%2.5e, z=%2.5e\n",xp[ijk],yp[ijk],zp[ijk]);
    }
dKhat1
=
oo2dx*(-Khat[-di + ijk] + Khat[di + ijk])
;

dKhat2
=
oo2dy*(-Khat[-dj + ijk] + Khat[dj + ijk])
;

dKhat3
=
oo2dz*(-Khat[-dk + ijk] + Khat[dk + ijk])
;

dA111
=
oo2dx*(-A11[-di + ijk] + A11[di + ijk])
;

dA112
=
oo2dx*(-A12[-di + ijk] + A12[di + ijk])
;

dA113
=
oo2dx*(-A13[-di + ijk] + A13[di + ijk])
;

dA122
=
oo2dx*(-A22[-di + ijk] + A22[di + ijk])
;

dA123
=
oo2dx*(-A23[-di + ijk] + A23[di + ijk])
;

dA133
=
oo2dx*(-A33[-di + ijk] + A33[di + ijk])
;

dA211
=
oo2dy*(-A11[-dj + ijk] + A11[dj + ijk])
;

dA212
=
oo2dy*(-A12[-dj + ijk] + A12[dj + ijk])
;

dA213
=
oo2dy*(-A13[-dj + ijk] + A13[dj + ijk])
;

dA222
=
oo2dy*(-A22[-dj + ijk] + A22[dj + ijk])
;

dA223
=
oo2dy*(-A23[-dj + ijk] + A23[dj + ijk])
;

dA233
=
oo2dy*(-A33[-dj + ijk] + A33[dj + ijk])
;

dA311
=
oo2dz*(-A11[-dk + ijk] + A11[dk + ijk])
;

dA312
=
oo2dz*(-A12[-dk + ijk] + A12[dk + ijk])
;

dA313
=
oo2dz*(-A13[-dk + ijk] + A13[dk + ijk])
;

dA322
=
oo2dz*(-A22[-dk + ijk] + A22[dk + ijk])
;

dA323
=
oo2dz*(-A23[-dk + ijk] + A23[dk + ijk])
;

dA333
=
oo2dz*(-A33[-dk + ijk] + A33[dk + ijk])
;

dG11
=
oo2dx*(-G1[-di + ijk] + G1[di + ijk])
;

dG12
=
oo2dx*(-G2[-di + ijk] + G2[di + ijk])
;

dG13
=
oo2dx*(-G3[-di + ijk] + G3[di + ijk])
;

dG21
=
oo2dy*(-G1[-dj + ijk] + G1[dj + ijk])
;

dG22
=
oo2dy*(-G2[-dj + ijk] + G2[dj + ijk])
;

dG23
=
oo2dy*(-G3[-dj + ijk] + G3[dj + ijk])
;

dG31
=
oo2dz*(-G1[-dk + ijk] + G1[dk + ijk])
;

dG32
=
oo2dz*(-G2[-dk + ijk] + G2[dk + ijk])
;

dG33
=
oo2dz*(-G3[-dk + ijk] + G3[dk + ijk])
;

dTheta1
=
oo2dx*(-Theta[-di + ijk] + Theta[di + ijk])
;

dTheta2
=
oo2dy*(-Theta[-dj + ijk] + Theta[dj + ijk])
;

dTheta3
=
oo2dz*(-Theta[-dk + ijk] + Theta[dk + ijk])
;

r
=
0
;

x
=
0
;

y
=
0
;

z
=
0
;

sup1
=
0
;

sup2
=
0
;

sup3
=
0
;


x=xp[ijk]; y=yp[ijk]; z=zp[ijk]; 


r=sqrt(x*x + y*y + z*z); 


sup1 = x/r; 


sup2 = y/r; 


sup3 = z/r; 

rTheta
=
-(dTheta1*sup1) - dTheta2*sup2 - dTheta3*sup3 - Theta[ijk]/r
;

rKhat
=
-((dKhat1*sup1 + dKhat2*sup2 + dKhat3*sup3 + Khat[ijk]/r)*sqrt(2.))
;

rG1
=
-(dG11*sup1) - dG21*sup2 - dG31*sup3 - G1[ijk]/r
;

rG2
=
-(dG12*sup1) - dG22*sup2 - dG32*sup3 - G2[ijk]/r
;

rG3
=
-(dG13*sup1) - dG23*sup2 - dG33*sup3 - G3[ijk]/r
;

rA11
=
-(dA111*sup1) - dA211*sup2 - dA311*sup3 - A11[ijk]/r
;

rA12
=
-(dA112*sup1) - dA212*sup2 - dA312*sup3 - A12[ijk]/r
;

rA13
=
-(dA113*sup1) - dA213*sup2 - dA313*sup3 - A13[ijk]/r
;

rA22
=
-(dA122*sup1) - dA222*sup2 - dA322*sup3 - A22[ijk]/r
;

rA23
=
-(dA123*sup1) - dA223*sup2 - dA323*sup3 - A23[ijk]/r
;

rA33
=
-(dA133*sup1) - dA233*sup2 - dA333*sup3 - A33[ijk]/r
;


if (CheckForNANandINF(11,                                          
    rA11,rA12,rA13,rA22,rA23,rA33, 
    rG1,rG2,rG3,rKhat,rTheta)) {
    printf("nans in RHS in z4_boundary_box_rad.m\n");
    printf("x=%2.5e, y=%2.5e, z=%2.5e\n",xp[ijk],yp[ijk],zp[ijk]);
  }

/* conditional */
if (addlinear) {

nKhat[ijk]
=
c*rKhat + pKhat[ijk]
;

nA11[ijk]
=
c*rA11 + pA11[ijk]
;

nA12[ijk]
=
c*rA12 + pA12[ijk]
;

nA13[ijk]
=
c*rA13 + pA13[ijk]
;

nA22[ijk]
=
c*rA22 + pA22[ijk]
;

nA23[ijk]
=
c*rA23 + pA23[ijk]
;

nA33[ijk]
=
c*rA33 + pA33[ijk]
;

nG1[ijk]
=
c*rG1 + pG1[ijk]
;

nG2[ijk]
=
c*rG2 + pG2[ijk]
;

nG3[ijk]
=
c*rG3 + pG3[ijk]
;

nTheta[ijk]
=
c*rTheta + pTheta[ijk]
;


} else { /* if (!addlinear) */

nKhat[ijk]
=
rKhat
;

nA11[ijk]
=
rA11
;

nA12[ijk]
=
rA12
;

nA13[ijk]
=
rA13
;

nA22[ijk]
=
rA22
;

nA23[ijk]
=
rA23
;

nA33[ijk]
=
rA33
;

nG1[ijk]
=
rG1
;

nG2[ijk]
=
rG2
;

nG3[ijk]
=
rG3
;

nTheta[ijk]
=
rTheta
;

}
/* if (addlinear) */




} endfor_ijk_openmp; /* loop i, j, k */

}  /* function */

/* z4_boundary_box_rad.c */
/* nvars = 78, nauxs = 51, n* = 199,  n/ = 42,  n+ = 302, n = 543, O = 0 */
