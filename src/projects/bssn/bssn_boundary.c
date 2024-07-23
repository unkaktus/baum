/* bssn_boundary.c */
/* Copyright (C) 1998 Bernd Bruegmann, 28.6.2012 */
/* Produced with Mathematica */

#include "bam.h"
#include "bssn.h"

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



void bssn_boundary(tVarList *unew, tVarList *upre, double c, tVarList *ucur)
{

tL *level = ucur->level;
int addlinear = (c != 0.0l);

double *g11 = vldataptr(ucur, METRIC_bssn_INDX_VAR);
double *g12 = vldataptr(ucur, 1 + METRIC_bssn_INDX_VAR);
double *g13 = vldataptr(ucur, 2 + METRIC_bssn_INDX_VAR);
double *g22 = vldataptr(ucur, 3 + METRIC_bssn_INDX_VAR);
double *g23 = vldataptr(ucur, 4 + METRIC_bssn_INDX_VAR);
double *g33 = vldataptr(ucur, 5 + METRIC_bssn_INDX_VAR);
double *chi = vldataptr(ucur, 6 + METRIC_bssn_INDX_VAR);
double *A11 = vldataptr(ucur, 7 + METRIC_bssn_INDX_VAR);
double *A12 = vldataptr(ucur, 8 + METRIC_bssn_INDX_VAR);
double *A13 = vldataptr(ucur, 9 + METRIC_bssn_INDX_VAR);
double *A22 = vldataptr(ucur, 10 + METRIC_bssn_INDX_VAR);
double *A23 = vldataptr(ucur, 11 + METRIC_bssn_INDX_VAR);
double *A33 = vldataptr(ucur, 12 + METRIC_bssn_INDX_VAR);
double *K = vldataptr(ucur, 13 + METRIC_bssn_INDX_VAR);
double *G1 = vldataptr(ucur, 14 + METRIC_bssn_INDX_VAR);
double *G2 = vldataptr(ucur, 15 + METRIC_bssn_INDX_VAR);
double *G3 = vldataptr(ucur, 16 + METRIC_bssn_INDX_VAR);
double *alpha = vldataptr(ucur, 17 + METRIC_bssn_INDX_VAR);
double *beta1 = vldataptr(ucur, 18 + METRIC_bssn_INDX_VAR);
double *beta2 = vldataptr(ucur, 19 + METRIC_bssn_INDX_VAR);
double *beta3 = vldataptr(ucur, 20 + METRIC_bssn_INDX_VAR);
double *B1 = vldataptr(ucur, 21 + METRIC_bssn_INDX_VAR);
double *B2 = vldataptr(ucur, 22 + METRIC_bssn_INDX_VAR);
double *B3 = vldataptr(ucur, 23 + METRIC_bssn_INDX_VAR);
double *ng11 = vldataptr(unew, METRIC_bssn_INDX_VAR);
double *ng12 = vldataptr(unew, 1 + METRIC_bssn_INDX_VAR);
double *ng13 = vldataptr(unew, 2 + METRIC_bssn_INDX_VAR);
double *ng22 = vldataptr(unew, 3 + METRIC_bssn_INDX_VAR);
double *ng23 = vldataptr(unew, 4 + METRIC_bssn_INDX_VAR);
double *ng33 = vldataptr(unew, 5 + METRIC_bssn_INDX_VAR);
double *nchi = vldataptr(unew, 6 + METRIC_bssn_INDX_VAR);
double *nA11 = vldataptr(unew, 7 + METRIC_bssn_INDX_VAR);
double *nA12 = vldataptr(unew, 8 + METRIC_bssn_INDX_VAR);
double *nA13 = vldataptr(unew, 9 + METRIC_bssn_INDX_VAR);
double *nA22 = vldataptr(unew, 10 + METRIC_bssn_INDX_VAR);
double *nA23 = vldataptr(unew, 11 + METRIC_bssn_INDX_VAR);
double *nA33 = vldataptr(unew, 12 + METRIC_bssn_INDX_VAR);
double *nK = vldataptr(unew, 13 + METRIC_bssn_INDX_VAR);
double *nG1 = vldataptr(unew, 14 + METRIC_bssn_INDX_VAR);
double *nG2 = vldataptr(unew, 15 + METRIC_bssn_INDX_VAR);
double *nG3 = vldataptr(unew, 16 + METRIC_bssn_INDX_VAR);
double *nalpha = vldataptr(unew, 17 + METRIC_bssn_INDX_VAR);
double *nbeta1 = vldataptr(unew, 18 + METRIC_bssn_INDX_VAR);
double *nbeta2 = vldataptr(unew, 19 + METRIC_bssn_INDX_VAR);
double *nbeta3 = vldataptr(unew, 20 + METRIC_bssn_INDX_VAR);
double *nB1 = vldataptr(unew, 21 + METRIC_bssn_INDX_VAR);
double *nB2 = vldataptr(unew, 22 + METRIC_bssn_INDX_VAR);
double *nB3 = vldataptr(unew, 23 + METRIC_bssn_INDX_VAR);
double *pg11 = vldataptr(upre, METRIC_bssn_INDX_VAR);
double *pg12 = vldataptr(upre, 1 + METRIC_bssn_INDX_VAR);
double *pg13 = vldataptr(upre, 2 + METRIC_bssn_INDX_VAR);
double *pg22 = vldataptr(upre, 3 + METRIC_bssn_INDX_VAR);
double *pg23 = vldataptr(upre, 4 + METRIC_bssn_INDX_VAR);
double *pg33 = vldataptr(upre, 5 + METRIC_bssn_INDX_VAR);
double *pchi = vldataptr(upre, 6 + METRIC_bssn_INDX_VAR);
double *pA11 = vldataptr(upre, 7 + METRIC_bssn_INDX_VAR);
double *pA12 = vldataptr(upre, 8 + METRIC_bssn_INDX_VAR);
double *pA13 = vldataptr(upre, 9 + METRIC_bssn_INDX_VAR);
double *pA22 = vldataptr(upre, 10 + METRIC_bssn_INDX_VAR);
double *pA23 = vldataptr(upre, 11 + METRIC_bssn_INDX_VAR);
double *pA33 = vldataptr(upre, 12 + METRIC_bssn_INDX_VAR);
double *pK = vldataptr(upre, 13 + METRIC_bssn_INDX_VAR);
double *pG1 = vldataptr(upre, 14 + METRIC_bssn_INDX_VAR);
double *pG2 = vldataptr(upre, 15 + METRIC_bssn_INDX_VAR);
double *pG3 = vldataptr(upre, 16 + METRIC_bssn_INDX_VAR);
double *palpha = vldataptr(upre, 17 + METRIC_bssn_INDX_VAR);
double *pbeta1 = vldataptr(upre, 18 + METRIC_bssn_INDX_VAR);
double *pbeta2 = vldataptr(upre, 19 + METRIC_bssn_INDX_VAR);
double *pbeta3 = vldataptr(upre, 20 + METRIC_bssn_INDX_VAR);
double *pB1 = vldataptr(upre, 21 + METRIC_bssn_INDX_VAR);
double *pB2 = vldataptr(upre, 22 + METRIC_bssn_INDX_VAR);
double *pB3 = vldataptr(upre, 23 + METRIC_bssn_INDX_VAR);

const int order               = Geti("order_centered");
const int N                   = Geti("boundary_N_extrapolate");
const int order_dissipation   = Geti("order_dissipation");
const double dissfactor       = get_dissipation_factor(level);
double shiftgammacoeff = Getd("bssn_shiftgammacoeff");
const double *xp = level->v[Ind("x")];
const double *yp = level->v[Ind("y")];
const double *zp = level->v[Ind("z")];
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

double ACABTF11 = 0.;
double ACABTF12 = 0.;
double ACABTF13 = 0.;
double ACABTF21 = 0.;
double ACABTF22 = 0.;
double ACABTF23 = 0.;
double ACABTF31 = 0.;
double ACABTF32 = 0.;
double ACABTF33 = 0.;
double ACsA1 = 0.;
double ACsA2 = 0.;
double ACsA3 = 0.;
double ACss = 0.;
double ADMginv11 = 0.;
double ADMginv12 = 0.;
double ADMginv13 = 0.;
double ADMginv22 = 0.;
double ADMginv23 = 0.;
double ADMginv33 = 0.;
double Ainv11 = 0.;
double Ainv12 = 0.;
double Ainv13 = 0.;
double Ainv22 = 0.;
double Ainv23 = 0.;
double Ainv33 = 0.;
double BA1 = 0.;
double BA2 = 0.;
double BA3 = 0.;
double Bs = 0.;
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
double DACABTF11 = 0.;
double DACABTF12 = 0.;
double DACABTF13 = 0.;
double DACABTF21 = 0.;
double DACABTF22 = 0.;
double DACABTF23 = 0.;
double DACABTF31 = 0.;
double DACABTF32 = 0.;
double DACABTF33 = 0.;
double DACsA1 = 0.;
double DACsA2 = 0.;
double DACsA3 = 0.;
double DACss = 0.;
double dASST111 = 0.;
double dASST112 = 0.;
double dASST113 = 0.;
double dASST122 = 0.;
double dASST123 = 0.;
double dASST133 = 0.;
double dASST211 = 0.;
double dASST212 = 0.;
double dASST213 = 0.;
double dASST222 = 0.;
double dASST223 = 0.;
double dASST233 = 0.;
double dASST311 = 0.;
double dASST312 = 0.;
double dASST313 = 0.;
double dASST322 = 0.;
double dASST323 = 0.;
double dASST333 = 0.;
double db11 = 0.;
double dB11 = 0.;
double db12 = 0.;
double dB12 = 0.;
double db13 = 0.;
double dB13 = 0.;
double db21 = 0.;
double dB21 = 0.;
double db22 = 0.;
double dB22 = 0.;
double db23 = 0.;
double dB23 = 0.;
double db31 = 0.;
double dB31 = 0.;
double db32 = 0.;
double dB32 = 0.;
double db33 = 0.;
double dB33 = 0.;
double DBA1 = 0.;
double DBA2 = 0.;
double DBA3 = 0.;
double DBs = 0.;
double dbSST11 = 0.;
double dBSST11 = 0.;
double dbSST12 = 0.;
double dBSST12 = 0.;
double dbSST13 = 0.;
double dBSST13 = 0.;
double dbSST21 = 0.;
double dBSST21 = 0.;
double dbSST22 = 0.;
double dBSST22 = 0.;
double dbSST23 = 0.;
double dBSST23 = 0.;
double dbSST31 = 0.;
double dBSST31 = 0.;
double dbSST32 = 0.;
double dBSST32 = 0.;
double dbSST33 = 0.;
double dBSST33 = 0.;
double ddRdr = 0.;
double detginv = 0.;
double dG11 = 0.;
double dg111 = 0.;
double dg112 = 0.;
double dg113 = 0.;
double dG12 = 0.;
double dg122 = 0.;
double dg123 = 0.;
double dG13 = 0.;
double dg133 = 0.;
double dG21 = 0.;
double dg211 = 0.;
double dg212 = 0.;
double dg213 = 0.;
double dG22 = 0.;
double dg222 = 0.;
double dg223 = 0.;
double dG23 = 0.;
double dg233 = 0.;
double dG31 = 0.;
double dg311 = 0.;
double dg312 = 0.;
double dg313 = 0.;
double dG32 = 0.;
double dg322 = 0.;
double dg323 = 0.;
double dG33 = 0.;
double dg333 = 0.;
double DGamA1 = 0.;
double DGamA2 = 0.;
double DGamA3 = 0.;
double DGams = 0.;
double dGSST11 = 0.;
double dgSST111 = 0.;
double dgSST112 = 0.;
double dgSST113 = 0.;
double dGSST12 = 0.;
double dgSST122 = 0.;
double dgSST123 = 0.;
double dGSST13 = 0.;
double dgSST133 = 0.;
double dGSST21 = 0.;
double dgSST211 = 0.;
double dgSST212 = 0.;
double dgSST213 = 0.;
double dGSST22 = 0.;
double dgSST222 = 0.;
double dgSST223 = 0.;
double dGSST23 = 0.;
double dgSST233 = 0.;
double dGSST31 = 0.;
double dgSST311 = 0.;
double dgSST312 = 0.;
double dgSST313 = 0.;
double dGSST32 = 0.;
double dgSST322 = 0.;
double dgSST323 = 0.;
double dGSST33 = 0.;
double dgSST333 = 0.;
double divbeta = 0.;
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
double DK = 0.;
double dK1 = 0.;
double dK2 = 0.;
double dK3 = 0.;
double dKSST1 = 0.;
double dKSST2 = 0.;
double dKSST3 = 0.;
double dRdr = 0.;
double GamA1 = 0.;
double GamA2 = 0.;
double GamA3 = 0.;
double Gams = 0.;
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
double lieg11 = 0.;
double lieg12 = 0.;
double lieg13 = 0.;
double lieg22 = 0.;
double lieg23 = 0.;
double lieg33 = 0.;
double modshatARG = 0.;
double muL = 0.;
double muS = 0.;
double oomodshat = 0.;
double qdd11 = 0.;
double qdd12 = 0.;
double qdd13 = 0.;
double qdd22 = 0.;
double qdd23 = 0.;
double qdd33 = 0.;
double qPhysuudd1111 = 0.;
double qPhysuudd1112 = 0.;
double qPhysuudd1113 = 0.;
double qPhysuudd1122 = 0.;
double qPhysuudd1123 = 0.;
double qPhysuudd1133 = 0.;
double qPhysuudd1211 = 0.;
double qPhysuudd1212 = 0.;
double qPhysuudd1213 = 0.;
double qPhysuudd1222 = 0.;
double qPhysuudd1223 = 0.;
double qPhysuudd1233 = 0.;
double qPhysuudd1311 = 0.;
double qPhysuudd1312 = 0.;
double qPhysuudd1313 = 0.;
double qPhysuudd1322 = 0.;
double qPhysuudd1323 = 0.;
double qPhysuudd1333 = 0.;
double qPhysuudd2211 = 0.;
double qPhysuudd2212 = 0.;
double qPhysuudd2213 = 0.;
double qPhysuudd2222 = 0.;
double qPhysuudd2223 = 0.;
double qPhysuudd2233 = 0.;
double qPhysuudd2311 = 0.;
double qPhysuudd2312 = 0.;
double qPhysuudd2313 = 0.;
double qPhysuudd2322 = 0.;
double qPhysuudd2323 = 0.;
double qPhysuudd2333 = 0.;
double qPhysuudd3311 = 0.;
double qPhysuudd3312 = 0.;
double qPhysuudd3313 = 0.;
double qPhysuudd3322 = 0.;
double qPhysuudd3323 = 0.;
double qPhysuudd3333 = 0.;
double qud11 = 0.;
double qud12 = 0.;
double qud13 = 0.;
double qud21 = 0.;
double qud22 = 0.;
double qud23 = 0.;
double qud31 = 0.;
double qud32 = 0.;
double qud33 = 0.;
double quu11 = 0.;
double quu12 = 0.;
double quu13 = 0.;
double quu22 = 0.;
double quu23 = 0.;
double quu33 = 0.;
double r = 0.;
double rA11 = 0.;
double rA12 = 0.;
double rA13 = 0.;
double rA22 = 0.;
double rA23 = 0.;
double rA33 = 0.;
double rACABTF11 = 0.;
double rACABTF12 = 0.;
double rACABTF13 = 0.;
double rACABTF22 = 0.;
double rACABTF23 = 0.;
double rACABTF33 = 0.;
double rACqq = 0.;
double rACsA1 = 0.;
double rACsA2 = 0.;
double rACsA3 = 0.;
double rACss = 0.;
double rB1 = 0.;
double rB2 = 0.;
double rB3 = 0.;
double rBA1 = 0.;
double rBA2 = 0.;
double rBA3 = 0.;
double rBs = 0.;
double rG1 = 0.;
double rG2 = 0.;
double rG3 = 0.;
double rGamA1 = 0.;
double rGamA2 = 0.;
double rGamA3 = 0.;
double rGams = 0.;
double rK = 0.;
double sdown1 = 0.;
double sdown2 = 0.;
double sdown3 = 0.;
double shat1 = 0.;
double shat2 = 0.;
double shat3 = 0.;
double sup1 = 0.;
double sup2 = 0.;
double sup3 = 0.;
double totdivbeta = 0.;
double vbetaA = 0.;
double vbetas = 0.;



forallpoints_ijk_openmp(level) {
if (!dequal(rp[ijk],level->bbox[1]-N*level->dx)) continue;
if (j==0 || j==box->n-1 || k==0 || k==box->o-1) continue;




if (CheckForNANandINF(24,                                                   
    g11[ijk],g12[ijk],g13[ijk],g22[ijk],g23[ijk],g33[ijk],
    A11[ijk],A12[ijk],A13[ijk],A22[ijk],A23[ijk],A33[ijk], 
    G1[ijk],G2[ijk],G3[ijk], K[ijk],chi[ijk],
    alpha[ijk],beta1[ijk],beta2[ijk],beta3[ijk],B1[ijk],B2[ijk],B3[ijk])) {
    printf("problem with vars in bssn_boundary.m\n");
    printf("x=%2.5e, y=%2.5e, z=%2.5e\n",xp[ijk],yp[ijk],zp[ijk]);}
if (order == 2 || boundaryNaway(1)) { 

db11
=
oo2dx*(-beta1[-di + ijk] + beta1[di + ijk])
;

db12
=
oo2dx*(-beta2[-di + ijk] + beta2[di + ijk])
;

db13
=
oo2dx*(-beta3[-di + ijk] + beta3[di + ijk])
;

db21
=
oo2dy*(-beta1[-dj + ijk] + beta1[dj + ijk])
;

db22
=
oo2dy*(-beta2[-dj + ijk] + beta2[dj + ijk])
;

db23
=
oo2dy*(-beta3[-dj + ijk] + beta3[dj + ijk])
;

db31
=
oo2dz*(-beta1[-dk + ijk] + beta1[dk + ijk])
;

db32
=
oo2dz*(-beta2[-dk + ijk] + beta2[dk + ijk])
;

db33
=
oo2dz*(-beta3[-dk + ijk] + beta3[dk + ijk])
;

dB11
=
oo2dx*(-B1[-di + ijk] + B1[di + ijk])
;

dB12
=
oo2dx*(-B2[-di + ijk] + B2[di + ijk])
;

dB13
=
oo2dx*(-B3[-di + ijk] + B3[di + ijk])
;

dB21
=
oo2dy*(-B1[-dj + ijk] + B1[dj + ijk])
;

dB22
=
oo2dy*(-B2[-dj + ijk] + B2[dj + ijk])
;

dB23
=
oo2dy*(-B3[-dj + ijk] + B3[dj + ijk])
;

dB31
=
oo2dz*(-B1[-dk + ijk] + B1[dk + ijk])
;

dB32
=
oo2dz*(-B2[-dk + ijk] + B2[dk + ijk])
;

dB33
=
oo2dz*(-B3[-dk + ijk] + B3[dk + ijk])
;

dg111
=
oo2dx*(-g11[-di + ijk] + g11[di + ijk])
;

dg112
=
oo2dx*(-g12[-di + ijk] + g12[di + ijk])
;

dg113
=
oo2dx*(-g13[-di + ijk] + g13[di + ijk])
;

dg122
=
oo2dx*(-g22[-di + ijk] + g22[di + ijk])
;

dg123
=
oo2dx*(-g23[-di + ijk] + g23[di + ijk])
;

dg133
=
oo2dx*(-g33[-di + ijk] + g33[di + ijk])
;

dg211
=
oo2dy*(-g11[-dj + ijk] + g11[dj + ijk])
;

dg212
=
oo2dy*(-g12[-dj + ijk] + g12[dj + ijk])
;

dg213
=
oo2dy*(-g13[-dj + ijk] + g13[dj + ijk])
;

dg222
=
oo2dy*(-g22[-dj + ijk] + g22[dj + ijk])
;

dg223
=
oo2dy*(-g23[-dj + ijk] + g23[dj + ijk])
;

dg233
=
oo2dy*(-g33[-dj + ijk] + g33[dj + ijk])
;

dg311
=
oo2dz*(-g11[-dk + ijk] + g11[dk + ijk])
;

dg312
=
oo2dz*(-g12[-dk + ijk] + g12[dk + ijk])
;

dg313
=
oo2dz*(-g13[-dk + ijk] + g13[dk + ijk])
;

dg322
=
oo2dz*(-g22[-dk + ijk] + g22[dk + ijk])
;

dg323
=
oo2dz*(-g23[-dk + ijk] + g23[dk + ijk])
;

dg333
=
oo2dz*(-g33[-dk + ijk] + g33[dk + ijk])
;

dK1
=
oo2dx*(-K[-di + ijk] + K[di + ijk])
;

dK2
=
oo2dy*(-K[-dj + ijk] + K[dj + ijk])
;

dK3
=
oo2dz*(-K[-dk + ijk] + K[dk + ijk])
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


} else if (order == 4 || boundaryNaway(2)) { 

db11
=
0.16666666666666666667*oo2dx*(beta1[-2*di + ijk] + 
    8.*(-beta1[-di + ijk] + beta1[di + ijk]) - beta1[2*di + ijk])
;

db12
=
0.16666666666666666667*oo2dx*(beta2[-2*di + ijk] + 
    8.*(-beta2[-di + ijk] + beta2[di + ijk]) - beta2[2*di + ijk])
;

db13
=
0.16666666666666666667*oo2dx*(beta3[-2*di + ijk] + 
    8.*(-beta3[-di + ijk] + beta3[di + ijk]) - beta3[2*di + ijk])
;

db21
=
0.16666666666666666667*oo2dy*(beta1[-2*dj + ijk] + 
    8.*(-beta1[-dj + ijk] + beta1[dj + ijk]) - beta1[2*dj + ijk])
;

db22
=
0.16666666666666666667*oo2dy*(beta2[-2*dj + ijk] + 
    8.*(-beta2[-dj + ijk] + beta2[dj + ijk]) - beta2[2*dj + ijk])
;

db23
=
0.16666666666666666667*oo2dy*(beta3[-2*dj + ijk] + 
    8.*(-beta3[-dj + ijk] + beta3[dj + ijk]) - beta3[2*dj + ijk])
;

db31
=
0.16666666666666666667*oo2dz*(beta1[-2*dk + ijk] + 
    8.*(-beta1[-dk + ijk] + beta1[dk + ijk]) - beta1[2*dk + ijk])
;

db32
=
0.16666666666666666667*oo2dz*(beta2[-2*dk + ijk] + 
    8.*(-beta2[-dk + ijk] + beta2[dk + ijk]) - beta2[2*dk + ijk])
;

db33
=
0.16666666666666666667*oo2dz*(beta3[-2*dk + ijk] + 
    8.*(-beta3[-dk + ijk] + beta3[dk + ijk]) - beta3[2*dk + ijk])
;

dB11
=
0.16666666666666666667*oo2dx*(B1[-2*di + ijk] + 
    8.*(-B1[-di + ijk] + B1[di + ijk]) - B1[2*di + ijk])
;

dB12
=
0.16666666666666666667*oo2dx*(B2[-2*di + ijk] + 
    8.*(-B2[-di + ijk] + B2[di + ijk]) - B2[2*di + ijk])
;

dB13
=
0.16666666666666666667*oo2dx*(B3[-2*di + ijk] + 
    8.*(-B3[-di + ijk] + B3[di + ijk]) - B3[2*di + ijk])
;

dB21
=
0.16666666666666666667*oo2dy*(B1[-2*dj + ijk] + 
    8.*(-B1[-dj + ijk] + B1[dj + ijk]) - B1[2*dj + ijk])
;

dB22
=
0.16666666666666666667*oo2dy*(B2[-2*dj + ijk] + 
    8.*(-B2[-dj + ijk] + B2[dj + ijk]) - B2[2*dj + ijk])
;

dB23
=
0.16666666666666666667*oo2dy*(B3[-2*dj + ijk] + 
    8.*(-B3[-dj + ijk] + B3[dj + ijk]) - B3[2*dj + ijk])
;

dB31
=
0.16666666666666666667*oo2dz*(B1[-2*dk + ijk] + 
    8.*(-B1[-dk + ijk] + B1[dk + ijk]) - B1[2*dk + ijk])
;

dB32
=
0.16666666666666666667*oo2dz*(B2[-2*dk + ijk] + 
    8.*(-B2[-dk + ijk] + B2[dk + ijk]) - B2[2*dk + ijk])
;

dB33
=
0.16666666666666666667*oo2dz*(B3[-2*dk + ijk] + 
    8.*(-B3[-dk + ijk] + B3[dk + ijk]) - B3[2*dk + ijk])
;

dg111
=
0.16666666666666666667*oo2dx*(g11[-2*di + ijk] + 
    8.*(-g11[-di + ijk] + g11[di + ijk]) - g11[2*di + ijk])
;

dg112
=
0.16666666666666666667*oo2dx*(g12[-2*di + ijk] + 
    8.*(-g12[-di + ijk] + g12[di + ijk]) - g12[2*di + ijk])
;

dg113
=
0.16666666666666666667*oo2dx*(g13[-2*di + ijk] + 
    8.*(-g13[-di + ijk] + g13[di + ijk]) - g13[2*di + ijk])
;

dg122
=
0.16666666666666666667*oo2dx*(g22[-2*di + ijk] + 
    8.*(-g22[-di + ijk] + g22[di + ijk]) - g22[2*di + ijk])
;

dg123
=
0.16666666666666666667*oo2dx*(g23[-2*di + ijk] + 
    8.*(-g23[-di + ijk] + g23[di + ijk]) - g23[2*di + ijk])
;

dg133
=
0.16666666666666666667*oo2dx*(g33[-2*di + ijk] + 
    8.*(-g33[-di + ijk] + g33[di + ijk]) - g33[2*di + ijk])
;

dg211
=
0.16666666666666666667*oo2dy*(g11[-2*dj + ijk] + 
    8.*(-g11[-dj + ijk] + g11[dj + ijk]) - g11[2*dj + ijk])
;

dg212
=
0.16666666666666666667*oo2dy*(g12[-2*dj + ijk] + 
    8.*(-g12[-dj + ijk] + g12[dj + ijk]) - g12[2*dj + ijk])
;

dg213
=
0.16666666666666666667*oo2dy*(g13[-2*dj + ijk] + 
    8.*(-g13[-dj + ijk] + g13[dj + ijk]) - g13[2*dj + ijk])
;

dg222
=
0.16666666666666666667*oo2dy*(g22[-2*dj + ijk] + 
    8.*(-g22[-dj + ijk] + g22[dj + ijk]) - g22[2*dj + ijk])
;

dg223
=
0.16666666666666666667*oo2dy*(g23[-2*dj + ijk] + 
    8.*(-g23[-dj + ijk] + g23[dj + ijk]) - g23[2*dj + ijk])
;

dg233
=
0.16666666666666666667*oo2dy*(g33[-2*dj + ijk] + 
    8.*(-g33[-dj + ijk] + g33[dj + ijk]) - g33[2*dj + ijk])
;

dg311
=
0.16666666666666666667*oo2dz*(g11[-2*dk + ijk] + 
    8.*(-g11[-dk + ijk] + g11[dk + ijk]) - g11[2*dk + ijk])
;

dg312
=
0.16666666666666666667*oo2dz*(g12[-2*dk + ijk] + 
    8.*(-g12[-dk + ijk] + g12[dk + ijk]) - g12[2*dk + ijk])
;

dg313
=
0.16666666666666666667*oo2dz*(g13[-2*dk + ijk] + 
    8.*(-g13[-dk + ijk] + g13[dk + ijk]) - g13[2*dk + ijk])
;

dg322
=
0.16666666666666666667*oo2dz*(g22[-2*dk + ijk] + 
    8.*(-g22[-dk + ijk] + g22[dk + ijk]) - g22[2*dk + ijk])
;

dg323
=
0.16666666666666666667*oo2dz*(g23[-2*dk + ijk] + 
    8.*(-g23[-dk + ijk] + g23[dk + ijk]) - g23[2*dk + ijk])
;

dg333
=
0.16666666666666666667*oo2dz*(g33[-2*dk + ijk] + 
    8.*(-g33[-dk + ijk] + g33[dk + ijk]) - g33[2*dk + ijk])
;

dK1
=
0.16666666666666666667*oo2dx*(K[-2*di + ijk] + 
    8.*(-K[-di + ijk] + K[di + ijk]) - K[2*di + ijk])
;

dK2
=
0.16666666666666666667*oo2dy*(K[-2*dj + ijk] + 
    8.*(-K[-dj + ijk] + K[dj + ijk]) - K[2*dj + ijk])
;

dK3
=
0.16666666666666666667*oo2dz*(K[-2*dk + ijk] + 
    8.*(-K[-dk + ijk] + K[dk + ijk]) - K[2*dk + ijk])
;

dA111
=
0.16666666666666666667*oo2dx*(A11[-2*di + ijk] + 
    8.*(-A11[-di + ijk] + A11[di + ijk]) - A11[2*di + ijk])
;

dA112
=
0.16666666666666666667*oo2dx*(A12[-2*di + ijk] + 
    8.*(-A12[-di + ijk] + A12[di + ijk]) - A12[2*di + ijk])
;

dA113
=
0.16666666666666666667*oo2dx*(A13[-2*di + ijk] + 
    8.*(-A13[-di + ijk] + A13[di + ijk]) - A13[2*di + ijk])
;

dA122
=
0.16666666666666666667*oo2dx*(A22[-2*di + ijk] + 
    8.*(-A22[-di + ijk] + A22[di + ijk]) - A22[2*di + ijk])
;

dA123
=
0.16666666666666666667*oo2dx*(A23[-2*di + ijk] + 
    8.*(-A23[-di + ijk] + A23[di + ijk]) - A23[2*di + ijk])
;

dA133
=
0.16666666666666666667*oo2dx*(A33[-2*di + ijk] + 
    8.*(-A33[-di + ijk] + A33[di + ijk]) - A33[2*di + ijk])
;

dA211
=
0.16666666666666666667*oo2dy*(A11[-2*dj + ijk] + 
    8.*(-A11[-dj + ijk] + A11[dj + ijk]) - A11[2*dj + ijk])
;

dA212
=
0.16666666666666666667*oo2dy*(A12[-2*dj + ijk] + 
    8.*(-A12[-dj + ijk] + A12[dj + ijk]) - A12[2*dj + ijk])
;

dA213
=
0.16666666666666666667*oo2dy*(A13[-2*dj + ijk] + 
    8.*(-A13[-dj + ijk] + A13[dj + ijk]) - A13[2*dj + ijk])
;

dA222
=
0.16666666666666666667*oo2dy*(A22[-2*dj + ijk] + 
    8.*(-A22[-dj + ijk] + A22[dj + ijk]) - A22[2*dj + ijk])
;

dA223
=
0.16666666666666666667*oo2dy*(A23[-2*dj + ijk] + 
    8.*(-A23[-dj + ijk] + A23[dj + ijk]) - A23[2*dj + ijk])
;

dA233
=
0.16666666666666666667*oo2dy*(A33[-2*dj + ijk] + 
    8.*(-A33[-dj + ijk] + A33[dj + ijk]) - A33[2*dj + ijk])
;

dA311
=
0.16666666666666666667*oo2dz*(A11[-2*dk + ijk] + 
    8.*(-A11[-dk + ijk] + A11[dk + ijk]) - A11[2*dk + ijk])
;

dA312
=
0.16666666666666666667*oo2dz*(A12[-2*dk + ijk] + 
    8.*(-A12[-dk + ijk] + A12[dk + ijk]) - A12[2*dk + ijk])
;

dA313
=
0.16666666666666666667*oo2dz*(A13[-2*dk + ijk] + 
    8.*(-A13[-dk + ijk] + A13[dk + ijk]) - A13[2*dk + ijk])
;

dA322
=
0.16666666666666666667*oo2dz*(A22[-2*dk + ijk] + 
    8.*(-A22[-dk + ijk] + A22[dk + ijk]) - A22[2*dk + ijk])
;

dA323
=
0.16666666666666666667*oo2dz*(A23[-2*dk + ijk] + 
    8.*(-A23[-dk + ijk] + A23[dk + ijk]) - A23[2*dk + ijk])
;

dA333
=
0.16666666666666666667*oo2dz*(A33[-2*dk + ijk] + 
    8.*(-A33[-dk + ijk] + A33[dk + ijk]) - A33[2*dk + ijk])
;

dG11
=
0.16666666666666666667*oo2dx*(G1[-2*di + ijk] + 
    8.*(-G1[-di + ijk] + G1[di + ijk]) - G1[2*di + ijk])
;

dG12
=
0.16666666666666666667*oo2dx*(G2[-2*di + ijk] + 
    8.*(-G2[-di + ijk] + G2[di + ijk]) - G2[2*di + ijk])
;

dG13
=
0.16666666666666666667*oo2dx*(G3[-2*di + ijk] + 
    8.*(-G3[-di + ijk] + G3[di + ijk]) - G3[2*di + ijk])
;

dG21
=
0.16666666666666666667*oo2dy*(G1[-2*dj + ijk] + 
    8.*(-G1[-dj + ijk] + G1[dj + ijk]) - G1[2*dj + ijk])
;

dG22
=
0.16666666666666666667*oo2dy*(G2[-2*dj + ijk] + 
    8.*(-G2[-dj + ijk] + G2[dj + ijk]) - G2[2*dj + ijk])
;

dG23
=
0.16666666666666666667*oo2dy*(G3[-2*dj + ijk] + 
    8.*(-G3[-dj + ijk] + G3[dj + ijk]) - G3[2*dj + ijk])
;

dG31
=
0.16666666666666666667*oo2dz*(G1[-2*dk + ijk] + 
    8.*(-G1[-dk + ijk] + G1[dk + ijk]) - G1[2*dk + ijk])
;

dG32
=
0.16666666666666666667*oo2dz*(G2[-2*dk + ijk] + 
    8.*(-G2[-dk + ijk] + G2[dk + ijk]) - G2[2*dk + ijk])
;

dG33
=
0.16666666666666666667*oo2dz*(G3[-2*dk + ijk] + 
    8.*(-G3[-dk + ijk] + G3[dk + ijk]) - G3[2*dk + ijk])
;


} else if (order == 6 || boundaryNaway(3)) { 

db11
=
0.033333333333333333333*oo2dx*(-beta1[-3*di + ijk] + 
    45.*(-beta1[-di + ijk] + beta1[di + ijk]) + 
    9.*(beta1[-2*di + ijk] - beta1[2*di + ijk]) + beta1[3*di + ijk])
;

db12
=
0.033333333333333333333*oo2dx*(-beta2[-3*di + ijk] + 
    45.*(-beta2[-di + ijk] + beta2[di + ijk]) + 
    9.*(beta2[-2*di + ijk] - beta2[2*di + ijk]) + beta2[3*di + ijk])
;

db13
=
0.033333333333333333333*oo2dx*(-beta3[-3*di + ijk] + 
    45.*(-beta3[-di + ijk] + beta3[di + ijk]) + 
    9.*(beta3[-2*di + ijk] - beta3[2*di + ijk]) + beta3[3*di + ijk])
;

db21
=
0.033333333333333333333*oo2dy*(-beta1[-3*dj + ijk] + 
    45.*(-beta1[-dj + ijk] + beta1[dj + ijk]) + 
    9.*(beta1[-2*dj + ijk] - beta1[2*dj + ijk]) + beta1[3*dj + ijk])
;

db22
=
0.033333333333333333333*oo2dy*(-beta2[-3*dj + ijk] + 
    45.*(-beta2[-dj + ijk] + beta2[dj + ijk]) + 
    9.*(beta2[-2*dj + ijk] - beta2[2*dj + ijk]) + beta2[3*dj + ijk])
;

db23
=
0.033333333333333333333*oo2dy*(-beta3[-3*dj + ijk] + 
    45.*(-beta3[-dj + ijk] + beta3[dj + ijk]) + 
    9.*(beta3[-2*dj + ijk] - beta3[2*dj + ijk]) + beta3[3*dj + ijk])
;

db31
=
0.033333333333333333333*oo2dz*(-beta1[-3*dk + ijk] + 
    45.*(-beta1[-dk + ijk] + beta1[dk + ijk]) + 
    9.*(beta1[-2*dk + ijk] - beta1[2*dk + ijk]) + beta1[3*dk + ijk])
;

db32
=
0.033333333333333333333*oo2dz*(-beta2[-3*dk + ijk] + 
    45.*(-beta2[-dk + ijk] + beta2[dk + ijk]) + 
    9.*(beta2[-2*dk + ijk] - beta2[2*dk + ijk]) + beta2[3*dk + ijk])
;

db33
=
0.033333333333333333333*oo2dz*(-beta3[-3*dk + ijk] + 
    45.*(-beta3[-dk + ijk] + beta3[dk + ijk]) + 
    9.*(beta3[-2*dk + ijk] - beta3[2*dk + ijk]) + beta3[3*dk + ijk])
;

dB11
=
0.033333333333333333333*oo2dx*(-B1[-3*di + ijk] + 
    45.*(-B1[-di + ijk] + B1[di + ijk]) + 
    9.*(B1[-2*di + ijk] - B1[2*di + ijk]) + B1[3*di + ijk])
;

dB12
=
0.033333333333333333333*oo2dx*(-B2[-3*di + ijk] + 
    45.*(-B2[-di + ijk] + B2[di + ijk]) + 
    9.*(B2[-2*di + ijk] - B2[2*di + ijk]) + B2[3*di + ijk])
;

dB13
=
0.033333333333333333333*oo2dx*(-B3[-3*di + ijk] + 
    45.*(-B3[-di + ijk] + B3[di + ijk]) + 
    9.*(B3[-2*di + ijk] - B3[2*di + ijk]) + B3[3*di + ijk])
;

dB21
=
0.033333333333333333333*oo2dy*(-B1[-3*dj + ijk] + 
    45.*(-B1[-dj + ijk] + B1[dj + ijk]) + 
    9.*(B1[-2*dj + ijk] - B1[2*dj + ijk]) + B1[3*dj + ijk])
;

dB22
=
0.033333333333333333333*oo2dy*(-B2[-3*dj + ijk] + 
    45.*(-B2[-dj + ijk] + B2[dj + ijk]) + 
    9.*(B2[-2*dj + ijk] - B2[2*dj + ijk]) + B2[3*dj + ijk])
;

dB23
=
0.033333333333333333333*oo2dy*(-B3[-3*dj + ijk] + 
    45.*(-B3[-dj + ijk] + B3[dj + ijk]) + 
    9.*(B3[-2*dj + ijk] - B3[2*dj + ijk]) + B3[3*dj + ijk])
;

dB31
=
0.033333333333333333333*oo2dz*(-B1[-3*dk + ijk] + 
    45.*(-B1[-dk + ijk] + B1[dk + ijk]) + 
    9.*(B1[-2*dk + ijk] - B1[2*dk + ijk]) + B1[3*dk + ijk])
;

dB32
=
0.033333333333333333333*oo2dz*(-B2[-3*dk + ijk] + 
    45.*(-B2[-dk + ijk] + B2[dk + ijk]) + 
    9.*(B2[-2*dk + ijk] - B2[2*dk + ijk]) + B2[3*dk + ijk])
;

dB33
=
0.033333333333333333333*oo2dz*(-B3[-3*dk + ijk] + 
    45.*(-B3[-dk + ijk] + B3[dk + ijk]) + 
    9.*(B3[-2*dk + ijk] - B3[2*dk + ijk]) + B3[3*dk + ijk])
;

dg111
=
0.033333333333333333333*oo2dx*(-g11[-3*di + ijk] + 
    45.*(-g11[-di + ijk] + g11[di + ijk]) + 
    9.*(g11[-2*di + ijk] - g11[2*di + ijk]) + g11[3*di + ijk])
;

dg112
=
0.033333333333333333333*oo2dx*(-g12[-3*di + ijk] + 
    45.*(-g12[-di + ijk] + g12[di + ijk]) + 
    9.*(g12[-2*di + ijk] - g12[2*di + ijk]) + g12[3*di + ijk])
;

dg113
=
0.033333333333333333333*oo2dx*(-g13[-3*di + ijk] + 
    45.*(-g13[-di + ijk] + g13[di + ijk]) + 
    9.*(g13[-2*di + ijk] - g13[2*di + ijk]) + g13[3*di + ijk])
;

dg122
=
0.033333333333333333333*oo2dx*(-g22[-3*di + ijk] + 
    45.*(-g22[-di + ijk] + g22[di + ijk]) + 
    9.*(g22[-2*di + ijk] - g22[2*di + ijk]) + g22[3*di + ijk])
;

dg123
=
0.033333333333333333333*oo2dx*(-g23[-3*di + ijk] + 
    45.*(-g23[-di + ijk] + g23[di + ijk]) + 
    9.*(g23[-2*di + ijk] - g23[2*di + ijk]) + g23[3*di + ijk])
;

dg133
=
0.033333333333333333333*oo2dx*(-g33[-3*di + ijk] + 
    45.*(-g33[-di + ijk] + g33[di + ijk]) + 
    9.*(g33[-2*di + ijk] - g33[2*di + ijk]) + g33[3*di + ijk])
;

dg211
=
0.033333333333333333333*oo2dy*(-g11[-3*dj + ijk] + 
    45.*(-g11[-dj + ijk] + g11[dj + ijk]) + 
    9.*(g11[-2*dj + ijk] - g11[2*dj + ijk]) + g11[3*dj + ijk])
;

dg212
=
0.033333333333333333333*oo2dy*(-g12[-3*dj + ijk] + 
    45.*(-g12[-dj + ijk] + g12[dj + ijk]) + 
    9.*(g12[-2*dj + ijk] - g12[2*dj + ijk]) + g12[3*dj + ijk])
;

dg213
=
0.033333333333333333333*oo2dy*(-g13[-3*dj + ijk] + 
    45.*(-g13[-dj + ijk] + g13[dj + ijk]) + 
    9.*(g13[-2*dj + ijk] - g13[2*dj + ijk]) + g13[3*dj + ijk])
;

dg222
=
0.033333333333333333333*oo2dy*(-g22[-3*dj + ijk] + 
    45.*(-g22[-dj + ijk] + g22[dj + ijk]) + 
    9.*(g22[-2*dj + ijk] - g22[2*dj + ijk]) + g22[3*dj + ijk])
;

dg223
=
0.033333333333333333333*oo2dy*(-g23[-3*dj + ijk] + 
    45.*(-g23[-dj + ijk] + g23[dj + ijk]) + 
    9.*(g23[-2*dj + ijk] - g23[2*dj + ijk]) + g23[3*dj + ijk])
;

dg233
=
0.033333333333333333333*oo2dy*(-g33[-3*dj + ijk] + 
    45.*(-g33[-dj + ijk] + g33[dj + ijk]) + 
    9.*(g33[-2*dj + ijk] - g33[2*dj + ijk]) + g33[3*dj + ijk])
;

dg311
=
0.033333333333333333333*oo2dz*(-g11[-3*dk + ijk] + 
    45.*(-g11[-dk + ijk] + g11[dk + ijk]) + 
    9.*(g11[-2*dk + ijk] - g11[2*dk + ijk]) + g11[3*dk + ijk])
;

dg312
=
0.033333333333333333333*oo2dz*(-g12[-3*dk + ijk] + 
    45.*(-g12[-dk + ijk] + g12[dk + ijk]) + 
    9.*(g12[-2*dk + ijk] - g12[2*dk + ijk]) + g12[3*dk + ijk])
;

dg313
=
0.033333333333333333333*oo2dz*(-g13[-3*dk + ijk] + 
    45.*(-g13[-dk + ijk] + g13[dk + ijk]) + 
    9.*(g13[-2*dk + ijk] - g13[2*dk + ijk]) + g13[3*dk + ijk])
;

dg322
=
0.033333333333333333333*oo2dz*(-g22[-3*dk + ijk] + 
    45.*(-g22[-dk + ijk] + g22[dk + ijk]) + 
    9.*(g22[-2*dk + ijk] - g22[2*dk + ijk]) + g22[3*dk + ijk])
;

dg323
=
0.033333333333333333333*oo2dz*(-g23[-3*dk + ijk] + 
    45.*(-g23[-dk + ijk] + g23[dk + ijk]) + 
    9.*(g23[-2*dk + ijk] - g23[2*dk + ijk]) + g23[3*dk + ijk])
;

dg333
=
0.033333333333333333333*oo2dz*(-g33[-3*dk + ijk] + 
    45.*(-g33[-dk + ijk] + g33[dk + ijk]) + 
    9.*(g33[-2*dk + ijk] - g33[2*dk + ijk]) + g33[3*dk + ijk])
;

dK1
=
0.033333333333333333333*oo2dx*(-K[-3*di + ijk] + 
    45.*(-K[-di + ijk] + K[di + ijk]) + 
    9.*(K[-2*di + ijk] - K[2*di + ijk]) + K[3*di + ijk])
;

dK2
=
0.033333333333333333333*oo2dy*(-K[-3*dj + ijk] + 
    45.*(-K[-dj + ijk] + K[dj + ijk]) + 
    9.*(K[-2*dj + ijk] - K[2*dj + ijk]) + K[3*dj + ijk])
;

dK3
=
0.033333333333333333333*oo2dz*(-K[-3*dk + ijk] + 
    45.*(-K[-dk + ijk] + K[dk + ijk]) + 
    9.*(K[-2*dk + ijk] - K[2*dk + ijk]) + K[3*dk + ijk])
;

dA111
=
0.033333333333333333333*oo2dx*(-A11[-3*di + ijk] + 
    45.*(-A11[-di + ijk] + A11[di + ijk]) + 
    9.*(A11[-2*di + ijk] - A11[2*di + ijk]) + A11[3*di + ijk])
;

dA112
=
0.033333333333333333333*oo2dx*(-A12[-3*di + ijk] + 
    45.*(-A12[-di + ijk] + A12[di + ijk]) + 
    9.*(A12[-2*di + ijk] - A12[2*di + ijk]) + A12[3*di + ijk])
;

dA113
=
0.033333333333333333333*oo2dx*(-A13[-3*di + ijk] + 
    45.*(-A13[-di + ijk] + A13[di + ijk]) + 
    9.*(A13[-2*di + ijk] - A13[2*di + ijk]) + A13[3*di + ijk])
;

dA122
=
0.033333333333333333333*oo2dx*(-A22[-3*di + ijk] + 
    45.*(-A22[-di + ijk] + A22[di + ijk]) + 
    9.*(A22[-2*di + ijk] - A22[2*di + ijk]) + A22[3*di + ijk])
;

dA123
=
0.033333333333333333333*oo2dx*(-A23[-3*di + ijk] + 
    45.*(-A23[-di + ijk] + A23[di + ijk]) + 
    9.*(A23[-2*di + ijk] - A23[2*di + ijk]) + A23[3*di + ijk])
;

dA133
=
0.033333333333333333333*oo2dx*(-A33[-3*di + ijk] + 
    45.*(-A33[-di + ijk] + A33[di + ijk]) + 
    9.*(A33[-2*di + ijk] - A33[2*di + ijk]) + A33[3*di + ijk])
;

dA211
=
0.033333333333333333333*oo2dy*(-A11[-3*dj + ijk] + 
    45.*(-A11[-dj + ijk] + A11[dj + ijk]) + 
    9.*(A11[-2*dj + ijk] - A11[2*dj + ijk]) + A11[3*dj + ijk])
;

dA212
=
0.033333333333333333333*oo2dy*(-A12[-3*dj + ijk] + 
    45.*(-A12[-dj + ijk] + A12[dj + ijk]) + 
    9.*(A12[-2*dj + ijk] - A12[2*dj + ijk]) + A12[3*dj + ijk])
;

dA213
=
0.033333333333333333333*oo2dy*(-A13[-3*dj + ijk] + 
    45.*(-A13[-dj + ijk] + A13[dj + ijk]) + 
    9.*(A13[-2*dj + ijk] - A13[2*dj + ijk]) + A13[3*dj + ijk])
;

dA222
=
0.033333333333333333333*oo2dy*(-A22[-3*dj + ijk] + 
    45.*(-A22[-dj + ijk] + A22[dj + ijk]) + 
    9.*(A22[-2*dj + ijk] - A22[2*dj + ijk]) + A22[3*dj + ijk])
;

dA223
=
0.033333333333333333333*oo2dy*(-A23[-3*dj + ijk] + 
    45.*(-A23[-dj + ijk] + A23[dj + ijk]) + 
    9.*(A23[-2*dj + ijk] - A23[2*dj + ijk]) + A23[3*dj + ijk])
;

dA233
=
0.033333333333333333333*oo2dy*(-A33[-3*dj + ijk] + 
    45.*(-A33[-dj + ijk] + A33[dj + ijk]) + 
    9.*(A33[-2*dj + ijk] - A33[2*dj + ijk]) + A33[3*dj + ijk])
;

dA311
=
0.033333333333333333333*oo2dz*(-A11[-3*dk + ijk] + 
    45.*(-A11[-dk + ijk] + A11[dk + ijk]) + 
    9.*(A11[-2*dk + ijk] - A11[2*dk + ijk]) + A11[3*dk + ijk])
;

dA312
=
0.033333333333333333333*oo2dz*(-A12[-3*dk + ijk] + 
    45.*(-A12[-dk + ijk] + A12[dk + ijk]) + 
    9.*(A12[-2*dk + ijk] - A12[2*dk + ijk]) + A12[3*dk + ijk])
;

dA313
=
0.033333333333333333333*oo2dz*(-A13[-3*dk + ijk] + 
    45.*(-A13[-dk + ijk] + A13[dk + ijk]) + 
    9.*(A13[-2*dk + ijk] - A13[2*dk + ijk]) + A13[3*dk + ijk])
;

dA322
=
0.033333333333333333333*oo2dz*(-A22[-3*dk + ijk] + 
    45.*(-A22[-dk + ijk] + A22[dk + ijk]) + 
    9.*(A22[-2*dk + ijk] - A22[2*dk + ijk]) + A22[3*dk + ijk])
;

dA323
=
0.033333333333333333333*oo2dz*(-A23[-3*dk + ijk] + 
    45.*(-A23[-dk + ijk] + A23[dk + ijk]) + 
    9.*(A23[-2*dk + ijk] - A23[2*dk + ijk]) + A23[3*dk + ijk])
;

dA333
=
0.033333333333333333333*oo2dz*(-A33[-3*dk + ijk] + 
    45.*(-A33[-dk + ijk] + A33[dk + ijk]) + 
    9.*(A33[-2*dk + ijk] - A33[2*dk + ijk]) + A33[3*dk + ijk])
;

dG11
=
0.033333333333333333333*oo2dx*(-G1[-3*di + ijk] + 
    45.*(-G1[-di + ijk] + G1[di + ijk]) + 
    9.*(G1[-2*di + ijk] - G1[2*di + ijk]) + G1[3*di + ijk])
;

dG12
=
0.033333333333333333333*oo2dx*(-G2[-3*di + ijk] + 
    45.*(-G2[-di + ijk] + G2[di + ijk]) + 
    9.*(G2[-2*di + ijk] - G2[2*di + ijk]) + G2[3*di + ijk])
;

dG13
=
0.033333333333333333333*oo2dx*(-G3[-3*di + ijk] + 
    45.*(-G3[-di + ijk] + G3[di + ijk]) + 
    9.*(G3[-2*di + ijk] - G3[2*di + ijk]) + G3[3*di + ijk])
;

dG21
=
0.033333333333333333333*oo2dy*(-G1[-3*dj + ijk] + 
    45.*(-G1[-dj + ijk] + G1[dj + ijk]) + 
    9.*(G1[-2*dj + ijk] - G1[2*dj + ijk]) + G1[3*dj + ijk])
;

dG22
=
0.033333333333333333333*oo2dy*(-G2[-3*dj + ijk] + 
    45.*(-G2[-dj + ijk] + G2[dj + ijk]) + 
    9.*(G2[-2*dj + ijk] - G2[2*dj + ijk]) + G2[3*dj + ijk])
;

dG23
=
0.033333333333333333333*oo2dy*(-G3[-3*dj + ijk] + 
    45.*(-G3[-dj + ijk] + G3[dj + ijk]) + 
    9.*(G3[-2*dj + ijk] - G3[2*dj + ijk]) + G3[3*dj + ijk])
;

dG31
=
0.033333333333333333333*oo2dz*(-G1[-3*dk + ijk] + 
    45.*(-G1[-dk + ijk] + G1[dk + ijk]) + 
    9.*(G1[-2*dk + ijk] - G1[2*dk + ijk]) + G1[3*dk + ijk])
;

dG32
=
0.033333333333333333333*oo2dz*(-G2[-3*dk + ijk] + 
    45.*(-G2[-dk + ijk] + G2[dk + ijk]) + 
    9.*(G2[-2*dk + ijk] - G2[2*dk + ijk]) + G2[3*dk + ijk])
;

dG33
=
0.033333333333333333333*oo2dz*(-G3[-3*dk + ijk] + 
    45.*(-G3[-dk + ijk] + G3[dk + ijk]) + 
    9.*(G3[-2*dk + ijk] - G3[2*dk + ijk]) + G3[3*dk + ijk])
;


} else if (order == 8 || boundaryNaway(4)) { 

db11
=
0.0023809523809523809524*oo2dx*
  (3.*beta1[-4*di + ijk] + 168.*beta1[-2*di + ijk] + 
    672.*(-beta1[-di + ijk] + beta1[di + ijk]) - 168.*beta1[2*di + ijk] + 
    32.*(-beta1[-3*di + ijk] + beta1[3*di + ijk]) - 3.*beta1[4*di + ijk])
;

db12
=
0.0023809523809523809524*oo2dx*
  (3.*beta2[-4*di + ijk] + 168.*beta2[-2*di + ijk] + 
    672.*(-beta2[-di + ijk] + beta2[di + ijk]) - 168.*beta2[2*di + ijk] + 
    32.*(-beta2[-3*di + ijk] + beta2[3*di + ijk]) - 3.*beta2[4*di + ijk])
;

db13
=
0.0023809523809523809524*oo2dx*
  (3.*beta3[-4*di + ijk] + 168.*beta3[-2*di + ijk] + 
    672.*(-beta3[-di + ijk] + beta3[di + ijk]) - 168.*beta3[2*di + ijk] + 
    32.*(-beta3[-3*di + ijk] + beta3[3*di + ijk]) - 3.*beta3[4*di + ijk])
;

db21
=
0.0023809523809523809524*oo2dy*
  (3.*beta1[-4*dj + ijk] + 168.*beta1[-2*dj + ijk] + 
    672.*(-beta1[-dj + ijk] + beta1[dj + ijk]) - 168.*beta1[2*dj + ijk] + 
    32.*(-beta1[-3*dj + ijk] + beta1[3*dj + ijk]) - 3.*beta1[4*dj + ijk])
;

db22
=
0.0023809523809523809524*oo2dy*
  (3.*beta2[-4*dj + ijk] + 168.*beta2[-2*dj + ijk] + 
    672.*(-beta2[-dj + ijk] + beta2[dj + ijk]) - 168.*beta2[2*dj + ijk] + 
    32.*(-beta2[-3*dj + ijk] + beta2[3*dj + ijk]) - 3.*beta2[4*dj + ijk])
;

db23
=
0.0023809523809523809524*oo2dy*
  (3.*beta3[-4*dj + ijk] + 168.*beta3[-2*dj + ijk] + 
    672.*(-beta3[-dj + ijk] + beta3[dj + ijk]) - 168.*beta3[2*dj + ijk] + 
    32.*(-beta3[-3*dj + ijk] + beta3[3*dj + ijk]) - 3.*beta3[4*dj + ijk])
;

db31
=
0.0023809523809523809524*oo2dz*
  (3.*beta1[-4*dk + ijk] + 168.*beta1[-2*dk + ijk] + 
    672.*(-beta1[-dk + ijk] + beta1[dk + ijk]) - 168.*beta1[2*dk + ijk] + 
    32.*(-beta1[-3*dk + ijk] + beta1[3*dk + ijk]) - 3.*beta1[4*dk + ijk])
;

db32
=
0.0023809523809523809524*oo2dz*
  (3.*beta2[-4*dk + ijk] + 168.*beta2[-2*dk + ijk] + 
    672.*(-beta2[-dk + ijk] + beta2[dk + ijk]) - 168.*beta2[2*dk + ijk] + 
    32.*(-beta2[-3*dk + ijk] + beta2[3*dk + ijk]) - 3.*beta2[4*dk + ijk])
;

db33
=
0.0023809523809523809524*oo2dz*
  (3.*beta3[-4*dk + ijk] + 168.*beta3[-2*dk + ijk] + 
    672.*(-beta3[-dk + ijk] + beta3[dk + ijk]) - 168.*beta3[2*dk + ijk] + 
    32.*(-beta3[-3*dk + ijk] + beta3[3*dk + ijk]) - 3.*beta3[4*dk + ijk])
;

dB11
=
0.0023809523809523809524*oo2dx*
  (3.*B1[-4*di + ijk] + 168.*B1[-2*di + ijk] + 
    672.*(-B1[-di + ijk] + B1[di + ijk]) - 168.*B1[2*di + ijk] + 
    32.*(-B1[-3*di + ijk] + B1[3*di + ijk]) - 3.*B1[4*di + ijk])
;

dB12
=
0.0023809523809523809524*oo2dx*
  (3.*B2[-4*di + ijk] + 168.*B2[-2*di + ijk] + 
    672.*(-B2[-di + ijk] + B2[di + ijk]) - 168.*B2[2*di + ijk] + 
    32.*(-B2[-3*di + ijk] + B2[3*di + ijk]) - 3.*B2[4*di + ijk])
;

dB13
=
0.0023809523809523809524*oo2dx*
  (3.*B3[-4*di + ijk] + 168.*B3[-2*di + ijk] + 
    672.*(-B3[-di + ijk] + B3[di + ijk]) - 168.*B3[2*di + ijk] + 
    32.*(-B3[-3*di + ijk] + B3[3*di + ijk]) - 3.*B3[4*di + ijk])
;

dB21
=
0.0023809523809523809524*oo2dy*
  (3.*B1[-4*dj + ijk] + 168.*B1[-2*dj + ijk] + 
    672.*(-B1[-dj + ijk] + B1[dj + ijk]) - 168.*B1[2*dj + ijk] + 
    32.*(-B1[-3*dj + ijk] + B1[3*dj + ijk]) - 3.*B1[4*dj + ijk])
;

dB22
=
0.0023809523809523809524*oo2dy*
  (3.*B2[-4*dj + ijk] + 168.*B2[-2*dj + ijk] + 
    672.*(-B2[-dj + ijk] + B2[dj + ijk]) - 168.*B2[2*dj + ijk] + 
    32.*(-B2[-3*dj + ijk] + B2[3*dj + ijk]) - 3.*B2[4*dj + ijk])
;

dB23
=
0.0023809523809523809524*oo2dy*
  (3.*B3[-4*dj + ijk] + 168.*B3[-2*dj + ijk] + 
    672.*(-B3[-dj + ijk] + B3[dj + ijk]) - 168.*B3[2*dj + ijk] + 
    32.*(-B3[-3*dj + ijk] + B3[3*dj + ijk]) - 3.*B3[4*dj + ijk])
;

dB31
=
0.0023809523809523809524*oo2dz*
  (3.*B1[-4*dk + ijk] + 168.*B1[-2*dk + ijk] + 
    672.*(-B1[-dk + ijk] + B1[dk + ijk]) - 168.*B1[2*dk + ijk] + 
    32.*(-B1[-3*dk + ijk] + B1[3*dk + ijk]) - 3.*B1[4*dk + ijk])
;

dB32
=
0.0023809523809523809524*oo2dz*
  (3.*B2[-4*dk + ijk] + 168.*B2[-2*dk + ijk] + 
    672.*(-B2[-dk + ijk] + B2[dk + ijk]) - 168.*B2[2*dk + ijk] + 
    32.*(-B2[-3*dk + ijk] + B2[3*dk + ijk]) - 3.*B2[4*dk + ijk])
;

dB33
=
0.0023809523809523809524*oo2dz*
  (3.*B3[-4*dk + ijk] + 168.*B3[-2*dk + ijk] + 
    672.*(-B3[-dk + ijk] + B3[dk + ijk]) - 168.*B3[2*dk + ijk] + 
    32.*(-B3[-3*dk + ijk] + B3[3*dk + ijk]) - 3.*B3[4*dk + ijk])
;

dg111
=
0.0023809523809523809524*oo2dx*
  (3.*g11[-4*di + ijk] + 168.*g11[-2*di + ijk] + 
    672.*(-g11[-di + ijk] + g11[di + ijk]) - 168.*g11[2*di + ijk] + 
    32.*(-g11[-3*di + ijk] + g11[3*di + ijk]) - 3.*g11[4*di + ijk])
;

dg112
=
0.0023809523809523809524*oo2dx*
  (3.*g12[-4*di + ijk] + 168.*g12[-2*di + ijk] + 
    672.*(-g12[-di + ijk] + g12[di + ijk]) - 168.*g12[2*di + ijk] + 
    32.*(-g12[-3*di + ijk] + g12[3*di + ijk]) - 3.*g12[4*di + ijk])
;

dg113
=
0.0023809523809523809524*oo2dx*
  (3.*g13[-4*di + ijk] + 168.*g13[-2*di + ijk] + 
    672.*(-g13[-di + ijk] + g13[di + ijk]) - 168.*g13[2*di + ijk] + 
    32.*(-g13[-3*di + ijk] + g13[3*di + ijk]) - 3.*g13[4*di + ijk])
;

dg122
=
0.0023809523809523809524*oo2dx*
  (3.*g22[-4*di + ijk] + 168.*g22[-2*di + ijk] + 
    672.*(-g22[-di + ijk] + g22[di + ijk]) - 168.*g22[2*di + ijk] + 
    32.*(-g22[-3*di + ijk] + g22[3*di + ijk]) - 3.*g22[4*di + ijk])
;

dg123
=
0.0023809523809523809524*oo2dx*
  (3.*g23[-4*di + ijk] + 168.*g23[-2*di + ijk] + 
    672.*(-g23[-di + ijk] + g23[di + ijk]) - 168.*g23[2*di + ijk] + 
    32.*(-g23[-3*di + ijk] + g23[3*di + ijk]) - 3.*g23[4*di + ijk])
;

dg133
=
0.0023809523809523809524*oo2dx*
  (3.*g33[-4*di + ijk] + 168.*g33[-2*di + ijk] + 
    672.*(-g33[-di + ijk] + g33[di + ijk]) - 168.*g33[2*di + ijk] + 
    32.*(-g33[-3*di + ijk] + g33[3*di + ijk]) - 3.*g33[4*di + ijk])
;

dg211
=
0.0023809523809523809524*oo2dy*
  (3.*g11[-4*dj + ijk] + 168.*g11[-2*dj + ijk] + 
    672.*(-g11[-dj + ijk] + g11[dj + ijk]) - 168.*g11[2*dj + ijk] + 
    32.*(-g11[-3*dj + ijk] + g11[3*dj + ijk]) - 3.*g11[4*dj + ijk])
;

dg212
=
0.0023809523809523809524*oo2dy*
  (3.*g12[-4*dj + ijk] + 168.*g12[-2*dj + ijk] + 
    672.*(-g12[-dj + ijk] + g12[dj + ijk]) - 168.*g12[2*dj + ijk] + 
    32.*(-g12[-3*dj + ijk] + g12[3*dj + ijk]) - 3.*g12[4*dj + ijk])
;

dg213
=
0.0023809523809523809524*oo2dy*
  (3.*g13[-4*dj + ijk] + 168.*g13[-2*dj + ijk] + 
    672.*(-g13[-dj + ijk] + g13[dj + ijk]) - 168.*g13[2*dj + ijk] + 
    32.*(-g13[-3*dj + ijk] + g13[3*dj + ijk]) - 3.*g13[4*dj + ijk])
;

dg222
=
0.0023809523809523809524*oo2dy*
  (3.*g22[-4*dj + ijk] + 168.*g22[-2*dj + ijk] + 
    672.*(-g22[-dj + ijk] + g22[dj + ijk]) - 168.*g22[2*dj + ijk] + 
    32.*(-g22[-3*dj + ijk] + g22[3*dj + ijk]) - 3.*g22[4*dj + ijk])
;

dg223
=
0.0023809523809523809524*oo2dy*
  (3.*g23[-4*dj + ijk] + 168.*g23[-2*dj + ijk] + 
    672.*(-g23[-dj + ijk] + g23[dj + ijk]) - 168.*g23[2*dj + ijk] + 
    32.*(-g23[-3*dj + ijk] + g23[3*dj + ijk]) - 3.*g23[4*dj + ijk])
;

dg233
=
0.0023809523809523809524*oo2dy*
  (3.*g33[-4*dj + ijk] + 168.*g33[-2*dj + ijk] + 
    672.*(-g33[-dj + ijk] + g33[dj + ijk]) - 168.*g33[2*dj + ijk] + 
    32.*(-g33[-3*dj + ijk] + g33[3*dj + ijk]) - 3.*g33[4*dj + ijk])
;

dg311
=
0.0023809523809523809524*oo2dz*
  (3.*g11[-4*dk + ijk] + 168.*g11[-2*dk + ijk] + 
    672.*(-g11[-dk + ijk] + g11[dk + ijk]) - 168.*g11[2*dk + ijk] + 
    32.*(-g11[-3*dk + ijk] + g11[3*dk + ijk]) - 3.*g11[4*dk + ijk])
;

dg312
=
0.0023809523809523809524*oo2dz*
  (3.*g12[-4*dk + ijk] + 168.*g12[-2*dk + ijk] + 
    672.*(-g12[-dk + ijk] + g12[dk + ijk]) - 168.*g12[2*dk + ijk] + 
    32.*(-g12[-3*dk + ijk] + g12[3*dk + ijk]) - 3.*g12[4*dk + ijk])
;

dg313
=
0.0023809523809523809524*oo2dz*
  (3.*g13[-4*dk + ijk] + 168.*g13[-2*dk + ijk] + 
    672.*(-g13[-dk + ijk] + g13[dk + ijk]) - 168.*g13[2*dk + ijk] + 
    32.*(-g13[-3*dk + ijk] + g13[3*dk + ijk]) - 3.*g13[4*dk + ijk])
;

dg322
=
0.0023809523809523809524*oo2dz*
  (3.*g22[-4*dk + ijk] + 168.*g22[-2*dk + ijk] + 
    672.*(-g22[-dk + ijk] + g22[dk + ijk]) - 168.*g22[2*dk + ijk] + 
    32.*(-g22[-3*dk + ijk] + g22[3*dk + ijk]) - 3.*g22[4*dk + ijk])
;

dg323
=
0.0023809523809523809524*oo2dz*
  (3.*g23[-4*dk + ijk] + 168.*g23[-2*dk + ijk] + 
    672.*(-g23[-dk + ijk] + g23[dk + ijk]) - 168.*g23[2*dk + ijk] + 
    32.*(-g23[-3*dk + ijk] + g23[3*dk + ijk]) - 3.*g23[4*dk + ijk])
;

dg333
=
0.0023809523809523809524*oo2dz*
  (3.*g33[-4*dk + ijk] + 168.*g33[-2*dk + ijk] + 
    672.*(-g33[-dk + ijk] + g33[dk + ijk]) - 168.*g33[2*dk + ijk] + 
    32.*(-g33[-3*dk + ijk] + g33[3*dk + ijk]) - 3.*g33[4*dk + ijk])
;

dK1
=
0.0023809523809523809524*oo2dx*
  (3.*K[-4*di + ijk] + 168.*K[-2*di + ijk] + 
    672.*(-K[-di + ijk] + K[di + ijk]) - 168.*K[2*di + ijk] + 
    32.*(-K[-3*di + ijk] + K[3*di + ijk]) - 3.*K[4*di + ijk])
;

dK2
=
0.0023809523809523809524*oo2dy*
  (3.*K[-4*dj + ijk] + 168.*K[-2*dj + ijk] + 
    672.*(-K[-dj + ijk] + K[dj + ijk]) - 168.*K[2*dj + ijk] + 
    32.*(-K[-3*dj + ijk] + K[3*dj + ijk]) - 3.*K[4*dj + ijk])
;

dK3
=
0.0023809523809523809524*oo2dz*
  (3.*K[-4*dk + ijk] + 168.*K[-2*dk + ijk] + 
    672.*(-K[-dk + ijk] + K[dk + ijk]) - 168.*K[2*dk + ijk] + 
    32.*(-K[-3*dk + ijk] + K[3*dk + ijk]) - 3.*K[4*dk + ijk])
;

dA111
=
0.0023809523809523809524*oo2dx*
  (3.*A11[-4*di + ijk] + 168.*A11[-2*di + ijk] + 
    672.*(-A11[-di + ijk] + A11[di + ijk]) - 168.*A11[2*di + ijk] + 
    32.*(-A11[-3*di + ijk] + A11[3*di + ijk]) - 3.*A11[4*di + ijk])
;

dA112
=
0.0023809523809523809524*oo2dx*
  (3.*A12[-4*di + ijk] + 168.*A12[-2*di + ijk] + 
    672.*(-A12[-di + ijk] + A12[di + ijk]) - 168.*A12[2*di + ijk] + 
    32.*(-A12[-3*di + ijk] + A12[3*di + ijk]) - 3.*A12[4*di + ijk])
;

dA113
=
0.0023809523809523809524*oo2dx*
  (3.*A13[-4*di + ijk] + 168.*A13[-2*di + ijk] + 
    672.*(-A13[-di + ijk] + A13[di + ijk]) - 168.*A13[2*di + ijk] + 
    32.*(-A13[-3*di + ijk] + A13[3*di + ijk]) - 3.*A13[4*di + ijk])
;

dA122
=
0.0023809523809523809524*oo2dx*
  (3.*A22[-4*di + ijk] + 168.*A22[-2*di + ijk] + 
    672.*(-A22[-di + ijk] + A22[di + ijk]) - 168.*A22[2*di + ijk] + 
    32.*(-A22[-3*di + ijk] + A22[3*di + ijk]) - 3.*A22[4*di + ijk])
;

dA123
=
0.0023809523809523809524*oo2dx*
  (3.*A23[-4*di + ijk] + 168.*A23[-2*di + ijk] + 
    672.*(-A23[-di + ijk] + A23[di + ijk]) - 168.*A23[2*di + ijk] + 
    32.*(-A23[-3*di + ijk] + A23[3*di + ijk]) - 3.*A23[4*di + ijk])
;

dA133
=
0.0023809523809523809524*oo2dx*
  (3.*A33[-4*di + ijk] + 168.*A33[-2*di + ijk] + 
    672.*(-A33[-di + ijk] + A33[di + ijk]) - 168.*A33[2*di + ijk] + 
    32.*(-A33[-3*di + ijk] + A33[3*di + ijk]) - 3.*A33[4*di + ijk])
;

dA211
=
0.0023809523809523809524*oo2dy*
  (3.*A11[-4*dj + ijk] + 168.*A11[-2*dj + ijk] + 
    672.*(-A11[-dj + ijk] + A11[dj + ijk]) - 168.*A11[2*dj + ijk] + 
    32.*(-A11[-3*dj + ijk] + A11[3*dj + ijk]) - 3.*A11[4*dj + ijk])
;

dA212
=
0.0023809523809523809524*oo2dy*
  (3.*A12[-4*dj + ijk] + 168.*A12[-2*dj + ijk] + 
    672.*(-A12[-dj + ijk] + A12[dj + ijk]) - 168.*A12[2*dj + ijk] + 
    32.*(-A12[-3*dj + ijk] + A12[3*dj + ijk]) - 3.*A12[4*dj + ijk])
;

dA213
=
0.0023809523809523809524*oo2dy*
  (3.*A13[-4*dj + ijk] + 168.*A13[-2*dj + ijk] + 
    672.*(-A13[-dj + ijk] + A13[dj + ijk]) - 168.*A13[2*dj + ijk] + 
    32.*(-A13[-3*dj + ijk] + A13[3*dj + ijk]) - 3.*A13[4*dj + ijk])
;

dA222
=
0.0023809523809523809524*oo2dy*
  (3.*A22[-4*dj + ijk] + 168.*A22[-2*dj + ijk] + 
    672.*(-A22[-dj + ijk] + A22[dj + ijk]) - 168.*A22[2*dj + ijk] + 
    32.*(-A22[-3*dj + ijk] + A22[3*dj + ijk]) - 3.*A22[4*dj + ijk])
;

dA223
=
0.0023809523809523809524*oo2dy*
  (3.*A23[-4*dj + ijk] + 168.*A23[-2*dj + ijk] + 
    672.*(-A23[-dj + ijk] + A23[dj + ijk]) - 168.*A23[2*dj + ijk] + 
    32.*(-A23[-3*dj + ijk] + A23[3*dj + ijk]) - 3.*A23[4*dj + ijk])
;

dA233
=
0.0023809523809523809524*oo2dy*
  (3.*A33[-4*dj + ijk] + 168.*A33[-2*dj + ijk] + 
    672.*(-A33[-dj + ijk] + A33[dj + ijk]) - 168.*A33[2*dj + ijk] + 
    32.*(-A33[-3*dj + ijk] + A33[3*dj + ijk]) - 3.*A33[4*dj + ijk])
;

dA311
=
0.0023809523809523809524*oo2dz*
  (3.*A11[-4*dk + ijk] + 168.*A11[-2*dk + ijk] + 
    672.*(-A11[-dk + ijk] + A11[dk + ijk]) - 168.*A11[2*dk + ijk] + 
    32.*(-A11[-3*dk + ijk] + A11[3*dk + ijk]) - 3.*A11[4*dk + ijk])
;

dA312
=
0.0023809523809523809524*oo2dz*
  (3.*A12[-4*dk + ijk] + 168.*A12[-2*dk + ijk] + 
    672.*(-A12[-dk + ijk] + A12[dk + ijk]) - 168.*A12[2*dk + ijk] + 
    32.*(-A12[-3*dk + ijk] + A12[3*dk + ijk]) - 3.*A12[4*dk + ijk])
;

dA313
=
0.0023809523809523809524*oo2dz*
  (3.*A13[-4*dk + ijk] + 168.*A13[-2*dk + ijk] + 
    672.*(-A13[-dk + ijk] + A13[dk + ijk]) - 168.*A13[2*dk + ijk] + 
    32.*(-A13[-3*dk + ijk] + A13[3*dk + ijk]) - 3.*A13[4*dk + ijk])
;

dA322
=
0.0023809523809523809524*oo2dz*
  (3.*A22[-4*dk + ijk] + 168.*A22[-2*dk + ijk] + 
    672.*(-A22[-dk + ijk] + A22[dk + ijk]) - 168.*A22[2*dk + ijk] + 
    32.*(-A22[-3*dk + ijk] + A22[3*dk + ijk]) - 3.*A22[4*dk + ijk])
;

dA323
=
0.0023809523809523809524*oo2dz*
  (3.*A23[-4*dk + ijk] + 168.*A23[-2*dk + ijk] + 
    672.*(-A23[-dk + ijk] + A23[dk + ijk]) - 168.*A23[2*dk + ijk] + 
    32.*(-A23[-3*dk + ijk] + A23[3*dk + ijk]) - 3.*A23[4*dk + ijk])
;

dA333
=
0.0023809523809523809524*oo2dz*
  (3.*A33[-4*dk + ijk] + 168.*A33[-2*dk + ijk] + 
    672.*(-A33[-dk + ijk] + A33[dk + ijk]) - 168.*A33[2*dk + ijk] + 
    32.*(-A33[-3*dk + ijk] + A33[3*dk + ijk]) - 3.*A33[4*dk + ijk])
;

dG11
=
0.0023809523809523809524*oo2dx*
  (3.*G1[-4*di + ijk] + 168.*G1[-2*di + ijk] + 
    672.*(-G1[-di + ijk] + G1[di + ijk]) - 168.*G1[2*di + ijk] + 
    32.*(-G1[-3*di + ijk] + G1[3*di + ijk]) - 3.*G1[4*di + ijk])
;

dG12
=
0.0023809523809523809524*oo2dx*
  (3.*G2[-4*di + ijk] + 168.*G2[-2*di + ijk] + 
    672.*(-G2[-di + ijk] + G2[di + ijk]) - 168.*G2[2*di + ijk] + 
    32.*(-G2[-3*di + ijk] + G2[3*di + ijk]) - 3.*G2[4*di + ijk])
;

dG13
=
0.0023809523809523809524*oo2dx*
  (3.*G3[-4*di + ijk] + 168.*G3[-2*di + ijk] + 
    672.*(-G3[-di + ijk] + G3[di + ijk]) - 168.*G3[2*di + ijk] + 
    32.*(-G3[-3*di + ijk] + G3[3*di + ijk]) - 3.*G3[4*di + ijk])
;

dG21
=
0.0023809523809523809524*oo2dy*
  (3.*G1[-4*dj + ijk] + 168.*G1[-2*dj + ijk] + 
    672.*(-G1[-dj + ijk] + G1[dj + ijk]) - 168.*G1[2*dj + ijk] + 
    32.*(-G1[-3*dj + ijk] + G1[3*dj + ijk]) - 3.*G1[4*dj + ijk])
;

dG22
=
0.0023809523809523809524*oo2dy*
  (3.*G2[-4*dj + ijk] + 168.*G2[-2*dj + ijk] + 
    672.*(-G2[-dj + ijk] + G2[dj + ijk]) - 168.*G2[2*dj + ijk] + 
    32.*(-G2[-3*dj + ijk] + G2[3*dj + ijk]) - 3.*G2[4*dj + ijk])
;

dG23
=
0.0023809523809523809524*oo2dy*
  (3.*G3[-4*dj + ijk] + 168.*G3[-2*dj + ijk] + 
    672.*(-G3[-dj + ijk] + G3[dj + ijk]) - 168.*G3[2*dj + ijk] + 
    32.*(-G3[-3*dj + ijk] + G3[3*dj + ijk]) - 3.*G3[4*dj + ijk])
;

dG31
=
0.0023809523809523809524*oo2dz*
  (3.*G1[-4*dk + ijk] + 168.*G1[-2*dk + ijk] + 
    672.*(-G1[-dk + ijk] + G1[dk + ijk]) - 168.*G1[2*dk + ijk] + 
    32.*(-G1[-3*dk + ijk] + G1[3*dk + ijk]) - 3.*G1[4*dk + ijk])
;

dG32
=
0.0023809523809523809524*oo2dz*
  (3.*G2[-4*dk + ijk] + 168.*G2[-2*dk + ijk] + 
    672.*(-G2[-dk + ijk] + G2[dk + ijk]) - 168.*G2[2*dk + ijk] + 
    32.*(-G2[-3*dk + ijk] + G2[3*dk + ijk]) - 3.*G2[4*dk + ijk])
;

dG33
=
0.0023809523809523809524*oo2dz*
  (3.*G3[-4*dk + ijk] + 168.*G3[-2*dk + ijk] + 
    672.*(-G3[-dk + ijk] + G3[dk + ijk]) - 168.*G3[2*dk + ijk] + 
    32.*(-G3[-3*dk + ijk] + G3[3*dk + ijk]) - 3.*G3[4*dk + ijk])
;


} else if (order == 10 || boundaryNaway(5)) { 

db11
=
0.00079365079365079365079*oo2dx*
  (-2.*beta1[-5*di + ijk] + 25.*beta1[-4*di + ijk] + 
    600.*beta1[-2*di + ijk] + 2100.*(-beta1[-di + ijk] + beta1[di + ijk]) - 
    600.*beta1[2*di + ijk] + 150.*
     (-beta1[-3*di + ijk] + beta1[3*di + ijk]) - 25.*beta1[4*di + ijk] + 
    2.*beta1[5*di + ijk])
;

db12
=
0.00079365079365079365079*oo2dx*
  (-2.*beta2[-5*di + ijk] + 25.*beta2[-4*di + ijk] + 
    600.*beta2[-2*di + ijk] + 2100.*(-beta2[-di + ijk] + beta2[di + ijk]) - 
    600.*beta2[2*di + ijk] + 150.*
     (-beta2[-3*di + ijk] + beta2[3*di + ijk]) - 25.*beta2[4*di + ijk] + 
    2.*beta2[5*di + ijk])
;

db13
=
0.00079365079365079365079*oo2dx*
  (-2.*beta3[-5*di + ijk] + 25.*beta3[-4*di + ijk] + 
    600.*beta3[-2*di + ijk] + 2100.*(-beta3[-di + ijk] + beta3[di + ijk]) - 
    600.*beta3[2*di + ijk] + 150.*
     (-beta3[-3*di + ijk] + beta3[3*di + ijk]) - 25.*beta3[4*di + ijk] + 
    2.*beta3[5*di + ijk])
;

db21
=
0.00079365079365079365079*oo2dy*
  (-2.*beta1[-5*dj + ijk] + 25.*beta1[-4*dj + ijk] + 
    600.*beta1[-2*dj + ijk] + 2100.*(-beta1[-dj + ijk] + beta1[dj + ijk]) - 
    600.*beta1[2*dj + ijk] + 150.*
     (-beta1[-3*dj + ijk] + beta1[3*dj + ijk]) - 25.*beta1[4*dj + ijk] + 
    2.*beta1[5*dj + ijk])
;

db22
=
0.00079365079365079365079*oo2dy*
  (-2.*beta2[-5*dj + ijk] + 25.*beta2[-4*dj + ijk] + 
    600.*beta2[-2*dj + ijk] + 2100.*(-beta2[-dj + ijk] + beta2[dj + ijk]) - 
    600.*beta2[2*dj + ijk] + 150.*
     (-beta2[-3*dj + ijk] + beta2[3*dj + ijk]) - 25.*beta2[4*dj + ijk] + 
    2.*beta2[5*dj + ijk])
;

db23
=
0.00079365079365079365079*oo2dy*
  (-2.*beta3[-5*dj + ijk] + 25.*beta3[-4*dj + ijk] + 
    600.*beta3[-2*dj + ijk] + 2100.*(-beta3[-dj + ijk] + beta3[dj + ijk]) - 
    600.*beta3[2*dj + ijk] + 150.*
     (-beta3[-3*dj + ijk] + beta3[3*dj + ijk]) - 25.*beta3[4*dj + ijk] + 
    2.*beta3[5*dj + ijk])
;

db31
=
0.00079365079365079365079*oo2dz*
  (-2.*beta1[-5*dk + ijk] + 25.*beta1[-4*dk + ijk] + 
    600.*beta1[-2*dk + ijk] + 2100.*(-beta1[-dk + ijk] + beta1[dk + ijk]) - 
    600.*beta1[2*dk + ijk] + 150.*
     (-beta1[-3*dk + ijk] + beta1[3*dk + ijk]) - 25.*beta1[4*dk + ijk] + 
    2.*beta1[5*dk + ijk])
;

db32
=
0.00079365079365079365079*oo2dz*
  (-2.*beta2[-5*dk + ijk] + 25.*beta2[-4*dk + ijk] + 
    600.*beta2[-2*dk + ijk] + 2100.*(-beta2[-dk + ijk] + beta2[dk + ijk]) - 
    600.*beta2[2*dk + ijk] + 150.*
     (-beta2[-3*dk + ijk] + beta2[3*dk + ijk]) - 25.*beta2[4*dk + ijk] + 
    2.*beta2[5*dk + ijk])
;

db33
=
0.00079365079365079365079*oo2dz*
  (-2.*beta3[-5*dk + ijk] + 25.*beta3[-4*dk + ijk] + 
    600.*beta3[-2*dk + ijk] + 2100.*(-beta3[-dk + ijk] + beta3[dk + ijk]) - 
    600.*beta3[2*dk + ijk] + 150.*
     (-beta3[-3*dk + ijk] + beta3[3*dk + ijk]) - 25.*beta3[4*dk + ijk] + 
    2.*beta3[5*dk + ijk])
;

dB11
=
0.00079365079365079365079*oo2dx*
  (-2.*B1[-5*di + ijk] + 25.*B1[-4*di + ijk] + 600.*B1[-2*di + ijk] + 
    2100.*(-B1[-di + ijk] + B1[di + ijk]) - 600.*B1[2*di + ijk] + 
    150.*(-B1[-3*di + ijk] + B1[3*di + ijk]) - 25.*B1[4*di + ijk] + 
    2.*B1[5*di + ijk])
;

dB12
=
0.00079365079365079365079*oo2dx*
  (-2.*B2[-5*di + ijk] + 25.*B2[-4*di + ijk] + 600.*B2[-2*di + ijk] + 
    2100.*(-B2[-di + ijk] + B2[di + ijk]) - 600.*B2[2*di + ijk] + 
    150.*(-B2[-3*di + ijk] + B2[3*di + ijk]) - 25.*B2[4*di + ijk] + 
    2.*B2[5*di + ijk])
;

dB13
=
0.00079365079365079365079*oo2dx*
  (-2.*B3[-5*di + ijk] + 25.*B3[-4*di + ijk] + 600.*B3[-2*di + ijk] + 
    2100.*(-B3[-di + ijk] + B3[di + ijk]) - 600.*B3[2*di + ijk] + 
    150.*(-B3[-3*di + ijk] + B3[3*di + ijk]) - 25.*B3[4*di + ijk] + 
    2.*B3[5*di + ijk])
;

dB21
=
0.00079365079365079365079*oo2dy*
  (-2.*B1[-5*dj + ijk] + 25.*B1[-4*dj + ijk] + 600.*B1[-2*dj + ijk] + 
    2100.*(-B1[-dj + ijk] + B1[dj + ijk]) - 600.*B1[2*dj + ijk] + 
    150.*(-B1[-3*dj + ijk] + B1[3*dj + ijk]) - 25.*B1[4*dj + ijk] + 
    2.*B1[5*dj + ijk])
;

dB22
=
0.00079365079365079365079*oo2dy*
  (-2.*B2[-5*dj + ijk] + 25.*B2[-4*dj + ijk] + 600.*B2[-2*dj + ijk] + 
    2100.*(-B2[-dj + ijk] + B2[dj + ijk]) - 600.*B2[2*dj + ijk] + 
    150.*(-B2[-3*dj + ijk] + B2[3*dj + ijk]) - 25.*B2[4*dj + ijk] + 
    2.*B2[5*dj + ijk])
;

dB23
=
0.00079365079365079365079*oo2dy*
  (-2.*B3[-5*dj + ijk] + 25.*B3[-4*dj + ijk] + 600.*B3[-2*dj + ijk] + 
    2100.*(-B3[-dj + ijk] + B3[dj + ijk]) - 600.*B3[2*dj + ijk] + 
    150.*(-B3[-3*dj + ijk] + B3[3*dj + ijk]) - 25.*B3[4*dj + ijk] + 
    2.*B3[5*dj + ijk])
;

dB31
=
0.00079365079365079365079*oo2dz*
  (-2.*B1[-5*dk + ijk] + 25.*B1[-4*dk + ijk] + 600.*B1[-2*dk + ijk] + 
    2100.*(-B1[-dk + ijk] + B1[dk + ijk]) - 600.*B1[2*dk + ijk] + 
    150.*(-B1[-3*dk + ijk] + B1[3*dk + ijk]) - 25.*B1[4*dk + ijk] + 
    2.*B1[5*dk + ijk])
;

dB32
=
0.00079365079365079365079*oo2dz*
  (-2.*B2[-5*dk + ijk] + 25.*B2[-4*dk + ijk] + 600.*B2[-2*dk + ijk] + 
    2100.*(-B2[-dk + ijk] + B2[dk + ijk]) - 600.*B2[2*dk + ijk] + 
    150.*(-B2[-3*dk + ijk] + B2[3*dk + ijk]) - 25.*B2[4*dk + ijk] + 
    2.*B2[5*dk + ijk])
;

dB33
=
0.00079365079365079365079*oo2dz*
  (-2.*B3[-5*dk + ijk] + 25.*B3[-4*dk + ijk] + 600.*B3[-2*dk + ijk] + 
    2100.*(-B3[-dk + ijk] + B3[dk + ijk]) - 600.*B3[2*dk + ijk] + 
    150.*(-B3[-3*dk + ijk] + B3[3*dk + ijk]) - 25.*B3[4*dk + ijk] + 
    2.*B3[5*dk + ijk])
;

dg111
=
0.00079365079365079365079*oo2dx*
  (-2.*g11[-5*di + ijk] + 25.*g11[-4*di + ijk] + 600.*g11[-2*di + ijk] + 
    2100.*(-g11[-di + ijk] + g11[di + ijk]) - 600.*g11[2*di + ijk] + 
    150.*(-g11[-3*di + ijk] + g11[3*di + ijk]) - 25.*g11[4*di + ijk] + 
    2.*g11[5*di + ijk])
;

dg112
=
0.00079365079365079365079*oo2dx*
  (-2.*g12[-5*di + ijk] + 25.*g12[-4*di + ijk] + 600.*g12[-2*di + ijk] + 
    2100.*(-g12[-di + ijk] + g12[di + ijk]) - 600.*g12[2*di + ijk] + 
    150.*(-g12[-3*di + ijk] + g12[3*di + ijk]) - 25.*g12[4*di + ijk] + 
    2.*g12[5*di + ijk])
;

dg113
=
0.00079365079365079365079*oo2dx*
  (-2.*g13[-5*di + ijk] + 25.*g13[-4*di + ijk] + 600.*g13[-2*di + ijk] + 
    2100.*(-g13[-di + ijk] + g13[di + ijk]) - 600.*g13[2*di + ijk] + 
    150.*(-g13[-3*di + ijk] + g13[3*di + ijk]) - 25.*g13[4*di + ijk] + 
    2.*g13[5*di + ijk])
;

dg122
=
0.00079365079365079365079*oo2dx*
  (-2.*g22[-5*di + ijk] + 25.*g22[-4*di + ijk] + 600.*g22[-2*di + ijk] + 
    2100.*(-g22[-di + ijk] + g22[di + ijk]) - 600.*g22[2*di + ijk] + 
    150.*(-g22[-3*di + ijk] + g22[3*di + ijk]) - 25.*g22[4*di + ijk] + 
    2.*g22[5*di + ijk])
;

dg123
=
0.00079365079365079365079*oo2dx*
  (-2.*g23[-5*di + ijk] + 25.*g23[-4*di + ijk] + 600.*g23[-2*di + ijk] + 
    2100.*(-g23[-di + ijk] + g23[di + ijk]) - 600.*g23[2*di + ijk] + 
    150.*(-g23[-3*di + ijk] + g23[3*di + ijk]) - 25.*g23[4*di + ijk] + 
    2.*g23[5*di + ijk])
;

dg133
=
0.00079365079365079365079*oo2dx*
  (-2.*g33[-5*di + ijk] + 25.*g33[-4*di + ijk] + 600.*g33[-2*di + ijk] + 
    2100.*(-g33[-di + ijk] + g33[di + ijk]) - 600.*g33[2*di + ijk] + 
    150.*(-g33[-3*di + ijk] + g33[3*di + ijk]) - 25.*g33[4*di + ijk] + 
    2.*g33[5*di + ijk])
;

dg211
=
0.00079365079365079365079*oo2dy*
  (-2.*g11[-5*dj + ijk] + 25.*g11[-4*dj + ijk] + 600.*g11[-2*dj + ijk] + 
    2100.*(-g11[-dj + ijk] + g11[dj + ijk]) - 600.*g11[2*dj + ijk] + 
    150.*(-g11[-3*dj + ijk] + g11[3*dj + ijk]) - 25.*g11[4*dj + ijk] + 
    2.*g11[5*dj + ijk])
;

dg212
=
0.00079365079365079365079*oo2dy*
  (-2.*g12[-5*dj + ijk] + 25.*g12[-4*dj + ijk] + 600.*g12[-2*dj + ijk] + 
    2100.*(-g12[-dj + ijk] + g12[dj + ijk]) - 600.*g12[2*dj + ijk] + 
    150.*(-g12[-3*dj + ijk] + g12[3*dj + ijk]) - 25.*g12[4*dj + ijk] + 
    2.*g12[5*dj + ijk])
;

dg213
=
0.00079365079365079365079*oo2dy*
  (-2.*g13[-5*dj + ijk] + 25.*g13[-4*dj + ijk] + 600.*g13[-2*dj + ijk] + 
    2100.*(-g13[-dj + ijk] + g13[dj + ijk]) - 600.*g13[2*dj + ijk] + 
    150.*(-g13[-3*dj + ijk] + g13[3*dj + ijk]) - 25.*g13[4*dj + ijk] + 
    2.*g13[5*dj + ijk])
;

dg222
=
0.00079365079365079365079*oo2dy*
  (-2.*g22[-5*dj + ijk] + 25.*g22[-4*dj + ijk] + 600.*g22[-2*dj + ijk] + 
    2100.*(-g22[-dj + ijk] + g22[dj + ijk]) - 600.*g22[2*dj + ijk] + 
    150.*(-g22[-3*dj + ijk] + g22[3*dj + ijk]) - 25.*g22[4*dj + ijk] + 
    2.*g22[5*dj + ijk])
;

dg223
=
0.00079365079365079365079*oo2dy*
  (-2.*g23[-5*dj + ijk] + 25.*g23[-4*dj + ijk] + 600.*g23[-2*dj + ijk] + 
    2100.*(-g23[-dj + ijk] + g23[dj + ijk]) - 600.*g23[2*dj + ijk] + 
    150.*(-g23[-3*dj + ijk] + g23[3*dj + ijk]) - 25.*g23[4*dj + ijk] + 
    2.*g23[5*dj + ijk])
;

dg233
=
0.00079365079365079365079*oo2dy*
  (-2.*g33[-5*dj + ijk] + 25.*g33[-4*dj + ijk] + 600.*g33[-2*dj + ijk] + 
    2100.*(-g33[-dj + ijk] + g33[dj + ijk]) - 600.*g33[2*dj + ijk] + 
    150.*(-g33[-3*dj + ijk] + g33[3*dj + ijk]) - 25.*g33[4*dj + ijk] + 
    2.*g33[5*dj + ijk])
;

dg311
=
0.00079365079365079365079*oo2dz*
  (-2.*g11[-5*dk + ijk] + 25.*g11[-4*dk + ijk] + 600.*g11[-2*dk + ijk] + 
    2100.*(-g11[-dk + ijk] + g11[dk + ijk]) - 600.*g11[2*dk + ijk] + 
    150.*(-g11[-3*dk + ijk] + g11[3*dk + ijk]) - 25.*g11[4*dk + ijk] + 
    2.*g11[5*dk + ijk])
;

dg312
=
0.00079365079365079365079*oo2dz*
  (-2.*g12[-5*dk + ijk] + 25.*g12[-4*dk + ijk] + 600.*g12[-2*dk + ijk] + 
    2100.*(-g12[-dk + ijk] + g12[dk + ijk]) - 600.*g12[2*dk + ijk] + 
    150.*(-g12[-3*dk + ijk] + g12[3*dk + ijk]) - 25.*g12[4*dk + ijk] + 
    2.*g12[5*dk + ijk])
;

dg313
=
0.00079365079365079365079*oo2dz*
  (-2.*g13[-5*dk + ijk] + 25.*g13[-4*dk + ijk] + 600.*g13[-2*dk + ijk] + 
    2100.*(-g13[-dk + ijk] + g13[dk + ijk]) - 600.*g13[2*dk + ijk] + 
    150.*(-g13[-3*dk + ijk] + g13[3*dk + ijk]) - 25.*g13[4*dk + ijk] + 
    2.*g13[5*dk + ijk])
;

dg322
=
0.00079365079365079365079*oo2dz*
  (-2.*g22[-5*dk + ijk] + 25.*g22[-4*dk + ijk] + 600.*g22[-2*dk + ijk] + 
    2100.*(-g22[-dk + ijk] + g22[dk + ijk]) - 600.*g22[2*dk + ijk] + 
    150.*(-g22[-3*dk + ijk] + g22[3*dk + ijk]) - 25.*g22[4*dk + ijk] + 
    2.*g22[5*dk + ijk])
;

dg323
=
0.00079365079365079365079*oo2dz*
  (-2.*g23[-5*dk + ijk] + 25.*g23[-4*dk + ijk] + 600.*g23[-2*dk + ijk] + 
    2100.*(-g23[-dk + ijk] + g23[dk + ijk]) - 600.*g23[2*dk + ijk] + 
    150.*(-g23[-3*dk + ijk] + g23[3*dk + ijk]) - 25.*g23[4*dk + ijk] + 
    2.*g23[5*dk + ijk])
;

dg333
=
0.00079365079365079365079*oo2dz*
  (-2.*g33[-5*dk + ijk] + 25.*g33[-4*dk + ijk] + 600.*g33[-2*dk + ijk] + 
    2100.*(-g33[-dk + ijk] + g33[dk + ijk]) - 600.*g33[2*dk + ijk] + 
    150.*(-g33[-3*dk + ijk] + g33[3*dk + ijk]) - 25.*g33[4*dk + ijk] + 
    2.*g33[5*dk + ijk])
;

dK1
=
0.00079365079365079365079*oo2dx*
  (-2.*K[-5*di + ijk] + 25.*K[-4*di + ijk] + 600.*K[-2*di + ijk] + 
    2100.*(-K[-di + ijk] + K[di + ijk]) - 600.*K[2*di + ijk] + 
    150.*(-K[-3*di + ijk] + K[3*di + ijk]) - 25.*K[4*di + ijk] + 
    2.*K[5*di + ijk])
;

dK2
=
0.00079365079365079365079*oo2dy*
  (-2.*K[-5*dj + ijk] + 25.*K[-4*dj + ijk] + 600.*K[-2*dj + ijk] + 
    2100.*(-K[-dj + ijk] + K[dj + ijk]) - 600.*K[2*dj + ijk] + 
    150.*(-K[-3*dj + ijk] + K[3*dj + ijk]) - 25.*K[4*dj + ijk] + 
    2.*K[5*dj + ijk])
;

dK3
=
0.00079365079365079365079*oo2dz*
  (-2.*K[-5*dk + ijk] + 25.*K[-4*dk + ijk] + 600.*K[-2*dk + ijk] + 
    2100.*(-K[-dk + ijk] + K[dk + ijk]) - 600.*K[2*dk + ijk] + 
    150.*(-K[-3*dk + ijk] + K[3*dk + ijk]) - 25.*K[4*dk + ijk] + 
    2.*K[5*dk + ijk])
;

dA111
=
0.00079365079365079365079*oo2dx*
  (-2.*A11[-5*di + ijk] + 25.*A11[-4*di + ijk] + 600.*A11[-2*di + ijk] + 
    2100.*(-A11[-di + ijk] + A11[di + ijk]) - 600.*A11[2*di + ijk] + 
    150.*(-A11[-3*di + ijk] + A11[3*di + ijk]) - 25.*A11[4*di + ijk] + 
    2.*A11[5*di + ijk])
;

dA112
=
0.00079365079365079365079*oo2dx*
  (-2.*A12[-5*di + ijk] + 25.*A12[-4*di + ijk] + 600.*A12[-2*di + ijk] + 
    2100.*(-A12[-di + ijk] + A12[di + ijk]) - 600.*A12[2*di + ijk] + 
    150.*(-A12[-3*di + ijk] + A12[3*di + ijk]) - 25.*A12[4*di + ijk] + 
    2.*A12[5*di + ijk])
;

dA113
=
0.00079365079365079365079*oo2dx*
  (-2.*A13[-5*di + ijk] + 25.*A13[-4*di + ijk] + 600.*A13[-2*di + ijk] + 
    2100.*(-A13[-di + ijk] + A13[di + ijk]) - 600.*A13[2*di + ijk] + 
    150.*(-A13[-3*di + ijk] + A13[3*di + ijk]) - 25.*A13[4*di + ijk] + 
    2.*A13[5*di + ijk])
;

dA122
=
0.00079365079365079365079*oo2dx*
  (-2.*A22[-5*di + ijk] + 25.*A22[-4*di + ijk] + 600.*A22[-2*di + ijk] + 
    2100.*(-A22[-di + ijk] + A22[di + ijk]) - 600.*A22[2*di + ijk] + 
    150.*(-A22[-3*di + ijk] + A22[3*di + ijk]) - 25.*A22[4*di + ijk] + 
    2.*A22[5*di + ijk])
;

dA123
=
0.00079365079365079365079*oo2dx*
  (-2.*A23[-5*di + ijk] + 25.*A23[-4*di + ijk] + 600.*A23[-2*di + ijk] + 
    2100.*(-A23[-di + ijk] + A23[di + ijk]) - 600.*A23[2*di + ijk] + 
    150.*(-A23[-3*di + ijk] + A23[3*di + ijk]) - 25.*A23[4*di + ijk] + 
    2.*A23[5*di + ijk])
;

dA133
=
0.00079365079365079365079*oo2dx*
  (-2.*A33[-5*di + ijk] + 25.*A33[-4*di + ijk] + 600.*A33[-2*di + ijk] + 
    2100.*(-A33[-di + ijk] + A33[di + ijk]) - 600.*A33[2*di + ijk] + 
    150.*(-A33[-3*di + ijk] + A33[3*di + ijk]) - 25.*A33[4*di + ijk] + 
    2.*A33[5*di + ijk])
;

dA211
=
0.00079365079365079365079*oo2dy*
  (-2.*A11[-5*dj + ijk] + 25.*A11[-4*dj + ijk] + 600.*A11[-2*dj + ijk] + 
    2100.*(-A11[-dj + ijk] + A11[dj + ijk]) - 600.*A11[2*dj + ijk] + 
    150.*(-A11[-3*dj + ijk] + A11[3*dj + ijk]) - 25.*A11[4*dj + ijk] + 
    2.*A11[5*dj + ijk])
;

dA212
=
0.00079365079365079365079*oo2dy*
  (-2.*A12[-5*dj + ijk] + 25.*A12[-4*dj + ijk] + 600.*A12[-2*dj + ijk] + 
    2100.*(-A12[-dj + ijk] + A12[dj + ijk]) - 600.*A12[2*dj + ijk] + 
    150.*(-A12[-3*dj + ijk] + A12[3*dj + ijk]) - 25.*A12[4*dj + ijk] + 
    2.*A12[5*dj + ijk])
;

dA213
=
0.00079365079365079365079*oo2dy*
  (-2.*A13[-5*dj + ijk] + 25.*A13[-4*dj + ijk] + 600.*A13[-2*dj + ijk] + 
    2100.*(-A13[-dj + ijk] + A13[dj + ijk]) - 600.*A13[2*dj + ijk] + 
    150.*(-A13[-3*dj + ijk] + A13[3*dj + ijk]) - 25.*A13[4*dj + ijk] + 
    2.*A13[5*dj + ijk])
;

dA222
=
0.00079365079365079365079*oo2dy*
  (-2.*A22[-5*dj + ijk] + 25.*A22[-4*dj + ijk] + 600.*A22[-2*dj + ijk] + 
    2100.*(-A22[-dj + ijk] + A22[dj + ijk]) - 600.*A22[2*dj + ijk] + 
    150.*(-A22[-3*dj + ijk] + A22[3*dj + ijk]) - 25.*A22[4*dj + ijk] + 
    2.*A22[5*dj + ijk])
;

dA223
=
0.00079365079365079365079*oo2dy*
  (-2.*A23[-5*dj + ijk] + 25.*A23[-4*dj + ijk] + 600.*A23[-2*dj + ijk] + 
    2100.*(-A23[-dj + ijk] + A23[dj + ijk]) - 600.*A23[2*dj + ijk] + 
    150.*(-A23[-3*dj + ijk] + A23[3*dj + ijk]) - 25.*A23[4*dj + ijk] + 
    2.*A23[5*dj + ijk])
;

dA233
=
0.00079365079365079365079*oo2dy*
  (-2.*A33[-5*dj + ijk] + 25.*A33[-4*dj + ijk] + 600.*A33[-2*dj + ijk] + 
    2100.*(-A33[-dj + ijk] + A33[dj + ijk]) - 600.*A33[2*dj + ijk] + 
    150.*(-A33[-3*dj + ijk] + A33[3*dj + ijk]) - 25.*A33[4*dj + ijk] + 
    2.*A33[5*dj + ijk])
;

dA311
=
0.00079365079365079365079*oo2dz*
  (-2.*A11[-5*dk + ijk] + 25.*A11[-4*dk + ijk] + 600.*A11[-2*dk + ijk] + 
    2100.*(-A11[-dk + ijk] + A11[dk + ijk]) - 600.*A11[2*dk + ijk] + 
    150.*(-A11[-3*dk + ijk] + A11[3*dk + ijk]) - 25.*A11[4*dk + ijk] + 
    2.*A11[5*dk + ijk])
;

dA312
=
0.00079365079365079365079*oo2dz*
  (-2.*A12[-5*dk + ijk] + 25.*A12[-4*dk + ijk] + 600.*A12[-2*dk + ijk] + 
    2100.*(-A12[-dk + ijk] + A12[dk + ijk]) - 600.*A12[2*dk + ijk] + 
    150.*(-A12[-3*dk + ijk] + A12[3*dk + ijk]) - 25.*A12[4*dk + ijk] + 
    2.*A12[5*dk + ijk])
;

dA313
=
0.00079365079365079365079*oo2dz*
  (-2.*A13[-5*dk + ijk] + 25.*A13[-4*dk + ijk] + 600.*A13[-2*dk + ijk] + 
    2100.*(-A13[-dk + ijk] + A13[dk + ijk]) - 600.*A13[2*dk + ijk] + 
    150.*(-A13[-3*dk + ijk] + A13[3*dk + ijk]) - 25.*A13[4*dk + ijk] + 
    2.*A13[5*dk + ijk])
;

dA322
=
0.00079365079365079365079*oo2dz*
  (-2.*A22[-5*dk + ijk] + 25.*A22[-4*dk + ijk] + 600.*A22[-2*dk + ijk] + 
    2100.*(-A22[-dk + ijk] + A22[dk + ijk]) - 600.*A22[2*dk + ijk] + 
    150.*(-A22[-3*dk + ijk] + A22[3*dk + ijk]) - 25.*A22[4*dk + ijk] + 
    2.*A22[5*dk + ijk])
;

dA323
=
0.00079365079365079365079*oo2dz*
  (-2.*A23[-5*dk + ijk] + 25.*A23[-4*dk + ijk] + 600.*A23[-2*dk + ijk] + 
    2100.*(-A23[-dk + ijk] + A23[dk + ijk]) - 600.*A23[2*dk + ijk] + 
    150.*(-A23[-3*dk + ijk] + A23[3*dk + ijk]) - 25.*A23[4*dk + ijk] + 
    2.*A23[5*dk + ijk])
;

dA333
=
0.00079365079365079365079*oo2dz*
  (-2.*A33[-5*dk + ijk] + 25.*A33[-4*dk + ijk] + 600.*A33[-2*dk + ijk] + 
    2100.*(-A33[-dk + ijk] + A33[dk + ijk]) - 600.*A33[2*dk + ijk] + 
    150.*(-A33[-3*dk + ijk] + A33[3*dk + ijk]) - 25.*A33[4*dk + ijk] + 
    2.*A33[5*dk + ijk])
;

dG11
=
0.00079365079365079365079*oo2dx*
  (-2.*G1[-5*di + ijk] + 25.*G1[-4*di + ijk] + 600.*G1[-2*di + ijk] + 
    2100.*(-G1[-di + ijk] + G1[di + ijk]) - 600.*G1[2*di + ijk] + 
    150.*(-G1[-3*di + ijk] + G1[3*di + ijk]) - 25.*G1[4*di + ijk] + 
    2.*G1[5*di + ijk])
;

dG12
=
0.00079365079365079365079*oo2dx*
  (-2.*G2[-5*di + ijk] + 25.*G2[-4*di + ijk] + 600.*G2[-2*di + ijk] + 
    2100.*(-G2[-di + ijk] + G2[di + ijk]) - 600.*G2[2*di + ijk] + 
    150.*(-G2[-3*di + ijk] + G2[3*di + ijk]) - 25.*G2[4*di + ijk] + 
    2.*G2[5*di + ijk])
;

dG13
=
0.00079365079365079365079*oo2dx*
  (-2.*G3[-5*di + ijk] + 25.*G3[-4*di + ijk] + 600.*G3[-2*di + ijk] + 
    2100.*(-G3[-di + ijk] + G3[di + ijk]) - 600.*G3[2*di + ijk] + 
    150.*(-G3[-3*di + ijk] + G3[3*di + ijk]) - 25.*G3[4*di + ijk] + 
    2.*G3[5*di + ijk])
;

dG21
=
0.00079365079365079365079*oo2dy*
  (-2.*G1[-5*dj + ijk] + 25.*G1[-4*dj + ijk] + 600.*G1[-2*dj + ijk] + 
    2100.*(-G1[-dj + ijk] + G1[dj + ijk]) - 600.*G1[2*dj + ijk] + 
    150.*(-G1[-3*dj + ijk] + G1[3*dj + ijk]) - 25.*G1[4*dj + ijk] + 
    2.*G1[5*dj + ijk])
;

dG22
=
0.00079365079365079365079*oo2dy*
  (-2.*G2[-5*dj + ijk] + 25.*G2[-4*dj + ijk] + 600.*G2[-2*dj + ijk] + 
    2100.*(-G2[-dj + ijk] + G2[dj + ijk]) - 600.*G2[2*dj + ijk] + 
    150.*(-G2[-3*dj + ijk] + G2[3*dj + ijk]) - 25.*G2[4*dj + ijk] + 
    2.*G2[5*dj + ijk])
;

dG23
=
0.00079365079365079365079*oo2dy*
  (-2.*G3[-5*dj + ijk] + 25.*G3[-4*dj + ijk] + 600.*G3[-2*dj + ijk] + 
    2100.*(-G3[-dj + ijk] + G3[dj + ijk]) - 600.*G3[2*dj + ijk] + 
    150.*(-G3[-3*dj + ijk] + G3[3*dj + ijk]) - 25.*G3[4*dj + ijk] + 
    2.*G3[5*dj + ijk])
;

dG31
=
0.00079365079365079365079*oo2dz*
  (-2.*G1[-5*dk + ijk] + 25.*G1[-4*dk + ijk] + 600.*G1[-2*dk + ijk] + 
    2100.*(-G1[-dk + ijk] + G1[dk + ijk]) - 600.*G1[2*dk + ijk] + 
    150.*(-G1[-3*dk + ijk] + G1[3*dk + ijk]) - 25.*G1[4*dk + ijk] + 
    2.*G1[5*dk + ijk])
;

dG32
=
0.00079365079365079365079*oo2dz*
  (-2.*G2[-5*dk + ijk] + 25.*G2[-4*dk + ijk] + 600.*G2[-2*dk + ijk] + 
    2100.*(-G2[-dk + ijk] + G2[dk + ijk]) - 600.*G2[2*dk + ijk] + 
    150.*(-G2[-3*dk + ijk] + G2[3*dk + ijk]) - 25.*G2[4*dk + ijk] + 
    2.*G2[5*dk + ijk])
;

dG33
=
0.00079365079365079365079*oo2dz*
  (-2.*G3[-5*dk + ijk] + 25.*G3[-4*dk + ijk] + 600.*G3[-2*dk + ijk] + 
    2100.*(-G3[-dk + ijk] + G3[dk + ijk]) - 600.*G3[2*dk + ijk] + 
    150.*(-G3[-3*dk + ijk] + G3[3*dk + ijk]) - 25.*G3[4*dk + ijk] + 
    2.*G3[5*dk + ijk])
;


} else errorexit("order not implemented");  

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

dbSST11
=
db11*Jac11 + db21*Jac21 + db31*Jac31
;

dbSST12
=
db12*Jac11 + db22*Jac21 + db32*Jac31
;

dbSST13
=
db13*Jac11 + db23*Jac21 + db33*Jac31
;

dbSST21
=
db11*Jac12 + db21*Jac22 + db31*Jac32
;

dbSST22
=
db12*Jac12 + db22*Jac22 + db32*Jac32
;

dbSST23
=
db13*Jac12 + db23*Jac22 + db33*Jac32
;

dbSST31
=
db11*Jac13 + db21*Jac23 + db31*Jac33
;

dbSST32
=
db12*Jac13 + db22*Jac23 + db32*Jac33
;

dbSST33
=
db13*Jac13 + db23*Jac23 + db33*Jac33
;

dBSST11
=
dB11*Jac11 + dB21*Jac21 + dB31*Jac31
;

dBSST12
=
dB12*Jac11 + dB22*Jac21 + dB32*Jac31
;

dBSST13
=
dB13*Jac11 + dB23*Jac21 + dB33*Jac31
;

dBSST21
=
dB11*Jac12 + dB21*Jac22 + dB31*Jac32
;

dBSST22
=
dB12*Jac12 + dB22*Jac22 + dB32*Jac32
;

dBSST23
=
dB13*Jac12 + dB23*Jac22 + dB33*Jac32
;

dBSST31
=
dB11*Jac13 + dB21*Jac23 + dB31*Jac33
;

dBSST32
=
dB12*Jac13 + dB22*Jac23 + dB32*Jac33
;

dBSST33
=
dB13*Jac13 + dB23*Jac23 + dB33*Jac33
;

dgSST111
=
dg111*Jac11 + dg211*Jac21 + dg311*Jac31
;

dgSST112
=
dg112*Jac11 + dg212*Jac21 + dg312*Jac31
;

dgSST113
=
dg113*Jac11 + dg213*Jac21 + dg313*Jac31
;

dgSST122
=
dg122*Jac11 + dg222*Jac21 + dg322*Jac31
;

dgSST123
=
dg123*Jac11 + dg223*Jac21 + dg323*Jac31
;

dgSST133
=
dg133*Jac11 + dg233*Jac21 + dg333*Jac31
;

dgSST211
=
dg111*Jac12 + dg211*Jac22 + dg311*Jac32
;

dgSST212
=
dg112*Jac12 + dg212*Jac22 + dg312*Jac32
;

dgSST213
=
dg113*Jac12 + dg213*Jac22 + dg313*Jac32
;

dgSST222
=
dg122*Jac12 + dg222*Jac22 + dg322*Jac32
;

dgSST223
=
dg123*Jac12 + dg223*Jac22 + dg323*Jac32
;

dgSST233
=
dg133*Jac12 + dg233*Jac22 + dg333*Jac32
;

dgSST311
=
dg111*Jac13 + dg211*Jac23 + dg311*Jac33
;

dgSST312
=
dg112*Jac13 + dg212*Jac23 + dg312*Jac33
;

dgSST313
=
dg113*Jac13 + dg213*Jac23 + dg313*Jac33
;

dgSST322
=
dg122*Jac13 + dg222*Jac23 + dg322*Jac33
;

dgSST323
=
dg123*Jac13 + dg223*Jac23 + dg323*Jac33
;

dgSST333
=
dg133*Jac13 + dg233*Jac23 + dg333*Jac33
;

dGSST11
=
dG11*Jac11 + dG21*Jac21 + dG31*Jac31
;

dGSST12
=
dG12*Jac11 + dG22*Jac21 + dG32*Jac31
;

dGSST13
=
dG13*Jac11 + dG23*Jac21 + dG33*Jac31
;

dGSST21
=
dG11*Jac12 + dG21*Jac22 + dG31*Jac32
;

dGSST22
=
dG12*Jac12 + dG22*Jac22 + dG32*Jac32
;

dGSST23
=
dG13*Jac12 + dG23*Jac22 + dG33*Jac32
;

dGSST31
=
dG11*Jac13 + dG21*Jac23 + dG31*Jac33
;

dGSST32
=
dG12*Jac13 + dG22*Jac23 + dG32*Jac33
;

dGSST33
=
dG13*Jac13 + dG23*Jac23 + dG33*Jac33
;

dKSST1
=
dK1*Jac11 + dK2*Jac21 + dK3*Jac31
;

dKSST2
=
dK1*Jac12 + dK2*Jac22 + dK3*Jac32
;

dKSST3
=
dK1*Jac13 + dK2*Jac23 + dK3*Jac33
;

dASST111
=
dA111*Jac11 + dA211*Jac21 + dA311*Jac31
;

dASST112
=
dA112*Jac11 + dA212*Jac21 + dA312*Jac31
;

dASST113
=
dA113*Jac11 + dA213*Jac21 + dA313*Jac31
;

dASST122
=
dA122*Jac11 + dA222*Jac21 + dA322*Jac31
;

dASST123
=
dA123*Jac11 + dA223*Jac21 + dA323*Jac31
;

dASST133
=
dA133*Jac11 + dA233*Jac21 + dA333*Jac31
;

dASST211
=
dA111*Jac12 + dA211*Jac22 + dA311*Jac32
;

dASST212
=
dA112*Jac12 + dA212*Jac22 + dA312*Jac32
;

dASST213
=
dA113*Jac12 + dA213*Jac22 + dA313*Jac32
;

dASST222
=
dA122*Jac12 + dA222*Jac22 + dA322*Jac32
;

dASST223
=
dA123*Jac12 + dA223*Jac22 + dA323*Jac32
;

dASST233
=
dA133*Jac12 + dA233*Jac22 + dA333*Jac32
;

dASST311
=
dA111*Jac13 + dA211*Jac23 + dA311*Jac33
;

dASST312
=
dA112*Jac13 + dA212*Jac23 + dA312*Jac33
;

dASST313
=
dA113*Jac13 + dA213*Jac23 + dA313*Jac33
;

dASST322
=
dA122*Jac13 + dA222*Jac23 + dA322*Jac33
;

dASST323
=
dA123*Jac13 + dA223*Jac23 + dA323*Jac33
;

dASST333
=
dA133*Jac13 + dA233*Jac23 + dA333*Jac33
;

db11
=
dbSST11
;

db12
=
dbSST12
;

db13
=
dbSST13
;

db21
=
dbSST21
;

db22
=
dbSST22
;

db23
=
dbSST23
;

db31
=
dbSST31
;

db32
=
dbSST32
;

db33
=
dbSST33
;

dB11
=
dBSST11
;

dB12
=
dBSST12
;

dB13
=
dBSST13
;

dB21
=
dBSST21
;

dB22
=
dBSST22
;

dB23
=
dBSST23
;

dB31
=
dBSST31
;

dB32
=
dBSST32
;

dB33
=
dBSST33
;

dg111
=
dgSST111
;

dg112
=
dgSST112
;

dg113
=
dgSST113
;

dg122
=
dgSST122
;

dg123
=
dgSST123
;

dg133
=
dgSST133
;

dg211
=
dgSST211
;

dg212
=
dgSST212
;

dg213
=
dgSST213
;

dg222
=
dgSST222
;

dg223
=
dgSST223
;

dg233
=
dgSST233
;

dg311
=
dgSST311
;

dg312
=
dgSST312
;

dg313
=
dgSST313
;

dg322
=
dgSST322
;

dg323
=
dgSST323
;

dg333
=
dgSST333
;

dG11
=
dGSST11
;

dG12
=
dGSST12
;

dG13
=
dGSST13
;

dG21
=
dGSST21
;

dG22
=
dGSST22
;

dG23
=
dGSST23
;

dG31
=
dGSST31
;

dG32
=
dGSST32
;

dG33
=
dGSST33
;

dK1
=
dKSST1
;

dK2
=
dKSST2
;

dK3
=
dKSST3
;

dA111
=
dASST111
;

dA112
=
dASST112
;

dA113
=
dASST113
;

dA122
=
dASST122
;

dA123
=
dASST123
;

dA133
=
dASST133
;

dA211
=
dASST211
;

dA212
=
dASST212
;

dA213
=
dASST213
;

dA222
=
dASST222
;

dA223
=
dASST223
;

dA233
=
dASST233
;

dA311
=
dASST311
;

dA312
=
dASST312
;

dA313
=
dASST313
;

dA322
=
dASST322
;

dA323
=
dASST323
;

dA333
=
dASST333
;

r
=
0
;

shat1
=
0
;

shat2
=
0
;

shat3
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


r=rp[ijk]; 


shat1=xp[ijk]/r;shat2=yp[ijk]/r;shat3=zp[ijk]/r; 

detginv
=
1/(2.*g12[ijk]*g13[ijk]*g23[ijk] + g11[ijk]*g22[ijk]*g33[ijk] - 
    g33[ijk]*pow2(g12[ijk]) - g22[ijk]*pow2(g13[ijk]) - 
    g11[ijk]*pow2(g23[ijk]))
;

ginv11
=
detginv*(g22[ijk]*g33[ijk] - pow2(g23[ijk]))
;

ginv12
=
detginv*(g13[ijk]*g23[ijk] - g12[ijk]*g33[ijk])
;

ginv13
=
detginv*(-(g13[ijk]*g22[ijk]) + g12[ijk]*g23[ijk])
;

ginv22
=
detginv*(g11[ijk]*g33[ijk] - pow2(g13[ijk]))
;

ginv23
=
detginv*(g12[ijk]*g13[ijk] - g11[ijk]*g23[ijk])
;

ginv33
=
detginv*(g11[ijk]*g22[ijk] - pow2(g12[ijk]))
;

ADMginv11
=
ginv11*chi[ijk]
;

ADMginv12
=
ginv12*chi[ijk]
;

ADMginv13
=
ginv13*chi[ijk]
;

ADMginv22
=
ginv22*chi[ijk]
;

ADMginv23
=
ginv23*chi[ijk]
;

ADMginv33
=
ginv33*chi[ijk]
;

modshatARG
=
2.*(ADMginv23*shat2*shat3 + shat1*(ADMginv12*shat2 + ADMginv13*shat3)) + 
  ADMginv11*pow2(shat1) + ADMginv22*pow2(shat2) + ADMginv33*pow2(shat3)
;


if (modshatARG<0.00001) {                           
      printf("modshat is wrong (%e)\n",modshatARG);
      modshatARG = 0.00001;
    }oomodshat
=
1/sqrt(modshatARG)
;

sdown1
=
oomodshat*shat1
;

sdown2
=
oomodshat*shat2
;

sdown3
=
oomodshat*shat3
;

sup1
=
ADMginv11*sdown1 + ADMginv12*sdown2 + ADMginv13*sdown3
;

sup2
=
ADMginv12*sdown1 + ADMginv22*sdown2 + ADMginv23*sdown3
;

sup3
=
ADMginv13*sdown1 + ADMginv23*sdown2 + ADMginv33*sdown3
;

qud11
=
1. - sdown1*sup1
;

qud12
=
-(sdown2*sup1)
;

qud13
=
-(sdown3*sup1)
;

qud21
=
-(sdown1*sup2)
;

qud22
=
1. - sdown2*sup2
;

qud23
=
-(sdown3*sup2)
;

qud31
=
-(sdown1*sup3)
;

qud32
=
-(sdown2*sup3)
;

qud33
=
1. - sdown3*sup3
;

qdd11
=
g11[ijk]/chi[ijk] - pow2(sdown1)
;

qdd12
=
-(sdown1*sdown2) + g12[ijk]/chi[ijk]
;

qdd13
=
-(sdown1*sdown3) + g13[ijk]/chi[ijk]
;

qdd22
=
g22[ijk]/chi[ijk] - pow2(sdown2)
;

qdd23
=
-(sdown2*sdown3) + g23[ijk]/chi[ijk]
;

qdd33
=
g33[ijk]/chi[ijk] - pow2(sdown3)
;

quu11
=
ADMginv11 - pow2(sup1)
;

quu12
=
ADMginv12 - sup1*sup2
;

quu13
=
ADMginv13 - sup1*sup3
;

quu22
=
ADMginv22 - pow2(sup2)
;

quu23
=
ADMginv23 - sup2*sup3
;

quu33
=
ADMginv33 - pow2(sup3)
;

qPhysuudd1111
=
-0.5*qdd11*quu11 + pow2(qud11)
;

qPhysuudd1112
=
qud11*qud12 - 0.5*qdd12*quu11
;

qPhysuudd1113
=
qud11*qud13 - 0.5*qdd13*quu11
;

qPhysuudd1122
=
-0.5*qdd22*quu11 + pow2(qud12)
;

qPhysuudd1123
=
qud12*qud13 - 0.5*qdd23*quu11
;

qPhysuudd1133
=
-0.5*qdd33*quu11 + pow2(qud13)
;

qPhysuudd1211
=
qud11*qud21 - 0.5*qdd11*quu12
;

qPhysuudd1212
=
0.5*(qud12*qud21 + qud11*qud22 - qdd12*quu12)
;

qPhysuudd1213
=
0.5*(qud13*qud21 + qud11*qud23 - qdd13*quu12)
;

qPhysuudd1222
=
qud12*qud22 - 0.5*qdd22*quu12
;

qPhysuudd1223
=
0.5*(qud13*qud22 + qud12*qud23 - qdd23*quu12)
;

qPhysuudd1233
=
qud13*qud23 - 0.5*qdd33*quu12
;

qPhysuudd1311
=
qud11*qud31 - 0.5*qdd11*quu13
;

qPhysuudd1312
=
0.5*(qud12*qud31 + qud11*qud32 - qdd12*quu13)
;

qPhysuudd1313
=
0.5*(qud13*qud31 + qud11*qud33 - qdd13*quu13)
;

qPhysuudd1322
=
qud12*qud32 - 0.5*qdd22*quu13
;

qPhysuudd1323
=
0.5*(qud13*qud32 + qud12*qud33 - qdd23*quu13)
;

qPhysuudd1333
=
qud13*qud33 - 0.5*qdd33*quu13
;

qPhysuudd2211
=
-0.5*qdd11*quu22 + pow2(qud21)
;

qPhysuudd2212
=
qud21*qud22 - 0.5*qdd12*quu22
;

qPhysuudd2213
=
qud21*qud23 - 0.5*qdd13*quu22
;

qPhysuudd2222
=
-0.5*qdd22*quu22 + pow2(qud22)
;

qPhysuudd2223
=
qud22*qud23 - 0.5*qdd23*quu22
;

qPhysuudd2233
=
-0.5*qdd33*quu22 + pow2(qud23)
;

qPhysuudd2311
=
qud21*qud31 - 0.5*qdd11*quu23
;

qPhysuudd2312
=
0.5*(qud22*qud31 + qud21*qud32 - qdd12*quu23)
;

qPhysuudd2313
=
0.5*(qud23*qud31 + qud21*qud33 - qdd13*quu23)
;

qPhysuudd2322
=
qud22*qud32 - 0.5*qdd22*quu23
;

qPhysuudd2323
=
0.5*(qud23*qud32 + qud22*qud33 - qdd23*quu23)
;

qPhysuudd2333
=
qud23*qud33 - 0.5*qdd33*quu23
;

qPhysuudd3311
=
-0.5*qdd11*quu33 + pow2(qud31)
;

qPhysuudd3312
=
qud31*qud32 - 0.5*qdd12*quu33
;

qPhysuudd3313
=
qud31*qud33 - 0.5*qdd13*quu33
;

qPhysuudd3322
=
-0.5*qdd22*quu33 + pow2(qud32)
;

qPhysuudd3323
=
qud32*qud33 - 0.5*qdd23*quu33
;

qPhysuudd3333
=
-0.5*qdd33*quu33 + pow2(qud33)
;

muL
=
2./alpha[ijk]
;

muS
=
shiftgammacoeff
;

vbetas
=
2.*sqrt((0.33333333333333333333*muS)/chi[ijk])
;

vbetaA
=
sqrt(muS/chi[ijk])
;

Ainv11
=
2.*(ginv11*(ginv12*A12[ijk] + ginv13*A13[ijk]) + ginv12*ginv13*A23[ijk]) + 
  A11[ijk]*pow2(ginv11) + A22[ijk]*pow2(ginv12) + A33[ijk]*pow2(ginv13)
;

Ainv12
=
ginv11*(ginv12*A11[ijk] + ginv22*A12[ijk] + ginv23*A13[ijk]) + 
  ginv12*(ginv13*A13[ijk] + ginv22*A22[ijk] + ginv23*A23[ijk]) + 
  ginv13*(ginv22*A23[ijk] + ginv23*A33[ijk]) + A12[ijk]*pow2(ginv12)
;

Ainv13
=
ginv11*(ginv13*A11[ijk] + ginv23*A12[ijk] + ginv33*A13[ijk]) + 
  ginv12*(ginv13*A12[ijk] + ginv23*A22[ijk] + ginv33*A23[ijk]) + 
  ginv13*(ginv23*A23[ijk] + ginv33*A33[ijk]) + A13[ijk]*pow2(ginv13)
;

Ainv22
=
2.*(ginv12*(ginv22*A12[ijk] + ginv23*A13[ijk]) + ginv22*ginv23*A23[ijk]) + 
  A11[ijk]*pow2(ginv12) + A22[ijk]*pow2(ginv22) + A33[ijk]*pow2(ginv23)
;

Ainv23
=
ginv13*(ginv22*A12[ijk] + ginv23*A13[ijk]) + 
  ginv12*(ginv13*A11[ijk] + ginv23*A12[ijk] + ginv33*A13[ijk]) + 
  ginv22*(ginv23*A22[ijk] + ginv33*A23[ijk]) + ginv23*ginv33*A33[ijk] + 
  A23[ijk]*pow2(ginv23)
;

Ainv33
=
2.*(ginv13*(ginv23*A12[ijk] + ginv33*A13[ijk]) + ginv23*ginv33*A23[ijk]) + 
  A11[ijk]*pow2(ginv13) + A22[ijk]*pow2(ginv23) + A33[ijk]*pow2(ginv33)
;

divbeta
=
db11 + db22 + db33
;

totdivbeta
=
0.66666666666666666667*divbeta
;

lieg11
=
dg111*beta1[ijk] + dg211*beta2[ijk] + dg311*beta3[ijk] + 
  (2.*db11 - totdivbeta)*g11[ijk] + 2.*(db12*g12[ijk] + db13*g13[ijk])
;

lieg12
=
dg112*beta1[ijk] + dg212*beta2[ijk] + dg312*beta3[ijk] + db21*g11[ijk] + 
  (db11 + db22 - totdivbeta)*g12[ijk] + db23*g13[ijk] + db12*g22[ijk] + 
  db13*g23[ijk]
;

lieg13
=
dg113*beta1[ijk] + dg213*beta2[ijk] + dg313*beta3[ijk] + db31*g11[ijk] + 
  db32*g12[ijk] + (db11 + db33 - totdivbeta)*g13[ijk] + db12*g23[ijk] + 
  db13*g33[ijk]
;

lieg22
=
dg122*beta1[ijk] + dg222*beta2[ijk] + dg322*beta3[ijk] - 
  totdivbeta*g22[ijk] + 2.*(db21*g12[ijk] + db22*g22[ijk] + db23*g23[ijk])
;

lieg23
=
dg123*beta1[ijk] + dg223*beta2[ijk] + dg323*beta3[ijk] + db31*g12[ijk] + 
  db21*g13[ijk] + db32*g22[ijk] + (db22 + db33 - totdivbeta)*g23[ijk] + 
  db23*g33[ijk]
;

lieg33
=
dg133*beta1[ijk] + dg233*beta2[ijk] + dg333*beta3[ijk] - 
  totdivbeta*g33[ijk] + 2.*(db31*g13[ijk] + db32*g23[ijk] + db33*g33[ijk])
;

ACss
=
2.*(sup1*(sup2*A12[ijk] + sup3*A13[ijk]) + sup2*sup3*A23[ijk]) + 
  A11[ijk]*pow2(sup1) + A22[ijk]*pow2(sup2) + A33[ijk]*pow2(sup3)
;

DACss
=
((2.*dA112 + dA211)*sup2 + 2.*dA113*sup3)*pow2(sup1) + 
  sup3*(dA311*pow2(sup1) + dA322*pow2(sup2) + 
     2.*((dA213 + dA312)*sup1*sup2 + dA223*pow2(sup2))) + 
  (dA233 + 2.*dA323)*sup2*pow2(sup3) + 
  sup1*((dA122 + 2.*dA212)*pow2(sup2) + dA133*pow2(sup3) + 
     2.*(dA123*sup2*sup3 + dA313*pow2(sup3))) + dA111*pow3(sup1) + 
  dA222*pow3(sup2) + dA333*pow3(sup3)
;

DK
=
dK1*sup1 + dK2*sup2 + dK3*sup3
;

Bs
=
sdown1*B1[ijk] + sdown2*B2[ijk] + sdown3*B3[ijk]
;

DBs
=
(dB11*sdown1 + dB12*sdown2 + dB13*sdown3)*sup1 + 
  (dB21*sdown1 + dB22*sdown2 + dB23*sdown3)*sup2 + 
  (dB31*sdown1 + dB32*sdown2 + dB33*sdown3)*sup3
;

Gams
=
sdown1*G1[ijk] + sdown2*G2[ijk] + sdown3*G3[ijk]
;

DGams
=
(dG11*sdown1 + dG12*sdown2 + dG13*sdown3)*sup1 + 
  (dG21*sdown1 + dG22*sdown2 + dG23*sdown3)*sup2 + 
  (dG31*sdown1 + dG32*sdown2 + dG33*sdown3)*sup3
;

BA1
=
qud11*B1[ijk] + qud12*B2[ijk] + qud13*B3[ijk]
;

BA2
=
qud21*B1[ijk] + qud22*B2[ijk] + qud23*B3[ijk]
;

BA3
=
qud31*B1[ijk] + qud32*B2[ijk] + qud33*B3[ijk]
;

DBA1
=
(dB11*qud11 + dB12*qud12 + dB13*qud13)*sup1 + 
  (dB21*qud11 + dB22*qud12 + dB23*qud13)*sup2 + 
  (dB31*qud11 + dB32*qud12 + dB33*qud13)*sup3
;

DBA2
=
(dB11*qud21 + dB12*qud22 + dB13*qud23)*sup1 + 
  (dB21*qud21 + dB22*qud22 + dB23*qud23)*sup2 + 
  (dB31*qud21 + dB32*qud22 + dB33*qud23)*sup3
;

DBA3
=
(dB11*qud31 + dB12*qud32 + dB13*qud33)*sup1 + 
  (dB21*qud31 + dB22*qud32 + dB23*qud33)*sup2 + 
  (dB31*qud31 + dB32*qud32 + dB33*qud33)*sup3
;

GamA1
=
qud11*G1[ijk] + qud12*G2[ijk] + qud13*G3[ijk]
;

GamA2
=
qud21*G1[ijk] + qud22*G2[ijk] + qud23*G3[ijk]
;

GamA3
=
qud31*G1[ijk] + qud32*G2[ijk] + qud33*G3[ijk]
;

DGamA1
=
(dG11*qud11 + dG12*qud12 + dG13*qud13)*sup1 + 
  (dG21*qud11 + dG22*qud12 + dG23*qud13)*sup2 + 
  (dG31*qud11 + dG32*qud12 + dG33*qud13)*sup3
;

DGamA2
=
(dG11*qud21 + dG12*qud22 + dG13*qud23)*sup1 + 
  (dG21*qud21 + dG22*qud22 + dG23*qud23)*sup2 + 
  (dG31*qud21 + dG32*qud22 + dG33*qud23)*sup3
;

DGamA3
=
(dG11*qud31 + dG12*qud32 + dG13*qud33)*sup1 + 
  (dG21*qud31 + dG22*qud32 + dG23*qud33)*sup2 + 
  (dG31*qud31 + dG32*qud32 + dG33*qud33)*sup3
;

ACsA1
=
sup1*(qud11*A11[ijk] + qud21*A12[ijk] + qud31*A13[ijk]) + 
  qud11*(sup2*A12[ijk] + sup3*A13[ijk]) + 
  sup2*(qud21*A22[ijk] + qud31*A23[ijk]) + 
  sup3*(qud21*A23[ijk] + qud31*A33[ijk])
;

ACsA2
=
sup1*(qud12*A11[ijk] + qud22*A12[ijk] + qud32*A13[ijk]) + 
  qud12*(sup2*A12[ijk] + sup3*A13[ijk]) + 
  sup2*(qud22*A22[ijk] + qud32*A23[ijk]) + 
  sup3*(qud22*A23[ijk] + qud32*A33[ijk])
;

ACsA3
=
sup1*(qud13*A11[ijk] + qud23*A12[ijk] + qud33*A13[ijk]) + 
  qud13*(sup2*A12[ijk] + sup3*A13[ijk]) + 
  sup2*(qud23*A22[ijk] + qud33*A23[ijk]) + 
  sup3*(qud23*A23[ijk] + qud33*A33[ijk])
;

DACsA1
=
((dA311*qud11 + dA312*qud21)*sup1 + 
     ((dA213 + dA312)*qud11 + (dA223 + dA322)*qud21 + dA233*qud31)*sup2 + 
     qud31*(dA313*sup1 + dA323*sup2))*sup3 + 
  sup1*(qud11*((dA112 + dA211)*sup2 + dA113*sup3) + 
     qud21*((dA122 + dA212)*sup2 + dA123*sup3) + 
     qud31*((dA123 + dA213)*sup2 + dA133*sup3)) + 
  (dA111*qud11 + dA112*qud21 + dA113*qud31)*pow2(sup1) + 
  (dA212*qud11 + dA222*qud21 + dA223*qud31)*pow2(sup2) + 
  (dA313*qud11 + dA323*qud21 + dA333*qud31)*pow2(sup3)
;

DACsA2
=
((dA311*qud12 + dA312*qud22)*sup1 + 
     ((dA213 + dA312)*qud12 + (dA223 + dA322)*qud22 + dA233*qud32)*sup2 + 
     qud32*(dA313*sup1 + dA323*sup2))*sup3 + 
  sup1*(qud12*((dA112 + dA211)*sup2 + dA113*sup3) + 
     qud22*((dA122 + dA212)*sup2 + dA123*sup3) + 
     qud32*((dA123 + dA213)*sup2 + dA133*sup3)) + 
  (dA111*qud12 + dA112*qud22 + dA113*qud32)*pow2(sup1) + 
  (dA212*qud12 + dA222*qud22 + dA223*qud32)*pow2(sup2) + 
  (dA313*qud12 + dA323*qud22 + dA333*qud32)*pow2(sup3)
;

DACsA3
=
((dA311*qud13 + dA312*qud23)*sup1 + 
     ((dA213 + dA312)*qud13 + (dA223 + dA322)*qud23 + dA233*qud33)*sup2 + 
     qud33*(dA313*sup1 + dA323*sup2))*sup3 + 
  sup1*(qud13*((dA112 + dA211)*sup2 + dA113*sup3) + 
     qud23*((dA122 + dA212)*sup2 + dA123*sup3) + 
     qud33*((dA123 + dA213)*sup2 + dA133*sup3)) + 
  (dA111*qud13 + dA112*qud23 + dA113*qud33)*pow2(sup1) + 
  (dA212*qud13 + dA222*qud23 + dA223*qud33)*pow2(sup2) + 
  (dA313*qud13 + dA323*qud23 + dA333*qud33)*pow2(sup3)
;

ACABTF11
=
qPhysuudd1111*A11[ijk] + qPhysuudd2211*A22[ijk] + 
  2.*(qPhysuudd1211*A12[ijk] + qPhysuudd1311*A13[ijk] + 
     qPhysuudd2311*A23[ijk]) + qPhysuudd3311*A33[ijk]
;

ACABTF12
=
qPhysuudd1112*A11[ijk] + qPhysuudd2212*A22[ijk] + 
  2.*(qPhysuudd1212*A12[ijk] + qPhysuudd1312*A13[ijk] + 
     qPhysuudd2312*A23[ijk]) + qPhysuudd3312*A33[ijk]
;

ACABTF13
=
qPhysuudd1113*A11[ijk] + qPhysuudd2213*A22[ijk] + 
  2.*(qPhysuudd1213*A12[ijk] + qPhysuudd1313*A13[ijk] + 
     qPhysuudd2313*A23[ijk]) + qPhysuudd3313*A33[ijk]
;

ACABTF21
=
qPhysuudd1112*A11[ijk] + qPhysuudd2212*A22[ijk] + 
  2.*(qPhysuudd1212*A12[ijk] + qPhysuudd1312*A13[ijk] + 
     qPhysuudd2312*A23[ijk]) + qPhysuudd3312*A33[ijk]
;

ACABTF22
=
qPhysuudd1122*A11[ijk] + qPhysuudd2222*A22[ijk] + 
  2.*(qPhysuudd1222*A12[ijk] + qPhysuudd1322*A13[ijk] + 
     qPhysuudd2322*A23[ijk]) + qPhysuudd3322*A33[ijk]
;

ACABTF23
=
qPhysuudd1123*A11[ijk] + qPhysuudd2223*A22[ijk] + 
  2.*(qPhysuudd1223*A12[ijk] + qPhysuudd1323*A13[ijk] + 
     qPhysuudd2323*A23[ijk]) + qPhysuudd3323*A33[ijk]
;

ACABTF31
=
qPhysuudd1113*A11[ijk] + qPhysuudd2213*A22[ijk] + 
  2.*(qPhysuudd1213*A12[ijk] + qPhysuudd1313*A13[ijk] + 
     qPhysuudd2313*A23[ijk]) + qPhysuudd3313*A33[ijk]
;

ACABTF32
=
qPhysuudd1123*A11[ijk] + qPhysuudd2223*A22[ijk] + 
  2.*(qPhysuudd1223*A12[ijk] + qPhysuudd1323*A13[ijk] + 
     qPhysuudd2323*A23[ijk]) + qPhysuudd3323*A33[ijk]
;

ACABTF33
=
qPhysuudd1133*A11[ijk] + qPhysuudd2233*A22[ijk] + 
  2.*(qPhysuudd1233*A12[ijk] + qPhysuudd1333*A13[ijk] + 
     qPhysuudd2333*A23[ijk]) + qPhysuudd3333*A33[ijk]
;

DACABTF11
=
(dA111*qPhysuudd1111 + dA122*qPhysuudd2211 + 
     2.*(dA112*qPhysuudd1211 + dA113*qPhysuudd1311 + 
        dA123*qPhysuudd2311) + dA133*qPhysuudd3311)*sup1 + 
  (dA211*qPhysuudd1111 + dA222*qPhysuudd2211 + 
     2.*(dA212*qPhysuudd1211 + dA213*qPhysuudd1311 + 
        dA223*qPhysuudd2311) + dA233*qPhysuudd3311)*sup2 + 
  (dA311*qPhysuudd1111 + dA322*qPhysuudd2211 + 
     2.*(dA312*qPhysuudd1211 + dA313*qPhysuudd1311 + dA323*qPhysuudd2311) + 
     dA333*qPhysuudd3311)*sup3
;

DACABTF12
=
(dA111*qPhysuudd1112 + dA122*qPhysuudd2212 + 
     2.*(dA112*qPhysuudd1212 + dA113*qPhysuudd1312 + 
        dA123*qPhysuudd2312) + dA133*qPhysuudd3312)*sup1 + 
  (dA211*qPhysuudd1112 + dA222*qPhysuudd2212 + 
     2.*(dA212*qPhysuudd1212 + dA213*qPhysuudd1312 + 
        dA223*qPhysuudd2312) + dA233*qPhysuudd3312)*sup2 + 
  (dA311*qPhysuudd1112 + dA322*qPhysuudd2212 + 
     2.*(dA312*qPhysuudd1212 + dA313*qPhysuudd1312 + dA323*qPhysuudd2312) + 
     dA333*qPhysuudd3312)*sup3
;

DACABTF13
=
(dA111*qPhysuudd1113 + dA122*qPhysuudd2213 + 
     2.*(dA112*qPhysuudd1213 + dA113*qPhysuudd1313 + 
        dA123*qPhysuudd2313) + dA133*qPhysuudd3313)*sup1 + 
  (dA211*qPhysuudd1113 + dA222*qPhysuudd2213 + 
     2.*(dA212*qPhysuudd1213 + dA213*qPhysuudd1313 + 
        dA223*qPhysuudd2313) + dA233*qPhysuudd3313)*sup2 + 
  (dA311*qPhysuudd1113 + dA322*qPhysuudd2213 + 
     2.*(dA312*qPhysuudd1213 + dA313*qPhysuudd1313 + dA323*qPhysuudd2313) + 
     dA333*qPhysuudd3313)*sup3
;

DACABTF21
=
(dA111*qPhysuudd1112 + dA122*qPhysuudd2212 + 
     2.*(dA112*qPhysuudd1212 + dA113*qPhysuudd1312 + 
        dA123*qPhysuudd2312) + dA133*qPhysuudd3312)*sup1 + 
  (dA211*qPhysuudd1112 + dA222*qPhysuudd2212 + 
     2.*(dA212*qPhysuudd1212 + dA213*qPhysuudd1312 + 
        dA223*qPhysuudd2312) + dA233*qPhysuudd3312)*sup2 + 
  (dA311*qPhysuudd1112 + dA322*qPhysuudd2212 + 
     2.*(dA312*qPhysuudd1212 + dA313*qPhysuudd1312 + dA323*qPhysuudd2312) + 
     dA333*qPhysuudd3312)*sup3
;

DACABTF22
=
(dA111*qPhysuudd1122 + dA122*qPhysuudd2222 + 
     2.*(dA112*qPhysuudd1222 + dA113*qPhysuudd1322 + 
        dA123*qPhysuudd2322) + dA133*qPhysuudd3322)*sup1 + 
  (dA211*qPhysuudd1122 + dA222*qPhysuudd2222 + 
     2.*(dA212*qPhysuudd1222 + dA213*qPhysuudd1322 + 
        dA223*qPhysuudd2322) + dA233*qPhysuudd3322)*sup2 + 
  (dA311*qPhysuudd1122 + dA322*qPhysuudd2222 + 
     2.*(dA312*qPhysuudd1222 + dA313*qPhysuudd1322 + dA323*qPhysuudd2322) + 
     dA333*qPhysuudd3322)*sup3
;

DACABTF23
=
(dA111*qPhysuudd1123 + dA122*qPhysuudd2223 + 
     2.*(dA112*qPhysuudd1223 + dA113*qPhysuudd1323 + 
        dA123*qPhysuudd2323) + dA133*qPhysuudd3323)*sup1 + 
  (dA211*qPhysuudd1123 + dA222*qPhysuudd2223 + 
     2.*(dA212*qPhysuudd1223 + dA213*qPhysuudd1323 + 
        dA223*qPhysuudd2323) + dA233*qPhysuudd3323)*sup2 + 
  (dA311*qPhysuudd1123 + dA322*qPhysuudd2223 + 
     2.*(dA312*qPhysuudd1223 + dA313*qPhysuudd1323 + dA323*qPhysuudd2323) + 
     dA333*qPhysuudd3323)*sup3
;

DACABTF31
=
(dA111*qPhysuudd1113 + dA122*qPhysuudd2213 + 
     2.*(dA112*qPhysuudd1213 + dA113*qPhysuudd1313 + 
        dA123*qPhysuudd2313) + dA133*qPhysuudd3313)*sup1 + 
  (dA211*qPhysuudd1113 + dA222*qPhysuudd2213 + 
     2.*(dA212*qPhysuudd1213 + dA213*qPhysuudd1313 + 
        dA223*qPhysuudd2313) + dA233*qPhysuudd3313)*sup2 + 
  (dA311*qPhysuudd1113 + dA322*qPhysuudd2213 + 
     2.*(dA312*qPhysuudd1213 + dA313*qPhysuudd1313 + dA323*qPhysuudd2313) + 
     dA333*qPhysuudd3313)*sup3
;

DACABTF32
=
(dA111*qPhysuudd1123 + dA122*qPhysuudd2223 + 
     2.*(dA112*qPhysuudd1223 + dA113*qPhysuudd1323 + 
        dA123*qPhysuudd2323) + dA133*qPhysuudd3323)*sup1 + 
  (dA211*qPhysuudd1123 + dA222*qPhysuudd2223 + 
     2.*(dA212*qPhysuudd1223 + dA213*qPhysuudd1323 + 
        dA223*qPhysuudd2323) + dA233*qPhysuudd3323)*sup2 + 
  (dA311*qPhysuudd1123 + dA322*qPhysuudd2223 + 
     2.*(dA312*qPhysuudd1223 + dA313*qPhysuudd1323 + dA323*qPhysuudd2323) + 
     dA333*qPhysuudd3323)*sup3
;

DACABTF33
=
(dA111*qPhysuudd1133 + dA122*qPhysuudd2233 + 
     2.*(dA112*qPhysuudd1233 + dA113*qPhysuudd1333 + 
        dA123*qPhysuudd2333) + dA133*qPhysuudd3333)*sup1 + 
  (dA211*qPhysuudd1133 + dA222*qPhysuudd2233 + 
     2.*(dA212*qPhysuudd1233 + dA213*qPhysuudd1333 + 
        dA223*qPhysuudd2333) + dA233*qPhysuudd3333)*sup2 + 
  (dA311*qPhysuudd1133 + dA322*qPhysuudd2233 + 
     2.*(dA312*qPhysuudd1233 + dA313*qPhysuudd1333 + dA323*qPhysuudd2333) + 
     dA333*qPhysuudd3333)*sup3
;

rK
=
dK1*beta1[ijk] + dK2*beta2[ijk] + dK3*beta3[ijk] - 
  alpha[ijk]*(DK + K[ijk]/r)*sqrt(muL)
;

rGams
=
-((DGams + Gams/r)*vbetas) + (dG11*sdown1 + dG12*sdown2 + dG13*sdown3)*
   beta1[ijk] + (dG21*sdown1 + dG22*sdown2 + dG23*sdown3)*beta2[ijk] + 
  (dG31*sdown1 + dG32*sdown2 + dG33*sdown3)*beta3[ijk]
;

rBs
=
-((DBs + Bs/r)*vbetas) + (dB11*sdown1 + dB12*sdown2 + dB13*sdown3)*
   beta1[ijk] + (dB21*sdown1 + dB22*sdown2 + dB23*sdown3)*beta2[ijk] + 
  (dB31*sdown1 + dB32*sdown2 + dB33*sdown3)*beta3[ijk]
;

rACss
=
-DACss - ACss/r + beta1[ijk]*(2.*
      (dA123*sup2*sup3 + sup1*(dA112*sup2 + dA113*sup3)) + 
     dA111*pow2(sup1) + dA122*pow2(sup2) + dA133*pow2(sup3)) + 
  beta2[ijk]*(2.*(dA223*sup2*sup3 + sup1*(dA212*sup2 + dA213*sup3)) + 
     dA211*pow2(sup1) + dA222*pow2(sup2) + dA233*pow2(sup3)) + 
  beta3[ijk]*(2.*(dA323*sup2*sup3 + sup1*(dA312*sup2 + dA313*sup3)) + 
     dA311*pow2(sup1) + dA322*pow2(sup2) + dA333*pow2(sup3))
;

rACqq
=
-rACss + (Ainv22*lieg22 + 2.*(Ainv12*lieg12 + Ainv13*lieg13 + 
        Ainv23*lieg23) - (2.*Ainv22*A22[ijk] + 
        4.*(Ainv12*A12[ijk] + Ainv13*A13[ijk] + Ainv23*A23[ijk]))*
      alpha[ijk] + Ainv11*(lieg11 - 2.*A11[ijk]*alpha[ijk]) + 
     Ainv33*(lieg33 - 2.*A33[ijk]*alpha[ijk]))*chi[ijk]
;

rGamA1
=
-((DGamA1 + GamA1/r)*vbetaA) + 
  (dG11*qud11 + dG12*qud12 + dG13*qud13)*beta1[ijk] + 
  (dG21*qud11 + dG22*qud12 + dG23*qud13)*beta2[ijk] + 
  (dG31*qud11 + dG32*qud12 + dG33*qud13)*beta3[ijk]
;

rGamA2
=
-((DGamA2 + GamA2/r)*vbetaA) + 
  (dG11*qud21 + dG12*qud22 + dG13*qud23)*beta1[ijk] + 
  (dG21*qud21 + dG22*qud22 + dG23*qud23)*beta2[ijk] + 
  (dG31*qud21 + dG32*qud22 + dG33*qud23)*beta3[ijk]
;

rGamA3
=
-((DGamA3 + GamA3/r)*vbetaA) + 
  (dG11*qud31 + dG12*qud32 + dG13*qud33)*beta1[ijk] + 
  (dG21*qud31 + dG22*qud32 + dG23*qud33)*beta2[ijk] + 
  (dG31*qud31 + dG32*qud32 + dG33*qud33)*beta3[ijk]
;

rBA1
=
-((DBA1 + BA1/r)*vbetaA) + (dB11*qud11 + dB12*qud12 + dB13*qud13)*
   beta1[ijk] + (dB21*qud11 + dB22*qud12 + dB23*qud13)*beta2[ijk] + 
  (dB31*qud11 + dB32*qud12 + dB33*qud13)*beta3[ijk]
;

rBA2
=
-((DBA2 + BA2/r)*vbetaA) + (dB11*qud21 + dB12*qud22 + dB13*qud23)*
   beta1[ijk] + (dB21*qud21 + dB22*qud22 + dB23*qud23)*beta2[ijk] + 
  (dB31*qud21 + dB32*qud22 + dB33*qud23)*beta3[ijk]
;

rBA3
=
-((DBA3 + BA3/r)*vbetaA) + (dB11*qud31 + dB12*qud32 + dB13*qud33)*
   beta1[ijk] + (dB21*qud31 + dB22*qud32 + dB23*qud33)*beta2[ijk] + 
  (dB31*qud31 + dB32*qud32 + dB33*qud33)*beta3[ijk]
;

rACsA1
=
-((DACsA1 + ACsA1/r)*alpha[ijk]) + 
  (dA112*qud11*sup2 + dA113*(qud31*sup1 + qud11*sup3))*beta1[ijk] + 
  sup1*((dA111*qud11 + dA112*qud21)*beta1[ijk] + 
     (dA211*qud11 + dA212*qud21 + dA213*qud31)*beta2[ijk] + 
     (dA311*qud11 + dA312*qud21 + dA313*qud31)*beta3[ijk]) + 
  sup2*((dA122*qud21 + dA123*qud31)*beta1[ijk] + 
     (dA212*qud11 + dA222*qud21 + dA223*qud31)*beta2[ijk] + 
     (dA312*qud11 + dA322*qud21 + dA323*qud31)*beta3[ijk]) + 
  sup3*((dA123*qud21 + dA133*qud31)*beta1[ijk] + 
     (dA213*qud11 + dA223*qud21 + dA233*qud31)*beta2[ijk] + 
     (dA313*qud11 + dA323*qud21 + dA333*qud31)*beta3[ijk])
;

rACsA2
=
-((DACsA2 + ACsA2/r)*alpha[ijk]) + 
  (dA112*qud12*sup2 + dA113*(qud32*sup1 + qud12*sup3))*beta1[ijk] + 
  sup1*((dA111*qud12 + dA112*qud22)*beta1[ijk] + 
     (dA211*qud12 + dA212*qud22 + dA213*qud32)*beta2[ijk] + 
     (dA311*qud12 + dA312*qud22 + dA313*qud32)*beta3[ijk]) + 
  sup2*((dA122*qud22 + dA123*qud32)*beta1[ijk] + 
     (dA212*qud12 + dA222*qud22 + dA223*qud32)*beta2[ijk] + 
     (dA312*qud12 + dA322*qud22 + dA323*qud32)*beta3[ijk]) + 
  sup3*((dA123*qud22 + dA133*qud32)*beta1[ijk] + 
     (dA213*qud12 + dA223*qud22 + dA233*qud32)*beta2[ijk] + 
     (dA313*qud12 + dA323*qud22 + dA333*qud32)*beta3[ijk])
;

rACsA3
=
-((DACsA3 + ACsA3/r)*alpha[ijk]) + 
  (dA112*qud13*sup2 + dA113*(qud33*sup1 + qud13*sup3))*beta1[ijk] + 
  sup1*((dA111*qud13 + dA112*qud23)*beta1[ijk] + 
     (dA211*qud13 + dA212*qud23 + dA213*qud33)*beta2[ijk] + 
     (dA311*qud13 + dA312*qud23 + dA313*qud33)*beta3[ijk]) + 
  sup2*((dA122*qud23 + dA123*qud33)*beta1[ijk] + 
     (dA212*qud13 + dA222*qud23 + dA223*qud33)*beta2[ijk] + 
     (dA312*qud13 + dA322*qud23 + dA323*qud33)*beta3[ijk]) + 
  sup3*((dA123*qud23 + dA133*qud33)*beta1[ijk] + 
     (dA213*qud13 + dA223*qud23 + dA233*qud33)*beta2[ijk] + 
     (dA313*qud13 + dA323*qud23 + dA333*qud33)*beta3[ijk])
;

rACABTF11
=
-((DACABTF11 + ACABTF11/r)*alpha[ijk]) + 
  (dA111*qPhysuudd1111 + dA122*qPhysuudd2211 + 
     2.*(dA112*qPhysuudd1211 + dA113*qPhysuudd1311 + 
        dA123*qPhysuudd2311) + dA133*qPhysuudd3311)*beta1[ijk] + 
  (dA211*qPhysuudd1111 + dA222*qPhysuudd2211 + 
     2.*(dA212*qPhysuudd1211 + dA213*qPhysuudd1311 + 
        dA223*qPhysuudd2311) + dA233*qPhysuudd3311)*beta2[ijk] + 
  (dA311*qPhysuudd1111 + dA322*qPhysuudd2211 + 
     2.*(dA312*qPhysuudd1211 + dA313*qPhysuudd1311 + dA323*qPhysuudd2311) + 
     dA333*qPhysuudd3311)*beta3[ijk]
;

rACABTF12
=
-((DACABTF21 + ACABTF21/r)*alpha[ijk]) + 
  (dA111*qPhysuudd1112 + dA122*qPhysuudd2212 + 
     2.*(dA112*qPhysuudd1212 + dA113*qPhysuudd1312 + 
        dA123*qPhysuudd2312) + dA133*qPhysuudd3312)*beta1[ijk] + 
  (dA211*qPhysuudd1112 + dA222*qPhysuudd2212 + 
     2.*(dA212*qPhysuudd1212 + dA213*qPhysuudd1312 + 
        dA223*qPhysuudd2312) + dA233*qPhysuudd3312)*beta2[ijk] + 
  (dA311*qPhysuudd1112 + dA322*qPhysuudd2212 + 
     2.*(dA312*qPhysuudd1212 + dA313*qPhysuudd1312 + dA323*qPhysuudd2312) + 
     dA333*qPhysuudd3312)*beta3[ijk]
;

rACABTF13
=
-((DACABTF31 + ACABTF31/r)*alpha[ijk]) + 
  (dA111*qPhysuudd1113 + dA122*qPhysuudd2213 + 
     2.*(dA112*qPhysuudd1213 + dA113*qPhysuudd1313 + 
        dA123*qPhysuudd2313) + dA133*qPhysuudd3313)*beta1[ijk] + 
  (dA211*qPhysuudd1113 + dA222*qPhysuudd2213 + 
     2.*(dA212*qPhysuudd1213 + dA213*qPhysuudd1313 + 
        dA223*qPhysuudd2313) + dA233*qPhysuudd3313)*beta2[ijk] + 
  (dA311*qPhysuudd1113 + dA322*qPhysuudd2213 + 
     2.*(dA312*qPhysuudd1213 + dA313*qPhysuudd1313 + dA323*qPhysuudd2313) + 
     dA333*qPhysuudd3313)*beta3[ijk]
;

rACABTF22
=
-((DACABTF22 + ACABTF22/r)*alpha[ijk]) + 
  (dA111*qPhysuudd1122 + dA122*qPhysuudd2222 + 
     2.*(dA112*qPhysuudd1222 + dA113*qPhysuudd1322 + 
        dA123*qPhysuudd2322) + dA133*qPhysuudd3322)*beta1[ijk] + 
  (dA211*qPhysuudd1122 + dA222*qPhysuudd2222 + 
     2.*(dA212*qPhysuudd1222 + dA213*qPhysuudd1322 + 
        dA223*qPhysuudd2322) + dA233*qPhysuudd3322)*beta2[ijk] + 
  (dA311*qPhysuudd1122 + dA322*qPhysuudd2222 + 
     2.*(dA312*qPhysuudd1222 + dA313*qPhysuudd1322 + dA323*qPhysuudd2322) + 
     dA333*qPhysuudd3322)*beta3[ijk]
;

rACABTF23
=
-((DACABTF32 + ACABTF32/r)*alpha[ijk]) + 
  (dA111*qPhysuudd1123 + dA122*qPhysuudd2223 + 
     2.*(dA112*qPhysuudd1223 + dA113*qPhysuudd1323 + 
        dA123*qPhysuudd2323) + dA133*qPhysuudd3323)*beta1[ijk] + 
  (dA211*qPhysuudd1123 + dA222*qPhysuudd2223 + 
     2.*(dA212*qPhysuudd1223 + dA213*qPhysuudd1323 + 
        dA223*qPhysuudd2323) + dA233*qPhysuudd3323)*beta2[ijk] + 
  (dA311*qPhysuudd1123 + dA322*qPhysuudd2223 + 
     2.*(dA312*qPhysuudd1223 + dA313*qPhysuudd1323 + dA323*qPhysuudd2323) + 
     dA333*qPhysuudd3323)*beta3[ijk]
;

rACABTF33
=
-((DACABTF33 + ACABTF33/r)*alpha[ijk]) + 
  (dA111*qPhysuudd1133 + dA122*qPhysuudd2233 + 
     2.*(dA112*qPhysuudd1233 + dA113*qPhysuudd1333 + 
        dA123*qPhysuudd2333) + dA133*qPhysuudd3333)*beta1[ijk] + 
  (dA211*qPhysuudd1133 + dA222*qPhysuudd2233 + 
     2.*(dA212*qPhysuudd1233 + dA213*qPhysuudd1333 + 
        dA223*qPhysuudd2333) + dA233*qPhysuudd3333)*beta2[ijk] + 
  (dA311*qPhysuudd1133 + dA322*qPhysuudd2233 + 
     2.*(dA312*qPhysuudd1233 + dA313*qPhysuudd1333 + dA323*qPhysuudd2333) + 
     dA333*qPhysuudd3333)*beta3[ijk]
;

rA11
=
rACABTF11 + 0.5*qdd11*rACqq + 2.*
   (qud11*rACsA1 + qud21*rACsA2 + qud31*rACsA3)*sdown1 + rACss*pow2(sdown1)
;

rA12
=
rACABTF12 + 0.5*qdd12*rACqq + (qud11*rACsA1 + qud21*rACsA2 + qud31*rACsA3)*
   sdown2 + sdown1*(qud12*rACsA1 + qud22*rACsA2 + qud32*rACsA3 + 
     rACss*sdown2)
;

rA13
=
rACABTF13 + 0.5*qdd13*rACqq + (qud11*rACsA1 + qud21*rACsA2 + qud31*rACsA3)*
   sdown3 + sdown1*(qud13*rACsA1 + qud23*rACsA2 + qud33*rACsA3 + 
     rACss*sdown3)
;

rA22
=
rACABTF22 + 0.5*qdd22*rACqq + 2.*
   (qud12*rACsA1 + qud22*rACsA2 + qud32*rACsA3)*sdown2 + rACss*pow2(sdown2)
;

rA23
=
rACABTF23 + 0.5*qdd23*rACqq + (qud12*rACsA1 + qud22*rACsA2 + qud32*rACsA3)*
   sdown3 + sdown2*(qud13*rACsA1 + qud23*rACsA2 + qud33*rACsA3 + 
     rACss*sdown3)
;

rA33
=
rACABTF33 + 0.5*qdd33*rACqq + 2.*
   (qud13*rACsA1 + qud23*rACsA2 + qud33*rACsA3)*sdown3 + rACss*pow2(sdown3)
;

rG1
=
qud11*rGamA1 + qud12*rGamA2 + qud13*rGamA3 + rGams*sup1
;

rG2
=
qud21*rGamA1 + qud22*rGamA2 + qud23*rGamA3 + rGams*sup2
;

rG3
=
qud31*rGamA1 + qud32*rGamA2 + qud33*rGamA3 + rGams*sup3
;

rB1
=
qud11*rBA1 + qud12*rBA2 + qud13*rBA3 + rBs*sup1
;

rB2
=
qud21*rBA1 + qud22*rBA2 + qud23*rBA3 + rBs*sup2
;

rB3
=
qud31*rBA1 + qud32*rBA2 + qud33*rBA3 + rBs*sup3
;


if (CheckForNANandINF(12, rA11,rA12,rA13,rA22,rA23,rA33,           
                       rB1,rB2,rB3,rG1,rG2,rK)) {
    printf("problem with RHS in bssn_boundary.m\n");
    printf("x=%2.5e, y=%2.5e, z=%2.5e\n",xp[ijk],yp[ijk],zp[ijk]);
  }
if (order_dissipation == 4 && boundaryNormore(2)) { 

rA11
=
rA11 - 2.*dissfactor*(oo2dx*(6.*A11[ijk] + A11[-2*di + ijk] - 
        4.*(A11[-di + ijk] + A11[di + ijk]) + A11[2*di + ijk]) + 
     oo2dy*(6.*A11[ijk] + A11[-2*dj + ijk] - 
        4.*(A11[-dj + ijk] + A11[dj + ijk]) + A11[2*dj + ijk]) + 
     oo2dz*(6.*A11[ijk] + A11[-2*dk + ijk] - 
        4.*(A11[-dk + ijk] + A11[dk + ijk]) + A11[2*dk + ijk]))
;

rA12
=
rA12 - 2.*dissfactor*(oo2dx*(6.*A12[ijk] + A12[-2*di + ijk] - 
        4.*(A12[-di + ijk] + A12[di + ijk]) + A12[2*di + ijk]) + 
     oo2dy*(6.*A12[ijk] + A12[-2*dj + ijk] - 
        4.*(A12[-dj + ijk] + A12[dj + ijk]) + A12[2*dj + ijk]) + 
     oo2dz*(6.*A12[ijk] + A12[-2*dk + ijk] - 
        4.*(A12[-dk + ijk] + A12[dk + ijk]) + A12[2*dk + ijk]))
;

rA13
=
rA13 - 2.*dissfactor*(oo2dx*(6.*A13[ijk] + A13[-2*di + ijk] - 
        4.*(A13[-di + ijk] + A13[di + ijk]) + A13[2*di + ijk]) + 
     oo2dy*(6.*A13[ijk] + A13[-2*dj + ijk] - 
        4.*(A13[-dj + ijk] + A13[dj + ijk]) + A13[2*dj + ijk]) + 
     oo2dz*(6.*A13[ijk] + A13[-2*dk + ijk] - 
        4.*(A13[-dk + ijk] + A13[dk + ijk]) + A13[2*dk + ijk]))
;

rA22
=
rA22 - 2.*dissfactor*(oo2dx*(6.*A22[ijk] + A22[-2*di + ijk] - 
        4.*(A22[-di + ijk] + A22[di + ijk]) + A22[2*di + ijk]) + 
     oo2dy*(6.*A22[ijk] + A22[-2*dj + ijk] - 
        4.*(A22[-dj + ijk] + A22[dj + ijk]) + A22[2*dj + ijk]) + 
     oo2dz*(6.*A22[ijk] + A22[-2*dk + ijk] - 
        4.*(A22[-dk + ijk] + A22[dk + ijk]) + A22[2*dk + ijk]))
;

rA23
=
rA23 - 2.*dissfactor*(oo2dx*(6.*A23[ijk] + A23[-2*di + ijk] - 
        4.*(A23[-di + ijk] + A23[di + ijk]) + A23[2*di + ijk]) + 
     oo2dy*(6.*A23[ijk] + A23[-2*dj + ijk] - 
        4.*(A23[-dj + ijk] + A23[dj + ijk]) + A23[2*dj + ijk]) + 
     oo2dz*(6.*A23[ijk] + A23[-2*dk + ijk] - 
        4.*(A23[-dk + ijk] + A23[dk + ijk]) + A23[2*dk + ijk]))
;

rA33
=
rA33 - 2.*dissfactor*(oo2dx*(6.*A33[ijk] + A33[-2*di + ijk] - 
        4.*(A33[-di + ijk] + A33[di + ijk]) + A33[2*di + ijk]) + 
     oo2dy*(6.*A33[ijk] + A33[-2*dj + ijk] - 
        4.*(A33[-dj + ijk] + A33[dj + ijk]) + A33[2*dj + ijk]) + 
     oo2dz*(6.*A33[ijk] + A33[-2*dk + ijk] - 
        4.*(A33[-dk + ijk] + A33[dk + ijk]) + A33[2*dk + ijk]))
;

rG1
=
rG1 - 2.*dissfactor*(oo2dx*(6.*G1[ijk] + G1[-2*di + ijk] - 
        4.*(G1[-di + ijk] + G1[di + ijk]) + G1[2*di + ijk]) + 
     oo2dy*(6.*G1[ijk] + G1[-2*dj + ijk] - 
        4.*(G1[-dj + ijk] + G1[dj + ijk]) + G1[2*dj + ijk]) + 
     oo2dz*(6.*G1[ijk] + G1[-2*dk + ijk] - 
        4.*(G1[-dk + ijk] + G1[dk + ijk]) + G1[2*dk + ijk]))
;

rG2
=
rG2 - 2.*dissfactor*(oo2dx*(6.*G2[ijk] + G2[-2*di + ijk] - 
        4.*(G2[-di + ijk] + G2[di + ijk]) + G2[2*di + ijk]) + 
     oo2dy*(6.*G2[ijk] + G2[-2*dj + ijk] - 
        4.*(G2[-dj + ijk] + G2[dj + ijk]) + G2[2*dj + ijk]) + 
     oo2dz*(6.*G2[ijk] + G2[-2*dk + ijk] - 
        4.*(G2[-dk + ijk] + G2[dk + ijk]) + G2[2*dk + ijk]))
;

rG3
=
rG3 - 2.*dissfactor*(oo2dx*(6.*G3[ijk] + G3[-2*di + ijk] - 
        4.*(G3[-di + ijk] + G3[di + ijk]) + G3[2*di + ijk]) + 
     oo2dy*(6.*G3[ijk] + G3[-2*dj + ijk] - 
        4.*(G3[-dj + ijk] + G3[dj + ijk]) + G3[2*dj + ijk]) + 
     oo2dz*(6.*G3[ijk] + G3[-2*dk + ijk] - 
        4.*(G3[-dk + ijk] + G3[dk + ijk]) + G3[2*dk + ijk]))
;

rK
=
rK - 2.*dissfactor*(oo2dx*(6.*K[ijk] + K[-2*di + ijk] - 
        4.*(K[-di + ijk] + K[di + ijk]) + K[2*di + ijk]) + 
     oo2dy*(6.*K[ijk] + K[-2*dj + ijk] - 4.*(K[-dj + ijk] + K[dj + ijk]) + 
        K[2*dj + ijk]) + oo2dz*
      (6.*K[ijk] + K[-2*dk + ijk] - 4.*(K[-dk + ijk] + K[dk + ijk]) + 
        K[2*dk + ijk]))
;

rB1
=
rB1 - 2.*dissfactor*(oo2dx*(6.*B1[ijk] + B1[-2*di + ijk] - 
        4.*(B1[-di + ijk] + B1[di + ijk]) + B1[2*di + ijk]) + 
     oo2dy*(6.*B1[ijk] + B1[-2*dj + ijk] - 
        4.*(B1[-dj + ijk] + B1[dj + ijk]) + B1[2*dj + ijk]) + 
     oo2dz*(6.*B1[ijk] + B1[-2*dk + ijk] - 
        4.*(B1[-dk + ijk] + B1[dk + ijk]) + B1[2*dk + ijk]))
;

rB2
=
rB2 - 2.*dissfactor*(oo2dx*(6.*B2[ijk] + B2[-2*di + ijk] - 
        4.*(B2[-di + ijk] + B2[di + ijk]) + B2[2*di + ijk]) + 
     oo2dy*(6.*B2[ijk] + B2[-2*dj + ijk] - 
        4.*(B2[-dj + ijk] + B2[dj + ijk]) + B2[2*dj + ijk]) + 
     oo2dz*(6.*B2[ijk] + B2[-2*dk + ijk] - 
        4.*(B2[-dk + ijk] + B2[dk + ijk]) + B2[2*dk + ijk]))
;

rB3
=
rB3 - 2.*dissfactor*(oo2dx*(6.*B3[ijk] + B3[-2*di + ijk] - 
        4.*(B3[-di + ijk] + B3[di + ijk]) + B3[2*di + ijk]) + 
     oo2dy*(6.*B3[ijk] + B3[-2*dj + ijk] - 
        4.*(B3[-dj + ijk] + B3[dj + ijk]) + B3[2*dj + ijk]) + 
     oo2dz*(6.*B3[ijk] + B3[-2*dk + ijk] - 
        4.*(B3[-dk + ijk] + B3[dk + ijk]) + B3[2*dk + ijk]))
;


} 


if (order_dissipation == 6 && boundaryNormore(3)) { 

rA11
=
rA11 + 2.*dissfactor*(oo2dx*(-20.*A11[ijk] + A11[-3*di + ijk] + 
        15.*(A11[-di + ijk] + A11[di + ijk]) - 
        6.*(A11[-2*di + ijk] + A11[2*di + ijk]) + A11[3*di + ijk]) + 
     oo2dy*(-20.*A11[ijk] + A11[-3*dj + ijk] + 
        15.*(A11[-dj + ijk] + A11[dj + ijk]) - 
        6.*(A11[-2*dj + ijk] + A11[2*dj + ijk]) + A11[3*dj + ijk]) + 
     oo2dz*(-20.*A11[ijk] + A11[-3*dk + ijk] + 
        15.*(A11[-dk + ijk] + A11[dk + ijk]) - 
        6.*(A11[-2*dk + ijk] + A11[2*dk + ijk]) + A11[3*dk + ijk]))
;

rA12
=
rA12 + 2.*dissfactor*(oo2dx*(-20.*A12[ijk] + A12[-3*di + ijk] + 
        15.*(A12[-di + ijk] + A12[di + ijk]) - 
        6.*(A12[-2*di + ijk] + A12[2*di + ijk]) + A12[3*di + ijk]) + 
     oo2dy*(-20.*A12[ijk] + A12[-3*dj + ijk] + 
        15.*(A12[-dj + ijk] + A12[dj + ijk]) - 
        6.*(A12[-2*dj + ijk] + A12[2*dj + ijk]) + A12[3*dj + ijk]) + 
     oo2dz*(-20.*A12[ijk] + A12[-3*dk + ijk] + 
        15.*(A12[-dk + ijk] + A12[dk + ijk]) - 
        6.*(A12[-2*dk + ijk] + A12[2*dk + ijk]) + A12[3*dk + ijk]))
;

rA13
=
rA13 + 2.*dissfactor*(oo2dx*(-20.*A13[ijk] + A13[-3*di + ijk] + 
        15.*(A13[-di + ijk] + A13[di + ijk]) - 
        6.*(A13[-2*di + ijk] + A13[2*di + ijk]) + A13[3*di + ijk]) + 
     oo2dy*(-20.*A13[ijk] + A13[-3*dj + ijk] + 
        15.*(A13[-dj + ijk] + A13[dj + ijk]) - 
        6.*(A13[-2*dj + ijk] + A13[2*dj + ijk]) + A13[3*dj + ijk]) + 
     oo2dz*(-20.*A13[ijk] + A13[-3*dk + ijk] + 
        15.*(A13[-dk + ijk] + A13[dk + ijk]) - 
        6.*(A13[-2*dk + ijk] + A13[2*dk + ijk]) + A13[3*dk + ijk]))
;

rA22
=
rA22 + 2.*dissfactor*(oo2dx*(-20.*A22[ijk] + A22[-3*di + ijk] + 
        15.*(A22[-di + ijk] + A22[di + ijk]) - 
        6.*(A22[-2*di + ijk] + A22[2*di + ijk]) + A22[3*di + ijk]) + 
     oo2dy*(-20.*A22[ijk] + A22[-3*dj + ijk] + 
        15.*(A22[-dj + ijk] + A22[dj + ijk]) - 
        6.*(A22[-2*dj + ijk] + A22[2*dj + ijk]) + A22[3*dj + ijk]) + 
     oo2dz*(-20.*A22[ijk] + A22[-3*dk + ijk] + 
        15.*(A22[-dk + ijk] + A22[dk + ijk]) - 
        6.*(A22[-2*dk + ijk] + A22[2*dk + ijk]) + A22[3*dk + ijk]))
;

rA23
=
rA23 + 2.*dissfactor*(oo2dx*(-20.*A23[ijk] + A23[-3*di + ijk] + 
        15.*(A23[-di + ijk] + A23[di + ijk]) - 
        6.*(A23[-2*di + ijk] + A23[2*di + ijk]) + A23[3*di + ijk]) + 
     oo2dy*(-20.*A23[ijk] + A23[-3*dj + ijk] + 
        15.*(A23[-dj + ijk] + A23[dj + ijk]) - 
        6.*(A23[-2*dj + ijk] + A23[2*dj + ijk]) + A23[3*dj + ijk]) + 
     oo2dz*(-20.*A23[ijk] + A23[-3*dk + ijk] + 
        15.*(A23[-dk + ijk] + A23[dk + ijk]) - 
        6.*(A23[-2*dk + ijk] + A23[2*dk + ijk]) + A23[3*dk + ijk]))
;

rA33
=
rA33 + 2.*dissfactor*(oo2dx*(-20.*A33[ijk] + A33[-3*di + ijk] + 
        15.*(A33[-di + ijk] + A33[di + ijk]) - 
        6.*(A33[-2*di + ijk] + A33[2*di + ijk]) + A33[3*di + ijk]) + 
     oo2dy*(-20.*A33[ijk] + A33[-3*dj + ijk] + 
        15.*(A33[-dj + ijk] + A33[dj + ijk]) - 
        6.*(A33[-2*dj + ijk] + A33[2*dj + ijk]) + A33[3*dj + ijk]) + 
     oo2dz*(-20.*A33[ijk] + A33[-3*dk + ijk] + 
        15.*(A33[-dk + ijk] + A33[dk + ijk]) - 
        6.*(A33[-2*dk + ijk] + A33[2*dk + ijk]) + A33[3*dk + ijk]))
;

rG1
=
rG1 + 2.*dissfactor*(oo2dx*(-20.*G1[ijk] + G1[-3*di + ijk] + 
        15.*(G1[-di + ijk] + G1[di + ijk]) - 
        6.*(G1[-2*di + ijk] + G1[2*di + ijk]) + G1[3*di + ijk]) + 
     oo2dy*(-20.*G1[ijk] + G1[-3*dj + ijk] + 
        15.*(G1[-dj + ijk] + G1[dj + ijk]) - 
        6.*(G1[-2*dj + ijk] + G1[2*dj + ijk]) + G1[3*dj + ijk]) + 
     oo2dz*(-20.*G1[ijk] + G1[-3*dk + ijk] + 
        15.*(G1[-dk + ijk] + G1[dk + ijk]) - 
        6.*(G1[-2*dk + ijk] + G1[2*dk + ijk]) + G1[3*dk + ijk]))
;

rG2
=
rG2 + 2.*dissfactor*(oo2dx*(-20.*G2[ijk] + G2[-3*di + ijk] + 
        15.*(G2[-di + ijk] + G2[di + ijk]) - 
        6.*(G2[-2*di + ijk] + G2[2*di + ijk]) + G2[3*di + ijk]) + 
     oo2dy*(-20.*G2[ijk] + G2[-3*dj + ijk] + 
        15.*(G2[-dj + ijk] + G2[dj + ijk]) - 
        6.*(G2[-2*dj + ijk] + G2[2*dj + ijk]) + G2[3*dj + ijk]) + 
     oo2dz*(-20.*G2[ijk] + G2[-3*dk + ijk] + 
        15.*(G2[-dk + ijk] + G2[dk + ijk]) - 
        6.*(G2[-2*dk + ijk] + G2[2*dk + ijk]) + G2[3*dk + ijk]))
;

rG3
=
rG3 + 2.*dissfactor*(oo2dx*(-20.*G3[ijk] + G3[-3*di + ijk] + 
        15.*(G3[-di + ijk] + G3[di + ijk]) - 
        6.*(G3[-2*di + ijk] + G3[2*di + ijk]) + G3[3*di + ijk]) + 
     oo2dy*(-20.*G3[ijk] + G3[-3*dj + ijk] + 
        15.*(G3[-dj + ijk] + G3[dj + ijk]) - 
        6.*(G3[-2*dj + ijk] + G3[2*dj + ijk]) + G3[3*dj + ijk]) + 
     oo2dz*(-20.*G3[ijk] + G3[-3*dk + ijk] + 
        15.*(G3[-dk + ijk] + G3[dk + ijk]) - 
        6.*(G3[-2*dk + ijk] + G3[2*dk + ijk]) + G3[3*dk + ijk]))
;

rK
=
rK + 2.*dissfactor*(oo2dx*(-20.*K[ijk] + K[-3*di + ijk] + 
        15.*(K[-di + ijk] + K[di + ijk]) - 
        6.*(K[-2*di + ijk] + K[2*di + ijk]) + K[3*di + ijk]) + 
     oo2dy*(-20.*K[ijk] + K[-3*dj + ijk] + 
        15.*(K[-dj + ijk] + K[dj + ijk]) - 
        6.*(K[-2*dj + ijk] + K[2*dj + ijk]) + K[3*dj + ijk]) + 
     oo2dz*(-20.*K[ijk] + K[-3*dk + ijk] + 
        15.*(K[-dk + ijk] + K[dk + ijk]) - 
        6.*(K[-2*dk + ijk] + K[2*dk + ijk]) + K[3*dk + ijk]))
;

rB1
=
rB1 + 2.*dissfactor*(oo2dx*(-20.*B1[ijk] + B1[-3*di + ijk] + 
        15.*(B1[-di + ijk] + B1[di + ijk]) - 
        6.*(B1[-2*di + ijk] + B1[2*di + ijk]) + B1[3*di + ijk]) + 
     oo2dy*(-20.*B1[ijk] + B1[-3*dj + ijk] + 
        15.*(B1[-dj + ijk] + B1[dj + ijk]) - 
        6.*(B1[-2*dj + ijk] + B1[2*dj + ijk]) + B1[3*dj + ijk]) + 
     oo2dz*(-20.*B1[ijk] + B1[-3*dk + ijk] + 
        15.*(B1[-dk + ijk] + B1[dk + ijk]) - 
        6.*(B1[-2*dk + ijk] + B1[2*dk + ijk]) + B1[3*dk + ijk]))
;

rB2
=
rB2 + 2.*dissfactor*(oo2dx*(-20.*B2[ijk] + B2[-3*di + ijk] + 
        15.*(B2[-di + ijk] + B2[di + ijk]) - 
        6.*(B2[-2*di + ijk] + B2[2*di + ijk]) + B2[3*di + ijk]) + 
     oo2dy*(-20.*B2[ijk] + B2[-3*dj + ijk] + 
        15.*(B2[-dj + ijk] + B2[dj + ijk]) - 
        6.*(B2[-2*dj + ijk] + B2[2*dj + ijk]) + B2[3*dj + ijk]) + 
     oo2dz*(-20.*B2[ijk] + B2[-3*dk + ijk] + 
        15.*(B2[-dk + ijk] + B2[dk + ijk]) - 
        6.*(B2[-2*dk + ijk] + B2[2*dk + ijk]) + B2[3*dk + ijk]))
;

rB3
=
rB3 + 2.*dissfactor*(oo2dx*(-20.*B3[ijk] + B3[-3*di + ijk] + 
        15.*(B3[-di + ijk] + B3[di + ijk]) - 
        6.*(B3[-2*di + ijk] + B3[2*di + ijk]) + B3[3*di + ijk]) + 
     oo2dy*(-20.*B3[ijk] + B3[-3*dj + ijk] + 
        15.*(B3[-dj + ijk] + B3[dj + ijk]) - 
        6.*(B3[-2*dj + ijk] + B3[2*dj + ijk]) + B3[3*dj + ijk]) + 
     oo2dz*(-20.*B3[ijk] + B3[-3*dk + ijk] + 
        15.*(B3[-dk + ijk] + B3[dk + ijk]) - 
        6.*(B3[-2*dk + ijk] + B3[2*dk + ijk]) + B3[3*dk + ijk]))
;


} 


if (order_dissipation == 8 && boundaryNormore(4)) { 

rA11
=
rA11 - 2.*dissfactor*(oo2dx*(70.*A11[ijk] + A11[-4*di + ijk] - 
        56.*(A11[-di + ijk] + A11[di + ijk]) + 
        28.*(A11[-2*di + ijk] + A11[2*di + ijk]) - 
        8.*(A11[-3*di + ijk] + A11[3*di + ijk]) + A11[4*di + ijk]) + 
     oo2dy*(70.*A11[ijk] + A11[-4*dj + ijk] - 
        56.*(A11[-dj + ijk] + A11[dj + ijk]) + 
        28.*(A11[-2*dj + ijk] + A11[2*dj + ijk]) - 
        8.*(A11[-3*dj + ijk] + A11[3*dj + ijk]) + A11[4*dj + ijk]) + 
     oo2dz*(70.*A11[ijk] + A11[-4*dk + ijk] - 
        56.*(A11[-dk + ijk] + A11[dk + ijk]) + 
        28.*(A11[-2*dk + ijk] + A11[2*dk + ijk]) - 
        8.*(A11[-3*dk + ijk] + A11[3*dk + ijk]) + A11[4*dk + ijk]))
;

rA12
=
rA12 - 2.*dissfactor*(oo2dx*(70.*A12[ijk] + A12[-4*di + ijk] - 
        56.*(A12[-di + ijk] + A12[di + ijk]) + 
        28.*(A12[-2*di + ijk] + A12[2*di + ijk]) - 
        8.*(A12[-3*di + ijk] + A12[3*di + ijk]) + A12[4*di + ijk]) + 
     oo2dy*(70.*A12[ijk] + A12[-4*dj + ijk] - 
        56.*(A12[-dj + ijk] + A12[dj + ijk]) + 
        28.*(A12[-2*dj + ijk] + A12[2*dj + ijk]) - 
        8.*(A12[-3*dj + ijk] + A12[3*dj + ijk]) + A12[4*dj + ijk]) + 
     oo2dz*(70.*A12[ijk] + A12[-4*dk + ijk] - 
        56.*(A12[-dk + ijk] + A12[dk + ijk]) + 
        28.*(A12[-2*dk + ijk] + A12[2*dk + ijk]) - 
        8.*(A12[-3*dk + ijk] + A12[3*dk + ijk]) + A12[4*dk + ijk]))
;

rA13
=
rA13 - 2.*dissfactor*(oo2dx*(70.*A13[ijk] + A13[-4*di + ijk] - 
        56.*(A13[-di + ijk] + A13[di + ijk]) + 
        28.*(A13[-2*di + ijk] + A13[2*di + ijk]) - 
        8.*(A13[-3*di + ijk] + A13[3*di + ijk]) + A13[4*di + ijk]) + 
     oo2dy*(70.*A13[ijk] + A13[-4*dj + ijk] - 
        56.*(A13[-dj + ijk] + A13[dj + ijk]) + 
        28.*(A13[-2*dj + ijk] + A13[2*dj + ijk]) - 
        8.*(A13[-3*dj + ijk] + A13[3*dj + ijk]) + A13[4*dj + ijk]) + 
     oo2dz*(70.*A13[ijk] + A13[-4*dk + ijk] - 
        56.*(A13[-dk + ijk] + A13[dk + ijk]) + 
        28.*(A13[-2*dk + ijk] + A13[2*dk + ijk]) - 
        8.*(A13[-3*dk + ijk] + A13[3*dk + ijk]) + A13[4*dk + ijk]))
;

rA22
=
rA22 - 2.*dissfactor*(oo2dx*(70.*A22[ijk] + A22[-4*di + ijk] - 
        56.*(A22[-di + ijk] + A22[di + ijk]) + 
        28.*(A22[-2*di + ijk] + A22[2*di + ijk]) - 
        8.*(A22[-3*di + ijk] + A22[3*di + ijk]) + A22[4*di + ijk]) + 
     oo2dy*(70.*A22[ijk] + A22[-4*dj + ijk] - 
        56.*(A22[-dj + ijk] + A22[dj + ijk]) + 
        28.*(A22[-2*dj + ijk] + A22[2*dj + ijk]) - 
        8.*(A22[-3*dj + ijk] + A22[3*dj + ijk]) + A22[4*dj + ijk]) + 
     oo2dz*(70.*A22[ijk] + A22[-4*dk + ijk] - 
        56.*(A22[-dk + ijk] + A22[dk + ijk]) + 
        28.*(A22[-2*dk + ijk] + A22[2*dk + ijk]) - 
        8.*(A22[-3*dk + ijk] + A22[3*dk + ijk]) + A22[4*dk + ijk]))
;

rA23
=
rA23 - 2.*dissfactor*(oo2dx*(70.*A23[ijk] + A23[-4*di + ijk] - 
        56.*(A23[-di + ijk] + A23[di + ijk]) + 
        28.*(A23[-2*di + ijk] + A23[2*di + ijk]) - 
        8.*(A23[-3*di + ijk] + A23[3*di + ijk]) + A23[4*di + ijk]) + 
     oo2dy*(70.*A23[ijk] + A23[-4*dj + ijk] - 
        56.*(A23[-dj + ijk] + A23[dj + ijk]) + 
        28.*(A23[-2*dj + ijk] + A23[2*dj + ijk]) - 
        8.*(A23[-3*dj + ijk] + A23[3*dj + ijk]) + A23[4*dj + ijk]) + 
     oo2dz*(70.*A23[ijk] + A23[-4*dk + ijk] - 
        56.*(A23[-dk + ijk] + A23[dk + ijk]) + 
        28.*(A23[-2*dk + ijk] + A23[2*dk + ijk]) - 
        8.*(A23[-3*dk + ijk] + A23[3*dk + ijk]) + A23[4*dk + ijk]))
;

rA33
=
rA33 - 2.*dissfactor*(oo2dx*(70.*A33[ijk] + A33[-4*di + ijk] - 
        56.*(A33[-di + ijk] + A33[di + ijk]) + 
        28.*(A33[-2*di + ijk] + A33[2*di + ijk]) - 
        8.*(A33[-3*di + ijk] + A33[3*di + ijk]) + A33[4*di + ijk]) + 
     oo2dy*(70.*A33[ijk] + A33[-4*dj + ijk] - 
        56.*(A33[-dj + ijk] + A33[dj + ijk]) + 
        28.*(A33[-2*dj + ijk] + A33[2*dj + ijk]) - 
        8.*(A33[-3*dj + ijk] + A33[3*dj + ijk]) + A33[4*dj + ijk]) + 
     oo2dz*(70.*A33[ijk] + A33[-4*dk + ijk] - 
        56.*(A33[-dk + ijk] + A33[dk + ijk]) + 
        28.*(A33[-2*dk + ijk] + A33[2*dk + ijk]) - 
        8.*(A33[-3*dk + ijk] + A33[3*dk + ijk]) + A33[4*dk + ijk]))
;

rG1
=
rG1 - 2.*dissfactor*(oo2dx*(70.*G1[ijk] + G1[-4*di + ijk] - 
        56.*(G1[-di + ijk] + G1[di + ijk]) + 
        28.*(G1[-2*di + ijk] + G1[2*di + ijk]) - 
        8.*(G1[-3*di + ijk] + G1[3*di + ijk]) + G1[4*di + ijk]) + 
     oo2dy*(70.*G1[ijk] + G1[-4*dj + ijk] - 
        56.*(G1[-dj + ijk] + G1[dj + ijk]) + 
        28.*(G1[-2*dj + ijk] + G1[2*dj + ijk]) - 
        8.*(G1[-3*dj + ijk] + G1[3*dj + ijk]) + G1[4*dj + ijk]) + 
     oo2dz*(70.*G1[ijk] + G1[-4*dk + ijk] - 
        56.*(G1[-dk + ijk] + G1[dk + ijk]) + 
        28.*(G1[-2*dk + ijk] + G1[2*dk + ijk]) - 
        8.*(G1[-3*dk + ijk] + G1[3*dk + ijk]) + G1[4*dk + ijk]))
;

rG2
=
rG2 - 2.*dissfactor*(oo2dx*(70.*G2[ijk] + G2[-4*di + ijk] - 
        56.*(G2[-di + ijk] + G2[di + ijk]) + 
        28.*(G2[-2*di + ijk] + G2[2*di + ijk]) - 
        8.*(G2[-3*di + ijk] + G2[3*di + ijk]) + G2[4*di + ijk]) + 
     oo2dy*(70.*G2[ijk] + G2[-4*dj + ijk] - 
        56.*(G2[-dj + ijk] + G2[dj + ijk]) + 
        28.*(G2[-2*dj + ijk] + G2[2*dj + ijk]) - 
        8.*(G2[-3*dj + ijk] + G2[3*dj + ijk]) + G2[4*dj + ijk]) + 
     oo2dz*(70.*G2[ijk] + G2[-4*dk + ijk] - 
        56.*(G2[-dk + ijk] + G2[dk + ijk]) + 
        28.*(G2[-2*dk + ijk] + G2[2*dk + ijk]) - 
        8.*(G2[-3*dk + ijk] + G2[3*dk + ijk]) + G2[4*dk + ijk]))
;

rG3
=
rG3 - 2.*dissfactor*(oo2dx*(70.*G3[ijk] + G3[-4*di + ijk] - 
        56.*(G3[-di + ijk] + G3[di + ijk]) + 
        28.*(G3[-2*di + ijk] + G3[2*di + ijk]) - 
        8.*(G3[-3*di + ijk] + G3[3*di + ijk]) + G3[4*di + ijk]) + 
     oo2dy*(70.*G3[ijk] + G3[-4*dj + ijk] - 
        56.*(G3[-dj + ijk] + G3[dj + ijk]) + 
        28.*(G3[-2*dj + ijk] + G3[2*dj + ijk]) - 
        8.*(G3[-3*dj + ijk] + G3[3*dj + ijk]) + G3[4*dj + ijk]) + 
     oo2dz*(70.*G3[ijk] + G3[-4*dk + ijk] - 
        56.*(G3[-dk + ijk] + G3[dk + ijk]) + 
        28.*(G3[-2*dk + ijk] + G3[2*dk + ijk]) - 
        8.*(G3[-3*dk + ijk] + G3[3*dk + ijk]) + G3[4*dk + ijk]))
;

rK
=
rK - 2.*dissfactor*(oo2dx*(70.*K[ijk] + K[-4*di + ijk] - 
        56.*(K[-di + ijk] + K[di + ijk]) + 
        28.*(K[-2*di + ijk] + K[2*di + ijk]) - 
        8.*(K[-3*di + ijk] + K[3*di + ijk]) + K[4*di + ijk]) + 
     oo2dy*(70.*K[ijk] + K[-4*dj + ijk] - 
        56.*(K[-dj + ijk] + K[dj + ijk]) + 
        28.*(K[-2*dj + ijk] + K[2*dj + ijk]) - 
        8.*(K[-3*dj + ijk] + K[3*dj + ijk]) + K[4*dj + ijk]) + 
     oo2dz*(70.*K[ijk] + K[-4*dk + ijk] - 56.*(K[-dk + ijk] + K[dk + ijk]) + 
        28.*(K[-2*dk + ijk] + K[2*dk + ijk]) - 
        8.*(K[-3*dk + ijk] + K[3*dk + ijk]) + K[4*dk + ijk]))
;

rB1
=
rB1 - 2.*dissfactor*(oo2dx*(70.*B1[ijk] + B1[-4*di + ijk] - 
        56.*(B1[-di + ijk] + B1[di + ijk]) + 
        28.*(B1[-2*di + ijk] + B1[2*di + ijk]) - 
        8.*(B1[-3*di + ijk] + B1[3*di + ijk]) + B1[4*di + ijk]) + 
     oo2dy*(70.*B1[ijk] + B1[-4*dj + ijk] - 
        56.*(B1[-dj + ijk] + B1[dj + ijk]) + 
        28.*(B1[-2*dj + ijk] + B1[2*dj + ijk]) - 
        8.*(B1[-3*dj + ijk] + B1[3*dj + ijk]) + B1[4*dj + ijk]) + 
     oo2dz*(70.*B1[ijk] + B1[-4*dk + ijk] - 
        56.*(B1[-dk + ijk] + B1[dk + ijk]) + 
        28.*(B1[-2*dk + ijk] + B1[2*dk + ijk]) - 
        8.*(B1[-3*dk + ijk] + B1[3*dk + ijk]) + B1[4*dk + ijk]))
;

rB2
=
rB2 - 2.*dissfactor*(oo2dx*(70.*B2[ijk] + B2[-4*di + ijk] - 
        56.*(B2[-di + ijk] + B2[di + ijk]) + 
        28.*(B2[-2*di + ijk] + B2[2*di + ijk]) - 
        8.*(B2[-3*di + ijk] + B2[3*di + ijk]) + B2[4*di + ijk]) + 
     oo2dy*(70.*B2[ijk] + B2[-4*dj + ijk] - 
        56.*(B2[-dj + ijk] + B2[dj + ijk]) + 
        28.*(B2[-2*dj + ijk] + B2[2*dj + ijk]) - 
        8.*(B2[-3*dj + ijk] + B2[3*dj + ijk]) + B2[4*dj + ijk]) + 
     oo2dz*(70.*B2[ijk] + B2[-4*dk + ijk] - 
        56.*(B2[-dk + ijk] + B2[dk + ijk]) + 
        28.*(B2[-2*dk + ijk] + B2[2*dk + ijk]) - 
        8.*(B2[-3*dk + ijk] + B2[3*dk + ijk]) + B2[4*dk + ijk]))
;

rB3
=
rB3 - 2.*dissfactor*(oo2dx*(70.*B3[ijk] + B3[-4*di + ijk] - 
        56.*(B3[-di + ijk] + B3[di + ijk]) + 
        28.*(B3[-2*di + ijk] + B3[2*di + ijk]) - 
        8.*(B3[-3*di + ijk] + B3[3*di + ijk]) + B3[4*di + ijk]) + 
     oo2dy*(70.*B3[ijk] + B3[-4*dj + ijk] - 
        56.*(B3[-dj + ijk] + B3[dj + ijk]) + 
        28.*(B3[-2*dj + ijk] + B3[2*dj + ijk]) - 
        8.*(B3[-3*dj + ijk] + B3[3*dj + ijk]) + B3[4*dj + ijk]) + 
     oo2dz*(70.*B3[ijk] + B3[-4*dk + ijk] - 
        56.*(B3[-dk + ijk] + B3[dk + ijk]) + 
        28.*(B3[-2*dk + ijk] + B3[2*dk + ijk]) - 
        8.*(B3[-3*dk + ijk] + B3[3*dk + ijk]) + B3[4*dk + ijk]))
;


} 


if (order_dissipation == 10 && boundaryNormore(5)) { 

rA11
=
rA11 + 2.*dissfactor*(oo2dx*(-252.*A11[ijk] + A11[-5*di + ijk] + 
        210.*(A11[-di + ijk] + A11[di + ijk]) - 
        120.*(A11[-2*di + ijk] + A11[2*di + ijk]) + 
        45.*(A11[-3*di + ijk] + A11[3*di + ijk]) - 
        10.*(A11[-4*di + ijk] + A11[4*di + ijk]) + A11[5*di + ijk]) + 
     oo2dy*(-252.*A11[ijk] + A11[-5*dj + ijk] + 
        210.*(A11[-dj + ijk] + A11[dj + ijk]) - 
        120.*(A11[-2*dj + ijk] + A11[2*dj + ijk]) + 
        45.*(A11[-3*dj + ijk] + A11[3*dj + ijk]) - 
        10.*(A11[-4*dj + ijk] + A11[4*dj + ijk]) + A11[5*dj + ijk]) + 
     oo2dz*(-252.*A11[ijk] + A11[-5*dk + ijk] + 
        210.*(A11[-dk + ijk] + A11[dk + ijk]) - 
        120.*(A11[-2*dk + ijk] + A11[2*dk + ijk]) + 
        45.*(A11[-3*dk + ijk] + A11[3*dk + ijk]) - 
        10.*(A11[-4*dk + ijk] + A11[4*dk + ijk]) + A11[5*dk + ijk]))
;

rA12
=
rA12 + 2.*dissfactor*(oo2dx*(-252.*A12[ijk] + A12[-5*di + ijk] + 
        210.*(A12[-di + ijk] + A12[di + ijk]) - 
        120.*(A12[-2*di + ijk] + A12[2*di + ijk]) + 
        45.*(A12[-3*di + ijk] + A12[3*di + ijk]) - 
        10.*(A12[-4*di + ijk] + A12[4*di + ijk]) + A12[5*di + ijk]) + 
     oo2dy*(-252.*A12[ijk] + A12[-5*dj + ijk] + 
        210.*(A12[-dj + ijk] + A12[dj + ijk]) - 
        120.*(A12[-2*dj + ijk] + A12[2*dj + ijk]) + 
        45.*(A12[-3*dj + ijk] + A12[3*dj + ijk]) - 
        10.*(A12[-4*dj + ijk] + A12[4*dj + ijk]) + A12[5*dj + ijk]) + 
     oo2dz*(-252.*A12[ijk] + A12[-5*dk + ijk] + 
        210.*(A12[-dk + ijk] + A12[dk + ijk]) - 
        120.*(A12[-2*dk + ijk] + A12[2*dk + ijk]) + 
        45.*(A12[-3*dk + ijk] + A12[3*dk + ijk]) - 
        10.*(A12[-4*dk + ijk] + A12[4*dk + ijk]) + A12[5*dk + ijk]))
;

rA13
=
rA13 + 2.*dissfactor*(oo2dx*(-252.*A13[ijk] + A13[-5*di + ijk] + 
        210.*(A13[-di + ijk] + A13[di + ijk]) - 
        120.*(A13[-2*di + ijk] + A13[2*di + ijk]) + 
        45.*(A13[-3*di + ijk] + A13[3*di + ijk]) - 
        10.*(A13[-4*di + ijk] + A13[4*di + ijk]) + A13[5*di + ijk]) + 
     oo2dy*(-252.*A13[ijk] + A13[-5*dj + ijk] + 
        210.*(A13[-dj + ijk] + A13[dj + ijk]) - 
        120.*(A13[-2*dj + ijk] + A13[2*dj + ijk]) + 
        45.*(A13[-3*dj + ijk] + A13[3*dj + ijk]) - 
        10.*(A13[-4*dj + ijk] + A13[4*dj + ijk]) + A13[5*dj + ijk]) + 
     oo2dz*(-252.*A13[ijk] + A13[-5*dk + ijk] + 
        210.*(A13[-dk + ijk] + A13[dk + ijk]) - 
        120.*(A13[-2*dk + ijk] + A13[2*dk + ijk]) + 
        45.*(A13[-3*dk + ijk] + A13[3*dk + ijk]) - 
        10.*(A13[-4*dk + ijk] + A13[4*dk + ijk]) + A13[5*dk + ijk]))
;

rA22
=
rA22 + 2.*dissfactor*(oo2dx*(-252.*A22[ijk] + A22[-5*di + ijk] + 
        210.*(A22[-di + ijk] + A22[di + ijk]) - 
        120.*(A22[-2*di + ijk] + A22[2*di + ijk]) + 
        45.*(A22[-3*di + ijk] + A22[3*di + ijk]) - 
        10.*(A22[-4*di + ijk] + A22[4*di + ijk]) + A22[5*di + ijk]) + 
     oo2dy*(-252.*A22[ijk] + A22[-5*dj + ijk] + 
        210.*(A22[-dj + ijk] + A22[dj + ijk]) - 
        120.*(A22[-2*dj + ijk] + A22[2*dj + ijk]) + 
        45.*(A22[-3*dj + ijk] + A22[3*dj + ijk]) - 
        10.*(A22[-4*dj + ijk] + A22[4*dj + ijk]) + A22[5*dj + ijk]) + 
     oo2dz*(-252.*A22[ijk] + A22[-5*dk + ijk] + 
        210.*(A22[-dk + ijk] + A22[dk + ijk]) - 
        120.*(A22[-2*dk + ijk] + A22[2*dk + ijk]) + 
        45.*(A22[-3*dk + ijk] + A22[3*dk + ijk]) - 
        10.*(A22[-4*dk + ijk] + A22[4*dk + ijk]) + A22[5*dk + ijk]))
;

rA23
=
rA23 + 2.*dissfactor*(oo2dx*(-252.*A23[ijk] + A23[-5*di + ijk] + 
        210.*(A23[-di + ijk] + A23[di + ijk]) - 
        120.*(A23[-2*di + ijk] + A23[2*di + ijk]) + 
        45.*(A23[-3*di + ijk] + A23[3*di + ijk]) - 
        10.*(A23[-4*di + ijk] + A23[4*di + ijk]) + A23[5*di + ijk]) + 
     oo2dy*(-252.*A23[ijk] + A23[-5*dj + ijk] + 
        210.*(A23[-dj + ijk] + A23[dj + ijk]) - 
        120.*(A23[-2*dj + ijk] + A23[2*dj + ijk]) + 
        45.*(A23[-3*dj + ijk] + A23[3*dj + ijk]) - 
        10.*(A23[-4*dj + ijk] + A23[4*dj + ijk]) + A23[5*dj + ijk]) + 
     oo2dz*(-252.*A23[ijk] + A23[-5*dk + ijk] + 
        210.*(A23[-dk + ijk] + A23[dk + ijk]) - 
        120.*(A23[-2*dk + ijk] + A23[2*dk + ijk]) + 
        45.*(A23[-3*dk + ijk] + A23[3*dk + ijk]) - 
        10.*(A23[-4*dk + ijk] + A23[4*dk + ijk]) + A23[5*dk + ijk]))
;

rA33
=
rA33 + 2.*dissfactor*(oo2dx*(-252.*A33[ijk] + A33[-5*di + ijk] + 
        210.*(A33[-di + ijk] + A33[di + ijk]) - 
        120.*(A33[-2*di + ijk] + A33[2*di + ijk]) + 
        45.*(A33[-3*di + ijk] + A33[3*di + ijk]) - 
        10.*(A33[-4*di + ijk] + A33[4*di + ijk]) + A33[5*di + ijk]) + 
     oo2dy*(-252.*A33[ijk] + A33[-5*dj + ijk] + 
        210.*(A33[-dj + ijk] + A33[dj + ijk]) - 
        120.*(A33[-2*dj + ijk] + A33[2*dj + ijk]) + 
        45.*(A33[-3*dj + ijk] + A33[3*dj + ijk]) - 
        10.*(A33[-4*dj + ijk] + A33[4*dj + ijk]) + A33[5*dj + ijk]) + 
     oo2dz*(-252.*A33[ijk] + A33[-5*dk + ijk] + 
        210.*(A33[-dk + ijk] + A33[dk + ijk]) - 
        120.*(A33[-2*dk + ijk] + A33[2*dk + ijk]) + 
        45.*(A33[-3*dk + ijk] + A33[3*dk + ijk]) - 
        10.*(A33[-4*dk + ijk] + A33[4*dk + ijk]) + A33[5*dk + ijk]))
;

rG1
=
rG1 + 2.*dissfactor*(oo2dx*(-252.*G1[ijk] + G1[-5*di + ijk] + 
        210.*(G1[-di + ijk] + G1[di + ijk]) - 
        120.*(G1[-2*di + ijk] + G1[2*di + ijk]) + 
        45.*(G1[-3*di + ijk] + G1[3*di + ijk]) - 
        10.*(G1[-4*di + ijk] + G1[4*di + ijk]) + G1[5*di + ijk]) + 
     oo2dy*(-252.*G1[ijk] + G1[-5*dj + ijk] + 
        210.*(G1[-dj + ijk] + G1[dj + ijk]) - 
        120.*(G1[-2*dj + ijk] + G1[2*dj + ijk]) + 
        45.*(G1[-3*dj + ijk] + G1[3*dj + ijk]) - 
        10.*(G1[-4*dj + ijk] + G1[4*dj + ijk]) + G1[5*dj + ijk]) + 
     oo2dz*(-252.*G1[ijk] + G1[-5*dk + ijk] + 
        210.*(G1[-dk + ijk] + G1[dk + ijk]) - 
        120.*(G1[-2*dk + ijk] + G1[2*dk + ijk]) + 
        45.*(G1[-3*dk + ijk] + G1[3*dk + ijk]) - 
        10.*(G1[-4*dk + ijk] + G1[4*dk + ijk]) + G1[5*dk + ijk]))
;

rG2
=
rG2 + 2.*dissfactor*(oo2dx*(-252.*G2[ijk] + G2[-5*di + ijk] + 
        210.*(G2[-di + ijk] + G2[di + ijk]) - 
        120.*(G2[-2*di + ijk] + G2[2*di + ijk]) + 
        45.*(G2[-3*di + ijk] + G2[3*di + ijk]) - 
        10.*(G2[-4*di + ijk] + G2[4*di + ijk]) + G2[5*di + ijk]) + 
     oo2dy*(-252.*G2[ijk] + G2[-5*dj + ijk] + 
        210.*(G2[-dj + ijk] + G2[dj + ijk]) - 
        120.*(G2[-2*dj + ijk] + G2[2*dj + ijk]) + 
        45.*(G2[-3*dj + ijk] + G2[3*dj + ijk]) - 
        10.*(G2[-4*dj + ijk] + G2[4*dj + ijk]) + G2[5*dj + ijk]) + 
     oo2dz*(-252.*G2[ijk] + G2[-5*dk + ijk] + 
        210.*(G2[-dk + ijk] + G2[dk + ijk]) - 
        120.*(G2[-2*dk + ijk] + G2[2*dk + ijk]) + 
        45.*(G2[-3*dk + ijk] + G2[3*dk + ijk]) - 
        10.*(G2[-4*dk + ijk] + G2[4*dk + ijk]) + G2[5*dk + ijk]))
;

rG3
=
rG3 + 2.*dissfactor*(oo2dx*(-252.*G3[ijk] + G3[-5*di + ijk] + 
        210.*(G3[-di + ijk] + G3[di + ijk]) - 
        120.*(G3[-2*di + ijk] + G3[2*di + ijk]) + 
        45.*(G3[-3*di + ijk] + G3[3*di + ijk]) - 
        10.*(G3[-4*di + ijk] + G3[4*di + ijk]) + G3[5*di + ijk]) + 
     oo2dy*(-252.*G3[ijk] + G3[-5*dj + ijk] + 
        210.*(G3[-dj + ijk] + G3[dj + ijk]) - 
        120.*(G3[-2*dj + ijk] + G3[2*dj + ijk]) + 
        45.*(G3[-3*dj + ijk] + G3[3*dj + ijk]) - 
        10.*(G3[-4*dj + ijk] + G3[4*dj + ijk]) + G3[5*dj + ijk]) + 
     oo2dz*(-252.*G3[ijk] + G3[-5*dk + ijk] + 
        210.*(G3[-dk + ijk] + G3[dk + ijk]) - 
        120.*(G3[-2*dk + ijk] + G3[2*dk + ijk]) + 
        45.*(G3[-3*dk + ijk] + G3[3*dk + ijk]) - 
        10.*(G3[-4*dk + ijk] + G3[4*dk + ijk]) + G3[5*dk + ijk]))
;

rK
=
rK + 2.*dissfactor*(oo2dx*(-252.*K[ijk] + K[-5*di + ijk] + 
        210.*(K[-di + ijk] + K[di + ijk]) - 
        120.*(K[-2*di + ijk] + K[2*di + ijk]) + 
        45.*(K[-3*di + ijk] + K[3*di + ijk]) - 
        10.*(K[-4*di + ijk] + K[4*di + ijk]) + K[5*di + ijk]) + 
     oo2dy*(-252.*K[ijk] + K[-5*dj + ijk] + 
        210.*(K[-dj + ijk] + K[dj + ijk]) - 
        120.*(K[-2*dj + ijk] + K[2*dj + ijk]) + 
        45.*(K[-3*dj + ijk] + K[3*dj + ijk]) - 
        10.*(K[-4*dj + ijk] + K[4*dj + ijk]) + K[5*dj + ijk]) + 
     oo2dz*(-252.*K[ijk] + K[-5*dk + ijk] + 
        210.*(K[-dk + ijk] + K[dk + ijk]) - 
        120.*(K[-2*dk + ijk] + K[2*dk + ijk]) + 
        45.*(K[-3*dk + ijk] + K[3*dk + ijk]) - 
        10.*(K[-4*dk + ijk] + K[4*dk + ijk]) + K[5*dk + ijk]))
;

rB1
=
rB1 + 2.*dissfactor*(oo2dx*(-252.*B1[ijk] + B1[-5*di + ijk] + 
        210.*(B1[-di + ijk] + B1[di + ijk]) - 
        120.*(B1[-2*di + ijk] + B1[2*di + ijk]) + 
        45.*(B1[-3*di + ijk] + B1[3*di + ijk]) - 
        10.*(B1[-4*di + ijk] + B1[4*di + ijk]) + B1[5*di + ijk]) + 
     oo2dy*(-252.*B1[ijk] + B1[-5*dj + ijk] + 
        210.*(B1[-dj + ijk] + B1[dj + ijk]) - 
        120.*(B1[-2*dj + ijk] + B1[2*dj + ijk]) + 
        45.*(B1[-3*dj + ijk] + B1[3*dj + ijk]) - 
        10.*(B1[-4*dj + ijk] + B1[4*dj + ijk]) + B1[5*dj + ijk]) + 
     oo2dz*(-252.*B1[ijk] + B1[-5*dk + ijk] + 
        210.*(B1[-dk + ijk] + B1[dk + ijk]) - 
        120.*(B1[-2*dk + ijk] + B1[2*dk + ijk]) + 
        45.*(B1[-3*dk + ijk] + B1[3*dk + ijk]) - 
        10.*(B1[-4*dk + ijk] + B1[4*dk + ijk]) + B1[5*dk + ijk]))
;

rB2
=
rB2 + 2.*dissfactor*(oo2dx*(-252.*B2[ijk] + B2[-5*di + ijk] + 
        210.*(B2[-di + ijk] + B2[di + ijk]) - 
        120.*(B2[-2*di + ijk] + B2[2*di + ijk]) + 
        45.*(B2[-3*di + ijk] + B2[3*di + ijk]) - 
        10.*(B2[-4*di + ijk] + B2[4*di + ijk]) + B2[5*di + ijk]) + 
     oo2dy*(-252.*B2[ijk] + B2[-5*dj + ijk] + 
        210.*(B2[-dj + ijk] + B2[dj + ijk]) - 
        120.*(B2[-2*dj + ijk] + B2[2*dj + ijk]) + 
        45.*(B2[-3*dj + ijk] + B2[3*dj + ijk]) - 
        10.*(B2[-4*dj + ijk] + B2[4*dj + ijk]) + B2[5*dj + ijk]) + 
     oo2dz*(-252.*B2[ijk] + B2[-5*dk + ijk] + 
        210.*(B2[-dk + ijk] + B2[dk + ijk]) - 
        120.*(B2[-2*dk + ijk] + B2[2*dk + ijk]) + 
        45.*(B2[-3*dk + ijk] + B2[3*dk + ijk]) - 
        10.*(B2[-4*dk + ijk] + B2[4*dk + ijk]) + B2[5*dk + ijk]))
;

rB3
=
rB3 + 2.*dissfactor*(oo2dx*(-252.*B3[ijk] + B3[-5*di + ijk] + 
        210.*(B3[-di + ijk] + B3[di + ijk]) - 
        120.*(B3[-2*di + ijk] + B3[2*di + ijk]) + 
        45.*(B3[-3*di + ijk] + B3[3*di + ijk]) - 
        10.*(B3[-4*di + ijk] + B3[4*di + ijk]) + B3[5*di + ijk]) + 
     oo2dy*(-252.*B3[ijk] + B3[-5*dj + ijk] + 
        210.*(B3[-dj + ijk] + B3[dj + ijk]) - 
        120.*(B3[-2*dj + ijk] + B3[2*dj + ijk]) + 
        45.*(B3[-3*dj + ijk] + B3[3*dj + ijk]) - 
        10.*(B3[-4*dj + ijk] + B3[4*dj + ijk]) + B3[5*dj + ijk]) + 
     oo2dz*(-252.*B3[ijk] + B3[-5*dk + ijk] + 
        210.*(B3[-dk + ijk] + B3[dk + ijk]) - 
        120.*(B3[-2*dk + ijk] + B3[2*dk + ijk]) + 
        45.*(B3[-3*dk + ijk] + B3[3*dk + ijk]) - 
        10.*(B3[-4*dk + ijk] + B3[4*dk + ijk]) + B3[5*dk + ijk]))
;


} 


if (order_dissipation == 12 && boundaryNormore(6)) { 

rA11
=
rA11 - 2.*dissfactor*(oo2dx*(924.*A11[ijk] + A11[-6*di + ijk] - 
        792.*(A11[-di + ijk] + A11[di + ijk]) + 
        495.*(A11[-2*di + ijk] + A11[2*di + ijk]) - 
        220.*(A11[-3*di + ijk] + A11[3*di + ijk]) + 
        66.*(A11[-4*di + ijk] + A11[4*di + ijk]) - 
        12.*(A11[-5*di + ijk] + A11[5*di + ijk]) + A11[6*di + ijk]) + 
     oo2dy*(924.*A11[ijk] + A11[-6*dj + ijk] - 
        792.*(A11[-dj + ijk] + A11[dj + ijk]) + 
        495.*(A11[-2*dj + ijk] + A11[2*dj + ijk]) - 
        220.*(A11[-3*dj + ijk] + A11[3*dj + ijk]) + 
        66.*(A11[-4*dj + ijk] + A11[4*dj + ijk]) - 
        12.*(A11[-5*dj + ijk] + A11[5*dj + ijk]) + A11[6*dj + ijk]) + 
     oo2dz*(924.*A11[ijk] + A11[-6*dk + ijk] - 
        792.*(A11[-dk + ijk] + A11[dk + ijk]) + 
        495.*(A11[-2*dk + ijk] + A11[2*dk + ijk]) - 
        220.*(A11[-3*dk + ijk] + A11[3*dk + ijk]) + 
        66.*(A11[-4*dk + ijk] + A11[4*dk + ijk]) - 
        12.*(A11[-5*dk + ijk] + A11[5*dk + ijk]) + A11[6*dk + ijk]))
;

rA12
=
rA12 - 2.*dissfactor*(oo2dx*(924.*A12[ijk] + A12[-6*di + ijk] - 
        792.*(A12[-di + ijk] + A12[di + ijk]) + 
        495.*(A12[-2*di + ijk] + A12[2*di + ijk]) - 
        220.*(A12[-3*di + ijk] + A12[3*di + ijk]) + 
        66.*(A12[-4*di + ijk] + A12[4*di + ijk]) - 
        12.*(A12[-5*di + ijk] + A12[5*di + ijk]) + A12[6*di + ijk]) + 
     oo2dy*(924.*A12[ijk] + A12[-6*dj + ijk] - 
        792.*(A12[-dj + ijk] + A12[dj + ijk]) + 
        495.*(A12[-2*dj + ijk] + A12[2*dj + ijk]) - 
        220.*(A12[-3*dj + ijk] + A12[3*dj + ijk]) + 
        66.*(A12[-4*dj + ijk] + A12[4*dj + ijk]) - 
        12.*(A12[-5*dj + ijk] + A12[5*dj + ijk]) + A12[6*dj + ijk]) + 
     oo2dz*(924.*A12[ijk] + A12[-6*dk + ijk] - 
        792.*(A12[-dk + ijk] + A12[dk + ijk]) + 
        495.*(A12[-2*dk + ijk] + A12[2*dk + ijk]) - 
        220.*(A12[-3*dk + ijk] + A12[3*dk + ijk]) + 
        66.*(A12[-4*dk + ijk] + A12[4*dk + ijk]) - 
        12.*(A12[-5*dk + ijk] + A12[5*dk + ijk]) + A12[6*dk + ijk]))
;

rA13
=
rA13 - 2.*dissfactor*(oo2dx*(924.*A13[ijk] + A13[-6*di + ijk] - 
        792.*(A13[-di + ijk] + A13[di + ijk]) + 
        495.*(A13[-2*di + ijk] + A13[2*di + ijk]) - 
        220.*(A13[-3*di + ijk] + A13[3*di + ijk]) + 
        66.*(A13[-4*di + ijk] + A13[4*di + ijk]) - 
        12.*(A13[-5*di + ijk] + A13[5*di + ijk]) + A13[6*di + ijk]) + 
     oo2dy*(924.*A13[ijk] + A13[-6*dj + ijk] - 
        792.*(A13[-dj + ijk] + A13[dj + ijk]) + 
        495.*(A13[-2*dj + ijk] + A13[2*dj + ijk]) - 
        220.*(A13[-3*dj + ijk] + A13[3*dj + ijk]) + 
        66.*(A13[-4*dj + ijk] + A13[4*dj + ijk]) - 
        12.*(A13[-5*dj + ijk] + A13[5*dj + ijk]) + A13[6*dj + ijk]) + 
     oo2dz*(924.*A13[ijk] + A13[-6*dk + ijk] - 
        792.*(A13[-dk + ijk] + A13[dk + ijk]) + 
        495.*(A13[-2*dk + ijk] + A13[2*dk + ijk]) - 
        220.*(A13[-3*dk + ijk] + A13[3*dk + ijk]) + 
        66.*(A13[-4*dk + ijk] + A13[4*dk + ijk]) - 
        12.*(A13[-5*dk + ijk] + A13[5*dk + ijk]) + A13[6*dk + ijk]))
;

rA22
=
rA22 - 2.*dissfactor*(oo2dx*(924.*A22[ijk] + A22[-6*di + ijk] - 
        792.*(A22[-di + ijk] + A22[di + ijk]) + 
        495.*(A22[-2*di + ijk] + A22[2*di + ijk]) - 
        220.*(A22[-3*di + ijk] + A22[3*di + ijk]) + 
        66.*(A22[-4*di + ijk] + A22[4*di + ijk]) - 
        12.*(A22[-5*di + ijk] + A22[5*di + ijk]) + A22[6*di + ijk]) + 
     oo2dy*(924.*A22[ijk] + A22[-6*dj + ijk] - 
        792.*(A22[-dj + ijk] + A22[dj + ijk]) + 
        495.*(A22[-2*dj + ijk] + A22[2*dj + ijk]) - 
        220.*(A22[-3*dj + ijk] + A22[3*dj + ijk]) + 
        66.*(A22[-4*dj + ijk] + A22[4*dj + ijk]) - 
        12.*(A22[-5*dj + ijk] + A22[5*dj + ijk]) + A22[6*dj + ijk]) + 
     oo2dz*(924.*A22[ijk] + A22[-6*dk + ijk] - 
        792.*(A22[-dk + ijk] + A22[dk + ijk]) + 
        495.*(A22[-2*dk + ijk] + A22[2*dk + ijk]) - 
        220.*(A22[-3*dk + ijk] + A22[3*dk + ijk]) + 
        66.*(A22[-4*dk + ijk] + A22[4*dk + ijk]) - 
        12.*(A22[-5*dk + ijk] + A22[5*dk + ijk]) + A22[6*dk + ijk]))
;

rA23
=
rA23 - 2.*dissfactor*(oo2dx*(924.*A23[ijk] + A23[-6*di + ijk] - 
        792.*(A23[-di + ijk] + A23[di + ijk]) + 
        495.*(A23[-2*di + ijk] + A23[2*di + ijk]) - 
        220.*(A23[-3*di + ijk] + A23[3*di + ijk]) + 
        66.*(A23[-4*di + ijk] + A23[4*di + ijk]) - 
        12.*(A23[-5*di + ijk] + A23[5*di + ijk]) + A23[6*di + ijk]) + 
     oo2dy*(924.*A23[ijk] + A23[-6*dj + ijk] - 
        792.*(A23[-dj + ijk] + A23[dj + ijk]) + 
        495.*(A23[-2*dj + ijk] + A23[2*dj + ijk]) - 
        220.*(A23[-3*dj + ijk] + A23[3*dj + ijk]) + 
        66.*(A23[-4*dj + ijk] + A23[4*dj + ijk]) - 
        12.*(A23[-5*dj + ijk] + A23[5*dj + ijk]) + A23[6*dj + ijk]) + 
     oo2dz*(924.*A23[ijk] + A23[-6*dk + ijk] - 
        792.*(A23[-dk + ijk] + A23[dk + ijk]) + 
        495.*(A23[-2*dk + ijk] + A23[2*dk + ijk]) - 
        220.*(A23[-3*dk + ijk] + A23[3*dk + ijk]) + 
        66.*(A23[-4*dk + ijk] + A23[4*dk + ijk]) - 
        12.*(A23[-5*dk + ijk] + A23[5*dk + ijk]) + A23[6*dk + ijk]))
;

rA33
=
rA33 - 2.*dissfactor*(oo2dx*(924.*A33[ijk] + A33[-6*di + ijk] - 
        792.*(A33[-di + ijk] + A33[di + ijk]) + 
        495.*(A33[-2*di + ijk] + A33[2*di + ijk]) - 
        220.*(A33[-3*di + ijk] + A33[3*di + ijk]) + 
        66.*(A33[-4*di + ijk] + A33[4*di + ijk]) - 
        12.*(A33[-5*di + ijk] + A33[5*di + ijk]) + A33[6*di + ijk]) + 
     oo2dy*(924.*A33[ijk] + A33[-6*dj + ijk] - 
        792.*(A33[-dj + ijk] + A33[dj + ijk]) + 
        495.*(A33[-2*dj + ijk] + A33[2*dj + ijk]) - 
        220.*(A33[-3*dj + ijk] + A33[3*dj + ijk]) + 
        66.*(A33[-4*dj + ijk] + A33[4*dj + ijk]) - 
        12.*(A33[-5*dj + ijk] + A33[5*dj + ijk]) + A33[6*dj + ijk]) + 
     oo2dz*(924.*A33[ijk] + A33[-6*dk + ijk] - 
        792.*(A33[-dk + ijk] + A33[dk + ijk]) + 
        495.*(A33[-2*dk + ijk] + A33[2*dk + ijk]) - 
        220.*(A33[-3*dk + ijk] + A33[3*dk + ijk]) + 
        66.*(A33[-4*dk + ijk] + A33[4*dk + ijk]) - 
        12.*(A33[-5*dk + ijk] + A33[5*dk + ijk]) + A33[6*dk + ijk]))
;

rG1
=
rG1 - 2.*dissfactor*(oo2dx*(924.*G1[ijk] + G1[-6*di + ijk] - 
        792.*(G1[-di + ijk] + G1[di + ijk]) + 
        495.*(G1[-2*di + ijk] + G1[2*di + ijk]) - 
        220.*(G1[-3*di + ijk] + G1[3*di + ijk]) + 
        66.*(G1[-4*di + ijk] + G1[4*di + ijk]) - 
        12.*(G1[-5*di + ijk] + G1[5*di + ijk]) + G1[6*di + ijk]) + 
     oo2dy*(924.*G1[ijk] + G1[-6*dj + ijk] - 
        792.*(G1[-dj + ijk] + G1[dj + ijk]) + 
        495.*(G1[-2*dj + ijk] + G1[2*dj + ijk]) - 
        220.*(G1[-3*dj + ijk] + G1[3*dj + ijk]) + 
        66.*(G1[-4*dj + ijk] + G1[4*dj + ijk]) - 
        12.*(G1[-5*dj + ijk] + G1[5*dj + ijk]) + G1[6*dj + ijk]) + 
     oo2dz*(924.*G1[ijk] + G1[-6*dk + ijk] - 
        792.*(G1[-dk + ijk] + G1[dk + ijk]) + 
        495.*(G1[-2*dk + ijk] + G1[2*dk + ijk]) - 
        220.*(G1[-3*dk + ijk] + G1[3*dk + ijk]) + 
        66.*(G1[-4*dk + ijk] + G1[4*dk + ijk]) - 
        12.*(G1[-5*dk + ijk] + G1[5*dk + ijk]) + G1[6*dk + ijk]))
;

rG2
=
rG2 - 2.*dissfactor*(oo2dx*(924.*G2[ijk] + G2[-6*di + ijk] - 
        792.*(G2[-di + ijk] + G2[di + ijk]) + 
        495.*(G2[-2*di + ijk] + G2[2*di + ijk]) - 
        220.*(G2[-3*di + ijk] + G2[3*di + ijk]) + 
        66.*(G2[-4*di + ijk] + G2[4*di + ijk]) - 
        12.*(G2[-5*di + ijk] + G2[5*di + ijk]) + G2[6*di + ijk]) + 
     oo2dy*(924.*G2[ijk] + G2[-6*dj + ijk] - 
        792.*(G2[-dj + ijk] + G2[dj + ijk]) + 
        495.*(G2[-2*dj + ijk] + G2[2*dj + ijk]) - 
        220.*(G2[-3*dj + ijk] + G2[3*dj + ijk]) + 
        66.*(G2[-4*dj + ijk] + G2[4*dj + ijk]) - 
        12.*(G2[-5*dj + ijk] + G2[5*dj + ijk]) + G2[6*dj + ijk]) + 
     oo2dz*(924.*G2[ijk] + G2[-6*dk + ijk] - 
        792.*(G2[-dk + ijk] + G2[dk + ijk]) + 
        495.*(G2[-2*dk + ijk] + G2[2*dk + ijk]) - 
        220.*(G2[-3*dk + ijk] + G2[3*dk + ijk]) + 
        66.*(G2[-4*dk + ijk] + G2[4*dk + ijk]) - 
        12.*(G2[-5*dk + ijk] + G2[5*dk + ijk]) + G2[6*dk + ijk]))
;

rG3
=
rG3 - 2.*dissfactor*(oo2dx*(924.*G3[ijk] + G3[-6*di + ijk] - 
        792.*(G3[-di + ijk] + G3[di + ijk]) + 
        495.*(G3[-2*di + ijk] + G3[2*di + ijk]) - 
        220.*(G3[-3*di + ijk] + G3[3*di + ijk]) + 
        66.*(G3[-4*di + ijk] + G3[4*di + ijk]) - 
        12.*(G3[-5*di + ijk] + G3[5*di + ijk]) + G3[6*di + ijk]) + 
     oo2dy*(924.*G3[ijk] + G3[-6*dj + ijk] - 
        792.*(G3[-dj + ijk] + G3[dj + ijk]) + 
        495.*(G3[-2*dj + ijk] + G3[2*dj + ijk]) - 
        220.*(G3[-3*dj + ijk] + G3[3*dj + ijk]) + 
        66.*(G3[-4*dj + ijk] + G3[4*dj + ijk]) - 
        12.*(G3[-5*dj + ijk] + G3[5*dj + ijk]) + G3[6*dj + ijk]) + 
     oo2dz*(924.*G3[ijk] + G3[-6*dk + ijk] - 
        792.*(G3[-dk + ijk] + G3[dk + ijk]) + 
        495.*(G3[-2*dk + ijk] + G3[2*dk + ijk]) - 
        220.*(G3[-3*dk + ijk] + G3[3*dk + ijk]) + 
        66.*(G3[-4*dk + ijk] + G3[4*dk + ijk]) - 
        12.*(G3[-5*dk + ijk] + G3[5*dk + ijk]) + G3[6*dk + ijk]))
;

rK
=
rK - 2.*dissfactor*(oo2dx*(924.*K[ijk] + K[-6*di + ijk] - 
        792.*(K[-di + ijk] + K[di + ijk]) + 
        495.*(K[-2*di + ijk] + K[2*di + ijk]) - 
        220.*(K[-3*di + ijk] + K[3*di + ijk]) + 
        66.*(K[-4*di + ijk] + K[4*di + ijk]) - 
        12.*(K[-5*di + ijk] + K[5*di + ijk]) + K[6*di + ijk]) + 
     oo2dy*(924.*K[ijk] + K[-6*dj + ijk] - 
        792.*(K[-dj + ijk] + K[dj + ijk]) + 
        495.*(K[-2*dj + ijk] + K[2*dj + ijk]) - 
        220.*(K[-3*dj + ijk] + K[3*dj + ijk]) + 
        66.*(K[-4*dj + ijk] + K[4*dj + ijk]) - 
        12.*(K[-5*dj + ijk] + K[5*dj + ijk]) + K[6*dj + ijk]) + 
     oo2dz*(924.*K[ijk] + K[-6*dk + ijk] - 
        792.*(K[-dk + ijk] + K[dk + ijk]) + 
        495.*(K[-2*dk + ijk] + K[2*dk + ijk]) - 
        220.*(K[-3*dk + ijk] + K[3*dk + ijk]) + 
        66.*(K[-4*dk + ijk] + K[4*dk + ijk]) - 
        12.*(K[-5*dk + ijk] + K[5*dk + ijk]) + K[6*dk + ijk]))
;

rB1
=
rB1 - 2.*dissfactor*(oo2dx*(924.*B1[ijk] + B1[-6*di + ijk] - 
        792.*(B1[-di + ijk] + B1[di + ijk]) + 
        495.*(B1[-2*di + ijk] + B1[2*di + ijk]) - 
        220.*(B1[-3*di + ijk] + B1[3*di + ijk]) + 
        66.*(B1[-4*di + ijk] + B1[4*di + ijk]) - 
        12.*(B1[-5*di + ijk] + B1[5*di + ijk]) + B1[6*di + ijk]) + 
     oo2dy*(924.*B1[ijk] + B1[-6*dj + ijk] - 
        792.*(B1[-dj + ijk] + B1[dj + ijk]) + 
        495.*(B1[-2*dj + ijk] + B1[2*dj + ijk]) - 
        220.*(B1[-3*dj + ijk] + B1[3*dj + ijk]) + 
        66.*(B1[-4*dj + ijk] + B1[4*dj + ijk]) - 
        12.*(B1[-5*dj + ijk] + B1[5*dj + ijk]) + B1[6*dj + ijk]) + 
     oo2dz*(924.*B1[ijk] + B1[-6*dk + ijk] - 
        792.*(B1[-dk + ijk] + B1[dk + ijk]) + 
        495.*(B1[-2*dk + ijk] + B1[2*dk + ijk]) - 
        220.*(B1[-3*dk + ijk] + B1[3*dk + ijk]) + 
        66.*(B1[-4*dk + ijk] + B1[4*dk + ijk]) - 
        12.*(B1[-5*dk + ijk] + B1[5*dk + ijk]) + B1[6*dk + ijk]))
;

rB2
=
rB2 - 2.*dissfactor*(oo2dx*(924.*B2[ijk] + B2[-6*di + ijk] - 
        792.*(B2[-di + ijk] + B2[di + ijk]) + 
        495.*(B2[-2*di + ijk] + B2[2*di + ijk]) - 
        220.*(B2[-3*di + ijk] + B2[3*di + ijk]) + 
        66.*(B2[-4*di + ijk] + B2[4*di + ijk]) - 
        12.*(B2[-5*di + ijk] + B2[5*di + ijk]) + B2[6*di + ijk]) + 
     oo2dy*(924.*B2[ijk] + B2[-6*dj + ijk] - 
        792.*(B2[-dj + ijk] + B2[dj + ijk]) + 
        495.*(B2[-2*dj + ijk] + B2[2*dj + ijk]) - 
        220.*(B2[-3*dj + ijk] + B2[3*dj + ijk]) + 
        66.*(B2[-4*dj + ijk] + B2[4*dj + ijk]) - 
        12.*(B2[-5*dj + ijk] + B2[5*dj + ijk]) + B2[6*dj + ijk]) + 
     oo2dz*(924.*B2[ijk] + B2[-6*dk + ijk] - 
        792.*(B2[-dk + ijk] + B2[dk + ijk]) + 
        495.*(B2[-2*dk + ijk] + B2[2*dk + ijk]) - 
        220.*(B2[-3*dk + ijk] + B2[3*dk + ijk]) + 
        66.*(B2[-4*dk + ijk] + B2[4*dk + ijk]) - 
        12.*(B2[-5*dk + ijk] + B2[5*dk + ijk]) + B2[6*dk + ijk]))
;

rB3
=
rB3 - 2.*dissfactor*(oo2dx*(924.*B3[ijk] + B3[-6*di + ijk] - 
        792.*(B3[-di + ijk] + B3[di + ijk]) + 
        495.*(B3[-2*di + ijk] + B3[2*di + ijk]) - 
        220.*(B3[-3*di + ijk] + B3[3*di + ijk]) + 
        66.*(B3[-4*di + ijk] + B3[4*di + ijk]) - 
        12.*(B3[-5*di + ijk] + B3[5*di + ijk]) + B3[6*di + ijk]) + 
     oo2dy*(924.*B3[ijk] + B3[-6*dj + ijk] - 
        792.*(B3[-dj + ijk] + B3[dj + ijk]) + 
        495.*(B3[-2*dj + ijk] + B3[2*dj + ijk]) - 
        220.*(B3[-3*dj + ijk] + B3[3*dj + ijk]) + 
        66.*(B3[-4*dj + ijk] + B3[4*dj + ijk]) - 
        12.*(B3[-5*dj + ijk] + B3[5*dj + ijk]) + B3[6*dj + ijk]) + 
     oo2dz*(924.*B3[ijk] + B3[-6*dk + ijk] - 
        792.*(B3[-dk + ijk] + B3[dk + ijk]) + 
        495.*(B3[-2*dk + ijk] + B3[2*dk + ijk]) - 
        220.*(B3[-3*dk + ijk] + B3[3*dk + ijk]) + 
        66.*(B3[-4*dk + ijk] + B3[4*dk + ijk]) - 
        12.*(B3[-5*dk + ijk] + B3[5*dk + ijk]) + B3[6*dk + ijk]))
;


} 



/* conditional */
if (addlinear) {

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

nK[ijk]
=
c*rK + pK[ijk]
;

nB1[ijk]
=
c*rB1 + pB1[ijk]
;

nB2[ijk]
=
c*rB2 + pB2[ijk]
;

nB3[ijk]
=
c*rB3 + pB3[ijk]
;


} else { /* if (!addlinear) */

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

nK[ijk]
=
rK
;

nB1[ijk]
=
rB1
;

nB2[ijk]
=
rB2
;

nB3[ijk]
=
rB3
;

}
/* if (addlinear) */




} endfor_ijk_openmp; /* loop i, j, k */



bampi_openmp_stop


}  /* function */

/* bssn_boundary.c */
/* nvars = 92, nauxs = 336, n* = 7151,  n/ = 105,  n+ = 10698, n = 17954, O = 0 */
