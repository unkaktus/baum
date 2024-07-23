/* z4_boundary_box.c */
/* Copyright (C) 1998 Bernd Bruegmann, 11.5.2015 */
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



void z4_boundary_box(tVarList *unew, tVarList *upre, double c, tVarList *ucur)
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

double AA11 = 0.;
double AA12 = 0.;
double AA13 = 0.;
double AA21 = 0.;
double AA22 = 0.;
double AA23 = 0.;
double AA31 = 0.;
double AA32 = 0.;
double AA33 = 0.;
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
double betaA1 = 0.;
double betaA2 = 0.;
double betaA3 = 0.;
double betas = 0.;
double cdA111 = 0.;
double cdA112 = 0.;
double cdA113 = 0.;
double cdA122 = 0.;
double cdA123 = 0.;
double cdA133 = 0.;
double cdA211 = 0.;
double cdA212 = 0.;
double cdA213 = 0.;
double cdA222 = 0.;
double cdA223 = 0.;
double cdA233 = 0.;
double cdA311 = 0.;
double cdA312 = 0.;
double cdA313 = 0.;
double cdA322 = 0.;
double cdA323 = 0.;
double cdA333 = 0.;
double cdda11 = 0.;
double cdda12 = 0.;
double cdda13 = 0.;
double cdda22 = 0.;
double cdda23 = 0.;
double cdda33 = 0.;
double cddf11 = 0.;
double cddf12 = 0.;
double cddf13 = 0.;
double cddf22 = 0.;
double cddf23 = 0.;
double cddf33 = 0.;
double chiguard = 0.;
double chiguarded = 0.;
double da1 = 0.;
double dA111 = 0.;
double dA112 = 0.;
double dA113 = 0.;
double dA122 = 0.;
double dA123 = 0.;
double dA133 = 0.;
double da2 = 0.;
double dA211 = 0.;
double dA212 = 0.;
double dA213 = 0.;
double dA222 = 0.;
double dA223 = 0.;
double dA233 = 0.;
double da3 = 0.;
double dA311 = 0.;
double dA312 = 0.;
double dA313 = 0.;
double dA322 = 0.;
double dA323 = 0.;
double dA333 = 0.;
double Dalpha = 0.;
double db11 = 0.;
double db12 = 0.;
double db13 = 0.;
double db21 = 0.;
double db22 = 0.;
double db23 = 0.;
double db31 = 0.;
double db32 = 0.;
double db33 = 0.;
double DbetaA1 = 0.;
double DbetaA2 = 0.;
double DbetaA3 = 0.;
double Dbetas = 0.;
double dchi1 = 0.;
double dchi2 = 0.;
double dchi3 = 0.;
double dda11 = 0.;
double dda12 = 0.;
double dda13 = 0.;
double dda22 = 0.;
double dda23 = 0.;
double dda33 = 0.;
double ddb111 = 0.;
double ddb112 = 0.;
double ddb113 = 0.;
double ddb121 = 0.;
double ddb122 = 0.;
double ddb123 = 0.;
double ddb131 = 0.;
double ddb132 = 0.;
double ddb133 = 0.;
double ddb221 = 0.;
double ddb222 = 0.;
double ddb223 = 0.;
double ddb231 = 0.;
double ddb232 = 0.;
double ddb233 = 0.;
double ddb331 = 0.;
double ddb332 = 0.;
double ddb333 = 0.;
double ddchi11 = 0.;
double ddchi12 = 0.;
double ddchi13 = 0.;
double ddchi22 = 0.;
double ddchi23 = 0.;
double ddchi33 = 0.;
double ddf11 = 0.;
double ddf12 = 0.;
double ddf13 = 0.;
double ddf22 = 0.;
double ddf23 = 0.;
double ddf33 = 0.;
double ddg1111 = 0.;
double ddg1112 = 0.;
double ddg1113 = 0.;
double ddg1122 = 0.;
double ddg1123 = 0.;
double ddg1133 = 0.;
double ddg1211 = 0.;
double ddg1212 = 0.;
double ddg1213 = 0.;
double ddg1222 = 0.;
double ddg1223 = 0.;
double ddg1233 = 0.;
double ddg1311 = 0.;
double ddg1312 = 0.;
double ddg1313 = 0.;
double ddg1322 = 0.;
double ddg1323 = 0.;
double ddg1333 = 0.;
double ddg2211 = 0.;
double ddg2212 = 0.;
double ddg2213 = 0.;
double ddg2222 = 0.;
double ddg2223 = 0.;
double ddg2233 = 0.;
double ddg2311 = 0.;
double ddg2312 = 0.;
double ddg2313 = 0.;
double ddg2322 = 0.;
double ddg2323 = 0.;
double ddg2333 = 0.;
double ddg3311 = 0.;
double ddg3312 = 0.;
double ddg3313 = 0.;
double ddg3322 = 0.;
double ddg3323 = 0.;
double ddg3333 = 0.;
double detginv = 0.;
double df1 = 0.;
double df2 = 0.;
double df3 = 0.;
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
double dGfromgdu11 = 0.;
double dGfromgdu12 = 0.;
double dGfromgdu13 = 0.;
double dGfromgdu21 = 0.;
double dGfromgdu22 = 0.;
double dGfromgdu23 = 0.;
double dGfromgdu31 = 0.;
double dGfromgdu32 = 0.;
double dGfromgdu33 = 0.;
double dginv111 = 0.;
double dginv112 = 0.;
double dginv113 = 0.;
double dginv122 = 0.;
double dginv123 = 0.;
double dginv133 = 0.;
double dginv211 = 0.;
double dginv212 = 0.;
double dginv213 = 0.;
double dginv222 = 0.;
double dginv223 = 0.;
double dginv233 = 0.;
double dginv311 = 0.;
double dginv312 = 0.;
double dginv313 = 0.;
double dginv322 = 0.;
double dginv323 = 0.;
double dginv333 = 0.;
double divbeta = 0.;
double DK = 0.;
double dK1 = 0.;
double dK2 = 0.;
double dK3 = 0.;
double DKhat = 0.;
double dKhat1 = 0.;
double dKhat2 = 0.;
double dKhat3 = 0.;
double DTheta = 0.;
double dTheta1 = 0.;
double dTheta2 = 0.;
double dTheta3 = 0.;
double f = 0.;
double ff = 0.;
double GamA1 = 0.;
double GamA2 = 0.;
double GamA3 = 0.;
double gamma111 = 0.;
double gamma112 = 0.;
double gamma113 = 0.;
double gamma122 = 0.;
double gamma123 = 0.;
double gamma133 = 0.;
double gamma211 = 0.;
double gamma212 = 0.;
double gamma213 = 0.;
double gamma222 = 0.;
double gamma223 = 0.;
double gamma233 = 0.;
double gamma311 = 0.;
double gamma312 = 0.;
double gamma313 = 0.;
double gamma322 = 0.;
double gamma323 = 0.;
double gamma333 = 0.;
double gammado111 = 0.;
double gammado112 = 0.;
double gammado113 = 0.;
double gammado122 = 0.;
double gammado123 = 0.;
double gammado133 = 0.;
double gammado211 = 0.;
double gammado212 = 0.;
double gammado213 = 0.;
double gammado222 = 0.;
double gammado223 = 0.;
double gammado233 = 0.;
double gammado311 = 0.;
double gammado312 = 0.;
double gammado313 = 0.;
double gammado322 = 0.;
double gammado323 = 0.;
double gammado333 = 0.;
double Gams = 0.;
double Gfromg1 = 0.;
double Gfromg2 = 0.;
double Gfromg3 = 0.;
double ginv11 = 0.;
double ginv12 = 0.;
double ginv13 = 0.;
double ginv22 = 0.;
double ginv23 = 0.;
double ginv33 = 0.;
double K = 0.;
double lieA11 = 0.;
double lieA12 = 0.;
double lieA13 = 0.;
double lieA22 = 0.;
double lieA23 = 0.;
double lieA33 = 0.;
double lieg11 = 0.;
double lieg12 = 0.;
double lieg13 = 0.;
double lieg22 = 0.;
double lieg23 = 0.;
double lieg33 = 0.;
double lienK = 0.;
double lienKhat = 0.;
double lienTheta = 0.;
double modshatARG = 0.;
double muL = 0.;
double muStilde = 0.;
double oochipsipower = 0.;
double oomodshat = 0.;
double psim4 = 0.;
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
double R11 = 0.;
double R12 = 0.;
double R13 = 0.;
double R22 = 0.;
double R23 = 0.;
double R33 = 0.;
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
double Rf11 = 0.;
double Rf12 = 0.;
double Rf13 = 0.;
double Rf22 = 0.;
double Rf23 = 0.;
double Rf33 = 0.;
double rG1 = 0.;
double rG2 = 0.;
double rG3 = 0.;
double rGamA1 = 0.;
double rGamA2 = 0.;
double rGamA3 = 0.;
double rGams = 0.;
double Rhat = 0.;
double rKhat = 0.;
double Rphi11 = 0.;
double Rphi12 = 0.;
double Rphi13 = 0.;
double Rphi22 = 0.;
double Rphi23 = 0.;
double Rphi33 = 0.;
double rTheta = 0.;
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
double trcdda = 0.;
double trcddf = 0.;
double vbetaA = 0.;
double vbetas = 0.;
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
    printf("x=%2.5e, y=%2.5e, z=%2.5e\n",xp[ijk],yp[ijk],zp[ijk]);}da1
=
oo2dx*(-alpha[-di + ijk] + alpha[di + ijk])
;

da2
=
oo2dy*(-alpha[-dj + ijk] + alpha[dj + ijk])
;

da3
=
oo2dz*(-alpha[-dk + ijk] + alpha[dk + ijk])
;

dda11
=
oodx2*(-2.*alpha[ijk] + alpha[-di + ijk] + alpha[di + ijk])
;

dda12
=
oo4dxdy*(alpha[-di - dj + ijk] - alpha[di - dj + ijk] - 
    alpha[-di + dj + ijk] + alpha[di + dj + ijk])
;

dda13
=
oo4dxdz*(alpha[-di - dk + ijk] - alpha[di - dk + ijk] - 
    alpha[-di + dk + ijk] + alpha[di + dk + ijk])
;

dda22
=
oody2*(-2.*alpha[ijk] + alpha[-dj + ijk] + alpha[dj + ijk])
;

dda23
=
oo4dydz*(alpha[-dj - dk + ijk] - alpha[dj - dk + ijk] - 
    alpha[-dj + dk + ijk] + alpha[dj + dk + ijk])
;

dda33
=
oodz2*(-2.*alpha[ijk] + alpha[-dk + ijk] + alpha[dk + ijk])
;

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

ddb111
=
oodx2*(-2.*beta1[ijk] + beta1[-di + ijk] + beta1[di + ijk])
;

ddb112
=
oodx2*(-2.*beta2[ijk] + beta2[-di + ijk] + beta2[di + ijk])
;

ddb113
=
oodx2*(-2.*beta3[ijk] + beta3[-di + ijk] + beta3[di + ijk])
;

ddb121
=
oo4dxdy*(beta1[-di - dj + ijk] - beta1[di - dj + ijk] - 
    beta1[-di + dj + ijk] + beta1[di + dj + ijk])
;

ddb122
=
oo4dxdy*(beta2[-di - dj + ijk] - beta2[di - dj + ijk] - 
    beta2[-di + dj + ijk] + beta2[di + dj + ijk])
;

ddb123
=
oo4dxdy*(beta3[-di - dj + ijk] - beta3[di - dj + ijk] - 
    beta3[-di + dj + ijk] + beta3[di + dj + ijk])
;

ddb131
=
oo4dxdz*(beta1[-di - dk + ijk] - beta1[di - dk + ijk] - 
    beta1[-di + dk + ijk] + beta1[di + dk + ijk])
;

ddb132
=
oo4dxdz*(beta2[-di - dk + ijk] - beta2[di - dk + ijk] - 
    beta2[-di + dk + ijk] + beta2[di + dk + ijk])
;

ddb133
=
oo4dxdz*(beta3[-di - dk + ijk] - beta3[di - dk + ijk] - 
    beta3[-di + dk + ijk] + beta3[di + dk + ijk])
;

ddb221
=
oody2*(-2.*beta1[ijk] + beta1[-dj + ijk] + beta1[dj + ijk])
;

ddb222
=
oody2*(-2.*beta2[ijk] + beta2[-dj + ijk] + beta2[dj + ijk])
;

ddb223
=
oody2*(-2.*beta3[ijk] + beta3[-dj + ijk] + beta3[dj + ijk])
;

ddb231
=
oo4dydz*(beta1[-dj - dk + ijk] - beta1[dj - dk + ijk] - 
    beta1[-dj + dk + ijk] + beta1[dj + dk + ijk])
;

ddb232
=
oo4dydz*(beta2[-dj - dk + ijk] - beta2[dj - dk + ijk] - 
    beta2[-dj + dk + ijk] + beta2[dj + dk + ijk])
;

ddb233
=
oo4dydz*(beta3[-dj - dk + ijk] - beta3[dj - dk + ijk] - 
    beta3[-dj + dk + ijk] + beta3[dj + dk + ijk])
;

ddb331
=
oodz2*(-2.*beta1[ijk] + beta1[-dk + ijk] + beta1[dk + ijk])
;

ddb332
=
oodz2*(-2.*beta2[ijk] + beta2[-dk + ijk] + beta2[dk + ijk])
;

ddb333
=
oodz2*(-2.*beta3[ijk] + beta3[-dk + ijk] + beta3[dk + ijk])
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

ddchi11
=
oodx2*(-2.*chi[ijk] + chi[-di + ijk] + chi[di + ijk])
;

ddchi12
=
oo4dxdy*(chi[-di - dj + ijk] - chi[di - dj + ijk] - chi[-di + dj + ijk] + 
    chi[di + dj + ijk])
;

ddchi13
=
oo4dxdz*(chi[-di - dk + ijk] - chi[di - dk + ijk] - chi[-di + dk + ijk] + 
    chi[di + dk + ijk])
;

ddchi22
=
oody2*(-2.*chi[ijk] + chi[-dj + ijk] + chi[dj + ijk])
;

ddchi23
=
oo4dydz*(chi[-dj - dk + ijk] - chi[dj - dk + ijk] - chi[-dj + dk + ijk] + 
    chi[dj + dk + ijk])
;

ddchi33
=
oodz2*(-2.*chi[ijk] + chi[-dk + ijk] + chi[dk + ijk])
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

ddg1111
=
oodx2*(-2.*g11[ijk] + g11[-di + ijk] + g11[di + ijk])
;

ddg1112
=
oodx2*(-2.*g12[ijk] + g12[-di + ijk] + g12[di + ijk])
;

ddg1113
=
oodx2*(-2.*g13[ijk] + g13[-di + ijk] + g13[di + ijk])
;

ddg1122
=
oodx2*(-2.*g22[ijk] + g22[-di + ijk] + g22[di + ijk])
;

ddg1123
=
oodx2*(-2.*g23[ijk] + g23[-di + ijk] + g23[di + ijk])
;

ddg1133
=
oodx2*(-2.*g33[ijk] + g33[-di + ijk] + g33[di + ijk])
;

ddg1211
=
oo4dxdy*(g11[-di - dj + ijk] - g11[di - dj + ijk] - g11[-di + dj + ijk] + 
    g11[di + dj + ijk])
;

ddg1212
=
oo4dxdy*(g12[-di - dj + ijk] - g12[di - dj + ijk] - g12[-di + dj + ijk] + 
    g12[di + dj + ijk])
;

ddg1213
=
oo4dxdy*(g13[-di - dj + ijk] - g13[di - dj + ijk] - g13[-di + dj + ijk] + 
    g13[di + dj + ijk])
;

ddg1222
=
oo4dxdy*(g22[-di - dj + ijk] - g22[di - dj + ijk] - g22[-di + dj + ijk] + 
    g22[di + dj + ijk])
;

ddg1223
=
oo4dxdy*(g23[-di - dj + ijk] - g23[di - dj + ijk] - g23[-di + dj + ijk] + 
    g23[di + dj + ijk])
;

ddg1233
=
oo4dxdy*(g33[-di - dj + ijk] - g33[di - dj + ijk] - g33[-di + dj + ijk] + 
    g33[di + dj + ijk])
;

ddg1311
=
oo4dxdz*(g11[-di - dk + ijk] - g11[di - dk + ijk] - g11[-di + dk + ijk] + 
    g11[di + dk + ijk])
;

ddg1312
=
oo4dxdz*(g12[-di - dk + ijk] - g12[di - dk + ijk] - g12[-di + dk + ijk] + 
    g12[di + dk + ijk])
;

ddg1313
=
oo4dxdz*(g13[-di - dk + ijk] - g13[di - dk + ijk] - g13[-di + dk + ijk] + 
    g13[di + dk + ijk])
;

ddg1322
=
oo4dxdz*(g22[-di - dk + ijk] - g22[di - dk + ijk] - g22[-di + dk + ijk] + 
    g22[di + dk + ijk])
;

ddg1323
=
oo4dxdz*(g23[-di - dk + ijk] - g23[di - dk + ijk] - g23[-di + dk + ijk] + 
    g23[di + dk + ijk])
;

ddg1333
=
oo4dxdz*(g33[-di - dk + ijk] - g33[di - dk + ijk] - g33[-di + dk + ijk] + 
    g33[di + dk + ijk])
;

ddg2211
=
oody2*(-2.*g11[ijk] + g11[-dj + ijk] + g11[dj + ijk])
;

ddg2212
=
oody2*(-2.*g12[ijk] + g12[-dj + ijk] + g12[dj + ijk])
;

ddg2213
=
oody2*(-2.*g13[ijk] + g13[-dj + ijk] + g13[dj + ijk])
;

ddg2222
=
oody2*(-2.*g22[ijk] + g22[-dj + ijk] + g22[dj + ijk])
;

ddg2223
=
oody2*(-2.*g23[ijk] + g23[-dj + ijk] + g23[dj + ijk])
;

ddg2233
=
oody2*(-2.*g33[ijk] + g33[-dj + ijk] + g33[dj + ijk])
;

ddg2311
=
oo4dydz*(g11[-dj - dk + ijk] - g11[dj - dk + ijk] - g11[-dj + dk + ijk] + 
    g11[dj + dk + ijk])
;

ddg2312
=
oo4dydz*(g12[-dj - dk + ijk] - g12[dj - dk + ijk] - g12[-dj + dk + ijk] + 
    g12[dj + dk + ijk])
;

ddg2313
=
oo4dydz*(g13[-dj - dk + ijk] - g13[dj - dk + ijk] - g13[-dj + dk + ijk] + 
    g13[dj + dk + ijk])
;

ddg2322
=
oo4dydz*(g22[-dj - dk + ijk] - g22[dj - dk + ijk] - g22[-dj + dk + ijk] + 
    g22[dj + dk + ijk])
;

ddg2323
=
oo4dydz*(g23[-dj - dk + ijk] - g23[dj - dk + ijk] - g23[-dj + dk + ijk] + 
    g23[dj + dk + ijk])
;

ddg2333
=
oo4dydz*(g33[-dj - dk + ijk] - g33[dj - dk + ijk] - g33[-dj + dk + ijk] + 
    g33[dj + dk + ijk])
;

ddg3311
=
oodz2*(-2.*g11[ijk] + g11[-dk + ijk] + g11[dk + ijk])
;

ddg3312
=
oodz2*(-2.*g12[ijk] + g12[-dk + ijk] + g12[dk + ijk])
;

ddg3313
=
oodz2*(-2.*g13[ijk] + g13[-dk + ijk] + g13[dk + ijk])
;

ddg3322
=
oodz2*(-2.*g22[ijk] + g22[-dk + ijk] + g22[dk + ijk])
;

ddg3323
=
oodz2*(-2.*g23[ijk] + g23[-dk + ijk] + g23[dk + ijk])
;

ddg3333
=
oodz2*(-2.*g33[ijk] + g33[-dk + ijk] + g33[dk + ijk])
;

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


x=xp[ijk]; y=yp[ijk]; z=zp[ijk]; 


r=sqrt(x*x + y*y + z*z); 


if (!(give_cube_normal (i,j,k, imax,jmax,kmax,N_points_away, &shat1,&shat2,&shat3))) 
     errorexit("not at boundary");
chiguarded
=
chi[ijk]
;

chiguard
=
chiDivFloor
;


if (chiguarded<chiguard) {                    
      printf("chi is wrong (%e)\n",chi[ijk]);
      chiguarded = chiguard;
    }detginv
=
1/(2.*g12[ijk]*g13[ijk]*g23[ijk] + g11[ijk]*g22[ijk]*g33[ijk] - 
    g33[ijk]*pow2(g12[ijk]) - g22[ijk]*pow2(g13[ijk]) - 
    g11[ijk]*pow2(g23[ijk]))
;


if (detginv<0.00001) {                           
      printf("detginv is wrong (%e)\n",detginv);
      detginv = 1.0;
    }ginv11
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
chiguarded*ginv11
;

ADMginv12
=
chiguarded*ginv12
;

ADMginv13
=
chiguarded*ginv13
;

ADMginv22
=
chiguarded*ginv22
;

ADMginv23
=
chiguarded*ginv23
;

ADMginv33
=
chiguarded*ginv33
;

modshatARG
=
2.*(ADMginv23*shat2*shat3 + shat1*(ADMginv12*shat2 + ADMginv13*shat3)) + 
  ADMginv11*pow2(shat1) + ADMginv22*pow2(shat2) + ADMginv33*pow2(shat3)
;


if (modshatARG<0.00001) {                                   
      printf("modshat is wrong (%e)\n",modshatARG);
      modshatARG = pow2(shat1) + pow2(shat2) + pow2(shat3);
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
g11[ijk]/chiguarded - pow2(sdown1)
;

qdd12
=
-(sdown1*sdown2) + g12[ijk]/chiguarded
;

qdd13
=
-(sdown1*sdown3) + g13[ijk]/chiguarded
;

qdd22
=
g22[ijk]/chiguarded - pow2(sdown2)
;

qdd23
=
-(sdown2*sdown3) + g23[ijk]/chiguarded
;

qdd33
=
g33[ijk]/chiguarded - pow2(sdown3)
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

muStilde
=
1/chiguarded
;

vbetas
=
2.*sqrt(0.33333333333333333333*muStilde)
;

vbetaA
=
sqrt(muStilde)
;

K
=
Khat[ijk] + 2.*Theta[ijk]
;

dK1
=
dKhat1 + 2.*dTheta1
;

dK2
=
dKhat2 + 2.*dTheta2
;

dK3
=
dKhat3 + 2.*dTheta3
;

dginv111
=
-2.*(dg123*ginv12*ginv13 + ginv11*(dg112*ginv12 + dg113*ginv13)) - 
  dg111*pow2(ginv11) - dg122*pow2(ginv12) - dg133*pow2(ginv13)
;

dginv112
=
-(ginv11*(dg111*ginv12 + dg112*ginv22 + dg113*ginv23)) - 
  ginv12*(dg113*ginv13 + dg122*ginv22 + dg123*ginv23) - 
  ginv13*(dg123*ginv22 + dg133*ginv23) - dg112*pow2(ginv12)
;

dginv113
=
-(ginv11*(dg111*ginv13 + dg112*ginv23 + dg113*ginv33)) - 
  ginv12*(dg112*ginv13 + dg122*ginv23 + dg123*ginv33) - 
  ginv13*(dg123*ginv23 + dg133*ginv33) - dg113*pow2(ginv13)
;

dginv122
=
-2.*(dg123*ginv22*ginv23 + ginv12*(dg112*ginv22 + dg113*ginv23)) - 
  dg111*pow2(ginv12) - dg122*pow2(ginv22) - dg133*pow2(ginv23)
;

dginv123
=
-(ginv13*(dg112*ginv22 + dg113*ginv23)) - dg133*ginv23*ginv33 - 
  ginv12*(dg111*ginv13 + dg112*ginv23 + dg113*ginv33) - 
  ginv22*(dg122*ginv23 + dg123*ginv33) - dg123*pow2(ginv23)
;

dginv133
=
-2.*(dg123*ginv23*ginv33 + ginv13*(dg112*ginv23 + dg113*ginv33)) - 
  dg111*pow2(ginv13) - dg122*pow2(ginv23) - dg133*pow2(ginv33)
;

dginv211
=
-2.*(dg223*ginv12*ginv13 + ginv11*(dg212*ginv12 + dg213*ginv13)) - 
  dg211*pow2(ginv11) - dg222*pow2(ginv12) - dg233*pow2(ginv13)
;

dginv212
=
-(ginv11*(dg211*ginv12 + dg212*ginv22 + dg213*ginv23)) - 
  ginv12*(dg213*ginv13 + dg222*ginv22 + dg223*ginv23) - 
  ginv13*(dg223*ginv22 + dg233*ginv23) - dg212*pow2(ginv12)
;

dginv213
=
-(ginv11*(dg211*ginv13 + dg212*ginv23 + dg213*ginv33)) - 
  ginv12*(dg212*ginv13 + dg222*ginv23 + dg223*ginv33) - 
  ginv13*(dg223*ginv23 + dg233*ginv33) - dg213*pow2(ginv13)
;

dginv222
=
-2.*(dg223*ginv22*ginv23 + ginv12*(dg212*ginv22 + dg213*ginv23)) - 
  dg211*pow2(ginv12) - dg222*pow2(ginv22) - dg233*pow2(ginv23)
;

dginv223
=
-(ginv13*(dg212*ginv22 + dg213*ginv23)) - dg233*ginv23*ginv33 - 
  ginv12*(dg211*ginv13 + dg212*ginv23 + dg213*ginv33) - 
  ginv22*(dg222*ginv23 + dg223*ginv33) - dg223*pow2(ginv23)
;

dginv233
=
-2.*(dg223*ginv23*ginv33 + ginv13*(dg212*ginv23 + dg213*ginv33)) - 
  dg211*pow2(ginv13) - dg222*pow2(ginv23) - dg233*pow2(ginv33)
;

dginv311
=
-2.*(dg323*ginv12*ginv13 + ginv11*(dg312*ginv12 + dg313*ginv13)) - 
  dg311*pow2(ginv11) - dg322*pow2(ginv12) - dg333*pow2(ginv13)
;

dginv312
=
-(ginv11*(dg311*ginv12 + dg312*ginv22 + dg313*ginv23)) - 
  ginv12*(dg313*ginv13 + dg322*ginv22 + dg323*ginv23) - 
  ginv13*(dg323*ginv22 + dg333*ginv23) - dg312*pow2(ginv12)
;

dginv313
=
-(ginv11*(dg311*ginv13 + dg312*ginv23 + dg313*ginv33)) - 
  ginv12*(dg312*ginv13 + dg322*ginv23 + dg323*ginv33) - 
  ginv13*(dg323*ginv23 + dg333*ginv33) - dg313*pow2(ginv13)
;

dginv322
=
-2.*(dg323*ginv22*ginv23 + ginv12*(dg312*ginv22 + dg313*ginv23)) - 
  dg311*pow2(ginv12) - dg322*pow2(ginv22) - dg333*pow2(ginv23)
;

dginv323
=
-(ginv13*(dg312*ginv22 + dg313*ginv23)) - dg333*ginv23*ginv33 - 
  ginv12*(dg311*ginv13 + dg312*ginv23 + dg313*ginv33) - 
  ginv22*(dg322*ginv23 + dg323*ginv33) - dg323*pow2(ginv23)
;

dginv333
=
-2.*(dg323*ginv23*ginv33 + ginv13*(dg312*ginv23 + dg313*ginv33)) - 
  dg311*pow2(ginv13) - dg322*pow2(ginv23) - dg333*pow2(ginv33)
;

gammado111
=
0.5*dg111
;

gammado112
=
0.5*dg211
;

gammado113
=
0.5*dg311
;

gammado122
=
-0.5*dg122 + dg212
;

gammado123
=
0.5*(-dg123 + dg213 + dg312)
;

gammado133
=
-0.5*dg133 + dg313
;

gammado211
=
dg112 - 0.5*dg211
;

gammado212
=
0.5*dg122
;

gammado213
=
0.5*(dg123 - dg213 + dg312)
;

gammado222
=
0.5*dg222
;

gammado223
=
0.5*dg322
;

gammado233
=
-0.5*dg233 + dg323
;

gammado311
=
dg113 - 0.5*dg311
;

gammado312
=
0.5*(dg123 + dg213 - dg312)
;

gammado313
=
0.5*dg133
;

gammado322
=
dg223 - 0.5*dg322
;

gammado323
=
0.5*dg233
;

gammado333
=
0.5*dg333
;

gamma111
=
gammado111*ginv11 + gammado211*ginv12 + gammado311*ginv13
;

gamma112
=
gammado112*ginv11 + gammado212*ginv12 + gammado312*ginv13
;

gamma113
=
gammado113*ginv11 + gammado213*ginv12 + gammado313*ginv13
;

gamma122
=
gammado122*ginv11 + gammado222*ginv12 + gammado322*ginv13
;

gamma123
=
gammado123*ginv11 + gammado223*ginv12 + gammado323*ginv13
;

gamma133
=
gammado133*ginv11 + gammado233*ginv12 + gammado333*ginv13
;

gamma211
=
gammado111*ginv12 + gammado211*ginv22 + gammado311*ginv23
;

gamma212
=
gammado112*ginv12 + gammado212*ginv22 + gammado312*ginv23
;

gamma213
=
gammado113*ginv12 + gammado213*ginv22 + gammado313*ginv23
;

gamma222
=
gammado122*ginv12 + gammado222*ginv22 + gammado322*ginv23
;

gamma223
=
gammado123*ginv12 + gammado223*ginv22 + gammado323*ginv23
;

gamma233
=
gammado133*ginv12 + gammado233*ginv22 + gammado333*ginv23
;

gamma311
=
gammado111*ginv13 + gammado211*ginv23 + gammado311*ginv33
;

gamma312
=
gammado112*ginv13 + gammado212*ginv23 + gammado312*ginv33
;

gamma313
=
gammado113*ginv13 + gammado213*ginv23 + gammado313*ginv33
;

gamma322
=
gammado122*ginv13 + gammado222*ginv23 + gammado322*ginv33
;

gamma323
=
gammado123*ginv13 + gammado223*ginv23 + gammado323*ginv33
;

gamma333
=
gammado133*ginv13 + gammado233*ginv23 + gammado333*ginv33
;

Gfromg1
=
gamma111*ginv11 + gamma122*ginv22 + 
  2.*(gamma112*ginv12 + gamma113*ginv13 + gamma123*ginv23) + gamma133*ginv33
;

Gfromg2
=
gamma211*ginv11 + gamma222*ginv22 + 
  2.*(gamma212*ginv12 + gamma213*ginv13 + gamma223*ginv23) + gamma233*ginv33
;

Gfromg3
=
gamma311*ginv11 + gamma322*ginv22 + 
  2.*(gamma312*ginv12 + gamma313*ginv13 + gamma323*ginv23) + gamma333*ginv33
;

dGfromgdu11
=
(ddg1111 - dg111*((8.*dg112 + 2.*dg211)*ginv12 + 
        (8.*dg113 + 2.*dg311)*ginv13) - 
     (dg113*(4.*dg112 + dg211) + dg112*dg311 + dg111*(dg213 + dg312))*
      ginv23 - ginv22*(dg112*dg211 + dg111*dg212 + 2.*pow2(dg112)) - 
     ginv33*(dg113*dg311 + dg111*dg313 + 2.*pow2(dg113)))*pow2(ginv11) + 
  (ddg1122 + ddg1212 - (dg123*(8.*dg112 + 2.*dg211) + 
        dg113*(4.*dg122 + 2.*dg212) + dg122*dg311 + 
        2.*(dg111*dg223 + dg112*(dg213 + dg312)) + dg111*dg322)*ginv13 - 
     (dg123*(4.*dg122 + 2.*dg212) + 
        2.*(dg113*dg222 + dg122*(dg213 + dg312) + dg112*(dg223 + dg322)))*
      ginv23 - ginv22*(3.*(dg122*dg212 + dg112*dg222) + 2.*pow2(dg122)) - 
     ginv33*(dg123*(dg213 + dg312) + dg122*dg313 + dg113*(dg223 + dg322) + 
        dg112*dg323 + 2.*pow2(dg123)))*pow2(ginv12) + 
  (ddg1133 + ddg1313 - (dg133*(4.*dg123 + 2.*(dg213 + dg312)) + 
        2.*(dg123*dg313 + dg113*(dg233 + dg323) + dg112*dg333))*ginv23 - 
     ginv22*(dg133*dg212 + dg113*dg223 + dg123*(dg213 + dg312) + 
        dg112*(dg233 + dg323) + 2.*pow2(dg123)) - 
     ginv33*(3.*(dg133*dg313 + dg113*dg333) + 2.*pow2(dg133)))*pow2(ginv13) \
+ ginv13*(ddg1333*ginv33 + ginv22*
      (ddg1223 - (dg133*dg222 + dg123*(4.*dg223 + dg322) + 
           dg122*(dg233 + dg323))*ginv23 - 
        (dg133*dg223 + dg123*(dg233 + 2.*dg323))*ginv33) + 
     ginv23*(ddg1233 + ddg1323 - 
        (dg133*(2.*dg233 + 3.*dg323) + 3.*dg123*dg333)*ginv33) - 
     (dg123*dg222 + dg122*dg223)*pow2(ginv22) - 
     (dg133*dg322 + 2.*(dg133*dg223 + dg123*(dg233 + dg323)) + 
        dg122*dg333)*pow2(ginv23) - 2.*dg133*dg333*pow2(ginv33)) + 
  ginv11*(ddg1313*ginv33 + ginv12*
      (2.*ddg1112 + ddg1211 - (dg113*(12.*dg112 + 3.*dg211) + 
           3.*dg112*dg311 + dg111*(8.*dg123 + 3.*(dg213 + dg312)))*ginv13 \
- (dg122*(4.*dg112 + dg211) + 6.*dg112*dg212 + dg111*dg222)*ginv22 - 
        (dg123*dg211 + dg122*dg311 + 
           4.*(dg113*(dg122 + dg212) + dg112*(dg123 + dg213 + dg312)) + 
           dg111*(dg223 + dg322))*ginv23 - 
        (dg123*dg311 + dg113*(4.*dg123 + 2.*(dg213 + dg312)) + 
           2.*dg112*dg313 + dg111*dg323)*ginv33) + 
     ginv22*(ddg1212 - (dg113*dg222 + 2.*(dg123*dg212 + dg112*dg223) + 
           dg122*(dg213 + dg312) + dg112*dg322)*ginv23 - 
        (dg113*dg223 + dg123*(dg213 + dg312) + dg112*dg323)*ginv33) + 
     ginv13*(2.*ddg1113 + ddg1311 - 
        (dg123*(4.*dg112 + dg211) + dg111*dg223 + 
           2.*(dg113*dg212 + dg112*(dg213 + dg312)))*ginv22 - 
        (dg133*dg211 + dg123*dg311 + 
           4.*(dg113*(dg123 + dg213 + dg312) + dg112*(dg133 + dg313)) + 
           dg111*(dg233 + dg323))*ginv23 - 
        (dg133*(4.*dg113 + dg311) + 6.*dg113*dg313 + dg111*dg333)*ginv33) + 
     ginv23*(ddg1213 + ddg1312 - 
        (dg133*(dg213 + dg312) + 2.*dg123*dg313 + 
           dg113*(dg233 + 2.*dg323) + dg112*dg333)*ginv33) - 
     (3.*dg112*dg211 + dg111*(4.*dg122 + 3.*dg212) + 6.*pow2(dg112))*
      pow2(ginv12) - (3.*dg113*dg311 + dg111*(4.*dg133 + 3.*dg313) + 
        6.*pow2(dg113))*pow2(ginv13) - 
     (dg122*dg212 + dg112*dg222)*pow2(ginv22) - 
     (dg133*dg212 + dg123*(dg213 + dg312) + dg122*dg313 + 
        dg113*(dg223 + dg322) + dg112*(dg233 + dg323))*pow2(ginv23) - 
     (dg133*dg313 + dg113*dg333)*pow2(ginv33)) + 
  ginv12*(ddg1323*ginv33 + ginv22*
      (ddg1222 - (3.*(dg123*dg222 + dg122*dg223) + 2.*dg122*dg322)*
         ginv23 - (dg123*(2.*dg223 + dg322) + dg122*dg323)*ginv33) + 
     ginv23*(ddg1223 + ddg1322 - 
        (dg133*(dg223 + dg322) + dg123*(dg233 + 4.*dg323) + dg122*dg333)*
         ginv33) + ginv13*(2.*ddg1123 + ddg1213 + ddg1312 - 
        (dg113*dg222 + 4.*(dg123*(dg122 + dg212) + dg112*dg223) + 
           dg122*(dg213 + dg312) + dg112*dg322)*ginv22 - 
        (dg133*(4.*dg123 + dg213 + dg312) + 4.*dg123*dg313 + 
           dg113*(dg233 + 4.*dg323) + dg112*dg333)*ginv33 - 
        ginv23*(2.*(dg133*dg212 + dg112*dg233 + dg122*dg313 + 
              dg113*dg322) + 4.*
            (dg122*dg133 + dg113*dg223 + dg123*(dg213 + dg312) + 
              dg112*dg323 + pow2(dg123)))) - 
     (dg133*(4.*dg112 + dg211) + dg113*(8.*dg123 + 2.*(dg213 + dg312)) + 
        2.*(dg123*dg311 + dg112*dg313) + dg111*(dg233 + 2.*dg323))*
      pow2(ginv13) - 2.*dg122*dg222*pow2(ginv22) - 
     (dg133*dg222 + 2.*dg123*(dg223 + dg322) + dg122*(dg233 + 2.*dg323))*
      pow2(ginv23) - (dg133*dg323 + dg123*dg333)*pow2(ginv33)) - 
  2.*pow2(dg111)*pow3(ginv11) - 
  (dg122*(4.*dg112 + dg211) + 2.*dg112*dg212 + dg111*dg222)*pow3(ginv12) - 
  (dg133*(4.*dg113 + dg311) + 2.*dg113*dg313 + dg111*dg333)*pow3(ginv13)
;

dGfromgdu12
=
(ddg1112 + ddg1211 - (4.*(dg112*dg113 + dg111*dg123) + 
        2.*(dg113*dg211 + dg112*dg311 + dg111*(dg213 + dg312)))*ginv13 - 
     (dg122*(6.*dg112 + 2.*dg211) + 6.*dg112*dg212 + 2.*dg111*dg222)*
      ginv22 - (4.*(dg113*(dg122 + dg212) + dg112*(dg123 + dg213)) + 
        dg122*dg311 + 2.*(dg123*dg211 + dg111*dg223 + dg112*dg312) + 
        dg111*dg322)*ginv23 - (dg123*dg311 + 
        dg113*(2.*(dg123 + dg213) + dg312) + dg112*dg313 + dg111*dg323)*
      ginv33)*pow2(ginv12) - ((2.*(dg113*dg123 + dg112*dg133) + 
        dg123*dg311 + dg113*dg312 + dg112*dg313 + dg111*dg323)*ginv22 + 
     (dg133*(4.*dg113 + dg311) + 2.*dg113*dg313 + dg111*dg333)*ginv23)*
   pow2(ginv13) + (ddg1222 - (4.*(dg123*dg222 + dg122*dg223) + 
        2.*dg122*dg322)*ginv23 - 
     (dg123*(2.*dg223 + dg322) + dg122*dg323)*ginv33)*pow2(ginv22) + 
  (ddg1233 + ddg1323 - (dg133*(2.*dg233 + 3.*dg323) + 3.*dg123*dg333)*
      ginv33)*pow2(ginv23) + ginv11*
   (ginv23*(ddg1113 - 2.*dg113*(dg133 + dg313)*ginv33) + 
     ginv22*(ddg1112 - (dg112*(4.*dg123 + 2.*dg213) + 
           2.*(dg113*(dg122 + dg212) + dg112*dg312))*ginv23 - 
        (dg113*(2.*dg123 + dg312) + dg112*dg313)*ginv33) + 
     ginv12*(ddg1111 - dg111*(6.*dg113 + 2.*dg311)*ginv13 - 
        (dg113*(8.*dg112 + 2.*dg211) + dg112*dg311 + 
           dg111*(2.*(dg123 + dg213) + dg312))*ginv23 - 
        ginv22*(2.*(dg112*dg211 + dg111*(dg122 + dg212)) + 
           6.*pow2(dg112)) - ginv33*
         (dg113*dg311 + dg111*dg313 + 2.*pow2(dg113))) - 
     ginv13*((dg112*(4.*dg113 + dg311) + dg111*(2.*dg123 + dg312))*
         ginv22 + ginv23*(dg113*dg311 + dg111*(2.*dg133 + dg313) + 
           4.*pow2(dg113))) - dg111*(6.*dg112 + 2.*dg211)*pow2(ginv12) - 
     2.*dg112*(dg122 + dg212)*pow2(ginv22) - 
     (2.*(dg112*dg133 + dg113*(dg123 + dg213)) + dg113*dg312 + dg112*dg313)*
      pow2(ginv23)) + ginv13*(ginv22*
      (ddg1123 + ddg1312 - (dg133*(2.*dg123 + dg312) + 
           2.*(dg123*dg313 + dg113*dg323) + dg112*dg333)*ginv33 - 
        ginv23*(2.*(dg133*(dg122 + dg212) + dg123*dg213 + dg113*dg223 + 
              dg112*dg233) + dg122*dg313 + dg113*dg322 + 
           4.*(dg123*dg312 + dg112*dg323 + pow2(dg123)))) + 
     ginv23*(ddg1133 + ddg1313 - 
        ginv33*(3.*(dg133*dg313 + dg113*dg333) + 2.*pow2(dg133))) - 
     (2.*(dg123*(dg122 + dg212) + dg112*dg223) + dg122*dg312 + 
        dg112*dg322)*pow2(ginv22) - 
     (dg133*(4.*dg123 + 2.*(dg213 + dg312)) + 
        2.*(dg123*dg313 + dg113*(dg233 + dg323) + dg112*dg333))*pow2(ginv23)\
) + ginv23*(ddg1333*ginv33 - 2.*dg133*dg333*pow2(ginv33)) + 
  ginv12*(ddg1313*ginv33 + ginv13*
      (ddg1113 + ddg1311 - (2.*
            (dg123*dg211 + dg113*(dg122 + dg212) + dg111*dg223) + 
           dg122*dg311 + dg112*(8.*dg123 + 2.*dg213 + 4.*dg312) + 
           dg111*dg322)*ginv22 - 
        (dg133*(4.*dg112 + 2.*dg211) + 
           dg113*(8.*dg123 + 4.*(dg213 + dg312)) + 4.*dg112*dg313 + 
           2.*(dg123*dg311 + dg111*(dg233 + dg323)))*ginv23 - 
        (dg133*(2.*dg113 + dg311) + 4.*dg113*dg313 + dg111*dg333)*ginv33) + 
     ginv23*(ddg1123 + 2.*ddg1213 + ddg1312 - 
        (2.*(dg133*(dg123 + dg213) + dg113*dg233) + dg133*dg312 + 
           4.*(dg123*dg313 + dg113*dg323) + dg112*dg333)*ginv33) + 
     ginv22*(ddg1122 + 2.*ddg1212 - 
        (4.*(dg122*dg213 + dg113*dg222) + 
           6.*(dg123*(dg122 + dg212) + dg112*dg223) + 
           3.*(dg122*dg312 + dg112*dg322))*ginv23 - 
        ginv33*(dg122*dg313 + dg113*dg322 + 
           2.*(dg113*dg223 + dg123*(dg213 + dg312) + dg112*dg323 + 
              pow2(dg123)))) - 
     2.*(dg113*dg311 + dg111*(dg133 + dg313) + pow2(dg113))*pow2(ginv13) - 
     (4.*(dg122*dg212 + dg112*dg222) + 2.*pow2(dg122))*pow2(ginv22) - 
     (4.*(dg123*dg213 + dg113*dg223) + 
        2.*(dg133*(dg122 + dg212) + dg123*dg312 + dg122*dg313 + 
           dg113*dg322 + dg112*(dg233 + dg323) + pow2(dg123)))*pow2(ginv23) \
- (dg133*dg313 + dg113*dg333)*pow2(ginv33)) + 
  ginv22*(ddg1323*ginv33 + ginv23*
      (2.*ddg1223 + ddg1322 - (2.*(dg133*dg223 + dg123*dg233) + 
           dg133*dg322 + 6.*dg123*dg323 + dg122*dg333)*ginv33) - 
     (2.*(dg133*dg222 + dg122*dg233) + dg123*(6.*dg223 + 3.*dg322) + 
        3.*dg122*dg323)*pow2(ginv23) - 
     (dg133*dg323 + dg123*dg333)*pow2(ginv33)) - 
  2.*((dg111*(dg112*ginv22 + dg113*ginv23) + ginv12*pow2(dg111))*
      pow2(ginv11) + (dg112*dg211 + dg111*(dg122 + dg212) + pow2(dg112))*
      pow3(ginv12) + dg122*dg222*pow3(ginv22)) - 
  (dg133*dg322 + 2.*(dg133*dg223 + dg123*(dg233 + dg323)) + dg122*dg333)*
   pow3(ginv23)
;

dGfromgdu13
=
-(((dg122*(4.*dg112 + dg211) + 2.*dg112*dg212 + dg111*dg222)*ginv23 + 
       (2.*(dg113*dg122 + dg112*dg123) + dg123*dg211 + dg113*dg212 + 
          dg112*dg213 + dg111*dg223)*ginv33 + 
       2.*ginv13*(dg112*dg211 + dg111*(dg122 + dg212) + pow2(dg112)))*
     pow2(ginv12)) + (ddg1113 + ddg1311 - 
     (dg123*(2.*dg112 + dg211) + dg113*dg212 + dg111*dg223 + 
        dg112*(dg213 + 2.*dg312))*ginv22 - 
     (dg133*dg211 + 2.*(dg113*dg213 + dg123*dg311) + 
        4.*(dg113*(dg123 + dg312) + dg112*(dg133 + dg313)) + 
        dg111*(dg233 + 2.*dg323))*ginv23 - 
     (dg133*(6.*dg113 + 2.*dg311) + 6.*dg113*dg313 + 2.*dg111*dg333)*ginv33\
)*pow2(ginv13) - (2.*dg122*dg222*ginv23 + 
     (dg123*dg222 + dg122*dg223)*ginv33)*pow2(ginv22) + 
  (ddg1223 + ddg1322 - (3.*(dg133*dg223 + dg123*dg233) + 6.*dg123*dg323 + 
        2.*(dg133*dg322 + dg122*dg333))*ginv33)*pow2(ginv23) + 
  ddg1333*pow2(ginv33) + ginv11*
   (ddg1113*ginv33 - ginv22*(2.*dg112*(dg122 + dg212)*ginv23 + 
        (dg113*dg212 + dg112*(2.*dg123 + dg213))*ginv33) + 
     ginv23*(ddg1112 - (dg113*(4.*dg123 + 2.*dg213) + 
           2.*(dg113*dg312 + dg112*(dg133 + dg313)))*ginv33) - 
     ginv12*(dg111*(6.*dg112 + 2.*dg211)*ginv13 + 
        (dg113*(4.*dg112 + dg211) + dg111*(2.*dg123 + dg213))*ginv33 + 
        ginv23*(dg112*dg211 + dg111*(2.*dg122 + dg212) + 4.*pow2(dg112))) + 
     ginv13*(ddg1111 - (dg113*(8.*dg112 + dg211) + 2.*dg112*dg311 + 
           dg111*(dg213 + 2.*(dg123 + dg312)))*ginv23 - 
        ginv22*(dg112*dg211 + dg111*dg212 + 2.*pow2(dg112)) - 
        ginv33*(2.*(dg113*dg311 + dg111*(dg133 + dg313)) + 6.*pow2(dg113))) \
- dg111*(6.*dg113 + 2.*dg311)*pow2(ginv13) - 
     (dg113*dg212 + dg112*dg213 + 
        2.*(dg113*dg122 + dg112*(dg123 + dg312)))*pow2(ginv23) - 
     2.*dg113*(dg133 + dg313)*pow2(ginv33)) + 
  ginv12*((ddg1123 + ddg1213)*ginv33 + 
     ginv13*(ddg1112 + ddg1211 - 
        (dg122*(2.*dg112 + dg211) + 4.*dg112*dg212 + dg111*dg222)*ginv22 - 
        (dg123*(8.*dg112 + 2.*dg211) + 
           4.*(dg113*(dg122 + dg212) + dg112*(dg213 + dg312)) + 
           2.*(dg122*dg311 + dg111*(dg223 + dg322)))*ginv23 - 
        (dg133*(2.*dg112 + dg211) + 
           dg113*(8.*dg123 + 4.*dg213 + 2.*dg312) + 
           2.*(dg123*dg311 + dg112*dg313) + dg111*(dg233 + 2.*dg323))*
         ginv33) - ginv22*((dg122*dg213 + dg113*dg222 + 
           2.*(dg123*(dg122 + dg212) + dg112*dg223))*ginv33 + 
        ginv23*(3.*(dg122*dg212 + dg112*dg222) + 2.*pow2(dg122))) + 
     ginv23*(ddg1122 + ddg1212 - 
        ginv33*(dg133*(2.*dg122 + dg212) + 
           2.*(dg123*dg312 + dg122*dg313 + dg113*dg322) + 
           dg112*(dg233 + 2.*dg323) + 
           4.*(dg123*dg213 + dg113*dg223 + pow2(dg123)))) - 
     (4.*(dg112*dg113 + dg111*dg123) + 
        2.*(dg113*dg211 + dg112*dg311 + dg111*(dg213 + dg312)))*
      pow2(ginv13) - (dg123*(4.*dg122 + 2.*dg212) + 
        2.*(dg113*dg222 + dg122*(dg213 + dg312) + dg112*(dg223 + dg322)))*
      pow2(ginv23) - (dg133*(2.*dg123 + dg213) + 2.*dg123*dg313 + 
        dg113*(dg233 + 2.*dg323))*pow2(ginv33)) + 
  ginv22*(ddg1223*ginv33 + ginv23*
      (ddg1222 - (dg133*dg222 + dg123*(6.*dg223 + 2.*dg322) + 
           dg122*(dg233 + 2.*dg323))*ginv33) - 
     (3.*(dg123*dg222 + dg122*dg223) + 2.*dg122*dg322)*pow2(ginv23) - 
     (dg133*dg223 + dg123*(dg233 + 2.*dg323))*pow2(ginv33)) + 
  ginv23*((ddg1233 + 2.*ddg1323)*ginv33 - 
     (dg133*(2.*dg233 + 4.*dg323) + 4.*dg123*dg333)*pow2(ginv33)) + 
  ginv13*((ddg1133 + 2.*ddg1313)*ginv33 + 
     ginv23*(ddg1123 + ddg1213 + 2.*ddg1312 - 
        (dg133*(6.*dg123 + 3.*dg213 + 4.*dg312) + 6.*dg123*dg313 + 
           dg113*(3.*dg233 + 6.*dg323) + 4.*dg112*dg333)*ginv33) + 
     ginv22*(ddg1212 - (dg123*(2.*dg122 + 4.*dg212) + dg113*dg222 + 
           dg122*(dg213 + 2.*dg312) + dg112*(4.*dg223 + 2.*dg322))*ginv23 \
- ginv33*(dg133*dg212 + dg112*(dg233 + 2.*dg323) + 
           2.*(dg113*dg223 + dg123*(dg213 + dg312) + pow2(dg123)))) - 
     (dg122*dg212 + dg112*dg222)*pow2(ginv22) - 
     (4.*(dg123*dg312 + dg112*dg323) + 
        2.*(dg133*(dg122 + dg212) + dg123*dg213 + dg112*dg233 + 
           dg122*dg313 + dg113*(dg223 + dg322) + pow2(dg123)))*pow2(ginv23) \
- (4.*(dg133*dg313 + dg113*dg333) + 2.*pow2(dg133))*pow2(ginv33)) - 
  (dg133*dg222 + 2.*dg123*(dg223 + dg322) + dg122*(dg233 + 2.*dg323))*
   pow3(ginv23) - 2.*((dg111*(dg112*ginv23 + dg113*ginv33) + 
        ginv13*pow2(dg111))*pow2(ginv11) + 
     (dg113*dg311 + dg111*(dg133 + dg313) + pow2(dg113))*pow3(ginv13) + 
     dg133*dg333*pow3(ginv33))
;

dGfromgdu21
=
(ddg1211 - (4.*(dg113*dg211 + dg111*dg213) + 2.*dg211*dg311)*ginv13 - 
     2.*(dg112 + dg211)*dg212*ginv22 - 
     (2.*(dg113*dg212 + (dg112 + dg211)*dg213) + dg212*dg311 + 
        dg211*dg312)*ginv23 - (dg213*(2.*dg113 + dg311) + dg211*dg313)*
      ginv33 - ginv12*(4.*(dg112*dg211 + dg111*dg212) + 2.*pow2(dg211)))*
   pow2(ginv11) + (ddg1222 + ddg2212 - 
     (4.*(dg212*(dg123 + dg213) + (dg112 + dg211)*dg223) + dg222*dg311 + 
        2.*(dg122*dg213 + dg113*dg222 + dg212*dg312) + dg211*dg322)*ginv13 \
- (2.*dg122 + 6.*dg212)*dg222*ginv22 - 
     ((2.*dg122 + 4.*dg212)*dg223 + 
        dg222*(4.*dg213 + 2.*(dg123 + dg312)) + 2.*dg212*dg322)*ginv23 - 
     (dg223*(2.*(dg123 + dg213) + dg312) + dg222*dg313 + dg213*dg322 + 
        dg212*dg323)*ginv33)*pow2(ginv12) + 
  (ddg1233 + ddg2313 - (2.*((dg123 + dg213)*dg223 + dg212*dg233) + 
        dg223*dg312 + dg212*dg323)*ginv22 - 
     (dg233*(4.*dg213 + 2.*dg312) + 
        2.*(dg123*dg233 + dg223*(dg133 + dg313) + dg213*dg323 + 
           dg212*dg333))*ginv23 - 
     (dg233*(2.*dg133 + 3.*dg313) + 3.*dg213*dg333)*ginv33)*pow2(ginv13) + 
  ginv11*(ddg2313*ginv33 + ginv22*
      (ddg2212 - (dg222*(2.*dg213 + dg312) + dg212*(4.*dg223 + dg322))*
         ginv23 - (dg223*(2.*dg213 + dg312) + dg212*dg323)*ginv33) + 
     ginv23*(ddg2213 + ddg2312 - 
        (dg233*(2.*dg213 + dg312) + 2.*(dg223*dg313 + dg213*dg323) + 
           dg212*dg333)*ginv33) + 
     ginv13*(2.*ddg1213 + ddg2311 - 
        (2.*(dg112 + dg211)*dg223 + 
           dg212*(4.*dg213 + 2.*(dg123 + dg312)))*ginv22 - 
        (2.*(dg133*dg213 + dg113*dg233) + dg233*dg311 + 6.*dg213*dg313 + 
           dg211*dg333)*ginv33 - 
        ginv23*(2.*(dg133*dg212 + dg123*dg213 + dg113*dg223 + 
              (dg112 + dg211)*dg233) + dg223*dg311 + dg211*dg323 + 
           4.*(dg213*dg312 + dg212*dg313 + pow2(dg213)))) + 
     ginv12*(2.*ddg1212 + ddg2211 - 
        (6.*(dg113*dg212 + dg112*dg213) + 4.*dg111*dg223 + 
           3.*dg212*dg311 + dg211*(4.*dg123 + 6.*dg213 + 3.*dg312))*ginv13 \
- (2.*(dg123*dg212 + dg122*dg213 + dg113*dg222 + 
              (dg112 + dg211)*dg223) + dg222*dg311 + 
           dg212*(8.*dg213 + 4.*dg312) + dg211*dg322)*ginv23 - 
        ginv22*(2.*(dg122*dg212 + (dg112 + dg211)*dg222) + 
           6.*pow2(dg212)) - ginv33*
         (dg223*dg311 + dg211*dg323 + 
           2.*(dg113*dg223 + dg213*(dg123 + dg312) + dg212*dg313 + 
              pow2(dg213)))) - 
     (6.*dg112*dg212 + dg211*(2.*dg122 + 6.*dg212) + 2.*dg111*dg222)*
      pow2(ginv12) - (2.*(dg133*dg211 + dg111*dg233) + 
        dg213*(6.*dg113 + 3.*dg311) + 3.*dg211*dg313)*pow2(ginv13) - 
     2.*dg212*dg222*pow2(ginv22) - 
     (2.*(dg213*dg223 + dg212*dg233) + dg223*dg312 + dg222*dg313 + 
        dg213*dg322 + dg212*dg323)*pow2(ginv23) - 
     (dg233*dg313 + dg213*dg333)*pow2(ginv33)) + 
  ginv12*(ddg2323*ginv33 + ginv13*
      (2.*ddg1223 + ddg2213 + ddg2312 - 
        (2.*((dg123 + dg213)*dg222 + dg122*dg223) + dg222*dg312 + 
           dg212*(8.*dg223 + dg322))*ginv22 - 
        (dg223*(8.*dg213 + 4.*(dg123 + dg312)) + 
           2.*(dg122*dg233 + dg222*(dg133 + dg313) + dg213*dg322) + 
           4.*dg212*(dg233 + dg323))*ginv23 - 
        (2.*(dg133*dg223 + (dg123 + dg213)*dg233) + dg233*dg312 + 
           4.*(dg223*dg313 + dg213*dg323) + dg212*dg333)*ginv33) + 
     ginv23*(ddg2223 + ddg2322 - 
        (dg233*(2.*dg223 + dg322) + 4.*dg223*dg323 + dg222*dg333)*ginv33) + 
     ginv22*(ddg2222 - dg222*(6.*dg223 + 2.*dg322)*ginv23 - 
        ginv33*(dg223*dg322 + dg222*dg323 + 2.*pow2(dg223))) - 
     (4.*(dg123*dg213 + dg113*dg223) + 
        2.*((dg112 + dg211)*dg233 + dg223*dg311 + dg213*dg312 + 
           dg212*(dg133 + dg313) + dg211*dg323 + pow2(dg213)))*pow2(ginv13) \
- 2.*(pow2(dg222)*pow2(ginv22) + 
        (dg223*dg322 + dg222*(dg233 + dg323) + pow2(dg223))*pow2(ginv23)) - 
     (dg233*dg323 + dg223*dg333)*pow2(ginv33)) + 
  ginv13*(ddg2333*ginv33 + ginv22*
      (ddg2223 - 2.*dg223*(dg233 + dg323)*ginv33 - 
        ginv23*(dg223*dg322 + dg222*(2.*dg233 + dg323) + 4.*pow2(dg223))) + 
     ginv23*(ddg2233 + ddg2323 - 
        ginv33*(3.*(dg233*dg323 + dg223*dg333) + 2.*pow2(dg233))) - 
     (dg233*(4.*dg223 + dg322) + 2.*dg223*dg323 + dg222*dg333)*
      pow2(ginv23) - 2.*(dg222*dg223*pow2(ginv22) + 
        dg233*dg333*pow2(ginv33))) - 
  2.*(dg111*dg211*pow3(ginv11) + 
     (dg122*dg212 + (dg112 + dg211)*dg222 + pow2(dg212))*pow3(ginv12)) - 
  (dg233*dg311 + 2.*(dg113*dg233 + dg213*(dg133 + dg313)) + dg211*dg333)*
   pow3(ginv13)
;

dGfromgdu22
=
-((2.*dg111*dg211*ginv12 + (dg112*dg211 + dg111*dg212)*ginv22 + 
       (dg113*dg211 + dg111*dg213)*ginv23)*pow2(ginv11)) + 
  (ddg1212 + ddg2211 - (2.*(dg123*dg211 + dg112*dg213 + dg111*dg223 + 
           dg212*(dg113 + dg311)) + dg211*(4.*dg213 + 2.*dg312))*ginv13 - 
     (2.*(dg123*dg212 + dg122*dg213 + dg113*dg222 + dg112*dg223) + 
        dg222*dg311 + dg212*(8.*dg213 + 2.*dg312) + 
        dg211*(4.*dg223 + dg322))*ginv23 - 
     ginv22*(4.*dg211*dg222 + 3.*(dg122*dg212 + dg112*dg222) + 
        6.*pow2(dg212)) - ginv33*
      (dg223*(dg113 + dg311) + dg213*(dg123 + dg312) + dg212*dg313 + 
        dg211*dg323 + 2.*pow2(dg213)))*pow2(ginv12) - 
  ((dg112*dg233 + dg223*(dg113 + dg311) + dg213*(dg123 + dg312) + 
        dg212*(dg133 + dg313) + dg211*dg323)*ginv22 + 
     (dg233*dg311 + 2.*(dg113*dg233 + dg213*(dg133 + dg313)) + 
        dg211*dg333)*ginv23)*pow2(ginv13) + 
  (ddg2222 - dg222*(8.*dg223 + 2.*dg322)*ginv23 - 
     ginv33*(dg223*dg322 + dg222*dg323 + 2.*pow2(dg223)))*pow2(ginv22) + 
  (ddg2233 + ddg2323 - ginv33*(3.*(dg233*dg323 + dg223*dg333) + 
        2.*pow2(dg233)))*pow2(ginv23) + 
  ginv13*(ginv22*(ddg1223 + ddg2312 - 
        (dg122*dg233 + dg222*(dg133 + dg313) + dg213*dg322 + 
           4.*(dg223*(dg123 + dg213 + dg312) + dg212*(dg233 + dg323)))*
         ginv23 - (dg233*(dg123 + dg312) + dg223*(dg133 + 2.*dg313) + 
           2.*dg213*dg323 + dg212*dg333)*ginv33) + 
     ginv23*(ddg1233 + ddg2313 - 
        (dg233*(2.*dg133 + 3.*dg313) + 3.*dg213*dg333)*ginv33) - 
     ((dg122 + 4.*dg212)*dg223 + dg222*(dg123 + dg312) + dg212*dg322)*
      pow2(ginv22) - (dg233*(4.*dg213 + 2.*dg312) + 
        2.*(dg123*dg233 + dg223*(dg133 + dg313) + dg213*dg323 + 
           dg212*dg333))*pow2(ginv23)) + 
  ginv11*(-(ginv13*((2.*(dg113*dg212 + dg112*dg213) + dg111*dg223 + 
             dg212*dg311 + dg211*(dg123 + dg312))*ginv22 + 
          (dg111*dg233 + dg213*(4.*dg113 + dg311) + dg211*(dg133 + dg313))*
           ginv23)) + ginv12*(ddg1211 - 
        (3.*(dg113*dg211 + dg111*dg213) + 2.*dg211*dg311)*ginv13 - 
        (6.*dg112*dg212 + dg211*(dg122 + 4.*dg212) + dg111*dg222)*ginv22 - 
        (4.*(dg113*dg212 + dg112*dg213) + dg111*dg223 + dg212*dg311 + 
           dg211*(dg123 + 4.*dg213 + dg312))*ginv23 - 
        (dg213*(2.*dg113 + dg311) + dg211*dg313)*ginv33) + 
     ginv22*(ddg1212 - (dg122*dg213 + dg113*dg222 + 2.*dg112*dg223 + 
           dg212*(4.*dg213 + 2.*(dg123 + dg312)))*ginv23 - 
        (dg113*dg223 + dg213*(dg123 + dg312) + dg212*dg313)*ginv33) + 
     ginv23*(ddg1213 - (dg113*dg233 + dg213*(dg133 + 2.*dg313))*ginv33) - 
     (3.*(dg112*dg211 + dg111*dg212) + 2.*pow2(dg211))*pow2(ginv12) - 
     (dg122*dg212 + dg112*dg222 + 2.*pow2(dg212))*pow2(ginv22) - 
     (dg113*dg223 + dg112*dg233 + dg213*(dg123 + dg312) + 
        dg212*(dg133 + dg313) + 2.*pow2(dg213))*pow2(ginv23)) + 
  ginv23*(ddg2333*ginv33 - 2.*dg233*dg333*pow2(ginv33)) + 
  ginv12*(ddg2313*ginv33 + ginv22*
      (ddg1222 + 2.*ddg2212 - ((3.*dg122 + 12.*dg212)*dg223 + 
           dg222*(8.*dg213 + 3.*(dg123 + dg312)) + 3.*dg212*dg322)*ginv23 \
- (dg223*(4.*dg213 + 2.*(dg123 + dg312)) + dg222*dg313 + dg213*dg322 + 
           2.*dg212*dg323)*ginv33) + 
     ginv23*(ddg1223 + 2.*ddg2213 + ddg2312 - 
        (dg233*(dg123 + 4.*dg213 + dg312) + dg223*(dg133 + 4.*dg313) + 
           4.*dg213*dg323 + dg212*dg333)*ginv33) + 
     ginv13*(ddg1213 + ddg2311 - 
        (dg122*dg213 + dg222*(dg113 + dg311) + 
           4.*((dg112 + dg211)*dg223 + dg212*(dg123 + dg213 + dg312)) + 
           dg211*dg322)*ginv22 - 
        (dg233*(dg113 + dg311) + dg213*(dg133 + 4.*dg313) + dg211*dg333)*
         ginv33 - ginv23*(2.*(dg133*dg212 + dg112*dg233 + dg223*dg311 + 
              dg211*dg323) + 4.*
            (dg113*dg223 + dg211*dg233 + dg213*(dg123 + dg312) + 
              dg212*dg313 + pow2(dg213)))) - 
     (dg111*dg233 + 2.*dg213*(dg113 + dg311) + dg211*(dg133 + 2.*dg313))*
      pow2(ginv13) - (2.*dg122 + 8.*dg212)*dg222*pow2(ginv22) - 
     ((dg122 + 4.*dg212)*dg233 + dg223*(8.*dg213 + 2.*(dg123 + dg312)) + 
        dg222*(dg133 + 2.*dg313) + 2.*(dg213*dg322 + dg212*dg323))*
      pow2(ginv23) - (dg233*dg313 + dg213*dg333)*pow2(ginv33)) + 
  ginv22*(ddg2323*ginv33 + ginv23*
      (2.*ddg2223 + ddg2322 - (dg233*(4.*dg223 + dg322) + 
           6.*dg223*dg323 + dg222*dg333)*ginv33) - 
     (3.*dg223*dg322 + dg222*(4.*dg233 + 3.*dg323) + 6.*pow2(dg223))*
      pow2(ginv23) - (dg233*dg323 + dg223*dg333)*pow2(ginv33)) - 
  (2.*dg112*dg212 + dg211*(dg122 + 4.*dg212) + dg111*dg222)*pow3(ginv12) - 
  2.*pow2(dg222)*pow3(ginv22) - 
  (dg233*(4.*dg223 + dg322) + 2.*dg223*dg323 + dg222*dg333)*pow3(ginv23)
;

dGfromgdu23
=
-((2.*dg111*dg211*ginv13 + (dg112*dg211 + dg111*dg212)*ginv23 + 
       (dg113*dg211 + dg111*dg213)*ginv33)*pow2(ginv11)) - 
  ((2.*dg112*dg212 + dg211*(dg122 + 4.*dg212) + dg111*dg222)*ginv13 + 
     (dg122*dg213 + dg212*(dg123 + 2.*dg213) + dg113*dg222 + 
        (dg112 + 2.*dg211)*dg223)*ginv33 + 
     2.*ginv23*(dg122*dg212 + (dg112 + dg211)*dg222 + pow2(dg212)))*
   pow2(ginv12) + (ddg1213 + ddg2311 - 
     ((dg112 + 2.*dg211)*dg223 + dg212*(dg123 + 2.*(dg213 + dg312)))*
      ginv22 - (3.*(dg133*dg213 + dg113*dg233) + 6.*dg213*dg313 + 
        2.*(dg233*dg311 + dg211*dg333))*ginv33 - 
     ginv23*(4.*(dg213*dg312 + dg212*dg313) + 
        2.*(dg133*dg212 + dg123*dg213 + (dg112 + dg211)*dg233 + 
           dg223*(dg113 + dg311) + dg211*dg323 + pow2(dg213))))*pow2(ginv13) \
+ (ddg2223 + ddg2322 - (dg233*(6.*dg223 + 2.*dg322) + 6.*dg223*dg323 + 
        2.*dg222*dg333)*ginv33)*pow2(ginv23) + ddg2333*pow2(ginv33) + 
  ginv11*(ddg1213*ginv33 + ginv13*
      (ddg1211 - 2.*(dg112 + dg211)*dg212*ginv22 - 
        (4.*(dg113*dg212 + dg112*dg213) + dg111*dg223 + 2.*dg212*dg311 + 
           dg211*(dg123 + 2.*(dg213 + dg312)))*ginv23 - 
        (dg111*dg233 + dg213*(6.*dg113 + 2.*dg311) + 
           dg211*(dg133 + 2.*dg313))*ginv33) - 
     ginv12*((4.*dg112*dg212 + dg211*(dg122 + 2.*dg212) + dg111*dg222)*
         ginv23 + (dg211*(dg123 + 2.*dg213) + 
           2.*(dg113*dg212 + dg112*dg213) + dg111*dg223)*ginv33 + 
        ginv13*(3.*(dg112*dg211 + dg111*dg212) + 2.*pow2(dg211))) - 
     ginv22*((dg212*(dg123 + 2.*dg213) + dg112*dg223)*ginv33 + 
        ginv23*(dg122*dg212 + dg112*dg222 + 2.*pow2(dg212))) + 
     ginv23*(ddg1212 - ginv33*(dg112*dg233 + dg212*(dg133 + 2.*dg313) + 
           2.*(dg113*dg223 + dg213*(dg123 + dg312) + pow2(dg213)))) - 
     (3.*(dg113*dg211 + dg111*dg213) + 2.*dg211*dg311)*pow2(ginv13) - 
     (dg122*dg213 + dg113*dg222 + dg112*dg223 + 
        dg212*(dg123 + 2.*(dg213 + dg312)))*pow2(ginv23) - 
     (dg113*dg233 + dg213*(dg133 + 2.*dg313))*pow2(ginv33)) + 
  ginv22*(ddg2223*ginv33 + ginv23*
      (ddg2222 - ginv33*(2.*(dg223*dg322 + dg222*(dg233 + dg323)) + 
           6.*pow2(dg223))) - dg222*(6.*dg223 + 2.*dg322)*pow2(ginv23) - 
     2.*dg223*(dg233 + dg323)*pow2(ginv33)) + 
  ginv12*((ddg1223 + ddg2213)*ginv33 - 
     ginv22*((2.*dg122 + 6.*dg212)*dg222*ginv23 + 
        ((dg123 + 2.*dg213)*dg222 + (dg122 + 4.*dg212)*dg223)*ginv33) + 
     ginv23*(ddg1222 + ddg2212 - 
        ((dg122 + 2.*dg212)*dg233 + 
           dg223*(4.*dg123 + 8.*dg213 + 2.*dg312) + 
           dg222*(dg133 + 2.*dg313) + 2.*(dg213*dg322 + dg212*dg323))*
         ginv33) + ginv13*(ddg1212 + ddg2211 - 
        (4.*(dg112 + dg211)*dg223 + 
           dg212*(8.*dg213 + 4.*(dg123 + dg312)) + 
           2.*(dg122*dg213 + dg222*(dg113 + dg311) + dg211*dg322))*ginv23 \
- ginv22*(dg122*dg212 + (dg112 + 2.*dg211)*dg222 + 4.*pow2(dg212)) - 
        ginv33*((dg112 + 2.*dg211)*dg233 + dg212*(dg133 + 2.*dg313) + 
           2.*(dg223*dg311 + dg213*dg312 + dg211*dg323) + 
           4.*(dg123*dg213 + dg113*dg223 + pow2(dg213)))) - 
     (2.*(dg123*dg211 + dg112*dg213 + dg111*dg223 + 
           dg212*(dg113 + dg311)) + dg211*(4.*dg213 + 2.*dg312))*
      pow2(ginv13) - ((2.*dg122 + 4.*dg212)*dg223 + 
        dg222*(4.*dg213 + 2.*(dg123 + dg312)) + 2.*dg212*dg322)*
      pow2(ginv23) - ((dg123 + 2.*dg213)*dg233 + 
        dg223*(dg133 + 2.*dg313) + 2.*dg213*dg323)*pow2(ginv33)) + 
  ginv13*((ddg1233 + 2.*ddg2313)*ginv33 + 
     ginv22*(ddg2212 - ((dg122 + 8.*dg212)*dg223 + 
           dg222*(dg123 + 2.*(dg213 + dg312)) + 2.*dg212*dg322)*ginv23 - 
        (dg223*(4.*dg213 + 2.*(dg123 + dg312)) + 2.*dg212*(dg233 + dg323))*
         ginv33) + ginv23*(ddg1223 + ddg2213 + 2.*ddg2312 - 
        (3.*(dg133*dg223 + dg123*dg233) + dg233*(6.*dg213 + 4.*dg312) + 
           6.*(dg223*dg313 + dg213*dg323) + 4.*dg212*dg333)*ginv33) - 
     2.*dg212*dg222*pow2(ginv22) - 
     ((dg122 + 4.*dg212)*dg233 + dg223*(2.*dg123 + 4.*(dg213 + dg312)) + 
        dg222*(dg133 + 2.*dg313) + 2.*dg213*dg322 + 4.*dg212*dg323)*
      pow2(ginv23) - (dg233*(2.*dg133 + 4.*dg313) + 4.*dg213*dg333)*
      pow2(ginv33)) + ginv23*((ddg2233 + 2.*ddg2323)*ginv33 - 
     (4.*(dg233*dg323 + dg223*dg333) + 2.*pow2(dg233))*pow2(ginv33)) - 
  (dg111*dg233 + 2.*dg213*(dg113 + dg311) + dg211*(dg133 + 2.*dg313))*
   pow3(ginv13) - 2.*((dg222*dg223*ginv33 + ginv23*pow2(dg222))*
      pow2(ginv22) + (dg223*dg322 + dg222*(dg233 + dg323) + pow2(dg223))*
      pow3(ginv23) + dg233*dg333*pow3(ginv33))
;

dGfromgdu31
=
(ddg1311 - ((4.*dg112 + 2.*dg211)*dg311 + 4.*dg111*dg312)*ginv12 - 
     (dg212*dg311 + (2.*dg112 + dg211)*dg312)*ginv22 - 
     (dg311*(dg213 + 2.*dg312) + dg211*dg313 + 
        2.*(dg113*dg312 + dg112*dg313))*ginv23 - 
     2.*(dg113 + dg311)*dg313*ginv33 - 
     ginv13*(4.*(dg113*dg311 + dg111*dg313) + 2.*pow2(dg311)))*pow2(ginv11) \
+ (ddg1322 + ddg2312 - (2.*dg122*dg322 + 3.*(dg222*dg312 + dg212*dg322))*
      ginv22 - ((2.*dg213 + 4.*dg312)*dg322 + 
        2.*(dg223*dg312 + dg222*dg313 + dg123*dg322 + 
           (dg122 + dg212)*dg323))*ginv23 - 
     (dg313*(dg223 + 2.*dg322) + (dg213 + 2.*(dg123 + dg312))*dg323)*
      ginv33 - ginv13*(4.*(dg123*dg312 + dg112*dg323) + 
        2.*(dg213*dg312 + (dg122 + dg212)*dg313 + dg113*dg322 + 
           dg311*(dg223 + dg322) + dg211*dg323 + pow2(dg312))))*pow2(ginv12) \
+ (ddg1333 + ddg3313 - (dg233*dg312 + dg223*dg313 + 
        (dg213 + 2.*(dg123 + dg312))*dg323 + dg212*dg333)*ginv22 - 
     (2.*(dg233*dg313 + dg133*dg323 + (dg123 + dg213)*dg333) + 
        4.*(dg313*dg323 + dg312*dg333))*ginv23 - 
     (2.*dg133 + 6.*dg313)*dg333*ginv33)*pow2(ginv13) + 
  ginv11*(ddg3313*ginv33 + ginv22*
      (ddg2312 - (dg222*dg313 + dg213*dg322 + 
           2.*(dg312*(dg223 + dg322) + dg212*dg323))*ginv23 - 
        (dg223*dg313 + (dg213 + 2.*dg312)*dg323)*ginv33) + 
     ginv23*(ddg2313 + ddg3312 - 
        (dg313*(dg233 + 4.*dg323) + (dg213 + 2.*dg312)*dg333)*ginv33) + 
     ginv12*(2.*ddg1312 + ddg2311 - 
        (dg311*(4.*dg123 + 3.*dg213 + 6.*dg312) + 3.*dg211*dg313 + 
           6.*(dg113*dg312 + dg112*dg313) + 4.*dg111*dg323)*ginv13 - 
        (dg222*dg311 + (2.*dg122 + 6.*dg212)*dg312 + 
           (2.*dg112 + dg211)*dg322)*ginv22 - 
        (4.*dg312*dg313 + 2.*((dg123 + dg213)*dg313 + 
              (dg113 + dg311)*dg323))*ginv33 - 
        ginv23*((2.*dg123 + 4.*dg213)*dg312 + dg311*(dg223 + 2.*dg322) + 
           dg211*dg323 + 2.*(dg122*dg313 + dg113*dg322 + dg112*dg323) + 
           4.*(dg212*dg313 + pow2(dg312)))) + 
     ginv13*(2.*ddg1313 + ddg3311 - 
        ((4.*dg213 + 8.*dg312)*dg313 + dg311*(dg233 + 2.*dg323) + 
           dg211*dg333 + 2.*(dg133*dg312 + dg123*dg313 + dg113*dg323 + 
              dg112*dg333))*ginv23 - 
        ginv22*(dg223*dg311 + dg211*dg323 + 
           2.*((dg123 + dg213)*dg312 + dg212*dg313 + dg112*dg323 + 
              pow2(dg312))) - ginv33*
         (2.*(dg133*dg313 + (dg113 + dg311)*dg333) + 6.*pow2(dg313))) - 
     ((2.*dg122 + 3.*dg212)*dg311 + (6.*dg112 + 3.*dg211)*dg312 + 
        2.*dg111*dg322)*pow2(ginv12) - 
     (6.*dg113*dg313 + dg311*(2.*dg133 + 6.*dg313) + 2.*dg111*dg333)*
      pow2(ginv13) - (dg222*dg312 + dg212*dg322)*pow2(ginv22) - 
     (dg313*(dg223 + 2.*dg322) + dg213*dg323 + dg312*(dg233 + 2.*dg323) + 
        dg212*dg333)*pow2(ginv23) - 2.*dg313*dg333*pow2(ginv33)) + 
  ginv12*(ddg3323*ginv33 + ginv13*
      (2.*ddg1323 + ddg2313 + ddg3312 - 
        (dg222*dg313 + (2.*dg123 + dg213)*dg322 + 
           dg312*(4.*dg223 + 2.*dg322) + (2.*dg122 + 4.*dg212)*dg323)*
         ginv22 - ((4.*dg213 + 8.*dg312)*dg323 + 
           4.*(dg313*(dg223 + dg322) + dg123*dg323) + 
           2.*(dg233*dg312 + dg133*dg322 + (dg122 + dg212)*dg333))*ginv23 \
- (dg313*(dg233 + 8.*dg323) + (dg213 + 2.*dg312)*dg333 + 
           2.*(dg133*dg323 + dg123*dg333))*ginv33) + 
     ginv22*(ddg2322 - 2.*(dg223 + dg322)*dg323*ginv33 - 
        ginv23*(3.*(dg223*dg322 + dg222*dg323) + 2.*pow2(dg322))) + 
     ginv23*(ddg2323 + ddg3322 - 
        ginv33*(dg233*dg323 + (dg223 + 2.*dg322)*dg333 + 4.*pow2(dg323))) - 
     (dg311*(dg233 + 4.*dg323) + 
        4.*((dg123 + dg312)*dg313 + dg113*dg323) + dg211*dg333 + 
        2.*(dg133*dg312 + dg213*dg313 + dg112*dg333))*pow2(ginv13) - 
     (2.*dg223*dg323 + dg322*(dg233 + 4.*dg323) + dg222*dg333)*
      pow2(ginv23) - 2.*(dg222*dg322*pow2(ginv22) + 
        dg323*dg333*pow2(ginv33))) + 
  ginv13*(ddg3333*ginv33 + ginv23*
      (ddg2333 + ddg3323 - (2.*dg233 + 6.*dg323)*dg333*ginv33) + 
     ginv22*(ddg2323 - (4.*dg223*dg323 + dg322*(dg233 + 2.*dg323) + 
           dg222*dg333)*ginv23 - 
        ginv33*(dg233*dg323 + dg223*dg333 + 2.*pow2(dg323))) - 
     (dg223*dg322 + dg222*dg323)*pow2(ginv22) - 
     2.*((dg233*dg323 + (dg223 + dg322)*dg333 + pow2(dg323))*pow2(ginv23) + 
        pow2(dg333)*pow2(ginv33))) - 
  (dg222*dg311 + dg211*dg322 + 2.*((dg122 + dg212)*dg312 + dg112*dg322))*
   pow3(ginv12) - 2.*(dg111*dg311*pow3(ginv11) + 
     (dg133*dg313 + (dg113 + dg311)*dg333 + pow2(dg313))*pow3(ginv13))
;

dGfromgdu32
=
-((2.*dg111*dg311*ginv12 + (dg112*dg311 + dg111*dg312)*ginv22 + 
       (dg113*dg311 + dg111*dg313)*ginv23)*pow2(ginv11)) + 
  (ddg1312 + ddg2311 - (4.*dg311*dg312 + 
        2.*((dg123 + dg213)*dg311 + dg113*dg312 + 
           (dg112 + dg211)*dg313 + dg111*dg323))*ginv13 - 
     ((3.*dg122 + 6.*dg212)*dg312 + 3.*dg112*dg322 + 
        2.*(dg222*dg311 + dg211*dg322))*ginv22 - 
     ((dg123 + 2.*(dg213 + dg312))*dg313 + (dg113 + 2.*dg311)*dg323)*
      ginv33 - ginv23*(4.*(dg213*dg312 + dg212*dg313) + 
        2.*(dg123*dg312 + dg122*dg313 + dg113*dg322 + 
           dg311*(dg223 + dg322) + (dg112 + dg211)*dg323 + pow2(dg312))))*
   pow2(ginv12) - ((dg123*dg313 + dg312*(dg133 + 2.*dg313) + 
        (dg113 + 2.*dg311)*dg323 + dg112*dg333)*ginv22 + 
     2.*ginv23*(dg133*dg313 + (dg113 + dg311)*dg333 + pow2(dg313)))*
   pow2(ginv13) + (ddg2322 - 2.*(dg223 + dg322)*dg323*ginv33 - 
     ginv23*(4.*(dg223*dg322 + dg222*dg323) + 2.*pow2(dg322)))*pow2(ginv22) \
+ (ddg2333 + ddg3323 - (2.*dg233 + 6.*dg323)*dg333*ginv33)*pow2(ginv23) + 
  ginv11*(-(ginv13*((dg311*(dg123 + 2.*dg312) + 
             2.*(dg113*dg312 + dg112*dg313) + dg111*dg323)*ginv22 + 
          (4.*dg113*dg313 + dg311*(dg133 + 2.*dg313) + dg111*dg333)*ginv23)\
) + ginv12*(ddg1311 - ((dg122 + 2.*dg212)*dg311 + 
           (6.*dg112 + 2.*dg211)*dg312 + dg111*dg322)*ginv22 - 
        (dg311*(dg123 + 2.*(dg213 + dg312)) + 2.*dg211*dg313 + 
           4.*(dg113*dg312 + dg112*dg313) + dg111*dg323)*ginv23 - 
        2.*(dg113 + dg311)*dg313*ginv33 - 
        ginv13*(3.*(dg113*dg311 + dg111*dg313) + 2.*pow2(dg311))) + 
     ginv22*(ddg1312 - ((dg123 + 2.*dg312)*dg313 + dg113*dg323)*ginv33 - 
        ginv23*(dg122*dg313 + dg113*dg322 + 
           2.*((dg123 + dg213)*dg312 + dg212*dg313 + dg112*dg323 + 
              pow2(dg312)))) + 
     ginv23*(ddg1313 - ginv33*
         (dg133*dg313 + dg113*dg333 + 2.*pow2(dg313))) - 
     ((3.*dg112 + 2.*dg211)*dg311 + 3.*dg111*dg312)*pow2(ginv12) - 
     ((dg122 + 2.*dg212)*dg312 + dg112*dg322)*pow2(ginv22) - 
     (dg133*dg312 + (dg123 + 2.*(dg213 + dg312))*dg313 + dg113*dg323 + 
        dg112*dg333)*pow2(ginv23)) + 
  ginv13*(ginv23*(ddg1333 + ddg3313 - (2.*dg133 + 6.*dg313)*dg333*ginv33) + 
     ginv22*(ddg1323 + ddg3312 - 
        (dg133*dg322 + (4.*dg123 + 2.*dg213 + 8.*dg312)*dg323 + 
           dg122*dg333 + 2.*(dg233*dg312 + dg313*(dg223 + dg322) + 
              dg212*dg333))*ginv23 - 
        ((dg133 + 4.*dg313)*dg323 + (dg123 + 2.*dg312)*dg333)*ginv33) - 
     (dg123*dg322 + dg122*dg323 + 
        2.*(dg312*(dg223 + dg322) + dg212*dg323))*pow2(ginv22) - 
     (2.*(dg233*dg313 + dg133*dg323 + (dg123 + dg213)*dg333) + 
        4.*(dg313*dg323 + dg312*dg333))*pow2(ginv23)) + 
  ginv12*(ddg3313*ginv33 + ginv22*
      (ddg1322 + 2.*ddg2312 - (4.*(dg222*dg313 + dg213*dg322) + 
           3.*(dg123*dg322 + dg122*dg323) + 
           6.*(dg312*(dg223 + dg322) + dg212*dg323))*ginv23 - 
        ((2.*dg213 + 4.*dg312)*dg323 + 
           2.*(dg313*(dg223 + dg322) + dg123*dg323))*ginv33) + 
     ginv23*(ddg1323 + 2.*ddg2313 + ddg3312 - 
        (dg133*dg323 + dg313*(2.*dg233 + 8.*dg323) + 
           (dg123 + 2.*(dg213 + dg312))*dg333)*ginv33) + 
     ginv13*(ddg1313 + ddg3311 - 
        (8.*dg312*dg313 + 4.*
            ((dg123 + dg213)*dg313 + (dg113 + dg311)*dg323) + 
           2.*(dg233*dg311 + dg133*dg312 + (dg112 + dg211)*dg333))*ginv23 \
- ginv22*(dg122*dg313 + dg113*dg322 + 
           2.*(dg213*dg312 + dg212*dg313 + dg311*(dg223 + dg322) + 
              dg211*dg323) + 4.*(dg123*dg312 + dg112*dg323 + pow2(dg312))) \
- ginv33*(dg133*dg313 + (dg113 + 2.*dg311)*dg333 + 4.*pow2(dg313))) - 
     (2.*dg113*dg313 + dg311*(dg133 + 4.*dg313) + dg111*dg333)*
      pow2(ginv13) - (2.*dg122*dg322 + 4.*(dg222*dg312 + dg212*dg322))*
      pow2(ginv22) - (dg133*dg322 + 
        4.*(dg313*(dg223 + dg322) + (dg213 + dg312)*dg323) + 
        dg122*dg333 + 2.*(dg233*dg312 + dg123*dg323 + dg212*dg333))*
      pow2(ginv23) - 2.*dg313*dg333*pow2(ginv33)) + 
  ginv22*(ddg3323*ginv33 + ginv23*
      (2.*ddg2323 + ddg3322 - ginv33*
         (2.*(dg233*dg323 + (dg223 + dg322)*dg333) + 6.*pow2(dg323))) - 
     (6.*dg223*dg323 + dg322*(2.*dg233 + 6.*dg323) + 2.*dg222*dg333)*
      pow2(ginv23) - 2.*dg323*dg333*pow2(ginv33)) + 
  ginv23*(ddg3333*ginv33 - 2.*pow2(dg333)*pow2(ginv33)) - 
  ((dg122 + 2.*dg212)*dg311 + 2.*(dg112 + dg211)*dg312 + dg111*dg322)*
   pow3(ginv12) - 2.*(dg222*dg322*pow3(ginv22) + 
     (dg233*dg323 + (dg223 + dg322)*dg333 + pow2(dg323))*pow3(ginv23))
;

dGfromgdu33
=
-((2.*dg111*dg311*ginv13 + (dg112*dg311 + dg111*dg312)*ginv23 + 
       (dg113*dg311 + dg111*dg313)*ginv33)*pow2(ginv11)) - 
  (((dg122 + 2.*dg212)*dg311 + 2.*(dg112 + dg211)*dg312 + dg111*dg322)*
      ginv13 + (dg222*dg311 + dg211*dg322 + 
        2.*((dg122 + dg212)*dg312 + dg112*dg322))*ginv23 + 
     (dg223*dg311 + (dg123 + dg213)*dg312 + (dg122 + dg212)*dg313 + 
        dg113*dg322 + (dg112 + dg211)*dg323)*ginv33)*pow2(ginv12) + 
  (ddg1313 + ddg3311 - ((2.*dg213 + 8.*dg312)*dg313 + 
        dg311*(dg233 + 4.*dg323) + dg211*dg333 + 
        2.*(dg133*dg312 + dg123*dg313 + dg113*dg323 + dg112*dg333))*ginv23 \
- ginv22*(dg223*dg311 + (dg123 + dg213)*dg312 + dg212*dg313 + 
        (dg112 + dg211)*dg323 + 2.*pow2(dg312)) - 
     ginv33*(4.*dg311*dg333 + 3.*(dg133*dg313 + dg113*dg333) + 
        6.*pow2(dg313)))*pow2(ginv13) - 
  (2.*dg222*dg322*ginv23 + (dg223*dg322 + dg222*dg323)*ginv33)*
   pow2(ginv22) + (ddg2323 + ddg3322 - 
     ginv33*(4.*dg322*dg333 + 3.*(dg233*dg323 + dg223*dg333) + 
        6.*pow2(dg323)))*pow2(ginv23) + ddg3333*pow2(ginv33) + 
  ginv13*((ddg1333 + 2.*ddg3313)*ginv33 + 
     ginv22*(ddg2312 - (dg222*dg313 + (dg123 + dg213)*dg322 + 
           dg122*dg323 + 4.*(dg312*(dg223 + dg322) + dg212*dg323))*ginv23 \
- (dg312*(dg233 + 4.*dg323) + 2.*(dg223*dg313 + (dg123 + dg213)*dg323) + 
           dg212*dg333)*ginv33) + 
     ginv23*(ddg1323 + ddg2313 + 2.*ddg3312 - 
        (12.*dg313*dg323 + (3.*dg213 + 8.*dg312)*dg333 + 
           3.*(dg233*dg313 + dg133*dg323 + dg123*dg333))*ginv33) - 
     (dg222*dg312 + dg212*dg322)*pow2(ginv22) - 
     ((dg133 + 4.*dg313)*dg322 + (2.*dg213 + 8.*dg312)*dg323 + 
        dg122*dg333 + 2.*(dg233*dg312 + dg223*dg313 + dg123*dg323 + 
           dg212*dg333))*pow2(ginv23) - 
     (2.*dg133 + 8.*dg313)*dg333*pow2(ginv33)) + 
  ginv23*((ddg2333 + 2.*ddg3323)*ginv33 - 
     (2.*dg233 + 8.*dg323)*dg333*pow2(ginv33)) + 
  ginv12*((ddg1323 + ddg2313)*ginv33 - 
     ginv22*((2.*dg122*dg322 + 3.*(dg222*dg312 + dg212*dg322))*ginv23 + 
        (dg222*dg313 + (dg123 + dg213)*dg322 + dg122*dg323 + 
           2.*(dg223*dg312 + dg212*dg323))*ginv33) + 
     ginv23*(ddg1322 + ddg2312 - 
        (dg233*dg312 + dg133*dg322 + 
           4.*(dg313*(dg223 + dg322) + (dg123 + dg213 + dg312)*dg323) + 
           (dg122 + dg212)*dg333)*ginv33) + 
     ginv13*(ddg1312 + ddg2311 - 
        (dg222*dg311 + (dg122 + 4.*dg212)*dg312 + (dg112 + dg211)*dg322)*
         ginv22 - (dg133*dg312 + dg311*(dg233 + 4.*dg323) + 
           4.*((dg123 + dg213 + dg312)*dg313 + dg113*dg323) + 
           (dg112 + dg211)*dg333)*ginv33 - 
        ginv23*(2.*(dg223*dg311 + dg122*dg313 + dg113*dg322 + 
              dg211*dg323) + 4.*
            ((dg123 + dg213)*dg312 + dg212*dg313 + dg311*dg322 + 
              dg112*dg323 + pow2(dg312)))) - 
     (4.*dg311*dg312 + 2.*((dg123 + dg213)*dg311 + dg113*dg312 + 
           (dg112 + dg211)*dg313 + dg111*dg323))*pow2(ginv13) - 
     ((2.*dg213 + 4.*dg312)*dg322 + 
        2.*(dg223*dg312 + dg222*dg313 + dg123*dg322 + 
           (dg122 + dg212)*dg323))*pow2(ginv23) - 
     (dg133*dg323 + dg313*(dg233 + 4.*dg323) + (dg123 + dg213)*dg333)*
      pow2(ginv33)) + ginv11*(ddg1313*ginv33 - 
     ginv12*(((3.*dg112 + 2.*dg211)*dg311 + 3.*dg111*dg312)*ginv13 + 
        ((dg122 + dg212)*dg311 + (4.*dg112 + dg211)*dg312 + dg111*dg322)*
         ginv23 + ((dg123 + dg213)*dg311 + dg211*dg313 + 
           2.*(dg113*dg312 + dg112*dg313) + dg111*dg323)*ginv33) - 
     ginv22*(((dg122 + 2.*dg212)*dg312 + dg112*dg322)*ginv23 + 
        ((dg123 + dg213)*dg312 + dg212*dg313 + dg112*dg323)*ginv33) + 
     ginv13*(ddg1311 - (dg212*dg311 + (2.*dg112 + dg211)*dg312)*ginv22 - 
        ((dg123 + dg213)*dg311 + 4.*(dg113 + dg311)*dg312 + 
           (4.*dg112 + dg211)*dg313 + dg111*dg323)*ginv23 - 
        (6.*dg113*dg313 + dg311*(dg133 + 4.*dg313) + dg111*dg333)*ginv33) + 
     ginv23*(ddg1312 - (dg312*(dg133 + 4.*dg313) + 
           2.*((dg123 + dg213)*dg313 + dg113*dg323) + dg112*dg333)*ginv33) \
- (3.*(dg113*dg311 + dg111*dg313) + 2.*pow2(dg311))*pow2(ginv13) - 
     ((dg123 + dg213)*dg312 + (dg122 + dg212)*dg313 + dg113*dg322 + 
        dg112*dg323 + 2.*pow2(dg312))*pow2(ginv23) - 
     (dg133*dg313 + dg113*dg333 + 2.*pow2(dg313))*pow2(ginv33)) + 
  ginv22*(ddg2323*ginv33 + ginv23*
      (ddg2322 - (6.*dg223*dg323 + dg322*(dg233 + 4.*dg323) + dg222*dg333)*
         ginv33) - (3.*(dg223*dg322 + dg222*dg323) + 2.*pow2(dg322))*
      pow2(ginv23) - (dg233*dg323 + dg223*dg333 + 2.*pow2(dg323))*
      pow2(ginv33)) - (2.*dg113*dg313 + dg311*(dg133 + 4.*dg313) + 
     dg111*dg333)*pow3(ginv13) - 
  (2.*dg223*dg323 + dg322*(dg233 + 4.*dg323) + dg222*dg333)*pow3(ginv23) - 
  2.*pow2(dg333)*pow3(ginv33)
;

R11
=
gammado111*Gfromg1 + gammado112*Gfromg2 + gammado113*Gfromg3 + 
  (-0.5*ddg1111 + 3.*gamma111*gammado111 + 
     2.*(gamma211*gammado112 + gamma311*gammado113) + 
     gamma211*gammado211 + gamma311*gammado311)*ginv11 + 
  (-ddg1211 + 3.*(gamma112*gammado111 + gamma111*gammado112) + 
     2.*(gamma212*gammado112 + gamma312*gammado113 + 
        gamma211*gammado122 + gamma311*gammado123) + gamma212*gammado211 + 
     gamma211*gammado212 + gamma312*gammado311 + gamma311*gammado312)*ginv12 \
+ (-ddg1311 + 3.*(gamma113*gammado111 + gamma111*gammado113) + 
     2.*(gamma213*gammado112 + gamma313*gammado113 + 
        gamma211*gammado123 + gamma311*gammado133) + gamma213*gammado211 + 
     gamma211*gammado213 + gamma313*gammado311 + gamma311*gammado313)*ginv13 \
+ (-0.5*ddg2211 + 3.*gamma112*gammado112 + 
     2.*(gamma212*gammado122 + gamma312*gammado123) + 
     gamma212*gammado212 + gamma312*gammado312)*ginv22 + 
  (-ddg2311 + 3.*(gamma113*gammado112 + gamma112*gammado113) + 
     2.*(gamma213*gammado122 + (gamma212 + gamma313)*gammado123 + 
        gamma312*gammado133) + gamma213*gammado212 + gamma212*gammado213 + 
     gamma313*gammado312 + gamma312*gammado313)*ginv23 + 
  (-0.5*ddg3311 + 3.*gamma113*gammado113 + 
     2.*(gamma213*gammado123 + gamma313*gammado133) + 
     gamma213*gammado213 + gamma313*gammado313)*ginv33 + dG11*g11[ijk] + 
  dG12*g12[ijk] + dG13*g13[ijk]
;

R12
=
(-0.5*ddg1112 + gamma112*gammado111 + (gamma111 + gamma212)*gammado112 + 
     gamma312*gammado113 + gamma111*gammado211 + 2.*gamma211*gammado212 + 
     gamma311*(gammado213 + gammado312))*ginv11 + 
  (-ddg1212 + gamma122*gammado111 + (2.*gamma112 + gamma222)*gammado112 + 
     gamma322*gammado113 + (gamma111 + gamma212)*gammado122 + 
     gamma112*gammado211 + (gamma111 + 2.*gamma212)*gammado212 + 
     2.*gamma211*gammado222 + gamma312*
      (gammado123 + gammado213 + gammado312) + 
     gamma311*(gammado223 + gammado322))*ginv12 + 
  (-ddg1312 + gamma123*gammado111 + (gamma113 + gamma223)*gammado112 + 
     (gamma112 + gamma323)*gammado113 + (gamma111 + gamma212)*gammado123 + 
     gamma312*gammado133 + gamma113*gammado211 + 
     (gamma111 + gamma313)*gammado213 + 
     2.*(gamma213*gammado212 + gamma211*gammado223) + 
     gamma313*gammado312 + gamma311*(gammado233 + gammado323))*ginv13 + 
  (-0.5*ddg2212 + gamma122*gammado112 + (gamma112 + gamma222)*gammado122 + 
     gamma322*gammado123 + gamma112*gammado212 + 2.*gamma212*gammado222 + 
     gamma312*(gammado223 + gammado322))*ginv22 + 
  (-ddg2312 + gamma123*gammado112 + gamma122*gammado113 + 
     (gamma113 + gamma223)*gammado122 + 
     (gamma112 + gamma222 + gamma323)*gammado123 + gamma322*gammado133 + 
     gamma113*gammado212 + gamma112*gammado213 + 
     2.*(gamma213*gammado222 + gamma212*gammado223) + 
     gamma313*(gammado223 + gammado322) + 
     gamma312*(gammado233 + gammado323))*ginv23 + 
  (-0.5*ddg3312 + gamma123*gammado113 + (gamma113 + gamma223)*gammado123 + 
     gamma323*gammado133 + gamma113*gammado213 + 2.*gamma213*gammado223 + 
     gamma313*(gammado233 + gammado323))*ginv33 + 
  0.5*((gammado112 + gammado211)*Gfromg1 + 
     (gammado122 + gammado212)*Gfromg2 + (gammado123 + gammado213)*Gfromg3 + 
     dG21*g11[ijk] + (dG11 + dG22)*g12[ijk] + dG23*g13[ijk] + 
     dG12*g22[ijk] + dG13*g23[ijk])
;

R13
=
(-0.5*ddg1113 + gamma113*gammado111 + gamma213*gammado112 + 
     (gamma111 + gamma313)*gammado113 + gamma111*gammado311 + 
     gamma211*(gammado213 + gammado312) + 2.*gamma311*gammado313)*ginv11 + 
  (-ddg1213 + gamma123*gammado111 + (gamma113 + gamma223)*gammado112 + 
     (gamma112 + gamma323)*gammado113 + gamma213*gammado122 + 
     (gamma111 + gamma313)*gammado123 + gamma112*gammado311 + 
     gamma111*gammado312 + gamma212*(gammado213 + gammado312) + 
     gamma211*(gammado223 + gammado322) + 
     2.*(gamma312*gammado313 + gamma311*gammado323))*ginv12 + 
  (-ddg1313 + gamma133*gammado111 + gamma233*gammado112 + 
     (2.*gamma113 + gamma333)*gammado113 + 
     (gamma111 + gamma313)*gammado133 + gamma113*gammado311 + 
     gamma213*(gammado123 + gammado213 + gammado312) + 
     (gamma111 + 2.*gamma313)*gammado313 + 
     gamma211*(gammado233 + gammado323) + 2.*gamma311*gammado333)*ginv13 + 
  (-0.5*ddg2213 + gamma123*gammado112 + gamma223*gammado122 + 
     (gamma112 + gamma323)*gammado123 + gamma112*gammado312 + 
     gamma212*(gammado223 + gammado322) + 2.*gamma312*gammado323)*ginv22 + 
  (-ddg2313 + gamma133*gammado112 + gamma123*gammado113 + 
     gamma233*gammado122 + (gamma113 + gamma223 + gamma333)*gammado123 + 
     (gamma112 + gamma323)*gammado133 + gamma113*gammado312 + 
     gamma112*gammado313 + gamma213*(gammado223 + gammado322) + 
     gamma212*(gammado233 + gammado323) + 
     2.*(gamma313*gammado323 + gamma312*gammado333))*ginv23 + 
  (-0.5*ddg3313 + gamma133*gammado113 + gamma233*gammado123 + 
     (gamma113 + gamma333)*gammado133 + gamma113*gammado313 + 
     gamma213*(gammado233 + gammado323) + 2.*gamma313*gammado333)*ginv33 + 
  0.5*((gammado113 + gammado311)*Gfromg1 + 
     (gammado123 + gammado312)*Gfromg2 + (gammado133 + gammado313)*Gfromg3 + 
     dG31*g11[ijk] + dG32*g12[ijk] + (dG11 + dG33)*g13[ijk] + 
     dG12*g23[ijk] + dG13*g33[ijk])
;

R22
=
gammado212*Gfromg1 + gammado222*Gfromg2 + gammado223*Gfromg3 + 
  (-0.5*ddg1122 + gamma112*(gammado112 + 2.*gammado211) + 
     3.*gamma212*gammado212 + gamma312*(2.*gammado213 + gammado312))*ginv11 \
+ (-ddg1222 + gamma122*(gammado112 + 2.*gammado211) + 
     gamma112*(gammado122 + 2.*gammado212) + 
     3.*(gamma222*gammado212 + gamma212*gammado222) + 
     2.*(gamma322*gammado213 + gamma312*gammado223) + 
     gamma322*gammado312 + gamma312*gammado322)*ginv12 + 
  (-ddg1322 + gamma123*(gammado112 + 2.*gammado211) + 
     gamma112*(gammado123 + 2.*gammado213) + 
     3.*(gamma223*gammado212 + gamma212*gammado223) + 
     2.*(gamma323*gammado213 + gamma312*gammado233) + 
     gamma323*gammado312 + gamma312*gammado323)*ginv13 + 
  (-0.5*ddg2222 + gamma122*(gammado122 + 2.*gammado212) + 
     3.*gamma222*gammado222 + gamma322*(2.*gammado223 + gammado322))*ginv22 \
+ (-ddg2322 + gamma123*(gammado122 + 2.*gammado212) + 
     gamma122*(gammado123 + 2.*gammado213) + 
     3.*(gamma223*gammado222 + gamma222*gammado223) + 
     2.*(gamma323*gammado223 + gamma322*gammado233) + 
     gamma323*gammado322 + gamma322*gammado323)*ginv23 + 
  (-0.5*ddg3322 + gamma123*(gammado123 + 2.*gammado213) + 
     3.*gamma223*gammado223 + gamma323*(2.*gammado233 + gammado323))*ginv33 \
+ dG21*g12[ijk] + dG22*g22[ijk] + dG23*g23[ijk]
;

R23
=
(-0.5*ddg1123 + gamma113*gammado211 + gamma213*gammado212 + 
     (gamma212 + gamma313)*gammado213 + 
     gamma112*(gammado113 + gammado311) + gamma212*gammado312 + 
     2.*gamma312*gammado313)*ginv11 + 
  (-ddg1223 + gamma123*gammado211 + (gamma113 + gamma223)*gammado212 + 
     (gamma222 + gamma323)*gammado213 + gamma213*gammado222 + 
     (gamma212 + gamma313)*gammado223 + 
     gamma122*(gammado113 + gammado311) + gamma222*gammado312 + 
     gamma112*(gammado123 + gammado312) + gamma212*gammado322 + 
     2.*(gamma322*gammado313 + gamma312*gammado323))*ginv12 + 
  (-ddg1323 + gamma133*gammado211 + gamma233*gammado212 + 
     (gamma113 + gamma223 + gamma333)*gammado213 + gamma213*gammado223 + 
     (gamma212 + gamma313)*gammado233 + 
     gamma123*(gammado113 + gammado311) + gamma223*gammado312 + 
     gamma112*(gammado133 + gammado313) + gamma212*gammado323 + 
     2.*(gamma323*gammado313 + gamma312*gammado333))*ginv13 + 
  (-0.5*ddg2223 + gamma123*gammado212 + gamma223*gammado222 + 
     (gamma222 + gamma323)*gammado223 + 
     gamma122*(gammado123 + gammado312) + gamma222*gammado322 + 
     2.*gamma322*gammado323)*ginv22 + 
  (-ddg2323 + gamma133*gammado212 + gamma233*gammado222 + 
     (2.*gamma223 + gamma333)*gammado223 + 
     (gamma222 + gamma323)*gammado233 + 
     gamma123*(gammado123 + gammado213 + gammado312) + 
     gamma122*(gammado133 + gammado313) + gamma223*gammado322 + 
     (gamma222 + 2.*gamma323)*gammado323 + 2.*gamma322*gammado333)*ginv23 + 
  (-0.5*ddg3323 + gamma133*gammado213 + gamma233*gammado223 + 
     (gamma223 + gamma333)*gammado233 + 
     gamma123*(gammado133 + gammado313) + gamma223*gammado323 + 
     2.*gamma323*gammado333)*ginv33 + 
  0.5*((gammado213 + gammado312)*Gfromg1 + 
     (gammado223 + gammado322)*Gfromg2 + (gammado233 + gammado323)*Gfromg3 + 
     dG31*g12[ijk] + dG21*g13[ijk] + dG32*g22[ijk] + 
     (dG22 + dG33)*g23[ijk] + dG23*g33[ijk])
;

R33
=
gammado313*Gfromg1 + gammado323*Gfromg2 + gammado333*Gfromg3 + 
  (-0.5*ddg1133 + gamma113*(gammado113 + 2.*gammado311) + 
     gamma213*(gammado213 + 2.*gammado312) + 3.*gamma313*gammado313)*ginv11 \
+ (-ddg1233 + gamma123*(gammado113 + 2.*gammado311) + 
     gamma113*(gammado123 + 2.*gammado312) + 
     gamma223*(gammado213 + 2.*gammado312) + 
     gamma213*(gammado223 + 2.*gammado322) + 
     3.*(gamma323*gammado313 + gamma313*gammado323))*ginv12 + 
  (-ddg1333 + gamma133*(gammado113 + 2.*gammado311) + 
     gamma233*(gammado213 + 2.*gammado312) + 
     gamma113*(gammado133 + 2.*gammado313) + 
     gamma213*(gammado233 + 2.*gammado323) + 
     3.*(gamma333*gammado313 + gamma313*gammado333))*ginv13 + 
  (-0.5*ddg2233 + gamma123*(gammado123 + 2.*gammado312) + 
     gamma223*(gammado223 + 2.*gammado322) + 3.*gamma323*gammado323)*ginv22 \
+ (-ddg2333 + gamma133*(gammado123 + 2.*gammado312) + 
     gamma123*(gammado133 + 2.*gammado313) + 
     gamma233*(gammado223 + 2.*gammado322) + 
     gamma223*(gammado233 + 2.*gammado323) + 
     3.*(gamma333*gammado323 + gamma323*gammado333))*ginv23 + 
  (-0.5*ddg3333 + gamma133*(gammado133 + 2.*gammado313) + 
     gamma233*(gammado233 + 2.*gammado323) + 3.*gamma333*gammado333)*ginv33 \
+ dG31*g13[ijk] + dG32*g23[ijk] + dG33*g33[ijk]
;

ff
=
chiguarded
;

oochipsipower
=
1/chipsipower
;

f
=
oochipsipower*log(ff)
;

psim4
=
exp(-4.*f)
;

df1
=
(dchi1*oochipsipower)/chiguarded
;

df2
=
(dchi2*oochipsipower)/chiguarded
;

df3
=
(dchi3*oochipsipower)/chiguarded
;

ddf11
=
(ddchi11*oochipsipower)/chiguarded - chipsipower*pow2(df1)
;

ddf12
=
-(chipsipower*df1*df2) + (ddchi12*oochipsipower)/chiguarded
;

ddf13
=
-(chipsipower*df1*df3) + (ddchi13*oochipsipower)/chiguarded
;

ddf22
=
(ddchi22*oochipsipower)/chiguarded - chipsipower*pow2(df2)
;

ddf23
=
-(chipsipower*df2*df3) + (ddchi23*oochipsipower)/chiguarded
;

ddf33
=
(ddchi33*oochipsipower)/chiguarded - chipsipower*pow2(df3)
;

cddf11
=
ddf11 - df1*gamma111 - df2*gamma211 - df3*gamma311
;

cddf12
=
ddf12 - df1*gamma112 - df2*gamma212 - df3*gamma312
;

cddf13
=
ddf13 - df1*gamma113 - df2*gamma213 - df3*gamma313
;

cddf22
=
ddf22 - df1*gamma122 - df2*gamma222 - df3*gamma322
;

cddf23
=
ddf23 - df1*gamma123 - df2*gamma223 - df3*gamma323
;

cddf33
=
ddf33 - df1*gamma133 - df2*gamma233 - df3*gamma333
;

trcddf
=
cddf11*ginv11 + cddf22*ginv22 + 
  2.*(cddf12*ginv12 + cddf13*ginv13 + cddf23*ginv23) + cddf33*ginv33
;

Rphi11
=
-2.*(cddf11 + trcddf*g11[ijk]) + (4. - 4.*ginv11*g11[ijk])*pow2(df1) - 
  g11[ijk]*(8.*(df1*(df2*ginv12 + df3*ginv13) + df2*df3*ginv23) + 
     4.*(ginv22*pow2(df2) + ginv33*pow2(df3)))
;

Rphi12
=
df1*df2*(4. - 8.*ginv12*g12[ijk]) - 2.*(cddf12 + trcddf*g12[ijk]) - 
  g12[ijk]*(8.*df3*(df1*ginv13 + df2*ginv23) + 
     4.*(ginv11*pow2(df1) + ginv22*pow2(df2) + ginv33*pow2(df3)))
;

Rphi13
=
df1*(4.*df3 - 8.*df2*ginv12*g13[ijk]) - 2.*(cddf13 + trcddf*g13[ijk]) - 
  g13[ijk]*(8.*df3*(df1*ginv13 + df2*ginv23) + 
     4.*(ginv11*pow2(df1) + ginv22*pow2(df2) + ginv33*pow2(df3)))
;

Rphi22
=
-2.*(cddf22 + trcddf*g22[ijk]) + (4. - 4.*ginv22*g22[ijk])*pow2(df2) - 
  g22[ijk]*(8.*(df1*(df2*ginv12 + df3*ginv13) + df2*df3*ginv23) + 
     4.*(ginv11*pow2(df1) + ginv33*pow2(df3)))
;

Rphi23
=
df2*(4.*df3 - 8.*df1*ginv12*g23[ijk]) - 2.*(cddf23 + trcddf*g23[ijk]) - 
  g23[ijk]*(8.*df3*(df1*ginv13 + df2*ginv23) + 
     4.*(ginv11*pow2(df1) + ginv22*pow2(df2) + ginv33*pow2(df3)))
;

Rphi33
=
-2.*(cddf33 + trcddf*g33[ijk]) - 
  g33[ijk]*(8.*(df1*(df2*ginv12 + df3*ginv13) + df2*df3*ginv23) + 
     4.*(ginv11*pow2(df1) + ginv22*pow2(df2))) + 
  (4. - 4.*ginv33*g33[ijk])*pow2(df3)
;

Rf11
=
R11 + Rphi11
;

Rf12
=
R12 + Rphi12
;

Rf13
=
R13 + Rphi13
;

Rf22
=
R22 + Rphi22
;

Rf23
=
R23 + Rphi23
;

Rf33
=
R33 + Rphi33
;

Rhat
=
psim4*(ginv11*Rf11 + ginv22*Rf22 + 
    2.*(ginv12*Rf12 + ginv13*Rf13 + ginv23*Rf23) + ginv33*Rf33)
;

cdda11
=
dda11 - da2*gamma211 - da3*gamma311 + 
  2.*((da2*df1 + da1*df2)*ginv12 + (da3*df1 + da1*df3)*ginv13 + 
     da2*df2*ginv22 + (da3*df2 + da2*df3)*ginv23 + da3*df3*ginv33)*g11[ijk] \
+ da1*(-4.*df1 - gamma111 + 2.*df1*ginv11*g11[ijk])
;

cdda12
=
dda12 - 2.*(da2*df1 + da1*df2) - da1*gamma112 - da2*gamma212 - 
  da3*gamma312 + 2.*(da1*df1*ginv11 + (da2*df1 + da1*df2)*ginv12 + 
     (da3*df1 + da1*df3)*ginv13 + da2*df2*ginv22 + 
     (da3*df2 + da2*df3)*ginv23 + da3*df3*ginv33)*g12[ijk]
;

cdda13
=
dda13 - 2.*(da3*df1 + da1*df3) - da1*gamma113 - da2*gamma213 - 
  da3*gamma313 + 2.*(da1*df1*ginv11 + (da2*df1 + da1*df2)*ginv12 + 
     (da3*df1 + da1*df3)*ginv13 + da2*df2*ginv22 + 
     (da3*df2 + da2*df3)*ginv23 + da3*df3*ginv33)*g13[ijk]
;

cdda22
=
dda22 - da1*gamma122 - da3*gamma322 + 
  2.*(da1*df1*ginv11 + (da2*df1 + da1*df2)*ginv12 + 
     (da3*df1 + da1*df3)*ginv13 + (da3*df2 + da2*df3)*ginv23 + 
     da3*df3*ginv33)*g22[ijk] + 
  da2*(-4.*df2 - gamma222 + 2.*df2*ginv22*g22[ijk])
;

cdda23
=
dda23 - 2.*(da3*df2 + da2*df3) - da1*gamma123 - da2*gamma223 - 
  da3*gamma323 + 2.*(da1*df1*ginv11 + (da2*df1 + da1*df2)*ginv12 + 
     (da3*df1 + da1*df3)*ginv13 + da2*df2*ginv22 + 
     (da3*df2 + da2*df3)*ginv23 + da3*df3*ginv33)*g23[ijk]
;

cdda33
=
dda33 - da1*gamma133 - da2*gamma233 + 
  2.*(da1*df1*ginv11 + (da2*df1 + da1*df2)*ginv12 + 
     (da3*df1 + da1*df3)*ginv13 + da2*df2*ginv22 + 
     (da3*df2 + da2*df3)*ginv23)*g33[ijk] + 
  da3*(-4.*df3 - gamma333 + 2.*df3*ginv33*g33[ijk])
;

dda12
=
dda12 - 2.*(da2*df1 + da1*df2) - da1*gamma112 - da2*gamma212 - 
  da3*gamma312 + 2.*(da1*df1*ginv11 + (da2*df1 + da1*df2)*ginv12 + 
     (da3*df1 + da1*df3)*ginv13 + da2*df2*ginv22 + 
     (da3*df2 + da2*df3)*ginv23 + da3*df3*ginv33)*g12[ijk]
;

dda13
=
dda13 - 2.*(da3*df1 + da1*df3) - da1*gamma113 - da2*gamma213 - 
  da3*gamma313 + 2.*(da1*df1*ginv11 + (da2*df1 + da1*df2)*ginv12 + 
     (da3*df1 + da1*df3)*ginv13 + da2*df2*ginv22 + 
     (da3*df2 + da2*df3)*ginv23 + da3*df3*ginv33)*g13[ijk]
;

dda23
=
dda23 - 2.*(da3*df2 + da2*df3) - da1*gamma123 - da2*gamma223 - 
  da3*gamma323 + 2.*(da1*df1*ginv11 + (da2*df1 + da1*df2)*ginv12 + 
     (da3*df1 + da1*df3)*ginv13 + da2*df2*ginv22 + 
     (da3*df2 + da2*df3)*ginv23 + da3*df3*ginv33)*g23[ijk]
;

trcdda
=
(cdda11*ginv11 + (cdda12 + dda12)*ginv12 + (cdda13 + dda13)*ginv13 + 
    cdda22*ginv22 + (cdda23 + dda23)*ginv23 + cdda33*ginv33)*psim4
;

AA11
=
2.*(ginv23*A12[ijk]*A13[ijk] + 
     A11[ijk]*(ginv12*A12[ijk] + ginv13*A13[ijk])) + ginv11*pow2(A11[ijk]) + 
  ginv22*pow2(A12[ijk]) + ginv33*pow2(A13[ijk])
;

AA12
=
A12[ijk]*(ginv11*A11[ijk] + ginv22*A22[ijk]) + ginv33*A13[ijk]*A23[ijk] + 
  ginv13*(A12[ijk]*A13[ijk] + A11[ijk]*A23[ijk]) + 
  ginv23*(A13[ijk]*A22[ijk] + A12[ijk]*A23[ijk]) + 
  ginv12*(A11[ijk]*A22[ijk] + pow2(A12[ijk]))
;

AA13
=
ginv22*A12[ijk]*A23[ijk] + ginv12*(A12[ijk]*A13[ijk] + A11[ijk]*A23[ijk]) + 
  A13[ijk]*(ginv11*A11[ijk] + ginv33*A33[ijk]) + 
  ginv23*(A13[ijk]*A23[ijk] + A12[ijk]*A33[ijk]) + 
  ginv13*(A11[ijk]*A33[ijk] + pow2(A13[ijk]))
;

AA21
=
A12[ijk]*(ginv11*A11[ijk] + ginv22*A22[ijk]) + ginv33*A13[ijk]*A23[ijk] + 
  ginv13*(A12[ijk]*A13[ijk] + A11[ijk]*A23[ijk]) + 
  ginv23*(A13[ijk]*A22[ijk] + A12[ijk]*A23[ijk]) + 
  ginv12*(A11[ijk]*A22[ijk] + pow2(A12[ijk]))
;

AA22
=
2.*(ginv23*A22[ijk]*A23[ijk] + 
     A12[ijk]*(ginv12*A22[ijk] + ginv13*A23[ijk])) + ginv11*pow2(A12[ijk]) + 
  ginv22*pow2(A22[ijk]) + ginv33*pow2(A23[ijk])
;

AA23
=
ginv11*A12[ijk]*A13[ijk] + ginv12*(A13[ijk]*A22[ijk] + A12[ijk]*A23[ijk]) + 
  A23[ijk]*(ginv22*A22[ijk] + ginv33*A33[ijk]) + 
  ginv13*(A13[ijk]*A23[ijk] + A12[ijk]*A33[ijk]) + 
  ginv23*(A22[ijk]*A33[ijk] + pow2(A23[ijk]))
;

AA31
=
ginv22*A12[ijk]*A23[ijk] + ginv12*(A12[ijk]*A13[ijk] + A11[ijk]*A23[ijk]) + 
  A13[ijk]*(ginv11*A11[ijk] + ginv33*A33[ijk]) + 
  ginv23*(A13[ijk]*A23[ijk] + A12[ijk]*A33[ijk]) + 
  ginv13*(A11[ijk]*A33[ijk] + pow2(A13[ijk]))
;

AA32
=
ginv11*A12[ijk]*A13[ijk] + ginv12*(A13[ijk]*A22[ijk] + A12[ijk]*A23[ijk]) + 
  A23[ijk]*(ginv22*A22[ijk] + ginv33*A33[ijk]) + 
  ginv13*(A13[ijk]*A23[ijk] + A12[ijk]*A33[ijk]) + 
  ginv23*(A22[ijk]*A33[ijk] + pow2(A23[ijk]))
;

AA33
=
2.*(ginv23*A23[ijk]*A33[ijk] + 
     A13[ijk]*(ginv12*A23[ijk] + ginv13*A33[ijk])) + ginv11*pow2(A13[ijk]) + 
  ginv22*pow2(A23[ijk]) + ginv33*pow2(A33[ijk])
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

cdA111
=
dA111 - 2.*(gamma111*A11[ijk] + gamma211*A12[ijk] + gamma311*A13[ijk])
;

cdA112
=
dA112 - gamma112*A11[ijk] - (gamma111 + gamma212)*A12[ijk] - 
  gamma312*A13[ijk] - gamma211*A22[ijk] - gamma311*A23[ijk]
;

cdA113
=
dA113 - gamma113*A11[ijk] - gamma213*A12[ijk] - 
  (gamma111 + gamma313)*A13[ijk] - gamma211*A23[ijk] - gamma311*A33[ijk]
;

cdA122
=
dA122 - 2.*(gamma112*A12[ijk] + gamma212*A22[ijk] + gamma312*A23[ijk])
;

cdA123
=
dA123 - gamma113*A12[ijk] - gamma112*A13[ijk] - gamma213*A22[ijk] - 
  (gamma212 + gamma313)*A23[ijk] - gamma312*A33[ijk]
;

cdA133
=
dA133 - 2.*(gamma113*A13[ijk] + gamma213*A23[ijk] + gamma313*A33[ijk])
;

cdA211
=
dA211 - 2.*(gamma112*A11[ijk] + gamma212*A12[ijk] + gamma312*A13[ijk])
;

cdA212
=
dA212 - gamma122*A11[ijk] - (gamma112 + gamma222)*A12[ijk] - 
  gamma322*A13[ijk] - gamma212*A22[ijk] - gamma312*A23[ijk]
;

cdA213
=
dA213 - gamma123*A11[ijk] - gamma223*A12[ijk] - 
  (gamma112 + gamma323)*A13[ijk] - gamma212*A23[ijk] - gamma312*A33[ijk]
;

cdA222
=
dA222 - 2.*(gamma122*A12[ijk] + gamma222*A22[ijk] + gamma322*A23[ijk])
;

cdA223
=
dA223 - gamma123*A12[ijk] - gamma122*A13[ijk] - gamma223*A22[ijk] - 
  (gamma222 + gamma323)*A23[ijk] - gamma322*A33[ijk]
;

cdA233
=
dA233 - 2.*(gamma123*A13[ijk] + gamma223*A23[ijk] + gamma323*A33[ijk])
;

cdA311
=
dA311 - 2.*(gamma113*A11[ijk] + gamma213*A12[ijk] + gamma313*A13[ijk])
;

cdA312
=
dA312 - gamma123*A11[ijk] - (gamma113 + gamma223)*A12[ijk] - 
  gamma323*A13[ijk] - gamma213*A22[ijk] - gamma313*A23[ijk]
;

cdA313
=
dA313 - gamma133*A11[ijk] - gamma233*A12[ijk] - 
  (gamma113 + gamma333)*A13[ijk] - gamma213*A23[ijk] - gamma313*A33[ijk]
;

cdA322
=
dA322 - 2.*(gamma123*A12[ijk] + gamma223*A22[ijk] + gamma323*A23[ijk])
;

cdA323
=
dA323 - gamma133*A12[ijk] - gamma123*A13[ijk] - gamma233*A22[ijk] - 
  (gamma223 + gamma333)*A23[ijk] - gamma323*A33[ijk]
;

cdA333
=
dA333 - 2.*(gamma133*A13[ijk] + gamma233*A23[ijk] + gamma333*A33[ijk])
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

lieA11
=
(2.*db11 - totdivbeta)*A11[ijk] + 2.*(db12*A12[ijk] + db13*A13[ijk]) + 
  dA111*beta1[ijk] + dA211*beta2[ijk] + dA311*beta3[ijk]
;

lieA12
=
db21*A11[ijk] + (db11 + db22 - totdivbeta)*A12[ijk] + db23*A13[ijk] + 
  db12*A22[ijk] + db13*A23[ijk] + dA112*beta1[ijk] + dA212*beta2[ijk] + 
  dA312*beta3[ijk]
;

lieA13
=
db31*A11[ijk] + db32*A12[ijk] + (db11 + db33 - totdivbeta)*A13[ijk] + 
  db12*A23[ijk] + db13*A33[ijk] + dA113*beta1[ijk] + dA213*beta2[ijk] + 
  dA313*beta3[ijk]
;

lieA22
=
-(totdivbeta*A22[ijk]) + 2.*(db21*A12[ijk] + db22*A22[ijk] + 
     db23*A23[ijk]) + dA122*beta1[ijk] + dA222*beta2[ijk] + dA322*beta3[ijk]
;

lieA23
=
db31*A12[ijk] + db21*A13[ijk] + db32*A22[ijk] + 
  (db22 + db33 - totdivbeta)*A23[ijk] + db23*A33[ijk] + dA123*beta1[ijk] + 
  dA223*beta2[ijk] + dA323*beta3[ijk]
;

lieA33
=
-(totdivbeta*A33[ijk]) + 2.*(db31*A13[ijk] + db32*A23[ijk] + 
     db33*A33[ijk]) + dA133*beta1[ijk] + dA233*beta2[ijk] + dA333*beta3[ijk]
;

betas
=
sdown1*beta1[ijk] + sdown2*beta2[ijk] + sdown3*beta3[ijk]
;

Dbetas
=
(db11*sdown1 + db12*sdown2 + db13*sdown3)*sup1 + 
  (db21*sdown1 + db22*sdown2 + db23*sdown3)*sup2 + 
  (db31*sdown1 + db32*sdown2 + db33*sdown3)*sup3
;

Dalpha
=
da1*sup1 + da2*sup2 + da3*sup3
;

DKhat
=
dKhat1*sup1 + dKhat2*sup2 + dKhat3*sup3
;

DK
=
dK1*sup1 + dK2*sup2 + dK3*sup3
;

DTheta
=
dTheta1*sup1 + dTheta2*sup2 + dTheta3*sup3
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

betaA1
=
qud11*beta1[ijk] + qud12*beta2[ijk] + qud13*beta3[ijk]
;

betaA2
=
qud21*beta1[ijk] + qud22*beta2[ijk] + qud23*beta3[ijk]
;

betaA3
=
qud31*beta1[ijk] + qud32*beta2[ijk] + qud33*beta3[ijk]
;

DbetaA1
=
(db11*qud11 + db12*qud12 + db13*qud13)*sup1 + 
  (db21*qud11 + db22*qud12 + db23*qud13)*sup2 + 
  (db31*qud11 + db32*qud12 + db33*qud13)*sup3
;

DbetaA2
=
(db11*qud21 + db12*qud22 + db13*qud23)*sup1 + 
  (db21*qud21 + db22*qud22 + db23*qud23)*sup2 + 
  (db31*qud21 + db32*qud22 + db33*qud23)*sup3
;

DbetaA3
=
(db11*qud31 + db12*qud32 + db13*qud33)*sup1 + 
  (db21*qud31 + db22*qud32 + db23*qud33)*sup2 + 
  (db31*qud31 + db32*qud32 + db33*qud33)*sup3
;

lienKhat
=
-((DKhat + Khat[ijk]/r)*sqrt(muL))
;

lienTheta
=
-DTheta - (kappa1*(2. + kappa2) + 1/r)*Theta[ijk]
;

lienK
=
lienKhat + 2.*lienTheta
;

rKhat
=
lienKhat*alpha[ijk] + dKhat1*beta1[ijk] + dKhat2*beta2[ijk] + 
  dKhat3*beta3[ijk]
;

rGams
=
((ddb111*quu11 + ddb221*quu22 + 
        2.*(ddb121*quu12 + ddb131*quu13 + ddb231*quu23) + ddb331*quu33)*
      sdown1 + (ddb112*quu11 + ddb222*quu22 + 2.*ddb232*quu23 + 
        ddb332*quu33)*sdown2 + 
     (ddb113*quu11 + ddb223*quu22 + 2.*(ddb133*quu13 + ddb233*quu23) + 
        ddb333*quu33)*sdown3 - 
     (ddb111*qud11 + ddb112*qud12 + ddb113*qud13 + ddb121*qud21 + 
        ddb122*qud22 + ddb123*qud23 + ddb131*qud31 + ddb132*qud32 + 
        ddb133*qud33)*sup1 - (ddb121*qud11 + ddb122*qud12 + 
        ddb123*qud13 + ddb221*qud21 + ddb222*qud22 + ddb223*qud23 + 
        ddb231*qud31 + ddb232*qud32 + ddb233*qud33)*sup2 - 
     (ddb131*qud11 + ddb132*qud12 + ddb133*qud13 + ddb231*qud21 + 
        ddb232*qud22 + ddb233*qud23 + ddb331*qud31 + ddb332*qud32 + 
        ddb333*qud33)*sup3)/chiguarded - (dG11 + dG22 + dG33)*vbetas + 
  (dG11*sdown1 + dG12*sdown2 + dG13*sdown3)*beta1[ijk] + 
  (dG21*sdown1 + dG22*sdown2 + dG23*sdown3)*beta2[ijk] + 
  (dG31*sdown1 + dG32*sdown2 + dG33*sdown3)*beta3[ijk] + 
  2.*(((ddb122*quu12 + ddb132*quu13)*sdown2 + ddb123*quu12*sdown3)/
      chiguarded + (0.33333333333333333333*
        (dTheta1*sup1 + dTheta2*sup2 + dTheta3*sup3)*alpha[ijk])/
      (chiguarded + chiguarded*vbetas) + 
     ((db11 + db22 + db33)*shiftdriver)/(vbetaA*sqrt(3.))) + 
  (1.3333333333333333333*(dKhat1*sup1 + dKhat2*sup2 + dKhat3*sup3)*
     alpha[ijk]*sqrt(muL))/(chiguarded*(vbetas + sqrt(muL)))
;

rTheta
=
lienTheta*alpha[ijk] + dTheta1*beta1[ijk] + dTheta2*beta2[ijk] + 
  dTheta3*beta3[ijk]
;

rACss
=
sup1*(2.*lieA13*sup3 + 1.3333333333333333333*chiguarded*dK1*alpha[ijk] + 
     sup2*(2.*lieA12 - cdda12*psim4 + 
        0.66666666666666666667*trcdda*g12[ijk])) + 
  sup3*(1.3333333333333333333*chiguarded*dK3*alpha[ijk] + 
     sup1*(-(dda13*psim4) + 0.66666666666666666667*trcdda*g13[ijk]) + 
     sup2*(2.*lieA23 - 0.66666666666666666667*Rhat*alpha[ijk]*g23[ijk])) + 
  0.66666666666666666667*(sup2*
      (-(Rhat*sup1*alpha[ijk]*g12[ijk]) + sup3*trcdda*g23[ijk]) + 
     ((-(dG11*sdown1) - dG12*sdown2 + dGfromgdu12*sdown2 + 
           dGfromgdu13*sdown3)*sup1 + 
        ((-dG21 + dGfromgdu21)*sdown1 + (-dG22 + dGfromgdu22)*sdown2 + 
           (-dG23 + dGfromgdu23)*sdown3)*sup2 + 
        (-(dG31*sdown1) - dG32*sdown2 + dGfromgdu32*sdown2 + 
           dGfromgdu33*sdown3)*sup3)*alpha[ijk]*pow2(chiguarded)) - 
  alpha[ijk]*(2.*(chiguarded*DTheta + AA13*sup1*sup3) + 
     0.66666666666666666667*sup3*
      (Rhat*sup1*g13[ijk] + dG33*sdown3*pow2(chiguarded))) + 
  (lieA11 - cdda11*psim4 - 2.*AA11*alpha[ijk] + 
     0.33333333333333333333*trcdda*g11[ijk])*pow2(sup1) + 
  (lieA22 - 2.*AA22*alpha[ijk] + 0.33333333333333333333*trcdda*g22[ijk])*
   pow2(sup2) - psim4*((cdda23 + dda23)*sup2*sup3 + 
     sup1*(dda12*sup2 + cdda13*sup3) + cdda22*pow2(sup2)) + 
  (lieA33 - cdda33*psim4 - 2.*AA33*alpha[ijk] + 
     0.33333333333333333333*trcdda*g33[ijk])*pow2(sup3) + 
  alpha[ijk]*(2.*((-AA31 + psim4*Rf13)*sup1*sup3 - 
        sup2*((AA12 + AA21)*sup1 + (AA23 + AA32)*sup3)) + 
     ginv11*(-2.*chiguarded*(cdA111*sup1 + cdA112*sup2) + 
        3.*dchi1*(sup1*A11[ijk] + sup2*A12[ijk]) + 
        sup3*(-2.*cdA113*chiguarded + 3.*dchi1*A13[ijk])) + 
     ginv22*(-2.*chiguarded*(cdA212*sup1 + cdA222*sup2) + 
        3.*dchi2*(sup1*A12[ijk] + sup2*A22[ijk]) + 
        sup3*(-2.*cdA223*chiguarded + 3.*dchi2*A23[ijk])) + 
     sup2*(1.3333333333333333333*chiguarded*dK2 + 
        2.*K*(sup1*A12[ijk] + sup3*A23[ijk])) + 
     ginv12*(-2.*chiguarded*((cdA112 + cdA211)*sup1 + 
           (cdA122 + cdA212)*sup2 + (cdA123 + cdA213)*sup3) + 
        3.*(sup1*(dchi2*A11[ijk] + dchi1*A12[ijk]) + 
           dchi2*(sup2*A12[ijk] + sup3*A13[ijk]) + 
           dchi1*(sup2*A22[ijk] + sup3*A23[ijk]))) + 
     ginv33*(-2.*chiguarded*(cdA313*sup1 + cdA323*sup2) + 
        3.*dchi3*(sup1*A13[ijk] + sup2*A23[ijk]) + 
        sup3*(-2.*cdA333*chiguarded + 3.*dchi3*A33[ijk])) + 
     ginv23*(-2.*chiguarded*((cdA213 + cdA312)*sup1 + 
           (cdA223 + cdA322)*sup2 + (cdA233 + cdA323)*sup3) + 
        3.*(sup1*(dchi3*A12[ijk] + dchi2*A13[ijk]) + 
           sup2*(dchi3*A22[ijk] + dchi2*A23[ijk]) + 
           sup3*(dchi3*A23[ijk] + dchi2*A33[ijk]))) + 
     ginv13*(-2.*chiguarded*((cdA113 + cdA311)*sup1 + 
           (cdA123 + cdA312)*sup2 + (cdA133 + cdA313)*sup3) + 
        3.*(dchi3*(sup1*A11[ijk] + sup2*A12[ijk]) + 
           (dchi1*sup1 + dchi3*sup3)*A13[ijk] + 
           dchi1*(sup2*A23[ijk] + sup3*A33[ijk]))) + 
     (0.66666666666666666667*((dGfromgdu11*sdown1 - dG13*sdown3)*sup1 + 
           dGfromgdu31*sdown1*sup3) + 
        kappa1*(-(Gfromg1*sdown1) - Gfromg2*sdown2 - Gfromg3*sdown3 + 
           sdown1*G1[ijk] + sdown2*G2[ijk] + sdown3*G3[ijk]))*
      pow2(chiguarded) + (psim4*Rf11 + K*A11[ijk])*pow2(sup1) + 
     psim4*(2.*sup2*(Rf12*sup1 + Rf23*sup3) + Rf22*pow2(sup2) + 
        Rf33*pow2(sup3)) + K*(2.*sup1*sup3*A13[ijk] + A22[ijk]*pow2(sup2) + 
        A33[ijk]*pow2(sup3)) + 
     0.33333333333333333333*(((dG11 - dGfromgdu11)*qud11 + 
           (dG12 - dGfromgdu12)*qud12 + (dG13 - dGfromgdu13)*qud13 + 
           (dG21 - dGfromgdu21)*qud21 + (dG22 - dGfromgdu22)*qud22 + 
           (dG23 - dGfromgdu23)*qud23 + (dG31 - dGfromgdu31)*qud31 + 
           (dG32 - dGfromgdu32)*qud32 + (dG33 - dGfromgdu33)*qud33)*
         pow2(chiguarded) - Rhat*
         (g11[ijk]*pow2(sup1) + g22[ijk]*pow2(sup2) + g33[ijk]*pow2(sup3))))
;

rACqq
=
-rACss + chiguarded*(2.*(Ainv12*lieg12 + Ainv13*lieg13 + Ainv23*lieg23) - 
     4.*(Ainv12*A12[ijk] + Ainv13*A13[ijk] + Ainv23*A23[ijk])*alpha[ijk] + 
     Ainv11*(lieg11 - 2.*A11[ijk]*alpha[ijk]) + 
     Ainv22*(lieg22 - 2.*A22[ijk]*alpha[ijk]) + 
     Ainv33*(lieg33 - 2.*A33[ijk]*alpha[ijk]))
;

rGamA1
=
-(((dG11*qud11 + dG12*qud12 + dG13*qud13)*sup1 + 
       (dG22*qud12 + dG23*qud13)*sup2 + (dG32*qud12 + dG33*qud13)*sup3 + 
       qud11*(dG21*sup2 + dG31*sup3))*vbetaA) - 
  (0.66666666666666666667*(dTheta1*quu11 + dTheta3*quu13)*alpha[ijk])/
   chiguarded + (2.3333333333333333333*
      ((ddb122*qud12 + ddb123*qud13)*quu12 + 
        (ddb131*qud11 + ddb132*qud12 + ddb133*qud13)*quu13) + 
     0.33333333333333333333*((ddb122*qud22 + ddb123*qud23 + 
           ddb131*qud31 + ddb132*qud32)*quu11 + 
        (ddb221*qud21 + ddb222*qud22 + ddb223*qud23 + ddb231*qud31 + 
           ddb232*qud32)*quu12 + 
        (ddb231*qud21 + ddb233*qud23 + ddb331*qud31 + ddb332*qud32 + 
           ddb333*qud33)*quu13) + 
     (ddb221*qud11 + ddb222*qud12 + ddb223*qud13)*quu22 + 
     2.*(ddb231*qud11 + ddb232*qud12 + ddb233*qud13)*quu23 + 
     (ddb331*qud11 + ddb332*qud12 + ddb333*qud13)*quu33 - 
     (((db11*quu11 + db21*quu12 + db31*quu13)*sdown1 + 
          (db12*quu11 + db22*quu12 + db32*quu13)*sdown2 + 
          (db13*quu11 + db23*quu12)*sdown3)*shiftdriver)/vbetaA + 
     (shiftdriver*(-(db33*quu13*sdown3) + 
          (db12*qud12 + db13*qud13)*sup1 + 
          (db22*qud12 + db23*qud13)*sup2 + 
          (db32*qud12 + db33*qud13)*sup3 + 
          qud11*(db11*sup1 + db21*sup2 + db31*sup3)))/vbetaA + 
     ((dG12*quu11 + dG22*quu12 + dG32*quu13)*sdown2 + 
        (dG13*quu11 + dG23*quu12 + dG33*quu13)*sdown3)*vbetaA + 
     quu11*(0.33333333333333333333*(ddb121*qud21 + ddb133*qud33) + 
        1.3333333333333333333*(ddb113*qud13 + ddb121*sdown1*sup2) + 
        dG11*sdown1*vbetaA) - 0.66666666666666666667*dTheta2*quu12*
      alpha[ijk] + quu12*(0.33333333333333333333*ddb233*qud33 + 
        ddb121*(2.3333333333333333333*qud11 + 
           1.3333333333333333333*sdown1*sup1) + 
        sdown1*(1.3333333333333333333*ddb231*sup3 + dG21*vbetaA) - 
        1.3333333333333333333*dKhat2*alpha[ijk]) + 
     1.3333333333333333333*((ddb131*quu13*sdown1 + 
           quu11*(ddb112*sdown2 + ddb113*sdown3))*sup1 + 
        sdown1*((ddb221*quu12 + ddb231*quu13)*sup2 + ddb331*quu13*sup3) + 
        sdown2*((ddb122*quu12 + ddb132*quu13)*sup1 + 
           (ddb122*quu11 + ddb222*quu12 + ddb232*quu13)*sup2 + 
           (ddb132*quu11 + ddb232*quu12 + ddb332*quu13)*sup3) + 
        sdown3*((ddb123*quu12 + ddb133*quu13)*sup1 + 
           (ddb123*quu11 + ddb223*quu12)*sup2 + 
           (ddb133*quu11 + ddb233*quu12 + ddb333*quu13)*sup3) + 
        quu11*(ddb112*qud12 + ddb111*(qud11 + sdown1*sup1) + 
           ddb131*sdown1*sup3 - dKhat1*alpha[ijk])) + 
     quu13*(0.33333333333333333333*ddb232*qud22 + dG31*sdown1*vbetaA + 
        1.3333333333333333333*(ddb233*sdown3*sup2 - dKhat3*alpha[ijk])))/
   chiguarded + (dG11*qud11 + dG12*qud12 + dG13*qud13)*beta1[ijk] + 
  (dG21*qud11 + dG22*qud12 + dG23*qud13)*beta2[ijk] + 
  (dG31*qud11 + dG32*qud12 + dG33*qud13)*beta3[ijk]
;

rGamA2
=
-(((dG11*qud21 + dG12*qud22 + dG13*qud23)*sup1 + 
       (dG22*qud22 + dG23*qud23)*sup2 + (dG32*qud22 + dG33*qud23)*sup3 + 
       qud21*(dG21*sup2 + dG31*sup3))*vbetaA) - 
  (0.66666666666666666667*(dTheta1*quu12 + dTheta3*quu23)*alpha[ijk])/
   chiguarded + ((ddb111*qud21 + ddb112*qud22 + ddb113*qud23)*quu11 + 
     2.*(ddb131*qud21 + ddb132*qud22 + ddb133*qud23)*quu13 + 
     2.3333333333333333333*((ddb121*qud21 + ddb122*qud22)*quu12 + 
        ddb231*qud21*quu23) + 0.33333333333333333333*
      ((ddb112*qud12 + ddb113*qud13 + ddb132*qud32 + ddb133*qud33)*
         quu12 + (ddb121*qud11 + ddb122*qud12 + ddb123*qud13 + 
           ddb231*qud31 + ddb232*qud32 + ddb233*qud33)*quu22 + 
        (ddb131*qud11 + ddb132*qud12 + ddb332*qud32 + ddb333*qud33)*quu23) \
+ (ddb332*qud22 + ddb333*qud23)*quu33 + 
     qud21*(1.3333333333333333333*ddb221*quu22 + ddb331*quu33 + 
        (db11*shiftdriver*sup1)/vbetaA) - 
     (((db11*quu12 + db21*quu22 + db31*quu23)*sdown1 + 
          (db12*quu12 + db22*quu22 + db32*quu23)*sdown2 + 
          (db13*quu12 + db23*quu22)*sdown3)*shiftdriver)/vbetaA + 
     (shiftdriver*(-(db33*quu23*sdown3) + 
          (db12*qud22 + db13*qud23)*sup1 + 
          (db21*qud21 + db22*qud22 + db23*qud23)*sup2 + 
          (db31*qud21 + db32*qud22 + db33*qud23)*sup3))/vbetaA + 
     (dG21*quu22*sdown1 + (dG22*quu22 + dG32*quu23)*sdown2 + 
        (dG23*quu22 + dG33*quu23)*sdown3 + 
        quu12*(dG12*sdown2 + dG13*sdown3))*vbetaA + 
     quu12*(2.3333333333333333333*ddb123*qud23 + 
        0.33333333333333333333*(ddb111*qud11 + ddb131*qud31) + 
        sdown1*(1.3333333333333333333*ddb121*sup2 + dG11*vbetaA)) - 
     0.66666666666666666667*dTheta2*quu22*alpha[ijk] + 
     1.3333333333333333333*((ddb132*quu23*sdown2 + ddb113*quu12*sdown3)*
         sup1 + sdown1*((ddb111*quu12 + ddb131*quu23)*sup1 + 
           ddb221*quu22*sup2 + ddb131*quu12*sup3) + 
        sdown3*((ddb123*quu22 + ddb133*quu23)*sup1 + 
           (ddb123*quu12 + ddb223*quu22)*sup2 + 
           (ddb133*quu12 + ddb233*quu22 + ddb333*quu23)*sup3) + 
        quu23*((ddb231*sdown1 + ddb232*sdown2)*sup2 + 
           ddb331*sdown1*sup3) + 
        sdown2*((ddb112*quu12 + ddb122*quu22)*sup1 + 
           (ddb232*quu22 + ddb332*quu23)*sup3 + 
           quu12*(ddb122*sup2 + ddb132*sup3)) - dKhat1*quu12*alpha[ijk] + 
        quu22*(ddb223*qud23 + ddb222*(qud22 + sdown2*sup2) + 
           sdown1*(ddb121*sup1 + ddb231*sup3) - dKhat2*alpha[ijk])) + 
     quu23*(2.3333333333333333333*(ddb232*qud22 + ddb233*qud23) + 
        0.33333333333333333333*(ddb133*qud13 + ddb331*qud31) + 
        dG31*sdown1*vbetaA + 1.3333333333333333333*
         (ddb233*sdown3*sup2 - dKhat3*alpha[ijk])))/chiguarded + 
  (dG11*qud21 + dG12*qud22 + dG13*qud23)*beta1[ijk] + 
  (dG21*qud21 + dG22*qud22 + dG23*qud23)*beta2[ijk] + 
  (dG31*qud21 + dG32*qud22 + dG33*qud23)*beta3[ijk]
;

rGamA3
=
-(((dG11*qud31 + dG12*qud32 + dG13*qud33)*sup1 + 
       (dG22*qud32 + dG23*qud33)*sup2 + (dG32*qud32 + dG33*qud33)*sup3 + 
       qud31*(dG21*sup2 + dG31*sup3))*vbetaA) - 
  (0.66666666666666666667*(dTheta1*quu13 + dTheta3*quu33)*alpha[ijk])/
   chiguarded + ((ddb111*qud31 + ddb112*qud32 + ddb113*qud33)*quu11 + 
     2.*(ddb121*qud31 + ddb122*qud32 + ddb123*qud33)*quu12 + 
     (ddb222*qud32 + ddb223*qud33)*quu22 + 
     2.3333333333333333333*((ddb132*qud32 + ddb133*qud33)*quu13 + 
        (ddb231*qud31 + ddb232*qud32 + ddb233*qud33)*quu23) + 
     0.33333333333333333333*((ddb111*qud11 + ddb112*qud12 + 
           ddb113*qud13 + ddb121*qud21 + ddb122*qud22 + ddb123*qud23)*
         quu13 + (ddb121*qud11 + ddb122*qud12 + ddb123*qud13 + 
           ddb221*qud21 + ddb222*qud22 + ddb223*qud23)*quu23 + 
        (ddb131*qud11 + ddb132*qud12 + ddb133*qud13 + ddb231*qud21 + 
           ddb232*qud22 + ddb233*qud23)*quu33) + 
     qud31*(2.3333333333333333333*ddb131*quu13 + ddb221*quu22 + 
        (db11*shiftdriver*sup1)/vbetaA) - 
     (((db11*quu13 + db21*quu23 + db31*quu33)*sdown1 + 
          (db12*quu13 + db22*quu23 + db32*quu33)*sdown2 + 
          (db13*quu13 + db23*quu23)*sdown3)*shiftdriver)/vbetaA + 
     (shiftdriver*(-(db33*quu33*sdown3) + 
          (db12*qud32 + db13*qud33)*sup1 + 
          (db21*qud31 + db22*qud32 + db23*qud33)*sup2 + 
          (db31*qud31 + db32*qud32 + db33*qud33)*sup3))/vbetaA + 
     ((dG11*quu13 + dG21*quu23 + dG31*quu33)*sdown1 + 
        (dG22*quu23 + dG32*quu33)*sdown2 + 
        (dG23*quu23 + dG33*quu33)*sdown3 + 
        quu13*(dG12*sdown2 + dG13*sdown3))*vbetaA - 
     0.66666666666666666667*dTheta2*quu23*alpha[ijk] + 
     1.3333333333333333333*((ddb132*quu33*sdown2 + ddb113*quu13*sdown3)*
         sup1 + (ddb223*quu23*sdown3 + 
           quu33*(ddb232*sdown2 + ddb233*sdown3))*sup2 + 
        sdown1*((ddb111*quu13 + ddb121*quu23)*sup1 + 
           quu13*(ddb121*sup2 + ddb131*sup3)) + 
        sdown2*((ddb112*quu13 + ddb122*quu23)*sup1 + 
           (ddb232*quu23 + ddb332*quu33)*sup3 + 
           quu13*(ddb122*sup2 + ddb132*sup3)) + 
        sdown3*((ddb123*quu23 + ddb133*quu33)*sup1 + 
           (ddb233*quu23 + ddb333*quu33)*sup3 + 
           quu13*(ddb123*sup2 + ddb133*sup3)) - dKhat1*quu13*alpha[ijk] + 
        quu23*((ddb221*sdown1 + ddb222*sdown2)*sup2 + 
           ddb231*sdown1*sup3 - dKhat2*alpha[ijk]) + 
        quu33*(ddb332*qud32 + ddb333*qud33 + 
           sdown1*(ddb131*sup1 + ddb231*sup2) + 
           ddb331*(qud31 + sdown1*sup3) - dKhat3*alpha[ijk])))/chiguarded + 
  (dG11*qud31 + dG12*qud32 + dG13*qud33)*beta1[ijk] + 
  (dG21*qud31 + dG22*qud32 + dG23*qud33)*beta2[ijk] + 
  (dG31*qud31 + dG32*qud32 + dG33*qud33)*beta3[ijk]
;

rACsA1
=
qud11*((-(cdda11*chiguarded) + lieA11)*sup1 + 
     chiguarded*(-(dda12*sup2) - dda13*sup3 + 
        0.66666666666666666667*dK1*alpha[ijk])) + 
  qud21*((-(cdda12*chiguarded) + lieA12)*sup1 + 
     chiguarded*(-(cdda22*sup2) - dda23*sup3 + 
        0.66666666666666666667*dK2*alpha[ijk])) + 
  qud31*((-(cdda13*chiguarded) + lieA13)*sup1 + 
     chiguarded*(-(cdda23*sup2) - cdda33*sup3 + 
        0.66666666666666666667*dK3*alpha[ijk])) + 
  sup2*(lieA22*qud21 + qud31*(lieA23 - 2.*AA23*alpha[ijk]) + 
     qud11*(lieA12 + chiguarded*Rf12*alpha[ijk]) - 
     0.5*dG21*qdd11*alpha[ijk]*pow2(chiguarded)) + 
  sup3*(lieA23*qud21 + qud31*(lieA33 - 2.*AA33*alpha[ijk]) + 
     qud11*(lieA13 + chiguarded*Rf13*alpha[ijk]) - 
     0.5*dG31*qdd11*alpha[ijk]*pow2(chiguarded)) + 
  alpha[ijk]*(-2.*((AA12*qud21 + AA13*qud31)*sup1 + 
        (AA21*qud11 + AA22*qud21)*sup2 + (AA31*qud11 + AA32*qud21)*sup3) + 
     K*(qud11*(sup2*A12[ijk] + sup3*A13[ijk]) + 
        qud21*(sup2*A22[ijk] + sup3*A23[ijk]) + 
        qud31*(sup1*A13[ijk] + sup2*A23[ijk] + sup3*A33[ijk])) + 
     chiguarded*(-(dTheta1*qud11) - dTheta2*qud21 - dTheta3*qud31 + 
        (qud11*Rf11 + qud21*Rf12 + qud31*Rf13)*sup1 + 
        (qud21*Rf22 + qud31*Rf23)*sup2 + (qud21*Rf23 + qud31*Rf33)*sup3 + 
        ginv11*(-(cdA111*qud11) - cdA112*qud21 - cdA113*qud31 + 
           (1.5*dchi1*(qud11*A11[ijk] + qud21*A12[ijk] + qud31*A13[ijk]))/
            chi[ijk]) + ginv22*
         (-(cdA212*qud11) - cdA222*qud21 - cdA223*qud31 + 
           (1.5*dchi2*(qud11*A12[ijk] + qud21*A22[ijk] + qud31*A23[ijk]))/
            chi[ijk]) + ginv12*
         (-((cdA112 + cdA211)*qud11) - (cdA122 + cdA212)*qud21 - 
           (cdA123 + cdA213)*qud31 + 
           (1.5*(qud11*(dchi2*A11[ijk] + dchi1*A12[ijk]) + 
                dchi2*(qud21*A12[ijk] + qud31*A13[ijk]) + 
                dchi1*(qud21*A22[ijk] + qud31*A23[ijk])))/chi[ijk]) + 
        ginv33*(-(cdA313*qud11) - cdA323*qud21 - cdA333*qud31 + 
           (1.5*dchi3*(qud11*A13[ijk] + qud21*A23[ijk] + qud31*A33[ijk]))/
            chi[ijk]) + ginv23*
         (-((cdA213 + cdA312)*qud11) - (cdA223 + cdA322)*qud21 - 
           (cdA233 + cdA323)*qud31 + 
           (1.5*(qud11*(dchi3*A12[ijk] + dchi2*A13[ijk]) + 
                qud21*(dchi3*A22[ijk] + dchi2*A23[ijk]) + 
                qud31*(dchi3*A23[ijk] + dchi2*A33[ijk])))/chi[ijk]) + 
        ginv13*(-((cdA113 + cdA311)*qud11) - (cdA123 + cdA312)*qud21 - 
           (cdA133 + cdA313)*qud31 + 
           (1.5*(dchi3*(qud11*A11[ijk] + qud21*A12[ijk]) + 
                (dchi1*qud11 + dchi3*qud31)*A13[ijk] + 
                dchi1*(qud21*A23[ijk] + qud31*A33[ijk])))/chi[ijk])) + 
     0.5*(-(Gfromg3*kappa1*qdd13) - dG11*qdd11*sup1 + 
        (dGfromgdu11*qdd11 + (-dG12 + dGfromgdu12)*qdd12 - dG13*qdd13)*
         sup1 + qdd12*((-dG22 + dGfromgdu22)*sup2 + 
           (-dG32 + dGfromgdu32)*sup3) + 
        qdd11*(dGfromgdu21*sup2 + dGfromgdu31*sup3 + kappa1*G1[ijk]) + 
        kappa1*(-(Gfromg1*qdd11) - Gfromg2*qdd12 + qdd12*G2[ijk]) + 
        qdd13*((-dG23 + dGfromgdu23)*sup2 + (-dG33 + dGfromgdu33)*sup3 + 
           kappa1*G3[ijk]))*pow2(chiguarded) + 
     sup1*(qud11*(-2.*AA11 + K*A11[ijk]) + K*qud21*A12[ijk] + 
        0.5*dGfromgdu13*qdd13*pow2(chiguarded)))
;

rACsA2
=
qud12*((-(cdda11*chiguarded) + lieA11)*sup1 + 
     chiguarded*(-(dda12*sup2) - dda13*sup3 + 
        0.66666666666666666667*dK1*alpha[ijk])) + 
  qud22*((-(cdda12*chiguarded) + lieA12)*sup1 + 
     chiguarded*(-(cdda22*sup2) - dda23*sup3 + 
        0.66666666666666666667*dK2*alpha[ijk])) + 
  qud32*((-(cdda13*chiguarded) + lieA13)*sup1 + 
     chiguarded*(-(cdda23*sup2) - cdda33*sup3 + 
        0.66666666666666666667*dK3*alpha[ijk])) + 
  sup2*(lieA22*qud22 + qud32*(lieA23 - 2.*AA23*alpha[ijk]) + 
     qud12*(lieA12 + chiguarded*Rf12*alpha[ijk]) - 
     0.5*dG21*qdd12*alpha[ijk]*pow2(chiguarded)) + 
  sup3*(lieA23*qud22 + qud32*(lieA33 - 2.*AA33*alpha[ijk]) + 
     qud12*(lieA13 + chiguarded*Rf13*alpha[ijk]) - 
     0.5*dG31*qdd12*alpha[ijk]*pow2(chiguarded)) + 
  alpha[ijk]*(-2.*((AA12*qud22 + AA13*qud32)*sup1 + 
        (AA21*qud12 + AA22*qud22)*sup2 + (AA31*qud12 + AA32*qud22)*sup3) + 
     K*(qud12*(sup2*A12[ijk] + sup3*A13[ijk]) + 
        qud22*(sup2*A22[ijk] + sup3*A23[ijk]) + 
        qud32*(sup1*A13[ijk] + sup2*A23[ijk] + sup3*A33[ijk])) + 
     chiguarded*(-(dTheta1*qud12) - dTheta2*qud22 - dTheta3*qud32 + 
        (qud12*Rf11 + qud22*Rf12 + qud32*Rf13)*sup1 + 
        (qud22*Rf22 + qud32*Rf23)*sup2 + (qud22*Rf23 + qud32*Rf33)*sup3 + 
        ginv11*(-(cdA111*qud12) - cdA112*qud22 - cdA113*qud32 + 
           (1.5*dchi1*(qud12*A11[ijk] + qud22*A12[ijk] + qud32*A13[ijk]))/
            chi[ijk]) + ginv22*
         (-(cdA212*qud12) - cdA222*qud22 - cdA223*qud32 + 
           (1.5*dchi2*(qud12*A12[ijk] + qud22*A22[ijk] + qud32*A23[ijk]))/
            chi[ijk]) + ginv12*
         (-((cdA112 + cdA211)*qud12) - (cdA122 + cdA212)*qud22 - 
           (cdA123 + cdA213)*qud32 + 
           (1.5*(qud12*(dchi2*A11[ijk] + dchi1*A12[ijk]) + 
                dchi2*(qud22*A12[ijk] + qud32*A13[ijk]) + 
                dchi1*(qud22*A22[ijk] + qud32*A23[ijk])))/chi[ijk]) + 
        ginv33*(-(cdA313*qud12) - cdA323*qud22 - cdA333*qud32 + 
           (1.5*dchi3*(qud12*A13[ijk] + qud22*A23[ijk] + qud32*A33[ijk]))/
            chi[ijk]) + ginv23*
         (-((cdA213 + cdA312)*qud12) - (cdA223 + cdA322)*qud22 - 
           (cdA233 + cdA323)*qud32 + 
           (1.5*(qud12*(dchi3*A12[ijk] + dchi2*A13[ijk]) + 
                qud22*(dchi3*A22[ijk] + dchi2*A23[ijk]) + 
                qud32*(dchi3*A23[ijk] + dchi2*A33[ijk])))/chi[ijk]) + 
        ginv13*(-((cdA113 + cdA311)*qud12) - (cdA123 + cdA312)*qud22 - 
           (cdA133 + cdA313)*qud32 + 
           (1.5*(dchi3*(qud12*A11[ijk] + qud22*A12[ijk]) + 
                (dchi1*qud12 + dchi3*qud32)*A13[ijk] + 
                dchi1*(qud22*A23[ijk] + qud32*A33[ijk])))/chi[ijk])) + 
     0.5*(-(Gfromg3*kappa1*qdd23) - dG11*qdd12*sup1 + 
        (dGfromgdu11*qdd12 + (-dG12 + dGfromgdu12)*qdd22 - dG13*qdd23)*
         sup1 + qdd22*((-dG22 + dGfromgdu22)*sup2 + 
           (-dG32 + dGfromgdu32)*sup3) + 
        qdd12*(dGfromgdu21*sup2 + dGfromgdu31*sup3 + kappa1*G1[ijk]) + 
        kappa1*(-(Gfromg1*qdd12) - Gfromg2*qdd22 + qdd22*G2[ijk]) + 
        qdd23*((-dG23 + dGfromgdu23)*sup2 + (-dG33 + dGfromgdu33)*sup3 + 
           kappa1*G3[ijk]))*pow2(chiguarded) + 
     sup1*(qud12*(-2.*AA11 + K*A11[ijk]) + K*qud22*A12[ijk] + 
        0.5*dGfromgdu13*qdd23*pow2(chiguarded)))
;

rACsA3
=
qud13*((-(cdda11*chiguarded) + lieA11)*sup1 + 
     chiguarded*(-(dda12*sup2) - dda13*sup3 + 
        0.66666666666666666667*dK1*alpha[ijk])) + 
  qud23*((-(cdda12*chiguarded) + lieA12)*sup1 + 
     chiguarded*(-(cdda22*sup2) - dda23*sup3 + 
        0.66666666666666666667*dK2*alpha[ijk])) + 
  qud33*((-(cdda13*chiguarded) + lieA13)*sup1 + 
     chiguarded*(-(cdda23*sup2) - cdda33*sup3 + 
        0.66666666666666666667*dK3*alpha[ijk])) + 
  sup2*(lieA22*qud23 + qud33*(lieA23 - 2.*AA23*alpha[ijk]) + 
     qud13*(lieA12 + chiguarded*Rf12*alpha[ijk]) - 
     0.5*dG21*qdd13*alpha[ijk]*pow2(chiguarded)) + 
  sup3*(lieA23*qud23 + qud33*(lieA33 - 2.*AA33*alpha[ijk]) + 
     qud13*(lieA13 + chiguarded*Rf13*alpha[ijk]) - 
     0.5*dG31*qdd13*alpha[ijk]*pow2(chiguarded)) + 
  alpha[ijk]*(-2.*((AA12*qud23 + AA13*qud33)*sup1 + 
        (AA21*qud13 + AA22*qud23)*sup2 + (AA31*qud13 + AA32*qud23)*sup3) + 
     K*(qud13*(sup2*A12[ijk] + sup3*A13[ijk]) + 
        qud23*(sup2*A22[ijk] + sup3*A23[ijk]) + 
        qud33*(sup1*A13[ijk] + sup2*A23[ijk] + sup3*A33[ijk])) + 
     chiguarded*(-(dTheta1*qud13) - dTheta2*qud23 - dTheta3*qud33 + 
        (qud13*Rf11 + qud23*Rf12 + qud33*Rf13)*sup1 + 
        (qud23*Rf22 + qud33*Rf23)*sup2 + (qud23*Rf23 + qud33*Rf33)*sup3 + 
        ginv11*(-(cdA111*qud13) - cdA112*qud23 - cdA113*qud33 + 
           (1.5*dchi1*(qud13*A11[ijk] + qud23*A12[ijk] + qud33*A13[ijk]))/
            chi[ijk]) + ginv22*
         (-(cdA212*qud13) - cdA222*qud23 - cdA223*qud33 + 
           (1.5*dchi2*(qud13*A12[ijk] + qud23*A22[ijk] + qud33*A23[ijk]))/
            chi[ijk]) + ginv12*
         (-((cdA112 + cdA211)*qud13) - (cdA122 + cdA212)*qud23 - 
           (cdA123 + cdA213)*qud33 + 
           (1.5*(qud13*(dchi2*A11[ijk] + dchi1*A12[ijk]) + 
                dchi2*(qud23*A12[ijk] + qud33*A13[ijk]) + 
                dchi1*(qud23*A22[ijk] + qud33*A23[ijk])))/chi[ijk]) + 
        ginv33*(-(cdA313*qud13) - cdA323*qud23 - cdA333*qud33 + 
           (1.5*dchi3*(qud13*A13[ijk] + qud23*A23[ijk] + qud33*A33[ijk]))/
            chi[ijk]) + ginv23*
         (-((cdA213 + cdA312)*qud13) - (cdA223 + cdA322)*qud23 - 
           (cdA233 + cdA323)*qud33 + 
           (1.5*(qud13*(dchi3*A12[ijk] + dchi2*A13[ijk]) + 
                qud23*(dchi3*A22[ijk] + dchi2*A23[ijk]) + 
                qud33*(dchi3*A23[ijk] + dchi2*A33[ijk])))/chi[ijk]) + 
        ginv13*(-((cdA113 + cdA311)*qud13) - (cdA123 + cdA312)*qud23 - 
           (cdA133 + cdA313)*qud33 + 
           (1.5*(dchi3*(qud13*A11[ijk] + qud23*A12[ijk]) + 
                (dchi1*qud13 + dchi3*qud33)*A13[ijk] + 
                dchi1*(qud23*A23[ijk] + qud33*A33[ijk])))/chi[ijk])) + 
     0.5*(-(Gfromg3*kappa1*qdd33) - dG11*qdd13*sup1 + 
        (dGfromgdu11*qdd13 + (-dG12 + dGfromgdu12)*qdd23 - dG13*qdd33)*
         sup1 + qdd23*((-dG22 + dGfromgdu22)*sup2 + 
           (-dG32 + dGfromgdu32)*sup3) + 
        qdd13*(dGfromgdu21*sup2 + dGfromgdu31*sup3 + kappa1*G1[ijk]) + 
        kappa1*(-(Gfromg1*qdd13) - Gfromg2*qdd23 + qdd23*G2[ijk]) + 
        qdd33*((-dG23 + dGfromgdu23)*sup2 + (-dG33 + dGfromgdu33)*sup3 + 
           kappa1*G3[ijk]))*pow2(chiguarded) + 
     sup1*(qud13*(-2.*AA11 + K*A11[ijk]) + K*qud23*A12[ijk] + 
        0.5*dGfromgdu13*qdd33*pow2(chiguarded)))
;

rACABTF11
=
(-(AA21*qPhysuudd1211) - AA31*qPhysuudd1311 + 
     (-(cdA112*qPhysuudd1211) - cdA113*qPhysuudd1311 + 
        cdA311*qPhysuudd1311 - cdA122*qPhysuudd2211 + 
        cdA212*qPhysuudd2211 + cdA213*qPhysuudd2311)*sup1 + 
     (-(cdA211*qPhysuudd1111) + cdA122*qPhysuudd1211 - 
        cdA212*qPhysuudd1211 + cdA123*qPhysuudd1311)*sup2 + 
     qPhysuudd1111*(cdA112*sup2 + cdA113*sup3 + 
        0.66666666666666666667*K*A11[ijk]) + 
     qPhysuudd1211*(cdA211*sup1 + 1.3333333333333333333*K*A12[ijk]) + 
     qPhysuudd2311*(-AA32 + 1.3333333333333333333*K*A23[ijk]) + 
     sup3*(-(cdA311*qPhysuudd1111) + 
        qPhysuudd1211*(cdA123 + (dchi3*A12[ijk])/chiguarded) + 
        (0.5*dchi3*(qPhysuudd1311*A13[ijk] + qPhysuudd2311*A23[ijk]))/
         chiguarded) + K*(1.3333333333333333333*qPhysuudd1311*A13[ijk] + 
        0.66666666666666666667*
         (qPhysuudd2211*A22[ijk] + qPhysuudd3311*A33[ijk])) + 
     ((-0.5*dchi3*qPhysuudd3311*sup1 + dchi2*qPhysuudd1311*sup2)*
         A13[ijk] + dchi1*(qPhysuudd2311*sup1 - 0.5*qPhysuudd1311*sup2)*
         A23[ijk] + 0.5*(-((dchi2*qPhysuudd1211 + dchi3*qPhysuudd1311)*
              sup1*A11[ijk]) + 
           qPhysuudd1111*(dchi2*sup2 + dchi3*sup3)*A11[ijk] - 
           (dchi3*qPhysuudd2311*sup1 + dchi1*qPhysuudd1111*sup2)*
            A12[ijk] + sup1*((dchi1*qPhysuudd1211 - 
                 dchi2*qPhysuudd2211)*A12[ijk] + 
              (dchi1*qPhysuudd1311 - dchi2*qPhysuudd2311)*A13[ijk]) + 
           (dchi1*(qPhysuudd2211*sup1 - qPhysuudd1211*sup2) + 
              dchi3*(-(qPhysuudd2311*sup2) + qPhysuudd2211*sup3))*
            A22[ijk] + sup2*((dchi2*qPhysuudd1211 - 
                 dchi3*qPhysuudd1311)*A12[ijk] + 
              (dchi2*qPhysuudd2311 - dchi3*qPhysuudd3311)*A23[ijk]) + 
           qPhysuudd3311*(dchi1*sup1 + dchi2*sup2)*A33[ijk] - 
           sup3*((dchi1*qPhysuudd1111 + dchi2*qPhysuudd1211)*A13[ijk] + 
              (dchi1*qPhysuudd1211 + dchi2*qPhysuudd2211)*A23[ijk] + 
              (dchi1*qPhysuudd1311 + dchi2*qPhysuudd2311)*A33[ijk])))/
      chiguarded)*alpha[ijk] + 
  qPhysuudd1111*(-(cdda11*chiguarded) + lieA11 - AA11*alpha[ijk]) + 
  qPhysuudd3311*(-(cdda33*chiguarded) + lieA33 - AA33*alpha[ijk] + 
     ((-cdA133 + cdA313)*sup1 + (-cdA233 + cdA323)*sup2)*alpha[ijk]) + 
  qPhysuudd2211*(-(cdda22*chiguarded) + lieA22 - AA22*alpha[ijk] + 
     (cdA223 - cdA322)*sup3*alpha[ijk]) + 
  qPhysuudd2311*(-(chiguarded*(cdda23 + dda23)) + 
     (-AA23 + cdA312*sup1 - cdA223*sup2 + cdA322*sup2 + cdA233*sup3 - 
        cdA323*sup3)*alpha[ijk] + 2.*(lieA23 - cdA123*sup1*alpha[ijk])) + 
  qPhysuudd1311*(-(chiguarded*(cdda13 + dda13)) + 
     (-AA13 + cdA312*sup2 + cdA133*sup3 - cdA313*sup3)*alpha[ijk] + 
     2.*(lieA13 - cdA213*sup2*alpha[ijk])) + 
  qPhysuudd1211*(-(chiguarded*(cdda12 + dda12)) + 
     (-AA12 + cdA213*sup3)*alpha[ijk] + 2.*(lieA12 - cdA312*sup3*alpha[ijk]))
;

rACABTF12
=
(-(AA21*qPhysuudd1212) - AA31*qPhysuudd1312 + 
     (-(cdA112*qPhysuudd1212) - cdA113*qPhysuudd1312 + 
        cdA311*qPhysuudd1312 - cdA122*qPhysuudd2212 + 
        cdA212*qPhysuudd2212 + cdA213*qPhysuudd2312)*sup1 + 
     (-(cdA211*qPhysuudd1112) + cdA122*qPhysuudd1212 - 
        cdA212*qPhysuudd1212 + cdA123*qPhysuudd1312)*sup2 + 
     qPhysuudd1112*(cdA112*sup2 + cdA113*sup3 + 
        0.66666666666666666667*K*A11[ijk]) + 
     qPhysuudd1212*(cdA211*sup1 + 1.3333333333333333333*K*A12[ijk]) + 
     qPhysuudd2312*(-AA32 + 1.3333333333333333333*K*A23[ijk]) + 
     sup3*(-(cdA311*qPhysuudd1112) + 
        qPhysuudd1212*(cdA123 + (dchi3*A12[ijk])/chiguarded) + 
        (0.5*dchi3*(qPhysuudd1312*A13[ijk] + qPhysuudd2312*A23[ijk]))/
         chiguarded) + K*(1.3333333333333333333*qPhysuudd1312*A13[ijk] + 
        0.66666666666666666667*
         (qPhysuudd2212*A22[ijk] + qPhysuudd3312*A33[ijk])) + 
     ((-0.5*dchi3*qPhysuudd3312*sup1 + dchi2*qPhysuudd1312*sup2)*
         A13[ijk] + dchi1*(qPhysuudd2312*sup1 - 0.5*qPhysuudd1312*sup2)*
         A23[ijk] + 0.5*(-((dchi2*qPhysuudd1212 + dchi3*qPhysuudd1312)*
              sup1*A11[ijk]) + 
           qPhysuudd1112*(dchi2*sup2 + dchi3*sup3)*A11[ijk] - 
           (dchi3*qPhysuudd2312*sup1 + dchi1*qPhysuudd1112*sup2)*
            A12[ijk] + sup1*((dchi1*qPhysuudd1212 - 
                 dchi2*qPhysuudd2212)*A12[ijk] + 
              (dchi1*qPhysuudd1312 - dchi2*qPhysuudd2312)*A13[ijk]) + 
           (dchi1*(qPhysuudd2212*sup1 - qPhysuudd1212*sup2) + 
              dchi3*(-(qPhysuudd2312*sup2) + qPhysuudd2212*sup3))*
            A22[ijk] + sup2*((dchi2*qPhysuudd1212 - 
                 dchi3*qPhysuudd1312)*A12[ijk] + 
              (dchi2*qPhysuudd2312 - dchi3*qPhysuudd3312)*A23[ijk]) + 
           qPhysuudd3312*(dchi1*sup1 + dchi2*sup2)*A33[ijk] - 
           sup3*((dchi1*qPhysuudd1112 + dchi2*qPhysuudd1212)*A13[ijk] + 
              (dchi1*qPhysuudd1212 + dchi2*qPhysuudd2212)*A23[ijk] + 
              (dchi1*qPhysuudd1312 + dchi2*qPhysuudd2312)*A33[ijk])))/
      chiguarded)*alpha[ijk] + 
  qPhysuudd1112*(-(cdda11*chiguarded) + lieA11 - AA11*alpha[ijk]) + 
  qPhysuudd3312*(-(cdda33*chiguarded) + lieA33 - AA33*alpha[ijk] + 
     ((-cdA133 + cdA313)*sup1 + (-cdA233 + cdA323)*sup2)*alpha[ijk]) + 
  qPhysuudd2212*(-(cdda22*chiguarded) + lieA22 - AA22*alpha[ijk] + 
     (cdA223 - cdA322)*sup3*alpha[ijk]) + 
  qPhysuudd2312*(-(chiguarded*(cdda23 + dda23)) + 
     (-AA23 + cdA312*sup1 - cdA223*sup2 + cdA322*sup2 + cdA233*sup3 - 
        cdA323*sup3)*alpha[ijk] + 2.*(lieA23 - cdA123*sup1*alpha[ijk])) + 
  qPhysuudd1312*(-(chiguarded*(cdda13 + dda13)) + 
     (-AA13 + cdA312*sup2 + cdA133*sup3 - cdA313*sup3)*alpha[ijk] + 
     2.*(lieA13 - cdA213*sup2*alpha[ijk])) + 
  qPhysuudd1212*(-(chiguarded*(cdda12 + dda12)) + 
     (-AA12 + cdA213*sup3)*alpha[ijk] + 2.*(lieA12 - cdA312*sup3*alpha[ijk]))
;

rACABTF13
=
(-(AA21*qPhysuudd1213) - AA31*qPhysuudd1313 + 
     (-(cdA112*qPhysuudd1213) - cdA113*qPhysuudd1313 + 
        cdA311*qPhysuudd1313 - cdA122*qPhysuudd2213 + 
        cdA212*qPhysuudd2213 + cdA213*qPhysuudd2313)*sup1 + 
     (-(cdA211*qPhysuudd1113) + cdA122*qPhysuudd1213 - 
        cdA212*qPhysuudd1213 + cdA123*qPhysuudd1313)*sup2 + 
     qPhysuudd1113*(cdA112*sup2 + cdA113*sup3 + 
        0.66666666666666666667*K*A11[ijk]) + 
     qPhysuudd1213*(cdA211*sup1 + 1.3333333333333333333*K*A12[ijk]) + 
     qPhysuudd2313*(-AA32 + 1.3333333333333333333*K*A23[ijk]) + 
     sup3*(-(cdA311*qPhysuudd1113) + 
        qPhysuudd1213*(cdA123 + (dchi3*A12[ijk])/chiguarded) + 
        (0.5*dchi3*(qPhysuudd1313*A13[ijk] + qPhysuudd2313*A23[ijk]))/
         chiguarded) + K*(1.3333333333333333333*qPhysuudd1313*A13[ijk] + 
        0.66666666666666666667*
         (qPhysuudd2213*A22[ijk] + qPhysuudd3313*A33[ijk])) + 
     ((-0.5*dchi3*qPhysuudd3313*sup1 + dchi2*qPhysuudd1313*sup2)*
         A13[ijk] + dchi1*(qPhysuudd2313*sup1 - 0.5*qPhysuudd1313*sup2)*
         A23[ijk] + 0.5*(-((dchi2*qPhysuudd1213 + dchi3*qPhysuudd1313)*
              sup1*A11[ijk]) + 
           qPhysuudd1113*(dchi2*sup2 + dchi3*sup3)*A11[ijk] - 
           (dchi3*qPhysuudd2313*sup1 + dchi1*qPhysuudd1113*sup2)*
            A12[ijk] + sup1*((dchi1*qPhysuudd1213 - 
                 dchi2*qPhysuudd2213)*A12[ijk] + 
              (dchi1*qPhysuudd1313 - dchi2*qPhysuudd2313)*A13[ijk]) + 
           (dchi1*(qPhysuudd2213*sup1 - qPhysuudd1213*sup2) + 
              dchi3*(-(qPhysuudd2313*sup2) + qPhysuudd2213*sup3))*
            A22[ijk] + sup2*((dchi2*qPhysuudd1213 - 
                 dchi3*qPhysuudd1313)*A12[ijk] + 
              (dchi2*qPhysuudd2313 - dchi3*qPhysuudd3313)*A23[ijk]) + 
           qPhysuudd3313*(dchi1*sup1 + dchi2*sup2)*A33[ijk] - 
           sup3*((dchi1*qPhysuudd1113 + dchi2*qPhysuudd1213)*A13[ijk] + 
              (dchi1*qPhysuudd1213 + dchi2*qPhysuudd2213)*A23[ijk] + 
              (dchi1*qPhysuudd1313 + dchi2*qPhysuudd2313)*A33[ijk])))/
      chiguarded)*alpha[ijk] + 
  qPhysuudd1113*(-(cdda11*chiguarded) + lieA11 - AA11*alpha[ijk]) + 
  qPhysuudd3313*(-(cdda33*chiguarded) + lieA33 - AA33*alpha[ijk] + 
     ((-cdA133 + cdA313)*sup1 + (-cdA233 + cdA323)*sup2)*alpha[ijk]) + 
  qPhysuudd2213*(-(cdda22*chiguarded) + lieA22 - AA22*alpha[ijk] + 
     (cdA223 - cdA322)*sup3*alpha[ijk]) + 
  qPhysuudd2313*(-(chiguarded*(cdda23 + dda23)) + 
     (-AA23 + cdA312*sup1 - cdA223*sup2 + cdA322*sup2 + cdA233*sup3 - 
        cdA323*sup3)*alpha[ijk] + 2.*(lieA23 - cdA123*sup1*alpha[ijk])) + 
  qPhysuudd1313*(-(chiguarded*(cdda13 + dda13)) + 
     (-AA13 + cdA312*sup2 + cdA133*sup3 - cdA313*sup3)*alpha[ijk] + 
     2.*(lieA13 - cdA213*sup2*alpha[ijk])) + 
  qPhysuudd1213*(-(chiguarded*(cdda12 + dda12)) + 
     (-AA12 + cdA213*sup3)*alpha[ijk] + 2.*(lieA12 - cdA312*sup3*alpha[ijk]))
;

rACABTF22
=
(-(AA21*qPhysuudd1222) - AA31*qPhysuudd1322 + 
     (-(cdA112*qPhysuudd1222) - cdA113*qPhysuudd1322 + 
        cdA311*qPhysuudd1322 - cdA122*qPhysuudd2222 + 
        cdA212*qPhysuudd2222 + cdA213*qPhysuudd2322)*sup1 + 
     (-(cdA211*qPhysuudd1122) + cdA122*qPhysuudd1222 - 
        cdA212*qPhysuudd1222 + cdA123*qPhysuudd1322)*sup2 + 
     qPhysuudd1122*(cdA112*sup2 + cdA113*sup3 + 
        0.66666666666666666667*K*A11[ijk]) + 
     qPhysuudd1222*(cdA211*sup1 + 1.3333333333333333333*K*A12[ijk]) + 
     qPhysuudd2322*(-AA32 + 1.3333333333333333333*K*A23[ijk]) + 
     sup3*(-(cdA311*qPhysuudd1122) + 
        qPhysuudd1222*(cdA123 + (dchi3*A12[ijk])/chiguarded) + 
        (0.5*dchi3*(qPhysuudd1322*A13[ijk] + qPhysuudd2322*A23[ijk]))/
         chiguarded) + K*(1.3333333333333333333*qPhysuudd1322*A13[ijk] + 
        0.66666666666666666667*
         (qPhysuudd2222*A22[ijk] + qPhysuudd3322*A33[ijk])) + 
     ((-0.5*dchi3*qPhysuudd3322*sup1 + dchi2*qPhysuudd1322*sup2)*
         A13[ijk] + dchi1*(qPhysuudd2322*sup1 - 0.5*qPhysuudd1322*sup2)*
         A23[ijk] + 0.5*(-((dchi2*qPhysuudd1222 + dchi3*qPhysuudd1322)*
              sup1*A11[ijk]) + 
           qPhysuudd1122*(dchi2*sup2 + dchi3*sup3)*A11[ijk] - 
           (dchi3*qPhysuudd2322*sup1 + dchi1*qPhysuudd1122*sup2)*
            A12[ijk] + sup1*((dchi1*qPhysuudd1222 - 
                 dchi2*qPhysuudd2222)*A12[ijk] + 
              (dchi1*qPhysuudd1322 - dchi2*qPhysuudd2322)*A13[ijk]) + 
           (dchi1*(qPhysuudd2222*sup1 - qPhysuudd1222*sup2) + 
              dchi3*(-(qPhysuudd2322*sup2) + qPhysuudd2222*sup3))*
            A22[ijk] + sup2*((dchi2*qPhysuudd1222 - 
                 dchi3*qPhysuudd1322)*A12[ijk] + 
              (dchi2*qPhysuudd2322 - dchi3*qPhysuudd3322)*A23[ijk]) + 
           qPhysuudd3322*(dchi1*sup1 + dchi2*sup2)*A33[ijk] - 
           sup3*((dchi1*qPhysuudd1122 + dchi2*qPhysuudd1222)*A13[ijk] + 
              (dchi1*qPhysuudd1222 + dchi2*qPhysuudd2222)*A23[ijk] + 
              (dchi1*qPhysuudd1322 + dchi2*qPhysuudd2322)*A33[ijk])))/
      chiguarded)*alpha[ijk] + 
  qPhysuudd1122*(-(cdda11*chiguarded) + lieA11 - AA11*alpha[ijk]) + 
  qPhysuudd3322*(-(cdda33*chiguarded) + lieA33 - AA33*alpha[ijk] + 
     ((-cdA133 + cdA313)*sup1 + (-cdA233 + cdA323)*sup2)*alpha[ijk]) + 
  qPhysuudd2222*(-(cdda22*chiguarded) + lieA22 - AA22*alpha[ijk] + 
     (cdA223 - cdA322)*sup3*alpha[ijk]) + 
  qPhysuudd2322*(-(chiguarded*(cdda23 + dda23)) + 
     (-AA23 + cdA312*sup1 - cdA223*sup2 + cdA322*sup2 + cdA233*sup3 - 
        cdA323*sup3)*alpha[ijk] + 2.*(lieA23 - cdA123*sup1*alpha[ijk])) + 
  qPhysuudd1322*(-(chiguarded*(cdda13 + dda13)) + 
     (-AA13 + cdA312*sup2 + cdA133*sup3 - cdA313*sup3)*alpha[ijk] + 
     2.*(lieA13 - cdA213*sup2*alpha[ijk])) + 
  qPhysuudd1222*(-(chiguarded*(cdda12 + dda12)) + 
     (-AA12 + cdA213*sup3)*alpha[ijk] + 2.*(lieA12 - cdA312*sup3*alpha[ijk]))
;

rACABTF23
=
(-(AA21*qPhysuudd1223) - AA31*qPhysuudd1323 + 
     (-(cdA112*qPhysuudd1223) - cdA113*qPhysuudd1323 + 
        cdA311*qPhysuudd1323 - cdA122*qPhysuudd2223 + 
        cdA212*qPhysuudd2223 + cdA213*qPhysuudd2323)*sup1 + 
     (-(cdA211*qPhysuudd1123) + cdA122*qPhysuudd1223 - 
        cdA212*qPhysuudd1223 + cdA123*qPhysuudd1323)*sup2 + 
     qPhysuudd1123*(cdA112*sup2 + cdA113*sup3 + 
        0.66666666666666666667*K*A11[ijk]) + 
     qPhysuudd1223*(cdA211*sup1 + 1.3333333333333333333*K*A12[ijk]) + 
     qPhysuudd2323*(-AA32 + 1.3333333333333333333*K*A23[ijk]) + 
     sup3*(-(cdA311*qPhysuudd1123) + 
        qPhysuudd1223*(cdA123 + (dchi3*A12[ijk])/chiguarded) + 
        (0.5*dchi3*(qPhysuudd1323*A13[ijk] + qPhysuudd2323*A23[ijk]))/
         chiguarded) + K*(1.3333333333333333333*qPhysuudd1323*A13[ijk] + 
        0.66666666666666666667*
         (qPhysuudd2223*A22[ijk] + qPhysuudd3323*A33[ijk])) + 
     ((-0.5*dchi3*qPhysuudd3323*sup1 + dchi2*qPhysuudd1323*sup2)*
         A13[ijk] + dchi1*(qPhysuudd2323*sup1 - 0.5*qPhysuudd1323*sup2)*
         A23[ijk] + 0.5*(-((dchi2*qPhysuudd1223 + dchi3*qPhysuudd1323)*
              sup1*A11[ijk]) + 
           qPhysuudd1123*(dchi2*sup2 + dchi3*sup3)*A11[ijk] - 
           (dchi3*qPhysuudd2323*sup1 + dchi1*qPhysuudd1123*sup2)*
            A12[ijk] + sup1*((dchi1*qPhysuudd1223 - 
                 dchi2*qPhysuudd2223)*A12[ijk] + 
              (dchi1*qPhysuudd1323 - dchi2*qPhysuudd2323)*A13[ijk]) + 
           (dchi1*(qPhysuudd2223*sup1 - qPhysuudd1223*sup2) + 
              dchi3*(-(qPhysuudd2323*sup2) + qPhysuudd2223*sup3))*
            A22[ijk] + sup2*((dchi2*qPhysuudd1223 - 
                 dchi3*qPhysuudd1323)*A12[ijk] + 
              (dchi2*qPhysuudd2323 - dchi3*qPhysuudd3323)*A23[ijk]) + 
           qPhysuudd3323*(dchi1*sup1 + dchi2*sup2)*A33[ijk] - 
           sup3*((dchi1*qPhysuudd1123 + dchi2*qPhysuudd1223)*A13[ijk] + 
              (dchi1*qPhysuudd1223 + dchi2*qPhysuudd2223)*A23[ijk] + 
              (dchi1*qPhysuudd1323 + dchi2*qPhysuudd2323)*A33[ijk])))/
      chiguarded)*alpha[ijk] + 
  qPhysuudd1123*(-(cdda11*chiguarded) + lieA11 - AA11*alpha[ijk]) + 
  qPhysuudd3323*(-(cdda33*chiguarded) + lieA33 - AA33*alpha[ijk] + 
     ((-cdA133 + cdA313)*sup1 + (-cdA233 + cdA323)*sup2)*alpha[ijk]) + 
  qPhysuudd2223*(-(cdda22*chiguarded) + lieA22 - AA22*alpha[ijk] + 
     (cdA223 - cdA322)*sup3*alpha[ijk]) + 
  qPhysuudd2323*(-(chiguarded*(cdda23 + dda23)) + 
     (-AA23 + cdA312*sup1 - cdA223*sup2 + cdA322*sup2 + cdA233*sup3 - 
        cdA323*sup3)*alpha[ijk] + 2.*(lieA23 - cdA123*sup1*alpha[ijk])) + 
  qPhysuudd1323*(-(chiguarded*(cdda13 + dda13)) + 
     (-AA13 + cdA312*sup2 + cdA133*sup3 - cdA313*sup3)*alpha[ijk] + 
     2.*(lieA13 - cdA213*sup2*alpha[ijk])) + 
  qPhysuudd1223*(-(chiguarded*(cdda12 + dda12)) + 
     (-AA12 + cdA213*sup3)*alpha[ijk] + 2.*(lieA12 - cdA312*sup3*alpha[ijk]))
;

rACABTF33
=
(-(AA21*qPhysuudd1233) - AA31*qPhysuudd1333 + 
     (-(cdA112*qPhysuudd1233) - cdA113*qPhysuudd1333 + 
        cdA311*qPhysuudd1333 - cdA122*qPhysuudd2233 + 
        cdA212*qPhysuudd2233 + cdA213*qPhysuudd2333)*sup1 + 
     (-(cdA211*qPhysuudd1133) + cdA122*qPhysuudd1233 - 
        cdA212*qPhysuudd1233 + cdA123*qPhysuudd1333)*sup2 + 
     qPhysuudd1133*(cdA112*sup2 + cdA113*sup3 + 
        0.66666666666666666667*K*A11[ijk]) + 
     qPhysuudd1233*(cdA211*sup1 + 1.3333333333333333333*K*A12[ijk]) + 
     qPhysuudd2333*(-AA32 + 1.3333333333333333333*K*A23[ijk]) + 
     sup3*(-(cdA311*qPhysuudd1133) + 
        qPhysuudd1233*(cdA123 + (dchi3*A12[ijk])/chiguarded) + 
        (0.5*dchi3*(qPhysuudd1333*A13[ijk] + qPhysuudd2333*A23[ijk]))/
         chiguarded) + K*(1.3333333333333333333*qPhysuudd1333*A13[ijk] + 
        0.66666666666666666667*
         (qPhysuudd2233*A22[ijk] + qPhysuudd3333*A33[ijk])) + 
     ((-0.5*dchi3*qPhysuudd3333*sup1 + dchi2*qPhysuudd1333*sup2)*
         A13[ijk] + dchi1*(qPhysuudd2333*sup1 - 0.5*qPhysuudd1333*sup2)*
         A23[ijk] + 0.5*(-((dchi2*qPhysuudd1233 + dchi3*qPhysuudd1333)*
              sup1*A11[ijk]) + 
           qPhysuudd1133*(dchi2*sup2 + dchi3*sup3)*A11[ijk] - 
           (dchi3*qPhysuudd2333*sup1 + dchi1*qPhysuudd1133*sup2)*
            A12[ijk] + sup1*((dchi1*qPhysuudd1233 - 
                 dchi2*qPhysuudd2233)*A12[ijk] + 
              (dchi1*qPhysuudd1333 - dchi2*qPhysuudd2333)*A13[ijk]) + 
           (dchi1*(qPhysuudd2233*sup1 - qPhysuudd1233*sup2) + 
              dchi3*(-(qPhysuudd2333*sup2) + qPhysuudd2233*sup3))*
            A22[ijk] + sup2*((dchi2*qPhysuudd1233 - 
                 dchi3*qPhysuudd1333)*A12[ijk] + 
              (dchi2*qPhysuudd2333 - dchi3*qPhysuudd3333)*A23[ijk]) + 
           qPhysuudd3333*(dchi1*sup1 + dchi2*sup2)*A33[ijk] - 
           sup3*((dchi1*qPhysuudd1133 + dchi2*qPhysuudd1233)*A13[ijk] + 
              (dchi1*qPhysuudd1233 + dchi2*qPhysuudd2233)*A23[ijk] + 
              (dchi1*qPhysuudd1333 + dchi2*qPhysuudd2333)*A33[ijk])))/
      chiguarded)*alpha[ijk] + 
  qPhysuudd1133*(-(cdda11*chiguarded) + lieA11 - AA11*alpha[ijk]) + 
  qPhysuudd3333*(-(cdda33*chiguarded) + lieA33 - AA33*alpha[ijk] + 
     ((-cdA133 + cdA313)*sup1 + (-cdA233 + cdA323)*sup2)*alpha[ijk]) + 
  qPhysuudd2233*(-(cdda22*chiguarded) + lieA22 - AA22*alpha[ijk] + 
     (cdA223 - cdA322)*sup3*alpha[ijk]) + 
  qPhysuudd2333*(-(chiguarded*(cdda23 + dda23)) + 
     (-AA23 + cdA312*sup1 - cdA223*sup2 + cdA322*sup2 + cdA233*sup3 - 
        cdA323*sup3)*alpha[ijk] + 2.*(lieA23 - cdA123*sup1*alpha[ijk])) + 
  qPhysuudd1333*(-(chiguarded*(cdda13 + dda13)) + 
     (-AA13 + cdA312*sup2 + cdA133*sup3 - cdA313*sup3)*alpha[ijk] + 
     2.*(lieA13 - cdA213*sup2*alpha[ijk])) + 
  qPhysuudd1233*(-(chiguarded*(cdda12 + dda12)) + 
     (-AA12 + cdA213*sup3)*alpha[ijk] + 2.*(lieA12 - cdA312*sup3*alpha[ijk]))
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


if (CheckForNANandINF(15,                                                     
    detginv, sdown1,sdown2,sdown3, sup1,sup2,sup3,
    qdd11,qdd12,qdd13,qdd22,qdd23,qdd33,
    muL,lienKhat)) {
    printf("problem with divisions in Z4_boundary.m\n");
    printf("x=%2.5e, y=%2.5e, z=%2.5e, r=%2.5e\n",xp[ijk],yp[ijk],zp[ijk],r);
  }
if (CheckForNANandINF(11,                                          
    rA11,rA12,rA13,rA22,rA23,rA33, 
    rG1,rG2,rG3,rKhat,rTheta)) {
    printf("nans in RHS in z4_boundary_box.m\n");
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

/* z4_boundary_box.c */
/* nvars = 78, nauxs = 435, n* = 7097,  n/ = 105,  n+ = 6297, n = 13499, O = 0 */
