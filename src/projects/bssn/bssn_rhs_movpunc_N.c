/* bssn_rhs_movpunc_N.c */
/* Copyright (C) 1998 Bernd Bruegmann, 6.3.2017 */
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
#define Tanh(x)    tanh(x)
#define Sech(x)    (1/cosh(x))
#define Cal(x,y,z) ((x)?(y):(z))
#define PROBLEM    printf("  %d pts away from boundary\n",boundaryaway(6)); \
                   printf("    %e %e %e \n    detg=%e \n", \
                   Ptr(level,"x")[ijk],Ptr(level,"y")[ijk],Ptr(level,"z")[ijk],detg); \
                   if (ijkinsidefinerlevel(box,ijk)==0) printf("    point is NOT inside finer box NOR in some symmetry area\n"); \
                   else printf("    point is inside finer box/ in symmetry \n"); \
                   printf("    gtilde = %e %e %e %e %e %e \n"  ,g11[ijk],g12[ijk],g13[ijk],g22[ijk], g23[ijk],g33[ijk]);\



void bssn_rhs_movpunc_N(tVarList *unew, tVarList *upre, double c, tVarList *ucur)
{

tL *level = ucur->level;

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

int index_bssn_eta = Ind("bssn_eta");
double *shiftDr = level->v[index_bssn_eta + 0];

const int addlinear           = (c != 0.0l);
const int order_centered      = Geti("order_centered");
const int order_advection     = Geti("order_advection");
const int advectionlopsided   = Geti("advection_lopsided");
const int advectionlopsided6  = Geti("advection_lopsided6");
const int advectionlopsided8  = Geti("advection_lopsided8");
const int advectionlopsided10 = Geti("advection_lopsided10");
const int order_dissipation = Geti("order_dissipation");
double dissfactor = get_dissipation_factor(level);
const double forceKzerofactor = Getv("bssn_forceKzero", "no");
const int subtractA      = Getv("bssn_subtractA", "yes");
const int normalizedetg  = Getv("bssn_normalizedetg", "yes");
const int constantlapse  = Getv("bssn_lapse", "constant");
const int oploglapse     = Getv("bssn_lapse", "1+log");
const int oploglapse2    = Getv("bssn_lapse", "1+log2");
const int oplogwithshift = Getv("bssn_lapse", "withshift");
const int harmoniclapse  = Getv("bssn_lapse", "harmonic");
const double gamma0factor    = Getv("bssn_shift", "gamma0");
const double gamma2factor    = Getv("bssn_shift", "gamma2");
const double withGadv        = Getv("bssn_shift", "withGadv");
const double withShiftadv    = Getv("bssn_shift", "withShiftadv");
const double withBadv        = Getv("bssn_shift", "withBadv");
const double withB           = !Getv("bssn_shift", "withoutB");
const double lapseharmonicf  = Getd("bssn_lapseharmonicf");
const double shiftalphapower = Getd("bssn_shiftalphapower");
const double shiftgammacoeff = Getd("bssn_shiftgammacoeff");
const double shiftGalphapower= Getd("bssn_shiftgammaalphapower");
const double shiftdriver     = Getd("bssn_shiftdriver");
const double chiDivFloor     = Getd("bssn_chi_div_floor");
const double chipsipower     = Getd("bssn_chi_psipower");
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
const int setRHSto0 = Getv("bssnRHSto0","yes");




int_advectionN_stencil;
bampi_openmp_start
int_advectionN_vars;


double dx = level->dx;
double dy = level->dy;
double dz = level->dz;
double oo2dx = 1/(2*dx), oo2dy = 1/(2*dy), oo2dz = 1/(2*dz);
double oodx2 = 1/(dx*dx), oody2 = 1/(dy*dy), oodz2 = 1/(dz*dz);
double oo4dxdy = 1/(4*dx*dy);
double oo4dxdz = 1/(4*dx*dz);
double oo4dydz = 1/(4*dy*dz);

double AA = 0.;
double AA11 = 0.;
double AA12 = 0.;
double AA13 = 0.;
double AA22 = 0.;
double AA23 = 0.;
double AA33 = 0.;
double absdpsim2 = 0.;
double advB1 = 0.;
double advB2 = 0.;
double advB3 = 0.;
double advbeta1 = 0.;
double advbeta2 = 0.;
double advbeta3 = 0.;
double advG1 = 0.;
double advG2 = 0.;
double advG3 = 0.;
double Ainv11 = 0.;
double Ainv12 = 0.;
double Ainv13 = 0.;
double Ainv22 = 0.;
double Ainv23 = 0.;
double Ainv33 = 0.;
double aux = 0.;
double betaF = 0.;
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
double daSST1 = 0.;
double dASST111 = 0.;
double dASST112 = 0.;
double dASST113 = 0.;
double dASST122 = 0.;
double dASST123 = 0.;
double dASST133 = 0.;
double daSST2 = 0.;
double dASST211 = 0.;
double dASST212 = 0.;
double dASST213 = 0.;
double dASST222 = 0.;
double dASST223 = 0.;
double dASST233 = 0.;
double daSST3 = 0.;
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
double dchi1 = 0.;
double dchi2 = 0.;
double dchi3 = 0.;
double dchiSST1 = 0.;
double dchiSST2 = 0.;
double dchiSST3 = 0.;
double dda11 = 0.;
double dda12 = 0.;
double dda13 = 0.;
double dda22 = 0.;
double dda23 = 0.;
double dda33 = 0.;
double ddaSST11 = 0.;
double ddaSST12 = 0.;
double ddaSST13 = 0.;
double ddaSST22 = 0.;
double ddaSST23 = 0.;
double ddaSST33 = 0.;
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
double ddbSST111 = 0.;
double ddbSST112 = 0.;
double ddbSST113 = 0.;
double ddbSST121 = 0.;
double ddbSST122 = 0.;
double ddbSST123 = 0.;
double ddbSST131 = 0.;
double ddbSST132 = 0.;
double ddbSST133 = 0.;
double ddbSST221 = 0.;
double ddbSST222 = 0.;
double ddbSST223 = 0.;
double ddbSST231 = 0.;
double ddbSST232 = 0.;
double ddbSST233 = 0.;
double ddbSST331 = 0.;
double ddbSST332 = 0.;
double ddbSST333 = 0.;
double ddchi11 = 0.;
double ddchi12 = 0.;
double ddchi13 = 0.;
double ddchi22 = 0.;
double ddchi23 = 0.;
double ddchi33 = 0.;
double ddchiSST11 = 0.;
double ddchiSST12 = 0.;
double ddchiSST13 = 0.;
double ddchiSST22 = 0.;
double ddchiSST23 = 0.;
double ddchiSST33 = 0.;
double ddf11 = 0.;
double ddf12 = 0.;
double ddf13 = 0.;
double ddf22 = 0.;
double ddf23 = 0.;
double ddf33 = 0.;
double ddRdr = 0.;
double deldelg1111 = 0.;
double deldelg1112 = 0.;
double deldelg1113 = 0.;
double deldelg1122 = 0.;
double deldelg1123 = 0.;
double deldelg1133 = 0.;
double deldelg1211 = 0.;
double deldelg1212 = 0.;
double deldelg1213 = 0.;
double deldelg1222 = 0.;
double deldelg1223 = 0.;
double deldelg1233 = 0.;
double deldelg1311 = 0.;
double deldelg1312 = 0.;
double deldelg1313 = 0.;
double deldelg1322 = 0.;
double deldelg1323 = 0.;
double deldelg1333 = 0.;
double deldelg2211 = 0.;
double deldelg2212 = 0.;
double deldelg2213 = 0.;
double deldelg2222 = 0.;
double deldelg2223 = 0.;
double deldelg2233 = 0.;
double deldelg2311 = 0.;
double deldelg2312 = 0.;
double deldelg2313 = 0.;
double deldelg2322 = 0.;
double deldelg2323 = 0.;
double deldelg2333 = 0.;
double deldelg3311 = 0.;
double deldelg3312 = 0.;
double deldelg3313 = 0.;
double deldelg3322 = 0.;
double deldelg3323 = 0.;
double deldelg3333 = 0.;
double deldelgSST1111 = 0.;
double deldelgSST1112 = 0.;
double deldelgSST1113 = 0.;
double deldelgSST1122 = 0.;
double deldelgSST1123 = 0.;
double deldelgSST1133 = 0.;
double deldelgSST1211 = 0.;
double deldelgSST1212 = 0.;
double deldelgSST1213 = 0.;
double deldelgSST1222 = 0.;
double deldelgSST1223 = 0.;
double deldelgSST1233 = 0.;
double deldelgSST1311 = 0.;
double deldelgSST1312 = 0.;
double deldelgSST1313 = 0.;
double deldelgSST1322 = 0.;
double deldelgSST1323 = 0.;
double deldelgSST1333 = 0.;
double deldelgSST2211 = 0.;
double deldelgSST2212 = 0.;
double deldelgSST2213 = 0.;
double deldelgSST2222 = 0.;
double deldelgSST2223 = 0.;
double deldelgSST2233 = 0.;
double deldelgSST2311 = 0.;
double deldelgSST2312 = 0.;
double deldelgSST2313 = 0.;
double deldelgSST2322 = 0.;
double deldelgSST2323 = 0.;
double deldelgSST2333 = 0.;
double deldelgSST3311 = 0.;
double deldelgSST3312 = 0.;
double deldelgSST3313 = 0.;
double deldelgSST3322 = 0.;
double deldelgSST3323 = 0.;
double deldelgSST3333 = 0.;
double delG11 = 0.;
double delg111 = 0.;
double delg112 = 0.;
double delg113 = 0.;
double delG12 = 0.;
double delg122 = 0.;
double delg123 = 0.;
double delG13 = 0.;
double delg133 = 0.;
double delG21 = 0.;
double delg211 = 0.;
double delg212 = 0.;
double delg213 = 0.;
double delG22 = 0.;
double delg222 = 0.;
double delg223 = 0.;
double delG23 = 0.;
double delg233 = 0.;
double delG31 = 0.;
double delg311 = 0.;
double delg312 = 0.;
double delg313 = 0.;
double delG32 = 0.;
double delg322 = 0.;
double delg323 = 0.;
double delG33 = 0.;
double delg333 = 0.;
double delGSST11 = 0.;
double delgSST111 = 0.;
double delgSST112 = 0.;
double delgSST113 = 0.;
double delGSST12 = 0.;
double delgSST122 = 0.;
double delgSST123 = 0.;
double delGSST13 = 0.;
double delgSST133 = 0.;
double delGSST21 = 0.;
double delgSST211 = 0.;
double delgSST212 = 0.;
double delgSST213 = 0.;
double delGSST22 = 0.;
double delgSST222 = 0.;
double delgSST223 = 0.;
double delGSST23 = 0.;
double delgSST233 = 0.;
double delGSST31 = 0.;
double delgSST311 = 0.;
double delgSST312 = 0.;
double delgSST313 = 0.;
double delGSST32 = 0.;
double delgSST322 = 0.;
double delgSST323 = 0.;
double delGSST33 = 0.;
double delgSST333 = 0.;
double detg = 0.;
double detginv = 0.;
double detnginv = 0.;
double df1 = 0.;
double df2 = 0.;
double df3 = 0.;
double divAinv1 = 0.;
double divAinv2 = 0.;
double divAinv3 = 0.;
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
double dK1 = 0.;
double dK2 = 0.;
double dK3 = 0.;
double dKSST1 = 0.;
double dKSST2 = 0.;
double dKSST3 = 0.;
double dpsim21 = 0.;
double dpsim22 = 0.;
double dpsim23 = 0.;
double dRdr = 0.;
double f = 0.;
double falpha = 0.;
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
double Gfromg1 = 0.;
double Gfromg2 = 0.;
double Gfromg3 = 0.;
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
double lieA11 = 0.;
double lieA12 = 0.;
double lieA13 = 0.;
double lieA22 = 0.;
double lieA23 = 0.;
double lieA33 = 0.;
double liealpha = 0.;
double liechi = 0.;
double lieg11 = 0.;
double lieg12 = 0.;
double lieg13 = 0.;
double lieg22 = 0.;
double lieg23 = 0.;
double lieg33 = 0.;
double lieK = 0.;
double oochipsipower = 0.;
double ootddivbeta1 = 0.;
double ootddivbeta2 = 0.;
double ootddivbeta3 = 0.;
double pseudolieG1 = 0.;
double pseudolieG2 = 0.;
double pseudolieG3 = 0.;
double psim2 = 0.;
double psim4 = 0.;
double R = 0.;
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
double ralpha = 0.;
double rB1 = 0.;
double rB2 = 0.;
double rB3 = 0.;
double rbeta1 = 0.;
double rbeta2 = 0.;
double rbeta3 = 0.;
double rchi = 0.;
double rG1 = 0.;
double rg11 = 0.;
double rg12 = 0.;
double rg13 = 0.;
double rG2 = 0.;
double rg22 = 0.;
double rg23 = 0.;
double rG3 = 0.;
double rg33 = 0.;
double rK = 0.;
double Rphi11 = 0.;
double Rphi12 = 0.;
double Rphi13 = 0.;
double Rphi22 = 0.;
double Rphi23 = 0.;
double Rphi33 = 0.;
double totdivbeta = 0.;
double traceA = 0.;
double trcdda = 0.;
double trcddf = 0.;
double w1 = 0.;
double w2 = 0.;
double w3 = 0.;



forinnerpoints_ijk_openmp(level) {



 if (alpha[ijk]<=0.) { printf("ATTENTION:  alpha<=0 at %d\n",ijk); 


 printf("  (%e %e  %e %e)\n",alpha[ijk],palpha[ijk],K[ijk],pK[ijk]); 


 printf("  x=%e y=%e z=%e\n",Ptr(level,"x")[ijk],Ptr(level,"y")[ijk],Ptr(level,"z")[ijk]);  


 if (ijkinsidefinerlevel(box,ijk)==0) printf("  point is NOT inside finer box NOR in some symmetry area\n");  


 else printf("  point is inside finer box/ in symmetry \n");}  


if (order_centered == 2 || boundaryNaway(1)) { 

da1
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

deldelg1111
=
oodx2*(-2.*g11[ijk] + g11[-di + ijk] + g11[di + ijk])
;

deldelg1112
=
oodx2*(-2.*g12[ijk] + g12[-di + ijk] + g12[di + ijk])
;

deldelg1113
=
oodx2*(-2.*g13[ijk] + g13[-di + ijk] + g13[di + ijk])
;

deldelg1122
=
oodx2*(-2.*g22[ijk] + g22[-di + ijk] + g22[di + ijk])
;

deldelg1123
=
oodx2*(-2.*g23[ijk] + g23[-di + ijk] + g23[di + ijk])
;

deldelg1133
=
oodx2*(-2.*g33[ijk] + g33[-di + ijk] + g33[di + ijk])
;

deldelg1211
=
oo4dxdy*(g11[-di - dj + ijk] - g11[di - dj + ijk] - g11[-di + dj + ijk] + 
    g11[di + dj + ijk])
;

deldelg1212
=
oo4dxdy*(g12[-di - dj + ijk] - g12[di - dj + ijk] - g12[-di + dj + ijk] + 
    g12[di + dj + ijk])
;

deldelg1213
=
oo4dxdy*(g13[-di - dj + ijk] - g13[di - dj + ijk] - g13[-di + dj + ijk] + 
    g13[di + dj + ijk])
;

deldelg1222
=
oo4dxdy*(g22[-di - dj + ijk] - g22[di - dj + ijk] - g22[-di + dj + ijk] + 
    g22[di + dj + ijk])
;

deldelg1223
=
oo4dxdy*(g23[-di - dj + ijk] - g23[di - dj + ijk] - g23[-di + dj + ijk] + 
    g23[di + dj + ijk])
;

deldelg1233
=
oo4dxdy*(g33[-di - dj + ijk] - g33[di - dj + ijk] - g33[-di + dj + ijk] + 
    g33[di + dj + ijk])
;

deldelg1311
=
oo4dxdz*(g11[-di - dk + ijk] - g11[di - dk + ijk] - g11[-di + dk + ijk] + 
    g11[di + dk + ijk])
;

deldelg1312
=
oo4dxdz*(g12[-di - dk + ijk] - g12[di - dk + ijk] - g12[-di + dk + ijk] + 
    g12[di + dk + ijk])
;

deldelg1313
=
oo4dxdz*(g13[-di - dk + ijk] - g13[di - dk + ijk] - g13[-di + dk + ijk] + 
    g13[di + dk + ijk])
;

deldelg1322
=
oo4dxdz*(g22[-di - dk + ijk] - g22[di - dk + ijk] - g22[-di + dk + ijk] + 
    g22[di + dk + ijk])
;

deldelg1323
=
oo4dxdz*(g23[-di - dk + ijk] - g23[di - dk + ijk] - g23[-di + dk + ijk] + 
    g23[di + dk + ijk])
;

deldelg1333
=
oo4dxdz*(g33[-di - dk + ijk] - g33[di - dk + ijk] - g33[-di + dk + ijk] + 
    g33[di + dk + ijk])
;

deldelg2211
=
oody2*(-2.*g11[ijk] + g11[-dj + ijk] + g11[dj + ijk])
;

deldelg2212
=
oody2*(-2.*g12[ijk] + g12[-dj + ijk] + g12[dj + ijk])
;

deldelg2213
=
oody2*(-2.*g13[ijk] + g13[-dj + ijk] + g13[dj + ijk])
;

deldelg2222
=
oody2*(-2.*g22[ijk] + g22[-dj + ijk] + g22[dj + ijk])
;

deldelg2223
=
oody2*(-2.*g23[ijk] + g23[-dj + ijk] + g23[dj + ijk])
;

deldelg2233
=
oody2*(-2.*g33[ijk] + g33[-dj + ijk] + g33[dj + ijk])
;

deldelg2311
=
oo4dydz*(g11[-dj - dk + ijk] - g11[dj - dk + ijk] - g11[-dj + dk + ijk] + 
    g11[dj + dk + ijk])
;

deldelg2312
=
oo4dydz*(g12[-dj - dk + ijk] - g12[dj - dk + ijk] - g12[-dj + dk + ijk] + 
    g12[dj + dk + ijk])
;

deldelg2313
=
oo4dydz*(g13[-dj - dk + ijk] - g13[dj - dk + ijk] - g13[-dj + dk + ijk] + 
    g13[dj + dk + ijk])
;

deldelg2322
=
oo4dydz*(g22[-dj - dk + ijk] - g22[dj - dk + ijk] - g22[-dj + dk + ijk] + 
    g22[dj + dk + ijk])
;

deldelg2323
=
oo4dydz*(g23[-dj - dk + ijk] - g23[dj - dk + ijk] - g23[-dj + dk + ijk] + 
    g23[dj + dk + ijk])
;

deldelg2333
=
oo4dydz*(g33[-dj - dk + ijk] - g33[dj - dk + ijk] - g33[-dj + dk + ijk] + 
    g33[dj + dk + ijk])
;

deldelg3311
=
oodz2*(-2.*g11[ijk] + g11[-dk + ijk] + g11[dk + ijk])
;

deldelg3312
=
oodz2*(-2.*g12[ijk] + g12[-dk + ijk] + g12[dk + ijk])
;

deldelg3313
=
oodz2*(-2.*g13[ijk] + g13[-dk + ijk] + g13[dk + ijk])
;

deldelg3322
=
oodz2*(-2.*g22[ijk] + g22[-dk + ijk] + g22[dk + ijk])
;

deldelg3323
=
oodz2*(-2.*g23[ijk] + g23[-dk + ijk] + g23[dk + ijk])
;

deldelg3333
=
oodz2*(-2.*g33[ijk] + g33[-dk + ijk] + g33[dk + ijk])
;

delG11
=
oo2dx*(-G1[-di + ijk] + G1[di + ijk])
;

delG12
=
oo2dx*(-G2[-di + ijk] + G2[di + ijk])
;

delG13
=
oo2dx*(-G3[-di + ijk] + G3[di + ijk])
;

delG21
=
oo2dy*(-G1[-dj + ijk] + G1[dj + ijk])
;

delG22
=
oo2dy*(-G2[-dj + ijk] + G2[dj + ijk])
;

delG23
=
oo2dy*(-G3[-dj + ijk] + G3[dj + ijk])
;

delG31
=
oo2dz*(-G1[-dk + ijk] + G1[dk + ijk])
;

delG32
=
oo2dz*(-G2[-dk + ijk] + G2[dk + ijk])
;

delG33
=
oo2dz*(-G3[-dk + ijk] + G3[dk + ijk])
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


} else if (order_centered == 4 || boundaryNaway(2)) { 


#ifdef REDUCEORDERTO2                           
  errorexit("CompilerFlag reduced order to 2");
 #else
da1
=
0.16666666666666666667*oo2dx*(alpha[-2*di + ijk] + 
    8.*(-alpha[-di + ijk] + alpha[di + ijk]) - alpha[2*di + ijk])
;

da2
=
0.16666666666666666667*oo2dy*(alpha[-2*dj + ijk] + 
    8.*(-alpha[-dj + ijk] + alpha[dj + ijk]) - alpha[2*dj + ijk])
;

da3
=
0.16666666666666666667*oo2dz*(alpha[-2*dk + ijk] + 
    8.*(-alpha[-dk + ijk] + alpha[dk + ijk]) - alpha[2*dk + ijk])
;

dda11
=
0.083333333333333333333*oodx2*(-30.*alpha[ijk] - alpha[-2*di + ijk] + 
    16.*(alpha[-di + ijk] + alpha[di + ijk]) - alpha[2*di + ijk])
;

dda12
=
0.027777777777777777778*oo2dx*oo2dy*
  (-alpha[2*(di - dj) + ijk] + 
    64.*(alpha[-di - dj + ijk] - alpha[di - dj + ijk] - 
       alpha[-di + dj + ijk] + alpha[di + dj + ijk]) + 
    8.*(-alpha[-di - 2*dj + ijk] + alpha[di - 2*dj + ijk] - 
       alpha[-2*di - dj + ijk] + alpha[2*di - dj + ijk] + 
       alpha[-2*di + dj + ijk] - alpha[2*di + dj + ijk] + 
       alpha[-di + 2*dj + ijk] - alpha[di + 2*dj + ijk]) - 
    alpha[2*(-di + dj) + ijk] + alpha[-2*(di + dj) + ijk] + 
    alpha[2*(di + dj) + ijk])
;

dda13
=
0.027777777777777777778*oo2dx*oo2dz*
  (-alpha[2*(di - dk) + ijk] + 
    64.*(alpha[-di - dk + ijk] - alpha[di - dk + ijk] - 
       alpha[-di + dk + ijk] + alpha[di + dk + ijk]) + 
    8.*(-alpha[-di - 2*dk + ijk] + alpha[di - 2*dk + ijk] - 
       alpha[-2*di - dk + ijk] + alpha[2*di - dk + ijk] + 
       alpha[-2*di + dk + ijk] - alpha[2*di + dk + ijk] + 
       alpha[-di + 2*dk + ijk] - alpha[di + 2*dk + ijk]) - 
    alpha[2*(-di + dk) + ijk] + alpha[-2*(di + dk) + ijk] + 
    alpha[2*(di + dk) + ijk])
;

dda22
=
0.083333333333333333333*oody2*(-30.*alpha[ijk] - alpha[-2*dj + ijk] + 
    16.*(alpha[-dj + ijk] + alpha[dj + ijk]) - alpha[2*dj + ijk])
;

dda23
=
0.027777777777777777778*oo2dy*oo2dz*
  (-alpha[2*(dj - dk) + ijk] + 
    64.*(alpha[-dj - dk + ijk] - alpha[dj - dk + ijk] - 
       alpha[-dj + dk + ijk] + alpha[dj + dk + ijk]) + 
    8.*(-alpha[-dj - 2*dk + ijk] + alpha[dj - 2*dk + ijk] - 
       alpha[-2*dj - dk + ijk] + alpha[2*dj - dk + ijk] + 
       alpha[-2*dj + dk + ijk] - alpha[2*dj + dk + ijk] + 
       alpha[-dj + 2*dk + ijk] - alpha[dj + 2*dk + ijk]) - 
    alpha[2*(-dj + dk) + ijk] + alpha[-2*(dj + dk) + ijk] + 
    alpha[2*(dj + dk) + ijk])
;

dda33
=
0.083333333333333333333*oodz2*(-30.*alpha[ijk] - alpha[-2*dk + ijk] + 
    16.*(alpha[-dk + ijk] + alpha[dk + ijk]) - alpha[2*dk + ijk])
;

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

ddb111
=
0.083333333333333333333*oodx2*(-30.*beta1[ijk] - beta1[-2*di + ijk] + 
    16.*(beta1[-di + ijk] + beta1[di + ijk]) - beta1[2*di + ijk])
;

ddb112
=
0.083333333333333333333*oodx2*(-30.*beta2[ijk] - beta2[-2*di + ijk] + 
    16.*(beta2[-di + ijk] + beta2[di + ijk]) - beta2[2*di + ijk])
;

ddb113
=
0.083333333333333333333*oodx2*(-30.*beta3[ijk] - beta3[-2*di + ijk] + 
    16.*(beta3[-di + ijk] + beta3[di + ijk]) - beta3[2*di + ijk])
;

ddb121
=
0.027777777777777777778*oo2dx*oo2dy*
  (-beta1[2*(di - dj) + ijk] + 
    64.*(beta1[-di - dj + ijk] - beta1[di - dj + ijk] - 
       beta1[-di + dj + ijk] + beta1[di + dj + ijk]) + 
    8.*(-beta1[-di - 2*dj + ijk] + beta1[di - 2*dj + ijk] - 
       beta1[-2*di - dj + ijk] + beta1[2*di - dj + ijk] + 
       beta1[-2*di + dj + ijk] - beta1[2*di + dj + ijk] + 
       beta1[-di + 2*dj + ijk] - beta1[di + 2*dj + ijk]) - 
    beta1[2*(-di + dj) + ijk] + beta1[-2*(di + dj) + ijk] + 
    beta1[2*(di + dj) + ijk])
;

ddb122
=
0.027777777777777777778*oo2dx*oo2dy*
  (-beta2[2*(di - dj) + ijk] + 
    64.*(beta2[-di - dj + ijk] - beta2[di - dj + ijk] - 
       beta2[-di + dj + ijk] + beta2[di + dj + ijk]) + 
    8.*(-beta2[-di - 2*dj + ijk] + beta2[di - 2*dj + ijk] - 
       beta2[-2*di - dj + ijk] + beta2[2*di - dj + ijk] + 
       beta2[-2*di + dj + ijk] - beta2[2*di + dj + ijk] + 
       beta2[-di + 2*dj + ijk] - beta2[di + 2*dj + ijk]) - 
    beta2[2*(-di + dj) + ijk] + beta2[-2*(di + dj) + ijk] + 
    beta2[2*(di + dj) + ijk])
;

ddb123
=
0.027777777777777777778*oo2dx*oo2dy*
  (-beta3[2*(di - dj) + ijk] + 
    64.*(beta3[-di - dj + ijk] - beta3[di - dj + ijk] - 
       beta3[-di + dj + ijk] + beta3[di + dj + ijk]) + 
    8.*(-beta3[-di - 2*dj + ijk] + beta3[di - 2*dj + ijk] - 
       beta3[-2*di - dj + ijk] + beta3[2*di - dj + ijk] + 
       beta3[-2*di + dj + ijk] - beta3[2*di + dj + ijk] + 
       beta3[-di + 2*dj + ijk] - beta3[di + 2*dj + ijk]) - 
    beta3[2*(-di + dj) + ijk] + beta3[-2*(di + dj) + ijk] + 
    beta3[2*(di + dj) + ijk])
;

ddb131
=
0.027777777777777777778*oo2dx*oo2dz*
  (-beta1[2*(di - dk) + ijk] + 
    64.*(beta1[-di - dk + ijk] - beta1[di - dk + ijk] - 
       beta1[-di + dk + ijk] + beta1[di + dk + ijk]) + 
    8.*(-beta1[-di - 2*dk + ijk] + beta1[di - 2*dk + ijk] - 
       beta1[-2*di - dk + ijk] + beta1[2*di - dk + ijk] + 
       beta1[-2*di + dk + ijk] - beta1[2*di + dk + ijk] + 
       beta1[-di + 2*dk + ijk] - beta1[di + 2*dk + ijk]) - 
    beta1[2*(-di + dk) + ijk] + beta1[-2*(di + dk) + ijk] + 
    beta1[2*(di + dk) + ijk])
;

ddb132
=
0.027777777777777777778*oo2dx*oo2dz*
  (-beta2[2*(di - dk) + ijk] + 
    64.*(beta2[-di - dk + ijk] - beta2[di - dk + ijk] - 
       beta2[-di + dk + ijk] + beta2[di + dk + ijk]) + 
    8.*(-beta2[-di - 2*dk + ijk] + beta2[di - 2*dk + ijk] - 
       beta2[-2*di - dk + ijk] + beta2[2*di - dk + ijk] + 
       beta2[-2*di + dk + ijk] - beta2[2*di + dk + ijk] + 
       beta2[-di + 2*dk + ijk] - beta2[di + 2*dk + ijk]) - 
    beta2[2*(-di + dk) + ijk] + beta2[-2*(di + dk) + ijk] + 
    beta2[2*(di + dk) + ijk])
;

ddb133
=
0.027777777777777777778*oo2dx*oo2dz*
  (-beta3[2*(di - dk) + ijk] + 
    64.*(beta3[-di - dk + ijk] - beta3[di - dk + ijk] - 
       beta3[-di + dk + ijk] + beta3[di + dk + ijk]) + 
    8.*(-beta3[-di - 2*dk + ijk] + beta3[di - 2*dk + ijk] - 
       beta3[-2*di - dk + ijk] + beta3[2*di - dk + ijk] + 
       beta3[-2*di + dk + ijk] - beta3[2*di + dk + ijk] + 
       beta3[-di + 2*dk + ijk] - beta3[di + 2*dk + ijk]) - 
    beta3[2*(-di + dk) + ijk] + beta3[-2*(di + dk) + ijk] + 
    beta3[2*(di + dk) + ijk])
;

ddb221
=
0.083333333333333333333*oody2*(-30.*beta1[ijk] - beta1[-2*dj + ijk] + 
    16.*(beta1[-dj + ijk] + beta1[dj + ijk]) - beta1[2*dj + ijk])
;

ddb222
=
0.083333333333333333333*oody2*(-30.*beta2[ijk] - beta2[-2*dj + ijk] + 
    16.*(beta2[-dj + ijk] + beta2[dj + ijk]) - beta2[2*dj + ijk])
;

ddb223
=
0.083333333333333333333*oody2*(-30.*beta3[ijk] - beta3[-2*dj + ijk] + 
    16.*(beta3[-dj + ijk] + beta3[dj + ijk]) - beta3[2*dj + ijk])
;

ddb231
=
0.027777777777777777778*oo2dy*oo2dz*
  (-beta1[2*(dj - dk) + ijk] + 
    64.*(beta1[-dj - dk + ijk] - beta1[dj - dk + ijk] - 
       beta1[-dj + dk + ijk] + beta1[dj + dk + ijk]) + 
    8.*(-beta1[-dj - 2*dk + ijk] + beta1[dj - 2*dk + ijk] - 
       beta1[-2*dj - dk + ijk] + beta1[2*dj - dk + ijk] + 
       beta1[-2*dj + dk + ijk] - beta1[2*dj + dk + ijk] + 
       beta1[-dj + 2*dk + ijk] - beta1[dj + 2*dk + ijk]) - 
    beta1[2*(-dj + dk) + ijk] + beta1[-2*(dj + dk) + ijk] + 
    beta1[2*(dj + dk) + ijk])
;

ddb232
=
0.027777777777777777778*oo2dy*oo2dz*
  (-beta2[2*(dj - dk) + ijk] + 
    64.*(beta2[-dj - dk + ijk] - beta2[dj - dk + ijk] - 
       beta2[-dj + dk + ijk] + beta2[dj + dk + ijk]) + 
    8.*(-beta2[-dj - 2*dk + ijk] + beta2[dj - 2*dk + ijk] - 
       beta2[-2*dj - dk + ijk] + beta2[2*dj - dk + ijk] + 
       beta2[-2*dj + dk + ijk] - beta2[2*dj + dk + ijk] + 
       beta2[-dj + 2*dk + ijk] - beta2[dj + 2*dk + ijk]) - 
    beta2[2*(-dj + dk) + ijk] + beta2[-2*(dj + dk) + ijk] + 
    beta2[2*(dj + dk) + ijk])
;

ddb233
=
0.027777777777777777778*oo2dy*oo2dz*
  (-beta3[2*(dj - dk) + ijk] + 
    64.*(beta3[-dj - dk + ijk] - beta3[dj - dk + ijk] - 
       beta3[-dj + dk + ijk] + beta3[dj + dk + ijk]) + 
    8.*(-beta3[-dj - 2*dk + ijk] + beta3[dj - 2*dk + ijk] - 
       beta3[-2*dj - dk + ijk] + beta3[2*dj - dk + ijk] + 
       beta3[-2*dj + dk + ijk] - beta3[2*dj + dk + ijk] + 
       beta3[-dj + 2*dk + ijk] - beta3[dj + 2*dk + ijk]) - 
    beta3[2*(-dj + dk) + ijk] + beta3[-2*(dj + dk) + ijk] + 
    beta3[2*(dj + dk) + ijk])
;

ddb331
=
0.083333333333333333333*oodz2*(-30.*beta1[ijk] - beta1[-2*dk + ijk] + 
    16.*(beta1[-dk + ijk] + beta1[dk + ijk]) - beta1[2*dk + ijk])
;

ddb332
=
0.083333333333333333333*oodz2*(-30.*beta2[ijk] - beta2[-2*dk + ijk] + 
    16.*(beta2[-dk + ijk] + beta2[dk + ijk]) - beta2[2*dk + ijk])
;

ddb333
=
0.083333333333333333333*oodz2*(-30.*beta3[ijk] - beta3[-2*dk + ijk] + 
    16.*(beta3[-dk + ijk] + beta3[dk + ijk]) - beta3[2*dk + ijk])
;

delg111
=
0.16666666666666666667*oo2dx*(g11[-2*di + ijk] + 
    8.*(-g11[-di + ijk] + g11[di + ijk]) - g11[2*di + ijk])
;

delg112
=
0.16666666666666666667*oo2dx*(g12[-2*di + ijk] + 
    8.*(-g12[-di + ijk] + g12[di + ijk]) - g12[2*di + ijk])
;

delg113
=
0.16666666666666666667*oo2dx*(g13[-2*di + ijk] + 
    8.*(-g13[-di + ijk] + g13[di + ijk]) - g13[2*di + ijk])
;

delg122
=
0.16666666666666666667*oo2dx*(g22[-2*di + ijk] + 
    8.*(-g22[-di + ijk] + g22[di + ijk]) - g22[2*di + ijk])
;

delg123
=
0.16666666666666666667*oo2dx*(g23[-2*di + ijk] + 
    8.*(-g23[-di + ijk] + g23[di + ijk]) - g23[2*di + ijk])
;

delg133
=
0.16666666666666666667*oo2dx*(g33[-2*di + ijk] + 
    8.*(-g33[-di + ijk] + g33[di + ijk]) - g33[2*di + ijk])
;

delg211
=
0.16666666666666666667*oo2dy*(g11[-2*dj + ijk] + 
    8.*(-g11[-dj + ijk] + g11[dj + ijk]) - g11[2*dj + ijk])
;

delg212
=
0.16666666666666666667*oo2dy*(g12[-2*dj + ijk] + 
    8.*(-g12[-dj + ijk] + g12[dj + ijk]) - g12[2*dj + ijk])
;

delg213
=
0.16666666666666666667*oo2dy*(g13[-2*dj + ijk] + 
    8.*(-g13[-dj + ijk] + g13[dj + ijk]) - g13[2*dj + ijk])
;

delg222
=
0.16666666666666666667*oo2dy*(g22[-2*dj + ijk] + 
    8.*(-g22[-dj + ijk] + g22[dj + ijk]) - g22[2*dj + ijk])
;

delg223
=
0.16666666666666666667*oo2dy*(g23[-2*dj + ijk] + 
    8.*(-g23[-dj + ijk] + g23[dj + ijk]) - g23[2*dj + ijk])
;

delg233
=
0.16666666666666666667*oo2dy*(g33[-2*dj + ijk] + 
    8.*(-g33[-dj + ijk] + g33[dj + ijk]) - g33[2*dj + ijk])
;

delg311
=
0.16666666666666666667*oo2dz*(g11[-2*dk + ijk] + 
    8.*(-g11[-dk + ijk] + g11[dk + ijk]) - g11[2*dk + ijk])
;

delg312
=
0.16666666666666666667*oo2dz*(g12[-2*dk + ijk] + 
    8.*(-g12[-dk + ijk] + g12[dk + ijk]) - g12[2*dk + ijk])
;

delg313
=
0.16666666666666666667*oo2dz*(g13[-2*dk + ijk] + 
    8.*(-g13[-dk + ijk] + g13[dk + ijk]) - g13[2*dk + ijk])
;

delg322
=
0.16666666666666666667*oo2dz*(g22[-2*dk + ijk] + 
    8.*(-g22[-dk + ijk] + g22[dk + ijk]) - g22[2*dk + ijk])
;

delg323
=
0.16666666666666666667*oo2dz*(g23[-2*dk + ijk] + 
    8.*(-g23[-dk + ijk] + g23[dk + ijk]) - g23[2*dk + ijk])
;

delg333
=
0.16666666666666666667*oo2dz*(g33[-2*dk + ijk] + 
    8.*(-g33[-dk + ijk] + g33[dk + ijk]) - g33[2*dk + ijk])
;

deldelg1111
=
0.083333333333333333333*oodx2*(-30.*g11[ijk] - g11[-2*di + ijk] + 
    16.*(g11[-di + ijk] + g11[di + ijk]) - g11[2*di + ijk])
;

deldelg1112
=
0.083333333333333333333*oodx2*(-30.*g12[ijk] - g12[-2*di + ijk] + 
    16.*(g12[-di + ijk] + g12[di + ijk]) - g12[2*di + ijk])
;

deldelg1113
=
0.083333333333333333333*oodx2*(-30.*g13[ijk] - g13[-2*di + ijk] + 
    16.*(g13[-di + ijk] + g13[di + ijk]) - g13[2*di + ijk])
;

deldelg1122
=
0.083333333333333333333*oodx2*(-30.*g22[ijk] - g22[-2*di + ijk] + 
    16.*(g22[-di + ijk] + g22[di + ijk]) - g22[2*di + ijk])
;

deldelg1123
=
0.083333333333333333333*oodx2*(-30.*g23[ijk] - g23[-2*di + ijk] + 
    16.*(g23[-di + ijk] + g23[di + ijk]) - g23[2*di + ijk])
;

deldelg1133
=
0.083333333333333333333*oodx2*(-30.*g33[ijk] - g33[-2*di + ijk] + 
    16.*(g33[-di + ijk] + g33[di + ijk]) - g33[2*di + ijk])
;

deldelg1211
=
0.027777777777777777778*oo2dx*oo2dy*
  (-g11[2*(di - dj) + ijk] + 64.*
     (g11[-di - dj + ijk] - g11[di - dj + ijk] - g11[-di + dj + ijk] + 
       g11[di + dj + ijk]) + 8.*
     (-g11[-di - 2*dj + ijk] + g11[di - 2*dj + ijk] - 
       g11[-2*di - dj + ijk] + g11[2*di - dj + ijk] + 
       g11[-2*di + dj + ijk] - g11[2*di + dj + ijk] + 
       g11[-di + 2*dj + ijk] - g11[di + 2*dj + ijk]) - 
    g11[2*(-di + dj) + ijk] + g11[-2*(di + dj) + ijk] + 
    g11[2*(di + dj) + ijk])
;

deldelg1212
=
0.027777777777777777778*oo2dx*oo2dy*
  (-g12[2*(di - dj) + ijk] + 64.*
     (g12[-di - dj + ijk] - g12[di - dj + ijk] - g12[-di + dj + ijk] + 
       g12[di + dj + ijk]) + 8.*
     (-g12[-di - 2*dj + ijk] + g12[di - 2*dj + ijk] - 
       g12[-2*di - dj + ijk] + g12[2*di - dj + ijk] + 
       g12[-2*di + dj + ijk] - g12[2*di + dj + ijk] + 
       g12[-di + 2*dj + ijk] - g12[di + 2*dj + ijk]) - 
    g12[2*(-di + dj) + ijk] + g12[-2*(di + dj) + ijk] + 
    g12[2*(di + dj) + ijk])
;

deldelg1213
=
0.027777777777777777778*oo2dx*oo2dy*
  (-g13[2*(di - dj) + ijk] + 64.*
     (g13[-di - dj + ijk] - g13[di - dj + ijk] - g13[-di + dj + ijk] + 
       g13[di + dj + ijk]) + 8.*
     (-g13[-di - 2*dj + ijk] + g13[di - 2*dj + ijk] - 
       g13[-2*di - dj + ijk] + g13[2*di - dj + ijk] + 
       g13[-2*di + dj + ijk] - g13[2*di + dj + ijk] + 
       g13[-di + 2*dj + ijk] - g13[di + 2*dj + ijk]) - 
    g13[2*(-di + dj) + ijk] + g13[-2*(di + dj) + ijk] + 
    g13[2*(di + dj) + ijk])
;

deldelg1222
=
0.027777777777777777778*oo2dx*oo2dy*
  (-g22[2*(di - dj) + ijk] + 64.*
     (g22[-di - dj + ijk] - g22[di - dj + ijk] - g22[-di + dj + ijk] + 
       g22[di + dj + ijk]) + 8.*
     (-g22[-di - 2*dj + ijk] + g22[di - 2*dj + ijk] - 
       g22[-2*di - dj + ijk] + g22[2*di - dj + ijk] + 
       g22[-2*di + dj + ijk] - g22[2*di + dj + ijk] + 
       g22[-di + 2*dj + ijk] - g22[di + 2*dj + ijk]) - 
    g22[2*(-di + dj) + ijk] + g22[-2*(di + dj) + ijk] + 
    g22[2*(di + dj) + ijk])
;

deldelg1223
=
0.027777777777777777778*oo2dx*oo2dy*
  (-g23[2*(di - dj) + ijk] + 64.*
     (g23[-di - dj + ijk] - g23[di - dj + ijk] - g23[-di + dj + ijk] + 
       g23[di + dj + ijk]) + 8.*
     (-g23[-di - 2*dj + ijk] + g23[di - 2*dj + ijk] - 
       g23[-2*di - dj + ijk] + g23[2*di - dj + ijk] + 
       g23[-2*di + dj + ijk] - g23[2*di + dj + ijk] + 
       g23[-di + 2*dj + ijk] - g23[di + 2*dj + ijk]) - 
    g23[2*(-di + dj) + ijk] + g23[-2*(di + dj) + ijk] + 
    g23[2*(di + dj) + ijk])
;

deldelg1233
=
0.027777777777777777778*oo2dx*oo2dy*
  (-g33[2*(di - dj) + ijk] + 64.*
     (g33[-di - dj + ijk] - g33[di - dj + ijk] - g33[-di + dj + ijk] + 
       g33[di + dj + ijk]) + 8.*
     (-g33[-di - 2*dj + ijk] + g33[di - 2*dj + ijk] - 
       g33[-2*di - dj + ijk] + g33[2*di - dj + ijk] + 
       g33[-2*di + dj + ijk] - g33[2*di + dj + ijk] + 
       g33[-di + 2*dj + ijk] - g33[di + 2*dj + ijk]) - 
    g33[2*(-di + dj) + ijk] + g33[-2*(di + dj) + ijk] + 
    g33[2*(di + dj) + ijk])
;

deldelg1311
=
0.027777777777777777778*oo2dx*oo2dz*
  (-g11[2*(di - dk) + ijk] + 64.*
     (g11[-di - dk + ijk] - g11[di - dk + ijk] - g11[-di + dk + ijk] + 
       g11[di + dk + ijk]) + 8.*
     (-g11[-di - 2*dk + ijk] + g11[di - 2*dk + ijk] - 
       g11[-2*di - dk + ijk] + g11[2*di - dk + ijk] + 
       g11[-2*di + dk + ijk] - g11[2*di + dk + ijk] + 
       g11[-di + 2*dk + ijk] - g11[di + 2*dk + ijk]) - 
    g11[2*(-di + dk) + ijk] + g11[-2*(di + dk) + ijk] + 
    g11[2*(di + dk) + ijk])
;

deldelg1312
=
0.027777777777777777778*oo2dx*oo2dz*
  (-g12[2*(di - dk) + ijk] + 64.*
     (g12[-di - dk + ijk] - g12[di - dk + ijk] - g12[-di + dk + ijk] + 
       g12[di + dk + ijk]) + 8.*
     (-g12[-di - 2*dk + ijk] + g12[di - 2*dk + ijk] - 
       g12[-2*di - dk + ijk] + g12[2*di - dk + ijk] + 
       g12[-2*di + dk + ijk] - g12[2*di + dk + ijk] + 
       g12[-di + 2*dk + ijk] - g12[di + 2*dk + ijk]) - 
    g12[2*(-di + dk) + ijk] + g12[-2*(di + dk) + ijk] + 
    g12[2*(di + dk) + ijk])
;

deldelg1313
=
0.027777777777777777778*oo2dx*oo2dz*
  (-g13[2*(di - dk) + ijk] + 64.*
     (g13[-di - dk + ijk] - g13[di - dk + ijk] - g13[-di + dk + ijk] + 
       g13[di + dk + ijk]) + 8.*
     (-g13[-di - 2*dk + ijk] + g13[di - 2*dk + ijk] - 
       g13[-2*di - dk + ijk] + g13[2*di - dk + ijk] + 
       g13[-2*di + dk + ijk] - g13[2*di + dk + ijk] + 
       g13[-di + 2*dk + ijk] - g13[di + 2*dk + ijk]) - 
    g13[2*(-di + dk) + ijk] + g13[-2*(di + dk) + ijk] + 
    g13[2*(di + dk) + ijk])
;

deldelg1322
=
0.027777777777777777778*oo2dx*oo2dz*
  (-g22[2*(di - dk) + ijk] + 64.*
     (g22[-di - dk + ijk] - g22[di - dk + ijk] - g22[-di + dk + ijk] + 
       g22[di + dk + ijk]) + 8.*
     (-g22[-di - 2*dk + ijk] + g22[di - 2*dk + ijk] - 
       g22[-2*di - dk + ijk] + g22[2*di - dk + ijk] + 
       g22[-2*di + dk + ijk] - g22[2*di + dk + ijk] + 
       g22[-di + 2*dk + ijk] - g22[di + 2*dk + ijk]) - 
    g22[2*(-di + dk) + ijk] + g22[-2*(di + dk) + ijk] + 
    g22[2*(di + dk) + ijk])
;

deldelg1323
=
0.027777777777777777778*oo2dx*oo2dz*
  (-g23[2*(di - dk) + ijk] + 64.*
     (g23[-di - dk + ijk] - g23[di - dk + ijk] - g23[-di + dk + ijk] + 
       g23[di + dk + ijk]) + 8.*
     (-g23[-di - 2*dk + ijk] + g23[di - 2*dk + ijk] - 
       g23[-2*di - dk + ijk] + g23[2*di - dk + ijk] + 
       g23[-2*di + dk + ijk] - g23[2*di + dk + ijk] + 
       g23[-di + 2*dk + ijk] - g23[di + 2*dk + ijk]) - 
    g23[2*(-di + dk) + ijk] + g23[-2*(di + dk) + ijk] + 
    g23[2*(di + dk) + ijk])
;

deldelg1333
=
0.027777777777777777778*oo2dx*oo2dz*
  (-g33[2*(di - dk) + ijk] + 64.*
     (g33[-di - dk + ijk] - g33[di - dk + ijk] - g33[-di + dk + ijk] + 
       g33[di + dk + ijk]) + 8.*
     (-g33[-di - 2*dk + ijk] + g33[di - 2*dk + ijk] - 
       g33[-2*di - dk + ijk] + g33[2*di - dk + ijk] + 
       g33[-2*di + dk + ijk] - g33[2*di + dk + ijk] + 
       g33[-di + 2*dk + ijk] - g33[di + 2*dk + ijk]) - 
    g33[2*(-di + dk) + ijk] + g33[-2*(di + dk) + ijk] + 
    g33[2*(di + dk) + ijk])
;

deldelg2211
=
0.083333333333333333333*oody2*(-30.*g11[ijk] - g11[-2*dj + ijk] + 
    16.*(g11[-dj + ijk] + g11[dj + ijk]) - g11[2*dj + ijk])
;

deldelg2212
=
0.083333333333333333333*oody2*(-30.*g12[ijk] - g12[-2*dj + ijk] + 
    16.*(g12[-dj + ijk] + g12[dj + ijk]) - g12[2*dj + ijk])
;

deldelg2213
=
0.083333333333333333333*oody2*(-30.*g13[ijk] - g13[-2*dj + ijk] + 
    16.*(g13[-dj + ijk] + g13[dj + ijk]) - g13[2*dj + ijk])
;

deldelg2222
=
0.083333333333333333333*oody2*(-30.*g22[ijk] - g22[-2*dj + ijk] + 
    16.*(g22[-dj + ijk] + g22[dj + ijk]) - g22[2*dj + ijk])
;

deldelg2223
=
0.083333333333333333333*oody2*(-30.*g23[ijk] - g23[-2*dj + ijk] + 
    16.*(g23[-dj + ijk] + g23[dj + ijk]) - g23[2*dj + ijk])
;

deldelg2233
=
0.083333333333333333333*oody2*(-30.*g33[ijk] - g33[-2*dj + ijk] + 
    16.*(g33[-dj + ijk] + g33[dj + ijk]) - g33[2*dj + ijk])
;

deldelg2311
=
0.027777777777777777778*oo2dy*oo2dz*
  (-g11[2*(dj - dk) + ijk] + 64.*
     (g11[-dj - dk + ijk] - g11[dj - dk + ijk] - g11[-dj + dk + ijk] + 
       g11[dj + dk + ijk]) + 8.*
     (-g11[-dj - 2*dk + ijk] + g11[dj - 2*dk + ijk] - 
       g11[-2*dj - dk + ijk] + g11[2*dj - dk + ijk] + 
       g11[-2*dj + dk + ijk] - g11[2*dj + dk + ijk] + 
       g11[-dj + 2*dk + ijk] - g11[dj + 2*dk + ijk]) - 
    g11[2*(-dj + dk) + ijk] + g11[-2*(dj + dk) + ijk] + 
    g11[2*(dj + dk) + ijk])
;

deldelg2312
=
0.027777777777777777778*oo2dy*oo2dz*
  (-g12[2*(dj - dk) + ijk] + 64.*
     (g12[-dj - dk + ijk] - g12[dj - dk + ijk] - g12[-dj + dk + ijk] + 
       g12[dj + dk + ijk]) + 8.*
     (-g12[-dj - 2*dk + ijk] + g12[dj - 2*dk + ijk] - 
       g12[-2*dj - dk + ijk] + g12[2*dj - dk + ijk] + 
       g12[-2*dj + dk + ijk] - g12[2*dj + dk + ijk] + 
       g12[-dj + 2*dk + ijk] - g12[dj + 2*dk + ijk]) - 
    g12[2*(-dj + dk) + ijk] + g12[-2*(dj + dk) + ijk] + 
    g12[2*(dj + dk) + ijk])
;

deldelg2313
=
0.027777777777777777778*oo2dy*oo2dz*
  (-g13[2*(dj - dk) + ijk] + 64.*
     (g13[-dj - dk + ijk] - g13[dj - dk + ijk] - g13[-dj + dk + ijk] + 
       g13[dj + dk + ijk]) + 8.*
     (-g13[-dj - 2*dk + ijk] + g13[dj - 2*dk + ijk] - 
       g13[-2*dj - dk + ijk] + g13[2*dj - dk + ijk] + 
       g13[-2*dj + dk + ijk] - g13[2*dj + dk + ijk] + 
       g13[-dj + 2*dk + ijk] - g13[dj + 2*dk + ijk]) - 
    g13[2*(-dj + dk) + ijk] + g13[-2*(dj + dk) + ijk] + 
    g13[2*(dj + dk) + ijk])
;

deldelg2322
=
0.027777777777777777778*oo2dy*oo2dz*
  (-g22[2*(dj - dk) + ijk] + 64.*
     (g22[-dj - dk + ijk] - g22[dj - dk + ijk] - g22[-dj + dk + ijk] + 
       g22[dj + dk + ijk]) + 8.*
     (-g22[-dj - 2*dk + ijk] + g22[dj - 2*dk + ijk] - 
       g22[-2*dj - dk + ijk] + g22[2*dj - dk + ijk] + 
       g22[-2*dj + dk + ijk] - g22[2*dj + dk + ijk] + 
       g22[-dj + 2*dk + ijk] - g22[dj + 2*dk + ijk]) - 
    g22[2*(-dj + dk) + ijk] + g22[-2*(dj + dk) + ijk] + 
    g22[2*(dj + dk) + ijk])
;

deldelg2323
=
0.027777777777777777778*oo2dy*oo2dz*
  (-g23[2*(dj - dk) + ijk] + 64.*
     (g23[-dj - dk + ijk] - g23[dj - dk + ijk] - g23[-dj + dk + ijk] + 
       g23[dj + dk + ijk]) + 8.*
     (-g23[-dj - 2*dk + ijk] + g23[dj - 2*dk + ijk] - 
       g23[-2*dj - dk + ijk] + g23[2*dj - dk + ijk] + 
       g23[-2*dj + dk + ijk] - g23[2*dj + dk + ijk] + 
       g23[-dj + 2*dk + ijk] - g23[dj + 2*dk + ijk]) - 
    g23[2*(-dj + dk) + ijk] + g23[-2*(dj + dk) + ijk] + 
    g23[2*(dj + dk) + ijk])
;

deldelg2333
=
0.027777777777777777778*oo2dy*oo2dz*
  (-g33[2*(dj - dk) + ijk] + 64.*
     (g33[-dj - dk + ijk] - g33[dj - dk + ijk] - g33[-dj + dk + ijk] + 
       g33[dj + dk + ijk]) + 8.*
     (-g33[-dj - 2*dk + ijk] + g33[dj - 2*dk + ijk] - 
       g33[-2*dj - dk + ijk] + g33[2*dj - dk + ijk] + 
       g33[-2*dj + dk + ijk] - g33[2*dj + dk + ijk] + 
       g33[-dj + 2*dk + ijk] - g33[dj + 2*dk + ijk]) - 
    g33[2*(-dj + dk) + ijk] + g33[-2*(dj + dk) + ijk] + 
    g33[2*(dj + dk) + ijk])
;

deldelg3311
=
0.083333333333333333333*oodz2*(-30.*g11[ijk] - g11[-2*dk + ijk] + 
    16.*(g11[-dk + ijk] + g11[dk + ijk]) - g11[2*dk + ijk])
;

deldelg3312
=
0.083333333333333333333*oodz2*(-30.*g12[ijk] - g12[-2*dk + ijk] + 
    16.*(g12[-dk + ijk] + g12[dk + ijk]) - g12[2*dk + ijk])
;

deldelg3313
=
0.083333333333333333333*oodz2*(-30.*g13[ijk] - g13[-2*dk + ijk] + 
    16.*(g13[-dk + ijk] + g13[dk + ijk]) - g13[2*dk + ijk])
;

deldelg3322
=
0.083333333333333333333*oodz2*(-30.*g22[ijk] - g22[-2*dk + ijk] + 
    16.*(g22[-dk + ijk] + g22[dk + ijk]) - g22[2*dk + ijk])
;

deldelg3323
=
0.083333333333333333333*oodz2*(-30.*g23[ijk] - g23[-2*dk + ijk] + 
    16.*(g23[-dk + ijk] + g23[dk + ijk]) - g23[2*dk + ijk])
;

deldelg3333
=
0.083333333333333333333*oodz2*(-30.*g33[ijk] - g33[-2*dk + ijk] + 
    16.*(g33[-dk + ijk] + g33[dk + ijk]) - g33[2*dk + ijk])
;

delG11
=
0.16666666666666666667*oo2dx*(G1[-2*di + ijk] + 
    8.*(-G1[-di + ijk] + G1[di + ijk]) - G1[2*di + ijk])
;

delG12
=
0.16666666666666666667*oo2dx*(G2[-2*di + ijk] + 
    8.*(-G2[-di + ijk] + G2[di + ijk]) - G2[2*di + ijk])
;

delG13
=
0.16666666666666666667*oo2dx*(G3[-2*di + ijk] + 
    8.*(-G3[-di + ijk] + G3[di + ijk]) - G3[2*di + ijk])
;

delG21
=
0.16666666666666666667*oo2dy*(G1[-2*dj + ijk] + 
    8.*(-G1[-dj + ijk] + G1[dj + ijk]) - G1[2*dj + ijk])
;

delG22
=
0.16666666666666666667*oo2dy*(G2[-2*dj + ijk] + 
    8.*(-G2[-dj + ijk] + G2[dj + ijk]) - G2[2*dj + ijk])
;

delG23
=
0.16666666666666666667*oo2dy*(G3[-2*dj + ijk] + 
    8.*(-G3[-dj + ijk] + G3[dj + ijk]) - G3[2*dj + ijk])
;

delG31
=
0.16666666666666666667*oo2dz*(G1[-2*dk + ijk] + 
    8.*(-G1[-dk + ijk] + G1[dk + ijk]) - G1[2*dk + ijk])
;

delG32
=
0.16666666666666666667*oo2dz*(G2[-2*dk + ijk] + 
    8.*(-G2[-dk + ijk] + G2[dk + ijk]) - G2[2*dk + ijk])
;

delG33
=
0.16666666666666666667*oo2dz*(G3[-2*dk + ijk] + 
    8.*(-G3[-dk + ijk] + G3[dk + ijk]) - G3[2*dk + ijk])
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

dchi1
=
0.16666666666666666667*oo2dx*(chi[-2*di + ijk] + 
    8.*(-chi[-di + ijk] + chi[di + ijk]) - chi[2*di + ijk])
;

dchi2
=
0.16666666666666666667*oo2dy*(chi[-2*dj + ijk] + 
    8.*(-chi[-dj + ijk] + chi[dj + ijk]) - chi[2*dj + ijk])
;

dchi3
=
0.16666666666666666667*oo2dz*(chi[-2*dk + ijk] + 
    8.*(-chi[-dk + ijk] + chi[dk + ijk]) - chi[2*dk + ijk])
;

ddchi11
=
0.083333333333333333333*oodx2*(-30.*chi[ijk] - chi[-2*di + ijk] + 
    16.*(chi[-di + ijk] + chi[di + ijk]) - chi[2*di + ijk])
;

ddchi12
=
0.027777777777777777778*oo2dx*oo2dy*
  (-chi[2*(di - dj) + ijk] + 64.*
     (chi[-di - dj + ijk] - chi[di - dj + ijk] - chi[-di + dj + ijk] + 
       chi[di + dj + ijk]) + 8.*
     (-chi[-di - 2*dj + ijk] + chi[di - 2*dj + ijk] - 
       chi[-2*di - dj + ijk] + chi[2*di - dj + ijk] + 
       chi[-2*di + dj + ijk] - chi[2*di + dj + ijk] + 
       chi[-di + 2*dj + ijk] - chi[di + 2*dj + ijk]) - 
    chi[2*(-di + dj) + ijk] + chi[-2*(di + dj) + ijk] + 
    chi[2*(di + dj) + ijk])
;

ddchi13
=
0.027777777777777777778*oo2dx*oo2dz*
  (-chi[2*(di - dk) + ijk] + 64.*
     (chi[-di - dk + ijk] - chi[di - dk + ijk] - chi[-di + dk + ijk] + 
       chi[di + dk + ijk]) + 8.*
     (-chi[-di - 2*dk + ijk] + chi[di - 2*dk + ijk] - 
       chi[-2*di - dk + ijk] + chi[2*di - dk + ijk] + 
       chi[-2*di + dk + ijk] - chi[2*di + dk + ijk] + 
       chi[-di + 2*dk + ijk] - chi[di + 2*dk + ijk]) - 
    chi[2*(-di + dk) + ijk] + chi[-2*(di + dk) + ijk] + 
    chi[2*(di + dk) + ijk])
;

ddchi22
=
0.083333333333333333333*oody2*(-30.*chi[ijk] - chi[-2*dj + ijk] + 
    16.*(chi[-dj + ijk] + chi[dj + ijk]) - chi[2*dj + ijk])
;

ddchi23
=
0.027777777777777777778*oo2dy*oo2dz*
  (-chi[2*(dj - dk) + ijk] + 64.*
     (chi[-dj - dk + ijk] - chi[dj - dk + ijk] - chi[-dj + dk + ijk] + 
       chi[dj + dk + ijk]) + 8.*
     (-chi[-dj - 2*dk + ijk] + chi[dj - 2*dk + ijk] - 
       chi[-2*dj - dk + ijk] + chi[2*dj - dk + ijk] + 
       chi[-2*dj + dk + ijk] - chi[2*dj + dk + ijk] + 
       chi[-dj + 2*dk + ijk] - chi[dj + 2*dk + ijk]) - 
    chi[2*(-dj + dk) + ijk] + chi[-2*(dj + dk) + ijk] + 
    chi[2*(dj + dk) + ijk])
;

ddchi33
=
0.083333333333333333333*oodz2*(-30.*chi[ijk] - chi[-2*dk + ijk] + 
    16.*(chi[-dk + ijk] + chi[dk + ijk]) - chi[2*dk + ijk])
;


#endif 


} else if (order_centered == 6 || boundaryNaway(3)) { 


#ifdef REDUCEORDERTO4                           
  errorexit("CompilerFlag reduced order to 4");
 #else
da1
=
0.033333333333333333333*oo2dx*(-alpha[-3*di + ijk] + 
    45.*(-alpha[-di + ijk] + alpha[di + ijk]) + 
    9.*(alpha[-2*di + ijk] - alpha[2*di + ijk]) + alpha[3*di + ijk])
;

da2
=
0.033333333333333333333*oo2dy*(-alpha[-3*dj + ijk] + 
    45.*(-alpha[-dj + ijk] + alpha[dj + ijk]) + 
    9.*(alpha[-2*dj + ijk] - alpha[2*dj + ijk]) + alpha[3*dj + ijk])
;

da3
=
0.033333333333333333333*oo2dz*(-alpha[-3*dk + ijk] + 
    45.*(-alpha[-dk + ijk] + alpha[dk + ijk]) + 
    9.*(alpha[-2*dk + ijk] - alpha[2*dk + ijk]) + alpha[3*dk + ijk])
;

dda11
=
0.0055555555555555555556*oodx2*
  (-490.*alpha[ijk] + 270.*(alpha[-di + ijk] + alpha[di + ijk]) - 
    27.*(alpha[-2*di + ijk] + alpha[2*di + ijk]) + 
    2.*(alpha[-3*di + ijk] + alpha[3*di + ijk]))
;

dda12
=
0.0011111111111111111111*oo2dx*oo2dy*
  (-alpha[3*(di - dj) + ijk] - 
    405.*(alpha[-di - 2*dj + ijk] + alpha[-2*di - dj + ijk]) + 
    2025.*(alpha[-di - dj + ijk] - alpha[di - dj + ijk] - 
       alpha[-di + dj + ijk] + alpha[di + dj + ijk]) + 
    405.*(alpha[di - 2*dj + ijk] + alpha[2*di - dj + ijk] + 
       alpha[-2*di + dj + ijk] - alpha[2*di + dj + ijk] + 
       alpha[-di + 2*dj + ijk] - alpha[di + 2*dj + ijk]) + 
    45.*(alpha[-di - 3*dj + ijk] - alpha[di - 3*dj + ijk] + 
       alpha[-3*di - dj + ijk] - alpha[3*di - dj + ijk] - 
       alpha[-3*di + dj + ijk] + alpha[3*di + dj + ijk] - 
       alpha[-di + 3*dj + ijk] + alpha[di + 3*dj + ijk]) + 
    9.*(-alpha[-2*di - 3*dj + ijk] + alpha[2*di - 3*dj + ijk] - 
       alpha[-3*di - 2*dj + ijk] + alpha[3*di - 2*dj + ijk] + 
       alpha[-3*di + 2*dj + ijk] - alpha[3*di + 2*dj + ijk] + 
       alpha[-2*di + 3*dj + ijk] - alpha[2*di + 3*dj + ijk]) - 
    alpha[3*(-di + dj) + ijk] + alpha[-3*(di + dj) + ijk] + 
    81.*(-alpha[2*(di - dj) + ijk] - alpha[2*(-di + dj) + ijk] + 
       alpha[-2*(di + dj) + ijk] + alpha[2*(di + dj) + ijk]) + 
    alpha[3*(di + dj) + ijk])
;

dda13
=
0.0011111111111111111111*oo2dx*oo2dz*
  (-alpha[3*(di - dk) + ijk] - 
    405.*(alpha[-di - 2*dk + ijk] + alpha[-2*di - dk + ijk]) + 
    2025.*(alpha[-di - dk + ijk] - alpha[di - dk + ijk] - 
       alpha[-di + dk + ijk] + alpha[di + dk + ijk]) + 
    405.*(alpha[di - 2*dk + ijk] + alpha[2*di - dk + ijk] + 
       alpha[-2*di + dk + ijk] - alpha[2*di + dk + ijk] + 
       alpha[-di + 2*dk + ijk] - alpha[di + 2*dk + ijk]) + 
    45.*(alpha[-di - 3*dk + ijk] - alpha[di - 3*dk + ijk] + 
       alpha[-3*di - dk + ijk] - alpha[3*di - dk + ijk] - 
       alpha[-3*di + dk + ijk] + alpha[3*di + dk + ijk] - 
       alpha[-di + 3*dk + ijk] + alpha[di + 3*dk + ijk]) + 
    9.*(-alpha[-2*di - 3*dk + ijk] + alpha[2*di - 3*dk + ijk] - 
       alpha[-3*di - 2*dk + ijk] + alpha[3*di - 2*dk + ijk] + 
       alpha[-3*di + 2*dk + ijk] - alpha[3*di + 2*dk + ijk] + 
       alpha[-2*di + 3*dk + ijk] - alpha[2*di + 3*dk + ijk]) - 
    alpha[3*(-di + dk) + ijk] + alpha[-3*(di + dk) + ijk] + 
    81.*(-alpha[2*(di - dk) + ijk] - alpha[2*(-di + dk) + ijk] + 
       alpha[-2*(di + dk) + ijk] + alpha[2*(di + dk) + ijk]) + 
    alpha[3*(di + dk) + ijk])
;

dda22
=
0.0055555555555555555556*oody2*
  (-490.*alpha[ijk] + 270.*(alpha[-dj + ijk] + alpha[dj + ijk]) - 
    27.*(alpha[-2*dj + ijk] + alpha[2*dj + ijk]) + 
    2.*(alpha[-3*dj + ijk] + alpha[3*dj + ijk]))
;

dda23
=
0.0011111111111111111111*oo2dy*oo2dz*
  (-alpha[3*(dj - dk) + ijk] - 
    405.*(alpha[-dj - 2*dk + ijk] + alpha[-2*dj - dk + ijk]) + 
    2025.*(alpha[-dj - dk + ijk] - alpha[dj - dk + ijk] - 
       alpha[-dj + dk + ijk] + alpha[dj + dk + ijk]) + 
    405.*(alpha[dj - 2*dk + ijk] + alpha[2*dj - dk + ijk] + 
       alpha[-2*dj + dk + ijk] - alpha[2*dj + dk + ijk] + 
       alpha[-dj + 2*dk + ijk] - alpha[dj + 2*dk + ijk]) + 
    45.*(alpha[-dj - 3*dk + ijk] - alpha[dj - 3*dk + ijk] + 
       alpha[-3*dj - dk + ijk] - alpha[3*dj - dk + ijk] - 
       alpha[-3*dj + dk + ijk] + alpha[3*dj + dk + ijk] - 
       alpha[-dj + 3*dk + ijk] + alpha[dj + 3*dk + ijk]) + 
    9.*(-alpha[-2*dj - 3*dk + ijk] + alpha[2*dj - 3*dk + ijk] - 
       alpha[-3*dj - 2*dk + ijk] + alpha[3*dj - 2*dk + ijk] + 
       alpha[-3*dj + 2*dk + ijk] - alpha[3*dj + 2*dk + ijk] + 
       alpha[-2*dj + 3*dk + ijk] - alpha[2*dj + 3*dk + ijk]) - 
    alpha[3*(-dj + dk) + ijk] + alpha[-3*(dj + dk) + ijk] + 
    81.*(-alpha[2*(dj - dk) + ijk] - alpha[2*(-dj + dk) + ijk] + 
       alpha[-2*(dj + dk) + ijk] + alpha[2*(dj + dk) + ijk]) + 
    alpha[3*(dj + dk) + ijk])
;

dda33
=
0.0055555555555555555556*oodz2*
  (-490.*alpha[ijk] + 270.*(alpha[-dk + ijk] + alpha[dk + ijk]) - 
    27.*(alpha[-2*dk + ijk] + alpha[2*dk + ijk]) + 
    2.*(alpha[-3*dk + ijk] + alpha[3*dk + ijk]))
;

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

ddb111
=
0.0055555555555555555556*oodx2*
  (-490.*beta1[ijk] + 270.*(beta1[-di + ijk] + beta1[di + ijk]) - 
    27.*(beta1[-2*di + ijk] + beta1[2*di + ijk]) + 
    2.*(beta1[-3*di + ijk] + beta1[3*di + ijk]))
;

ddb112
=
0.0055555555555555555556*oodx2*
  (-490.*beta2[ijk] + 270.*(beta2[-di + ijk] + beta2[di + ijk]) - 
    27.*(beta2[-2*di + ijk] + beta2[2*di + ijk]) + 
    2.*(beta2[-3*di + ijk] + beta2[3*di + ijk]))
;

ddb113
=
0.0055555555555555555556*oodx2*
  (-490.*beta3[ijk] + 270.*(beta3[-di + ijk] + beta3[di + ijk]) - 
    27.*(beta3[-2*di + ijk] + beta3[2*di + ijk]) + 
    2.*(beta3[-3*di + ijk] + beta3[3*di + ijk]))
;

ddb121
=
0.0011111111111111111111*oo2dx*oo2dy*
  (-beta1[3*(di - dj) + ijk] - 
    405.*(beta1[-di - 2*dj + ijk] + beta1[-2*di - dj + ijk]) + 
    2025.*(beta1[-di - dj + ijk] - beta1[di - dj + ijk] - 
       beta1[-di + dj + ijk] + beta1[di + dj + ijk]) + 
    405.*(beta1[di - 2*dj + ijk] + beta1[2*di - dj + ijk] + 
       beta1[-2*di + dj + ijk] - beta1[2*di + dj + ijk] + 
       beta1[-di + 2*dj + ijk] - beta1[di + 2*dj + ijk]) + 
    45.*(beta1[-di - 3*dj + ijk] - beta1[di - 3*dj + ijk] + 
       beta1[-3*di - dj + ijk] - beta1[3*di - dj + ijk] - 
       beta1[-3*di + dj + ijk] + beta1[3*di + dj + ijk] - 
       beta1[-di + 3*dj + ijk] + beta1[di + 3*dj + ijk]) + 
    9.*(-beta1[-2*di - 3*dj + ijk] + beta1[2*di - 3*dj + ijk] - 
       beta1[-3*di - 2*dj + ijk] + beta1[3*di - 2*dj + ijk] + 
       beta1[-3*di + 2*dj + ijk] - beta1[3*di + 2*dj + ijk] + 
       beta1[-2*di + 3*dj + ijk] - beta1[2*di + 3*dj + ijk]) - 
    beta1[3*(-di + dj) + ijk] + beta1[-3*(di + dj) + ijk] + 
    81.*(-beta1[2*(di - dj) + ijk] - beta1[2*(-di + dj) + ijk] + 
       beta1[-2*(di + dj) + ijk] + beta1[2*(di + dj) + ijk]) + 
    beta1[3*(di + dj) + ijk])
;

ddb122
=
0.0011111111111111111111*oo2dx*oo2dy*
  (-beta2[3*(di - dj) + ijk] - 
    405.*(beta2[-di - 2*dj + ijk] + beta2[-2*di - dj + ijk]) + 
    2025.*(beta2[-di - dj + ijk] - beta2[di - dj + ijk] - 
       beta2[-di + dj + ijk] + beta2[di + dj + ijk]) + 
    405.*(beta2[di - 2*dj + ijk] + beta2[2*di - dj + ijk] + 
       beta2[-2*di + dj + ijk] - beta2[2*di + dj + ijk] + 
       beta2[-di + 2*dj + ijk] - beta2[di + 2*dj + ijk]) + 
    45.*(beta2[-di - 3*dj + ijk] - beta2[di - 3*dj + ijk] + 
       beta2[-3*di - dj + ijk] - beta2[3*di - dj + ijk] - 
       beta2[-3*di + dj + ijk] + beta2[3*di + dj + ijk] - 
       beta2[-di + 3*dj + ijk] + beta2[di + 3*dj + ijk]) + 
    9.*(-beta2[-2*di - 3*dj + ijk] + beta2[2*di - 3*dj + ijk] - 
       beta2[-3*di - 2*dj + ijk] + beta2[3*di - 2*dj + ijk] + 
       beta2[-3*di + 2*dj + ijk] - beta2[3*di + 2*dj + ijk] + 
       beta2[-2*di + 3*dj + ijk] - beta2[2*di + 3*dj + ijk]) - 
    beta2[3*(-di + dj) + ijk] + beta2[-3*(di + dj) + ijk] + 
    81.*(-beta2[2*(di - dj) + ijk] - beta2[2*(-di + dj) + ijk] + 
       beta2[-2*(di + dj) + ijk] + beta2[2*(di + dj) + ijk]) + 
    beta2[3*(di + dj) + ijk])
;

ddb123
=
0.0011111111111111111111*oo2dx*oo2dy*
  (-beta3[3*(di - dj) + ijk] - 
    405.*(beta3[-di - 2*dj + ijk] + beta3[-2*di - dj + ijk]) + 
    2025.*(beta3[-di - dj + ijk] - beta3[di - dj + ijk] - 
       beta3[-di + dj + ijk] + beta3[di + dj + ijk]) + 
    405.*(beta3[di - 2*dj + ijk] + beta3[2*di - dj + ijk] + 
       beta3[-2*di + dj + ijk] - beta3[2*di + dj + ijk] + 
       beta3[-di + 2*dj + ijk] - beta3[di + 2*dj + ijk]) + 
    45.*(beta3[-di - 3*dj + ijk] - beta3[di - 3*dj + ijk] + 
       beta3[-3*di - dj + ijk] - beta3[3*di - dj + ijk] - 
       beta3[-3*di + dj + ijk] + beta3[3*di + dj + ijk] - 
       beta3[-di + 3*dj + ijk] + beta3[di + 3*dj + ijk]) + 
    9.*(-beta3[-2*di - 3*dj + ijk] + beta3[2*di - 3*dj + ijk] - 
       beta3[-3*di - 2*dj + ijk] + beta3[3*di - 2*dj + ijk] + 
       beta3[-3*di + 2*dj + ijk] - beta3[3*di + 2*dj + ijk] + 
       beta3[-2*di + 3*dj + ijk] - beta3[2*di + 3*dj + ijk]) - 
    beta3[3*(-di + dj) + ijk] + beta3[-3*(di + dj) + ijk] + 
    81.*(-beta3[2*(di - dj) + ijk] - beta3[2*(-di + dj) + ijk] + 
       beta3[-2*(di + dj) + ijk] + beta3[2*(di + dj) + ijk]) + 
    beta3[3*(di + dj) + ijk])
;

ddb131
=
0.0011111111111111111111*oo2dx*oo2dz*
  (-beta1[3*(di - dk) + ijk] - 
    405.*(beta1[-di - 2*dk + ijk] + beta1[-2*di - dk + ijk]) + 
    2025.*(beta1[-di - dk + ijk] - beta1[di - dk + ijk] - 
       beta1[-di + dk + ijk] + beta1[di + dk + ijk]) + 
    405.*(beta1[di - 2*dk + ijk] + beta1[2*di - dk + ijk] + 
       beta1[-2*di + dk + ijk] - beta1[2*di + dk + ijk] + 
       beta1[-di + 2*dk + ijk] - beta1[di + 2*dk + ijk]) + 
    45.*(beta1[-di - 3*dk + ijk] - beta1[di - 3*dk + ijk] + 
       beta1[-3*di - dk + ijk] - beta1[3*di - dk + ijk] - 
       beta1[-3*di + dk + ijk] + beta1[3*di + dk + ijk] - 
       beta1[-di + 3*dk + ijk] + beta1[di + 3*dk + ijk]) + 
    9.*(-beta1[-2*di - 3*dk + ijk] + beta1[2*di - 3*dk + ijk] - 
       beta1[-3*di - 2*dk + ijk] + beta1[3*di - 2*dk + ijk] + 
       beta1[-3*di + 2*dk + ijk] - beta1[3*di + 2*dk + ijk] + 
       beta1[-2*di + 3*dk + ijk] - beta1[2*di + 3*dk + ijk]) - 
    beta1[3*(-di + dk) + ijk] + beta1[-3*(di + dk) + ijk] + 
    81.*(-beta1[2*(di - dk) + ijk] - beta1[2*(-di + dk) + ijk] + 
       beta1[-2*(di + dk) + ijk] + beta1[2*(di + dk) + ijk]) + 
    beta1[3*(di + dk) + ijk])
;

ddb132
=
0.0011111111111111111111*oo2dx*oo2dz*
  (-beta2[3*(di - dk) + ijk] - 
    405.*(beta2[-di - 2*dk + ijk] + beta2[-2*di - dk + ijk]) + 
    2025.*(beta2[-di - dk + ijk] - beta2[di - dk + ijk] - 
       beta2[-di + dk + ijk] + beta2[di + dk + ijk]) + 
    405.*(beta2[di - 2*dk + ijk] + beta2[2*di - dk + ijk] + 
       beta2[-2*di + dk + ijk] - beta2[2*di + dk + ijk] + 
       beta2[-di + 2*dk + ijk] - beta2[di + 2*dk + ijk]) + 
    45.*(beta2[-di - 3*dk + ijk] - beta2[di - 3*dk + ijk] + 
       beta2[-3*di - dk + ijk] - beta2[3*di - dk + ijk] - 
       beta2[-3*di + dk + ijk] + beta2[3*di + dk + ijk] - 
       beta2[-di + 3*dk + ijk] + beta2[di + 3*dk + ijk]) + 
    9.*(-beta2[-2*di - 3*dk + ijk] + beta2[2*di - 3*dk + ijk] - 
       beta2[-3*di - 2*dk + ijk] + beta2[3*di - 2*dk + ijk] + 
       beta2[-3*di + 2*dk + ijk] - beta2[3*di + 2*dk + ijk] + 
       beta2[-2*di + 3*dk + ijk] - beta2[2*di + 3*dk + ijk]) - 
    beta2[3*(-di + dk) + ijk] + beta2[-3*(di + dk) + ijk] + 
    81.*(-beta2[2*(di - dk) + ijk] - beta2[2*(-di + dk) + ijk] + 
       beta2[-2*(di + dk) + ijk] + beta2[2*(di + dk) + ijk]) + 
    beta2[3*(di + dk) + ijk])
;

ddb133
=
0.0011111111111111111111*oo2dx*oo2dz*
  (-beta3[3*(di - dk) + ijk] - 
    405.*(beta3[-di - 2*dk + ijk] + beta3[-2*di - dk + ijk]) + 
    2025.*(beta3[-di - dk + ijk] - beta3[di - dk + ijk] - 
       beta3[-di + dk + ijk] + beta3[di + dk + ijk]) + 
    405.*(beta3[di - 2*dk + ijk] + beta3[2*di - dk + ijk] + 
       beta3[-2*di + dk + ijk] - beta3[2*di + dk + ijk] + 
       beta3[-di + 2*dk + ijk] - beta3[di + 2*dk + ijk]) + 
    45.*(beta3[-di - 3*dk + ijk] - beta3[di - 3*dk + ijk] + 
       beta3[-3*di - dk + ijk] - beta3[3*di - dk + ijk] - 
       beta3[-3*di + dk + ijk] + beta3[3*di + dk + ijk] - 
       beta3[-di + 3*dk + ijk] + beta3[di + 3*dk + ijk]) + 
    9.*(-beta3[-2*di - 3*dk + ijk] + beta3[2*di - 3*dk + ijk] - 
       beta3[-3*di - 2*dk + ijk] + beta3[3*di - 2*dk + ijk] + 
       beta3[-3*di + 2*dk + ijk] - beta3[3*di + 2*dk + ijk] + 
       beta3[-2*di + 3*dk + ijk] - beta3[2*di + 3*dk + ijk]) - 
    beta3[3*(-di + dk) + ijk] + beta3[-3*(di + dk) + ijk] + 
    81.*(-beta3[2*(di - dk) + ijk] - beta3[2*(-di + dk) + ijk] + 
       beta3[-2*(di + dk) + ijk] + beta3[2*(di + dk) + ijk]) + 
    beta3[3*(di + dk) + ijk])
;

ddb221
=
0.0055555555555555555556*oody2*
  (-490.*beta1[ijk] + 270.*(beta1[-dj + ijk] + beta1[dj + ijk]) - 
    27.*(beta1[-2*dj + ijk] + beta1[2*dj + ijk]) + 
    2.*(beta1[-3*dj + ijk] + beta1[3*dj + ijk]))
;

ddb222
=
0.0055555555555555555556*oody2*
  (-490.*beta2[ijk] + 270.*(beta2[-dj + ijk] + beta2[dj + ijk]) - 
    27.*(beta2[-2*dj + ijk] + beta2[2*dj + ijk]) + 
    2.*(beta2[-3*dj + ijk] + beta2[3*dj + ijk]))
;

ddb223
=
0.0055555555555555555556*oody2*
  (-490.*beta3[ijk] + 270.*(beta3[-dj + ijk] + beta3[dj + ijk]) - 
    27.*(beta3[-2*dj + ijk] + beta3[2*dj + ijk]) + 
    2.*(beta3[-3*dj + ijk] + beta3[3*dj + ijk]))
;

ddb231
=
0.0011111111111111111111*oo2dy*oo2dz*
  (-beta1[3*(dj - dk) + ijk] - 
    405.*(beta1[-dj - 2*dk + ijk] + beta1[-2*dj - dk + ijk]) + 
    2025.*(beta1[-dj - dk + ijk] - beta1[dj - dk + ijk] - 
       beta1[-dj + dk + ijk] + beta1[dj + dk + ijk]) + 
    405.*(beta1[dj - 2*dk + ijk] + beta1[2*dj - dk + ijk] + 
       beta1[-2*dj + dk + ijk] - beta1[2*dj + dk + ijk] + 
       beta1[-dj + 2*dk + ijk] - beta1[dj + 2*dk + ijk]) + 
    45.*(beta1[-dj - 3*dk + ijk] - beta1[dj - 3*dk + ijk] + 
       beta1[-3*dj - dk + ijk] - beta1[3*dj - dk + ijk] - 
       beta1[-3*dj + dk + ijk] + beta1[3*dj + dk + ijk] - 
       beta1[-dj + 3*dk + ijk] + beta1[dj + 3*dk + ijk]) + 
    9.*(-beta1[-2*dj - 3*dk + ijk] + beta1[2*dj - 3*dk + ijk] - 
       beta1[-3*dj - 2*dk + ijk] + beta1[3*dj - 2*dk + ijk] + 
       beta1[-3*dj + 2*dk + ijk] - beta1[3*dj + 2*dk + ijk] + 
       beta1[-2*dj + 3*dk + ijk] - beta1[2*dj + 3*dk + ijk]) - 
    beta1[3*(-dj + dk) + ijk] + beta1[-3*(dj + dk) + ijk] + 
    81.*(-beta1[2*(dj - dk) + ijk] - beta1[2*(-dj + dk) + ijk] + 
       beta1[-2*(dj + dk) + ijk] + beta1[2*(dj + dk) + ijk]) + 
    beta1[3*(dj + dk) + ijk])
;

ddb232
=
0.0011111111111111111111*oo2dy*oo2dz*
  (-beta2[3*(dj - dk) + ijk] - 
    405.*(beta2[-dj - 2*dk + ijk] + beta2[-2*dj - dk + ijk]) + 
    2025.*(beta2[-dj - dk + ijk] - beta2[dj - dk + ijk] - 
       beta2[-dj + dk + ijk] + beta2[dj + dk + ijk]) + 
    405.*(beta2[dj - 2*dk + ijk] + beta2[2*dj - dk + ijk] + 
       beta2[-2*dj + dk + ijk] - beta2[2*dj + dk + ijk] + 
       beta2[-dj + 2*dk + ijk] - beta2[dj + 2*dk + ijk]) + 
    45.*(beta2[-dj - 3*dk + ijk] - beta2[dj - 3*dk + ijk] + 
       beta2[-3*dj - dk + ijk] - beta2[3*dj - dk + ijk] - 
       beta2[-3*dj + dk + ijk] + beta2[3*dj + dk + ijk] - 
       beta2[-dj + 3*dk + ijk] + beta2[dj + 3*dk + ijk]) + 
    9.*(-beta2[-2*dj - 3*dk + ijk] + beta2[2*dj - 3*dk + ijk] - 
       beta2[-3*dj - 2*dk + ijk] + beta2[3*dj - 2*dk + ijk] + 
       beta2[-3*dj + 2*dk + ijk] - beta2[3*dj + 2*dk + ijk] + 
       beta2[-2*dj + 3*dk + ijk] - beta2[2*dj + 3*dk + ijk]) - 
    beta2[3*(-dj + dk) + ijk] + beta2[-3*(dj + dk) + ijk] + 
    81.*(-beta2[2*(dj - dk) + ijk] - beta2[2*(-dj + dk) + ijk] + 
       beta2[-2*(dj + dk) + ijk] + beta2[2*(dj + dk) + ijk]) + 
    beta2[3*(dj + dk) + ijk])
;

ddb233
=
0.0011111111111111111111*oo2dy*oo2dz*
  (-beta3[3*(dj - dk) + ijk] - 
    405.*(beta3[-dj - 2*dk + ijk] + beta3[-2*dj - dk + ijk]) + 
    2025.*(beta3[-dj - dk + ijk] - beta3[dj - dk + ijk] - 
       beta3[-dj + dk + ijk] + beta3[dj + dk + ijk]) + 
    405.*(beta3[dj - 2*dk + ijk] + beta3[2*dj - dk + ijk] + 
       beta3[-2*dj + dk + ijk] - beta3[2*dj + dk + ijk] + 
       beta3[-dj + 2*dk + ijk] - beta3[dj + 2*dk + ijk]) + 
    45.*(beta3[-dj - 3*dk + ijk] - beta3[dj - 3*dk + ijk] + 
       beta3[-3*dj - dk + ijk] - beta3[3*dj - dk + ijk] - 
       beta3[-3*dj + dk + ijk] + beta3[3*dj + dk + ijk] - 
       beta3[-dj + 3*dk + ijk] + beta3[dj + 3*dk + ijk]) + 
    9.*(-beta3[-2*dj - 3*dk + ijk] + beta3[2*dj - 3*dk + ijk] - 
       beta3[-3*dj - 2*dk + ijk] + beta3[3*dj - 2*dk + ijk] + 
       beta3[-3*dj + 2*dk + ijk] - beta3[3*dj + 2*dk + ijk] + 
       beta3[-2*dj + 3*dk + ijk] - beta3[2*dj + 3*dk + ijk]) - 
    beta3[3*(-dj + dk) + ijk] + beta3[-3*(dj + dk) + ijk] + 
    81.*(-beta3[2*(dj - dk) + ijk] - beta3[2*(-dj + dk) + ijk] + 
       beta3[-2*(dj + dk) + ijk] + beta3[2*(dj + dk) + ijk]) + 
    beta3[3*(dj + dk) + ijk])
;

ddb331
=
0.0055555555555555555556*oodz2*
  (-490.*beta1[ijk] + 270.*(beta1[-dk + ijk] + beta1[dk + ijk]) - 
    27.*(beta1[-2*dk + ijk] + beta1[2*dk + ijk]) + 
    2.*(beta1[-3*dk + ijk] + beta1[3*dk + ijk]))
;

ddb332
=
0.0055555555555555555556*oodz2*
  (-490.*beta2[ijk] + 270.*(beta2[-dk + ijk] + beta2[dk + ijk]) - 
    27.*(beta2[-2*dk + ijk] + beta2[2*dk + ijk]) + 
    2.*(beta2[-3*dk + ijk] + beta2[3*dk + ijk]))
;

ddb333
=
0.0055555555555555555556*oodz2*
  (-490.*beta3[ijk] + 270.*(beta3[-dk + ijk] + beta3[dk + ijk]) - 
    27.*(beta3[-2*dk + ijk] + beta3[2*dk + ijk]) + 
    2.*(beta3[-3*dk + ijk] + beta3[3*dk + ijk]))
;

delg111
=
0.033333333333333333333*oo2dx*(-g11[-3*di + ijk] + 
    45.*(-g11[-di + ijk] + g11[di + ijk]) + 
    9.*(g11[-2*di + ijk] - g11[2*di + ijk]) + g11[3*di + ijk])
;

delg112
=
0.033333333333333333333*oo2dx*(-g12[-3*di + ijk] + 
    45.*(-g12[-di + ijk] + g12[di + ijk]) + 
    9.*(g12[-2*di + ijk] - g12[2*di + ijk]) + g12[3*di + ijk])
;

delg113
=
0.033333333333333333333*oo2dx*(-g13[-3*di + ijk] + 
    45.*(-g13[-di + ijk] + g13[di + ijk]) + 
    9.*(g13[-2*di + ijk] - g13[2*di + ijk]) + g13[3*di + ijk])
;

delg122
=
0.033333333333333333333*oo2dx*(-g22[-3*di + ijk] + 
    45.*(-g22[-di + ijk] + g22[di + ijk]) + 
    9.*(g22[-2*di + ijk] - g22[2*di + ijk]) + g22[3*di + ijk])
;

delg123
=
0.033333333333333333333*oo2dx*(-g23[-3*di + ijk] + 
    45.*(-g23[-di + ijk] + g23[di + ijk]) + 
    9.*(g23[-2*di + ijk] - g23[2*di + ijk]) + g23[3*di + ijk])
;

delg133
=
0.033333333333333333333*oo2dx*(-g33[-3*di + ijk] + 
    45.*(-g33[-di + ijk] + g33[di + ijk]) + 
    9.*(g33[-2*di + ijk] - g33[2*di + ijk]) + g33[3*di + ijk])
;

delg211
=
0.033333333333333333333*oo2dy*(-g11[-3*dj + ijk] + 
    45.*(-g11[-dj + ijk] + g11[dj + ijk]) + 
    9.*(g11[-2*dj + ijk] - g11[2*dj + ijk]) + g11[3*dj + ijk])
;

delg212
=
0.033333333333333333333*oo2dy*(-g12[-3*dj + ijk] + 
    45.*(-g12[-dj + ijk] + g12[dj + ijk]) + 
    9.*(g12[-2*dj + ijk] - g12[2*dj + ijk]) + g12[3*dj + ijk])
;

delg213
=
0.033333333333333333333*oo2dy*(-g13[-3*dj + ijk] + 
    45.*(-g13[-dj + ijk] + g13[dj + ijk]) + 
    9.*(g13[-2*dj + ijk] - g13[2*dj + ijk]) + g13[3*dj + ijk])
;

delg222
=
0.033333333333333333333*oo2dy*(-g22[-3*dj + ijk] + 
    45.*(-g22[-dj + ijk] + g22[dj + ijk]) + 
    9.*(g22[-2*dj + ijk] - g22[2*dj + ijk]) + g22[3*dj + ijk])
;

delg223
=
0.033333333333333333333*oo2dy*(-g23[-3*dj + ijk] + 
    45.*(-g23[-dj + ijk] + g23[dj + ijk]) + 
    9.*(g23[-2*dj + ijk] - g23[2*dj + ijk]) + g23[3*dj + ijk])
;

delg233
=
0.033333333333333333333*oo2dy*(-g33[-3*dj + ijk] + 
    45.*(-g33[-dj + ijk] + g33[dj + ijk]) + 
    9.*(g33[-2*dj + ijk] - g33[2*dj + ijk]) + g33[3*dj + ijk])
;

delg311
=
0.033333333333333333333*oo2dz*(-g11[-3*dk + ijk] + 
    45.*(-g11[-dk + ijk] + g11[dk + ijk]) + 
    9.*(g11[-2*dk + ijk] - g11[2*dk + ijk]) + g11[3*dk + ijk])
;

delg312
=
0.033333333333333333333*oo2dz*(-g12[-3*dk + ijk] + 
    45.*(-g12[-dk + ijk] + g12[dk + ijk]) + 
    9.*(g12[-2*dk + ijk] - g12[2*dk + ijk]) + g12[3*dk + ijk])
;

delg313
=
0.033333333333333333333*oo2dz*(-g13[-3*dk + ijk] + 
    45.*(-g13[-dk + ijk] + g13[dk + ijk]) + 
    9.*(g13[-2*dk + ijk] - g13[2*dk + ijk]) + g13[3*dk + ijk])
;

delg322
=
0.033333333333333333333*oo2dz*(-g22[-3*dk + ijk] + 
    45.*(-g22[-dk + ijk] + g22[dk + ijk]) + 
    9.*(g22[-2*dk + ijk] - g22[2*dk + ijk]) + g22[3*dk + ijk])
;

delg323
=
0.033333333333333333333*oo2dz*(-g23[-3*dk + ijk] + 
    45.*(-g23[-dk + ijk] + g23[dk + ijk]) + 
    9.*(g23[-2*dk + ijk] - g23[2*dk + ijk]) + g23[3*dk + ijk])
;

delg333
=
0.033333333333333333333*oo2dz*(-g33[-3*dk + ijk] + 
    45.*(-g33[-dk + ijk] + g33[dk + ijk]) + 
    9.*(g33[-2*dk + ijk] - g33[2*dk + ijk]) + g33[3*dk + ijk])
;

deldelg1111
=
0.0055555555555555555556*oodx2*
  (-490.*g11[ijk] + 270.*(g11[-di + ijk] + g11[di + ijk]) - 
    27.*(g11[-2*di + ijk] + g11[2*di + ijk]) + 
    2.*(g11[-3*di + ijk] + g11[3*di + ijk]))
;

deldelg1112
=
0.0055555555555555555556*oodx2*
  (-490.*g12[ijk] + 270.*(g12[-di + ijk] + g12[di + ijk]) - 
    27.*(g12[-2*di + ijk] + g12[2*di + ijk]) + 
    2.*(g12[-3*di + ijk] + g12[3*di + ijk]))
;

deldelg1113
=
0.0055555555555555555556*oodx2*
  (-490.*g13[ijk] + 270.*(g13[-di + ijk] + g13[di + ijk]) - 
    27.*(g13[-2*di + ijk] + g13[2*di + ijk]) + 
    2.*(g13[-3*di + ijk] + g13[3*di + ijk]))
;

deldelg1122
=
0.0055555555555555555556*oodx2*
  (-490.*g22[ijk] + 270.*(g22[-di + ijk] + g22[di + ijk]) - 
    27.*(g22[-2*di + ijk] + g22[2*di + ijk]) + 
    2.*(g22[-3*di + ijk] + g22[3*di + ijk]))
;

deldelg1123
=
0.0055555555555555555556*oodx2*
  (-490.*g23[ijk] + 270.*(g23[-di + ijk] + g23[di + ijk]) - 
    27.*(g23[-2*di + ijk] + g23[2*di + ijk]) + 
    2.*(g23[-3*di + ijk] + g23[3*di + ijk]))
;

deldelg1133
=
0.0055555555555555555556*oodx2*
  (-490.*g33[ijk] + 270.*(g33[-di + ijk] + g33[di + ijk]) - 
    27.*(g33[-2*di + ijk] + g33[2*di + ijk]) + 
    2.*(g33[-3*di + ijk] + g33[3*di + ijk]))
;

deldelg1211
=
0.0011111111111111111111*oo2dx*oo2dy*
  (-g11[3*(di - dj) + ijk] - 405.*
     (g11[-di - 2*dj + ijk] + g11[-2*di - dj + ijk]) + 
    2025.*(g11[-di - dj + ijk] - g11[di - dj + ijk] - g11[-di + dj + ijk] + 
       g11[di + dj + ijk]) + 405.*
     (g11[di - 2*dj + ijk] + g11[2*di - dj + ijk] + g11[-2*di + dj + ijk] - 
       g11[2*di + dj + ijk] + g11[-di + 2*dj + ijk] - g11[di + 2*dj + ijk]) \
+ 45.*(g11[-di - 3*dj + ijk] - g11[di - 3*dj + ijk] + 
       g11[-3*di - dj + ijk] - g11[3*di - dj + ijk] - 
       g11[-3*di + dj + ijk] + g11[3*di + dj + ijk] - 
       g11[-di + 3*dj + ijk] + g11[di + 3*dj + ijk]) + 
    9.*(-g11[-2*di - 3*dj + ijk] + g11[2*di - 3*dj + ijk] - 
       g11[-3*di - 2*dj + ijk] + g11[3*di - 2*dj + ijk] + 
       g11[-3*di + 2*dj + ijk] - g11[3*di + 2*dj + ijk] + 
       g11[-2*di + 3*dj + ijk] - g11[2*di + 3*dj + ijk]) - 
    g11[3*(-di + dj) + ijk] + g11[-3*(di + dj) + ijk] + 
    81.*(-g11[2*(di - dj) + ijk] - g11[2*(-di + dj) + ijk] + 
       g11[-2*(di + dj) + ijk] + g11[2*(di + dj) + ijk]) + 
    g11[3*(di + dj) + ijk])
;

deldelg1212
=
0.0011111111111111111111*oo2dx*oo2dy*
  (-g12[3*(di - dj) + ijk] - 405.*
     (g12[-di - 2*dj + ijk] + g12[-2*di - dj + ijk]) + 
    2025.*(g12[-di - dj + ijk] - g12[di - dj + ijk] - g12[-di + dj + ijk] + 
       g12[di + dj + ijk]) + 405.*
     (g12[di - 2*dj + ijk] + g12[2*di - dj + ijk] + g12[-2*di + dj + ijk] - 
       g12[2*di + dj + ijk] + g12[-di + 2*dj + ijk] - g12[di + 2*dj + ijk]) \
+ 45.*(g12[-di - 3*dj + ijk] - g12[di - 3*dj + ijk] + 
       g12[-3*di - dj + ijk] - g12[3*di - dj + ijk] - 
       g12[-3*di + dj + ijk] + g12[3*di + dj + ijk] - 
       g12[-di + 3*dj + ijk] + g12[di + 3*dj + ijk]) + 
    9.*(-g12[-2*di - 3*dj + ijk] + g12[2*di - 3*dj + ijk] - 
       g12[-3*di - 2*dj + ijk] + g12[3*di - 2*dj + ijk] + 
       g12[-3*di + 2*dj + ijk] - g12[3*di + 2*dj + ijk] + 
       g12[-2*di + 3*dj + ijk] - g12[2*di + 3*dj + ijk]) - 
    g12[3*(-di + dj) + ijk] + g12[-3*(di + dj) + ijk] + 
    81.*(-g12[2*(di - dj) + ijk] - g12[2*(-di + dj) + ijk] + 
       g12[-2*(di + dj) + ijk] + g12[2*(di + dj) + ijk]) + 
    g12[3*(di + dj) + ijk])
;

deldelg1213
=
0.0011111111111111111111*oo2dx*oo2dy*
  (-g13[3*(di - dj) + ijk] - 405.*
     (g13[-di - 2*dj + ijk] + g13[-2*di - dj + ijk]) + 
    2025.*(g13[-di - dj + ijk] - g13[di - dj + ijk] - g13[-di + dj + ijk] + 
       g13[di + dj + ijk]) + 405.*
     (g13[di - 2*dj + ijk] + g13[2*di - dj + ijk] + g13[-2*di + dj + ijk] - 
       g13[2*di + dj + ijk] + g13[-di + 2*dj + ijk] - g13[di + 2*dj + ijk]) \
+ 45.*(g13[-di - 3*dj + ijk] - g13[di - 3*dj + ijk] + 
       g13[-3*di - dj + ijk] - g13[3*di - dj + ijk] - 
       g13[-3*di + dj + ijk] + g13[3*di + dj + ijk] - 
       g13[-di + 3*dj + ijk] + g13[di + 3*dj + ijk]) + 
    9.*(-g13[-2*di - 3*dj + ijk] + g13[2*di - 3*dj + ijk] - 
       g13[-3*di - 2*dj + ijk] + g13[3*di - 2*dj + ijk] + 
       g13[-3*di + 2*dj + ijk] - g13[3*di + 2*dj + ijk] + 
       g13[-2*di + 3*dj + ijk] - g13[2*di + 3*dj + ijk]) - 
    g13[3*(-di + dj) + ijk] + g13[-3*(di + dj) + ijk] + 
    81.*(-g13[2*(di - dj) + ijk] - g13[2*(-di + dj) + ijk] + 
       g13[-2*(di + dj) + ijk] + g13[2*(di + dj) + ijk]) + 
    g13[3*(di + dj) + ijk])
;

deldelg1222
=
0.0011111111111111111111*oo2dx*oo2dy*
  (-g22[3*(di - dj) + ijk] - 405.*
     (g22[-di - 2*dj + ijk] + g22[-2*di - dj + ijk]) + 
    2025.*(g22[-di - dj + ijk] - g22[di - dj + ijk] - g22[-di + dj + ijk] + 
       g22[di + dj + ijk]) + 405.*
     (g22[di - 2*dj + ijk] + g22[2*di - dj + ijk] + g22[-2*di + dj + ijk] - 
       g22[2*di + dj + ijk] + g22[-di + 2*dj + ijk] - g22[di + 2*dj + ijk]) \
+ 45.*(g22[-di - 3*dj + ijk] - g22[di - 3*dj + ijk] + 
       g22[-3*di - dj + ijk] - g22[3*di - dj + ijk] - 
       g22[-3*di + dj + ijk] + g22[3*di + dj + ijk] - 
       g22[-di + 3*dj + ijk] + g22[di + 3*dj + ijk]) + 
    9.*(-g22[-2*di - 3*dj + ijk] + g22[2*di - 3*dj + ijk] - 
       g22[-3*di - 2*dj + ijk] + g22[3*di - 2*dj + ijk] + 
       g22[-3*di + 2*dj + ijk] - g22[3*di + 2*dj + ijk] + 
       g22[-2*di + 3*dj + ijk] - g22[2*di + 3*dj + ijk]) - 
    g22[3*(-di + dj) + ijk] + g22[-3*(di + dj) + ijk] + 
    81.*(-g22[2*(di - dj) + ijk] - g22[2*(-di + dj) + ijk] + 
       g22[-2*(di + dj) + ijk] + g22[2*(di + dj) + ijk]) + 
    g22[3*(di + dj) + ijk])
;

deldelg1223
=
0.0011111111111111111111*oo2dx*oo2dy*
  (-g23[3*(di - dj) + ijk] - 405.*
     (g23[-di - 2*dj + ijk] + g23[-2*di - dj + ijk]) + 
    2025.*(g23[-di - dj + ijk] - g23[di - dj + ijk] - g23[-di + dj + ijk] + 
       g23[di + dj + ijk]) + 405.*
     (g23[di - 2*dj + ijk] + g23[2*di - dj + ijk] + g23[-2*di + dj + ijk] - 
       g23[2*di + dj + ijk] + g23[-di + 2*dj + ijk] - g23[di + 2*dj + ijk]) \
+ 45.*(g23[-di - 3*dj + ijk] - g23[di - 3*dj + ijk] + 
       g23[-3*di - dj + ijk] - g23[3*di - dj + ijk] - 
       g23[-3*di + dj + ijk] + g23[3*di + dj + ijk] - 
       g23[-di + 3*dj + ijk] + g23[di + 3*dj + ijk]) + 
    9.*(-g23[-2*di - 3*dj + ijk] + g23[2*di - 3*dj + ijk] - 
       g23[-3*di - 2*dj + ijk] + g23[3*di - 2*dj + ijk] + 
       g23[-3*di + 2*dj + ijk] - g23[3*di + 2*dj + ijk] + 
       g23[-2*di + 3*dj + ijk] - g23[2*di + 3*dj + ijk]) - 
    g23[3*(-di + dj) + ijk] + g23[-3*(di + dj) + ijk] + 
    81.*(-g23[2*(di - dj) + ijk] - g23[2*(-di + dj) + ijk] + 
       g23[-2*(di + dj) + ijk] + g23[2*(di + dj) + ijk]) + 
    g23[3*(di + dj) + ijk])
;

deldelg1233
=
0.0011111111111111111111*oo2dx*oo2dy*
  (-g33[3*(di - dj) + ijk] - 405.*
     (g33[-di - 2*dj + ijk] + g33[-2*di - dj + ijk]) + 
    2025.*(g33[-di - dj + ijk] - g33[di - dj + ijk] - g33[-di + dj + ijk] + 
       g33[di + dj + ijk]) + 405.*
     (g33[di - 2*dj + ijk] + g33[2*di - dj + ijk] + g33[-2*di + dj + ijk] - 
       g33[2*di + dj + ijk] + g33[-di + 2*dj + ijk] - g33[di + 2*dj + ijk]) \
+ 45.*(g33[-di - 3*dj + ijk] - g33[di - 3*dj + ijk] + 
       g33[-3*di - dj + ijk] - g33[3*di - dj + ijk] - 
       g33[-3*di + dj + ijk] + g33[3*di + dj + ijk] - 
       g33[-di + 3*dj + ijk] + g33[di + 3*dj + ijk]) + 
    9.*(-g33[-2*di - 3*dj + ijk] + g33[2*di - 3*dj + ijk] - 
       g33[-3*di - 2*dj + ijk] + g33[3*di - 2*dj + ijk] + 
       g33[-3*di + 2*dj + ijk] - g33[3*di + 2*dj + ijk] + 
       g33[-2*di + 3*dj + ijk] - g33[2*di + 3*dj + ijk]) - 
    g33[3*(-di + dj) + ijk] + g33[-3*(di + dj) + ijk] + 
    81.*(-g33[2*(di - dj) + ijk] - g33[2*(-di + dj) + ijk] + 
       g33[-2*(di + dj) + ijk] + g33[2*(di + dj) + ijk]) + 
    g33[3*(di + dj) + ijk])
;

deldelg1311
=
0.0011111111111111111111*oo2dx*oo2dz*
  (-g11[3*(di - dk) + ijk] - 405.*
     (g11[-di - 2*dk + ijk] + g11[-2*di - dk + ijk]) + 
    2025.*(g11[-di - dk + ijk] - g11[di - dk + ijk] - g11[-di + dk + ijk] + 
       g11[di + dk + ijk]) + 405.*
     (g11[di - 2*dk + ijk] + g11[2*di - dk + ijk] + g11[-2*di + dk + ijk] - 
       g11[2*di + dk + ijk] + g11[-di + 2*dk + ijk] - g11[di + 2*dk + ijk]) \
+ 45.*(g11[-di - 3*dk + ijk] - g11[di - 3*dk + ijk] + 
       g11[-3*di - dk + ijk] - g11[3*di - dk + ijk] - 
       g11[-3*di + dk + ijk] + g11[3*di + dk + ijk] - 
       g11[-di + 3*dk + ijk] + g11[di + 3*dk + ijk]) + 
    9.*(-g11[-2*di - 3*dk + ijk] + g11[2*di - 3*dk + ijk] - 
       g11[-3*di - 2*dk + ijk] + g11[3*di - 2*dk + ijk] + 
       g11[-3*di + 2*dk + ijk] - g11[3*di + 2*dk + ijk] + 
       g11[-2*di + 3*dk + ijk] - g11[2*di + 3*dk + ijk]) - 
    g11[3*(-di + dk) + ijk] + g11[-3*(di + dk) + ijk] + 
    81.*(-g11[2*(di - dk) + ijk] - g11[2*(-di + dk) + ijk] + 
       g11[-2*(di + dk) + ijk] + g11[2*(di + dk) + ijk]) + 
    g11[3*(di + dk) + ijk])
;

deldelg1312
=
0.0011111111111111111111*oo2dx*oo2dz*
  (-g12[3*(di - dk) + ijk] - 405.*
     (g12[-di - 2*dk + ijk] + g12[-2*di - dk + ijk]) + 
    2025.*(g12[-di - dk + ijk] - g12[di - dk + ijk] - g12[-di + dk + ijk] + 
       g12[di + dk + ijk]) + 405.*
     (g12[di - 2*dk + ijk] + g12[2*di - dk + ijk] + g12[-2*di + dk + ijk] - 
       g12[2*di + dk + ijk] + g12[-di + 2*dk + ijk] - g12[di + 2*dk + ijk]) \
+ 45.*(g12[-di - 3*dk + ijk] - g12[di - 3*dk + ijk] + 
       g12[-3*di - dk + ijk] - g12[3*di - dk + ijk] - 
       g12[-3*di + dk + ijk] + g12[3*di + dk + ijk] - 
       g12[-di + 3*dk + ijk] + g12[di + 3*dk + ijk]) + 
    9.*(-g12[-2*di - 3*dk + ijk] + g12[2*di - 3*dk + ijk] - 
       g12[-3*di - 2*dk + ijk] + g12[3*di - 2*dk + ijk] + 
       g12[-3*di + 2*dk + ijk] - g12[3*di + 2*dk + ijk] + 
       g12[-2*di + 3*dk + ijk] - g12[2*di + 3*dk + ijk]) - 
    g12[3*(-di + dk) + ijk] + g12[-3*(di + dk) + ijk] + 
    81.*(-g12[2*(di - dk) + ijk] - g12[2*(-di + dk) + ijk] + 
       g12[-2*(di + dk) + ijk] + g12[2*(di + dk) + ijk]) + 
    g12[3*(di + dk) + ijk])
;

deldelg1313
=
0.0011111111111111111111*oo2dx*oo2dz*
  (-g13[3*(di - dk) + ijk] - 405.*
     (g13[-di - 2*dk + ijk] + g13[-2*di - dk + ijk]) + 
    2025.*(g13[-di - dk + ijk] - g13[di - dk + ijk] - g13[-di + dk + ijk] + 
       g13[di + dk + ijk]) + 405.*
     (g13[di - 2*dk + ijk] + g13[2*di - dk + ijk] + g13[-2*di + dk + ijk] - 
       g13[2*di + dk + ijk] + g13[-di + 2*dk + ijk] - g13[di + 2*dk + ijk]) \
+ 45.*(g13[-di - 3*dk + ijk] - g13[di - 3*dk + ijk] + 
       g13[-3*di - dk + ijk] - g13[3*di - dk + ijk] - 
       g13[-3*di + dk + ijk] + g13[3*di + dk + ijk] - 
       g13[-di + 3*dk + ijk] + g13[di + 3*dk + ijk]) + 
    9.*(-g13[-2*di - 3*dk + ijk] + g13[2*di - 3*dk + ijk] - 
       g13[-3*di - 2*dk + ijk] + g13[3*di - 2*dk + ijk] + 
       g13[-3*di + 2*dk + ijk] - g13[3*di + 2*dk + ijk] + 
       g13[-2*di + 3*dk + ijk] - g13[2*di + 3*dk + ijk]) - 
    g13[3*(-di + dk) + ijk] + g13[-3*(di + dk) + ijk] + 
    81.*(-g13[2*(di - dk) + ijk] - g13[2*(-di + dk) + ijk] + 
       g13[-2*(di + dk) + ijk] + g13[2*(di + dk) + ijk]) + 
    g13[3*(di + dk) + ijk])
;

deldelg1322
=
0.0011111111111111111111*oo2dx*oo2dz*
  (-g22[3*(di - dk) + ijk] - 405.*
     (g22[-di - 2*dk + ijk] + g22[-2*di - dk + ijk]) + 
    2025.*(g22[-di - dk + ijk] - g22[di - dk + ijk] - g22[-di + dk + ijk] + 
       g22[di + dk + ijk]) + 405.*
     (g22[di - 2*dk + ijk] + g22[2*di - dk + ijk] + g22[-2*di + dk + ijk] - 
       g22[2*di + dk + ijk] + g22[-di + 2*dk + ijk] - g22[di + 2*dk + ijk]) \
+ 45.*(g22[-di - 3*dk + ijk] - g22[di - 3*dk + ijk] + 
       g22[-3*di - dk + ijk] - g22[3*di - dk + ijk] - 
       g22[-3*di + dk + ijk] + g22[3*di + dk + ijk] - 
       g22[-di + 3*dk + ijk] + g22[di + 3*dk + ijk]) + 
    9.*(-g22[-2*di - 3*dk + ijk] + g22[2*di - 3*dk + ijk] - 
       g22[-3*di - 2*dk + ijk] + g22[3*di - 2*dk + ijk] + 
       g22[-3*di + 2*dk + ijk] - g22[3*di + 2*dk + ijk] + 
       g22[-2*di + 3*dk + ijk] - g22[2*di + 3*dk + ijk]) - 
    g22[3*(-di + dk) + ijk] + g22[-3*(di + dk) + ijk] + 
    81.*(-g22[2*(di - dk) + ijk] - g22[2*(-di + dk) + ijk] + 
       g22[-2*(di + dk) + ijk] + g22[2*(di + dk) + ijk]) + 
    g22[3*(di + dk) + ijk])
;

deldelg1323
=
0.0011111111111111111111*oo2dx*oo2dz*
  (-g23[3*(di - dk) + ijk] - 405.*
     (g23[-di - 2*dk + ijk] + g23[-2*di - dk + ijk]) + 
    2025.*(g23[-di - dk + ijk] - g23[di - dk + ijk] - g23[-di + dk + ijk] + 
       g23[di + dk + ijk]) + 405.*
     (g23[di - 2*dk + ijk] + g23[2*di - dk + ijk] + g23[-2*di + dk + ijk] - 
       g23[2*di + dk + ijk] + g23[-di + 2*dk + ijk] - g23[di + 2*dk + ijk]) \
+ 45.*(g23[-di - 3*dk + ijk] - g23[di - 3*dk + ijk] + 
       g23[-3*di - dk + ijk] - g23[3*di - dk + ijk] - 
       g23[-3*di + dk + ijk] + g23[3*di + dk + ijk] - 
       g23[-di + 3*dk + ijk] + g23[di + 3*dk + ijk]) + 
    9.*(-g23[-2*di - 3*dk + ijk] + g23[2*di - 3*dk + ijk] - 
       g23[-3*di - 2*dk + ijk] + g23[3*di - 2*dk + ijk] + 
       g23[-3*di + 2*dk + ijk] - g23[3*di + 2*dk + ijk] + 
       g23[-2*di + 3*dk + ijk] - g23[2*di + 3*dk + ijk]) - 
    g23[3*(-di + dk) + ijk] + g23[-3*(di + dk) + ijk] + 
    81.*(-g23[2*(di - dk) + ijk] - g23[2*(-di + dk) + ijk] + 
       g23[-2*(di + dk) + ijk] + g23[2*(di + dk) + ijk]) + 
    g23[3*(di + dk) + ijk])
;

deldelg1333
=
0.0011111111111111111111*oo2dx*oo2dz*
  (-g33[3*(di - dk) + ijk] - 405.*
     (g33[-di - 2*dk + ijk] + g33[-2*di - dk + ijk]) + 
    2025.*(g33[-di - dk + ijk] - g33[di - dk + ijk] - g33[-di + dk + ijk] + 
       g33[di + dk + ijk]) + 405.*
     (g33[di - 2*dk + ijk] + g33[2*di - dk + ijk] + g33[-2*di + dk + ijk] - 
       g33[2*di + dk + ijk] + g33[-di + 2*dk + ijk] - g33[di + 2*dk + ijk]) \
+ 45.*(g33[-di - 3*dk + ijk] - g33[di - 3*dk + ijk] + 
       g33[-3*di - dk + ijk] - g33[3*di - dk + ijk] - 
       g33[-3*di + dk + ijk] + g33[3*di + dk + ijk] - 
       g33[-di + 3*dk + ijk] + g33[di + 3*dk + ijk]) + 
    9.*(-g33[-2*di - 3*dk + ijk] + g33[2*di - 3*dk + ijk] - 
       g33[-3*di - 2*dk + ijk] + g33[3*di - 2*dk + ijk] + 
       g33[-3*di + 2*dk + ijk] - g33[3*di + 2*dk + ijk] + 
       g33[-2*di + 3*dk + ijk] - g33[2*di + 3*dk + ijk]) - 
    g33[3*(-di + dk) + ijk] + g33[-3*(di + dk) + ijk] + 
    81.*(-g33[2*(di - dk) + ijk] - g33[2*(-di + dk) + ijk] + 
       g33[-2*(di + dk) + ijk] + g33[2*(di + dk) + ijk]) + 
    g33[3*(di + dk) + ijk])
;

deldelg2211
=
0.0055555555555555555556*oody2*
  (-490.*g11[ijk] + 270.*(g11[-dj + ijk] + g11[dj + ijk]) - 
    27.*(g11[-2*dj + ijk] + g11[2*dj + ijk]) + 
    2.*(g11[-3*dj + ijk] + g11[3*dj + ijk]))
;

deldelg2212
=
0.0055555555555555555556*oody2*
  (-490.*g12[ijk] + 270.*(g12[-dj + ijk] + g12[dj + ijk]) - 
    27.*(g12[-2*dj + ijk] + g12[2*dj + ijk]) + 
    2.*(g12[-3*dj + ijk] + g12[3*dj + ijk]))
;

deldelg2213
=
0.0055555555555555555556*oody2*
  (-490.*g13[ijk] + 270.*(g13[-dj + ijk] + g13[dj + ijk]) - 
    27.*(g13[-2*dj + ijk] + g13[2*dj + ijk]) + 
    2.*(g13[-3*dj + ijk] + g13[3*dj + ijk]))
;

deldelg2222
=
0.0055555555555555555556*oody2*
  (-490.*g22[ijk] + 270.*(g22[-dj + ijk] + g22[dj + ijk]) - 
    27.*(g22[-2*dj + ijk] + g22[2*dj + ijk]) + 
    2.*(g22[-3*dj + ijk] + g22[3*dj + ijk]))
;

deldelg2223
=
0.0055555555555555555556*oody2*
  (-490.*g23[ijk] + 270.*(g23[-dj + ijk] + g23[dj + ijk]) - 
    27.*(g23[-2*dj + ijk] + g23[2*dj + ijk]) + 
    2.*(g23[-3*dj + ijk] + g23[3*dj + ijk]))
;

deldelg2233
=
0.0055555555555555555556*oody2*
  (-490.*g33[ijk] + 270.*(g33[-dj + ijk] + g33[dj + ijk]) - 
    27.*(g33[-2*dj + ijk] + g33[2*dj + ijk]) + 
    2.*(g33[-3*dj + ijk] + g33[3*dj + ijk]))
;

deldelg2311
=
0.0011111111111111111111*oo2dy*oo2dz*
  (-g11[3*(dj - dk) + ijk] - 405.*
     (g11[-dj - 2*dk + ijk] + g11[-2*dj - dk + ijk]) + 
    2025.*(g11[-dj - dk + ijk] - g11[dj - dk + ijk] - g11[-dj + dk + ijk] + 
       g11[dj + dk + ijk]) + 405.*
     (g11[dj - 2*dk + ijk] + g11[2*dj - dk + ijk] + g11[-2*dj + dk + ijk] - 
       g11[2*dj + dk + ijk] + g11[-dj + 2*dk + ijk] - g11[dj + 2*dk + ijk]) \
+ 45.*(g11[-dj - 3*dk + ijk] - g11[dj - 3*dk + ijk] + 
       g11[-3*dj - dk + ijk] - g11[3*dj - dk + ijk] - 
       g11[-3*dj + dk + ijk] + g11[3*dj + dk + ijk] - 
       g11[-dj + 3*dk + ijk] + g11[dj + 3*dk + ijk]) + 
    9.*(-g11[-2*dj - 3*dk + ijk] + g11[2*dj - 3*dk + ijk] - 
       g11[-3*dj - 2*dk + ijk] + g11[3*dj - 2*dk + ijk] + 
       g11[-3*dj + 2*dk + ijk] - g11[3*dj + 2*dk + ijk] + 
       g11[-2*dj + 3*dk + ijk] - g11[2*dj + 3*dk + ijk]) - 
    g11[3*(-dj + dk) + ijk] + g11[-3*(dj + dk) + ijk] + 
    81.*(-g11[2*(dj - dk) + ijk] - g11[2*(-dj + dk) + ijk] + 
       g11[-2*(dj + dk) + ijk] + g11[2*(dj + dk) + ijk]) + 
    g11[3*(dj + dk) + ijk])
;

deldelg2312
=
0.0011111111111111111111*oo2dy*oo2dz*
  (-g12[3*(dj - dk) + ijk] - 405.*
     (g12[-dj - 2*dk + ijk] + g12[-2*dj - dk + ijk]) + 
    2025.*(g12[-dj - dk + ijk] - g12[dj - dk + ijk] - g12[-dj + dk + ijk] + 
       g12[dj + dk + ijk]) + 405.*
     (g12[dj - 2*dk + ijk] + g12[2*dj - dk + ijk] + g12[-2*dj + dk + ijk] - 
       g12[2*dj + dk + ijk] + g12[-dj + 2*dk + ijk] - g12[dj + 2*dk + ijk]) \
+ 45.*(g12[-dj - 3*dk + ijk] - g12[dj - 3*dk + ijk] + 
       g12[-3*dj - dk + ijk] - g12[3*dj - dk + ijk] - 
       g12[-3*dj + dk + ijk] + g12[3*dj + dk + ijk] - 
       g12[-dj + 3*dk + ijk] + g12[dj + 3*dk + ijk]) + 
    9.*(-g12[-2*dj - 3*dk + ijk] + g12[2*dj - 3*dk + ijk] - 
       g12[-3*dj - 2*dk + ijk] + g12[3*dj - 2*dk + ijk] + 
       g12[-3*dj + 2*dk + ijk] - g12[3*dj + 2*dk + ijk] + 
       g12[-2*dj + 3*dk + ijk] - g12[2*dj + 3*dk + ijk]) - 
    g12[3*(-dj + dk) + ijk] + g12[-3*(dj + dk) + ijk] + 
    81.*(-g12[2*(dj - dk) + ijk] - g12[2*(-dj + dk) + ijk] + 
       g12[-2*(dj + dk) + ijk] + g12[2*(dj + dk) + ijk]) + 
    g12[3*(dj + dk) + ijk])
;

deldelg2313
=
0.0011111111111111111111*oo2dy*oo2dz*
  (-g13[3*(dj - dk) + ijk] - 405.*
     (g13[-dj - 2*dk + ijk] + g13[-2*dj - dk + ijk]) + 
    2025.*(g13[-dj - dk + ijk] - g13[dj - dk + ijk] - g13[-dj + dk + ijk] + 
       g13[dj + dk + ijk]) + 405.*
     (g13[dj - 2*dk + ijk] + g13[2*dj - dk + ijk] + g13[-2*dj + dk + ijk] - 
       g13[2*dj + dk + ijk] + g13[-dj + 2*dk + ijk] - g13[dj + 2*dk + ijk]) \
+ 45.*(g13[-dj - 3*dk + ijk] - g13[dj - 3*dk + ijk] + 
       g13[-3*dj - dk + ijk] - g13[3*dj - dk + ijk] - 
       g13[-3*dj + dk + ijk] + g13[3*dj + dk + ijk] - 
       g13[-dj + 3*dk + ijk] + g13[dj + 3*dk + ijk]) + 
    9.*(-g13[-2*dj - 3*dk + ijk] + g13[2*dj - 3*dk + ijk] - 
       g13[-3*dj - 2*dk + ijk] + g13[3*dj - 2*dk + ijk] + 
       g13[-3*dj + 2*dk + ijk] - g13[3*dj + 2*dk + ijk] + 
       g13[-2*dj + 3*dk + ijk] - g13[2*dj + 3*dk + ijk]) - 
    g13[3*(-dj + dk) + ijk] + g13[-3*(dj + dk) + ijk] + 
    81.*(-g13[2*(dj - dk) + ijk] - g13[2*(-dj + dk) + ijk] + 
       g13[-2*(dj + dk) + ijk] + g13[2*(dj + dk) + ijk]) + 
    g13[3*(dj + dk) + ijk])
;

deldelg2322
=
0.0011111111111111111111*oo2dy*oo2dz*
  (-g22[3*(dj - dk) + ijk] - 405.*
     (g22[-dj - 2*dk + ijk] + g22[-2*dj - dk + ijk]) + 
    2025.*(g22[-dj - dk + ijk] - g22[dj - dk + ijk] - g22[-dj + dk + ijk] + 
       g22[dj + dk + ijk]) + 405.*
     (g22[dj - 2*dk + ijk] + g22[2*dj - dk + ijk] + g22[-2*dj + dk + ijk] - 
       g22[2*dj + dk + ijk] + g22[-dj + 2*dk + ijk] - g22[dj + 2*dk + ijk]) \
+ 45.*(g22[-dj - 3*dk + ijk] - g22[dj - 3*dk + ijk] + 
       g22[-3*dj - dk + ijk] - g22[3*dj - dk + ijk] - 
       g22[-3*dj + dk + ijk] + g22[3*dj + dk + ijk] - 
       g22[-dj + 3*dk + ijk] + g22[dj + 3*dk + ijk]) + 
    9.*(-g22[-2*dj - 3*dk + ijk] + g22[2*dj - 3*dk + ijk] - 
       g22[-3*dj - 2*dk + ijk] + g22[3*dj - 2*dk + ijk] + 
       g22[-3*dj + 2*dk + ijk] - g22[3*dj + 2*dk + ijk] + 
       g22[-2*dj + 3*dk + ijk] - g22[2*dj + 3*dk + ijk]) - 
    g22[3*(-dj + dk) + ijk] + g22[-3*(dj + dk) + ijk] + 
    81.*(-g22[2*(dj - dk) + ijk] - g22[2*(-dj + dk) + ijk] + 
       g22[-2*(dj + dk) + ijk] + g22[2*(dj + dk) + ijk]) + 
    g22[3*(dj + dk) + ijk])
;

deldelg2323
=
0.0011111111111111111111*oo2dy*oo2dz*
  (-g23[3*(dj - dk) + ijk] - 405.*
     (g23[-dj - 2*dk + ijk] + g23[-2*dj - dk + ijk]) + 
    2025.*(g23[-dj - dk + ijk] - g23[dj - dk + ijk] - g23[-dj + dk + ijk] + 
       g23[dj + dk + ijk]) + 405.*
     (g23[dj - 2*dk + ijk] + g23[2*dj - dk + ijk] + g23[-2*dj + dk + ijk] - 
       g23[2*dj + dk + ijk] + g23[-dj + 2*dk + ijk] - g23[dj + 2*dk + ijk]) \
+ 45.*(g23[-dj - 3*dk + ijk] - g23[dj - 3*dk + ijk] + 
       g23[-3*dj - dk + ijk] - g23[3*dj - dk + ijk] - 
       g23[-3*dj + dk + ijk] + g23[3*dj + dk + ijk] - 
       g23[-dj + 3*dk + ijk] + g23[dj + 3*dk + ijk]) + 
    9.*(-g23[-2*dj - 3*dk + ijk] + g23[2*dj - 3*dk + ijk] - 
       g23[-3*dj - 2*dk + ijk] + g23[3*dj - 2*dk + ijk] + 
       g23[-3*dj + 2*dk + ijk] - g23[3*dj + 2*dk + ijk] + 
       g23[-2*dj + 3*dk + ijk] - g23[2*dj + 3*dk + ijk]) - 
    g23[3*(-dj + dk) + ijk] + g23[-3*(dj + dk) + ijk] + 
    81.*(-g23[2*(dj - dk) + ijk] - g23[2*(-dj + dk) + ijk] + 
       g23[-2*(dj + dk) + ijk] + g23[2*(dj + dk) + ijk]) + 
    g23[3*(dj + dk) + ijk])
;

deldelg2333
=
0.0011111111111111111111*oo2dy*oo2dz*
  (-g33[3*(dj - dk) + ijk] - 405.*
     (g33[-dj - 2*dk + ijk] + g33[-2*dj - dk + ijk]) + 
    2025.*(g33[-dj - dk + ijk] - g33[dj - dk + ijk] - g33[-dj + dk + ijk] + 
       g33[dj + dk + ijk]) + 405.*
     (g33[dj - 2*dk + ijk] + g33[2*dj - dk + ijk] + g33[-2*dj + dk + ijk] - 
       g33[2*dj + dk + ijk] + g33[-dj + 2*dk + ijk] - g33[dj + 2*dk + ijk]) \
+ 45.*(g33[-dj - 3*dk + ijk] - g33[dj - 3*dk + ijk] + 
       g33[-3*dj - dk + ijk] - g33[3*dj - dk + ijk] - 
       g33[-3*dj + dk + ijk] + g33[3*dj + dk + ijk] - 
       g33[-dj + 3*dk + ijk] + g33[dj + 3*dk + ijk]) + 
    9.*(-g33[-2*dj - 3*dk + ijk] + g33[2*dj - 3*dk + ijk] - 
       g33[-3*dj - 2*dk + ijk] + g33[3*dj - 2*dk + ijk] + 
       g33[-3*dj + 2*dk + ijk] - g33[3*dj + 2*dk + ijk] + 
       g33[-2*dj + 3*dk + ijk] - g33[2*dj + 3*dk + ijk]) - 
    g33[3*(-dj + dk) + ijk] + g33[-3*(dj + dk) + ijk] + 
    81.*(-g33[2*(dj - dk) + ijk] - g33[2*(-dj + dk) + ijk] + 
       g33[-2*(dj + dk) + ijk] + g33[2*(dj + dk) + ijk]) + 
    g33[3*(dj + dk) + ijk])
;

deldelg3311
=
0.0055555555555555555556*oodz2*
  (-490.*g11[ijk] + 270.*(g11[-dk + ijk] + g11[dk + ijk]) - 
    27.*(g11[-2*dk + ijk] + g11[2*dk + ijk]) + 
    2.*(g11[-3*dk + ijk] + g11[3*dk + ijk]))
;

deldelg3312
=
0.0055555555555555555556*oodz2*
  (-490.*g12[ijk] + 270.*(g12[-dk + ijk] + g12[dk + ijk]) - 
    27.*(g12[-2*dk + ijk] + g12[2*dk + ijk]) + 
    2.*(g12[-3*dk + ijk] + g12[3*dk + ijk]))
;

deldelg3313
=
0.0055555555555555555556*oodz2*
  (-490.*g13[ijk] + 270.*(g13[-dk + ijk] + g13[dk + ijk]) - 
    27.*(g13[-2*dk + ijk] + g13[2*dk + ijk]) + 
    2.*(g13[-3*dk + ijk] + g13[3*dk + ijk]))
;

deldelg3322
=
0.0055555555555555555556*oodz2*
  (-490.*g22[ijk] + 270.*(g22[-dk + ijk] + g22[dk + ijk]) - 
    27.*(g22[-2*dk + ijk] + g22[2*dk + ijk]) + 
    2.*(g22[-3*dk + ijk] + g22[3*dk + ijk]))
;

deldelg3323
=
0.0055555555555555555556*oodz2*
  (-490.*g23[ijk] + 270.*(g23[-dk + ijk] + g23[dk + ijk]) - 
    27.*(g23[-2*dk + ijk] + g23[2*dk + ijk]) + 
    2.*(g23[-3*dk + ijk] + g23[3*dk + ijk]))
;

deldelg3333
=
0.0055555555555555555556*oodz2*
  (-490.*g33[ijk] + 270.*(g33[-dk + ijk] + g33[dk + ijk]) - 
    27.*(g33[-2*dk + ijk] + g33[2*dk + ijk]) + 
    2.*(g33[-3*dk + ijk] + g33[3*dk + ijk]))
;

delG11
=
0.033333333333333333333*oo2dx*(-G1[-3*di + ijk] + 
    45.*(-G1[-di + ijk] + G1[di + ijk]) + 
    9.*(G1[-2*di + ijk] - G1[2*di + ijk]) + G1[3*di + ijk])
;

delG12
=
0.033333333333333333333*oo2dx*(-G2[-3*di + ijk] + 
    45.*(-G2[-di + ijk] + G2[di + ijk]) + 
    9.*(G2[-2*di + ijk] - G2[2*di + ijk]) + G2[3*di + ijk])
;

delG13
=
0.033333333333333333333*oo2dx*(-G3[-3*di + ijk] + 
    45.*(-G3[-di + ijk] + G3[di + ijk]) + 
    9.*(G3[-2*di + ijk] - G3[2*di + ijk]) + G3[3*di + ijk])
;

delG21
=
0.033333333333333333333*oo2dy*(-G1[-3*dj + ijk] + 
    45.*(-G1[-dj + ijk] + G1[dj + ijk]) + 
    9.*(G1[-2*dj + ijk] - G1[2*dj + ijk]) + G1[3*dj + ijk])
;

delG22
=
0.033333333333333333333*oo2dy*(-G2[-3*dj + ijk] + 
    45.*(-G2[-dj + ijk] + G2[dj + ijk]) + 
    9.*(G2[-2*dj + ijk] - G2[2*dj + ijk]) + G2[3*dj + ijk])
;

delG23
=
0.033333333333333333333*oo2dy*(-G3[-3*dj + ijk] + 
    45.*(-G3[-dj + ijk] + G3[dj + ijk]) + 
    9.*(G3[-2*dj + ijk] - G3[2*dj + ijk]) + G3[3*dj + ijk])
;

delG31
=
0.033333333333333333333*oo2dz*(-G1[-3*dk + ijk] + 
    45.*(-G1[-dk + ijk] + G1[dk + ijk]) + 
    9.*(G1[-2*dk + ijk] - G1[2*dk + ijk]) + G1[3*dk + ijk])
;

delG32
=
0.033333333333333333333*oo2dz*(-G2[-3*dk + ijk] + 
    45.*(-G2[-dk + ijk] + G2[dk + ijk]) + 
    9.*(G2[-2*dk + ijk] - G2[2*dk + ijk]) + G2[3*dk + ijk])
;

delG33
=
0.033333333333333333333*oo2dz*(-G3[-3*dk + ijk] + 
    45.*(-G3[-dk + ijk] + G3[dk + ijk]) + 
    9.*(G3[-2*dk + ijk] - G3[2*dk + ijk]) + G3[3*dk + ijk])
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

dchi1
=
0.033333333333333333333*oo2dx*(-chi[-3*di + ijk] + 
    45.*(-chi[-di + ijk] + chi[di + ijk]) + 
    9.*(chi[-2*di + ijk] - chi[2*di + ijk]) + chi[3*di + ijk])
;

dchi2
=
0.033333333333333333333*oo2dy*(-chi[-3*dj + ijk] + 
    45.*(-chi[-dj + ijk] + chi[dj + ijk]) + 
    9.*(chi[-2*dj + ijk] - chi[2*dj + ijk]) + chi[3*dj + ijk])
;

dchi3
=
0.033333333333333333333*oo2dz*(-chi[-3*dk + ijk] + 
    45.*(-chi[-dk + ijk] + chi[dk + ijk]) + 
    9.*(chi[-2*dk + ijk] - chi[2*dk + ijk]) + chi[3*dk + ijk])
;

ddchi11
=
0.0055555555555555555556*oodx2*
  (-490.*chi[ijk] + 270.*(chi[-di + ijk] + chi[di + ijk]) - 
    27.*(chi[-2*di + ijk] + chi[2*di + ijk]) + 
    2.*(chi[-3*di + ijk] + chi[3*di + ijk]))
;

ddchi12
=
0.0011111111111111111111*oo2dx*oo2dy*
  (-chi[3*(di - dj) + ijk] - 405.*
     (chi[-di - 2*dj + ijk] + chi[-2*di - dj + ijk]) + 
    2025.*(chi[-di - dj + ijk] - chi[di - dj + ijk] - chi[-di + dj + ijk] + 
       chi[di + dj + ijk]) + 405.*
     (chi[di - 2*dj + ijk] + chi[2*di - dj + ijk] + chi[-2*di + dj + ijk] - 
       chi[2*di + dj + ijk] + chi[-di + 2*dj + ijk] - chi[di + 2*dj + ijk]) \
+ 45.*(chi[-di - 3*dj + ijk] - chi[di - 3*dj + ijk] + 
       chi[-3*di - dj + ijk] - chi[3*di - dj + ijk] - 
       chi[-3*di + dj + ijk] + chi[3*di + dj + ijk] - 
       chi[-di + 3*dj + ijk] + chi[di + 3*dj + ijk]) + 
    9.*(-chi[-2*di - 3*dj + ijk] + chi[2*di - 3*dj + ijk] - 
       chi[-3*di - 2*dj + ijk] + chi[3*di - 2*dj + ijk] + 
       chi[-3*di + 2*dj + ijk] - chi[3*di + 2*dj + ijk] + 
       chi[-2*di + 3*dj + ijk] - chi[2*di + 3*dj + ijk]) - 
    chi[3*(-di + dj) + ijk] + chi[-3*(di + dj) + ijk] + 
    81.*(-chi[2*(di - dj) + ijk] - chi[2*(-di + dj) + ijk] + 
       chi[-2*(di + dj) + ijk] + chi[2*(di + dj) + ijk]) + 
    chi[3*(di + dj) + ijk])
;

ddchi13
=
0.0011111111111111111111*oo2dx*oo2dz*
  (-chi[3*(di - dk) + ijk] - 405.*
     (chi[-di - 2*dk + ijk] + chi[-2*di - dk + ijk]) + 
    2025.*(chi[-di - dk + ijk] - chi[di - dk + ijk] - chi[-di + dk + ijk] + 
       chi[di + dk + ijk]) + 405.*
     (chi[di - 2*dk + ijk] + chi[2*di - dk + ijk] + chi[-2*di + dk + ijk] - 
       chi[2*di + dk + ijk] + chi[-di + 2*dk + ijk] - chi[di + 2*dk + ijk]) \
+ 45.*(chi[-di - 3*dk + ijk] - chi[di - 3*dk + ijk] + 
       chi[-3*di - dk + ijk] - chi[3*di - dk + ijk] - 
       chi[-3*di + dk + ijk] + chi[3*di + dk + ijk] - 
       chi[-di + 3*dk + ijk] + chi[di + 3*dk + ijk]) + 
    9.*(-chi[-2*di - 3*dk + ijk] + chi[2*di - 3*dk + ijk] - 
       chi[-3*di - 2*dk + ijk] + chi[3*di - 2*dk + ijk] + 
       chi[-3*di + 2*dk + ijk] - chi[3*di + 2*dk + ijk] + 
       chi[-2*di + 3*dk + ijk] - chi[2*di + 3*dk + ijk]) - 
    chi[3*(-di + dk) + ijk] + chi[-3*(di + dk) + ijk] + 
    81.*(-chi[2*(di - dk) + ijk] - chi[2*(-di + dk) + ijk] + 
       chi[-2*(di + dk) + ijk] + chi[2*(di + dk) + ijk]) + 
    chi[3*(di + dk) + ijk])
;

ddchi22
=
0.0055555555555555555556*oody2*
  (-490.*chi[ijk] + 270.*(chi[-dj + ijk] + chi[dj + ijk]) - 
    27.*(chi[-2*dj + ijk] + chi[2*dj + ijk]) + 
    2.*(chi[-3*dj + ijk] + chi[3*dj + ijk]))
;

ddchi23
=
0.0011111111111111111111*oo2dy*oo2dz*
  (-chi[3*(dj - dk) + ijk] - 405.*
     (chi[-dj - 2*dk + ijk] + chi[-2*dj - dk + ijk]) + 
    2025.*(chi[-dj - dk + ijk] - chi[dj - dk + ijk] - chi[-dj + dk + ijk] + 
       chi[dj + dk + ijk]) + 405.*
     (chi[dj - 2*dk + ijk] + chi[2*dj - dk + ijk] + chi[-2*dj + dk + ijk] - 
       chi[2*dj + dk + ijk] + chi[-dj + 2*dk + ijk] - chi[dj + 2*dk + ijk]) \
+ 45.*(chi[-dj - 3*dk + ijk] - chi[dj - 3*dk + ijk] + 
       chi[-3*dj - dk + ijk] - chi[3*dj - dk + ijk] - 
       chi[-3*dj + dk + ijk] + chi[3*dj + dk + ijk] - 
       chi[-dj + 3*dk + ijk] + chi[dj + 3*dk + ijk]) + 
    9.*(-chi[-2*dj - 3*dk + ijk] + chi[2*dj - 3*dk + ijk] - 
       chi[-3*dj - 2*dk + ijk] + chi[3*dj - 2*dk + ijk] + 
       chi[-3*dj + 2*dk + ijk] - chi[3*dj + 2*dk + ijk] + 
       chi[-2*dj + 3*dk + ijk] - chi[2*dj + 3*dk + ijk]) - 
    chi[3*(-dj + dk) + ijk] + chi[-3*(dj + dk) + ijk] + 
    81.*(-chi[2*(dj - dk) + ijk] - chi[2*(-dj + dk) + ijk] + 
       chi[-2*(dj + dk) + ijk] + chi[2*(dj + dk) + ijk]) + 
    chi[3*(dj + dk) + ijk])
;

ddchi33
=
0.0055555555555555555556*oodz2*
  (-490.*chi[ijk] + 270.*(chi[-dk + ijk] + chi[dk + ijk]) - 
    27.*(chi[-2*dk + ijk] + chi[2*dk + ijk]) + 
    2.*(chi[-3*dk + ijk] + chi[3*dk + ijk]))
;


#endif 


} else if (order_centered == 8 || boundaryNaway(4)) { 


#ifdef REDUCEORDERTO6                           
  errorexit("CompilerFlag reduced order to 6");
 #else
da1
=
0.0023809523809523809524*oo2dx*
  (3.*alpha[-4*di + ijk] + 168.*alpha[-2*di + ijk] + 
    672.*(-alpha[-di + ijk] + alpha[di + ijk]) - 168.*alpha[2*di + ijk] + 
    32.*(-alpha[-3*di + ijk] + alpha[3*di + ijk]) - 3.*alpha[4*di + ijk])
;

da2
=
0.0023809523809523809524*oo2dy*
  (3.*alpha[-4*dj + ijk] + 168.*alpha[-2*dj + ijk] + 
    672.*(-alpha[-dj + ijk] + alpha[dj + ijk]) - 168.*alpha[2*dj + ijk] + 
    32.*(-alpha[-3*dj + ijk] + alpha[3*dj + ijk]) - 3.*alpha[4*dj + ijk])
;

da3
=
0.0023809523809523809524*oo2dz*
  (3.*alpha[-4*dk + ijk] + 168.*alpha[-2*dk + ijk] + 
    672.*(-alpha[-dk + ijk] + alpha[dk + ijk]) - 168.*alpha[2*dk + ijk] + 
    32.*(-alpha[-3*dk + ijk] + alpha[3*dk + ijk]) - 3.*alpha[4*dk + ijk])
;

dda11
=
0.0001984126984126984127*oodx2*
  (-14350.*alpha[ijk] + 8064.*(alpha[-di + ijk] + alpha[di + ijk]) - 
    1008.*(alpha[-2*di + ijk] + alpha[2*di + ijk]) + 
    128.*(alpha[-3*di + ijk] + alpha[3*di + ijk]) - 
    9.*(alpha[-4*di + ijk] + alpha[4*di + ijk]))
;

dda12
=
5.6689342403628117914e-6*oo2dx*oo2dy*
  (451584.*(alpha[-di - dj + ijk] - alpha[di - dj + ijk] - 
       alpha[-di + dj + ijk] + alpha[di + dj + ijk]) + 
    112896.*(-alpha[-di - 2*dj + ijk] + alpha[di - 2*dj + ijk] - 
       alpha[-2*di - dj + ijk] + alpha[2*di - dj + ijk] + 
       alpha[-2*di + dj + ijk] - alpha[2*di + dj + ijk] + 
       alpha[-di + 2*dj + ijk] - alpha[di + 2*dj + ijk]) - 
    21504.*(alpha[di - 3*dj + ijk] + alpha[3*di - dj + ijk] + 
       alpha[-3*di + dj + ijk] + alpha[-di + 3*dj + ijk]) + 
    21504.*(alpha[-di - 3*dj + ijk] + alpha[-3*di - dj + ijk] + 
       alpha[3*di + dj + ijk] + alpha[di + 3*dj + ijk]) + 
    5376.*(-alpha[-2*di - 3*dj + ijk] + alpha[2*di - 3*dj + ijk] - 
       alpha[-3*di - 2*dj + ijk] + alpha[3*di - 2*dj + ijk] + 
       alpha[-3*di + 2*dj + ijk] - alpha[3*di + 2*dj + ijk] + 
       alpha[-2*di + 3*dj + ijk] - alpha[2*di + 3*dj + ijk]) - 
    504.*(alpha[2*di - 4*dj + ijk] + alpha[4*di - 2*dj + ijk] + 
       alpha[-4*di + 2*dj + ijk] + alpha[-2*di + 4*dj + ijk]) + 
    2016.*(-alpha[-di - 4*dj + ijk] + alpha[di - 4*dj + ijk] - 
       alpha[-4*di - dj + ijk] + alpha[4*di - dj + ijk] + 
       alpha[-4*di + dj + ijk] - alpha[4*di + dj + ijk] + 
       alpha[-di + 4*dj + ijk] - alpha[di + 4*dj + ijk]) + 
    504.*(alpha[-2*di - 4*dj + ijk] + alpha[-4*di - 2*dj + ijk] + 
       alpha[4*di + 2*dj + ijk] + alpha[2*di + 4*dj + ijk]) + 
    96.*(-alpha[-3*di - 4*dj + ijk] + alpha[3*di - 4*dj + ijk] - 
       alpha[-4*di - 3*dj + ijk] + alpha[4*di - 3*dj + ijk] + 
       alpha[-4*di + 3*dj + ijk] - alpha[4*di + 3*dj + ijk] + 
       alpha[-3*di + 4*dj + ijk] - alpha[3*di + 4*dj + ijk]) - 
    28224.*(alpha[2*(di - dj) + ijk] + alpha[2*(-di + dj) + ijk]) - 
    1024.*(alpha[3*(di - dj) + ijk] + alpha[3*(-di + dj) + ijk]) - 
    9.*(alpha[4*(di - dj) + ijk] + alpha[4*(-di + dj) + ijk]) + 
    28224.*(alpha[-2*(di + dj) + ijk] + alpha[2*(di + dj) + ijk]) + 
    1024.*(alpha[-3*(di + dj) + ijk] + alpha[3*(di + dj) + ijk]) + 
    9.*(alpha[-4*(di + dj) + ijk] + alpha[4*(di + dj) + ijk]))
;

dda13
=
5.6689342403628117914e-6*oo2dx*oo2dz*
  (451584.*(alpha[-di - dk + ijk] - alpha[di - dk + ijk] - 
       alpha[-di + dk + ijk] + alpha[di + dk + ijk]) + 
    112896.*(-alpha[-di - 2*dk + ijk] + alpha[di - 2*dk + ijk] - 
       alpha[-2*di - dk + ijk] + alpha[2*di - dk + ijk] + 
       alpha[-2*di + dk + ijk] - alpha[2*di + dk + ijk] + 
       alpha[-di + 2*dk + ijk] - alpha[di + 2*dk + ijk]) - 
    21504.*(alpha[di - 3*dk + ijk] + alpha[3*di - dk + ijk] + 
       alpha[-3*di + dk + ijk] + alpha[-di + 3*dk + ijk]) + 
    21504.*(alpha[-di - 3*dk + ijk] + alpha[-3*di - dk + ijk] + 
       alpha[3*di + dk + ijk] + alpha[di + 3*dk + ijk]) + 
    5376.*(-alpha[-2*di - 3*dk + ijk] + alpha[2*di - 3*dk + ijk] - 
       alpha[-3*di - 2*dk + ijk] + alpha[3*di - 2*dk + ijk] + 
       alpha[-3*di + 2*dk + ijk] - alpha[3*di + 2*dk + ijk] + 
       alpha[-2*di + 3*dk + ijk] - alpha[2*di + 3*dk + ijk]) - 
    504.*(alpha[2*di - 4*dk + ijk] + alpha[4*di - 2*dk + ijk] + 
       alpha[-4*di + 2*dk + ijk] + alpha[-2*di + 4*dk + ijk]) + 
    2016.*(-alpha[-di - 4*dk + ijk] + alpha[di - 4*dk + ijk] - 
       alpha[-4*di - dk + ijk] + alpha[4*di - dk + ijk] + 
       alpha[-4*di + dk + ijk] - alpha[4*di + dk + ijk] + 
       alpha[-di + 4*dk + ijk] - alpha[di + 4*dk + ijk]) + 
    504.*(alpha[-2*di - 4*dk + ijk] + alpha[-4*di - 2*dk + ijk] + 
       alpha[4*di + 2*dk + ijk] + alpha[2*di + 4*dk + ijk]) + 
    96.*(-alpha[-3*di - 4*dk + ijk] + alpha[3*di - 4*dk + ijk] - 
       alpha[-4*di - 3*dk + ijk] + alpha[4*di - 3*dk + ijk] + 
       alpha[-4*di + 3*dk + ijk] - alpha[4*di + 3*dk + ijk] + 
       alpha[-3*di + 4*dk + ijk] - alpha[3*di + 4*dk + ijk]) - 
    28224.*(alpha[2*(di - dk) + ijk] + alpha[2*(-di + dk) + ijk]) - 
    1024.*(alpha[3*(di - dk) + ijk] + alpha[3*(-di + dk) + ijk]) - 
    9.*(alpha[4*(di - dk) + ijk] + alpha[4*(-di + dk) + ijk]) + 
    28224.*(alpha[-2*(di + dk) + ijk] + alpha[2*(di + dk) + ijk]) + 
    1024.*(alpha[-3*(di + dk) + ijk] + alpha[3*(di + dk) + ijk]) + 
    9.*(alpha[-4*(di + dk) + ijk] + alpha[4*(di + dk) + ijk]))
;

dda22
=
0.0001984126984126984127*oody2*
  (-14350.*alpha[ijk] + 8064.*(alpha[-dj + ijk] + alpha[dj + ijk]) - 
    1008.*(alpha[-2*dj + ijk] + alpha[2*dj + ijk]) + 
    128.*(alpha[-3*dj + ijk] + alpha[3*dj + ijk]) - 
    9.*(alpha[-4*dj + ijk] + alpha[4*dj + ijk]))
;

dda23
=
5.6689342403628117914e-6*oo2dy*oo2dz*
  (451584.*(alpha[-dj - dk + ijk] - alpha[dj - dk + ijk] - 
       alpha[-dj + dk + ijk] + alpha[dj + dk + ijk]) + 
    112896.*(-alpha[-dj - 2*dk + ijk] + alpha[dj - 2*dk + ijk] - 
       alpha[-2*dj - dk + ijk] + alpha[2*dj - dk + ijk] + 
       alpha[-2*dj + dk + ijk] - alpha[2*dj + dk + ijk] + 
       alpha[-dj + 2*dk + ijk] - alpha[dj + 2*dk + ijk]) - 
    21504.*(alpha[dj - 3*dk + ijk] + alpha[3*dj - dk + ijk] + 
       alpha[-3*dj + dk + ijk] + alpha[-dj + 3*dk + ijk]) + 
    21504.*(alpha[-dj - 3*dk + ijk] + alpha[-3*dj - dk + ijk] + 
       alpha[3*dj + dk + ijk] + alpha[dj + 3*dk + ijk]) + 
    5376.*(-alpha[-2*dj - 3*dk + ijk] + alpha[2*dj - 3*dk + ijk] - 
       alpha[-3*dj - 2*dk + ijk] + alpha[3*dj - 2*dk + ijk] + 
       alpha[-3*dj + 2*dk + ijk] - alpha[3*dj + 2*dk + ijk] + 
       alpha[-2*dj + 3*dk + ijk] - alpha[2*dj + 3*dk + ijk]) - 
    504.*(alpha[2*dj - 4*dk + ijk] + alpha[4*dj - 2*dk + ijk] + 
       alpha[-4*dj + 2*dk + ijk] + alpha[-2*dj + 4*dk + ijk]) + 
    2016.*(-alpha[-dj - 4*dk + ijk] + alpha[dj - 4*dk + ijk] - 
       alpha[-4*dj - dk + ijk] + alpha[4*dj - dk + ijk] + 
       alpha[-4*dj + dk + ijk] - alpha[4*dj + dk + ijk] + 
       alpha[-dj + 4*dk + ijk] - alpha[dj + 4*dk + ijk]) + 
    504.*(alpha[-2*dj - 4*dk + ijk] + alpha[-4*dj - 2*dk + ijk] + 
       alpha[4*dj + 2*dk + ijk] + alpha[2*dj + 4*dk + ijk]) + 
    96.*(-alpha[-3*dj - 4*dk + ijk] + alpha[3*dj - 4*dk + ijk] - 
       alpha[-4*dj - 3*dk + ijk] + alpha[4*dj - 3*dk + ijk] + 
       alpha[-4*dj + 3*dk + ijk] - alpha[4*dj + 3*dk + ijk] + 
       alpha[-3*dj + 4*dk + ijk] - alpha[3*dj + 4*dk + ijk]) - 
    28224.*(alpha[2*(dj - dk) + ijk] + alpha[2*(-dj + dk) + ijk]) - 
    1024.*(alpha[3*(dj - dk) + ijk] + alpha[3*(-dj + dk) + ijk]) - 
    9.*(alpha[4*(dj - dk) + ijk] + alpha[4*(-dj + dk) + ijk]) + 
    28224.*(alpha[-2*(dj + dk) + ijk] + alpha[2*(dj + dk) + ijk]) + 
    1024.*(alpha[-3*(dj + dk) + ijk] + alpha[3*(dj + dk) + ijk]) + 
    9.*(alpha[-4*(dj + dk) + ijk] + alpha[4*(dj + dk) + ijk]))
;

dda33
=
0.0001984126984126984127*oodz2*
  (-14350.*alpha[ijk] + 8064.*(alpha[-dk + ijk] + alpha[dk + ijk]) - 
    1008.*(alpha[-2*dk + ijk] + alpha[2*dk + ijk]) + 
    128.*(alpha[-3*dk + ijk] + alpha[3*dk + ijk]) - 
    9.*(alpha[-4*dk + ijk] + alpha[4*dk + ijk]))
;

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

ddb111
=
0.0001984126984126984127*oodx2*
  (-14350.*beta1[ijk] + 8064.*(beta1[-di + ijk] + beta1[di + ijk]) - 
    1008.*(beta1[-2*di + ijk] + beta1[2*di + ijk]) + 
    128.*(beta1[-3*di + ijk] + beta1[3*di + ijk]) - 
    9.*(beta1[-4*di + ijk] + beta1[4*di + ijk]))
;

ddb112
=
0.0001984126984126984127*oodx2*
  (-14350.*beta2[ijk] + 8064.*(beta2[-di + ijk] + beta2[di + ijk]) - 
    1008.*(beta2[-2*di + ijk] + beta2[2*di + ijk]) + 
    128.*(beta2[-3*di + ijk] + beta2[3*di + ijk]) - 
    9.*(beta2[-4*di + ijk] + beta2[4*di + ijk]))
;

ddb113
=
0.0001984126984126984127*oodx2*
  (-14350.*beta3[ijk] + 8064.*(beta3[-di + ijk] + beta3[di + ijk]) - 
    1008.*(beta3[-2*di + ijk] + beta3[2*di + ijk]) + 
    128.*(beta3[-3*di + ijk] + beta3[3*di + ijk]) - 
    9.*(beta3[-4*di + ijk] + beta3[4*di + ijk]))
;

ddb121
=
5.6689342403628117914e-6*oo2dx*oo2dy*
  (451584.*(beta1[-di - dj + ijk] - beta1[di - dj + ijk] - 
       beta1[-di + dj + ijk] + beta1[di + dj + ijk]) + 
    112896.*(-beta1[-di - 2*dj + ijk] + beta1[di - 2*dj + ijk] - 
       beta1[-2*di - dj + ijk] + beta1[2*di - dj + ijk] + 
       beta1[-2*di + dj + ijk] - beta1[2*di + dj + ijk] + 
       beta1[-di + 2*dj + ijk] - beta1[di + 2*dj + ijk]) - 
    21504.*(beta1[di - 3*dj + ijk] + beta1[3*di - dj + ijk] + 
       beta1[-3*di + dj + ijk] + beta1[-di + 3*dj + ijk]) + 
    21504.*(beta1[-di - 3*dj + ijk] + beta1[-3*di - dj + ijk] + 
       beta1[3*di + dj + ijk] + beta1[di + 3*dj + ijk]) + 
    5376.*(-beta1[-2*di - 3*dj + ijk] + beta1[2*di - 3*dj + ijk] - 
       beta1[-3*di - 2*dj + ijk] + beta1[3*di - 2*dj + ijk] + 
       beta1[-3*di + 2*dj + ijk] - beta1[3*di + 2*dj + ijk] + 
       beta1[-2*di + 3*dj + ijk] - beta1[2*di + 3*dj + ijk]) - 
    504.*(beta1[2*di - 4*dj + ijk] + beta1[4*di - 2*dj + ijk] + 
       beta1[-4*di + 2*dj + ijk] + beta1[-2*di + 4*dj + ijk]) + 
    2016.*(-beta1[-di - 4*dj + ijk] + beta1[di - 4*dj + ijk] - 
       beta1[-4*di - dj + ijk] + beta1[4*di - dj + ijk] + 
       beta1[-4*di + dj + ijk] - beta1[4*di + dj + ijk] + 
       beta1[-di + 4*dj + ijk] - beta1[di + 4*dj + ijk]) + 
    504.*(beta1[-2*di - 4*dj + ijk] + beta1[-4*di - 2*dj + ijk] + 
       beta1[4*di + 2*dj + ijk] + beta1[2*di + 4*dj + ijk]) + 
    96.*(-beta1[-3*di - 4*dj + ijk] + beta1[3*di - 4*dj + ijk] - 
       beta1[-4*di - 3*dj + ijk] + beta1[4*di - 3*dj + ijk] + 
       beta1[-4*di + 3*dj + ijk] - beta1[4*di + 3*dj + ijk] + 
       beta1[-3*di + 4*dj + ijk] - beta1[3*di + 4*dj + ijk]) - 
    28224.*(beta1[2*(di - dj) + ijk] + beta1[2*(-di + dj) + ijk]) - 
    1024.*(beta1[3*(di - dj) + ijk] + beta1[3*(-di + dj) + ijk]) - 
    9.*(beta1[4*(di - dj) + ijk] + beta1[4*(-di + dj) + ijk]) + 
    28224.*(beta1[-2*(di + dj) + ijk] + beta1[2*(di + dj) + ijk]) + 
    1024.*(beta1[-3*(di + dj) + ijk] + beta1[3*(di + dj) + ijk]) + 
    9.*(beta1[-4*(di + dj) + ijk] + beta1[4*(di + dj) + ijk]))
;

ddb122
=
5.6689342403628117914e-6*oo2dx*oo2dy*
  (451584.*(beta2[-di - dj + ijk] - beta2[di - dj + ijk] - 
       beta2[-di + dj + ijk] + beta2[di + dj + ijk]) + 
    112896.*(-beta2[-di - 2*dj + ijk] + beta2[di - 2*dj + ijk] - 
       beta2[-2*di - dj + ijk] + beta2[2*di - dj + ijk] + 
       beta2[-2*di + dj + ijk] - beta2[2*di + dj + ijk] + 
       beta2[-di + 2*dj + ijk] - beta2[di + 2*dj + ijk]) - 
    21504.*(beta2[di - 3*dj + ijk] + beta2[3*di - dj + ijk] + 
       beta2[-3*di + dj + ijk] + beta2[-di + 3*dj + ijk]) + 
    21504.*(beta2[-di - 3*dj + ijk] + beta2[-3*di - dj + ijk] + 
       beta2[3*di + dj + ijk] + beta2[di + 3*dj + ijk]) + 
    5376.*(-beta2[-2*di - 3*dj + ijk] + beta2[2*di - 3*dj + ijk] - 
       beta2[-3*di - 2*dj + ijk] + beta2[3*di - 2*dj + ijk] + 
       beta2[-3*di + 2*dj + ijk] - beta2[3*di + 2*dj + ijk] + 
       beta2[-2*di + 3*dj + ijk] - beta2[2*di + 3*dj + ijk]) - 
    504.*(beta2[2*di - 4*dj + ijk] + beta2[4*di - 2*dj + ijk] + 
       beta2[-4*di + 2*dj + ijk] + beta2[-2*di + 4*dj + ijk]) + 
    2016.*(-beta2[-di - 4*dj + ijk] + beta2[di - 4*dj + ijk] - 
       beta2[-4*di - dj + ijk] + beta2[4*di - dj + ijk] + 
       beta2[-4*di + dj + ijk] - beta2[4*di + dj + ijk] + 
       beta2[-di + 4*dj + ijk] - beta2[di + 4*dj + ijk]) + 
    504.*(beta2[-2*di - 4*dj + ijk] + beta2[-4*di - 2*dj + ijk] + 
       beta2[4*di + 2*dj + ijk] + beta2[2*di + 4*dj + ijk]) + 
    96.*(-beta2[-3*di - 4*dj + ijk] + beta2[3*di - 4*dj + ijk] - 
       beta2[-4*di - 3*dj + ijk] + beta2[4*di - 3*dj + ijk] + 
       beta2[-4*di + 3*dj + ijk] - beta2[4*di + 3*dj + ijk] + 
       beta2[-3*di + 4*dj + ijk] - beta2[3*di + 4*dj + ijk]) - 
    28224.*(beta2[2*(di - dj) + ijk] + beta2[2*(-di + dj) + ijk]) - 
    1024.*(beta2[3*(di - dj) + ijk] + beta2[3*(-di + dj) + ijk]) - 
    9.*(beta2[4*(di - dj) + ijk] + beta2[4*(-di + dj) + ijk]) + 
    28224.*(beta2[-2*(di + dj) + ijk] + beta2[2*(di + dj) + ijk]) + 
    1024.*(beta2[-3*(di + dj) + ijk] + beta2[3*(di + dj) + ijk]) + 
    9.*(beta2[-4*(di + dj) + ijk] + beta2[4*(di + dj) + ijk]))
;

ddb123
=
5.6689342403628117914e-6*oo2dx*oo2dy*
  (451584.*(beta3[-di - dj + ijk] - beta3[di - dj + ijk] - 
       beta3[-di + dj + ijk] + beta3[di + dj + ijk]) + 
    112896.*(-beta3[-di - 2*dj + ijk] + beta3[di - 2*dj + ijk] - 
       beta3[-2*di - dj + ijk] + beta3[2*di - dj + ijk] + 
       beta3[-2*di + dj + ijk] - beta3[2*di + dj + ijk] + 
       beta3[-di + 2*dj + ijk] - beta3[di + 2*dj + ijk]) - 
    21504.*(beta3[di - 3*dj + ijk] + beta3[3*di - dj + ijk] + 
       beta3[-3*di + dj + ijk] + beta3[-di + 3*dj + ijk]) + 
    21504.*(beta3[-di - 3*dj + ijk] + beta3[-3*di - dj + ijk] + 
       beta3[3*di + dj + ijk] + beta3[di + 3*dj + ijk]) + 
    5376.*(-beta3[-2*di - 3*dj + ijk] + beta3[2*di - 3*dj + ijk] - 
       beta3[-3*di - 2*dj + ijk] + beta3[3*di - 2*dj + ijk] + 
       beta3[-3*di + 2*dj + ijk] - beta3[3*di + 2*dj + ijk] + 
       beta3[-2*di + 3*dj + ijk] - beta3[2*di + 3*dj + ijk]) - 
    504.*(beta3[2*di - 4*dj + ijk] + beta3[4*di - 2*dj + ijk] + 
       beta3[-4*di + 2*dj + ijk] + beta3[-2*di + 4*dj + ijk]) + 
    2016.*(-beta3[-di - 4*dj + ijk] + beta3[di - 4*dj + ijk] - 
       beta3[-4*di - dj + ijk] + beta3[4*di - dj + ijk] + 
       beta3[-4*di + dj + ijk] - beta3[4*di + dj + ijk] + 
       beta3[-di + 4*dj + ijk] - beta3[di + 4*dj + ijk]) + 
    504.*(beta3[-2*di - 4*dj + ijk] + beta3[-4*di - 2*dj + ijk] + 
       beta3[4*di + 2*dj + ijk] + beta3[2*di + 4*dj + ijk]) + 
    96.*(-beta3[-3*di - 4*dj + ijk] + beta3[3*di - 4*dj + ijk] - 
       beta3[-4*di - 3*dj + ijk] + beta3[4*di - 3*dj + ijk] + 
       beta3[-4*di + 3*dj + ijk] - beta3[4*di + 3*dj + ijk] + 
       beta3[-3*di + 4*dj + ijk] - beta3[3*di + 4*dj + ijk]) - 
    28224.*(beta3[2*(di - dj) + ijk] + beta3[2*(-di + dj) + ijk]) - 
    1024.*(beta3[3*(di - dj) + ijk] + beta3[3*(-di + dj) + ijk]) - 
    9.*(beta3[4*(di - dj) + ijk] + beta3[4*(-di + dj) + ijk]) + 
    28224.*(beta3[-2*(di + dj) + ijk] + beta3[2*(di + dj) + ijk]) + 
    1024.*(beta3[-3*(di + dj) + ijk] + beta3[3*(di + dj) + ijk]) + 
    9.*(beta3[-4*(di + dj) + ijk] + beta3[4*(di + dj) + ijk]))
;

ddb131
=
5.6689342403628117914e-6*oo2dx*oo2dz*
  (451584.*(beta1[-di - dk + ijk] - beta1[di - dk + ijk] - 
       beta1[-di + dk + ijk] + beta1[di + dk + ijk]) + 
    112896.*(-beta1[-di - 2*dk + ijk] + beta1[di - 2*dk + ijk] - 
       beta1[-2*di - dk + ijk] + beta1[2*di - dk + ijk] + 
       beta1[-2*di + dk + ijk] - beta1[2*di + dk + ijk] + 
       beta1[-di + 2*dk + ijk] - beta1[di + 2*dk + ijk]) - 
    21504.*(beta1[di - 3*dk + ijk] + beta1[3*di - dk + ijk] + 
       beta1[-3*di + dk + ijk] + beta1[-di + 3*dk + ijk]) + 
    21504.*(beta1[-di - 3*dk + ijk] + beta1[-3*di - dk + ijk] + 
       beta1[3*di + dk + ijk] + beta1[di + 3*dk + ijk]) + 
    5376.*(-beta1[-2*di - 3*dk + ijk] + beta1[2*di - 3*dk + ijk] - 
       beta1[-3*di - 2*dk + ijk] + beta1[3*di - 2*dk + ijk] + 
       beta1[-3*di + 2*dk + ijk] - beta1[3*di + 2*dk + ijk] + 
       beta1[-2*di + 3*dk + ijk] - beta1[2*di + 3*dk + ijk]) - 
    504.*(beta1[2*di - 4*dk + ijk] + beta1[4*di - 2*dk + ijk] + 
       beta1[-4*di + 2*dk + ijk] + beta1[-2*di + 4*dk + ijk]) + 
    2016.*(-beta1[-di - 4*dk + ijk] + beta1[di - 4*dk + ijk] - 
       beta1[-4*di - dk + ijk] + beta1[4*di - dk + ijk] + 
       beta1[-4*di + dk + ijk] - beta1[4*di + dk + ijk] + 
       beta1[-di + 4*dk + ijk] - beta1[di + 4*dk + ijk]) + 
    504.*(beta1[-2*di - 4*dk + ijk] + beta1[-4*di - 2*dk + ijk] + 
       beta1[4*di + 2*dk + ijk] + beta1[2*di + 4*dk + ijk]) + 
    96.*(-beta1[-3*di - 4*dk + ijk] + beta1[3*di - 4*dk + ijk] - 
       beta1[-4*di - 3*dk + ijk] + beta1[4*di - 3*dk + ijk] + 
       beta1[-4*di + 3*dk + ijk] - beta1[4*di + 3*dk + ijk] + 
       beta1[-3*di + 4*dk + ijk] - beta1[3*di + 4*dk + ijk]) - 
    28224.*(beta1[2*(di - dk) + ijk] + beta1[2*(-di + dk) + ijk]) - 
    1024.*(beta1[3*(di - dk) + ijk] + beta1[3*(-di + dk) + ijk]) - 
    9.*(beta1[4*(di - dk) + ijk] + beta1[4*(-di + dk) + ijk]) + 
    28224.*(beta1[-2*(di + dk) + ijk] + beta1[2*(di + dk) + ijk]) + 
    1024.*(beta1[-3*(di + dk) + ijk] + beta1[3*(di + dk) + ijk]) + 
    9.*(beta1[-4*(di + dk) + ijk] + beta1[4*(di + dk) + ijk]))
;

ddb132
=
5.6689342403628117914e-6*oo2dx*oo2dz*
  (451584.*(beta2[-di - dk + ijk] - beta2[di - dk + ijk] - 
       beta2[-di + dk + ijk] + beta2[di + dk + ijk]) + 
    112896.*(-beta2[-di - 2*dk + ijk] + beta2[di - 2*dk + ijk] - 
       beta2[-2*di - dk + ijk] + beta2[2*di - dk + ijk] + 
       beta2[-2*di + dk + ijk] - beta2[2*di + dk + ijk] + 
       beta2[-di + 2*dk + ijk] - beta2[di + 2*dk + ijk]) - 
    21504.*(beta2[di - 3*dk + ijk] + beta2[3*di - dk + ijk] + 
       beta2[-3*di + dk + ijk] + beta2[-di + 3*dk + ijk]) + 
    21504.*(beta2[-di - 3*dk + ijk] + beta2[-3*di - dk + ijk] + 
       beta2[3*di + dk + ijk] + beta2[di + 3*dk + ijk]) + 
    5376.*(-beta2[-2*di - 3*dk + ijk] + beta2[2*di - 3*dk + ijk] - 
       beta2[-3*di - 2*dk + ijk] + beta2[3*di - 2*dk + ijk] + 
       beta2[-3*di + 2*dk + ijk] - beta2[3*di + 2*dk + ijk] + 
       beta2[-2*di + 3*dk + ijk] - beta2[2*di + 3*dk + ijk]) - 
    504.*(beta2[2*di - 4*dk + ijk] + beta2[4*di - 2*dk + ijk] + 
       beta2[-4*di + 2*dk + ijk] + beta2[-2*di + 4*dk + ijk]) + 
    2016.*(-beta2[-di - 4*dk + ijk] + beta2[di - 4*dk + ijk] - 
       beta2[-4*di - dk + ijk] + beta2[4*di - dk + ijk] + 
       beta2[-4*di + dk + ijk] - beta2[4*di + dk + ijk] + 
       beta2[-di + 4*dk + ijk] - beta2[di + 4*dk + ijk]) + 
    504.*(beta2[-2*di - 4*dk + ijk] + beta2[-4*di - 2*dk + ijk] + 
       beta2[4*di + 2*dk + ijk] + beta2[2*di + 4*dk + ijk]) + 
    96.*(-beta2[-3*di - 4*dk + ijk] + beta2[3*di - 4*dk + ijk] - 
       beta2[-4*di - 3*dk + ijk] + beta2[4*di - 3*dk + ijk] + 
       beta2[-4*di + 3*dk + ijk] - beta2[4*di + 3*dk + ijk] + 
       beta2[-3*di + 4*dk + ijk] - beta2[3*di + 4*dk + ijk]) - 
    28224.*(beta2[2*(di - dk) + ijk] + beta2[2*(-di + dk) + ijk]) - 
    1024.*(beta2[3*(di - dk) + ijk] + beta2[3*(-di + dk) + ijk]) - 
    9.*(beta2[4*(di - dk) + ijk] + beta2[4*(-di + dk) + ijk]) + 
    28224.*(beta2[-2*(di + dk) + ijk] + beta2[2*(di + dk) + ijk]) + 
    1024.*(beta2[-3*(di + dk) + ijk] + beta2[3*(di + dk) + ijk]) + 
    9.*(beta2[-4*(di + dk) + ijk] + beta2[4*(di + dk) + ijk]))
;

ddb133
=
5.6689342403628117914e-6*oo2dx*oo2dz*
  (451584.*(beta3[-di - dk + ijk] - beta3[di - dk + ijk] - 
       beta3[-di + dk + ijk] + beta3[di + dk + ijk]) + 
    112896.*(-beta3[-di - 2*dk + ijk] + beta3[di - 2*dk + ijk] - 
       beta3[-2*di - dk + ijk] + beta3[2*di - dk + ijk] + 
       beta3[-2*di + dk + ijk] - beta3[2*di + dk + ijk] + 
       beta3[-di + 2*dk + ijk] - beta3[di + 2*dk + ijk]) - 
    21504.*(beta3[di - 3*dk + ijk] + beta3[3*di - dk + ijk] + 
       beta3[-3*di + dk + ijk] + beta3[-di + 3*dk + ijk]) + 
    21504.*(beta3[-di - 3*dk + ijk] + beta3[-3*di - dk + ijk] + 
       beta3[3*di + dk + ijk] + beta3[di + 3*dk + ijk]) + 
    5376.*(-beta3[-2*di - 3*dk + ijk] + beta3[2*di - 3*dk + ijk] - 
       beta3[-3*di - 2*dk + ijk] + beta3[3*di - 2*dk + ijk] + 
       beta3[-3*di + 2*dk + ijk] - beta3[3*di + 2*dk + ijk] + 
       beta3[-2*di + 3*dk + ijk] - beta3[2*di + 3*dk + ijk]) - 
    504.*(beta3[2*di - 4*dk + ijk] + beta3[4*di - 2*dk + ijk] + 
       beta3[-4*di + 2*dk + ijk] + beta3[-2*di + 4*dk + ijk]) + 
    2016.*(-beta3[-di - 4*dk + ijk] + beta3[di - 4*dk + ijk] - 
       beta3[-4*di - dk + ijk] + beta3[4*di - dk + ijk] + 
       beta3[-4*di + dk + ijk] - beta3[4*di + dk + ijk] + 
       beta3[-di + 4*dk + ijk] - beta3[di + 4*dk + ijk]) + 
    504.*(beta3[-2*di - 4*dk + ijk] + beta3[-4*di - 2*dk + ijk] + 
       beta3[4*di + 2*dk + ijk] + beta3[2*di + 4*dk + ijk]) + 
    96.*(-beta3[-3*di - 4*dk + ijk] + beta3[3*di - 4*dk + ijk] - 
       beta3[-4*di - 3*dk + ijk] + beta3[4*di - 3*dk + ijk] + 
       beta3[-4*di + 3*dk + ijk] - beta3[4*di + 3*dk + ijk] + 
       beta3[-3*di + 4*dk + ijk] - beta3[3*di + 4*dk + ijk]) - 
    28224.*(beta3[2*(di - dk) + ijk] + beta3[2*(-di + dk) + ijk]) - 
    1024.*(beta3[3*(di - dk) + ijk] + beta3[3*(-di + dk) + ijk]) - 
    9.*(beta3[4*(di - dk) + ijk] + beta3[4*(-di + dk) + ijk]) + 
    28224.*(beta3[-2*(di + dk) + ijk] + beta3[2*(di + dk) + ijk]) + 
    1024.*(beta3[-3*(di + dk) + ijk] + beta3[3*(di + dk) + ijk]) + 
    9.*(beta3[-4*(di + dk) + ijk] + beta3[4*(di + dk) + ijk]))
;

ddb221
=
0.0001984126984126984127*oody2*
  (-14350.*beta1[ijk] + 8064.*(beta1[-dj + ijk] + beta1[dj + ijk]) - 
    1008.*(beta1[-2*dj + ijk] + beta1[2*dj + ijk]) + 
    128.*(beta1[-3*dj + ijk] + beta1[3*dj + ijk]) - 
    9.*(beta1[-4*dj + ijk] + beta1[4*dj + ijk]))
;

ddb222
=
0.0001984126984126984127*oody2*
  (-14350.*beta2[ijk] + 8064.*(beta2[-dj + ijk] + beta2[dj + ijk]) - 
    1008.*(beta2[-2*dj + ijk] + beta2[2*dj + ijk]) + 
    128.*(beta2[-3*dj + ijk] + beta2[3*dj + ijk]) - 
    9.*(beta2[-4*dj + ijk] + beta2[4*dj + ijk]))
;

ddb223
=
0.0001984126984126984127*oody2*
  (-14350.*beta3[ijk] + 8064.*(beta3[-dj + ijk] + beta3[dj + ijk]) - 
    1008.*(beta3[-2*dj + ijk] + beta3[2*dj + ijk]) + 
    128.*(beta3[-3*dj + ijk] + beta3[3*dj + ijk]) - 
    9.*(beta3[-4*dj + ijk] + beta3[4*dj + ijk]))
;

ddb231
=
5.6689342403628117914e-6*oo2dy*oo2dz*
  (451584.*(beta1[-dj - dk + ijk] - beta1[dj - dk + ijk] - 
       beta1[-dj + dk + ijk] + beta1[dj + dk + ijk]) + 
    112896.*(-beta1[-dj - 2*dk + ijk] + beta1[dj - 2*dk + ijk] - 
       beta1[-2*dj - dk + ijk] + beta1[2*dj - dk + ijk] + 
       beta1[-2*dj + dk + ijk] - beta1[2*dj + dk + ijk] + 
       beta1[-dj + 2*dk + ijk] - beta1[dj + 2*dk + ijk]) - 
    21504.*(beta1[dj - 3*dk + ijk] + beta1[3*dj - dk + ijk] + 
       beta1[-3*dj + dk + ijk] + beta1[-dj + 3*dk + ijk]) + 
    21504.*(beta1[-dj - 3*dk + ijk] + beta1[-3*dj - dk + ijk] + 
       beta1[3*dj + dk + ijk] + beta1[dj + 3*dk + ijk]) + 
    5376.*(-beta1[-2*dj - 3*dk + ijk] + beta1[2*dj - 3*dk + ijk] - 
       beta1[-3*dj - 2*dk + ijk] + beta1[3*dj - 2*dk + ijk] + 
       beta1[-3*dj + 2*dk + ijk] - beta1[3*dj + 2*dk + ijk] + 
       beta1[-2*dj + 3*dk + ijk] - beta1[2*dj + 3*dk + ijk]) - 
    504.*(beta1[2*dj - 4*dk + ijk] + beta1[4*dj - 2*dk + ijk] + 
       beta1[-4*dj + 2*dk + ijk] + beta1[-2*dj + 4*dk + ijk]) + 
    2016.*(-beta1[-dj - 4*dk + ijk] + beta1[dj - 4*dk + ijk] - 
       beta1[-4*dj - dk + ijk] + beta1[4*dj - dk + ijk] + 
       beta1[-4*dj + dk + ijk] - beta1[4*dj + dk + ijk] + 
       beta1[-dj + 4*dk + ijk] - beta1[dj + 4*dk + ijk]) + 
    504.*(beta1[-2*dj - 4*dk + ijk] + beta1[-4*dj - 2*dk + ijk] + 
       beta1[4*dj + 2*dk + ijk] + beta1[2*dj + 4*dk + ijk]) + 
    96.*(-beta1[-3*dj - 4*dk + ijk] + beta1[3*dj - 4*dk + ijk] - 
       beta1[-4*dj - 3*dk + ijk] + beta1[4*dj - 3*dk + ijk] + 
       beta1[-4*dj + 3*dk + ijk] - beta1[4*dj + 3*dk + ijk] + 
       beta1[-3*dj + 4*dk + ijk] - beta1[3*dj + 4*dk + ijk]) - 
    28224.*(beta1[2*(dj - dk) + ijk] + beta1[2*(-dj + dk) + ijk]) - 
    1024.*(beta1[3*(dj - dk) + ijk] + beta1[3*(-dj + dk) + ijk]) - 
    9.*(beta1[4*(dj - dk) + ijk] + beta1[4*(-dj + dk) + ijk]) + 
    28224.*(beta1[-2*(dj + dk) + ijk] + beta1[2*(dj + dk) + ijk]) + 
    1024.*(beta1[-3*(dj + dk) + ijk] + beta1[3*(dj + dk) + ijk]) + 
    9.*(beta1[-4*(dj + dk) + ijk] + beta1[4*(dj + dk) + ijk]))
;

ddb232
=
5.6689342403628117914e-6*oo2dy*oo2dz*
  (451584.*(beta2[-dj - dk + ijk] - beta2[dj - dk + ijk] - 
       beta2[-dj + dk + ijk] + beta2[dj + dk + ijk]) + 
    112896.*(-beta2[-dj - 2*dk + ijk] + beta2[dj - 2*dk + ijk] - 
       beta2[-2*dj - dk + ijk] + beta2[2*dj - dk + ijk] + 
       beta2[-2*dj + dk + ijk] - beta2[2*dj + dk + ijk] + 
       beta2[-dj + 2*dk + ijk] - beta2[dj + 2*dk + ijk]) - 
    21504.*(beta2[dj - 3*dk + ijk] + beta2[3*dj - dk + ijk] + 
       beta2[-3*dj + dk + ijk] + beta2[-dj + 3*dk + ijk]) + 
    21504.*(beta2[-dj - 3*dk + ijk] + beta2[-3*dj - dk + ijk] + 
       beta2[3*dj + dk + ijk] + beta2[dj + 3*dk + ijk]) + 
    5376.*(-beta2[-2*dj - 3*dk + ijk] + beta2[2*dj - 3*dk + ijk] - 
       beta2[-3*dj - 2*dk + ijk] + beta2[3*dj - 2*dk + ijk] + 
       beta2[-3*dj + 2*dk + ijk] - beta2[3*dj + 2*dk + ijk] + 
       beta2[-2*dj + 3*dk + ijk] - beta2[2*dj + 3*dk + ijk]) - 
    504.*(beta2[2*dj - 4*dk + ijk] + beta2[4*dj - 2*dk + ijk] + 
       beta2[-4*dj + 2*dk + ijk] + beta2[-2*dj + 4*dk + ijk]) + 
    2016.*(-beta2[-dj - 4*dk + ijk] + beta2[dj - 4*dk + ijk] - 
       beta2[-4*dj - dk + ijk] + beta2[4*dj - dk + ijk] + 
       beta2[-4*dj + dk + ijk] - beta2[4*dj + dk + ijk] + 
       beta2[-dj + 4*dk + ijk] - beta2[dj + 4*dk + ijk]) + 
    504.*(beta2[-2*dj - 4*dk + ijk] + beta2[-4*dj - 2*dk + ijk] + 
       beta2[4*dj + 2*dk + ijk] + beta2[2*dj + 4*dk + ijk]) + 
    96.*(-beta2[-3*dj - 4*dk + ijk] + beta2[3*dj - 4*dk + ijk] - 
       beta2[-4*dj - 3*dk + ijk] + beta2[4*dj - 3*dk + ijk] + 
       beta2[-4*dj + 3*dk + ijk] - beta2[4*dj + 3*dk + ijk] + 
       beta2[-3*dj + 4*dk + ijk] - beta2[3*dj + 4*dk + ijk]) - 
    28224.*(beta2[2*(dj - dk) + ijk] + beta2[2*(-dj + dk) + ijk]) - 
    1024.*(beta2[3*(dj - dk) + ijk] + beta2[3*(-dj + dk) + ijk]) - 
    9.*(beta2[4*(dj - dk) + ijk] + beta2[4*(-dj + dk) + ijk]) + 
    28224.*(beta2[-2*(dj + dk) + ijk] + beta2[2*(dj + dk) + ijk]) + 
    1024.*(beta2[-3*(dj + dk) + ijk] + beta2[3*(dj + dk) + ijk]) + 
    9.*(beta2[-4*(dj + dk) + ijk] + beta2[4*(dj + dk) + ijk]))
;

ddb233
=
5.6689342403628117914e-6*oo2dy*oo2dz*
  (451584.*(beta3[-dj - dk + ijk] - beta3[dj - dk + ijk] - 
       beta3[-dj + dk + ijk] + beta3[dj + dk + ijk]) + 
    112896.*(-beta3[-dj - 2*dk + ijk] + beta3[dj - 2*dk + ijk] - 
       beta3[-2*dj - dk + ijk] + beta3[2*dj - dk + ijk] + 
       beta3[-2*dj + dk + ijk] - beta3[2*dj + dk + ijk] + 
       beta3[-dj + 2*dk + ijk] - beta3[dj + 2*dk + ijk]) - 
    21504.*(beta3[dj - 3*dk + ijk] + beta3[3*dj - dk + ijk] + 
       beta3[-3*dj + dk + ijk] + beta3[-dj + 3*dk + ijk]) + 
    21504.*(beta3[-dj - 3*dk + ijk] + beta3[-3*dj - dk + ijk] + 
       beta3[3*dj + dk + ijk] + beta3[dj + 3*dk + ijk]) + 
    5376.*(-beta3[-2*dj - 3*dk + ijk] + beta3[2*dj - 3*dk + ijk] - 
       beta3[-3*dj - 2*dk + ijk] + beta3[3*dj - 2*dk + ijk] + 
       beta3[-3*dj + 2*dk + ijk] - beta3[3*dj + 2*dk + ijk] + 
       beta3[-2*dj + 3*dk + ijk] - beta3[2*dj + 3*dk + ijk]) - 
    504.*(beta3[2*dj - 4*dk + ijk] + beta3[4*dj - 2*dk + ijk] + 
       beta3[-4*dj + 2*dk + ijk] + beta3[-2*dj + 4*dk + ijk]) + 
    2016.*(-beta3[-dj - 4*dk + ijk] + beta3[dj - 4*dk + ijk] - 
       beta3[-4*dj - dk + ijk] + beta3[4*dj - dk + ijk] + 
       beta3[-4*dj + dk + ijk] - beta3[4*dj + dk + ijk] + 
       beta3[-dj + 4*dk + ijk] - beta3[dj + 4*dk + ijk]) + 
    504.*(beta3[-2*dj - 4*dk + ijk] + beta3[-4*dj - 2*dk + ijk] + 
       beta3[4*dj + 2*dk + ijk] + beta3[2*dj + 4*dk + ijk]) + 
    96.*(-beta3[-3*dj - 4*dk + ijk] + beta3[3*dj - 4*dk + ijk] - 
       beta3[-4*dj - 3*dk + ijk] + beta3[4*dj - 3*dk + ijk] + 
       beta3[-4*dj + 3*dk + ijk] - beta3[4*dj + 3*dk + ijk] + 
       beta3[-3*dj + 4*dk + ijk] - beta3[3*dj + 4*dk + ijk]) - 
    28224.*(beta3[2*(dj - dk) + ijk] + beta3[2*(-dj + dk) + ijk]) - 
    1024.*(beta3[3*(dj - dk) + ijk] + beta3[3*(-dj + dk) + ijk]) - 
    9.*(beta3[4*(dj - dk) + ijk] + beta3[4*(-dj + dk) + ijk]) + 
    28224.*(beta3[-2*(dj + dk) + ijk] + beta3[2*(dj + dk) + ijk]) + 
    1024.*(beta3[-3*(dj + dk) + ijk] + beta3[3*(dj + dk) + ijk]) + 
    9.*(beta3[-4*(dj + dk) + ijk] + beta3[4*(dj + dk) + ijk]))
;

ddb331
=
0.0001984126984126984127*oodz2*
  (-14350.*beta1[ijk] + 8064.*(beta1[-dk + ijk] + beta1[dk + ijk]) - 
    1008.*(beta1[-2*dk + ijk] + beta1[2*dk + ijk]) + 
    128.*(beta1[-3*dk + ijk] + beta1[3*dk + ijk]) - 
    9.*(beta1[-4*dk + ijk] + beta1[4*dk + ijk]))
;

ddb332
=
0.0001984126984126984127*oodz2*
  (-14350.*beta2[ijk] + 8064.*(beta2[-dk + ijk] + beta2[dk + ijk]) - 
    1008.*(beta2[-2*dk + ijk] + beta2[2*dk + ijk]) + 
    128.*(beta2[-3*dk + ijk] + beta2[3*dk + ijk]) - 
    9.*(beta2[-4*dk + ijk] + beta2[4*dk + ijk]))
;

ddb333
=
0.0001984126984126984127*oodz2*
  (-14350.*beta3[ijk] + 8064.*(beta3[-dk + ijk] + beta3[dk + ijk]) - 
    1008.*(beta3[-2*dk + ijk] + beta3[2*dk + ijk]) + 
    128.*(beta3[-3*dk + ijk] + beta3[3*dk + ijk]) - 
    9.*(beta3[-4*dk + ijk] + beta3[4*dk + ijk]))
;

delg111
=
0.0023809523809523809524*oo2dx*
  (3.*g11[-4*di + ijk] + 168.*g11[-2*di + ijk] + 
    672.*(-g11[-di + ijk] + g11[di + ijk]) - 168.*g11[2*di + ijk] + 
    32.*(-g11[-3*di + ijk] + g11[3*di + ijk]) - 3.*g11[4*di + ijk])
;

delg112
=
0.0023809523809523809524*oo2dx*
  (3.*g12[-4*di + ijk] + 168.*g12[-2*di + ijk] + 
    672.*(-g12[-di + ijk] + g12[di + ijk]) - 168.*g12[2*di + ijk] + 
    32.*(-g12[-3*di + ijk] + g12[3*di + ijk]) - 3.*g12[4*di + ijk])
;

delg113
=
0.0023809523809523809524*oo2dx*
  (3.*g13[-4*di + ijk] + 168.*g13[-2*di + ijk] + 
    672.*(-g13[-di + ijk] + g13[di + ijk]) - 168.*g13[2*di + ijk] + 
    32.*(-g13[-3*di + ijk] + g13[3*di + ijk]) - 3.*g13[4*di + ijk])
;

delg122
=
0.0023809523809523809524*oo2dx*
  (3.*g22[-4*di + ijk] + 168.*g22[-2*di + ijk] + 
    672.*(-g22[-di + ijk] + g22[di + ijk]) - 168.*g22[2*di + ijk] + 
    32.*(-g22[-3*di + ijk] + g22[3*di + ijk]) - 3.*g22[4*di + ijk])
;

delg123
=
0.0023809523809523809524*oo2dx*
  (3.*g23[-4*di + ijk] + 168.*g23[-2*di + ijk] + 
    672.*(-g23[-di + ijk] + g23[di + ijk]) - 168.*g23[2*di + ijk] + 
    32.*(-g23[-3*di + ijk] + g23[3*di + ijk]) - 3.*g23[4*di + ijk])
;

delg133
=
0.0023809523809523809524*oo2dx*
  (3.*g33[-4*di + ijk] + 168.*g33[-2*di + ijk] + 
    672.*(-g33[-di + ijk] + g33[di + ijk]) - 168.*g33[2*di + ijk] + 
    32.*(-g33[-3*di + ijk] + g33[3*di + ijk]) - 3.*g33[4*di + ijk])
;

delg211
=
0.0023809523809523809524*oo2dy*
  (3.*g11[-4*dj + ijk] + 168.*g11[-2*dj + ijk] + 
    672.*(-g11[-dj + ijk] + g11[dj + ijk]) - 168.*g11[2*dj + ijk] + 
    32.*(-g11[-3*dj + ijk] + g11[3*dj + ijk]) - 3.*g11[4*dj + ijk])
;

delg212
=
0.0023809523809523809524*oo2dy*
  (3.*g12[-4*dj + ijk] + 168.*g12[-2*dj + ijk] + 
    672.*(-g12[-dj + ijk] + g12[dj + ijk]) - 168.*g12[2*dj + ijk] + 
    32.*(-g12[-3*dj + ijk] + g12[3*dj + ijk]) - 3.*g12[4*dj + ijk])
;

delg213
=
0.0023809523809523809524*oo2dy*
  (3.*g13[-4*dj + ijk] + 168.*g13[-2*dj + ijk] + 
    672.*(-g13[-dj + ijk] + g13[dj + ijk]) - 168.*g13[2*dj + ijk] + 
    32.*(-g13[-3*dj + ijk] + g13[3*dj + ijk]) - 3.*g13[4*dj + ijk])
;

delg222
=
0.0023809523809523809524*oo2dy*
  (3.*g22[-4*dj + ijk] + 168.*g22[-2*dj + ijk] + 
    672.*(-g22[-dj + ijk] + g22[dj + ijk]) - 168.*g22[2*dj + ijk] + 
    32.*(-g22[-3*dj + ijk] + g22[3*dj + ijk]) - 3.*g22[4*dj + ijk])
;

delg223
=
0.0023809523809523809524*oo2dy*
  (3.*g23[-4*dj + ijk] + 168.*g23[-2*dj + ijk] + 
    672.*(-g23[-dj + ijk] + g23[dj + ijk]) - 168.*g23[2*dj + ijk] + 
    32.*(-g23[-3*dj + ijk] + g23[3*dj + ijk]) - 3.*g23[4*dj + ijk])
;

delg233
=
0.0023809523809523809524*oo2dy*
  (3.*g33[-4*dj + ijk] + 168.*g33[-2*dj + ijk] + 
    672.*(-g33[-dj + ijk] + g33[dj + ijk]) - 168.*g33[2*dj + ijk] + 
    32.*(-g33[-3*dj + ijk] + g33[3*dj + ijk]) - 3.*g33[4*dj + ijk])
;

delg311
=
0.0023809523809523809524*oo2dz*
  (3.*g11[-4*dk + ijk] + 168.*g11[-2*dk + ijk] + 
    672.*(-g11[-dk + ijk] + g11[dk + ijk]) - 168.*g11[2*dk + ijk] + 
    32.*(-g11[-3*dk + ijk] + g11[3*dk + ijk]) - 3.*g11[4*dk + ijk])
;

delg312
=
0.0023809523809523809524*oo2dz*
  (3.*g12[-4*dk + ijk] + 168.*g12[-2*dk + ijk] + 
    672.*(-g12[-dk + ijk] + g12[dk + ijk]) - 168.*g12[2*dk + ijk] + 
    32.*(-g12[-3*dk + ijk] + g12[3*dk + ijk]) - 3.*g12[4*dk + ijk])
;

delg313
=
0.0023809523809523809524*oo2dz*
  (3.*g13[-4*dk + ijk] + 168.*g13[-2*dk + ijk] + 
    672.*(-g13[-dk + ijk] + g13[dk + ijk]) - 168.*g13[2*dk + ijk] + 
    32.*(-g13[-3*dk + ijk] + g13[3*dk + ijk]) - 3.*g13[4*dk + ijk])
;

delg322
=
0.0023809523809523809524*oo2dz*
  (3.*g22[-4*dk + ijk] + 168.*g22[-2*dk + ijk] + 
    672.*(-g22[-dk + ijk] + g22[dk + ijk]) - 168.*g22[2*dk + ijk] + 
    32.*(-g22[-3*dk + ijk] + g22[3*dk + ijk]) - 3.*g22[4*dk + ijk])
;

delg323
=
0.0023809523809523809524*oo2dz*
  (3.*g23[-4*dk + ijk] + 168.*g23[-2*dk + ijk] + 
    672.*(-g23[-dk + ijk] + g23[dk + ijk]) - 168.*g23[2*dk + ijk] + 
    32.*(-g23[-3*dk + ijk] + g23[3*dk + ijk]) - 3.*g23[4*dk + ijk])
;

delg333
=
0.0023809523809523809524*oo2dz*
  (3.*g33[-4*dk + ijk] + 168.*g33[-2*dk + ijk] + 
    672.*(-g33[-dk + ijk] + g33[dk + ijk]) - 168.*g33[2*dk + ijk] + 
    32.*(-g33[-3*dk + ijk] + g33[3*dk + ijk]) - 3.*g33[4*dk + ijk])
;

deldelg1111
=
0.0001984126984126984127*oodx2*
  (-14350.*g11[ijk] + 8064.*(g11[-di + ijk] + g11[di + ijk]) - 
    1008.*(g11[-2*di + ijk] + g11[2*di + ijk]) + 
    128.*(g11[-3*di + ijk] + g11[3*di + ijk]) - 
    9.*(g11[-4*di + ijk] + g11[4*di + ijk]))
;

deldelg1112
=
0.0001984126984126984127*oodx2*
  (-14350.*g12[ijk] + 8064.*(g12[-di + ijk] + g12[di + ijk]) - 
    1008.*(g12[-2*di + ijk] + g12[2*di + ijk]) + 
    128.*(g12[-3*di + ijk] + g12[3*di + ijk]) - 
    9.*(g12[-4*di + ijk] + g12[4*di + ijk]))
;

deldelg1113
=
0.0001984126984126984127*oodx2*
  (-14350.*g13[ijk] + 8064.*(g13[-di + ijk] + g13[di + ijk]) - 
    1008.*(g13[-2*di + ijk] + g13[2*di + ijk]) + 
    128.*(g13[-3*di + ijk] + g13[3*di + ijk]) - 
    9.*(g13[-4*di + ijk] + g13[4*di + ijk]))
;

deldelg1122
=
0.0001984126984126984127*oodx2*
  (-14350.*g22[ijk] + 8064.*(g22[-di + ijk] + g22[di + ijk]) - 
    1008.*(g22[-2*di + ijk] + g22[2*di + ijk]) + 
    128.*(g22[-3*di + ijk] + g22[3*di + ijk]) - 
    9.*(g22[-4*di + ijk] + g22[4*di + ijk]))
;

deldelg1123
=
0.0001984126984126984127*oodx2*
  (-14350.*g23[ijk] + 8064.*(g23[-di + ijk] + g23[di + ijk]) - 
    1008.*(g23[-2*di + ijk] + g23[2*di + ijk]) + 
    128.*(g23[-3*di + ijk] + g23[3*di + ijk]) - 
    9.*(g23[-4*di + ijk] + g23[4*di + ijk]))
;

deldelg1133
=
0.0001984126984126984127*oodx2*
  (-14350.*g33[ijk] + 8064.*(g33[-di + ijk] + g33[di + ijk]) - 
    1008.*(g33[-2*di + ijk] + g33[2*di + ijk]) + 
    128.*(g33[-3*di + ijk] + g33[3*di + ijk]) - 
    9.*(g33[-4*di + ijk] + g33[4*di + ijk]))
;

deldelg1211
=
5.6689342403628117914e-6*oo2dx*oo2dy*
  (451584.*(g11[-di - dj + ijk] - g11[di - dj + ijk] - 
       g11[-di + dj + ijk] + g11[di + dj + ijk]) + 
    112896.*(-g11[-di - 2*dj + ijk] + g11[di - 2*dj + ijk] - 
       g11[-2*di - dj + ijk] + g11[2*di - dj + ijk] + 
       g11[-2*di + dj + ijk] - g11[2*di + dj + ijk] + 
       g11[-di + 2*dj + ijk] - g11[di + 2*dj + ijk]) - 
    21504.*(g11[di - 3*dj + ijk] + g11[3*di - dj + ijk] + 
       g11[-3*di + dj + ijk] + g11[-di + 3*dj + ijk]) + 
    21504.*(g11[-di - 3*dj + ijk] + g11[-3*di - dj + ijk] + 
       g11[3*di + dj + ijk] + g11[di + 3*dj + ijk]) + 
    5376.*(-g11[-2*di - 3*dj + ijk] + g11[2*di - 3*dj + ijk] - 
       g11[-3*di - 2*dj + ijk] + g11[3*di - 2*dj + ijk] + 
       g11[-3*di + 2*dj + ijk] - g11[3*di + 2*dj + ijk] + 
       g11[-2*di + 3*dj + ijk] - g11[2*di + 3*dj + ijk]) - 
    504.*(g11[2*di - 4*dj + ijk] + g11[4*di - 2*dj + ijk] + 
       g11[-4*di + 2*dj + ijk] + g11[-2*di + 4*dj + ijk]) + 
    2016.*(-g11[-di - 4*dj + ijk] + g11[di - 4*dj + ijk] - 
       g11[-4*di - dj + ijk] + g11[4*di - dj + ijk] + 
       g11[-4*di + dj + ijk] - g11[4*di + dj + ijk] + 
       g11[-di + 4*dj + ijk] - g11[di + 4*dj + ijk]) + 
    504.*(g11[-2*di - 4*dj + ijk] + g11[-4*di - 2*dj + ijk] + 
       g11[4*di + 2*dj + ijk] + g11[2*di + 4*dj + ijk]) + 
    96.*(-g11[-3*di - 4*dj + ijk] + g11[3*di - 4*dj + ijk] - 
       g11[-4*di - 3*dj + ijk] + g11[4*di - 3*dj + ijk] + 
       g11[-4*di + 3*dj + ijk] - g11[4*di + 3*dj + ijk] + 
       g11[-3*di + 4*dj + ijk] - g11[3*di + 4*dj + ijk]) - 
    28224.*(g11[2*(di - dj) + ijk] + g11[2*(-di + dj) + ijk]) - 
    1024.*(g11[3*(di - dj) + ijk] + g11[3*(-di + dj) + ijk]) - 
    9.*(g11[4*(di - dj) + ijk] + g11[4*(-di + dj) + ijk]) + 
    28224.*(g11[-2*(di + dj) + ijk] + g11[2*(di + dj) + ijk]) + 
    1024.*(g11[-3*(di + dj) + ijk] + g11[3*(di + dj) + ijk]) + 
    9.*(g11[-4*(di + dj) + ijk] + g11[4*(di + dj) + ijk]))
;

deldelg1212
=
5.6689342403628117914e-6*oo2dx*oo2dy*
  (451584.*(g12[-di - dj + ijk] - g12[di - dj + ijk] - 
       g12[-di + dj + ijk] + g12[di + dj + ijk]) + 
    112896.*(-g12[-di - 2*dj + ijk] + g12[di - 2*dj + ijk] - 
       g12[-2*di - dj + ijk] + g12[2*di - dj + ijk] + 
       g12[-2*di + dj + ijk] - g12[2*di + dj + ijk] + 
       g12[-di + 2*dj + ijk] - g12[di + 2*dj + ijk]) - 
    21504.*(g12[di - 3*dj + ijk] + g12[3*di - dj + ijk] + 
       g12[-3*di + dj + ijk] + g12[-di + 3*dj + ijk]) + 
    21504.*(g12[-di - 3*dj + ijk] + g12[-3*di - dj + ijk] + 
       g12[3*di + dj + ijk] + g12[di + 3*dj + ijk]) + 
    5376.*(-g12[-2*di - 3*dj + ijk] + g12[2*di - 3*dj + ijk] - 
       g12[-3*di - 2*dj + ijk] + g12[3*di - 2*dj + ijk] + 
       g12[-3*di + 2*dj + ijk] - g12[3*di + 2*dj + ijk] + 
       g12[-2*di + 3*dj + ijk] - g12[2*di + 3*dj + ijk]) - 
    504.*(g12[2*di - 4*dj + ijk] + g12[4*di - 2*dj + ijk] + 
       g12[-4*di + 2*dj + ijk] + g12[-2*di + 4*dj + ijk]) + 
    2016.*(-g12[-di - 4*dj + ijk] + g12[di - 4*dj + ijk] - 
       g12[-4*di - dj + ijk] + g12[4*di - dj + ijk] + 
       g12[-4*di + dj + ijk] - g12[4*di + dj + ijk] + 
       g12[-di + 4*dj + ijk] - g12[di + 4*dj + ijk]) + 
    504.*(g12[-2*di - 4*dj + ijk] + g12[-4*di - 2*dj + ijk] + 
       g12[4*di + 2*dj + ijk] + g12[2*di + 4*dj + ijk]) + 
    96.*(-g12[-3*di - 4*dj + ijk] + g12[3*di - 4*dj + ijk] - 
       g12[-4*di - 3*dj + ijk] + g12[4*di - 3*dj + ijk] + 
       g12[-4*di + 3*dj + ijk] - g12[4*di + 3*dj + ijk] + 
       g12[-3*di + 4*dj + ijk] - g12[3*di + 4*dj + ijk]) - 
    28224.*(g12[2*(di - dj) + ijk] + g12[2*(-di + dj) + ijk]) - 
    1024.*(g12[3*(di - dj) + ijk] + g12[3*(-di + dj) + ijk]) - 
    9.*(g12[4*(di - dj) + ijk] + g12[4*(-di + dj) + ijk]) + 
    28224.*(g12[-2*(di + dj) + ijk] + g12[2*(di + dj) + ijk]) + 
    1024.*(g12[-3*(di + dj) + ijk] + g12[3*(di + dj) + ijk]) + 
    9.*(g12[-4*(di + dj) + ijk] + g12[4*(di + dj) + ijk]))
;

deldelg1213
=
5.6689342403628117914e-6*oo2dx*oo2dy*
  (451584.*(g13[-di - dj + ijk] - g13[di - dj + ijk] - 
       g13[-di + dj + ijk] + g13[di + dj + ijk]) + 
    112896.*(-g13[-di - 2*dj + ijk] + g13[di - 2*dj + ijk] - 
       g13[-2*di - dj + ijk] + g13[2*di - dj + ijk] + 
       g13[-2*di + dj + ijk] - g13[2*di + dj + ijk] + 
       g13[-di + 2*dj + ijk] - g13[di + 2*dj + ijk]) - 
    21504.*(g13[di - 3*dj + ijk] + g13[3*di - dj + ijk] + 
       g13[-3*di + dj + ijk] + g13[-di + 3*dj + ijk]) + 
    21504.*(g13[-di - 3*dj + ijk] + g13[-3*di - dj + ijk] + 
       g13[3*di + dj + ijk] + g13[di + 3*dj + ijk]) + 
    5376.*(-g13[-2*di - 3*dj + ijk] + g13[2*di - 3*dj + ijk] - 
       g13[-3*di - 2*dj + ijk] + g13[3*di - 2*dj + ijk] + 
       g13[-3*di + 2*dj + ijk] - g13[3*di + 2*dj + ijk] + 
       g13[-2*di + 3*dj + ijk] - g13[2*di + 3*dj + ijk]) - 
    504.*(g13[2*di - 4*dj + ijk] + g13[4*di - 2*dj + ijk] + 
       g13[-4*di + 2*dj + ijk] + g13[-2*di + 4*dj + ijk]) + 
    2016.*(-g13[-di - 4*dj + ijk] + g13[di - 4*dj + ijk] - 
       g13[-4*di - dj + ijk] + g13[4*di - dj + ijk] + 
       g13[-4*di + dj + ijk] - g13[4*di + dj + ijk] + 
       g13[-di + 4*dj + ijk] - g13[di + 4*dj + ijk]) + 
    504.*(g13[-2*di - 4*dj + ijk] + g13[-4*di - 2*dj + ijk] + 
       g13[4*di + 2*dj + ijk] + g13[2*di + 4*dj + ijk]) + 
    96.*(-g13[-3*di - 4*dj + ijk] + g13[3*di - 4*dj + ijk] - 
       g13[-4*di - 3*dj + ijk] + g13[4*di - 3*dj + ijk] + 
       g13[-4*di + 3*dj + ijk] - g13[4*di + 3*dj + ijk] + 
       g13[-3*di + 4*dj + ijk] - g13[3*di + 4*dj + ijk]) - 
    28224.*(g13[2*(di - dj) + ijk] + g13[2*(-di + dj) + ijk]) - 
    1024.*(g13[3*(di - dj) + ijk] + g13[3*(-di + dj) + ijk]) - 
    9.*(g13[4*(di - dj) + ijk] + g13[4*(-di + dj) + ijk]) + 
    28224.*(g13[-2*(di + dj) + ijk] + g13[2*(di + dj) + ijk]) + 
    1024.*(g13[-3*(di + dj) + ijk] + g13[3*(di + dj) + ijk]) + 
    9.*(g13[-4*(di + dj) + ijk] + g13[4*(di + dj) + ijk]))
;

deldelg1222
=
5.6689342403628117914e-6*oo2dx*oo2dy*
  (451584.*(g22[-di - dj + ijk] - g22[di - dj + ijk] - 
       g22[-di + dj + ijk] + g22[di + dj + ijk]) + 
    112896.*(-g22[-di - 2*dj + ijk] + g22[di - 2*dj + ijk] - 
       g22[-2*di - dj + ijk] + g22[2*di - dj + ijk] + 
       g22[-2*di + dj + ijk] - g22[2*di + dj + ijk] + 
       g22[-di + 2*dj + ijk] - g22[di + 2*dj + ijk]) - 
    21504.*(g22[di - 3*dj + ijk] + g22[3*di - dj + ijk] + 
       g22[-3*di + dj + ijk] + g22[-di + 3*dj + ijk]) + 
    21504.*(g22[-di - 3*dj + ijk] + g22[-3*di - dj + ijk] + 
       g22[3*di + dj + ijk] + g22[di + 3*dj + ijk]) + 
    5376.*(-g22[-2*di - 3*dj + ijk] + g22[2*di - 3*dj + ijk] - 
       g22[-3*di - 2*dj + ijk] + g22[3*di - 2*dj + ijk] + 
       g22[-3*di + 2*dj + ijk] - g22[3*di + 2*dj + ijk] + 
       g22[-2*di + 3*dj + ijk] - g22[2*di + 3*dj + ijk]) - 
    504.*(g22[2*di - 4*dj + ijk] + g22[4*di - 2*dj + ijk] + 
       g22[-4*di + 2*dj + ijk] + g22[-2*di + 4*dj + ijk]) + 
    2016.*(-g22[-di - 4*dj + ijk] + g22[di - 4*dj + ijk] - 
       g22[-4*di - dj + ijk] + g22[4*di - dj + ijk] + 
       g22[-4*di + dj + ijk] - g22[4*di + dj + ijk] + 
       g22[-di + 4*dj + ijk] - g22[di + 4*dj + ijk]) + 
    504.*(g22[-2*di - 4*dj + ijk] + g22[-4*di - 2*dj + ijk] + 
       g22[4*di + 2*dj + ijk] + g22[2*di + 4*dj + ijk]) + 
    96.*(-g22[-3*di - 4*dj + ijk] + g22[3*di - 4*dj + ijk] - 
       g22[-4*di - 3*dj + ijk] + g22[4*di - 3*dj + ijk] + 
       g22[-4*di + 3*dj + ijk] - g22[4*di + 3*dj + ijk] + 
       g22[-3*di + 4*dj + ijk] - g22[3*di + 4*dj + ijk]) - 
    28224.*(g22[2*(di - dj) + ijk] + g22[2*(-di + dj) + ijk]) - 
    1024.*(g22[3*(di - dj) + ijk] + g22[3*(-di + dj) + ijk]) - 
    9.*(g22[4*(di - dj) + ijk] + g22[4*(-di + dj) + ijk]) + 
    28224.*(g22[-2*(di + dj) + ijk] + g22[2*(di + dj) + ijk]) + 
    1024.*(g22[-3*(di + dj) + ijk] + g22[3*(di + dj) + ijk]) + 
    9.*(g22[-4*(di + dj) + ijk] + g22[4*(di + dj) + ijk]))
;

deldelg1223
=
5.6689342403628117914e-6*oo2dx*oo2dy*
  (451584.*(g23[-di - dj + ijk] - g23[di - dj + ijk] - 
       g23[-di + dj + ijk] + g23[di + dj + ijk]) + 
    112896.*(-g23[-di - 2*dj + ijk] + g23[di - 2*dj + ijk] - 
       g23[-2*di - dj + ijk] + g23[2*di - dj + ijk] + 
       g23[-2*di + dj + ijk] - g23[2*di + dj + ijk] + 
       g23[-di + 2*dj + ijk] - g23[di + 2*dj + ijk]) - 
    21504.*(g23[di - 3*dj + ijk] + g23[3*di - dj + ijk] + 
       g23[-3*di + dj + ijk] + g23[-di + 3*dj + ijk]) + 
    21504.*(g23[-di - 3*dj + ijk] + g23[-3*di - dj + ijk] + 
       g23[3*di + dj + ijk] + g23[di + 3*dj + ijk]) + 
    5376.*(-g23[-2*di - 3*dj + ijk] + g23[2*di - 3*dj + ijk] - 
       g23[-3*di - 2*dj + ijk] + g23[3*di - 2*dj + ijk] + 
       g23[-3*di + 2*dj + ijk] - g23[3*di + 2*dj + ijk] + 
       g23[-2*di + 3*dj + ijk] - g23[2*di + 3*dj + ijk]) - 
    504.*(g23[2*di - 4*dj + ijk] + g23[4*di - 2*dj + ijk] + 
       g23[-4*di + 2*dj + ijk] + g23[-2*di + 4*dj + ijk]) + 
    2016.*(-g23[-di - 4*dj + ijk] + g23[di - 4*dj + ijk] - 
       g23[-4*di - dj + ijk] + g23[4*di - dj + ijk] + 
       g23[-4*di + dj + ijk] - g23[4*di + dj + ijk] + 
       g23[-di + 4*dj + ijk] - g23[di + 4*dj + ijk]) + 
    504.*(g23[-2*di - 4*dj + ijk] + g23[-4*di - 2*dj + ijk] + 
       g23[4*di + 2*dj + ijk] + g23[2*di + 4*dj + ijk]) + 
    96.*(-g23[-3*di - 4*dj + ijk] + g23[3*di - 4*dj + ijk] - 
       g23[-4*di - 3*dj + ijk] + g23[4*di - 3*dj + ijk] + 
       g23[-4*di + 3*dj + ijk] - g23[4*di + 3*dj + ijk] + 
       g23[-3*di + 4*dj + ijk] - g23[3*di + 4*dj + ijk]) - 
    28224.*(g23[2*(di - dj) + ijk] + g23[2*(-di + dj) + ijk]) - 
    1024.*(g23[3*(di - dj) + ijk] + g23[3*(-di + dj) + ijk]) - 
    9.*(g23[4*(di - dj) + ijk] + g23[4*(-di + dj) + ijk]) + 
    28224.*(g23[-2*(di + dj) + ijk] + g23[2*(di + dj) + ijk]) + 
    1024.*(g23[-3*(di + dj) + ijk] + g23[3*(di + dj) + ijk]) + 
    9.*(g23[-4*(di + dj) + ijk] + g23[4*(di + dj) + ijk]))
;

deldelg1233
=
5.6689342403628117914e-6*oo2dx*oo2dy*
  (451584.*(g33[-di - dj + ijk] - g33[di - dj + ijk] - 
       g33[-di + dj + ijk] + g33[di + dj + ijk]) + 
    112896.*(-g33[-di - 2*dj + ijk] + g33[di - 2*dj + ijk] - 
       g33[-2*di - dj + ijk] + g33[2*di - dj + ijk] + 
       g33[-2*di + dj + ijk] - g33[2*di + dj + ijk] + 
       g33[-di + 2*dj + ijk] - g33[di + 2*dj + ijk]) - 
    21504.*(g33[di - 3*dj + ijk] + g33[3*di - dj + ijk] + 
       g33[-3*di + dj + ijk] + g33[-di + 3*dj + ijk]) + 
    21504.*(g33[-di - 3*dj + ijk] + g33[-3*di - dj + ijk] + 
       g33[3*di + dj + ijk] + g33[di + 3*dj + ijk]) + 
    5376.*(-g33[-2*di - 3*dj + ijk] + g33[2*di - 3*dj + ijk] - 
       g33[-3*di - 2*dj + ijk] + g33[3*di - 2*dj + ijk] + 
       g33[-3*di + 2*dj + ijk] - g33[3*di + 2*dj + ijk] + 
       g33[-2*di + 3*dj + ijk] - g33[2*di + 3*dj + ijk]) - 
    504.*(g33[2*di - 4*dj + ijk] + g33[4*di - 2*dj + ijk] + 
       g33[-4*di + 2*dj + ijk] + g33[-2*di + 4*dj + ijk]) + 
    2016.*(-g33[-di - 4*dj + ijk] + g33[di - 4*dj + ijk] - 
       g33[-4*di - dj + ijk] + g33[4*di - dj + ijk] + 
       g33[-4*di + dj + ijk] - g33[4*di + dj + ijk] + 
       g33[-di + 4*dj + ijk] - g33[di + 4*dj + ijk]) + 
    504.*(g33[-2*di - 4*dj + ijk] + g33[-4*di - 2*dj + ijk] + 
       g33[4*di + 2*dj + ijk] + g33[2*di + 4*dj + ijk]) + 
    96.*(-g33[-3*di - 4*dj + ijk] + g33[3*di - 4*dj + ijk] - 
       g33[-4*di - 3*dj + ijk] + g33[4*di - 3*dj + ijk] + 
       g33[-4*di + 3*dj + ijk] - g33[4*di + 3*dj + ijk] + 
       g33[-3*di + 4*dj + ijk] - g33[3*di + 4*dj + ijk]) - 
    28224.*(g33[2*(di - dj) + ijk] + g33[2*(-di + dj) + ijk]) - 
    1024.*(g33[3*(di - dj) + ijk] + g33[3*(-di + dj) + ijk]) - 
    9.*(g33[4*(di - dj) + ijk] + g33[4*(-di + dj) + ijk]) + 
    28224.*(g33[-2*(di + dj) + ijk] + g33[2*(di + dj) + ijk]) + 
    1024.*(g33[-3*(di + dj) + ijk] + g33[3*(di + dj) + ijk]) + 
    9.*(g33[-4*(di + dj) + ijk] + g33[4*(di + dj) + ijk]))
;

deldelg1311
=
5.6689342403628117914e-6*oo2dx*oo2dz*
  (451584.*(g11[-di - dk + ijk] - g11[di - dk + ijk] - 
       g11[-di + dk + ijk] + g11[di + dk + ijk]) + 
    112896.*(-g11[-di - 2*dk + ijk] + g11[di - 2*dk + ijk] - 
       g11[-2*di - dk + ijk] + g11[2*di - dk + ijk] + 
       g11[-2*di + dk + ijk] - g11[2*di + dk + ijk] + 
       g11[-di + 2*dk + ijk] - g11[di + 2*dk + ijk]) - 
    21504.*(g11[di - 3*dk + ijk] + g11[3*di - dk + ijk] + 
       g11[-3*di + dk + ijk] + g11[-di + 3*dk + ijk]) + 
    21504.*(g11[-di - 3*dk + ijk] + g11[-3*di - dk + ijk] + 
       g11[3*di + dk + ijk] + g11[di + 3*dk + ijk]) + 
    5376.*(-g11[-2*di - 3*dk + ijk] + g11[2*di - 3*dk + ijk] - 
       g11[-3*di - 2*dk + ijk] + g11[3*di - 2*dk + ijk] + 
       g11[-3*di + 2*dk + ijk] - g11[3*di + 2*dk + ijk] + 
       g11[-2*di + 3*dk + ijk] - g11[2*di + 3*dk + ijk]) - 
    504.*(g11[2*di - 4*dk + ijk] + g11[4*di - 2*dk + ijk] + 
       g11[-4*di + 2*dk + ijk] + g11[-2*di + 4*dk + ijk]) + 
    2016.*(-g11[-di - 4*dk + ijk] + g11[di - 4*dk + ijk] - 
       g11[-4*di - dk + ijk] + g11[4*di - dk + ijk] + 
       g11[-4*di + dk + ijk] - g11[4*di + dk + ijk] + 
       g11[-di + 4*dk + ijk] - g11[di + 4*dk + ijk]) + 
    504.*(g11[-2*di - 4*dk + ijk] + g11[-4*di - 2*dk + ijk] + 
       g11[4*di + 2*dk + ijk] + g11[2*di + 4*dk + ijk]) + 
    96.*(-g11[-3*di - 4*dk + ijk] + g11[3*di - 4*dk + ijk] - 
       g11[-4*di - 3*dk + ijk] + g11[4*di - 3*dk + ijk] + 
       g11[-4*di + 3*dk + ijk] - g11[4*di + 3*dk + ijk] + 
       g11[-3*di + 4*dk + ijk] - g11[3*di + 4*dk + ijk]) - 
    28224.*(g11[2*(di - dk) + ijk] + g11[2*(-di + dk) + ijk]) - 
    1024.*(g11[3*(di - dk) + ijk] + g11[3*(-di + dk) + ijk]) - 
    9.*(g11[4*(di - dk) + ijk] + g11[4*(-di + dk) + ijk]) + 
    28224.*(g11[-2*(di + dk) + ijk] + g11[2*(di + dk) + ijk]) + 
    1024.*(g11[-3*(di + dk) + ijk] + g11[3*(di + dk) + ijk]) + 
    9.*(g11[-4*(di + dk) + ijk] + g11[4*(di + dk) + ijk]))
;

deldelg1312
=
5.6689342403628117914e-6*oo2dx*oo2dz*
  (451584.*(g12[-di - dk + ijk] - g12[di - dk + ijk] - 
       g12[-di + dk + ijk] + g12[di + dk + ijk]) + 
    112896.*(-g12[-di - 2*dk + ijk] + g12[di - 2*dk + ijk] - 
       g12[-2*di - dk + ijk] + g12[2*di - dk + ijk] + 
       g12[-2*di + dk + ijk] - g12[2*di + dk + ijk] + 
       g12[-di + 2*dk + ijk] - g12[di + 2*dk + ijk]) - 
    21504.*(g12[di - 3*dk + ijk] + g12[3*di - dk + ijk] + 
       g12[-3*di + dk + ijk] + g12[-di + 3*dk + ijk]) + 
    21504.*(g12[-di - 3*dk + ijk] + g12[-3*di - dk + ijk] + 
       g12[3*di + dk + ijk] + g12[di + 3*dk + ijk]) + 
    5376.*(-g12[-2*di - 3*dk + ijk] + g12[2*di - 3*dk + ijk] - 
       g12[-3*di - 2*dk + ijk] + g12[3*di - 2*dk + ijk] + 
       g12[-3*di + 2*dk + ijk] - g12[3*di + 2*dk + ijk] + 
       g12[-2*di + 3*dk + ijk] - g12[2*di + 3*dk + ijk]) - 
    504.*(g12[2*di - 4*dk + ijk] + g12[4*di - 2*dk + ijk] + 
       g12[-4*di + 2*dk + ijk] + g12[-2*di + 4*dk + ijk]) + 
    2016.*(-g12[-di - 4*dk + ijk] + g12[di - 4*dk + ijk] - 
       g12[-4*di - dk + ijk] + g12[4*di - dk + ijk] + 
       g12[-4*di + dk + ijk] - g12[4*di + dk + ijk] + 
       g12[-di + 4*dk + ijk] - g12[di + 4*dk + ijk]) + 
    504.*(g12[-2*di - 4*dk + ijk] + g12[-4*di - 2*dk + ijk] + 
       g12[4*di + 2*dk + ijk] + g12[2*di + 4*dk + ijk]) + 
    96.*(-g12[-3*di - 4*dk + ijk] + g12[3*di - 4*dk + ijk] - 
       g12[-4*di - 3*dk + ijk] + g12[4*di - 3*dk + ijk] + 
       g12[-4*di + 3*dk + ijk] - g12[4*di + 3*dk + ijk] + 
       g12[-3*di + 4*dk + ijk] - g12[3*di + 4*dk + ijk]) - 
    28224.*(g12[2*(di - dk) + ijk] + g12[2*(-di + dk) + ijk]) - 
    1024.*(g12[3*(di - dk) + ijk] + g12[3*(-di + dk) + ijk]) - 
    9.*(g12[4*(di - dk) + ijk] + g12[4*(-di + dk) + ijk]) + 
    28224.*(g12[-2*(di + dk) + ijk] + g12[2*(di + dk) + ijk]) + 
    1024.*(g12[-3*(di + dk) + ijk] + g12[3*(di + dk) + ijk]) + 
    9.*(g12[-4*(di + dk) + ijk] + g12[4*(di + dk) + ijk]))
;

deldelg1313
=
5.6689342403628117914e-6*oo2dx*oo2dz*
  (451584.*(g13[-di - dk + ijk] - g13[di - dk + ijk] - 
       g13[-di + dk + ijk] + g13[di + dk + ijk]) + 
    112896.*(-g13[-di - 2*dk + ijk] + g13[di - 2*dk + ijk] - 
       g13[-2*di - dk + ijk] + g13[2*di - dk + ijk] + 
       g13[-2*di + dk + ijk] - g13[2*di + dk + ijk] + 
       g13[-di + 2*dk + ijk] - g13[di + 2*dk + ijk]) - 
    21504.*(g13[di - 3*dk + ijk] + g13[3*di - dk + ijk] + 
       g13[-3*di + dk + ijk] + g13[-di + 3*dk + ijk]) + 
    21504.*(g13[-di - 3*dk + ijk] + g13[-3*di - dk + ijk] + 
       g13[3*di + dk + ijk] + g13[di + 3*dk + ijk]) + 
    5376.*(-g13[-2*di - 3*dk + ijk] + g13[2*di - 3*dk + ijk] - 
       g13[-3*di - 2*dk + ijk] + g13[3*di - 2*dk + ijk] + 
       g13[-3*di + 2*dk + ijk] - g13[3*di + 2*dk + ijk] + 
       g13[-2*di + 3*dk + ijk] - g13[2*di + 3*dk + ijk]) - 
    504.*(g13[2*di - 4*dk + ijk] + g13[4*di - 2*dk + ijk] + 
       g13[-4*di + 2*dk + ijk] + g13[-2*di + 4*dk + ijk]) + 
    2016.*(-g13[-di - 4*dk + ijk] + g13[di - 4*dk + ijk] - 
       g13[-4*di - dk + ijk] + g13[4*di - dk + ijk] + 
       g13[-4*di + dk + ijk] - g13[4*di + dk + ijk] + 
       g13[-di + 4*dk + ijk] - g13[di + 4*dk + ijk]) + 
    504.*(g13[-2*di - 4*dk + ijk] + g13[-4*di - 2*dk + ijk] + 
       g13[4*di + 2*dk + ijk] + g13[2*di + 4*dk + ijk]) + 
    96.*(-g13[-3*di - 4*dk + ijk] + g13[3*di - 4*dk + ijk] - 
       g13[-4*di - 3*dk + ijk] + g13[4*di - 3*dk + ijk] + 
       g13[-4*di + 3*dk + ijk] - g13[4*di + 3*dk + ijk] + 
       g13[-3*di + 4*dk + ijk] - g13[3*di + 4*dk + ijk]) - 
    28224.*(g13[2*(di - dk) + ijk] + g13[2*(-di + dk) + ijk]) - 
    1024.*(g13[3*(di - dk) + ijk] + g13[3*(-di + dk) + ijk]) - 
    9.*(g13[4*(di - dk) + ijk] + g13[4*(-di + dk) + ijk]) + 
    28224.*(g13[-2*(di + dk) + ijk] + g13[2*(di + dk) + ijk]) + 
    1024.*(g13[-3*(di + dk) + ijk] + g13[3*(di + dk) + ijk]) + 
    9.*(g13[-4*(di + dk) + ijk] + g13[4*(di + dk) + ijk]))
;

deldelg1322
=
5.6689342403628117914e-6*oo2dx*oo2dz*
  (451584.*(g22[-di - dk + ijk] - g22[di - dk + ijk] - 
       g22[-di + dk + ijk] + g22[di + dk + ijk]) + 
    112896.*(-g22[-di - 2*dk + ijk] + g22[di - 2*dk + ijk] - 
       g22[-2*di - dk + ijk] + g22[2*di - dk + ijk] + 
       g22[-2*di + dk + ijk] - g22[2*di + dk + ijk] + 
       g22[-di + 2*dk + ijk] - g22[di + 2*dk + ijk]) - 
    21504.*(g22[di - 3*dk + ijk] + g22[3*di - dk + ijk] + 
       g22[-3*di + dk + ijk] + g22[-di + 3*dk + ijk]) + 
    21504.*(g22[-di - 3*dk + ijk] + g22[-3*di - dk + ijk] + 
       g22[3*di + dk + ijk] + g22[di + 3*dk + ijk]) + 
    5376.*(-g22[-2*di - 3*dk + ijk] + g22[2*di - 3*dk + ijk] - 
       g22[-3*di - 2*dk + ijk] + g22[3*di - 2*dk + ijk] + 
       g22[-3*di + 2*dk + ijk] - g22[3*di + 2*dk + ijk] + 
       g22[-2*di + 3*dk + ijk] - g22[2*di + 3*dk + ijk]) - 
    504.*(g22[2*di - 4*dk + ijk] + g22[4*di - 2*dk + ijk] + 
       g22[-4*di + 2*dk + ijk] + g22[-2*di + 4*dk + ijk]) + 
    2016.*(-g22[-di - 4*dk + ijk] + g22[di - 4*dk + ijk] - 
       g22[-4*di - dk + ijk] + g22[4*di - dk + ijk] + 
       g22[-4*di + dk + ijk] - g22[4*di + dk + ijk] + 
       g22[-di + 4*dk + ijk] - g22[di + 4*dk + ijk]) + 
    504.*(g22[-2*di - 4*dk + ijk] + g22[-4*di - 2*dk + ijk] + 
       g22[4*di + 2*dk + ijk] + g22[2*di + 4*dk + ijk]) + 
    96.*(-g22[-3*di - 4*dk + ijk] + g22[3*di - 4*dk + ijk] - 
       g22[-4*di - 3*dk + ijk] + g22[4*di - 3*dk + ijk] + 
       g22[-4*di + 3*dk + ijk] - g22[4*di + 3*dk + ijk] + 
       g22[-3*di + 4*dk + ijk] - g22[3*di + 4*dk + ijk]) - 
    28224.*(g22[2*(di - dk) + ijk] + g22[2*(-di + dk) + ijk]) - 
    1024.*(g22[3*(di - dk) + ijk] + g22[3*(-di + dk) + ijk]) - 
    9.*(g22[4*(di - dk) + ijk] + g22[4*(-di + dk) + ijk]) + 
    28224.*(g22[-2*(di + dk) + ijk] + g22[2*(di + dk) + ijk]) + 
    1024.*(g22[-3*(di + dk) + ijk] + g22[3*(di + dk) + ijk]) + 
    9.*(g22[-4*(di + dk) + ijk] + g22[4*(di + dk) + ijk]))
;

deldelg1323
=
5.6689342403628117914e-6*oo2dx*oo2dz*
  (451584.*(g23[-di - dk + ijk] - g23[di - dk + ijk] - 
       g23[-di + dk + ijk] + g23[di + dk + ijk]) + 
    112896.*(-g23[-di - 2*dk + ijk] + g23[di - 2*dk + ijk] - 
       g23[-2*di - dk + ijk] + g23[2*di - dk + ijk] + 
       g23[-2*di + dk + ijk] - g23[2*di + dk + ijk] + 
       g23[-di + 2*dk + ijk] - g23[di + 2*dk + ijk]) - 
    21504.*(g23[di - 3*dk + ijk] + g23[3*di - dk + ijk] + 
       g23[-3*di + dk + ijk] + g23[-di + 3*dk + ijk]) + 
    21504.*(g23[-di - 3*dk + ijk] + g23[-3*di - dk + ijk] + 
       g23[3*di + dk + ijk] + g23[di + 3*dk + ijk]) + 
    5376.*(-g23[-2*di - 3*dk + ijk] + g23[2*di - 3*dk + ijk] - 
       g23[-3*di - 2*dk + ijk] + g23[3*di - 2*dk + ijk] + 
       g23[-3*di + 2*dk + ijk] - g23[3*di + 2*dk + ijk] + 
       g23[-2*di + 3*dk + ijk] - g23[2*di + 3*dk + ijk]) - 
    504.*(g23[2*di - 4*dk + ijk] + g23[4*di - 2*dk + ijk] + 
       g23[-4*di + 2*dk + ijk] + g23[-2*di + 4*dk + ijk]) + 
    2016.*(-g23[-di - 4*dk + ijk] + g23[di - 4*dk + ijk] - 
       g23[-4*di - dk + ijk] + g23[4*di - dk + ijk] + 
       g23[-4*di + dk + ijk] - g23[4*di + dk + ijk] + 
       g23[-di + 4*dk + ijk] - g23[di + 4*dk + ijk]) + 
    504.*(g23[-2*di - 4*dk + ijk] + g23[-4*di - 2*dk + ijk] + 
       g23[4*di + 2*dk + ijk] + g23[2*di + 4*dk + ijk]) + 
    96.*(-g23[-3*di - 4*dk + ijk] + g23[3*di - 4*dk + ijk] - 
       g23[-4*di - 3*dk + ijk] + g23[4*di - 3*dk + ijk] + 
       g23[-4*di + 3*dk + ijk] - g23[4*di + 3*dk + ijk] + 
       g23[-3*di + 4*dk + ijk] - g23[3*di + 4*dk + ijk]) - 
    28224.*(g23[2*(di - dk) + ijk] + g23[2*(-di + dk) + ijk]) - 
    1024.*(g23[3*(di - dk) + ijk] + g23[3*(-di + dk) + ijk]) - 
    9.*(g23[4*(di - dk) + ijk] + g23[4*(-di + dk) + ijk]) + 
    28224.*(g23[-2*(di + dk) + ijk] + g23[2*(di + dk) + ijk]) + 
    1024.*(g23[-3*(di + dk) + ijk] + g23[3*(di + dk) + ijk]) + 
    9.*(g23[-4*(di + dk) + ijk] + g23[4*(di + dk) + ijk]))
;

deldelg1333
=
5.6689342403628117914e-6*oo2dx*oo2dz*
  (451584.*(g33[-di - dk + ijk] - g33[di - dk + ijk] - 
       g33[-di + dk + ijk] + g33[di + dk + ijk]) + 
    112896.*(-g33[-di - 2*dk + ijk] + g33[di - 2*dk + ijk] - 
       g33[-2*di - dk + ijk] + g33[2*di - dk + ijk] + 
       g33[-2*di + dk + ijk] - g33[2*di + dk + ijk] + 
       g33[-di + 2*dk + ijk] - g33[di + 2*dk + ijk]) - 
    21504.*(g33[di - 3*dk + ijk] + g33[3*di - dk + ijk] + 
       g33[-3*di + dk + ijk] + g33[-di + 3*dk + ijk]) + 
    21504.*(g33[-di - 3*dk + ijk] + g33[-3*di - dk + ijk] + 
       g33[3*di + dk + ijk] + g33[di + 3*dk + ijk]) + 
    5376.*(-g33[-2*di - 3*dk + ijk] + g33[2*di - 3*dk + ijk] - 
       g33[-3*di - 2*dk + ijk] + g33[3*di - 2*dk + ijk] + 
       g33[-3*di + 2*dk + ijk] - g33[3*di + 2*dk + ijk] + 
       g33[-2*di + 3*dk + ijk] - g33[2*di + 3*dk + ijk]) - 
    504.*(g33[2*di - 4*dk + ijk] + g33[4*di - 2*dk + ijk] + 
       g33[-4*di + 2*dk + ijk] + g33[-2*di + 4*dk + ijk]) + 
    2016.*(-g33[-di - 4*dk + ijk] + g33[di - 4*dk + ijk] - 
       g33[-4*di - dk + ijk] + g33[4*di - dk + ijk] + 
       g33[-4*di + dk + ijk] - g33[4*di + dk + ijk] + 
       g33[-di + 4*dk + ijk] - g33[di + 4*dk + ijk]) + 
    504.*(g33[-2*di - 4*dk + ijk] + g33[-4*di - 2*dk + ijk] + 
       g33[4*di + 2*dk + ijk] + g33[2*di + 4*dk + ijk]) + 
    96.*(-g33[-3*di - 4*dk + ijk] + g33[3*di - 4*dk + ijk] - 
       g33[-4*di - 3*dk + ijk] + g33[4*di - 3*dk + ijk] + 
       g33[-4*di + 3*dk + ijk] - g33[4*di + 3*dk + ijk] + 
       g33[-3*di + 4*dk + ijk] - g33[3*di + 4*dk + ijk]) - 
    28224.*(g33[2*(di - dk) + ijk] + g33[2*(-di + dk) + ijk]) - 
    1024.*(g33[3*(di - dk) + ijk] + g33[3*(-di + dk) + ijk]) - 
    9.*(g33[4*(di - dk) + ijk] + g33[4*(-di + dk) + ijk]) + 
    28224.*(g33[-2*(di + dk) + ijk] + g33[2*(di + dk) + ijk]) + 
    1024.*(g33[-3*(di + dk) + ijk] + g33[3*(di + dk) + ijk]) + 
    9.*(g33[-4*(di + dk) + ijk] + g33[4*(di + dk) + ijk]))
;

deldelg2211
=
0.0001984126984126984127*oody2*
  (-14350.*g11[ijk] + 8064.*(g11[-dj + ijk] + g11[dj + ijk]) - 
    1008.*(g11[-2*dj + ijk] + g11[2*dj + ijk]) + 
    128.*(g11[-3*dj + ijk] + g11[3*dj + ijk]) - 
    9.*(g11[-4*dj + ijk] + g11[4*dj + ijk]))
;

deldelg2212
=
0.0001984126984126984127*oody2*
  (-14350.*g12[ijk] + 8064.*(g12[-dj + ijk] + g12[dj + ijk]) - 
    1008.*(g12[-2*dj + ijk] + g12[2*dj + ijk]) + 
    128.*(g12[-3*dj + ijk] + g12[3*dj + ijk]) - 
    9.*(g12[-4*dj + ijk] + g12[4*dj + ijk]))
;

deldelg2213
=
0.0001984126984126984127*oody2*
  (-14350.*g13[ijk] + 8064.*(g13[-dj + ijk] + g13[dj + ijk]) - 
    1008.*(g13[-2*dj + ijk] + g13[2*dj + ijk]) + 
    128.*(g13[-3*dj + ijk] + g13[3*dj + ijk]) - 
    9.*(g13[-4*dj + ijk] + g13[4*dj + ijk]))
;

deldelg2222
=
0.0001984126984126984127*oody2*
  (-14350.*g22[ijk] + 8064.*(g22[-dj + ijk] + g22[dj + ijk]) - 
    1008.*(g22[-2*dj + ijk] + g22[2*dj + ijk]) + 
    128.*(g22[-3*dj + ijk] + g22[3*dj + ijk]) - 
    9.*(g22[-4*dj + ijk] + g22[4*dj + ijk]))
;

deldelg2223
=
0.0001984126984126984127*oody2*
  (-14350.*g23[ijk] + 8064.*(g23[-dj + ijk] + g23[dj + ijk]) - 
    1008.*(g23[-2*dj + ijk] + g23[2*dj + ijk]) + 
    128.*(g23[-3*dj + ijk] + g23[3*dj + ijk]) - 
    9.*(g23[-4*dj + ijk] + g23[4*dj + ijk]))
;

deldelg2233
=
0.0001984126984126984127*oody2*
  (-14350.*g33[ijk] + 8064.*(g33[-dj + ijk] + g33[dj + ijk]) - 
    1008.*(g33[-2*dj + ijk] + g33[2*dj + ijk]) + 
    128.*(g33[-3*dj + ijk] + g33[3*dj + ijk]) - 
    9.*(g33[-4*dj + ijk] + g33[4*dj + ijk]))
;

deldelg2311
=
5.6689342403628117914e-6*oo2dy*oo2dz*
  (451584.*(g11[-dj - dk + ijk] - g11[dj - dk + ijk] - 
       g11[-dj + dk + ijk] + g11[dj + dk + ijk]) + 
    112896.*(-g11[-dj - 2*dk + ijk] + g11[dj - 2*dk + ijk] - 
       g11[-2*dj - dk + ijk] + g11[2*dj - dk + ijk] + 
       g11[-2*dj + dk + ijk] - g11[2*dj + dk + ijk] + 
       g11[-dj + 2*dk + ijk] - g11[dj + 2*dk + ijk]) - 
    21504.*(g11[dj - 3*dk + ijk] + g11[3*dj - dk + ijk] + 
       g11[-3*dj + dk + ijk] + g11[-dj + 3*dk + ijk]) + 
    21504.*(g11[-dj - 3*dk + ijk] + g11[-3*dj - dk + ijk] + 
       g11[3*dj + dk + ijk] + g11[dj + 3*dk + ijk]) + 
    5376.*(-g11[-2*dj - 3*dk + ijk] + g11[2*dj - 3*dk + ijk] - 
       g11[-3*dj - 2*dk + ijk] + g11[3*dj - 2*dk + ijk] + 
       g11[-3*dj + 2*dk + ijk] - g11[3*dj + 2*dk + ijk] + 
       g11[-2*dj + 3*dk + ijk] - g11[2*dj + 3*dk + ijk]) - 
    504.*(g11[2*dj - 4*dk + ijk] + g11[4*dj - 2*dk + ijk] + 
       g11[-4*dj + 2*dk + ijk] + g11[-2*dj + 4*dk + ijk]) + 
    2016.*(-g11[-dj - 4*dk + ijk] + g11[dj - 4*dk + ijk] - 
       g11[-4*dj - dk + ijk] + g11[4*dj - dk + ijk] + 
       g11[-4*dj + dk + ijk] - g11[4*dj + dk + ijk] + 
       g11[-dj + 4*dk + ijk] - g11[dj + 4*dk + ijk]) + 
    504.*(g11[-2*dj - 4*dk + ijk] + g11[-4*dj - 2*dk + ijk] + 
       g11[4*dj + 2*dk + ijk] + g11[2*dj + 4*dk + ijk]) + 
    96.*(-g11[-3*dj - 4*dk + ijk] + g11[3*dj - 4*dk + ijk] - 
       g11[-4*dj - 3*dk + ijk] + g11[4*dj - 3*dk + ijk] + 
       g11[-4*dj + 3*dk + ijk] - g11[4*dj + 3*dk + ijk] + 
       g11[-3*dj + 4*dk + ijk] - g11[3*dj + 4*dk + ijk]) - 
    28224.*(g11[2*(dj - dk) + ijk] + g11[2*(-dj + dk) + ijk]) - 
    1024.*(g11[3*(dj - dk) + ijk] + g11[3*(-dj + dk) + ijk]) - 
    9.*(g11[4*(dj - dk) + ijk] + g11[4*(-dj + dk) + ijk]) + 
    28224.*(g11[-2*(dj + dk) + ijk] + g11[2*(dj + dk) + ijk]) + 
    1024.*(g11[-3*(dj + dk) + ijk] + g11[3*(dj + dk) + ijk]) + 
    9.*(g11[-4*(dj + dk) + ijk] + g11[4*(dj + dk) + ijk]))
;

deldelg2312
=
5.6689342403628117914e-6*oo2dy*oo2dz*
  (451584.*(g12[-dj - dk + ijk] - g12[dj - dk + ijk] - 
       g12[-dj + dk + ijk] + g12[dj + dk + ijk]) + 
    112896.*(-g12[-dj - 2*dk + ijk] + g12[dj - 2*dk + ijk] - 
       g12[-2*dj - dk + ijk] + g12[2*dj - dk + ijk] + 
       g12[-2*dj + dk + ijk] - g12[2*dj + dk + ijk] + 
       g12[-dj + 2*dk + ijk] - g12[dj + 2*dk + ijk]) - 
    21504.*(g12[dj - 3*dk + ijk] + g12[3*dj - dk + ijk] + 
       g12[-3*dj + dk + ijk] + g12[-dj + 3*dk + ijk]) + 
    21504.*(g12[-dj - 3*dk + ijk] + g12[-3*dj - dk + ijk] + 
       g12[3*dj + dk + ijk] + g12[dj + 3*dk + ijk]) + 
    5376.*(-g12[-2*dj - 3*dk + ijk] + g12[2*dj - 3*dk + ijk] - 
       g12[-3*dj - 2*dk + ijk] + g12[3*dj - 2*dk + ijk] + 
       g12[-3*dj + 2*dk + ijk] - g12[3*dj + 2*dk + ijk] + 
       g12[-2*dj + 3*dk + ijk] - g12[2*dj + 3*dk + ijk]) - 
    504.*(g12[2*dj - 4*dk + ijk] + g12[4*dj - 2*dk + ijk] + 
       g12[-4*dj + 2*dk + ijk] + g12[-2*dj + 4*dk + ijk]) + 
    2016.*(-g12[-dj - 4*dk + ijk] + g12[dj - 4*dk + ijk] - 
       g12[-4*dj - dk + ijk] + g12[4*dj - dk + ijk] + 
       g12[-4*dj + dk + ijk] - g12[4*dj + dk + ijk] + 
       g12[-dj + 4*dk + ijk] - g12[dj + 4*dk + ijk]) + 
    504.*(g12[-2*dj - 4*dk + ijk] + g12[-4*dj - 2*dk + ijk] + 
       g12[4*dj + 2*dk + ijk] + g12[2*dj + 4*dk + ijk]) + 
    96.*(-g12[-3*dj - 4*dk + ijk] + g12[3*dj - 4*dk + ijk] - 
       g12[-4*dj - 3*dk + ijk] + g12[4*dj - 3*dk + ijk] + 
       g12[-4*dj + 3*dk + ijk] - g12[4*dj + 3*dk + ijk] + 
       g12[-3*dj + 4*dk + ijk] - g12[3*dj + 4*dk + ijk]) - 
    28224.*(g12[2*(dj - dk) + ijk] + g12[2*(-dj + dk) + ijk]) - 
    1024.*(g12[3*(dj - dk) + ijk] + g12[3*(-dj + dk) + ijk]) - 
    9.*(g12[4*(dj - dk) + ijk] + g12[4*(-dj + dk) + ijk]) + 
    28224.*(g12[-2*(dj + dk) + ijk] + g12[2*(dj + dk) + ijk]) + 
    1024.*(g12[-3*(dj + dk) + ijk] + g12[3*(dj + dk) + ijk]) + 
    9.*(g12[-4*(dj + dk) + ijk] + g12[4*(dj + dk) + ijk]))
;

deldelg2313
=
5.6689342403628117914e-6*oo2dy*oo2dz*
  (451584.*(g13[-dj - dk + ijk] - g13[dj - dk + ijk] - 
       g13[-dj + dk + ijk] + g13[dj + dk + ijk]) + 
    112896.*(-g13[-dj - 2*dk + ijk] + g13[dj - 2*dk + ijk] - 
       g13[-2*dj - dk + ijk] + g13[2*dj - dk + ijk] + 
       g13[-2*dj + dk + ijk] - g13[2*dj + dk + ijk] + 
       g13[-dj + 2*dk + ijk] - g13[dj + 2*dk + ijk]) - 
    21504.*(g13[dj - 3*dk + ijk] + g13[3*dj - dk + ijk] + 
       g13[-3*dj + dk + ijk] + g13[-dj + 3*dk + ijk]) + 
    21504.*(g13[-dj - 3*dk + ijk] + g13[-3*dj - dk + ijk] + 
       g13[3*dj + dk + ijk] + g13[dj + 3*dk + ijk]) + 
    5376.*(-g13[-2*dj - 3*dk + ijk] + g13[2*dj - 3*dk + ijk] - 
       g13[-3*dj - 2*dk + ijk] + g13[3*dj - 2*dk + ijk] + 
       g13[-3*dj + 2*dk + ijk] - g13[3*dj + 2*dk + ijk] + 
       g13[-2*dj + 3*dk + ijk] - g13[2*dj + 3*dk + ijk]) - 
    504.*(g13[2*dj - 4*dk + ijk] + g13[4*dj - 2*dk + ijk] + 
       g13[-4*dj + 2*dk + ijk] + g13[-2*dj + 4*dk + ijk]) + 
    2016.*(-g13[-dj - 4*dk + ijk] + g13[dj - 4*dk + ijk] - 
       g13[-4*dj - dk + ijk] + g13[4*dj - dk + ijk] + 
       g13[-4*dj + dk + ijk] - g13[4*dj + dk + ijk] + 
       g13[-dj + 4*dk + ijk] - g13[dj + 4*dk + ijk]) + 
    504.*(g13[-2*dj - 4*dk + ijk] + g13[-4*dj - 2*dk + ijk] + 
       g13[4*dj + 2*dk + ijk] + g13[2*dj + 4*dk + ijk]) + 
    96.*(-g13[-3*dj - 4*dk + ijk] + g13[3*dj - 4*dk + ijk] - 
       g13[-4*dj - 3*dk + ijk] + g13[4*dj - 3*dk + ijk] + 
       g13[-4*dj + 3*dk + ijk] - g13[4*dj + 3*dk + ijk] + 
       g13[-3*dj + 4*dk + ijk] - g13[3*dj + 4*dk + ijk]) - 
    28224.*(g13[2*(dj - dk) + ijk] + g13[2*(-dj + dk) + ijk]) - 
    1024.*(g13[3*(dj - dk) + ijk] + g13[3*(-dj + dk) + ijk]) - 
    9.*(g13[4*(dj - dk) + ijk] + g13[4*(-dj + dk) + ijk]) + 
    28224.*(g13[-2*(dj + dk) + ijk] + g13[2*(dj + dk) + ijk]) + 
    1024.*(g13[-3*(dj + dk) + ijk] + g13[3*(dj + dk) + ijk]) + 
    9.*(g13[-4*(dj + dk) + ijk] + g13[4*(dj + dk) + ijk]))
;

deldelg2322
=
5.6689342403628117914e-6*oo2dy*oo2dz*
  (451584.*(g22[-dj - dk + ijk] - g22[dj - dk + ijk] - 
       g22[-dj + dk + ijk] + g22[dj + dk + ijk]) + 
    112896.*(-g22[-dj - 2*dk + ijk] + g22[dj - 2*dk + ijk] - 
       g22[-2*dj - dk + ijk] + g22[2*dj - dk + ijk] + 
       g22[-2*dj + dk + ijk] - g22[2*dj + dk + ijk] + 
       g22[-dj + 2*dk + ijk] - g22[dj + 2*dk + ijk]) - 
    21504.*(g22[dj - 3*dk + ijk] + g22[3*dj - dk + ijk] + 
       g22[-3*dj + dk + ijk] + g22[-dj + 3*dk + ijk]) + 
    21504.*(g22[-dj - 3*dk + ijk] + g22[-3*dj - dk + ijk] + 
       g22[3*dj + dk + ijk] + g22[dj + 3*dk + ijk]) + 
    5376.*(-g22[-2*dj - 3*dk + ijk] + g22[2*dj - 3*dk + ijk] - 
       g22[-3*dj - 2*dk + ijk] + g22[3*dj - 2*dk + ijk] + 
       g22[-3*dj + 2*dk + ijk] - g22[3*dj + 2*dk + ijk] + 
       g22[-2*dj + 3*dk + ijk] - g22[2*dj + 3*dk + ijk]) - 
    504.*(g22[2*dj - 4*dk + ijk] + g22[4*dj - 2*dk + ijk] + 
       g22[-4*dj + 2*dk + ijk] + g22[-2*dj + 4*dk + ijk]) + 
    2016.*(-g22[-dj - 4*dk + ijk] + g22[dj - 4*dk + ijk] - 
       g22[-4*dj - dk + ijk] + g22[4*dj - dk + ijk] + 
       g22[-4*dj + dk + ijk] - g22[4*dj + dk + ijk] + 
       g22[-dj + 4*dk + ijk] - g22[dj + 4*dk + ijk]) + 
    504.*(g22[-2*dj - 4*dk + ijk] + g22[-4*dj - 2*dk + ijk] + 
       g22[4*dj + 2*dk + ijk] + g22[2*dj + 4*dk + ijk]) + 
    96.*(-g22[-3*dj - 4*dk + ijk] + g22[3*dj - 4*dk + ijk] - 
       g22[-4*dj - 3*dk + ijk] + g22[4*dj - 3*dk + ijk] + 
       g22[-4*dj + 3*dk + ijk] - g22[4*dj + 3*dk + ijk] + 
       g22[-3*dj + 4*dk + ijk] - g22[3*dj + 4*dk + ijk]) - 
    28224.*(g22[2*(dj - dk) + ijk] + g22[2*(-dj + dk) + ijk]) - 
    1024.*(g22[3*(dj - dk) + ijk] + g22[3*(-dj + dk) + ijk]) - 
    9.*(g22[4*(dj - dk) + ijk] + g22[4*(-dj + dk) + ijk]) + 
    28224.*(g22[-2*(dj + dk) + ijk] + g22[2*(dj + dk) + ijk]) + 
    1024.*(g22[-3*(dj + dk) + ijk] + g22[3*(dj + dk) + ijk]) + 
    9.*(g22[-4*(dj + dk) + ijk] + g22[4*(dj + dk) + ijk]))
;

deldelg2323
=
5.6689342403628117914e-6*oo2dy*oo2dz*
  (451584.*(g23[-dj - dk + ijk] - g23[dj - dk + ijk] - 
       g23[-dj + dk + ijk] + g23[dj + dk + ijk]) + 
    112896.*(-g23[-dj - 2*dk + ijk] + g23[dj - 2*dk + ijk] - 
       g23[-2*dj - dk + ijk] + g23[2*dj - dk + ijk] + 
       g23[-2*dj + dk + ijk] - g23[2*dj + dk + ijk] + 
       g23[-dj + 2*dk + ijk] - g23[dj + 2*dk + ijk]) - 
    21504.*(g23[dj - 3*dk + ijk] + g23[3*dj - dk + ijk] + 
       g23[-3*dj + dk + ijk] + g23[-dj + 3*dk + ijk]) + 
    21504.*(g23[-dj - 3*dk + ijk] + g23[-3*dj - dk + ijk] + 
       g23[3*dj + dk + ijk] + g23[dj + 3*dk + ijk]) + 
    5376.*(-g23[-2*dj - 3*dk + ijk] + g23[2*dj - 3*dk + ijk] - 
       g23[-3*dj - 2*dk + ijk] + g23[3*dj - 2*dk + ijk] + 
       g23[-3*dj + 2*dk + ijk] - g23[3*dj + 2*dk + ijk] + 
       g23[-2*dj + 3*dk + ijk] - g23[2*dj + 3*dk + ijk]) - 
    504.*(g23[2*dj - 4*dk + ijk] + g23[4*dj - 2*dk + ijk] + 
       g23[-4*dj + 2*dk + ijk] + g23[-2*dj + 4*dk + ijk]) + 
    2016.*(-g23[-dj - 4*dk + ijk] + g23[dj - 4*dk + ijk] - 
       g23[-4*dj - dk + ijk] + g23[4*dj - dk + ijk] + 
       g23[-4*dj + dk + ijk] - g23[4*dj + dk + ijk] + 
       g23[-dj + 4*dk + ijk] - g23[dj + 4*dk + ijk]) + 
    504.*(g23[-2*dj - 4*dk + ijk] + g23[-4*dj - 2*dk + ijk] + 
       g23[4*dj + 2*dk + ijk] + g23[2*dj + 4*dk + ijk]) + 
    96.*(-g23[-3*dj - 4*dk + ijk] + g23[3*dj - 4*dk + ijk] - 
       g23[-4*dj - 3*dk + ijk] + g23[4*dj - 3*dk + ijk] + 
       g23[-4*dj + 3*dk + ijk] - g23[4*dj + 3*dk + ijk] + 
       g23[-3*dj + 4*dk + ijk] - g23[3*dj + 4*dk + ijk]) - 
    28224.*(g23[2*(dj - dk) + ijk] + g23[2*(-dj + dk) + ijk]) - 
    1024.*(g23[3*(dj - dk) + ijk] + g23[3*(-dj + dk) + ijk]) - 
    9.*(g23[4*(dj - dk) + ijk] + g23[4*(-dj + dk) + ijk]) + 
    28224.*(g23[-2*(dj + dk) + ijk] + g23[2*(dj + dk) + ijk]) + 
    1024.*(g23[-3*(dj + dk) + ijk] + g23[3*(dj + dk) + ijk]) + 
    9.*(g23[-4*(dj + dk) + ijk] + g23[4*(dj + dk) + ijk]))
;

deldelg2333
=
5.6689342403628117914e-6*oo2dy*oo2dz*
  (451584.*(g33[-dj - dk + ijk] - g33[dj - dk + ijk] - 
       g33[-dj + dk + ijk] + g33[dj + dk + ijk]) + 
    112896.*(-g33[-dj - 2*dk + ijk] + g33[dj - 2*dk + ijk] - 
       g33[-2*dj - dk + ijk] + g33[2*dj - dk + ijk] + 
       g33[-2*dj + dk + ijk] - g33[2*dj + dk + ijk] + 
       g33[-dj + 2*dk + ijk] - g33[dj + 2*dk + ijk]) - 
    21504.*(g33[dj - 3*dk + ijk] + g33[3*dj - dk + ijk] + 
       g33[-3*dj + dk + ijk] + g33[-dj + 3*dk + ijk]) + 
    21504.*(g33[-dj - 3*dk + ijk] + g33[-3*dj - dk + ijk] + 
       g33[3*dj + dk + ijk] + g33[dj + 3*dk + ijk]) + 
    5376.*(-g33[-2*dj - 3*dk + ijk] + g33[2*dj - 3*dk + ijk] - 
       g33[-3*dj - 2*dk + ijk] + g33[3*dj - 2*dk + ijk] + 
       g33[-3*dj + 2*dk + ijk] - g33[3*dj + 2*dk + ijk] + 
       g33[-2*dj + 3*dk + ijk] - g33[2*dj + 3*dk + ijk]) - 
    504.*(g33[2*dj - 4*dk + ijk] + g33[4*dj - 2*dk + ijk] + 
       g33[-4*dj + 2*dk + ijk] + g33[-2*dj + 4*dk + ijk]) + 
    2016.*(-g33[-dj - 4*dk + ijk] + g33[dj - 4*dk + ijk] - 
       g33[-4*dj - dk + ijk] + g33[4*dj - dk + ijk] + 
       g33[-4*dj + dk + ijk] - g33[4*dj + dk + ijk] + 
       g33[-dj + 4*dk + ijk] - g33[dj + 4*dk + ijk]) + 
    504.*(g33[-2*dj - 4*dk + ijk] + g33[-4*dj - 2*dk + ijk] + 
       g33[4*dj + 2*dk + ijk] + g33[2*dj + 4*dk + ijk]) + 
    96.*(-g33[-3*dj - 4*dk + ijk] + g33[3*dj - 4*dk + ijk] - 
       g33[-4*dj - 3*dk + ijk] + g33[4*dj - 3*dk + ijk] + 
       g33[-4*dj + 3*dk + ijk] - g33[4*dj + 3*dk + ijk] + 
       g33[-3*dj + 4*dk + ijk] - g33[3*dj + 4*dk + ijk]) - 
    28224.*(g33[2*(dj - dk) + ijk] + g33[2*(-dj + dk) + ijk]) - 
    1024.*(g33[3*(dj - dk) + ijk] + g33[3*(-dj + dk) + ijk]) - 
    9.*(g33[4*(dj - dk) + ijk] + g33[4*(-dj + dk) + ijk]) + 
    28224.*(g33[-2*(dj + dk) + ijk] + g33[2*(dj + dk) + ijk]) + 
    1024.*(g33[-3*(dj + dk) + ijk] + g33[3*(dj + dk) + ijk]) + 
    9.*(g33[-4*(dj + dk) + ijk] + g33[4*(dj + dk) + ijk]))
;

deldelg3311
=
0.0001984126984126984127*oodz2*
  (-14350.*g11[ijk] + 8064.*(g11[-dk + ijk] + g11[dk + ijk]) - 
    1008.*(g11[-2*dk + ijk] + g11[2*dk + ijk]) + 
    128.*(g11[-3*dk + ijk] + g11[3*dk + ijk]) - 
    9.*(g11[-4*dk + ijk] + g11[4*dk + ijk]))
;

deldelg3312
=
0.0001984126984126984127*oodz2*
  (-14350.*g12[ijk] + 8064.*(g12[-dk + ijk] + g12[dk + ijk]) - 
    1008.*(g12[-2*dk + ijk] + g12[2*dk + ijk]) + 
    128.*(g12[-3*dk + ijk] + g12[3*dk + ijk]) - 
    9.*(g12[-4*dk + ijk] + g12[4*dk + ijk]))
;

deldelg3313
=
0.0001984126984126984127*oodz2*
  (-14350.*g13[ijk] + 8064.*(g13[-dk + ijk] + g13[dk + ijk]) - 
    1008.*(g13[-2*dk + ijk] + g13[2*dk + ijk]) + 
    128.*(g13[-3*dk + ijk] + g13[3*dk + ijk]) - 
    9.*(g13[-4*dk + ijk] + g13[4*dk + ijk]))
;

deldelg3322
=
0.0001984126984126984127*oodz2*
  (-14350.*g22[ijk] + 8064.*(g22[-dk + ijk] + g22[dk + ijk]) - 
    1008.*(g22[-2*dk + ijk] + g22[2*dk + ijk]) + 
    128.*(g22[-3*dk + ijk] + g22[3*dk + ijk]) - 
    9.*(g22[-4*dk + ijk] + g22[4*dk + ijk]))
;

deldelg3323
=
0.0001984126984126984127*oodz2*
  (-14350.*g23[ijk] + 8064.*(g23[-dk + ijk] + g23[dk + ijk]) - 
    1008.*(g23[-2*dk + ijk] + g23[2*dk + ijk]) + 
    128.*(g23[-3*dk + ijk] + g23[3*dk + ijk]) - 
    9.*(g23[-4*dk + ijk] + g23[4*dk + ijk]))
;

deldelg3333
=
0.0001984126984126984127*oodz2*
  (-14350.*g33[ijk] + 8064.*(g33[-dk + ijk] + g33[dk + ijk]) - 
    1008.*(g33[-2*dk + ijk] + g33[2*dk + ijk]) + 
    128.*(g33[-3*dk + ijk] + g33[3*dk + ijk]) - 
    9.*(g33[-4*dk + ijk] + g33[4*dk + ijk]))
;

delG11
=
0.0023809523809523809524*oo2dx*
  (3.*G1[-4*di + ijk] + 168.*G1[-2*di + ijk] + 
    672.*(-G1[-di + ijk] + G1[di + ijk]) - 168.*G1[2*di + ijk] + 
    32.*(-G1[-3*di + ijk] + G1[3*di + ijk]) - 3.*G1[4*di + ijk])
;

delG12
=
0.0023809523809523809524*oo2dx*
  (3.*G2[-4*di + ijk] + 168.*G2[-2*di + ijk] + 
    672.*(-G2[-di + ijk] + G2[di + ijk]) - 168.*G2[2*di + ijk] + 
    32.*(-G2[-3*di + ijk] + G2[3*di + ijk]) - 3.*G2[4*di + ijk])
;

delG13
=
0.0023809523809523809524*oo2dx*
  (3.*G3[-4*di + ijk] + 168.*G3[-2*di + ijk] + 
    672.*(-G3[-di + ijk] + G3[di + ijk]) - 168.*G3[2*di + ijk] + 
    32.*(-G3[-3*di + ijk] + G3[3*di + ijk]) - 3.*G3[4*di + ijk])
;

delG21
=
0.0023809523809523809524*oo2dy*
  (3.*G1[-4*dj + ijk] + 168.*G1[-2*dj + ijk] + 
    672.*(-G1[-dj + ijk] + G1[dj + ijk]) - 168.*G1[2*dj + ijk] + 
    32.*(-G1[-3*dj + ijk] + G1[3*dj + ijk]) - 3.*G1[4*dj + ijk])
;

delG22
=
0.0023809523809523809524*oo2dy*
  (3.*G2[-4*dj + ijk] + 168.*G2[-2*dj + ijk] + 
    672.*(-G2[-dj + ijk] + G2[dj + ijk]) - 168.*G2[2*dj + ijk] + 
    32.*(-G2[-3*dj + ijk] + G2[3*dj + ijk]) - 3.*G2[4*dj + ijk])
;

delG23
=
0.0023809523809523809524*oo2dy*
  (3.*G3[-4*dj + ijk] + 168.*G3[-2*dj + ijk] + 
    672.*(-G3[-dj + ijk] + G3[dj + ijk]) - 168.*G3[2*dj + ijk] + 
    32.*(-G3[-3*dj + ijk] + G3[3*dj + ijk]) - 3.*G3[4*dj + ijk])
;

delG31
=
0.0023809523809523809524*oo2dz*
  (3.*G1[-4*dk + ijk] + 168.*G1[-2*dk + ijk] + 
    672.*(-G1[-dk + ijk] + G1[dk + ijk]) - 168.*G1[2*dk + ijk] + 
    32.*(-G1[-3*dk + ijk] + G1[3*dk + ijk]) - 3.*G1[4*dk + ijk])
;

delG32
=
0.0023809523809523809524*oo2dz*
  (3.*G2[-4*dk + ijk] + 168.*G2[-2*dk + ijk] + 
    672.*(-G2[-dk + ijk] + G2[dk + ijk]) - 168.*G2[2*dk + ijk] + 
    32.*(-G2[-3*dk + ijk] + G2[3*dk + ijk]) - 3.*G2[4*dk + ijk])
;

delG33
=
0.0023809523809523809524*oo2dz*
  (3.*G3[-4*dk + ijk] + 168.*G3[-2*dk + ijk] + 
    672.*(-G3[-dk + ijk] + G3[dk + ijk]) - 168.*G3[2*dk + ijk] + 
    32.*(-G3[-3*dk + ijk] + G3[3*dk + ijk]) - 3.*G3[4*dk + ijk])
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

dchi1
=
0.0023809523809523809524*oo2dx*
  (3.*chi[-4*di + ijk] + 168.*chi[-2*di + ijk] + 
    672.*(-chi[-di + ijk] + chi[di + ijk]) - 168.*chi[2*di + ijk] + 
    32.*(-chi[-3*di + ijk] + chi[3*di + ijk]) - 3.*chi[4*di + ijk])
;

dchi2
=
0.0023809523809523809524*oo2dy*
  (3.*chi[-4*dj + ijk] + 168.*chi[-2*dj + ijk] + 
    672.*(-chi[-dj + ijk] + chi[dj + ijk]) - 168.*chi[2*dj + ijk] + 
    32.*(-chi[-3*dj + ijk] + chi[3*dj + ijk]) - 3.*chi[4*dj + ijk])
;

dchi3
=
0.0023809523809523809524*oo2dz*
  (3.*chi[-4*dk + ijk] + 168.*chi[-2*dk + ijk] + 
    672.*(-chi[-dk + ijk] + chi[dk + ijk]) - 168.*chi[2*dk + ijk] + 
    32.*(-chi[-3*dk + ijk] + chi[3*dk + ijk]) - 3.*chi[4*dk + ijk])
;

ddchi11
=
0.0001984126984126984127*oodx2*
  (-14350.*chi[ijk] + 8064.*(chi[-di + ijk] + chi[di + ijk]) - 
    1008.*(chi[-2*di + ijk] + chi[2*di + ijk]) + 
    128.*(chi[-3*di + ijk] + chi[3*di + ijk]) - 
    9.*(chi[-4*di + ijk] + chi[4*di + ijk]))
;

ddchi12
=
5.6689342403628117914e-6*oo2dx*oo2dy*
  (451584.*(chi[-di - dj + ijk] - chi[di - dj + ijk] - 
       chi[-di + dj + ijk] + chi[di + dj + ijk]) + 
    112896.*(-chi[-di - 2*dj + ijk] + chi[di - 2*dj + ijk] - 
       chi[-2*di - dj + ijk] + chi[2*di - dj + ijk] + 
       chi[-2*di + dj + ijk] - chi[2*di + dj + ijk] + 
       chi[-di + 2*dj + ijk] - chi[di + 2*dj + ijk]) - 
    21504.*(chi[di - 3*dj + ijk] + chi[3*di - dj + ijk] + 
       chi[-3*di + dj + ijk] + chi[-di + 3*dj + ijk]) + 
    21504.*(chi[-di - 3*dj + ijk] + chi[-3*di - dj + ijk] + 
       chi[3*di + dj + ijk] + chi[di + 3*dj + ijk]) + 
    5376.*(-chi[-2*di - 3*dj + ijk] + chi[2*di - 3*dj + ijk] - 
       chi[-3*di - 2*dj + ijk] + chi[3*di - 2*dj + ijk] + 
       chi[-3*di + 2*dj + ijk] - chi[3*di + 2*dj + ijk] + 
       chi[-2*di + 3*dj + ijk] - chi[2*di + 3*dj + ijk]) - 
    504.*(chi[2*di - 4*dj + ijk] + chi[4*di - 2*dj + ijk] + 
       chi[-4*di + 2*dj + ijk] + chi[-2*di + 4*dj + ijk]) + 
    2016.*(-chi[-di - 4*dj + ijk] + chi[di - 4*dj + ijk] - 
       chi[-4*di - dj + ijk] + chi[4*di - dj + ijk] + 
       chi[-4*di + dj + ijk] - chi[4*di + dj + ijk] + 
       chi[-di + 4*dj + ijk] - chi[di + 4*dj + ijk]) + 
    504.*(chi[-2*di - 4*dj + ijk] + chi[-4*di - 2*dj + ijk] + 
       chi[4*di + 2*dj + ijk] + chi[2*di + 4*dj + ijk]) + 
    96.*(-chi[-3*di - 4*dj + ijk] + chi[3*di - 4*dj + ijk] - 
       chi[-4*di - 3*dj + ijk] + chi[4*di - 3*dj + ijk] + 
       chi[-4*di + 3*dj + ijk] - chi[4*di + 3*dj + ijk] + 
       chi[-3*di + 4*dj + ijk] - chi[3*di + 4*dj + ijk]) - 
    28224.*(chi[2*(di - dj) + ijk] + chi[2*(-di + dj) + ijk]) - 
    1024.*(chi[3*(di - dj) + ijk] + chi[3*(-di + dj) + ijk]) - 
    9.*(chi[4*(di - dj) + ijk] + chi[4*(-di + dj) + ijk]) + 
    28224.*(chi[-2*(di + dj) + ijk] + chi[2*(di + dj) + ijk]) + 
    1024.*(chi[-3*(di + dj) + ijk] + chi[3*(di + dj) + ijk]) + 
    9.*(chi[-4*(di + dj) + ijk] + chi[4*(di + dj) + ijk]))
;

ddchi13
=
5.6689342403628117914e-6*oo2dx*oo2dz*
  (451584.*(chi[-di - dk + ijk] - chi[di - dk + ijk] - 
       chi[-di + dk + ijk] + chi[di + dk + ijk]) + 
    112896.*(-chi[-di - 2*dk + ijk] + chi[di - 2*dk + ijk] - 
       chi[-2*di - dk + ijk] + chi[2*di - dk + ijk] + 
       chi[-2*di + dk + ijk] - chi[2*di + dk + ijk] + 
       chi[-di + 2*dk + ijk] - chi[di + 2*dk + ijk]) - 
    21504.*(chi[di - 3*dk + ijk] + chi[3*di - dk + ijk] + 
       chi[-3*di + dk + ijk] + chi[-di + 3*dk + ijk]) + 
    21504.*(chi[-di - 3*dk + ijk] + chi[-3*di - dk + ijk] + 
       chi[3*di + dk + ijk] + chi[di + 3*dk + ijk]) + 
    5376.*(-chi[-2*di - 3*dk + ijk] + chi[2*di - 3*dk + ijk] - 
       chi[-3*di - 2*dk + ijk] + chi[3*di - 2*dk + ijk] + 
       chi[-3*di + 2*dk + ijk] - chi[3*di + 2*dk + ijk] + 
       chi[-2*di + 3*dk + ijk] - chi[2*di + 3*dk + ijk]) - 
    504.*(chi[2*di - 4*dk + ijk] + chi[4*di - 2*dk + ijk] + 
       chi[-4*di + 2*dk + ijk] + chi[-2*di + 4*dk + ijk]) + 
    2016.*(-chi[-di - 4*dk + ijk] + chi[di - 4*dk + ijk] - 
       chi[-4*di - dk + ijk] + chi[4*di - dk + ijk] + 
       chi[-4*di + dk + ijk] - chi[4*di + dk + ijk] + 
       chi[-di + 4*dk + ijk] - chi[di + 4*dk + ijk]) + 
    504.*(chi[-2*di - 4*dk + ijk] + chi[-4*di - 2*dk + ijk] + 
       chi[4*di + 2*dk + ijk] + chi[2*di + 4*dk + ijk]) + 
    96.*(-chi[-3*di - 4*dk + ijk] + chi[3*di - 4*dk + ijk] - 
       chi[-4*di - 3*dk + ijk] + chi[4*di - 3*dk + ijk] + 
       chi[-4*di + 3*dk + ijk] - chi[4*di + 3*dk + ijk] + 
       chi[-3*di + 4*dk + ijk] - chi[3*di + 4*dk + ijk]) - 
    28224.*(chi[2*(di - dk) + ijk] + chi[2*(-di + dk) + ijk]) - 
    1024.*(chi[3*(di - dk) + ijk] + chi[3*(-di + dk) + ijk]) - 
    9.*(chi[4*(di - dk) + ijk] + chi[4*(-di + dk) + ijk]) + 
    28224.*(chi[-2*(di + dk) + ijk] + chi[2*(di + dk) + ijk]) + 
    1024.*(chi[-3*(di + dk) + ijk] + chi[3*(di + dk) + ijk]) + 
    9.*(chi[-4*(di + dk) + ijk] + chi[4*(di + dk) + ijk]))
;

ddchi22
=
0.0001984126984126984127*oody2*
  (-14350.*chi[ijk] + 8064.*(chi[-dj + ijk] + chi[dj + ijk]) - 
    1008.*(chi[-2*dj + ijk] + chi[2*dj + ijk]) + 
    128.*(chi[-3*dj + ijk] + chi[3*dj + ijk]) - 
    9.*(chi[-4*dj + ijk] + chi[4*dj + ijk]))
;

ddchi23
=
5.6689342403628117914e-6*oo2dy*oo2dz*
  (451584.*(chi[-dj - dk + ijk] - chi[dj - dk + ijk] - 
       chi[-dj + dk + ijk] + chi[dj + dk + ijk]) + 
    112896.*(-chi[-dj - 2*dk + ijk] + chi[dj - 2*dk + ijk] - 
       chi[-2*dj - dk + ijk] + chi[2*dj - dk + ijk] + 
       chi[-2*dj + dk + ijk] - chi[2*dj + dk + ijk] + 
       chi[-dj + 2*dk + ijk] - chi[dj + 2*dk + ijk]) - 
    21504.*(chi[dj - 3*dk + ijk] + chi[3*dj - dk + ijk] + 
       chi[-3*dj + dk + ijk] + chi[-dj + 3*dk + ijk]) + 
    21504.*(chi[-dj - 3*dk + ijk] + chi[-3*dj - dk + ijk] + 
       chi[3*dj + dk + ijk] + chi[dj + 3*dk + ijk]) + 
    5376.*(-chi[-2*dj - 3*dk + ijk] + chi[2*dj - 3*dk + ijk] - 
       chi[-3*dj - 2*dk + ijk] + chi[3*dj - 2*dk + ijk] + 
       chi[-3*dj + 2*dk + ijk] - chi[3*dj + 2*dk + ijk] + 
       chi[-2*dj + 3*dk + ijk] - chi[2*dj + 3*dk + ijk]) - 
    504.*(chi[2*dj - 4*dk + ijk] + chi[4*dj - 2*dk + ijk] + 
       chi[-4*dj + 2*dk + ijk] + chi[-2*dj + 4*dk + ijk]) + 
    2016.*(-chi[-dj - 4*dk + ijk] + chi[dj - 4*dk + ijk] - 
       chi[-4*dj - dk + ijk] + chi[4*dj - dk + ijk] + 
       chi[-4*dj + dk + ijk] - chi[4*dj + dk + ijk] + 
       chi[-dj + 4*dk + ijk] - chi[dj + 4*dk + ijk]) + 
    504.*(chi[-2*dj - 4*dk + ijk] + chi[-4*dj - 2*dk + ijk] + 
       chi[4*dj + 2*dk + ijk] + chi[2*dj + 4*dk + ijk]) + 
    96.*(-chi[-3*dj - 4*dk + ijk] + chi[3*dj - 4*dk + ijk] - 
       chi[-4*dj - 3*dk + ijk] + chi[4*dj - 3*dk + ijk] + 
       chi[-4*dj + 3*dk + ijk] - chi[4*dj + 3*dk + ijk] + 
       chi[-3*dj + 4*dk + ijk] - chi[3*dj + 4*dk + ijk]) - 
    28224.*(chi[2*(dj - dk) + ijk] + chi[2*(-dj + dk) + ijk]) - 
    1024.*(chi[3*(dj - dk) + ijk] + chi[3*(-dj + dk) + ijk]) - 
    9.*(chi[4*(dj - dk) + ijk] + chi[4*(-dj + dk) + ijk]) + 
    28224.*(chi[-2*(dj + dk) + ijk] + chi[2*(dj + dk) + ijk]) + 
    1024.*(chi[-3*(dj + dk) + ijk] + chi[3*(dj + dk) + ijk]) + 
    9.*(chi[-4*(dj + dk) + ijk] + chi[4*(dj + dk) + ijk]))
;

ddchi33
=
0.0001984126984126984127*oodz2*
  (-14350.*chi[ijk] + 8064.*(chi[-dk + ijk] + chi[dk + ijk]) - 
    1008.*(chi[-2*dk + ijk] + chi[2*dk + ijk]) + 
    128.*(chi[-3*dk + ijk] + chi[3*dk + ijk]) - 
    9.*(chi[-4*dk + ijk] + chi[4*dk + ijk]))
;


#endif 


} else { 


#ifdef REDUCEORDERTO8                           
  errorexit("CompilerFlag reduced order to 8");
 #else
da1
=
0.00079365079365079365079*oo2dx*
  (-2.*alpha[-5*di + ijk] + 25.*alpha[-4*di + ijk] + 
    600.*alpha[-2*di + ijk] + 2100.*(-alpha[-di + ijk] + alpha[di + ijk]) - 
    600.*alpha[2*di + ijk] + 150.*
     (-alpha[-3*di + ijk] + alpha[3*di + ijk]) - 25.*alpha[4*di + ijk] + 
    2.*alpha[5*di + ijk])
;

da2
=
0.00079365079365079365079*oo2dy*
  (-2.*alpha[-5*dj + ijk] + 25.*alpha[-4*dj + ijk] + 
    600.*alpha[-2*dj + ijk] + 2100.*(-alpha[-dj + ijk] + alpha[dj + ijk]) - 
    600.*alpha[2*dj + ijk] + 150.*
     (-alpha[-3*dj + ijk] + alpha[3*dj + ijk]) - 25.*alpha[4*dj + ijk] + 
    2.*alpha[5*dj + ijk])
;

da3
=
0.00079365079365079365079*oo2dz*
  (-2.*alpha[-5*dk + ijk] + 25.*alpha[-4*dk + ijk] + 
    600.*alpha[-2*dk + ijk] + 2100.*(-alpha[-dk + ijk] + alpha[dk + ijk]) - 
    600.*alpha[2*dk + ijk] + 150.*
     (-alpha[-3*dk + ijk] + alpha[3*dk + ijk]) - 25.*alpha[4*dk + ijk] + 
    2.*alpha[5*dk + ijk])
;

dda11
=
0.00003968253968253968254*oodx2*
  (-73766.*alpha[ijk] + 42000.*(alpha[-di + ijk] + alpha[di + ijk]) - 
    6000.*(alpha[-2*di + ijk] + alpha[2*di + ijk]) + 
    1000.*(alpha[-3*di + ijk] + alpha[3*di + ijk]) - 
    125.*(alpha[-4*di + ijk] + alpha[4*di + ijk]) + 
    8.*(alpha[-5*di + ijk] + alpha[5*di + ijk]))
;

dda12
=
6.2988158226253464349e-7*oo2dx*oo2dy*
  (4.41e6*(alpha[-di - dj + ijk] - alpha[di - dj + ijk] - 
       alpha[-di + dj + ijk] + alpha[di + dj + ijk]) + 
    1.26e6*(-alpha[-di - 2*dj + ijk] + alpha[di - 2*dj + ijk] - 
       alpha[-2*di - dj + ijk] + alpha[2*di - dj + ijk] + 
       alpha[-2*di + dj + ijk] - alpha[2*di + dj + ijk] + 
       alpha[-di + 2*dj + ijk] - alpha[di + 2*dj + ijk]) - 
    315000.*(alpha[di - 3*dj + ijk] + alpha[3*di - dj + ijk] + 
       alpha[-3*di + dj + ijk] + alpha[-di + 3*dj + ijk]) + 
    315000.*(alpha[-di - 3*dj + ijk] + alpha[-3*di - dj + ijk] + 
       alpha[3*di + dj + ijk] + alpha[di + 3*dj + ijk]) + 
    90000.*(-alpha[-2*di - 3*dj + ijk] + alpha[2*di - 3*dj + ijk] - 
       alpha[-3*di - 2*dj + ijk] + alpha[3*di - 2*dj + ijk] + 
       alpha[-3*di + 2*dj + ijk] - alpha[3*di + 2*dj + ijk] + 
       alpha[-2*di + 3*dj + ijk] - alpha[2*di + 3*dj + ijk]) - 
    15000.*(alpha[2*di - 4*dj + ijk] + alpha[4*di - 2*dj + ijk] + 
       alpha[-4*di + 2*dj + ijk] + alpha[-2*di + 4*dj + ijk]) + 
    52500.*(-alpha[-di - 4*dj + ijk] + alpha[di - 4*dj + ijk] - 
       alpha[-4*di - dj + ijk] + alpha[4*di - dj + ijk] + 
       alpha[-4*di + dj + ijk] - alpha[4*di + dj + ijk] + 
       alpha[-di + 4*dj + ijk] - alpha[di + 4*dj + ijk]) + 
    15000.*(alpha[-2*di - 4*dj + ijk] + alpha[-4*di - 2*dj + ijk] + 
       alpha[4*di + 2*dj + ijk] + alpha[2*di + 4*dj + ijk]) + 
    3750.*(-alpha[-3*di - 4*dj + ijk] + alpha[3*di - 4*dj + ijk] - 
       alpha[-4*di - 3*dj + ijk] + alpha[4*di - 3*dj + ijk] + 
       alpha[-4*di + 3*dj + ijk] - alpha[4*di + 3*dj + ijk] + 
       alpha[-3*di + 4*dj + ijk] - alpha[3*di + 4*dj + ijk]) - 
    300.*(alpha[3*di - 5*dj + ijk] + alpha[5*di - 3*dj + ijk] + 
       alpha[-5*di + 3*dj + ijk] + alpha[-3*di + 5*dj + ijk]) - 
    4200.*(alpha[di - 5*dj + ijk] + alpha[5*di - dj + ijk] + 
       alpha[-5*di + dj + ijk] + alpha[-di + 5*dj + ijk]) + 
    4200.*(alpha[-di - 5*dj + ijk] + alpha[-5*di - dj + ijk] + 
       alpha[5*di + dj + ijk] + alpha[di + 5*dj + ijk]) + 
    1200.*(-alpha[-2*di - 5*dj + ijk] + alpha[2*di - 5*dj + ijk] - 
       alpha[-5*di - 2*dj + ijk] + alpha[5*di - 2*dj + ijk] + 
       alpha[-5*di + 2*dj + ijk] - alpha[5*di + 2*dj + ijk] + 
       alpha[-2*di + 5*dj + ijk] - alpha[2*di + 5*dj + ijk]) + 
    300.*(alpha[-3*di - 5*dj + ijk] + alpha[-5*di - 3*dj + ijk] + 
       alpha[5*di + 3*dj + ijk] + alpha[3*di + 5*dj + ijk]) + 
    50.*(-alpha[-4*di - 5*dj + ijk] + alpha[4*di - 5*dj + ijk] - 
       alpha[-5*di - 4*dj + ijk] + alpha[5*di - 4*dj + ijk] + 
       alpha[-5*di + 4*dj + ijk] - alpha[5*di + 4*dj + ijk] + 
       alpha[-4*di + 5*dj + ijk] - alpha[4*di + 5*dj + ijk]) - 
    360000.*(alpha[2*(di - dj) + ijk] + alpha[2*(-di + dj) + ijk]) - 
    22500.*(alpha[3*(di - dj) + ijk] + alpha[3*(-di + dj) + ijk]) - 
    625.*(alpha[4*(di - dj) + ijk] + alpha[4*(-di + dj) + ijk]) - 
    4.*(alpha[5*(di - dj) + ijk] + alpha[5*(-di + dj) + ijk]) + 
    360000.*(alpha[-2*(di + dj) + ijk] + alpha[2*(di + dj) + ijk]) + 
    22500.*(alpha[-3*(di + dj) + ijk] + alpha[3*(di + dj) + ijk]) + 
    625.*(alpha[-4*(di + dj) + ijk] + alpha[4*(di + dj) + ijk]) + 
    4.*(alpha[-5*(di + dj) + ijk] + alpha[5*(di + dj) + ijk]))
;

dda13
=
6.2988158226253464349e-7*oo2dx*oo2dz*
  (4.41e6*(alpha[-di - dk + ijk] - alpha[di - dk + ijk] - 
       alpha[-di + dk + ijk] + alpha[di + dk + ijk]) + 
    1.26e6*(-alpha[-di - 2*dk + ijk] + alpha[di - 2*dk + ijk] - 
       alpha[-2*di - dk + ijk] + alpha[2*di - dk + ijk] + 
       alpha[-2*di + dk + ijk] - alpha[2*di + dk + ijk] + 
       alpha[-di + 2*dk + ijk] - alpha[di + 2*dk + ijk]) - 
    315000.*(alpha[di - 3*dk + ijk] + alpha[3*di - dk + ijk] + 
       alpha[-3*di + dk + ijk] + alpha[-di + 3*dk + ijk]) + 
    315000.*(alpha[-di - 3*dk + ijk] + alpha[-3*di - dk + ijk] + 
       alpha[3*di + dk + ijk] + alpha[di + 3*dk + ijk]) + 
    90000.*(-alpha[-2*di - 3*dk + ijk] + alpha[2*di - 3*dk + ijk] - 
       alpha[-3*di - 2*dk + ijk] + alpha[3*di - 2*dk + ijk] + 
       alpha[-3*di + 2*dk + ijk] - alpha[3*di + 2*dk + ijk] + 
       alpha[-2*di + 3*dk + ijk] - alpha[2*di + 3*dk + ijk]) - 
    15000.*(alpha[2*di - 4*dk + ijk] + alpha[4*di - 2*dk + ijk] + 
       alpha[-4*di + 2*dk + ijk] + alpha[-2*di + 4*dk + ijk]) + 
    52500.*(-alpha[-di - 4*dk + ijk] + alpha[di - 4*dk + ijk] - 
       alpha[-4*di - dk + ijk] + alpha[4*di - dk + ijk] + 
       alpha[-4*di + dk + ijk] - alpha[4*di + dk + ijk] + 
       alpha[-di + 4*dk + ijk] - alpha[di + 4*dk + ijk]) + 
    15000.*(alpha[-2*di - 4*dk + ijk] + alpha[-4*di - 2*dk + ijk] + 
       alpha[4*di + 2*dk + ijk] + alpha[2*di + 4*dk + ijk]) + 
    3750.*(-alpha[-3*di - 4*dk + ijk] + alpha[3*di - 4*dk + ijk] - 
       alpha[-4*di - 3*dk + ijk] + alpha[4*di - 3*dk + ijk] + 
       alpha[-4*di + 3*dk + ijk] - alpha[4*di + 3*dk + ijk] + 
       alpha[-3*di + 4*dk + ijk] - alpha[3*di + 4*dk + ijk]) - 
    300.*(alpha[3*di - 5*dk + ijk] + alpha[5*di - 3*dk + ijk] + 
       alpha[-5*di + 3*dk + ijk] + alpha[-3*di + 5*dk + ijk]) - 
    4200.*(alpha[di - 5*dk + ijk] + alpha[5*di - dk + ijk] + 
       alpha[-5*di + dk + ijk] + alpha[-di + 5*dk + ijk]) + 
    4200.*(alpha[-di - 5*dk + ijk] + alpha[-5*di - dk + ijk] + 
       alpha[5*di + dk + ijk] + alpha[di + 5*dk + ijk]) + 
    1200.*(-alpha[-2*di - 5*dk + ijk] + alpha[2*di - 5*dk + ijk] - 
       alpha[-5*di - 2*dk + ijk] + alpha[5*di - 2*dk + ijk] + 
       alpha[-5*di + 2*dk + ijk] - alpha[5*di + 2*dk + ijk] + 
       alpha[-2*di + 5*dk + ijk] - alpha[2*di + 5*dk + ijk]) + 
    300.*(alpha[-3*di - 5*dk + ijk] + alpha[-5*di - 3*dk + ijk] + 
       alpha[5*di + 3*dk + ijk] + alpha[3*di + 5*dk + ijk]) + 
    50.*(-alpha[-4*di - 5*dk + ijk] + alpha[4*di - 5*dk + ijk] - 
       alpha[-5*di - 4*dk + ijk] + alpha[5*di - 4*dk + ijk] + 
       alpha[-5*di + 4*dk + ijk] - alpha[5*di + 4*dk + ijk] + 
       alpha[-4*di + 5*dk + ijk] - alpha[4*di + 5*dk + ijk]) - 
    360000.*(alpha[2*(di - dk) + ijk] + alpha[2*(-di + dk) + ijk]) - 
    22500.*(alpha[3*(di - dk) + ijk] + alpha[3*(-di + dk) + ijk]) - 
    625.*(alpha[4*(di - dk) + ijk] + alpha[4*(-di + dk) + ijk]) - 
    4.*(alpha[5*(di - dk) + ijk] + alpha[5*(-di + dk) + ijk]) + 
    360000.*(alpha[-2*(di + dk) + ijk] + alpha[2*(di + dk) + ijk]) + 
    22500.*(alpha[-3*(di + dk) + ijk] + alpha[3*(di + dk) + ijk]) + 
    625.*(alpha[-4*(di + dk) + ijk] + alpha[4*(di + dk) + ijk]) + 
    4.*(alpha[-5*(di + dk) + ijk] + alpha[5*(di + dk) + ijk]))
;

dda22
=
0.00003968253968253968254*oody2*
  (-73766.*alpha[ijk] + 42000.*(alpha[-dj + ijk] + alpha[dj + ijk]) - 
    6000.*(alpha[-2*dj + ijk] + alpha[2*dj + ijk]) + 
    1000.*(alpha[-3*dj + ijk] + alpha[3*dj + ijk]) - 
    125.*(alpha[-4*dj + ijk] + alpha[4*dj + ijk]) + 
    8.*(alpha[-5*dj + ijk] + alpha[5*dj + ijk]))
;

dda23
=
6.2988158226253464349e-7*oo2dy*oo2dz*
  (4.41e6*(alpha[-dj - dk + ijk] - alpha[dj - dk + ijk] - 
       alpha[-dj + dk + ijk] + alpha[dj + dk + ijk]) + 
    1.26e6*(-alpha[-dj - 2*dk + ijk] + alpha[dj - 2*dk + ijk] - 
       alpha[-2*dj - dk + ijk] + alpha[2*dj - dk + ijk] + 
       alpha[-2*dj + dk + ijk] - alpha[2*dj + dk + ijk] + 
       alpha[-dj + 2*dk + ijk] - alpha[dj + 2*dk + ijk]) - 
    315000.*(alpha[dj - 3*dk + ijk] + alpha[3*dj - dk + ijk] + 
       alpha[-3*dj + dk + ijk] + alpha[-dj + 3*dk + ijk]) + 
    315000.*(alpha[-dj - 3*dk + ijk] + alpha[-3*dj - dk + ijk] + 
       alpha[3*dj + dk + ijk] + alpha[dj + 3*dk + ijk]) + 
    90000.*(-alpha[-2*dj - 3*dk + ijk] + alpha[2*dj - 3*dk + ijk] - 
       alpha[-3*dj - 2*dk + ijk] + alpha[3*dj - 2*dk + ijk] + 
       alpha[-3*dj + 2*dk + ijk] - alpha[3*dj + 2*dk + ijk] + 
       alpha[-2*dj + 3*dk + ijk] - alpha[2*dj + 3*dk + ijk]) - 
    15000.*(alpha[2*dj - 4*dk + ijk] + alpha[4*dj - 2*dk + ijk] + 
       alpha[-4*dj + 2*dk + ijk] + alpha[-2*dj + 4*dk + ijk]) + 
    52500.*(-alpha[-dj - 4*dk + ijk] + alpha[dj - 4*dk + ijk] - 
       alpha[-4*dj - dk + ijk] + alpha[4*dj - dk + ijk] + 
       alpha[-4*dj + dk + ijk] - alpha[4*dj + dk + ijk] + 
       alpha[-dj + 4*dk + ijk] - alpha[dj + 4*dk + ijk]) + 
    15000.*(alpha[-2*dj - 4*dk + ijk] + alpha[-4*dj - 2*dk + ijk] + 
       alpha[4*dj + 2*dk + ijk] + alpha[2*dj + 4*dk + ijk]) + 
    3750.*(-alpha[-3*dj - 4*dk + ijk] + alpha[3*dj - 4*dk + ijk] - 
       alpha[-4*dj - 3*dk + ijk] + alpha[4*dj - 3*dk + ijk] + 
       alpha[-4*dj + 3*dk + ijk] - alpha[4*dj + 3*dk + ijk] + 
       alpha[-3*dj + 4*dk + ijk] - alpha[3*dj + 4*dk + ijk]) - 
    300.*(alpha[3*dj - 5*dk + ijk] + alpha[5*dj - 3*dk + ijk] + 
       alpha[-5*dj + 3*dk + ijk] + alpha[-3*dj + 5*dk + ijk]) - 
    4200.*(alpha[dj - 5*dk + ijk] + alpha[5*dj - dk + ijk] + 
       alpha[-5*dj + dk + ijk] + alpha[-dj + 5*dk + ijk]) + 
    4200.*(alpha[-dj - 5*dk + ijk] + alpha[-5*dj - dk + ijk] + 
       alpha[5*dj + dk + ijk] + alpha[dj + 5*dk + ijk]) + 
    1200.*(-alpha[-2*dj - 5*dk + ijk] + alpha[2*dj - 5*dk + ijk] - 
       alpha[-5*dj - 2*dk + ijk] + alpha[5*dj - 2*dk + ijk] + 
       alpha[-5*dj + 2*dk + ijk] - alpha[5*dj + 2*dk + ijk] + 
       alpha[-2*dj + 5*dk + ijk] - alpha[2*dj + 5*dk + ijk]) + 
    300.*(alpha[-3*dj - 5*dk + ijk] + alpha[-5*dj - 3*dk + ijk] + 
       alpha[5*dj + 3*dk + ijk] + alpha[3*dj + 5*dk + ijk]) + 
    50.*(-alpha[-4*dj - 5*dk + ijk] + alpha[4*dj - 5*dk + ijk] - 
       alpha[-5*dj - 4*dk + ijk] + alpha[5*dj - 4*dk + ijk] + 
       alpha[-5*dj + 4*dk + ijk] - alpha[5*dj + 4*dk + ijk] + 
       alpha[-4*dj + 5*dk + ijk] - alpha[4*dj + 5*dk + ijk]) - 
    360000.*(alpha[2*(dj - dk) + ijk] + alpha[2*(-dj + dk) + ijk]) - 
    22500.*(alpha[3*(dj - dk) + ijk] + alpha[3*(-dj + dk) + ijk]) - 
    625.*(alpha[4*(dj - dk) + ijk] + alpha[4*(-dj + dk) + ijk]) - 
    4.*(alpha[5*(dj - dk) + ijk] + alpha[5*(-dj + dk) + ijk]) + 
    360000.*(alpha[-2*(dj + dk) + ijk] + alpha[2*(dj + dk) + ijk]) + 
    22500.*(alpha[-3*(dj + dk) + ijk] + alpha[3*(dj + dk) + ijk]) + 
    625.*(alpha[-4*(dj + dk) + ijk] + alpha[4*(dj + dk) + ijk]) + 
    4.*(alpha[-5*(dj + dk) + ijk] + alpha[5*(dj + dk) + ijk]))
;

dda33
=
0.00003968253968253968254*oodz2*
  (-73766.*alpha[ijk] + 42000.*(alpha[-dk + ijk] + alpha[dk + ijk]) - 
    6000.*(alpha[-2*dk + ijk] + alpha[2*dk + ijk]) + 
    1000.*(alpha[-3*dk + ijk] + alpha[3*dk + ijk]) - 
    125.*(alpha[-4*dk + ijk] + alpha[4*dk + ijk]) + 
    8.*(alpha[-5*dk + ijk] + alpha[5*dk + ijk]))
;

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

ddb111
=
0.00003968253968253968254*oodx2*
  (-73766.*beta1[ijk] + 42000.*(beta1[-di + ijk] + beta1[di + ijk]) - 
    6000.*(beta1[-2*di + ijk] + beta1[2*di + ijk]) + 
    1000.*(beta1[-3*di + ijk] + beta1[3*di + ijk]) - 
    125.*(beta1[-4*di + ijk] + beta1[4*di + ijk]) + 
    8.*(beta1[-5*di + ijk] + beta1[5*di + ijk]))
;

ddb112
=
0.00003968253968253968254*oodx2*
  (-73766.*beta2[ijk] + 42000.*(beta2[-di + ijk] + beta2[di + ijk]) - 
    6000.*(beta2[-2*di + ijk] + beta2[2*di + ijk]) + 
    1000.*(beta2[-3*di + ijk] + beta2[3*di + ijk]) - 
    125.*(beta2[-4*di + ijk] + beta2[4*di + ijk]) + 
    8.*(beta2[-5*di + ijk] + beta2[5*di + ijk]))
;

ddb113
=
0.00003968253968253968254*oodx2*
  (-73766.*beta3[ijk] + 42000.*(beta3[-di + ijk] + beta3[di + ijk]) - 
    6000.*(beta3[-2*di + ijk] + beta3[2*di + ijk]) + 
    1000.*(beta3[-3*di + ijk] + beta3[3*di + ijk]) - 
    125.*(beta3[-4*di + ijk] + beta3[4*di + ijk]) + 
    8.*(beta3[-5*di + ijk] + beta3[5*di + ijk]))
;

ddb121
=
6.2988158226253464349e-7*oo2dx*oo2dy*
  (4.41e6*(beta1[-di - dj + ijk] - beta1[di - dj + ijk] - 
       beta1[-di + dj + ijk] + beta1[di + dj + ijk]) + 
    1.26e6*(-beta1[-di - 2*dj + ijk] + beta1[di - 2*dj + ijk] - 
       beta1[-2*di - dj + ijk] + beta1[2*di - dj + ijk] + 
       beta1[-2*di + dj + ijk] - beta1[2*di + dj + ijk] + 
       beta1[-di + 2*dj + ijk] - beta1[di + 2*dj + ijk]) - 
    315000.*(beta1[di - 3*dj + ijk] + beta1[3*di - dj + ijk] + 
       beta1[-3*di + dj + ijk] + beta1[-di + 3*dj + ijk]) + 
    315000.*(beta1[-di - 3*dj + ijk] + beta1[-3*di - dj + ijk] + 
       beta1[3*di + dj + ijk] + beta1[di + 3*dj + ijk]) + 
    90000.*(-beta1[-2*di - 3*dj + ijk] + beta1[2*di - 3*dj + ijk] - 
       beta1[-3*di - 2*dj + ijk] + beta1[3*di - 2*dj + ijk] + 
       beta1[-3*di + 2*dj + ijk] - beta1[3*di + 2*dj + ijk] + 
       beta1[-2*di + 3*dj + ijk] - beta1[2*di + 3*dj + ijk]) - 
    15000.*(beta1[2*di - 4*dj + ijk] + beta1[4*di - 2*dj + ijk] + 
       beta1[-4*di + 2*dj + ijk] + beta1[-2*di + 4*dj + ijk]) + 
    52500.*(-beta1[-di - 4*dj + ijk] + beta1[di - 4*dj + ijk] - 
       beta1[-4*di - dj + ijk] + beta1[4*di - dj + ijk] + 
       beta1[-4*di + dj + ijk] - beta1[4*di + dj + ijk] + 
       beta1[-di + 4*dj + ijk] - beta1[di + 4*dj + ijk]) + 
    15000.*(beta1[-2*di - 4*dj + ijk] + beta1[-4*di - 2*dj + ijk] + 
       beta1[4*di + 2*dj + ijk] + beta1[2*di + 4*dj + ijk]) + 
    3750.*(-beta1[-3*di - 4*dj + ijk] + beta1[3*di - 4*dj + ijk] - 
       beta1[-4*di - 3*dj + ijk] + beta1[4*di - 3*dj + ijk] + 
       beta1[-4*di + 3*dj + ijk] - beta1[4*di + 3*dj + ijk] + 
       beta1[-3*di + 4*dj + ijk] - beta1[3*di + 4*dj + ijk]) - 
    300.*(beta1[3*di - 5*dj + ijk] + beta1[5*di - 3*dj + ijk] + 
       beta1[-5*di + 3*dj + ijk] + beta1[-3*di + 5*dj + ijk]) - 
    4200.*(beta1[di - 5*dj + ijk] + beta1[5*di - dj + ijk] + 
       beta1[-5*di + dj + ijk] + beta1[-di + 5*dj + ijk]) + 
    4200.*(beta1[-di - 5*dj + ijk] + beta1[-5*di - dj + ijk] + 
       beta1[5*di + dj + ijk] + beta1[di + 5*dj + ijk]) + 
    1200.*(-beta1[-2*di - 5*dj + ijk] + beta1[2*di - 5*dj + ijk] - 
       beta1[-5*di - 2*dj + ijk] + beta1[5*di - 2*dj + ijk] + 
       beta1[-5*di + 2*dj + ijk] - beta1[5*di + 2*dj + ijk] + 
       beta1[-2*di + 5*dj + ijk] - beta1[2*di + 5*dj + ijk]) + 
    300.*(beta1[-3*di - 5*dj + ijk] + beta1[-5*di - 3*dj + ijk] + 
       beta1[5*di + 3*dj + ijk] + beta1[3*di + 5*dj + ijk]) + 
    50.*(-beta1[-4*di - 5*dj + ijk] + beta1[4*di - 5*dj + ijk] - 
       beta1[-5*di - 4*dj + ijk] + beta1[5*di - 4*dj + ijk] + 
       beta1[-5*di + 4*dj + ijk] - beta1[5*di + 4*dj + ijk] + 
       beta1[-4*di + 5*dj + ijk] - beta1[4*di + 5*dj + ijk]) - 
    360000.*(beta1[2*(di - dj) + ijk] + beta1[2*(-di + dj) + ijk]) - 
    22500.*(beta1[3*(di - dj) + ijk] + beta1[3*(-di + dj) + ijk]) - 
    625.*(beta1[4*(di - dj) + ijk] + beta1[4*(-di + dj) + ijk]) - 
    4.*(beta1[5*(di - dj) + ijk] + beta1[5*(-di + dj) + ijk]) + 
    360000.*(beta1[-2*(di + dj) + ijk] + beta1[2*(di + dj) + ijk]) + 
    22500.*(beta1[-3*(di + dj) + ijk] + beta1[3*(di + dj) + ijk]) + 
    625.*(beta1[-4*(di + dj) + ijk] + beta1[4*(di + dj) + ijk]) + 
    4.*(beta1[-5*(di + dj) + ijk] + beta1[5*(di + dj) + ijk]))
;

ddb122
=
6.2988158226253464349e-7*oo2dx*oo2dy*
  (4.41e6*(beta2[-di - dj + ijk] - beta2[di - dj + ijk] - 
       beta2[-di + dj + ijk] + beta2[di + dj + ijk]) + 
    1.26e6*(-beta2[-di - 2*dj + ijk] + beta2[di - 2*dj + ijk] - 
       beta2[-2*di - dj + ijk] + beta2[2*di - dj + ijk] + 
       beta2[-2*di + dj + ijk] - beta2[2*di + dj + ijk] + 
       beta2[-di + 2*dj + ijk] - beta2[di + 2*dj + ijk]) - 
    315000.*(beta2[di - 3*dj + ijk] + beta2[3*di - dj + ijk] + 
       beta2[-3*di + dj + ijk] + beta2[-di + 3*dj + ijk]) + 
    315000.*(beta2[-di - 3*dj + ijk] + beta2[-3*di - dj + ijk] + 
       beta2[3*di + dj + ijk] + beta2[di + 3*dj + ijk]) + 
    90000.*(-beta2[-2*di - 3*dj + ijk] + beta2[2*di - 3*dj + ijk] - 
       beta2[-3*di - 2*dj + ijk] + beta2[3*di - 2*dj + ijk] + 
       beta2[-3*di + 2*dj + ijk] - beta2[3*di + 2*dj + ijk] + 
       beta2[-2*di + 3*dj + ijk] - beta2[2*di + 3*dj + ijk]) - 
    15000.*(beta2[2*di - 4*dj + ijk] + beta2[4*di - 2*dj + ijk] + 
       beta2[-4*di + 2*dj + ijk] + beta2[-2*di + 4*dj + ijk]) + 
    52500.*(-beta2[-di - 4*dj + ijk] + beta2[di - 4*dj + ijk] - 
       beta2[-4*di - dj + ijk] + beta2[4*di - dj + ijk] + 
       beta2[-4*di + dj + ijk] - beta2[4*di + dj + ijk] + 
       beta2[-di + 4*dj + ijk] - beta2[di + 4*dj + ijk]) + 
    15000.*(beta2[-2*di - 4*dj + ijk] + beta2[-4*di - 2*dj + ijk] + 
       beta2[4*di + 2*dj + ijk] + beta2[2*di + 4*dj + ijk]) + 
    3750.*(-beta2[-3*di - 4*dj + ijk] + beta2[3*di - 4*dj + ijk] - 
       beta2[-4*di - 3*dj + ijk] + beta2[4*di - 3*dj + ijk] + 
       beta2[-4*di + 3*dj + ijk] - beta2[4*di + 3*dj + ijk] + 
       beta2[-3*di + 4*dj + ijk] - beta2[3*di + 4*dj + ijk]) - 
    300.*(beta2[3*di - 5*dj + ijk] + beta2[5*di - 3*dj + ijk] + 
       beta2[-5*di + 3*dj + ijk] + beta2[-3*di + 5*dj + ijk]) - 
    4200.*(beta2[di - 5*dj + ijk] + beta2[5*di - dj + ijk] + 
       beta2[-5*di + dj + ijk] + beta2[-di + 5*dj + ijk]) + 
    4200.*(beta2[-di - 5*dj + ijk] + beta2[-5*di - dj + ijk] + 
       beta2[5*di + dj + ijk] + beta2[di + 5*dj + ijk]) + 
    1200.*(-beta2[-2*di - 5*dj + ijk] + beta2[2*di - 5*dj + ijk] - 
       beta2[-5*di - 2*dj + ijk] + beta2[5*di - 2*dj + ijk] + 
       beta2[-5*di + 2*dj + ijk] - beta2[5*di + 2*dj + ijk] + 
       beta2[-2*di + 5*dj + ijk] - beta2[2*di + 5*dj + ijk]) + 
    300.*(beta2[-3*di - 5*dj + ijk] + beta2[-5*di - 3*dj + ijk] + 
       beta2[5*di + 3*dj + ijk] + beta2[3*di + 5*dj + ijk]) + 
    50.*(-beta2[-4*di - 5*dj + ijk] + beta2[4*di - 5*dj + ijk] - 
       beta2[-5*di - 4*dj + ijk] + beta2[5*di - 4*dj + ijk] + 
       beta2[-5*di + 4*dj + ijk] - beta2[5*di + 4*dj + ijk] + 
       beta2[-4*di + 5*dj + ijk] - beta2[4*di + 5*dj + ijk]) - 
    360000.*(beta2[2*(di - dj) + ijk] + beta2[2*(-di + dj) + ijk]) - 
    22500.*(beta2[3*(di - dj) + ijk] + beta2[3*(-di + dj) + ijk]) - 
    625.*(beta2[4*(di - dj) + ijk] + beta2[4*(-di + dj) + ijk]) - 
    4.*(beta2[5*(di - dj) + ijk] + beta2[5*(-di + dj) + ijk]) + 
    360000.*(beta2[-2*(di + dj) + ijk] + beta2[2*(di + dj) + ijk]) + 
    22500.*(beta2[-3*(di + dj) + ijk] + beta2[3*(di + dj) + ijk]) + 
    625.*(beta2[-4*(di + dj) + ijk] + beta2[4*(di + dj) + ijk]) + 
    4.*(beta2[-5*(di + dj) + ijk] + beta2[5*(di + dj) + ijk]))
;

ddb123
=
6.2988158226253464349e-7*oo2dx*oo2dy*
  (4.41e6*(beta3[-di - dj + ijk] - beta3[di - dj + ijk] - 
       beta3[-di + dj + ijk] + beta3[di + dj + ijk]) + 
    1.26e6*(-beta3[-di - 2*dj + ijk] + beta3[di - 2*dj + ijk] - 
       beta3[-2*di - dj + ijk] + beta3[2*di - dj + ijk] + 
       beta3[-2*di + dj + ijk] - beta3[2*di + dj + ijk] + 
       beta3[-di + 2*dj + ijk] - beta3[di + 2*dj + ijk]) - 
    315000.*(beta3[di - 3*dj + ijk] + beta3[3*di - dj + ijk] + 
       beta3[-3*di + dj + ijk] + beta3[-di + 3*dj + ijk]) + 
    315000.*(beta3[-di - 3*dj + ijk] + beta3[-3*di - dj + ijk] + 
       beta3[3*di + dj + ijk] + beta3[di + 3*dj + ijk]) + 
    90000.*(-beta3[-2*di - 3*dj + ijk] + beta3[2*di - 3*dj + ijk] - 
       beta3[-3*di - 2*dj + ijk] + beta3[3*di - 2*dj + ijk] + 
       beta3[-3*di + 2*dj + ijk] - beta3[3*di + 2*dj + ijk] + 
       beta3[-2*di + 3*dj + ijk] - beta3[2*di + 3*dj + ijk]) - 
    15000.*(beta3[2*di - 4*dj + ijk] + beta3[4*di - 2*dj + ijk] + 
       beta3[-4*di + 2*dj + ijk] + beta3[-2*di + 4*dj + ijk]) + 
    52500.*(-beta3[-di - 4*dj + ijk] + beta3[di - 4*dj + ijk] - 
       beta3[-4*di - dj + ijk] + beta3[4*di - dj + ijk] + 
       beta3[-4*di + dj + ijk] - beta3[4*di + dj + ijk] + 
       beta3[-di + 4*dj + ijk] - beta3[di + 4*dj + ijk]) + 
    15000.*(beta3[-2*di - 4*dj + ijk] + beta3[-4*di - 2*dj + ijk] + 
       beta3[4*di + 2*dj + ijk] + beta3[2*di + 4*dj + ijk]) + 
    3750.*(-beta3[-3*di - 4*dj + ijk] + beta3[3*di - 4*dj + ijk] - 
       beta3[-4*di - 3*dj + ijk] + beta3[4*di - 3*dj + ijk] + 
       beta3[-4*di + 3*dj + ijk] - beta3[4*di + 3*dj + ijk] + 
       beta3[-3*di + 4*dj + ijk] - beta3[3*di + 4*dj + ijk]) - 
    300.*(beta3[3*di - 5*dj + ijk] + beta3[5*di - 3*dj + ijk] + 
       beta3[-5*di + 3*dj + ijk] + beta3[-3*di + 5*dj + ijk]) - 
    4200.*(beta3[di - 5*dj + ijk] + beta3[5*di - dj + ijk] + 
       beta3[-5*di + dj + ijk] + beta3[-di + 5*dj + ijk]) + 
    4200.*(beta3[-di - 5*dj + ijk] + beta3[-5*di - dj + ijk] + 
       beta3[5*di + dj + ijk] + beta3[di + 5*dj + ijk]) + 
    1200.*(-beta3[-2*di - 5*dj + ijk] + beta3[2*di - 5*dj + ijk] - 
       beta3[-5*di - 2*dj + ijk] + beta3[5*di - 2*dj + ijk] + 
       beta3[-5*di + 2*dj + ijk] - beta3[5*di + 2*dj + ijk] + 
       beta3[-2*di + 5*dj + ijk] - beta3[2*di + 5*dj + ijk]) + 
    300.*(beta3[-3*di - 5*dj + ijk] + beta3[-5*di - 3*dj + ijk] + 
       beta3[5*di + 3*dj + ijk] + beta3[3*di + 5*dj + ijk]) + 
    50.*(-beta3[-4*di - 5*dj + ijk] + beta3[4*di - 5*dj + ijk] - 
       beta3[-5*di - 4*dj + ijk] + beta3[5*di - 4*dj + ijk] + 
       beta3[-5*di + 4*dj + ijk] - beta3[5*di + 4*dj + ijk] + 
       beta3[-4*di + 5*dj + ijk] - beta3[4*di + 5*dj + ijk]) - 
    360000.*(beta3[2*(di - dj) + ijk] + beta3[2*(-di + dj) + ijk]) - 
    22500.*(beta3[3*(di - dj) + ijk] + beta3[3*(-di + dj) + ijk]) - 
    625.*(beta3[4*(di - dj) + ijk] + beta3[4*(-di + dj) + ijk]) - 
    4.*(beta3[5*(di - dj) + ijk] + beta3[5*(-di + dj) + ijk]) + 
    360000.*(beta3[-2*(di + dj) + ijk] + beta3[2*(di + dj) + ijk]) + 
    22500.*(beta3[-3*(di + dj) + ijk] + beta3[3*(di + dj) + ijk]) + 
    625.*(beta3[-4*(di + dj) + ijk] + beta3[4*(di + dj) + ijk]) + 
    4.*(beta3[-5*(di + dj) + ijk] + beta3[5*(di + dj) + ijk]))
;

ddb131
=
6.2988158226253464349e-7*oo2dx*oo2dz*
  (4.41e6*(beta1[-di - dk + ijk] - beta1[di - dk + ijk] - 
       beta1[-di + dk + ijk] + beta1[di + dk + ijk]) + 
    1.26e6*(-beta1[-di - 2*dk + ijk] + beta1[di - 2*dk + ijk] - 
       beta1[-2*di - dk + ijk] + beta1[2*di - dk + ijk] + 
       beta1[-2*di + dk + ijk] - beta1[2*di + dk + ijk] + 
       beta1[-di + 2*dk + ijk] - beta1[di + 2*dk + ijk]) - 
    315000.*(beta1[di - 3*dk + ijk] + beta1[3*di - dk + ijk] + 
       beta1[-3*di + dk + ijk] + beta1[-di + 3*dk + ijk]) + 
    315000.*(beta1[-di - 3*dk + ijk] + beta1[-3*di - dk + ijk] + 
       beta1[3*di + dk + ijk] + beta1[di + 3*dk + ijk]) + 
    90000.*(-beta1[-2*di - 3*dk + ijk] + beta1[2*di - 3*dk + ijk] - 
       beta1[-3*di - 2*dk + ijk] + beta1[3*di - 2*dk + ijk] + 
       beta1[-3*di + 2*dk + ijk] - beta1[3*di + 2*dk + ijk] + 
       beta1[-2*di + 3*dk + ijk] - beta1[2*di + 3*dk + ijk]) - 
    15000.*(beta1[2*di - 4*dk + ijk] + beta1[4*di - 2*dk + ijk] + 
       beta1[-4*di + 2*dk + ijk] + beta1[-2*di + 4*dk + ijk]) + 
    52500.*(-beta1[-di - 4*dk + ijk] + beta1[di - 4*dk + ijk] - 
       beta1[-4*di - dk + ijk] + beta1[4*di - dk + ijk] + 
       beta1[-4*di + dk + ijk] - beta1[4*di + dk + ijk] + 
       beta1[-di + 4*dk + ijk] - beta1[di + 4*dk + ijk]) + 
    15000.*(beta1[-2*di - 4*dk + ijk] + beta1[-4*di - 2*dk + ijk] + 
       beta1[4*di + 2*dk + ijk] + beta1[2*di + 4*dk + ijk]) + 
    3750.*(-beta1[-3*di - 4*dk + ijk] + beta1[3*di - 4*dk + ijk] - 
       beta1[-4*di - 3*dk + ijk] + beta1[4*di - 3*dk + ijk] + 
       beta1[-4*di + 3*dk + ijk] - beta1[4*di + 3*dk + ijk] + 
       beta1[-3*di + 4*dk + ijk] - beta1[3*di + 4*dk + ijk]) - 
    300.*(beta1[3*di - 5*dk + ijk] + beta1[5*di - 3*dk + ijk] + 
       beta1[-5*di + 3*dk + ijk] + beta1[-3*di + 5*dk + ijk]) - 
    4200.*(beta1[di - 5*dk + ijk] + beta1[5*di - dk + ijk] + 
       beta1[-5*di + dk + ijk] + beta1[-di + 5*dk + ijk]) + 
    4200.*(beta1[-di - 5*dk + ijk] + beta1[-5*di - dk + ijk] + 
       beta1[5*di + dk + ijk] + beta1[di + 5*dk + ijk]) + 
    1200.*(-beta1[-2*di - 5*dk + ijk] + beta1[2*di - 5*dk + ijk] - 
       beta1[-5*di - 2*dk + ijk] + beta1[5*di - 2*dk + ijk] + 
       beta1[-5*di + 2*dk + ijk] - beta1[5*di + 2*dk + ijk] + 
       beta1[-2*di + 5*dk + ijk] - beta1[2*di + 5*dk + ijk]) + 
    300.*(beta1[-3*di - 5*dk + ijk] + beta1[-5*di - 3*dk + ijk] + 
       beta1[5*di + 3*dk + ijk] + beta1[3*di + 5*dk + ijk]) + 
    50.*(-beta1[-4*di - 5*dk + ijk] + beta1[4*di - 5*dk + ijk] - 
       beta1[-5*di - 4*dk + ijk] + beta1[5*di - 4*dk + ijk] + 
       beta1[-5*di + 4*dk + ijk] - beta1[5*di + 4*dk + ijk] + 
       beta1[-4*di + 5*dk + ijk] - beta1[4*di + 5*dk + ijk]) - 
    360000.*(beta1[2*(di - dk) + ijk] + beta1[2*(-di + dk) + ijk]) - 
    22500.*(beta1[3*(di - dk) + ijk] + beta1[3*(-di + dk) + ijk]) - 
    625.*(beta1[4*(di - dk) + ijk] + beta1[4*(-di + dk) + ijk]) - 
    4.*(beta1[5*(di - dk) + ijk] + beta1[5*(-di + dk) + ijk]) + 
    360000.*(beta1[-2*(di + dk) + ijk] + beta1[2*(di + dk) + ijk]) + 
    22500.*(beta1[-3*(di + dk) + ijk] + beta1[3*(di + dk) + ijk]) + 
    625.*(beta1[-4*(di + dk) + ijk] + beta1[4*(di + dk) + ijk]) + 
    4.*(beta1[-5*(di + dk) + ijk] + beta1[5*(di + dk) + ijk]))
;

ddb132
=
6.2988158226253464349e-7*oo2dx*oo2dz*
  (4.41e6*(beta2[-di - dk + ijk] - beta2[di - dk + ijk] - 
       beta2[-di + dk + ijk] + beta2[di + dk + ijk]) + 
    1.26e6*(-beta2[-di - 2*dk + ijk] + beta2[di - 2*dk + ijk] - 
       beta2[-2*di - dk + ijk] + beta2[2*di - dk + ijk] + 
       beta2[-2*di + dk + ijk] - beta2[2*di + dk + ijk] + 
       beta2[-di + 2*dk + ijk] - beta2[di + 2*dk + ijk]) - 
    315000.*(beta2[di - 3*dk + ijk] + beta2[3*di - dk + ijk] + 
       beta2[-3*di + dk + ijk] + beta2[-di + 3*dk + ijk]) + 
    315000.*(beta2[-di - 3*dk + ijk] + beta2[-3*di - dk + ijk] + 
       beta2[3*di + dk + ijk] + beta2[di + 3*dk + ijk]) + 
    90000.*(-beta2[-2*di - 3*dk + ijk] + beta2[2*di - 3*dk + ijk] - 
       beta2[-3*di - 2*dk + ijk] + beta2[3*di - 2*dk + ijk] + 
       beta2[-3*di + 2*dk + ijk] - beta2[3*di + 2*dk + ijk] + 
       beta2[-2*di + 3*dk + ijk] - beta2[2*di + 3*dk + ijk]) - 
    15000.*(beta2[2*di - 4*dk + ijk] + beta2[4*di - 2*dk + ijk] + 
       beta2[-4*di + 2*dk + ijk] + beta2[-2*di + 4*dk + ijk]) + 
    52500.*(-beta2[-di - 4*dk + ijk] + beta2[di - 4*dk + ijk] - 
       beta2[-4*di - dk + ijk] + beta2[4*di - dk + ijk] + 
       beta2[-4*di + dk + ijk] - beta2[4*di + dk + ijk] + 
       beta2[-di + 4*dk + ijk] - beta2[di + 4*dk + ijk]) + 
    15000.*(beta2[-2*di - 4*dk + ijk] + beta2[-4*di - 2*dk + ijk] + 
       beta2[4*di + 2*dk + ijk] + beta2[2*di + 4*dk + ijk]) + 
    3750.*(-beta2[-3*di - 4*dk + ijk] + beta2[3*di - 4*dk + ijk] - 
       beta2[-4*di - 3*dk + ijk] + beta2[4*di - 3*dk + ijk] + 
       beta2[-4*di + 3*dk + ijk] - beta2[4*di + 3*dk + ijk] + 
       beta2[-3*di + 4*dk + ijk] - beta2[3*di + 4*dk + ijk]) - 
    300.*(beta2[3*di - 5*dk + ijk] + beta2[5*di - 3*dk + ijk] + 
       beta2[-5*di + 3*dk + ijk] + beta2[-3*di + 5*dk + ijk]) - 
    4200.*(beta2[di - 5*dk + ijk] + beta2[5*di - dk + ijk] + 
       beta2[-5*di + dk + ijk] + beta2[-di + 5*dk + ijk]) + 
    4200.*(beta2[-di - 5*dk + ijk] + beta2[-5*di - dk + ijk] + 
       beta2[5*di + dk + ijk] + beta2[di + 5*dk + ijk]) + 
    1200.*(-beta2[-2*di - 5*dk + ijk] + beta2[2*di - 5*dk + ijk] - 
       beta2[-5*di - 2*dk + ijk] + beta2[5*di - 2*dk + ijk] + 
       beta2[-5*di + 2*dk + ijk] - beta2[5*di + 2*dk + ijk] + 
       beta2[-2*di + 5*dk + ijk] - beta2[2*di + 5*dk + ijk]) + 
    300.*(beta2[-3*di - 5*dk + ijk] + beta2[-5*di - 3*dk + ijk] + 
       beta2[5*di + 3*dk + ijk] + beta2[3*di + 5*dk + ijk]) + 
    50.*(-beta2[-4*di - 5*dk + ijk] + beta2[4*di - 5*dk + ijk] - 
       beta2[-5*di - 4*dk + ijk] + beta2[5*di - 4*dk + ijk] + 
       beta2[-5*di + 4*dk + ijk] - beta2[5*di + 4*dk + ijk] + 
       beta2[-4*di + 5*dk + ijk] - beta2[4*di + 5*dk + ijk]) - 
    360000.*(beta2[2*(di - dk) + ijk] + beta2[2*(-di + dk) + ijk]) - 
    22500.*(beta2[3*(di - dk) + ijk] + beta2[3*(-di + dk) + ijk]) - 
    625.*(beta2[4*(di - dk) + ijk] + beta2[4*(-di + dk) + ijk]) - 
    4.*(beta2[5*(di - dk) + ijk] + beta2[5*(-di + dk) + ijk]) + 
    360000.*(beta2[-2*(di + dk) + ijk] + beta2[2*(di + dk) + ijk]) + 
    22500.*(beta2[-3*(di + dk) + ijk] + beta2[3*(di + dk) + ijk]) + 
    625.*(beta2[-4*(di + dk) + ijk] + beta2[4*(di + dk) + ijk]) + 
    4.*(beta2[-5*(di + dk) + ijk] + beta2[5*(di + dk) + ijk]))
;

ddb133
=
6.2988158226253464349e-7*oo2dx*oo2dz*
  (4.41e6*(beta3[-di - dk + ijk] - beta3[di - dk + ijk] - 
       beta3[-di + dk + ijk] + beta3[di + dk + ijk]) + 
    1.26e6*(-beta3[-di - 2*dk + ijk] + beta3[di - 2*dk + ijk] - 
       beta3[-2*di - dk + ijk] + beta3[2*di - dk + ijk] + 
       beta3[-2*di + dk + ijk] - beta3[2*di + dk + ijk] + 
       beta3[-di + 2*dk + ijk] - beta3[di + 2*dk + ijk]) - 
    315000.*(beta3[di - 3*dk + ijk] + beta3[3*di - dk + ijk] + 
       beta3[-3*di + dk + ijk] + beta3[-di + 3*dk + ijk]) + 
    315000.*(beta3[-di - 3*dk + ijk] + beta3[-3*di - dk + ijk] + 
       beta3[3*di + dk + ijk] + beta3[di + 3*dk + ijk]) + 
    90000.*(-beta3[-2*di - 3*dk + ijk] + beta3[2*di - 3*dk + ijk] - 
       beta3[-3*di - 2*dk + ijk] + beta3[3*di - 2*dk + ijk] + 
       beta3[-3*di + 2*dk + ijk] - beta3[3*di + 2*dk + ijk] + 
       beta3[-2*di + 3*dk + ijk] - beta3[2*di + 3*dk + ijk]) - 
    15000.*(beta3[2*di - 4*dk + ijk] + beta3[4*di - 2*dk + ijk] + 
       beta3[-4*di + 2*dk + ijk] + beta3[-2*di + 4*dk + ijk]) + 
    52500.*(-beta3[-di - 4*dk + ijk] + beta3[di - 4*dk + ijk] - 
       beta3[-4*di - dk + ijk] + beta3[4*di - dk + ijk] + 
       beta3[-4*di + dk + ijk] - beta3[4*di + dk + ijk] + 
       beta3[-di + 4*dk + ijk] - beta3[di + 4*dk + ijk]) + 
    15000.*(beta3[-2*di - 4*dk + ijk] + beta3[-4*di - 2*dk + ijk] + 
       beta3[4*di + 2*dk + ijk] + beta3[2*di + 4*dk + ijk]) + 
    3750.*(-beta3[-3*di - 4*dk + ijk] + beta3[3*di - 4*dk + ijk] - 
       beta3[-4*di - 3*dk + ijk] + beta3[4*di - 3*dk + ijk] + 
       beta3[-4*di + 3*dk + ijk] - beta3[4*di + 3*dk + ijk] + 
       beta3[-3*di + 4*dk + ijk] - beta3[3*di + 4*dk + ijk]) - 
    300.*(beta3[3*di - 5*dk + ijk] + beta3[5*di - 3*dk + ijk] + 
       beta3[-5*di + 3*dk + ijk] + beta3[-3*di + 5*dk + ijk]) - 
    4200.*(beta3[di - 5*dk + ijk] + beta3[5*di - dk + ijk] + 
       beta3[-5*di + dk + ijk] + beta3[-di + 5*dk + ijk]) + 
    4200.*(beta3[-di - 5*dk + ijk] + beta3[-5*di - dk + ijk] + 
       beta3[5*di + dk + ijk] + beta3[di + 5*dk + ijk]) + 
    1200.*(-beta3[-2*di - 5*dk + ijk] + beta3[2*di - 5*dk + ijk] - 
       beta3[-5*di - 2*dk + ijk] + beta3[5*di - 2*dk + ijk] + 
       beta3[-5*di + 2*dk + ijk] - beta3[5*di + 2*dk + ijk] + 
       beta3[-2*di + 5*dk + ijk] - beta3[2*di + 5*dk + ijk]) + 
    300.*(beta3[-3*di - 5*dk + ijk] + beta3[-5*di - 3*dk + ijk] + 
       beta3[5*di + 3*dk + ijk] + beta3[3*di + 5*dk + ijk]) + 
    50.*(-beta3[-4*di - 5*dk + ijk] + beta3[4*di - 5*dk + ijk] - 
       beta3[-5*di - 4*dk + ijk] + beta3[5*di - 4*dk + ijk] + 
       beta3[-5*di + 4*dk + ijk] - beta3[5*di + 4*dk + ijk] + 
       beta3[-4*di + 5*dk + ijk] - beta3[4*di + 5*dk + ijk]) - 
    360000.*(beta3[2*(di - dk) + ijk] + beta3[2*(-di + dk) + ijk]) - 
    22500.*(beta3[3*(di - dk) + ijk] + beta3[3*(-di + dk) + ijk]) - 
    625.*(beta3[4*(di - dk) + ijk] + beta3[4*(-di + dk) + ijk]) - 
    4.*(beta3[5*(di - dk) + ijk] + beta3[5*(-di + dk) + ijk]) + 
    360000.*(beta3[-2*(di + dk) + ijk] + beta3[2*(di + dk) + ijk]) + 
    22500.*(beta3[-3*(di + dk) + ijk] + beta3[3*(di + dk) + ijk]) + 
    625.*(beta3[-4*(di + dk) + ijk] + beta3[4*(di + dk) + ijk]) + 
    4.*(beta3[-5*(di + dk) + ijk] + beta3[5*(di + dk) + ijk]))
;

ddb221
=
0.00003968253968253968254*oody2*
  (-73766.*beta1[ijk] + 42000.*(beta1[-dj + ijk] + beta1[dj + ijk]) - 
    6000.*(beta1[-2*dj + ijk] + beta1[2*dj + ijk]) + 
    1000.*(beta1[-3*dj + ijk] + beta1[3*dj + ijk]) - 
    125.*(beta1[-4*dj + ijk] + beta1[4*dj + ijk]) + 
    8.*(beta1[-5*dj + ijk] + beta1[5*dj + ijk]))
;

ddb222
=
0.00003968253968253968254*oody2*
  (-73766.*beta2[ijk] + 42000.*(beta2[-dj + ijk] + beta2[dj + ijk]) - 
    6000.*(beta2[-2*dj + ijk] + beta2[2*dj + ijk]) + 
    1000.*(beta2[-3*dj + ijk] + beta2[3*dj + ijk]) - 
    125.*(beta2[-4*dj + ijk] + beta2[4*dj + ijk]) + 
    8.*(beta2[-5*dj + ijk] + beta2[5*dj + ijk]))
;

ddb223
=
0.00003968253968253968254*oody2*
  (-73766.*beta3[ijk] + 42000.*(beta3[-dj + ijk] + beta3[dj + ijk]) - 
    6000.*(beta3[-2*dj + ijk] + beta3[2*dj + ijk]) + 
    1000.*(beta3[-3*dj + ijk] + beta3[3*dj + ijk]) - 
    125.*(beta3[-4*dj + ijk] + beta3[4*dj + ijk]) + 
    8.*(beta3[-5*dj + ijk] + beta3[5*dj + ijk]))
;

ddb231
=
6.2988158226253464349e-7*oo2dy*oo2dz*
  (4.41e6*(beta1[-dj - dk + ijk] - beta1[dj - dk + ijk] - 
       beta1[-dj + dk + ijk] + beta1[dj + dk + ijk]) + 
    1.26e6*(-beta1[-dj - 2*dk + ijk] + beta1[dj - 2*dk + ijk] - 
       beta1[-2*dj - dk + ijk] + beta1[2*dj - dk + ijk] + 
       beta1[-2*dj + dk + ijk] - beta1[2*dj + dk + ijk] + 
       beta1[-dj + 2*dk + ijk] - beta1[dj + 2*dk + ijk]) - 
    315000.*(beta1[dj - 3*dk + ijk] + beta1[3*dj - dk + ijk] + 
       beta1[-3*dj + dk + ijk] + beta1[-dj + 3*dk + ijk]) + 
    315000.*(beta1[-dj - 3*dk + ijk] + beta1[-3*dj - dk + ijk] + 
       beta1[3*dj + dk + ijk] + beta1[dj + 3*dk + ijk]) + 
    90000.*(-beta1[-2*dj - 3*dk + ijk] + beta1[2*dj - 3*dk + ijk] - 
       beta1[-3*dj - 2*dk + ijk] + beta1[3*dj - 2*dk + ijk] + 
       beta1[-3*dj + 2*dk + ijk] - beta1[3*dj + 2*dk + ijk] + 
       beta1[-2*dj + 3*dk + ijk] - beta1[2*dj + 3*dk + ijk]) - 
    15000.*(beta1[2*dj - 4*dk + ijk] + beta1[4*dj - 2*dk + ijk] + 
       beta1[-4*dj + 2*dk + ijk] + beta1[-2*dj + 4*dk + ijk]) + 
    52500.*(-beta1[-dj - 4*dk + ijk] + beta1[dj - 4*dk + ijk] - 
       beta1[-4*dj - dk + ijk] + beta1[4*dj - dk + ijk] + 
       beta1[-4*dj + dk + ijk] - beta1[4*dj + dk + ijk] + 
       beta1[-dj + 4*dk + ijk] - beta1[dj + 4*dk + ijk]) + 
    15000.*(beta1[-2*dj - 4*dk + ijk] + beta1[-4*dj - 2*dk + ijk] + 
       beta1[4*dj + 2*dk + ijk] + beta1[2*dj + 4*dk + ijk]) + 
    3750.*(-beta1[-3*dj - 4*dk + ijk] + beta1[3*dj - 4*dk + ijk] - 
       beta1[-4*dj - 3*dk + ijk] + beta1[4*dj - 3*dk + ijk] + 
       beta1[-4*dj + 3*dk + ijk] - beta1[4*dj + 3*dk + ijk] + 
       beta1[-3*dj + 4*dk + ijk] - beta1[3*dj + 4*dk + ijk]) - 
    300.*(beta1[3*dj - 5*dk + ijk] + beta1[5*dj - 3*dk + ijk] + 
       beta1[-5*dj + 3*dk + ijk] + beta1[-3*dj + 5*dk + ijk]) - 
    4200.*(beta1[dj - 5*dk + ijk] + beta1[5*dj - dk + ijk] + 
       beta1[-5*dj + dk + ijk] + beta1[-dj + 5*dk + ijk]) + 
    4200.*(beta1[-dj - 5*dk + ijk] + beta1[-5*dj - dk + ijk] + 
       beta1[5*dj + dk + ijk] + beta1[dj + 5*dk + ijk]) + 
    1200.*(-beta1[-2*dj - 5*dk + ijk] + beta1[2*dj - 5*dk + ijk] - 
       beta1[-5*dj - 2*dk + ijk] + beta1[5*dj - 2*dk + ijk] + 
       beta1[-5*dj + 2*dk + ijk] - beta1[5*dj + 2*dk + ijk] + 
       beta1[-2*dj + 5*dk + ijk] - beta1[2*dj + 5*dk + ijk]) + 
    300.*(beta1[-3*dj - 5*dk + ijk] + beta1[-5*dj - 3*dk + ijk] + 
       beta1[5*dj + 3*dk + ijk] + beta1[3*dj + 5*dk + ijk]) + 
    50.*(-beta1[-4*dj - 5*dk + ijk] + beta1[4*dj - 5*dk + ijk] - 
       beta1[-5*dj - 4*dk + ijk] + beta1[5*dj - 4*dk + ijk] + 
       beta1[-5*dj + 4*dk + ijk] - beta1[5*dj + 4*dk + ijk] + 
       beta1[-4*dj + 5*dk + ijk] - beta1[4*dj + 5*dk + ijk]) - 
    360000.*(beta1[2*(dj - dk) + ijk] + beta1[2*(-dj + dk) + ijk]) - 
    22500.*(beta1[3*(dj - dk) + ijk] + beta1[3*(-dj + dk) + ijk]) - 
    625.*(beta1[4*(dj - dk) + ijk] + beta1[4*(-dj + dk) + ijk]) - 
    4.*(beta1[5*(dj - dk) + ijk] + beta1[5*(-dj + dk) + ijk]) + 
    360000.*(beta1[-2*(dj + dk) + ijk] + beta1[2*(dj + dk) + ijk]) + 
    22500.*(beta1[-3*(dj + dk) + ijk] + beta1[3*(dj + dk) + ijk]) + 
    625.*(beta1[-4*(dj + dk) + ijk] + beta1[4*(dj + dk) + ijk]) + 
    4.*(beta1[-5*(dj + dk) + ijk] + beta1[5*(dj + dk) + ijk]))
;

ddb232
=
6.2988158226253464349e-7*oo2dy*oo2dz*
  (4.41e6*(beta2[-dj - dk + ijk] - beta2[dj - dk + ijk] - 
       beta2[-dj + dk + ijk] + beta2[dj + dk + ijk]) + 
    1.26e6*(-beta2[-dj - 2*dk + ijk] + beta2[dj - 2*dk + ijk] - 
       beta2[-2*dj - dk + ijk] + beta2[2*dj - dk + ijk] + 
       beta2[-2*dj + dk + ijk] - beta2[2*dj + dk + ijk] + 
       beta2[-dj + 2*dk + ijk] - beta2[dj + 2*dk + ijk]) - 
    315000.*(beta2[dj - 3*dk + ijk] + beta2[3*dj - dk + ijk] + 
       beta2[-3*dj + dk + ijk] + beta2[-dj + 3*dk + ijk]) + 
    315000.*(beta2[-dj - 3*dk + ijk] + beta2[-3*dj - dk + ijk] + 
       beta2[3*dj + dk + ijk] + beta2[dj + 3*dk + ijk]) + 
    90000.*(-beta2[-2*dj - 3*dk + ijk] + beta2[2*dj - 3*dk + ijk] - 
       beta2[-3*dj - 2*dk + ijk] + beta2[3*dj - 2*dk + ijk] + 
       beta2[-3*dj + 2*dk + ijk] - beta2[3*dj + 2*dk + ijk] + 
       beta2[-2*dj + 3*dk + ijk] - beta2[2*dj + 3*dk + ijk]) - 
    15000.*(beta2[2*dj - 4*dk + ijk] + beta2[4*dj - 2*dk + ijk] + 
       beta2[-4*dj + 2*dk + ijk] + beta2[-2*dj + 4*dk + ijk]) + 
    52500.*(-beta2[-dj - 4*dk + ijk] + beta2[dj - 4*dk + ijk] - 
       beta2[-4*dj - dk + ijk] + beta2[4*dj - dk + ijk] + 
       beta2[-4*dj + dk + ijk] - beta2[4*dj + dk + ijk] + 
       beta2[-dj + 4*dk + ijk] - beta2[dj + 4*dk + ijk]) + 
    15000.*(beta2[-2*dj - 4*dk + ijk] + beta2[-4*dj - 2*dk + ijk] + 
       beta2[4*dj + 2*dk + ijk] + beta2[2*dj + 4*dk + ijk]) + 
    3750.*(-beta2[-3*dj - 4*dk + ijk] + beta2[3*dj - 4*dk + ijk] - 
       beta2[-4*dj - 3*dk + ijk] + beta2[4*dj - 3*dk + ijk] + 
       beta2[-4*dj + 3*dk + ijk] - beta2[4*dj + 3*dk + ijk] + 
       beta2[-3*dj + 4*dk + ijk] - beta2[3*dj + 4*dk + ijk]) - 
    300.*(beta2[3*dj - 5*dk + ijk] + beta2[5*dj - 3*dk + ijk] + 
       beta2[-5*dj + 3*dk + ijk] + beta2[-3*dj + 5*dk + ijk]) - 
    4200.*(beta2[dj - 5*dk + ijk] + beta2[5*dj - dk + ijk] + 
       beta2[-5*dj + dk + ijk] + beta2[-dj + 5*dk + ijk]) + 
    4200.*(beta2[-dj - 5*dk + ijk] + beta2[-5*dj - dk + ijk] + 
       beta2[5*dj + dk + ijk] + beta2[dj + 5*dk + ijk]) + 
    1200.*(-beta2[-2*dj - 5*dk + ijk] + beta2[2*dj - 5*dk + ijk] - 
       beta2[-5*dj - 2*dk + ijk] + beta2[5*dj - 2*dk + ijk] + 
       beta2[-5*dj + 2*dk + ijk] - beta2[5*dj + 2*dk + ijk] + 
       beta2[-2*dj + 5*dk + ijk] - beta2[2*dj + 5*dk + ijk]) + 
    300.*(beta2[-3*dj - 5*dk + ijk] + beta2[-5*dj - 3*dk + ijk] + 
       beta2[5*dj + 3*dk + ijk] + beta2[3*dj + 5*dk + ijk]) + 
    50.*(-beta2[-4*dj - 5*dk + ijk] + beta2[4*dj - 5*dk + ijk] - 
       beta2[-5*dj - 4*dk + ijk] + beta2[5*dj - 4*dk + ijk] + 
       beta2[-5*dj + 4*dk + ijk] - beta2[5*dj + 4*dk + ijk] + 
       beta2[-4*dj + 5*dk + ijk] - beta2[4*dj + 5*dk + ijk]) - 
    360000.*(beta2[2*(dj - dk) + ijk] + beta2[2*(-dj + dk) + ijk]) - 
    22500.*(beta2[3*(dj - dk) + ijk] + beta2[3*(-dj + dk) + ijk]) - 
    625.*(beta2[4*(dj - dk) + ijk] + beta2[4*(-dj + dk) + ijk]) - 
    4.*(beta2[5*(dj - dk) + ijk] + beta2[5*(-dj + dk) + ijk]) + 
    360000.*(beta2[-2*(dj + dk) + ijk] + beta2[2*(dj + dk) + ijk]) + 
    22500.*(beta2[-3*(dj + dk) + ijk] + beta2[3*(dj + dk) + ijk]) + 
    625.*(beta2[-4*(dj + dk) + ijk] + beta2[4*(dj + dk) + ijk]) + 
    4.*(beta2[-5*(dj + dk) + ijk] + beta2[5*(dj + dk) + ijk]))
;

ddb233
=
6.2988158226253464349e-7*oo2dy*oo2dz*
  (4.41e6*(beta3[-dj - dk + ijk] - beta3[dj - dk + ijk] - 
       beta3[-dj + dk + ijk] + beta3[dj + dk + ijk]) + 
    1.26e6*(-beta3[-dj - 2*dk + ijk] + beta3[dj - 2*dk + ijk] - 
       beta3[-2*dj - dk + ijk] + beta3[2*dj - dk + ijk] + 
       beta3[-2*dj + dk + ijk] - beta3[2*dj + dk + ijk] + 
       beta3[-dj + 2*dk + ijk] - beta3[dj + 2*dk + ijk]) - 
    315000.*(beta3[dj - 3*dk + ijk] + beta3[3*dj - dk + ijk] + 
       beta3[-3*dj + dk + ijk] + beta3[-dj + 3*dk + ijk]) + 
    315000.*(beta3[-dj - 3*dk + ijk] + beta3[-3*dj - dk + ijk] + 
       beta3[3*dj + dk + ijk] + beta3[dj + 3*dk + ijk]) + 
    90000.*(-beta3[-2*dj - 3*dk + ijk] + beta3[2*dj - 3*dk + ijk] - 
       beta3[-3*dj - 2*dk + ijk] + beta3[3*dj - 2*dk + ijk] + 
       beta3[-3*dj + 2*dk + ijk] - beta3[3*dj + 2*dk + ijk] + 
       beta3[-2*dj + 3*dk + ijk] - beta3[2*dj + 3*dk + ijk]) - 
    15000.*(beta3[2*dj - 4*dk + ijk] + beta3[4*dj - 2*dk + ijk] + 
       beta3[-4*dj + 2*dk + ijk] + beta3[-2*dj + 4*dk + ijk]) + 
    52500.*(-beta3[-dj - 4*dk + ijk] + beta3[dj - 4*dk + ijk] - 
       beta3[-4*dj - dk + ijk] + beta3[4*dj - dk + ijk] + 
       beta3[-4*dj + dk + ijk] - beta3[4*dj + dk + ijk] + 
       beta3[-dj + 4*dk + ijk] - beta3[dj + 4*dk + ijk]) + 
    15000.*(beta3[-2*dj - 4*dk + ijk] + beta3[-4*dj - 2*dk + ijk] + 
       beta3[4*dj + 2*dk + ijk] + beta3[2*dj + 4*dk + ijk]) + 
    3750.*(-beta3[-3*dj - 4*dk + ijk] + beta3[3*dj - 4*dk + ijk] - 
       beta3[-4*dj - 3*dk + ijk] + beta3[4*dj - 3*dk + ijk] + 
       beta3[-4*dj + 3*dk + ijk] - beta3[4*dj + 3*dk + ijk] + 
       beta3[-3*dj + 4*dk + ijk] - beta3[3*dj + 4*dk + ijk]) - 
    300.*(beta3[3*dj - 5*dk + ijk] + beta3[5*dj - 3*dk + ijk] + 
       beta3[-5*dj + 3*dk + ijk] + beta3[-3*dj + 5*dk + ijk]) - 
    4200.*(beta3[dj - 5*dk + ijk] + beta3[5*dj - dk + ijk] + 
       beta3[-5*dj + dk + ijk] + beta3[-dj + 5*dk + ijk]) + 
    4200.*(beta3[-dj - 5*dk + ijk] + beta3[-5*dj - dk + ijk] + 
       beta3[5*dj + dk + ijk] + beta3[dj + 5*dk + ijk]) + 
    1200.*(-beta3[-2*dj - 5*dk + ijk] + beta3[2*dj - 5*dk + ijk] - 
       beta3[-5*dj - 2*dk + ijk] + beta3[5*dj - 2*dk + ijk] + 
       beta3[-5*dj + 2*dk + ijk] - beta3[5*dj + 2*dk + ijk] + 
       beta3[-2*dj + 5*dk + ijk] - beta3[2*dj + 5*dk + ijk]) + 
    300.*(beta3[-3*dj - 5*dk + ijk] + beta3[-5*dj - 3*dk + ijk] + 
       beta3[5*dj + 3*dk + ijk] + beta3[3*dj + 5*dk + ijk]) + 
    50.*(-beta3[-4*dj - 5*dk + ijk] + beta3[4*dj - 5*dk + ijk] - 
       beta3[-5*dj - 4*dk + ijk] + beta3[5*dj - 4*dk + ijk] + 
       beta3[-5*dj + 4*dk + ijk] - beta3[5*dj + 4*dk + ijk] + 
       beta3[-4*dj + 5*dk + ijk] - beta3[4*dj + 5*dk + ijk]) - 
    360000.*(beta3[2*(dj - dk) + ijk] + beta3[2*(-dj + dk) + ijk]) - 
    22500.*(beta3[3*(dj - dk) + ijk] + beta3[3*(-dj + dk) + ijk]) - 
    625.*(beta3[4*(dj - dk) + ijk] + beta3[4*(-dj + dk) + ijk]) - 
    4.*(beta3[5*(dj - dk) + ijk] + beta3[5*(-dj + dk) + ijk]) + 
    360000.*(beta3[-2*(dj + dk) + ijk] + beta3[2*(dj + dk) + ijk]) + 
    22500.*(beta3[-3*(dj + dk) + ijk] + beta3[3*(dj + dk) + ijk]) + 
    625.*(beta3[-4*(dj + dk) + ijk] + beta3[4*(dj + dk) + ijk]) + 
    4.*(beta3[-5*(dj + dk) + ijk] + beta3[5*(dj + dk) + ijk]))
;

ddb331
=
0.00003968253968253968254*oodz2*
  (-73766.*beta1[ijk] + 42000.*(beta1[-dk + ijk] + beta1[dk + ijk]) - 
    6000.*(beta1[-2*dk + ijk] + beta1[2*dk + ijk]) + 
    1000.*(beta1[-3*dk + ijk] + beta1[3*dk + ijk]) - 
    125.*(beta1[-4*dk + ijk] + beta1[4*dk + ijk]) + 
    8.*(beta1[-5*dk + ijk] + beta1[5*dk + ijk]))
;

ddb332
=
0.00003968253968253968254*oodz2*
  (-73766.*beta2[ijk] + 42000.*(beta2[-dk + ijk] + beta2[dk + ijk]) - 
    6000.*(beta2[-2*dk + ijk] + beta2[2*dk + ijk]) + 
    1000.*(beta2[-3*dk + ijk] + beta2[3*dk + ijk]) - 
    125.*(beta2[-4*dk + ijk] + beta2[4*dk + ijk]) + 
    8.*(beta2[-5*dk + ijk] + beta2[5*dk + ijk]))
;

ddb333
=
0.00003968253968253968254*oodz2*
  (-73766.*beta3[ijk] + 42000.*(beta3[-dk + ijk] + beta3[dk + ijk]) - 
    6000.*(beta3[-2*dk + ijk] + beta3[2*dk + ijk]) + 
    1000.*(beta3[-3*dk + ijk] + beta3[3*dk + ijk]) - 
    125.*(beta3[-4*dk + ijk] + beta3[4*dk + ijk]) + 
    8.*(beta3[-5*dk + ijk] + beta3[5*dk + ijk]))
;

delg111
=
0.00079365079365079365079*oo2dx*
  (-2.*g11[-5*di + ijk] + 25.*g11[-4*di + ijk] + 600.*g11[-2*di + ijk] + 
    2100.*(-g11[-di + ijk] + g11[di + ijk]) - 600.*g11[2*di + ijk] + 
    150.*(-g11[-3*di + ijk] + g11[3*di + ijk]) - 25.*g11[4*di + ijk] + 
    2.*g11[5*di + ijk])
;

delg112
=
0.00079365079365079365079*oo2dx*
  (-2.*g12[-5*di + ijk] + 25.*g12[-4*di + ijk] + 600.*g12[-2*di + ijk] + 
    2100.*(-g12[-di + ijk] + g12[di + ijk]) - 600.*g12[2*di + ijk] + 
    150.*(-g12[-3*di + ijk] + g12[3*di + ijk]) - 25.*g12[4*di + ijk] + 
    2.*g12[5*di + ijk])
;

delg113
=
0.00079365079365079365079*oo2dx*
  (-2.*g13[-5*di + ijk] + 25.*g13[-4*di + ijk] + 600.*g13[-2*di + ijk] + 
    2100.*(-g13[-di + ijk] + g13[di + ijk]) - 600.*g13[2*di + ijk] + 
    150.*(-g13[-3*di + ijk] + g13[3*di + ijk]) - 25.*g13[4*di + ijk] + 
    2.*g13[5*di + ijk])
;

delg122
=
0.00079365079365079365079*oo2dx*
  (-2.*g22[-5*di + ijk] + 25.*g22[-4*di + ijk] + 600.*g22[-2*di + ijk] + 
    2100.*(-g22[-di + ijk] + g22[di + ijk]) - 600.*g22[2*di + ijk] + 
    150.*(-g22[-3*di + ijk] + g22[3*di + ijk]) - 25.*g22[4*di + ijk] + 
    2.*g22[5*di + ijk])
;

delg123
=
0.00079365079365079365079*oo2dx*
  (-2.*g23[-5*di + ijk] + 25.*g23[-4*di + ijk] + 600.*g23[-2*di + ijk] + 
    2100.*(-g23[-di + ijk] + g23[di + ijk]) - 600.*g23[2*di + ijk] + 
    150.*(-g23[-3*di + ijk] + g23[3*di + ijk]) - 25.*g23[4*di + ijk] + 
    2.*g23[5*di + ijk])
;

delg133
=
0.00079365079365079365079*oo2dx*
  (-2.*g33[-5*di + ijk] + 25.*g33[-4*di + ijk] + 600.*g33[-2*di + ijk] + 
    2100.*(-g33[-di + ijk] + g33[di + ijk]) - 600.*g33[2*di + ijk] + 
    150.*(-g33[-3*di + ijk] + g33[3*di + ijk]) - 25.*g33[4*di + ijk] + 
    2.*g33[5*di + ijk])
;

delg211
=
0.00079365079365079365079*oo2dy*
  (-2.*g11[-5*dj + ijk] + 25.*g11[-4*dj + ijk] + 600.*g11[-2*dj + ijk] + 
    2100.*(-g11[-dj + ijk] + g11[dj + ijk]) - 600.*g11[2*dj + ijk] + 
    150.*(-g11[-3*dj + ijk] + g11[3*dj + ijk]) - 25.*g11[4*dj + ijk] + 
    2.*g11[5*dj + ijk])
;

delg212
=
0.00079365079365079365079*oo2dy*
  (-2.*g12[-5*dj + ijk] + 25.*g12[-4*dj + ijk] + 600.*g12[-2*dj + ijk] + 
    2100.*(-g12[-dj + ijk] + g12[dj + ijk]) - 600.*g12[2*dj + ijk] + 
    150.*(-g12[-3*dj + ijk] + g12[3*dj + ijk]) - 25.*g12[4*dj + ijk] + 
    2.*g12[5*dj + ijk])
;

delg213
=
0.00079365079365079365079*oo2dy*
  (-2.*g13[-5*dj + ijk] + 25.*g13[-4*dj + ijk] + 600.*g13[-2*dj + ijk] + 
    2100.*(-g13[-dj + ijk] + g13[dj + ijk]) - 600.*g13[2*dj + ijk] + 
    150.*(-g13[-3*dj + ijk] + g13[3*dj + ijk]) - 25.*g13[4*dj + ijk] + 
    2.*g13[5*dj + ijk])
;

delg222
=
0.00079365079365079365079*oo2dy*
  (-2.*g22[-5*dj + ijk] + 25.*g22[-4*dj + ijk] + 600.*g22[-2*dj + ijk] + 
    2100.*(-g22[-dj + ijk] + g22[dj + ijk]) - 600.*g22[2*dj + ijk] + 
    150.*(-g22[-3*dj + ijk] + g22[3*dj + ijk]) - 25.*g22[4*dj + ijk] + 
    2.*g22[5*dj + ijk])
;

delg223
=
0.00079365079365079365079*oo2dy*
  (-2.*g23[-5*dj + ijk] + 25.*g23[-4*dj + ijk] + 600.*g23[-2*dj + ijk] + 
    2100.*(-g23[-dj + ijk] + g23[dj + ijk]) - 600.*g23[2*dj + ijk] + 
    150.*(-g23[-3*dj + ijk] + g23[3*dj + ijk]) - 25.*g23[4*dj + ijk] + 
    2.*g23[5*dj + ijk])
;

delg233
=
0.00079365079365079365079*oo2dy*
  (-2.*g33[-5*dj + ijk] + 25.*g33[-4*dj + ijk] + 600.*g33[-2*dj + ijk] + 
    2100.*(-g33[-dj + ijk] + g33[dj + ijk]) - 600.*g33[2*dj + ijk] + 
    150.*(-g33[-3*dj + ijk] + g33[3*dj + ijk]) - 25.*g33[4*dj + ijk] + 
    2.*g33[5*dj + ijk])
;

delg311
=
0.00079365079365079365079*oo2dz*
  (-2.*g11[-5*dk + ijk] + 25.*g11[-4*dk + ijk] + 600.*g11[-2*dk + ijk] + 
    2100.*(-g11[-dk + ijk] + g11[dk + ijk]) - 600.*g11[2*dk + ijk] + 
    150.*(-g11[-3*dk + ijk] + g11[3*dk + ijk]) - 25.*g11[4*dk + ijk] + 
    2.*g11[5*dk + ijk])
;

delg312
=
0.00079365079365079365079*oo2dz*
  (-2.*g12[-5*dk + ijk] + 25.*g12[-4*dk + ijk] + 600.*g12[-2*dk + ijk] + 
    2100.*(-g12[-dk + ijk] + g12[dk + ijk]) - 600.*g12[2*dk + ijk] + 
    150.*(-g12[-3*dk + ijk] + g12[3*dk + ijk]) - 25.*g12[4*dk + ijk] + 
    2.*g12[5*dk + ijk])
;

delg313
=
0.00079365079365079365079*oo2dz*
  (-2.*g13[-5*dk + ijk] + 25.*g13[-4*dk + ijk] + 600.*g13[-2*dk + ijk] + 
    2100.*(-g13[-dk + ijk] + g13[dk + ijk]) - 600.*g13[2*dk + ijk] + 
    150.*(-g13[-3*dk + ijk] + g13[3*dk + ijk]) - 25.*g13[4*dk + ijk] + 
    2.*g13[5*dk + ijk])
;

delg322
=
0.00079365079365079365079*oo2dz*
  (-2.*g22[-5*dk + ijk] + 25.*g22[-4*dk + ijk] + 600.*g22[-2*dk + ijk] + 
    2100.*(-g22[-dk + ijk] + g22[dk + ijk]) - 600.*g22[2*dk + ijk] + 
    150.*(-g22[-3*dk + ijk] + g22[3*dk + ijk]) - 25.*g22[4*dk + ijk] + 
    2.*g22[5*dk + ijk])
;

delg323
=
0.00079365079365079365079*oo2dz*
  (-2.*g23[-5*dk + ijk] + 25.*g23[-4*dk + ijk] + 600.*g23[-2*dk + ijk] + 
    2100.*(-g23[-dk + ijk] + g23[dk + ijk]) - 600.*g23[2*dk + ijk] + 
    150.*(-g23[-3*dk + ijk] + g23[3*dk + ijk]) - 25.*g23[4*dk + ijk] + 
    2.*g23[5*dk + ijk])
;

delg333
=
0.00079365079365079365079*oo2dz*
  (-2.*g33[-5*dk + ijk] + 25.*g33[-4*dk + ijk] + 600.*g33[-2*dk + ijk] + 
    2100.*(-g33[-dk + ijk] + g33[dk + ijk]) - 600.*g33[2*dk + ijk] + 
    150.*(-g33[-3*dk + ijk] + g33[3*dk + ijk]) - 25.*g33[4*dk + ijk] + 
    2.*g33[5*dk + ijk])
;

deldelg1111
=
0.00003968253968253968254*oodx2*
  (-73766.*g11[ijk] + 42000.*(g11[-di + ijk] + g11[di + ijk]) - 
    6000.*(g11[-2*di + ijk] + g11[2*di + ijk]) + 
    1000.*(g11[-3*di + ijk] + g11[3*di + ijk]) - 
    125.*(g11[-4*di + ijk] + g11[4*di + ijk]) + 
    8.*(g11[-5*di + ijk] + g11[5*di + ijk]))
;

deldelg1112
=
0.00003968253968253968254*oodx2*
  (-73766.*g12[ijk] + 42000.*(g12[-di + ijk] + g12[di + ijk]) - 
    6000.*(g12[-2*di + ijk] + g12[2*di + ijk]) + 
    1000.*(g12[-3*di + ijk] + g12[3*di + ijk]) - 
    125.*(g12[-4*di + ijk] + g12[4*di + ijk]) + 
    8.*(g12[-5*di + ijk] + g12[5*di + ijk]))
;

deldelg1113
=
0.00003968253968253968254*oodx2*
  (-73766.*g13[ijk] + 42000.*(g13[-di + ijk] + g13[di + ijk]) - 
    6000.*(g13[-2*di + ijk] + g13[2*di + ijk]) + 
    1000.*(g13[-3*di + ijk] + g13[3*di + ijk]) - 
    125.*(g13[-4*di + ijk] + g13[4*di + ijk]) + 
    8.*(g13[-5*di + ijk] + g13[5*di + ijk]))
;

deldelg1122
=
0.00003968253968253968254*oodx2*
  (-73766.*g22[ijk] + 42000.*(g22[-di + ijk] + g22[di + ijk]) - 
    6000.*(g22[-2*di + ijk] + g22[2*di + ijk]) + 
    1000.*(g22[-3*di + ijk] + g22[3*di + ijk]) - 
    125.*(g22[-4*di + ijk] + g22[4*di + ijk]) + 
    8.*(g22[-5*di + ijk] + g22[5*di + ijk]))
;

deldelg1123
=
0.00003968253968253968254*oodx2*
  (-73766.*g23[ijk] + 42000.*(g23[-di + ijk] + g23[di + ijk]) - 
    6000.*(g23[-2*di + ijk] + g23[2*di + ijk]) + 
    1000.*(g23[-3*di + ijk] + g23[3*di + ijk]) - 
    125.*(g23[-4*di + ijk] + g23[4*di + ijk]) + 
    8.*(g23[-5*di + ijk] + g23[5*di + ijk]))
;

deldelg1133
=
0.00003968253968253968254*oodx2*
  (-73766.*g33[ijk] + 42000.*(g33[-di + ijk] + g33[di + ijk]) - 
    6000.*(g33[-2*di + ijk] + g33[2*di + ijk]) + 
    1000.*(g33[-3*di + ijk] + g33[3*di + ijk]) - 
    125.*(g33[-4*di + ijk] + g33[4*di + ijk]) + 
    8.*(g33[-5*di + ijk] + g33[5*di + ijk]))
;

deldelg1211
=
6.2988158226253464349e-7*oo2dx*oo2dy*
  (4.41e6*(g11[-di - dj + ijk] - g11[di - dj + ijk] - g11[-di + dj + ijk] + 
       g11[di + dj + ijk]) + 1.26e6*
     (-g11[-di - 2*dj + ijk] + g11[di - 2*dj + ijk] - 
       g11[-2*di - dj + ijk] + g11[2*di - dj + ijk] + 
       g11[-2*di + dj + ijk] - g11[2*di + dj + ijk] + 
       g11[-di + 2*dj + ijk] - g11[di + 2*dj + ijk]) - 
    315000.*(g11[di - 3*dj + ijk] + g11[3*di - dj + ijk] + 
       g11[-3*di + dj + ijk] + g11[-di + 3*dj + ijk]) + 
    315000.*(g11[-di - 3*dj + ijk] + g11[-3*di - dj + ijk] + 
       g11[3*di + dj + ijk] + g11[di + 3*dj + ijk]) + 
    90000.*(-g11[-2*di - 3*dj + ijk] + g11[2*di - 3*dj + ijk] - 
       g11[-3*di - 2*dj + ijk] + g11[3*di - 2*dj + ijk] + 
       g11[-3*di + 2*dj + ijk] - g11[3*di + 2*dj + ijk] + 
       g11[-2*di + 3*dj + ijk] - g11[2*di + 3*dj + ijk]) - 
    15000.*(g11[2*di - 4*dj + ijk] + g11[4*di - 2*dj + ijk] + 
       g11[-4*di + 2*dj + ijk] + g11[-2*di + 4*dj + ijk]) + 
    52500.*(-g11[-di - 4*dj + ijk] + g11[di - 4*dj + ijk] - 
       g11[-4*di - dj + ijk] + g11[4*di - dj + ijk] + 
       g11[-4*di + dj + ijk] - g11[4*di + dj + ijk] + 
       g11[-di + 4*dj + ijk] - g11[di + 4*dj + ijk]) + 
    15000.*(g11[-2*di - 4*dj + ijk] + g11[-4*di - 2*dj + ijk] + 
       g11[4*di + 2*dj + ijk] + g11[2*di + 4*dj + ijk]) + 
    3750.*(-g11[-3*di - 4*dj + ijk] + g11[3*di - 4*dj + ijk] - 
       g11[-4*di - 3*dj + ijk] + g11[4*di - 3*dj + ijk] + 
       g11[-4*di + 3*dj + ijk] - g11[4*di + 3*dj + ijk] + 
       g11[-3*di + 4*dj + ijk] - g11[3*di + 4*dj + ijk]) - 
    300.*(g11[3*di - 5*dj + ijk] + g11[5*di - 3*dj + ijk] + 
       g11[-5*di + 3*dj + ijk] + g11[-3*di + 5*dj + ijk]) - 
    4200.*(g11[di - 5*dj + ijk] + g11[5*di - dj + ijk] + 
       g11[-5*di + dj + ijk] + g11[-di + 5*dj + ijk]) + 
    4200.*(g11[-di - 5*dj + ijk] + g11[-5*di - dj + ijk] + 
       g11[5*di + dj + ijk] + g11[di + 5*dj + ijk]) + 
    1200.*(-g11[-2*di - 5*dj + ijk] + g11[2*di - 5*dj + ijk] - 
       g11[-5*di - 2*dj + ijk] + g11[5*di - 2*dj + ijk] + 
       g11[-5*di + 2*dj + ijk] - g11[5*di + 2*dj + ijk] + 
       g11[-2*di + 5*dj + ijk] - g11[2*di + 5*dj + ijk]) + 
    300.*(g11[-3*di - 5*dj + ijk] + g11[-5*di - 3*dj + ijk] + 
       g11[5*di + 3*dj + ijk] + g11[3*di + 5*dj + ijk]) + 
    50.*(-g11[-4*di - 5*dj + ijk] + g11[4*di - 5*dj + ijk] - 
       g11[-5*di - 4*dj + ijk] + g11[5*di - 4*dj + ijk] + 
       g11[-5*di + 4*dj + ijk] - g11[5*di + 4*dj + ijk] + 
       g11[-4*di + 5*dj + ijk] - g11[4*di + 5*dj + ijk]) - 
    360000.*(g11[2*(di - dj) + ijk] + g11[2*(-di + dj) + ijk]) - 
    22500.*(g11[3*(di - dj) + ijk] + g11[3*(-di + dj) + ijk]) - 
    625.*(g11[4*(di - dj) + ijk] + g11[4*(-di + dj) + ijk]) - 
    4.*(g11[5*(di - dj) + ijk] + g11[5*(-di + dj) + ijk]) + 
    360000.*(g11[-2*(di + dj) + ijk] + g11[2*(di + dj) + ijk]) + 
    22500.*(g11[-3*(di + dj) + ijk] + g11[3*(di + dj) + ijk]) + 
    625.*(g11[-4*(di + dj) + ijk] + g11[4*(di + dj) + ijk]) + 
    4.*(g11[-5*(di + dj) + ijk] + g11[5*(di + dj) + ijk]))
;

deldelg1212
=
6.2988158226253464349e-7*oo2dx*oo2dy*
  (4.41e6*(g12[-di - dj + ijk] - g12[di - dj + ijk] - g12[-di + dj + ijk] + 
       g12[di + dj + ijk]) + 1.26e6*
     (-g12[-di - 2*dj + ijk] + g12[di - 2*dj + ijk] - 
       g12[-2*di - dj + ijk] + g12[2*di - dj + ijk] + 
       g12[-2*di + dj + ijk] - g12[2*di + dj + ijk] + 
       g12[-di + 2*dj + ijk] - g12[di + 2*dj + ijk]) - 
    315000.*(g12[di - 3*dj + ijk] + g12[3*di - dj + ijk] + 
       g12[-3*di + dj + ijk] + g12[-di + 3*dj + ijk]) + 
    315000.*(g12[-di - 3*dj + ijk] + g12[-3*di - dj + ijk] + 
       g12[3*di + dj + ijk] + g12[di + 3*dj + ijk]) + 
    90000.*(-g12[-2*di - 3*dj + ijk] + g12[2*di - 3*dj + ijk] - 
       g12[-3*di - 2*dj + ijk] + g12[3*di - 2*dj + ijk] + 
       g12[-3*di + 2*dj + ijk] - g12[3*di + 2*dj + ijk] + 
       g12[-2*di + 3*dj + ijk] - g12[2*di + 3*dj + ijk]) - 
    15000.*(g12[2*di - 4*dj + ijk] + g12[4*di - 2*dj + ijk] + 
       g12[-4*di + 2*dj + ijk] + g12[-2*di + 4*dj + ijk]) + 
    52500.*(-g12[-di - 4*dj + ijk] + g12[di - 4*dj + ijk] - 
       g12[-4*di - dj + ijk] + g12[4*di - dj + ijk] + 
       g12[-4*di + dj + ijk] - g12[4*di + dj + ijk] + 
       g12[-di + 4*dj + ijk] - g12[di + 4*dj + ijk]) + 
    15000.*(g12[-2*di - 4*dj + ijk] + g12[-4*di - 2*dj + ijk] + 
       g12[4*di + 2*dj + ijk] + g12[2*di + 4*dj + ijk]) + 
    3750.*(-g12[-3*di - 4*dj + ijk] + g12[3*di - 4*dj + ijk] - 
       g12[-4*di - 3*dj + ijk] + g12[4*di - 3*dj + ijk] + 
       g12[-4*di + 3*dj + ijk] - g12[4*di + 3*dj + ijk] + 
       g12[-3*di + 4*dj + ijk] - g12[3*di + 4*dj + ijk]) - 
    300.*(g12[3*di - 5*dj + ijk] + g12[5*di - 3*dj + ijk] + 
       g12[-5*di + 3*dj + ijk] + g12[-3*di + 5*dj + ijk]) - 
    4200.*(g12[di - 5*dj + ijk] + g12[5*di - dj + ijk] + 
       g12[-5*di + dj + ijk] + g12[-di + 5*dj + ijk]) + 
    4200.*(g12[-di - 5*dj + ijk] + g12[-5*di - dj + ijk] + 
       g12[5*di + dj + ijk] + g12[di + 5*dj + ijk]) + 
    1200.*(-g12[-2*di - 5*dj + ijk] + g12[2*di - 5*dj + ijk] - 
       g12[-5*di - 2*dj + ijk] + g12[5*di - 2*dj + ijk] + 
       g12[-5*di + 2*dj + ijk] - g12[5*di + 2*dj + ijk] + 
       g12[-2*di + 5*dj + ijk] - g12[2*di + 5*dj + ijk]) + 
    300.*(g12[-3*di - 5*dj + ijk] + g12[-5*di - 3*dj + ijk] + 
       g12[5*di + 3*dj + ijk] + g12[3*di + 5*dj + ijk]) + 
    50.*(-g12[-4*di - 5*dj + ijk] + g12[4*di - 5*dj + ijk] - 
       g12[-5*di - 4*dj + ijk] + g12[5*di - 4*dj + ijk] + 
       g12[-5*di + 4*dj + ijk] - g12[5*di + 4*dj + ijk] + 
       g12[-4*di + 5*dj + ijk] - g12[4*di + 5*dj + ijk]) - 
    360000.*(g12[2*(di - dj) + ijk] + g12[2*(-di + dj) + ijk]) - 
    22500.*(g12[3*(di - dj) + ijk] + g12[3*(-di + dj) + ijk]) - 
    625.*(g12[4*(di - dj) + ijk] + g12[4*(-di + dj) + ijk]) - 
    4.*(g12[5*(di - dj) + ijk] + g12[5*(-di + dj) + ijk]) + 
    360000.*(g12[-2*(di + dj) + ijk] + g12[2*(di + dj) + ijk]) + 
    22500.*(g12[-3*(di + dj) + ijk] + g12[3*(di + dj) + ijk]) + 
    625.*(g12[-4*(di + dj) + ijk] + g12[4*(di + dj) + ijk]) + 
    4.*(g12[-5*(di + dj) + ijk] + g12[5*(di + dj) + ijk]))
;

deldelg1213
=
6.2988158226253464349e-7*oo2dx*oo2dy*
  (4.41e6*(g13[-di - dj + ijk] - g13[di - dj + ijk] - g13[-di + dj + ijk] + 
       g13[di + dj + ijk]) + 1.26e6*
     (-g13[-di - 2*dj + ijk] + g13[di - 2*dj + ijk] - 
       g13[-2*di - dj + ijk] + g13[2*di - dj + ijk] + 
       g13[-2*di + dj + ijk] - g13[2*di + dj + ijk] + 
       g13[-di + 2*dj + ijk] - g13[di + 2*dj + ijk]) - 
    315000.*(g13[di - 3*dj + ijk] + g13[3*di - dj + ijk] + 
       g13[-3*di + dj + ijk] + g13[-di + 3*dj + ijk]) + 
    315000.*(g13[-di - 3*dj + ijk] + g13[-3*di - dj + ijk] + 
       g13[3*di + dj + ijk] + g13[di + 3*dj + ijk]) + 
    90000.*(-g13[-2*di - 3*dj + ijk] + g13[2*di - 3*dj + ijk] - 
       g13[-3*di - 2*dj + ijk] + g13[3*di - 2*dj + ijk] + 
       g13[-3*di + 2*dj + ijk] - g13[3*di + 2*dj + ijk] + 
       g13[-2*di + 3*dj + ijk] - g13[2*di + 3*dj + ijk]) - 
    15000.*(g13[2*di - 4*dj + ijk] + g13[4*di - 2*dj + ijk] + 
       g13[-4*di + 2*dj + ijk] + g13[-2*di + 4*dj + ijk]) + 
    52500.*(-g13[-di - 4*dj + ijk] + g13[di - 4*dj + ijk] - 
       g13[-4*di - dj + ijk] + g13[4*di - dj + ijk] + 
       g13[-4*di + dj + ijk] - g13[4*di + dj + ijk] + 
       g13[-di + 4*dj + ijk] - g13[di + 4*dj + ijk]) + 
    15000.*(g13[-2*di - 4*dj + ijk] + g13[-4*di - 2*dj + ijk] + 
       g13[4*di + 2*dj + ijk] + g13[2*di + 4*dj + ijk]) + 
    3750.*(-g13[-3*di - 4*dj + ijk] + g13[3*di - 4*dj + ijk] - 
       g13[-4*di - 3*dj + ijk] + g13[4*di - 3*dj + ijk] + 
       g13[-4*di + 3*dj + ijk] - g13[4*di + 3*dj + ijk] + 
       g13[-3*di + 4*dj + ijk] - g13[3*di + 4*dj + ijk]) - 
    300.*(g13[3*di - 5*dj + ijk] + g13[5*di - 3*dj + ijk] + 
       g13[-5*di + 3*dj + ijk] + g13[-3*di + 5*dj + ijk]) - 
    4200.*(g13[di - 5*dj + ijk] + g13[5*di - dj + ijk] + 
       g13[-5*di + dj + ijk] + g13[-di + 5*dj + ijk]) + 
    4200.*(g13[-di - 5*dj + ijk] + g13[-5*di - dj + ijk] + 
       g13[5*di + dj + ijk] + g13[di + 5*dj + ijk]) + 
    1200.*(-g13[-2*di - 5*dj + ijk] + g13[2*di - 5*dj + ijk] - 
       g13[-5*di - 2*dj + ijk] + g13[5*di - 2*dj + ijk] + 
       g13[-5*di + 2*dj + ijk] - g13[5*di + 2*dj + ijk] + 
       g13[-2*di + 5*dj + ijk] - g13[2*di + 5*dj + ijk]) + 
    300.*(g13[-3*di - 5*dj + ijk] + g13[-5*di - 3*dj + ijk] + 
       g13[5*di + 3*dj + ijk] + g13[3*di + 5*dj + ijk]) + 
    50.*(-g13[-4*di - 5*dj + ijk] + g13[4*di - 5*dj + ijk] - 
       g13[-5*di - 4*dj + ijk] + g13[5*di - 4*dj + ijk] + 
       g13[-5*di + 4*dj + ijk] - g13[5*di + 4*dj + ijk] + 
       g13[-4*di + 5*dj + ijk] - g13[4*di + 5*dj + ijk]) - 
    360000.*(g13[2*(di - dj) + ijk] + g13[2*(-di + dj) + ijk]) - 
    22500.*(g13[3*(di - dj) + ijk] + g13[3*(-di + dj) + ijk]) - 
    625.*(g13[4*(di - dj) + ijk] + g13[4*(-di + dj) + ijk]) - 
    4.*(g13[5*(di - dj) + ijk] + g13[5*(-di + dj) + ijk]) + 
    360000.*(g13[-2*(di + dj) + ijk] + g13[2*(di + dj) + ijk]) + 
    22500.*(g13[-3*(di + dj) + ijk] + g13[3*(di + dj) + ijk]) + 
    625.*(g13[-4*(di + dj) + ijk] + g13[4*(di + dj) + ijk]) + 
    4.*(g13[-5*(di + dj) + ijk] + g13[5*(di + dj) + ijk]))
;

deldelg1222
=
6.2988158226253464349e-7*oo2dx*oo2dy*
  (4.41e6*(g22[-di - dj + ijk] - g22[di - dj + ijk] - g22[-di + dj + ijk] + 
       g22[di + dj + ijk]) + 1.26e6*
     (-g22[-di - 2*dj + ijk] + g22[di - 2*dj + ijk] - 
       g22[-2*di - dj + ijk] + g22[2*di - dj + ijk] + 
       g22[-2*di + dj + ijk] - g22[2*di + dj + ijk] + 
       g22[-di + 2*dj + ijk] - g22[di + 2*dj + ijk]) - 
    315000.*(g22[di - 3*dj + ijk] + g22[3*di - dj + ijk] + 
       g22[-3*di + dj + ijk] + g22[-di + 3*dj + ijk]) + 
    315000.*(g22[-di - 3*dj + ijk] + g22[-3*di - dj + ijk] + 
       g22[3*di + dj + ijk] + g22[di + 3*dj + ijk]) + 
    90000.*(-g22[-2*di - 3*dj + ijk] + g22[2*di - 3*dj + ijk] - 
       g22[-3*di - 2*dj + ijk] + g22[3*di - 2*dj + ijk] + 
       g22[-3*di + 2*dj + ijk] - g22[3*di + 2*dj + ijk] + 
       g22[-2*di + 3*dj + ijk] - g22[2*di + 3*dj + ijk]) - 
    15000.*(g22[2*di - 4*dj + ijk] + g22[4*di - 2*dj + ijk] + 
       g22[-4*di + 2*dj + ijk] + g22[-2*di + 4*dj + ijk]) + 
    52500.*(-g22[-di - 4*dj + ijk] + g22[di - 4*dj + ijk] - 
       g22[-4*di - dj + ijk] + g22[4*di - dj + ijk] + 
       g22[-4*di + dj + ijk] - g22[4*di + dj + ijk] + 
       g22[-di + 4*dj + ijk] - g22[di + 4*dj + ijk]) + 
    15000.*(g22[-2*di - 4*dj + ijk] + g22[-4*di - 2*dj + ijk] + 
       g22[4*di + 2*dj + ijk] + g22[2*di + 4*dj + ijk]) + 
    3750.*(-g22[-3*di - 4*dj + ijk] + g22[3*di - 4*dj + ijk] - 
       g22[-4*di - 3*dj + ijk] + g22[4*di - 3*dj + ijk] + 
       g22[-4*di + 3*dj + ijk] - g22[4*di + 3*dj + ijk] + 
       g22[-3*di + 4*dj + ijk] - g22[3*di + 4*dj + ijk]) - 
    300.*(g22[3*di - 5*dj + ijk] + g22[5*di - 3*dj + ijk] + 
       g22[-5*di + 3*dj + ijk] + g22[-3*di + 5*dj + ijk]) - 
    4200.*(g22[di - 5*dj + ijk] + g22[5*di - dj + ijk] + 
       g22[-5*di + dj + ijk] + g22[-di + 5*dj + ijk]) + 
    4200.*(g22[-di - 5*dj + ijk] + g22[-5*di - dj + ijk] + 
       g22[5*di + dj + ijk] + g22[di + 5*dj + ijk]) + 
    1200.*(-g22[-2*di - 5*dj + ijk] + g22[2*di - 5*dj + ijk] - 
       g22[-5*di - 2*dj + ijk] + g22[5*di - 2*dj + ijk] + 
       g22[-5*di + 2*dj + ijk] - g22[5*di + 2*dj + ijk] + 
       g22[-2*di + 5*dj + ijk] - g22[2*di + 5*dj + ijk]) + 
    300.*(g22[-3*di - 5*dj + ijk] + g22[-5*di - 3*dj + ijk] + 
       g22[5*di + 3*dj + ijk] + g22[3*di + 5*dj + ijk]) + 
    50.*(-g22[-4*di - 5*dj + ijk] + g22[4*di - 5*dj + ijk] - 
       g22[-5*di - 4*dj + ijk] + g22[5*di - 4*dj + ijk] + 
       g22[-5*di + 4*dj + ijk] - g22[5*di + 4*dj + ijk] + 
       g22[-4*di + 5*dj + ijk] - g22[4*di + 5*dj + ijk]) - 
    360000.*(g22[2*(di - dj) + ijk] + g22[2*(-di + dj) + ijk]) - 
    22500.*(g22[3*(di - dj) + ijk] + g22[3*(-di + dj) + ijk]) - 
    625.*(g22[4*(di - dj) + ijk] + g22[4*(-di + dj) + ijk]) - 
    4.*(g22[5*(di - dj) + ijk] + g22[5*(-di + dj) + ijk]) + 
    360000.*(g22[-2*(di + dj) + ijk] + g22[2*(di + dj) + ijk]) + 
    22500.*(g22[-3*(di + dj) + ijk] + g22[3*(di + dj) + ijk]) + 
    625.*(g22[-4*(di + dj) + ijk] + g22[4*(di + dj) + ijk]) + 
    4.*(g22[-5*(di + dj) + ijk] + g22[5*(di + dj) + ijk]))
;

deldelg1223
=
6.2988158226253464349e-7*oo2dx*oo2dy*
  (4.41e6*(g23[-di - dj + ijk] - g23[di - dj + ijk] - g23[-di + dj + ijk] + 
       g23[di + dj + ijk]) + 1.26e6*
     (-g23[-di - 2*dj + ijk] + g23[di - 2*dj + ijk] - 
       g23[-2*di - dj + ijk] + g23[2*di - dj + ijk] + 
       g23[-2*di + dj + ijk] - g23[2*di + dj + ijk] + 
       g23[-di + 2*dj + ijk] - g23[di + 2*dj + ijk]) - 
    315000.*(g23[di - 3*dj + ijk] + g23[3*di - dj + ijk] + 
       g23[-3*di + dj + ijk] + g23[-di + 3*dj + ijk]) + 
    315000.*(g23[-di - 3*dj + ijk] + g23[-3*di - dj + ijk] + 
       g23[3*di + dj + ijk] + g23[di + 3*dj + ijk]) + 
    90000.*(-g23[-2*di - 3*dj + ijk] + g23[2*di - 3*dj + ijk] - 
       g23[-3*di - 2*dj + ijk] + g23[3*di - 2*dj + ijk] + 
       g23[-3*di + 2*dj + ijk] - g23[3*di + 2*dj + ijk] + 
       g23[-2*di + 3*dj + ijk] - g23[2*di + 3*dj + ijk]) - 
    15000.*(g23[2*di - 4*dj + ijk] + g23[4*di - 2*dj + ijk] + 
       g23[-4*di + 2*dj + ijk] + g23[-2*di + 4*dj + ijk]) + 
    52500.*(-g23[-di - 4*dj + ijk] + g23[di - 4*dj + ijk] - 
       g23[-4*di - dj + ijk] + g23[4*di - dj + ijk] + 
       g23[-4*di + dj + ijk] - g23[4*di + dj + ijk] + 
       g23[-di + 4*dj + ijk] - g23[di + 4*dj + ijk]) + 
    15000.*(g23[-2*di - 4*dj + ijk] + g23[-4*di - 2*dj + ijk] + 
       g23[4*di + 2*dj + ijk] + g23[2*di + 4*dj + ijk]) + 
    3750.*(-g23[-3*di - 4*dj + ijk] + g23[3*di - 4*dj + ijk] - 
       g23[-4*di - 3*dj + ijk] + g23[4*di - 3*dj + ijk] + 
       g23[-4*di + 3*dj + ijk] - g23[4*di + 3*dj + ijk] + 
       g23[-3*di + 4*dj + ijk] - g23[3*di + 4*dj + ijk]) - 
    300.*(g23[3*di - 5*dj + ijk] + g23[5*di - 3*dj + ijk] + 
       g23[-5*di + 3*dj + ijk] + g23[-3*di + 5*dj + ijk]) - 
    4200.*(g23[di - 5*dj + ijk] + g23[5*di - dj + ijk] + 
       g23[-5*di + dj + ijk] + g23[-di + 5*dj + ijk]) + 
    4200.*(g23[-di - 5*dj + ijk] + g23[-5*di - dj + ijk] + 
       g23[5*di + dj + ijk] + g23[di + 5*dj + ijk]) + 
    1200.*(-g23[-2*di - 5*dj + ijk] + g23[2*di - 5*dj + ijk] - 
       g23[-5*di - 2*dj + ijk] + g23[5*di - 2*dj + ijk] + 
       g23[-5*di + 2*dj + ijk] - g23[5*di + 2*dj + ijk] + 
       g23[-2*di + 5*dj + ijk] - g23[2*di + 5*dj + ijk]) + 
    300.*(g23[-3*di - 5*dj + ijk] + g23[-5*di - 3*dj + ijk] + 
       g23[5*di + 3*dj + ijk] + g23[3*di + 5*dj + ijk]) + 
    50.*(-g23[-4*di - 5*dj + ijk] + g23[4*di - 5*dj + ijk] - 
       g23[-5*di - 4*dj + ijk] + g23[5*di - 4*dj + ijk] + 
       g23[-5*di + 4*dj + ijk] - g23[5*di + 4*dj + ijk] + 
       g23[-4*di + 5*dj + ijk] - g23[4*di + 5*dj + ijk]) - 
    360000.*(g23[2*(di - dj) + ijk] + g23[2*(-di + dj) + ijk]) - 
    22500.*(g23[3*(di - dj) + ijk] + g23[3*(-di + dj) + ijk]) - 
    625.*(g23[4*(di - dj) + ijk] + g23[4*(-di + dj) + ijk]) - 
    4.*(g23[5*(di - dj) + ijk] + g23[5*(-di + dj) + ijk]) + 
    360000.*(g23[-2*(di + dj) + ijk] + g23[2*(di + dj) + ijk]) + 
    22500.*(g23[-3*(di + dj) + ijk] + g23[3*(di + dj) + ijk]) + 
    625.*(g23[-4*(di + dj) + ijk] + g23[4*(di + dj) + ijk]) + 
    4.*(g23[-5*(di + dj) + ijk] + g23[5*(di + dj) + ijk]))
;

deldelg1233
=
6.2988158226253464349e-7*oo2dx*oo2dy*
  (4.41e6*(g33[-di - dj + ijk] - g33[di - dj + ijk] - g33[-di + dj + ijk] + 
       g33[di + dj + ijk]) + 1.26e6*
     (-g33[-di - 2*dj + ijk] + g33[di - 2*dj + ijk] - 
       g33[-2*di - dj + ijk] + g33[2*di - dj + ijk] + 
       g33[-2*di + dj + ijk] - g33[2*di + dj + ijk] + 
       g33[-di + 2*dj + ijk] - g33[di + 2*dj + ijk]) - 
    315000.*(g33[di - 3*dj + ijk] + g33[3*di - dj + ijk] + 
       g33[-3*di + dj + ijk] + g33[-di + 3*dj + ijk]) + 
    315000.*(g33[-di - 3*dj + ijk] + g33[-3*di - dj + ijk] + 
       g33[3*di + dj + ijk] + g33[di + 3*dj + ijk]) + 
    90000.*(-g33[-2*di - 3*dj + ijk] + g33[2*di - 3*dj + ijk] - 
       g33[-3*di - 2*dj + ijk] + g33[3*di - 2*dj + ijk] + 
       g33[-3*di + 2*dj + ijk] - g33[3*di + 2*dj + ijk] + 
       g33[-2*di + 3*dj + ijk] - g33[2*di + 3*dj + ijk]) - 
    15000.*(g33[2*di - 4*dj + ijk] + g33[4*di - 2*dj + ijk] + 
       g33[-4*di + 2*dj + ijk] + g33[-2*di + 4*dj + ijk]) + 
    52500.*(-g33[-di - 4*dj + ijk] + g33[di - 4*dj + ijk] - 
       g33[-4*di - dj + ijk] + g33[4*di - dj + ijk] + 
       g33[-4*di + dj + ijk] - g33[4*di + dj + ijk] + 
       g33[-di + 4*dj + ijk] - g33[di + 4*dj + ijk]) + 
    15000.*(g33[-2*di - 4*dj + ijk] + g33[-4*di - 2*dj + ijk] + 
       g33[4*di + 2*dj + ijk] + g33[2*di + 4*dj + ijk]) + 
    3750.*(-g33[-3*di - 4*dj + ijk] + g33[3*di - 4*dj + ijk] - 
       g33[-4*di - 3*dj + ijk] + g33[4*di - 3*dj + ijk] + 
       g33[-4*di + 3*dj + ijk] - g33[4*di + 3*dj + ijk] + 
       g33[-3*di + 4*dj + ijk] - g33[3*di + 4*dj + ijk]) - 
    300.*(g33[3*di - 5*dj + ijk] + g33[5*di - 3*dj + ijk] + 
       g33[-5*di + 3*dj + ijk] + g33[-3*di + 5*dj + ijk]) - 
    4200.*(g33[di - 5*dj + ijk] + g33[5*di - dj + ijk] + 
       g33[-5*di + dj + ijk] + g33[-di + 5*dj + ijk]) + 
    4200.*(g33[-di - 5*dj + ijk] + g33[-5*di - dj + ijk] + 
       g33[5*di + dj + ijk] + g33[di + 5*dj + ijk]) + 
    1200.*(-g33[-2*di - 5*dj + ijk] + g33[2*di - 5*dj + ijk] - 
       g33[-5*di - 2*dj + ijk] + g33[5*di - 2*dj + ijk] + 
       g33[-5*di + 2*dj + ijk] - g33[5*di + 2*dj + ijk] + 
       g33[-2*di + 5*dj + ijk] - g33[2*di + 5*dj + ijk]) + 
    300.*(g33[-3*di - 5*dj + ijk] + g33[-5*di - 3*dj + ijk] + 
       g33[5*di + 3*dj + ijk] + g33[3*di + 5*dj + ijk]) + 
    50.*(-g33[-4*di - 5*dj + ijk] + g33[4*di - 5*dj + ijk] - 
       g33[-5*di - 4*dj + ijk] + g33[5*di - 4*dj + ijk] + 
       g33[-5*di + 4*dj + ijk] - g33[5*di + 4*dj + ijk] + 
       g33[-4*di + 5*dj + ijk] - g33[4*di + 5*dj + ijk]) - 
    360000.*(g33[2*(di - dj) + ijk] + g33[2*(-di + dj) + ijk]) - 
    22500.*(g33[3*(di - dj) + ijk] + g33[3*(-di + dj) + ijk]) - 
    625.*(g33[4*(di - dj) + ijk] + g33[4*(-di + dj) + ijk]) - 
    4.*(g33[5*(di - dj) + ijk] + g33[5*(-di + dj) + ijk]) + 
    360000.*(g33[-2*(di + dj) + ijk] + g33[2*(di + dj) + ijk]) + 
    22500.*(g33[-3*(di + dj) + ijk] + g33[3*(di + dj) + ijk]) + 
    625.*(g33[-4*(di + dj) + ijk] + g33[4*(di + dj) + ijk]) + 
    4.*(g33[-5*(di + dj) + ijk] + g33[5*(di + dj) + ijk]))
;

deldelg1311
=
6.2988158226253464349e-7*oo2dx*oo2dz*
  (4.41e6*(g11[-di - dk + ijk] - g11[di - dk + ijk] - g11[-di + dk + ijk] + 
       g11[di + dk + ijk]) + 1.26e6*
     (-g11[-di - 2*dk + ijk] + g11[di - 2*dk + ijk] - 
       g11[-2*di - dk + ijk] + g11[2*di - dk + ijk] + 
       g11[-2*di + dk + ijk] - g11[2*di + dk + ijk] + 
       g11[-di + 2*dk + ijk] - g11[di + 2*dk + ijk]) - 
    315000.*(g11[di - 3*dk + ijk] + g11[3*di - dk + ijk] + 
       g11[-3*di + dk + ijk] + g11[-di + 3*dk + ijk]) + 
    315000.*(g11[-di - 3*dk + ijk] + g11[-3*di - dk + ijk] + 
       g11[3*di + dk + ijk] + g11[di + 3*dk + ijk]) + 
    90000.*(-g11[-2*di - 3*dk + ijk] + g11[2*di - 3*dk + ijk] - 
       g11[-3*di - 2*dk + ijk] + g11[3*di - 2*dk + ijk] + 
       g11[-3*di + 2*dk + ijk] - g11[3*di + 2*dk + ijk] + 
       g11[-2*di + 3*dk + ijk] - g11[2*di + 3*dk + ijk]) - 
    15000.*(g11[2*di - 4*dk + ijk] + g11[4*di - 2*dk + ijk] + 
       g11[-4*di + 2*dk + ijk] + g11[-2*di + 4*dk + ijk]) + 
    52500.*(-g11[-di - 4*dk + ijk] + g11[di - 4*dk + ijk] - 
       g11[-4*di - dk + ijk] + g11[4*di - dk + ijk] + 
       g11[-4*di + dk + ijk] - g11[4*di + dk + ijk] + 
       g11[-di + 4*dk + ijk] - g11[di + 4*dk + ijk]) + 
    15000.*(g11[-2*di - 4*dk + ijk] + g11[-4*di - 2*dk + ijk] + 
       g11[4*di + 2*dk + ijk] + g11[2*di + 4*dk + ijk]) + 
    3750.*(-g11[-3*di - 4*dk + ijk] + g11[3*di - 4*dk + ijk] - 
       g11[-4*di - 3*dk + ijk] + g11[4*di - 3*dk + ijk] + 
       g11[-4*di + 3*dk + ijk] - g11[4*di + 3*dk + ijk] + 
       g11[-3*di + 4*dk + ijk] - g11[3*di + 4*dk + ijk]) - 
    300.*(g11[3*di - 5*dk + ijk] + g11[5*di - 3*dk + ijk] + 
       g11[-5*di + 3*dk + ijk] + g11[-3*di + 5*dk + ijk]) - 
    4200.*(g11[di - 5*dk + ijk] + g11[5*di - dk + ijk] + 
       g11[-5*di + dk + ijk] + g11[-di + 5*dk + ijk]) + 
    4200.*(g11[-di - 5*dk + ijk] + g11[-5*di - dk + ijk] + 
       g11[5*di + dk + ijk] + g11[di + 5*dk + ijk]) + 
    1200.*(-g11[-2*di - 5*dk + ijk] + g11[2*di - 5*dk + ijk] - 
       g11[-5*di - 2*dk + ijk] + g11[5*di - 2*dk + ijk] + 
       g11[-5*di + 2*dk + ijk] - g11[5*di + 2*dk + ijk] + 
       g11[-2*di + 5*dk + ijk] - g11[2*di + 5*dk + ijk]) + 
    300.*(g11[-3*di - 5*dk + ijk] + g11[-5*di - 3*dk + ijk] + 
       g11[5*di + 3*dk + ijk] + g11[3*di + 5*dk + ijk]) + 
    50.*(-g11[-4*di - 5*dk + ijk] + g11[4*di - 5*dk + ijk] - 
       g11[-5*di - 4*dk + ijk] + g11[5*di - 4*dk + ijk] + 
       g11[-5*di + 4*dk + ijk] - g11[5*di + 4*dk + ijk] + 
       g11[-4*di + 5*dk + ijk] - g11[4*di + 5*dk + ijk]) - 
    360000.*(g11[2*(di - dk) + ijk] + g11[2*(-di + dk) + ijk]) - 
    22500.*(g11[3*(di - dk) + ijk] + g11[3*(-di + dk) + ijk]) - 
    625.*(g11[4*(di - dk) + ijk] + g11[4*(-di + dk) + ijk]) - 
    4.*(g11[5*(di - dk) + ijk] + g11[5*(-di + dk) + ijk]) + 
    360000.*(g11[-2*(di + dk) + ijk] + g11[2*(di + dk) + ijk]) + 
    22500.*(g11[-3*(di + dk) + ijk] + g11[3*(di + dk) + ijk]) + 
    625.*(g11[-4*(di + dk) + ijk] + g11[4*(di + dk) + ijk]) + 
    4.*(g11[-5*(di + dk) + ijk] + g11[5*(di + dk) + ijk]))
;

deldelg1312
=
6.2988158226253464349e-7*oo2dx*oo2dz*
  (4.41e6*(g12[-di - dk + ijk] - g12[di - dk + ijk] - g12[-di + dk + ijk] + 
       g12[di + dk + ijk]) + 1.26e6*
     (-g12[-di - 2*dk + ijk] + g12[di - 2*dk + ijk] - 
       g12[-2*di - dk + ijk] + g12[2*di - dk + ijk] + 
       g12[-2*di + dk + ijk] - g12[2*di + dk + ijk] + 
       g12[-di + 2*dk + ijk] - g12[di + 2*dk + ijk]) - 
    315000.*(g12[di - 3*dk + ijk] + g12[3*di - dk + ijk] + 
       g12[-3*di + dk + ijk] + g12[-di + 3*dk + ijk]) + 
    315000.*(g12[-di - 3*dk + ijk] + g12[-3*di - dk + ijk] + 
       g12[3*di + dk + ijk] + g12[di + 3*dk + ijk]) + 
    90000.*(-g12[-2*di - 3*dk + ijk] + g12[2*di - 3*dk + ijk] - 
       g12[-3*di - 2*dk + ijk] + g12[3*di - 2*dk + ijk] + 
       g12[-3*di + 2*dk + ijk] - g12[3*di + 2*dk + ijk] + 
       g12[-2*di + 3*dk + ijk] - g12[2*di + 3*dk + ijk]) - 
    15000.*(g12[2*di - 4*dk + ijk] + g12[4*di - 2*dk + ijk] + 
       g12[-4*di + 2*dk + ijk] + g12[-2*di + 4*dk + ijk]) + 
    52500.*(-g12[-di - 4*dk + ijk] + g12[di - 4*dk + ijk] - 
       g12[-4*di - dk + ijk] + g12[4*di - dk + ijk] + 
       g12[-4*di + dk + ijk] - g12[4*di + dk + ijk] + 
       g12[-di + 4*dk + ijk] - g12[di + 4*dk + ijk]) + 
    15000.*(g12[-2*di - 4*dk + ijk] + g12[-4*di - 2*dk + ijk] + 
       g12[4*di + 2*dk + ijk] + g12[2*di + 4*dk + ijk]) + 
    3750.*(-g12[-3*di - 4*dk + ijk] + g12[3*di - 4*dk + ijk] - 
       g12[-4*di - 3*dk + ijk] + g12[4*di - 3*dk + ijk] + 
       g12[-4*di + 3*dk + ijk] - g12[4*di + 3*dk + ijk] + 
       g12[-3*di + 4*dk + ijk] - g12[3*di + 4*dk + ijk]) - 
    300.*(g12[3*di - 5*dk + ijk] + g12[5*di - 3*dk + ijk] + 
       g12[-5*di + 3*dk + ijk] + g12[-3*di + 5*dk + ijk]) - 
    4200.*(g12[di - 5*dk + ijk] + g12[5*di - dk + ijk] + 
       g12[-5*di + dk + ijk] + g12[-di + 5*dk + ijk]) + 
    4200.*(g12[-di - 5*dk + ijk] + g12[-5*di - dk + ijk] + 
       g12[5*di + dk + ijk] + g12[di + 5*dk + ijk]) + 
    1200.*(-g12[-2*di - 5*dk + ijk] + g12[2*di - 5*dk + ijk] - 
       g12[-5*di - 2*dk + ijk] + g12[5*di - 2*dk + ijk] + 
       g12[-5*di + 2*dk + ijk] - g12[5*di + 2*dk + ijk] + 
       g12[-2*di + 5*dk + ijk] - g12[2*di + 5*dk + ijk]) + 
    300.*(g12[-3*di - 5*dk + ijk] + g12[-5*di - 3*dk + ijk] + 
       g12[5*di + 3*dk + ijk] + g12[3*di + 5*dk + ijk]) + 
    50.*(-g12[-4*di - 5*dk + ijk] + g12[4*di - 5*dk + ijk] - 
       g12[-5*di - 4*dk + ijk] + g12[5*di - 4*dk + ijk] + 
       g12[-5*di + 4*dk + ijk] - g12[5*di + 4*dk + ijk] + 
       g12[-4*di + 5*dk + ijk] - g12[4*di + 5*dk + ijk]) - 
    360000.*(g12[2*(di - dk) + ijk] + g12[2*(-di + dk) + ijk]) - 
    22500.*(g12[3*(di - dk) + ijk] + g12[3*(-di + dk) + ijk]) - 
    625.*(g12[4*(di - dk) + ijk] + g12[4*(-di + dk) + ijk]) - 
    4.*(g12[5*(di - dk) + ijk] + g12[5*(-di + dk) + ijk]) + 
    360000.*(g12[-2*(di + dk) + ijk] + g12[2*(di + dk) + ijk]) + 
    22500.*(g12[-3*(di + dk) + ijk] + g12[3*(di + dk) + ijk]) + 
    625.*(g12[-4*(di + dk) + ijk] + g12[4*(di + dk) + ijk]) + 
    4.*(g12[-5*(di + dk) + ijk] + g12[5*(di + dk) + ijk]))
;

deldelg1313
=
6.2988158226253464349e-7*oo2dx*oo2dz*
  (4.41e6*(g13[-di - dk + ijk] - g13[di - dk + ijk] - g13[-di + dk + ijk] + 
       g13[di + dk + ijk]) + 1.26e6*
     (-g13[-di - 2*dk + ijk] + g13[di - 2*dk + ijk] - 
       g13[-2*di - dk + ijk] + g13[2*di - dk + ijk] + 
       g13[-2*di + dk + ijk] - g13[2*di + dk + ijk] + 
       g13[-di + 2*dk + ijk] - g13[di + 2*dk + ijk]) - 
    315000.*(g13[di - 3*dk + ijk] + g13[3*di - dk + ijk] + 
       g13[-3*di + dk + ijk] + g13[-di + 3*dk + ijk]) + 
    315000.*(g13[-di - 3*dk + ijk] + g13[-3*di - dk + ijk] + 
       g13[3*di + dk + ijk] + g13[di + 3*dk + ijk]) + 
    90000.*(-g13[-2*di - 3*dk + ijk] + g13[2*di - 3*dk + ijk] - 
       g13[-3*di - 2*dk + ijk] + g13[3*di - 2*dk + ijk] + 
       g13[-3*di + 2*dk + ijk] - g13[3*di + 2*dk + ijk] + 
       g13[-2*di + 3*dk + ijk] - g13[2*di + 3*dk + ijk]) - 
    15000.*(g13[2*di - 4*dk + ijk] + g13[4*di - 2*dk + ijk] + 
       g13[-4*di + 2*dk + ijk] + g13[-2*di + 4*dk + ijk]) + 
    52500.*(-g13[-di - 4*dk + ijk] + g13[di - 4*dk + ijk] - 
       g13[-4*di - dk + ijk] + g13[4*di - dk + ijk] + 
       g13[-4*di + dk + ijk] - g13[4*di + dk + ijk] + 
       g13[-di + 4*dk + ijk] - g13[di + 4*dk + ijk]) + 
    15000.*(g13[-2*di - 4*dk + ijk] + g13[-4*di - 2*dk + ijk] + 
       g13[4*di + 2*dk + ijk] + g13[2*di + 4*dk + ijk]) + 
    3750.*(-g13[-3*di - 4*dk + ijk] + g13[3*di - 4*dk + ijk] - 
       g13[-4*di - 3*dk + ijk] + g13[4*di - 3*dk + ijk] + 
       g13[-4*di + 3*dk + ijk] - g13[4*di + 3*dk + ijk] + 
       g13[-3*di + 4*dk + ijk] - g13[3*di + 4*dk + ijk]) - 
    300.*(g13[3*di - 5*dk + ijk] + g13[5*di - 3*dk + ijk] + 
       g13[-5*di + 3*dk + ijk] + g13[-3*di + 5*dk + ijk]) - 
    4200.*(g13[di - 5*dk + ijk] + g13[5*di - dk + ijk] + 
       g13[-5*di + dk + ijk] + g13[-di + 5*dk + ijk]) + 
    4200.*(g13[-di - 5*dk + ijk] + g13[-5*di - dk + ijk] + 
       g13[5*di + dk + ijk] + g13[di + 5*dk + ijk]) + 
    1200.*(-g13[-2*di - 5*dk + ijk] + g13[2*di - 5*dk + ijk] - 
       g13[-5*di - 2*dk + ijk] + g13[5*di - 2*dk + ijk] + 
       g13[-5*di + 2*dk + ijk] - g13[5*di + 2*dk + ijk] + 
       g13[-2*di + 5*dk + ijk] - g13[2*di + 5*dk + ijk]) + 
    300.*(g13[-3*di - 5*dk + ijk] + g13[-5*di - 3*dk + ijk] + 
       g13[5*di + 3*dk + ijk] + g13[3*di + 5*dk + ijk]) + 
    50.*(-g13[-4*di - 5*dk + ijk] + g13[4*di - 5*dk + ijk] - 
       g13[-5*di - 4*dk + ijk] + g13[5*di - 4*dk + ijk] + 
       g13[-5*di + 4*dk + ijk] - g13[5*di + 4*dk + ijk] + 
       g13[-4*di + 5*dk + ijk] - g13[4*di + 5*dk + ijk]) - 
    360000.*(g13[2*(di - dk) + ijk] + g13[2*(-di + dk) + ijk]) - 
    22500.*(g13[3*(di - dk) + ijk] + g13[3*(-di + dk) + ijk]) - 
    625.*(g13[4*(di - dk) + ijk] + g13[4*(-di + dk) + ijk]) - 
    4.*(g13[5*(di - dk) + ijk] + g13[5*(-di + dk) + ijk]) + 
    360000.*(g13[-2*(di + dk) + ijk] + g13[2*(di + dk) + ijk]) + 
    22500.*(g13[-3*(di + dk) + ijk] + g13[3*(di + dk) + ijk]) + 
    625.*(g13[-4*(di + dk) + ijk] + g13[4*(di + dk) + ijk]) + 
    4.*(g13[-5*(di + dk) + ijk] + g13[5*(di + dk) + ijk]))
;

deldelg1322
=
6.2988158226253464349e-7*oo2dx*oo2dz*
  (4.41e6*(g22[-di - dk + ijk] - g22[di - dk + ijk] - g22[-di + dk + ijk] + 
       g22[di + dk + ijk]) + 1.26e6*
     (-g22[-di - 2*dk + ijk] + g22[di - 2*dk + ijk] - 
       g22[-2*di - dk + ijk] + g22[2*di - dk + ijk] + 
       g22[-2*di + dk + ijk] - g22[2*di + dk + ijk] + 
       g22[-di + 2*dk + ijk] - g22[di + 2*dk + ijk]) - 
    315000.*(g22[di - 3*dk + ijk] + g22[3*di - dk + ijk] + 
       g22[-3*di + dk + ijk] + g22[-di + 3*dk + ijk]) + 
    315000.*(g22[-di - 3*dk + ijk] + g22[-3*di - dk + ijk] + 
       g22[3*di + dk + ijk] + g22[di + 3*dk + ijk]) + 
    90000.*(-g22[-2*di - 3*dk + ijk] + g22[2*di - 3*dk + ijk] - 
       g22[-3*di - 2*dk + ijk] + g22[3*di - 2*dk + ijk] + 
       g22[-3*di + 2*dk + ijk] - g22[3*di + 2*dk + ijk] + 
       g22[-2*di + 3*dk + ijk] - g22[2*di + 3*dk + ijk]) - 
    15000.*(g22[2*di - 4*dk + ijk] + g22[4*di - 2*dk + ijk] + 
       g22[-4*di + 2*dk + ijk] + g22[-2*di + 4*dk + ijk]) + 
    52500.*(-g22[-di - 4*dk + ijk] + g22[di - 4*dk + ijk] - 
       g22[-4*di - dk + ijk] + g22[4*di - dk + ijk] + 
       g22[-4*di + dk + ijk] - g22[4*di + dk + ijk] + 
       g22[-di + 4*dk + ijk] - g22[di + 4*dk + ijk]) + 
    15000.*(g22[-2*di - 4*dk + ijk] + g22[-4*di - 2*dk + ijk] + 
       g22[4*di + 2*dk + ijk] + g22[2*di + 4*dk + ijk]) + 
    3750.*(-g22[-3*di - 4*dk + ijk] + g22[3*di - 4*dk + ijk] - 
       g22[-4*di - 3*dk + ijk] + g22[4*di - 3*dk + ijk] + 
       g22[-4*di + 3*dk + ijk] - g22[4*di + 3*dk + ijk] + 
       g22[-3*di + 4*dk + ijk] - g22[3*di + 4*dk + ijk]) - 
    300.*(g22[3*di - 5*dk + ijk] + g22[5*di - 3*dk + ijk] + 
       g22[-5*di + 3*dk + ijk] + g22[-3*di + 5*dk + ijk]) - 
    4200.*(g22[di - 5*dk + ijk] + g22[5*di - dk + ijk] + 
       g22[-5*di + dk + ijk] + g22[-di + 5*dk + ijk]) + 
    4200.*(g22[-di - 5*dk + ijk] + g22[-5*di - dk + ijk] + 
       g22[5*di + dk + ijk] + g22[di + 5*dk + ijk]) + 
    1200.*(-g22[-2*di - 5*dk + ijk] + g22[2*di - 5*dk + ijk] - 
       g22[-5*di - 2*dk + ijk] + g22[5*di - 2*dk + ijk] + 
       g22[-5*di + 2*dk + ijk] - g22[5*di + 2*dk + ijk] + 
       g22[-2*di + 5*dk + ijk] - g22[2*di + 5*dk + ijk]) + 
    300.*(g22[-3*di - 5*dk + ijk] + g22[-5*di - 3*dk + ijk] + 
       g22[5*di + 3*dk + ijk] + g22[3*di + 5*dk + ijk]) + 
    50.*(-g22[-4*di - 5*dk + ijk] + g22[4*di - 5*dk + ijk] - 
       g22[-5*di - 4*dk + ijk] + g22[5*di - 4*dk + ijk] + 
       g22[-5*di + 4*dk + ijk] - g22[5*di + 4*dk + ijk] + 
       g22[-4*di + 5*dk + ijk] - g22[4*di + 5*dk + ijk]) - 
    360000.*(g22[2*(di - dk) + ijk] + g22[2*(-di + dk) + ijk]) - 
    22500.*(g22[3*(di - dk) + ijk] + g22[3*(-di + dk) + ijk]) - 
    625.*(g22[4*(di - dk) + ijk] + g22[4*(-di + dk) + ijk]) - 
    4.*(g22[5*(di - dk) + ijk] + g22[5*(-di + dk) + ijk]) + 
    360000.*(g22[-2*(di + dk) + ijk] + g22[2*(di + dk) + ijk]) + 
    22500.*(g22[-3*(di + dk) + ijk] + g22[3*(di + dk) + ijk]) + 
    625.*(g22[-4*(di + dk) + ijk] + g22[4*(di + dk) + ijk]) + 
    4.*(g22[-5*(di + dk) + ijk] + g22[5*(di + dk) + ijk]))
;

deldelg1323
=
6.2988158226253464349e-7*oo2dx*oo2dz*
  (4.41e6*(g23[-di - dk + ijk] - g23[di - dk + ijk] - g23[-di + dk + ijk] + 
       g23[di + dk + ijk]) + 1.26e6*
     (-g23[-di - 2*dk + ijk] + g23[di - 2*dk + ijk] - 
       g23[-2*di - dk + ijk] + g23[2*di - dk + ijk] + 
       g23[-2*di + dk + ijk] - g23[2*di + dk + ijk] + 
       g23[-di + 2*dk + ijk] - g23[di + 2*dk + ijk]) - 
    315000.*(g23[di - 3*dk + ijk] + g23[3*di - dk + ijk] + 
       g23[-3*di + dk + ijk] + g23[-di + 3*dk + ijk]) + 
    315000.*(g23[-di - 3*dk + ijk] + g23[-3*di - dk + ijk] + 
       g23[3*di + dk + ijk] + g23[di + 3*dk + ijk]) + 
    90000.*(-g23[-2*di - 3*dk + ijk] + g23[2*di - 3*dk + ijk] - 
       g23[-3*di - 2*dk + ijk] + g23[3*di - 2*dk + ijk] + 
       g23[-3*di + 2*dk + ijk] - g23[3*di + 2*dk + ijk] + 
       g23[-2*di + 3*dk + ijk] - g23[2*di + 3*dk + ijk]) - 
    15000.*(g23[2*di - 4*dk + ijk] + g23[4*di - 2*dk + ijk] + 
       g23[-4*di + 2*dk + ijk] + g23[-2*di + 4*dk + ijk]) + 
    52500.*(-g23[-di - 4*dk + ijk] + g23[di - 4*dk + ijk] - 
       g23[-4*di - dk + ijk] + g23[4*di - dk + ijk] + 
       g23[-4*di + dk + ijk] - g23[4*di + dk + ijk] + 
       g23[-di + 4*dk + ijk] - g23[di + 4*dk + ijk]) + 
    15000.*(g23[-2*di - 4*dk + ijk] + g23[-4*di - 2*dk + ijk] + 
       g23[4*di + 2*dk + ijk] + g23[2*di + 4*dk + ijk]) + 
    3750.*(-g23[-3*di - 4*dk + ijk] + g23[3*di - 4*dk + ijk] - 
       g23[-4*di - 3*dk + ijk] + g23[4*di - 3*dk + ijk] + 
       g23[-4*di + 3*dk + ijk] - g23[4*di + 3*dk + ijk] + 
       g23[-3*di + 4*dk + ijk] - g23[3*di + 4*dk + ijk]) - 
    300.*(g23[3*di - 5*dk + ijk] + g23[5*di - 3*dk + ijk] + 
       g23[-5*di + 3*dk + ijk] + g23[-3*di + 5*dk + ijk]) - 
    4200.*(g23[di - 5*dk + ijk] + g23[5*di - dk + ijk] + 
       g23[-5*di + dk + ijk] + g23[-di + 5*dk + ijk]) + 
    4200.*(g23[-di - 5*dk + ijk] + g23[-5*di - dk + ijk] + 
       g23[5*di + dk + ijk] + g23[di + 5*dk + ijk]) + 
    1200.*(-g23[-2*di - 5*dk + ijk] + g23[2*di - 5*dk + ijk] - 
       g23[-5*di - 2*dk + ijk] + g23[5*di - 2*dk + ijk] + 
       g23[-5*di + 2*dk + ijk] - g23[5*di + 2*dk + ijk] + 
       g23[-2*di + 5*dk + ijk] - g23[2*di + 5*dk + ijk]) + 
    300.*(g23[-3*di - 5*dk + ijk] + g23[-5*di - 3*dk + ijk] + 
       g23[5*di + 3*dk + ijk] + g23[3*di + 5*dk + ijk]) + 
    50.*(-g23[-4*di - 5*dk + ijk] + g23[4*di - 5*dk + ijk] - 
       g23[-5*di - 4*dk + ijk] + g23[5*di - 4*dk + ijk] + 
       g23[-5*di + 4*dk + ijk] - g23[5*di + 4*dk + ijk] + 
       g23[-4*di + 5*dk + ijk] - g23[4*di + 5*dk + ijk]) - 
    360000.*(g23[2*(di - dk) + ijk] + g23[2*(-di + dk) + ijk]) - 
    22500.*(g23[3*(di - dk) + ijk] + g23[3*(-di + dk) + ijk]) - 
    625.*(g23[4*(di - dk) + ijk] + g23[4*(-di + dk) + ijk]) - 
    4.*(g23[5*(di - dk) + ijk] + g23[5*(-di + dk) + ijk]) + 
    360000.*(g23[-2*(di + dk) + ijk] + g23[2*(di + dk) + ijk]) + 
    22500.*(g23[-3*(di + dk) + ijk] + g23[3*(di + dk) + ijk]) + 
    625.*(g23[-4*(di + dk) + ijk] + g23[4*(di + dk) + ijk]) + 
    4.*(g23[-5*(di + dk) + ijk] + g23[5*(di + dk) + ijk]))
;

deldelg1333
=
6.2988158226253464349e-7*oo2dx*oo2dz*
  (4.41e6*(g33[-di - dk + ijk] - g33[di - dk + ijk] - g33[-di + dk + ijk] + 
       g33[di + dk + ijk]) + 1.26e6*
     (-g33[-di - 2*dk + ijk] + g33[di - 2*dk + ijk] - 
       g33[-2*di - dk + ijk] + g33[2*di - dk + ijk] + 
       g33[-2*di + dk + ijk] - g33[2*di + dk + ijk] + 
       g33[-di + 2*dk + ijk] - g33[di + 2*dk + ijk]) - 
    315000.*(g33[di - 3*dk + ijk] + g33[3*di - dk + ijk] + 
       g33[-3*di + dk + ijk] + g33[-di + 3*dk + ijk]) + 
    315000.*(g33[-di - 3*dk + ijk] + g33[-3*di - dk + ijk] + 
       g33[3*di + dk + ijk] + g33[di + 3*dk + ijk]) + 
    90000.*(-g33[-2*di - 3*dk + ijk] + g33[2*di - 3*dk + ijk] - 
       g33[-3*di - 2*dk + ijk] + g33[3*di - 2*dk + ijk] + 
       g33[-3*di + 2*dk + ijk] - g33[3*di + 2*dk + ijk] + 
       g33[-2*di + 3*dk + ijk] - g33[2*di + 3*dk + ijk]) - 
    15000.*(g33[2*di - 4*dk + ijk] + g33[4*di - 2*dk + ijk] + 
       g33[-4*di + 2*dk + ijk] + g33[-2*di + 4*dk + ijk]) + 
    52500.*(-g33[-di - 4*dk + ijk] + g33[di - 4*dk + ijk] - 
       g33[-4*di - dk + ijk] + g33[4*di - dk + ijk] + 
       g33[-4*di + dk + ijk] - g33[4*di + dk + ijk] + 
       g33[-di + 4*dk + ijk] - g33[di + 4*dk + ijk]) + 
    15000.*(g33[-2*di - 4*dk + ijk] + g33[-4*di - 2*dk + ijk] + 
       g33[4*di + 2*dk + ijk] + g33[2*di + 4*dk + ijk]) + 
    3750.*(-g33[-3*di - 4*dk + ijk] + g33[3*di - 4*dk + ijk] - 
       g33[-4*di - 3*dk + ijk] + g33[4*di - 3*dk + ijk] + 
       g33[-4*di + 3*dk + ijk] - g33[4*di + 3*dk + ijk] + 
       g33[-3*di + 4*dk + ijk] - g33[3*di + 4*dk + ijk]) - 
    300.*(g33[3*di - 5*dk + ijk] + g33[5*di - 3*dk + ijk] + 
       g33[-5*di + 3*dk + ijk] + g33[-3*di + 5*dk + ijk]) - 
    4200.*(g33[di - 5*dk + ijk] + g33[5*di - dk + ijk] + 
       g33[-5*di + dk + ijk] + g33[-di + 5*dk + ijk]) + 
    4200.*(g33[-di - 5*dk + ijk] + g33[-5*di - dk + ijk] + 
       g33[5*di + dk + ijk] + g33[di + 5*dk + ijk]) + 
    1200.*(-g33[-2*di - 5*dk + ijk] + g33[2*di - 5*dk + ijk] - 
       g33[-5*di - 2*dk + ijk] + g33[5*di - 2*dk + ijk] + 
       g33[-5*di + 2*dk + ijk] - g33[5*di + 2*dk + ijk] + 
       g33[-2*di + 5*dk + ijk] - g33[2*di + 5*dk + ijk]) + 
    300.*(g33[-3*di - 5*dk + ijk] + g33[-5*di - 3*dk + ijk] + 
       g33[5*di + 3*dk + ijk] + g33[3*di + 5*dk + ijk]) + 
    50.*(-g33[-4*di - 5*dk + ijk] + g33[4*di - 5*dk + ijk] - 
       g33[-5*di - 4*dk + ijk] + g33[5*di - 4*dk + ijk] + 
       g33[-5*di + 4*dk + ijk] - g33[5*di + 4*dk + ijk] + 
       g33[-4*di + 5*dk + ijk] - g33[4*di + 5*dk + ijk]) - 
    360000.*(g33[2*(di - dk) + ijk] + g33[2*(-di + dk) + ijk]) - 
    22500.*(g33[3*(di - dk) + ijk] + g33[3*(-di + dk) + ijk]) - 
    625.*(g33[4*(di - dk) + ijk] + g33[4*(-di + dk) + ijk]) - 
    4.*(g33[5*(di - dk) + ijk] + g33[5*(-di + dk) + ijk]) + 
    360000.*(g33[-2*(di + dk) + ijk] + g33[2*(di + dk) + ijk]) + 
    22500.*(g33[-3*(di + dk) + ijk] + g33[3*(di + dk) + ijk]) + 
    625.*(g33[-4*(di + dk) + ijk] + g33[4*(di + dk) + ijk]) + 
    4.*(g33[-5*(di + dk) + ijk] + g33[5*(di + dk) + ijk]))
;

deldelg2211
=
0.00003968253968253968254*oody2*
  (-73766.*g11[ijk] + 42000.*(g11[-dj + ijk] + g11[dj + ijk]) - 
    6000.*(g11[-2*dj + ijk] + g11[2*dj + ijk]) + 
    1000.*(g11[-3*dj + ijk] + g11[3*dj + ijk]) - 
    125.*(g11[-4*dj + ijk] + g11[4*dj + ijk]) + 
    8.*(g11[-5*dj + ijk] + g11[5*dj + ijk]))
;

deldelg2212
=
0.00003968253968253968254*oody2*
  (-73766.*g12[ijk] + 42000.*(g12[-dj + ijk] + g12[dj + ijk]) - 
    6000.*(g12[-2*dj + ijk] + g12[2*dj + ijk]) + 
    1000.*(g12[-3*dj + ijk] + g12[3*dj + ijk]) - 
    125.*(g12[-4*dj + ijk] + g12[4*dj + ijk]) + 
    8.*(g12[-5*dj + ijk] + g12[5*dj + ijk]))
;

deldelg2213
=
0.00003968253968253968254*oody2*
  (-73766.*g13[ijk] + 42000.*(g13[-dj + ijk] + g13[dj + ijk]) - 
    6000.*(g13[-2*dj + ijk] + g13[2*dj + ijk]) + 
    1000.*(g13[-3*dj + ijk] + g13[3*dj + ijk]) - 
    125.*(g13[-4*dj + ijk] + g13[4*dj + ijk]) + 
    8.*(g13[-5*dj + ijk] + g13[5*dj + ijk]))
;

deldelg2222
=
0.00003968253968253968254*oody2*
  (-73766.*g22[ijk] + 42000.*(g22[-dj + ijk] + g22[dj + ijk]) - 
    6000.*(g22[-2*dj + ijk] + g22[2*dj + ijk]) + 
    1000.*(g22[-3*dj + ijk] + g22[3*dj + ijk]) - 
    125.*(g22[-4*dj + ijk] + g22[4*dj + ijk]) + 
    8.*(g22[-5*dj + ijk] + g22[5*dj + ijk]))
;

deldelg2223
=
0.00003968253968253968254*oody2*
  (-73766.*g23[ijk] + 42000.*(g23[-dj + ijk] + g23[dj + ijk]) - 
    6000.*(g23[-2*dj + ijk] + g23[2*dj + ijk]) + 
    1000.*(g23[-3*dj + ijk] + g23[3*dj + ijk]) - 
    125.*(g23[-4*dj + ijk] + g23[4*dj + ijk]) + 
    8.*(g23[-5*dj + ijk] + g23[5*dj + ijk]))
;

deldelg2233
=
0.00003968253968253968254*oody2*
  (-73766.*g33[ijk] + 42000.*(g33[-dj + ijk] + g33[dj + ijk]) - 
    6000.*(g33[-2*dj + ijk] + g33[2*dj + ijk]) + 
    1000.*(g33[-3*dj + ijk] + g33[3*dj + ijk]) - 
    125.*(g33[-4*dj + ijk] + g33[4*dj + ijk]) + 
    8.*(g33[-5*dj + ijk] + g33[5*dj + ijk]))
;

deldelg2311
=
6.2988158226253464349e-7*oo2dy*oo2dz*
  (4.41e6*(g11[-dj - dk + ijk] - g11[dj - dk + ijk] - g11[-dj + dk + ijk] + 
       g11[dj + dk + ijk]) + 1.26e6*
     (-g11[-dj - 2*dk + ijk] + g11[dj - 2*dk + ijk] - 
       g11[-2*dj - dk + ijk] + g11[2*dj - dk + ijk] + 
       g11[-2*dj + dk + ijk] - g11[2*dj + dk + ijk] + 
       g11[-dj + 2*dk + ijk] - g11[dj + 2*dk + ijk]) - 
    315000.*(g11[dj - 3*dk + ijk] + g11[3*dj - dk + ijk] + 
       g11[-3*dj + dk + ijk] + g11[-dj + 3*dk + ijk]) + 
    315000.*(g11[-dj - 3*dk + ijk] + g11[-3*dj - dk + ijk] + 
       g11[3*dj + dk + ijk] + g11[dj + 3*dk + ijk]) + 
    90000.*(-g11[-2*dj - 3*dk + ijk] + g11[2*dj - 3*dk + ijk] - 
       g11[-3*dj - 2*dk + ijk] + g11[3*dj - 2*dk + ijk] + 
       g11[-3*dj + 2*dk + ijk] - g11[3*dj + 2*dk + ijk] + 
       g11[-2*dj + 3*dk + ijk] - g11[2*dj + 3*dk + ijk]) - 
    15000.*(g11[2*dj - 4*dk + ijk] + g11[4*dj - 2*dk + ijk] + 
       g11[-4*dj + 2*dk + ijk] + g11[-2*dj + 4*dk + ijk]) + 
    52500.*(-g11[-dj - 4*dk + ijk] + g11[dj - 4*dk + ijk] - 
       g11[-4*dj - dk + ijk] + g11[4*dj - dk + ijk] + 
       g11[-4*dj + dk + ijk] - g11[4*dj + dk + ijk] + 
       g11[-dj + 4*dk + ijk] - g11[dj + 4*dk + ijk]) + 
    15000.*(g11[-2*dj - 4*dk + ijk] + g11[-4*dj - 2*dk + ijk] + 
       g11[4*dj + 2*dk + ijk] + g11[2*dj + 4*dk + ijk]) + 
    3750.*(-g11[-3*dj - 4*dk + ijk] + g11[3*dj - 4*dk + ijk] - 
       g11[-4*dj - 3*dk + ijk] + g11[4*dj - 3*dk + ijk] + 
       g11[-4*dj + 3*dk + ijk] - g11[4*dj + 3*dk + ijk] + 
       g11[-3*dj + 4*dk + ijk] - g11[3*dj + 4*dk + ijk]) - 
    300.*(g11[3*dj - 5*dk + ijk] + g11[5*dj - 3*dk + ijk] + 
       g11[-5*dj + 3*dk + ijk] + g11[-3*dj + 5*dk + ijk]) - 
    4200.*(g11[dj - 5*dk + ijk] + g11[5*dj - dk + ijk] + 
       g11[-5*dj + dk + ijk] + g11[-dj + 5*dk + ijk]) + 
    4200.*(g11[-dj - 5*dk + ijk] + g11[-5*dj - dk + ijk] + 
       g11[5*dj + dk + ijk] + g11[dj + 5*dk + ijk]) + 
    1200.*(-g11[-2*dj - 5*dk + ijk] + g11[2*dj - 5*dk + ijk] - 
       g11[-5*dj - 2*dk + ijk] + g11[5*dj - 2*dk + ijk] + 
       g11[-5*dj + 2*dk + ijk] - g11[5*dj + 2*dk + ijk] + 
       g11[-2*dj + 5*dk + ijk] - g11[2*dj + 5*dk + ijk]) + 
    300.*(g11[-3*dj - 5*dk + ijk] + g11[-5*dj - 3*dk + ijk] + 
       g11[5*dj + 3*dk + ijk] + g11[3*dj + 5*dk + ijk]) + 
    50.*(-g11[-4*dj - 5*dk + ijk] + g11[4*dj - 5*dk + ijk] - 
       g11[-5*dj - 4*dk + ijk] + g11[5*dj - 4*dk + ijk] + 
       g11[-5*dj + 4*dk + ijk] - g11[5*dj + 4*dk + ijk] + 
       g11[-4*dj + 5*dk + ijk] - g11[4*dj + 5*dk + ijk]) - 
    360000.*(g11[2*(dj - dk) + ijk] + g11[2*(-dj + dk) + ijk]) - 
    22500.*(g11[3*(dj - dk) + ijk] + g11[3*(-dj + dk) + ijk]) - 
    625.*(g11[4*(dj - dk) + ijk] + g11[4*(-dj + dk) + ijk]) - 
    4.*(g11[5*(dj - dk) + ijk] + g11[5*(-dj + dk) + ijk]) + 
    360000.*(g11[-2*(dj + dk) + ijk] + g11[2*(dj + dk) + ijk]) + 
    22500.*(g11[-3*(dj + dk) + ijk] + g11[3*(dj + dk) + ijk]) + 
    625.*(g11[-4*(dj + dk) + ijk] + g11[4*(dj + dk) + ijk]) + 
    4.*(g11[-5*(dj + dk) + ijk] + g11[5*(dj + dk) + ijk]))
;

deldelg2312
=
6.2988158226253464349e-7*oo2dy*oo2dz*
  (4.41e6*(g12[-dj - dk + ijk] - g12[dj - dk + ijk] - g12[-dj + dk + ijk] + 
       g12[dj + dk + ijk]) + 1.26e6*
     (-g12[-dj - 2*dk + ijk] + g12[dj - 2*dk + ijk] - 
       g12[-2*dj - dk + ijk] + g12[2*dj - dk + ijk] + 
       g12[-2*dj + dk + ijk] - g12[2*dj + dk + ijk] + 
       g12[-dj + 2*dk + ijk] - g12[dj + 2*dk + ijk]) - 
    315000.*(g12[dj - 3*dk + ijk] + g12[3*dj - dk + ijk] + 
       g12[-3*dj + dk + ijk] + g12[-dj + 3*dk + ijk]) + 
    315000.*(g12[-dj - 3*dk + ijk] + g12[-3*dj - dk + ijk] + 
       g12[3*dj + dk + ijk] + g12[dj + 3*dk + ijk]) + 
    90000.*(-g12[-2*dj - 3*dk + ijk] + g12[2*dj - 3*dk + ijk] - 
       g12[-3*dj - 2*dk + ijk] + g12[3*dj - 2*dk + ijk] + 
       g12[-3*dj + 2*dk + ijk] - g12[3*dj + 2*dk + ijk] + 
       g12[-2*dj + 3*dk + ijk] - g12[2*dj + 3*dk + ijk]) - 
    15000.*(g12[2*dj - 4*dk + ijk] + g12[4*dj - 2*dk + ijk] + 
       g12[-4*dj + 2*dk + ijk] + g12[-2*dj + 4*dk + ijk]) + 
    52500.*(-g12[-dj - 4*dk + ijk] + g12[dj - 4*dk + ijk] - 
       g12[-4*dj - dk + ijk] + g12[4*dj - dk + ijk] + 
       g12[-4*dj + dk + ijk] - g12[4*dj + dk + ijk] + 
       g12[-dj + 4*dk + ijk] - g12[dj + 4*dk + ijk]) + 
    15000.*(g12[-2*dj - 4*dk + ijk] + g12[-4*dj - 2*dk + ijk] + 
       g12[4*dj + 2*dk + ijk] + g12[2*dj + 4*dk + ijk]) + 
    3750.*(-g12[-3*dj - 4*dk + ijk] + g12[3*dj - 4*dk + ijk] - 
       g12[-4*dj - 3*dk + ijk] + g12[4*dj - 3*dk + ijk] + 
       g12[-4*dj + 3*dk + ijk] - g12[4*dj + 3*dk + ijk] + 
       g12[-3*dj + 4*dk + ijk] - g12[3*dj + 4*dk + ijk]) - 
    300.*(g12[3*dj - 5*dk + ijk] + g12[5*dj - 3*dk + ijk] + 
       g12[-5*dj + 3*dk + ijk] + g12[-3*dj + 5*dk + ijk]) - 
    4200.*(g12[dj - 5*dk + ijk] + g12[5*dj - dk + ijk] + 
       g12[-5*dj + dk + ijk] + g12[-dj + 5*dk + ijk]) + 
    4200.*(g12[-dj - 5*dk + ijk] + g12[-5*dj - dk + ijk] + 
       g12[5*dj + dk + ijk] + g12[dj + 5*dk + ijk]) + 
    1200.*(-g12[-2*dj - 5*dk + ijk] + g12[2*dj - 5*dk + ijk] - 
       g12[-5*dj - 2*dk + ijk] + g12[5*dj - 2*dk + ijk] + 
       g12[-5*dj + 2*dk + ijk] - g12[5*dj + 2*dk + ijk] + 
       g12[-2*dj + 5*dk + ijk] - g12[2*dj + 5*dk + ijk]) + 
    300.*(g12[-3*dj - 5*dk + ijk] + g12[-5*dj - 3*dk + ijk] + 
       g12[5*dj + 3*dk + ijk] + g12[3*dj + 5*dk + ijk]) + 
    50.*(-g12[-4*dj - 5*dk + ijk] + g12[4*dj - 5*dk + ijk] - 
       g12[-5*dj - 4*dk + ijk] + g12[5*dj - 4*dk + ijk] + 
       g12[-5*dj + 4*dk + ijk] - g12[5*dj + 4*dk + ijk] + 
       g12[-4*dj + 5*dk + ijk] - g12[4*dj + 5*dk + ijk]) - 
    360000.*(g12[2*(dj - dk) + ijk] + g12[2*(-dj + dk) + ijk]) - 
    22500.*(g12[3*(dj - dk) + ijk] + g12[3*(-dj + dk) + ijk]) - 
    625.*(g12[4*(dj - dk) + ijk] + g12[4*(-dj + dk) + ijk]) - 
    4.*(g12[5*(dj - dk) + ijk] + g12[5*(-dj + dk) + ijk]) + 
    360000.*(g12[-2*(dj + dk) + ijk] + g12[2*(dj + dk) + ijk]) + 
    22500.*(g12[-3*(dj + dk) + ijk] + g12[3*(dj + dk) + ijk]) + 
    625.*(g12[-4*(dj + dk) + ijk] + g12[4*(dj + dk) + ijk]) + 
    4.*(g12[-5*(dj + dk) + ijk] + g12[5*(dj + dk) + ijk]))
;

deldelg2313
=
6.2988158226253464349e-7*oo2dy*oo2dz*
  (4.41e6*(g13[-dj - dk + ijk] - g13[dj - dk + ijk] - g13[-dj + dk + ijk] + 
       g13[dj + dk + ijk]) + 1.26e6*
     (-g13[-dj - 2*dk + ijk] + g13[dj - 2*dk + ijk] - 
       g13[-2*dj - dk + ijk] + g13[2*dj - dk + ijk] + 
       g13[-2*dj + dk + ijk] - g13[2*dj + dk + ijk] + 
       g13[-dj + 2*dk + ijk] - g13[dj + 2*dk + ijk]) - 
    315000.*(g13[dj - 3*dk + ijk] + g13[3*dj - dk + ijk] + 
       g13[-3*dj + dk + ijk] + g13[-dj + 3*dk + ijk]) + 
    315000.*(g13[-dj - 3*dk + ijk] + g13[-3*dj - dk + ijk] + 
       g13[3*dj + dk + ijk] + g13[dj + 3*dk + ijk]) + 
    90000.*(-g13[-2*dj - 3*dk + ijk] + g13[2*dj - 3*dk + ijk] - 
       g13[-3*dj - 2*dk + ijk] + g13[3*dj - 2*dk + ijk] + 
       g13[-3*dj + 2*dk + ijk] - g13[3*dj + 2*dk + ijk] + 
       g13[-2*dj + 3*dk + ijk] - g13[2*dj + 3*dk + ijk]) - 
    15000.*(g13[2*dj - 4*dk + ijk] + g13[4*dj - 2*dk + ijk] + 
       g13[-4*dj + 2*dk + ijk] + g13[-2*dj + 4*dk + ijk]) + 
    52500.*(-g13[-dj - 4*dk + ijk] + g13[dj - 4*dk + ijk] - 
       g13[-4*dj - dk + ijk] + g13[4*dj - dk + ijk] + 
       g13[-4*dj + dk + ijk] - g13[4*dj + dk + ijk] + 
       g13[-dj + 4*dk + ijk] - g13[dj + 4*dk + ijk]) + 
    15000.*(g13[-2*dj - 4*dk + ijk] + g13[-4*dj - 2*dk + ijk] + 
       g13[4*dj + 2*dk + ijk] + g13[2*dj + 4*dk + ijk]) + 
    3750.*(-g13[-3*dj - 4*dk + ijk] + g13[3*dj - 4*dk + ijk] - 
       g13[-4*dj - 3*dk + ijk] + g13[4*dj - 3*dk + ijk] + 
       g13[-4*dj + 3*dk + ijk] - g13[4*dj + 3*dk + ijk] + 
       g13[-3*dj + 4*dk + ijk] - g13[3*dj + 4*dk + ijk]) - 
    300.*(g13[3*dj - 5*dk + ijk] + g13[5*dj - 3*dk + ijk] + 
       g13[-5*dj + 3*dk + ijk] + g13[-3*dj + 5*dk + ijk]) - 
    4200.*(g13[dj - 5*dk + ijk] + g13[5*dj - dk + ijk] + 
       g13[-5*dj + dk + ijk] + g13[-dj + 5*dk + ijk]) + 
    4200.*(g13[-dj - 5*dk + ijk] + g13[-5*dj - dk + ijk] + 
       g13[5*dj + dk + ijk] + g13[dj + 5*dk + ijk]) + 
    1200.*(-g13[-2*dj - 5*dk + ijk] + g13[2*dj - 5*dk + ijk] - 
       g13[-5*dj - 2*dk + ijk] + g13[5*dj - 2*dk + ijk] + 
       g13[-5*dj + 2*dk + ijk] - g13[5*dj + 2*dk + ijk] + 
       g13[-2*dj + 5*dk + ijk] - g13[2*dj + 5*dk + ijk]) + 
    300.*(g13[-3*dj - 5*dk + ijk] + g13[-5*dj - 3*dk + ijk] + 
       g13[5*dj + 3*dk + ijk] + g13[3*dj + 5*dk + ijk]) + 
    50.*(-g13[-4*dj - 5*dk + ijk] + g13[4*dj - 5*dk + ijk] - 
       g13[-5*dj - 4*dk + ijk] + g13[5*dj - 4*dk + ijk] + 
       g13[-5*dj + 4*dk + ijk] - g13[5*dj + 4*dk + ijk] + 
       g13[-4*dj + 5*dk + ijk] - g13[4*dj + 5*dk + ijk]) - 
    360000.*(g13[2*(dj - dk) + ijk] + g13[2*(-dj + dk) + ijk]) - 
    22500.*(g13[3*(dj - dk) + ijk] + g13[3*(-dj + dk) + ijk]) - 
    625.*(g13[4*(dj - dk) + ijk] + g13[4*(-dj + dk) + ijk]) - 
    4.*(g13[5*(dj - dk) + ijk] + g13[5*(-dj + dk) + ijk]) + 
    360000.*(g13[-2*(dj + dk) + ijk] + g13[2*(dj + dk) + ijk]) + 
    22500.*(g13[-3*(dj + dk) + ijk] + g13[3*(dj + dk) + ijk]) + 
    625.*(g13[-4*(dj + dk) + ijk] + g13[4*(dj + dk) + ijk]) + 
    4.*(g13[-5*(dj + dk) + ijk] + g13[5*(dj + dk) + ijk]))
;

deldelg2322
=
6.2988158226253464349e-7*oo2dy*oo2dz*
  (4.41e6*(g22[-dj - dk + ijk] - g22[dj - dk + ijk] - g22[-dj + dk + ijk] + 
       g22[dj + dk + ijk]) + 1.26e6*
     (-g22[-dj - 2*dk + ijk] + g22[dj - 2*dk + ijk] - 
       g22[-2*dj - dk + ijk] + g22[2*dj - dk + ijk] + 
       g22[-2*dj + dk + ijk] - g22[2*dj + dk + ijk] + 
       g22[-dj + 2*dk + ijk] - g22[dj + 2*dk + ijk]) - 
    315000.*(g22[dj - 3*dk + ijk] + g22[3*dj - dk + ijk] + 
       g22[-3*dj + dk + ijk] + g22[-dj + 3*dk + ijk]) + 
    315000.*(g22[-dj - 3*dk + ijk] + g22[-3*dj - dk + ijk] + 
       g22[3*dj + dk + ijk] + g22[dj + 3*dk + ijk]) + 
    90000.*(-g22[-2*dj - 3*dk + ijk] + g22[2*dj - 3*dk + ijk] - 
       g22[-3*dj - 2*dk + ijk] + g22[3*dj - 2*dk + ijk] + 
       g22[-3*dj + 2*dk + ijk] - g22[3*dj + 2*dk + ijk] + 
       g22[-2*dj + 3*dk + ijk] - g22[2*dj + 3*dk + ijk]) - 
    15000.*(g22[2*dj - 4*dk + ijk] + g22[4*dj - 2*dk + ijk] + 
       g22[-4*dj + 2*dk + ijk] + g22[-2*dj + 4*dk + ijk]) + 
    52500.*(-g22[-dj - 4*dk + ijk] + g22[dj - 4*dk + ijk] - 
       g22[-4*dj - dk + ijk] + g22[4*dj - dk + ijk] + 
       g22[-4*dj + dk + ijk] - g22[4*dj + dk + ijk] + 
       g22[-dj + 4*dk + ijk] - g22[dj + 4*dk + ijk]) + 
    15000.*(g22[-2*dj - 4*dk + ijk] + g22[-4*dj - 2*dk + ijk] + 
       g22[4*dj + 2*dk + ijk] + g22[2*dj + 4*dk + ijk]) + 
    3750.*(-g22[-3*dj - 4*dk + ijk] + g22[3*dj - 4*dk + ijk] - 
       g22[-4*dj - 3*dk + ijk] + g22[4*dj - 3*dk + ijk] + 
       g22[-4*dj + 3*dk + ijk] - g22[4*dj + 3*dk + ijk] + 
       g22[-3*dj + 4*dk + ijk] - g22[3*dj + 4*dk + ijk]) - 
    300.*(g22[3*dj - 5*dk + ijk] + g22[5*dj - 3*dk + ijk] + 
       g22[-5*dj + 3*dk + ijk] + g22[-3*dj + 5*dk + ijk]) - 
    4200.*(g22[dj - 5*dk + ijk] + g22[5*dj - dk + ijk] + 
       g22[-5*dj + dk + ijk] + g22[-dj + 5*dk + ijk]) + 
    4200.*(g22[-dj - 5*dk + ijk] + g22[-5*dj - dk + ijk] + 
       g22[5*dj + dk + ijk] + g22[dj + 5*dk + ijk]) + 
    1200.*(-g22[-2*dj - 5*dk + ijk] + g22[2*dj - 5*dk + ijk] - 
       g22[-5*dj - 2*dk + ijk] + g22[5*dj - 2*dk + ijk] + 
       g22[-5*dj + 2*dk + ijk] - g22[5*dj + 2*dk + ijk] + 
       g22[-2*dj + 5*dk + ijk] - g22[2*dj + 5*dk + ijk]) + 
    300.*(g22[-3*dj - 5*dk + ijk] + g22[-5*dj - 3*dk + ijk] + 
       g22[5*dj + 3*dk + ijk] + g22[3*dj + 5*dk + ijk]) + 
    50.*(-g22[-4*dj - 5*dk + ijk] + g22[4*dj - 5*dk + ijk] - 
       g22[-5*dj - 4*dk + ijk] + g22[5*dj - 4*dk + ijk] + 
       g22[-5*dj + 4*dk + ijk] - g22[5*dj + 4*dk + ijk] + 
       g22[-4*dj + 5*dk + ijk] - g22[4*dj + 5*dk + ijk]) - 
    360000.*(g22[2*(dj - dk) + ijk] + g22[2*(-dj + dk) + ijk]) - 
    22500.*(g22[3*(dj - dk) + ijk] + g22[3*(-dj + dk) + ijk]) - 
    625.*(g22[4*(dj - dk) + ijk] + g22[4*(-dj + dk) + ijk]) - 
    4.*(g22[5*(dj - dk) + ijk] + g22[5*(-dj + dk) + ijk]) + 
    360000.*(g22[-2*(dj + dk) + ijk] + g22[2*(dj + dk) + ijk]) + 
    22500.*(g22[-3*(dj + dk) + ijk] + g22[3*(dj + dk) + ijk]) + 
    625.*(g22[-4*(dj + dk) + ijk] + g22[4*(dj + dk) + ijk]) + 
    4.*(g22[-5*(dj + dk) + ijk] + g22[5*(dj + dk) + ijk]))
;

deldelg2323
=
6.2988158226253464349e-7*oo2dy*oo2dz*
  (4.41e6*(g23[-dj - dk + ijk] - g23[dj - dk + ijk] - g23[-dj + dk + ijk] + 
       g23[dj + dk + ijk]) + 1.26e6*
     (-g23[-dj - 2*dk + ijk] + g23[dj - 2*dk + ijk] - 
       g23[-2*dj - dk + ijk] + g23[2*dj - dk + ijk] + 
       g23[-2*dj + dk + ijk] - g23[2*dj + dk + ijk] + 
       g23[-dj + 2*dk + ijk] - g23[dj + 2*dk + ijk]) - 
    315000.*(g23[dj - 3*dk + ijk] + g23[3*dj - dk + ijk] + 
       g23[-3*dj + dk + ijk] + g23[-dj + 3*dk + ijk]) + 
    315000.*(g23[-dj - 3*dk + ijk] + g23[-3*dj - dk + ijk] + 
       g23[3*dj + dk + ijk] + g23[dj + 3*dk + ijk]) + 
    90000.*(-g23[-2*dj - 3*dk + ijk] + g23[2*dj - 3*dk + ijk] - 
       g23[-3*dj - 2*dk + ijk] + g23[3*dj - 2*dk + ijk] + 
       g23[-3*dj + 2*dk + ijk] - g23[3*dj + 2*dk + ijk] + 
       g23[-2*dj + 3*dk + ijk] - g23[2*dj + 3*dk + ijk]) - 
    15000.*(g23[2*dj - 4*dk + ijk] + g23[4*dj - 2*dk + ijk] + 
       g23[-4*dj + 2*dk + ijk] + g23[-2*dj + 4*dk + ijk]) + 
    52500.*(-g23[-dj - 4*dk + ijk] + g23[dj - 4*dk + ijk] - 
       g23[-4*dj - dk + ijk] + g23[4*dj - dk + ijk] + 
       g23[-4*dj + dk + ijk] - g23[4*dj + dk + ijk] + 
       g23[-dj + 4*dk + ijk] - g23[dj + 4*dk + ijk]) + 
    15000.*(g23[-2*dj - 4*dk + ijk] + g23[-4*dj - 2*dk + ijk] + 
       g23[4*dj + 2*dk + ijk] + g23[2*dj + 4*dk + ijk]) + 
    3750.*(-g23[-3*dj - 4*dk + ijk] + g23[3*dj - 4*dk + ijk] - 
       g23[-4*dj - 3*dk + ijk] + g23[4*dj - 3*dk + ijk] + 
       g23[-4*dj + 3*dk + ijk] - g23[4*dj + 3*dk + ijk] + 
       g23[-3*dj + 4*dk + ijk] - g23[3*dj + 4*dk + ijk]) - 
    300.*(g23[3*dj - 5*dk + ijk] + g23[5*dj - 3*dk + ijk] + 
       g23[-5*dj + 3*dk + ijk] + g23[-3*dj + 5*dk + ijk]) - 
    4200.*(g23[dj - 5*dk + ijk] + g23[5*dj - dk + ijk] + 
       g23[-5*dj + dk + ijk] + g23[-dj + 5*dk + ijk]) + 
    4200.*(g23[-dj - 5*dk + ijk] + g23[-5*dj - dk + ijk] + 
       g23[5*dj + dk + ijk] + g23[dj + 5*dk + ijk]) + 
    1200.*(-g23[-2*dj - 5*dk + ijk] + g23[2*dj - 5*dk + ijk] - 
       g23[-5*dj - 2*dk + ijk] + g23[5*dj - 2*dk + ijk] + 
       g23[-5*dj + 2*dk + ijk] - g23[5*dj + 2*dk + ijk] + 
       g23[-2*dj + 5*dk + ijk] - g23[2*dj + 5*dk + ijk]) + 
    300.*(g23[-3*dj - 5*dk + ijk] + g23[-5*dj - 3*dk + ijk] + 
       g23[5*dj + 3*dk + ijk] + g23[3*dj + 5*dk + ijk]) + 
    50.*(-g23[-4*dj - 5*dk + ijk] + g23[4*dj - 5*dk + ijk] - 
       g23[-5*dj - 4*dk + ijk] + g23[5*dj - 4*dk + ijk] + 
       g23[-5*dj + 4*dk + ijk] - g23[5*dj + 4*dk + ijk] + 
       g23[-4*dj + 5*dk + ijk] - g23[4*dj + 5*dk + ijk]) - 
    360000.*(g23[2*(dj - dk) + ijk] + g23[2*(-dj + dk) + ijk]) - 
    22500.*(g23[3*(dj - dk) + ijk] + g23[3*(-dj + dk) + ijk]) - 
    625.*(g23[4*(dj - dk) + ijk] + g23[4*(-dj + dk) + ijk]) - 
    4.*(g23[5*(dj - dk) + ijk] + g23[5*(-dj + dk) + ijk]) + 
    360000.*(g23[-2*(dj + dk) + ijk] + g23[2*(dj + dk) + ijk]) + 
    22500.*(g23[-3*(dj + dk) + ijk] + g23[3*(dj + dk) + ijk]) + 
    625.*(g23[-4*(dj + dk) + ijk] + g23[4*(dj + dk) + ijk]) + 
    4.*(g23[-5*(dj + dk) + ijk] + g23[5*(dj + dk) + ijk]))
;

deldelg2333
=
6.2988158226253464349e-7*oo2dy*oo2dz*
  (4.41e6*(g33[-dj - dk + ijk] - g33[dj - dk + ijk] - g33[-dj + dk + ijk] + 
       g33[dj + dk + ijk]) + 1.26e6*
     (-g33[-dj - 2*dk + ijk] + g33[dj - 2*dk + ijk] - 
       g33[-2*dj - dk + ijk] + g33[2*dj - dk + ijk] + 
       g33[-2*dj + dk + ijk] - g33[2*dj + dk + ijk] + 
       g33[-dj + 2*dk + ijk] - g33[dj + 2*dk + ijk]) - 
    315000.*(g33[dj - 3*dk + ijk] + g33[3*dj - dk + ijk] + 
       g33[-3*dj + dk + ijk] + g33[-dj + 3*dk + ijk]) + 
    315000.*(g33[-dj - 3*dk + ijk] + g33[-3*dj - dk + ijk] + 
       g33[3*dj + dk + ijk] + g33[dj + 3*dk + ijk]) + 
    90000.*(-g33[-2*dj - 3*dk + ijk] + g33[2*dj - 3*dk + ijk] - 
       g33[-3*dj - 2*dk + ijk] + g33[3*dj - 2*dk + ijk] + 
       g33[-3*dj + 2*dk + ijk] - g33[3*dj + 2*dk + ijk] + 
       g33[-2*dj + 3*dk + ijk] - g33[2*dj + 3*dk + ijk]) - 
    15000.*(g33[2*dj - 4*dk + ijk] + g33[4*dj - 2*dk + ijk] + 
       g33[-4*dj + 2*dk + ijk] + g33[-2*dj + 4*dk + ijk]) + 
    52500.*(-g33[-dj - 4*dk + ijk] + g33[dj - 4*dk + ijk] - 
       g33[-4*dj - dk + ijk] + g33[4*dj - dk + ijk] + 
       g33[-4*dj + dk + ijk] - g33[4*dj + dk + ijk] + 
       g33[-dj + 4*dk + ijk] - g33[dj + 4*dk + ijk]) + 
    15000.*(g33[-2*dj - 4*dk + ijk] + g33[-4*dj - 2*dk + ijk] + 
       g33[4*dj + 2*dk + ijk] + g33[2*dj + 4*dk + ijk]) + 
    3750.*(-g33[-3*dj - 4*dk + ijk] + g33[3*dj - 4*dk + ijk] - 
       g33[-4*dj - 3*dk + ijk] + g33[4*dj - 3*dk + ijk] + 
       g33[-4*dj + 3*dk + ijk] - g33[4*dj + 3*dk + ijk] + 
       g33[-3*dj + 4*dk + ijk] - g33[3*dj + 4*dk + ijk]) - 
    300.*(g33[3*dj - 5*dk + ijk] + g33[5*dj - 3*dk + ijk] + 
       g33[-5*dj + 3*dk + ijk] + g33[-3*dj + 5*dk + ijk]) - 
    4200.*(g33[dj - 5*dk + ijk] + g33[5*dj - dk + ijk] + 
       g33[-5*dj + dk + ijk] + g33[-dj + 5*dk + ijk]) + 
    4200.*(g33[-dj - 5*dk + ijk] + g33[-5*dj - dk + ijk] + 
       g33[5*dj + dk + ijk] + g33[dj + 5*dk + ijk]) + 
    1200.*(-g33[-2*dj - 5*dk + ijk] + g33[2*dj - 5*dk + ijk] - 
       g33[-5*dj - 2*dk + ijk] + g33[5*dj - 2*dk + ijk] + 
       g33[-5*dj + 2*dk + ijk] - g33[5*dj + 2*dk + ijk] + 
       g33[-2*dj + 5*dk + ijk] - g33[2*dj + 5*dk + ijk]) + 
    300.*(g33[-3*dj - 5*dk + ijk] + g33[-5*dj - 3*dk + ijk] + 
       g33[5*dj + 3*dk + ijk] + g33[3*dj + 5*dk + ijk]) + 
    50.*(-g33[-4*dj - 5*dk + ijk] + g33[4*dj - 5*dk + ijk] - 
       g33[-5*dj - 4*dk + ijk] + g33[5*dj - 4*dk + ijk] + 
       g33[-5*dj + 4*dk + ijk] - g33[5*dj + 4*dk + ijk] + 
       g33[-4*dj + 5*dk + ijk] - g33[4*dj + 5*dk + ijk]) - 
    360000.*(g33[2*(dj - dk) + ijk] + g33[2*(-dj + dk) + ijk]) - 
    22500.*(g33[3*(dj - dk) + ijk] + g33[3*(-dj + dk) + ijk]) - 
    625.*(g33[4*(dj - dk) + ijk] + g33[4*(-dj + dk) + ijk]) - 
    4.*(g33[5*(dj - dk) + ijk] + g33[5*(-dj + dk) + ijk]) + 
    360000.*(g33[-2*(dj + dk) + ijk] + g33[2*(dj + dk) + ijk]) + 
    22500.*(g33[-3*(dj + dk) + ijk] + g33[3*(dj + dk) + ijk]) + 
    625.*(g33[-4*(dj + dk) + ijk] + g33[4*(dj + dk) + ijk]) + 
    4.*(g33[-5*(dj + dk) + ijk] + g33[5*(dj + dk) + ijk]))
;

deldelg3311
=
0.00003968253968253968254*oodz2*
  (-73766.*g11[ijk] + 42000.*(g11[-dk + ijk] + g11[dk + ijk]) - 
    6000.*(g11[-2*dk + ijk] + g11[2*dk + ijk]) + 
    1000.*(g11[-3*dk + ijk] + g11[3*dk + ijk]) - 
    125.*(g11[-4*dk + ijk] + g11[4*dk + ijk]) + 
    8.*(g11[-5*dk + ijk] + g11[5*dk + ijk]))
;

deldelg3312
=
0.00003968253968253968254*oodz2*
  (-73766.*g12[ijk] + 42000.*(g12[-dk + ijk] + g12[dk + ijk]) - 
    6000.*(g12[-2*dk + ijk] + g12[2*dk + ijk]) + 
    1000.*(g12[-3*dk + ijk] + g12[3*dk + ijk]) - 
    125.*(g12[-4*dk + ijk] + g12[4*dk + ijk]) + 
    8.*(g12[-5*dk + ijk] + g12[5*dk + ijk]))
;

deldelg3313
=
0.00003968253968253968254*oodz2*
  (-73766.*g13[ijk] + 42000.*(g13[-dk + ijk] + g13[dk + ijk]) - 
    6000.*(g13[-2*dk + ijk] + g13[2*dk + ijk]) + 
    1000.*(g13[-3*dk + ijk] + g13[3*dk + ijk]) - 
    125.*(g13[-4*dk + ijk] + g13[4*dk + ijk]) + 
    8.*(g13[-5*dk + ijk] + g13[5*dk + ijk]))
;

deldelg3322
=
0.00003968253968253968254*oodz2*
  (-73766.*g22[ijk] + 42000.*(g22[-dk + ijk] + g22[dk + ijk]) - 
    6000.*(g22[-2*dk + ijk] + g22[2*dk + ijk]) + 
    1000.*(g22[-3*dk + ijk] + g22[3*dk + ijk]) - 
    125.*(g22[-4*dk + ijk] + g22[4*dk + ijk]) + 
    8.*(g22[-5*dk + ijk] + g22[5*dk + ijk]))
;

deldelg3323
=
0.00003968253968253968254*oodz2*
  (-73766.*g23[ijk] + 42000.*(g23[-dk + ijk] + g23[dk + ijk]) - 
    6000.*(g23[-2*dk + ijk] + g23[2*dk + ijk]) + 
    1000.*(g23[-3*dk + ijk] + g23[3*dk + ijk]) - 
    125.*(g23[-4*dk + ijk] + g23[4*dk + ijk]) + 
    8.*(g23[-5*dk + ijk] + g23[5*dk + ijk]))
;

deldelg3333
=
0.00003968253968253968254*oodz2*
  (-73766.*g33[ijk] + 42000.*(g33[-dk + ijk] + g33[dk + ijk]) - 
    6000.*(g33[-2*dk + ijk] + g33[2*dk + ijk]) + 
    1000.*(g33[-3*dk + ijk] + g33[3*dk + ijk]) - 
    125.*(g33[-4*dk + ijk] + g33[4*dk + ijk]) + 
    8.*(g33[-5*dk + ijk] + g33[5*dk + ijk]))
;

delG11
=
0.00079365079365079365079*oo2dx*
  (-2.*G1[-5*di + ijk] + 25.*G1[-4*di + ijk] + 600.*G1[-2*di + ijk] + 
    2100.*(-G1[-di + ijk] + G1[di + ijk]) - 600.*G1[2*di + ijk] + 
    150.*(-G1[-3*di + ijk] + G1[3*di + ijk]) - 25.*G1[4*di + ijk] + 
    2.*G1[5*di + ijk])
;

delG12
=
0.00079365079365079365079*oo2dx*
  (-2.*G2[-5*di + ijk] + 25.*G2[-4*di + ijk] + 600.*G2[-2*di + ijk] + 
    2100.*(-G2[-di + ijk] + G2[di + ijk]) - 600.*G2[2*di + ijk] + 
    150.*(-G2[-3*di + ijk] + G2[3*di + ijk]) - 25.*G2[4*di + ijk] + 
    2.*G2[5*di + ijk])
;

delG13
=
0.00079365079365079365079*oo2dx*
  (-2.*G3[-5*di + ijk] + 25.*G3[-4*di + ijk] + 600.*G3[-2*di + ijk] + 
    2100.*(-G3[-di + ijk] + G3[di + ijk]) - 600.*G3[2*di + ijk] + 
    150.*(-G3[-3*di + ijk] + G3[3*di + ijk]) - 25.*G3[4*di + ijk] + 
    2.*G3[5*di + ijk])
;

delG21
=
0.00079365079365079365079*oo2dy*
  (-2.*G1[-5*dj + ijk] + 25.*G1[-4*dj + ijk] + 600.*G1[-2*dj + ijk] + 
    2100.*(-G1[-dj + ijk] + G1[dj + ijk]) - 600.*G1[2*dj + ijk] + 
    150.*(-G1[-3*dj + ijk] + G1[3*dj + ijk]) - 25.*G1[4*dj + ijk] + 
    2.*G1[5*dj + ijk])
;

delG22
=
0.00079365079365079365079*oo2dy*
  (-2.*G2[-5*dj + ijk] + 25.*G2[-4*dj + ijk] + 600.*G2[-2*dj + ijk] + 
    2100.*(-G2[-dj + ijk] + G2[dj + ijk]) - 600.*G2[2*dj + ijk] + 
    150.*(-G2[-3*dj + ijk] + G2[3*dj + ijk]) - 25.*G2[4*dj + ijk] + 
    2.*G2[5*dj + ijk])
;

delG23
=
0.00079365079365079365079*oo2dy*
  (-2.*G3[-5*dj + ijk] + 25.*G3[-4*dj + ijk] + 600.*G3[-2*dj + ijk] + 
    2100.*(-G3[-dj + ijk] + G3[dj + ijk]) - 600.*G3[2*dj + ijk] + 
    150.*(-G3[-3*dj + ijk] + G3[3*dj + ijk]) - 25.*G3[4*dj + ijk] + 
    2.*G3[5*dj + ijk])
;

delG31
=
0.00079365079365079365079*oo2dz*
  (-2.*G1[-5*dk + ijk] + 25.*G1[-4*dk + ijk] + 600.*G1[-2*dk + ijk] + 
    2100.*(-G1[-dk + ijk] + G1[dk + ijk]) - 600.*G1[2*dk + ijk] + 
    150.*(-G1[-3*dk + ijk] + G1[3*dk + ijk]) - 25.*G1[4*dk + ijk] + 
    2.*G1[5*dk + ijk])
;

delG32
=
0.00079365079365079365079*oo2dz*
  (-2.*G2[-5*dk + ijk] + 25.*G2[-4*dk + ijk] + 600.*G2[-2*dk + ijk] + 
    2100.*(-G2[-dk + ijk] + G2[dk + ijk]) - 600.*G2[2*dk + ijk] + 
    150.*(-G2[-3*dk + ijk] + G2[3*dk + ijk]) - 25.*G2[4*dk + ijk] + 
    2.*G2[5*dk + ijk])
;

delG33
=
0.00079365079365079365079*oo2dz*
  (-2.*G3[-5*dk + ijk] + 25.*G3[-4*dk + ijk] + 600.*G3[-2*dk + ijk] + 
    2100.*(-G3[-dk + ijk] + G3[dk + ijk]) - 600.*G3[2*dk + ijk] + 
    150.*(-G3[-3*dk + ijk] + G3[3*dk + ijk]) - 25.*G3[4*dk + ijk] + 
    2.*G3[5*dk + ijk])
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

dchi1
=
0.00079365079365079365079*oo2dx*
  (-2.*chi[-5*di + ijk] + 25.*chi[-4*di + ijk] + 600.*chi[-2*di + ijk] + 
    2100.*(-chi[-di + ijk] + chi[di + ijk]) - 600.*chi[2*di + ijk] + 
    150.*(-chi[-3*di + ijk] + chi[3*di + ijk]) - 25.*chi[4*di + ijk] + 
    2.*chi[5*di + ijk])
;

dchi2
=
0.00079365079365079365079*oo2dy*
  (-2.*chi[-5*dj + ijk] + 25.*chi[-4*dj + ijk] + 600.*chi[-2*dj + ijk] + 
    2100.*(-chi[-dj + ijk] + chi[dj + ijk]) - 600.*chi[2*dj + ijk] + 
    150.*(-chi[-3*dj + ijk] + chi[3*dj + ijk]) - 25.*chi[4*dj + ijk] + 
    2.*chi[5*dj + ijk])
;

dchi3
=
0.00079365079365079365079*oo2dz*
  (-2.*chi[-5*dk + ijk] + 25.*chi[-4*dk + ijk] + 600.*chi[-2*dk + ijk] + 
    2100.*(-chi[-dk + ijk] + chi[dk + ijk]) - 600.*chi[2*dk + ijk] + 
    150.*(-chi[-3*dk + ijk] + chi[3*dk + ijk]) - 25.*chi[4*dk + ijk] + 
    2.*chi[5*dk + ijk])
;

ddchi11
=
0.00003968253968253968254*oodx2*
  (-73766.*chi[ijk] + 42000.*(chi[-di + ijk] + chi[di + ijk]) - 
    6000.*(chi[-2*di + ijk] + chi[2*di + ijk]) + 
    1000.*(chi[-3*di + ijk] + chi[3*di + ijk]) - 
    125.*(chi[-4*di + ijk] + chi[4*di + ijk]) + 
    8.*(chi[-5*di + ijk] + chi[5*di + ijk]))
;

ddchi12
=
6.2988158226253464349e-7*oo2dx*oo2dy*
  (4.41e6*(chi[-di - dj + ijk] - chi[di - dj + ijk] - chi[-di + dj + ijk] + 
       chi[di + dj + ijk]) + 1.26e6*
     (-chi[-di - 2*dj + ijk] + chi[di - 2*dj + ijk] - 
       chi[-2*di - dj + ijk] + chi[2*di - dj + ijk] + 
       chi[-2*di + dj + ijk] - chi[2*di + dj + ijk] + 
       chi[-di + 2*dj + ijk] - chi[di + 2*dj + ijk]) - 
    315000.*(chi[di - 3*dj + ijk] + chi[3*di - dj + ijk] + 
       chi[-3*di + dj + ijk] + chi[-di + 3*dj + ijk]) + 
    315000.*(chi[-di - 3*dj + ijk] + chi[-3*di - dj + ijk] + 
       chi[3*di + dj + ijk] + chi[di + 3*dj + ijk]) + 
    90000.*(-chi[-2*di - 3*dj + ijk] + chi[2*di - 3*dj + ijk] - 
       chi[-3*di - 2*dj + ijk] + chi[3*di - 2*dj + ijk] + 
       chi[-3*di + 2*dj + ijk] - chi[3*di + 2*dj + ijk] + 
       chi[-2*di + 3*dj + ijk] - chi[2*di + 3*dj + ijk]) - 
    15000.*(chi[2*di - 4*dj + ijk] + chi[4*di - 2*dj + ijk] + 
       chi[-4*di + 2*dj + ijk] + chi[-2*di + 4*dj + ijk]) + 
    52500.*(-chi[-di - 4*dj + ijk] + chi[di - 4*dj + ijk] - 
       chi[-4*di - dj + ijk] + chi[4*di - dj + ijk] + 
       chi[-4*di + dj + ijk] - chi[4*di + dj + ijk] + 
       chi[-di + 4*dj + ijk] - chi[di + 4*dj + ijk]) + 
    15000.*(chi[-2*di - 4*dj + ijk] + chi[-4*di - 2*dj + ijk] + 
       chi[4*di + 2*dj + ijk] + chi[2*di + 4*dj + ijk]) + 
    3750.*(-chi[-3*di - 4*dj + ijk] + chi[3*di - 4*dj + ijk] - 
       chi[-4*di - 3*dj + ijk] + chi[4*di - 3*dj + ijk] + 
       chi[-4*di + 3*dj + ijk] - chi[4*di + 3*dj + ijk] + 
       chi[-3*di + 4*dj + ijk] - chi[3*di + 4*dj + ijk]) - 
    300.*(chi[3*di - 5*dj + ijk] + chi[5*di - 3*dj + ijk] + 
       chi[-5*di + 3*dj + ijk] + chi[-3*di + 5*dj + ijk]) - 
    4200.*(chi[di - 5*dj + ijk] + chi[5*di - dj + ijk] + 
       chi[-5*di + dj + ijk] + chi[-di + 5*dj + ijk]) + 
    4200.*(chi[-di - 5*dj + ijk] + chi[-5*di - dj + ijk] + 
       chi[5*di + dj + ijk] + chi[di + 5*dj + ijk]) + 
    1200.*(-chi[-2*di - 5*dj + ijk] + chi[2*di - 5*dj + ijk] - 
       chi[-5*di - 2*dj + ijk] + chi[5*di - 2*dj + ijk] + 
       chi[-5*di + 2*dj + ijk] - chi[5*di + 2*dj + ijk] + 
       chi[-2*di + 5*dj + ijk] - chi[2*di + 5*dj + ijk]) + 
    300.*(chi[-3*di - 5*dj + ijk] + chi[-5*di - 3*dj + ijk] + 
       chi[5*di + 3*dj + ijk] + chi[3*di + 5*dj + ijk]) + 
    50.*(-chi[-4*di - 5*dj + ijk] + chi[4*di - 5*dj + ijk] - 
       chi[-5*di - 4*dj + ijk] + chi[5*di - 4*dj + ijk] + 
       chi[-5*di + 4*dj + ijk] - chi[5*di + 4*dj + ijk] + 
       chi[-4*di + 5*dj + ijk] - chi[4*di + 5*dj + ijk]) - 
    360000.*(chi[2*(di - dj) + ijk] + chi[2*(-di + dj) + ijk]) - 
    22500.*(chi[3*(di - dj) + ijk] + chi[3*(-di + dj) + ijk]) - 
    625.*(chi[4*(di - dj) + ijk] + chi[4*(-di + dj) + ijk]) - 
    4.*(chi[5*(di - dj) + ijk] + chi[5*(-di + dj) + ijk]) + 
    360000.*(chi[-2*(di + dj) + ijk] + chi[2*(di + dj) + ijk]) + 
    22500.*(chi[-3*(di + dj) + ijk] + chi[3*(di + dj) + ijk]) + 
    625.*(chi[-4*(di + dj) + ijk] + chi[4*(di + dj) + ijk]) + 
    4.*(chi[-5*(di + dj) + ijk] + chi[5*(di + dj) + ijk]))
;

ddchi13
=
6.2988158226253464349e-7*oo2dx*oo2dz*
  (4.41e6*(chi[-di - dk + ijk] - chi[di - dk + ijk] - chi[-di + dk + ijk] + 
       chi[di + dk + ijk]) + 1.26e6*
     (-chi[-di - 2*dk + ijk] + chi[di - 2*dk + ijk] - 
       chi[-2*di - dk + ijk] + chi[2*di - dk + ijk] + 
       chi[-2*di + dk + ijk] - chi[2*di + dk + ijk] + 
       chi[-di + 2*dk + ijk] - chi[di + 2*dk + ijk]) - 
    315000.*(chi[di - 3*dk + ijk] + chi[3*di - dk + ijk] + 
       chi[-3*di + dk + ijk] + chi[-di + 3*dk + ijk]) + 
    315000.*(chi[-di - 3*dk + ijk] + chi[-3*di - dk + ijk] + 
       chi[3*di + dk + ijk] + chi[di + 3*dk + ijk]) + 
    90000.*(-chi[-2*di - 3*dk + ijk] + chi[2*di - 3*dk + ijk] - 
       chi[-3*di - 2*dk + ijk] + chi[3*di - 2*dk + ijk] + 
       chi[-3*di + 2*dk + ijk] - chi[3*di + 2*dk + ijk] + 
       chi[-2*di + 3*dk + ijk] - chi[2*di + 3*dk + ijk]) - 
    15000.*(chi[2*di - 4*dk + ijk] + chi[4*di - 2*dk + ijk] + 
       chi[-4*di + 2*dk + ijk] + chi[-2*di + 4*dk + ijk]) + 
    52500.*(-chi[-di - 4*dk + ijk] + chi[di - 4*dk + ijk] - 
       chi[-4*di - dk + ijk] + chi[4*di - dk + ijk] + 
       chi[-4*di + dk + ijk] - chi[4*di + dk + ijk] + 
       chi[-di + 4*dk + ijk] - chi[di + 4*dk + ijk]) + 
    15000.*(chi[-2*di - 4*dk + ijk] + chi[-4*di - 2*dk + ijk] + 
       chi[4*di + 2*dk + ijk] + chi[2*di + 4*dk + ijk]) + 
    3750.*(-chi[-3*di - 4*dk + ijk] + chi[3*di - 4*dk + ijk] - 
       chi[-4*di - 3*dk + ijk] + chi[4*di - 3*dk + ijk] + 
       chi[-4*di + 3*dk + ijk] - chi[4*di + 3*dk + ijk] + 
       chi[-3*di + 4*dk + ijk] - chi[3*di + 4*dk + ijk]) - 
    300.*(chi[3*di - 5*dk + ijk] + chi[5*di - 3*dk + ijk] + 
       chi[-5*di + 3*dk + ijk] + chi[-3*di + 5*dk + ijk]) - 
    4200.*(chi[di - 5*dk + ijk] + chi[5*di - dk + ijk] + 
       chi[-5*di + dk + ijk] + chi[-di + 5*dk + ijk]) + 
    4200.*(chi[-di - 5*dk + ijk] + chi[-5*di - dk + ijk] + 
       chi[5*di + dk + ijk] + chi[di + 5*dk + ijk]) + 
    1200.*(-chi[-2*di - 5*dk + ijk] + chi[2*di - 5*dk + ijk] - 
       chi[-5*di - 2*dk + ijk] + chi[5*di - 2*dk + ijk] + 
       chi[-5*di + 2*dk + ijk] - chi[5*di + 2*dk + ijk] + 
       chi[-2*di + 5*dk + ijk] - chi[2*di + 5*dk + ijk]) + 
    300.*(chi[-3*di - 5*dk + ijk] + chi[-5*di - 3*dk + ijk] + 
       chi[5*di + 3*dk + ijk] + chi[3*di + 5*dk + ijk]) + 
    50.*(-chi[-4*di - 5*dk + ijk] + chi[4*di - 5*dk + ijk] - 
       chi[-5*di - 4*dk + ijk] + chi[5*di - 4*dk + ijk] + 
       chi[-5*di + 4*dk + ijk] - chi[5*di + 4*dk + ijk] + 
       chi[-4*di + 5*dk + ijk] - chi[4*di + 5*dk + ijk]) - 
    360000.*(chi[2*(di - dk) + ijk] + chi[2*(-di + dk) + ijk]) - 
    22500.*(chi[3*(di - dk) + ijk] + chi[3*(-di + dk) + ijk]) - 
    625.*(chi[4*(di - dk) + ijk] + chi[4*(-di + dk) + ijk]) - 
    4.*(chi[5*(di - dk) + ijk] + chi[5*(-di + dk) + ijk]) + 
    360000.*(chi[-2*(di + dk) + ijk] + chi[2*(di + dk) + ijk]) + 
    22500.*(chi[-3*(di + dk) + ijk] + chi[3*(di + dk) + ijk]) + 
    625.*(chi[-4*(di + dk) + ijk] + chi[4*(di + dk) + ijk]) + 
    4.*(chi[-5*(di + dk) + ijk] + chi[5*(di + dk) + ijk]))
;

ddchi22
=
0.00003968253968253968254*oody2*
  (-73766.*chi[ijk] + 42000.*(chi[-dj + ijk] + chi[dj + ijk]) - 
    6000.*(chi[-2*dj + ijk] + chi[2*dj + ijk]) + 
    1000.*(chi[-3*dj + ijk] + chi[3*dj + ijk]) - 
    125.*(chi[-4*dj + ijk] + chi[4*dj + ijk]) + 
    8.*(chi[-5*dj + ijk] + chi[5*dj + ijk]))
;

ddchi23
=
6.2988158226253464349e-7*oo2dy*oo2dz*
  (4.41e6*(chi[-dj - dk + ijk] - chi[dj - dk + ijk] - chi[-dj + dk + ijk] + 
       chi[dj + dk + ijk]) + 1.26e6*
     (-chi[-dj - 2*dk + ijk] + chi[dj - 2*dk + ijk] - 
       chi[-2*dj - dk + ijk] + chi[2*dj - dk + ijk] + 
       chi[-2*dj + dk + ijk] - chi[2*dj + dk + ijk] + 
       chi[-dj + 2*dk + ijk] - chi[dj + 2*dk + ijk]) - 
    315000.*(chi[dj - 3*dk + ijk] + chi[3*dj - dk + ijk] + 
       chi[-3*dj + dk + ijk] + chi[-dj + 3*dk + ijk]) + 
    315000.*(chi[-dj - 3*dk + ijk] + chi[-3*dj - dk + ijk] + 
       chi[3*dj + dk + ijk] + chi[dj + 3*dk + ijk]) + 
    90000.*(-chi[-2*dj - 3*dk + ijk] + chi[2*dj - 3*dk + ijk] - 
       chi[-3*dj - 2*dk + ijk] + chi[3*dj - 2*dk + ijk] + 
       chi[-3*dj + 2*dk + ijk] - chi[3*dj + 2*dk + ijk] + 
       chi[-2*dj + 3*dk + ijk] - chi[2*dj + 3*dk + ijk]) - 
    15000.*(chi[2*dj - 4*dk + ijk] + chi[4*dj - 2*dk + ijk] + 
       chi[-4*dj + 2*dk + ijk] + chi[-2*dj + 4*dk + ijk]) + 
    52500.*(-chi[-dj - 4*dk + ijk] + chi[dj - 4*dk + ijk] - 
       chi[-4*dj - dk + ijk] + chi[4*dj - dk + ijk] + 
       chi[-4*dj + dk + ijk] - chi[4*dj + dk + ijk] + 
       chi[-dj + 4*dk + ijk] - chi[dj + 4*dk + ijk]) + 
    15000.*(chi[-2*dj - 4*dk + ijk] + chi[-4*dj - 2*dk + ijk] + 
       chi[4*dj + 2*dk + ijk] + chi[2*dj + 4*dk + ijk]) + 
    3750.*(-chi[-3*dj - 4*dk + ijk] + chi[3*dj - 4*dk + ijk] - 
       chi[-4*dj - 3*dk + ijk] + chi[4*dj - 3*dk + ijk] + 
       chi[-4*dj + 3*dk + ijk] - chi[4*dj + 3*dk + ijk] + 
       chi[-3*dj + 4*dk + ijk] - chi[3*dj + 4*dk + ijk]) - 
    300.*(chi[3*dj - 5*dk + ijk] + chi[5*dj - 3*dk + ijk] + 
       chi[-5*dj + 3*dk + ijk] + chi[-3*dj + 5*dk + ijk]) - 
    4200.*(chi[dj - 5*dk + ijk] + chi[5*dj - dk + ijk] + 
       chi[-5*dj + dk + ijk] + chi[-dj + 5*dk + ijk]) + 
    4200.*(chi[-dj - 5*dk + ijk] + chi[-5*dj - dk + ijk] + 
       chi[5*dj + dk + ijk] + chi[dj + 5*dk + ijk]) + 
    1200.*(-chi[-2*dj - 5*dk + ijk] + chi[2*dj - 5*dk + ijk] - 
       chi[-5*dj - 2*dk + ijk] + chi[5*dj - 2*dk + ijk] + 
       chi[-5*dj + 2*dk + ijk] - chi[5*dj + 2*dk + ijk] + 
       chi[-2*dj + 5*dk + ijk] - chi[2*dj + 5*dk + ijk]) + 
    300.*(chi[-3*dj - 5*dk + ijk] + chi[-5*dj - 3*dk + ijk] + 
       chi[5*dj + 3*dk + ijk] + chi[3*dj + 5*dk + ijk]) + 
    50.*(-chi[-4*dj - 5*dk + ijk] + chi[4*dj - 5*dk + ijk] - 
       chi[-5*dj - 4*dk + ijk] + chi[5*dj - 4*dk + ijk] + 
       chi[-5*dj + 4*dk + ijk] - chi[5*dj + 4*dk + ijk] + 
       chi[-4*dj + 5*dk + ijk] - chi[4*dj + 5*dk + ijk]) - 
    360000.*(chi[2*(dj - dk) + ijk] + chi[2*(-dj + dk) + ijk]) - 
    22500.*(chi[3*(dj - dk) + ijk] + chi[3*(-dj + dk) + ijk]) - 
    625.*(chi[4*(dj - dk) + ijk] + chi[4*(-dj + dk) + ijk]) - 
    4.*(chi[5*(dj - dk) + ijk] + chi[5*(-dj + dk) + ijk]) + 
    360000.*(chi[-2*(dj + dk) + ijk] + chi[2*(dj + dk) + ijk]) + 
    22500.*(chi[-3*(dj + dk) + ijk] + chi[3*(dj + dk) + ijk]) + 
    625.*(chi[-4*(dj + dk) + ijk] + chi[4*(dj + dk) + ijk]) + 
    4.*(chi[-5*(dj + dk) + ijk] + chi[5*(dj + dk) + ijk]))
;

ddchi33
=
0.00003968253968253968254*oodz2*
  (-73766.*chi[ijk] + 42000.*(chi[-dk + ijk] + chi[dk + ijk]) - 
    6000.*(chi[-2*dk + ijk] + chi[2*dk + ijk]) + 
    1000.*(chi[-3*dk + ijk] + chi[3*dk + ijk]) - 
    125.*(chi[-4*dk + ijk] + chi[4*dk + ijk]) + 
    8.*(chi[-5*dk + ijk] + chi[5*dk + ijk]))
;


#endif 


} 

w1
=
2.*oo2dx*beta1[ijk]
;

w2
=
2.*oo2dy*beta2[ijk]
;

w3
=
2.*oo2dz*beta3[ijk]
;


if (order_advection == 2 || boundaryNaway(1)) { 


set_advection2(w1, w2, w3); 

lieg11
=
advu0*g11[advi0] + advu1*g11[advi1] + advu2*g11[advi2] + advv0*g11[advj0] + 
  advv1*g11[advj1] + advv2*g11[advj2] + advw0*g11[advk0] + 
  advw1*g11[advk1] + advw2*g11[advk2]
;

lieg12
=
advu0*g12[advi0] + advu1*g12[advi1] + advu2*g12[advi2] + advv0*g12[advj0] + 
  advv1*g12[advj1] + advv2*g12[advj2] + advw0*g12[advk0] + 
  advw1*g12[advk1] + advw2*g12[advk2]
;

lieg13
=
advu0*g13[advi0] + advu1*g13[advi1] + advu2*g13[advi2] + advv0*g13[advj0] + 
  advv1*g13[advj1] + advv2*g13[advj2] + advw0*g13[advk0] + 
  advw1*g13[advk1] + advw2*g13[advk2]
;

lieg22
=
advu0*g22[advi0] + advu1*g22[advi1] + advu2*g22[advi2] + advv0*g22[advj0] + 
  advv1*g22[advj1] + advv2*g22[advj2] + advw0*g22[advk0] + 
  advw1*g22[advk1] + advw2*g22[advk2]
;

lieg23
=
advu0*g23[advi0] + advu1*g23[advi1] + advu2*g23[advi2] + advv0*g23[advj0] + 
  advv1*g23[advj1] + advv2*g23[advj2] + advw0*g23[advk0] + 
  advw1*g23[advk1] + advw2*g23[advk2]
;

lieg33
=
advu0*g33[advi0] + advu1*g33[advi1] + advu2*g33[advi2] + advv0*g33[advj0] + 
  advv1*g33[advj1] + advv2*g33[advj2] + advw0*g33[advk0] + 
  advw1*g33[advk1] + advw2*g33[advk2]
;

lieA11
=
advu0*A11[advi0] + advu1*A11[advi1] + advu2*A11[advi2] + advv0*A11[advj0] + 
  advv1*A11[advj1] + advv2*A11[advj2] + advw0*A11[advk0] + 
  advw1*A11[advk1] + advw2*A11[advk2]
;

lieA12
=
advu0*A12[advi0] + advu1*A12[advi1] + advu2*A12[advi2] + advv0*A12[advj0] + 
  advv1*A12[advj1] + advv2*A12[advj2] + advw0*A12[advk0] + 
  advw1*A12[advk1] + advw2*A12[advk2]
;

lieA13
=
advu0*A13[advi0] + advu1*A13[advi1] + advu2*A13[advi2] + advv0*A13[advj0] + 
  advv1*A13[advj1] + advv2*A13[advj2] + advw0*A13[advk0] + 
  advw1*A13[advk1] + advw2*A13[advk2]
;

lieA22
=
advu0*A22[advi0] + advu1*A22[advi1] + advu2*A22[advi2] + advv0*A22[advj0] + 
  advv1*A22[advj1] + advv2*A22[advj2] + advw0*A22[advk0] + 
  advw1*A22[advk1] + advw2*A22[advk2]
;

lieA23
=
advu0*A23[advi0] + advu1*A23[advi1] + advu2*A23[advi2] + advv0*A23[advj0] + 
  advv1*A23[advj1] + advv2*A23[advj2] + advw0*A23[advk0] + 
  advw1*A23[advk1] + advw2*A23[advk2]
;

lieA33
=
advu0*A33[advi0] + advu1*A33[advi1] + advu2*A33[advi2] + advv0*A33[advj0] + 
  advv1*A33[advj1] + advv2*A33[advj2] + advw0*A33[advk0] + 
  advw1*A33[advk1] + advw2*A33[advk2]
;

lieK
=
advu0*K[advi0] + advu1*K[advi1] + advu2*K[advi2] + advv0*K[advj0] + 
  advv1*K[advj1] + advv2*K[advj2] + advw0*K[advk0] + advw1*K[advk1] + 
  advw2*K[advk2]
;

liechi
=
advu0*chi[advi0] + advu1*chi[advi1] + advu2*chi[advi2] + advv0*chi[advj0] + 
  advv1*chi[advj1] + advv2*chi[advj2] + advw0*chi[advk0] + 
  advw1*chi[advk1] + advw2*chi[advk2]
;

liealpha
=
advu0*alpha[advi0] + advu1*alpha[advi1] + advu2*alpha[advi2] + 
  advv0*alpha[advj0] + advv1*alpha[advj1] + advv2*alpha[advj2] + 
  advw0*alpha[advk0] + advw1*alpha[advk1] + advw2*alpha[advk2]
;

advG1
=
advu0*G1[advi0] + advu1*G1[advi1] + advu2*G1[advi2] + advv0*G1[advj0] + 
  advv1*G1[advj1] + advv2*G1[advj2] + advw0*G1[advk0] + advw1*G1[advk1] + 
  advw2*G1[advk2]
;

advG2
=
advu0*G2[advi0] + advu1*G2[advi1] + advu2*G2[advi2] + advv0*G2[advj0] + 
  advv1*G2[advj1] + advv2*G2[advj2] + advw0*G2[advk0] + advw1*G2[advk1] + 
  advw2*G2[advk2]
;

advG3
=
advu0*G3[advi0] + advu1*G3[advi1] + advu2*G3[advi2] + advv0*G3[advj0] + 
  advv1*G3[advj1] + advv2*G3[advj2] + advw0*G3[advk0] + advw1*G3[advk1] + 
  advw2*G3[advk2]
;

advbeta1
=
advu0*beta1[advi0] + advu1*beta1[advi1] + advu2*beta1[advi2] + 
  advv0*beta1[advj0] + advv1*beta1[advj1] + advv2*beta1[advj2] + 
  advw0*beta1[advk0] + advw1*beta1[advk1] + advw2*beta1[advk2]
;

advbeta2
=
advu0*beta2[advi0] + advu1*beta2[advi1] + advu2*beta2[advi2] + 
  advv0*beta2[advj0] + advv1*beta2[advj1] + advv2*beta2[advj2] + 
  advw0*beta2[advk0] + advw1*beta2[advk1] + advw2*beta2[advk2]
;

advbeta3
=
advu0*beta3[advi0] + advu1*beta3[advi1] + advu2*beta3[advi2] + 
  advv0*beta3[advj0] + advv1*beta3[advj1] + advv2*beta3[advj2] + 
  advw0*beta3[advk0] + advw1*beta3[advk1] + advw2*beta3[advk2]
;

advB1
=
advu0*B1[advi0] + advu1*B1[advi1] + advu2*B1[advi2] + advv0*B1[advj0] + 
  advv1*B1[advj1] + advv2*B1[advj2] + advw0*B1[advk0] + advw1*B1[advk1] + 
  advw2*B1[advk2]
;

advB2
=
advu0*B2[advi0] + advu1*B2[advi1] + advu2*B2[advi2] + advv0*B2[advj0] + 
  advv1*B2[advj1] + advv2*B2[advj2] + advw0*B2[advk0] + advw1*B2[advk1] + 
  advw2*B2[advk2]
;

advB3
=
advu0*B3[advi0] + advu1*B3[advi1] + advu2*B3[advi2] + advv0*B3[advj0] + 
  advv1*B3[advj1] + advv2*B3[advj2] + advw0*B3[advk0] + advw1*B3[advk1] + 
  advw2*B3[advk2]
;


} else if (order_advection == 4 || boundaryNaway(2)) { 


#ifdef REDUCEORDERTO2                           
  errorexit("CompilerFlag reduced order to 2");
 #else

set_advection4(w1, w2, w3); 

lieg11
=
advu0*g11[advi0] + advu1*g11[advi1] + advu2*g11[advi2] + advu3*g11[advi3] + 
  advu4*g11[advi4] + advv0*g11[advj0] + advv1*g11[advj1] + 
  advv2*g11[advj2] + advv3*g11[advj3] + advv4*g11[advj4] + 
  advw0*g11[advk0] + advw1*g11[advk1] + advw2*g11[advk2] + 
  advw3*g11[advk3] + advw4*g11[advk4]
;

lieg12
=
advu0*g12[advi0] + advu1*g12[advi1] + advu2*g12[advi2] + advu3*g12[advi3] + 
  advu4*g12[advi4] + advv0*g12[advj0] + advv1*g12[advj1] + 
  advv2*g12[advj2] + advv3*g12[advj3] + advv4*g12[advj4] + 
  advw0*g12[advk0] + advw1*g12[advk1] + advw2*g12[advk2] + 
  advw3*g12[advk3] + advw4*g12[advk4]
;

lieg13
=
advu0*g13[advi0] + advu1*g13[advi1] + advu2*g13[advi2] + advu3*g13[advi3] + 
  advu4*g13[advi4] + advv0*g13[advj0] + advv1*g13[advj1] + 
  advv2*g13[advj2] + advv3*g13[advj3] + advv4*g13[advj4] + 
  advw0*g13[advk0] + advw1*g13[advk1] + advw2*g13[advk2] + 
  advw3*g13[advk3] + advw4*g13[advk4]
;

lieg22
=
advu0*g22[advi0] + advu1*g22[advi1] + advu2*g22[advi2] + advu3*g22[advi3] + 
  advu4*g22[advi4] + advv0*g22[advj0] + advv1*g22[advj1] + 
  advv2*g22[advj2] + advv3*g22[advj3] + advv4*g22[advj4] + 
  advw0*g22[advk0] + advw1*g22[advk1] + advw2*g22[advk2] + 
  advw3*g22[advk3] + advw4*g22[advk4]
;

lieg23
=
advu0*g23[advi0] + advu1*g23[advi1] + advu2*g23[advi2] + advu3*g23[advi3] + 
  advu4*g23[advi4] + advv0*g23[advj0] + advv1*g23[advj1] + 
  advv2*g23[advj2] + advv3*g23[advj3] + advv4*g23[advj4] + 
  advw0*g23[advk0] + advw1*g23[advk1] + advw2*g23[advk2] + 
  advw3*g23[advk3] + advw4*g23[advk4]
;

lieg33
=
advu0*g33[advi0] + advu1*g33[advi1] + advu2*g33[advi2] + advu3*g33[advi3] + 
  advu4*g33[advi4] + advv0*g33[advj0] + advv1*g33[advj1] + 
  advv2*g33[advj2] + advv3*g33[advj3] + advv4*g33[advj4] + 
  advw0*g33[advk0] + advw1*g33[advk1] + advw2*g33[advk2] + 
  advw3*g33[advk3] + advw4*g33[advk4]
;

lieA11
=
advu0*A11[advi0] + advu1*A11[advi1] + advu2*A11[advi2] + advu3*A11[advi3] + 
  advu4*A11[advi4] + advv0*A11[advj0] + advv1*A11[advj1] + 
  advv2*A11[advj2] + advv3*A11[advj3] + advv4*A11[advj4] + 
  advw0*A11[advk0] + advw1*A11[advk1] + advw2*A11[advk2] + 
  advw3*A11[advk3] + advw4*A11[advk4]
;

lieA12
=
advu0*A12[advi0] + advu1*A12[advi1] + advu2*A12[advi2] + advu3*A12[advi3] + 
  advu4*A12[advi4] + advv0*A12[advj0] + advv1*A12[advj1] + 
  advv2*A12[advj2] + advv3*A12[advj3] + advv4*A12[advj4] + 
  advw0*A12[advk0] + advw1*A12[advk1] + advw2*A12[advk2] + 
  advw3*A12[advk3] + advw4*A12[advk4]
;

lieA13
=
advu0*A13[advi0] + advu1*A13[advi1] + advu2*A13[advi2] + advu3*A13[advi3] + 
  advu4*A13[advi4] + advv0*A13[advj0] + advv1*A13[advj1] + 
  advv2*A13[advj2] + advv3*A13[advj3] + advv4*A13[advj4] + 
  advw0*A13[advk0] + advw1*A13[advk1] + advw2*A13[advk2] + 
  advw3*A13[advk3] + advw4*A13[advk4]
;

lieA22
=
advu0*A22[advi0] + advu1*A22[advi1] + advu2*A22[advi2] + advu3*A22[advi3] + 
  advu4*A22[advi4] + advv0*A22[advj0] + advv1*A22[advj1] + 
  advv2*A22[advj2] + advv3*A22[advj3] + advv4*A22[advj4] + 
  advw0*A22[advk0] + advw1*A22[advk1] + advw2*A22[advk2] + 
  advw3*A22[advk3] + advw4*A22[advk4]
;

lieA23
=
advu0*A23[advi0] + advu1*A23[advi1] + advu2*A23[advi2] + advu3*A23[advi3] + 
  advu4*A23[advi4] + advv0*A23[advj0] + advv1*A23[advj1] + 
  advv2*A23[advj2] + advv3*A23[advj3] + advv4*A23[advj4] + 
  advw0*A23[advk0] + advw1*A23[advk1] + advw2*A23[advk2] + 
  advw3*A23[advk3] + advw4*A23[advk4]
;

lieA33
=
advu0*A33[advi0] + advu1*A33[advi1] + advu2*A33[advi2] + advu3*A33[advi3] + 
  advu4*A33[advi4] + advv0*A33[advj0] + advv1*A33[advj1] + 
  advv2*A33[advj2] + advv3*A33[advj3] + advv4*A33[advj4] + 
  advw0*A33[advk0] + advw1*A33[advk1] + advw2*A33[advk2] + 
  advw3*A33[advk3] + advw4*A33[advk4]
;

lieK
=
advu0*K[advi0] + advu1*K[advi1] + advu2*K[advi2] + advu3*K[advi3] + 
  advu4*K[advi4] + advv0*K[advj0] + advv1*K[advj1] + advv2*K[advj2] + 
  advv3*K[advj3] + advv4*K[advj4] + advw0*K[advk0] + advw1*K[advk1] + 
  advw2*K[advk2] + advw3*K[advk3] + advw4*K[advk4]
;

liechi
=
advu0*chi[advi0] + advu1*chi[advi1] + advu2*chi[advi2] + advu3*chi[advi3] + 
  advu4*chi[advi4] + advv0*chi[advj0] + advv1*chi[advj1] + 
  advv2*chi[advj2] + advv3*chi[advj3] + advv4*chi[advj4] + 
  advw0*chi[advk0] + advw1*chi[advk1] + advw2*chi[advk2] + 
  advw3*chi[advk3] + advw4*chi[advk4]
;

liealpha
=
advu0*alpha[advi0] + advu1*alpha[advi1] + advu2*alpha[advi2] + 
  advu3*alpha[advi3] + advu4*alpha[advi4] + advv0*alpha[advj0] + 
  advv1*alpha[advj1] + advv2*alpha[advj2] + advv3*alpha[advj3] + 
  advv4*alpha[advj4] + advw0*alpha[advk0] + advw1*alpha[advk1] + 
  advw2*alpha[advk2] + advw3*alpha[advk3] + advw4*alpha[advk4]
;

advG1
=
advu0*G1[advi0] + advu1*G1[advi1] + advu2*G1[advi2] + advu3*G1[advi3] + 
  advu4*G1[advi4] + advv0*G1[advj0] + advv1*G1[advj1] + advv2*G1[advj2] + 
  advv3*G1[advj3] + advv4*G1[advj4] + advw0*G1[advk0] + advw1*G1[advk1] + 
  advw2*G1[advk2] + advw3*G1[advk3] + advw4*G1[advk4]
;

advG2
=
advu0*G2[advi0] + advu1*G2[advi1] + advu2*G2[advi2] + advu3*G2[advi3] + 
  advu4*G2[advi4] + advv0*G2[advj0] + advv1*G2[advj1] + advv2*G2[advj2] + 
  advv3*G2[advj3] + advv4*G2[advj4] + advw0*G2[advk0] + advw1*G2[advk1] + 
  advw2*G2[advk2] + advw3*G2[advk3] + advw4*G2[advk4]
;

advG3
=
advu0*G3[advi0] + advu1*G3[advi1] + advu2*G3[advi2] + advu3*G3[advi3] + 
  advu4*G3[advi4] + advv0*G3[advj0] + advv1*G3[advj1] + advv2*G3[advj2] + 
  advv3*G3[advj3] + advv4*G3[advj4] + advw0*G3[advk0] + advw1*G3[advk1] + 
  advw2*G3[advk2] + advw3*G3[advk3] + advw4*G3[advk4]
;

advbeta1
=
advu0*beta1[advi0] + advu1*beta1[advi1] + advu2*beta1[advi2] + 
  advu3*beta1[advi3] + advu4*beta1[advi4] + advv0*beta1[advj0] + 
  advv1*beta1[advj1] + advv2*beta1[advj2] + advv3*beta1[advj3] + 
  advv4*beta1[advj4] + advw0*beta1[advk0] + advw1*beta1[advk1] + 
  advw2*beta1[advk2] + advw3*beta1[advk3] + advw4*beta1[advk4]
;

advbeta2
=
advu0*beta2[advi0] + advu1*beta2[advi1] + advu2*beta2[advi2] + 
  advu3*beta2[advi3] + advu4*beta2[advi4] + advv0*beta2[advj0] + 
  advv1*beta2[advj1] + advv2*beta2[advj2] + advv3*beta2[advj3] + 
  advv4*beta2[advj4] + advw0*beta2[advk0] + advw1*beta2[advk1] + 
  advw2*beta2[advk2] + advw3*beta2[advk3] + advw4*beta2[advk4]
;

advbeta3
=
advu0*beta3[advi0] + advu1*beta3[advi1] + advu2*beta3[advi2] + 
  advu3*beta3[advi3] + advu4*beta3[advi4] + advv0*beta3[advj0] + 
  advv1*beta3[advj1] + advv2*beta3[advj2] + advv3*beta3[advj3] + 
  advv4*beta3[advj4] + advw0*beta3[advk0] + advw1*beta3[advk1] + 
  advw2*beta3[advk2] + advw3*beta3[advk3] + advw4*beta3[advk4]
;

advB1
=
advu0*B1[advi0] + advu1*B1[advi1] + advu2*B1[advi2] + advu3*B1[advi3] + 
  advu4*B1[advi4] + advv0*B1[advj0] + advv1*B1[advj1] + advv2*B1[advj2] + 
  advv3*B1[advj3] + advv4*B1[advj4] + advw0*B1[advk0] + advw1*B1[advk1] + 
  advw2*B1[advk2] + advw3*B1[advk3] + advw4*B1[advk4]
;

advB2
=
advu0*B2[advi0] + advu1*B2[advi1] + advu2*B2[advi2] + advu3*B2[advi3] + 
  advu4*B2[advi4] + advv0*B2[advj0] + advv1*B2[advj1] + advv2*B2[advj2] + 
  advv3*B2[advj3] + advv4*B2[advj4] + advw0*B2[advk0] + advw1*B2[advk1] + 
  advw2*B2[advk2] + advw3*B2[advk3] + advw4*B2[advk4]
;

advB3
=
advu0*B3[advi0] + advu1*B3[advi1] + advu2*B3[advi2] + advu3*B3[advi3] + 
  advu4*B3[advi4] + advv0*B3[advj0] + advv1*B3[advj1] + advv2*B3[advj2] + 
  advv3*B3[advj3] + advv4*B3[advj4] + advw0*B3[advk0] + advw1*B3[advk1] + 
  advw2*B3[advk2] + advw3*B3[advk3] + advw4*B3[advk4]
;


#endif 


} else if (order_advection == 6 || boundaryNaway(3)) { 


#ifdef REDUCEORDERTO4                           
  errorexit("CompilerFlag reduced order to 4");
 #else

set_advection6(w1, w2, w3); 

lieg11
=
advu0*g11[advi0] + advu1*g11[advi1] + advu2*g11[advi2] + advu3*g11[advi3] + 
  advu4*g11[advi4] + advu5*g11[advi5] + advu6*g11[advi6] + 
  advv0*g11[advj0] + advv1*g11[advj1] + advv2*g11[advj2] + 
  advv3*g11[advj3] + advv4*g11[advj4] + advv5*g11[advj5] + 
  advv6*g11[advj6] + advw0*g11[advk0] + advw1*g11[advk1] + 
  advw2*g11[advk2] + advw3*g11[advk3] + advw4*g11[advk4] + 
  advw5*g11[advk5] + advw6*g11[advk6]
;

lieg12
=
advu0*g12[advi0] + advu1*g12[advi1] + advu2*g12[advi2] + advu3*g12[advi3] + 
  advu4*g12[advi4] + advu5*g12[advi5] + advu6*g12[advi6] + 
  advv0*g12[advj0] + advv1*g12[advj1] + advv2*g12[advj2] + 
  advv3*g12[advj3] + advv4*g12[advj4] + advv5*g12[advj5] + 
  advv6*g12[advj6] + advw0*g12[advk0] + advw1*g12[advk1] + 
  advw2*g12[advk2] + advw3*g12[advk3] + advw4*g12[advk4] + 
  advw5*g12[advk5] + advw6*g12[advk6]
;

lieg13
=
advu0*g13[advi0] + advu1*g13[advi1] + advu2*g13[advi2] + advu3*g13[advi3] + 
  advu4*g13[advi4] + advu5*g13[advi5] + advu6*g13[advi6] + 
  advv0*g13[advj0] + advv1*g13[advj1] + advv2*g13[advj2] + 
  advv3*g13[advj3] + advv4*g13[advj4] + advv5*g13[advj5] + 
  advv6*g13[advj6] + advw0*g13[advk0] + advw1*g13[advk1] + 
  advw2*g13[advk2] + advw3*g13[advk3] + advw4*g13[advk4] + 
  advw5*g13[advk5] + advw6*g13[advk6]
;

lieg22
=
advu0*g22[advi0] + advu1*g22[advi1] + advu2*g22[advi2] + advu3*g22[advi3] + 
  advu4*g22[advi4] + advu5*g22[advi5] + advu6*g22[advi6] + 
  advv0*g22[advj0] + advv1*g22[advj1] + advv2*g22[advj2] + 
  advv3*g22[advj3] + advv4*g22[advj4] + advv5*g22[advj5] + 
  advv6*g22[advj6] + advw0*g22[advk0] + advw1*g22[advk1] + 
  advw2*g22[advk2] + advw3*g22[advk3] + advw4*g22[advk4] + 
  advw5*g22[advk5] + advw6*g22[advk6]
;

lieg23
=
advu0*g23[advi0] + advu1*g23[advi1] + advu2*g23[advi2] + advu3*g23[advi3] + 
  advu4*g23[advi4] + advu5*g23[advi5] + advu6*g23[advi6] + 
  advv0*g23[advj0] + advv1*g23[advj1] + advv2*g23[advj2] + 
  advv3*g23[advj3] + advv4*g23[advj4] + advv5*g23[advj5] + 
  advv6*g23[advj6] + advw0*g23[advk0] + advw1*g23[advk1] + 
  advw2*g23[advk2] + advw3*g23[advk3] + advw4*g23[advk4] + 
  advw5*g23[advk5] + advw6*g23[advk6]
;

lieg33
=
advu0*g33[advi0] + advu1*g33[advi1] + advu2*g33[advi2] + advu3*g33[advi3] + 
  advu4*g33[advi4] + advu5*g33[advi5] + advu6*g33[advi6] + 
  advv0*g33[advj0] + advv1*g33[advj1] + advv2*g33[advj2] + 
  advv3*g33[advj3] + advv4*g33[advj4] + advv5*g33[advj5] + 
  advv6*g33[advj6] + advw0*g33[advk0] + advw1*g33[advk1] + 
  advw2*g33[advk2] + advw3*g33[advk3] + advw4*g33[advk4] + 
  advw5*g33[advk5] + advw6*g33[advk6]
;

lieA11
=
advu0*A11[advi0] + advu1*A11[advi1] + advu2*A11[advi2] + advu3*A11[advi3] + 
  advu4*A11[advi4] + advu5*A11[advi5] + advu6*A11[advi6] + 
  advv0*A11[advj0] + advv1*A11[advj1] + advv2*A11[advj2] + 
  advv3*A11[advj3] + advv4*A11[advj4] + advv5*A11[advj5] + 
  advv6*A11[advj6] + advw0*A11[advk0] + advw1*A11[advk1] + 
  advw2*A11[advk2] + advw3*A11[advk3] + advw4*A11[advk4] + 
  advw5*A11[advk5] + advw6*A11[advk6]
;

lieA12
=
advu0*A12[advi0] + advu1*A12[advi1] + advu2*A12[advi2] + advu3*A12[advi3] + 
  advu4*A12[advi4] + advu5*A12[advi5] + advu6*A12[advi6] + 
  advv0*A12[advj0] + advv1*A12[advj1] + advv2*A12[advj2] + 
  advv3*A12[advj3] + advv4*A12[advj4] + advv5*A12[advj5] + 
  advv6*A12[advj6] + advw0*A12[advk0] + advw1*A12[advk1] + 
  advw2*A12[advk2] + advw3*A12[advk3] + advw4*A12[advk4] + 
  advw5*A12[advk5] + advw6*A12[advk6]
;

lieA13
=
advu0*A13[advi0] + advu1*A13[advi1] + advu2*A13[advi2] + advu3*A13[advi3] + 
  advu4*A13[advi4] + advu5*A13[advi5] + advu6*A13[advi6] + 
  advv0*A13[advj0] + advv1*A13[advj1] + advv2*A13[advj2] + 
  advv3*A13[advj3] + advv4*A13[advj4] + advv5*A13[advj5] + 
  advv6*A13[advj6] + advw0*A13[advk0] + advw1*A13[advk1] + 
  advw2*A13[advk2] + advw3*A13[advk3] + advw4*A13[advk4] + 
  advw5*A13[advk5] + advw6*A13[advk6]
;

lieA22
=
advu0*A22[advi0] + advu1*A22[advi1] + advu2*A22[advi2] + advu3*A22[advi3] + 
  advu4*A22[advi4] + advu5*A22[advi5] + advu6*A22[advi6] + 
  advv0*A22[advj0] + advv1*A22[advj1] + advv2*A22[advj2] + 
  advv3*A22[advj3] + advv4*A22[advj4] + advv5*A22[advj5] + 
  advv6*A22[advj6] + advw0*A22[advk0] + advw1*A22[advk1] + 
  advw2*A22[advk2] + advw3*A22[advk3] + advw4*A22[advk4] + 
  advw5*A22[advk5] + advw6*A22[advk6]
;

lieA23
=
advu0*A23[advi0] + advu1*A23[advi1] + advu2*A23[advi2] + advu3*A23[advi3] + 
  advu4*A23[advi4] + advu5*A23[advi5] + advu6*A23[advi6] + 
  advv0*A23[advj0] + advv1*A23[advj1] + advv2*A23[advj2] + 
  advv3*A23[advj3] + advv4*A23[advj4] + advv5*A23[advj5] + 
  advv6*A23[advj6] + advw0*A23[advk0] + advw1*A23[advk1] + 
  advw2*A23[advk2] + advw3*A23[advk3] + advw4*A23[advk4] + 
  advw5*A23[advk5] + advw6*A23[advk6]
;

lieA33
=
advu0*A33[advi0] + advu1*A33[advi1] + advu2*A33[advi2] + advu3*A33[advi3] + 
  advu4*A33[advi4] + advu5*A33[advi5] + advu6*A33[advi6] + 
  advv0*A33[advj0] + advv1*A33[advj1] + advv2*A33[advj2] + 
  advv3*A33[advj3] + advv4*A33[advj4] + advv5*A33[advj5] + 
  advv6*A33[advj6] + advw0*A33[advk0] + advw1*A33[advk1] + 
  advw2*A33[advk2] + advw3*A33[advk3] + advw4*A33[advk4] + 
  advw5*A33[advk5] + advw6*A33[advk6]
;

lieK
=
advu0*K[advi0] + advu1*K[advi1] + advu2*K[advi2] + advu3*K[advi3] + 
  advu4*K[advi4] + advu5*K[advi5] + advu6*K[advi6] + advv0*K[advj0] + 
  advv1*K[advj1] + advv2*K[advj2] + advv3*K[advj3] + advv4*K[advj4] + 
  advv5*K[advj5] + advv6*K[advj6] + advw0*K[advk0] + advw1*K[advk1] + 
  advw2*K[advk2] + advw3*K[advk3] + advw4*K[advk4] + advw5*K[advk5] + 
  advw6*K[advk6]
;

liechi
=
advu0*chi[advi0] + advu1*chi[advi1] + advu2*chi[advi2] + advu3*chi[advi3] + 
  advu4*chi[advi4] + advu5*chi[advi5] + advu6*chi[advi6] + 
  advv0*chi[advj0] + advv1*chi[advj1] + advv2*chi[advj2] + 
  advv3*chi[advj3] + advv4*chi[advj4] + advv5*chi[advj5] + 
  advv6*chi[advj6] + advw0*chi[advk0] + advw1*chi[advk1] + 
  advw2*chi[advk2] + advw3*chi[advk3] + advw4*chi[advk4] + 
  advw5*chi[advk5] + advw6*chi[advk6]
;

liealpha
=
advu0*alpha[advi0] + advu1*alpha[advi1] + advu2*alpha[advi2] + 
  advu3*alpha[advi3] + advu4*alpha[advi4] + advu5*alpha[advi5] + 
  advu6*alpha[advi6] + advv0*alpha[advj0] + advv1*alpha[advj1] + 
  advv2*alpha[advj2] + advv3*alpha[advj3] + advv4*alpha[advj4] + 
  advv5*alpha[advj5] + advv6*alpha[advj6] + advw0*alpha[advk0] + 
  advw1*alpha[advk1] + advw2*alpha[advk2] + advw3*alpha[advk3] + 
  advw4*alpha[advk4] + advw5*alpha[advk5] + advw6*alpha[advk6]
;

advG1
=
advu0*G1[advi0] + advu1*G1[advi1] + advu2*G1[advi2] + advu3*G1[advi3] + 
  advu4*G1[advi4] + advu5*G1[advi5] + advu6*G1[advi6] + advv0*G1[advj0] + 
  advv1*G1[advj1] + advv2*G1[advj2] + advv3*G1[advj3] + advv4*G1[advj4] + 
  advv5*G1[advj5] + advv6*G1[advj6] + advw0*G1[advk0] + advw1*G1[advk1] + 
  advw2*G1[advk2] + advw3*G1[advk3] + advw4*G1[advk4] + advw5*G1[advk5] + 
  advw6*G1[advk6]
;

advG2
=
advu0*G2[advi0] + advu1*G2[advi1] + advu2*G2[advi2] + advu3*G2[advi3] + 
  advu4*G2[advi4] + advu5*G2[advi5] + advu6*G2[advi6] + advv0*G2[advj0] + 
  advv1*G2[advj1] + advv2*G2[advj2] + advv3*G2[advj3] + advv4*G2[advj4] + 
  advv5*G2[advj5] + advv6*G2[advj6] + advw0*G2[advk0] + advw1*G2[advk1] + 
  advw2*G2[advk2] + advw3*G2[advk3] + advw4*G2[advk4] + advw5*G2[advk5] + 
  advw6*G2[advk6]
;

advG3
=
advu0*G3[advi0] + advu1*G3[advi1] + advu2*G3[advi2] + advu3*G3[advi3] + 
  advu4*G3[advi4] + advu5*G3[advi5] + advu6*G3[advi6] + advv0*G3[advj0] + 
  advv1*G3[advj1] + advv2*G3[advj2] + advv3*G3[advj3] + advv4*G3[advj4] + 
  advv5*G3[advj5] + advv6*G3[advj6] + advw0*G3[advk0] + advw1*G3[advk1] + 
  advw2*G3[advk2] + advw3*G3[advk3] + advw4*G3[advk4] + advw5*G3[advk5] + 
  advw6*G3[advk6]
;

advbeta1
=
advu0*beta1[advi0] + advu1*beta1[advi1] + advu2*beta1[advi2] + 
  advu3*beta1[advi3] + advu4*beta1[advi4] + advu5*beta1[advi5] + 
  advu6*beta1[advi6] + advv0*beta1[advj0] + advv1*beta1[advj1] + 
  advv2*beta1[advj2] + advv3*beta1[advj3] + advv4*beta1[advj4] + 
  advv5*beta1[advj5] + advv6*beta1[advj6] + advw0*beta1[advk0] + 
  advw1*beta1[advk1] + advw2*beta1[advk2] + advw3*beta1[advk3] + 
  advw4*beta1[advk4] + advw5*beta1[advk5] + advw6*beta1[advk6]
;

advbeta2
=
advu0*beta2[advi0] + advu1*beta2[advi1] + advu2*beta2[advi2] + 
  advu3*beta2[advi3] + advu4*beta2[advi4] + advu5*beta2[advi5] + 
  advu6*beta2[advi6] + advv0*beta2[advj0] + advv1*beta2[advj1] + 
  advv2*beta2[advj2] + advv3*beta2[advj3] + advv4*beta2[advj4] + 
  advv5*beta2[advj5] + advv6*beta2[advj6] + advw0*beta2[advk0] + 
  advw1*beta2[advk1] + advw2*beta2[advk2] + advw3*beta2[advk3] + 
  advw4*beta2[advk4] + advw5*beta2[advk5] + advw6*beta2[advk6]
;

advbeta3
=
advu0*beta3[advi0] + advu1*beta3[advi1] + advu2*beta3[advi2] + 
  advu3*beta3[advi3] + advu4*beta3[advi4] + advu5*beta3[advi5] + 
  advu6*beta3[advi6] + advv0*beta3[advj0] + advv1*beta3[advj1] + 
  advv2*beta3[advj2] + advv3*beta3[advj3] + advv4*beta3[advj4] + 
  advv5*beta3[advj5] + advv6*beta3[advj6] + advw0*beta3[advk0] + 
  advw1*beta3[advk1] + advw2*beta3[advk2] + advw3*beta3[advk3] + 
  advw4*beta3[advk4] + advw5*beta3[advk5] + advw6*beta3[advk6]
;

advB1
=
advu0*B1[advi0] + advu1*B1[advi1] + advu2*B1[advi2] + advu3*B1[advi3] + 
  advu4*B1[advi4] + advu5*B1[advi5] + advu6*B1[advi6] + advv0*B1[advj0] + 
  advv1*B1[advj1] + advv2*B1[advj2] + advv3*B1[advj3] + advv4*B1[advj4] + 
  advv5*B1[advj5] + advv6*B1[advj6] + advw0*B1[advk0] + advw1*B1[advk1] + 
  advw2*B1[advk2] + advw3*B1[advk3] + advw4*B1[advk4] + advw5*B1[advk5] + 
  advw6*B1[advk6]
;

advB2
=
advu0*B2[advi0] + advu1*B2[advi1] + advu2*B2[advi2] + advu3*B2[advi3] + 
  advu4*B2[advi4] + advu5*B2[advi5] + advu6*B2[advi6] + advv0*B2[advj0] + 
  advv1*B2[advj1] + advv2*B2[advj2] + advv3*B2[advj3] + advv4*B2[advj4] + 
  advv5*B2[advj5] + advv6*B2[advj6] + advw0*B2[advk0] + advw1*B2[advk1] + 
  advw2*B2[advk2] + advw3*B2[advk3] + advw4*B2[advk4] + advw5*B2[advk5] + 
  advw6*B2[advk6]
;

advB3
=
advu0*B3[advi0] + advu1*B3[advi1] + advu2*B3[advi2] + advu3*B3[advi3] + 
  advu4*B3[advi4] + advu5*B3[advi5] + advu6*B3[advi6] + advv0*B3[advj0] + 
  advv1*B3[advj1] + advv2*B3[advj2] + advv3*B3[advj3] + advv4*B3[advj4] + 
  advv5*B3[advj5] + advv6*B3[advj6] + advw0*B3[advk0] + advw1*B3[advk1] + 
  advw2*B3[advk2] + advw3*B3[advk3] + advw4*B3[advk4] + advw5*B3[advk5] + 
  advw6*B3[advk6]
;


#endif 


} else if (order_advection == 8 || boundaryNaway(4)) { 


#ifdef REDUCEORDERTO6                           
  errorexit("CompilerFlag reduced order to 6");
 #else

set_advection8(w1, w2, w3); 

lieg11
=
advu0*g11[advi0] + advu1*g11[advi1] + advu2*g11[advi2] + advu3*g11[advi3] + 
  advu4*g11[advi4] + advu5*g11[advi5] + advu6*g11[advi6] + 
  advu7*g11[advi7] + advu8*g11[advi8] + advv0*g11[advj0] + 
  advv1*g11[advj1] + advv2*g11[advj2] + advv3*g11[advj3] + 
  advv4*g11[advj4] + advv5*g11[advj5] + advv6*g11[advj6] + 
  advv7*g11[advj7] + advv8*g11[advj8] + advw0*g11[advk0] + 
  advw1*g11[advk1] + advw2*g11[advk2] + advw3*g11[advk3] + 
  advw4*g11[advk4] + advw5*g11[advk5] + advw6*g11[advk6] + 
  advw7*g11[advk7] + advw8*g11[advk8]
;

lieg12
=
advu0*g12[advi0] + advu1*g12[advi1] + advu2*g12[advi2] + advu3*g12[advi3] + 
  advu4*g12[advi4] + advu5*g12[advi5] + advu6*g12[advi6] + 
  advu7*g12[advi7] + advu8*g12[advi8] + advv0*g12[advj0] + 
  advv1*g12[advj1] + advv2*g12[advj2] + advv3*g12[advj3] + 
  advv4*g12[advj4] + advv5*g12[advj5] + advv6*g12[advj6] + 
  advv7*g12[advj7] + advv8*g12[advj8] + advw0*g12[advk0] + 
  advw1*g12[advk1] + advw2*g12[advk2] + advw3*g12[advk3] + 
  advw4*g12[advk4] + advw5*g12[advk5] + advw6*g12[advk6] + 
  advw7*g12[advk7] + advw8*g12[advk8]
;

lieg13
=
advu0*g13[advi0] + advu1*g13[advi1] + advu2*g13[advi2] + advu3*g13[advi3] + 
  advu4*g13[advi4] + advu5*g13[advi5] + advu6*g13[advi6] + 
  advu7*g13[advi7] + advu8*g13[advi8] + advv0*g13[advj0] + 
  advv1*g13[advj1] + advv2*g13[advj2] + advv3*g13[advj3] + 
  advv4*g13[advj4] + advv5*g13[advj5] + advv6*g13[advj6] + 
  advv7*g13[advj7] + advv8*g13[advj8] + advw0*g13[advk0] + 
  advw1*g13[advk1] + advw2*g13[advk2] + advw3*g13[advk3] + 
  advw4*g13[advk4] + advw5*g13[advk5] + advw6*g13[advk6] + 
  advw7*g13[advk7] + advw8*g13[advk8]
;

lieg22
=
advu0*g22[advi0] + advu1*g22[advi1] + advu2*g22[advi2] + advu3*g22[advi3] + 
  advu4*g22[advi4] + advu5*g22[advi5] + advu6*g22[advi6] + 
  advu7*g22[advi7] + advu8*g22[advi8] + advv0*g22[advj0] + 
  advv1*g22[advj1] + advv2*g22[advj2] + advv3*g22[advj3] + 
  advv4*g22[advj4] + advv5*g22[advj5] + advv6*g22[advj6] + 
  advv7*g22[advj7] + advv8*g22[advj8] + advw0*g22[advk0] + 
  advw1*g22[advk1] + advw2*g22[advk2] + advw3*g22[advk3] + 
  advw4*g22[advk4] + advw5*g22[advk5] + advw6*g22[advk6] + 
  advw7*g22[advk7] + advw8*g22[advk8]
;

lieg23
=
advu0*g23[advi0] + advu1*g23[advi1] + advu2*g23[advi2] + advu3*g23[advi3] + 
  advu4*g23[advi4] + advu5*g23[advi5] + advu6*g23[advi6] + 
  advu7*g23[advi7] + advu8*g23[advi8] + advv0*g23[advj0] + 
  advv1*g23[advj1] + advv2*g23[advj2] + advv3*g23[advj3] + 
  advv4*g23[advj4] + advv5*g23[advj5] + advv6*g23[advj6] + 
  advv7*g23[advj7] + advv8*g23[advj8] + advw0*g23[advk0] + 
  advw1*g23[advk1] + advw2*g23[advk2] + advw3*g23[advk3] + 
  advw4*g23[advk4] + advw5*g23[advk5] + advw6*g23[advk6] + 
  advw7*g23[advk7] + advw8*g23[advk8]
;

lieg33
=
advu0*g33[advi0] + advu1*g33[advi1] + advu2*g33[advi2] + advu3*g33[advi3] + 
  advu4*g33[advi4] + advu5*g33[advi5] + advu6*g33[advi6] + 
  advu7*g33[advi7] + advu8*g33[advi8] + advv0*g33[advj0] + 
  advv1*g33[advj1] + advv2*g33[advj2] + advv3*g33[advj3] + 
  advv4*g33[advj4] + advv5*g33[advj5] + advv6*g33[advj6] + 
  advv7*g33[advj7] + advv8*g33[advj8] + advw0*g33[advk0] + 
  advw1*g33[advk1] + advw2*g33[advk2] + advw3*g33[advk3] + 
  advw4*g33[advk4] + advw5*g33[advk5] + advw6*g33[advk6] + 
  advw7*g33[advk7] + advw8*g33[advk8]
;

lieA11
=
advu0*A11[advi0] + advu1*A11[advi1] + advu2*A11[advi2] + advu3*A11[advi3] + 
  advu4*A11[advi4] + advu5*A11[advi5] + advu6*A11[advi6] + 
  advu7*A11[advi7] + advu8*A11[advi8] + advv0*A11[advj0] + 
  advv1*A11[advj1] + advv2*A11[advj2] + advv3*A11[advj3] + 
  advv4*A11[advj4] + advv5*A11[advj5] + advv6*A11[advj6] + 
  advv7*A11[advj7] + advv8*A11[advj8] + advw0*A11[advk0] + 
  advw1*A11[advk1] + advw2*A11[advk2] + advw3*A11[advk3] + 
  advw4*A11[advk4] + advw5*A11[advk5] + advw6*A11[advk6] + 
  advw7*A11[advk7] + advw8*A11[advk8]
;

lieA12
=
advu0*A12[advi0] + advu1*A12[advi1] + advu2*A12[advi2] + advu3*A12[advi3] + 
  advu4*A12[advi4] + advu5*A12[advi5] + advu6*A12[advi6] + 
  advu7*A12[advi7] + advu8*A12[advi8] + advv0*A12[advj0] + 
  advv1*A12[advj1] + advv2*A12[advj2] + advv3*A12[advj3] + 
  advv4*A12[advj4] + advv5*A12[advj5] + advv6*A12[advj6] + 
  advv7*A12[advj7] + advv8*A12[advj8] + advw0*A12[advk0] + 
  advw1*A12[advk1] + advw2*A12[advk2] + advw3*A12[advk3] + 
  advw4*A12[advk4] + advw5*A12[advk5] + advw6*A12[advk6] + 
  advw7*A12[advk7] + advw8*A12[advk8]
;

lieA13
=
advu0*A13[advi0] + advu1*A13[advi1] + advu2*A13[advi2] + advu3*A13[advi3] + 
  advu4*A13[advi4] + advu5*A13[advi5] + advu6*A13[advi6] + 
  advu7*A13[advi7] + advu8*A13[advi8] + advv0*A13[advj0] + 
  advv1*A13[advj1] + advv2*A13[advj2] + advv3*A13[advj3] + 
  advv4*A13[advj4] + advv5*A13[advj5] + advv6*A13[advj6] + 
  advv7*A13[advj7] + advv8*A13[advj8] + advw0*A13[advk0] + 
  advw1*A13[advk1] + advw2*A13[advk2] + advw3*A13[advk3] + 
  advw4*A13[advk4] + advw5*A13[advk5] + advw6*A13[advk6] + 
  advw7*A13[advk7] + advw8*A13[advk8]
;

lieA22
=
advu0*A22[advi0] + advu1*A22[advi1] + advu2*A22[advi2] + advu3*A22[advi3] + 
  advu4*A22[advi4] + advu5*A22[advi5] + advu6*A22[advi6] + 
  advu7*A22[advi7] + advu8*A22[advi8] + advv0*A22[advj0] + 
  advv1*A22[advj1] + advv2*A22[advj2] + advv3*A22[advj3] + 
  advv4*A22[advj4] + advv5*A22[advj5] + advv6*A22[advj6] + 
  advv7*A22[advj7] + advv8*A22[advj8] + advw0*A22[advk0] + 
  advw1*A22[advk1] + advw2*A22[advk2] + advw3*A22[advk3] + 
  advw4*A22[advk4] + advw5*A22[advk5] + advw6*A22[advk6] + 
  advw7*A22[advk7] + advw8*A22[advk8]
;

lieA23
=
advu0*A23[advi0] + advu1*A23[advi1] + advu2*A23[advi2] + advu3*A23[advi3] + 
  advu4*A23[advi4] + advu5*A23[advi5] + advu6*A23[advi6] + 
  advu7*A23[advi7] + advu8*A23[advi8] + advv0*A23[advj0] + 
  advv1*A23[advj1] + advv2*A23[advj2] + advv3*A23[advj3] + 
  advv4*A23[advj4] + advv5*A23[advj5] + advv6*A23[advj6] + 
  advv7*A23[advj7] + advv8*A23[advj8] + advw0*A23[advk0] + 
  advw1*A23[advk1] + advw2*A23[advk2] + advw3*A23[advk3] + 
  advw4*A23[advk4] + advw5*A23[advk5] + advw6*A23[advk6] + 
  advw7*A23[advk7] + advw8*A23[advk8]
;

lieA33
=
advu0*A33[advi0] + advu1*A33[advi1] + advu2*A33[advi2] + advu3*A33[advi3] + 
  advu4*A33[advi4] + advu5*A33[advi5] + advu6*A33[advi6] + 
  advu7*A33[advi7] + advu8*A33[advi8] + advv0*A33[advj0] + 
  advv1*A33[advj1] + advv2*A33[advj2] + advv3*A33[advj3] + 
  advv4*A33[advj4] + advv5*A33[advj5] + advv6*A33[advj6] + 
  advv7*A33[advj7] + advv8*A33[advj8] + advw0*A33[advk0] + 
  advw1*A33[advk1] + advw2*A33[advk2] + advw3*A33[advk3] + 
  advw4*A33[advk4] + advw5*A33[advk5] + advw6*A33[advk6] + 
  advw7*A33[advk7] + advw8*A33[advk8]
;

lieK
=
advu0*K[advi0] + advu1*K[advi1] + advu2*K[advi2] + advu3*K[advi3] + 
  advu4*K[advi4] + advu5*K[advi5] + advu6*K[advi6] + advu7*K[advi7] + 
  advu8*K[advi8] + advv0*K[advj0] + advv1*K[advj1] + advv2*K[advj2] + 
  advv3*K[advj3] + advv4*K[advj4] + advv5*K[advj5] + advv6*K[advj6] + 
  advv7*K[advj7] + advv8*K[advj8] + advw0*K[advk0] + advw1*K[advk1] + 
  advw2*K[advk2] + advw3*K[advk3] + advw4*K[advk4] + advw5*K[advk5] + 
  advw6*K[advk6] + advw7*K[advk7] + advw8*K[advk8]
;

liechi
=
advu0*chi[advi0] + advu1*chi[advi1] + advu2*chi[advi2] + advu3*chi[advi3] + 
  advu4*chi[advi4] + advu5*chi[advi5] + advu6*chi[advi6] + 
  advu7*chi[advi7] + advu8*chi[advi8] + advv0*chi[advj0] + 
  advv1*chi[advj1] + advv2*chi[advj2] + advv3*chi[advj3] + 
  advv4*chi[advj4] + advv5*chi[advj5] + advv6*chi[advj6] + 
  advv7*chi[advj7] + advv8*chi[advj8] + advw0*chi[advk0] + 
  advw1*chi[advk1] + advw2*chi[advk2] + advw3*chi[advk3] + 
  advw4*chi[advk4] + advw5*chi[advk5] + advw6*chi[advk6] + 
  advw7*chi[advk7] + advw8*chi[advk8]
;

liealpha
=
advu0*alpha[advi0] + advu1*alpha[advi1] + advu2*alpha[advi2] + 
  advu3*alpha[advi3] + advu4*alpha[advi4] + advu5*alpha[advi5] + 
  advu6*alpha[advi6] + advu7*alpha[advi7] + advu8*alpha[advi8] + 
  advv0*alpha[advj0] + advv1*alpha[advj1] + advv2*alpha[advj2] + 
  advv3*alpha[advj3] + advv4*alpha[advj4] + advv5*alpha[advj5] + 
  advv6*alpha[advj6] + advv7*alpha[advj7] + advv8*alpha[advj8] + 
  advw0*alpha[advk0] + advw1*alpha[advk1] + advw2*alpha[advk2] + 
  advw3*alpha[advk3] + advw4*alpha[advk4] + advw5*alpha[advk5] + 
  advw6*alpha[advk6] + advw7*alpha[advk7] + advw8*alpha[advk8]
;

advG1
=
advu0*G1[advi0] + advu1*G1[advi1] + advu2*G1[advi2] + advu3*G1[advi3] + 
  advu4*G1[advi4] + advu5*G1[advi5] + advu6*G1[advi6] + advu7*G1[advi7] + 
  advu8*G1[advi8] + advv0*G1[advj0] + advv1*G1[advj1] + advv2*G1[advj2] + 
  advv3*G1[advj3] + advv4*G1[advj4] + advv5*G1[advj5] + advv6*G1[advj6] + 
  advv7*G1[advj7] + advv8*G1[advj8] + advw0*G1[advk0] + advw1*G1[advk1] + 
  advw2*G1[advk2] + advw3*G1[advk3] + advw4*G1[advk4] + advw5*G1[advk5] + 
  advw6*G1[advk6] + advw7*G1[advk7] + advw8*G1[advk8]
;

advG2
=
advu0*G2[advi0] + advu1*G2[advi1] + advu2*G2[advi2] + advu3*G2[advi3] + 
  advu4*G2[advi4] + advu5*G2[advi5] + advu6*G2[advi6] + advu7*G2[advi7] + 
  advu8*G2[advi8] + advv0*G2[advj0] + advv1*G2[advj1] + advv2*G2[advj2] + 
  advv3*G2[advj3] + advv4*G2[advj4] + advv5*G2[advj5] + advv6*G2[advj6] + 
  advv7*G2[advj7] + advv8*G2[advj8] + advw0*G2[advk0] + advw1*G2[advk1] + 
  advw2*G2[advk2] + advw3*G2[advk3] + advw4*G2[advk4] + advw5*G2[advk5] + 
  advw6*G2[advk6] + advw7*G2[advk7] + advw8*G2[advk8]
;

advG3
=
advu0*G3[advi0] + advu1*G3[advi1] + advu2*G3[advi2] + advu3*G3[advi3] + 
  advu4*G3[advi4] + advu5*G3[advi5] + advu6*G3[advi6] + advu7*G3[advi7] + 
  advu8*G3[advi8] + advv0*G3[advj0] + advv1*G3[advj1] + advv2*G3[advj2] + 
  advv3*G3[advj3] + advv4*G3[advj4] + advv5*G3[advj5] + advv6*G3[advj6] + 
  advv7*G3[advj7] + advv8*G3[advj8] + advw0*G3[advk0] + advw1*G3[advk1] + 
  advw2*G3[advk2] + advw3*G3[advk3] + advw4*G3[advk4] + advw5*G3[advk5] + 
  advw6*G3[advk6] + advw7*G3[advk7] + advw8*G3[advk8]
;

advbeta1
=
advu0*beta1[advi0] + advu1*beta1[advi1] + advu2*beta1[advi2] + 
  advu3*beta1[advi3] + advu4*beta1[advi4] + advu5*beta1[advi5] + 
  advu6*beta1[advi6] + advu7*beta1[advi7] + advu8*beta1[advi8] + 
  advv0*beta1[advj0] + advv1*beta1[advj1] + advv2*beta1[advj2] + 
  advv3*beta1[advj3] + advv4*beta1[advj4] + advv5*beta1[advj5] + 
  advv6*beta1[advj6] + advv7*beta1[advj7] + advv8*beta1[advj8] + 
  advw0*beta1[advk0] + advw1*beta1[advk1] + advw2*beta1[advk2] + 
  advw3*beta1[advk3] + advw4*beta1[advk4] + advw5*beta1[advk5] + 
  advw6*beta1[advk6] + advw7*beta1[advk7] + advw8*beta1[advk8]
;

advbeta2
=
advu0*beta2[advi0] + advu1*beta2[advi1] + advu2*beta2[advi2] + 
  advu3*beta2[advi3] + advu4*beta2[advi4] + advu5*beta2[advi5] + 
  advu6*beta2[advi6] + advu7*beta2[advi7] + advu8*beta2[advi8] + 
  advv0*beta2[advj0] + advv1*beta2[advj1] + advv2*beta2[advj2] + 
  advv3*beta2[advj3] + advv4*beta2[advj4] + advv5*beta2[advj5] + 
  advv6*beta2[advj6] + advv7*beta2[advj7] + advv8*beta2[advj8] + 
  advw0*beta2[advk0] + advw1*beta2[advk1] + advw2*beta2[advk2] + 
  advw3*beta2[advk3] + advw4*beta2[advk4] + advw5*beta2[advk5] + 
  advw6*beta2[advk6] + advw7*beta2[advk7] + advw8*beta2[advk8]
;

advbeta3
=
advu0*beta3[advi0] + advu1*beta3[advi1] + advu2*beta3[advi2] + 
  advu3*beta3[advi3] + advu4*beta3[advi4] + advu5*beta3[advi5] + 
  advu6*beta3[advi6] + advu7*beta3[advi7] + advu8*beta3[advi8] + 
  advv0*beta3[advj0] + advv1*beta3[advj1] + advv2*beta3[advj2] + 
  advv3*beta3[advj3] + advv4*beta3[advj4] + advv5*beta3[advj5] + 
  advv6*beta3[advj6] + advv7*beta3[advj7] + advv8*beta3[advj8] + 
  advw0*beta3[advk0] + advw1*beta3[advk1] + advw2*beta3[advk2] + 
  advw3*beta3[advk3] + advw4*beta3[advk4] + advw5*beta3[advk5] + 
  advw6*beta3[advk6] + advw7*beta3[advk7] + advw8*beta3[advk8]
;

advB1
=
advu0*B1[advi0] + advu1*B1[advi1] + advu2*B1[advi2] + advu3*B1[advi3] + 
  advu4*B1[advi4] + advu5*B1[advi5] + advu6*B1[advi6] + advu7*B1[advi7] + 
  advu8*B1[advi8] + advv0*B1[advj0] + advv1*B1[advj1] + advv2*B1[advj2] + 
  advv3*B1[advj3] + advv4*B1[advj4] + advv5*B1[advj5] + advv6*B1[advj6] + 
  advv7*B1[advj7] + advv8*B1[advj8] + advw0*B1[advk0] + advw1*B1[advk1] + 
  advw2*B1[advk2] + advw3*B1[advk3] + advw4*B1[advk4] + advw5*B1[advk5] + 
  advw6*B1[advk6] + advw7*B1[advk7] + advw8*B1[advk8]
;

advB2
=
advu0*B2[advi0] + advu1*B2[advi1] + advu2*B2[advi2] + advu3*B2[advi3] + 
  advu4*B2[advi4] + advu5*B2[advi5] + advu6*B2[advi6] + advu7*B2[advi7] + 
  advu8*B2[advi8] + advv0*B2[advj0] + advv1*B2[advj1] + advv2*B2[advj2] + 
  advv3*B2[advj3] + advv4*B2[advj4] + advv5*B2[advj5] + advv6*B2[advj6] + 
  advv7*B2[advj7] + advv8*B2[advj8] + advw0*B2[advk0] + advw1*B2[advk1] + 
  advw2*B2[advk2] + advw3*B2[advk3] + advw4*B2[advk4] + advw5*B2[advk5] + 
  advw6*B2[advk6] + advw7*B2[advk7] + advw8*B2[advk8]
;

advB3
=
advu0*B3[advi0] + advu1*B3[advi1] + advu2*B3[advi2] + advu3*B3[advi3] + 
  advu4*B3[advi4] + advu5*B3[advi5] + advu6*B3[advi6] + advu7*B3[advi7] + 
  advu8*B3[advi8] + advv0*B3[advj0] + advv1*B3[advj1] + advv2*B3[advj2] + 
  advv3*B3[advj3] + advv4*B3[advj4] + advv5*B3[advj5] + advv6*B3[advj6] + 
  advv7*B3[advj7] + advv8*B3[advj8] + advw0*B3[advk0] + advw1*B3[advk1] + 
  advw2*B3[advk2] + advw3*B3[advk3] + advw4*B3[advk4] + advw5*B3[advk5] + 
  advw6*B3[advk6] + advw7*B3[advk7] + advw8*B3[advk8]
;


#endif 


} else if (order_advection == 10 || boundaryNaway(5)) { 


#ifdef REDUCEORDERTO8                           
  errorexit("CompilerFlag reduced order to 8");
 #else

set_advection10(w1, w2, w3); 

lieg11
=
advu0*g11[advi0] + advu1*g11[advi1] + advu10*g11[advi10] + 
  advu2*g11[advi2] + advu3*g11[advi3] + advu4*g11[advi4] + 
  advu5*g11[advi5] + advu6*g11[advi6] + advu7*g11[advi7] + 
  advu8*g11[advi8] + advu9*g11[advi9] + advv0*g11[advj0] + 
  advv1*g11[advj1] + advv10*g11[advj10] + advv2*g11[advj2] + 
  advv3*g11[advj3] + advv4*g11[advj4] + advv5*g11[advj5] + 
  advv6*g11[advj6] + advv7*g11[advj7] + advv8*g11[advj8] + 
  advv9*g11[advj9] + advw0*g11[advk0] + advw1*g11[advk1] + 
  advw10*g11[advk10] + advw2*g11[advk2] + advw3*g11[advk3] + 
  advw4*g11[advk4] + advw5*g11[advk5] + advw6*g11[advk6] + 
  advw7*g11[advk7] + advw8*g11[advk8] + advw9*g11[advk9]
;

lieg12
=
advu0*g12[advi0] + advu1*g12[advi1] + advu10*g12[advi10] + 
  advu2*g12[advi2] + advu3*g12[advi3] + advu4*g12[advi4] + 
  advu5*g12[advi5] + advu6*g12[advi6] + advu7*g12[advi7] + 
  advu8*g12[advi8] + advu9*g12[advi9] + advv0*g12[advj0] + 
  advv1*g12[advj1] + advv10*g12[advj10] + advv2*g12[advj2] + 
  advv3*g12[advj3] + advv4*g12[advj4] + advv5*g12[advj5] + 
  advv6*g12[advj6] + advv7*g12[advj7] + advv8*g12[advj8] + 
  advv9*g12[advj9] + advw0*g12[advk0] + advw1*g12[advk1] + 
  advw10*g12[advk10] + advw2*g12[advk2] + advw3*g12[advk3] + 
  advw4*g12[advk4] + advw5*g12[advk5] + advw6*g12[advk6] + 
  advw7*g12[advk7] + advw8*g12[advk8] + advw9*g12[advk9]
;

lieg13
=
advu0*g13[advi0] + advu1*g13[advi1] + advu10*g13[advi10] + 
  advu2*g13[advi2] + advu3*g13[advi3] + advu4*g13[advi4] + 
  advu5*g13[advi5] + advu6*g13[advi6] + advu7*g13[advi7] + 
  advu8*g13[advi8] + advu9*g13[advi9] + advv0*g13[advj0] + 
  advv1*g13[advj1] + advv10*g13[advj10] + advv2*g13[advj2] + 
  advv3*g13[advj3] + advv4*g13[advj4] + advv5*g13[advj5] + 
  advv6*g13[advj6] + advv7*g13[advj7] + advv8*g13[advj8] + 
  advv9*g13[advj9] + advw0*g13[advk0] + advw1*g13[advk1] + 
  advw10*g13[advk10] + advw2*g13[advk2] + advw3*g13[advk3] + 
  advw4*g13[advk4] + advw5*g13[advk5] + advw6*g13[advk6] + 
  advw7*g13[advk7] + advw8*g13[advk8] + advw9*g13[advk9]
;

lieg22
=
advu0*g22[advi0] + advu1*g22[advi1] + advu10*g22[advi10] + 
  advu2*g22[advi2] + advu3*g22[advi3] + advu4*g22[advi4] + 
  advu5*g22[advi5] + advu6*g22[advi6] + advu7*g22[advi7] + 
  advu8*g22[advi8] + advu9*g22[advi9] + advv0*g22[advj0] + 
  advv1*g22[advj1] + advv10*g22[advj10] + advv2*g22[advj2] + 
  advv3*g22[advj3] + advv4*g22[advj4] + advv5*g22[advj5] + 
  advv6*g22[advj6] + advv7*g22[advj7] + advv8*g22[advj8] + 
  advv9*g22[advj9] + advw0*g22[advk0] + advw1*g22[advk1] + 
  advw10*g22[advk10] + advw2*g22[advk2] + advw3*g22[advk3] + 
  advw4*g22[advk4] + advw5*g22[advk5] + advw6*g22[advk6] + 
  advw7*g22[advk7] + advw8*g22[advk8] + advw9*g22[advk9]
;

lieg23
=
advu0*g23[advi0] + advu1*g23[advi1] + advu10*g23[advi10] + 
  advu2*g23[advi2] + advu3*g23[advi3] + advu4*g23[advi4] + 
  advu5*g23[advi5] + advu6*g23[advi6] + advu7*g23[advi7] + 
  advu8*g23[advi8] + advu9*g23[advi9] + advv0*g23[advj0] + 
  advv1*g23[advj1] + advv10*g23[advj10] + advv2*g23[advj2] + 
  advv3*g23[advj3] + advv4*g23[advj4] + advv5*g23[advj5] + 
  advv6*g23[advj6] + advv7*g23[advj7] + advv8*g23[advj8] + 
  advv9*g23[advj9] + advw0*g23[advk0] + advw1*g23[advk1] + 
  advw10*g23[advk10] + advw2*g23[advk2] + advw3*g23[advk3] + 
  advw4*g23[advk4] + advw5*g23[advk5] + advw6*g23[advk6] + 
  advw7*g23[advk7] + advw8*g23[advk8] + advw9*g23[advk9]
;

lieg33
=
advu0*g33[advi0] + advu1*g33[advi1] + advu10*g33[advi10] + 
  advu2*g33[advi2] + advu3*g33[advi3] + advu4*g33[advi4] + 
  advu5*g33[advi5] + advu6*g33[advi6] + advu7*g33[advi7] + 
  advu8*g33[advi8] + advu9*g33[advi9] + advv0*g33[advj0] + 
  advv1*g33[advj1] + advv10*g33[advj10] + advv2*g33[advj2] + 
  advv3*g33[advj3] + advv4*g33[advj4] + advv5*g33[advj5] + 
  advv6*g33[advj6] + advv7*g33[advj7] + advv8*g33[advj8] + 
  advv9*g33[advj9] + advw0*g33[advk0] + advw1*g33[advk1] + 
  advw10*g33[advk10] + advw2*g33[advk2] + advw3*g33[advk3] + 
  advw4*g33[advk4] + advw5*g33[advk5] + advw6*g33[advk6] + 
  advw7*g33[advk7] + advw8*g33[advk8] + advw9*g33[advk9]
;

lieA11
=
advu0*A11[advi0] + advu1*A11[advi1] + advu10*A11[advi10] + 
  advu2*A11[advi2] + advu3*A11[advi3] + advu4*A11[advi4] + 
  advu5*A11[advi5] + advu6*A11[advi6] + advu7*A11[advi7] + 
  advu8*A11[advi8] + advu9*A11[advi9] + advv0*A11[advj0] + 
  advv1*A11[advj1] + advv10*A11[advj10] + advv2*A11[advj2] + 
  advv3*A11[advj3] + advv4*A11[advj4] + advv5*A11[advj5] + 
  advv6*A11[advj6] + advv7*A11[advj7] + advv8*A11[advj8] + 
  advv9*A11[advj9] + advw0*A11[advk0] + advw1*A11[advk1] + 
  advw10*A11[advk10] + advw2*A11[advk2] + advw3*A11[advk3] + 
  advw4*A11[advk4] + advw5*A11[advk5] + advw6*A11[advk6] + 
  advw7*A11[advk7] + advw8*A11[advk8] + advw9*A11[advk9]
;

lieA12
=
advu0*A12[advi0] + advu1*A12[advi1] + advu10*A12[advi10] + 
  advu2*A12[advi2] + advu3*A12[advi3] + advu4*A12[advi4] + 
  advu5*A12[advi5] + advu6*A12[advi6] + advu7*A12[advi7] + 
  advu8*A12[advi8] + advu9*A12[advi9] + advv0*A12[advj0] + 
  advv1*A12[advj1] + advv10*A12[advj10] + advv2*A12[advj2] + 
  advv3*A12[advj3] + advv4*A12[advj4] + advv5*A12[advj5] + 
  advv6*A12[advj6] + advv7*A12[advj7] + advv8*A12[advj8] + 
  advv9*A12[advj9] + advw0*A12[advk0] + advw1*A12[advk1] + 
  advw10*A12[advk10] + advw2*A12[advk2] + advw3*A12[advk3] + 
  advw4*A12[advk4] + advw5*A12[advk5] + advw6*A12[advk6] + 
  advw7*A12[advk7] + advw8*A12[advk8] + advw9*A12[advk9]
;

lieA13
=
advu0*A13[advi0] + advu1*A13[advi1] + advu10*A13[advi10] + 
  advu2*A13[advi2] + advu3*A13[advi3] + advu4*A13[advi4] + 
  advu5*A13[advi5] + advu6*A13[advi6] + advu7*A13[advi7] + 
  advu8*A13[advi8] + advu9*A13[advi9] + advv0*A13[advj0] + 
  advv1*A13[advj1] + advv10*A13[advj10] + advv2*A13[advj2] + 
  advv3*A13[advj3] + advv4*A13[advj4] + advv5*A13[advj5] + 
  advv6*A13[advj6] + advv7*A13[advj7] + advv8*A13[advj8] + 
  advv9*A13[advj9] + advw0*A13[advk0] + advw1*A13[advk1] + 
  advw10*A13[advk10] + advw2*A13[advk2] + advw3*A13[advk3] + 
  advw4*A13[advk4] + advw5*A13[advk5] + advw6*A13[advk6] + 
  advw7*A13[advk7] + advw8*A13[advk8] + advw9*A13[advk9]
;

lieA22
=
advu0*A22[advi0] + advu1*A22[advi1] + advu10*A22[advi10] + 
  advu2*A22[advi2] + advu3*A22[advi3] + advu4*A22[advi4] + 
  advu5*A22[advi5] + advu6*A22[advi6] + advu7*A22[advi7] + 
  advu8*A22[advi8] + advu9*A22[advi9] + advv0*A22[advj0] + 
  advv1*A22[advj1] + advv10*A22[advj10] + advv2*A22[advj2] + 
  advv3*A22[advj3] + advv4*A22[advj4] + advv5*A22[advj5] + 
  advv6*A22[advj6] + advv7*A22[advj7] + advv8*A22[advj8] + 
  advv9*A22[advj9] + advw0*A22[advk0] + advw1*A22[advk1] + 
  advw10*A22[advk10] + advw2*A22[advk2] + advw3*A22[advk3] + 
  advw4*A22[advk4] + advw5*A22[advk5] + advw6*A22[advk6] + 
  advw7*A22[advk7] + advw8*A22[advk8] + advw9*A22[advk9]
;

lieA23
=
advu0*A23[advi0] + advu1*A23[advi1] + advu10*A23[advi10] + 
  advu2*A23[advi2] + advu3*A23[advi3] + advu4*A23[advi4] + 
  advu5*A23[advi5] + advu6*A23[advi6] + advu7*A23[advi7] + 
  advu8*A23[advi8] + advu9*A23[advi9] + advv0*A23[advj0] + 
  advv1*A23[advj1] + advv10*A23[advj10] + advv2*A23[advj2] + 
  advv3*A23[advj3] + advv4*A23[advj4] + advv5*A23[advj5] + 
  advv6*A23[advj6] + advv7*A23[advj7] + advv8*A23[advj8] + 
  advv9*A23[advj9] + advw0*A23[advk0] + advw1*A23[advk1] + 
  advw10*A23[advk10] + advw2*A23[advk2] + advw3*A23[advk3] + 
  advw4*A23[advk4] + advw5*A23[advk5] + advw6*A23[advk6] + 
  advw7*A23[advk7] + advw8*A23[advk8] + advw9*A23[advk9]
;

lieA33
=
advu0*A33[advi0] + advu1*A33[advi1] + advu10*A33[advi10] + 
  advu2*A33[advi2] + advu3*A33[advi3] + advu4*A33[advi4] + 
  advu5*A33[advi5] + advu6*A33[advi6] + advu7*A33[advi7] + 
  advu8*A33[advi8] + advu9*A33[advi9] + advv0*A33[advj0] + 
  advv1*A33[advj1] + advv10*A33[advj10] + advv2*A33[advj2] + 
  advv3*A33[advj3] + advv4*A33[advj4] + advv5*A33[advj5] + 
  advv6*A33[advj6] + advv7*A33[advj7] + advv8*A33[advj8] + 
  advv9*A33[advj9] + advw0*A33[advk0] + advw1*A33[advk1] + 
  advw10*A33[advk10] + advw2*A33[advk2] + advw3*A33[advk3] + 
  advw4*A33[advk4] + advw5*A33[advk5] + advw6*A33[advk6] + 
  advw7*A33[advk7] + advw8*A33[advk8] + advw9*A33[advk9]
;

lieK
=
advu0*K[advi0] + advu1*K[advi1] + advu10*K[advi10] + advu2*K[advi2] + 
  advu3*K[advi3] + advu4*K[advi4] + advu5*K[advi5] + advu6*K[advi6] + 
  advu7*K[advi7] + advu8*K[advi8] + advu9*K[advi9] + advv0*K[advj0] + 
  advv1*K[advj1] + advv10*K[advj10] + advv2*K[advj2] + advv3*K[advj3] + 
  advv4*K[advj4] + advv5*K[advj5] + advv6*K[advj6] + advv7*K[advj7] + 
  advv8*K[advj8] + advv9*K[advj9] + advw0*K[advk0] + advw1*K[advk1] + 
  advw10*K[advk10] + advw2*K[advk2] + advw3*K[advk3] + advw4*K[advk4] + 
  advw5*K[advk5] + advw6*K[advk6] + advw7*K[advk7] + advw8*K[advk8] + 
  advw9*K[advk9]
;

liechi
=
advu0*chi[advi0] + advu1*chi[advi1] + advu10*chi[advi10] + 
  advu2*chi[advi2] + advu3*chi[advi3] + advu4*chi[advi4] + 
  advu5*chi[advi5] + advu6*chi[advi6] + advu7*chi[advi7] + 
  advu8*chi[advi8] + advu9*chi[advi9] + advv0*chi[advj0] + 
  advv1*chi[advj1] + advv10*chi[advj10] + advv2*chi[advj2] + 
  advv3*chi[advj3] + advv4*chi[advj4] + advv5*chi[advj5] + 
  advv6*chi[advj6] + advv7*chi[advj7] + advv8*chi[advj8] + 
  advv9*chi[advj9] + advw0*chi[advk0] + advw1*chi[advk1] + 
  advw10*chi[advk10] + advw2*chi[advk2] + advw3*chi[advk3] + 
  advw4*chi[advk4] + advw5*chi[advk5] + advw6*chi[advk6] + 
  advw7*chi[advk7] + advw8*chi[advk8] + advw9*chi[advk9]
;

liealpha
=
advu0*alpha[advi0] + advu1*alpha[advi1] + advu10*alpha[advi10] + 
  advu2*alpha[advi2] + advu3*alpha[advi3] + advu4*alpha[advi4] + 
  advu5*alpha[advi5] + advu6*alpha[advi6] + advu7*alpha[advi7] + 
  advu8*alpha[advi8] + advu9*alpha[advi9] + advv0*alpha[advj0] + 
  advv1*alpha[advj1] + advv10*alpha[advj10] + advv2*alpha[advj2] + 
  advv3*alpha[advj3] + advv4*alpha[advj4] + advv5*alpha[advj5] + 
  advv6*alpha[advj6] + advv7*alpha[advj7] + advv8*alpha[advj8] + 
  advv9*alpha[advj9] + advw0*alpha[advk0] + advw1*alpha[advk1] + 
  advw10*alpha[advk10] + advw2*alpha[advk2] + advw3*alpha[advk3] + 
  advw4*alpha[advk4] + advw5*alpha[advk5] + advw6*alpha[advk6] + 
  advw7*alpha[advk7] + advw8*alpha[advk8] + advw9*alpha[advk9]
;

advG1
=
advu0*G1[advi0] + advu1*G1[advi1] + advu10*G1[advi10] + advu2*G1[advi2] + 
  advu3*G1[advi3] + advu4*G1[advi4] + advu5*G1[advi5] + advu6*G1[advi6] + 
  advu7*G1[advi7] + advu8*G1[advi8] + advu9*G1[advi9] + advv0*G1[advj0] + 
  advv1*G1[advj1] + advv10*G1[advj10] + advv2*G1[advj2] + advv3*G1[advj3] + 
  advv4*G1[advj4] + advv5*G1[advj5] + advv6*G1[advj6] + advv7*G1[advj7] + 
  advv8*G1[advj8] + advv9*G1[advj9] + advw0*G1[advk0] + advw1*G1[advk1] + 
  advw10*G1[advk10] + advw2*G1[advk2] + advw3*G1[advk3] + advw4*G1[advk4] + 
  advw5*G1[advk5] + advw6*G1[advk6] + advw7*G1[advk7] + advw8*G1[advk8] + 
  advw9*G1[advk9]
;

advG2
=
advu0*G2[advi0] + advu1*G2[advi1] + advu10*G2[advi10] + advu2*G2[advi2] + 
  advu3*G2[advi3] + advu4*G2[advi4] + advu5*G2[advi5] + advu6*G2[advi6] + 
  advu7*G2[advi7] + advu8*G2[advi8] + advu9*G2[advi9] + advv0*G2[advj0] + 
  advv1*G2[advj1] + advv10*G2[advj10] + advv2*G2[advj2] + advv3*G2[advj3] + 
  advv4*G2[advj4] + advv5*G2[advj5] + advv6*G2[advj6] + advv7*G2[advj7] + 
  advv8*G2[advj8] + advv9*G2[advj9] + advw0*G2[advk0] + advw1*G2[advk1] + 
  advw10*G2[advk10] + advw2*G2[advk2] + advw3*G2[advk3] + advw4*G2[advk4] + 
  advw5*G2[advk5] + advw6*G2[advk6] + advw7*G2[advk7] + advw8*G2[advk8] + 
  advw9*G2[advk9]
;

advG3
=
advu0*G3[advi0] + advu1*G3[advi1] + advu10*G3[advi10] + advu2*G3[advi2] + 
  advu3*G3[advi3] + advu4*G3[advi4] + advu5*G3[advi5] + advu6*G3[advi6] + 
  advu7*G3[advi7] + advu8*G3[advi8] + advu9*G3[advi9] + advv0*G3[advj0] + 
  advv1*G3[advj1] + advv10*G3[advj10] + advv2*G3[advj2] + advv3*G3[advj3] + 
  advv4*G3[advj4] + advv5*G3[advj5] + advv6*G3[advj6] + advv7*G3[advj7] + 
  advv8*G3[advj8] + advv9*G3[advj9] + advw0*G3[advk0] + advw1*G3[advk1] + 
  advw10*G3[advk10] + advw2*G3[advk2] + advw3*G3[advk3] + advw4*G3[advk4] + 
  advw5*G3[advk5] + advw6*G3[advk6] + advw7*G3[advk7] + advw8*G3[advk8] + 
  advw9*G3[advk9]
;

advbeta1
=
advu0*beta1[advi0] + advu1*beta1[advi1] + advu10*beta1[advi10] + 
  advu2*beta1[advi2] + advu3*beta1[advi3] + advu4*beta1[advi4] + 
  advu5*beta1[advi5] + advu6*beta1[advi6] + advu7*beta1[advi7] + 
  advu8*beta1[advi8] + advu9*beta1[advi9] + advv0*beta1[advj0] + 
  advv1*beta1[advj1] + advv10*beta1[advj10] + advv2*beta1[advj2] + 
  advv3*beta1[advj3] + advv4*beta1[advj4] + advv5*beta1[advj5] + 
  advv6*beta1[advj6] + advv7*beta1[advj7] + advv8*beta1[advj8] + 
  advv9*beta1[advj9] + advw0*beta1[advk0] + advw1*beta1[advk1] + 
  advw10*beta1[advk10] + advw2*beta1[advk2] + advw3*beta1[advk3] + 
  advw4*beta1[advk4] + advw5*beta1[advk5] + advw6*beta1[advk6] + 
  advw7*beta1[advk7] + advw8*beta1[advk8] + advw9*beta1[advk9]
;

advbeta2
=
advu0*beta2[advi0] + advu1*beta2[advi1] + advu10*beta2[advi10] + 
  advu2*beta2[advi2] + advu3*beta2[advi3] + advu4*beta2[advi4] + 
  advu5*beta2[advi5] + advu6*beta2[advi6] + advu7*beta2[advi7] + 
  advu8*beta2[advi8] + advu9*beta2[advi9] + advv0*beta2[advj0] + 
  advv1*beta2[advj1] + advv10*beta2[advj10] + advv2*beta2[advj2] + 
  advv3*beta2[advj3] + advv4*beta2[advj4] + advv5*beta2[advj5] + 
  advv6*beta2[advj6] + advv7*beta2[advj7] + advv8*beta2[advj8] + 
  advv9*beta2[advj9] + advw0*beta2[advk0] + advw1*beta2[advk1] + 
  advw10*beta2[advk10] + advw2*beta2[advk2] + advw3*beta2[advk3] + 
  advw4*beta2[advk4] + advw5*beta2[advk5] + advw6*beta2[advk6] + 
  advw7*beta2[advk7] + advw8*beta2[advk8] + advw9*beta2[advk9]
;

advbeta3
=
advu0*beta3[advi0] + advu1*beta3[advi1] + advu10*beta3[advi10] + 
  advu2*beta3[advi2] + advu3*beta3[advi3] + advu4*beta3[advi4] + 
  advu5*beta3[advi5] + advu6*beta3[advi6] + advu7*beta3[advi7] + 
  advu8*beta3[advi8] + advu9*beta3[advi9] + advv0*beta3[advj0] + 
  advv1*beta3[advj1] + advv10*beta3[advj10] + advv2*beta3[advj2] + 
  advv3*beta3[advj3] + advv4*beta3[advj4] + advv5*beta3[advj5] + 
  advv6*beta3[advj6] + advv7*beta3[advj7] + advv8*beta3[advj8] + 
  advv9*beta3[advj9] + advw0*beta3[advk0] + advw1*beta3[advk1] + 
  advw10*beta3[advk10] + advw2*beta3[advk2] + advw3*beta3[advk3] + 
  advw4*beta3[advk4] + advw5*beta3[advk5] + advw6*beta3[advk6] + 
  advw7*beta3[advk7] + advw8*beta3[advk8] + advw9*beta3[advk9]
;

advB1
=
advu0*B1[advi0] + advu1*B1[advi1] + advu10*B1[advi10] + advu2*B1[advi2] + 
  advu3*B1[advi3] + advu4*B1[advi4] + advu5*B1[advi5] + advu6*B1[advi6] + 
  advu7*B1[advi7] + advu8*B1[advi8] + advu9*B1[advi9] + advv0*B1[advj0] + 
  advv1*B1[advj1] + advv10*B1[advj10] + advv2*B1[advj2] + advv3*B1[advj3] + 
  advv4*B1[advj4] + advv5*B1[advj5] + advv6*B1[advj6] + advv7*B1[advj7] + 
  advv8*B1[advj8] + advv9*B1[advj9] + advw0*B1[advk0] + advw1*B1[advk1] + 
  advw10*B1[advk10] + advw2*B1[advk2] + advw3*B1[advk3] + advw4*B1[advk4] + 
  advw5*B1[advk5] + advw6*B1[advk6] + advw7*B1[advk7] + advw8*B1[advk8] + 
  advw9*B1[advk9]
;

advB2
=
advu0*B2[advi0] + advu1*B2[advi1] + advu10*B2[advi10] + advu2*B2[advi2] + 
  advu3*B2[advi3] + advu4*B2[advi4] + advu5*B2[advi5] + advu6*B2[advi6] + 
  advu7*B2[advi7] + advu8*B2[advi8] + advu9*B2[advi9] + advv0*B2[advj0] + 
  advv1*B2[advj1] + advv10*B2[advj10] + advv2*B2[advj2] + advv3*B2[advj3] + 
  advv4*B2[advj4] + advv5*B2[advj5] + advv6*B2[advj6] + advv7*B2[advj7] + 
  advv8*B2[advj8] + advv9*B2[advj9] + advw0*B2[advk0] + advw1*B2[advk1] + 
  advw10*B2[advk10] + advw2*B2[advk2] + advw3*B2[advk3] + advw4*B2[advk4] + 
  advw5*B2[advk5] + advw6*B2[advk6] + advw7*B2[advk7] + advw8*B2[advk8] + 
  advw9*B2[advk9]
;

advB3
=
advu0*B3[advi0] + advu1*B3[advi1] + advu10*B3[advi10] + advu2*B3[advi2] + 
  advu3*B3[advi3] + advu4*B3[advi4] + advu5*B3[advi5] + advu6*B3[advi6] + 
  advu7*B3[advi7] + advu8*B3[advi8] + advu9*B3[advi9] + advv0*B3[advj0] + 
  advv1*B3[advj1] + advv10*B3[advj10] + advv2*B3[advj2] + advv3*B3[advj3] + 
  advv4*B3[advj4] + advv5*B3[advj5] + advv6*B3[advj6] + advv7*B3[advj7] + 
  advv8*B3[advj8] + advv9*B3[advj9] + advw0*B3[advk0] + advw1*B3[advk1] + 
  advw10*B3[advk10] + advw2*B3[advk2] + advw3*B3[advk3] + advw4*B3[advk4] + 
  advw5*B3[advk5] + advw6*B3[advk6] + advw7*B3[advk7] + advw8*B3[advk8] + 
  advw9*B3[advk9]
;


#endif 


} else errorexit("order not implemented");  


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


if (order_advection == 2 || boundaryNaway(1)) { 

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


} else if (order_advection == 4 || boundaryNaway(2)) { 

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


} else if (order_advection == 6 || boundaryNaway(3)) { 

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


} else if (order_centered == 8 || boundaryNaway(4)) { 

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


} else { 

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


} 

daSST1
=
da1*Jac11 + da2*Jac21 + da3*Jac31
;

daSST2
=
da1*Jac12 + da2*Jac22 + da3*Jac32
;

daSST3
=
da1*Jac13 + da2*Jac23 + da3*Jac33
;

ddaSST11
=
da1*DJac111 + da2*DJac211 + da3*DJac311 + 
  2.*(dda23*Jac21*Jac31 + Jac11*(dda12*Jac21 + dda13*Jac31)) + 
  dda11*pow2(Jac11) + dda22*pow2(Jac21) + dda33*pow2(Jac31)
;

ddaSST12
=
da1*DJac112 + da2*DJac212 + da3*DJac312 + 
  Jac12*(dda11*Jac11 + dda12*Jac21 + dda13*Jac31) + 
  Jac22*(dda12*Jac11 + dda22*Jac21 + dda23*Jac31) + 
  (dda13*Jac11 + dda23*Jac21 + dda33*Jac31)*Jac32
;

ddaSST13
=
da1*DJac113 + da2*DJac213 + da3*DJac313 + 
  Jac13*(dda11*Jac11 + dda12*Jac21 + dda13*Jac31) + 
  Jac23*(dda12*Jac11 + dda22*Jac21 + dda23*Jac31) + 
  (dda13*Jac11 + dda23*Jac21 + dda33*Jac31)*Jac33
;

ddaSST22
=
da1*DJac122 + da2*DJac222 + da3*DJac322 + 
  2.*(dda23*Jac22*Jac32 + Jac12*(dda12*Jac22 + dda13*Jac32)) + 
  dda11*pow2(Jac12) + dda22*pow2(Jac22) + dda33*pow2(Jac32)
;

ddaSST23
=
da1*DJac123 + da2*DJac223 + da3*DJac323 + 
  Jac13*(dda11*Jac12 + dda12*Jac22 + dda13*Jac32) + 
  Jac23*(dda12*Jac12 + dda22*Jac22 + dda23*Jac32) + 
  (dda13*Jac12 + dda23*Jac22 + dda33*Jac32)*Jac33
;

ddaSST33
=
da1*DJac133 + da2*DJac233 + da3*DJac333 + 
  2.*(dda23*Jac23*Jac33 + Jac13*(dda12*Jac23 + dda13*Jac33)) + 
  dda11*pow2(Jac13) + dda22*pow2(Jac23) + dda33*pow2(Jac33)
;

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

ddbSST111
=
db11*DJac111 + db21*DJac211 + db31*DJac311 + 
  2.*(ddb231*Jac21*Jac31 + Jac11*(ddb121*Jac21 + ddb131*Jac31)) + 
  ddb111*pow2(Jac11) + ddb221*pow2(Jac21) + ddb331*pow2(Jac31)
;

ddbSST112
=
db12*DJac111 + db22*DJac211 + db32*DJac311 + 
  2.*(ddb232*Jac21*Jac31 + Jac11*(ddb122*Jac21 + ddb132*Jac31)) + 
  ddb112*pow2(Jac11) + ddb222*pow2(Jac21) + ddb332*pow2(Jac31)
;

ddbSST113
=
db13*DJac111 + db23*DJac211 + db33*DJac311 + 
  2.*(ddb233*Jac21*Jac31 + Jac11*(ddb123*Jac21 + ddb133*Jac31)) + 
  ddb113*pow2(Jac11) + ddb223*pow2(Jac21) + ddb333*pow2(Jac31)
;

ddbSST121
=
db11*DJac112 + db21*DJac212 + db31*DJac312 + 
  Jac12*(ddb111*Jac11 + ddb121*Jac21 + ddb131*Jac31) + 
  Jac22*(ddb121*Jac11 + ddb221*Jac21 + ddb231*Jac31) + 
  (ddb131*Jac11 + ddb231*Jac21 + ddb331*Jac31)*Jac32
;

ddbSST122
=
db12*DJac112 + db22*DJac212 + db32*DJac312 + 
  Jac12*(ddb112*Jac11 + ddb122*Jac21 + ddb132*Jac31) + 
  Jac22*(ddb122*Jac11 + ddb222*Jac21 + ddb232*Jac31) + 
  (ddb132*Jac11 + ddb232*Jac21 + ddb332*Jac31)*Jac32
;

ddbSST123
=
db13*DJac112 + db23*DJac212 + db33*DJac312 + 
  Jac12*(ddb113*Jac11 + ddb123*Jac21 + ddb133*Jac31) + 
  Jac22*(ddb123*Jac11 + ddb223*Jac21 + ddb233*Jac31) + 
  (ddb133*Jac11 + ddb233*Jac21 + ddb333*Jac31)*Jac32
;

ddbSST131
=
db11*DJac113 + db21*DJac213 + db31*DJac313 + 
  Jac13*(ddb111*Jac11 + ddb121*Jac21 + ddb131*Jac31) + 
  Jac23*(ddb121*Jac11 + ddb221*Jac21 + ddb231*Jac31) + 
  (ddb131*Jac11 + ddb231*Jac21 + ddb331*Jac31)*Jac33
;

ddbSST132
=
db12*DJac113 + db22*DJac213 + db32*DJac313 + 
  Jac13*(ddb112*Jac11 + ddb122*Jac21 + ddb132*Jac31) + 
  Jac23*(ddb122*Jac11 + ddb222*Jac21 + ddb232*Jac31) + 
  (ddb132*Jac11 + ddb232*Jac21 + ddb332*Jac31)*Jac33
;

ddbSST133
=
db13*DJac113 + db23*DJac213 + db33*DJac313 + 
  Jac13*(ddb113*Jac11 + ddb123*Jac21 + ddb133*Jac31) + 
  Jac23*(ddb123*Jac11 + ddb223*Jac21 + ddb233*Jac31) + 
  (ddb133*Jac11 + ddb233*Jac21 + ddb333*Jac31)*Jac33
;

ddbSST221
=
db11*DJac122 + db21*DJac222 + db31*DJac322 + 
  2.*(ddb231*Jac22*Jac32 + Jac12*(ddb121*Jac22 + ddb131*Jac32)) + 
  ddb111*pow2(Jac12) + ddb221*pow2(Jac22) + ddb331*pow2(Jac32)
;

ddbSST222
=
db12*DJac122 + db22*DJac222 + db32*DJac322 + 
  2.*(ddb232*Jac22*Jac32 + Jac12*(ddb122*Jac22 + ddb132*Jac32)) + 
  ddb112*pow2(Jac12) + ddb222*pow2(Jac22) + ddb332*pow2(Jac32)
;

ddbSST223
=
db13*DJac122 + db23*DJac222 + db33*DJac322 + 
  2.*(ddb233*Jac22*Jac32 + Jac12*(ddb123*Jac22 + ddb133*Jac32)) + 
  ddb113*pow2(Jac12) + ddb223*pow2(Jac22) + ddb333*pow2(Jac32)
;

ddbSST231
=
db11*DJac123 + db21*DJac223 + db31*DJac323 + 
  Jac13*(ddb111*Jac12 + ddb121*Jac22 + ddb131*Jac32) + 
  Jac23*(ddb121*Jac12 + ddb221*Jac22 + ddb231*Jac32) + 
  (ddb131*Jac12 + ddb231*Jac22 + ddb331*Jac32)*Jac33
;

ddbSST232
=
db12*DJac123 + db22*DJac223 + db32*DJac323 + 
  Jac13*(ddb112*Jac12 + ddb122*Jac22 + ddb132*Jac32) + 
  Jac23*(ddb122*Jac12 + ddb222*Jac22 + ddb232*Jac32) + 
  (ddb132*Jac12 + ddb232*Jac22 + ddb332*Jac32)*Jac33
;

ddbSST233
=
db13*DJac123 + db23*DJac223 + db33*DJac323 + 
  Jac13*(ddb113*Jac12 + ddb123*Jac22 + ddb133*Jac32) + 
  Jac23*(ddb123*Jac12 + ddb223*Jac22 + ddb233*Jac32) + 
  (ddb133*Jac12 + ddb233*Jac22 + ddb333*Jac32)*Jac33
;

ddbSST331
=
db11*DJac133 + db21*DJac233 + db31*DJac333 + 
  2.*(ddb231*Jac23*Jac33 + Jac13*(ddb121*Jac23 + ddb131*Jac33)) + 
  ddb111*pow2(Jac13) + ddb221*pow2(Jac23) + ddb331*pow2(Jac33)
;

ddbSST332
=
db12*DJac133 + db22*DJac233 + db32*DJac333 + 
  2.*(ddb232*Jac23*Jac33 + Jac13*(ddb122*Jac23 + ddb132*Jac33)) + 
  ddb112*pow2(Jac13) + ddb222*pow2(Jac23) + ddb332*pow2(Jac33)
;

ddbSST333
=
db13*DJac133 + db23*DJac233 + db33*DJac333 + 
  2.*(ddb233*Jac23*Jac33 + Jac13*(ddb123*Jac23 + ddb133*Jac33)) + 
  ddb113*pow2(Jac13) + ddb223*pow2(Jac23) + ddb333*pow2(Jac33)
;

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

deldelgSST1111
=
delg111*DJac111 + delg211*DJac211 + delg311*DJac311 + 
  2.*(deldelg2311*Jac21*Jac31 + 
     Jac11*(deldelg1211*Jac21 + deldelg1311*Jac31)) + 
  deldelg1111*pow2(Jac11) + deldelg2211*pow2(Jac21) + deldelg3311*pow2(Jac31)
;

deldelgSST1112
=
delg112*DJac111 + delg212*DJac211 + delg312*DJac311 + 
  2.*(deldelg2312*Jac21*Jac31 + 
     Jac11*(deldelg1212*Jac21 + deldelg1312*Jac31)) + 
  deldelg1112*pow2(Jac11) + deldelg2212*pow2(Jac21) + deldelg3312*pow2(Jac31)
;

deldelgSST1113
=
delg113*DJac111 + delg213*DJac211 + delg313*DJac311 + 
  2.*(deldelg2313*Jac21*Jac31 + 
     Jac11*(deldelg1213*Jac21 + deldelg1313*Jac31)) + 
  deldelg1113*pow2(Jac11) + deldelg2213*pow2(Jac21) + deldelg3313*pow2(Jac31)
;

deldelgSST1122
=
delg122*DJac111 + delg222*DJac211 + delg322*DJac311 + 
  2.*(deldelg2322*Jac21*Jac31 + 
     Jac11*(deldelg1222*Jac21 + deldelg1322*Jac31)) + 
  deldelg1122*pow2(Jac11) + deldelg2222*pow2(Jac21) + deldelg3322*pow2(Jac31)
;

deldelgSST1123
=
delg123*DJac111 + delg223*DJac211 + delg323*DJac311 + 
  2.*(deldelg2323*Jac21*Jac31 + 
     Jac11*(deldelg1223*Jac21 + deldelg1323*Jac31)) + 
  deldelg1123*pow2(Jac11) + deldelg2223*pow2(Jac21) + deldelg3323*pow2(Jac31)
;

deldelgSST1133
=
delg133*DJac111 + delg233*DJac211 + delg333*DJac311 + 
  2.*(deldelg2333*Jac21*Jac31 + 
     Jac11*(deldelg1233*Jac21 + deldelg1333*Jac31)) + 
  deldelg1133*pow2(Jac11) + deldelg2233*pow2(Jac21) + deldelg3333*pow2(Jac31)
;

deldelgSST1211
=
delg111*DJac112 + delg211*DJac212 + delg311*DJac312 + 
  Jac12*(deldelg1111*Jac11 + deldelg1211*Jac21 + deldelg1311*Jac31) + 
  Jac22*(deldelg1211*Jac11 + deldelg2211*Jac21 + deldelg2311*Jac31) + 
  (deldelg1311*Jac11 + deldelg2311*Jac21 + deldelg3311*Jac31)*Jac32
;

deldelgSST1212
=
delg112*DJac112 + delg212*DJac212 + delg312*DJac312 + 
  Jac12*(deldelg1112*Jac11 + deldelg1212*Jac21 + deldelg1312*Jac31) + 
  Jac22*(deldelg1212*Jac11 + deldelg2212*Jac21 + deldelg2312*Jac31) + 
  (deldelg1312*Jac11 + deldelg2312*Jac21 + deldelg3312*Jac31)*Jac32
;

deldelgSST1213
=
delg113*DJac112 + delg213*DJac212 + delg313*DJac312 + 
  Jac12*(deldelg1113*Jac11 + deldelg1213*Jac21 + deldelg1313*Jac31) + 
  Jac22*(deldelg1213*Jac11 + deldelg2213*Jac21 + deldelg2313*Jac31) + 
  (deldelg1313*Jac11 + deldelg2313*Jac21 + deldelg3313*Jac31)*Jac32
;

deldelgSST1222
=
delg122*DJac112 + delg222*DJac212 + delg322*DJac312 + 
  Jac12*(deldelg1122*Jac11 + deldelg1222*Jac21 + deldelg1322*Jac31) + 
  Jac22*(deldelg1222*Jac11 + deldelg2222*Jac21 + deldelg2322*Jac31) + 
  (deldelg1322*Jac11 + deldelg2322*Jac21 + deldelg3322*Jac31)*Jac32
;

deldelgSST1223
=
delg123*DJac112 + delg223*DJac212 + delg323*DJac312 + 
  Jac12*(deldelg1123*Jac11 + deldelg1223*Jac21 + deldelg1323*Jac31) + 
  Jac22*(deldelg1223*Jac11 + deldelg2223*Jac21 + deldelg2323*Jac31) + 
  (deldelg1323*Jac11 + deldelg2323*Jac21 + deldelg3323*Jac31)*Jac32
;

deldelgSST1233
=
delg133*DJac112 + delg233*DJac212 + delg333*DJac312 + 
  Jac12*(deldelg1133*Jac11 + deldelg1233*Jac21 + deldelg1333*Jac31) + 
  Jac22*(deldelg1233*Jac11 + deldelg2233*Jac21 + deldelg2333*Jac31) + 
  (deldelg1333*Jac11 + deldelg2333*Jac21 + deldelg3333*Jac31)*Jac32
;

deldelgSST1311
=
delg111*DJac113 + delg211*DJac213 + delg311*DJac313 + 
  Jac13*(deldelg1111*Jac11 + deldelg1211*Jac21 + deldelg1311*Jac31) + 
  Jac23*(deldelg1211*Jac11 + deldelg2211*Jac21 + deldelg2311*Jac31) + 
  (deldelg1311*Jac11 + deldelg2311*Jac21 + deldelg3311*Jac31)*Jac33
;

deldelgSST1312
=
delg112*DJac113 + delg212*DJac213 + delg312*DJac313 + 
  Jac13*(deldelg1112*Jac11 + deldelg1212*Jac21 + deldelg1312*Jac31) + 
  Jac23*(deldelg1212*Jac11 + deldelg2212*Jac21 + deldelg2312*Jac31) + 
  (deldelg1312*Jac11 + deldelg2312*Jac21 + deldelg3312*Jac31)*Jac33
;

deldelgSST1313
=
delg113*DJac113 + delg213*DJac213 + delg313*DJac313 + 
  Jac13*(deldelg1113*Jac11 + deldelg1213*Jac21 + deldelg1313*Jac31) + 
  Jac23*(deldelg1213*Jac11 + deldelg2213*Jac21 + deldelg2313*Jac31) + 
  (deldelg1313*Jac11 + deldelg2313*Jac21 + deldelg3313*Jac31)*Jac33
;

deldelgSST1322
=
delg122*DJac113 + delg222*DJac213 + delg322*DJac313 + 
  Jac13*(deldelg1122*Jac11 + deldelg1222*Jac21 + deldelg1322*Jac31) + 
  Jac23*(deldelg1222*Jac11 + deldelg2222*Jac21 + deldelg2322*Jac31) + 
  (deldelg1322*Jac11 + deldelg2322*Jac21 + deldelg3322*Jac31)*Jac33
;

deldelgSST1323
=
delg123*DJac113 + delg223*DJac213 + delg323*DJac313 + 
  Jac13*(deldelg1123*Jac11 + deldelg1223*Jac21 + deldelg1323*Jac31) + 
  Jac23*(deldelg1223*Jac11 + deldelg2223*Jac21 + deldelg2323*Jac31) + 
  (deldelg1323*Jac11 + deldelg2323*Jac21 + deldelg3323*Jac31)*Jac33
;

deldelgSST1333
=
delg133*DJac113 + delg233*DJac213 + delg333*DJac313 + 
  Jac13*(deldelg1133*Jac11 + deldelg1233*Jac21 + deldelg1333*Jac31) + 
  Jac23*(deldelg1233*Jac11 + deldelg2233*Jac21 + deldelg2333*Jac31) + 
  (deldelg1333*Jac11 + deldelg2333*Jac21 + deldelg3333*Jac31)*Jac33
;

deldelgSST2211
=
delg111*DJac122 + delg211*DJac222 + delg311*DJac322 + 
  2.*(deldelg2311*Jac22*Jac32 + 
     Jac12*(deldelg1211*Jac22 + deldelg1311*Jac32)) + 
  deldelg1111*pow2(Jac12) + deldelg2211*pow2(Jac22) + deldelg3311*pow2(Jac32)
;

deldelgSST2212
=
delg112*DJac122 + delg212*DJac222 + delg312*DJac322 + 
  2.*(deldelg2312*Jac22*Jac32 + 
     Jac12*(deldelg1212*Jac22 + deldelg1312*Jac32)) + 
  deldelg1112*pow2(Jac12) + deldelg2212*pow2(Jac22) + deldelg3312*pow2(Jac32)
;

deldelgSST2213
=
delg113*DJac122 + delg213*DJac222 + delg313*DJac322 + 
  2.*(deldelg2313*Jac22*Jac32 + 
     Jac12*(deldelg1213*Jac22 + deldelg1313*Jac32)) + 
  deldelg1113*pow2(Jac12) + deldelg2213*pow2(Jac22) + deldelg3313*pow2(Jac32)
;

deldelgSST2222
=
delg122*DJac122 + delg222*DJac222 + delg322*DJac322 + 
  2.*(deldelg2322*Jac22*Jac32 + 
     Jac12*(deldelg1222*Jac22 + deldelg1322*Jac32)) + 
  deldelg1122*pow2(Jac12) + deldelg2222*pow2(Jac22) + deldelg3322*pow2(Jac32)
;

deldelgSST2223
=
delg123*DJac122 + delg223*DJac222 + delg323*DJac322 + 
  2.*(deldelg2323*Jac22*Jac32 + 
     Jac12*(deldelg1223*Jac22 + deldelg1323*Jac32)) + 
  deldelg1123*pow2(Jac12) + deldelg2223*pow2(Jac22) + deldelg3323*pow2(Jac32)
;

deldelgSST2233
=
delg133*DJac122 + delg233*DJac222 + delg333*DJac322 + 
  2.*(deldelg2333*Jac22*Jac32 + 
     Jac12*(deldelg1233*Jac22 + deldelg1333*Jac32)) + 
  deldelg1133*pow2(Jac12) + deldelg2233*pow2(Jac22) + deldelg3333*pow2(Jac32)
;

deldelgSST2311
=
delg111*DJac123 + delg211*DJac223 + delg311*DJac323 + 
  Jac13*(deldelg1111*Jac12 + deldelg1211*Jac22 + deldelg1311*Jac32) + 
  Jac23*(deldelg1211*Jac12 + deldelg2211*Jac22 + deldelg2311*Jac32) + 
  (deldelg1311*Jac12 + deldelg2311*Jac22 + deldelg3311*Jac32)*Jac33
;

deldelgSST2312
=
delg112*DJac123 + delg212*DJac223 + delg312*DJac323 + 
  Jac13*(deldelg1112*Jac12 + deldelg1212*Jac22 + deldelg1312*Jac32) + 
  Jac23*(deldelg1212*Jac12 + deldelg2212*Jac22 + deldelg2312*Jac32) + 
  (deldelg1312*Jac12 + deldelg2312*Jac22 + deldelg3312*Jac32)*Jac33
;

deldelgSST2313
=
delg113*DJac123 + delg213*DJac223 + delg313*DJac323 + 
  Jac13*(deldelg1113*Jac12 + deldelg1213*Jac22 + deldelg1313*Jac32) + 
  Jac23*(deldelg1213*Jac12 + deldelg2213*Jac22 + deldelg2313*Jac32) + 
  (deldelg1313*Jac12 + deldelg2313*Jac22 + deldelg3313*Jac32)*Jac33
;

deldelgSST2322
=
delg122*DJac123 + delg222*DJac223 + delg322*DJac323 + 
  Jac13*(deldelg1122*Jac12 + deldelg1222*Jac22 + deldelg1322*Jac32) + 
  Jac23*(deldelg1222*Jac12 + deldelg2222*Jac22 + deldelg2322*Jac32) + 
  (deldelg1322*Jac12 + deldelg2322*Jac22 + deldelg3322*Jac32)*Jac33
;

deldelgSST2323
=
delg123*DJac123 + delg223*DJac223 + delg323*DJac323 + 
  Jac13*(deldelg1123*Jac12 + deldelg1223*Jac22 + deldelg1323*Jac32) + 
  Jac23*(deldelg1223*Jac12 + deldelg2223*Jac22 + deldelg2323*Jac32) + 
  (deldelg1323*Jac12 + deldelg2323*Jac22 + deldelg3323*Jac32)*Jac33
;

deldelgSST2333
=
delg133*DJac123 + delg233*DJac223 + delg333*DJac323 + 
  Jac13*(deldelg1133*Jac12 + deldelg1233*Jac22 + deldelg1333*Jac32) + 
  Jac23*(deldelg1233*Jac12 + deldelg2233*Jac22 + deldelg2333*Jac32) + 
  (deldelg1333*Jac12 + deldelg2333*Jac22 + deldelg3333*Jac32)*Jac33
;

deldelgSST3311
=
delg111*DJac133 + delg211*DJac233 + delg311*DJac333 + 
  2.*(deldelg2311*Jac23*Jac33 + 
     Jac13*(deldelg1211*Jac23 + deldelg1311*Jac33)) + 
  deldelg1111*pow2(Jac13) + deldelg2211*pow2(Jac23) + deldelg3311*pow2(Jac33)
;

deldelgSST3312
=
delg112*DJac133 + delg212*DJac233 + delg312*DJac333 + 
  2.*(deldelg2312*Jac23*Jac33 + 
     Jac13*(deldelg1212*Jac23 + deldelg1312*Jac33)) + 
  deldelg1112*pow2(Jac13) + deldelg2212*pow2(Jac23) + deldelg3312*pow2(Jac33)
;

deldelgSST3313
=
delg113*DJac133 + delg213*DJac233 + delg313*DJac333 + 
  2.*(deldelg2313*Jac23*Jac33 + 
     Jac13*(deldelg1213*Jac23 + deldelg1313*Jac33)) + 
  deldelg1113*pow2(Jac13) + deldelg2213*pow2(Jac23) + deldelg3313*pow2(Jac33)
;

deldelgSST3322
=
delg122*DJac133 + delg222*DJac233 + delg322*DJac333 + 
  2.*(deldelg2322*Jac23*Jac33 + 
     Jac13*(deldelg1222*Jac23 + deldelg1322*Jac33)) + 
  deldelg1122*pow2(Jac13) + deldelg2222*pow2(Jac23) + deldelg3322*pow2(Jac33)
;

deldelgSST3323
=
delg123*DJac133 + delg223*DJac233 + delg323*DJac333 + 
  2.*(deldelg2323*Jac23*Jac33 + 
     Jac13*(deldelg1223*Jac23 + deldelg1323*Jac33)) + 
  deldelg1123*pow2(Jac13) + deldelg2223*pow2(Jac23) + deldelg3323*pow2(Jac33)
;

deldelgSST3333
=
delg133*DJac133 + delg233*DJac233 + delg333*DJac333 + 
  2.*(deldelg2333*Jac23*Jac33 + 
     Jac13*(deldelg1233*Jac23 + deldelg1333*Jac33)) + 
  deldelg1133*pow2(Jac13) + deldelg2233*pow2(Jac23) + deldelg3333*pow2(Jac33)
;

delGSST11
=
delG11*Jac11 + delG21*Jac21 + delG31*Jac31
;

delGSST12
=
delG12*Jac11 + delG22*Jac21 + delG32*Jac31
;

delGSST13
=
delG13*Jac11 + delG23*Jac21 + delG33*Jac31
;

delGSST21
=
delG11*Jac12 + delG21*Jac22 + delG31*Jac32
;

delGSST22
=
delG12*Jac12 + delG22*Jac22 + delG32*Jac32
;

delGSST23
=
delG13*Jac12 + delG23*Jac22 + delG33*Jac32
;

delGSST31
=
delG11*Jac13 + delG21*Jac23 + delG31*Jac33
;

delGSST32
=
delG12*Jac13 + delG22*Jac23 + delG32*Jac33
;

delGSST33
=
delG13*Jac13 + delG23*Jac23 + delG33*Jac33
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

ddchiSST11
=
dchi1*DJac111 + dchi2*DJac211 + dchi3*DJac311 + 
  2.*(ddchi23*Jac21*Jac31 + Jac11*(ddchi12*Jac21 + ddchi13*Jac31)) + 
  ddchi11*pow2(Jac11) + ddchi22*pow2(Jac21) + ddchi33*pow2(Jac31)
;

ddchiSST12
=
dchi1*DJac112 + dchi2*DJac212 + dchi3*DJac312 + 
  Jac12*(ddchi11*Jac11 + ddchi12*Jac21 + ddchi13*Jac31) + 
  Jac22*(ddchi12*Jac11 + ddchi22*Jac21 + ddchi23*Jac31) + 
  (ddchi13*Jac11 + ddchi23*Jac21 + ddchi33*Jac31)*Jac32
;

ddchiSST13
=
dchi1*DJac113 + dchi2*DJac213 + dchi3*DJac313 + 
  Jac13*(ddchi11*Jac11 + ddchi12*Jac21 + ddchi13*Jac31) + 
  Jac23*(ddchi12*Jac11 + ddchi22*Jac21 + ddchi23*Jac31) + 
  (ddchi13*Jac11 + ddchi23*Jac21 + ddchi33*Jac31)*Jac33
;

ddchiSST22
=
dchi1*DJac122 + dchi2*DJac222 + dchi3*DJac322 + 
  2.*(ddchi23*Jac22*Jac32 + Jac12*(ddchi12*Jac22 + ddchi13*Jac32)) + 
  ddchi11*pow2(Jac12) + ddchi22*pow2(Jac22) + ddchi33*pow2(Jac32)
;

ddchiSST23
=
dchi1*DJac123 + dchi2*DJac223 + dchi3*DJac323 + 
  Jac13*(ddchi11*Jac12 + ddchi12*Jac22 + ddchi13*Jac32) + 
  Jac23*(ddchi12*Jac12 + ddchi22*Jac22 + ddchi23*Jac32) + 
  (ddchi13*Jac12 + ddchi23*Jac22 + ddchi33*Jac32)*Jac33
;

ddchiSST33
=
dchi1*DJac133 + dchi2*DJac233 + dchi3*DJac333 + 
  2.*(ddchi23*Jac23*Jac33 + Jac13*(ddchi12*Jac23 + ddchi13*Jac33)) + 
  ddchi11*pow2(Jac13) + ddchi22*pow2(Jac23) + ddchi33*pow2(Jac33)
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

da1
=
daSST1
;

da2
=
daSST2
;

da3
=
daSST3
;

dda11
=
ddaSST11
;

dda12
=
ddaSST12
;

dda13
=
ddaSST13
;

dda22
=
ddaSST22
;

dda23
=
ddaSST23
;

dda33
=
ddaSST33
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

ddb111
=
ddbSST111
;

ddb112
=
ddbSST112
;

ddb113
=
ddbSST113
;

ddb121
=
ddbSST121
;

ddb122
=
ddbSST122
;

ddb123
=
ddbSST123
;

ddb131
=
ddbSST131
;

ddb132
=
ddbSST132
;

ddb133
=
ddbSST133
;

ddb221
=
ddbSST221
;

ddb222
=
ddbSST222
;

ddb223
=
ddbSST223
;

ddb231
=
ddbSST231
;

ddb232
=
ddbSST232
;

ddb233
=
ddbSST233
;

ddb331
=
ddbSST331
;

ddb332
=
ddbSST332
;

ddb333
=
ddbSST333
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

deldelg1111
=
deldelgSST1111
;

deldelg1112
=
deldelgSST1112
;

deldelg1113
=
deldelgSST1113
;

deldelg1122
=
deldelgSST1122
;

deldelg1123
=
deldelgSST1123
;

deldelg1133
=
deldelgSST1133
;

deldelg1211
=
deldelgSST1211
;

deldelg1212
=
deldelgSST1212
;

deldelg1213
=
deldelgSST1213
;

deldelg1222
=
deldelgSST1222
;

deldelg1223
=
deldelgSST1223
;

deldelg1233
=
deldelgSST1233
;

deldelg1311
=
deldelgSST1311
;

deldelg1312
=
deldelgSST1312
;

deldelg1313
=
deldelgSST1313
;

deldelg1322
=
deldelgSST1322
;

deldelg1323
=
deldelgSST1323
;

deldelg1333
=
deldelgSST1333
;

deldelg2211
=
deldelgSST2211
;

deldelg2212
=
deldelgSST2212
;

deldelg2213
=
deldelgSST2213
;

deldelg2222
=
deldelgSST2222
;

deldelg2223
=
deldelgSST2223
;

deldelg2233
=
deldelgSST2233
;

deldelg2311
=
deldelgSST2311
;

deldelg2312
=
deldelgSST2312
;

deldelg2313
=
deldelgSST2313
;

deldelg2322
=
deldelgSST2322
;

deldelg2323
=
deldelgSST2323
;

deldelg2333
=
deldelgSST2333
;

deldelg3311
=
deldelgSST3311
;

deldelg3312
=
deldelgSST3312
;

deldelg3313
=
deldelgSST3313
;

deldelg3322
=
deldelgSST3322
;

deldelg3323
=
deldelgSST3323
;

deldelg3333
=
deldelgSST3333
;

delG11
=
delGSST11
;

delG12
=
delGSST12
;

delG13
=
delGSST13
;

delG21
=
delGSST21
;

delG22
=
delGSST22
;

delG23
=
delGSST23
;

delG31
=
delGSST31
;

delG32
=
delGSST32
;

delG33
=
delGSST33
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

ddchi11
=
ddchiSST11
;

ddchi12
=
ddchiSST12
;

ddchi13
=
ddchiSST13
;

ddchi22
=
ddchiSST22
;

ddchi23
=
ddchiSST23
;

ddchi33
=
ddchiSST33
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

lieg11
=
delg111*beta1[ijk] + delg211*beta2[ijk] + delg311*beta3[ijk]
;

lieg12
=
delg112*beta1[ijk] + delg212*beta2[ijk] + delg312*beta3[ijk]
;

lieg13
=
delg113*beta1[ijk] + delg213*beta2[ijk] + delg313*beta3[ijk]
;

lieg22
=
delg122*beta1[ijk] + delg222*beta2[ijk] + delg322*beta3[ijk]
;

lieg23
=
delg123*beta1[ijk] + delg223*beta2[ijk] + delg323*beta3[ijk]
;

lieg33
=
delg133*beta1[ijk] + delg233*beta2[ijk] + delg333*beta3[ijk]
;

lieA11
=
dA111*beta1[ijk] + dA211*beta2[ijk] + dA311*beta3[ijk]
;

lieA12
=
dA112*beta1[ijk] + dA212*beta2[ijk] + dA312*beta3[ijk]
;

lieA13
=
dA113*beta1[ijk] + dA213*beta2[ijk] + dA313*beta3[ijk]
;

lieA22
=
dA122*beta1[ijk] + dA222*beta2[ijk] + dA322*beta3[ijk]
;

lieA23
=
dA123*beta1[ijk] + dA223*beta2[ijk] + dA323*beta3[ijk]
;

lieA33
=
dA133*beta1[ijk] + dA233*beta2[ijk] + dA333*beta3[ijk]
;

lieK
=
dK1*beta1[ijk] + dK2*beta2[ijk] + dK3*beta3[ijk]
;

liechi
=
dchi1*beta1[ijk] + dchi2*beta2[ijk] + dchi3*beta3[ijk]
;

liealpha
=
da1*beta1[ijk] + da2*beta2[ijk] + da3*beta3[ijk]
;

advG1
=
delG11*beta1[ijk] + delG21*beta2[ijk] + delG31*beta3[ijk]
;

advG2
=
delG12*beta1[ijk] + delG22*beta2[ijk] + delG32*beta3[ijk]
;

advG3
=
delG13*beta1[ijk] + delG23*beta2[ijk] + delG33*beta3[ijk]
;

advbeta1
=
db11*beta1[ijk] + db21*beta2[ijk] + db31*beta3[ijk]
;

advbeta2
=
db12*beta1[ijk] + db22*beta2[ijk] + db32*beta3[ijk]
;

advbeta3
=
db13*beta1[ijk] + db23*beta2[ijk] + db33*beta3[ijk]
;

advB1
=
dB11*beta1[ijk] + dB21*beta2[ijk] + dB31*beta3[ijk]
;

advB2
=
dB12*beta1[ijk] + dB22*beta2[ijk] + dB32*beta3[ijk]
;

advB3
=
dB13*beta1[ijk] + dB23*beta2[ijk] + dB33*beta3[ijk]
;


} 

detg
=
2.*g12[ijk]*g13[ijk]*g23[ijk] + g11[ijk]*g22[ijk]*g33[ijk] - 
  g33[ijk]*pow2(g12[ijk]) - g22[ijk]*pow2(g13[ijk]) - g11[ijk]*pow2(g23[ijk])
;


if (detg<=0.01) { PROBLEM; } 

detginv
=
1/detg
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

gammado111
=
0.5*delg111
;

gammado112
=
0.5*delg211
;

gammado113
=
0.5*delg311
;

gammado122
=
-0.5*delg122 + delg212
;

gammado123
=
0.5*(-delg123 + delg213 + delg312)
;

gammado133
=
-0.5*delg133 + delg313
;

gammado211
=
delg112 - 0.5*delg211
;

gammado212
=
0.5*delg122
;

gammado213
=
0.5*(delg123 - delg213 + delg312)
;

gammado222
=
0.5*delg222
;

gammado223
=
0.5*delg322
;

gammado233
=
-0.5*delg233 + delg323
;

gammado311
=
delg113 - 0.5*delg311
;

gammado312
=
0.5*(delg123 + delg213 - delg312)
;

gammado313
=
0.5*delg133
;

gammado322
=
delg223 - 0.5*delg322
;

gammado323
=
0.5*delg233
;

gammado333
=
0.5*delg333
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

R11
=
gammado111*Gfromg1 + gammado112*Gfromg2 + gammado113*Gfromg3 + 
  (-0.5*deldelg1111 + 3.*gamma111*gammado111 + 
     2.*(gamma211*gammado112 + gamma311*gammado113) + 
     gamma211*gammado211 + gamma311*gammado311)*ginv11 + 
  (-deldelg1211 + 3.*(gamma112*gammado111 + gamma111*gammado112) + 
     2.*(gamma212*gammado112 + gamma312*gammado113 + 
        gamma211*gammado122 + gamma311*gammado123) + gamma212*gammado211 + 
     gamma211*gammado212 + gamma312*gammado311 + gamma311*gammado312)*ginv12 \
+ (-deldelg1311 + 3.*(gamma113*gammado111 + gamma111*gammado113) + 
     2.*(gamma213*gammado112 + gamma313*gammado113 + 
        gamma211*gammado123 + gamma311*gammado133) + gamma213*gammado211 + 
     gamma211*gammado213 + gamma313*gammado311 + gamma311*gammado313)*ginv13 \
+ (-0.5*deldelg2211 + 3.*gamma112*gammado112 + 
     2.*(gamma212*gammado122 + gamma312*gammado123) + 
     gamma212*gammado212 + gamma312*gammado312)*ginv22 + 
  (-deldelg2311 + 3.*(gamma113*gammado112 + gamma112*gammado113) + 
     2.*(gamma213*gammado122 + (gamma212 + gamma313)*gammado123 + 
        gamma312*gammado133) + gamma213*gammado212 + gamma212*gammado213 + 
     gamma313*gammado312 + gamma312*gammado313)*ginv23 + 
  (-0.5*deldelg3311 + 3.*gamma113*gammado113 + 
     2.*(gamma213*gammado123 + gamma313*gammado133) + 
     gamma213*gammado213 + gamma313*gammado313)*ginv33 + delG11*g11[ijk] + 
  delG12*g12[ijk] + delG13*g13[ijk]
;

R12
=
(-0.5*deldelg1112 + gamma112*gammado111 + 
     (gamma111 + gamma212)*gammado112 + gamma312*gammado113 + 
     gamma111*gammado211 + 2.*gamma211*gammado212 + 
     gamma311*(gammado213 + gammado312))*ginv11 + 
  (-deldelg1212 + gamma122*gammado111 + 
     (2.*gamma112 + gamma222)*gammado112 + gamma322*gammado113 + 
     (gamma111 + gamma212)*gammado122 + gamma112*gammado211 + 
     (gamma111 + 2.*gamma212)*gammado212 + 2.*gamma211*gammado222 + 
     gamma312*(gammado123 + gammado213 + gammado312) + 
     gamma311*(gammado223 + gammado322))*ginv12 + 
  (-deldelg1312 + gamma123*gammado111 + (gamma113 + gamma223)*gammado112 + 
     (gamma112 + gamma323)*gammado113 + (gamma111 + gamma212)*gammado123 + 
     gamma312*gammado133 + gamma113*gammado211 + 
     (gamma111 + gamma313)*gammado213 + 
     2.*(gamma213*gammado212 + gamma211*gammado223) + 
     gamma313*gammado312 + gamma311*(gammado233 + gammado323))*ginv13 + 
  (-0.5*deldelg2212 + gamma122*gammado112 + 
     (gamma112 + gamma222)*gammado122 + gamma322*gammado123 + 
     gamma112*gammado212 + 2.*gamma212*gammado222 + 
     gamma312*(gammado223 + gammado322))*ginv22 + 
  (-deldelg2312 + gamma123*gammado112 + gamma122*gammado113 + 
     (gamma113 + gamma223)*gammado122 + 
     (gamma112 + gamma222 + gamma323)*gammado123 + gamma322*gammado133 + 
     gamma113*gammado212 + gamma112*gammado213 + 
     2.*(gamma213*gammado222 + gamma212*gammado223) + 
     gamma313*(gammado223 + gammado322) + 
     gamma312*(gammado233 + gammado323))*ginv23 + 
  (-0.5*deldelg3312 + gamma123*gammado113 + 
     (gamma113 + gamma223)*gammado123 + gamma323*gammado133 + 
     gamma113*gammado213 + 2.*gamma213*gammado223 + 
     gamma313*(gammado233 + gammado323))*ginv33 + 
  0.5*((gammado112 + gammado211)*Gfromg1 + 
     (gammado122 + gammado212)*Gfromg2 + (gammado123 + gammado213)*Gfromg3 + 
     delG21*g11[ijk] + (delG11 + delG22)*g12[ijk] + delG23*g13[ijk] + 
     delG12*g22[ijk] + delG13*g23[ijk])
;

R13
=
(-0.5*deldelg1113 + gamma113*gammado111 + gamma213*gammado112 + 
     (gamma111 + gamma313)*gammado113 + gamma111*gammado311 + 
     gamma211*(gammado213 + gammado312) + 2.*gamma311*gammado313)*ginv11 + 
  (-deldelg1213 + gamma123*gammado111 + (gamma113 + gamma223)*gammado112 + 
     (gamma112 + gamma323)*gammado113 + gamma213*gammado122 + 
     (gamma111 + gamma313)*gammado123 + gamma112*gammado311 + 
     gamma111*gammado312 + gamma212*(gammado213 + gammado312) + 
     gamma211*(gammado223 + gammado322) + 
     2.*(gamma312*gammado313 + gamma311*gammado323))*ginv12 + 
  (-deldelg1313 + gamma133*gammado111 + gamma233*gammado112 + 
     (2.*gamma113 + gamma333)*gammado113 + 
     (gamma111 + gamma313)*gammado133 + gamma113*gammado311 + 
     gamma213*(gammado123 + gammado213 + gammado312) + 
     (gamma111 + 2.*gamma313)*gammado313 + 
     gamma211*(gammado233 + gammado323) + 2.*gamma311*gammado333)*ginv13 + 
  (-0.5*deldelg2213 + gamma123*gammado112 + gamma223*gammado122 + 
     (gamma112 + gamma323)*gammado123 + gamma112*gammado312 + 
     gamma212*(gammado223 + gammado322) + 2.*gamma312*gammado323)*ginv22 + 
  (-deldelg2313 + gamma133*gammado112 + gamma123*gammado113 + 
     gamma233*gammado122 + (gamma113 + gamma223 + gamma333)*gammado123 + 
     (gamma112 + gamma323)*gammado133 + gamma113*gammado312 + 
     gamma112*gammado313 + gamma213*(gammado223 + gammado322) + 
     gamma212*(gammado233 + gammado323) + 
     2.*(gamma313*gammado323 + gamma312*gammado333))*ginv23 + 
  (-0.5*deldelg3313 + gamma133*gammado113 + gamma233*gammado123 + 
     (gamma113 + gamma333)*gammado133 + gamma113*gammado313 + 
     gamma213*(gammado233 + gammado323) + 2.*gamma313*gammado333)*ginv33 + 
  0.5*((gammado113 + gammado311)*Gfromg1 + 
     (gammado123 + gammado312)*Gfromg2 + (gammado133 + gammado313)*Gfromg3 + 
     delG31*g11[ijk] + delG32*g12[ijk] + (delG11 + delG33)*g13[ijk] + 
     delG12*g23[ijk] + delG13*g33[ijk])
;

R22
=
gammado212*Gfromg1 + gammado222*Gfromg2 + gammado223*Gfromg3 + 
  (-0.5*deldelg1122 + gamma112*(gammado112 + 2.*gammado211) + 
     3.*gamma212*gammado212 + gamma312*(2.*gammado213 + gammado312))*ginv11 \
+ (-deldelg1222 + gamma122*(gammado112 + 2.*gammado211) + 
     gamma112*(gammado122 + 2.*gammado212) + 
     3.*(gamma222*gammado212 + gamma212*gammado222) + 
     2.*(gamma322*gammado213 + gamma312*gammado223) + 
     gamma322*gammado312 + gamma312*gammado322)*ginv12 + 
  (-deldelg1322 + gamma123*(gammado112 + 2.*gammado211) + 
     gamma112*(gammado123 + 2.*gammado213) + 
     3.*(gamma223*gammado212 + gamma212*gammado223) + 
     2.*(gamma323*gammado213 + gamma312*gammado233) + 
     gamma323*gammado312 + gamma312*gammado323)*ginv13 + 
  (-0.5*deldelg2222 + gamma122*(gammado122 + 2.*gammado212) + 
     3.*gamma222*gammado222 + gamma322*(2.*gammado223 + gammado322))*ginv22 \
+ (-deldelg2322 + gamma123*(gammado122 + 2.*gammado212) + 
     gamma122*(gammado123 + 2.*gammado213) + 
     3.*(gamma223*gammado222 + gamma222*gammado223) + 
     2.*(gamma323*gammado223 + gamma322*gammado233) + 
     gamma323*gammado322 + gamma322*gammado323)*ginv23 + 
  (-0.5*deldelg3322 + gamma123*(gammado123 + 2.*gammado213) + 
     3.*gamma223*gammado223 + gamma323*(2.*gammado233 + gammado323))*ginv33 \
+ delG21*g12[ijk] + delG22*g22[ijk] + delG23*g23[ijk]
;

R23
=
(-0.5*deldelg1123 + gamma113*gammado211 + gamma213*gammado212 + 
     (gamma212 + gamma313)*gammado213 + 
     gamma112*(gammado113 + gammado311) + gamma212*gammado312 + 
     2.*gamma312*gammado313)*ginv11 + 
  (-deldelg1223 + gamma123*gammado211 + (gamma113 + gamma223)*gammado212 + 
     (gamma222 + gamma323)*gammado213 + gamma213*gammado222 + 
     (gamma212 + gamma313)*gammado223 + 
     gamma122*(gammado113 + gammado311) + gamma222*gammado312 + 
     gamma112*(gammado123 + gammado312) + gamma212*gammado322 + 
     2.*(gamma322*gammado313 + gamma312*gammado323))*ginv12 + 
  (-deldelg1323 + gamma133*gammado211 + gamma233*gammado212 + 
     (gamma113 + gamma223 + gamma333)*gammado213 + gamma213*gammado223 + 
     (gamma212 + gamma313)*gammado233 + 
     gamma123*(gammado113 + gammado311) + gamma223*gammado312 + 
     gamma112*(gammado133 + gammado313) + gamma212*gammado323 + 
     2.*(gamma323*gammado313 + gamma312*gammado333))*ginv13 + 
  (-0.5*deldelg2223 + gamma123*gammado212 + gamma223*gammado222 + 
     (gamma222 + gamma323)*gammado223 + 
     gamma122*(gammado123 + gammado312) + gamma222*gammado322 + 
     2.*gamma322*gammado323)*ginv22 + 
  (-deldelg2323 + gamma133*gammado212 + gamma233*gammado222 + 
     (2.*gamma223 + gamma333)*gammado223 + 
     (gamma222 + gamma323)*gammado233 + 
     gamma123*(gammado123 + gammado213 + gammado312) + 
     gamma122*(gammado133 + gammado313) + gamma223*gammado322 + 
     (gamma222 + 2.*gamma323)*gammado323 + 2.*gamma322*gammado333)*ginv23 + 
  (-0.5*deldelg3323 + gamma133*gammado213 + gamma233*gammado223 + 
     (gamma223 + gamma333)*gammado233 + 
     gamma123*(gammado133 + gammado313) + gamma223*gammado323 + 
     2.*gamma323*gammado333)*ginv33 + 
  0.5*((gammado213 + gammado312)*Gfromg1 + 
     (gammado223 + gammado322)*Gfromg2 + (gammado233 + gammado323)*Gfromg3 + 
     delG31*g12[ijk] + delG21*g13[ijk] + delG32*g22[ijk] + 
     (delG22 + delG33)*g23[ijk] + delG23*g33[ijk])
;

R33
=
gammado313*Gfromg1 + gammado323*Gfromg2 + gammado333*Gfromg3 + 
  (-0.5*deldelg1133 + gamma113*(gammado113 + 2.*gammado311) + 
     gamma213*(gammado213 + 2.*gammado312) + 3.*gamma313*gammado313)*ginv11 \
+ (-deldelg1233 + gamma123*(gammado113 + 2.*gammado311) + 
     gamma113*(gammado123 + 2.*gammado312) + 
     gamma223*(gammado213 + 2.*gammado312) + 
     gamma213*(gammado223 + 2.*gammado322) + 
     3.*(gamma323*gammado313 + gamma313*gammado323))*ginv12 + 
  (-deldelg1333 + gamma133*(gammado113 + 2.*gammado311) + 
     gamma233*(gammado213 + 2.*gammado312) + 
     gamma113*(gammado133 + 2.*gammado313) + 
     gamma213*(gammado233 + 2.*gammado323) + 
     3.*(gamma333*gammado313 + gamma313*gammado333))*ginv13 + 
  (-0.5*deldelg2233 + gamma123*(gammado123 + 2.*gammado312) + 
     gamma223*(gammado223 + 2.*gammado322) + 3.*gamma323*gammado323)*ginv22 \
+ (-deldelg2333 + gamma133*(gammado123 + 2.*gammado312) + 
     gamma123*(gammado133 + 2.*gammado313) + 
     gamma233*(gammado223 + 2.*gammado322) + 
     gamma223*(gammado233 + 2.*gammado323) + 
     3.*(gamma333*gammado323 + gamma323*gammado333))*ginv23 + 
  (-0.5*deldelg3333 + gamma133*(gammado133 + 2.*gammado313) + 
     gamma233*(gammado233 + 2.*gammado323) + 3.*gamma333*gammado333)*ginv33 \
+ delG31*g13[ijk] + delG32*g23[ijk] + delG33*g33[ijk]
;

chiguard
=
chiDivFloor
;

chiguarded
=
chi[ijk]
;


if (chiguarded < chiguard) chiguarded = chiguard; 

oochipsipower
=
1/chipsipower
;

f
=
oochipsipower*log(chiguarded)
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

psim2
=
exp(-2.*f)
;

psim4
=
exp(-4.*f)
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

trcdda
=
(cdda11*ginv11 + cdda22*ginv22 + 
    2.*(cdda12*ginv12 + cdda13*ginv13 + cdda23*ginv23) + cdda33*ginv33)*psim4
;

K[ijk]
=
forceKzerofactor*K[ijk]
;

dK1
=
dK1*forceKzerofactor
;

dK2
=
dK2*forceKzerofactor
;

dK3
=
dK3*forceKzerofactor
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

AA33
=
2.*(ginv23*A23[ijk]*A33[ijk] + 
     A13[ijk]*(ginv12*A23[ijk] + ginv13*A33[ijk])) + ginv11*pow2(A13[ijk]) + 
  ginv22*pow2(A23[ijk]) + ginv33*pow2(A33[ijk])
;

AA
=
AA11*ginv11 + AA22*ginv22 + 2.*(AA12*ginv12 + AA13*ginv13 + AA23*ginv23) + 
  AA33*ginv33
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

divAinv1
=
-6.*(Ainv11*df1 + Ainv12*df2 + Ainv13*df3) - Ainv11*gamma111 - 
  Ainv22*gamma122 - 2.*(Ainv12*gamma112 + Ainv13*gamma113 + 
     Ainv23*gamma123) - Ainv33*gamma133 + 
  0.66666666666666666667*(dK1*ginv11 + dK2*ginv12 + dK3*ginv13)
;

divAinv2
=
-6.*(Ainv12*df1 + Ainv22*df2 + Ainv23*df3) - Ainv11*gamma211 - 
  Ainv22*gamma222 - 2.*(Ainv12*gamma212 + Ainv13*gamma213 + 
     Ainv23*gamma223) - Ainv33*gamma233 + 
  0.66666666666666666667*(dK1*ginv12 + dK2*ginv22 + dK3*ginv23)
;

divAinv3
=
-6.*(Ainv13*df1 + Ainv23*df2 + Ainv33*df3) - Ainv11*gamma311 - 
  Ainv22*gamma322 - 2.*(Ainv12*gamma312 + Ainv13*gamma313 + 
     Ainv23*gamma323) - Ainv33*gamma333 + 
  0.66666666666666666667*(dK1*ginv13 + dK2*ginv23 + dK3*ginv33)
;

R
=
AA - 0.66666666666666666667*pow2(K[ijk])
;

divbeta
=
db11 + db22 + db33
;

totdivbeta
=
0.66666666666666666667*divbeta
;

ootddivbeta1
=
0.33333333333333333333*(ddb111 + ddb122 + ddb133)
;

ootddivbeta2
=
0.33333333333333333333*(ddb121 + ddb222 + ddb233)
;

ootddivbeta3
=
0.33333333333333333333*(ddb131 + ddb232 + ddb333)
;

lieg11
=
lieg11 + (2.*db11 - totdivbeta)*g11[ijk] + 2.*(db12*g12[ijk] + db13*g13[ijk])
;

lieg12
=
lieg12 + db21*g11[ijk] + (db11 + db22 - totdivbeta)*g12[ijk] + 
  db23*g13[ijk] + db12*g22[ijk] + db13*g23[ijk]
;

lieg13
=
lieg13 + db31*g11[ijk] + db32*g12[ijk] + 
  (db11 + db33 - totdivbeta)*g13[ijk] + db12*g23[ijk] + db13*g33[ijk]
;

lieg22
=
lieg22 - totdivbeta*g22[ijk] + 
  2.*(db21*g12[ijk] + db22*g22[ijk] + db23*g23[ijk])
;

lieg23
=
lieg23 + db31*g12[ijk] + db21*g13[ijk] + db32*g22[ijk] + 
  (db22 + db33 - totdivbeta)*g23[ijk] + db23*g33[ijk]
;

lieg33
=
lieg33 - totdivbeta*g33[ijk] + 
  2.*(db31*g13[ijk] + db32*g23[ijk] + db33*g33[ijk])
;

lieA11
=
lieA11 + (2.*db11 - totdivbeta)*A11[ijk] + 2.*(db12*A12[ijk] + db13*A13[ijk])
;

lieA12
=
lieA12 + db21*A11[ijk] + (db11 + db22 - totdivbeta)*A12[ijk] + 
  db23*A13[ijk] + db12*A22[ijk] + db13*A23[ijk]
;

lieA13
=
lieA13 + db31*A11[ijk] + db32*A12[ijk] + 
  (db11 + db33 - totdivbeta)*A13[ijk] + db12*A23[ijk] + db13*A33[ijk]
;

lieA22
=
lieA22 - totdivbeta*A22[ijk] + 
  2.*(db21*A12[ijk] + db22*A22[ijk] + db23*A23[ijk])
;

lieA23
=
lieA23 + db31*A12[ijk] + db21*A13[ijk] + db32*A22[ijk] + 
  (db22 + db33 - totdivbeta)*A23[ijk] + db23*A33[ijk]
;

lieA33
=
lieA33 - totdivbeta*A33[ijk] + 
  2.*(db31*A13[ijk] + db32*A23[ijk] + db33*A33[ijk])
;

liechi
=
0.16666666666666666667*chiguarded*chipsipower*divbeta + liechi
;

pseudolieG1
=
advG1 - db11*Gfromg1 - db21*Gfromg2 - db31*Gfromg3 + ddb221*ginv22 + 
  2.*ddb231*ginv23 + ddb331*ginv33 + ginv11*(ddb111 + ootddivbeta1) + 
  ginv12*(2.*ddb121 + ootddivbeta2) + ginv13*(2.*ddb131 + ootddivbeta3) + 
  Gfromg1*totdivbeta
;

pseudolieG2
=
advG2 - db12*Gfromg1 - db22*Gfromg2 - db32*Gfromg3 + ddb112*ginv11 + 
  2.*ddb132*ginv13 + ddb332*ginv33 + ginv12*(2.*ddb122 + ootddivbeta1) + 
  ginv22*(ddb222 + ootddivbeta2) + ginv23*(2.*ddb232 + ootddivbeta3) + 
  Gfromg2*totdivbeta
;

pseudolieG3
=
advG3 - db13*Gfromg1 - db23*Gfromg2 - db33*Gfromg3 + ddb113*ginv11 + 
  2.*ddb123*ginv12 + ddb223*ginv22 + ginv13*(2.*ddb133 + ootddivbeta1) + 
  ginv23*(2.*ddb233 + ootddivbeta2) + ginv33*(ddb333 + ootddivbeta3) + 
  Gfromg3*totdivbeta
;

rg11
=
lieg11 - 2.*A11[ijk]*alpha[ijk]
;

rg12
=
lieg12 - 2.*A12[ijk]*alpha[ijk]
;

rg13
=
lieg13 - 2.*A13[ijk]*alpha[ijk]
;

rg22
=
lieg22 - 2.*A22[ijk]*alpha[ijk]
;

rg23
=
lieg23 - 2.*A23[ijk]*alpha[ijk]
;

rg33
=
lieg33 - 2.*A33[ijk]*alpha[ijk]
;

rA11
=
lieA11 + psim4*(-cdda11 + R11*alpha[ijk]) + 
  0.33333333333333333333*trcdda*g11[ijk] + 
  alpha[ijk]*(-2.*AA11 + psim4*Rphi11 - 0.33333333333333333333*R*g11[ijk] + 
     A11[ijk]*K[ijk])
;

rA12
=
lieA12 + psim4*(-cdda12 + R12*alpha[ijk]) + 
  0.33333333333333333333*trcdda*g12[ijk] + 
  alpha[ijk]*(-2.*AA12 + psim4*Rphi12 - 0.33333333333333333333*R*g12[ijk] + 
     A12[ijk]*K[ijk])
;

rA13
=
lieA13 + psim4*(-cdda13 + R13*alpha[ijk]) + 
  0.33333333333333333333*trcdda*g13[ijk] + 
  alpha[ijk]*(-2.*AA13 + psim4*Rphi13 - 0.33333333333333333333*R*g13[ijk] + 
     A13[ijk]*K[ijk])
;

rA22
=
lieA22 + psim4*(-cdda22 + R22*alpha[ijk]) + 
  0.33333333333333333333*trcdda*g22[ijk] + 
  alpha[ijk]*(-2.*AA22 + psim4*Rphi22 - 0.33333333333333333333*R*g22[ijk] + 
     A22[ijk]*K[ijk])
;

rA23
=
lieA23 + psim4*(-cdda23 + R23*alpha[ijk]) + 
  0.33333333333333333333*trcdda*g23[ijk] + 
  alpha[ijk]*(-2.*AA23 + psim4*Rphi23 - 0.33333333333333333333*R*g23[ijk] + 
     A23[ijk]*K[ijk])
;

rA33
=
lieA33 + psim4*(-cdda33 + R33*alpha[ijk]) + 
  0.33333333333333333333*trcdda*g33[ijk] + 
  alpha[ijk]*(-2.*AA33 + psim4*Rphi33 - 0.33333333333333333333*R*g33[ijk] + 
     A33[ijk]*K[ijk])
;

rG1
=
pseudolieG1 - 2.*(Ainv11*da1 + Ainv12*da2 + Ainv13*da3 + divAinv1*alpha[ijk])
;

rG2
=
pseudolieG2 - 2.*(Ainv12*da1 + Ainv22*da2 + Ainv23*da3 + divAinv2*alpha[ijk])
;

rG3
=
pseudolieG3 - 2.*(Ainv13*da1 + Ainv23*da2 + Ainv33*da3 + divAinv3*alpha[ijk])
;

rK
=
lieK - trcdda + alpha[ijk]*(AA + 0.33333333333333333333*pow2(K[ijk]))
;

rchi
=
liechi - 0.16666666666666666667*chiguarded*chipsipower*alpha[ijk]*K[ijk]
;

falpha
=
lapseharmonicf*oploglapse + (8.*oploglapse2)/(9. - 3.*alpha[ijk]) + 
  harmoniclapse*alpha[ijk]
;

ralpha
=
liealpha*oplogwithshift - falpha*alpha[ijk]*K[ijk]
;

betaF
=
shiftgammacoeff*Power(alpha[ijk],shiftalphapower)
;

rbeta1
=
advbeta1*withShiftadv - shiftdriver*beta1[ijk] + 
  withB*((betaF*gamma0factor + gamma2factor)*B1[ijk] + 
     shiftdriver*beta1[ijk]) + 
  shiftgammacoeff*Power(alpha[ijk],shiftGalphapower)*G1[ijk] - 
  shiftgammacoeff*withB*Power(alpha[ijk],shiftGalphapower)*G1[ijk]
;

rbeta2
=
advbeta2*withShiftadv - shiftdriver*beta2[ijk] + 
  withB*((betaF*gamma0factor + gamma2factor)*B2[ijk] + 
     shiftdriver*beta2[ijk]) + 
  shiftgammacoeff*Power(alpha[ijk],shiftGalphapower)*G2[ijk] - 
  shiftgammacoeff*withB*Power(alpha[ijk],shiftGalphapower)*G2[ijk]
;

rbeta3
=
advbeta3*withShiftadv - shiftdriver*beta3[ijk] + 
  withB*((betaF*gamma0factor + gamma2factor)*B3[ijk] + 
     shiftdriver*beta3[ijk]) + 
  shiftgammacoeff*Power(alpha[ijk],shiftGalphapower)*G3[ijk] - 
  shiftgammacoeff*withB*Power(alpha[ijk],shiftGalphapower)*G3[ijk]
;


if (use_eta) { 

dpsim21
=
-2.*df1*psim2
;

dpsim22
=
-2.*df2*psim2
;

dpsim23
=
-2.*df3*psim2
;

absdpsim2
=
Sqrt(2.*(dpsim21*(dpsim22*ginv12 + dpsim23*ginv13) + 
      dpsim22*dpsim23*ginv23) + ginv11*pow2(dpsim21) + 
   ginv22*pow2(dpsim22) + ginv33*pow2(dpsim23))
;


shiftDr[ijk] = bssn_eta_set(xp[ijk],yp[ijk],zp[ijk], absdpsim2,psim2); 

rB1
=
(gamma0factor + betaF*gamma2factor)*rG1 + advB1*withBadv - 
  advG1*(gamma0factor + betaF*gamma2factor)*withGadv - B1[ijk]*shiftDr[ijk]
;

rB2
=
(gamma0factor + betaF*gamma2factor)*rG2 + advB2*withBadv - 
  advG2*(gamma0factor + betaF*gamma2factor)*withGadv - B2[ijk]*shiftDr[ijk]
;

rB3
=
(gamma0factor + betaF*gamma2factor)*rG3 + advB3*withBadv - 
  advG3*(gamma0factor + betaF*gamma2factor)*withGadv - B3[ijk]*shiftDr[ijk]
;


} else { 

rB1
=
withB*((gamma0factor + betaF*gamma2factor)*rG1 + advB1*withBadv - 
    advG1*(gamma0factor + betaF*gamma2factor)*withGadv - shiftdriver*B1[ijk])
;

rB2
=
withB*((gamma0factor + betaF*gamma2factor)*rG2 + advB2*withBadv - 
    advG2*(gamma0factor + betaF*gamma2factor)*withGadv - shiftdriver*B2[ijk])
;

rB3
=
withB*((gamma0factor + betaF*gamma2factor)*rG3 + advB3*withBadv - 
    advG3*(gamma0factor + betaF*gamma2factor)*withGadv - shiftdriver*B3[ijk])
;


} 


if (order_dissipation == 4 && boundaryNormore(2)) { 

rg11
=
rg11 - 2.*dissfactor*(oo2dx*(6.*g11[ijk] + g11[-2*di + ijk] - 
        4.*(g11[-di + ijk] + g11[di + ijk]) + g11[2*di + ijk]) + 
     oo2dy*(6.*g11[ijk] + g11[-2*dj + ijk] - 
        4.*(g11[-dj + ijk] + g11[dj + ijk]) + g11[2*dj + ijk]) + 
     oo2dz*(6.*g11[ijk] + g11[-2*dk + ijk] - 
        4.*(g11[-dk + ijk] + g11[dk + ijk]) + g11[2*dk + ijk]))
;

rg12
=
rg12 - 2.*dissfactor*(oo2dx*(6.*g12[ijk] + g12[-2*di + ijk] - 
        4.*(g12[-di + ijk] + g12[di + ijk]) + g12[2*di + ijk]) + 
     oo2dy*(6.*g12[ijk] + g12[-2*dj + ijk] - 
        4.*(g12[-dj + ijk] + g12[dj + ijk]) + g12[2*dj + ijk]) + 
     oo2dz*(6.*g12[ijk] + g12[-2*dk + ijk] - 
        4.*(g12[-dk + ijk] + g12[dk + ijk]) + g12[2*dk + ijk]))
;

rg13
=
rg13 - 2.*dissfactor*(oo2dx*(6.*g13[ijk] + g13[-2*di + ijk] - 
        4.*(g13[-di + ijk] + g13[di + ijk]) + g13[2*di + ijk]) + 
     oo2dy*(6.*g13[ijk] + g13[-2*dj + ijk] - 
        4.*(g13[-dj + ijk] + g13[dj + ijk]) + g13[2*dj + ijk]) + 
     oo2dz*(6.*g13[ijk] + g13[-2*dk + ijk] - 
        4.*(g13[-dk + ijk] + g13[dk + ijk]) + g13[2*dk + ijk]))
;

rg22
=
rg22 - 2.*dissfactor*(oo2dx*(6.*g22[ijk] + g22[-2*di + ijk] - 
        4.*(g22[-di + ijk] + g22[di + ijk]) + g22[2*di + ijk]) + 
     oo2dy*(6.*g22[ijk] + g22[-2*dj + ijk] - 
        4.*(g22[-dj + ijk] + g22[dj + ijk]) + g22[2*dj + ijk]) + 
     oo2dz*(6.*g22[ijk] + g22[-2*dk + ijk] - 
        4.*(g22[-dk + ijk] + g22[dk + ijk]) + g22[2*dk + ijk]))
;

rg23
=
rg23 - 2.*dissfactor*(oo2dx*(6.*g23[ijk] + g23[-2*di + ijk] - 
        4.*(g23[-di + ijk] + g23[di + ijk]) + g23[2*di + ijk]) + 
     oo2dy*(6.*g23[ijk] + g23[-2*dj + ijk] - 
        4.*(g23[-dj + ijk] + g23[dj + ijk]) + g23[2*dj + ijk]) + 
     oo2dz*(6.*g23[ijk] + g23[-2*dk + ijk] - 
        4.*(g23[-dk + ijk] + g23[dk + ijk]) + g23[2*dk + ijk]))
;

rg33
=
rg33 - 2.*dissfactor*(oo2dx*(6.*g33[ijk] + g33[-2*di + ijk] - 
        4.*(g33[-di + ijk] + g33[di + ijk]) + g33[2*di + ijk]) + 
     oo2dy*(6.*g33[ijk] + g33[-2*dj + ijk] - 
        4.*(g33[-dj + ijk] + g33[dj + ijk]) + g33[2*dj + ijk]) + 
     oo2dz*(6.*g33[ijk] + g33[-2*dk + ijk] - 
        4.*(g33[-dk + ijk] + g33[dk + ijk]) + g33[2*dk + ijk]))
;

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

rchi
=
rchi - 2.*dissfactor*(oo2dx*(6.*chi[ijk] + chi[-2*di + ijk] - 
        4.*(chi[-di + ijk] + chi[di + ijk]) + chi[2*di + ijk]) + 
     oo2dy*(6.*chi[ijk] + chi[-2*dj + ijk] - 
        4.*(chi[-dj + ijk] + chi[dj + ijk]) + chi[2*dj + ijk]) + 
     oo2dz*(6.*chi[ijk] + chi[-2*dk + ijk] - 
        4.*(chi[-dk + ijk] + chi[dk + ijk]) + chi[2*dk + ijk]))
;


} 


if (order_dissipation == 6 && boundaryNormore(3)) { 


#ifdef REDUCEORDERTO2                           
  errorexit("CompilerFlag reduced order to 2");
 #else
rg11
=
rg11 + 2.*dissfactor*(oo2dx*(-20.*g11[ijk] + g11[-3*di + ijk] + 
        15.*(g11[-di + ijk] + g11[di + ijk]) - 
        6.*(g11[-2*di + ijk] + g11[2*di + ijk]) + g11[3*di + ijk]) + 
     oo2dy*(-20.*g11[ijk] + g11[-3*dj + ijk] + 
        15.*(g11[-dj + ijk] + g11[dj + ijk]) - 
        6.*(g11[-2*dj + ijk] + g11[2*dj + ijk]) + g11[3*dj + ijk]) + 
     oo2dz*(-20.*g11[ijk] + g11[-3*dk + ijk] + 
        15.*(g11[-dk + ijk] + g11[dk + ijk]) - 
        6.*(g11[-2*dk + ijk] + g11[2*dk + ijk]) + g11[3*dk + ijk]))
;

rg12
=
rg12 + 2.*dissfactor*(oo2dx*(-20.*g12[ijk] + g12[-3*di + ijk] + 
        15.*(g12[-di + ijk] + g12[di + ijk]) - 
        6.*(g12[-2*di + ijk] + g12[2*di + ijk]) + g12[3*di + ijk]) + 
     oo2dy*(-20.*g12[ijk] + g12[-3*dj + ijk] + 
        15.*(g12[-dj + ijk] + g12[dj + ijk]) - 
        6.*(g12[-2*dj + ijk] + g12[2*dj + ijk]) + g12[3*dj + ijk]) + 
     oo2dz*(-20.*g12[ijk] + g12[-3*dk + ijk] + 
        15.*(g12[-dk + ijk] + g12[dk + ijk]) - 
        6.*(g12[-2*dk + ijk] + g12[2*dk + ijk]) + g12[3*dk + ijk]))
;

rg13
=
rg13 + 2.*dissfactor*(oo2dx*(-20.*g13[ijk] + g13[-3*di + ijk] + 
        15.*(g13[-di + ijk] + g13[di + ijk]) - 
        6.*(g13[-2*di + ijk] + g13[2*di + ijk]) + g13[3*di + ijk]) + 
     oo2dy*(-20.*g13[ijk] + g13[-3*dj + ijk] + 
        15.*(g13[-dj + ijk] + g13[dj + ijk]) - 
        6.*(g13[-2*dj + ijk] + g13[2*dj + ijk]) + g13[3*dj + ijk]) + 
     oo2dz*(-20.*g13[ijk] + g13[-3*dk + ijk] + 
        15.*(g13[-dk + ijk] + g13[dk + ijk]) - 
        6.*(g13[-2*dk + ijk] + g13[2*dk + ijk]) + g13[3*dk + ijk]))
;

rg22
=
rg22 + 2.*dissfactor*(oo2dx*(-20.*g22[ijk] + g22[-3*di + ijk] + 
        15.*(g22[-di + ijk] + g22[di + ijk]) - 
        6.*(g22[-2*di + ijk] + g22[2*di + ijk]) + g22[3*di + ijk]) + 
     oo2dy*(-20.*g22[ijk] + g22[-3*dj + ijk] + 
        15.*(g22[-dj + ijk] + g22[dj + ijk]) - 
        6.*(g22[-2*dj + ijk] + g22[2*dj + ijk]) + g22[3*dj + ijk]) + 
     oo2dz*(-20.*g22[ijk] + g22[-3*dk + ijk] + 
        15.*(g22[-dk + ijk] + g22[dk + ijk]) - 
        6.*(g22[-2*dk + ijk] + g22[2*dk + ijk]) + g22[3*dk + ijk]))
;

rg23
=
rg23 + 2.*dissfactor*(oo2dx*(-20.*g23[ijk] + g23[-3*di + ijk] + 
        15.*(g23[-di + ijk] + g23[di + ijk]) - 
        6.*(g23[-2*di + ijk] + g23[2*di + ijk]) + g23[3*di + ijk]) + 
     oo2dy*(-20.*g23[ijk] + g23[-3*dj + ijk] + 
        15.*(g23[-dj + ijk] + g23[dj + ijk]) - 
        6.*(g23[-2*dj + ijk] + g23[2*dj + ijk]) + g23[3*dj + ijk]) + 
     oo2dz*(-20.*g23[ijk] + g23[-3*dk + ijk] + 
        15.*(g23[-dk + ijk] + g23[dk + ijk]) - 
        6.*(g23[-2*dk + ijk] + g23[2*dk + ijk]) + g23[3*dk + ijk]))
;

rg33
=
rg33 + 2.*dissfactor*(oo2dx*(-20.*g33[ijk] + g33[-3*di + ijk] + 
        15.*(g33[-di + ijk] + g33[di + ijk]) - 
        6.*(g33[-2*di + ijk] + g33[2*di + ijk]) + g33[3*di + ijk]) + 
     oo2dy*(-20.*g33[ijk] + g33[-3*dj + ijk] + 
        15.*(g33[-dj + ijk] + g33[dj + ijk]) - 
        6.*(g33[-2*dj + ijk] + g33[2*dj + ijk]) + g33[3*dj + ijk]) + 
     oo2dz*(-20.*g33[ijk] + g33[-3*dk + ijk] + 
        15.*(g33[-dk + ijk] + g33[dk + ijk]) - 
        6.*(g33[-2*dk + ijk] + g33[2*dk + ijk]) + g33[3*dk + ijk]))
;

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

rchi
=
rchi + 2.*dissfactor*(oo2dx*(-20.*chi[ijk] + chi[-3*di + ijk] + 
        15.*(chi[-di + ijk] + chi[di + ijk]) - 
        6.*(chi[-2*di + ijk] + chi[2*di + ijk]) + chi[3*di + ijk]) + 
     oo2dy*(-20.*chi[ijk] + chi[-3*dj + ijk] + 
        15.*(chi[-dj + ijk] + chi[dj + ijk]) - 
        6.*(chi[-2*dj + ijk] + chi[2*dj + ijk]) + chi[3*dj + ijk]) + 
     oo2dz*(-20.*chi[ijk] + chi[-3*dk + ijk] + 
        15.*(chi[-dk + ijk] + chi[dk + ijk]) - 
        6.*(chi[-2*dk + ijk] + chi[2*dk + ijk]) + chi[3*dk + ijk]))
;


#endif 


} 


if (order_dissipation == 8 && boundaryNormore(4)) { 


#ifdef REDUCEORDERTO4                           
  errorexit("CompilerFlag reduced order to 4");
 #else
rg11
=
rg11 - 2.*dissfactor*(oo2dx*(70.*g11[ijk] + g11[-4*di + ijk] - 
        56.*(g11[-di + ijk] + g11[di + ijk]) + 
        28.*(g11[-2*di + ijk] + g11[2*di + ijk]) - 
        8.*(g11[-3*di + ijk] + g11[3*di + ijk]) + g11[4*di + ijk]) + 
     oo2dy*(70.*g11[ijk] + g11[-4*dj + ijk] - 
        56.*(g11[-dj + ijk] + g11[dj + ijk]) + 
        28.*(g11[-2*dj + ijk] + g11[2*dj + ijk]) - 
        8.*(g11[-3*dj + ijk] + g11[3*dj + ijk]) + g11[4*dj + ijk]) + 
     oo2dz*(70.*g11[ijk] + g11[-4*dk + ijk] - 
        56.*(g11[-dk + ijk] + g11[dk + ijk]) + 
        28.*(g11[-2*dk + ijk] + g11[2*dk + ijk]) - 
        8.*(g11[-3*dk + ijk] + g11[3*dk + ijk]) + g11[4*dk + ijk]))
;

rg12
=
rg12 - 2.*dissfactor*(oo2dx*(70.*g12[ijk] + g12[-4*di + ijk] - 
        56.*(g12[-di + ijk] + g12[di + ijk]) + 
        28.*(g12[-2*di + ijk] + g12[2*di + ijk]) - 
        8.*(g12[-3*di + ijk] + g12[3*di + ijk]) + g12[4*di + ijk]) + 
     oo2dy*(70.*g12[ijk] + g12[-4*dj + ijk] - 
        56.*(g12[-dj + ijk] + g12[dj + ijk]) + 
        28.*(g12[-2*dj + ijk] + g12[2*dj + ijk]) - 
        8.*(g12[-3*dj + ijk] + g12[3*dj + ijk]) + g12[4*dj + ijk]) + 
     oo2dz*(70.*g12[ijk] + g12[-4*dk + ijk] - 
        56.*(g12[-dk + ijk] + g12[dk + ijk]) + 
        28.*(g12[-2*dk + ijk] + g12[2*dk + ijk]) - 
        8.*(g12[-3*dk + ijk] + g12[3*dk + ijk]) + g12[4*dk + ijk]))
;

rg13
=
rg13 - 2.*dissfactor*(oo2dx*(70.*g13[ijk] + g13[-4*di + ijk] - 
        56.*(g13[-di + ijk] + g13[di + ijk]) + 
        28.*(g13[-2*di + ijk] + g13[2*di + ijk]) - 
        8.*(g13[-3*di + ijk] + g13[3*di + ijk]) + g13[4*di + ijk]) + 
     oo2dy*(70.*g13[ijk] + g13[-4*dj + ijk] - 
        56.*(g13[-dj + ijk] + g13[dj + ijk]) + 
        28.*(g13[-2*dj + ijk] + g13[2*dj + ijk]) - 
        8.*(g13[-3*dj + ijk] + g13[3*dj + ijk]) + g13[4*dj + ijk]) + 
     oo2dz*(70.*g13[ijk] + g13[-4*dk + ijk] - 
        56.*(g13[-dk + ijk] + g13[dk + ijk]) + 
        28.*(g13[-2*dk + ijk] + g13[2*dk + ijk]) - 
        8.*(g13[-3*dk + ijk] + g13[3*dk + ijk]) + g13[4*dk + ijk]))
;

rg22
=
rg22 - 2.*dissfactor*(oo2dx*(70.*g22[ijk] + g22[-4*di + ijk] - 
        56.*(g22[-di + ijk] + g22[di + ijk]) + 
        28.*(g22[-2*di + ijk] + g22[2*di + ijk]) - 
        8.*(g22[-3*di + ijk] + g22[3*di + ijk]) + g22[4*di + ijk]) + 
     oo2dy*(70.*g22[ijk] + g22[-4*dj + ijk] - 
        56.*(g22[-dj + ijk] + g22[dj + ijk]) + 
        28.*(g22[-2*dj + ijk] + g22[2*dj + ijk]) - 
        8.*(g22[-3*dj + ijk] + g22[3*dj + ijk]) + g22[4*dj + ijk]) + 
     oo2dz*(70.*g22[ijk] + g22[-4*dk + ijk] - 
        56.*(g22[-dk + ijk] + g22[dk + ijk]) + 
        28.*(g22[-2*dk + ijk] + g22[2*dk + ijk]) - 
        8.*(g22[-3*dk + ijk] + g22[3*dk + ijk]) + g22[4*dk + ijk]))
;

rg23
=
rg23 - 2.*dissfactor*(oo2dx*(70.*g23[ijk] + g23[-4*di + ijk] - 
        56.*(g23[-di + ijk] + g23[di + ijk]) + 
        28.*(g23[-2*di + ijk] + g23[2*di + ijk]) - 
        8.*(g23[-3*di + ijk] + g23[3*di + ijk]) + g23[4*di + ijk]) + 
     oo2dy*(70.*g23[ijk] + g23[-4*dj + ijk] - 
        56.*(g23[-dj + ijk] + g23[dj + ijk]) + 
        28.*(g23[-2*dj + ijk] + g23[2*dj + ijk]) - 
        8.*(g23[-3*dj + ijk] + g23[3*dj + ijk]) + g23[4*dj + ijk]) + 
     oo2dz*(70.*g23[ijk] + g23[-4*dk + ijk] - 
        56.*(g23[-dk + ijk] + g23[dk + ijk]) + 
        28.*(g23[-2*dk + ijk] + g23[2*dk + ijk]) - 
        8.*(g23[-3*dk + ijk] + g23[3*dk + ijk]) + g23[4*dk + ijk]))
;

rg33
=
rg33 - 2.*dissfactor*(oo2dx*(70.*g33[ijk] + g33[-4*di + ijk] - 
        56.*(g33[-di + ijk] + g33[di + ijk]) + 
        28.*(g33[-2*di + ijk] + g33[2*di + ijk]) - 
        8.*(g33[-3*di + ijk] + g33[3*di + ijk]) + g33[4*di + ijk]) + 
     oo2dy*(70.*g33[ijk] + g33[-4*dj + ijk] - 
        56.*(g33[-dj + ijk] + g33[dj + ijk]) + 
        28.*(g33[-2*dj + ijk] + g33[2*dj + ijk]) - 
        8.*(g33[-3*dj + ijk] + g33[3*dj + ijk]) + g33[4*dj + ijk]) + 
     oo2dz*(70.*g33[ijk] + g33[-4*dk + ijk] - 
        56.*(g33[-dk + ijk] + g33[dk + ijk]) + 
        28.*(g33[-2*dk + ijk] + g33[2*dk + ijk]) - 
        8.*(g33[-3*dk + ijk] + g33[3*dk + ijk]) + g33[4*dk + ijk]))
;

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

rchi
=
rchi - 2.*dissfactor*(oo2dx*(70.*chi[ijk] + chi[-4*di + ijk] - 
        56.*(chi[-di + ijk] + chi[di + ijk]) + 
        28.*(chi[-2*di + ijk] + chi[2*di + ijk]) - 
        8.*(chi[-3*di + ijk] + chi[3*di + ijk]) + chi[4*di + ijk]) + 
     oo2dy*(70.*chi[ijk] + chi[-4*dj + ijk] - 
        56.*(chi[-dj + ijk] + chi[dj + ijk]) + 
        28.*(chi[-2*dj + ijk] + chi[2*dj + ijk]) - 
        8.*(chi[-3*dj + ijk] + chi[3*dj + ijk]) + chi[4*dj + ijk]) + 
     oo2dz*(70.*chi[ijk] + chi[-4*dk + ijk] - 
        56.*(chi[-dk + ijk] + chi[dk + ijk]) + 
        28.*(chi[-2*dk + ijk] + chi[2*dk + ijk]) - 
        8.*(chi[-3*dk + ijk] + chi[3*dk + ijk]) + chi[4*dk + ijk]))
;


#endif 


} 


if (order_dissipation == 10 && boundaryNormore(5)) { 


#ifdef REDUCEORDERTO6                           
  errorexit("CompilerFlag reduced order to 6");
 #else
rg11
=
rg11 + 2.*dissfactor*(oo2dx*(-252.*g11[ijk] + g11[-5*di + ijk] + 
        210.*(g11[-di + ijk] + g11[di + ijk]) - 
        120.*(g11[-2*di + ijk] + g11[2*di + ijk]) + 
        45.*(g11[-3*di + ijk] + g11[3*di + ijk]) - 
        10.*(g11[-4*di + ijk] + g11[4*di + ijk]) + g11[5*di + ijk]) + 
     oo2dy*(-252.*g11[ijk] + g11[-5*dj + ijk] + 
        210.*(g11[-dj + ijk] + g11[dj + ijk]) - 
        120.*(g11[-2*dj + ijk] + g11[2*dj + ijk]) + 
        45.*(g11[-3*dj + ijk] + g11[3*dj + ijk]) - 
        10.*(g11[-4*dj + ijk] + g11[4*dj + ijk]) + g11[5*dj + ijk]) + 
     oo2dz*(-252.*g11[ijk] + g11[-5*dk + ijk] + 
        210.*(g11[-dk + ijk] + g11[dk + ijk]) - 
        120.*(g11[-2*dk + ijk] + g11[2*dk + ijk]) + 
        45.*(g11[-3*dk + ijk] + g11[3*dk + ijk]) - 
        10.*(g11[-4*dk + ijk] + g11[4*dk + ijk]) + g11[5*dk + ijk]))
;

rg12
=
rg12 + 2.*dissfactor*(oo2dx*(-252.*g12[ijk] + g12[-5*di + ijk] + 
        210.*(g12[-di + ijk] + g12[di + ijk]) - 
        120.*(g12[-2*di + ijk] + g12[2*di + ijk]) + 
        45.*(g12[-3*di + ijk] + g12[3*di + ijk]) - 
        10.*(g12[-4*di + ijk] + g12[4*di + ijk]) + g12[5*di + ijk]) + 
     oo2dy*(-252.*g12[ijk] + g12[-5*dj + ijk] + 
        210.*(g12[-dj + ijk] + g12[dj + ijk]) - 
        120.*(g12[-2*dj + ijk] + g12[2*dj + ijk]) + 
        45.*(g12[-3*dj + ijk] + g12[3*dj + ijk]) - 
        10.*(g12[-4*dj + ijk] + g12[4*dj + ijk]) + g12[5*dj + ijk]) + 
     oo2dz*(-252.*g12[ijk] + g12[-5*dk + ijk] + 
        210.*(g12[-dk + ijk] + g12[dk + ijk]) - 
        120.*(g12[-2*dk + ijk] + g12[2*dk + ijk]) + 
        45.*(g12[-3*dk + ijk] + g12[3*dk + ijk]) - 
        10.*(g12[-4*dk + ijk] + g12[4*dk + ijk]) + g12[5*dk + ijk]))
;

rg13
=
rg13 + 2.*dissfactor*(oo2dx*(-252.*g13[ijk] + g13[-5*di + ijk] + 
        210.*(g13[-di + ijk] + g13[di + ijk]) - 
        120.*(g13[-2*di + ijk] + g13[2*di + ijk]) + 
        45.*(g13[-3*di + ijk] + g13[3*di + ijk]) - 
        10.*(g13[-4*di + ijk] + g13[4*di + ijk]) + g13[5*di + ijk]) + 
     oo2dy*(-252.*g13[ijk] + g13[-5*dj + ijk] + 
        210.*(g13[-dj + ijk] + g13[dj + ijk]) - 
        120.*(g13[-2*dj + ijk] + g13[2*dj + ijk]) + 
        45.*(g13[-3*dj + ijk] + g13[3*dj + ijk]) - 
        10.*(g13[-4*dj + ijk] + g13[4*dj + ijk]) + g13[5*dj + ijk]) + 
     oo2dz*(-252.*g13[ijk] + g13[-5*dk + ijk] + 
        210.*(g13[-dk + ijk] + g13[dk + ijk]) - 
        120.*(g13[-2*dk + ijk] + g13[2*dk + ijk]) + 
        45.*(g13[-3*dk + ijk] + g13[3*dk + ijk]) - 
        10.*(g13[-4*dk + ijk] + g13[4*dk + ijk]) + g13[5*dk + ijk]))
;

rg22
=
rg22 + 2.*dissfactor*(oo2dx*(-252.*g22[ijk] + g22[-5*di + ijk] + 
        210.*(g22[-di + ijk] + g22[di + ijk]) - 
        120.*(g22[-2*di + ijk] + g22[2*di + ijk]) + 
        45.*(g22[-3*di + ijk] + g22[3*di + ijk]) - 
        10.*(g22[-4*di + ijk] + g22[4*di + ijk]) + g22[5*di + ijk]) + 
     oo2dy*(-252.*g22[ijk] + g22[-5*dj + ijk] + 
        210.*(g22[-dj + ijk] + g22[dj + ijk]) - 
        120.*(g22[-2*dj + ijk] + g22[2*dj + ijk]) + 
        45.*(g22[-3*dj + ijk] + g22[3*dj + ijk]) - 
        10.*(g22[-4*dj + ijk] + g22[4*dj + ijk]) + g22[5*dj + ijk]) + 
     oo2dz*(-252.*g22[ijk] + g22[-5*dk + ijk] + 
        210.*(g22[-dk + ijk] + g22[dk + ijk]) - 
        120.*(g22[-2*dk + ijk] + g22[2*dk + ijk]) + 
        45.*(g22[-3*dk + ijk] + g22[3*dk + ijk]) - 
        10.*(g22[-4*dk + ijk] + g22[4*dk + ijk]) + g22[5*dk + ijk]))
;

rg23
=
rg23 + 2.*dissfactor*(oo2dx*(-252.*g23[ijk] + g23[-5*di + ijk] + 
        210.*(g23[-di + ijk] + g23[di + ijk]) - 
        120.*(g23[-2*di + ijk] + g23[2*di + ijk]) + 
        45.*(g23[-3*di + ijk] + g23[3*di + ijk]) - 
        10.*(g23[-4*di + ijk] + g23[4*di + ijk]) + g23[5*di + ijk]) + 
     oo2dy*(-252.*g23[ijk] + g23[-5*dj + ijk] + 
        210.*(g23[-dj + ijk] + g23[dj + ijk]) - 
        120.*(g23[-2*dj + ijk] + g23[2*dj + ijk]) + 
        45.*(g23[-3*dj + ijk] + g23[3*dj + ijk]) - 
        10.*(g23[-4*dj + ijk] + g23[4*dj + ijk]) + g23[5*dj + ijk]) + 
     oo2dz*(-252.*g23[ijk] + g23[-5*dk + ijk] + 
        210.*(g23[-dk + ijk] + g23[dk + ijk]) - 
        120.*(g23[-2*dk + ijk] + g23[2*dk + ijk]) + 
        45.*(g23[-3*dk + ijk] + g23[3*dk + ijk]) - 
        10.*(g23[-4*dk + ijk] + g23[4*dk + ijk]) + g23[5*dk + ijk]))
;

rg33
=
rg33 + 2.*dissfactor*(oo2dx*(-252.*g33[ijk] + g33[-5*di + ijk] + 
        210.*(g33[-di + ijk] + g33[di + ijk]) - 
        120.*(g33[-2*di + ijk] + g33[2*di + ijk]) + 
        45.*(g33[-3*di + ijk] + g33[3*di + ijk]) - 
        10.*(g33[-4*di + ijk] + g33[4*di + ijk]) + g33[5*di + ijk]) + 
     oo2dy*(-252.*g33[ijk] + g33[-5*dj + ijk] + 
        210.*(g33[-dj + ijk] + g33[dj + ijk]) - 
        120.*(g33[-2*dj + ijk] + g33[2*dj + ijk]) + 
        45.*(g33[-3*dj + ijk] + g33[3*dj + ijk]) - 
        10.*(g33[-4*dj + ijk] + g33[4*dj + ijk]) + g33[5*dj + ijk]) + 
     oo2dz*(-252.*g33[ijk] + g33[-5*dk + ijk] + 
        210.*(g33[-dk + ijk] + g33[dk + ijk]) - 
        120.*(g33[-2*dk + ijk] + g33[2*dk + ijk]) + 
        45.*(g33[-3*dk + ijk] + g33[3*dk + ijk]) - 
        10.*(g33[-4*dk + ijk] + g33[4*dk + ijk]) + g33[5*dk + ijk]))
;

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

rchi
=
rchi + 2.*dissfactor*(oo2dx*(-252.*chi[ijk] + chi[-5*di + ijk] + 
        210.*(chi[-di + ijk] + chi[di + ijk]) - 
        120.*(chi[-2*di + ijk] + chi[2*di + ijk]) + 
        45.*(chi[-3*di + ijk] + chi[3*di + ijk]) - 
        10.*(chi[-4*di + ijk] + chi[4*di + ijk]) + chi[5*di + ijk]) + 
     oo2dy*(-252.*chi[ijk] + chi[-5*dj + ijk] + 
        210.*(chi[-dj + ijk] + chi[dj + ijk]) - 
        120.*(chi[-2*dj + ijk] + chi[2*dj + ijk]) + 
        45.*(chi[-3*dj + ijk] + chi[3*dj + ijk]) - 
        10.*(chi[-4*dj + ijk] + chi[4*dj + ijk]) + chi[5*dj + ijk]) + 
     oo2dz*(-252.*chi[ijk] + chi[-5*dk + ijk] + 
        210.*(chi[-dk + ijk] + chi[dk + ijk]) - 
        120.*(chi[-2*dk + ijk] + chi[2*dk + ijk]) + 
        45.*(chi[-3*dk + ijk] + chi[3*dk + ijk]) - 
        10.*(chi[-4*dk + ijk] + chi[4*dk + ijk]) + chi[5*dk + ijk]))
;


#endif 


} 


if (order_dissipation == 12 && boundaryNormore(6)) { 


#ifdef REDUCEORDERTO8                           
  errorexit("CompilerFlag reduced order to 8");
 #else
rg11
=
rg11 - 2.*dissfactor*(oo2dx*(924.*g11[ijk] + g11[-6*di + ijk] - 
        792.*(g11[-di + ijk] + g11[di + ijk]) + 
        495.*(g11[-2*di + ijk] + g11[2*di + ijk]) - 
        220.*(g11[-3*di + ijk] + g11[3*di + ijk]) + 
        66.*(g11[-4*di + ijk] + g11[4*di + ijk]) - 
        12.*(g11[-5*di + ijk] + g11[5*di + ijk]) + g11[6*di + ijk]) + 
     oo2dy*(924.*g11[ijk] + g11[-6*dj + ijk] - 
        792.*(g11[-dj + ijk] + g11[dj + ijk]) + 
        495.*(g11[-2*dj + ijk] + g11[2*dj + ijk]) - 
        220.*(g11[-3*dj + ijk] + g11[3*dj + ijk]) + 
        66.*(g11[-4*dj + ijk] + g11[4*dj + ijk]) - 
        12.*(g11[-5*dj + ijk] + g11[5*dj + ijk]) + g11[6*dj + ijk]) + 
     oo2dz*(924.*g11[ijk] + g11[-6*dk + ijk] - 
        792.*(g11[-dk + ijk] + g11[dk + ijk]) + 
        495.*(g11[-2*dk + ijk] + g11[2*dk + ijk]) - 
        220.*(g11[-3*dk + ijk] + g11[3*dk + ijk]) + 
        66.*(g11[-4*dk + ijk] + g11[4*dk + ijk]) - 
        12.*(g11[-5*dk + ijk] + g11[5*dk + ijk]) + g11[6*dk + ijk]))
;

rg12
=
rg12 - 2.*dissfactor*(oo2dx*(924.*g12[ijk] + g12[-6*di + ijk] - 
        792.*(g12[-di + ijk] + g12[di + ijk]) + 
        495.*(g12[-2*di + ijk] + g12[2*di + ijk]) - 
        220.*(g12[-3*di + ijk] + g12[3*di + ijk]) + 
        66.*(g12[-4*di + ijk] + g12[4*di + ijk]) - 
        12.*(g12[-5*di + ijk] + g12[5*di + ijk]) + g12[6*di + ijk]) + 
     oo2dy*(924.*g12[ijk] + g12[-6*dj + ijk] - 
        792.*(g12[-dj + ijk] + g12[dj + ijk]) + 
        495.*(g12[-2*dj + ijk] + g12[2*dj + ijk]) - 
        220.*(g12[-3*dj + ijk] + g12[3*dj + ijk]) + 
        66.*(g12[-4*dj + ijk] + g12[4*dj + ijk]) - 
        12.*(g12[-5*dj + ijk] + g12[5*dj + ijk]) + g12[6*dj + ijk]) + 
     oo2dz*(924.*g12[ijk] + g12[-6*dk + ijk] - 
        792.*(g12[-dk + ijk] + g12[dk + ijk]) + 
        495.*(g12[-2*dk + ijk] + g12[2*dk + ijk]) - 
        220.*(g12[-3*dk + ijk] + g12[3*dk + ijk]) + 
        66.*(g12[-4*dk + ijk] + g12[4*dk + ijk]) - 
        12.*(g12[-5*dk + ijk] + g12[5*dk + ijk]) + g12[6*dk + ijk]))
;

rg13
=
rg13 - 2.*dissfactor*(oo2dx*(924.*g13[ijk] + g13[-6*di + ijk] - 
        792.*(g13[-di + ijk] + g13[di + ijk]) + 
        495.*(g13[-2*di + ijk] + g13[2*di + ijk]) - 
        220.*(g13[-3*di + ijk] + g13[3*di + ijk]) + 
        66.*(g13[-4*di + ijk] + g13[4*di + ijk]) - 
        12.*(g13[-5*di + ijk] + g13[5*di + ijk]) + g13[6*di + ijk]) + 
     oo2dy*(924.*g13[ijk] + g13[-6*dj + ijk] - 
        792.*(g13[-dj + ijk] + g13[dj + ijk]) + 
        495.*(g13[-2*dj + ijk] + g13[2*dj + ijk]) - 
        220.*(g13[-3*dj + ijk] + g13[3*dj + ijk]) + 
        66.*(g13[-4*dj + ijk] + g13[4*dj + ijk]) - 
        12.*(g13[-5*dj + ijk] + g13[5*dj + ijk]) + g13[6*dj + ijk]) + 
     oo2dz*(924.*g13[ijk] + g13[-6*dk + ijk] - 
        792.*(g13[-dk + ijk] + g13[dk + ijk]) + 
        495.*(g13[-2*dk + ijk] + g13[2*dk + ijk]) - 
        220.*(g13[-3*dk + ijk] + g13[3*dk + ijk]) + 
        66.*(g13[-4*dk + ijk] + g13[4*dk + ijk]) - 
        12.*(g13[-5*dk + ijk] + g13[5*dk + ijk]) + g13[6*dk + ijk]))
;

rg22
=
rg22 - 2.*dissfactor*(oo2dx*(924.*g22[ijk] + g22[-6*di + ijk] - 
        792.*(g22[-di + ijk] + g22[di + ijk]) + 
        495.*(g22[-2*di + ijk] + g22[2*di + ijk]) - 
        220.*(g22[-3*di + ijk] + g22[3*di + ijk]) + 
        66.*(g22[-4*di + ijk] + g22[4*di + ijk]) - 
        12.*(g22[-5*di + ijk] + g22[5*di + ijk]) + g22[6*di + ijk]) + 
     oo2dy*(924.*g22[ijk] + g22[-6*dj + ijk] - 
        792.*(g22[-dj + ijk] + g22[dj + ijk]) + 
        495.*(g22[-2*dj + ijk] + g22[2*dj + ijk]) - 
        220.*(g22[-3*dj + ijk] + g22[3*dj + ijk]) + 
        66.*(g22[-4*dj + ijk] + g22[4*dj + ijk]) - 
        12.*(g22[-5*dj + ijk] + g22[5*dj + ijk]) + g22[6*dj + ijk]) + 
     oo2dz*(924.*g22[ijk] + g22[-6*dk + ijk] - 
        792.*(g22[-dk + ijk] + g22[dk + ijk]) + 
        495.*(g22[-2*dk + ijk] + g22[2*dk + ijk]) - 
        220.*(g22[-3*dk + ijk] + g22[3*dk + ijk]) + 
        66.*(g22[-4*dk + ijk] + g22[4*dk + ijk]) - 
        12.*(g22[-5*dk + ijk] + g22[5*dk + ijk]) + g22[6*dk + ijk]))
;

rg23
=
rg23 - 2.*dissfactor*(oo2dx*(924.*g23[ijk] + g23[-6*di + ijk] - 
        792.*(g23[-di + ijk] + g23[di + ijk]) + 
        495.*(g23[-2*di + ijk] + g23[2*di + ijk]) - 
        220.*(g23[-3*di + ijk] + g23[3*di + ijk]) + 
        66.*(g23[-4*di + ijk] + g23[4*di + ijk]) - 
        12.*(g23[-5*di + ijk] + g23[5*di + ijk]) + g23[6*di + ijk]) + 
     oo2dy*(924.*g23[ijk] + g23[-6*dj + ijk] - 
        792.*(g23[-dj + ijk] + g23[dj + ijk]) + 
        495.*(g23[-2*dj + ijk] + g23[2*dj + ijk]) - 
        220.*(g23[-3*dj + ijk] + g23[3*dj + ijk]) + 
        66.*(g23[-4*dj + ijk] + g23[4*dj + ijk]) - 
        12.*(g23[-5*dj + ijk] + g23[5*dj + ijk]) + g23[6*dj + ijk]) + 
     oo2dz*(924.*g23[ijk] + g23[-6*dk + ijk] - 
        792.*(g23[-dk + ijk] + g23[dk + ijk]) + 
        495.*(g23[-2*dk + ijk] + g23[2*dk + ijk]) - 
        220.*(g23[-3*dk + ijk] + g23[3*dk + ijk]) + 
        66.*(g23[-4*dk + ijk] + g23[4*dk + ijk]) - 
        12.*(g23[-5*dk + ijk] + g23[5*dk + ijk]) + g23[6*dk + ijk]))
;

rg33
=
rg33 - 2.*dissfactor*(oo2dx*(924.*g33[ijk] + g33[-6*di + ijk] - 
        792.*(g33[-di + ijk] + g33[di + ijk]) + 
        495.*(g33[-2*di + ijk] + g33[2*di + ijk]) - 
        220.*(g33[-3*di + ijk] + g33[3*di + ijk]) + 
        66.*(g33[-4*di + ijk] + g33[4*di + ijk]) - 
        12.*(g33[-5*di + ijk] + g33[5*di + ijk]) + g33[6*di + ijk]) + 
     oo2dy*(924.*g33[ijk] + g33[-6*dj + ijk] - 
        792.*(g33[-dj + ijk] + g33[dj + ijk]) + 
        495.*(g33[-2*dj + ijk] + g33[2*dj + ijk]) - 
        220.*(g33[-3*dj + ijk] + g33[3*dj + ijk]) + 
        66.*(g33[-4*dj + ijk] + g33[4*dj + ijk]) - 
        12.*(g33[-5*dj + ijk] + g33[5*dj + ijk]) + g33[6*dj + ijk]) + 
     oo2dz*(924.*g33[ijk] + g33[-6*dk + ijk] - 
        792.*(g33[-dk + ijk] + g33[dk + ijk]) + 
        495.*(g33[-2*dk + ijk] + g33[2*dk + ijk]) - 
        220.*(g33[-3*dk + ijk] + g33[3*dk + ijk]) + 
        66.*(g33[-4*dk + ijk] + g33[4*dk + ijk]) - 
        12.*(g33[-5*dk + ijk] + g33[5*dk + ijk]) + g33[6*dk + ijk]))
;

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

rchi
=
rchi - 2.*dissfactor*(oo2dx*(924.*chi[ijk] + chi[-6*di + ijk] - 
        792.*(chi[-di + ijk] + chi[di + ijk]) + 
        495.*(chi[-2*di + ijk] + chi[2*di + ijk]) - 
        220.*(chi[-3*di + ijk] + chi[3*di + ijk]) + 
        66.*(chi[-4*di + ijk] + chi[4*di + ijk]) - 
        12.*(chi[-5*di + ijk] + chi[5*di + ijk]) + chi[6*di + ijk]) + 
     oo2dy*(924.*chi[ijk] + chi[-6*dj + ijk] - 
        792.*(chi[-dj + ijk] + chi[dj + ijk]) + 
        495.*(chi[-2*dj + ijk] + chi[2*dj + ijk]) - 
        220.*(chi[-3*dj + ijk] + chi[3*dj + ijk]) + 
        66.*(chi[-4*dj + ijk] + chi[4*dj + ijk]) - 
        12.*(chi[-5*dj + ijk] + chi[5*dj + ijk]) + chi[6*dj + ijk]) + 
     oo2dz*(924.*chi[ijk] + chi[-6*dk + ijk] - 
        792.*(chi[-dk + ijk] + chi[dk + ijk]) + 
        495.*(chi[-2*dk + ijk] + chi[2*dk + ijk]) - 
        220.*(chi[-3*dk + ijk] + chi[3*dk + ijk]) + 
        66.*(chi[-4*dk + ijk] + chi[4*dk + ijk]) - 
        12.*(chi[-5*dk + ijk] + chi[5*dk + ijk]) + chi[6*dk + ijk]))
;


#endif 


} 


if (setRHSto0) { 

rg11
=
0
;

rg12
=
0
;

rg13
=
0
;

rg22
=
0
;

rg23
=
0
;

rg33
=
0
;

rA11
=
0
;

rA12
=
0
;

rA13
=
0
;

rA22
=
0
;

rA23
=
0
;

rA33
=
0
;

rG1
=
0
;

rG2
=
0
;

rG3
=
0
;

rK
=
0
;

rchi
=
0
;

rbeta1
=
0
;

rbeta2
=
0
;

rbeta3
=
0
;

rB1
=
0
;

rB2
=
0
;

rB3
=
0
;

ralpha
=
0
;


} 


CheckForNANandINF(24, rA11,rA12,rA13,rA22,rA23,rA33, ralpha,rB1,rB2,rB3,rbeta1,rbeta2,rbeta3,rchi, rG1,rg11,rg12,rg13,rG2,rg22,rg23,rG3,g33,rK); 



/* conditional */
if (addlinear) {

ng11[ijk]
=
c*rg11 + pg11[ijk]
;

ng12[ijk]
=
c*rg12 + pg12[ijk]
;

ng13[ijk]
=
c*rg13 + pg13[ijk]
;

ng22[ijk]
=
c*rg22 + pg22[ijk]
;

ng23[ijk]
=
c*rg23 + pg23[ijk]
;

ng33[ijk]
=
c*rg33 + pg33[ijk]
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

nK[ijk]
=
forceKzerofactor*(c*rK + pK[ijk])
;

nchi[ijk]
=
c*rchi + pchi[ijk]
;

nalpha[ijk]
=
c*ralpha + palpha[ijk]
;

nbeta1[ijk]
=
c*rbeta1 + pbeta1[ijk]
;

nbeta2[ijk]
=
c*rbeta2 + pbeta2[ijk]
;

nbeta3[ijk]
=
c*rbeta3 + pbeta3[ijk]
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

detnginv
=
1/(2.*ng12[ijk]*ng13[ijk]*ng23[ijk] + ng11[ijk]*ng22[ijk]*ng33[ijk] - 
    ng33[ijk]*pow2(ng12[ijk]) - ng22[ijk]*pow2(ng13[ijk]) - 
    ng11[ijk]*pow2(ng23[ijk]))
;



/* conditional */
if (subtractA) {

traceA
=
detnginv*(-2.*nA23[ijk]*ng11[ijk]*ng23[ijk] + 
    nA11[ijk]*ng22[ijk]*ng33[ijk] + 
    ng11[ijk]*(nA33[ijk]*ng22[ijk] + nA22[ijk]*ng33[ijk]) + 
    2.*(ng13[ijk]*(nA23[ijk]*ng12[ijk] - nA13[ijk]*ng22[ijk] + 
          nA12[ijk]*ng23[ijk]) + 
       ng12[ijk]*(nA13[ijk]*ng23[ijk] - nA12[ijk]*ng33[ijk])) - 
    nA33[ijk]*pow2(ng12[ijk]) - nA22[ijk]*pow2(ng13[ijk]) - 
    nA11[ijk]*pow2(ng23[ijk]))
;

aux
=
-0.33333333333333333333*traceA
;

nA11[ijk]
=
nA11[ijk] + aux*ng11[ijk]
;

nA12[ijk]
=
nA12[ijk] + aux*ng12[ijk]
;

nA13[ijk]
=
nA13[ijk] + aux*ng13[ijk]
;

nA22[ijk]
=
nA22[ijk] + aux*ng22[ijk]
;

nA23[ijk]
=
nA23[ijk] + aux*ng23[ijk]
;

nA33[ijk]
=
nA33[ijk] + aux*ng33[ijk]
;

}
/* if (subtractA) */




/* conditional */
if (normalizedetg) {

aux
=
Power(detnginv,0.3333333333333333)
;

ng11[ijk]
=
aux*ng11[ijk]
;

ng12[ijk]
=
aux*ng12[ijk]
;

ng13[ijk]
=
aux*ng13[ijk]
;

ng22[ijk]
=
aux*ng22[ijk]
;

ng23[ijk]
=
aux*ng23[ijk]
;

ng33[ijk]
=
aux*ng33[ijk]
;

}
/* if (normalizedetg) */



} else { /* if (!normalizedetg) */

ng11[ijk]
=
rg11
;

ng12[ijk]
=
rg12
;

ng13[ijk]
=
rg13
;

ng22[ijk]
=
rg22
;

ng23[ijk]
=
rg23
;

ng33[ijk]
=
rg33
;

nchi[ijk]
=
rchi
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

nK[ijk]
=
forceKzerofactor*rK
;

nalpha[ijk]
=
ralpha
;

nbeta1[ijk]
=
rbeta1
;

nbeta2[ijk]
=
rbeta2
;

nbeta3[ijk]
=
rbeta3
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
/* if (normalizedetg) */




} endfor_ijk_openmp; /* loop i, j, k */



bampi_openmp_stop


}  /* function */

/* bssn_rhs_movpunc_N.c */
/* nvars = 78, nauxs = 478, n* = 24368,  n/ = 94,  n+ = 44133, n = 68595, O = 0 */
