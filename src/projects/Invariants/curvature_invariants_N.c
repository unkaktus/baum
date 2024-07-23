/* curvature_invariants_N.c */
/* Copyright (C) 1998 Bernd Bruegmann, 28.2.2019 */
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




void curvature_invariants_N(tVarList *u)
{

double *g11 = vldataptr(u, 0);
double *g12 = vldataptr(u, 1);
double *g13 = vldataptr(u, 2);
double *g22 = vldataptr(u, 3);
double *g23 = vldataptr(u, 4);
double *g33 = vldataptr(u, 5);
double *K11 = vldataptr(u, 6);
double *K12 = vldataptr(u, 7);
double *K13 = vldataptr(u, 8);
double *K22 = vldataptr(u, 9);
double *K23 = vldataptr(u, 10);
double *K33 = vldataptr(u, 11);
double *alpha = vldataptr(u, 12);
double *beta1 = vldataptr(u, 13);
double *beta2 = vldataptr(u, 14);
double *beta3 = vldataptr(u, 15);
double *xp = vldataptr(u, 16);
double *yp = vldataptr(u, 17);
double *zp = vldataptr(u, 18);
double *rpsi4 = vldataptr(u, 19);
double *ipsi4 = vldataptr(u, 20);
double *Csqr = vldataptr(u, 21);
tL *level = u->level;

double oosqrt2 = 1.0/sqrt(2);
double sqrt2 = sqrt(2);

const int kinnersley = Getv("invariants_tetrad","kinnersley");
const int gausscodacciflat = Getv("gauss_codacci_mainardi", "flat");
const int order_centered = Geti("order_centered");
const double x0 = Getd("sphere_x0");
const double y0 = Getd("sphere_y0");
const double z0 = Getd("sphere_z0");

const int useShellsTransfo = (Getv("grid","shells") && level->l==0);
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

double B11 = 0.;
double B12 = 0.;
double B13 = 0.;
double B21 = 0.;
double B22 = 0.;
double B23 = 0.;
double B31 = 0.;
double B32 = 0.;
double B33 = 0.;
double cosphi = 0.;
double costheta = 0.;
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
double ddgSST1111 = 0.;
double ddgSST1112 = 0.;
double ddgSST1113 = 0.;
double ddgSST1121 = 0.;
double ddgSST1122 = 0.;
double ddgSST1123 = 0.;
double ddgSST1131 = 0.;
double ddgSST1132 = 0.;
double ddgSST1133 = 0.;
double ddgSST1211 = 0.;
double ddgSST1212 = 0.;
double ddgSST1213 = 0.;
double ddgSST1221 = 0.;
double ddgSST1222 = 0.;
double ddgSST1223 = 0.;
double ddgSST1231 = 0.;
double ddgSST1232 = 0.;
double ddgSST1233 = 0.;
double ddgSST1311 = 0.;
double ddgSST1312 = 0.;
double ddgSST1313 = 0.;
double ddgSST1321 = 0.;
double ddgSST1322 = 0.;
double ddgSST1323 = 0.;
double ddgSST1331 = 0.;
double ddgSST1332 = 0.;
double ddgSST1333 = 0.;
double ddgSST2111 = 0.;
double ddgSST2112 = 0.;
double ddgSST2113 = 0.;
double ddgSST2121 = 0.;
double ddgSST2122 = 0.;
double ddgSST2123 = 0.;
double ddgSST2131 = 0.;
double ddgSST2132 = 0.;
double ddgSST2133 = 0.;
double ddgSST2211 = 0.;
double ddgSST2212 = 0.;
double ddgSST2213 = 0.;
double ddgSST2221 = 0.;
double ddgSST2222 = 0.;
double ddgSST2223 = 0.;
double ddgSST2231 = 0.;
double ddgSST2232 = 0.;
double ddgSST2233 = 0.;
double ddgSST2311 = 0.;
double ddgSST2312 = 0.;
double ddgSST2313 = 0.;
double ddgSST2321 = 0.;
double ddgSST2322 = 0.;
double ddgSST2323 = 0.;
double ddgSST2331 = 0.;
double ddgSST2332 = 0.;
double ddgSST2333 = 0.;
double ddgSST3111 = 0.;
double ddgSST3112 = 0.;
double ddgSST3113 = 0.;
double ddgSST3121 = 0.;
double ddgSST3122 = 0.;
double ddgSST3123 = 0.;
double ddgSST3131 = 0.;
double ddgSST3132 = 0.;
double ddgSST3133 = 0.;
double ddgSST3211 = 0.;
double ddgSST3212 = 0.;
double ddgSST3213 = 0.;
double ddgSST3221 = 0.;
double ddgSST3222 = 0.;
double ddgSST3223 = 0.;
double ddgSST3231 = 0.;
double ddgSST3232 = 0.;
double ddgSST3233 = 0.;
double ddgSST3311 = 0.;
double ddgSST3312 = 0.;
double ddgSST3313 = 0.;
double ddgSST3321 = 0.;
double ddgSST3322 = 0.;
double ddgSST3323 = 0.;
double ddgSST3331 = 0.;
double ddgSST3332 = 0.;
double ddgSST3333 = 0.;
double ddRdr = 0.;
double detg = 0.;
double dg111 = 0.;
double dg112 = 0.;
double dg113 = 0.;
double dg122 = 0.;
double dg123 = 0.;
double dg133 = 0.;
double dg211 = 0.;
double dg212 = 0.;
double dg213 = 0.;
double dg222 = 0.;
double dg223 = 0.;
double dg233 = 0.;
double dg311 = 0.;
double dg312 = 0.;
double dg313 = 0.;
double dg322 = 0.;
double dg323 = 0.;
double dg333 = 0.;
double dgSST111 = 0.;
double dgSST112 = 0.;
double dgSST113 = 0.;
double dgSST121 = 0.;
double dgSST122 = 0.;
double dgSST123 = 0.;
double dgSST131 = 0.;
double dgSST132 = 0.;
double dgSST133 = 0.;
double dgSST211 = 0.;
double dgSST212 = 0.;
double dgSST213 = 0.;
double dgSST221 = 0.;
double dgSST222 = 0.;
double dgSST223 = 0.;
double dgSST231 = 0.;
double dgSST232 = 0.;
double dgSST233 = 0.;
double dgSST311 = 0.;
double dgSST312 = 0.;
double dgSST313 = 0.;
double dgSST321 = 0.;
double dgSST322 = 0.;
double dgSST323 = 0.;
double dgSST331 = 0.;
double dgSST332 = 0.;
double dgSST333 = 0.;
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
double dK111 = 0.;
double dK112 = 0.;
double dK113 = 0.;
double dK122 = 0.;
double dK123 = 0.;
double dK133 = 0.;
double dK211 = 0.;
double dK212 = 0.;
double dK213 = 0.;
double dK222 = 0.;
double dK223 = 0.;
double dK233 = 0.;
double dK311 = 0.;
double dK312 = 0.;
double dK313 = 0.;
double dK322 = 0.;
double dK323 = 0.;
double dK333 = 0.;
double dKSST111 = 0.;
double dKSST112 = 0.;
double dKSST113 = 0.;
double dKSST121 = 0.;
double dKSST122 = 0.;
double dKSST123 = 0.;
double dKSST131 = 0.;
double dKSST132 = 0.;
double dKSST133 = 0.;
double dKSST211 = 0.;
double dKSST212 = 0.;
double dKSST213 = 0.;
double dKSST221 = 0.;
double dKSST222 = 0.;
double dKSST223 = 0.;
double dKSST231 = 0.;
double dKSST232 = 0.;
double dKSST233 = 0.;
double dKSST311 = 0.;
double dKSST312 = 0.;
double dKSST313 = 0.;
double dKSST321 = 0.;
double dKSST322 = 0.;
double dKSST323 = 0.;
double dKSST331 = 0.;
double dKSST332 = 0.;
double dKSST333 = 0.;
double dRdr = 0.;
double E11 = 0.;
double E12 = 0.;
double E13 = 0.;
double E21 = 0.;
double E22 = 0.;
double E23 = 0.;
double E31 = 0.;
double E32 = 0.;
double E33 = 0.;
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
double ginv11 = 0.;
double ginv12 = 0.;
double ginv13 = 0.;
double ginv21 = 0.;
double ginv22 = 0.;
double ginv23 = 0.;
double ginv31 = 0.;
double ginv32 = 0.;
double ginv33 = 0.;
double ginvdet11 = 0.;
double ginvdet12 = 0.;
double ginvdet13 = 0.;
double ginvdet22 = 0.;
double ginvdet23 = 0.;
double ginvdet33 = 0.;
double im1 = 0.;
double im2 = 0.;
double im3 = 0.;
double imb1 = 0.;
double imb2 = 0.;
double imb3 = 0.;
double Jac11 = 0.;
double Jac12 = 0.;
double Jac13 = 0.;
double Jac21 = 0.;
double Jac22 = 0.;
double Jac23 = 0.;
double Jac31 = 0.;
double Jac32 = 0.;
double Jac33 = 0.;
double n0 = 0.;
double n1 = 0.;
double n2 = 0.;
double n3 = 0.;
double r = 0.;
double R = 0.;
double R11 = 0.;
double R12 = 0.;
double R1212 = 0.;
double R1213 = 0.;
double R1223 = 0.;
double R13 = 0.;
double R1313 = 0.;
double R1323 = 0.;
double R22 = 0.;
double R23 = 0.;
double R2323 = 0.;
double R33 = 0.;
double Riemm11 = 0.;
double Riemm112 = 0.;
double Riemm113 = 0.;
double Riemm12 = 0.;
double Riemm1212 = 0.;
double Riemm1213 = 0.;
double Riemm1223 = 0.;
double Riemm123 = 0.;
double Riemm13 = 0.;
double Riemm1313 = 0.;
double Riemm1323 = 0.;
double Riemm212 = 0.;
double Riemm213 = 0.;
double Riemm22 = 0.;
double Riemm223 = 0.;
double Riemm23 = 0.;
double Riemm2323 = 0.;
double Riemm312 = 0.;
double Riemm313 = 0.;
double Riemm323 = 0.;
double Riemm33 = 0.;
double rm1 = 0.;
double rm2 = 0.;
double rm3 = 0.;
double rmb1 = 0.;
double rmb2 = 0.;
double rmb3 = 0.;
double sinphi = 0.;
double sintheta = 0.;
double Sqrtdetg = 0.;
double tfactor = 0.;
double trK = 0.;
double v11 = 0.;
double v12 = 0.;
double v13 = 0.;
double v21 = 0.;
double v22 = 0.;
double v23 = 0.;
double v31 = 0.;
double v32 = 0.;
double v33 = 0.;
double vecp1 = 0.;
double vecp2 = 0.;
double vecp3 = 0.;
double vecr1 = 0.;
double vecr2 = 0.;
double vecr3 = 0.;
double vect1 = 0.;
double vect2 = 0.;
double vect3 = 0.;
double vecu0 = 0.;
double vecu1 = 0.;
double vecu2 = 0.;
double vecu3 = 0.;
double w11 = 0.;
double w12 = 0.;
double w13 = 0.;
double w22 = 0.;
double w23 = 0.;
double w33 = 0.;
double x = 0.;
double y = 0.;
double z = 0.;



forinnerpoints_ijk_openmp(level) {



if (order_centered == 2 || boundary1away) { 

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

dK111
=
oo2dx*(-K11[-di + ijk] + K11[di + ijk])
;

dK112
=
oo2dx*(-K12[-di + ijk] + K12[di + ijk])
;

dK113
=
oo2dx*(-K13[-di + ijk] + K13[di + ijk])
;

dK122
=
oo2dx*(-K22[-di + ijk] + K22[di + ijk])
;

dK123
=
oo2dx*(-K23[-di + ijk] + K23[di + ijk])
;

dK133
=
oo2dx*(-K33[-di + ijk] + K33[di + ijk])
;

dK211
=
oo2dy*(-K11[-dj + ijk] + K11[dj + ijk])
;

dK212
=
oo2dy*(-K12[-dj + ijk] + K12[dj + ijk])
;

dK213
=
oo2dy*(-K13[-dj + ijk] + K13[dj + ijk])
;

dK222
=
oo2dy*(-K22[-dj + ijk] + K22[dj + ijk])
;

dK223
=
oo2dy*(-K23[-dj + ijk] + K23[dj + ijk])
;

dK233
=
oo2dy*(-K33[-dj + ijk] + K33[dj + ijk])
;

dK311
=
oo2dz*(-K11[-dk + ijk] + K11[dk + ijk])
;

dK312
=
oo2dz*(-K12[-dk + ijk] + K12[dk + ijk])
;

dK313
=
oo2dz*(-K13[-dk + ijk] + K13[dk + ijk])
;

dK322
=
oo2dz*(-K22[-dk + ijk] + K22[dk + ijk])
;

dK323
=
oo2dz*(-K23[-dk + ijk] + K23[dk + ijk])
;

dK333
=
oo2dz*(-K33[-dk + ijk] + K33[dk + ijk])
;


} else if (order_centered == 4 || boundaryNaway(2)) { 

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

ddg1111
=
0.083333333333333333333*oodx2*(-30.*g11[ijk] - g11[-2*di + ijk] + 
    16.*(g11[-di + ijk] + g11[di + ijk]) - g11[2*di + ijk])
;

ddg1112
=
0.083333333333333333333*oodx2*(-30.*g12[ijk] - g12[-2*di + ijk] + 
    16.*(g12[-di + ijk] + g12[di + ijk]) - g12[2*di + ijk])
;

ddg1113
=
0.083333333333333333333*oodx2*(-30.*g13[ijk] - g13[-2*di + ijk] + 
    16.*(g13[-di + ijk] + g13[di + ijk]) - g13[2*di + ijk])
;

ddg1122
=
0.083333333333333333333*oodx2*(-30.*g22[ijk] - g22[-2*di + ijk] + 
    16.*(g22[-di + ijk] + g22[di + ijk]) - g22[2*di + ijk])
;

ddg1123
=
0.083333333333333333333*oodx2*(-30.*g23[ijk] - g23[-2*di + ijk] + 
    16.*(g23[-di + ijk] + g23[di + ijk]) - g23[2*di + ijk])
;

ddg1133
=
0.083333333333333333333*oodx2*(-30.*g33[ijk] - g33[-2*di + ijk] + 
    16.*(g33[-di + ijk] + g33[di + ijk]) - g33[2*di + ijk])
;

ddg1211
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

ddg1212
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

ddg1213
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

ddg1222
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

ddg1223
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

ddg1233
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

ddg1311
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

ddg1312
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

ddg1313
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

ddg1322
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

ddg1323
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

ddg1333
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

ddg2211
=
0.083333333333333333333*oody2*(-30.*g11[ijk] - g11[-2*dj + ijk] + 
    16.*(g11[-dj + ijk] + g11[dj + ijk]) - g11[2*dj + ijk])
;

ddg2212
=
0.083333333333333333333*oody2*(-30.*g12[ijk] - g12[-2*dj + ijk] + 
    16.*(g12[-dj + ijk] + g12[dj + ijk]) - g12[2*dj + ijk])
;

ddg2213
=
0.083333333333333333333*oody2*(-30.*g13[ijk] - g13[-2*dj + ijk] + 
    16.*(g13[-dj + ijk] + g13[dj + ijk]) - g13[2*dj + ijk])
;

ddg2222
=
0.083333333333333333333*oody2*(-30.*g22[ijk] - g22[-2*dj + ijk] + 
    16.*(g22[-dj + ijk] + g22[dj + ijk]) - g22[2*dj + ijk])
;

ddg2223
=
0.083333333333333333333*oody2*(-30.*g23[ijk] - g23[-2*dj + ijk] + 
    16.*(g23[-dj + ijk] + g23[dj + ijk]) - g23[2*dj + ijk])
;

ddg2233
=
0.083333333333333333333*oody2*(-30.*g33[ijk] - g33[-2*dj + ijk] + 
    16.*(g33[-dj + ijk] + g33[dj + ijk]) - g33[2*dj + ijk])
;

ddg2311
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

ddg2312
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

ddg2313
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

ddg2322
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

ddg2323
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

ddg2333
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

ddg3311
=
0.083333333333333333333*oodz2*(-30.*g11[ijk] - g11[-2*dk + ijk] + 
    16.*(g11[-dk + ijk] + g11[dk + ijk]) - g11[2*dk + ijk])
;

ddg3312
=
0.083333333333333333333*oodz2*(-30.*g12[ijk] - g12[-2*dk + ijk] + 
    16.*(g12[-dk + ijk] + g12[dk + ijk]) - g12[2*dk + ijk])
;

ddg3313
=
0.083333333333333333333*oodz2*(-30.*g13[ijk] - g13[-2*dk + ijk] + 
    16.*(g13[-dk + ijk] + g13[dk + ijk]) - g13[2*dk + ijk])
;

ddg3322
=
0.083333333333333333333*oodz2*(-30.*g22[ijk] - g22[-2*dk + ijk] + 
    16.*(g22[-dk + ijk] + g22[dk + ijk]) - g22[2*dk + ijk])
;

ddg3323
=
0.083333333333333333333*oodz2*(-30.*g23[ijk] - g23[-2*dk + ijk] + 
    16.*(g23[-dk + ijk] + g23[dk + ijk]) - g23[2*dk + ijk])
;

ddg3333
=
0.083333333333333333333*oodz2*(-30.*g33[ijk] - g33[-2*dk + ijk] + 
    16.*(g33[-dk + ijk] + g33[dk + ijk]) - g33[2*dk + ijk])
;

dK111
=
0.16666666666666666667*oo2dx*(K11[-2*di + ijk] + 
    8.*(-K11[-di + ijk] + K11[di + ijk]) - K11[2*di + ijk])
;

dK112
=
0.16666666666666666667*oo2dx*(K12[-2*di + ijk] + 
    8.*(-K12[-di + ijk] + K12[di + ijk]) - K12[2*di + ijk])
;

dK113
=
0.16666666666666666667*oo2dx*(K13[-2*di + ijk] + 
    8.*(-K13[-di + ijk] + K13[di + ijk]) - K13[2*di + ijk])
;

dK122
=
0.16666666666666666667*oo2dx*(K22[-2*di + ijk] + 
    8.*(-K22[-di + ijk] + K22[di + ijk]) - K22[2*di + ijk])
;

dK123
=
0.16666666666666666667*oo2dx*(K23[-2*di + ijk] + 
    8.*(-K23[-di + ijk] + K23[di + ijk]) - K23[2*di + ijk])
;

dK133
=
0.16666666666666666667*oo2dx*(K33[-2*di + ijk] + 
    8.*(-K33[-di + ijk] + K33[di + ijk]) - K33[2*di + ijk])
;

dK211
=
0.16666666666666666667*oo2dy*(K11[-2*dj + ijk] + 
    8.*(-K11[-dj + ijk] + K11[dj + ijk]) - K11[2*dj + ijk])
;

dK212
=
0.16666666666666666667*oo2dy*(K12[-2*dj + ijk] + 
    8.*(-K12[-dj + ijk] + K12[dj + ijk]) - K12[2*dj + ijk])
;

dK213
=
0.16666666666666666667*oo2dy*(K13[-2*dj + ijk] + 
    8.*(-K13[-dj + ijk] + K13[dj + ijk]) - K13[2*dj + ijk])
;

dK222
=
0.16666666666666666667*oo2dy*(K22[-2*dj + ijk] + 
    8.*(-K22[-dj + ijk] + K22[dj + ijk]) - K22[2*dj + ijk])
;

dK223
=
0.16666666666666666667*oo2dy*(K23[-2*dj + ijk] + 
    8.*(-K23[-dj + ijk] + K23[dj + ijk]) - K23[2*dj + ijk])
;

dK233
=
0.16666666666666666667*oo2dy*(K33[-2*dj + ijk] + 
    8.*(-K33[-dj + ijk] + K33[dj + ijk]) - K33[2*dj + ijk])
;

dK311
=
0.16666666666666666667*oo2dz*(K11[-2*dk + ijk] + 
    8.*(-K11[-dk + ijk] + K11[dk + ijk]) - K11[2*dk + ijk])
;

dK312
=
0.16666666666666666667*oo2dz*(K12[-2*dk + ijk] + 
    8.*(-K12[-dk + ijk] + K12[dk + ijk]) - K12[2*dk + ijk])
;

dK313
=
0.16666666666666666667*oo2dz*(K13[-2*dk + ijk] + 
    8.*(-K13[-dk + ijk] + K13[dk + ijk]) - K13[2*dk + ijk])
;

dK322
=
0.16666666666666666667*oo2dz*(K22[-2*dk + ijk] + 
    8.*(-K22[-dk + ijk] + K22[dk + ijk]) - K22[2*dk + ijk])
;

dK323
=
0.16666666666666666667*oo2dz*(K23[-2*dk + ijk] + 
    8.*(-K23[-dk + ijk] + K23[dk + ijk]) - K23[2*dk + ijk])
;

dK333
=
0.16666666666666666667*oo2dz*(K33[-2*dk + ijk] + 
    8.*(-K33[-dk + ijk] + K33[dk + ijk]) - K33[2*dk + ijk])
;


} else if (order_centered == 6 || boundaryNaway(3)) { 

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

ddg1111
=
0.0055555555555555555556*oodx2*
  (-490.*g11[ijk] + 270.*(g11[-di + ijk] + g11[di + ijk]) - 
    27.*(g11[-2*di + ijk] + g11[2*di + ijk]) + 
    2.*(g11[-3*di + ijk] + g11[3*di + ijk]))
;

ddg1112
=
0.0055555555555555555556*oodx2*
  (-490.*g12[ijk] + 270.*(g12[-di + ijk] + g12[di + ijk]) - 
    27.*(g12[-2*di + ijk] + g12[2*di + ijk]) + 
    2.*(g12[-3*di + ijk] + g12[3*di + ijk]))
;

ddg1113
=
0.0055555555555555555556*oodx2*
  (-490.*g13[ijk] + 270.*(g13[-di + ijk] + g13[di + ijk]) - 
    27.*(g13[-2*di + ijk] + g13[2*di + ijk]) + 
    2.*(g13[-3*di + ijk] + g13[3*di + ijk]))
;

ddg1122
=
0.0055555555555555555556*oodx2*
  (-490.*g22[ijk] + 270.*(g22[-di + ijk] + g22[di + ijk]) - 
    27.*(g22[-2*di + ijk] + g22[2*di + ijk]) + 
    2.*(g22[-3*di + ijk] + g22[3*di + ijk]))
;

ddg1123
=
0.0055555555555555555556*oodx2*
  (-490.*g23[ijk] + 270.*(g23[-di + ijk] + g23[di + ijk]) - 
    27.*(g23[-2*di + ijk] + g23[2*di + ijk]) + 
    2.*(g23[-3*di + ijk] + g23[3*di + ijk]))
;

ddg1133
=
0.0055555555555555555556*oodx2*
  (-490.*g33[ijk] + 270.*(g33[-di + ijk] + g33[di + ijk]) - 
    27.*(g33[-2*di + ijk] + g33[2*di + ijk]) + 
    2.*(g33[-3*di + ijk] + g33[3*di + ijk]))
;

ddg1211
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

ddg1212
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

ddg1213
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

ddg1222
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

ddg1223
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

ddg1233
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

ddg1311
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

ddg1312
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

ddg1313
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

ddg1322
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

ddg1323
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

ddg1333
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

ddg2211
=
0.0055555555555555555556*oody2*
  (-490.*g11[ijk] + 270.*(g11[-dj + ijk] + g11[dj + ijk]) - 
    27.*(g11[-2*dj + ijk] + g11[2*dj + ijk]) + 
    2.*(g11[-3*dj + ijk] + g11[3*dj + ijk]))
;

ddg2212
=
0.0055555555555555555556*oody2*
  (-490.*g12[ijk] + 270.*(g12[-dj + ijk] + g12[dj + ijk]) - 
    27.*(g12[-2*dj + ijk] + g12[2*dj + ijk]) + 
    2.*(g12[-3*dj + ijk] + g12[3*dj + ijk]))
;

ddg2213
=
0.0055555555555555555556*oody2*
  (-490.*g13[ijk] + 270.*(g13[-dj + ijk] + g13[dj + ijk]) - 
    27.*(g13[-2*dj + ijk] + g13[2*dj + ijk]) + 
    2.*(g13[-3*dj + ijk] + g13[3*dj + ijk]))
;

ddg2222
=
0.0055555555555555555556*oody2*
  (-490.*g22[ijk] + 270.*(g22[-dj + ijk] + g22[dj + ijk]) - 
    27.*(g22[-2*dj + ijk] + g22[2*dj + ijk]) + 
    2.*(g22[-3*dj + ijk] + g22[3*dj + ijk]))
;

ddg2223
=
0.0055555555555555555556*oody2*
  (-490.*g23[ijk] + 270.*(g23[-dj + ijk] + g23[dj + ijk]) - 
    27.*(g23[-2*dj + ijk] + g23[2*dj + ijk]) + 
    2.*(g23[-3*dj + ijk] + g23[3*dj + ijk]))
;

ddg2233
=
0.0055555555555555555556*oody2*
  (-490.*g33[ijk] + 270.*(g33[-dj + ijk] + g33[dj + ijk]) - 
    27.*(g33[-2*dj + ijk] + g33[2*dj + ijk]) + 
    2.*(g33[-3*dj + ijk] + g33[3*dj + ijk]))
;

ddg2311
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

ddg2312
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

ddg2313
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

ddg2322
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

ddg2323
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

ddg2333
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

ddg3311
=
0.0055555555555555555556*oodz2*
  (-490.*g11[ijk] + 270.*(g11[-dk + ijk] + g11[dk + ijk]) - 
    27.*(g11[-2*dk + ijk] + g11[2*dk + ijk]) + 
    2.*(g11[-3*dk + ijk] + g11[3*dk + ijk]))
;

ddg3312
=
0.0055555555555555555556*oodz2*
  (-490.*g12[ijk] + 270.*(g12[-dk + ijk] + g12[dk + ijk]) - 
    27.*(g12[-2*dk + ijk] + g12[2*dk + ijk]) + 
    2.*(g12[-3*dk + ijk] + g12[3*dk + ijk]))
;

ddg3313
=
0.0055555555555555555556*oodz2*
  (-490.*g13[ijk] + 270.*(g13[-dk + ijk] + g13[dk + ijk]) - 
    27.*(g13[-2*dk + ijk] + g13[2*dk + ijk]) + 
    2.*(g13[-3*dk + ijk] + g13[3*dk + ijk]))
;

ddg3322
=
0.0055555555555555555556*oodz2*
  (-490.*g22[ijk] + 270.*(g22[-dk + ijk] + g22[dk + ijk]) - 
    27.*(g22[-2*dk + ijk] + g22[2*dk + ijk]) + 
    2.*(g22[-3*dk + ijk] + g22[3*dk + ijk]))
;

ddg3323
=
0.0055555555555555555556*oodz2*
  (-490.*g23[ijk] + 270.*(g23[-dk + ijk] + g23[dk + ijk]) - 
    27.*(g23[-2*dk + ijk] + g23[2*dk + ijk]) + 
    2.*(g23[-3*dk + ijk] + g23[3*dk + ijk]))
;

ddg3333
=
0.0055555555555555555556*oodz2*
  (-490.*g33[ijk] + 270.*(g33[-dk + ijk] + g33[dk + ijk]) - 
    27.*(g33[-2*dk + ijk] + g33[2*dk + ijk]) + 
    2.*(g33[-3*dk + ijk] + g33[3*dk + ijk]))
;

dK111
=
0.033333333333333333333*oo2dx*(-K11[-3*di + ijk] + 
    45.*(-K11[-di + ijk] + K11[di + ijk]) + 
    9.*(K11[-2*di + ijk] - K11[2*di + ijk]) + K11[3*di + ijk])
;

dK112
=
0.033333333333333333333*oo2dx*(-K12[-3*di + ijk] + 
    45.*(-K12[-di + ijk] + K12[di + ijk]) + 
    9.*(K12[-2*di + ijk] - K12[2*di + ijk]) + K12[3*di + ijk])
;

dK113
=
0.033333333333333333333*oo2dx*(-K13[-3*di + ijk] + 
    45.*(-K13[-di + ijk] + K13[di + ijk]) + 
    9.*(K13[-2*di + ijk] - K13[2*di + ijk]) + K13[3*di + ijk])
;

dK122
=
0.033333333333333333333*oo2dx*(-K22[-3*di + ijk] + 
    45.*(-K22[-di + ijk] + K22[di + ijk]) + 
    9.*(K22[-2*di + ijk] - K22[2*di + ijk]) + K22[3*di + ijk])
;

dK123
=
0.033333333333333333333*oo2dx*(-K23[-3*di + ijk] + 
    45.*(-K23[-di + ijk] + K23[di + ijk]) + 
    9.*(K23[-2*di + ijk] - K23[2*di + ijk]) + K23[3*di + ijk])
;

dK133
=
0.033333333333333333333*oo2dx*(-K33[-3*di + ijk] + 
    45.*(-K33[-di + ijk] + K33[di + ijk]) + 
    9.*(K33[-2*di + ijk] - K33[2*di + ijk]) + K33[3*di + ijk])
;

dK211
=
0.033333333333333333333*oo2dy*(-K11[-3*dj + ijk] + 
    45.*(-K11[-dj + ijk] + K11[dj + ijk]) + 
    9.*(K11[-2*dj + ijk] - K11[2*dj + ijk]) + K11[3*dj + ijk])
;

dK212
=
0.033333333333333333333*oo2dy*(-K12[-3*dj + ijk] + 
    45.*(-K12[-dj + ijk] + K12[dj + ijk]) + 
    9.*(K12[-2*dj + ijk] - K12[2*dj + ijk]) + K12[3*dj + ijk])
;

dK213
=
0.033333333333333333333*oo2dy*(-K13[-3*dj + ijk] + 
    45.*(-K13[-dj + ijk] + K13[dj + ijk]) + 
    9.*(K13[-2*dj + ijk] - K13[2*dj + ijk]) + K13[3*dj + ijk])
;

dK222
=
0.033333333333333333333*oo2dy*(-K22[-3*dj + ijk] + 
    45.*(-K22[-dj + ijk] + K22[dj + ijk]) + 
    9.*(K22[-2*dj + ijk] - K22[2*dj + ijk]) + K22[3*dj + ijk])
;

dK223
=
0.033333333333333333333*oo2dy*(-K23[-3*dj + ijk] + 
    45.*(-K23[-dj + ijk] + K23[dj + ijk]) + 
    9.*(K23[-2*dj + ijk] - K23[2*dj + ijk]) + K23[3*dj + ijk])
;

dK233
=
0.033333333333333333333*oo2dy*(-K33[-3*dj + ijk] + 
    45.*(-K33[-dj + ijk] + K33[dj + ijk]) + 
    9.*(K33[-2*dj + ijk] - K33[2*dj + ijk]) + K33[3*dj + ijk])
;

dK311
=
0.033333333333333333333*oo2dz*(-K11[-3*dk + ijk] + 
    45.*(-K11[-dk + ijk] + K11[dk + ijk]) + 
    9.*(K11[-2*dk + ijk] - K11[2*dk + ijk]) + K11[3*dk + ijk])
;

dK312
=
0.033333333333333333333*oo2dz*(-K12[-3*dk + ijk] + 
    45.*(-K12[-dk + ijk] + K12[dk + ijk]) + 
    9.*(K12[-2*dk + ijk] - K12[2*dk + ijk]) + K12[3*dk + ijk])
;

dK313
=
0.033333333333333333333*oo2dz*(-K13[-3*dk + ijk] + 
    45.*(-K13[-dk + ijk] + K13[dk + ijk]) + 
    9.*(K13[-2*dk + ijk] - K13[2*dk + ijk]) + K13[3*dk + ijk])
;

dK322
=
0.033333333333333333333*oo2dz*(-K22[-3*dk + ijk] + 
    45.*(-K22[-dk + ijk] + K22[dk + ijk]) + 
    9.*(K22[-2*dk + ijk] - K22[2*dk + ijk]) + K22[3*dk + ijk])
;

dK323
=
0.033333333333333333333*oo2dz*(-K23[-3*dk + ijk] + 
    45.*(-K23[-dk + ijk] + K23[dk + ijk]) + 
    9.*(K23[-2*dk + ijk] - K23[2*dk + ijk]) + K23[3*dk + ijk])
;

dK333
=
0.033333333333333333333*oo2dz*(-K33[-3*dk + ijk] + 
    45.*(-K33[-dk + ijk] + K33[dk + ijk]) + 
    9.*(K33[-2*dk + ijk] - K33[2*dk + ijk]) + K33[3*dk + ijk])
;


} else errorexit("order is not implemented yet"); 


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

dgSST121
=
dg112*Jac11 + dg212*Jac21 + dg312*Jac31
;

dgSST122
=
dg122*Jac11 + dg222*Jac21 + dg322*Jac31
;

dgSST123
=
dg123*Jac11 + dg223*Jac21 + dg323*Jac31
;

dgSST131
=
dg113*Jac11 + dg213*Jac21 + dg313*Jac31
;

dgSST132
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

dgSST221
=
dg112*Jac12 + dg212*Jac22 + dg312*Jac32
;

dgSST222
=
dg122*Jac12 + dg222*Jac22 + dg322*Jac32
;

dgSST223
=
dg123*Jac12 + dg223*Jac22 + dg323*Jac32
;

dgSST231
=
dg113*Jac12 + dg213*Jac22 + dg313*Jac32
;

dgSST232
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

dgSST321
=
dg112*Jac13 + dg212*Jac23 + dg312*Jac33
;

dgSST322
=
dg122*Jac13 + dg222*Jac23 + dg322*Jac33
;

dgSST323
=
dg123*Jac13 + dg223*Jac23 + dg323*Jac33
;

dgSST331
=
dg113*Jac13 + dg213*Jac23 + dg313*Jac33
;

dgSST332
=
dg123*Jac13 + dg223*Jac23 + dg323*Jac33
;

dgSST333
=
dg133*Jac13 + dg233*Jac23 + dg333*Jac33
;

ddgSST1111
=
dg111*DJac111 + dg211*DJac211 + dg311*DJac311 + 
  2.*(ddg2311*Jac21*Jac31 + Jac11*(ddg1211*Jac21 + ddg1311*Jac31)) + 
  ddg1111*pow2(Jac11) + ddg2211*pow2(Jac21) + ddg3311*pow2(Jac31)
;

ddgSST1112
=
dg112*DJac111 + dg212*DJac211 + dg312*DJac311 + 
  2.*(ddg2312*Jac21*Jac31 + Jac11*(ddg1212*Jac21 + ddg1312*Jac31)) + 
  ddg1112*pow2(Jac11) + ddg2212*pow2(Jac21) + ddg3312*pow2(Jac31)
;

ddgSST1113
=
dg113*DJac111 + dg213*DJac211 + dg313*DJac311 + 
  2.*(ddg2313*Jac21*Jac31 + Jac11*(ddg1213*Jac21 + ddg1313*Jac31)) + 
  ddg1113*pow2(Jac11) + ddg2213*pow2(Jac21) + ddg3313*pow2(Jac31)
;

ddgSST1121
=
dg112*DJac111 + dg212*DJac211 + dg312*DJac311 + 
  2.*(ddg2312*Jac21*Jac31 + Jac11*(ddg1212*Jac21 + ddg1312*Jac31)) + 
  ddg1112*pow2(Jac11) + ddg2212*pow2(Jac21) + ddg3312*pow2(Jac31)
;

ddgSST1122
=
dg122*DJac111 + dg222*DJac211 + dg322*DJac311 + 
  2.*(ddg2322*Jac21*Jac31 + Jac11*(ddg1222*Jac21 + ddg1322*Jac31)) + 
  ddg1122*pow2(Jac11) + ddg2222*pow2(Jac21) + ddg3322*pow2(Jac31)
;

ddgSST1123
=
dg123*DJac111 + dg223*DJac211 + dg323*DJac311 + 
  2.*(ddg2323*Jac21*Jac31 + Jac11*(ddg1223*Jac21 + ddg1323*Jac31)) + 
  ddg1123*pow2(Jac11) + ddg2223*pow2(Jac21) + ddg3323*pow2(Jac31)
;

ddgSST1131
=
dg113*DJac111 + dg213*DJac211 + dg313*DJac311 + 
  2.*(ddg2313*Jac21*Jac31 + Jac11*(ddg1213*Jac21 + ddg1313*Jac31)) + 
  ddg1113*pow2(Jac11) + ddg2213*pow2(Jac21) + ddg3313*pow2(Jac31)
;

ddgSST1132
=
dg123*DJac111 + dg223*DJac211 + dg323*DJac311 + 
  2.*(ddg2323*Jac21*Jac31 + Jac11*(ddg1223*Jac21 + ddg1323*Jac31)) + 
  ddg1123*pow2(Jac11) + ddg2223*pow2(Jac21) + ddg3323*pow2(Jac31)
;

ddgSST1133
=
dg133*DJac111 + dg233*DJac211 + dg333*DJac311 + 
  2.*(ddg2333*Jac21*Jac31 + Jac11*(ddg1233*Jac21 + ddg1333*Jac31)) + 
  ddg1133*pow2(Jac11) + ddg2233*pow2(Jac21) + ddg3333*pow2(Jac31)
;

ddgSST1211
=
dg111*DJac112 + dg211*DJac212 + dg311*DJac312 + 
  Jac12*(ddg1111*Jac11 + ddg1211*Jac21 + ddg1311*Jac31) + 
  Jac22*(ddg1211*Jac11 + ddg2211*Jac21 + ddg2311*Jac31) + 
  (ddg1311*Jac11 + ddg2311*Jac21 + ddg3311*Jac31)*Jac32
;

ddgSST1212
=
dg112*DJac112 + dg212*DJac212 + dg312*DJac312 + 
  Jac12*(ddg1112*Jac11 + ddg1212*Jac21 + ddg1312*Jac31) + 
  Jac22*(ddg1212*Jac11 + ddg2212*Jac21 + ddg2312*Jac31) + 
  (ddg1312*Jac11 + ddg2312*Jac21 + ddg3312*Jac31)*Jac32
;

ddgSST1213
=
dg113*DJac112 + dg213*DJac212 + dg313*DJac312 + 
  Jac12*(ddg1113*Jac11 + ddg1213*Jac21 + ddg1313*Jac31) + 
  Jac22*(ddg1213*Jac11 + ddg2213*Jac21 + ddg2313*Jac31) + 
  (ddg1313*Jac11 + ddg2313*Jac21 + ddg3313*Jac31)*Jac32
;

ddgSST1221
=
dg112*DJac112 + dg212*DJac212 + dg312*DJac312 + 
  Jac12*(ddg1112*Jac11 + ddg1212*Jac21 + ddg1312*Jac31) + 
  Jac22*(ddg1212*Jac11 + ddg2212*Jac21 + ddg2312*Jac31) + 
  (ddg1312*Jac11 + ddg2312*Jac21 + ddg3312*Jac31)*Jac32
;

ddgSST1222
=
dg122*DJac112 + dg222*DJac212 + dg322*DJac312 + 
  Jac12*(ddg1122*Jac11 + ddg1222*Jac21 + ddg1322*Jac31) + 
  Jac22*(ddg1222*Jac11 + ddg2222*Jac21 + ddg2322*Jac31) + 
  (ddg1322*Jac11 + ddg2322*Jac21 + ddg3322*Jac31)*Jac32
;

ddgSST1223
=
dg123*DJac112 + dg223*DJac212 + dg323*DJac312 + 
  Jac12*(ddg1123*Jac11 + ddg1223*Jac21 + ddg1323*Jac31) + 
  Jac22*(ddg1223*Jac11 + ddg2223*Jac21 + ddg2323*Jac31) + 
  (ddg1323*Jac11 + ddg2323*Jac21 + ddg3323*Jac31)*Jac32
;

ddgSST1231
=
dg113*DJac112 + dg213*DJac212 + dg313*DJac312 + 
  Jac12*(ddg1113*Jac11 + ddg1213*Jac21 + ddg1313*Jac31) + 
  Jac22*(ddg1213*Jac11 + ddg2213*Jac21 + ddg2313*Jac31) + 
  (ddg1313*Jac11 + ddg2313*Jac21 + ddg3313*Jac31)*Jac32
;

ddgSST1232
=
dg123*DJac112 + dg223*DJac212 + dg323*DJac312 + 
  Jac12*(ddg1123*Jac11 + ddg1223*Jac21 + ddg1323*Jac31) + 
  Jac22*(ddg1223*Jac11 + ddg2223*Jac21 + ddg2323*Jac31) + 
  (ddg1323*Jac11 + ddg2323*Jac21 + ddg3323*Jac31)*Jac32
;

ddgSST1233
=
dg133*DJac112 + dg233*DJac212 + dg333*DJac312 + 
  Jac12*(ddg1133*Jac11 + ddg1233*Jac21 + ddg1333*Jac31) + 
  Jac22*(ddg1233*Jac11 + ddg2233*Jac21 + ddg2333*Jac31) + 
  (ddg1333*Jac11 + ddg2333*Jac21 + ddg3333*Jac31)*Jac32
;

ddgSST1311
=
dg111*DJac113 + dg211*DJac213 + dg311*DJac313 + 
  Jac13*(ddg1111*Jac11 + ddg1211*Jac21 + ddg1311*Jac31) + 
  Jac23*(ddg1211*Jac11 + ddg2211*Jac21 + ddg2311*Jac31) + 
  (ddg1311*Jac11 + ddg2311*Jac21 + ddg3311*Jac31)*Jac33
;

ddgSST1312
=
dg112*DJac113 + dg212*DJac213 + dg312*DJac313 + 
  Jac13*(ddg1112*Jac11 + ddg1212*Jac21 + ddg1312*Jac31) + 
  Jac23*(ddg1212*Jac11 + ddg2212*Jac21 + ddg2312*Jac31) + 
  (ddg1312*Jac11 + ddg2312*Jac21 + ddg3312*Jac31)*Jac33
;

ddgSST1313
=
dg113*DJac113 + dg213*DJac213 + dg313*DJac313 + 
  Jac13*(ddg1113*Jac11 + ddg1213*Jac21 + ddg1313*Jac31) + 
  Jac23*(ddg1213*Jac11 + ddg2213*Jac21 + ddg2313*Jac31) + 
  (ddg1313*Jac11 + ddg2313*Jac21 + ddg3313*Jac31)*Jac33
;

ddgSST1321
=
dg112*DJac113 + dg212*DJac213 + dg312*DJac313 + 
  Jac13*(ddg1112*Jac11 + ddg1212*Jac21 + ddg1312*Jac31) + 
  Jac23*(ddg1212*Jac11 + ddg2212*Jac21 + ddg2312*Jac31) + 
  (ddg1312*Jac11 + ddg2312*Jac21 + ddg3312*Jac31)*Jac33
;

ddgSST1322
=
dg122*DJac113 + dg222*DJac213 + dg322*DJac313 + 
  Jac13*(ddg1122*Jac11 + ddg1222*Jac21 + ddg1322*Jac31) + 
  Jac23*(ddg1222*Jac11 + ddg2222*Jac21 + ddg2322*Jac31) + 
  (ddg1322*Jac11 + ddg2322*Jac21 + ddg3322*Jac31)*Jac33
;

ddgSST1323
=
dg123*DJac113 + dg223*DJac213 + dg323*DJac313 + 
  Jac13*(ddg1123*Jac11 + ddg1223*Jac21 + ddg1323*Jac31) + 
  Jac23*(ddg1223*Jac11 + ddg2223*Jac21 + ddg2323*Jac31) + 
  (ddg1323*Jac11 + ddg2323*Jac21 + ddg3323*Jac31)*Jac33
;

ddgSST1331
=
dg113*DJac113 + dg213*DJac213 + dg313*DJac313 + 
  Jac13*(ddg1113*Jac11 + ddg1213*Jac21 + ddg1313*Jac31) + 
  Jac23*(ddg1213*Jac11 + ddg2213*Jac21 + ddg2313*Jac31) + 
  (ddg1313*Jac11 + ddg2313*Jac21 + ddg3313*Jac31)*Jac33
;

ddgSST1332
=
dg123*DJac113 + dg223*DJac213 + dg323*DJac313 + 
  Jac13*(ddg1123*Jac11 + ddg1223*Jac21 + ddg1323*Jac31) + 
  Jac23*(ddg1223*Jac11 + ddg2223*Jac21 + ddg2323*Jac31) + 
  (ddg1323*Jac11 + ddg2323*Jac21 + ddg3323*Jac31)*Jac33
;

ddgSST1333
=
dg133*DJac113 + dg233*DJac213 + dg333*DJac313 + 
  Jac13*(ddg1133*Jac11 + ddg1233*Jac21 + ddg1333*Jac31) + 
  Jac23*(ddg1233*Jac11 + ddg2233*Jac21 + ddg2333*Jac31) + 
  (ddg1333*Jac11 + ddg2333*Jac21 + ddg3333*Jac31)*Jac33
;

ddgSST2111
=
dg111*DJac112 + dg211*DJac212 + dg311*DJac312 + 
  Jac12*(ddg1111*Jac11 + ddg1211*Jac21 + ddg1311*Jac31) + 
  Jac22*(ddg1211*Jac11 + ddg2211*Jac21 + ddg2311*Jac31) + 
  (ddg1311*Jac11 + ddg2311*Jac21 + ddg3311*Jac31)*Jac32
;

ddgSST2112
=
dg112*DJac112 + dg212*DJac212 + dg312*DJac312 + 
  Jac12*(ddg1112*Jac11 + ddg1212*Jac21 + ddg1312*Jac31) + 
  Jac22*(ddg1212*Jac11 + ddg2212*Jac21 + ddg2312*Jac31) + 
  (ddg1312*Jac11 + ddg2312*Jac21 + ddg3312*Jac31)*Jac32
;

ddgSST2113
=
dg113*DJac112 + dg213*DJac212 + dg313*DJac312 + 
  Jac12*(ddg1113*Jac11 + ddg1213*Jac21 + ddg1313*Jac31) + 
  Jac22*(ddg1213*Jac11 + ddg2213*Jac21 + ddg2313*Jac31) + 
  (ddg1313*Jac11 + ddg2313*Jac21 + ddg3313*Jac31)*Jac32
;

ddgSST2121
=
dg112*DJac112 + dg212*DJac212 + dg312*DJac312 + 
  Jac12*(ddg1112*Jac11 + ddg1212*Jac21 + ddg1312*Jac31) + 
  Jac22*(ddg1212*Jac11 + ddg2212*Jac21 + ddg2312*Jac31) + 
  (ddg1312*Jac11 + ddg2312*Jac21 + ddg3312*Jac31)*Jac32
;

ddgSST2122
=
dg122*DJac112 + dg222*DJac212 + dg322*DJac312 + 
  Jac12*(ddg1122*Jac11 + ddg1222*Jac21 + ddg1322*Jac31) + 
  Jac22*(ddg1222*Jac11 + ddg2222*Jac21 + ddg2322*Jac31) + 
  (ddg1322*Jac11 + ddg2322*Jac21 + ddg3322*Jac31)*Jac32
;

ddgSST2123
=
dg123*DJac112 + dg223*DJac212 + dg323*DJac312 + 
  Jac12*(ddg1123*Jac11 + ddg1223*Jac21 + ddg1323*Jac31) + 
  Jac22*(ddg1223*Jac11 + ddg2223*Jac21 + ddg2323*Jac31) + 
  (ddg1323*Jac11 + ddg2323*Jac21 + ddg3323*Jac31)*Jac32
;

ddgSST2131
=
dg113*DJac112 + dg213*DJac212 + dg313*DJac312 + 
  Jac12*(ddg1113*Jac11 + ddg1213*Jac21 + ddg1313*Jac31) + 
  Jac22*(ddg1213*Jac11 + ddg2213*Jac21 + ddg2313*Jac31) + 
  (ddg1313*Jac11 + ddg2313*Jac21 + ddg3313*Jac31)*Jac32
;

ddgSST2132
=
dg123*DJac112 + dg223*DJac212 + dg323*DJac312 + 
  Jac12*(ddg1123*Jac11 + ddg1223*Jac21 + ddg1323*Jac31) + 
  Jac22*(ddg1223*Jac11 + ddg2223*Jac21 + ddg2323*Jac31) + 
  (ddg1323*Jac11 + ddg2323*Jac21 + ddg3323*Jac31)*Jac32
;

ddgSST2133
=
dg133*DJac112 + dg233*DJac212 + dg333*DJac312 + 
  Jac12*(ddg1133*Jac11 + ddg1233*Jac21 + ddg1333*Jac31) + 
  Jac22*(ddg1233*Jac11 + ddg2233*Jac21 + ddg2333*Jac31) + 
  (ddg1333*Jac11 + ddg2333*Jac21 + ddg3333*Jac31)*Jac32
;

ddgSST2211
=
dg111*DJac122 + dg211*DJac222 + dg311*DJac322 + 
  2.*(ddg2311*Jac22*Jac32 + Jac12*(ddg1211*Jac22 + ddg1311*Jac32)) + 
  ddg1111*pow2(Jac12) + ddg2211*pow2(Jac22) + ddg3311*pow2(Jac32)
;

ddgSST2212
=
dg112*DJac122 + dg212*DJac222 + dg312*DJac322 + 
  2.*(ddg2312*Jac22*Jac32 + Jac12*(ddg1212*Jac22 + ddg1312*Jac32)) + 
  ddg1112*pow2(Jac12) + ddg2212*pow2(Jac22) + ddg3312*pow2(Jac32)
;

ddgSST2213
=
dg113*DJac122 + dg213*DJac222 + dg313*DJac322 + 
  2.*(ddg2313*Jac22*Jac32 + Jac12*(ddg1213*Jac22 + ddg1313*Jac32)) + 
  ddg1113*pow2(Jac12) + ddg2213*pow2(Jac22) + ddg3313*pow2(Jac32)
;

ddgSST2221
=
dg112*DJac122 + dg212*DJac222 + dg312*DJac322 + 
  2.*(ddg2312*Jac22*Jac32 + Jac12*(ddg1212*Jac22 + ddg1312*Jac32)) + 
  ddg1112*pow2(Jac12) + ddg2212*pow2(Jac22) + ddg3312*pow2(Jac32)
;

ddgSST2222
=
dg122*DJac122 + dg222*DJac222 + dg322*DJac322 + 
  2.*(ddg2322*Jac22*Jac32 + Jac12*(ddg1222*Jac22 + ddg1322*Jac32)) + 
  ddg1122*pow2(Jac12) + ddg2222*pow2(Jac22) + ddg3322*pow2(Jac32)
;

ddgSST2223
=
dg123*DJac122 + dg223*DJac222 + dg323*DJac322 + 
  2.*(ddg2323*Jac22*Jac32 + Jac12*(ddg1223*Jac22 + ddg1323*Jac32)) + 
  ddg1123*pow2(Jac12) + ddg2223*pow2(Jac22) + ddg3323*pow2(Jac32)
;

ddgSST2231
=
dg113*DJac122 + dg213*DJac222 + dg313*DJac322 + 
  2.*(ddg2313*Jac22*Jac32 + Jac12*(ddg1213*Jac22 + ddg1313*Jac32)) + 
  ddg1113*pow2(Jac12) + ddg2213*pow2(Jac22) + ddg3313*pow2(Jac32)
;

ddgSST2232
=
dg123*DJac122 + dg223*DJac222 + dg323*DJac322 + 
  2.*(ddg2323*Jac22*Jac32 + Jac12*(ddg1223*Jac22 + ddg1323*Jac32)) + 
  ddg1123*pow2(Jac12) + ddg2223*pow2(Jac22) + ddg3323*pow2(Jac32)
;

ddgSST2233
=
dg133*DJac122 + dg233*DJac222 + dg333*DJac322 + 
  2.*(ddg2333*Jac22*Jac32 + Jac12*(ddg1233*Jac22 + ddg1333*Jac32)) + 
  ddg1133*pow2(Jac12) + ddg2233*pow2(Jac22) + ddg3333*pow2(Jac32)
;

ddgSST2311
=
dg111*DJac123 + dg211*DJac223 + dg311*DJac323 + 
  Jac13*(ddg1111*Jac12 + ddg1211*Jac22 + ddg1311*Jac32) + 
  Jac23*(ddg1211*Jac12 + ddg2211*Jac22 + ddg2311*Jac32) + 
  (ddg1311*Jac12 + ddg2311*Jac22 + ddg3311*Jac32)*Jac33
;

ddgSST2312
=
dg112*DJac123 + dg212*DJac223 + dg312*DJac323 + 
  Jac13*(ddg1112*Jac12 + ddg1212*Jac22 + ddg1312*Jac32) + 
  Jac23*(ddg1212*Jac12 + ddg2212*Jac22 + ddg2312*Jac32) + 
  (ddg1312*Jac12 + ddg2312*Jac22 + ddg3312*Jac32)*Jac33
;

ddgSST2313
=
dg113*DJac123 + dg213*DJac223 + dg313*DJac323 + 
  Jac13*(ddg1113*Jac12 + ddg1213*Jac22 + ddg1313*Jac32) + 
  Jac23*(ddg1213*Jac12 + ddg2213*Jac22 + ddg2313*Jac32) + 
  (ddg1313*Jac12 + ddg2313*Jac22 + ddg3313*Jac32)*Jac33
;

ddgSST2321
=
dg112*DJac123 + dg212*DJac223 + dg312*DJac323 + 
  Jac13*(ddg1112*Jac12 + ddg1212*Jac22 + ddg1312*Jac32) + 
  Jac23*(ddg1212*Jac12 + ddg2212*Jac22 + ddg2312*Jac32) + 
  (ddg1312*Jac12 + ddg2312*Jac22 + ddg3312*Jac32)*Jac33
;

ddgSST2322
=
dg122*DJac123 + dg222*DJac223 + dg322*DJac323 + 
  Jac13*(ddg1122*Jac12 + ddg1222*Jac22 + ddg1322*Jac32) + 
  Jac23*(ddg1222*Jac12 + ddg2222*Jac22 + ddg2322*Jac32) + 
  (ddg1322*Jac12 + ddg2322*Jac22 + ddg3322*Jac32)*Jac33
;

ddgSST2323
=
dg123*DJac123 + dg223*DJac223 + dg323*DJac323 + 
  Jac13*(ddg1123*Jac12 + ddg1223*Jac22 + ddg1323*Jac32) + 
  Jac23*(ddg1223*Jac12 + ddg2223*Jac22 + ddg2323*Jac32) + 
  (ddg1323*Jac12 + ddg2323*Jac22 + ddg3323*Jac32)*Jac33
;

ddgSST2331
=
dg113*DJac123 + dg213*DJac223 + dg313*DJac323 + 
  Jac13*(ddg1113*Jac12 + ddg1213*Jac22 + ddg1313*Jac32) + 
  Jac23*(ddg1213*Jac12 + ddg2213*Jac22 + ddg2313*Jac32) + 
  (ddg1313*Jac12 + ddg2313*Jac22 + ddg3313*Jac32)*Jac33
;

ddgSST2332
=
dg123*DJac123 + dg223*DJac223 + dg323*DJac323 + 
  Jac13*(ddg1123*Jac12 + ddg1223*Jac22 + ddg1323*Jac32) + 
  Jac23*(ddg1223*Jac12 + ddg2223*Jac22 + ddg2323*Jac32) + 
  (ddg1323*Jac12 + ddg2323*Jac22 + ddg3323*Jac32)*Jac33
;

ddgSST2333
=
dg133*DJac123 + dg233*DJac223 + dg333*DJac323 + 
  Jac13*(ddg1133*Jac12 + ddg1233*Jac22 + ddg1333*Jac32) + 
  Jac23*(ddg1233*Jac12 + ddg2233*Jac22 + ddg2333*Jac32) + 
  (ddg1333*Jac12 + ddg2333*Jac22 + ddg3333*Jac32)*Jac33
;

ddgSST3111
=
dg111*DJac113 + dg211*DJac213 + dg311*DJac313 + 
  Jac13*(ddg1111*Jac11 + ddg1211*Jac21 + ddg1311*Jac31) + 
  Jac23*(ddg1211*Jac11 + ddg2211*Jac21 + ddg2311*Jac31) + 
  (ddg1311*Jac11 + ddg2311*Jac21 + ddg3311*Jac31)*Jac33
;

ddgSST3112
=
dg112*DJac113 + dg212*DJac213 + dg312*DJac313 + 
  Jac13*(ddg1112*Jac11 + ddg1212*Jac21 + ddg1312*Jac31) + 
  Jac23*(ddg1212*Jac11 + ddg2212*Jac21 + ddg2312*Jac31) + 
  (ddg1312*Jac11 + ddg2312*Jac21 + ddg3312*Jac31)*Jac33
;

ddgSST3113
=
dg113*DJac113 + dg213*DJac213 + dg313*DJac313 + 
  Jac13*(ddg1113*Jac11 + ddg1213*Jac21 + ddg1313*Jac31) + 
  Jac23*(ddg1213*Jac11 + ddg2213*Jac21 + ddg2313*Jac31) + 
  (ddg1313*Jac11 + ddg2313*Jac21 + ddg3313*Jac31)*Jac33
;

ddgSST3121
=
dg112*DJac113 + dg212*DJac213 + dg312*DJac313 + 
  Jac13*(ddg1112*Jac11 + ddg1212*Jac21 + ddg1312*Jac31) + 
  Jac23*(ddg1212*Jac11 + ddg2212*Jac21 + ddg2312*Jac31) + 
  (ddg1312*Jac11 + ddg2312*Jac21 + ddg3312*Jac31)*Jac33
;

ddgSST3122
=
dg122*DJac113 + dg222*DJac213 + dg322*DJac313 + 
  Jac13*(ddg1122*Jac11 + ddg1222*Jac21 + ddg1322*Jac31) + 
  Jac23*(ddg1222*Jac11 + ddg2222*Jac21 + ddg2322*Jac31) + 
  (ddg1322*Jac11 + ddg2322*Jac21 + ddg3322*Jac31)*Jac33
;

ddgSST3123
=
dg123*DJac113 + dg223*DJac213 + dg323*DJac313 + 
  Jac13*(ddg1123*Jac11 + ddg1223*Jac21 + ddg1323*Jac31) + 
  Jac23*(ddg1223*Jac11 + ddg2223*Jac21 + ddg2323*Jac31) + 
  (ddg1323*Jac11 + ddg2323*Jac21 + ddg3323*Jac31)*Jac33
;

ddgSST3131
=
dg113*DJac113 + dg213*DJac213 + dg313*DJac313 + 
  Jac13*(ddg1113*Jac11 + ddg1213*Jac21 + ddg1313*Jac31) + 
  Jac23*(ddg1213*Jac11 + ddg2213*Jac21 + ddg2313*Jac31) + 
  (ddg1313*Jac11 + ddg2313*Jac21 + ddg3313*Jac31)*Jac33
;

ddgSST3132
=
dg123*DJac113 + dg223*DJac213 + dg323*DJac313 + 
  Jac13*(ddg1123*Jac11 + ddg1223*Jac21 + ddg1323*Jac31) + 
  Jac23*(ddg1223*Jac11 + ddg2223*Jac21 + ddg2323*Jac31) + 
  (ddg1323*Jac11 + ddg2323*Jac21 + ddg3323*Jac31)*Jac33
;

ddgSST3133
=
dg133*DJac113 + dg233*DJac213 + dg333*DJac313 + 
  Jac13*(ddg1133*Jac11 + ddg1233*Jac21 + ddg1333*Jac31) + 
  Jac23*(ddg1233*Jac11 + ddg2233*Jac21 + ddg2333*Jac31) + 
  (ddg1333*Jac11 + ddg2333*Jac21 + ddg3333*Jac31)*Jac33
;

ddgSST3211
=
dg111*DJac123 + dg211*DJac223 + dg311*DJac323 + 
  Jac13*(ddg1111*Jac12 + ddg1211*Jac22 + ddg1311*Jac32) + 
  Jac23*(ddg1211*Jac12 + ddg2211*Jac22 + ddg2311*Jac32) + 
  (ddg1311*Jac12 + ddg2311*Jac22 + ddg3311*Jac32)*Jac33
;

ddgSST3212
=
dg112*DJac123 + dg212*DJac223 + dg312*DJac323 + 
  Jac13*(ddg1112*Jac12 + ddg1212*Jac22 + ddg1312*Jac32) + 
  Jac23*(ddg1212*Jac12 + ddg2212*Jac22 + ddg2312*Jac32) + 
  (ddg1312*Jac12 + ddg2312*Jac22 + ddg3312*Jac32)*Jac33
;

ddgSST3213
=
dg113*DJac123 + dg213*DJac223 + dg313*DJac323 + 
  Jac13*(ddg1113*Jac12 + ddg1213*Jac22 + ddg1313*Jac32) + 
  Jac23*(ddg1213*Jac12 + ddg2213*Jac22 + ddg2313*Jac32) + 
  (ddg1313*Jac12 + ddg2313*Jac22 + ddg3313*Jac32)*Jac33
;

ddgSST3221
=
dg112*DJac123 + dg212*DJac223 + dg312*DJac323 + 
  Jac13*(ddg1112*Jac12 + ddg1212*Jac22 + ddg1312*Jac32) + 
  Jac23*(ddg1212*Jac12 + ddg2212*Jac22 + ddg2312*Jac32) + 
  (ddg1312*Jac12 + ddg2312*Jac22 + ddg3312*Jac32)*Jac33
;

ddgSST3222
=
dg122*DJac123 + dg222*DJac223 + dg322*DJac323 + 
  Jac13*(ddg1122*Jac12 + ddg1222*Jac22 + ddg1322*Jac32) + 
  Jac23*(ddg1222*Jac12 + ddg2222*Jac22 + ddg2322*Jac32) + 
  (ddg1322*Jac12 + ddg2322*Jac22 + ddg3322*Jac32)*Jac33
;

ddgSST3223
=
dg123*DJac123 + dg223*DJac223 + dg323*DJac323 + 
  Jac13*(ddg1123*Jac12 + ddg1223*Jac22 + ddg1323*Jac32) + 
  Jac23*(ddg1223*Jac12 + ddg2223*Jac22 + ddg2323*Jac32) + 
  (ddg1323*Jac12 + ddg2323*Jac22 + ddg3323*Jac32)*Jac33
;

ddgSST3231
=
dg113*DJac123 + dg213*DJac223 + dg313*DJac323 + 
  Jac13*(ddg1113*Jac12 + ddg1213*Jac22 + ddg1313*Jac32) + 
  Jac23*(ddg1213*Jac12 + ddg2213*Jac22 + ddg2313*Jac32) + 
  (ddg1313*Jac12 + ddg2313*Jac22 + ddg3313*Jac32)*Jac33
;

ddgSST3232
=
dg123*DJac123 + dg223*DJac223 + dg323*DJac323 + 
  Jac13*(ddg1123*Jac12 + ddg1223*Jac22 + ddg1323*Jac32) + 
  Jac23*(ddg1223*Jac12 + ddg2223*Jac22 + ddg2323*Jac32) + 
  (ddg1323*Jac12 + ddg2323*Jac22 + ddg3323*Jac32)*Jac33
;

ddgSST3233
=
dg133*DJac123 + dg233*DJac223 + dg333*DJac323 + 
  Jac13*(ddg1133*Jac12 + ddg1233*Jac22 + ddg1333*Jac32) + 
  Jac23*(ddg1233*Jac12 + ddg2233*Jac22 + ddg2333*Jac32) + 
  (ddg1333*Jac12 + ddg2333*Jac22 + ddg3333*Jac32)*Jac33
;

ddgSST3311
=
dg111*DJac133 + dg211*DJac233 + dg311*DJac333 + 
  2.*(ddg2311*Jac23*Jac33 + Jac13*(ddg1211*Jac23 + ddg1311*Jac33)) + 
  ddg1111*pow2(Jac13) + ddg2211*pow2(Jac23) + ddg3311*pow2(Jac33)
;

ddgSST3312
=
dg112*DJac133 + dg212*DJac233 + dg312*DJac333 + 
  2.*(ddg2312*Jac23*Jac33 + Jac13*(ddg1212*Jac23 + ddg1312*Jac33)) + 
  ddg1112*pow2(Jac13) + ddg2212*pow2(Jac23) + ddg3312*pow2(Jac33)
;

ddgSST3313
=
dg113*DJac133 + dg213*DJac233 + dg313*DJac333 + 
  2.*(ddg2313*Jac23*Jac33 + Jac13*(ddg1213*Jac23 + ddg1313*Jac33)) + 
  ddg1113*pow2(Jac13) + ddg2213*pow2(Jac23) + ddg3313*pow2(Jac33)
;

ddgSST3321
=
dg112*DJac133 + dg212*DJac233 + dg312*DJac333 + 
  2.*(ddg2312*Jac23*Jac33 + Jac13*(ddg1212*Jac23 + ddg1312*Jac33)) + 
  ddg1112*pow2(Jac13) + ddg2212*pow2(Jac23) + ddg3312*pow2(Jac33)
;

ddgSST3322
=
dg122*DJac133 + dg222*DJac233 + dg322*DJac333 + 
  2.*(ddg2322*Jac23*Jac33 + Jac13*(ddg1222*Jac23 + ddg1322*Jac33)) + 
  ddg1122*pow2(Jac13) + ddg2222*pow2(Jac23) + ddg3322*pow2(Jac33)
;

ddgSST3323
=
dg123*DJac133 + dg223*DJac233 + dg323*DJac333 + 
  2.*(ddg2323*Jac23*Jac33 + Jac13*(ddg1223*Jac23 + ddg1323*Jac33)) + 
  ddg1123*pow2(Jac13) + ddg2223*pow2(Jac23) + ddg3323*pow2(Jac33)
;

ddgSST3331
=
dg113*DJac133 + dg213*DJac233 + dg313*DJac333 + 
  2.*(ddg2313*Jac23*Jac33 + Jac13*(ddg1213*Jac23 + ddg1313*Jac33)) + 
  ddg1113*pow2(Jac13) + ddg2213*pow2(Jac23) + ddg3313*pow2(Jac33)
;

ddgSST3332
=
dg123*DJac133 + dg223*DJac233 + dg323*DJac333 + 
  2.*(ddg2323*Jac23*Jac33 + Jac13*(ddg1223*Jac23 + ddg1323*Jac33)) + 
  ddg1123*pow2(Jac13) + ddg2223*pow2(Jac23) + ddg3323*pow2(Jac33)
;

ddgSST3333
=
dg133*DJac133 + dg233*DJac233 + dg333*DJac333 + 
  2.*(ddg2333*Jac23*Jac33 + Jac13*(ddg1233*Jac23 + ddg1333*Jac33)) + 
  ddg1133*pow2(Jac13) + ddg2233*pow2(Jac23) + ddg3333*pow2(Jac33)
;

dKSST111
=
dK111*Jac11 + dK211*Jac21 + dK311*Jac31
;

dKSST112
=
dK112*Jac11 + dK212*Jac21 + dK312*Jac31
;

dKSST113
=
dK113*Jac11 + dK213*Jac21 + dK313*Jac31
;

dKSST121
=
dK112*Jac11 + dK212*Jac21 + dK312*Jac31
;

dKSST122
=
dK122*Jac11 + dK222*Jac21 + dK322*Jac31
;

dKSST123
=
dK123*Jac11 + dK223*Jac21 + dK323*Jac31
;

dKSST131
=
dK113*Jac11 + dK213*Jac21 + dK313*Jac31
;

dKSST132
=
dK123*Jac11 + dK223*Jac21 + dK323*Jac31
;

dKSST133
=
dK133*Jac11 + dK233*Jac21 + dK333*Jac31
;

dKSST211
=
dK111*Jac12 + dK211*Jac22 + dK311*Jac32
;

dKSST212
=
dK112*Jac12 + dK212*Jac22 + dK312*Jac32
;

dKSST213
=
dK113*Jac12 + dK213*Jac22 + dK313*Jac32
;

dKSST221
=
dK112*Jac12 + dK212*Jac22 + dK312*Jac32
;

dKSST222
=
dK122*Jac12 + dK222*Jac22 + dK322*Jac32
;

dKSST223
=
dK123*Jac12 + dK223*Jac22 + dK323*Jac32
;

dKSST231
=
dK113*Jac12 + dK213*Jac22 + dK313*Jac32
;

dKSST232
=
dK123*Jac12 + dK223*Jac22 + dK323*Jac32
;

dKSST233
=
dK133*Jac12 + dK233*Jac22 + dK333*Jac32
;

dKSST311
=
dK111*Jac13 + dK211*Jac23 + dK311*Jac33
;

dKSST312
=
dK112*Jac13 + dK212*Jac23 + dK312*Jac33
;

dKSST313
=
dK113*Jac13 + dK213*Jac23 + dK313*Jac33
;

dKSST321
=
dK112*Jac13 + dK212*Jac23 + dK312*Jac33
;

dKSST322
=
dK122*Jac13 + dK222*Jac23 + dK322*Jac33
;

dKSST323
=
dK123*Jac13 + dK223*Jac23 + dK323*Jac33
;

dKSST331
=
dK113*Jac13 + dK213*Jac23 + dK313*Jac33
;

dKSST332
=
dK123*Jac13 + dK223*Jac23 + dK323*Jac33
;

dKSST333
=
dK133*Jac13 + dK233*Jac23 + dK333*Jac33
;

dg111
=
dgSST111
;

dg112
=
dgSST121
;

dg113
=
dgSST131
;

dg122
=
dgSST122
;

dg123
=
dgSST132
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
dgSST221
;

dg213
=
dgSST231
;

dg222
=
dgSST222
;

dg223
=
dgSST232
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
dgSST321
;

dg313
=
dgSST331
;

dg322
=
dgSST322
;

dg323
=
dgSST332
;

dg333
=
dgSST333
;

ddg1111
=
ddgSST1111
;

ddg1112
=
ddgSST1121
;

ddg1113
=
ddgSST1131
;

ddg1122
=
ddgSST1122
;

ddg1123
=
ddgSST1132
;

ddg1133
=
ddgSST1133
;

ddg1211
=
ddgSST2111
;

ddg1212
=
ddgSST2121
;

ddg1213
=
ddgSST2131
;

ddg1222
=
ddgSST2122
;

ddg1223
=
ddgSST2132
;

ddg1233
=
ddgSST2133
;

ddg1311
=
ddgSST3111
;

ddg1312
=
ddgSST3121
;

ddg1313
=
ddgSST3131
;

ddg1322
=
ddgSST3122
;

ddg1323
=
ddgSST3132
;

ddg1333
=
ddgSST3133
;

ddg2211
=
ddgSST2211
;

ddg2212
=
ddgSST2221
;

ddg2213
=
ddgSST2231
;

ddg2222
=
ddgSST2222
;

ddg2223
=
ddgSST2232
;

ddg2233
=
ddgSST2233
;

ddg2311
=
ddgSST3211
;

ddg2312
=
ddgSST3221
;

ddg2313
=
ddgSST3231
;

ddg2322
=
ddgSST3222
;

ddg2323
=
ddgSST3232
;

ddg2333
=
ddgSST3233
;

ddg3311
=
ddgSST3311
;

ddg3312
=
ddgSST3321
;

ddg3313
=
ddgSST3331
;

ddg3322
=
ddgSST3322
;

ddg3323
=
ddgSST3332
;

ddg3333
=
ddgSST3333
;

dK111
=
dKSST111
;

dK112
=
dKSST121
;

dK113
=
dKSST131
;

dK122
=
dKSST122
;

dK123
=
dKSST132
;

dK133
=
dKSST133
;

dK211
=
dKSST211
;

dK212
=
dKSST221
;

dK213
=
dKSST231
;

dK222
=
dKSST222
;

dK223
=
dKSST232
;

dK233
=
dKSST233
;

dK311
=
dKSST311
;

dK312
=
dKSST321
;

dK313
=
dKSST331
;

dK322
=
dKSST322
;

dK323
=
dKSST332
;

dK333
=
dKSST333
;


} 

ginvdet11
=
g22[ijk]*g33[ijk] - pow2(g23[ijk])
;

ginvdet12
=
g13[ijk]*g23[ijk] - g12[ijk]*g33[ijk]
;

ginvdet13
=
-(g13[ijk]*g22[ijk]) + g12[ijk]*g23[ijk]
;

ginvdet22
=
g11[ijk]*g33[ijk] - pow2(g13[ijk])
;

ginvdet23
=
g12[ijk]*g13[ijk] - g11[ijk]*g23[ijk]
;

ginvdet33
=
g11[ijk]*g22[ijk] - pow2(g12[ijk])
;

detg
=
ginvdet11*g11[ijk] + ginvdet12*g12[ijk] + ginvdet13*g13[ijk]
;

ginv11
=
ginvdet11/detg
;

ginv12
=
ginvdet12/detg
;

ginv13
=
ginvdet13/detg
;

ginv21
=
ginvdet12/detg
;

ginv22
=
ginvdet22/detg
;

ginv23
=
ginvdet23/detg
;

ginv31
=
ginvdet13/detg
;

ginv32
=
ginvdet23/detg
;

ginv33
=
ginvdet33/detg
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
gammado111*ginv21 + gammado211*ginv22 + gammado311*ginv23
;

gamma212
=
gammado112*ginv21 + gammado212*ginv22 + gammado312*ginv23
;

gamma213
=
gammado113*ginv21 + gammado213*ginv22 + gammado313*ginv23
;

gamma222
=
gammado122*ginv21 + gammado222*ginv22 + gammado322*ginv23
;

gamma223
=
gammado123*ginv21 + gammado223*ginv22 + gammado323*ginv23
;

gamma233
=
gammado133*ginv21 + gammado233*ginv22 + gammado333*ginv23
;

gamma311
=
gammado111*ginv31 + gammado211*ginv32 + gammado311*ginv33
;

gamma312
=
gammado112*ginv31 + gammado212*ginv32 + gammado312*ginv33
;

gamma313
=
gammado113*ginv31 + gammado213*ginv32 + gammado313*ginv33
;

gamma322
=
gammado122*ginv31 + gammado222*ginv32 + gammado322*ginv33
;

gamma323
=
gammado123*ginv31 + gammado223*ginv32 + gammado323*ginv33
;

gamma333
=
gammado133*ginv31 + gammado233*ginv32 + gammado333*ginv33
;

R11
=
(ddg1212 + 0.5*(-ddg1122 - ddg2211) + gamma112*gammado112 - 
     gamma111*gammado122 + gamma212*gammado212 - gamma211*gammado222 + 
     gamma312*gammado312 - gamma311*gammado322)*ginv22 + 
  (ddg1213 + 0.5*(-ddg1123 - ddg2311) + gamma112*gammado113 - 
     gamma111*gammado123 + gamma212*gammado213 - gamma211*gammado223 + 
     gamma312*gammado313 - gamma311*gammado323)*ginv23 + 
  0.5*((ddg1112 - ddg1211)*ginv12 + (ddg1113 - ddg1311)*ginv13 + 
     (-ddg1112 + ddg1211)*ginv21 + (-ddg1113 + ddg1311)*ginv31) + 
  (ddg1312 + 0.5*(-ddg1123 - ddg2311) + gamma112*gammado113 - 
     gamma111*gammado123 + gamma212*gammado213 - gamma211*gammado223 + 
     gamma312*gammado313 - gamma311*gammado323)*ginv32 + 
  (ddg1313 + 0.5*(-ddg1133 - ddg3311) + gamma113*gammado113 - 
     gamma111*gammado133 + gamma213*gammado213 - gamma211*gammado233 + 
     gamma313*gammado313 - gamma311*gammado333)*ginv33
;

R12
=
(0.5*(ddg1122 - ddg1212) - gamma112*gammado112 + gamma111*gammado122 - 
     gamma212*gammado212 + gamma211*gammado222 - gamma312*gammado312 + 
     gamma311*gammado322)*ginv12 + 
  (0.5*(ddg1123 - ddg1312) - gamma112*gammado113 + gamma111*gammado123 - 
     gamma212*gammado213 + gamma211*gammado223 - gamma312*gammado313 + 
     gamma311*gammado323)*ginv13 + 
  0.5*((-ddg1212 + ddg2211)*ginv21 + (ddg2213 - ddg2312)*ginv23 + 
     (-ddg1213 + ddg2311)*ginv31) + 
  (0.5*(-ddg1223 + ddg1322) + gamma113*gammado122 - gamma112*gammado123 + 
     gamma213*gammado222 - gamma212*gammado223 + gamma313*gammado322 - 
     gamma312*gammado323)*ginv32 + 
  (0.5*(-ddg1233 + ddg1323 + ddg2313 - ddg3312) + gamma113*gammado123 - 
     gamma112*gammado133 + gamma213*gammado223 - gamma212*gammado233 + 
     gamma313*gammado323 - gamma312*gammado333)*ginv33
;

R13
=
(0.5*(-ddg1312 + ddg2311) - gamma112*gammado113 + gamma111*gammado123 - 
     gamma212*gammado213 + gamma211*gammado223 - gamma312*gammado313 + 
     gamma311*gammado323)*ginv21 + 
  (0.5*(ddg1223 - ddg1322 - ddg2213 + ddg2312) - gamma113*gammado122 + 
     gamma112*gammado123 - gamma213*gammado222 + gamma212*gammado223 - 
     gamma313*gammado322 + gamma312*gammado323)*ginv22 + 
  0.5*((ddg1123 - ddg1213)*ginv12 + (ddg1133 - ddg1313)*ginv13 + 
     (ddg1233 - ddg1323)*ginv23) + 
  (0.5*(-ddg1313 + ddg3311) - gamma113*gammado113 + gamma111*gammado133 - 
     gamma213*gammado213 + gamma211*gammado233 - gamma313*gammado313 + 
     gamma311*gammado333)*ginv31 + 
  (0.5*(-ddg2313 + ddg3312) - gamma113*gammado123 + gamma112*gammado133 - 
     gamma213*gammado223 + gamma212*gammado233 - gamma313*gammado323 + 
     gamma312*gammado333)*ginv32
;

R22
=
(ddg1212 + 0.5*(-ddg1122 - ddg2211) + gamma112*gammado112 - 
     gamma111*gammado122 + gamma212*gammado212 - gamma211*gammado222 + 
     gamma312*gammado312 - gamma311*gammado322)*ginv11 + 
  (ddg1223 + 0.5*(-ddg1322 - ddg2213) - gamma113*gammado122 + 
     gamma112*gammado123 - gamma213*gammado222 + gamma212*gammado223 - 
     gamma313*gammado322 + gamma312*gammado323)*ginv13 + 
  (0.5*(-ddg1322 - ddg2213) + ddg2312 - gamma113*gammado122 + 
     gamma112*gammado123 - gamma213*gammado222 + gamma212*gammado223 - 
     gamma313*gammado322 + gamma312*gammado323)*ginv31 + 
  0.5*((ddg1222 - ddg2212)*ginv12 + (-ddg1222 + ddg2212)*ginv21 + 
     (ddg2223 - ddg2322)*ginv23 + (-ddg2223 + ddg2322)*ginv32) + 
  (ddg2323 + 0.5*(-ddg2233 - ddg3322) + gamma123*gammado123 - 
     gamma122*gammado133 + gamma223*gammado223 - gamma222*gammado233 + 
     gamma323*gammado323 - gamma322*gammado333)*ginv33
;

R23
=
(0.5*(-ddg1123 + ddg1213 + ddg1312 - ddg2311) + gamma112*gammado113 - 
     gamma111*gammado123 + gamma212*gammado213 - gamma211*gammado223 + 
     gamma312*gammado313 - gamma311*gammado323)*ginv11 + 
  (0.5*(ddg1322 - ddg2312) + gamma113*gammado122 - gamma112*gammado123 + 
     gamma213*gammado222 - gamma212*gammado223 + gamma313*gammado322 - 
     gamma312*gammado323)*ginv12 + 
  0.5*((ddg1233 - ddg2313)*ginv13 + (-ddg1223 + ddg2213)*ginv21 + 
     (ddg2233 - ddg2323)*ginv23) + 
  (0.5*(-ddg1323 + ddg3312) - gamma113*gammado123 + gamma112*gammado133 - 
     gamma213*gammado223 + gamma212*gammado233 - gamma313*gammado323 + 
     gamma312*gammado333)*ginv31 + 
  (0.5*(-ddg2323 + ddg3322) - gamma123*gammado123 + gamma122*gammado133 - 
     gamma223*gammado223 + gamma222*gammado233 - gamma323*gammado323 + 
     gamma322*gammado333)*ginv32
;

R33
=
(ddg1313 + 0.5*(-ddg1133 - ddg3311) + gamma113*gammado113 - 
     gamma111*gammado133 + gamma213*gammado213 - gamma211*gammado233 + 
     gamma313*gammado313 - gamma311*gammado333)*ginv11 + 
  (ddg1323 + 0.5*(-ddg1233 - ddg3312) + gamma113*gammado123 - 
     gamma112*gammado133 + gamma213*gammado223 - gamma212*gammado233 + 
     gamma313*gammado323 - gamma312*gammado333)*ginv12 + 
  (ddg2313 + 0.5*(-ddg1233 - ddg3312) + gamma113*gammado123 - 
     gamma112*gammado133 + gamma213*gammado223 - gamma212*gammado233 + 
     gamma313*gammado323 - gamma312*gammado333)*ginv21 + 
  (ddg2323 + 0.5*(-ddg2233 - ddg3322) + gamma123*gammado123 - 
     gamma122*gammado133 + gamma223*gammado223 - gamma222*gammado233 + 
     gamma323*gammado323 - gamma322*gammado333)*ginv22 + 
  0.5*((ddg1333 - ddg3313)*ginv13 + (ddg2333 - ddg3323)*ginv23 + 
     (-ddg1333 + ddg3313)*ginv31 + (-ddg2333 + ddg3323)*ginv32)
;

R
=
ginv11*R11 + (ginv12 + ginv21)*R12 + (ginv13 + ginv31)*R13 + ginv22*R22 + 
  (ginv23 + ginv32)*R23 + ginv33*R33
;

R1212
=
-2.*R12*g12[ijk] + R11*g22[ijk] + g11[ijk]*(R22 - 0.5*R*g22[ijk]) + 
  0.5*R*pow2(g12[ijk])
;

R1213
=
-(R13*g12[ijk]) - R12*g13[ijk] + 0.5*R*g12[ijk]*g13[ijk] + R11*g23[ijk] + 
  g11[ijk]*(R23 - 0.5*R*g23[ijk])
;

R1223
=
-(R22*g13[ijk]) - R13*g22[ijk] + 0.5*R*g13[ijk]*g22[ijk] + R12*g23[ijk] + 
  g12[ijk]*(R23 - 0.5*R*g23[ijk])
;

R1313
=
-2.*R13*g13[ijk] + R11*g33[ijk] + g11[ijk]*(R33 - 0.5*R*g33[ijk]) + 
  0.5*R*pow2(g13[ijk])
;

R1323
=
-(R23*g13[ijk]) - R13*g23[ijk] + 0.5*R*g13[ijk]*g23[ijk] + R12*g33[ijk] + 
  g12[ijk]*(R33 - 0.5*R*g33[ijk])
;

R2323
=
-2.*R23*g23[ijk] + R22*g33[ijk] + g22[ijk]*(R33 - 0.5*R*g33[ijk]) + 
  0.5*R*pow2(g23[ijk])
;

trK
=
ginv11*K11[ijk] + (ginv12 + ginv21)*K12[ijk] + (ginv13 + ginv31)*K13[ijk] + 
  ginv22*K22[ijk] + (ginv23 + ginv32)*K23[ijk] + ginv33*K33[ijk]
;

Riemm1212
=
R1212 + K11[ijk]*K22[ijk] - pow2(K12[ijk])
;

Riemm1213
=
R1213 - K12[ijk]*K13[ijk] + K11[ijk]*K23[ijk]
;

Riemm1223
=
R1223 - K13[ijk]*K22[ijk] + K12[ijk]*K23[ijk]
;

Riemm1313
=
R1313 + K11[ijk]*K33[ijk] - pow2(K13[ijk])
;

Riemm1323
=
R1323 - K13[ijk]*K23[ijk] + K12[ijk]*K33[ijk]
;

Riemm2323
=
R2323 + K22[ijk]*K33[ijk] - pow2(K23[ijk])
;

Riemm112
=
dK112 - dK211 + gamma112*K11[ijk] - gamma111*K12[ijk] + gamma212*K12[ijk] + 
  gamma312*K13[ijk] - gamma211*K22[ijk] - gamma311*K23[ijk]
;

Riemm113
=
dK113 - dK311 + gamma113*K11[ijk] + gamma213*K12[ijk] - gamma111*K13[ijk] + 
  gamma313*K13[ijk] - gamma211*K23[ijk] - gamma311*K33[ijk]
;

Riemm123
=
dK213 - dK312 + gamma113*K12[ijk] - gamma112*K13[ijk] + gamma213*K22[ijk] - 
  gamma212*K23[ijk] + gamma313*K23[ijk] - gamma312*K33[ijk]
;

Riemm212
=
dK122 - dK212 + gamma122*K11[ijk] - gamma112*K12[ijk] + gamma222*K12[ijk] + 
  gamma322*K13[ijk] - gamma212*K22[ijk] - gamma312*K23[ijk]
;

Riemm213
=
dK123 - dK312 + gamma123*K11[ijk] + gamma223*K12[ijk] - gamma112*K13[ijk] + 
  gamma323*K13[ijk] - gamma212*K23[ijk] - gamma312*K33[ijk]
;

Riemm223
=
dK223 - dK322 + gamma123*K12[ijk] - gamma122*K13[ijk] + gamma223*K22[ijk] - 
  gamma222*K23[ijk] + gamma323*K23[ijk] - gamma322*K33[ijk]
;

Riemm312
=
dK123 - dK213 + gamma123*K11[ijk] - gamma113*K12[ijk] + gamma223*K12[ijk] + 
  gamma323*K13[ijk] - gamma213*K22[ijk] - gamma313*K23[ijk]
;

Riemm313
=
dK133 - dK313 + gamma133*K11[ijk] + gamma233*K12[ijk] - gamma113*K13[ijk] + 
  gamma333*K13[ijk] - gamma213*K23[ijk] - gamma313*K33[ijk]
;

Riemm323
=
dK233 - dK323 + gamma133*K12[ijk] - gamma123*K13[ijk] + gamma233*K22[ijk] - 
  gamma223*K23[ijk] + gamma333*K23[ijk] - gamma323*K33[ijk]
;

Riemm11
=
R11 - (ginv23 + ginv32)*K12[ijk]*K13[ijk] + 
  K11[ijk]*(trK - ginv12*K12[ijk] - ginv13*K13[ijk]) - 
  K11[ijk]*(ginv21*K12[ijk] + ginv31*K13[ijk]) - ginv11*pow2(K11[ijk]) - 
  ginv22*pow2(K12[ijk]) - ginv33*pow2(K13[ijk])
;

Riemm12
=
R12 - (ginv21*K11[ijk] + ginv23*K13[ijk])*K22[ijk] - 
  (ginv31*K11[ijk] + ginv33*K13[ijk])*K23[ijk] + 
  K12[ijk]*(trK - ginv11*K11[ijk] - ginv13*K13[ijk] - ginv22*K22[ijk] - 
     ginv32*K23[ijk]) - ginv12*pow2(K12[ijk])
;

Riemm13
=
R13 - (ginv21*K11[ijk] + ginv22*K12[ijk])*K23[ijk] - 
  (ginv31*K11[ijk] + ginv32*K12[ijk])*K33[ijk] + 
  K13[ijk]*(trK - ginv11*K11[ijk] - ginv12*K12[ijk] - ginv23*K23[ijk] - 
     ginv33*K33[ijk]) - ginv13*pow2(K13[ijk])
;

Riemm22
=
R22 - (ginv31*K12[ijk] + ginv32*K22[ijk])*K23[ijk] - 
  K12[ijk]*(ginv12*K22[ijk] + ginv13*K23[ijk]) + 
  K22[ijk]*(trK - ginv21*K12[ijk] - ginv23*K23[ijk]) - 
  ginv11*pow2(K12[ijk]) - ginv22*pow2(K22[ijk]) - ginv33*pow2(K23[ijk])
;

Riemm23
=
R23 - K13[ijk]*(ginv11*K12[ijk] + ginv12*K22[ijk]) - 
  ginv32*K22[ijk]*K33[ijk] - K12[ijk]*(ginv21*K23[ijk] + ginv31*K33[ijk]) + 
  K23[ijk]*(trK - ginv13*K13[ijk] - ginv22*K22[ijk] - ginv33*K33[ijk]) - 
  ginv23*pow2(K23[ijk])
;

Riemm33
=
R33 + (trK - ginv31*K13[ijk] - (ginv23 + ginv32)*K23[ijk])*K33[ijk] - 
  K13[ijk]*((ginv12 + ginv21)*K23[ijk] + ginv13*K33[ijk]) - 
  ginv11*pow2(K13[ijk]) - ginv22*pow2(K23[ijk]) - ginv33*pow2(K33[ijk])
;

tfactor
=
oosqrt2
;


if (Csqr) { 

Sqrtdetg
=
Sqrt(detg)
;

E11
=
Riemm11
;

E12
=
Riemm12
;

E13
=
Riemm13
;

E21
=
Riemm12
;

E22
=
Riemm22
;

E23
=
Riemm23
;

E31
=
Riemm13
;

E32
=
Riemm23
;

E33
=
Riemm33
;

B11
=
(-(Riemm123*g11[ijk]) + Riemm113*g12[ijk] - Riemm112*g13[ijk])/Sqrtdetg
;

B12
=
(-(Riemm223*g11[ijk]) + Riemm213*g12[ijk] - Riemm212*g13[ijk])/Sqrtdetg
;

B13
=
(-(Riemm323*g11[ijk]) + Riemm313*g12[ijk] - Riemm312*g13[ijk])/Sqrtdetg
;

B21
=
(-(Riemm123*g12[ijk]) + Riemm113*g22[ijk] - Riemm112*g23[ijk])/Sqrtdetg
;

B22
=
(-(Riemm223*g12[ijk]) + Riemm213*g22[ijk] - Riemm212*g23[ijk])/Sqrtdetg
;

B23
=
(-(Riemm323*g12[ijk]) + Riemm313*g22[ijk] - Riemm312*g23[ijk])/Sqrtdetg
;

B31
=
(-(Riemm123*g13[ijk]) + Riemm113*g23[ijk] - Riemm112*g33[ijk])/Sqrtdetg
;

B32
=
(-(Riemm223*g13[ijk]) + Riemm213*g23[ijk] - Riemm212*g33[ijk])/Sqrtdetg
;

B33
=
(-(Riemm323*g13[ijk]) + Riemm313*g23[ijk] - Riemm312*g33[ijk])/Sqrtdetg
;

Csqr[ijk]
=
(-((B23 + B32)*B33) + (E23 + E32)*E33)*ginv32*ginv33 + 
  ginv12*((-(B11*(B23 + B32)) + E11*(E23 + E32))*ginv13 + 
     2.*(-(B12*B21) + E12*E21)*ginv21 + 
     (-((B12 + B21)*B22) + (E12 + E21)*E22)*ginv22 + 
     (-(B12*B23) - B21*B32 + E12*E23 + E21*E32)*ginv23 + 
     (-(B13*B21) - B12*B31 + E13*E21 + E12*E31)*ginv31 + 
     (-(B22*(B13 + B31)) + E22*(E13 + E31))*ginv32 + 
     (-(B13*B23) - B31*B32 + E13*E23 + E31*E32)*ginv33) + 
  ginv21*((-((B12 + B21)*B22) + (E12 + E21)*E22)*ginv22 + 
     (-(B22*(B13 + B31)) + E22*(E13 + E31))*ginv23 + 
     (-(B11*(B23 + B32)) + E11*(E23 + E32))*ginv31 + 
     (-(B12*B23) - B21*B32 + E12*E23 + E21*E32)*ginv32 + 
     (-(B13*B23) - B31*B32 + E13*E23 + E31*E32)*ginv33) + 
  ginv13*((-(B13*B21) - B12*B31 + E13*E21 + E12*E31)*ginv21 + 
     (-(B21*B23) - B12*B32 + E21*E23 + E12*E32)*ginv22 + 
     (-((B12 + B21)*B33) + (E12 + E21)*E33)*ginv23 + 
     2.*(-(B13*B31) + E13*E31)*ginv31 + 
     (-(B23*B31) - B13*B32 + E23*E31 + E13*E32)*ginv32 + 
     (-((B13 + B31)*B33) + (E13 + E31)*E33)*ginv33) + 
  ginv31*((-((B12 + B21)*B33) + (E12 + E21)*E33)*ginv32 + 
     (-((B13 + B31)*B33) + (E13 + E31)*E33)*ginv33) + 
  ginv23*((-(B23*B31) - B13*B32 + E23*E31 + E13*E32)*ginv31 + 
     2.*(-(B23*B32) + E23*E32)*ginv32 + 
     (-((B23 + B32)*B33) + (E23 + E32)*E33)*ginv33) + 
  ginv11*((-(B11*(B12 + B21)) + E11*(E12 + E21))*(ginv12 + ginv21) + 
     (-(B11*(B13 + B31)) + E11*(E13 + E31))*(ginv13 + ginv31) + 
     (-(B12*B13) - B21*B31 + E12*E13 + E21*E31)*(ginv23 + ginv32) + 
     ginv22*(-pow2(B12) - pow2(B21) + pow2(E12) + pow2(E21)) + 
     ginv33*(-pow2(B13) - pow2(B31) + pow2(E13) + pow2(E31))) + 
  ginv22*((-(B21*B23) - B12*B32 + E21*E23 + E12*E32)*ginv31 + 
     (-(B22*(B23 + B32)) + E22*(E23 + E32))*(ginv23 + ginv32) + 
     ginv33*(-pow2(B23) - pow2(B32) + pow2(E23) + pow2(E32))) + 
  (-pow2(B11) + pow2(E11))*pow2(ginv11) + 
  (-(B11*B22) + E11*E22)*(pow2(ginv12) + pow2(ginv21)) + 
  (-pow2(B22) + pow2(E22))*pow2(ginv22) + 
  (-(B11*B33) + E11*E33)*(pow2(ginv13) + pow2(ginv31)) + 
  (-(B22*B33) + E22*E33)*(pow2(ginv23) + pow2(ginv32)) + 
  (-pow2(B33) + pow2(E33))*pow2(ginv33)
;


} 

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



/* conditional */
if (kinnersley) {

n0
=
0.5
;

n1
=
-0.5*cosphi*sintheta
;

n2
=
-0.5*sinphi*sintheta
;

n3
=
-0.5*costheta
;

rm1
=
cosphi*costheta*tfactor
;

rm2
=
costheta*sinphi*tfactor
;

rm3
=
-(sintheta*tfactor)
;

im1
=
-(sinphi*tfactor)
;

im2
=
cosphi*tfactor
;

im3
=
0
;


} else { /* if (!kinnersley) */

v11
=
-y
;

v12
=
x
;

v13
=
0
;

v21
=
x
;

v22
=
y
;

v23
=
z
;

v31
=
Sqrt(detg)*(ginv13*(-(v12*v21) + v11*v22) + ginv12*(v13*v21 - v11*v23) + 
    ginv11*(-(v13*v22) + v12*v23))
;

v32
=
Sqrt(detg)*(ginv23*(-(v12*v21) + v11*v22) + ginv22*(v13*v21 - v11*v23) + 
    ginv21*(-(v13*v22) + v12*v23))
;

v33
=
Sqrt(detg)*(ginv33*(-(v12*v21) + v11*v22) + ginv32*(v13*v21 - v11*v23) + 
    ginv31*(-(v13*v22) + v12*v23))
;

w11
=
2.*(v11*(v12*g12[ijk] + v13*g13[ijk]) + v12*v13*g23[ijk]) + 
  g11[ijk]*pow2(v11) + g22[ijk]*pow2(v12) + g33[ijk]*pow2(v13)
;

v11
=
v11/Sqrt(w11)
;

v12
=
v12/Sqrt(w11)
;

v13
=
v13/Sqrt(w11)
;

w12
=
v21*(v11*g11[ijk] + v12*g12[ijk] + v13*g13[ijk]) + 
  v11*(v22*g12[ijk] + v23*g13[ijk]) + v22*(v12*g22[ijk] + v13*g23[ijk]) + 
  v23*(v12*g23[ijk] + v13*g33[ijk])
;

v21
=
v21 - v11*w12
;

v22
=
v22 - v12*w12
;

v23
=
v23 - v13*w12
;

w22
=
2.*(v21*(v22*g12[ijk] + v23*g13[ijk]) + v22*v23*g23[ijk]) + 
  g11[ijk]*pow2(v21) + g22[ijk]*pow2(v22) + g33[ijk]*pow2(v23)
;

v21
=
v21/Sqrt(w22)
;

v22
=
v22/Sqrt(w22)
;

v23
=
v23/Sqrt(w22)
;

w13
=
v31*(v11*g11[ijk] + v12*g12[ijk] + v13*g13[ijk]) + 
  v11*(v32*g12[ijk] + v33*g13[ijk]) + v32*(v12*g22[ijk] + v13*g23[ijk]) + 
  v33*(v12*g23[ijk] + v13*g33[ijk])
;

w23
=
v31*(v21*g11[ijk] + v22*g12[ijk] + v23*g13[ijk]) + 
  v21*(v32*g12[ijk] + v33*g13[ijk]) + v32*(v22*g22[ijk] + v23*g23[ijk]) + 
  v33*(v22*g23[ijk] + v23*g33[ijk])
;

v31
=
v31 - v11*w13 - v21*w23
;

v32
=
v32 - v12*w13 - v22*w23
;

v33
=
v33 - v13*w13 - v23*w23
;

w33
=
2.*(v31*(v32*g12[ijk] + v33*g13[ijk]) + v32*v33*g23[ijk]) + 
  g11[ijk]*pow2(v31) + g22[ijk]*pow2(v32) + g33[ijk]*pow2(v33)
;

v31
=
v31/Sqrt(w33)
;

v32
=
v32/Sqrt(w33)
;

v33
=
v33/Sqrt(w33)
;

vecr1
=
v21
;

vecr2
=
v22
;

vecr3
=
v23
;

vect1
=
v31
;

vect2
=
v32
;

vect3
=
v33
;

vecp1
=
v11
;

vecp2
=
v12
;

vecp3
=
v13
;

vecu0
=
1/alpha[ijk]
;

vecu1
=
-(beta1[ijk]/alpha[ijk])
;

vecu2
=
-(beta2[ijk]/alpha[ijk])
;

vecu3
=
-(beta3[ijk]/alpha[ijk])
;

n0
=
tfactor*vecu0
;

n1
=
tfactor*(-vecr1 + vecu1)
;

n2
=
tfactor*(-vecr2 + vecu2)
;

n3
=
tfactor*(-vecr3 + vecu3)
;

rm1
=
tfactor*vect1
;

rm2
=
tfactor*vect2
;

rm3
=
tfactor*vect3
;

im1
=
tfactor*vecp1
;

im2
=
tfactor*vecp2
;

im3
=
tfactor*vecp3
;

}
/* if (kinnersley) */


rmb1
=
rm1
;

rmb2
=
rm2
;

rmb3
=
rm3
;

imb1
=
-im1
;

imb2
=
-im2
;

imb3
=
-im3
;



/* conditional */
if (gausscodacciflat) {

rpsi4[ijk]
=
pow2(imb3)*(2.*n0*n1*Riemm313 + Riemm1313*pow2(n1) + Riemm2323*pow2(n2)) + 
  pow2(imb1)*(2.*(-(n0*(n2*Riemm112 + n3*Riemm113)) + n2*n3*Riemm1213) + 
     Riemm11*pow2(n0) + Riemm1212*pow2(n2) + Riemm1313*pow2(n3)) + 
  pow2(n0)*(2.*(imb2*imb3*Riemm23 - rmb1*(Riemm12*rmb2 + Riemm13*rmb3)) + 
     Riemm22*pow2(imb2) + Riemm33*pow2(imb3) - Riemm11*pow2(rmb1)) + 
  pow2(n3)*(-2.*Riemm1323*rmb1*rmb2 + Riemm2323*pow2(imb2) - 
     Riemm1313*pow2(rmb1)) + Riemm1212*
   (2.*n1*n2*(-(imb1*imb2) + rmb1*rmb2) + pow2(imb2)*pow2(n1) - 
     pow2(n2)*pow2(rmb1)) + (-2.*n0*n1*Riemm212 + 2.*n0*n3*Riemm223 - 
     Riemm22*pow2(n0) - Riemm1212*pow2(n1) - Riemm2323*pow2(n3))*pow2(rmb2) \
+ 2.*(-(n0*(n1*(Riemm213*rmb2*rmb3 + 
             rmb1*(Riemm112*rmb2 + Riemm113*rmb3)) + 
          n3*(imb2*imb3*Riemm323 + Riemm223*pow2(imb2)))) + 
     imb1*(n0*(n1*(imb2*Riemm112 + imb3*Riemm113) + imb3*n2*Riemm123) + 
        (imb2*Riemm12 + imb3*Riemm13)*pow2(n0)) - 
     imb1*(n1*((imb3*n2 + imb2*n3)*Riemm1213 + imb3*n3*Riemm1313) + 
        imb2*n0*n2*Riemm212 + imb3*
         (n2*(n3*Riemm1323 + n0*Riemm312) + Riemm1223*pow2(n2))) + 
     rmb3*(rmb2*(n1*(-(n2*Riemm1223) + n3*Riemm1323) + n2*n3*Riemm2323 + 
           n0*(-(n1*Riemm312) + n3*Riemm323) - Riemm23*pow2(n0) - 
           Riemm1213*pow2(n1)) + 
        rmb1*(n0*(-(n2*Riemm123) + n3*Riemm313) + Riemm1223*pow2(n2))) + 
     imb2*(imb3*(n1*(n2*Riemm1223 + n0*Riemm213) + n0*n2*Riemm223 + 
           Riemm1213*pow2(n1)) + imb1*Riemm1323*pow2(n3)) + 
     n3*(-(imb2*(imb1*n0*(Riemm123 + Riemm213) + 
             imb3*(n1*Riemm1323 + n2*Riemm2323))) - 
        imb1*imb3*n0*Riemm313 + (n1*Riemm1313 + n2*Riemm1323)*rmb1*rmb3 + 
        Riemm1223*(imb1*imb2*n2 - n2*rmb1*rmb2 - n1*pow2(imb2)) + 
        n0*Riemm113*pow2(rmb1)) - 
     n2*(n0*Riemm223*rmb2*rmb3 + n3*Riemm1213*pow2(rmb1)) + 
     n0*(n3*(Riemm123 + Riemm213)*rmb1*rmb2 + 
        n2*(rmb1*(Riemm212*rmb2 + Riemm312*rmb3) + Riemm323*pow2(imb3) + 
           Riemm112*pow2(rmb1))) + 
     n1*(n0*(imb2*imb3*Riemm312 + Riemm212*pow2(imb2)) + 
        n2*(Riemm1213*rmb1*rmb3 + Riemm1323*pow2(imb3)) + 
        n3*(Riemm1213*rmb1*rmb2 + Riemm1223*pow2(rmb2)))) - 
  (2.*(n1*(n2*Riemm1323 + n0*Riemm313) + n0*n2*Riemm323) + 
     Riemm33*pow2(n0) + Riemm1313*pow2(n1) + Riemm2323*pow2(n2))*pow2(rmb3)
;

ipsi4[ijk]
=
n3*((4.*imb2*n1 - 2.*imb1*n2)*Riemm1223*rmb2 + 
     Riemm1213*((2.*imb2*n1 - 4.*imb1*n2)*rmb1 + 2.*imb1*n1*rmb2) + 
     2.*imb1*n1*Riemm1313*rmb3) + 
  n0*((4.*imb1*n2*Riemm112 - 2.*imb3*n1*Riemm113)*rmb1 + 
     2.*imb1*Riemm123*(n3*rmb2 - n2*rmb3)) - 
  2.*(rmb1*(imb2*n0*n1*Riemm112 + imb1*Riemm11*pow2(n0)) + 
     imb1*(rmb2*(n0*n1*Riemm112 + Riemm12*pow2(n0)) + 
        rmb3*(n0*n1*Riemm113 + Riemm13*pow2(n0))) + 
     (imb2*Riemm1212*rmb2 + imb3*Riemm1313*rmb3)*pow2(n1)) + 
  rmb3*(n1*((-4.*imb3*n2 + 2.*imb2*n3)*Riemm1323 - 
        n0*(2.*imb2*Riemm312 + 4.*imb3*Riemm313)) + 
     n2*(2.*n3*(imb1*Riemm1323 + imb2*Riemm2323) + 
        n0*(2.*(-(imb2*Riemm223) + imb1*Riemm312) - 4.*imb3*Riemm323)) + 
     2.*(n0*n3*(imb1*Riemm313 + imb2*Riemm323) - 
        imb2*(n0*n1*Riemm213 + Riemm23*pow2(n0)) + 
        Riemm1213*(imb1*n1*n2 - imb2*pow2(n1)) + 
        Riemm1223*(-(imb2*n1*n2) + imb1*pow2(n2)) - 
        imb3*(Riemm33*pow2(n0) + Riemm2323*pow2(n2)))) + 
  rmb1*(imb1*(4.*n0*n3*Riemm113 - 
        2.*(Riemm1212*pow2(n2) + Riemm1313*pow2(n3))) + 
     2.*(n3*(imb3*n2*Riemm1323 + n0*(imb2*Riemm123 + imb3*Riemm313)) + 
        imb3*(-(n0*n2*Riemm123) + n1*(n2*Riemm1213 + n3*Riemm1313) + 
           n0*n2*Riemm312 - Riemm13*pow2(n0) + Riemm1223*pow2(n2)) + 
        imb2*(-(n2*n3*Riemm1223) + n2*(n1*Riemm1212 + n0*Riemm212) + 
           n0*n3*Riemm213 - Riemm12*pow2(n0) - Riemm1323*pow2(n3)))) + 
  rmb2*(n1*(2.*n2*(imb1*Riemm1212 - imb3*Riemm1223) - 
        n0*(4.*imb2*Riemm212 + 2.*imb3*Riemm312)) + 
     2.*(imb3*(-(n0*(n1*Riemm213 + n2*Riemm223)) + 
           n3*(n1*Riemm1323 + n2*Riemm2323 + n0*Riemm323) - 
           Riemm23*pow2(n0) - Riemm1213*pow2(n1)) + 
        imb1*(n0*(n2*Riemm212 + n3*Riemm213) - Riemm1323*pow2(n3))) + 
     imb2*(4.*n0*n3*Riemm223 - 2.*(Riemm22*pow2(n0) + Riemm2323*pow2(n3))))
;


} else { /* if (!gausscodacciflat) */

rpsi4[ijk]
=
Riemm23*(imb2*imb3 - rmb2*rmb3) + 
  rmb1*((Riemm1212*rmb2 + Riemm1213*rmb3)*vecr1*vecr2 - 
     Riemm213*rmb2*vecr3) + imb1*
   (imb3*(Riemm13 - Riemm123*vecr2) + 
     imb2*(Riemm12 - Riemm112*vecr1 + Riemm123*vecr3)) + 
  rmb2*((Riemm213 + Riemm312)*rmb3*vecr1 - 
     rmb1*(Riemm212*vecr2 + Riemm123*vecr3)) - 
  rmb1*(Riemm13*rmb3 + rmb2*(Riemm12 + Riemm1223*vecr2*vecr3)) + 
  rmb3*((rmb1*(-Riemm313 + Riemm1313*vecr1) + 
        rmb2*(Riemm1323*vecr1 + Riemm2323*vecr2))*vecr3 + 
     vecr2*(Riemm223*rmb2 + rmb1*(-Riemm312 + Riemm1323*vecr3))) - 
  vecr1*(imb3*(imb1*Riemm113 + imb2*(Riemm213 + Riemm312)) + 
     Riemm1223*rmb2*rmb3*vecr2 + Riemm212*pow2(imb2) + Riemm313*pow2(imb3)) \
- vecr3*(imb3*((imb1*Riemm1313 + imb2*Riemm1323)*vecr1 + 
        (imb1*Riemm1323 + imb2*Riemm2323)*vecr2) + Riemm113*pow2(rmb1)) + 
  vecr3*(Riemm323*(imb2*imb3 - rmb2*rmb3) + 
     imb1*(imb3*Riemm313 + imb2*(Riemm213 + Riemm1223*vecr2)) + 
     Riemm113*pow2(imb1) - Riemm1223*vecr1*pow2(imb2) + 
     Riemm223*(pow2(imb2) - pow2(rmb2))) + 
  vecr1*(imb2*((-(imb1*Riemm1212) + imb3*Riemm1223)*vecr2 - 
        imb1*Riemm1213*vecr3) + 
     rmb1*(Riemm113*rmb3 + rmb2*(Riemm112 + Riemm1213*vecr3)) + 
     Riemm1323*vecr2*pow2(imb3) + (Riemm212 + Riemm1223*vecr3)*pow2(rmb2) + 
     Riemm313*pow2(rmb3)) - vecr2*
   (imb2*imb3*Riemm223 + Riemm1213*vecr3*pow2(rmb1) + 
     Riemm1323*vecr1*pow2(rmb3)) + 
  vecr2*(Riemm123*rmb1*rmb3 + imb1*
      (imb2*Riemm212 + imb3*(Riemm312 - Riemm1213*vecr1)) + 
     Riemm1213*vecr3*pow2(imb1) + Riemm112*(pow2(imb1) - pow2(rmb1)) + 
     Riemm323*(-pow2(imb3) + pow2(rmb3))) + 
  (Riemm1213*(imb2*imb3 - rmb2*rmb3) + 
     0.5*(-(Riemm1212*pow2(rmb2)) - Riemm1313*pow2(rmb3)))*pow2(vecr1) + 
  (-(imb1*imb3*Riemm1223) + Riemm1223*rmb1*rmb3 - 
     0.5*Riemm1212*pow2(rmb1) + 0.5*Riemm2323*(pow2(imb3) - pow2(rmb3)))*
   pow2(vecr2) + (Riemm1323*(imb1*imb2 - rmb1*rmb2) + 
     0.5*(-(Riemm1313*pow2(rmb1)) + Riemm2323*(pow2(imb2) - pow2(rmb2))))*
   pow2(vecr3) + 0.5*(Riemm11*(pow2(imb1) - pow2(rmb1)) + 
     Riemm22*(pow2(imb2) - pow2(rmb2)) + Riemm33*(pow2(imb3) - pow2(rmb3)) + 
     (Riemm1212*pow2(imb2) + Riemm1313*pow2(imb3))*pow2(vecr1) + 
     pow2(imb1)*(Riemm1212*pow2(vecr2) + Riemm1313*pow2(vecr3)))
;

ipsi4[ijk]
=
-(imb2*(Riemm213*rmb1 + 2.*Riemm223*rmb2)*vecr3) + 
  (rmb1*(imb3*Riemm1313*vecr1 - imb2*Riemm1223*vecr2) + 
     imb2*rmb3*(-Riemm323 + Riemm1323*vecr1 + Riemm2323*vecr2))*vecr3 + 
  vecr1*((imb3*(Riemm1213*rmb1 - Riemm1223*rmb2) + 
        imb1*(Riemm1212*rmb2 + Riemm1213*rmb3))*vecr2 + 
     (imb2*Riemm1213*rmb1 + imb3*Riemm1323*rmb2 + imb1*Riemm1313*rmb3)*vecr3\
) + imb3*(rmb2*((Riemm213 + Riemm312)*vecr1 + Riemm223*vecr2) + 
     rmb1*((Riemm123 - Riemm312)*vecr2 - Riemm313*vecr3)) + 
  vecr2*(((-(imb1*Riemm1223) + imb3*Riemm2323)*rmb2 + imb1*Riemm1323*rmb3)*
      vecr3 + rmb1*(imb2*Riemm1212*vecr1 + imb3*Riemm1323*vecr3) + 
     2.*(imb3*Riemm323*rmb3 - imb1*Riemm1213*rmb1*vecr3)) + 
  rmb2*(-(imb1*(Riemm212*vecr2 + Riemm123*vecr3)) + 
     vecr1*(imb1*(Riemm112 + Riemm1213*vecr3) + 
        2.*imb2*(Riemm212 + Riemm1223*vecr3))) - 
  rmb3*(imb1*Riemm13 + imb3*(Riemm33 + 2.*Riemm1323*vecr1*vecr2) + 
     imb2*(Riemm23 + Riemm1223*vecr1*vecr2 + Riemm1213*pow2(vecr1))) + 
  rmb3*((imb1*Riemm113 + imb2*(Riemm213 + Riemm312) + 2.*imb3*Riemm313)*
      vecr1 + (imb1*Riemm123 + imb2*Riemm223)*vecr2 - 
     imb1*(Riemm312*vecr2 + Riemm313*vecr3) - imb3*Riemm1313*pow2(vecr1) + 
     (imb1*Riemm1223 - imb3*Riemm2323)*pow2(vecr2)) + 
  rmb1*((imb2*Riemm112 + imb3*Riemm113)*vecr1 - imb2*Riemm123*vecr3 + 
     imb3*Riemm1223*pow2(vecr2) - 
     imb1*(2.*Riemm112*vecr2 + Riemm1212*pow2(vecr2))) - 
  rmb1*(imb3*Riemm13 + imb1*(Riemm11 + 2.*Riemm113*vecr3 + 
        Riemm1313*pow2(vecr3)) + 
     imb2*(Riemm12 + Riemm212*vecr2 + Riemm1323*pow2(vecr3))) - 
  rmb2*(imb3*(Riemm23 + Riemm323*vecr3 + Riemm1213*pow2(vecr1)) + 
     imb1*(Riemm12 + Riemm213*vecr3 + Riemm1323*pow2(vecr3)) + 
     imb2*(Riemm22 + Riemm1212*pow2(vecr1) + Riemm2323*pow2(vecr3)))
;

}
/* if (gausscodacciflat) */




} endfor_ijk_openmp; /* loop i, j, k */



bampi_openmp_stop


}  /* function */

/* curvature_invariants_N.c */
/* nvars = 38, nauxs = 395, n* = 4914,  n/ = 113,  n+ = 7468, n = 12495, O = 0 */
