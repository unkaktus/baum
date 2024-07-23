/* adm_constraints_N.c */
/* Copyright (C) 1998 Bernd Bruegmann, 18.4.2012 */
/* Produced with Mathematica */

#include "bam.h"
#include "adm.h"

#define Power(x,y) (pow((double) (x), (double) (y)))
#define Sqrt(x)    sqrt(x)
#define Log(x)     log((double) (x))
#define pow2(x)    ((x)*(x))
#define pow4(x)    ((x)*(x)*(x)*(x))
#define pow2inv(x) (1.0/((x)*(x)))

#define Tanh(x)    tanh(x)
#define Sech(x)    (1/cosh(x))
#define Cal(x,y,z) ((x)?(y):(z))




void adm_constraints_N(tVarList *ucon, tVarList *uadm)
{

tL *level = ucon->level;

double *ham = vldataptr(ucon, 0);
double *mom1 = vldataptr(ucon, 1);
double *mom2 = vldataptr(ucon, 2);
double *mom3 = vldataptr(ucon, 3);
double *normham = vldataptr(ucon, 4);
double *normmom1 = vldataptr(ucon, 5);
double *normmom2 = vldataptr(ucon, 6);
double *normmom3 = vldataptr(ucon, 7);
double *g11 = vldataptr(uadm, 0);
double *g12 = vldataptr(uadm, 1);
double *g13 = vldataptr(uadm, 2);
double *g22 = vldataptr(uadm, 3);
double *g23 = vldataptr(uadm, 4);
double *g33 = vldataptr(uadm, 5);
double *K11 = vldataptr(uadm, 6);
double *K12 = vldataptr(uadm, 7);
double *K13 = vldataptr(uadm, 8);
double *K22 = vldataptr(uadm, 9);
double *K23 = vldataptr(uadm, 10);
double *K33 = vldataptr(uadm, 11);
double *psi = vldataptr(uadm, 12);
double *dpop1 = vldataptr(uadm, 13);
double *dpop2 = vldataptr(uadm, 14);
double *dpop3 = vldataptr(uadm, 15);
double *ddpop11 = vldataptr(uadm, 16);
double *ddpop12 = vldataptr(uadm, 17);
double *ddpop13 = vldataptr(uadm, 18);
double *ddpop22 = vldataptr(uadm, 19);
double *ddpop23 = vldataptr(uadm, 20);
double *ddpop33 = vldataptr(uadm, 21);
double *rho = vldataptr(uadm, 22);
double *S1 = vldataptr(uadm, 23);
double *S2 = vldataptr(uadm, 24);
double *S3 = vldataptr(uadm, 25);

const int normConstr = Getv("adm_normalizedConstraints", "yes");
const int TermByTerm = Getv("adm_normalizedConstraints", "TermByTerm");
const int usepsi = 0;
const int order_centered = Geti("order_centered");
const int matter = ( (Getv("physics","matter")) && (!level->shells) );
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

double cdKdA1 = 0.;
double cdKdA2 = 0.;
double cdKdA3 = 0.;
double cdKdB1 = 0.;
double cdKdB2 = 0.;
double cdKdB3 = 0.;
double cdKdC1 = 0.;
double cdKdC2 = 0.;
double cdKdC3 = 0.;
double cdKudd111 = 0.;
double cdKudd112 = 0.;
double cdKudd113 = 0.;
double cdKudd122 = 0.;
double cdKudd123 = 0.;
double cdKudd133 = 0.;
double cdKudd211 = 0.;
double cdKudd212 = 0.;
double cdKudd213 = 0.;
double cdKudd222 = 0.;
double cdKudd223 = 0.;
double cdKudd233 = 0.;
double cdKudd311 = 0.;
double cdKudd312 = 0.;
double cdKudd313 = 0.;
double cdKudd322 = 0.;
double cdKudd323 = 0.;
double cdKudd333 = 0.;
double codelK111 = 0.;
double codelK112 = 0.;
double codelK113 = 0.;
double codelK122 = 0.;
double codelK123 = 0.;
double codelK133 = 0.;
double codelK211 = 0.;
double codelK212 = 0.;
double codelK213 = 0.;
double codelK222 = 0.;
double codelK223 = 0.;
double codelK233 = 0.;
double codelK311 = 0.;
double codelK312 = 0.;
double codelK313 = 0.;
double codelK322 = 0.;
double codelK323 = 0.;
double codelK333 = 0.;
double codelKA111 = 0.;
double codelKA112 = 0.;
double codelKA113 = 0.;
double codelKA122 = 0.;
double codelKA123 = 0.;
double codelKA133 = 0.;
double codelKA211 = 0.;
double codelKA212 = 0.;
double codelKA213 = 0.;
double codelKA222 = 0.;
double codelKA223 = 0.;
double codelKA233 = 0.;
double codelKA311 = 0.;
double codelKA312 = 0.;
double codelKA313 = 0.;
double codelKA322 = 0.;
double codelKA323 = 0.;
double codelKA333 = 0.;
double codelKB111 = 0.;
double codelKB112 = 0.;
double codelKB113 = 0.;
double codelKB122 = 0.;
double codelKB123 = 0.;
double codelKB133 = 0.;
double codelKB211 = 0.;
double codelKB212 = 0.;
double codelKB213 = 0.;
double codelKB222 = 0.;
double codelKB223 = 0.;
double codelKB233 = 0.;
double codelKB311 = 0.;
double codelKB312 = 0.;
double codelKB313 = 0.;
double codelKB322 = 0.;
double codelKB323 = 0.;
double codelKB333 = 0.;
double codelKC111 = 0.;
double codelKC112 = 0.;
double codelKC113 = 0.;
double codelKC122 = 0.;
double codelKC123 = 0.;
double codelKC133 = 0.;
double codelKC211 = 0.;
double codelKC212 = 0.;
double codelKC213 = 0.;
double codelKC222 = 0.;
double codelKC223 = 0.;
double codelKC233 = 0.;
double codelKC311 = 0.;
double codelKC312 = 0.;
double codelKC313 = 0.;
double codelKC322 = 0.;
double codelKC323 = 0.;
double codelKC333 = 0.;
double codelTrKA1 = 0.;
double codelTrKA2 = 0.;
double codelTrKA3 = 0.;
double codelTrKB1 = 0.;
double codelTrKB2 = 0.;
double codelTrKB3 = 0.;
double codelTrKC1 = 0.;
double codelTrKC2 = 0.;
double codelTrKC3 = 0.;
double ddRdr = 0.;
double deldelf11 = 0.;
double deldelf12 = 0.;
double deldelf13 = 0.;
double deldelf22 = 0.;
double deldelf23 = 0.;
double deldelf33 = 0.;
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
double delf1 = 0.;
double delf2 = 0.;
double delf3 = 0.;
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
double delK111 = 0.;
double delK112 = 0.;
double delK113 = 0.;
double delK122 = 0.;
double delK123 = 0.;
double delK133 = 0.;
double delK211 = 0.;
double delK212 = 0.;
double delK213 = 0.;
double delK222 = 0.;
double delK223 = 0.;
double delK233 = 0.;
double delK311 = 0.;
double delK312 = 0.;
double delK313 = 0.;
double delK322 = 0.;
double delK323 = 0.;
double delK333 = 0.;
double delKSST111 = 0.;
double delKSST112 = 0.;
double delKSST113 = 0.;
double delKSST121 = 0.;
double delKSST122 = 0.;
double delKSST123 = 0.;
double delKSST131 = 0.;
double delKSST132 = 0.;
double delKSST133 = 0.;
double delKSST211 = 0.;
double delKSST212 = 0.;
double delKSST213 = 0.;
double delKSST221 = 0.;
double delKSST222 = 0.;
double delKSST223 = 0.;
double delKSST231 = 0.;
double delKSST232 = 0.;
double delKSST233 = 0.;
double delKSST311 = 0.;
double delKSST312 = 0.;
double delKSST313 = 0.;
double delKSST321 = 0.;
double delKSST322 = 0.;
double delKSST323 = 0.;
double delKSST331 = 0.;
double delKSST332 = 0.;
double delKSST333 = 0.;
double denom = 0.;
double denom1 = 0.;
double denom2 = 0.;
double denom3 = 0.;
double detginvf = 0.;
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
double f = 0.;
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
double K = 0.;
double Kud11 = 0.;
double Kud12 = 0.;
double Kud13 = 0.;
double Kud21 = 0.;
double Kud22 = 0.;
double Kud23 = 0.;
double Kud31 = 0.;
double Kud32 = 0.;
double Kud33 = 0.;
double KudKud = 0.;
double momu1 = 0.;
double momu2 = 0.;
double momu3 = 0.;
double R = 0.;
double R11 = 0.;
double R12 = 0.;
double R13 = 0.;
double R22 = 0.;
double R23 = 0.;
double R33 = 0.;
double RA = 0.;
double RA11 = 0.;
double RA12 = 0.;
double RA13 = 0.;
double RA22 = 0.;
double RA23 = 0.;
double RA33 = 0.;
double RB = 0.;
double RB11 = 0.;
double RB12 = 0.;
double RB13 = 0.;
double RB22 = 0.;
double RB23 = 0.;
double RB33 = 0.;
double RC = 0.;
double RC11 = 0.;
double RC12 = 0.;
double RC13 = 0.;
double RC22 = 0.;
double RC23 = 0.;
double RC33 = 0.;
double RD = 0.;
double RD11 = 0.;
double RD12 = 0.;
double RD13 = 0.;
double RD22 = 0.;
double RD23 = 0.;
double RD33 = 0.;
double trK = 0.;



forinnerpoints_ijk_openmp(level) {



if (order_centered == 2 || boundary1away) { 

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

delK111
=
oo2dx*(-K11[-di + ijk] + K11[di + ijk])
;

delK112
=
oo2dx*(-K12[-di + ijk] + K12[di + ijk])
;

delK113
=
oo2dx*(-K13[-di + ijk] + K13[di + ijk])
;

delK122
=
oo2dx*(-K22[-di + ijk] + K22[di + ijk])
;

delK123
=
oo2dx*(-K23[-di + ijk] + K23[di + ijk])
;

delK133
=
oo2dx*(-K33[-di + ijk] + K33[di + ijk])
;

delK211
=
oo2dy*(-K11[-dj + ijk] + K11[dj + ijk])
;

delK212
=
oo2dy*(-K12[-dj + ijk] + K12[dj + ijk])
;

delK213
=
oo2dy*(-K13[-dj + ijk] + K13[dj + ijk])
;

delK222
=
oo2dy*(-K22[-dj + ijk] + K22[dj + ijk])
;

delK223
=
oo2dy*(-K23[-dj + ijk] + K23[dj + ijk])
;

delK233
=
oo2dy*(-K33[-dj + ijk] + K33[dj + ijk])
;

delK311
=
oo2dz*(-K11[-dk + ijk] + K11[dk + ijk])
;

delK312
=
oo2dz*(-K12[-dk + ijk] + K12[dk + ijk])
;

delK313
=
oo2dz*(-K13[-dk + ijk] + K13[dk + ijk])
;

delK322
=
oo2dz*(-K22[-dk + ijk] + K22[dk + ijk])
;

delK323
=
oo2dz*(-K23[-dk + ijk] + K23[dk + ijk])
;

delK333
=
oo2dz*(-K33[-dk + ijk] + K33[dk + ijk])
;


} else if (order_centered == 4 || boundaryNaway(2)) { 

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

delK111
=
0.16666666666666666667*oo2dx*(K11[-2*di + ijk] + 
    8.*(-K11[-di + ijk] + K11[di + ijk]) - K11[2*di + ijk])
;

delK112
=
0.16666666666666666667*oo2dx*(K12[-2*di + ijk] + 
    8.*(-K12[-di + ijk] + K12[di + ijk]) - K12[2*di + ijk])
;

delK113
=
0.16666666666666666667*oo2dx*(K13[-2*di + ijk] + 
    8.*(-K13[-di + ijk] + K13[di + ijk]) - K13[2*di + ijk])
;

delK122
=
0.16666666666666666667*oo2dx*(K22[-2*di + ijk] + 
    8.*(-K22[-di + ijk] + K22[di + ijk]) - K22[2*di + ijk])
;

delK123
=
0.16666666666666666667*oo2dx*(K23[-2*di + ijk] + 
    8.*(-K23[-di + ijk] + K23[di + ijk]) - K23[2*di + ijk])
;

delK133
=
0.16666666666666666667*oo2dx*(K33[-2*di + ijk] + 
    8.*(-K33[-di + ijk] + K33[di + ijk]) - K33[2*di + ijk])
;

delK211
=
0.16666666666666666667*oo2dy*(K11[-2*dj + ijk] + 
    8.*(-K11[-dj + ijk] + K11[dj + ijk]) - K11[2*dj + ijk])
;

delK212
=
0.16666666666666666667*oo2dy*(K12[-2*dj + ijk] + 
    8.*(-K12[-dj + ijk] + K12[dj + ijk]) - K12[2*dj + ijk])
;

delK213
=
0.16666666666666666667*oo2dy*(K13[-2*dj + ijk] + 
    8.*(-K13[-dj + ijk] + K13[dj + ijk]) - K13[2*dj + ijk])
;

delK222
=
0.16666666666666666667*oo2dy*(K22[-2*dj + ijk] + 
    8.*(-K22[-dj + ijk] + K22[dj + ijk]) - K22[2*dj + ijk])
;

delK223
=
0.16666666666666666667*oo2dy*(K23[-2*dj + ijk] + 
    8.*(-K23[-dj + ijk] + K23[dj + ijk]) - K23[2*dj + ijk])
;

delK233
=
0.16666666666666666667*oo2dy*(K33[-2*dj + ijk] + 
    8.*(-K33[-dj + ijk] + K33[dj + ijk]) - K33[2*dj + ijk])
;

delK311
=
0.16666666666666666667*oo2dz*(K11[-2*dk + ijk] + 
    8.*(-K11[-dk + ijk] + K11[dk + ijk]) - K11[2*dk + ijk])
;

delK312
=
0.16666666666666666667*oo2dz*(K12[-2*dk + ijk] + 
    8.*(-K12[-dk + ijk] + K12[dk + ijk]) - K12[2*dk + ijk])
;

delK313
=
0.16666666666666666667*oo2dz*(K13[-2*dk + ijk] + 
    8.*(-K13[-dk + ijk] + K13[dk + ijk]) - K13[2*dk + ijk])
;

delK322
=
0.16666666666666666667*oo2dz*(K22[-2*dk + ijk] + 
    8.*(-K22[-dk + ijk] + K22[dk + ijk]) - K22[2*dk + ijk])
;

delK323
=
0.16666666666666666667*oo2dz*(K23[-2*dk + ijk] + 
    8.*(-K23[-dk + ijk] + K23[dk + ijk]) - K23[2*dk + ijk])
;

delK333
=
0.16666666666666666667*oo2dz*(K33[-2*dk + ijk] + 
    8.*(-K33[-dk + ijk] + K33[dk + ijk]) - K33[2*dk + ijk])
;


} else if (order_centered == 6 || boundaryNaway(3)){ 

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

delK111
=
0.033333333333333333333*oo2dx*(-K11[-3*di + ijk] + 
    45.*(-K11[-di + ijk] + K11[di + ijk]) + 
    9.*(K11[-2*di + ijk] - K11[2*di + ijk]) + K11[3*di + ijk])
;

delK112
=
0.033333333333333333333*oo2dx*(-K12[-3*di + ijk] + 
    45.*(-K12[-di + ijk] + K12[di + ijk]) + 
    9.*(K12[-2*di + ijk] - K12[2*di + ijk]) + K12[3*di + ijk])
;

delK113
=
0.033333333333333333333*oo2dx*(-K13[-3*di + ijk] + 
    45.*(-K13[-di + ijk] + K13[di + ijk]) + 
    9.*(K13[-2*di + ijk] - K13[2*di + ijk]) + K13[3*di + ijk])
;

delK122
=
0.033333333333333333333*oo2dx*(-K22[-3*di + ijk] + 
    45.*(-K22[-di + ijk] + K22[di + ijk]) + 
    9.*(K22[-2*di + ijk] - K22[2*di + ijk]) + K22[3*di + ijk])
;

delK123
=
0.033333333333333333333*oo2dx*(-K23[-3*di + ijk] + 
    45.*(-K23[-di + ijk] + K23[di + ijk]) + 
    9.*(K23[-2*di + ijk] - K23[2*di + ijk]) + K23[3*di + ijk])
;

delK133
=
0.033333333333333333333*oo2dx*(-K33[-3*di + ijk] + 
    45.*(-K33[-di + ijk] + K33[di + ijk]) + 
    9.*(K33[-2*di + ijk] - K33[2*di + ijk]) + K33[3*di + ijk])
;

delK211
=
0.033333333333333333333*oo2dy*(-K11[-3*dj + ijk] + 
    45.*(-K11[-dj + ijk] + K11[dj + ijk]) + 
    9.*(K11[-2*dj + ijk] - K11[2*dj + ijk]) + K11[3*dj + ijk])
;

delK212
=
0.033333333333333333333*oo2dy*(-K12[-3*dj + ijk] + 
    45.*(-K12[-dj + ijk] + K12[dj + ijk]) + 
    9.*(K12[-2*dj + ijk] - K12[2*dj + ijk]) + K12[3*dj + ijk])
;

delK213
=
0.033333333333333333333*oo2dy*(-K13[-3*dj + ijk] + 
    45.*(-K13[-dj + ijk] + K13[dj + ijk]) + 
    9.*(K13[-2*dj + ijk] - K13[2*dj + ijk]) + K13[3*dj + ijk])
;

delK222
=
0.033333333333333333333*oo2dy*(-K22[-3*dj + ijk] + 
    45.*(-K22[-dj + ijk] + K22[dj + ijk]) + 
    9.*(K22[-2*dj + ijk] - K22[2*dj + ijk]) + K22[3*dj + ijk])
;

delK223
=
0.033333333333333333333*oo2dy*(-K23[-3*dj + ijk] + 
    45.*(-K23[-dj + ijk] + K23[dj + ijk]) + 
    9.*(K23[-2*dj + ijk] - K23[2*dj + ijk]) + K23[3*dj + ijk])
;

delK233
=
0.033333333333333333333*oo2dy*(-K33[-3*dj + ijk] + 
    45.*(-K33[-dj + ijk] + K33[dj + ijk]) + 
    9.*(K33[-2*dj + ijk] - K33[2*dj + ijk]) + K33[3*dj + ijk])
;

delK311
=
0.033333333333333333333*oo2dz*(-K11[-3*dk + ijk] + 
    45.*(-K11[-dk + ijk] + K11[dk + ijk]) + 
    9.*(K11[-2*dk + ijk] - K11[2*dk + ijk]) + K11[3*dk + ijk])
;

delK312
=
0.033333333333333333333*oo2dz*(-K12[-3*dk + ijk] + 
    45.*(-K12[-dk + ijk] + K12[dk + ijk]) + 
    9.*(K12[-2*dk + ijk] - K12[2*dk + ijk]) + K12[3*dk + ijk])
;

delK313
=
0.033333333333333333333*oo2dz*(-K13[-3*dk + ijk] + 
    45.*(-K13[-dk + ijk] + K13[dk + ijk]) + 
    9.*(K13[-2*dk + ijk] - K13[2*dk + ijk]) + K13[3*dk + ijk])
;

delK322
=
0.033333333333333333333*oo2dz*(-K22[-3*dk + ijk] + 
    45.*(-K22[-dk + ijk] + K22[dk + ijk]) + 
    9.*(K22[-2*dk + ijk] - K22[2*dk + ijk]) + K22[3*dk + ijk])
;

delK323
=
0.033333333333333333333*oo2dz*(-K23[-3*dk + ijk] + 
    45.*(-K23[-dk + ijk] + K23[dk + ijk]) + 
    9.*(K23[-2*dk + ijk] - K23[2*dk + ijk]) + K23[3*dk + ijk])
;

delK333
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

delKSST111
=
delK111*Jac11 + delK211*Jac21 + delK311*Jac31
;

delKSST112
=
delK112*Jac11 + delK212*Jac21 + delK312*Jac31
;

delKSST113
=
delK113*Jac11 + delK213*Jac21 + delK313*Jac31
;

delKSST121
=
delK112*Jac11 + delK212*Jac21 + delK312*Jac31
;

delKSST122
=
delK122*Jac11 + delK222*Jac21 + delK322*Jac31
;

delKSST123
=
delK123*Jac11 + delK223*Jac21 + delK323*Jac31
;

delKSST131
=
delK113*Jac11 + delK213*Jac21 + delK313*Jac31
;

delKSST132
=
delK123*Jac11 + delK223*Jac21 + delK323*Jac31
;

delKSST133
=
delK133*Jac11 + delK233*Jac21 + delK333*Jac31
;

delKSST211
=
delK111*Jac12 + delK211*Jac22 + delK311*Jac32
;

delKSST212
=
delK112*Jac12 + delK212*Jac22 + delK312*Jac32
;

delKSST213
=
delK113*Jac12 + delK213*Jac22 + delK313*Jac32
;

delKSST221
=
delK112*Jac12 + delK212*Jac22 + delK312*Jac32
;

delKSST222
=
delK122*Jac12 + delK222*Jac22 + delK322*Jac32
;

delKSST223
=
delK123*Jac12 + delK223*Jac22 + delK323*Jac32
;

delKSST231
=
delK113*Jac12 + delK213*Jac22 + delK313*Jac32
;

delKSST232
=
delK123*Jac12 + delK223*Jac22 + delK323*Jac32
;

delKSST233
=
delK133*Jac12 + delK233*Jac22 + delK333*Jac32
;

delKSST311
=
delK111*Jac13 + delK211*Jac23 + delK311*Jac33
;

delKSST312
=
delK112*Jac13 + delK212*Jac23 + delK312*Jac33
;

delKSST313
=
delK113*Jac13 + delK213*Jac23 + delK313*Jac33
;

delKSST321
=
delK112*Jac13 + delK212*Jac23 + delK312*Jac33
;

delKSST322
=
delK122*Jac13 + delK222*Jac23 + delK322*Jac33
;

delKSST323
=
delK123*Jac13 + delK223*Jac23 + delK323*Jac33
;

delKSST331
=
delK113*Jac13 + delK213*Jac23 + delK313*Jac33
;

delKSST332
=
delK123*Jac13 + delK223*Jac23 + delK323*Jac33
;

delKSST333
=
delK133*Jac13 + delK233*Jac23 + delK333*Jac33
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

delK111
=
delKSST111
;

delK112
=
delKSST121
;

delK113
=
delKSST131
;

delK122
=
delKSST122
;

delK123
=
delKSST132
;

delK133
=
delKSST133
;

delK211
=
delKSST211
;

delK212
=
delKSST221
;

delK213
=
delKSST231
;

delK222
=
delKSST222
;

delK223
=
delKSST232
;

delK233
=
delKSST233
;

delK311
=
delKSST311
;

delK312
=
delKSST321
;

delK313
=
delKSST331
;

delK322
=
delKSST322
;

delK323
=
delKSST332
;

delK333
=
delKSST333
;


} 

f
=
1.
;



/* conditional */
if (usepsi) {

f
=
pow4(psi[ijk])
;

delf1
=
4.*f*dpop1[ijk]
;

delf2
=
4.*f*dpop2[ijk]
;

delf3
=
4.*f*dpop3[ijk]
;

deldelf11
=
f*(4.*ddpop11[ijk] + 12.*pow2(dpop1[ijk]))
;

deldelf12
=
f*(4.*ddpop12[ijk] + 12.*dpop1[ijk]*dpop2[ijk])
;

deldelf13
=
f*(4.*ddpop13[ijk] + 12.*dpop1[ijk]*dpop3[ijk])
;

deldelf22
=
f*(4.*ddpop22[ijk] + 12.*pow2(dpop2[ijk]))
;

deldelf23
=
f*(4.*ddpop23[ijk] + 12.*dpop2[ijk]*dpop3[ijk])
;

deldelf33
=
f*(4.*ddpop33[ijk] + 12.*pow2(dpop3[ijk]))
;

deldelg1111
=
2.*delf1*delg111 + deldelg1111*f + deldelf11*g11[ijk]
;

deldelg1112
=
2.*delf1*delg112 + deldelg1112*f + deldelf11*g12[ijk]
;

deldelg1113
=
2.*delf1*delg113 + deldelg1113*f + deldelf11*g13[ijk]
;

deldelg1122
=
2.*delf1*delg122 + deldelg1122*f + deldelf11*g22[ijk]
;

deldelg1123
=
2.*delf1*delg123 + deldelg1123*f + deldelf11*g23[ijk]
;

deldelg1133
=
2.*delf1*delg133 + deldelg1133*f + deldelf11*g33[ijk]
;

deldelg1211
=
delf2*delg111 + delf1*delg211 + deldelg1211*f + deldelf12*g11[ijk]
;

deldelg1212
=
delf2*delg112 + delf1*delg212 + deldelg1212*f + deldelf12*g12[ijk]
;

deldelg1213
=
delf2*delg113 + delf1*delg213 + deldelg1213*f + deldelf12*g13[ijk]
;

deldelg1222
=
delf2*delg122 + delf1*delg222 + deldelg1222*f + deldelf12*g22[ijk]
;

deldelg1223
=
delf2*delg123 + delf1*delg223 + deldelg1223*f + deldelf12*g23[ijk]
;

deldelg1233
=
delf2*delg133 + delf1*delg233 + deldelg1233*f + deldelf12*g33[ijk]
;

deldelg1311
=
delf3*delg111 + delf1*delg311 + deldelg1311*f + deldelf13*g11[ijk]
;

deldelg1312
=
delf3*delg112 + delf1*delg312 + deldelg1312*f + deldelf13*g12[ijk]
;

deldelg1313
=
delf3*delg113 + delf1*delg313 + deldelg1313*f + deldelf13*g13[ijk]
;

deldelg1322
=
delf3*delg122 + delf1*delg322 + deldelg1322*f + deldelf13*g22[ijk]
;

deldelg1323
=
delf3*delg123 + delf1*delg323 + deldelg1323*f + deldelf13*g23[ijk]
;

deldelg1333
=
delf3*delg133 + delf1*delg333 + deldelg1333*f + deldelf13*g33[ijk]
;

deldelg2211
=
2.*delf2*delg211 + deldelg2211*f + deldelf22*g11[ijk]
;

deldelg2212
=
2.*delf2*delg212 + deldelg2212*f + deldelf22*g12[ijk]
;

deldelg2213
=
2.*delf2*delg213 + deldelg2213*f + deldelf22*g13[ijk]
;

deldelg2222
=
2.*delf2*delg222 + deldelg2222*f + deldelf22*g22[ijk]
;

deldelg2223
=
2.*delf2*delg223 + deldelg2223*f + deldelf22*g23[ijk]
;

deldelg2233
=
2.*delf2*delg233 + deldelg2233*f + deldelf22*g33[ijk]
;

deldelg2311
=
delf3*delg211 + delf2*delg311 + deldelg2311*f + deldelf23*g11[ijk]
;

deldelg2312
=
delf3*delg212 + delf2*delg312 + deldelg2312*f + deldelf23*g12[ijk]
;

deldelg2313
=
delf3*delg213 + delf2*delg313 + deldelg2313*f + deldelf23*g13[ijk]
;

deldelg2322
=
delf3*delg222 + delf2*delg322 + deldelg2322*f + deldelf23*g22[ijk]
;

deldelg2323
=
delf3*delg223 + delf2*delg323 + deldelg2323*f + deldelf23*g23[ijk]
;

deldelg2333
=
delf3*delg233 + delf2*delg333 + deldelg2333*f + deldelf23*g33[ijk]
;

deldelg3311
=
2.*delf3*delg311 + deldelg3311*f + deldelf33*g11[ijk]
;

deldelg3312
=
2.*delf3*delg312 + deldelg3312*f + deldelf33*g12[ijk]
;

deldelg3313
=
2.*delf3*delg313 + deldelg3313*f + deldelf33*g13[ijk]
;

deldelg3322
=
2.*delf3*delg322 + deldelg3322*f + deldelf33*g22[ijk]
;

deldelg3323
=
2.*delf3*delg323 + deldelg3323*f + deldelf33*g23[ijk]
;

deldelg3333
=
2.*delf3*delg333 + deldelg3333*f + deldelf33*g33[ijk]
;

delg111
=
delg111*f + delf1*g11[ijk]
;

delg112
=
delg112*f + delf1*g12[ijk]
;

delg113
=
delg113*f + delf1*g13[ijk]
;

delg122
=
delg122*f + delf1*g22[ijk]
;

delg123
=
delg123*f + delf1*g23[ijk]
;

delg133
=
delg133*f + delf1*g33[ijk]
;

delg211
=
delg211*f + delf2*g11[ijk]
;

delg212
=
delg212*f + delf2*g12[ijk]
;

delg213
=
delg213*f + delf2*g13[ijk]
;

delg222
=
delg222*f + delf2*g22[ijk]
;

delg223
=
delg223*f + delf2*g23[ijk]
;

delg233
=
delg233*f + delf2*g33[ijk]
;

delg311
=
delg311*f + delf3*g11[ijk]
;

delg312
=
delg312*f + delf3*g12[ijk]
;

delg313
=
delg313*f + delf3*g13[ijk]
;

delg322
=
delg322*f + delf3*g22[ijk]
;

delg323
=
delg323*f + delf3*g23[ijk]
;

delg333
=
delg333*f + delf3*g33[ijk]
;

}
/* if (usepsi) */


detginvf
=
1/(f*(2.*g12[ijk]*g13[ijk]*g23[ijk] + g11[ijk]*g22[ijk]*g33[ijk] - 
      g33[ijk]*pow2(g12[ijk]) - g22[ijk]*pow2(g13[ijk]) - 
      g11[ijk]*pow2(g23[ijk])))
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

R11
=
(gamma112*gammado111 - gamma111*gammado112 + gamma212*gammado211 - 
     gamma211*gammado212 + gamma312*gammado311 - gamma311*gammado312)*ginv12 \
+ (gamma113*gammado111 - gamma111*gammado113 + gamma213*gammado211 - 
     gamma211*gammado213 + gamma313*gammado311 - gamma311*gammado313)*ginv13 \
+ (deldelg1212 + 0.5*(-deldelg1122 - deldelg2211) + gamma112*gammado112 - 
     gamma111*gammado122 + gamma212*gammado212 - gamma211*gammado222 + 
     gamma312*gammado312 - gamma311*gammado322)*ginv22 + 
  (-deldelg1123 + deldelg1213 + deldelg1312 - deldelg2311 + 
     gamma113*gammado112 + gamma112*gammado113 + gamma213*gammado212 + 
     gamma212*gammado213 + gamma313*gammado312 + gamma312*gammado313 - 
     2.*(gamma111*gammado123 + gamma211*gammado223 + gamma311*gammado323))*
   ginv23 + (deldelg1313 + 0.5*(-deldelg1133 - deldelg3311) + 
     gamma113*gammado113 - gamma111*gammado133 + gamma213*gammado213 - 
     gamma211*gammado233 + gamma313*gammado313 - gamma311*gammado333)*ginv33
;

R12
=
(-(gamma112*gammado111) + gamma111*gammado112 - gamma212*gammado211 + 
     gamma211*gammado212 - gamma312*gammado311 + gamma311*gammado312)*ginv11 \
+ (-deldelg1212 + 0.5*(deldelg1122 + deldelg2211) - gamma112*gammado112 + 
     gamma111*gammado122 - gamma212*gammado212 + gamma211*gammado222 - 
     gamma312*gammado312 + gamma311*gammado322)*ginv12 + 
  (0.5*(deldelg1123 - deldelg1213 - deldelg1312 + deldelg2311) + 
     gamma113*gammado112 + gamma111*gammado123 + gamma213*gammado212 + 
     gamma211*gammado223 + gamma313*gammado312 - 
     2.*(gamma112*gammado113 + gamma212*gammado213 + 
        gamma312*gammado313) + gamma311*gammado323)*ginv13 + 
  (0.5*(-deldelg1223 + deldelg1322 + deldelg2213 - deldelg2312) + 
     gamma113*gammado122 - gamma112*gammado123 + gamma213*gammado222 - 
     gamma212*gammado223 + gamma313*gammado322 - gamma312*gammado323)*ginv23 \
+ (0.5*(-deldelg1233 + deldelg1323 + deldelg2313 - deldelg3312) + 
     gamma113*gammado123 - gamma112*gammado133 + gamma213*gammado223 - 
     gamma212*gammado233 + gamma313*gammado323 - gamma312*gammado333)*ginv33
;

R13
=
(0.5*(deldelg1123 - deldelg1213 - deldelg1312 + deldelg2311) + 
     gamma123*gammado111 - gamma113*gammado112 + gamma223*gammado211 - 
     gamma213*gammado212 + gamma323*gammado311 - gamma313*gammado312)*ginv12 \
+ (-deldelg1313 + 0.5*(deldelg1133 + deldelg3311) + gamma133*gammado111 - 
     gamma113*gammado113 + gamma233*gammado211 - gamma213*gammado213 + 
     gamma333*gammado311 - gamma313*gammado313)*ginv13 + 
  (0.5*(deldelg1223 - deldelg1322 - deldelg2213 + deldelg2312) + 
     gamma123*gammado112 - gamma113*gammado122 + gamma223*gammado212 - 
     gamma213*gammado222 + gamma323*gammado312 - gamma313*gammado322)*ginv22 \
+ (0.5*(deldelg1233 - deldelg1323 - deldelg2313 + deldelg3312) + 
     gamma133*gammado112 + gamma123*gammado113 + gamma233*gammado212 + 
     gamma223*gammado213 + gamma333*gammado312 + gamma323*gammado313 - 
     2.*(gamma113*gammado123 + gamma213*gammado223 + gamma313*gammado323))*
   ginv23 + (gamma133*gammado113 - gamma113*gammado133 + 
     gamma233*gammado213 - gamma213*gammado233 + gamma333*gammado313 - 
     gamma313*gammado333)*ginv33
;

R22
=
(deldelg1212 + 0.5*(-deldelg1122 - deldelg2211) - gamma122*gammado111 + 
     gamma112*gammado112 - gamma222*gammado211 + gamma212*gammado212 - 
     gamma322*gammado311 + gamma312*gammado312)*ginv11 + 
  (-(gamma122*gammado112) + gamma112*gammado122 - gamma222*gammado212 + 
     gamma212*gammado222 - gamma322*gammado312 + gamma312*gammado322)*ginv12 \
+ (deldelg1223 - deldelg1322 - deldelg2213 + deldelg2312 + 
     gamma123*gammado112 + gamma112*gammado123 + gamma223*gammado212 + 
     gamma212*gammado223 + gamma323*gammado312 - 
     2.*(gamma122*gammado113 + gamma222*gammado213 + 
        gamma322*gammado313) + gamma312*gammado323)*ginv13 + 
  (gamma123*gammado122 - gamma122*gammado123 + gamma223*gammado222 - 
     gamma222*gammado223 + gamma323*gammado322 - gamma322*gammado323)*ginv23 \
+ (deldelg2323 + 0.5*(-deldelg2233 - deldelg3322) + gamma123*gammado123 - 
     gamma122*gammado133 + gamma223*gammado223 - gamma222*gammado233 + 
     gamma323*gammado323 - gamma322*gammado333)*ginv33
;

R23
=
(0.5*(-deldelg1123 + deldelg1213 + deldelg1312 - deldelg2311) - 
     gamma123*gammado111 + gamma113*gammado112 - gamma223*gammado211 + 
     gamma213*gammado212 - gamma323*gammado311 + gamma313*gammado312)*ginv11 \
+ (0.5*(-deldelg1223 + deldelg1322 + deldelg2213 - deldelg2312) - 
     gamma123*gammado112 + gamma113*gammado122 - gamma223*gammado212 + 
     gamma213*gammado222 - gamma323*gammado312 + gamma313*gammado322)*ginv12 \
+ (0.5*(deldelg1233 - deldelg1323 - deldelg2313 + deldelg3312) + 
     gamma133*gammado112 + gamma113*gammado123 + gamma233*gammado212 + 
     gamma213*gammado223 + gamma333*gammado312 - 
     2.*(gamma123*gammado113 + gamma223*gammado213 + 
        gamma323*gammado313) + gamma313*gammado323)*ginv13 + 
  (-deldelg2323 + 0.5*(deldelg2233 + deldelg3322) + gamma133*gammado122 - 
     gamma123*gammado123 + gamma233*gammado222 - gamma223*gammado223 + 
     gamma333*gammado322 - gamma323*gammado323)*ginv23 + 
  (gamma133*gammado123 - gamma123*gammado133 + gamma233*gammado223 - 
     gamma223*gammado233 + gamma333*gammado323 - gamma323*gammado333)*ginv33
;

R33
=
(deldelg1313 + 0.5*(-deldelg1133 - deldelg3311) - gamma133*gammado111 + 
     gamma113*gammado113 - gamma233*gammado211 + gamma213*gammado213 - 
     gamma333*gammado311 + gamma313*gammado313)*ginv11 + 
  (-deldelg1233 + deldelg1323 + deldelg2313 - deldelg3312 + 
     gamma123*gammado113 + gamma113*gammado123 + gamma223*gammado213 + 
     gamma213*gammado223 - 2.*
      (gamma133*gammado112 + gamma233*gammado212 + gamma333*gammado312) + 
     gamma323*gammado313 + gamma313*gammado323)*ginv12 + 
  (-(gamma133*gammado113) + gamma113*gammado133 - gamma233*gammado213 + 
     gamma213*gammado233 - gamma333*gammado313 + gamma313*gammado333)*ginv13 \
+ (deldelg2323 + 0.5*(-deldelg2233 - deldelg3322) - gamma133*gammado122 + 
     gamma123*gammado123 - gamma233*gammado222 + gamma223*gammado223 - 
     gamma333*gammado322 + gamma323*gammado323)*ginv22 + 
  (-(gamma133*gammado123) + gamma123*gammado133 - gamma233*gammado223 + 
     gamma223*gammado233 - gamma333*gammado323 + gamma323*gammado333)*ginv23
;

R
=
ginv11*R11 + ginv22*R22 + 2.*(ginv12*R12 + ginv13*R13 + ginv23*R23) + 
  ginv33*R33
;

Kud11
=
ginv11*K11[ijk] + ginv12*K12[ijk] + ginv13*K13[ijk]
;

Kud12
=
ginv11*K12[ijk] + ginv12*K22[ijk] + ginv13*K23[ijk]
;

Kud13
=
ginv11*K13[ijk] + ginv12*K23[ijk] + ginv13*K33[ijk]
;

Kud21
=
ginv12*K11[ijk] + ginv22*K12[ijk] + ginv23*K13[ijk]
;

Kud22
=
ginv12*K12[ijk] + ginv22*K22[ijk] + ginv23*K23[ijk]
;

Kud23
=
ginv12*K13[ijk] + ginv22*K23[ijk] + ginv23*K33[ijk]
;

Kud31
=
ginv13*K11[ijk] + ginv23*K12[ijk] + ginv33*K13[ijk]
;

Kud32
=
ginv13*K12[ijk] + ginv23*K22[ijk] + ginv33*K23[ijk]
;

Kud33
=
ginv13*K13[ijk] + ginv23*K23[ijk] + ginv33*K33[ijk]
;

K
=
Kud11 + Kud22 + Kud33
;

KudKud
=
2.*(Kud12*Kud21 + Kud13*Kud31 + Kud23*Kud32) + pow2(Kud11) + pow2(Kud22) + 
  pow2(Kud33)
;

codelK111
=
delK111 - 2.*(gamma111*K11[ijk] + gamma211*K12[ijk] + gamma311*K13[ijk])
;

codelK112
=
delK112 - gamma112*K11[ijk] - (gamma111 + gamma212)*K12[ijk] - 
  gamma312*K13[ijk] - gamma211*K22[ijk] - gamma311*K23[ijk]
;

codelK113
=
delK113 - gamma113*K11[ijk] - gamma213*K12[ijk] - 
  (gamma111 + gamma313)*K13[ijk] - gamma211*K23[ijk] - gamma311*K33[ijk]
;

codelK122
=
delK122 - 2.*(gamma112*K12[ijk] + gamma212*K22[ijk] + gamma312*K23[ijk])
;

codelK123
=
delK123 - gamma113*K12[ijk] - gamma112*K13[ijk] - gamma213*K22[ijk] - 
  (gamma212 + gamma313)*K23[ijk] - gamma312*K33[ijk]
;

codelK133
=
delK133 - 2.*(gamma113*K13[ijk] + gamma213*K23[ijk] + gamma313*K33[ijk])
;

codelK211
=
delK211 - 2.*(gamma112*K11[ijk] + gamma212*K12[ijk] + gamma312*K13[ijk])
;

codelK212
=
delK212 - gamma122*K11[ijk] - (gamma112 + gamma222)*K12[ijk] - 
  gamma322*K13[ijk] - gamma212*K22[ijk] - gamma312*K23[ijk]
;

codelK213
=
delK213 - gamma123*K11[ijk] - gamma223*K12[ijk] - 
  (gamma112 + gamma323)*K13[ijk] - gamma212*K23[ijk] - gamma312*K33[ijk]
;

codelK222
=
delK222 - 2.*(gamma122*K12[ijk] + gamma222*K22[ijk] + gamma322*K23[ijk])
;

codelK223
=
delK223 - gamma123*K12[ijk] - gamma122*K13[ijk] - gamma223*K22[ijk] - 
  (gamma222 + gamma323)*K23[ijk] - gamma322*K33[ijk]
;

codelK233
=
delK233 - 2.*(gamma123*K13[ijk] + gamma223*K23[ijk] + gamma323*K33[ijk])
;

codelK311
=
delK311 - 2.*(gamma113*K11[ijk] + gamma213*K12[ijk] + gamma313*K13[ijk])
;

codelK312
=
delK312 - gamma123*K11[ijk] - (gamma113 + gamma223)*K12[ijk] - 
  gamma323*K13[ijk] - gamma213*K22[ijk] - gamma313*K23[ijk]
;

codelK313
=
delK313 - gamma133*K11[ijk] - gamma233*K12[ijk] - 
  (gamma113 + gamma333)*K13[ijk] - gamma213*K23[ijk] - gamma313*K33[ijk]
;

codelK322
=
delK322 - 2.*(gamma123*K12[ijk] + gamma223*K22[ijk] + gamma323*K23[ijk])
;

codelK323
=
delK323 - gamma133*K12[ijk] - gamma123*K13[ijk] - gamma233*K22[ijk] - 
  (gamma223 + gamma333)*K23[ijk] - gamma323*K33[ijk]
;

codelK333
=
delK333 - 2.*(gamma133*K13[ijk] + gamma233*K23[ijk] + gamma333*K33[ijk])
;

cdKudd111
=
codelK111*ginv11 + codelK211*ginv12 + codelK311*ginv13
;

cdKudd112
=
codelK112*ginv11 + codelK212*ginv12 + codelK312*ginv13
;

cdKudd113
=
codelK113*ginv11 + codelK213*ginv12 + codelK313*ginv13
;

cdKudd122
=
codelK122*ginv11 + codelK222*ginv12 + codelK322*ginv13
;

cdKudd123
=
codelK123*ginv11 + codelK223*ginv12 + codelK323*ginv13
;

cdKudd133
=
codelK133*ginv11 + codelK233*ginv12 + codelK333*ginv13
;

cdKudd211
=
codelK111*ginv12 + codelK211*ginv22 + codelK311*ginv23
;

cdKudd212
=
codelK112*ginv12 + codelK212*ginv22 + codelK312*ginv23
;

cdKudd213
=
codelK113*ginv12 + codelK213*ginv22 + codelK313*ginv23
;

cdKudd222
=
codelK122*ginv12 + codelK222*ginv22 + codelK322*ginv23
;

cdKudd223
=
codelK123*ginv12 + codelK223*ginv22 + codelK323*ginv23
;

cdKudd233
=
codelK133*ginv12 + codelK233*ginv22 + codelK333*ginv23
;

cdKudd311
=
codelK111*ginv13 + codelK211*ginv23 + codelK311*ginv33
;

cdKudd312
=
codelK112*ginv13 + codelK212*ginv23 + codelK312*ginv33
;

cdKudd313
=
codelK113*ginv13 + codelK213*ginv23 + codelK313*ginv33
;

cdKudd322
=
codelK122*ginv13 + codelK222*ginv23 + codelK322*ginv33
;

cdKudd323
=
codelK123*ginv13 + codelK223*ginv23 + codelK323*ginv33
;

cdKudd333
=
codelK133*ginv13 + codelK233*ginv23 + codelK333*ginv33
;



/* conditional */
if (matter) {

ham[ijk]
=
-KudKud + R - 16.*PI*rho[ijk] + pow2(K)
;


} else { /* if (!matter) */

ham[ijk]
=
-KudKud + R + pow2(K)
;

}
/* if (matter) */


momu1
=
(cdKudd212 + cdKudd313)*ginv11 + 
  (-cdKudd112 + cdKudd222 + cdKudd323)*ginv12 + 
  (-cdKudd113 + cdKudd223 + cdKudd333)*ginv13 - cdKudd122*ginv22 - 
  2.*cdKudd123*ginv23 - cdKudd133*ginv33
;

momu2
=
-(cdKudd211*ginv11) + (cdKudd111 - cdKudd212 + cdKudd313)*ginv12 - 
  2.*cdKudd213*ginv13 + (cdKudd112 + cdKudd323)*ginv22 + 
  (cdKudd113 - cdKudd223 + cdKudd333)*ginv23 - cdKudd233*ginv33
;

momu3
=
-(cdKudd311*ginv11) - 2.*cdKudd312*ginv12 + 
  (cdKudd111 + cdKudd212 - cdKudd313)*ginv13 - cdKudd322*ginv22 + 
  (cdKudd112 + cdKudd222 - cdKudd323)*ginv23 + (cdKudd113 + cdKudd223)*ginv33
;



/* conditional */
if (matter) {

mom1[ijk]
=
f*(g11[ijk]*(momu1 - 8.*PI*S1[ijk]) + g12[ijk]*(momu2 - 8.*PI*S2[ijk]) + 
    g13[ijk]*(momu3 - 8.*PI*S3[ijk]))
;

mom2[ijk]
=
f*(g12[ijk]*(momu1 - 8.*PI*S1[ijk]) + g22[ijk]*(momu2 - 8.*PI*S2[ijk]) + 
    g23[ijk]*(momu3 - 8.*PI*S3[ijk]))
;

mom3[ijk]
=
f*(g13[ijk]*(momu1 - 8.*PI*S1[ijk]) + g23[ijk]*(momu2 - 8.*PI*S2[ijk]) + 
    g33[ijk]*(momu3 - 8.*PI*S3[ijk]))
;


} else { /* if (!matter) */

mom1[ijk]
=
f*(momu1*g11[ijk] + momu2*g12[ijk] + momu3*g13[ijk])
;

mom2[ijk]
=
f*(momu1*g12[ijk] + momu2*g22[ijk] + momu3*g23[ijk])
;

mom3[ijk]
=
f*(momu1*g13[ijk] + momu2*g23[ijk] + momu3*g33[ijk])
;

}
/* if (matter) */


trK
=
K
;



/* conditional */
if (normConstr) {



/* conditional */
if (TermByTerm) {

RA11
=
-(deldelg1111*ginv11) - deldelg2211*ginv22 - 
  2.*(deldelg1211*ginv12 + deldelg1311*ginv13 + deldelg2311*ginv23) - 
  deldelg3311*ginv33
;

RA12
=
-(deldelg1112*ginv11) - deldelg2212*ginv22 - 
  2.*(deldelg1212*ginv12 + deldelg1312*ginv13 + deldelg2312*ginv23) - 
  deldelg3312*ginv33
;

RA13
=
-(deldelg1113*ginv11) - deldelg2213*ginv22 - 
  2.*(deldelg1213*ginv12 + deldelg1313*ginv13 + deldelg2313*ginv23) - 
  deldelg3313*ginv33
;

RA22
=
-(deldelg1122*ginv11) - deldelg2222*ginv22 - 
  2.*(deldelg1222*ginv12 + deldelg1322*ginv13 + deldelg2322*ginv23) - 
  deldelg3322*ginv33
;

RA23
=
-(deldelg1123*ginv11) - deldelg2223*ginv22 - 
  2.*(deldelg1223*ginv12 + deldelg1323*ginv13 + deldelg2323*ginv23) - 
  deldelg3323*ginv33
;

RA33
=
-(deldelg1133*ginv11) - deldelg2233*ginv22 - 
  2.*(deldelg1233*ginv12 + deldelg1333*ginv13 + deldelg2333*ginv23) - 
  deldelg3333*ginv33
;

RB11
=
deldelg1111*ginv11 + (deldelg1112 + deldelg1211)*ginv12 + 
  (deldelg1113 + deldelg1311)*ginv13 + deldelg1212*ginv22 + 
  (deldelg1213 + deldelg1312)*ginv23 + deldelg1313*ginv33
;

RB12
=
deldelg1211*ginv11 + (deldelg1212 + deldelg2211)*ginv12 + 
  (deldelg1213 + deldelg2311)*ginv13 + deldelg2212*ginv22 + 
  (deldelg2213 + deldelg2312)*ginv23 + deldelg2313*ginv33
;

RB13
=
deldelg1311*ginv11 + (deldelg1312 + deldelg2311)*ginv12 + 
  (deldelg1313 + deldelg3311)*ginv13 + deldelg2312*ginv22 + 
  (deldelg2313 + deldelg3312)*ginv23 + deldelg3313*ginv33
;

RB22
=
deldelg1212*ginv11 + (deldelg1222 + deldelg2212)*ginv12 + 
  (deldelg1223 + deldelg2312)*ginv13 + deldelg2222*ginv22 + 
  (deldelg2223 + deldelg2322)*ginv23 + deldelg2323*ginv33
;

RB23
=
deldelg1312*ginv11 + (deldelg1322 + deldelg2312)*ginv12 + 
  (deldelg1323 + deldelg3312)*ginv13 + deldelg2322*ginv22 + 
  (deldelg2323 + deldelg3322)*ginv23 + deldelg3323*ginv33
;

RB33
=
deldelg1313*ginv11 + (deldelg1323 + deldelg2313)*ginv12 + 
  (deldelg1333 + deldelg3313)*ginv13 + deldelg2323*ginv22 + 
  (deldelg2333 + deldelg3323)*ginv23 + deldelg3333*ginv33
;

RC11
=
(gamma111*gammado111 + gamma211*gammado211 + gamma311*gammado311)*ginv11 + 
  (gamma112*gammado111 + gamma111*gammado112 + gamma212*gammado211 + 
     gamma211*gammado212 + gamma312*gammado311 + gamma311*gammado312)*ginv12 \
+ (gamma113*gammado111 + gamma111*gammado113 + gamma213*gammado211 + 
     gamma211*gammado213 + gamma313*gammado311 + gamma311*gammado313)*ginv13 \
+ (gamma112*gammado112 + gamma212*gammado212 + gamma312*gammado312)*ginv22 + 
  (gamma113*gammado112 + gamma112*gammado113 + gamma213*gammado212 + 
     gamma212*gammado213 + gamma313*gammado312 + gamma312*gammado313)*ginv23 \
+ (gamma113*gammado113 + gamma213*gammado213 + gamma313*gammado313)*ginv33
;

RC12
=
(gamma111*gammado112 + gamma211*gammado212 + gamma311*gammado312)*ginv11 + 
  (gamma112*gammado112 + gamma111*gammado122 + gamma212*gammado212 + 
     gamma211*gammado222 + gamma312*gammado312 + gamma311*gammado322)*ginv12 \
+ (gamma113*gammado112 + gamma111*gammado123 + gamma213*gammado212 + 
     gamma211*gammado223 + gamma313*gammado312 + gamma311*gammado323)*ginv13 \
+ (gamma112*gammado122 + gamma212*gammado222 + gamma312*gammado322)*ginv22 + 
  (gamma113*gammado122 + gamma112*gammado123 + gamma213*gammado222 + 
     gamma212*gammado223 + gamma313*gammado322 + gamma312*gammado323)*ginv23 \
+ (gamma113*gammado123 + gamma213*gammado223 + gamma313*gammado323)*ginv33
;

RC13
=
(gamma111*gammado113 + gamma211*gammado213 + gamma311*gammado313)*ginv11 + 
  (gamma112*gammado113 + gamma111*gammado123 + gamma212*gammado213 + 
     gamma211*gammado223 + gamma312*gammado313 + gamma311*gammado323)*ginv12 \
+ (gamma113*gammado113 + gamma111*gammado133 + gamma213*gammado213 + 
     gamma211*gammado233 + gamma313*gammado313 + gamma311*gammado333)*ginv13 \
+ (gamma112*gammado123 + gamma212*gammado223 + gamma312*gammado323)*ginv22 + 
  (gamma113*gammado123 + gamma112*gammado133 + gamma213*gammado223 + 
     gamma212*gammado233 + gamma313*gammado323 + gamma312*gammado333)*ginv23 \
+ (gamma113*gammado133 + gamma213*gammado233 + gamma313*gammado333)*ginv33
;

RC22
=
(gamma112*gammado112 + gamma212*gammado212 + gamma312*gammado312)*ginv11 + 
  (gamma122*gammado112 + gamma112*gammado122 + gamma222*gammado212 + 
     gamma212*gammado222 + gamma322*gammado312 + gamma312*gammado322)*ginv12 \
+ (gamma123*gammado112 + gamma112*gammado123 + gamma223*gammado212 + 
     gamma212*gammado223 + gamma323*gammado312 + gamma312*gammado323)*ginv13 \
+ (gamma122*gammado122 + gamma222*gammado222 + gamma322*gammado322)*ginv22 + 
  (gamma123*gammado122 + gamma122*gammado123 + gamma223*gammado222 + 
     gamma222*gammado223 + gamma323*gammado322 + gamma322*gammado323)*ginv23 \
+ (gamma123*gammado123 + gamma223*gammado223 + gamma323*gammado323)*ginv33
;

RC23
=
(gamma112*gammado113 + gamma212*gammado213 + gamma312*gammado313)*ginv11 + 
  (gamma122*gammado113 + gamma112*gammado123 + gamma222*gammado213 + 
     gamma212*gammado223 + gamma322*gammado313 + gamma312*gammado323)*ginv12 \
+ (gamma123*gammado113 + gamma112*gammado133 + gamma223*gammado213 + 
     gamma212*gammado233 + gamma323*gammado313 + gamma312*gammado333)*ginv13 \
+ (gamma122*gammado123 + gamma222*gammado223 + gamma322*gammado323)*ginv22 + 
  (gamma123*gammado123 + gamma122*gammado133 + gamma223*gammado223 + 
     gamma222*gammado233 + gamma323*gammado323 + gamma322*gammado333)*ginv23 \
+ (gamma123*gammado133 + gamma223*gammado233 + gamma323*gammado333)*ginv33
;

RC33
=
(gamma113*gammado113 + gamma213*gammado213 + gamma313*gammado313)*ginv11 + 
  (gamma123*gammado113 + gamma113*gammado123 + gamma223*gammado213 + 
     gamma213*gammado223 + gamma323*gammado313 + gamma313*gammado323)*ginv12 \
+ (gamma133*gammado113 + gamma113*gammado133 + gamma233*gammado213 + 
     gamma213*gammado233 + gamma333*gammado313 + gamma313*gammado333)*ginv13 \
+ (gamma123*gammado123 + gamma223*gammado223 + gamma323*gammado323)*ginv22 + 
  (gamma133*gammado123 + gamma123*gammado133 + gamma233*gammado223 + 
     gamma223*gammado233 + gamma333*gammado323 + gamma323*gammado333)*ginv23 \
+ (gamma133*gammado133 + gamma233*gammado233 + gamma333*gammado333)*ginv33
;

RD11
=
-((gamma111*gammado111 + gamma211*gammado211 + gamma311*gammado311)*
     ginv11) - (gamma111*gammado122 + gamma211*gammado222 + 
     gamma311*gammado322)*ginv22 - 
  2.*((gamma111*gammado112 + gamma211*gammado212 + gamma311*gammado312)*
      ginv12 + (gamma111*gammado113 + gamma211*gammado213 + 
        gamma311*gammado313)*ginv13 + 
     (gamma111*gammado123 + gamma211*gammado223 + gamma311*gammado323)*
      ginv23) - (gamma111*gammado133 + gamma211*gammado233 + 
     gamma311*gammado333)*ginv33
;

RD12
=
-((gamma112*gammado111 + gamma212*gammado211 + gamma312*gammado311)*
     ginv11) - (gamma112*gammado122 + gamma212*gammado222 + 
     gamma312*gammado322)*ginv22 - 
  2.*((gamma112*gammado112 + gamma212*gammado212 + gamma312*gammado312)*
      ginv12 + (gamma112*gammado113 + gamma212*gammado213 + 
        gamma312*gammado313)*ginv13 + 
     (gamma112*gammado123 + gamma212*gammado223 + gamma312*gammado323)*
      ginv23) - (gamma112*gammado133 + gamma212*gammado233 + 
     gamma312*gammado333)*ginv33
;

RD13
=
-((gamma113*gammado111 + gamma213*gammado211 + gamma313*gammado311)*
     ginv11) - (gamma113*gammado122 + gamma213*gammado222 + 
     gamma313*gammado322)*ginv22 - 
  2.*((gamma113*gammado112 + gamma213*gammado212 + gamma313*gammado312)*
      ginv12 + (gamma113*gammado113 + gamma213*gammado213 + 
        gamma313*gammado313)*ginv13 + 
     (gamma113*gammado123 + gamma213*gammado223 + gamma313*gammado323)*
      ginv23) - (gamma113*gammado133 + gamma213*gammado233 + 
     gamma313*gammado333)*ginv33
;

RD22
=
-((gamma122*gammado111 + gamma222*gammado211 + gamma322*gammado311)*
     ginv11) - (gamma122*gammado122 + gamma222*gammado222 + 
     gamma322*gammado322)*ginv22 - 
  2.*((gamma122*gammado112 + gamma222*gammado212 + gamma322*gammado312)*
      ginv12 + (gamma122*gammado113 + gamma222*gammado213 + 
        gamma322*gammado313)*ginv13 + 
     (gamma122*gammado123 + gamma222*gammado223 + gamma322*gammado323)*
      ginv23) - (gamma122*gammado133 + gamma222*gammado233 + 
     gamma322*gammado333)*ginv33
;

RD23
=
-((gamma123*gammado111 + gamma223*gammado211 + gamma323*gammado311)*
     ginv11) - (gamma123*gammado122 + gamma223*gammado222 + 
     gamma323*gammado322)*ginv22 - 
  2.*((gamma123*gammado112 + gamma223*gammado212 + gamma323*gammado312)*
      ginv12 + (gamma123*gammado113 + gamma223*gammado213 + 
        gamma323*gammado313)*ginv13 + 
     (gamma123*gammado123 + gamma223*gammado223 + gamma323*gammado323)*
      ginv23) - (gamma123*gammado133 + gamma223*gammado233 + 
     gamma323*gammado333)*ginv33
;

RD33
=
-((gamma133*gammado111 + gamma233*gammado211 + gamma333*gammado311)*
     ginv11) - (gamma133*gammado122 + gamma233*gammado222 + 
     gamma333*gammado322)*ginv22 - 
  2.*((gamma133*gammado112 + gamma233*gammado212 + gamma333*gammado312)*
      ginv12 + (gamma133*gammado113 + gamma233*gammado213 + 
        gamma333*gammado313)*ginv13 + 
     (gamma133*gammado123 + gamma233*gammado223 + gamma333*gammado323)*
      ginv23) - (gamma133*gammado133 + gamma233*gammado233 + 
     gamma333*gammado333)*ginv33
;

RA
=
ginv11*RA11 + ginv22*RA22 + 2.*(ginv12*RA12 + ginv13*RA13 + ginv23*RA23) + 
  ginv33*RA33
;

RB
=
ginv11*RB11 + ginv22*RB22 + 2.*(ginv12*RB12 + ginv13*RB13 + ginv23*RB23) + 
  ginv33*RB33
;

RC
=
ginv11*RC11 + ginv22*RC22 + 2.*(ginv12*RC12 + ginv13*RC13 + ginv23*RC23) + 
  ginv33*RC33
;

RD
=
ginv11*RD11 + ginv22*RD22 + 2.*(ginv12*RD12 + ginv13*RD13 + ginv23*RD23) + 
  ginv33*RD33
;

denom
=
fabs(KudKud) + fabs(RA) + fabs(RB) + fabs(RC) + fabs(RD) + pow2(K)
;


} else { /* if (!TermByTerm) */

denom
=
fabs(KudKud) + fabs(R) + pow2(K)
;

}
/* if (TermByTerm) */


normham[ijk]
=
Cal(denom <= 0.,0.,ham[ijk]/denom)
;



/* conditional */
if (TermByTerm) {

codelKA111
=
delK111
;

codelKA112
=
delK112
;

codelKA113
=
delK113
;

codelKA122
=
delK122
;

codelKA123
=
delK123
;

codelKA133
=
delK133
;

codelKA211
=
delK211
;

codelKA212
=
delK212
;

codelKA213
=
delK213
;

codelKA222
=
delK222
;

codelKA223
=
delK223
;

codelKA233
=
delK233
;

codelKA311
=
delK311
;

codelKA312
=
delK312
;

codelKA313
=
delK313
;

codelKA322
=
delK322
;

codelKA323
=
delK323
;

codelKA333
=
delK333
;

codelKB111
=
-(gamma111*K11[ijk]) - gamma211*K12[ijk] - gamma311*K13[ijk]
;

codelKB112
=
-(gamma111*K12[ijk]) - gamma211*K22[ijk] - gamma311*K23[ijk]
;

codelKB113
=
-(gamma111*K13[ijk]) - gamma211*K23[ijk] - gamma311*K33[ijk]
;

codelKB122
=
-(gamma112*K12[ijk]) - gamma212*K22[ijk] - gamma312*K23[ijk]
;

codelKB123
=
-(gamma112*K13[ijk]) - gamma212*K23[ijk] - gamma312*K33[ijk]
;

codelKB133
=
-(gamma113*K13[ijk]) - gamma213*K23[ijk] - gamma313*K33[ijk]
;

codelKB211
=
-(gamma112*K11[ijk]) - gamma212*K12[ijk] - gamma312*K13[ijk]
;

codelKB212
=
-(gamma112*K12[ijk]) - gamma212*K22[ijk] - gamma312*K23[ijk]
;

codelKB213
=
-(gamma112*K13[ijk]) - gamma212*K23[ijk] - gamma312*K33[ijk]
;

codelKB222
=
-(gamma122*K12[ijk]) - gamma222*K22[ijk] - gamma322*K23[ijk]
;

codelKB223
=
-(gamma122*K13[ijk]) - gamma222*K23[ijk] - gamma322*K33[ijk]
;

codelKB233
=
-(gamma123*K13[ijk]) - gamma223*K23[ijk] - gamma323*K33[ijk]
;

codelKB311
=
-(gamma113*K11[ijk]) - gamma213*K12[ijk] - gamma313*K13[ijk]
;

codelKB312
=
-(gamma113*K12[ijk]) - gamma213*K22[ijk] - gamma313*K23[ijk]
;

codelKB313
=
-(gamma113*K13[ijk]) - gamma213*K23[ijk] - gamma313*K33[ijk]
;

codelKB322
=
-(gamma123*K12[ijk]) - gamma223*K22[ijk] - gamma323*K23[ijk]
;

codelKB323
=
-(gamma123*K13[ijk]) - gamma223*K23[ijk] - gamma323*K33[ijk]
;

codelKB333
=
-(gamma133*K13[ijk]) - gamma233*K23[ijk] - gamma333*K33[ijk]
;

codelKC111
=
-(gamma111*K11[ijk]) - gamma211*K12[ijk] - gamma311*K13[ijk]
;

codelKC112
=
-(gamma111*K12[ijk]) - gamma211*K22[ijk] - gamma311*K23[ijk]
;

codelKC113
=
-(gamma111*K13[ijk]) - gamma211*K23[ijk] - gamma311*K33[ijk]
;

codelKC122
=
-(gamma112*K12[ijk]) - gamma212*K22[ijk] - gamma312*K23[ijk]
;

codelKC123
=
-(gamma112*K13[ijk]) - gamma212*K23[ijk] - gamma312*K33[ijk]
;

codelKC133
=
-(gamma113*K13[ijk]) - gamma213*K23[ijk] - gamma313*K33[ijk]
;

codelKC211
=
-(gamma112*K11[ijk]) - gamma212*K12[ijk] - gamma312*K13[ijk]
;

codelKC212
=
-(gamma112*K12[ijk]) - gamma212*K22[ijk] - gamma312*K23[ijk]
;

codelKC213
=
-(gamma112*K13[ijk]) - gamma212*K23[ijk] - gamma312*K33[ijk]
;

codelKC222
=
-(gamma122*K12[ijk]) - gamma222*K22[ijk] - gamma322*K23[ijk]
;

codelKC223
=
-(gamma122*K13[ijk]) - gamma222*K23[ijk] - gamma322*K33[ijk]
;

codelKC233
=
-(gamma123*K13[ijk]) - gamma223*K23[ijk] - gamma323*K33[ijk]
;

codelKC311
=
-(gamma113*K11[ijk]) - gamma213*K12[ijk] - gamma313*K13[ijk]
;

codelKC312
=
-(gamma113*K12[ijk]) - gamma213*K22[ijk] - gamma313*K23[ijk]
;

codelKC313
=
-(gamma113*K13[ijk]) - gamma213*K23[ijk] - gamma313*K33[ijk]
;

codelKC322
=
-(gamma123*K12[ijk]) - gamma223*K22[ijk] - gamma323*K23[ijk]
;

codelKC323
=
-(gamma123*K13[ijk]) - gamma223*K23[ijk] - gamma323*K33[ijk]
;

codelKC333
=
-(gamma133*K13[ijk]) - gamma233*K23[ijk] - gamma333*K33[ijk]
;

cdKdA1
=
codelKA111*ginv11 + (codelKA112 + codelKA211)*ginv12 + 
  (codelKA113 + codelKA311)*ginv13 + codelKA212*ginv22 + 
  (codelKA213 + codelKA312)*ginv23 + codelKA313*ginv33
;

cdKdA2
=
codelKA112*ginv11 + (codelKA122 + codelKA212)*ginv12 + 
  (codelKA123 + codelKA312)*ginv13 + codelKA222*ginv22 + 
  (codelKA223 + codelKA322)*ginv23 + codelKA323*ginv33
;

cdKdA3
=
codelKA113*ginv11 + (codelKA123 + codelKA213)*ginv12 + 
  (codelKA133 + codelKA313)*ginv13 + codelKA223*ginv22 + 
  (codelKA233 + codelKA323)*ginv23 + codelKA333*ginv33
;

cdKdB1
=
codelKB111*ginv11 + (codelKB112 + codelKB211)*ginv12 + 
  (codelKB113 + codelKB311)*ginv13 + codelKB212*ginv22 + 
  (codelKB213 + codelKB312)*ginv23 + codelKB313*ginv33
;

cdKdB2
=
codelKB112*ginv11 + (codelKB122 + codelKB212)*ginv12 + 
  (codelKB123 + codelKB312)*ginv13 + codelKB222*ginv22 + 
  (codelKB223 + codelKB322)*ginv23 + codelKB323*ginv33
;

cdKdB3
=
codelKB113*ginv11 + (codelKB123 + codelKB213)*ginv12 + 
  (codelKB133 + codelKB313)*ginv13 + codelKB223*ginv22 + 
  (codelKB233 + codelKB323)*ginv23 + codelKB333*ginv33
;

cdKdC1
=
codelKC111*ginv11 + (codelKC112 + codelKC211)*ginv12 + 
  (codelKC113 + codelKC311)*ginv13 + codelKC212*ginv22 + 
  (codelKC213 + codelKC312)*ginv23 + codelKC313*ginv33
;

cdKdC2
=
codelKC112*ginv11 + (codelKC122 + codelKC212)*ginv12 + 
  (codelKC123 + codelKC312)*ginv13 + codelKC222*ginv22 + 
  (codelKC223 + codelKC322)*ginv23 + codelKC323*ginv33
;

cdKdC3
=
codelKC113*ginv11 + (codelKC123 + codelKC213)*ginv12 + 
  (codelKC133 + codelKC313)*ginv13 + codelKC223*ginv22 + 
  (codelKC233 + codelKC323)*ginv23 + codelKC333*ginv33
;

codelTrKA1
=
codelKA111*ginv11 + codelKA122*ginv22 + 
  2.*(codelKA112*ginv12 + codelKA113*ginv13 + codelKA123*ginv23) + 
  codelKA133*ginv33
;

codelTrKA2
=
codelKA211*ginv11 + codelKA222*ginv22 + 
  2.*(codelKA212*ginv12 + codelKA213*ginv13 + codelKA223*ginv23) + 
  codelKA233*ginv33
;

codelTrKA3
=
codelKA311*ginv11 + codelKA322*ginv22 + 
  2.*(codelKA312*ginv12 + codelKA313*ginv13 + codelKA323*ginv23) + 
  codelKA333*ginv33
;

codelTrKB1
=
codelKB111*ginv11 + codelKB122*ginv22 + 
  2.*(codelKB112*ginv12 + codelKB113*ginv13 + codelKB123*ginv23) + 
  codelKB133*ginv33
;

codelTrKB2
=
codelKB211*ginv11 + codelKB222*ginv22 + 
  2.*(codelKB212*ginv12 + codelKB213*ginv13 + codelKB223*ginv23) + 
  codelKB233*ginv33
;

codelTrKB3
=
codelKB311*ginv11 + codelKB322*ginv22 + 
  2.*(codelKB312*ginv12 + codelKB313*ginv13 + codelKB323*ginv23) + 
  codelKB333*ginv33
;

codelTrKC1
=
codelKC111*ginv11 + codelKC122*ginv22 + 
  2.*(codelKC112*ginv12 + codelKC113*ginv13 + codelKC123*ginv23) + 
  codelKC133*ginv33
;

codelTrKC2
=
codelKC211*ginv11 + codelKC222*ginv22 + 
  2.*(codelKC212*ginv12 + codelKC213*ginv13 + codelKC223*ginv23) + 
  codelKC233*ginv33
;

codelTrKC3
=
codelKC311*ginv11 + codelKC322*ginv22 + 
  2.*(codelKC312*ginv12 + codelKC313*ginv13 + codelKC323*ginv23) + 
  codelKC333*ginv33
;

denom1
=
fabs(cdKdA1) + fabs(cdKdB1) + fabs(cdKdC1) + fabs(codelTrKA1) + 
  fabs(codelTrKB1) + fabs(codelTrKC1)
;

denom2
=
fabs(cdKdA2) + fabs(cdKdB2) + fabs(cdKdC2) + fabs(codelTrKA2) + 
  fabs(codelTrKB2) + fabs(codelTrKC2)
;

denom3
=
fabs(cdKdA3) + fabs(cdKdB3) + fabs(cdKdC3) + fabs(codelTrKA3) + 
  fabs(codelTrKB3) + fabs(codelTrKC3)
;


} else { /* if (!TermByTerm) */

denom1
=
fabs(cdKudd111 + cdKudd212 + cdKudd313) + 
  fabs(codelK111*ginv11 + codelK122*ginv22 + 
    2.*(codelK112*ginv12 + codelK113*ginv13 + codelK123*ginv23) + 
    codelK133*ginv33)
;

denom2
=
fabs(cdKudd112 + cdKudd222 + cdKudd323) + 
  fabs(codelK211*ginv11 + codelK222*ginv22 + 
    2.*(codelK212*ginv12 + codelK213*ginv13 + codelK223*ginv23) + 
    codelK233*ginv33)
;

denom3
=
fabs(cdKudd113 + cdKudd223 + cdKudd333) + 
  fabs(codelK311*ginv11 + codelK322*ginv22 + 
    2.*(codelK312*ginv12 + codelK313*ginv13 + codelK323*ginv23) + 
    codelK333*ginv33)
;

}
/* if (TermByTerm) */


normmom1[ijk]
=
Cal(denom1 <= 0.,0.,mom1[ijk]/denom1)
;

normmom2[ijk]
=
Cal(denom2 <= 0.,0.,mom2[ijk]/denom2)
;

normmom3[ijk]
=
Cal(denom3 <= 0.,0.,mom3[ijk]/denom3)
;

}
/* if (TermByTerm) */




} endfor_ijk_openmp; /* loop i, j, k */



bampi_openmp_stop


}  /* function */

/* adm_constraints_N.c */
/* nvars = 39, nauxs = 397, n* = 4380,  n/ = 102,  n+ = 7130, n = 11612, O = 0 */
