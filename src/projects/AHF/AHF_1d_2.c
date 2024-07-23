/* AHF_1d_2.c */
/* Copyright (C) 1998 Bernd Bruegmann, 27.5.2009 */
/* Produced with Mathematica */

#include "bam.h"
#include "AHF.h"

#define Power(x,y) (pow((double) (x), (double) (y)))
#define Log(x)     log((double) (x)))
#define pow2(x)    ((x)*(x))
#define pow2inv(x) (1.0/((x)*(x)))

#define Cal(x,y,z) ((x)?(y):(z))




void AHF_1d_2(tVarList *vars, tVarList *ahfvars)
{

tL *level = ahfvars->level;
double *g11 = vldataptr(vars, 0);
double *g12 = vldataptr(vars, 1);
double *g13 = vldataptr(vars, 2);
double *g22 = vldataptr(vars, 3);
double *g23 = vldataptr(vars, 4);
double *g33 = vldataptr(vars, 5);
double *A11 = vldataptr(vars, 6);
double *A12 = vldataptr(vars, 7);
double *A13 = vldataptr(vars, 8);
double *A22 = vldataptr(vars, 9);
double *A23 = vldataptr(vars, 10);
double *A33 = vldataptr(vars, 11);
double *K = vldataptr(vars, 12);
double *phi = vldataptr(vars, 13);
double *s1 = vldataptr(ahfvars, 0);
double *s2 = vldataptr(ahfvars, 1);
double *s3 = vldataptr(ahfvars, 2);
double *H = vldataptr(ahfvars, 3);


double dx = level->dx;
double dy = level->dy;
double dz = level->dz;
double oo2dx = 1/(2*dx), oo2dy = 1/(2*dy), oo2dz = 1/(2*dz);
double oodx2 = 1/(dx*dx), oody2 = 1/(dy*dy), oodz2 = 1/(dz*dz);
double oo4dxdy = 1/(4*dx*dy);
double oo4dxdz = 1/(4*dx*dz);
double oo4dydz = 1/(4*dy*dz);
long ccc;

double deldelg1111;
double deldelg1112;
double deldelg1113;
double deldelg1122;
double deldelg1123;
double deldelg1133;
double deldelg1211;
double deldelg1212;
double deldelg1213;
double deldelg1222;
double deldelg1223;
double deldelg1233;
double deldelg1311;
double deldelg1312;
double deldelg1313;
double deldelg1322;
double deldelg1323;
double deldelg1333;
double deldelg2211;
double deldelg2212;
double deldelg2213;
double deldelg2222;
double deldelg2223;
double deldelg2233;
double deldelg2311;
double deldelg2312;
double deldelg2313;
double deldelg2322;
double deldelg2323;
double deldelg2333;
double deldelg3311;
double deldelg3312;
double deldelg3313;
double deldelg3322;
double deldelg3323;
double deldelg3333;
double delg111;
double delg112;
double delg113;
double delg122;
double delg123;
double delg133;
double delg211;
double delg212;
double delg213;
double delg222;
double delg223;
double delg233;
double delg311;
double delg312;
double delg313;
double delg322;
double delg323;
double delg333;
double detginv;
double dphi1;
double dphi2;
double dphi3;
double ds11;
double Ds11;
double ds12;
double Ds12;
double ds13;
double Ds13;
double ds22;
double Ds22;
double ds23;
double Ds23;
double ds33;
double Ds33;
double gamma111;
double gamma112;
double gamma113;
double gamma122;
double gamma123;
double gamma133;
double gamma211;
double gamma212;
double gamma213;
double gamma222;
double gamma223;
double gamma233;
double gamma311;
double gamma312;
double gamma313;
double gamma322;
double gamma323;
double gamma333;
double gammado111;
double gammado112;
double gammado113;
double gammado122;
double gammado123;
double gammado133;
double gammado211;
double gammado212;
double gammado213;
double gammado222;
double gammado223;
double gammado233;
double gammado311;
double gammado312;
double gammado313;
double gammado322;
double gammado323;
double gammado333;
double ggamma111;
double ggamma112;
double ggamma113;
double ggamma121;
double ggamma122;
double ggamma123;
double ggamma131;
double ggamma132;
double ggamma133;
double ggamma211;
double ggamma212;
double ggamma213;
double ggamma221;
double ggamma222;
double ggamma223;
double ggamma231;
double ggamma232;
double ggamma233;
double ggamma311;
double ggamma312;
double ggamma313;
double ggamma321;
double ggamma322;
double ggamma323;
double ggamma331;
double ggamma332;
double ggamma333;
double ginv11;
double ginv12;
double ginv13;
double ginv22;
double ginv23;
double ginv33;
double KK11;
double KK12;
double KK13;
double KK22;
double KK23;
double KK33;



forinner25(level,oo2dx,oo2dy,oo2dz) {


delg111
=
oo2dx*(-g11[mcc] + g11[pcc])
;

delg112
=
oo2dx*(-g12[mcc] + g12[pcc])
;

delg113
=
oo2dx*(-g13[mcc] + g13[pcc])
;

delg122
=
oo2dx*(-g22[mcc] + g22[pcc])
;

delg123
=
oo2dx*(-g23[mcc] + g23[pcc])
;

delg133
=
oo2dx*(-g33[mcc] + g33[pcc])
;

delg211
=
oo2dy*(-g11[cmc] + g11[cpc])
;

delg212
=
oo2dy*(-g12[cmc] + g12[cpc])
;

delg213
=
oo2dy*(-g13[cmc] + g13[cpc])
;

delg222
=
oo2dy*(-g22[cmc] + g22[cpc])
;

delg223
=
oo2dy*(-g23[cmc] + g23[cpc])
;

delg233
=
oo2dy*(-g33[cmc] + g33[cpc])
;

delg311
=
oo2dz*(-g11[ccm] + g11[ccp])
;

delg312
=
oo2dz*(-g12[ccm] + g12[ccp])
;

delg313
=
oo2dz*(-g13[ccm] + g13[ccp])
;

delg322
=
oo2dz*(-g22[ccm] + g22[ccp])
;

delg323
=
oo2dz*(-g23[ccm] + g23[ccp])
;

delg333
=
oo2dz*(-g33[ccm] + g33[ccp])
;

deldelg1111
=
oodx2*(-2.*g11[ccc] + g11[mcc] + g11[pcc])
;

deldelg1112
=
oodx2*(-2.*g12[ccc] + g12[mcc] + g12[pcc])
;

deldelg1113
=
oodx2*(-2.*g13[ccc] + g13[mcc] + g13[pcc])
;

deldelg1122
=
oodx2*(-2.*g22[ccc] + g22[mcc] + g22[pcc])
;

deldelg1123
=
oodx2*(-2.*g23[ccc] + g23[mcc] + g23[pcc])
;

deldelg1133
=
oodx2*(-2.*g33[ccc] + g33[mcc] + g33[pcc])
;

deldelg1211
=
oo4dxdy*(g11[mmc] - g11[mpc] - g11[pmc] + g11[ppc])
;

deldelg1212
=
oo4dxdy*(g12[mmc] - g12[mpc] - g12[pmc] + g12[ppc])
;

deldelg1213
=
oo4dxdy*(g13[mmc] - g13[mpc] - g13[pmc] + g13[ppc])
;

deldelg1222
=
oo4dxdy*(g22[mmc] - g22[mpc] - g22[pmc] + g22[ppc])
;

deldelg1223
=
oo4dxdy*(g23[mmc] - g23[mpc] - g23[pmc] + g23[ppc])
;

deldelg1233
=
oo4dxdy*(g33[mmc] - g33[mpc] - g33[pmc] + g33[ppc])
;

deldelg1311
=
oo4dxdz*(g11[mcm] - g11[mcp] - g11[pcm] + g11[pcp])
;

deldelg1312
=
oo4dxdz*(g12[mcm] - g12[mcp] - g12[pcm] + g12[pcp])
;

deldelg1313
=
oo4dxdz*(g13[mcm] - g13[mcp] - g13[pcm] + g13[pcp])
;

deldelg1322
=
oo4dxdz*(g22[mcm] - g22[mcp] - g22[pcm] + g22[pcp])
;

deldelg1323
=
oo4dxdz*(g23[mcm] - g23[mcp] - g23[pcm] + g23[pcp])
;

deldelg1333
=
oo4dxdz*(g33[mcm] - g33[mcp] - g33[pcm] + g33[pcp])
;

deldelg2211
=
oody2*(-2.*g11[ccc] + g11[cmc] + g11[cpc])
;

deldelg2212
=
oody2*(-2.*g12[ccc] + g12[cmc] + g12[cpc])
;

deldelg2213
=
oody2*(-2.*g13[ccc] + g13[cmc] + g13[cpc])
;

deldelg2222
=
oody2*(-2.*g22[ccc] + g22[cmc] + g22[cpc])
;

deldelg2223
=
oody2*(-2.*g23[ccc] + g23[cmc] + g23[cpc])
;

deldelg2233
=
oody2*(-2.*g33[ccc] + g33[cmc] + g33[cpc])
;

deldelg2311
=
oo4dydz*(g11[cmm] - g11[cmp] - g11[cpm] + g11[cpp])
;

deldelg2312
=
oo4dydz*(g12[cmm] - g12[cmp] - g12[cpm] + g12[cpp])
;

deldelg2313
=
oo4dydz*(g13[cmm] - g13[cmp] - g13[cpm] + g13[cpp])
;

deldelg2322
=
oo4dydz*(g22[cmm] - g22[cmp] - g22[cpm] + g22[cpp])
;

deldelg2323
=
oo4dydz*(g23[cmm] - g23[cmp] - g23[cpm] + g23[cpp])
;

deldelg2333
=
oo4dydz*(g33[cmm] - g33[cmp] - g33[cpm] + g33[cpp])
;

deldelg3311
=
oodz2*(-2.*g11[ccc] + g11[ccm] + g11[ccp])
;

deldelg3312
=
oodz2*(-2.*g12[ccc] + g12[ccm] + g12[ccp])
;

deldelg3313
=
oodz2*(-2.*g13[ccc] + g13[ccm] + g13[ccp])
;

deldelg3322
=
oodz2*(-2.*g22[ccc] + g22[ccm] + g22[ccp])
;

deldelg3323
=
oodz2*(-2.*g23[ccc] + g23[ccm] + g23[ccp])
;

deldelg3333
=
oodz2*(-2.*g33[ccc] + g33[ccm] + g33[ccp])
;

detginv
=
1/(2.*g12[ccc]*g13[ccc]*g23[ccc] + g11[ccc]*g22[ccc]*g33[ccc] - 
    g33[ccc]*pow2(g12[ccc]) - g22[ccc]*pow2(g13[ccc]) - 
    g11[ccc]*pow2(g23[ccc]))
;

ginv11
=
detginv*(g22[ccc]*g33[ccc] - pow2(g23[ccc]))
;

ginv12
=
detginv*(g13[ccc]*g23[ccc] - g12[ccc]*g33[ccc])
;

ginv13
=
detginv*(-(g13[ccc]*g22[ccc]) + g12[ccc]*g23[ccc])
;

ginv22
=
detginv*(g11[ccc]*g33[ccc] - pow2(g13[ccc]))
;

ginv23
=
detginv*(g12[ccc]*g13[ccc] - g11[ccc]*g23[ccc])
;

ginv33
=
detginv*(g11[ccc]*g22[ccc] - pow2(g12[ccc]))
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

ds11
=
oo2dx*(-s1[mcc] + s1[pcc])
;

ds12
=
oo2dy*(-s1[cmc] + s1[cpc])
;

ds13
=
oo2dz*(-s1[ccm] + s1[ccp])
;

ds22
=
oo2dy*(-s2[cmc] + s2[cpc])
;

ds23
=
oo2dz*(-s2[ccm] + s2[ccp])
;

ds33
=
oo2dz*(-s3[ccm] + s3[ccp])
;

dphi1
=
oo2dx*(-phi[mcc] + phi[pcc])
;

dphi2
=
oo2dy*(-phi[cmc] + phi[cpc])
;

dphi3
=
oo2dz*(-phi[ccm] + phi[ccp])
;

ggamma111
=
gamma111 - 2.*(dphi2*ginv12 + dphi3*ginv13)*g11[ccc] + 
  dphi1*(4. - 2.*ginv11*g11[ccc])
;

ggamma112
=
gamma112 - 2.*(dphi1*ginv11 + dphi3*ginv13)*g12[ccc] + 
  dphi2*(2. - 2.*ginv12*g12[ccc])
;

ggamma113
=
gamma113 - 2.*(dphi1*ginv11 + dphi2*ginv12)*g13[ccc] + 
  dphi3*(2. - 2.*ginv13*g13[ccc])
;

ggamma121
=
gamma112 - 2.*(dphi1*ginv11 + dphi3*ginv13)*g12[ccc] + 
  dphi2*(2. - 2.*ginv12*g12[ccc])
;

ggamma122
=
gamma122 - 2.*(dphi1*ginv11 + dphi2*ginv12 + dphi3*ginv13)*g22[ccc]
;

ggamma123
=
gamma123 - 2.*(dphi1*ginv11 + dphi2*ginv12 + dphi3*ginv13)*g23[ccc]
;

ggamma131
=
gamma113 - 2.*(dphi1*ginv11 + dphi2*ginv12)*g13[ccc] + 
  dphi3*(2. - 2.*ginv13*g13[ccc])
;

ggamma132
=
gamma123 - 2.*(dphi1*ginv11 + dphi2*ginv12 + dphi3*ginv13)*g23[ccc]
;

ggamma133
=
gamma133 - 2.*(dphi1*ginv11 + dphi2*ginv12 + dphi3*ginv13)*g33[ccc]
;

ggamma211
=
gamma211 - 2.*(dphi1*ginv12 + dphi2*ginv22 + dphi3*ginv23)*g11[ccc]
;

ggamma212
=
gamma212 - 2.*(dphi2*ginv22 + dphi3*ginv23)*g12[ccc] + 
  dphi1*(2. - 2.*ginv12*g12[ccc])
;

ggamma213
=
gamma213 - 2.*(dphi1*ginv12 + dphi2*ginv22 + dphi3*ginv23)*g13[ccc]
;

ggamma221
=
gamma212 - 2.*(dphi2*ginv22 + dphi3*ginv23)*g12[ccc] + 
  dphi1*(2. - 2.*ginv12*g12[ccc])
;

ggamma222
=
gamma222 - 2.*(dphi1*ginv12 + dphi3*ginv23)*g22[ccc] + 
  dphi2*(4. - 2.*ginv22*g22[ccc])
;

ggamma223
=
gamma223 - 2.*(dphi1*ginv12 + dphi2*ginv22)*g23[ccc] + 
  dphi3*(2. - 2.*ginv23*g23[ccc])
;

ggamma231
=
gamma213 - 2.*(dphi1*ginv12 + dphi2*ginv22 + dphi3*ginv23)*g13[ccc]
;

ggamma232
=
gamma223 - 2.*(dphi1*ginv12 + dphi2*ginv22)*g23[ccc] + 
  dphi3*(2. - 2.*ginv23*g23[ccc])
;

ggamma233
=
gamma233 - 2.*(dphi1*ginv12 + dphi2*ginv22 + dphi3*ginv23)*g33[ccc]
;

ggamma311
=
gamma311 - 2.*(dphi1*ginv13 + dphi2*ginv23 + dphi3*ginv33)*g11[ccc]
;

ggamma312
=
gamma312 - 2.*(dphi1*ginv13 + dphi2*ginv23 + dphi3*ginv33)*g12[ccc]
;

ggamma313
=
gamma313 - 2.*(dphi2*ginv23 + dphi3*ginv33)*g13[ccc] + 
  dphi1*(2. - 2.*ginv13*g13[ccc])
;

ggamma321
=
gamma312 - 2.*(dphi1*ginv13 + dphi2*ginv23 + dphi3*ginv33)*g12[ccc]
;

ggamma322
=
gamma322 - 2.*(dphi1*ginv13 + dphi2*ginv23 + dphi3*ginv33)*g22[ccc]
;

ggamma323
=
gamma323 - 2.*(dphi1*ginv13 + dphi3*ginv33)*g23[ccc] + 
  dphi2*(2. - 2.*ginv23*g23[ccc])
;

ggamma331
=
gamma313 - 2.*(dphi2*ginv23 + dphi3*ginv33)*g13[ccc] + 
  dphi1*(2. - 2.*ginv13*g13[ccc])
;

ggamma332
=
gamma323 - 2.*(dphi1*ginv13 + dphi3*ginv33)*g23[ccc] + 
  dphi2*(2. - 2.*ginv23*g23[ccc])
;

ggamma333
=
gamma333 - 2.*(dphi1*ginv13 + dphi2*ginv23)*g33[ccc] + 
  dphi3*(4. - 2.*ginv33*g33[ccc])
;

Ds11
=
ds11 + ggamma111*s1[ccc] + ggamma112*s2[ccc] + ggamma113*s3[ccc]
;

Ds12
=
ds12 + ggamma211*s1[ccc] + ggamma212*s2[ccc] + ggamma213*s3[ccc]
;

Ds13
=
ds13 + ggamma311*s1[ccc] + ggamma312*s2[ccc] + ggamma313*s3[ccc]
;

Ds22
=
ds22 + ggamma221*s1[ccc] + ggamma222*s2[ccc] + ggamma223*s3[ccc]
;

Ds23
=
ds23 + ggamma321*s1[ccc] + ggamma322*s2[ccc] + ggamma323*s3[ccc]
;

Ds33
=
ds33 + ggamma331*s1[ccc] + ggamma332*s2[ccc] + ggamma333*s3[ccc]
;

KK11
=
exp(4.*phi[ccc])*(A11[ccc] + 0.33333333333333333333*g11[ccc]*K[ccc])
;

KK12
=
exp(4.*phi[ccc])*(A12[ccc] + 0.33333333333333333333*g12[ccc]*K[ccc])
;

KK13
=
exp(4.*phi[ccc])*(A13[ccc] + 0.33333333333333333333*g13[ccc]*K[ccc])
;

KK22
=
exp(4.*phi[ccc])*(A22[ccc] + 0.33333333333333333333*g22[ccc]*K[ccc])
;

KK23
=
exp(4.*phi[ccc])*(A23[ccc] + 0.33333333333333333333*g23[ccc]*K[ccc])
;

KK33
=
exp(4.*phi[ccc])*(A33[ccc] + 0.33333333333333333333*g33[ccc]*K[ccc])
;

H[ccc]
=
Ds11 + Ds22 + Ds33 - K[ccc] + 2.*
   (KK23*s2[ccc]*s3[ccc] + s1[ccc]*(KK12*s2[ccc] + KK13*s3[ccc])) + 
  KK11*pow2(s1[ccc]) + KK22*pow2(s2[ccc]) + KK33*pow2(s3[ccc])
;



} endforinner; /* loop i, j, k */

}  /* function */

/* AHF_1d_2.c */
/* nvars = 18, nauxs = 145, n* = 437,  n/ = 21,  n+ = 358, n = 816, O = 0 */
