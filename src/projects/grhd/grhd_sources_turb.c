/* grhd_sources_turb.c */
/* Copyright (C) 1998 Bernd Bruegmann, 20.6.2023 */
/* Produced with Mathematica */

#include "bam.h"
#include "grhd.h"

#define Power(x,y) (pow((double) (x), (double) (y)))
#define Log(x)     log((double) (x))
#define pow2(x)    ((x)*(x))
#define pow2inv(x) (1.0/((x)*(x)))

#define Cal(x,y,z) ((x)?(y):(z))

#define SetSourceToZero 0




void grhd_sources_turb(tVarList *unew, tVarList *upre, double c, tVarList *ucur, tVarList *primVars, tVarList *otherVars, tVarList *matterADMVars)
{

tL *level = ucur->level;
int addlinear = (c != 0.0l);

const int order_centered = Geti("order_centered");
const int useprim = Getv("grhd_source_computation","prim");
double cooling = Getd("grhd_cooling");
double tcool = Getd("grhd_cooling_time");
double EOSK = EOS.K;
double EOSGmo = EOS.GAMMAMO;
double lmix = GRHD.turb_lmix;
double pres_eos;
double cs2;

double *MD = vldataptr(ucur, MATTER.INDX_VAR_q);
double *Mtau = vldataptr(ucur, 1 + MATTER.INDX_VAR_q);
double *MS1 = vldataptr(ucur, 2 + MATTER.INDX_VAR_q);
double *MS2 = vldataptr(ucur, 3 + MATTER.INDX_VAR_q);
double *MS3 = vldataptr(ucur, 4 + MATTER.INDX_VAR_q);
double *nMD = vldataptr(unew, MATTER.INDX_VAR_q);
double *nMtau = vldataptr(unew, 1 + MATTER.INDX_VAR_q);
double *nMS1 = vldataptr(unew, 2 + MATTER.INDX_VAR_q);
double *nMS2 = vldataptr(unew, 3 + MATTER.INDX_VAR_q);
double *nMS3 = vldataptr(unew, 4 + MATTER.INDX_VAR_q);
double *pMD = vldataptr(upre, MATTER.INDX_VAR_q);
double *pMtau = vldataptr(upre, 1 + MATTER.INDX_VAR_q);
double *pMS1 = vldataptr(upre, 2 + MATTER.INDX_VAR_q);
double *pMS2 = vldataptr(upre, 3 + MATTER.INDX_VAR_q);
double *pMS3 = vldataptr(upre, 4 + MATTER.INDX_VAR_q);
double *alpha = vldataptr(ucur, MATTER.INDX_VAR_alp);
double *beta1 = vldataptr(ucur, MATTER.INDX_VAR_bet);
double *beta2 = vldataptr(ucur, 1 + MATTER.INDX_VAR_bet);
double *beta3 = vldataptr(ucur, 2 + MATTER.INDX_VAR_bet);
int index_adm_gxx = Ind("adm_gxx");
double *g11 = level->v[index_adm_gxx + 0];
double *g12 = level->v[index_adm_gxx + 1];
double *g13 = level->v[index_adm_gxx + 2];
double *g22 = level->v[index_adm_gxx + 3];
double *g23 = level->v[index_adm_gxx + 4];
double *g33 = level->v[index_adm_gxx + 5];
int index_adm_Kxx = Ind("adm_Kxx");
double *K11 = level->v[index_adm_Kxx + 0];
double *K12 = level->v[index_adm_Kxx + 1];
double *K13 = level->v[index_adm_Kxx + 2];
double *K22 = level->v[index_adm_Kxx + 3];
double *K23 = level->v[index_adm_Kxx + 4];
double *K33 = level->v[index_adm_Kxx + 5];
int index_matter_mask = Ind("matter_mask");
double *mask = level->v[index_matter_mask + 0];
const double *xp = level->v[Ind("x")];
const double *yp = level->v[Ind("y")];
const double *zp = level->v[Ind("z")];
double *Mrho = vldataptr(primVars, 0);
double *Mepsl = vldataptr(primVars, 1);
double *Mv1 = vldataptr(primVars, 2);
double *Mv2 = vldataptr(primVars, 3);
double *Mv3 = vldataptr(primVars, 4);
double *Mp = vldataptr(primVars, 5);
double *Mvsqr = vldataptr(primVars, 6);
double *detmetric = vldataptr(primVars, 7);
double *MturbTau11 = vldataptr(otherVars, GRHD.IND_turbTau);
double *MturbTau12 = vldataptr(otherVars, 1 + GRHD.IND_turbTau);
double *MturbTau13 = vldataptr(otherVars, 2 + GRHD.IND_turbTau);
double *MturbTau22 = vldataptr(otherVars, 3 + GRHD.IND_turbTau);
double *MturbTau23 = vldataptr(otherVars, 4 + GRHD.IND_turbTau);
double *MturbTau33 = vldataptr(otherVars, 5 + GRHD.IND_turbTau);
double *Madmrho = vldataptr(matterADMVars, 0);
double *MadmS1 = vldataptr(matterADMVars, 1);
double *MadmS2 = vldataptr(matterADMVars, 2);
double *MadmS3 = vldataptr(matterADMVars, 3);
double *MadmSS11 = vldataptr(matterADMVars, 4);
double *MadmSS12 = vldataptr(matterADMVars, 5);
double *MadmSS13 = vldataptr(matterADMVars, 6);
double *MadmSS22 = vldataptr(matterADMVars, 7);
double *MadmSS23 = vldataptr(matterADMVars, 8);
double *MadmSS33 = vldataptr(matterADMVars, 9);
double *MadmST = vldataptr(matterADMVars, 10);




bampi_openmp_start


double dx = level->dx;
double dy = level->dy;
double dz = level->dz;
double oo2dx = 1/(2*dx), oo2dy = 1/(2*dy), oo2dz = 1/(2*dz);
double oodx2 = 1/(dx*dx), oody2 = 1/(dy*dy), oodz2 = 1/(dz*dz);
double oo4dxdy = 1/(4*dx*dy);
double oo4dxdz = 1/(4*dx*dz);
double oo4dydz = 1/(4*dy*dz);

double betadown1 = 0.;
double betadown2 = 0.;
double betadown3 = 0.;
double conD = 0.;
double conS1 = 0.;
double conS2 = 0.;
double conS3 = 0.;
double conSi1 = 0.;
double conSi2 = 0.;
double conSi3 = 0.;
double conT = 0.;
double dalpha1 = 0.;
double dalpha2 = 0.;
double dalpha3 = 0.;
double dbeta11 = 0.;
double dbeta12 = 0.;
double dbeta13 = 0.;
double dbeta21 = 0.;
double dbeta22 = 0.;
double dbeta23 = 0.;
double dbeta31 = 0.;
double dbeta32 = 0.;
double dbeta33 = 0.;
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
double delv11 = 0.;
double delv12 = 0.;
double delv13 = 0.;
double delv21 = 0.;
double delv22 = 0.;
double delv23 = 0.;
double delv31 = 0.;
double delv32 = 0.;
double delv33 = 0.;
double detginv = 0.;
double EOSKoGmo = 0.;
double epsth = 0.;
double ginv11 = 0.;
double ginv12 = 0.;
double ginv13 = 0.;
double ginv22 = 0.;
double ginv23 = 0.;
double ginv33 = 0.;
double Lambda = 0.;
double matrixdetg = 0.;
double sqD = 0.;
double sqrtdetgamma = 0.;
double sqrtmdetg = 0.;
double sqS1 = 0.;
double sqS2 = 0.;
double sqS3 = 0.;
double sqT = 0.;
double T00 = 0.;
double T0i1 = 0.;
double T0i2 = 0.;
double T0i3 = 0.;
double T0j1 = 0.;
double T0j2 = 0.;
double T0j3 = 0.;
double Tij11 = 0.;
double Tij12 = 0.;
double Tij13 = 0.;
double Tij22 = 0.;
double Tij23 = 0.;
double Tij33 = 0.;
double udown1 = 0.;
double udown2 = 0.;
double udown3 = 0.;
double uup1 = 0.;
double uup2 = 0.;
double uup3 = 0.;
double uupt = 0.;
double v1 = 0.;
double v2 = 0.;
double v3 = 0.;
double vshift1 = 0.;
double vshift2 = 0.;
double vshift3 = 0.;
double vsq = 0.;
double vup1 = 0.;
double vup2 = 0.;
double vup3 = 0.;
double W = 0.;
double W2hrho = 0.;



forinnerpoints_ijk_openmp(level) {


matrixdetg
=
2.*g12[ijk]*g13[ijk]*g23[ijk] + g11[ijk]*g22[ijk]*g33[ijk] - 
  g33[ijk]*pow2(g12[ijk]) - g22[ijk]*pow2(g13[ijk]) - g11[ijk]*pow2(g23[ijk])
;


 if (((MATTER.USEMASK) && (mask[ijk]>.99)) || matrixdetg<=0.) { 

sqD
=
0
;

sqS1
=
0
;

sqS2
=
0
;

sqS3
=
0
;

sqT
=
0
;

Madmrho[ijk]
=
0
;

MadmS1[ijk]
=
0
;

MadmS2[ijk]
=
0
;

MadmS3[ijk]
=
0
;

MadmSS11[ijk]
=
0
;

MadmSS12[ijk]
=
0
;

MadmSS13[ijk]
=
0
;

MadmSS22[ijk]
=
0
;

MadmSS23[ijk]
=
0
;

MadmSS33[ijk]
=
0
;

MadmST[ijk]
=
0
;


 } else { 


 if (alpha[ijk]<=ALPHAMIN){  


 printf("alpha = %e < alpha_min = %e\n",alpha[ijk],ALPHAMIN); 


 printf("  x=%e y=%e z=%e\n",Ptr(level,"x")[ijk],Ptr(level,"y")[ijk],Ptr(level,"z")[ijk]);  


 if (ijkinsidefinerlevel(box,ijk)==0) printf("  point is NOT inside finer box NOR in some symmetry area\n");  


 else printf("  point is inside finer box/ in symmetry \n");}  


                                                                                                                                      
        CheckForNANandINF(11,MD[ijk],MS1[ijk],MS2[ijk],MS3[ijk],Mtau[ijk], Mrho[ijk],Mepsl[ijk],Mp[ijk],Mv1[ijk],Mv2[ijk],Mv3[ijk]); 

if (order_centered == 2 || boundaryNaway(1)) { 

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

dalpha1
=
oo2dx*(-alpha[-di + ijk] + alpha[di + ijk])
;

dalpha2
=
oo2dy*(-alpha[-dj + ijk] + alpha[dj + ijk])
;

dalpha3
=
oo2dz*(-alpha[-dk + ijk] + alpha[dk + ijk])
;

dbeta11
=
oo2dx*(-beta1[-di + ijk] + beta1[di + ijk])
;

dbeta12
=
oo2dx*(-beta2[-di + ijk] + beta2[di + ijk])
;

dbeta13
=
oo2dx*(-beta3[-di + ijk] + beta3[di + ijk])
;

dbeta21
=
oo2dy*(-beta1[-dj + ijk] + beta1[dj + ijk])
;

dbeta22
=
oo2dy*(-beta2[-dj + ijk] + beta2[dj + ijk])
;

dbeta23
=
oo2dy*(-beta3[-dj + ijk] + beta3[dj + ijk])
;

dbeta31
=
oo2dz*(-beta1[-dk + ijk] + beta1[dk + ijk])
;

dbeta32
=
oo2dz*(-beta2[-dk + ijk] + beta2[dk + ijk])
;

dbeta33
=
oo2dz*(-beta3[-dk + ijk] + beta3[dk + ijk])
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

dalpha1
=
0.16666666666666666667*oo2dx*(alpha[-2*di + ijk] + 
    8.*(-alpha[-di + ijk] + alpha[di + ijk]) - alpha[2*di + ijk])
;

dalpha2
=
0.16666666666666666667*oo2dy*(alpha[-2*dj + ijk] + 
    8.*(-alpha[-dj + ijk] + alpha[dj + ijk]) - alpha[2*dj + ijk])
;

dalpha3
=
0.16666666666666666667*oo2dz*(alpha[-2*dk + ijk] + 
    8.*(-alpha[-dk + ijk] + alpha[dk + ijk]) - alpha[2*dk + ijk])
;

dbeta11
=
0.16666666666666666667*oo2dx*(beta1[-2*di + ijk] + 
    8.*(-beta1[-di + ijk] + beta1[di + ijk]) - beta1[2*di + ijk])
;

dbeta12
=
0.16666666666666666667*oo2dx*(beta2[-2*di + ijk] + 
    8.*(-beta2[-di + ijk] + beta2[di + ijk]) - beta2[2*di + ijk])
;

dbeta13
=
0.16666666666666666667*oo2dx*(beta3[-2*di + ijk] + 
    8.*(-beta3[-di + ijk] + beta3[di + ijk]) - beta3[2*di + ijk])
;

dbeta21
=
0.16666666666666666667*oo2dy*(beta1[-2*dj + ijk] + 
    8.*(-beta1[-dj + ijk] + beta1[dj + ijk]) - beta1[2*dj + ijk])
;

dbeta22
=
0.16666666666666666667*oo2dy*(beta2[-2*dj + ijk] + 
    8.*(-beta2[-dj + ijk] + beta2[dj + ijk]) - beta2[2*dj + ijk])
;

dbeta23
=
0.16666666666666666667*oo2dy*(beta3[-2*dj + ijk] + 
    8.*(-beta3[-dj + ijk] + beta3[dj + ijk]) - beta3[2*dj + ijk])
;

dbeta31
=
0.16666666666666666667*oo2dz*(beta1[-2*dk + ijk] + 
    8.*(-beta1[-dk + ijk] + beta1[dk + ijk]) - beta1[2*dk + ijk])
;

dbeta32
=
0.16666666666666666667*oo2dz*(beta2[-2*dk + ijk] + 
    8.*(-beta2[-dk + ijk] + beta2[dk + ijk]) - beta2[2*dk + ijk])
;

dbeta33
=
0.16666666666666666667*oo2dz*(beta3[-2*dk + ijk] + 
    8.*(-beta3[-dk + ijk] + beta3[dk + ijk]) - beta3[2*dk + ijk])
;


} else if (order_centered == 6 || boundaryNaway(3)) { 

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

dalpha1
=
0.033333333333333333333*oo2dx*(-alpha[-3*di + ijk] + 
    45.*(-alpha[-di + ijk] + alpha[di + ijk]) + 
    9.*(alpha[-2*di + ijk] - alpha[2*di + ijk]) + alpha[3*di + ijk])
;

dalpha2
=
0.033333333333333333333*oo2dy*(-alpha[-3*dj + ijk] + 
    45.*(-alpha[-dj + ijk] + alpha[dj + ijk]) + 
    9.*(alpha[-2*dj + ijk] - alpha[2*dj + ijk]) + alpha[3*dj + ijk])
;

dalpha3
=
0.033333333333333333333*oo2dz*(-alpha[-3*dk + ijk] + 
    45.*(-alpha[-dk + ijk] + alpha[dk + ijk]) + 
    9.*(alpha[-2*dk + ijk] - alpha[2*dk + ijk]) + alpha[3*dk + ijk])
;

dbeta11
=
0.033333333333333333333*oo2dx*(-beta1[-3*di + ijk] + 
    45.*(-beta1[-di + ijk] + beta1[di + ijk]) + 
    9.*(beta1[-2*di + ijk] - beta1[2*di + ijk]) + beta1[3*di + ijk])
;

dbeta12
=
0.033333333333333333333*oo2dx*(-beta2[-3*di + ijk] + 
    45.*(-beta2[-di + ijk] + beta2[di + ijk]) + 
    9.*(beta2[-2*di + ijk] - beta2[2*di + ijk]) + beta2[3*di + ijk])
;

dbeta13
=
0.033333333333333333333*oo2dx*(-beta3[-3*di + ijk] + 
    45.*(-beta3[-di + ijk] + beta3[di + ijk]) + 
    9.*(beta3[-2*di + ijk] - beta3[2*di + ijk]) + beta3[3*di + ijk])
;

dbeta21
=
0.033333333333333333333*oo2dy*(-beta1[-3*dj + ijk] + 
    45.*(-beta1[-dj + ijk] + beta1[dj + ijk]) + 
    9.*(beta1[-2*dj + ijk] - beta1[2*dj + ijk]) + beta1[3*dj + ijk])
;

dbeta22
=
0.033333333333333333333*oo2dy*(-beta2[-3*dj + ijk] + 
    45.*(-beta2[-dj + ijk] + beta2[dj + ijk]) + 
    9.*(beta2[-2*dj + ijk] - beta2[2*dj + ijk]) + beta2[3*dj + ijk])
;

dbeta23
=
0.033333333333333333333*oo2dy*(-beta3[-3*dj + ijk] + 
    45.*(-beta3[-dj + ijk] + beta3[dj + ijk]) + 
    9.*(beta3[-2*dj + ijk] - beta3[2*dj + ijk]) + beta3[3*dj + ijk])
;

dbeta31
=
0.033333333333333333333*oo2dz*(-beta1[-3*dk + ijk] + 
    45.*(-beta1[-dk + ijk] + beta1[dk + ijk]) + 
    9.*(beta1[-2*dk + ijk] - beta1[2*dk + ijk]) + beta1[3*dk + ijk])
;

dbeta32
=
0.033333333333333333333*oo2dz*(-beta2[-3*dk + ijk] + 
    45.*(-beta2[-dk + ijk] + beta2[dk + ijk]) + 
    9.*(beta2[-2*dk + ijk] - beta2[2*dk + ijk]) + beta2[3*dk + ijk])
;

dbeta33
=
0.033333333333333333333*oo2dz*(-beta3[-3*dk + ijk] + 
    45.*(-beta3[-dk + ijk] + beta3[dk + ijk]) + 
    9.*(beta3[-2*dk + ijk] - beta3[2*dk + ijk]) + beta3[3*dk + ijk])
;


} else errorexit("order is not implemented yet"); 

matrixdetg
=
2.*g12[ijk]*g13[ijk]*g23[ijk] + g11[ijk]*g22[ijk]*g33[ijk] - 
  g33[ijk]*pow2(g12[ijk]) - g22[ijk]*pow2(g13[ijk]) - g11[ijk]*pow2(g23[ijk])
;

detginv
=
1/matrixdetg
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

sqrtdetgamma
=
sqrt(matrixdetg)
;

sqrtmdetg
=
sqrtdetgamma*alpha[ijk]
;

vup1
=
Mv1[ijk]
;

vup2
=
Mv2[ijk]
;

vup3
=
Mv3[ijk]
;

v1
=
g11[ijk]*Mv1[ijk] + g12[ijk]*Mv2[ijk] + g13[ijk]*Mv3[ijk]
;

v2
=
g12[ijk]*Mv1[ijk] + g22[ijk]*Mv2[ijk] + g23[ijk]*Mv3[ijk]
;

v3
=
g13[ijk]*Mv1[ijk] + g23[ijk]*Mv2[ijk] + g33[ijk]*Mv3[ijk]
;

vsq
=
v1*vup1 + v2*vup2 + v3*vup3
;

W
=
1/sqrt(1. - vsq)
;


if (vsq>=1.) W = GRHD.HRSC_WlorMAX; 

W2hrho
=
(Mp[ijk] + Mrho[ijk] + Mepsl[ijk]*Mrho[ijk])*pow2(W)
;



/* conditional */
if (cooling) {

uupt
=
W/alpha[ijk]
;

uup1
=
W*(vup1 - beta1[ijk]/alpha[ijk])
;

uup2
=
W*(vup2 - beta2[ijk]/alpha[ijk])
;

uup3
=
W*(vup3 - beta3[ijk]/alpha[ijk])
;

betadown1
=
beta1[ijk]*g11[ijk] + beta2[ijk]*g12[ijk] + beta3[ijk]*g13[ijk]
;

betadown2
=
beta1[ijk]*g12[ijk] + beta2[ijk]*g22[ijk] + beta3[ijk]*g23[ijk]
;

betadown3
=
beta1[ijk]*g13[ijk] + beta2[ijk]*g23[ijk] + beta3[ijk]*g33[ijk]
;

udown1
=
betadown1*uupt + uup1*g11[ijk] + uup2*g12[ijk] + uup3*g13[ijk]
;

udown2
=
betadown2*uupt + uup1*g12[ijk] + uup2*g22[ijk] + uup3*g23[ijk]
;

udown3
=
betadown3*uupt + uup1*g13[ijk] + uup2*g23[ijk] + uup3*g33[ijk]
;

EOSKoGmo
=
EOSK/EOSGmo
;

epsth
=
Mepsl[ijk] - EOSKoGmo*Power(Mrho[ijk],EOSGmo)
;

Lambda
=
(epsth*Mrho[ijk])/tcool
;


} else { /* if (!cooling) */

Lambda
=
cooling
;

}
/* if (cooling) */



 if (alpha[ijk]>=0.2) { 


if (order_centered == 2 || boundaryNaway(1)) { 

delv11
=
oo2dx*(-Mv1[-di + ijk] + Mv1[di + ijk])
;

delv12
=
oo2dx*(-Mv2[-di + ijk] + Mv2[di + ijk])
;

delv13
=
oo2dx*(-Mv3[-di + ijk] + Mv3[di + ijk])
;

delv21
=
oo2dy*(-Mv1[-dj + ijk] + Mv1[dj + ijk])
;

delv22
=
oo2dy*(-Mv2[-dj + ijk] + Mv2[dj + ijk])
;

delv23
=
oo2dy*(-Mv3[-dj + ijk] + Mv3[dj + ijk])
;

delv31
=
oo2dz*(-Mv1[-dk + ijk] + Mv1[dk + ijk])
;

delv32
=
oo2dz*(-Mv2[-dk + ijk] + Mv2[dk + ijk])
;

delv33
=
oo2dz*(-Mv3[-dk + ijk] + Mv3[dk + ijk])
;


} else if (order_centered == 4 || boundaryNaway(2)) { 

delv11
=
0.16666666666666666667*oo2dx*(Mv1[-2*di + ijk] + 
    8.*(-Mv1[-di + ijk] + Mv1[di + ijk]) - Mv1[2*di + ijk])
;

delv12
=
0.16666666666666666667*oo2dx*(Mv2[-2*di + ijk] + 
    8.*(-Mv2[-di + ijk] + Mv2[di + ijk]) - Mv2[2*di + ijk])
;

delv13
=
0.16666666666666666667*oo2dx*(Mv3[-2*di + ijk] + 
    8.*(-Mv3[-di + ijk] + Mv3[di + ijk]) - Mv3[2*di + ijk])
;

delv21
=
0.16666666666666666667*oo2dy*(Mv1[-2*dj + ijk] + 
    8.*(-Mv1[-dj + ijk] + Mv1[dj + ijk]) - Mv1[2*dj + ijk])
;

delv22
=
0.16666666666666666667*oo2dy*(Mv2[-2*dj + ijk] + 
    8.*(-Mv2[-dj + ijk] + Mv2[dj + ijk]) - Mv2[2*dj + ijk])
;

delv23
=
0.16666666666666666667*oo2dy*(Mv3[-2*dj + ijk] + 
    8.*(-Mv3[-dj + ijk] + Mv3[dj + ijk]) - Mv3[2*dj + ijk])
;

delv31
=
0.16666666666666666667*oo2dz*(Mv1[-2*dk + ijk] + 
    8.*(-Mv1[-dk + ijk] + Mv1[dk + ijk]) - Mv1[2*dk + ijk])
;

delv32
=
0.16666666666666666667*oo2dz*(Mv2[-2*dk + ijk] + 
    8.*(-Mv2[-dk + ijk] + Mv2[dk + ijk]) - Mv2[2*dk + ijk])
;

delv33
=
0.16666666666666666667*oo2dz*(Mv3[-2*dk + ijk] + 
    8.*(-Mv3[-dk + ijk] + Mv3[dk + ijk]) - Mv3[2*dk + ijk])
;


} else if (order_centered == 6 || boundaryNaway(3)) { 

delv11
=
0.033333333333333333333*oo2dx*(-Mv1[-3*di + ijk] + 
    45.*(-Mv1[-di + ijk] + Mv1[di + ijk]) + 
    9.*(Mv1[-2*di + ijk] - Mv1[2*di + ijk]) + Mv1[3*di + ijk])
;

delv12
=
0.033333333333333333333*oo2dx*(-Mv2[-3*di + ijk] + 
    45.*(-Mv2[-di + ijk] + Mv2[di + ijk]) + 
    9.*(Mv2[-2*di + ijk] - Mv2[2*di + ijk]) + Mv2[3*di + ijk])
;

delv13
=
0.033333333333333333333*oo2dx*(-Mv3[-3*di + ijk] + 
    45.*(-Mv3[-di + ijk] + Mv3[di + ijk]) + 
    9.*(Mv3[-2*di + ijk] - Mv3[2*di + ijk]) + Mv3[3*di + ijk])
;

delv21
=
0.033333333333333333333*oo2dy*(-Mv1[-3*dj + ijk] + 
    45.*(-Mv1[-dj + ijk] + Mv1[dj + ijk]) + 
    9.*(Mv1[-2*dj + ijk] - Mv1[2*dj + ijk]) + Mv1[3*dj + ijk])
;

delv22
=
0.033333333333333333333*oo2dy*(-Mv2[-3*dj + ijk] + 
    45.*(-Mv2[-dj + ijk] + Mv2[dj + ijk]) + 
    9.*(Mv2[-2*dj + ijk] - Mv2[2*dj + ijk]) + Mv2[3*dj + ijk])
;

delv23
=
0.033333333333333333333*oo2dy*(-Mv3[-3*dj + ijk] + 
    45.*(-Mv3[-dj + ijk] + Mv3[dj + ijk]) + 
    9.*(Mv3[-2*dj + ijk] - Mv3[2*dj + ijk]) + Mv3[3*dj + ijk])
;

delv31
=
0.033333333333333333333*oo2dz*(-Mv1[-3*dk + ijk] + 
    45.*(-Mv1[-dk + ijk] + Mv1[dk + ijk]) + 
    9.*(Mv1[-2*dk + ijk] - Mv1[2*dk + ijk]) + Mv1[3*dk + ijk])
;

delv32
=
0.033333333333333333333*oo2dz*(-Mv2[-3*dk + ijk] + 
    45.*(-Mv2[-dk + ijk] + Mv2[dk + ijk]) + 
    9.*(Mv2[-2*dk + ijk] - Mv2[2*dk + ijk]) + Mv2[3*dk + ijk])
;

delv33
=
0.033333333333333333333*oo2dz*(-Mv3[-3*dk + ijk] + 
    45.*(-Mv3[-dk + ijk] + Mv3[dk + ijk]) + 
    9.*(Mv3[-2*dk + ijk] - Mv3[2*dk + ijk]) + Mv3[3*dk + ijk])
;


} else errorexit("order is not implemented yet"); 


EOS.comp("re","","","pc","","", Mrho[ijk],Mepsl[ijk], &pres_eos,&cs2); 

MturbTau11[ijk]
=
lmix*W2hrho*(-(delg111*vup1) - delg211*vup2 - delg311*vup3 - 
    2.*(delv11*g11[ijk] + delv12*g12[ijk] + delv13*g13[ijk]) + 
    g11[ijk]*(ginv22*(0.33333333333333333333*
           (delg122*vup1 + delg222*vup2 + delg322*vup3) + 
          0.66666666666666666667*
           (delv21*g12[ijk] + delv22*g22[ijk] + delv23*g23[ijk])) + 
       ginv33*(0.33333333333333333333*
           (delg133*vup1 + delg233*vup2 + delg333*vup3) + 
          0.66666666666666666667*
           (delv31*g13[ijk] + delv32*g23[ijk] + delv33*g33[ijk]))) + 
    ginv11*(g11[ijk]*(0.33333333333333333333*
           (delg111*vup1 + delg211*vup2 + delg311*vup3) + 
          0.66666666666666666667*(delv12*g12[ijk] + delv13*g13[ijk])) + 
       0.66666666666666666667*delv11*pow2(g11[ijk])) + 
    0.66666666666666666667*(ginv23*g11[ijk]*
        (delg123*vup1 + delg223*vup2 + delg323*vup3 + delv31*g12[ijk] + 
          delv21*g13[ijk] + delv32*g22[ijk] + (delv22 + delv33)*g23[ijk] + 
          delv23*g33[ijk]) + ginv12*
        (g11[ijk]*(delg112*vup1 + delg212*vup2 + delg312*vup3 + 
             (delv11 + delv22)*g12[ijk] + delv23*g13[ijk] + 
             delv12*g22[ijk] + delv13*g23[ijk]) + delv21*pow2(g11[ijk])) + 
       ginv13*(g11[ijk]*(delg113*vup1 + delg213*vup2 + delg313*vup3 + 
             delv32*g12[ijk] + (delv11 + delv33)*g13[ijk] + 
             delv12*g23[ijk] + delv13*g33[ijk]) + delv31*pow2(g11[ijk]))))*
  sqrt(cs2)
;

MturbTau12[ijk]
=
lmix*W2hrho*(-(delg112*vup1) - delg212*vup2 - delg312*vup3 - 
    delv21*g11[ijk] - (delv11 + delv22)*g12[ijk] - delv23*g13[ijk] - 
    delv12*g22[ijk] - delv13*g23[ijk] + 
    ginv33*g12[ijk]*(0.33333333333333333333*
        (delg133*vup1 + delg233*vup2 + delg333*vup3) + 
       0.66666666666666666667*(delv31*g13[ijk] + delv32*g23[ijk] + 
          delv33*g33[ijk])) + ginv22*
     (g12[ijk]*(0.33333333333333333333*
           (delg122*vup1 + delg222*vup2 + delg322*vup3) + 
          0.66666666666666666667*(delv22*g22[ijk] + delv23*g23[ijk])) + 
       0.66666666666666666667*delv21*pow2(g12[ijk])) + 
    ginv11*((0.33333333333333333333*
           (delg111*vup1 + delg211*vup2 + delg311*vup3) + 
          0.66666666666666666667*delv11*g11[ijk])*g12[ijk] + 
       0.66666666666666666667*(delv13*g12[ijk]*g13[ijk] + 
          delv12*pow2(g12[ijk]))) + 
    0.66666666666666666667*(ginv12*
        (g12[ijk]*(delg112*vup1 + delg212*vup2 + delg312*vup3 + 
             delv21*g11[ijk] + delv23*g13[ijk] + delv12*g22[ijk] + 
             delv13*g23[ijk]) + (delv11 + delv22)*pow2(g12[ijk])) + 
       ginv23*(g12[ijk]*(delg123*vup1 + delg223*vup2 + delg323*vup3 + 
             delv21*g13[ijk] + delv32*g22[ijk] + 
             (delv22 + delv33)*g23[ijk] + delv23*g33[ijk]) + 
          delv31*pow2(g12[ijk])) + 
       ginv13*(g12[ijk]*(delg113*vup1 + delg213*vup2 + delg313*vup3 + 
             delv31*g11[ijk] + (delv11 + delv33)*g13[ijk] + 
             delv12*g23[ijk] + delv13*g33[ijk]) + delv32*pow2(g12[ijk]))))*
  sqrt(cs2)
;

MturbTau13[ijk]
=
lmix*W2hrho*(-(delg113*vup1) - delg213*vup2 - delg313*vup3 - 
    delv31*g11[ijk] - delv32*g12[ijk] - (delv11 + delv33)*g13[ijk] - 
    delv12*g23[ijk] + ginv22*g13[ijk]*
     (0.33333333333333333333*
        (delg122*vup1 + delg222*vup2 + delg322*vup3) + 
       0.66666666666666666667*(delv21*g12[ijk] + delv22*g22[ijk] + 
          delv23*g23[ijk])) - delv13*g33[ijk] + 
    ginv33*(g13[ijk]*(0.33333333333333333333*
           (delg133*vup1 + delg233*vup2 + delg333*vup3) + 
          0.66666666666666666667*(delv32*g23[ijk] + delv33*g33[ijk])) + 
       0.66666666666666666667*delv31*pow2(g13[ijk])) + 
    ginv11*((0.33333333333333333333*
           (delg111*vup1 + delg211*vup2 + delg311*vup3) + 
          0.66666666666666666667*delv11*g11[ijk])*g13[ijk] + 
       0.66666666666666666667*(delv12*g12[ijk]*g13[ijk] + 
          delv13*pow2(g13[ijk]))) + 
    0.66666666666666666667*(ginv23*
        (g13[ijk]*(delg123*vup1 + delg223*vup2 + delg323*vup3 + 
             delv31*g12[ijk] + delv32*g22[ijk] + 
             (delv22 + delv33)*g23[ijk] + delv23*g33[ijk]) + 
          delv21*pow2(g13[ijk])) + 
       ginv12*(g13[ijk]*(delg112*vup1 + delg212*vup2 + delg312*vup3 + 
             delv21*g11[ijk] + (delv11 + delv22)*g12[ijk] + 
             delv12*g22[ijk] + delv13*g23[ijk]) + delv23*pow2(g13[ijk])) + 
       ginv13*(g13[ijk]*(delg113*vup1 + delg213*vup2 + delg313*vup3 + 
             delv31*g11[ijk] + delv32*g12[ijk] + delv12*g23[ijk] + 
             delv13*g33[ijk]) + (delv11 + delv33)*pow2(g13[ijk]))))*sqrt(cs2)
;

MturbTau22[ijk]
=
lmix*W2hrho*(-(delg122*vup1) - delg222*vup2 - delg322*vup3 - 
    2.*(delv21*g12[ijk] + delv22*g22[ijk] + delv23*g23[ijk]) + 
    g22[ijk]*(ginv11*(0.33333333333333333333*
           (delg111*vup1 + delg211*vup2 + delg311*vup3) + 
          0.66666666666666666667*
           (delv11*g11[ijk] + delv12*g12[ijk] + delv13*g13[ijk])) + 
       ginv33*(0.33333333333333333333*
           (delg133*vup1 + delg233*vup2 + delg333*vup3) + 
          0.66666666666666666667*
           (delv31*g13[ijk] + delv32*g23[ijk] + delv33*g33[ijk]))) + 
    ginv22*((0.33333333333333333333*
           (delg122*vup1 + delg222*vup2 + delg322*vup3) + 
          0.66666666666666666667*delv21*g12[ijk])*g22[ijk] + 
       0.66666666666666666667*(delv23*g22[ijk]*g23[ijk] + 
          delv22*pow2(g22[ijk]))) + 
    0.66666666666666666667*(ginv13*g22[ijk]*
        (delg113*vup1 + delg213*vup2 + delg313*vup3 + delv31*g11[ijk] + 
          delv32*g12[ijk] + (delv11 + delv33)*g13[ijk] + delv12*g23[ijk] + 
          delv13*g33[ijk]) + ginv12*
        (g22[ijk]*(delg112*vup1 + delg212*vup2 + delg312*vup3 + 
             delv21*g11[ijk] + (delv11 + delv22)*g12[ijk] + 
             delv23*g13[ijk] + delv13*g23[ijk]) + delv12*pow2(g22[ijk])) + 
       ginv23*(g22[ijk]*(delg123*vup1 + delg223*vup2 + delg323*vup3 + 
             delv31*g12[ijk] + delv21*g13[ijk] + 
             (delv22 + delv33)*g23[ijk] + delv23*g33[ijk]) + 
          delv32*pow2(g22[ijk]))))*sqrt(cs2)
;

MturbTau23[ijk]
=
lmix*W2hrho*(-(delg123*vup1) - delg223*vup2 - delg323*vup3 - 
    delv31*g12[ijk] - delv21*g13[ijk] - delv32*g22[ijk] - 
    (delv22 + delv33)*g23[ijk] + 
    ginv11*(0.33333333333333333333*
        (delg111*vup1 + delg211*vup2 + delg311*vup3) + 
       0.66666666666666666667*
        (delv11*g11[ijk] + delv12*g12[ijk] + delv13*g13[ijk]))*g23[ijk] - 
    delv23*g33[ijk] + ginv22*((0.33333333333333333333*
           (delg122*vup1 + delg222*vup2 + delg322*vup3) + 
          0.66666666666666666667*delv21*g12[ijk])*g23[ijk] + 
       0.66666666666666666667*(delv22*g22[ijk]*g23[ijk] + 
          delv23*pow2(g23[ijk]))) + 
    ginv33*((0.33333333333333333333*
           (delg133*vup1 + delg233*vup2 + delg333*vup3) + 
          0.66666666666666666667*delv31*g13[ijk])*g23[ijk] + 
       0.66666666666666666667*(delv33*g23[ijk]*g33[ijk] + 
          delv32*pow2(g23[ijk]))) + 
    0.66666666666666666667*(ginv13*
        (g23[ijk]*(delg113*vup1 + delg213*vup2 + delg313*vup3 + 
             delv31*g11[ijk] + delv32*g12[ijk] + 
             (delv11 + delv33)*g13[ijk] + delv13*g33[ijk]) + 
          delv12*pow2(g23[ijk])) + 
       ginv12*((delg112*vup1 + delg212*vup2 + delg312*vup3 + 
             delv21*g11[ijk] + (delv11 + delv22)*g12[ijk] + 
             delv23*g13[ijk] + delv12*g22[ijk])*g23[ijk] + 
          delv13*pow2(g23[ijk])) + 
       ginv23*(g23[ijk]*(delg123*vup1 + delg223*vup2 + delg323*vup3 + 
             delv31*g12[ijk] + delv21*g13[ijk] + delv32*g22[ijk] + 
             delv23*g33[ijk]) + (delv22 + delv33)*pow2(g23[ijk]))))*sqrt(cs2)
;

MturbTau33[ijk]
=
lmix*W2hrho*(-(delg133*vup1) - delg233*vup2 - delg333*vup3 + 
    (ginv11*(0.33333333333333333333*
           (delg111*vup1 + delg211*vup2 + delg311*vup3) + 
          0.66666666666666666667*
           (delv11*g11[ijk] + delv12*g12[ijk] + delv13*g13[ijk])) + 
       ginv22*(0.33333333333333333333*
           (delg122*vup1 + delg222*vup2 + delg322*vup3) + 
          0.66666666666666666667*
           (delv21*g12[ijk] + delv22*g22[ijk] + delv23*g23[ijk])))*g33[ijk] \
- 2.*(delv31*g13[ijk] + delv32*g23[ijk] + delv33*g33[ijk]) + 
    0.66666666666666666667*(ginv12*
        (delg112*vup1 + delg212*vup2 + delg312*vup3 + delv21*g11[ijk] + 
          (delv11 + delv22)*g12[ijk] + delv23*g13[ijk] + 
          delv12*g22[ijk] + delv13*g23[ijk])*g33[ijk] + 
       ginv13*((delg113*vup1 + delg213*vup2 + delg313*vup3 + 
             delv31*g11[ijk] + delv32*g12[ijk] + 
             (delv11 + delv33)*g13[ijk] + delv12*g23[ijk])*g33[ijk] + 
          delv13*pow2(g33[ijk])) + 
       ginv23*((delg123*vup1 + delg223*vup2 + delg323*vup3 + 
             delv31*g12[ijk] + delv21*g13[ijk] + delv32*g22[ijk] + 
             (delv22 + delv33)*g23[ijk])*g33[ijk] + delv23*pow2(g33[ijk]))) \
+ ginv33*((0.33333333333333333333*
           (delg133*vup1 + delg233*vup2 + delg333*vup3) + 
          0.66666666666666666667*delv31*g13[ijk])*g33[ijk] + 
       0.66666666666666666667*(delv32*g23[ijk]*g33[ijk] + 
          delv33*pow2(g33[ijk]))))*sqrt(cs2)
;


} else { 

MturbTau11[ijk]
=
0
;

MturbTau12[ijk]
=
0
;

MturbTau13[ijk]
=
0
;

MturbTau22[ijk]
=
0
;

MturbTau23[ijk]
=
0
;

MturbTau33[ijk]
=
0
;


} 


if (CheckForNANandINF(6, MturbTau11[ijk], MturbTau12[ijk], MturbTau13[ijk], MturbTau22[ijk], MturbTau23[ijk], MturbTau33[ijk]) ) { 

MturbTau11[ijk]
=
0
;

MturbTau12[ijk]
=
0
;

MturbTau13[ijk]
=
0
;

MturbTau22[ijk]
=
0
;

MturbTau23[ijk]
=
0
;

MturbTau33[ijk]
=
0
;


} 



/* conditional */
if (useprim) {

vshift1
=
vup1 - beta1[ijk]/alpha[ijk]
;

vshift2
=
vup2 - beta2[ijk]/alpha[ijk]
;

vshift3
=
vup3 - beta3[ijk]/alpha[ijk]
;

T00
=
(W2hrho - Mp[ijk])*pow2inv(alpha[ijk])
;

T0i1
=
(vshift1*W2hrho)/alpha[ijk] + beta1[ijk]*Mp[ijk]*pow2inv(alpha[ijk])
;

T0i2
=
(vshift2*W2hrho)/alpha[ijk] + beta2[ijk]*Mp[ijk]*pow2inv(alpha[ijk])
;

T0i3
=
(vshift3*W2hrho)/alpha[ijk] + beta3[ijk]*Mp[ijk]*pow2inv(alpha[ijk])
;

Tij11
=
MturbTau11[ijk] + W2hrho*pow2(vshift1) + 
  Mp[ijk]*(ginv11 - pow2(beta1[ijk])*pow2inv(alpha[ijk]))
;

Tij12
=
vshift1*vshift2*W2hrho + MturbTau12[ijk] + 
  Mp[ijk]*(ginv12 - beta1[ijk]*beta2[ijk]*pow2inv(alpha[ijk]))
;

Tij13
=
vshift1*vshift3*W2hrho + MturbTau13[ijk] + 
  Mp[ijk]*(ginv13 - beta1[ijk]*beta3[ijk]*pow2inv(alpha[ijk]))
;

Tij22
=
MturbTau22[ijk] + W2hrho*pow2(vshift2) + 
  Mp[ijk]*(ginv22 - pow2(beta2[ijk])*pow2inv(alpha[ijk]))
;

Tij23
=
vshift2*vshift3*W2hrho + MturbTau23[ijk] + 
  Mp[ijk]*(ginv23 - beta2[ijk]*beta3[ijk]*pow2inv(alpha[ijk]))
;

Tij33
=
MturbTau33[ijk] + W2hrho*pow2(vshift3) + 
  Mp[ijk]*(ginv33 - pow2(beta3[ijk])*pow2inv(alpha[ijk]))
;

T0j1
=
(v1*W2hrho)/alpha[ijk]
;

T0j2
=
(v2*W2hrho)/alpha[ijk]
;

T0j3
=
(v3*W2hrho)/alpha[ijk]
;


} else { /* if (!useprim) */

conD
=
MD[ijk]/sqrtdetgamma
;

conS1
=
MS1[ijk]/sqrtdetgamma
;

conS2
=
MS2[ijk]/sqrtdetgamma
;

conS3
=
MS3[ijk]/sqrtdetgamma
;

conSi1
=
conS1*ginv11 + conS2*ginv12 + conS3*ginv13
;

conSi2
=
conS1*ginv12 + conS2*ginv22 + conS3*ginv23
;

conSi3
=
conS1*ginv13 + conS2*ginv23 + conS3*ginv33
;

conT
=
Mtau[ijk]/sqrtdetgamma
;

T00
=
(conD + conT)*pow2inv(alpha[ijk])
;

T0i1
=
conSi1/alpha[ijk] - T00*beta1[ijk]
;

T0i2
=
conSi2/alpha[ijk] - T00*beta2[ijk]
;

T0i3
=
conSi3/alpha[ijk] - T00*beta3[ijk]
;

Tij11
=
conSi1*(vup1 - (2.*beta1[ijk])/alpha[ijk]) + ginv11*Mp[ijk] + 
  MturbTau11[ijk] + T00*pow2(beta1[ijk])
;

Tij12
=
beta1[ijk]*(-(conSi2/alpha[ijk]) + T00*beta2[ijk]) + 
  conSi1*(vup2 - beta2[ijk]/alpha[ijk]) + ginv12*Mp[ijk] + MturbTau12[ijk]
;

Tij13
=
beta1[ijk]*(-(conSi3/alpha[ijk]) + T00*beta3[ijk]) + 
  conSi1*(vup3 - beta3[ijk]/alpha[ijk]) + ginv13*Mp[ijk] + MturbTau13[ijk]
;

Tij22
=
conSi2*(vup2 - (2.*beta2[ijk])/alpha[ijk]) + ginv22*Mp[ijk] + 
  MturbTau22[ijk] + T00*pow2(beta2[ijk])
;

Tij23
=
beta2[ijk]*(-(conSi3/alpha[ijk]) + T00*beta3[ijk]) + 
  conSi2*(vup3 - beta3[ijk]/alpha[ijk]) + ginv23*Mp[ijk] + MturbTau23[ijk]
;

Tij33
=
conSi3*(vup3 - (2.*beta3[ijk])/alpha[ijk]) + ginv33*Mp[ijk] + 
  MturbTau33[ijk] + T00*pow2(beta3[ijk])
;

T0j1
=
conS1/alpha[ijk]
;

T0j2
=
conS2/alpha[ijk]
;

T0j3
=
conS3/alpha[ijk]
;

}
/* if (useprim) */


sqD
=
0.
;

sqT
=
sqrtmdetg*(-(dalpha1*(T0i1 + T00*beta1[ijk])) - 
    dalpha2*(T0i2 + T00*beta2[ijk]) - dalpha3*(T0i3 + T00*beta3[ijk]) + 
    2.*((Tij12 + T0i1*beta2[ijk] + beta1[ijk]*(T0i2 + T00*beta2[ijk]))*
        K12[ijk] + (Tij13 + (T0i1 + T00*beta1[ijk])*beta3[ijk])*K13[ijk] + 
       (Tij23 + T00*beta2[ijk]*beta3[ijk])*K23[ijk] + 
       T0i2*(beta2[ijk]*K22[ijk] + beta3[ijk]*K23[ijk]) + 
       T0i3*(beta1[ijk]*K13[ijk] + beta2[ijk]*K23[ijk] + 
          beta3[ijk]*K33[ijk])) - Lambda*uupt*pow2(alpha[ijk]) + 
    K11[ijk]*(Tij11 + 2.*T0i1*beta1[ijk] + T00*pow2(beta1[ijk])) + 
    K22[ijk]*(Tij22 + T00*pow2(beta2[ijk])) + 
    K33[ijk]*(Tij33 + T00*pow2(beta3[ijk])))
;

sqS1
=
sqrtmdetg*(dbeta11*T0j1 + dbeta12*T0j2 + dbeta13*T0j3 - 
    (dalpha1*T00 + Lambda*udown1)*alpha[ijk] + 
    (delg111*T0i1 + delg112*T0i2 + delg113*T0i3)*beta1[ijk] + 
    (delg112*T0i1 + delg122*T0i2 + delg123*T0i3)*beta2[ijk] + 
    (delg113*T0i1 + delg123*T0i2 + delg133*T0i3)*beta3[ijk] + 
    1.*(delg112*(Tij12 + T00*beta1[ijk]*beta2[ijk]) + 
       delg113*(Tij13 + T00*beta1[ijk]*beta3[ijk]) + 
       delg123*(Tij23 + T00*beta2[ijk]*beta3[ijk])) + 
    0.5*(delg111*(Tij11 + T00*pow2(beta1[ijk])) + 
       delg122*(Tij22 + T00*pow2(beta2[ijk])) + 
       delg133*(Tij33 + T00*pow2(beta3[ijk]))))
;

sqS2
=
sqrtmdetg*(dbeta21*T0j1 + dbeta22*T0j2 + dbeta23*T0j3 - 
    (dalpha2*T00 + Lambda*udown2)*alpha[ijk] + 
    (delg211*T0i1 + delg212*T0i2 + delg213*T0i3)*beta1[ijk] + 
    (delg212*T0i1 + delg222*T0i2 + delg223*T0i3)*beta2[ijk] + 
    (delg213*T0i1 + delg223*T0i2 + delg233*T0i3)*beta3[ijk] + 
    1.*(delg212*(Tij12 + T00*beta1[ijk]*beta2[ijk]) + 
       delg213*(Tij13 + T00*beta1[ijk]*beta3[ijk]) + 
       delg223*(Tij23 + T00*beta2[ijk]*beta3[ijk])) + 
    0.5*(delg211*(Tij11 + T00*pow2(beta1[ijk])) + 
       delg222*(Tij22 + T00*pow2(beta2[ijk])) + 
       delg233*(Tij33 + T00*pow2(beta3[ijk]))))
;

sqS3
=
sqrtmdetg*(dbeta31*T0j1 + dbeta32*T0j2 + dbeta33*T0j3 - 
    (dalpha3*T00 + Lambda*udown3)*alpha[ijk] + 
    (delg311*T0i1 + delg312*T0i2 + delg313*T0i3)*beta1[ijk] + 
    (delg312*T0i1 + delg322*T0i2 + delg323*T0i3)*beta2[ijk] + 
    (delg313*T0i1 + delg323*T0i2 + delg333*T0i3)*beta3[ijk] + 
    1.*(delg312*(Tij12 + T00*beta1[ijk]*beta2[ijk]) + 
       delg313*(Tij13 + T00*beta1[ijk]*beta3[ijk]) + 
       delg323*(Tij23 + T00*beta2[ijk]*beta3[ijk])) + 
    0.5*(delg311*(Tij11 + T00*pow2(beta1[ijk])) + 
       delg322*(Tij22 + T00*pow2(beta2[ijk])) + 
       delg333*(Tij33 + T00*pow2(beta3[ijk]))))
;

Madmrho[ijk]
=
W2hrho - Mp[ijk]
;

MadmS1[ijk]
=
vup1*W2hrho
;

MadmS2[ijk]
=
vup2*W2hrho
;

MadmS3[ijk]
=
vup3*W2hrho
;

MadmSS11[ijk]
=
g11[ijk]*Mp[ijk] + W2hrho*pow2(v1)
;

MadmSS12[ijk]
=
v1*v2*W2hrho + g12[ijk]*Mp[ijk]
;

MadmSS13[ijk]
=
v1*v3*W2hrho + g13[ijk]*Mp[ijk]
;

MadmSS22[ijk]
=
g22[ijk]*Mp[ijk] + W2hrho*pow2(v2)
;

MadmSS23[ijk]
=
v2*v3*W2hrho + g23[ijk]*Mp[ijk]
;

MadmSS33[ijk]
=
g33[ijk]*Mp[ijk] + W2hrho*pow2(v3)
;

MadmST[ijk]
=
(v1*vup1 + v2*vup2 + v3*vup3)*W2hrho + 3.*Mp[ijk]
;


if (SetSourceToZero) { 

sqD
=
0
;

sqT
=
0
;

sqS1
=
0
;

sqS2
=
0
;

sqS3
=
0
;


} 


if (Mrho[ijk]<GRHD.ATM_RHOATM) { 

sqS1
=
0
;

sqS2
=
0
;

sqS3
=
0
;


} 


if (EOS.COLD) { 

sqT
=
0
;


} 


if (CheckForNANandINF(16, Madmrho[ijk], MadmS1[ijk],MadmS2[ijk],MadmS3[ijk],   MadmSS11[ijk],MadmSS12[ijk],MadmSS13[ijk],MadmSS22[ijk],MadmSS23[ijk],MadmSS33[ijk], MadmST[ijk], sqD,sqT,sqS1,sqS2,sqS3) ) { 

sqD
=
0
;

sqT
=
0
;

sqS1
=
0
;

sqS2
=
0
;

sqS3
=
0
;

Madmrho[ijk]
=
0
;

MadmS1[ijk]
=
0
;

MadmS2[ijk]
=
0
;

MadmS3[ijk]
=
0
;

MadmSS11[ijk]
=
0
;

MadmSS12[ijk]
=
0
;

MadmSS13[ijk]
=
0
;

MadmSS22[ijk]
=
0
;

MadmSS23[ijk]
=
0
;

MadmSS33[ijk]
=
0
;

MadmST[ijk]
=
0
;


} 


} 



/* conditional */
if (addlinear) {

nMD[ijk]
=
c*sqD + nMD[ijk]
;

nMtau[ijk]
=
c*sqT + nMtau[ijk]
;

nMS1[ijk]
=
c*sqS1 + nMS1[ijk]
;

nMS2[ijk]
=
c*sqS2 + nMS2[ijk]
;

nMS3[ijk]
=
c*sqS3 + nMS3[ijk]
;


} else { /* if (!addlinear) */

nMD[ijk]
=
sqD + nMD[ijk]
;

nMtau[ijk]
=
sqT + nMtau[ijk]
;

nMS1[ijk]
=
sqS1 + nMS1[ijk]
;

nMS2[ijk]
=
sqS2 + nMS2[ijk]
;

nMS3[ijk]
=
sqS3 + nMS3[ijk]
;

}
/* if (addlinear) */



if (CheckForNANandINF(16, Madmrho[ijk], MadmS1[ijk],MadmS2[ijk],MadmS3[ijk], MadmSS11[ijk],MadmSS12[ijk],MadmSS13[ijk],MadmSS22[ijk],MadmSS23[ijk],MadmSS33[ijk], MadmST[ijk], sqD,sqT,sqS1,sqS2,sqS3)) {  


 printf("problem with nan's inside con2source\n");  


 printf("  x=%e y=%e z=%e\n",Ptr(level,"x")[ijk],Ptr(level,"y")[ijk],Ptr(level,"z")[ijk]);  


 printf("  sD=%e sT=%e sSx=%e sSy=%e sSz=%e\n",sqD,sqT,sqS1,sqS2,sqS3);  


 printf("  %e %e %e %e   %e %e %e  %e\n",v1,v2,v3,vsq,  W2hrho,g11[ijk],Mp[ijk],  MadmSS12[ijk]);  


 if (ijkinsidefinerlevel(box,ijk)==0) printf("  point is NOT inside finer box NOR in some symmetry area\n");  


 else printf("  point is inside finer box/ in symmetry \n");  


 }  



} endfor_ijk_openmp; /* loop i, j, k */



bampi_openmp_stop


}  /* function */

/* grhd_sources_turb.c */
/* nvars = 59, nauxs = 100, n* = 1413,  n/ = 77,  n+ = 1809, n = 3299, O = 0 */
