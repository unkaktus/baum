/* bssn_rhs_add_source.c */
/* Copyright (C) 1998 Bernd Bruegmann, 28.5.2024 */
/* Produced with Mathematica */

#include "bam.h"
#include "bssn.h"

#define Power(x,y) (pow((double) (x), (double) (y)))
#define Log(x)     log((double) (x)))
#define pow2(x)    ((x)*(x))
#define pow2inv(x) (1.0/((x)*(x)))

#define Cal(x,y,z) ((x)?(y):(z))




void bssn_rhs_add_source(tVarList *unew, tVarList *upre, double c, tVarList *ucur)
{

tL *level = ucur->level;
int addlinear = (c != 0.0l);

const double chiDivFloor = Getd("bssn_chi_div_floor");
const double chipsipower = Getd("bssn_chi_psipower");

int index_adm_rho = Ind("adm_rho");
double *Madmrho = level->v[index_adm_rho + 0];
int index_adm_Sx = Ind("adm_Sx");
double *MadmS1 = level->v[index_adm_Sx + 0];
double *MadmS2 = level->v[index_adm_Sx + 1];
double *MadmS3 = level->v[index_adm_Sx + 2];
int index_adm_SSxx = Ind("adm_SSxx");
double *MadmSS11 = level->v[index_adm_SSxx + 0];
double *MadmSS12 = level->v[index_adm_SSxx + 1];
double *MadmSS13 = level->v[index_adm_SSxx + 2];
double *MadmSS22 = level->v[index_adm_SSxx + 3];
double *MadmSS23 = level->v[index_adm_SSxx + 4];
double *MadmSS33 = level->v[index_adm_SSxx + 5];
int index_adm_ST = Ind("adm_ST");
double *MadmST = level->v[index_adm_ST + 0];
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
const double *xp = level->v[Ind("x")];
const double *yp = level->v[Ind("y")];
const double *zp = level->v[Ind("z")];
const double gamma0factor    = Getv("bssn_shift", "gamma0");
const double gamma2factor    = Getv("bssn_shift", "gamma2");
const double withGadv        = Getv("bssn_shift", "withGadv");
const double withShiftadv    = Getv("bssn_shift", "withShiftadv");
const double withBadv        = Getv("bssn_shift", "withBadv");
const double withB           = !Getv("bssn_shift", "withoutB");
const double lapseharmonicf  = Getd("bssn_lapseharmonicf");
const double shiftalphapower = Getd("bssn_shiftalphapower");
const double shiftgammacoeff = Getd("bssn_shiftgammacoeff");
const double shiftdriver     = Getd("bssn_shiftdriver");
const int setRHSto0 = Getv("bssnRHSto0","yes")+ Getv("bssnsourceRHSto0","yes");




bampi_openmp_start


double dx = level->dx;
double dy = level->dy;
double dz = level->dz;
double oo2dx = 1/(2*dx), oo2dy = 1/(2*dy), oo2dz = 1/(2*dz);
double oodx2 = 1/(dx*dx), oody2 = 1/(dy*dy), oodz2 = 1/(dz*dz);
double oo4dxdy = 1/(4*dx*dy);
double oo4dxdz = 1/(4*dx*dz);
double oo4dydz = 1/(4*dy*dz);

double betaF = 0.;
double chiguard = 0.;
double chiguarded = 0.;
double detg = 0.;
double detginv = 0.;
double f = 0.;
double ff = 0.;
double ginv11 = 0.;
double ginv12 = 0.;
double ginv13 = 0.;
double ginv22 = 0.;
double ginv23 = 0.;
double ginv33 = 0.;
double MadmSSTF11 = 0.;
double MadmSSTF12 = 0.;
double MadmSSTF13 = 0.;
double MadmSSTF22 = 0.;
double MadmSSTF23 = 0.;
double MadmSSTF33 = 0.;
double metric11 = 0.;
double metric12 = 0.;
double metric13 = 0.;
double metric22 = 0.;
double metric23 = 0.;
double metric33 = 0.;
double metricinv11 = 0.;
double metricinv12 = 0.;
double metricinv13 = 0.;
double metricinv22 = 0.;
double metricinv23 = 0.;
double metricinv33 = 0.;
double oochipsipower = 0.;
double psi4 = 0.;
double psim4 = 0.;
double r2A11 = 0.;
double r2A12 = 0.;
double r2A13 = 0.;
double r2A22 = 0.;
double r2A23 = 0.;
double r2A33 = 0.;
double rA11 = 0.;
double rA12 = 0.;
double rA13 = 0.;
double rA22 = 0.;
double rA23 = 0.;
double rA33 = 0.;
double rB1 = 0.;
double rB2 = 0.;
double rB3 = 0.;
double rG1 = 0.;
double rG2 = 0.;
double rG3 = 0.;
double rK = 0.;



forinnerpoints_ijk_openmp(level) {


detg
=
2.*g12[ijk]*g13[ijk]*g23[ijk] + g11[ijk]*g22[ijk]*g33[ijk] - 
  g33[ijk]*pow2(g12[ijk]) - g22[ijk]*pow2(g13[ijk]) - g11[ijk]*pow2(g23[ijk])
;

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

chiguard
=
chiDivFloor
;

chiguarded
=
chi[ijk]
;


if (chiguarded < chiguard) chiguarded = chiguard; 

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

psi4
=
exp(4.*f)
;

metric11
=
psi4*g11[ijk]
;

metric12
=
psi4*g12[ijk]
;

metric13
=
psi4*g13[ijk]
;

metric22
=
psi4*g22[ijk]
;

metric23
=
psi4*g23[ijk]
;

metric33
=
psi4*g33[ijk]
;

metricinv11
=
ginv11*psim4
;

metricinv12
=
ginv12*psim4
;

metricinv13
=
ginv13*psim4
;

metricinv22
=
ginv22*psim4
;

metricinv23
=
ginv23*psim4
;

metricinv33
=
ginv33*psim4
;


//if (detg<=0.) errorexit(" detg <= 0.");  


//if (psi4<=0.) errorexit(" psi^4 <= 0.");  

MadmSSTF11
=
MadmSS11[ijk] - 0.33333333333333333333*metric11*MadmST[ijk]
;

MadmSSTF12
=
MadmSS12[ijk] - 0.33333333333333333333*metric12*MadmST[ijk]
;

MadmSSTF13
=
MadmSS13[ijk] - 0.33333333333333333333*metric13*MadmST[ijk]
;

MadmSSTF22
=
MadmSS22[ijk] - 0.33333333333333333333*metric22*MadmST[ijk]
;

MadmSSTF23
=
MadmSS23[ijk] - 0.33333333333333333333*metric23*MadmST[ijk]
;

MadmSSTF33
=
MadmSS33[ijk] - 0.33333333333333333333*metric33*MadmST[ijk]
;

rA11
=
-8.*MadmSSTF11*PI*psim4*alpha[ijk]
;

rA12
=
-8.*MadmSSTF12*PI*psim4*alpha[ijk]
;

rA13
=
-8.*MadmSSTF13*PI*psim4*alpha[ijk]
;

rA22
=
-8.*MadmSSTF22*PI*psim4*alpha[ijk]
;

rA23
=
-8.*MadmSSTF23*PI*psim4*alpha[ijk]
;

rA33
=
-8.*MadmSSTF33*PI*psim4*alpha[ijk]
;

rK
=
4.*PI*alpha[ijk]*(Madmrho[ijk] + MadmST[ijk])
;

rG1
=
-16.*PI*psi4*alpha[ijk]*MadmS1[ijk]
;

rG2
=
-16.*PI*psi4*alpha[ijk]*MadmS2[ijk]
;

rG3
=
-16.*PI*psi4*alpha[ijk]*MadmS3[ijk]
;

r2A11
=
-5.3333333333333333333*PI*alpha[ijk]*g11[ijk]*Madmrho[ijk]
;

r2A12
=
-5.3333333333333333333*PI*alpha[ijk]*g12[ijk]*Madmrho[ijk]
;

r2A13
=
-5.3333333333333333333*PI*alpha[ijk]*g13[ijk]*Madmrho[ijk]
;

r2A22
=
-5.3333333333333333333*PI*alpha[ijk]*g22[ijk]*Madmrho[ijk]
;

r2A23
=
-5.3333333333333333333*PI*alpha[ijk]*g23[ijk]*Madmrho[ijk]
;

r2A33
=
-5.3333333333333333333*PI*alpha[ijk]*g33[ijk]*Madmrho[ijk]
;

betaF
=
shiftgammacoeff*Power(alpha[ijk],shiftalphapower)
;

rB1
=
(gamma0factor + betaF*gamma2factor)*rG1*withB
;

rB2
=
(gamma0factor + betaF*gamma2factor)*rG2*withB
;

rB3
=
(gamma0factor + betaF*gamma2factor)*rG3*withB
;


if (setRHSto0 || CheckForNANandINF(19, rA11,rA12,rA13,rA22,rA23,rA33,    
    r2A11,r2A12,r2A13,r2A22,r2A23,r2A33, rG1,rG2,rG3,rK, rB1,rB2,rB3)) {
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

rK
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

r2A11
=
0
;

r2A12
=
0
;

r2A13
=
0
;

r2A22
=
0
;

r2A23
=
0
;

r2A33
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


} 



/* conditional */
if (addlinear) {

nK[ijk]
=
c*rK + nK[ijk]
;

nA11[ijk]
=
c*(r2A11 + rA11) + nA11[ijk]
;

nA12[ijk]
=
c*(r2A12 + rA12) + nA12[ijk]
;

nA13[ijk]
=
c*(r2A13 + rA13) + nA13[ijk]
;

nA22[ijk]
=
c*(r2A22 + rA22) + nA22[ijk]
;

nA23[ijk]
=
c*(r2A23 + rA23) + nA23[ijk]
;

nA33[ijk]
=
c*(r2A33 + rA33) + nA33[ijk]
;

nG1[ijk]
=
c*rG1 + nG1[ijk]
;

nG2[ijk]
=
c*rG2 + nG2[ijk]
;

nG3[ijk]
=
c*rG3 + nG3[ijk]
;

nB1[ijk]
=
c*rB1 + nB1[ijk]
;

nB2[ijk]
=
c*rB2 + nB2[ijk]
;

nB3[ijk]
=
c*rB3 + nB3[ijk]
;


} else { /* if (!addlinear) */

nK[ijk]
=
rK + nK[ijk]
;

nA11[ijk]
=
r2A11 + rA11 + nA11[ijk]
;

nA12[ijk]
=
r2A12 + rA12 + nA12[ijk]
;

nA13[ijk]
=
r2A13 + rA13 + nA13[ijk]
;

nA22[ijk]
=
r2A22 + rA22 + nA22[ijk]
;

nA23[ijk]
=
r2A23 + rA23 + nA23[ijk]
;

nA33[ijk]
=
r2A33 + rA33 + nA33[ijk]
;

nG1[ijk]
=
rG1 + nG1[ijk]
;

nG2[ijk]
=
rG2 + nG2[ijk]
;

nG3[ijk]
=
rG3 + nG3[ijk]
;

nB1[ijk]
=
rB1 + nB1[ijk]
;

nB2[ijk]
=
rB2 + nB2[ijk]
;

nB3[ijk]
=
rB3 + nB3[ijk]
;

}
/* if (addlinear) */



if (CheckForNANandINF(13, nA11[ijk],nA12[ijk],nA13[ijk],nA22[ijk],nA23[ijk],nA33[ijk], nG1[ijk],nG2[ijk],nG3[ijk], nK[ijk], nB1[ijk],nB2[ijk],nB3[ijk])) 
    printf("the ultimative test if there are no nans coming from the matter part ... failed (nans inside bssn_rhs_add_source)\n");


} endfor_ijk_openmp; /* loop i, j, k */



bampi_openmp_stop


}  /* function */

/* bssn_rhs_add_source.c */
/* nvars = 62, nauxs = 53, n* = 232,  n/ = 32,  n+ = 151, n = 415, O = 0 */
