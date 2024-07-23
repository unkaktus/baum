/* grhd_p2a.c */
/* Copyright (C) 1998 Bernd Bruegmann, 15.8.2012 */
/* Produced with Mathematica */

#include "bam.h"
#include "grhd.h"

#define Power(x,y) (pow((double) (x), (double) (y)))
#define Log(x)     log((double) (x)))
#define pow2(x)    ((x)*(x))
#define pow2inv(x) (1.0/((x)*(x)))

#define Cal(x,y,z) ((x)?(y):(z))




void grhd_p2a(tL *level)
{

int index_adm_gxx = Ind("adm_gxx");
double *g11 = level->v[index_adm_gxx + 0];
double *g12 = level->v[index_adm_gxx + 1];
double *g13 = level->v[index_adm_gxx + 2];
double *g22 = level->v[index_adm_gxx + 3];
double *g23 = level->v[index_adm_gxx + 4];
double *g33 = level->v[index_adm_gxx + 5];
tVarList *vl_give_matterADMVars_vl = give_matterADMVars_vl(level);
double *Madmrho = vldataptr(vl_give_matterADMVars_vl, 0);
double *MadmS1 = vldataptr(vl_give_matterADMVars_vl, 1);
double *MadmS2 = vldataptr(vl_give_matterADMVars_vl, 2);
double *MadmS3 = vldataptr(vl_give_matterADMVars_vl, 3);
double *MadmSS11 = vldataptr(vl_give_matterADMVars_vl, 4);
double *MadmSS12 = vldataptr(vl_give_matterADMVars_vl, 5);
double *MadmSS13 = vldataptr(vl_give_matterADMVars_vl, 6);
double *MadmSS22 = vldataptr(vl_give_matterADMVars_vl, 7);
double *MadmSS23 = vldataptr(vl_give_matterADMVars_vl, 8);
double *MadmSS33 = vldataptr(vl_give_matterADMVars_vl, 9);
double *MadmST = vldataptr(vl_give_matterADMVars_vl, 10);
vlfree(vl_give_matterADMVars_vl);
tVarList *vl_give_primVars_vl = give_primVars_vl(level);
double *Mrho = vldataptr(vl_give_primVars_vl, 0);
double *Mepsl = vldataptr(vl_give_primVars_vl, 1);
double *Mv1 = vldataptr(vl_give_primVars_vl, 2);
double *Mv2 = vldataptr(vl_give_primVars_vl, 3);
double *Mv3 = vldataptr(vl_give_primVars_vl, 4);
double *Mp = vldataptr(vl_give_primVars_vl, 5);
double *Mvsqr = vldataptr(vl_give_primVars_vl, 6);
double *detmetric = vldataptr(vl_give_primVars_vl, 7);
vlfree(vl_give_primVars_vl);




bampi_openmp_start


double dx = level->dx;
double dy = level->dy;
double dz = level->dz;
double oo2dx = 1/(2*dx), oo2dy = 1/(2*dy), oo2dz = 1/(2*dz);
double oodx2 = 1/(dx*dx), oody2 = 1/(dy*dy), oodz2 = 1/(dz*dz);
double oo4dxdy = 1/(4*dx*dy);
double oo4dxdz = 1/(4*dx*dz);
double oo4dydz = 1/(4*dy*dz);

double v1 = 0.;
double v2 = 0.;
double v3 = 0.;
double vsq = 0.;
double vup1 = 0.;
double vup2 = 0.;
double vup3 = 0.;
double W = 0.;
double W2hrho = 0.;



forinnerpoints_ijk_openmp(level) {


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


if (vsq>=GRHD.HRSC_VMAX) W = GRHD.HRSC_WlorMAX; 

W2hrho
=
(Mp[ijk] + Mrho[ijk] + Mepsl[ijk]*Mrho[ijk])*pow2(W)
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



} endfor_ijk_openmp; /* loop i, j, k */



bampi_openmp_stop


}  /* function */

/* grhd_p2a.c */
/* nvars = 24, nauxs = 9, n* = 89,  n/ = 21,  n+ = 36, n = 146, O = 0 */
