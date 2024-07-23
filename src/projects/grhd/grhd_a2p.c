/* grhd_a2p.c */
/* Copyright (C) 1998 Bernd Bruegmann, 15.8.2012 */
/* Produced with Mathematica */

#include "bam.h"
#include "grhd.h"

#define Power(x,y) (pow((double) (x), (double) (y)))
#define Log(x)     log((double) (x)))
#define pow2(x)    ((x)*(x))
#define pow2inv(x) (1.0/((x)*(x)))

#define Cal(x,y,z) ((x)?(y):(z))




void grhd_a2p(tL* level)
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




forinnerpoints_ijk_openmp(level) {


Mv1[ijk]
=
0
;

Mv2[ijk]
=
0
;

Mv3[ijk]
=
0
;



} endfor_ijk_openmp; /* loop i, j, k */



bampi_openmp_stop


}  /* function */

/* grhd_a2p.c */
/* nvars = 24, nauxs = 0, n* = 52,  n/ = 20,  n+ = 15, n = 87, O = 0 */
