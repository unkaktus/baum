/* compute_magvars.c */
/* A Neuweiler 01/23 */

#include "bam.h"
#include "hydroanalysis.h"

#define PR 0

#define Power(x,y) pow((double) (x), (double) (y))
#define Sqrt(x)    sqrt((double) (x))
#define Log(x)     log((double) (x))
#define pow2(x)    ((x)*(x))
#define pow4(x)    ((x)*(x)*(x)*(x))
#define pow2inv(x) (1.0/((x)*(x)))

#define Tanh(x)    tanh(x)
#define Sech(x)    (1/cosh(x))
#define Cal(x,y,z) ((x)?(y):(z))


#define ALPHAMIN 1e-6


/* perform MHD analysis
   registered in ANALYZE */
void compute_magvars(tL *level)
{

  if (PR) printf(" \t in compute_magvars()\n");

  double *Bx   = Ptr(level, "grmhd_Bx");
  double *By   = Ptr(level, "grmhd_By");
  double *Bz   = Ptr(level, "grmhd_Bz");
  double *bx   = Ptr(level, "grmhd_sbx");
  double *by   = Ptr(level, "grmhd_sby");
  double *bz   = Ptr(level, "grmhd_sbz");
  double *b2   = Ptr(level, "grmhd_sb2");
  double *b0   = Ptr(level, "grmhd_sb0");
  double *vx   = Ptr(level, "grhd_vx");
  double *vy   = Ptr(level, "grhd_vy");
  double *vz   = Ptr(level, "grhd_vz");
  double *v2   = Ptr(level, "grhd_v2");

  double *Emag = Ptr(level, "hydroa_Emag");
  double *Pmag = Ptr(level, "hydroa_Pmag");
  double *Bnorm = Ptr(level, "hydroa_B");
  double *divB = Ptr(level, "hydroa_divB");
  double *Emagtor = Ptr(level, "hydroa_Emagtor");
  double *Emagpol = Ptr(level, "hydroa_Emagpol");
  double *Btor = Ptr(level, "hydroa_Btor");
  double *Bpol = Ptr(level, "hydroa_Bpol");
  double *Fpoy = Ptr(level, "hydroa_Fpoy");

  double *mask = Ptr(level, "matter_mask");
  double *gxx   = Ptr(level, "adm_gxx");
  double *gxy   = Ptr(level, "adm_gxy");
  double *gxz   = Ptr(level, "adm_gxz");
  double *gyy   = Ptr(level, "adm_gyy");
  double *gyz   = Ptr(level, "adm_gyz");
  double *gzz   = Ptr(level, "adm_gzz");
  double *alpha = Ptr(level, "alpha");
  double *betax = Ptr(level, "betax");
  double *betay = Ptr(level, "betay");
  double *betaz = Ptr(level, "betaz");

  double *xp = level->v[Ind("x")];
  double *yp = level->v[Ind("y")];
  double *zp = level->v[Ind("z")];

  int order = Geti("order_centered");

  double dx = level->dx;
  double dy = level->dy;
  double dz = level->dz;
  double oo2dx = 1/(2*dx), oo2dy = 1/(2*dy), oo2dz = 1/(2*dz);

  bampi_openmp_start

  double sqrtg, W, B2;

  double dBxdx, dBydy, dBzdz;
  double detgamma, sqrtgamma;

  double Bxtor, Bytor, Bztor, Bxpol, Bypol, Bzpol;
  double B2tor, B2pol, sb2tor, sb2pol;
  double vr, br, ut, bt;

  forinnerpoints_ijk_openmp(level) {
    detgamma        = detg(gxx[ijk],gxy[ijk],gxz[ijk],gyy[ijk],gyz[ijk],gzz[ijk]); 
    sqrtgamma       = sqrt(detgamma);
    sqrtg     = sqrtgamma*alpha[ijk];
    W         = 1/sqrt(1. - v2[ijk]);

    if ((MATTER.USEMASK) && (mask[ijk]>.99)) { 
      Emag[ijk] = 0.;
      Pmag[ijk] = 0.;
      Bnorm[ijk]= 0.;
      divB[ijk] = 0.;
      Emagtor[ijk] = 0.;
      Emagpol[ijk] = 0.;
      Btor[ijk] = 0;
      Bpol[ijk] = 0;
      Fpoy[ijk] = 0;
      continue; 
    } 


    B2 =    Bx[ijk]*(gxx[ijk]*Bx[ijk] + gxy[ijk]*By[ijk] + gxz[ijk]*Bz[ijk]) 
          + By[ijk]*(gxy[ijk]*Bx[ijk] + gyy[ijk]*By[ijk] + gyz[ijk]*Bz[ijk])
          + Bz[ijk]*(gxz[ijk]*Bx[ijk] + gyz[ijk]*By[ijk] + gzz[ijk]*Bz[ijk]);

    /* compute 2nd order derivatives */
    dBxdx = oo2dx*(-Bx[-di + ijk] + Bx[di + ijk]);
    dBydy = oo2dy*(-By[-dj + ijk] + By[dj + ijk]);
    dBzdz = oo2dz*(-Bz[-dk + ijk] + Bz[dk + ijk]);

    Emag[ijk]   = 0.5*W*sqrtgamma*b2[ijk] ; //1./(8.*PI)*W*sqrtgamma*b2[ijk] ;
    Pmag[ijk]   = 0.5*b2[ijk];
    Bnorm[ijk]  = sqrt(B2)/sqrtgamma;
    divB[ijk]   = -1./sqrtg*(dBxdx + dBydy + dBzdz);

    /* toroidal and poloidal components */
    Bxtor = - yp[ijk]/(xp[ijk]*xp[ijk]+yp[ijk]*yp[ijk])* (-yp[ijk]*Bx[ijk]+xp[ijk]*By[ijk]);
    Bytor =   xp[ijk]/(xp[ijk]*xp[ijk]+yp[ijk]*yp[ijk])* (-yp[ijk]*Bx[ijk]+xp[ijk]*By[ijk]);
    Bztor =   0.;

    B2tor = Bxtor*(gxx[ijk]*Bxtor + gxy[ijk]*Bytor + gxz[ijk]*Bztor) 
          + Bytor*(gxy[ijk]*Bxtor + gyy[ijk]*Bytor + gyz[ijk]*Bztor)
          + Bztor*(gxz[ijk]*Bxtor + gyz[ijk]*Bytor + gzz[ijk]*Bztor);

    Bxpol =   xp[ijk]/(xp[ijk]*xp[ijk]+yp[ijk]*yp[ijk])* ( xp[ijk]*Bx[ijk]+yp[ijk]*By[ijk]);
    Bypol =   yp[ijk]/(xp[ijk]*xp[ijk]+yp[ijk]*yp[ijk])* ( xp[ijk]*Bx[ijk]+yp[ijk]*By[ijk]);
    Bzpol =   Bz[ijk];

    B2pol = Bxpol*(gxx[ijk]*Bxpol + gxy[ijk]*Bypol + gxz[ijk]*Bzpol) 
          + Bypol*(gxy[ijk]*Bxpol + gyy[ijk]*Bypol + gyz[ijk]*Bzpol)
          + Bzpol*(gxz[ijk]*Bxpol + gyz[ijk]*Bypol + gzz[ijk]*Bzpol) ;

    sb2tor =  (B2tor/W/W  
              + pow((Bxtor*(gxx[ijk]*vx[ijk]+gxy[ijk]*vy[ijk]+gxz[ijk]*vz[ijk])
                  +  Bytor*(gxy[ijk]*vx[ijk]+gyy[ijk]*vy[ijk]+gyz[ijk]*vz[ijk])
                  +  Bztor*(gxz[ijk]*vx[ijk]+gyz[ijk]*vy[ijk]+gzz[ijk]*vz[ijk])),2.))/sqrtgamma/sqrtgamma;

    sb2pol =  (B2pol/W/W  
              + pow((Bxpol*(gxx[ijk]*vx[ijk]+gxy[ijk]*vy[ijk]+gxz[ijk]*vz[ijk])
                  +  Bypol*(gxy[ijk]*vx[ijk]+gyy[ijk]*vy[ijk]+gyz[ijk]*vz[ijk])
                  +  Bzpol*(gxz[ijk]*vx[ijk]+gyz[ijk]*vy[ijk]+gzz[ijk]*vz[ijk])),2.))/sqrtgamma/sqrtgamma;

    Btor[ijk] = sqrt(B2tor)/sqrtgamma;
    Bpol[ijk] = sqrt(B2pol)/sqrtgamma;

    Emagtor[ijk]   = 0.5*W*sqrtgamma*sb2tor ; // 1./(8.*PI)*W*sqrtgamma*sb2tor ;
    Emagpol[ijk]   = 0.5*W*sqrtgamma*sb2pol ; // 1./(8.*PI)*W*sqrtgamma*sb2pol ;

    vr = (vx[ijk]*xp[ijk]+vy[ijk]*yp[ijk]+vz[ijk]*zp[ijk])/sqrt(xp[ijk]*xp[ijk]+yp[ijk]*yp[ijk]+zp[ijk]*zp[ijk]);
    br = (bx[ijk]*xp[ijk]+by[ijk]*yp[ijk]+bz[ijk]*zp[ijk])/sqrt(xp[ijk]*xp[ijk]+yp[ijk]*yp[ijk]+zp[ijk]*zp[ijk]);
    ut = (-alpha[ijk] + (betax[ijk]*vx[ijk]+betay[ijk]*vy[ijk]+betaz[ijk]*vz[ijk]))*W;
    bt = (-alpha[ijk] + (betax[ijk]*vx[ijk]+betay[ijk]*vy[ijk]+betaz[ijk]*vz[ijk]))*b0[ijk]*alpha[ijk] 
          + (betax[ijk]*Bx[ijk] + betay[ijk]*By[ijk] + betaz[ijk]*Bz[ijk])/W;
    Fpoy[ijk] = - sqrtg*(b2[ijk]*vr*ut - br*bt);

  } endfor_ijk_openmp; /* loop i, j, k */
  
  bampi_openmp_stop

}
