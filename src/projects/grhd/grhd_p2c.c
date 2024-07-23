/* grhd_p2c 
   sbernuz 03/2012 */

#include "bam.h"
#include "grhd.h"

#define PR 0


// compute v_i and v^2
void grhd_compute_v2_pt(double vx, double vy, double vz,
			double gxx, double gxy, double gxz, double gyy, double gyz, double gzz,
			double *vlowx, double *vlowy, double *vlowz, double *v2 )
{
  
  double vlx = gxx*vx + gxy*vy + gxz*vz; 
  double vly = gxy*vx + gyy*vy + gyz*vz; 
  double vlz = gxz*vx + gyz*vy + gzz*vz;

  *vlowx = vlx;
  *vlowy = vly;
  *vlowz = vlz;
  *v2    = vlx*vx + vly*vy + vlz*vz;

}


// compute detg and diag of inverse metric (safe computation)
void grhd_compute_detg_invg_pt(double g11, double g12, double g13, 
			       double g22, double g23, double g33,
			       double *det,
			       double *i11, double *i12, double *i13, 
			       double *i22, double *i23, double *i33)
{
  
  double detg, oodetg, gginv11, gginv12, gginv13, gginv22, gginv23, gginv33;
  
  gginv11 = g22*g33 - g23*g23;
  gginv12 = g13*g23 - g12*g33;
  gginv13 = g12*g23 - g13*g22;
  gginv22 = g11*g33 - g13*g13;
  gginv23 = g12*g13 - g11*g23;
  gginv33 = g11*g22 - g12*g12;
  
  detg    = g11*gginv11 + g12*gginv12 + g13*gginv13;

  // make sure metric is always invertible (useful? testme)
  /*
    if (!finite(detg)) {
    if (PR) printf("metric inversion: detg not finite, set to one and go on.\n");
    detg = 1.; 
    }
    detg   = DMAX( detg , DETGMIN );
  */
  oodetg = 1./detg;

  *i11 = gginv11 * oodetg;
  *i12 = gginv12 * oodetg;
  *i13 = gginv13 * oodetg;
  *i22 = gginv22 * oodetg;
  *i23 = gginv23 * oodetg;
  *i33 = gginv33 * oodetg;
  *det = detg;
}


// conservatives
void grhd_compute_q_pt(double gxx, double gxy, double gxz, double gyy, double gyz, double gzz, double detg, 
		       double rho, double epsl, double p, double vlowx, double vlowy, double vlowz, double Wlor,
		       double *D, double *T, double *Sx, double *Sy, double *Sz)
{

  double sqrtdetgamma = sqrt(detg);
  /* double sqrtdetgamma = sqrt(fabs(detg)); */
  double W2hrho = Wlor*Wlor*( rho + rho*epsl + p );
  
  *D   = sqrtdetgamma * Wlor   * rho;
  *T   = sqrtdetgamma * (W2hrho - p - Wlor*rho);
  *Sx  = sqrtdetgamma * W2hrho * vlowx ;
  *Sy  = sqrtdetgamma * W2hrho * vlowy ;
  *Sz  = sqrtdetgamma * W2hrho * vlowz ;
  
  if (PR) printf("%2.6e %2.6e %2.6e %2.6e %2.6e + %2.6e %2.6e   ->  %2.6e %2.6e %2.6e %2.6e %2.6e\n",
      rho,epsl,vlowx,vlowy,vlowz, sqrtdetgamma,W2hrho, *D,*T,*Sx,*Sy,*Sz);

}


// phys fluxes, x direction
void grhd_compute_fx_pt(double alpha, double betax, double detg, 
			double pres, double vx, 
			double D, double T, double Sx, double Sy, double Sz,
			double* fD, double* fT, double* fSx, double* fSy, double* fSz)
{
  
  double sqrtg = sqrt(detg);
  /* double sqrtg = sqrt(fabs(detg)); */
  double avb   = alpha * vx - betax;
  
  *fD  = D  * avb;
  *fT  = T  * avb + alpha * sqrtg * pres * vx;    
  *fSx = Sx * avb + alpha * sqrtg * pres;  
  *fSy = Sy * avb;
  *fSz = Sz * avb;

}


// add gamma^1k*tau_ki to S-fluxes (for turbulence)
void grhd_addturb_fx_pt(double tTau11, double tTau12, double tTau13, double tTau22, double tTau23, double tTau33,
			double alpha, double detg, double gup11, double gup12, double gup13,
			double* fSx, double* fSy, double* fSz) {
	     	
	double sqrtg = sqrt(detg);

	double gtaux, gtauy, gtauz;
	
	
	gtaux = gup11* tTau11 + gup12* tTau12 + gup13* tTau13;
	gtauy = gup11* tTau12 + gup12* tTau22 + gup13* tTau23;
	gtauz = gup11* tTau13 + gup12* tTau23 + gup13* tTau33;
	
	*fSx += alpha* sqrtg* gtaux;
	*fSy += alpha* sqrtg* gtauy;
	*fSz += alpha* sqrtg* gtauz;
	
}
