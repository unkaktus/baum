/* grhd_eig.c 
   sbernuz 11/2011 */

#include "bam.h"
#include "grhd.h"


#define FLATSPACETIME 0
#define PR 0
#define DEBUG 0

#define ONE 1.0
#define ZERO 0.
#define TWO 2.
#define BIG 1e32
#define TINY 1e-32


// grhd eigenvalues in x-direction 
void grhd_eig( double alpha, double betax, double vx, double v2, double gxxup, double cs2,
	       double *lam0, double *lamm, double *lamp )
{
  
  double Csqrt = sqrt( cs2*( ( ONE - v2  ) * ( gxxup  * ( ONE - v2 * cs2 ) - vx * vx * ( ONE-cs2  ) ) ) );
  
  // lambda0 
  *lam0 = alpha  * vx  - betax ;  

  // lambda- 
  *lamm = alpha  *( vx  *( ONE - cs2  ) - Csqrt )/( ONE - v2  * cs2  ) - betax ;  

  // lambda+ 
  *lamp = alpha  *( vx  *( ONE - cs2  ) + Csqrt )/( ONE - v2  * cs2  ) - betax ;
  
  // cure eigenvalues 
  int cure = 0;  
  if (!finite(*lamm)) {
    cure++;
    *lamm = 0.;
  }
  if (!finite(*lamp)) {
    cure++;
    *lamp = 0.;
  }
  if (!finite(*lam0)) {
    cure++;
    *lam0 = 0.;
  }  
  
#if (FLATSPACETIME)
  // check, 1d flat spacetime:
  Csqrt  = sqrt(cs2);
  *lam0  = vx;
  *lamm  = (vx - Csqrt)/(ONE - vx*Csqrt); // lam-
  *lamp  = (vx + Csqrt)/(ONE + vx*Csqrt); // lam+
#endif
  
  if (PR) if (cure>0) printf("(PR=%d) cured %d eigenvalues \n",PR,cure);

}

/* WT: same as grhd_eig, but without NANs */
/* grhd eigenvalues in x-direction */
void grhd_eigenvals(double alpha, double betax, double vx, double v2,
                    double gxxup, double cs2,
	            double *lam0, double *lamm, double *lamp )
{
  int cure = 0;  
  double ONE_minus_v2cs2, C, Csqrt;

  // lambda0 
  *lam0 = alpha*vx - betax;

  /* terms in lamm, lamp */
  ONE_minus_v2cs2 = ONE - v2*cs2;
  C = cs2*( (ONE-v2) * ( gxxup*(ONE_minus_v2cs2) - vx*vx*(ONE-cs2) ) );

  /* check for div by 0 and sqrt of neg. numbers in lamm, lamp */
  if(ONE_minus_v2cs2==0.0) cure = cure | 1;      /* set bit0 */
  if(C<0.0) { C=fabs(cs2); cure = cure | 2; }    /* set bit1 */
  Csqrt = sqrt(C);

  // cure eigenvalues 
  if(cure)
    *lamm = *lamp = 0.0; /* we could also try *lamm = *lamp = *lam0; */
  else
  {
    // lambda-
    *lamm = alpha*( vx*(ONE-cs2) - Csqrt )/(ONE_minus_v2cs2) - betax;
    // lambda+ 
    *lamp = alpha*( vx*(ONE-cs2) + Csqrt )/(ONE_minus_v2cs2) - betax;
  }
  
#if (FLATSPACETIME)
  // check, 1d flat spacetime:
  Csqrt  = sqrt(cs2);
  *lam0  = vx;
  *lamm  = (vx - Csqrt)/(ONE - vx*Csqrt); // lam-
  *lamp  = (vx + Csqrt)/(ONE + vx*Csqrt); // lam+
#endif
  
  if(cure>0) printf("(PR=%d) grhd_eigenvals cure=%d\n", PR, cure);
}


// hydro left/right eigenvectors for flat spacetime
void grhd_eigLR_flat(double gxx, double gxy, double gxz, double gyy, double gyz, double gzz, 
		double alpha, double betax, double gxxup, double detg, 
		double h, double W, double rho, double cs2, double kappa, 
		double vx, double vy, double vz, double vlowx, double vlowy, double vlowz,
		double lam0, double lamm, double lamp, 
		double *L, double *R )
{


  // flat spacetime and ideal gas EoS
  // Donat et al. J Comp Phys 146 (1998)  

  double K  = h; 
  double Ap = (ONE - vx*vx )/(ONE - vx*lamp);
  double Am = (ONE - vx*vx )/(ONE - vx*lamm);
  double delta = h*h*h*W*(K-ONE)*(ONE-vx*vx)*(Ap*lamp-Am*lamm);
  double v2 = vx*vx + vy*vy + vz*vz;
  int i;  
  
  i = -1;
  
  // r01 
  R[++i] = K/(h*W);
  R[++i] = vx;
  R[++i] = vy;
  R[++i] = vz;
  R[++i] = ONE-K/(h*W);
  
  // r02
  R[++i] = W*vy;
  R[++i] = TWO*h*W*W*vx*vy;
  R[++i] = h*(ONE + TWO*W*W*vy*vy);
  R[++i] = TWO*h*W*W*vy*vz;
  R[++i] = TWO*h*W*W*vy-W*vy;
  
  // r03 
  R[++i] = W*vz;
  R[++i] = TWO*h*W*W*vx*vz;
  R[++i] = TWO*h*W*W*vy*vz;
  R[++i] = h*(ONE + TWO*W*W*vz*vz);
  R[++i] = TWO*h*W*W*vz-W*vz;
  
  // rm
  R[++i] = ONE;
  R[++i] = h*W*Am*lamm;
  R[++i] = h*W*vy;
  R[++i] = h*W*vz;
  R[++i] = h*W*Am - ONE; 

  // rp
  R[++i] = ONE;
  R[++i] = h*W*Ap*lamp;
  R[++i] = h*W*vy;
  R[++i] = h*W*vz;
  R[++i] = h*W*Ap - ONE;
  
  i = -1;

  // l01   
  L[++i] = W/(K-ONE)*(h-W);
  L[++i] = W/(K-ONE)*(W*vx);
  L[++i] = W/(K-ONE)*(W*vy);
  L[++i] = W/(K-ONE)*(W*vz);
  L[++i] = W/(K-ONE)*(-W);
    
  // l02 
  L[++i] = ONE/(h*(ONE-vx*vx))*(-vy);
  L[++i] = ONE/(h*(ONE-vx*vx))*(vx*vy);
  L[++i] = ONE/(h*(ONE-vx*vx))*(ONE-vx*vx);
  L[++i] = ONE/(h*(ONE-vx*vx))*(ZERO);
  L[++i] = ONE/(h*(ONE-vx*vx))*(-vy);
  
  // l03 
  L[++i] = ONE/(h*(ONE-vx*vx))*(-vz);
  L[++i] = ONE/(h*(ONE-vx*vx))*(vx*vz);
  L[++i] = ONE/(h*(ONE-vx*vx))*(ZERO);
  L[++i] = ONE/(h*(ONE-vx*vx))*(ONE-vx*vx);
  L[++i] = ONE/(h*(ONE-vx*vx))*(-vz);

  // lm 
  L[++i] = h*h/delta*( h*W*Ap*(vx-lamp)-vx-W*W*(v2-vx*vx)*(TWO*K-ONE)*(vx-Ap*lamp)+K*Ap*lamp );
  L[++i] = h*h/delta*( ONE+W*W*(v2-vx*vx)*(TWO*K-ONE)*(ONE-Ap)-K*Ap );
  L[++i] = h*h/delta*( W*W*vy*(TWO*K-ONE)*Ap*(vx-lamp) );
  L[++i] = h*h/delta*( W*W*vz*(TWO*K-ONE)*Ap*(vx-lamp) );
  L[++i] = h*h/delta*( -vx-W*W*(v2-vx*vx)*(TWO*K-ONE)*(vx-Ap*lamp)+K*Ap*lamp );

  // lp 
  L[++i] = -h*h/delta*( h*W*Am*(vx-lamm)-vx-W*W*(v2-vx*vx)*(TWO*K-ONE)*(vx-Am*lamm)+K*Am*lamm );
  L[++i] = -h*h/delta*( ONE+W*W*(v2-vx*vx)*(TWO*K-ONE)*(ONE-Am)-K*Am );
  L[++i] = -h*h/delta*( W*W*vy*(TWO*K-ONE)*Am*(vx-lamm) );
  L[++i] = -h*h/delta*( W*W*vz*(TWO*K-ONE)*Am*(vx-lamm) );
  L[++i] = -h*h/delta*( -vx-W*W*(v2-vx*vx)*(TWO*K-ONE)*(vx-Am*lamm)+K*Am*lamm );  
  
}


// hydro left/right eigenvectors 
void grhd_eigLR_gen(double gxx, double gxy, double gxz, double gyy, double gyz, double gzz, 
		    double alpha, double betax, double gxxup, double detg, 
		    double h, double W, double rho, double cs2, double kappa, 
		    double vx, double vy, double vz, double vlowx, double vlowy, double vlowz,
		    double lam0, double lamm, double lamp, 
		    double *L, double *R )
{
  
  // General expressions, Sec 6.2 of
  // http://relativity.livingreviews.org/Articles/lrr-2008-7/
   
  double W2  = W*W;
  double tW2 = TWO*W2;
  double h2  = h*h;
  double hW = h*W; 
  //if (hW<=ZERO) hW = TINY;
  
  double tmpLm = (alpha<=ZERO) ? BIG : (lamm+betax)/alpha;
  double tmpLp = (alpha<=ZERO) ? BIG : (lamp+betax)/alpha;

  double Vm = (vx-tmpLm)/(gxxup-vx*tmpLm);
  double Vp = (vx-tmpLp)/(gxxup-vx*tmpLp);

  double Cp = vlowx-Vp;
  double Cm = vlowx-Vm;

  double Am = (gxxup-vx*vx)/(gxxup-vx*tmpLm);
  double Ap = (gxxup-vx*vx)/(gxxup-vx*tmpLp);

  double kk = kappa/rho;  
  double K  = kk/(kk-cs2);
  //double K  = (kk==cs2)  ? BIG  : kk/(kk-cs2);
  //double K  = (kk==cs2)  ? ZERO  : kk/(kk-cs2); // assume k=kATM=0 (Cold ATM)
  //double K  = h;
  //printf(" in eig: kappa=%e rho=%e kk=%e K=%e hW=%e\n",kappa,rho,kk,K,hW);

  double tKmo  = (TWO * K - ONE);
  double omKAp = (ONE - K * Ap);
  double omKAm = (ONE - K * Am);
  double omK   = (ONE - K);

  double WoKmo = -W/omK; 
  //double WoKmo = (omK==ZERO) ? BIG : -W/omK;

  double W2oKmo = WoKmo*W;

  double cxx = gyy * gzz - gyz * gyz;
  double cxy = gxz * gyz - gxy * gzz;
  double cxz = gxy * gyz - gxz * gyy;

  double xsi    = cxx - detg * vx * vx;
  double W2xsi  = W2 * xsi;
  double hxsi   = h*xsi;
  //double oohxsi = ONE/hxsi; 
  double oohxsi = (hxsi==ZERO) ? BIG : ONE/hxsi;

  double delta    = h2 * hW * omK * (Cm - Cp) * xsi; 
  //double h2odelta = h2 / delta; 
  double h2odelta = (delta==ZERO) ? BIG : h2/delta;

  // note:    m <---------> p
  double tmpL5m = ( omK * (Vp * (W2xsi - cxx) - detg * vx) - K * W2xsi * Vp );
  double tmpL5p = ( omK * (Vm * (W2xsi - cxx) - detg * vx) - K * W2xsi * Vm );

  double tmpx = tKmo * (W2xsi * vx - cxx * vx);
  double tmpy = tKmo * (W2xsi * vy - cxy * vx);
  double tmpz = tKmo * (W2xsi * vz - cxz * vx);

  double vx2mo = (vx * vlowx - ONE);

  int i;  


  // right, Eqs 108-114
    
  i = -1;
  
  // r01 
  R[++i] = K/(hW); 
  R[++i] = vlowx;     
  R[++i] = vlowy;
  R[++i] = vlowz;
  R[++i] = ONE-K/(hW);
  
  // r02
  R[++i] = W * vlowy;
  R[++i] = h * (gxy + tW2 * vlowy * vlowx);
  R[++i] = h * (gyy + tW2 * vlowy * vlowy);
  R[++i] = h * (gyz + tW2 * vlowy * vlowz);
  R[++i] = W * vlowy * (TWO * hW - ONE);

  // r03 
  R[++i] = W * vlowz;
  R[++i] = h * (gxz + tW2 * vlowz * vlowx);
  R[++i] = h * (gyz + tW2 * vlowz * vlowy);
  R[++i] = h * (gzz + tW2 * vlowz * vlowz);
  R[++i] = W * vlowz * (TWO * hW - ONE);

  // rm
  R[++i] = ONE;
  R[++i] = hW * Cm;
  R[++i] = hW * vlowy;
  R[++i] = hW * vlowz;
  R[++i] = hW * Am - ONE; 

  // rp
  R[++i] = ONE;
  R[++i] = hW * Cp;
  R[++i] = hW * vlowy;
  R[++i] = hW * vlowz;
  R[++i] = hW * Ap - ONE;


  // left, Eqs 115-124

  i=-1;

  // l01   
  L[++i] =   WoKmo * (h - W);
  L[++i] =  W2oKmo * vx;
  L[++i] =  W2oKmo * vy;
  L[++i] =  W2oKmo * vz;
  L[++i] = -W2oKmo; 

  // l02 
  L[++i] = oohxsi * (gyz *   vlowz  - gzz * vlowy     );
  L[++i] = oohxsi * (gzz *   vlowy  - gyz * vlowz     ) * vx;
  L[++i] = oohxsi * (gzz * (-vx2mo) + gxz * vlowz * vx);
  L[++i] = oohxsi * (gyz *   vx2mo  - gxz * vlowy * vx); 
  L[++i] = oohxsi * (gyz *   vlowz  - gzz * vlowy     );
 
  // l03 
  L[++i] = oohxsi * (gyz *   vlowy  - gyy * vlowz     );
  L[++i] = oohxsi * (gyy *   vlowz  - gyz * vlowy     ) * vx;
  L[++i] = oohxsi * (gyz *   vx2mo  - gxy * vlowz * vx);
  L[++i] = oohxsi * (gyy * (-vx2mo) + gxy * vlowy * vx);
  L[++i] = oohxsi * (gyz *   vlowy  - gyy * vlowz     );

  // lm 
  L[++i] = h2odelta * ( hW  * Vp * xsi + tmpL5m );
  L[++i] = h2odelta * ( cxx * omKAp + Vp * tmpx );
  L[++i] = h2odelta * ( cxy * omKAp + Vp * tmpy );
  L[++i] = h2odelta * ( cxz * omKAp + Vp * tmpz );
  L[++i] = h2odelta * ( tmpL5m );

  // lp 
  L[++i] = - h2odelta * ( hW  * Vm * xsi + tmpL5p );
  L[++i] = - h2odelta * ( cxx * omKAm + Vm * tmpx ); 
  L[++i] = - h2odelta * ( cxy * omKAm + Vm * tmpy );
  L[++i] = - h2odelta * ( cxz * omKAm + Vm * tmpz );
  L[++i] = - h2odelta * ( tmpL5p );
  
}

/* WT: same as grhd_eigLR_gen, but without NANs */
/* hydro left/right eigenvectors */
void grhd_eigenvecLR_gen(double gxx, double gxy, double gxz,
                         double gyy, double gyz, double gzz,
                         double alpha, double betax, double gxxup, double detg,
                         double h, double W,
                         double rho, double cs2, double kappa,
                         double vx, double vy, double vz,
                         double vlowx, double vlowy, double vlowz,
                         double lam0, double lamm, double lamp,
                         double *L, double *R)
{
  // General expressions, Sec 6.2 of
  // http://relativity.livingreviews.org/Articles/lrr-2008-7/

  double W2, tW2, h2, hW, oohW;
  double tmpLm, tmpLp;
  double gxxup_vxtmpLm, gxxup_vxtmpLp;
  double Vm, Vp;
  double Cp, Cm;
  double Am, Ap;
  double kk, K;  
  double tKmo, omKAp, omKAm, omK;
  double WoKmo; 
  double W2oKmo;
  double cxx, cxy, cxz;
  double xsi, W2xsi, hxsi;
  double oohxsi;
  double delta, h2odelta;
  double tmpL5m, tmpL5p;
  double tmpx, tmpy, tmpz;
  double vx2mo;
  int i;

  W2  = W*W;
  tW2 = TWO*W2;
  h2  = h*h;
  hW = h*W; 
  //if (hW<=ZERO) hW = TINY;
  /* oohW = 1.0/hW; */
  if(hW!=0.0)  oohW = 1.0/hW;
  else         oohW = BIG;

  tmpLm = (alpha<=ZERO) ? BIG : (lamm+betax)/alpha;
  tmpLp = (alpha<=ZERO) ? BIG : (lamp+betax)/alpha;

  gxxup_vxtmpLm = gxxup-vx*tmpLm;
  gxxup_vxtmpLp = gxxup-vx*tmpLp;
  if(gxxup_vxtmpLm==0.0)  gxxup_vxtmpLm = TINY;
  if(gxxup_vxtmpLp==0.0)  gxxup_vxtmpLp = TINY;

  Vm = (vx-tmpLm)/(gxxup_vxtmpLm);
  Vp = (vx-tmpLp)/(gxxup_vxtmpLp);

  Cp = vlowx-Vp;
  Cm = vlowx-Vm;

  Am = (gxxup-vx*vx)/(gxxup_vxtmpLm);
  Ap = (gxxup-vx*vx)/(gxxup_vxtmpLp);

  /* kk = kappa/rho;  K  = kk/(kk-cs2); */
  if(rho!=0.0)  kk = kappa/rho;
  else          kk = BIG;
  if(kk!=cs2)   K  = kk/(kk-cs2);
  else          K  = BIG;
  //K  = (kk==cs2)  ? BIG  : kk/(kk-cs2);
  //K  = (kk==cs2)  ? ZERO : kk/(kk-cs2); // assume k=kATM=0 (Cold ATM)
  //K  = h;
  //printf(" in eig: kappa=%e rho=%e kk=%e K=%e hW=%e\n",kappa,rho,kk,K,hW);

  tKmo  = (TWO * K - ONE);
  omKAp = (ONE - K * Ap);
  omKAm = (ONE - K * Am);
  omK   = (ONE - K);

  /* WoKmo = -W/omK; */
  if(omK!=0.0)  WoKmo = -W/omK;
  else          WoKmo = BIG;

  W2oKmo = WoKmo*W;

  cxx = gyy * gzz - gyz * gyz;
  cxy = gxz * gyz - gxy * gzz;
  cxz = gxy * gyz - gxz * gyy;

  xsi    = cxx - detg * vx * vx;
  W2xsi  = W2 * xsi;
  hxsi   = h*xsi;
  //oohxsi = ONE/hxsi; 
  oohxsi = (hxsi==ZERO) ? BIG : ONE/hxsi;

  delta    = h2 * hW * omK * (Cm - Cp) * xsi; 
  //h2odelta = h2 / delta; 
  h2odelta = (delta==ZERO) ? BIG : h2/delta;

  // note:    m <---------> p
  tmpL5m = ( omK * (Vp * (W2xsi - cxx) - detg * vx) - K * W2xsi * Vp );
  tmpL5p = ( omK * (Vm * (W2xsi - cxx) - detg * vx) - K * W2xsi * Vm );

  tmpx = tKmo * (W2xsi * vx - cxx * vx);
  tmpy = tKmo * (W2xsi * vy - cxy * vx);
  tmpz = tKmo * (W2xsi * vz - cxz * vx);

  vx2mo = (vx * vlowx - ONE);

  // right, Eqs 108-114
  i = -1;

  // r01 
  R[++i] = K*oohW;
  R[++i] = vlowx;
  R[++i] = vlowy;
  R[++i] = vlowz;
  R[++i] = ONE-K*oohW;

  // r02
  R[++i] = W * vlowy;
  R[++i] = h * (gxy + tW2 * vlowy * vlowx);
  R[++i] = h * (gyy + tW2 * vlowy * vlowy);
  R[++i] = h * (gyz + tW2 * vlowy * vlowz);
  R[++i] = W * vlowy * (TWO * hW - ONE);

  // r03 
  R[++i] = W * vlowz;
  R[++i] = h * (gxz + tW2 * vlowz * vlowx);
  R[++i] = h * (gyz + tW2 * vlowz * vlowy);
  R[++i] = h * (gzz + tW2 * vlowz * vlowz);
  R[++i] = W * vlowz * (TWO * hW - ONE);

  // rm
  R[++i] = ONE;
  R[++i] = hW * Cm;
  R[++i] = hW * vlowy;
  R[++i] = hW * vlowz;
  R[++i] = hW * Am - ONE; 

  // rp
  R[++i] = ONE;
  R[++i] = hW * Cp;
  R[++i] = hW * vlowy;
  R[++i] = hW * vlowz;
  R[++i] = hW * Ap - ONE;


  // left, Eqs 115-124
  i=-1;

  // l01   
  L[++i] =   WoKmo * (h - W);
  L[++i] =  W2oKmo * vx;
  L[++i] =  W2oKmo * vy;
  L[++i] =  W2oKmo * vz;
  L[++i] = -W2oKmo; 

  // l02 
  L[++i] = oohxsi * (gyz *   vlowz  - gzz * vlowy     );
  L[++i] = oohxsi * (gzz *   vlowy  - gyz * vlowz     ) * vx;
  L[++i] = oohxsi * (gzz * (-vx2mo) + gxz * vlowz * vx);
  L[++i] = oohxsi * (gyz *   vx2mo  - gxz * vlowy * vx); 
  L[++i] = oohxsi * (gyz *   vlowz  - gzz * vlowy     );
 
  // l03 
  L[++i] = oohxsi * (gyz *   vlowy  - gyy * vlowz     );
  L[++i] = oohxsi * (gyy *   vlowz  - gyz * vlowy     ) * vx;
  L[++i] = oohxsi * (gyz *   vx2mo  - gxy * vlowz * vx);
  L[++i] = oohxsi * (gyy * (-vx2mo) + gxy * vlowy * vx);
  L[++i] = oohxsi * (gyz *   vlowy  - gyy * vlowz     );

  // lm 
  L[++i] = h2odelta * ( hW  * Vp * xsi + tmpL5m );
  L[++i] = h2odelta * ( cxx * omKAp + Vp * tmpx );
  L[++i] = h2odelta * ( cxy * omKAp + Vp * tmpy );
  L[++i] = h2odelta * ( cxz * omKAp + Vp * tmpz );
  L[++i] = h2odelta * ( tmpL5m );

  // lp 
  L[++i] = - h2odelta * ( hW  * Vm * xsi + tmpL5p );
  L[++i] = - h2odelta * ( cxx * omKAm + Vm * tmpx ); 
  L[++i] = - h2odelta * ( cxy * omKAm + Vm * tmpy );
  L[++i] = - h2odelta * ( cxz * omKAm + Vm * tmpz );
  L[++i] = - h2odelta * ( tmpL5p );
}


// hydro left/right eigenvectors (normalized)
void grhd_eigLR_norm(double gxx, double gxy, double gxz, double gyy, double gyz, double gzz, 
		     double alpha, double betax, double gxxup, double detg, 
		     double h, double W, double rho, double cs2, double kappa, 
		     double vx, double vy, double vz, double vlowx, double vlowy, double vlowz,
		     double lam0, double lamm, double lamp, 
		     double *L, double *R )
{

  // General expressions, Sec 6.2 of
  // http://relativity.livingreviews.org/Articles/lrr-2008-7/
  // normalized after
  // http://arxiv.org/abs/1206.5972

  double W2  = W*W;
  double tW2 = TWO*W2;
  double h2  = h*h;
  double hW = h*W; 
  //if (hW<=ZERO) hW = TINY;
  
  double tmpLm = (alpha<=ZERO) ? BIG : (lamm+betax)/alpha;
  double tmpLp = (alpha<=ZERO) ? BIG : (lamp+betax)/alpha;

  double Vm = (vx-tmpLm)/(gxxup-vx*tmpLm);
  double Vp = (vx-tmpLp)/(gxxup-vx*tmpLp);

  double Cp = vlowx-Vp;
  double Cm = vlowx-Vm;

  double Am = (gxxup-vx*vx)/(gxxup-vx*tmpLm);
  double Ap = (gxxup-vx*vx)/(gxxup-vx*tmpLp);

  //double kk = (rho<=RHOATM) ? kappa/RHOATM : kappa/rho;
  //
  double kk = kappa/rho;
  
  //
  double K  = kk/(kk-cs2);
  //double K  = (kk==cs2)  ? BIG  : kk/(kk-cs2);
  //double K  = (kk==cs2)  ? ZERO  : kk/(kk-cs2); // assume k=kATM=0 (Cold ATM)
  //double K  = h;

  //printf(" in eig: kappa=%e rho=%e kk=%e K=%e hW=%e\n",kappa,rho,kk,K,hW);

  double tKmo  = (TWO * K - ONE);
  double omKAp = (ONE - K * Ap);
  double omKAm = (ONE - K * Am);
  double omK   = (ONE - K);

  double WoKmo = -W/omK; 
  //double WoKmo = (omK==ZERO) ? BIG : -W/omK;

  double W2oKmo = WoKmo*W;

  double cxx = gyy * gzz - gyz * gyz;
  double cxy = gxz * gyz - gxy * gzz;
  double cxz = gxy * gyz - gxz * gyy;

  double xsi    = cxx - detg * vx * vx;
  double W2xsi  = W2 * xsi;
  double hxsi   = h*xsi;
  //double oohxsi = ONE/hxsi; 
  double oohxsi = (hxsi==ZERO) ? BIG : ONE/hxsi;

  double delta    = h2 * hW * omK * (Cm - Cp) * xsi; 
  //double h2odelta = h2 / delta; 
  double h2odelta = (delta==ZERO) ? BIG : h2/delta;

  // note:    m <---------> p
  double tmpL5m = ( omK * (Vp * (W2xsi - cxx) - detg * vx) - K * W2xsi * Vp );
  double tmpL5p = ( omK * (Vm * (W2xsi - cxx) - detg * vx) - K * W2xsi * Vm );

  double tmpx = tKmo * (W2xsi * vx - cxx * vx);
  double tmpy = tKmo * (W2xsi * vy - cxy * vx);
  double tmpz = tKmo * (W2xsi * vz - cxz * vx);

  double vx2mo = (vx * vlowx - ONE);

  double norm  = kk-cs2;
  double Wocs2 = W/cs2;

  int i;  


  // right, Eqs 108-114 r01 normalized
    
  i = -1;
  
  // r01 
  R[++i] = kk/(hW); 
  R[++i] = norm*vlowx;     
  R[++i] = norm*vlowy;
  R[++i] = norm*vlowz;
  R[++i] = norm-kk/(hW);
  
  // r02
  R[++i] = W * vlowy;
  R[++i] = h * (gxy + tW2 * vlowy * vlowx);
  R[++i] = h * (gyy + tW2 * vlowy * vlowy);
  R[++i] = h * (gyz + tW2 * vlowy * vlowz);
  R[++i] = W * vlowy * (TWO * hW - ONE);

  // r03 
  R[++i] = W * vlowz;
  R[++i] = h * (gxz + tW2 * vlowz * vlowx);
  R[++i] = h * (gyz + tW2 * vlowz * vlowy);
  R[++i] = h * (gzz + tW2 * vlowz * vlowz);
  R[++i] = W * vlowz * (TWO * hW - ONE);

  // rm
  R[++i] = ONE;
  R[++i] = hW * Cm;
  R[++i] = hW * vlowy;
  R[++i] = hW * vlowz;
  R[++i] = hW * Am - ONE; 

  // rp
  R[++i] = ONE;
  R[++i] = hW * Cp;
  R[++i] = hW * vlowy;
  R[++i] = hW * vlowz;
  R[++i] = hW * Ap - ONE;


  // left, Eqs 115-124 l01 normalized

  i=-1;

  // l01   
  /*
  L[++i] =   WoKmo * (h - W)/norm; 
  L[++i] =  W2oKmo * vx/norm; 
  L[++i] =  W2oKmo * vy/norm; 
  L[++i] =  W2oKmo * vz/norm; 
  L[++i] = -W2oKmo/norm; 
  */
  L[++i] =   Wocs2 * (h - W); 
  L[++i] =  W*Wocs2 * vx; 
  L[++i] =  W*Wocs2 * vy; 
  L[++i] =  W*Wocs2 * vz; 
  L[++i] = -W*Wocs2; 


  // l02 
  L[++i] = oohxsi * (gyz *   vlowz  - gzz * vlowy     );
  L[++i] = oohxsi * (gzz *   vlowy  - gyz * vlowz     ) * vx;
  L[++i] = oohxsi * (gzz * (-vx2mo) + gxz * vlowz * vx);
  L[++i] = oohxsi * (gyz *   vx2mo  - gxz * vlowy * vx); 
  L[++i] = oohxsi * (gyz *   vlowz  - gzz * vlowy     );
 
  // l03 
  L[++i] = oohxsi * (gyz *   vlowy  - gyy * vlowz     );
  L[++i] = oohxsi * (gyy *   vlowz  - gyz * vlowy     ) * vx;
  L[++i] = oohxsi * (gyz *   vx2mo  - gxy * vlowz * vx);
  L[++i] = oohxsi * (gyy * (-vx2mo) + gxy * vlowy * vx);
  L[++i] = oohxsi * (gyz *   vlowy  - gyy * vlowz     );

  // lm 
  L[++i] = h2odelta * ( hW  * Vp * xsi + tmpL5m );
  L[++i] = h2odelta * ( cxx * omKAp + Vp * tmpx );
  L[++i] = h2odelta * ( cxy * omKAp + Vp * tmpy );
  L[++i] = h2odelta * ( cxz * omKAp + Vp * tmpz );
  L[++i] = h2odelta * ( tmpL5m );

  // lp 
  L[++i] = - h2odelta * ( hW  * Vm * xsi + tmpL5p );
  L[++i] = - h2odelta * ( cxx * omKAm + Vm * tmpx ); 
  L[++i] = - h2odelta * ( cxy * omKAm + Vm * tmpy );
  L[++i] = - h2odelta * ( cxz * omKAm + Vm * tmpz );
  L[++i] = - h2odelta * ( tmpL5p );

}


// wrapper
void grhd_eigLR(double gxx, double gxy, double gxz, double gyy, double gyz, double gzz, 
		double alpha, double betax, double gxxup, double detg, 
		double h, double W, double rho, double cs2, double kappa, 
		double vx, double vy, double vz, double vlowx, double vlowy, double vlowz,
		double lam0, double lamm, double lamp, 
		double *L, double *R )

{

#if (1)
  grhd_eigLR_gen( gxx,  gxy,  gxz,  gyy,  gyz,  gzz, 
		  alpha,  betax,  gxxup,  detg, 
		  h,  W,  rho,  cs2,  kappa, 
		  vx,  vy,  vz,  vlowx,  vlowy,  vlowz,
		  lam0,  lamm,  lamp, 
		  L,  R );
#endif

#if (0)
  grhd_eigLR_norm( gxx,  gxy,  gxz,  gyy,  gyz,  gzz, 
		  alpha,  betax,  gxxup,  detg, 
		  h,  W,  rho,  cs2,  kappa, 
		  vx,  vy,  vz,  vlowx,  vlowy,  vlowz,
		  lam0,  lamm,  lamp, 
		  L,  R );
#endif

#if (FLATSPACETIME)
  grhd_eigLR_flat( gxx,  gxy,  gxz,  gyy,  gyz,  gzz, 
		  alpha,  betax,  gxxup,  detg, 
		  h,  W,  rho,  cs2,  kappa, 
		  vx,  vy,  vz,  vlowx,  vlowy,  vlowz,
		  lam0,  lamm,  lamp, 
		  L,  R );
#endif

}


#if (DEBUG)
main()
{

  double gxx = 1.;
  double gxy = 0.; 
  double gxz = 0.;
  double gyy = 1.;
  double gyz = 0.;
  double gzz = 1.;
  double alpha = 1.;
  double betax = 0.;
  double gupxx = 1.;
  double detg  = 1.; 
  
  double eos_gammamo = 0.666666666666;
  double rho  = 1e-3;
  double epsl = 100./(eos_gammamo*pow(rho,eos_gammamo));
  double pres =  eos_gammamo*rho*epsl;
  double vx   = 0.99;
  double vy   = 0.01;
  double vz   = 0.01;

  double vlowx, vlowy, vlowz;

  vlowx = vx;
  vlowy = vy;
  vlowz = vz;

  double v2 = vx*vlowx + vy*vlowy + vz*vlowz;
  double h  = 1.+epsl+pres/rho;
  double Wlor  = 1./sqrt(1.-v2);
  double kappa = eos_gammamo*rho;
  double chi   = eos_gammamo*epsl;
  double cs2   = (chi + pres/(rho*rho)*kappa)/h;
    

  double lamo, lamm, lamp;
  double L[25],R[25];

  grhd_eig( alpha, betax, vx, v2, gupxx, cs2,
	    &lamo, &lamm, &lamp );

  grhd_eigLR( gxx, gxy, gxz, gyy, gyz, gzz, 
	      alpha, betax,  gupxx, detg, 
	      h, Wlor, rho, cs2, kappa, 
	      vx, vy, vz, vlowx, vlowy,  vlowz,
	      lamo,  lamm,  lamp, 
	      L, R );
  
  int k,j;
  int ne = 5;
  
  printf("L\n");
  for( k = 0 ; k < ne ; k++ ) 
    {
      printf(" k=%d\n",k);
      for( j = 0 ; j < ne ; j++ ) printf(" %e", L[k*ne + j]);
      printf("\n"); 
    }
  printf("R\n");
  for( k = 0 ; k < ne ; k++ ) 
    {
      printf(" k=%d\n",k);
      for( j = 0 ; j < ne ; j++ ) printf(" %e", R[k*ne + j]);
      printf("\n"); 
    }
  printf("-----\n");

}
#endif
