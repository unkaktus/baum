/* 
   These routines implements analytic fits of cold EoS Sly, FPS and APR using the formulas derived in [2]. 
   The fit is performed in a thermodynamically consistent way as described in [1,2].
   
   G=c=Msun=1 units are used

   Refs:
   ================================================ 
   Shibata et al., 
   [1] Phys.Rev. D71 (2005) 084021
   arXiv:astro-ph/0503119
   [2] Phys.Rev. D73 (2006) 064027 
   arXiv:astro-ph/0603145 
   ================================================ 

   Author: S.Bernuzzi 
*/



#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>




/* ================================================  
   Analytic fits : epsl( rho ) , Eq.(25) of [2]
   ================================================ */

void Shib_EpslFit( const int which_eos,
		 double* Ecold_p,
		 double* DEcoldDrho_p,                  
		 double* D2EcoldD2rho_p, 
		 const double rho ) 
{
  
  /* see Tab.I of [1,2] for the coefs */
  static double  eos_SLy[16] = {0.0, 0.1037, 0.1956, 39264., 1.9503, 254.83, 1.3823, -1.234, 1.2e5, 9e5, 4., 0.75, 0.057, 0.138, 0.84, 0.338};
  static double  eos_FPS[16] = {0.0, 0.15806, 0.220, 5956.4, 1.633, 170.68, 1.1056, -0.703, 2e4, 9e5, 5., 0.75, 0.0627, 0.1387, 0.56, 0.308};
  static double  eos_APR[16] = {0.0, 0.0889, 0.1821, 5.945e5, 2.4265, 1600.0, 2.165, -6.96, 1.2e5, 9e5, 4., 0.75, 0.057, 0.138, 0.84, 0.338};
  
  static double* eos_index[3] = {eos_SLy,eos_FPS,eos_APR};
  
  double* eosP;
  
  double Ecold,DEcoldDrho,D2EcoldD2rho;
  double arg1, arg2;
  double tmp_exp;
  double Fermi_1,Fermi_m1, Fermi_2,Fermi_m2;
  double dFermi_1,dFermi_m1, dFermi_2,dFermi_m2;
  double d2Fermi_1,d2Fermi_m1, d2Fermi_2,d2Fermi_m2;
  double oorho, oorho2;
  double rhop2,rhop4,rhop6,rhop13,rhop15;
  double tmp1,tmp2,tmp3,tmp4a,tmp4b,tmp4,arg56;
  double Qr1, Qr2, Qr3;
  double dQr1, dQr2, dQr3;
  double d2Qr1, d2Qr2, d2Qr3;
  
  eosP = eos_index[which_eos];
  
  oorho  = 1./rho;
  oorho2 = oorho*oorho;

  arg1 = eosP[8]*rho - eosP[10];
  arg2 = eosP[9]*rho - eosP[11];

  /* ---------------------------------------
     Compute the Fermi functions separating 
     the three regions of the FIT
     ------------------------------------ */ 
    if (arg1 > 40.0 ) {
      Fermi_1  = 0.0;
      Fermi_m1 = 1.0;
    } else if (arg1 < -40.0) {
      Fermi_1  = 1.0;
      Fermi_m1 = 0.0;
    } else {
      tmp_exp = exp(arg1);
      Fermi_1  = 1.0 / (1.0 + tmp_exp);
      Fermi_m1 = tmp_exp * Fermi_1;
    }
    
    dFermi_m1 = eosP[8] * Fermi_1 * Fermi_m1;
    dFermi_1  = - dFermi_m1;

    d2Fermi_m1 =  eosP[8] * ( dFermi_1 * Fermi_m1 + Fermi_1 * dFermi_m1 );
    d2Fermi_1  = - d2Fermi_m1;

    if (arg2 > 40.0 ) {
      Fermi_2  = 0.0;
      Fermi_m2 = 1.0;
    } else if (arg2 < -40.0) {
      Fermi_2  = 1.0;
      Fermi_m2 = 0.0;
    } else {
      tmp_exp = exp(arg2);
      Fermi_2  = 1.0 / (1.0 + tmp_exp);
      Fermi_m2 = tmp_exp * Fermi_2;
    }

    dFermi_m2 = eosP[9] * Fermi_2 * Fermi_m2;
    dFermi_2  = - dFermi_m2;

    d2Fermi_m2 = eosP[9] * ( dFermi_2 * Fermi_m2 + Fermi_2 * dFermi_m2 );
    d2Fermi_2  = - d2Fermi_m2;

    rhop2  = pow(rho,eosP[2]);
    rhop4  = pow(rho,eosP[4]);
    rhop6  = pow(rho,eosP[6]);
    rhop13 = pow(rho,eosP[13]);
    rhop15 = pow(rho,eosP[15]);

    arg56 = 1.+eosP[5] * rhop6;
    tmp1  = ( 1. + eosP[1] * rhop2 + eosP[3] * rhop4 );
    tmp2  = pow( arg56, eosP[7]);

    tmp3  = eosP[1] * eosP[2] * rhop2 * oorho + eosP[3] * eosP[4] * rhop4 * oorho;
    tmp4a = eosP[7] * tmp2 / arg56;
    tmp4b = eosP[5] * eosP[6] * rhop6 * oorho;
    tmp4  = tmp4a * tmp4b;

    Qr1   = tmp1 * tmp2 - 1.;
    dQr1  = tmp3 * tmp2 + tmp1 * tmp4; 
    d2Qr1 = ( eosP[1] * eosP[2] * (eosP[2]-1.) * rhop2 * oorho2 + eosP[3] * eosP[4] * (eosP[4]-1.) * rhop4 * oorho2 ) * tmp2 
      + 2. * tmp3 * tmp4 
      + tmp1*( ( eosP[7] * (eosP[7]-1.) * tmp2 /(arg56*arg56) * tmp4b ) * tmp4b +
	       ( eosP[5] * eosP[6] * (eosP[6]-1.) * rhop6 * oorho2    ) * tmp4a ); 
    
    Qr2   = eosP[12] * rhop13;
    dQr2  = eosP[12] * eosP[13] * rhop13*oorho;
    d2Qr2 = eosP[12] * eosP[13] * (eosP[13]-1.) * rhop13*oorho2;
	    
    Qr3   = eosP[14] * rhop15;
    dQr3  = eosP[14] * eosP[15] * rhop15*oorho;
    d2Qr3 = eosP[14] * eosP[15] * (eosP[15]-1.) * rhop15*oorho2;

    Ecold        = Qr1 * Fermi_m1 + Qr2 * Fermi_1 * Fermi_m2 +  Qr3 * Fermi_2;

    DEcoldDrho   = dQr1 * Fermi_m1 + Qr1 * dFermi_m1 
      + dQr2 * Fermi_1 * Fermi_m2  + Qr2 * dFermi_1 * Fermi_m2 + Qr2 * Fermi_1 * dFermi_m2 
      + dQr3 * Fermi_2 + Qr3 * dFermi_2;
    
    D2EcoldD2rho = d2Qr1 * Fermi_m1 + 2.0 * dQr1 * dFermi_m1 + Qr1 * d2Fermi_m1 
      + d2Qr2 * Fermi_1  + 2.0 * dQr2 * dFermi_1  + Qr2 * d2Fermi_1; 
    
    *Ecold_p        = Ecold;
    *DEcoldDrho_p   = DEcoldDrho;
    *D2EcoldD2rho_p = D2EcoldD2rho;
    
}




/* ================================================  
   Analytic fits routine 
   Cold EOS, Thermodynamical consistent procedure
   ================================================ */

void Shib_EOSFit_Cold( const int which_eos,
		       const int n_elems,
		       const double* rho_p,
		       double* pressure_p,
		       double* DpDrho_p,
		       double* DpDepsl_p,
		       double* cs2_p,
		       double* ene_p,
		       double* epsl_p,
		       double* DepslDrho_p,
		       double* D2epslD2rho_p )
{
  
  /* tmp qts at point */
  double rho;
  double pressure,Pcold,DpDrho,DpDepsl;
  double Ecold, DEcoldDrho, D2EcoldD2rho;
  double cs2;

  int i;

  for (i=0;i<n_elems;i++) {
    
    rho = rho_p[i];
    
    /* ****************************************
       Compute the Cold EOS
       ************************************* */
    
    Shib_EpslFit(which_eos,&Ecold,&DEcoldDrho,&D2EcoldD2rho,rho);
    
    /* -------------------------------------------
       Termodinamically consistent relations @ T = 0                      
       1-par EOS: P(rho,epsl) = P(rho)       
       Pcold = rho*rho*DEcold_Drho         
       ------------------------------------------- */
    
    Pcold    = rho*rho*DEcoldDrho;
    pressure = Pcold;

    DpDrho   = 2.0/rho*Pcold + rho*rho*D2EcoldD2rho; 
    DpDepsl  = 0.0;
    
    cs2 = ( DpDrho + pressure*DpDepsl/(rho*rho) ) /
      ( 1.0 + Ecold + pressure/rho );
    
    /* ****************************************
       Output the computed EOS values
       ************************************* */
    
    pressure_p[i] =  pressure;	
    DpDrho_p[i]   =  DpDrho;	
    DpDepsl_p[i]  =  DpDepsl;	
    cs2_p[i]      =  cs2;	
    ene_p[i]      =  rho*(Ecold+1.);
    epsl_p[i]        = Ecold;
    DepslDrho_p[i]   = DEcoldDrho;
    D2epslD2rho_p[i] = D2EcoldD2rho;

  } /* end for */
  
}


