/* 
   These routines implements analytic fits of cold EoS Sly and FPS using the formulas derived in [1]. 
   The fit can be perfomed either in a non-thermodynamically consistent way [1] or in a 
   thermodynamically consistent one as described in [2]. 
   The former is however more accurate. 

   G=c=Msun=1 units are used in i/o 

   Refs:
   ================================================ 
   [1] Haensel, P.; Potekhin, A. Y. 
   Astronomy and Astrophysics, v.428, p.191-197 (2004) 
   arXiv:astro-ph/0408324 
   ================================================ 
   [2] Shibata et al., 
   Phys Rev D 71, 084021 (2005) 
   arXiv:astro-ph/0503119 
   ================================================ 

   Author: S.Bernuzzi 
   (Taken, adpted and extended from R.De Pietri implementation in Whisky)
*/


#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>



/* ================================================  
   Analytic fits : epsl( rho ) , Eq.(15) of [1]
   ================================================ */

void HP_EpslFit( const int which_eos,
		 double* Ecold_p,
		 double* DEcoldDrho_p,                  
		 double* D2EcoldD2rho_p, 
		 const double rho ) 
{
  
  static double  eos_SLy[8] = {0.0, 0.423, 2.42, 0.031, 0.78, 0.238, 0.912, 3.674};
  static double  eos_FPS[8] = {0.0, 0.320, 2.17, 0.173, 3.01, 0.540, 0.847, 3.681};
  
  static double* eos_index[2] = {eos_SLy,eos_FPS};
  
  double* eosP;
  
  double Ecold,DEcoldDrho,D2EcoldD2rho;
  double arg1;
  double tmp_exp;
  double Fermi_1,Fermi_m1;
  double dFermi_1,dFermi_m1;
  double d2Fermi_1,d2Fermi_m1;
  double tmp1,tmp2,tmp3,tmp4,tmp5;
  double Qr1,Qr2;
  double dQr1,dQr2;
  double d2Qr1,d2Qr2;
  
  double n,log10_n;

  /* ---------------------------------------
   Before eq (4) pag. 192:  n = \rho /m_0 and  m_0 = 1.66 e-24 g   
   n in fm^-3 -> 10^-39 * [(cactusM/cactusL^3))/1.66e-24 (in CGS)] 
   --------------------------------------- */
  
  const double Kd = 3.720718892340736e+02; /* 372.04;  */
  const double Vlog10 = 0.4342944819032518276511289; /* =  1/Log(10) *//* 0.43429448190325; */
    
  eosP = eos_index[which_eos];
  
  n       = Kd * rho;
  log10_n = Vlog10*log(n);
  
  arg1 = eosP[6]*(log10_n+eosP[7]);
  
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
    
    dFermi_m1 = ((eosP[6]*Vlog10*Kd)/n) * Fermi_1 * Fermi_m1;
    dFermi_1  = - dFermi_m1;
    d2Fermi_1 = ((eosP[6]*Vlog10*Kd)/n) * dFermi_1 * (Fermi_1 - Fermi_m1)
      - Kd*dFermi_1/n;
    d2Fermi_m1  = - d2Fermi_1;
    
    tmp1= eosP[1]*pow(n,eosP[2]);
    tmp2= eosP[3]*pow(n,eosP[4]);
    tmp3= (1.0+eosP[5]*n);
    Qr1   = (tmp1+tmp2)/(tmp3*tmp3);
    dQr1  = ( Kd*eosP[2]*(tmp1/n)+Kd*eosP[4]*(tmp2/n))/(tmp3*tmp3)
      - 2.0 * (tmp1+tmp2)/(tmp3*tmp3*tmp3) * Kd * eosP[5];
    d2Qr1 = Kd*Kd*
      (( eosP[2]*(eosP[2]-1.0)*(tmp1/(n*n)) 
	 +eosP[4]*(eosP[4]-1.0)*(tmp2/(n*n)))/(tmp3*tmp3)
       - 4.0 * (eosP[2]*(tmp1/n)+eosP[4]*(tmp2/n))/(tmp3*tmp3*tmp3) * eosP[5]
       + 6.0 * (tmp1+tmp2)/(tmp3*tmp3*tmp3*tmp3) * eosP[5]* eosP[5]
       );
    
    tmp4  = 2.1 *pow(n,0.585);
    tmp5  = (8e-6 + tmp4);
    Qr2   = n / tmp5;
    dQr2  = Kd/ tmp5 - Kd * 0.585 * tmp4/(tmp5*tmp5);
    d2Qr2 = Kd*Kd*tmp4*(
			- 0.585 * tmp5
			- 0.585 * 0.585 *tmp5
			+ 2.0*0.585* 0.585 * tmp4
			)/(n*tmp5*tmp5*tmp5); 
    
    Ecold        = Qr1 * Fermi_m1 + Qr2 * Fermi_1;
    DEcoldDrho   = dQr1 * Fermi_m1 + Qr1 * dFermi_m1
      + dQr2 * Fermi_1  + Qr2 * dFermi_1;
    D2EcoldD2rho = d2Qr1 * Fermi_m1 + 2.0 * dQr1 * dFermi_m1 + Qr1 * d2Fermi_m1
      + d2Qr2 * Fermi_1  + 2.0 * dQr2 * dFermi_1  + Qr2 * d2Fermi_1;
    
    *Ecold_p        = Ecold;
    *DEcoldDrho_p   = DEcoldDrho;
    *D2EcoldD2rho_p = D2EcoldD2rho;
    
}





/* ================================================  
   Analytic fits : press( ene ) , Eq.(14) of [1]
   ================================================ */

void HP_PresFit( const int which_eos,
		 double* pres_p,
		 double* DPDene_p,                  
		 const double ene ) 

{
    double press,dPdE,E;
    double arg1,arg2,arg3,arg4;
    double Q1,Q2,Q3,Q4;
    double dQ1,dQ2,dQ3,dQ4;
    double tmp_exp;
    double Fermi_1,Fermi_m1;
    double Fermi_2,Fermi_m2;
    double Fermi_3,Fermi_m3;
    double Fermi_4,Fermi_m4;
    double dFermi_1,dFermi_2,dFermi_3,dFermi_4;
    
    static double  eos_SLy[19] = {0.0, 6.22, 6.121, 0.005925, 0.16326, 6.48, 11.4971, 19.105, 0.8938, 6.54,
				  11.4950, -22.775, 1.5707, 4.3, 14.08, 27.80, -1.653, 1.50, 14.67};
    static double  eos_FPS[19] = {0.0, 6.22, 6.121, 0.006004, 0.16345, 6.50, 11.8440, 17.24, 1.065, 6.54, 
				  11.8421, -22.003, 1.5552, 9.3, 14.19, 23.73, -1.508, 1.79, 15.13};
    static double* eos_index[2] = {eos_SLy,eos_FPS};
    double* eosP;

    double xi,zeta,DzetaDxi;
    const double Kxi   = 6.175859799237819e+17;
    const double Kzeta = 1.801598783618508e-39;
    const double Vlog10 = 0.4342944819032518276511289; /* =  1/Log(10) */
    
    eosP = eos_index[which_eos];
    E = ene;
    xi = Vlog10 * log(E*Kxi);
    
    arg1 = eosP[5]*(xi-eosP[6]);
    arg2 = eosP[9]*(eosP[10]-xi);
    arg3 = eosP[13]*(eosP[14]-xi);
    arg4 = eosP[17]*(eosP[18]-xi);

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
    /* -----------------------------*/
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
    /* -----------------------------*/
    if (arg3 > 40.0 ) {
        Fermi_3  = 0.0;
        Fermi_m3 = 1.0;
    } else if (arg3 < -40.0) {
        Fermi_3  = 1.0;
        Fermi_m3 = 0.0;
    } else {
        tmp_exp = exp(arg3);
        Fermi_3  = 1.0 / (1.0 + tmp_exp);
        Fermi_m3 = tmp_exp * Fermi_3;
    }
    /* -----------------------------*/
    if (arg4 > 40.0 ) {
        Fermi_4  = 0.0;
        Fermi_m4 = 1.0;
    } else if (arg4 < -40.0) {
        Fermi_4  = 1.0;
        Fermi_m4 = 0.0;
    } else {
        tmp_exp = exp(arg4);
        Fermi_4  = 1.0 / (1.0 + tmp_exp);
        Fermi_m4 = tmp_exp * Fermi_4;
    }
    /* -----------------------------*/
    dFermi_1 = -eosP[5] * Fermi_1 * Fermi_m1;
    dFermi_2 = eosP[9]  * Fermi_2 * Fermi_m2;
    dFermi_3 = eosP[13] * Fermi_3 * Fermi_m3;
    dFermi_4 = eosP[17] * Fermi_4 * Fermi_m4;

    Q1 = (eosP[1]+eosP[2]*xi+eosP[3]*xi*xi*xi)/
         (1.0+eosP[4]*xi);
    Q2 = eosP[7]+ eosP[8]*xi;
    Q3 = eosP[11]+ eosP[12]*xi;
    Q4 = eosP[15]+ eosP[16]*xi;

    dQ1 = (  (eosP[2]+3.0*eosP[3]*xi*xi) 
           - Q1 * eosP[4]
          ) / (1.0+eosP[4]*xi);
    dQ2 = eosP[8];
    dQ3 = eosP[12];
    dQ4 = eosP[16];
    
    zeta     = Q1 * Fermi_1 + Q2* Fermi_2 + Q3* Fermi_3 + Q4* Fermi_4;
    DzetaDxi =dQ1 * Fermi_1 +dQ2* Fermi_2 +dQ3* Fermi_3 +dQ4* Fermi_4
             + Q1 *dFermi_1 + Q2*dFermi_2 + Q3*dFermi_3 + Q4*dFermi_4;
    
    press = Kzeta*exp(zeta/Vlog10); 
    dPdE = (press / Vlog10) * (DzetaDxi * Vlog10/E);
    
    *pres_p = press;
    *DPDene_p   = dPdE;

}





/* ================================================  
   Analytic fits routine : HP fit
   Cold EOS, Thermodynamical NOT consistent procedure
   ================================================ */

void HP_EOSFit_Cold( const int which_eos,
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
  double rho, epsl;
  double ene,pressure,Pcold,DPcoldDene,DpDrho,DpDepsl;
  double Ecold, DEcoldDrho, D2EcoldD2rho;
  double cs2;

  int i;

  for (i=0;i<n_elems;i++) {
    
    rho = rho_p[i];
    
    /* ****************************************
       Compute the Cold EOS
       ************************************* */
    
    HP_EpslFit(which_eos,&Ecold,&DEcoldDrho,&D2EcoldD2rho,rho);
    
    epsl = Ecold; 
    ene  = rho * (1.0 + epsl);
    
    HP_PresFit(which_eos,&Pcold,&DPcoldDene,ene);
    
    pressure = Pcold;
    
    DpDrho  = DPcoldDene * (1.0 + epsl + rho * DEcoldDrho);     
    DpDepsl = 0.0;
    
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
    epsl_p[i]        = epsl;
    DepslDrho_p[i]   = DEcoldDrho;
    D2epslD2rho_p[i] = D2EcoldD2rho;

  } /* end for */
  
}





/* ================================================  
   Analytic fits routine : HP fit al la Shibata 
   Cold EOS, Thermodynamical consistent procedure
   ================================================ */

void HPS_EOSFit_Cold( const int which_eos,
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
    
    HP_EpslFit(which_eos,&Ecold,&DEcoldDrho,&D2EcoldD2rho,rho);
    
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


