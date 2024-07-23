/* eos_analytic.c */
/* mth 06/12, sbernuz 08/12 */

#include "bam.h"
#include "eos.h"

#define PR 1


/*****************************************************************************/
/* general routine for sound speed for p,rho,epsl */
double eos_cs2_rep(double rho, double epsl, double p, double dpdrho, double dpdepsl) 
{
  if (rho==0.)
    return 0.;
  else
    return( ( dpdrho + p*dpdepsl/(rho*rho) )/( 1.0 + epsl + p/rho ) );
}




/*****************************************************************************/
/* no EoS, causes an error */
int eos_dummy() 
{
  errorexit("this eos function pointer is not set");
  return 0;
}




/*****************************************************************************/
/* dust EoS */
int eos_dust(double *p, double *cs2, double *dpdrho, double *dpdepsl, double *rho, double *epsl) 
{
  *p      = 0.;
  *dpdepsl= 0.;
  *dpdrho = 0.;
  *cs2    = 0.;
  *epsl   = 0.;
  
  return 0;
}




/*****************************************************************************/
/* polytropic EoS */
int eos_poly(double *p, double *cs2, double *dpdrho, double *dpdepsl, double *rho, double *epsl) 
{

  double krhotog1 = EOS.K * pow( *rho, EOS.GAMMA-1.); // compute less pow !

  *p      = *rho * krhotog1;
  *dpdepsl= 0.;
  *dpdrho = (EOS.GAMMA) * krhotog1;
  if (*epsl<=-1.) {
    // in this case rest epsl
    // (e.g. used when this eos is the cold part of ideal gas)
    *epsl = krhotog1 /(EOS.GAMMAMO);
  }
  *cs2    = eos_cs2_rep(*rho,*epsl, *p,*dpdrho,*dpdepsl);

  if (!finite(*rho) || !finite(*epsl) || !finite(*cs2)) return 1;
  return 0;
}




/*****************************************************************************/
/* polytropic EoS for P(H) */
int eos_poly_H(double *h, double *rho, double *p, double *epsl)
{

  double K = EOS.K;
  double G1 = EOS.GAMMA-1.;

  double H1 = 2.*sinh(0.5*(*h))*exp(0.5*(*h)); // exp(H)-1
  double n  = pow( G1/(EOS.GAMMA*K)*H1, 1.0/G1 ); // = rho
  double krhotog1 = K*pow(n,G1);

  *p    = n * krhotog1;
  *rho  = n;
  *epsl = krhotog1 / G1;
  
  if (!finite(*rho) || !finite(*p) || !finite(*epsl)) return 1;
  return 0;
}




/*****************************************************************************/
/* piecewise polytropic EoS */
/* Refs. 
   Read et al. Phys.Rev.D79:124032,2009 
   Hotokezaka et al Phys.Rev.D83:124008,2011 (2011) 1105.4370 
*/

double epsltreshold = 0.;

void eos_load_pwp()
{
  const int pr = 1;
  
  tU units;
  set_units(&units);
     
  char *fname = Gets("eos_tab_file");
  char sline[1024], spar[512],sval[512];
  
  // open file 
  if (pr) printf(" reading PWP pars from file:\n    %s\n", fname);
  FILE *ifp = fopen(fname, "r");
  if (!ifp) errorexits(" could not open file %s", fname);
  
  /* read the file */
  int i = 0;
  while (fgets(sline,1024,ifp)) {
    
    // skip comments
    if (sline[0] == '#') 
      continue;
    
    // read pars line by line
    if (strrchr(sline,'=')){
      sscanf(sline,"%s = %s\n",spar,sval);
      if         (strncmp(spar,"EOS_PWP_N",9)==0)      {
        EOS.PWP_N     = atoi(sval);     
	EOS.PWP_N++;
        EOS.PWP_RHO   = (double*) malloc ((EOS.PWP_N)*sizeof(double));
        EOS.PWP_K     = (double*) malloc ((EOS.PWP_N)*sizeof(double));
        EOS.PWP_GAMMA = (double*) malloc ((EOS.PWP_N)*sizeof(double));
        EOS.PWP_a     = (double*) malloc ((EOS.PWP_N)*sizeof(double));
        EOS.PWP_p     = (double*) malloc ((EOS.PWP_N)*sizeof(double));
        EOS.PWP_epsl  = (double*) malloc ((EOS.PWP_N)*sizeof(double));
        EOS.PWP_eta   = (double*) malloc ((EOS.PWP_N)*sizeof(double));
      } else if (strncmp(spar,"EOS_PWP_LG10P1",14)==0) {
	EOS.PWP_LG10P1 = atof(sval);
      } else if (strncmp(spar,"EOS_PWP_LG10R1",14)==0) {
	EOS.PWP_LG10R1 = atof(sval);
      } else if (strncmp(spar,"EOS_PWP_K0",10)==0) {
        EOS.PWP_K0 = atof(sval);
      } else if (strncmp(spar,"EOS_PWP_GAMMA0",14)==0) {
        EOS.PWP_GAMMA0  = atof(sval);
      } else {
        errorexit(" problem during reading PWP eos file");
      }
      continue;
    }
    
    // read table pars
    i++;
    sscanf(sline,"%lf %lf",&(EOS.PWP_RHO[i]),&(EOS.PWP_GAMMA[i]));
    if (i==EOS.PWP_N-1) break;
  }
  
  // close file
  fclose (ifp);
  

  /* compute parameters */
  
  EOS.PWP_GAMMA[0] = EOS.PWP_GAMMA0;
  
  double lg10rho0,lg10rho1, eps_imo,eps_i, krhotog1_imo,krhotog1_i;
  
  lg10rho1 = EOS.PWP_LG10R1; // log10 fiducial density "rho1"
  lg10rho0 = ( EOS.PWP_LG10P1 - lg10rho1*EOS.PWP_GAMMA[1] - log10(EOS.PWP_K0) ) 
    / ( EOS.PWP_GAMMA0-EOS.PWP_GAMMA[1] );

  if (pr) {
    
    printf(" lg10R0 %.16e lg10R1 %e\n",lg10rho0,lg10rho1);
    printf("     R0 %.16e     R1 %e\n",pow(10.,lg10rho0),pow(10.,lg10rho1));
    printf("     R0 %e     R1 %e [dimensionless]\n",pow(10.,lg10rho0)/units.Mdens_cgs,pow(10.,lg10rho1)/units.Mdens_cgs);
     
  }

  // dividing densities 
  // note that RHO[i] has different indexes convention than in the literature ! 
  EOS.PWP_RHO[0]  = 0.;
  EOS.PWP_RHO[1]  = pow( 10., lg10rho0 );
  EOS.PWP_RHO[2]  = pow( 10., lg10rho1 );

  // switch to dimensionless units to rho and K0
  // other quantities follow 
  double Kgcm3 = pow( units.Mdens_cgs, 1.-EOS.PWP_GAMMA0 );
  EOS.PWP_K0 /= Kgcm3;
  for (i=0; i<EOS.PWP_N; i++) EOS.PWP_RHO[i] /= units.Mdens_cgs;

  // polytropic constants 
  EOS.PWP_K[0] = EOS.PWP_K0;
  for (i=1; i<EOS.PWP_N; i++)
    EOS.PWP_K[i] = EOS.PWP_K[i-1]*pow( EOS.PWP_RHO[i],EOS.PWP_GAMMA[i-1]-EOS.PWP_GAMMA[i]);

  // compute other quantities at diving densities
  EOS.PWP_epsl[0]  = 0.;
  EOS.PWP_a[0]     = 0.;
  EOS.PWP_eta[0]   = 0.; 
  EOS.PWP_p[0]     = 0.;
  
  for (i=1; i<EOS.PWP_N; i++) {
    
    krhotog1_imo = EOS.PWP_K[i-1] * pow( EOS.PWP_RHO[i],EOS.PWP_GAMMA[i-1] -1);
    krhotog1_i   = EOS.PWP_K[i]   * pow( EOS.PWP_RHO[i],EOS.PWP_GAMMA[i]   -1);

    eps_imo = krhotog1_imo/(EOS.PWP_GAMMA[i-1]-1);
    eps_i   = krhotog1_i  /(EOS.PWP_GAMMA[i]  -1);
    
    // a_i    
    EOS.PWP_a[i] = EOS.PWP_a[i-1] + eps_imo - eps_i;
    
    // p_i, eps_i, eta_i
    EOS.PWP_p[i]    = EOS.PWP_RHO[i] * krhotog1_i;    
    EOS.PWP_epsl[i] = EOS.PWP_a[i] + eps_i;
    EOS.PWP_eta[i]  = EOS.PWP_a[i] + eps_i * EOS.PWP_GAMMA[i];

  }

  if (pr) {
    printf("  Piecewise polytropic parameters - Dimensionless units\n");
    printf("  conv fact [rho] %e [K] %e\n",units.Mdens_cgs,Kgcm3);
    printf("  i    rho             K             Gamma           a   ");
    printf("              p              epsl            eta\n");
    for (i=0; i<EOS.PWP_N; i++) 
      printf("  %d %.8e %.8e %.8e %.8e   %.8e %.8e %.8e\n",
             i,
             EOS.PWP_RHO[i],EOS.PWP_K[i],EOS.PWP_GAMMA[i],EOS.PWP_a[i],
             EOS.PWP_p[i],EOS.PWP_epsl[i],EOS.PWP_eta[i]);
  }


  // set "epsltreshold"
  double p,cs2,dpdrho,dpdepsl, rho;
  if (Getv("eos","pwphot")+GetvLax("grhd_use_atmosphere","yes")+GetvLax("grmhd_use_atmosphere","yes")) {
    epsltreshold = -10.;
    rho = GetdLax("grhd_atm_level"); // * GetdLax("grhd_atm_factor"); 
    if(Getv("physics", "grmhd")) rho = GetdLax("grmhd_atm_level"); // * GetdLax("grmhd_atm_factor"); 
    eos_pwp(&p,&cs2,&dpdrho,&dpdepsl, &rho, &epsltreshold);
    if (pr) printf("set epsltreshold = %4.3e (rho = %4.3e)\n",epsltreshold,rho);
  }


  double dummy, c; 
  rho = 1e-15;
  FILE *fp = fopen("eos_pwp", "w+");
  if(PR) {
        c = 1./999.*(log10(6e-03/1e-15));
        while(rho<6e-03){
                eos_pwp(&p, &dummy, &dummy, &dummy, &rho, &dummy);
                fprintf(fp, "%le %le\n", rho, p);
                rho = rho*pow(10,c);
        }
        fclose(fp);
  }
  // debug:
  //errorexit(" stop");
  
}

int eos_pwp(double *p, double *cs2, double *dpdrho, double *dpdepsl, double *rho, double *epsl) 
{
  // find segment
  int i = 0;
  while (i<EOS.PWP_N-1 && *rho>EOS.PWP_RHO[i+1]) {
    i++;
  }
    
  //for (i=0; i<EOS.PWP_N; i++) if (*rho <= EOS.PWP_RHO[i] ) break;
  
  double krhotog1 = EOS.PWP_K[i] * pow((*rho),EOS.PWP_GAMMA[i]-1);

  *p      = (*rho) * krhotog1;
  if (*epsl<=-1.) {
    *epsl   = EOS.PWP_a[i] + krhotog1/(EOS.PWP_GAMMA[i]-1);
  }
  *dpdrho = EOS.PWP_GAMMA[i] * krhotog1;
  *dpdepsl= 0.;
  *cs2    = eos_cs2_rep(*rho,*epsl, *p,*dpdrho,*dpdepsl);
  
  if (!finite(*rho) || !finite(*epsl) || !finite(*cs2)) return 1; 
  return 0;
}

int eos_pwp_H(double *H, double *rho, double *p, double *epsl) 
{
  int i;
  double eta = *H -1;
  
  // see eq A.8-10 Read et al.
  
  // routine needs better check

  // find segment
  for (i=0; i<EOS.PWP_N; i++) if (eta < EOS.PWP_eta[i]) break;
  
  double Ni = 1./(EOS.PWP_GAMMA[i]-1);

  double tmp  = (eta-EOS.PWP_a[i])/(EOS.PWP_K[i]*(Ni+1));
  double tmpn = pow( tmp, Ni ); 

  *p    = EOS.PWP_K[i] * tmp *tmpn;
  *rho  = tmpn;
  *epsl = (EOS.PWP_a[i] + Ni*eta)/(Ni+1);


  /*

  int pr = 0;
  
  int i,j;
  double pres, f,df;

  // find segment
  for (i=0; i<EOS.PWP_N; i++) if (*H < EOS.PWP_H[i+1]) break;

  if (pr) {
    printf("======================================\n");
    printf("  find p(H):\n");
    printf("  H=%e   i=%d   H_i-1=%e H_i=%e\n",*H,i,EOS.PWP_H[i],EOS.PWP_H[i+1]);
  }
  
  // find it in the first integral from p1=0
  pres = EOS.PWP_p[i];
  j = 0;
  do {
    
    f    = *H - (EOS.PWP_H[i] + comp_H_pwp( EOS.PWP_p[i],pres, EOS.PWP_GAMMA[i],EOS.PWP_K[i]));
    df   = 1./(pres + pow( pres/EOS.PWP_K[i] , 1./EOS.PWP_GAMMA[i]));
    
    if (pr) {
      printf("  NR %d:   p=%e    f=%e df=%e f/df=%e\n",j,pres,f,df, f/df);
    }
    
    pres = pres + f/df;
    j++;
  } while (fabs(f/df)>1e-12);
  
  // compute primitives with the stored values
  *rho   = pow( pres/EOS.PWP_K[i] , 1./EOS.PWP_GAMMA[i]);
  *p     = pres;
  *epsl  = (1.+EOS.PWP_a[i])*(*rho) + EOS.PWP_K[i]/(EOS.PWP_GAMMA[i]-1) * pow( *rho, EOS.PWP_GAMMA[i]);
  // *e     = (*rho)*(1.+(*epsl));
  
  if (pr) {
    printf("  -> p = %e   rho = %e\n", *p,*rho);
    printf("======================================\n");
    errorexit("");
  }
  
  
  if (!finite(*rho) || !finite(*p)) return 1;
  return 0;
  */

}

int eos_pwpHot(double *p, double *cs2, double *dpdrho, double *dpdepsl, double *rho, double *epsl) 
{
  
  double ec, pres, ccc, kkk, eh, krhotog1;

  //printf(1," here \n ");

  // find segment
  int i = 0;
  while (i<EOS.PWP_N-1 && *rho>EOS.PWP_RHO[i+1]) {
    i++;
  }
    
  krhotog1 = EOS.PWP_K[i] * pow((*rho),EOS.PWP_GAMMA[i]-1);

  pres = krhotog1*(*rho);
  ccc  = EOS.PWP_GAMMA[i] * krhotog1;
  ec   = EOS.PWP_a[i] +  krhotog1/(EOS.PWP_GAMMA[i]-1.); 

  eh    = *epsl - ec;

  // do not allow hot part if unphysical
  //if (eh<0.) eh = 0.;
  // alt., use atm level (need testing)
  if (eh<epsltreshold) eh = 0.; 

  ccc  += (EOS.GAMMAMO)*eh;
  kkk   = (EOS.GAMMAMO)*(*rho);
  pres += kkk*eh;

  *p      = pres;
  *dpdrho = ccc;
  *dpdepsl= kkk;
  *cs2    = eos_cs2_rep(*rho,*epsl, pres,ccc,kkk);

  if (!finite(*rho) || !finite(*epsl) || !finite(*cs2)) return 1; 
  return 0;

}

int eos_pwp_p (double *p, double *cs2, double *drhodp, double *dedp, double *rho, double *e)
{

 int i=0;
 double epsl, despdr;

 while( i<EOS.PWP_N-1 && *p > EOS.PWP_p[i+1] ) i++;

 *rho	 = pow( (*p/EOS.PWP_K[i]) , (1./EOS.PWP_GAMMA[i]) );
 *e  	 = (1.+ EOS.PWP_a[i])*(*rho) + *p/(EOS.PWP_GAMMA[i] - 1.);
 *drhodp = *rho/(EOS.PWP_GAMMA[i]*(*p));
 *dedp   = 1./(EOS.PWP_GAMMA[i]*(*p))*(1. + EOS.PWP_a[i])*pow( (*p/EOS.PWP_K[i]), (1./EOS.PWP_GAMMA[i]) ) + 1./(EOS.PWP_GAMMA[i] - 1.);

 return 0;

}


/*****************************************************************************/
/* ideal gas */
int eos_ideal(double *p, double *cs2, double *dpdrho, double *dpdepsl, double *rho, double *epsl) 
{
  
  *p      = EOS.GAMMAMO*(*rho)*(*epsl);
  *dpdepsl= EOS.GAMMAMO*(*rho);
  *dpdrho = EOS.GAMMAMO*(*epsl);
  *cs2    = eos_cs2_rep(*rho,*epsl, *p,*dpdrho,*dpdepsl);
  
  if (!finite(*rho) || !finite(*epsl) || !finite(*cs2)) return 1;
  return 0;
}

/**************************************************************************/
/* HG: Altered routine to load pwps and fix the energy gap in the ID */
void eos_load_pwp_p()
{
  const int pr = 1;
  
  tU units;
  set_units(&units);

  char *fname = Gets("eos_pwp_ID");
  char sline[1024], spar[512],sval[512];
  
  // open file 
  if (pr) printf(" reading PWP pars from file:\n    %s\n", fname);
  FILE *ifp = fopen(fname, "r");
  if (!ifp) errorexits(" could not open file %s", fname);
  
  /* read the file */
  int i = 0;
  while (fgets(sline,1024,ifp)) {
    
    // skip comments
    if (sline[0] == '#') 
      continue;
    
    // read pars line by line
    if (strrchr(sline,'=')){
      sscanf(sline,"%s = %s\n",spar,sval);
      if         (strncmp(spar,"EOS_PWP_N",9)==0)      {
        EOS.PWP_N     = atoi(sval);     
	EOS.PWP_N++;
        EOS.PWP_RHO   = (double*) malloc ((EOS.PWP_N)*sizeof(double));
        EOS.PWP_K     = (double*) malloc ((EOS.PWP_N)*sizeof(double));
        EOS.PWP_GAMMA = (double*) malloc ((EOS.PWP_N)*sizeof(double));
        EOS.PWP_a     = (double*) malloc ((EOS.PWP_N)*sizeof(double));
        EOS.PWP_p     = (double*) malloc ((EOS.PWP_N)*sizeof(double));
        EOS.PWP_epsl  = (double*) malloc ((EOS.PWP_N)*sizeof(double));
        EOS.PWP_eta   = (double*) malloc ((EOS.PWP_N)*sizeof(double));
      } else if (strncmp(spar,"EOS_PWP_LG10P1",14)==0) {
	EOS.PWP_LG10P1 = atof(sval);
      } else if (strncmp(spar,"EOS_PWP_LG10R1",14)==0) {
	EOS.PWP_LG10R1 = atof(sval);
      } else if (strncmp(spar,"EOS_PWP_K0",10)==0) {
        EOS.PWP_K0 = atof(sval);
      } else if (strncmp(spar,"EOS_PWP_GAMMA0",14)==0) {
        EOS.PWP_GAMMA0  = atof(sval);
      } else {
        errorexit(" problem during reading PWP eos file");
      }
      continue;
    }
    
    // read table pars
    i++;
    sscanf(sline,"%lf %lf",&(EOS.PWP_RHO[i]),&(EOS.PWP_GAMMA[i]));
    if (i==EOS.PWP_N-1) break;
  }
  
  // close file
  fclose (ifp);
  

  /* compute parameters */
  
  EOS.PWP_GAMMA[0] = EOS.PWP_GAMMA0;
  
  double lg10rho0,lg10rho1, eps_imo,eps_i, krhotog1_imo,krhotog1_i;
  
  lg10rho1 = EOS.PWP_LG10R1; // log10 fiducial density "rho1"
  lg10rho0 = ( EOS.PWP_LG10P1 - lg10rho1*EOS.PWP_GAMMA[1] - log10(EOS.PWP_K0) ) 
    / ( EOS.PWP_GAMMA0-EOS.PWP_GAMMA[1] );

  if (pr) {
    
    printf(" lg10R0 %e lg10R1 %e\n",lg10rho0,lg10rho1);
    printf("     R0 %e     R1 %e\n",pow(10.,lg10rho0),pow(10.,lg10rho1));
    printf("     R0 %e     R1 %e [dimensionless]\n",pow(10.,lg10rho0)/units.Mdens_cgs,pow(10.,lg10rho1)/units.Mdens_cgs);
     
  }

  // dividing densities 
  // note that RHO[i] has different indexes convention than in the literature ! 
  EOS.PWP_RHO[0]  = 0.;
  EOS.PWP_RHO[1]  = pow( 10., lg10rho0 );
  EOS.PWP_RHO[2]  = pow( 10., lg10rho1 );

  // switch to dimensionless units to rho and K0
  // other quantities follow 
  double Kgcm3 = pow( units.Mdens_cgs, 1.-EOS.PWP_GAMMA0 );
  EOS.PWP_K0 /= Kgcm3;
  for (i=0; i<EOS.PWP_N; i++) EOS.PWP_RHO[i] /= units.Mdens_cgs;

  // polytropic constants 
  EOS.PWP_K[0] = EOS.PWP_K0;
  for (i=1; i<EOS.PWP_N; i++)
    EOS.PWP_K[i] = EOS.PWP_K[i-1]*pow( EOS.PWP_RHO[i],EOS.PWP_GAMMA[i-1]-EOS.PWP_GAMMA[i]);

  // compute other quantities at diving densities
  EOS.PWP_epsl[0]  = 0.;
  EOS.PWP_a[0]     = 0.;
  EOS.PWP_eta[0]   = 0.; 
  EOS.PWP_p[0]     = 0.;
  
  for (i=1; i<EOS.PWP_N; i++) {
    
    krhotog1_imo = EOS.PWP_K[i-1] * pow( EOS.PWP_RHO[i],EOS.PWP_GAMMA[i-1] -1);
    krhotog1_i   = EOS.PWP_K[i]   * pow( EOS.PWP_RHO[i],EOS.PWP_GAMMA[i]   -1);

    eps_imo = krhotog1_imo/(EOS.PWP_GAMMA[i-1]-1);
    eps_i   = krhotog1_i  /(EOS.PWP_GAMMA[i]  -1);
    
    // a_i    
    EOS.PWP_a[i] = EOS.PWP_a[i-1] + eps_imo - eps_i;
    
    // p_i, eps_i, eta_i
    EOS.PWP_p[i]    = EOS.PWP_RHO[i] * krhotog1_i;    
    EOS.PWP_epsl[i] = EOS.PWP_a[i] + eps_i;
    EOS.PWP_eta[i]  = EOS.PWP_a[i] + eps_i * EOS.PWP_GAMMA[i];

  }

  if (pr) {
    printf("  Piecewise polytropic parameters - Dimensionless units\n");
    printf("  conv fact [rho] %e [K] %e\n",units.Mdens_cgs,Kgcm3);
    printf("  i    rho             K             Gamma           a   ");
    printf("              p              epsl            eta\n");
    for (i=0; i<EOS.PWP_N; i++) 
      printf("  %d %.8e %.8e %.8e %.8e   %.8e %.8e %.8e\n",
             i,
             EOS.PWP_RHO[i],EOS.PWP_K[i],EOS.PWP_GAMMA[i],EOS.PWP_a[i],
             EOS.PWP_p[i],EOS.PWP_epsl[i],EOS.PWP_eta[i]);
  }


  // set "epsltreshold"
  double p,cs2,dpdrho,dpdepsl, rho;
  if (Getv("eos","pwphot")+GetvLax("grhd_use_atmosphere","yes")+GetvLax("grmhd_use_atmosphere","yes")) {
    epsltreshold = -10.;
    rho = GetdLax("grhd_atm_level"); // * GetdLax("grhd_atm_factor"); 
    if(Getv("physics", "grmhd")) rho = GetdLax("grmhd_atm_level"); // * GetdLax("grmhd_atm_factor"); 
    eos_pwp(&p,&cs2,&dpdrho,&dpdepsl, &rho, &epsltreshold);
    if (pr) printf("set epsltreshold = %4.3e (rho = %4.3e)\n",epsltreshold,rho);
  }


  // debug:
  //errorexit(" stop");
  
}


// FS: functions for adiabatic index
double EOS_Gamma(double rho, double p, double epsl, double cs2) {
  if (EOS.type==PWPHOT || EOS.type==PWP) {
    double rhoh = rho*(1 + epsl) + p;
    return cs2*rhoh/p;
  }
  else if (EOS.type==IDEAL) return EOS.GAMMA;
  else errorexit("HLLC Riemann solver does not work with the specified EOS \n");
}

double EOS3D_Gamma(double rho, double Y, double T) {
  double p, dpdrho, dpdepsl, dummy;
  if (EOS.type==TAB3D) {
    EOS.use3DT(&rho,&dummy,&Y,&p,&T,&dummy,&dpdrho,&dpdepsl);
    return rho/p*dpdrho + dpdepsl/rho;
  }
  else errorexit("HLLC Riemann solver does not work with the specified EOS \n");
}