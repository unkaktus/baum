/* eos_tab1D.h */
/* mth 06/12 */

#include "bam.h"
#include "eos.h"


#define PR 0
#define ln10 2.302585092994046e+00


// arrays for EoS 1d tables
static double *eostab_Logrho  = NULL;
static double *eostab_Logepsl = NULL; 
static double *eostab_DLogepsl = NULL; 
//static double *eostab_Logpre  = NULL;
static int eostab_size;



/*****************************************************************************/
/* interpolation routines with derivatives */

/* Locate the index of the nearest xval element in the array x */
int locate1d( double *x, int Nx, double xval)
{
  int ju,jm,jl;
  int ascnd;
  
  jl=-1;
  ju=Nx;
  
  if (xval <= x[0]) { 
    if (xval < x[0]) if (PR) printf("  eos_table warning:  pt to locate is outside (xval<xx).\n"); 
    return 0; 
  } else if (xval >= x[Nx-1]) { 
    if (xval > x[Nx-1])if (PR)  printf("  eos_table warning:  pt to locate is outside (xval>xx).\n"); 
    return Nx-1; 
  }
  
  ascnd = (x[Nx-1] >= x[0]);
  
  while (ju-jl > 1) {
    
    jm = (ju+jl) >> 1;
    
    if ((xval >= x[jm]) == ascnd)
      jl=jm;
    else
      ju=jm;
    
  }
  
  return jl;
}

/* Hermite 2 pts interpolation with drvts */
void interp1d_her4( double *f, double *df, double *x, int Nx, double xv, double *fv_p, double *dfv_p, double *ddfv_p )
{
  /* Given the values in xv, it returns the interpolated values fv and 
  its 1st and 2nd derivatives dfv, ddfv of the fuction f(x) 
  Hermite 2 pts interpolation is used */
  
  if (Nx < 2) errorexit(" too few points for interpolation");
  
  int i = locate1d(x,Nx,xv);
  
  if( i < 0 ){ 
    if (PR) printf(" too few points on the left => interpolation maybe be inaccurate! (rho=%e)\n",pow(10.,xv));
    i = 1;
  }
  if( i > (Nx-1) ){ 
    if (PR) printf(" too few points on the right => interpolation maybe be inaccurate! (rho=%e)\n",pow(10.,xv));
    i = Nx-1;
  }
  
  double xi    =  x[i]; 
  double fi    =  f[i];  
  double dfi   = df[i];  
  
  double xipo  =  x[i+1];
  double fipo  =  f[i+1];
  double dfipo = df[i+1];
  
  double w   = (xv-xi)/(xipo-xi);
  double w2  = w*w;
  double w3  = w2*w;
  double h0w = 2.*w3 - 3.*w2 + 1.;
  double h1w = w3 - 2.*w2 + w;
  
  double mw   = 1. - w;
  double mw2  = mw*mw;
  double mw3  = mw2*mw;
  double h0mw = 2.*mw3 - 3.*mw2 + 1.;
  double h1mw = mw3 - 2.*mw2 + mw;
  
  double dw = 1./(xipo-xi);
  
  *fv_p   = fi*h0w + fipo*h0mw + (xipo-xi)*(dfi*h1w - dfipo*h1mw);
  
  *dfv_p  = 6.*( fi*(w2-w) - fipo*(mw2-mw) )*dw 
      + dfi*(3.*w2-4.*w+1.) + dfipo*(3.*mw2-4.*mw+1.);
  
  *ddfv_p = ( fi*(12.*w-6.) + fipo*(12.*mw-6.) )*dw*dw 
      + ( dfi*(6.*w-4.) - dfipo*(6.*mw-4.) )*dw;
  
}

/* Lagrange 4 pts interpolation with drvts */
void interp1d_lag4( double *f, double *df, double *x, int Nx, double xv, double *fv_p, double *dfv_p, double *ddfv_p )
{
  /* Given the values in xv, it returns the interpolated values fv and 
  its 1st and 2nd derivatives dfv, ddfv of the fuction f(x) 
  Lagrangian 4 pts interpolation is used */
  
  if (Nx < 4) errorexit(" too few points for interpolation");
  
  int i = locate1d(x,Nx,xv);
  
  if( i < 1 ){ 
    if (PR) printf(" too few points on the left => interpolation maybe be inaccurate! (rho=%e)\n",pow(10.,xv));
    i = 1;
  }
  if( i > (Nx-3) ){ 
    if (PR) printf(" too few points on the right => interpolation maybe be inaccurate! (rho=%e)\n",pow(10.,xv));
    i = Nx-3;
  }
  
  double ximo =  x[i-1];
  double xi   =  x[i];
  double xipo =  x[i+1]; 
  double xipt =  x[i+2]; 
  
  double C1   = (f[i] - f[i-1])/(xi - ximo);
  double C2   = (-f[i] + f[i+1])/(-xi + xipo);
  double C3   = (-f[i+1] + f[i+2])/(-xipo + xipt);
  double CC1  = (-C1 + C2)/(-ximo + xipo);
  double CC2  = (-C2 + C3)/(-xi + xipt);
  double CCC1 = (-CC1 + CC2)/(-ximo + xipt);
  
  *fv_p   = f[i-1] + (-ximo + xv)*(C1 + (-xi + xv)*(CC1 + CCC1*(-xipo + xv)));
  *dfv_p  = C1 - (CC1 - CCC1*(xi + xipo - 2.*xv))*(ximo - xv)
      + (-xi + xv)*(CC1 + CCC1*(-xipo + xv));
  *ddfv_p = 2.*(CC1 - CCC1*(xi + ximo + xipo - 3.*xv));
}

/* Linear interpolation with drvts */
void interp1d_lag1( double *f, double *df, double *x, int Nx, double xv, double *fv_p, double *dfv_p, double *ddfv_p )
{
  
  /* Given the values in xv, it returns the interpolated values fv and 
  its 1st and 2nd derivatives dfv, ddfv of the fuction f(x) 
  Linear interpolation is used (check and testing) */
  
  if (Nx < 2) errorexit(" too few points for interpolation");
  
  int i;
  i = locate1d(x,Nx,xv);
  
  if( i < 0 ){ 
    if (PR) printf(" too few points on the left => interpolation maybe be inaccurate! (rho=%e)\n",pow(10.,xv));
    i = 0;
  }
  if( i > (Nx-1) ){ 
    if (PR) printf(" too few points on the right => interpolation maybe be inaccurate! (rho=%e)\n",pow(10.,xv));
    i = Nx-1;
  }
  
  double xi   =  x[i];
  double xipo =  x[i+1]; 
  
  double fi   =  f[i];
  double fipo =  f[i+1]; 
  
  *dfv_p  = (fipo-fi)/(xipo-xi);
  *fv_p   = fi + *dfv_p*(xv - xi);    
  *ddfv_p = 0.;
}

// wrapper 
void interp_eos_1d(double Logrho, double *LogE, double *LogDEDrho, double *LogD2EDrho2) { 
  /* interp_eostab1d_lag4(eostab_Logepsl,eostab_Logrho,eostab_size,Logrho,LogE,LogDEDrho,LogD2EDrho2); */
  interp1d_her4(eostab_Logepsl,eostab_DLogepsl,eostab_Logrho,eostab_size,Logrho,LogE,LogDEDrho,LogD2EDrho2); 
  /* interp_eostab1d_lag1(eostab_Logepsl,eostab_Logrho,eostab_size,Logrho,LogE,LogDEDrho,LogD2EDrho2); */
}






/*****************************************************************************/
/*  load 1d tables */
void eos_load_tab1d()
{

  /* 
  1D tables 
  -----------------------------------------------------
     Tables must contain the line: 
     
  172  <-- Number of lines
     
     and 4 cols:
     
  index  n_B [fm^{-3}]  e [g/cm^3]   p [dyn/cm^2]
     
  n_B = baryon number density, 
  e = total energy density,
  p = pressure.
     
  lines starting with # are comments
     
  !!! IMPORTANT !!!
  to use routine "grhd_locate" below tabs have to be ordered in rho (n_B)
     
  */
  

  FILE *pf;
  
  if (GetsLax("eos_tab_file")==0) errorexit(" need eos_tab_file");
  
  char *fname = Gets("eos_tab_file");
  printf(" read in eos table file:     %s\n", fname);
  
  char s[1024];
  
  //  units  
  static double UTS_mdens_cgs = 6.1764e+17; // g cm^-3 
  static double UTS_pres_cgs  = 5.5511e+38; // dyn cm^-2 
  static double UTS_n_cgs     = 1e+39;      // cm^-3 
  static double mB            = 1.66e-24;   // g 
  
  int i  = 0;
  
  int tmpi;
  double nB, rho, ene, epsl, pres;
  
  // open file 
  pf = fopen(fname, "r");
  if (!pf) 
    errorexits(" could not open file %s", fname);
  
  // import tabs 
  int firstcall = 1;
  while (fgets(s, 1024, pf)) {
    
    // comments 
    if (s[0] == '#') continue; 
    
    // no entires 
    if (firstcall) {
      firstcall = 0;
      sscanf(s, "%d <-- Number of lines",&eostab_size);
      // allocate tab arrays 
      eostab_Logrho    =  malloc( eostab_size*sizeof(double) );
      eostab_Logepsl   =  malloc( eostab_size*sizeof(double) );
      eostab_DLogepsl  =  malloc( eostab_size*sizeof(double) );
      //eostab_Logpre  =  malloc( eostab_size*sizeof(double) );
      continue; 
    }
    
    // import data 
    if (fscanf(pf,"%d %lf %lf %lf",&tmpi,&nB,&ene,&pres)==0) errorexit("problem during reading");

    // switch to c = G = Msun = 1 units 
    rho  = nB  *UTS_n_cgs*mB/UTS_mdens_cgs;
    ene  = ene /UTS_mdens_cgs;
    pres = pres/UTS_pres_cgs;
    
    epsl = ene/rho - 1.;  

    eostab_Logrho[i] = log10(rho);
    eostab_Logepsl[i]= log10(epsl);
    
    // the following quantity is
    //
    // d Log epsl   rho  d epsl      P
    // ---------- = ---- ------ = ------- 
    // d Log rho    epsl d rho    epsl rho
    //
    // the last line uses 1st thermodynamics law @ T=0
    // P =  rho^2 epsl'(rho)
    //
    eostab_DLogepsl[i]= pres/(epsl*rho);
    
    //eostab_Logpre[i] = log10(pres);
    
    if (PR) printf("%d   %e %e %e  ->  %e %e %e\n",i,rho,epsl,pres,eostab_Logrho[i],eostab_Logepsl[i],eostab_DLogepsl[i]);//,eostab_Logpre[i]);
    
    i++;
    if (i == eostab_size) break;
    
  }
  
  // close file 
  fclose (pf);

  // check order of tab: locate routine assumes monotonic rho 
  for (i=1; i<eostab_size-1; i++) {
    if ((eostab_Logrho[i]-eostab_Logrho[i-1])<0.) 
      // for now we exit... 
      // todo: implement better method to locate tab position 
      errorexit(" EoS tab must be ordered with monotonic n_B");    
  }

}


// Extract limits of the table 
void eos_tab1d_Logrhobnd(double *Logrho_max,double *Logrho_min)
{
  *Logrho_max = eostab_Logrho[eostab_size];
  *Logrho_min = eostab_Logrho[0];
}









/*****************************************************************************/
/* interpolate 1d EoS Tab (Termodynamically consistent procedure) */
int eos_cold_tab1d(double *p, double *cs2, double *dpdrho, double *dpdepsl, double *rho, double *epsl) 
{
  double LogE, LogDEDrho, LogD2EDrho2; /* tmp Log eps and drvts */
  double E, DEDrho, D2EDrho2;          /* tmp eps and drvts */
  
  double Logrho = log10(*rho);
  if (!finite(Logrho)) return 1;
  
  interp_eos_1d(Logrho,&LogE,&LogDEDrho,&LogD2EDrho2);
  
  E        = pow(10., LogE);// this is epsl from tab
  DEDrho   = E*LogDEDrho/(*rho); 
  D2EDrho2 = (E*LogD2EDrho2/(ln10) + LogDEDrho*(DEDrho*(*rho)-E))/((*rho)*(*rho)); 
  
  *p      = (*rho)*(*rho)*DEDrho;
  *epsl   = E; // reset epsl here
  *dpdrho = 2.*(*rho)*DEDrho + (*rho)*(*rho)*D2EDrho2;
  *dpdepsl= 0.;    
  *cs2    = eos_cs2_rep(*rho,E, *p,*dpdepsl,*dpdrho);
  
  if (PR) printf("%e %e -> %e %e %e -> %e %e %e %e\n", *rho,*epsl, E,DEDrho,D2EDrho2, *p,*dpdrho,*dpdepsl,*cs2);
  if (!finite(*rho) || !finite(*epsl) || !finite(*cs2)) return 1;
  return 0;
}

































