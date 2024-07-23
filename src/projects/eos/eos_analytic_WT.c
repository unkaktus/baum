/* eos_analytic_WT.c */
/* Wolfgang Tichy 9/2015  */

#include "bam.h"
#include "eos.h"


/******************************************************************/    
/* In this file we collect EOS routines that never introduce NANs */
/******************************************************************/    


/**************************************************************************/
/* general routine that returns sound speed given p,rho,epsl and p derivs */
double eos_cs2_rep_WT(double rho, double epsl, double p, 
                      double dpdrho, double dpdepsl)
{
  double h;
  double rho2;

  /* There is no well defined sound speed if rho=0, so return 0. */
  if(rho==0.0) return 0.0;

  /* spec. enthalpy h. h>=0 for WEC and p>=0 */
  h = 1.0 + epsl + p/rho;
  if(h==0.0) return 0.0;

  rho2 = rho*rho;
  if(rho2==0.0) return 0.0;
  /* Note: rho^2 could be rounded to zero even if rho is not zero */

  return ( dpdrho + p*dpdepsl/(rho2) ) / h;
  // return ( dpdrho + porho*(dpdepsl/rho) ) / h;
}


/* EOS made from cold piecewise polytropic EOS 
   plus a hot part that obeys a Gamma-law
   Output: double *p, double *cs2, double *dpdrho, double *dpdepsl
   Input:  double *rho, double *epsl */
int eos_pwpHot_WT(double *p, double *cs2, double *dpdrho, double *dpdepsl,
                  double *rho, double *epsl)
{
  double ec, pres, ccc, kkk, eh, krhotog1;

  /* find polytrope piece */
  int i = 0;
  while(i<EOS.PWP_N-1 && *rho>EOS.PWP_RHO[i+1])  i++;

  /* cold polytrope */    
  krhotog1 = EOS.PWP_K[i] * pow(fabs(*rho), EOS.PWP_GAMMA[i]-1);
  /* We use fabs in pow above to avoid NANs if rho<0.
     If rho<0 the pressure below will be negative, which may be a nice
     continuous continuation of it. */
  pres = krhotog1*(*rho);
  ccc  = EOS.PWP_GAMMA[i] * krhotog1;
  ec   = EOS.PWP_a[i] +  krhotog1/(EOS.PWP_GAMMA[i]-1.); 

  /* find hot epsl */
  eh    = *epsl - ec;

  /* do not allow hot part if unphysical */
  if(eh<0.0) eh=0.0;

  /* compute hot dp/drho, dp/depsl = dp/deh, p */
  ccc  += (EOS.GAMMAMO)*eh;
  kkk   = (EOS.GAMMAMO)*(*rho);
  pres += kkk*eh;

  /* set output vars */
  *p      = pres;
  *dpdrho = ccc;
  *dpdepsl= kkk;
  *cs2    = eos_cs2_rep_WT(*rho,*epsl, pres,ccc,kkk);

  return 0;
}


/* ideal gas EOS
   Output: double *p, double *cs2, double *dpdrho, double *dpdepsl
   Input:  double *rho, double *epsl */
int eos_ideal_WT(double *p, double *cs2, double *dpdrho, double *dpdepsl,
                 double *rho, double *epsl) 
{
  
  *p      = EOS.GAMMAMO*(*rho)*(*epsl);
  *dpdepsl= EOS.GAMMAMO*(*rho);
  *dpdrho = EOS.GAMMAMO*(*epsl);
  *cs2    = eos_cs2_rep_WT(*rho,*epsl, *p,*dpdrho,*dpdepsl);
//printf("*rho=%g *epsl=%g *p=%g *dpdrho=%g *dpdepsl=%g *cs2=%g\n",
//*rho,*epsl, *p,*dpdrho,*dpdepsl,*cs2);
  return 0;
}
