/* hrsc_rec1d_ppm.c 
   PPM reconstruction routines 
   harcoded for uniform grids
   
   Refs: 
   [C&W] Colella & Woodward J.Comput.Phys. 54 (1984) 
   [M&M] Marti & Mueller J.Comput.Phys. 123 (1996) 
   [M&a] Mignone et al ApJ (2005) 
   
   sbernuz */



#include "bam.h"
#include "matter.h"




/* 1D interpolation routine 
   Quartic-polynomial MC-limited interpolation. 
   See (1.6) C&W, p.177  
   Here is the pointwise MC-VanLeer limited slopes */

// void PPM_MCslope ( const double *a, double *delta_a, const int i )
// {

//   double delta_a_sign, delta_a_tmp;    
//   double tmp1,tmp2;
    
//   delta_a [i]  = 0.0;
    
//   if ( (a[i+1]-a[i])*(a[i]-a[i-1]) > 0.0 ) {
	
// 	  delta_a_tmp  = (a [i+1] - a [i-1]);
	
// 	  delta_a_sign = (delta_a_tmp>0)?1.:-1.;
	
// 	  tmp1         = DMIN( fabs ( a[i+1] - a[i] ), fabs ( a[i] - a[i-1] ) );
// 	  tmp2         = DMIN ( 0.5 *fabs (delta_a_tmp ), 2.0*tmp1);
// 	  delta_a [i]  = tmp2* delta_a_sign;
	
//   } 	
    
// }


void PPM_MCslope ( const double *a, double *delta_a, const int i )
{
  
   double S, delta_a_lim;    
   double D2m, D2c, D2p, S2, D2lim;

   const double CVL = 1.25;

   double deltamm =  a[i-1] - a[i-2];
   double deltam  =  a[i]   - a[i-1];
   double deltac  = (a[i+1] - a[i-1])*0.5;
   double deltap  =  a[i+1] - a[i];
   double deltapp =  a[i+2] - a[i+1];

   if (DMIN(deltam*deltap,deltapp*deltamm)>=0) {

      delta_a_lim = 2*DMIN(fabs(deltam),fabs(deltap));
	   S = (deltac>0)?1.:-1.;
      delta_a[i] = DMIN(fabs(deltac),fabs(delta_a_lim))*S;
   
   }

   else {

      D2m = a[i]   - 2*a[i-1] + a[i-2];
      D2c = a[i+1] - 2*a[i]   + a[i-1];
      D2p = a[i+2] - 2*a[i+1] + a[i];

      S2 = (D2c>0)?1.:-1.;

      D2lim = DMIN(fabs(D2c),DMIN(DMAX(S2*D2m,0),DMAX(S2*D2p,0)));

      if (S2*deltac<0) delta_a_lim = DMIN(CVL*1.5*D2lim,2*fabs(deltam));
      else             delta_a_lim = DMIN(CVL*1.5*D2lim,2*fabs(deltap));

	   S = (deltac>0)?1.:-1.;

      delta_a[i] = DMIN(fabs(deltac),fabs(delta_a_lim))*S;

   }
    
}


double PPM_MCslope_pt ( const double *a, const int i )
{
  
   double S, delta_a_lim;    
   double D2m, D2c, D2p, S2, D2lim;

   const double CVL = 1.25;

   double deltamm =  a[i-1] - a[i-2];
   double deltam  =  a[i]   - a[i-1];
   double deltac  = (a[i+1] - a[i-1])*0.5;
   double deltap  =  a[i+1] - a[i];
   double deltapp =  a[i+2] - a[i+1];

   if (DMIN(deltam*deltap,deltapp*deltamm)>=0) {

      delta_a_lim = 2*DMIN(fabs(deltam),fabs(deltap));
	   S = (deltac>0)?1.:-1.;
      return DMIN(fabs(deltac),fabs(delta_a_lim))*S;
   
   }

   else {

      D2m = a[i]   - 2*a[i-1] + a[i-2];
      D2c = a[i+1] - 2*a[i]   + a[i-1];
      D2p = a[i+2] - 2*a[i+1] + a[i];

      S2 = (D2c>0)?1.:-1.;

      D2lim = DMIN(fabs(D2c),DMIN(DMAX(S2*D2m,0),DMAX(S2*D2p,0)));

      if (S2*deltac<0) delta_a_lim = DMIN(CVL*1.5*D2lim,2*fabs(deltam));
      else             delta_a_lim = DMIN(CVL*1.5*D2lim,2*fabs(deltap));

	   S = (deltac>0)?1.:-1.;

      return DMIN(fabs(deltac),fabs(delta_a_lim))*S;

   }
    
}



/* 1D interpolation routine 
   Quartic-polynomial MC-limited interpolation. 
   See (1.6) C&W, p.177  
   Here is the pointwise correction to the variables 
   Note that differently from the others pointwise routines this one refers to interfaces,
   i.e. reconstructs the point x_i+1/2 (not the plus/mins edges of the ith cell) 
   the reason is to avoid double computations */

void PPM_interp1D ( const double *a, double *amins, double *aplus, 
		    const int i, 
		    const double *delta_a ) 
{
    
  aplus [i]   = 0.5 * ( a [i] + a [i+1] ) +  ( delta_a [i] - delta_a [i+1] )/6 - ( 3*(delta_a[i+1] - delta_a[i]) - (delta_a[i+2] - delta_a[i-1]) )/30;
  amins [i+1] = aplus [i];
	
}




 
/* CD steepening.
   See (1.14-17) and Appendix A p. 197 of C&W  
   Here is the pointwise correction to the variables */

void PPM_cdsteep1D ( const double *a, const double *delta_a, double *amins, double *aplus, 
		     const int i, 
		     const double epsilon1, const double eta1, const double eta2 ) 
{

  double delta_a_1_cent, delta_a_2_mins, delta_a_2_plus;
  double criterion1, criterion2;
  double eta,eta_tilde;
  double tmp;

  const double half  = 0.5;
  const double one   = 1.0;
  const double two   = 2.0;
  const double six   = 6.0;

  delta_a_2_mins = a [i] - two * a [i - 1] + a [i - 2]; 
  delta_a_2_plus = a [i + 2] - two * a [i + 1] + a [i]; 
  delta_a_1_cent = a [i + 1] - a [i - 1];
  
  /* cure den */
  if( delta_a_1_cent == 0.0 ) delta_a_1_cent = 1e-10;
  
  criterion1 = delta_a_2_plus * delta_a_2_mins;
  criterion2 = fabs ( delta_a_1_cent ) - epsilon1 * DMIN (fabs ( a [i + 1] ), fabs ( a [i - 1] ) );
  
  eta_tilde = 0.0;
  
  if ( criterion1 < 0.0 && criterion2 > 0.0 )
    eta_tilde = - ( a [i + 2] - two*a [i + 1] + two * a[i - 1] - a[i - 2] )/( six * delta_a_1_cent );
  
  tmp = DMIN ( eta1 * ( eta_tilde - eta2 ), one );
  eta = DMAX ( 0.0, tmp );
  amins [i] = amins [i] * ( one - eta ) + eta * ( a[i-1] + half * delta_a [i-1] );
  aplus [i] = aplus [i] * ( one - eta ) + eta * ( a[i+1] - half * delta_a [i+1] );

} /* end PPM_cdsteep1D */





/* Flattening.
   See (70-71) of C&W and Appendix A of M&M 
   Here is the pointwise correction to the variables */

void PPM_flat1D ( const double *a, double *amins, double *aplus, 
		  const int i, 
		  const double omega_flat, const double om_omega_flat ) 
{
    
    amins [i] = a [i] * omega_flat + amins [i] * om_omega_flat;
    aplus [i] = a [i] * omega_flat + aplus [i] * om_omega_flat;
    
} /* end PPM_flat1D */





/* Monotonization.
   See (1.10) of C&W 
   Here is the pointwise correction to the variables */

void PPM_monoton1D ( const double *a, double *amins, double *aplus, const int i )
{
   const double CVL = 1.25;

   double D2m, D2c, D2p, S2, D2lim, delta_a_lim;

   double deltamm =  a[i-1] - a[i-2];
   double deltam  =  a[i]   - a[i-1];
   double deltac  = (a[i+1] - a[i-1])*0.5;
   double deltap  =  a[i+1] - a[i];
   double deltapp =  a[i+2] - a[i+1];

   int monotonize = 1;

   if (DMIN(deltam*deltap,deltapp*deltamm)<0) {

      D2m = a[i]   - 2*a[i-1] + a[i-2];
      D2c = a[i+1] - 2*a[i]   + a[i-1];
      D2p = a[i+2] - 2*a[i+1] + a[i];

      S2 = (D2c>0)?1.:-1.;

      D2lim = DMIN(fabs(D2c),DMIN(DMAX(S2*D2m,0),DMAX(S2*D2p,0)));

      if (S2*deltac<0) delta_a_lim = DMIN(CVL*1.5*D2lim,2*fabs(deltam));
      else             delta_a_lim = DMIN(CVL*1.5*D2lim,2*fabs(deltap));

      if (fabs(deltac)<fabs(delta_a_lim)) monotonize = 0;

   }

   if (monotonize) {

      if( ( ( aplus [i] - a [i] ) * ( amins [i] - a [i] ) ) >= 0.0 ) {
         amins [i] = a [i];
         aplus [i] = a [i];
      }

      if ( ( aplus [i] - amins [i] )*( amins [i] - 3 * a [i] + 2 * aplus [i] ) < 0.0 )
         amins [i] = 3 * a [i] - 2 * aplus [i];
            
      if ( ( aplus [i] - amins [i] )*( - aplus [i] + 3 * a [i] - 2 * amins [i] ) < 0.0 )
         aplus [i] = 3 * a [i] - 2 * amins [i];
         
   }

} /* end PPM_monoton1D */


void PPM_monoton1D_pt ( const double *a, double *amins, double *aplus, const int i )
{
   const double CVL = 1.25;

   double D2m, D2c, D2p, S2, D2lim, delta_a_lim;
   double res_plus, res_mins;

   double deltamm =  a[i-1] - a[i-2];
   double deltam  =  a[i]   - a[i-1];
   double deltac  = (a[i+1] - a[i-1])*0.5;
   double deltap  =  a[i+1] - a[i];
   double deltapp =  a[i+2] - a[i+1];

   int monotonize = 1;

   if (DMIN(deltam*deltap,deltapp*deltamm)<0) {

      D2m = a[i]   - 2*a[i-1] + a[i-2];
      D2c = a[i+1] - 2*a[i]   + a[i-1];
      D2p = a[i+2] - 2*a[i+1] + a[i];

      S2 = (D2c>0)?1.:-1.;

      D2lim = DMIN(fabs(D2c),DMIN(DMAX(S2*D2m,0),DMAX(S2*D2p,0)));

      if (S2*deltac<0) delta_a_lim = DMIN(CVL*1.5*D2lim,2*fabs(deltam));
      else             delta_a_lim = DMIN(CVL*1.5*D2lim,2*fabs(deltap));

      if (fabs(deltac)<fabs(delta_a_lim)) monotonize = 0;

   }

   if (monotonize) {

      if( ( ( *aplus - a [i] ) * ( *amins - a [i] ) ) >= 0.0 ) {
         *amins = a [i];
         *aplus = a [i];
      }

      if ( ( *aplus - *amins )*(   (*amins) - 3 * a [i] + 2 * (*aplus) ) < 0.0 )
         *amins = 3 * a [i] - 2 * (*aplus);
            
      if ( ( *aplus - *amins )*( - (*aplus) + 3 * a [i] - 2 * (*amins) ) < 0.0 )
         *aplus = 3 * a [i] - 2 * (*amins);
         
   }


} /* end PPM_monoton1D */