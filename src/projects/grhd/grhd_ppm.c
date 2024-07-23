/* Routine for 1D PPM reconstruction */

/* A note on size and indexes 

   I assume here 1D primitive vars of size: 
   ( "physical" pts ) + 2 x ( "ghost" pts ) = nx + 2 x 2 
   
   The index range is:
   
   -3  -2   -1    0    1    2    ...    nx   nx+1 nx+2 nx+3 nx+4
   x    x    x    x    o    o    ...    o    x    x    x    x 

   e.g. rho [ -3...nx+4 ] 

   plus/mins vectors require same convention for indexes
   they're filled only on these points: rhoplus [ -1...nx+1 ] */

#include "bam.h"
#include "grhd.h"

void Rec1D_PPM_grhd( 
		      double *rho, double *epsl, double *vx, double *vy, double *vz, double *pres, 
		      double *rhomins, double *epslmins, double *vxmins, double *vymins, double *vzmins,
		      double *rhoplus, double *epslplus, double *vxplus, double *vyplus, double *vzplus,
		      int nx)
{
  
  int Nghost = MATTER.HRSC_NGHOST;
  int Ntot = nx+2*Nghost;
  int off  = Nghost-1;

  /* contains all informations */
  double **buffer = (double**) malloc (6*sizeof(double*));
  for (int l=0; l<6; l++) 
    buffer[l] = (double*) malloc ((Ntot)*sizeof(double));

  double *delta_rho  = buffer[0] + off;
  double *delta_epsl = buffer[1] + off;
  double *delta_vx   = buffer[2] + off;
  double *delta_vy   = buffer[3] + off;
  double *delta_vz   = buffer[4] + off;
  double *omega      = buffer[5] + off;


  /* VARS FOR CD STEEPENING */
  double criterion_cdsteep_lhs, criterion_cdsteep_rhs;
  
  /* VARS FOR FLATTENING */
  double criterion_flat, q_flat, deltap_2, delta_omega, omega_flat, om_omega_flat;
  //static double *deltap, *deltav;
  double deltap, deltav;

  int i,l;
  double gamma_eos,cs2,h;
  double tmp;

  double tmp_aux, tmp_p, tmp_kappa, tmp_chi; // tmp for EOS call
  
  
  
  /* PARAMETERS */
  const double K0       = 1.0;   
  const double eta1     = 5.0;
  const double eta2     = 0.05;
  const double epsilon1 = 0.1;
  const double omega1   = 0.52; 
  const double omega2   = 10.0;
  const double epsilon2 = 0.5; 

  
  
  /* INTERPOLATION step 
     VanLeer limited slopes */
  for( i = -1 ; i <= nx+2 ; i++ ) {

      // NOTE These slopes are used in
      // i.)  interpolation step
      // ii.) in the CD-steepening step (only delta_rho for now)
      
      PPM_MCslope ( rho,  delta_rho,  i );
      PPM_MCslope ( vx,   delta_vx,   i );
      PPM_MCslope ( vy,   delta_vy,   i );
      PPM_MCslope ( vz,   delta_vz,   i );
      PPM_MCslope ( epsl, delta_epsl, i );
      
  }
  
  
  
  /* INTERPOLATION step  
     compute values at interfaces a_i+1/2 
     See (1.6) C&W, p.177 */
  for( i = -1 ; i <= nx+2 ; i++ ) {      

      // NOTE The interpolation reconstructs interfaces (not cells), 
      // since we need also the values
      //               amins[0]   aplus[nx+1]   
      // for the monotonization step
      // the usual loop on cells
      //         for( i = 0 ; i <= nx+1 ; i++ ) { 
      // has to be extended 
      //         for( i = -1 ; i <= nx+2 ; i++ ) { 
      
      PPM_interp1D ( rho,  rhomins,  rhoplus,  i, delta_rho  );
      PPM_interp1D ( vx,   vxmins,   vxplus,   i, delta_vx   );
      PPM_interp1D ( vy,   vymins,   vyplus,   i, delta_vy   );
      PPM_interp1D ( vz,   vzmins,   vzplus,   i, delta_vz   );
      PPM_interp1D ( epsl, epslmins, epslplus, i, delta_epsl );
      
  }
  

  
  
  /* CD-STEEPENING step */
  for( i = 0 ; i <= nx+1 ; i++ ) {
      
      EOS.comp("re","","","pc","","", rho[i],epsl[i], &tmp_aux,&cs2);

      h         = 1.+epsl[i]+ pres[i]/rho[i];
      gamma_eos = cs2*rho[i]*h/pres[i];
      
      criterion_cdsteep_lhs = 
	  gamma_eos * K0 * 
	  fabs ( rho [i + 1] - rho [i - 1] ) / 
	  ( DMIN (rho [i + 1], rho [i - 1] ) );
      
      criterion_cdsteep_rhs = fabs ( pres [i + 1] - pres [i - 1] ) / 
	  ( DMIN ( pres [i + 1], pres [i - 1] ) );
      
      if( criterion_cdsteep_lhs >= criterion_cdsteep_rhs ) {
	  	  
	  PPM_cdsteep1D ( rho, delta_rho, rhomins, rhoplus, i, epsilon1, eta1, eta2 ); 

	  // FIXME (Thu Nov 12 13:08:08 CET 2009): in principle this correction si applied only on rho var, so 
	  // (i)  PPM_cdsteep1D could be hardcoded here for rho 
	  // (ii) PPM_interpolate1D could be changed/eliminate/different output for MC-slopes vectors
	  
	  /*
	    PPM_cdsteep1D ( epsl, delta_epsl, epslmins, epslplus, i, epsilon1, eta1, eta2 );
	    PPM_cdsteep1D ( vx, delta_vx, vxmins, vxplus, i, epsilon1, eta1, eta2 );
	    PPM_cdsteep1D ( vy, delta_vy, vymins, vyplus, i, epsilon1, eta1, eta2 );
	    PPM_cdsteep1D ( vz, delta_vz, vzmins, vzplus, i, epsilon1, eta1, eta2 ); 
	  */

        }
      
  }
  
  
  
  /* FLATTENING step
     compute omega */
  for( i = -1 ; i <= nx+2 ; i++ ) {
      
      deltap = pres [i + 1] - pres [i - 1]; 
      deltav = vx [i + 1] - vx [i - 1];
      
      criterion_flat = epsilon2 * DMIN ( pres [i + 1], pres [i - 1] ) - fabs( deltap );
      
      if ((criterion_flat < 0.0) && (deltav < 0.0) ) q_flat = 1.0;
      else q_flat = 0.0;
    
      deltap_2 = pres [i + 2] - pres [i - 2];
      
      /* cure den */
      if ( deltap_2 == 0.0 ) {
	  delta_omega = - omega1;
	  if( deltap == 0.0 ) delta_omega = 1.0 - omega1;
      } else {
	  delta_omega = deltap / deltap_2 - omega1;
      }
      
      tmp = DMAX ( 0.0, delta_omega * omega2 );
      omega [i] = DMIN (1.0, q_flat * tmp );
  }

  
   /* FLATTENING step 
      correct values */
  for( i = 0 ; i <= nx+1 ; i++ ) {
      
      deltap = pres [i + 1] - pres [i - 1]; 
      deltav = vx [i + 1] - vx [i - 1];

      /* check for shock's direction  */
      if( deltap < 0.0 )    omega_flat = DMAX ( omega [i], omega [i + 1] );
      else                  omega_flat = DMAX ( omega [i], omega [i - 1] );
      
      om_omega_flat = 1.0 - omega_flat;

      PPM_flat1D ( rho,  rhomins,  rhoplus,  i, omega_flat, om_omega_flat );
      PPM_flat1D ( vx,   vxmins,   vxplus,   i, omega_flat, om_omega_flat );
      PPM_flat1D ( vy,   vymins,   vyplus,   i, omega_flat, om_omega_flat );
      PPM_flat1D ( vz,   vzmins,   vzplus,   i, omega_flat, om_omega_flat );
      PPM_flat1D ( epsl, epslmins, epslplus, i, omega_flat, om_omega_flat );
      
  }

  
  //   /* MONOTONIZATION step */
  for( i = 0 ; i <= nx+1 ; i++ ) {

      PPM_monoton1D ( rho,  rhomins,  rhoplus,  i );
      PPM_monoton1D ( vx,   vxmins,   vxplus,   i );
      PPM_monoton1D ( vy,   vymins,   vyplus,   i );
      PPM_monoton1D ( vz,   vzmins,   vzplus,   i );
      PPM_monoton1D ( epsl, epslmins, epslplus, i );
  
  }

  for (l=0; l<6; l++) 
    free(buffer[l]);
  free(buffer);


}