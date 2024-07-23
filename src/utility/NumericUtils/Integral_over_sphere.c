/* integral_over_sphere.c */
/* Jose Gonzalez 9/04 */
/* mth 08/11 */

#include "bam.h"






double integral_over_sphere_old(tL *level, int ncircles, int npoints,
				double r, int index)
{
      
  /* ncircles = number of circles for the sphere (without the poles) -> must be ODD  */
  /* npoints  = number of points for each circle of the sphere       -> must be EVEN */
  /* index    = index of the variable that we want to integrate                      */
  
  int i,j,l=0;
  double dtheta,dphi;
  double local_sum=0,integral=0;
  int np_proc,np_rem;
  /* WT: double coeff_theta[ncircles+1],coeff_phi[npoints+1]; should really be: */
  double *coeff_theta = dmalloc(ncircles+2);
  double *coeff_phi   = dmalloc(npoints+2);
  int np_total=(ncircles+1)*(npoints+1);
  int init_point,end_point;
  double x,y,z,theta,phi;
  /* WT: double matrix[ncircles+1][npoints+1]; should really be: */
  double **matrix = dmatrix(0, ncircles+1, 0, npoints+1);
  
  if (ncircles%2==0){
    errorexit("Call to integral over sphere with even number of circles");
  }
  if ((npoints+1)%2==0){
    errorexit("Call to integral over sphere with odd number of points");
  }

  if (!(Getv("grid", "box")))
      errorexit("Integral over sphere only tested in Box grid");

  if (Getv("grid", "bitant")) {
      /* z=0 reflection symmetry */
      dtheta = PI/(2*(ncircles+1));
      dphi = 2.0*PI/npoints;

  } else if (Getv("grid", "rotant")) {
      /* x=y=0 inversion symmetry */
      dtheta = PI/(ncircles+1);
      dphi = PI/npoints;

  } else if (Getv("grid", "quadrant")) {
      /* z=0 reflection and x=y=0 inversion symmetry */
      dtheta = PI/(2*(ncircles+1));
      dphi = PI/npoints;

  } else if (Getv("grid", "octant")) {
      /* x=0, y=0, and z=0 reflection symmetry */
      dtheta = PI/(2*(ncircles+1));
      dphi = PI/(2*npoints);
  } else if (Getv("grid", "qreflect")) {
      /* z=0 and y=0 reflection symmetry, without rotation */
      dtheta = PI/(2*(ncircles+1));
      dphi = PI/npoints;
  } else {
      /* Full grid */
      dtheta = PI/(ncircles+1);
      dphi = 2.0*PI/npoints;
     
  }
 
  /* Theta coefficients */
  for(i=1; i<=ncircles-2; i=i+2){
      coeff_theta[i] = 2;
      coeff_theta[i+1] = 1;
  }
  coeff_theta[ncircles] = 2;
  coeff_theta[ncircles+1] = 0.5;

  /* Phi coefficients */
  coeff_phi[1] = 0.5;
  for(i=2; i<=npoints; i=i+2){
      coeff_phi[i] = 2;
      coeff_phi[i+1] = 1;
  }
  coeff_phi[npoints+1] = 0.5;

  np_proc = np_total/bampi_size();
  np_rem = np_total%bampi_size();
  if(bampi_rank()<np_rem) np_proc = np_proc + 1;

  if(bampi_rank()<np_rem){
      init_point = bampi_rank()*np_proc + 1;
  }else{
      init_point = bampi_rank()*np_proc+np_rem+1;
  }

  end_point = init_point + np_proc - 1;

  for(i=1; i<=ncircles+1; i++){
      for(j=1; j<=npoints+1; j++){
	  
	  theta = i*dtheta;
	  phi = (j-1)*dphi;

	  x = r*sin(theta)*cos(phi);
	  y = r*sin(theta)*sin(phi);
	  z = r*cos(theta);

          errorexit("fixme: outdated call to interpolate");
          // matrix[i][j] = interpolate_xyz_scalar_box(level, x, y, z, index, 0)*sin(theta);
	  l = l+1;

	  if(!((l>=init_point)&&(l<=end_point))) matrix[i][j] = 0;
	 
	  local_sum = local_sum + 4*dtheta*dphi/9*coeff_theta[i]*matrix[i][j]*coeff_phi[j];

      }
  }

  bampi_allreduce_sum_vector(&local_sum, &integral, 1);

  if (Getv("grid", "box")) {
      if (Getv("grid", "bitant")){
	  /* z=0 reflection symmetry */
	  integral = 2*integral;
	  //printf("bitant\n");
      } else if (Getv("grid", "rotant")){
	  /* x=y=0 inversion symmetry */
	  integral = 2*integral;
	  //printf("rotant\n");
      } else if (Getv("grid", "quadrant")){
	  /* z=0 reflection and x=y=0 inversion symmetry */
	  integral = 4*integral;
	  //printf("quadrant\n");
      } else if (Getv("grid", "octant")){
	  /* x=0, y=0, and z=0 reflection symmetry */
	  integral = 8*integral;
	  //printf("octant\n");
      } else if (Getv("grid", "qreflect")) {
          /* z=0 reflection and x=y=0 inversion symmetry without rotation*/
	  integral = 4*integral;
      }  else {
	  /* Full grid */
	  integral = integral;
	  //printf("full\n");
      }
  }

  /* free vars and the dmatrix */
  free(coeff_theta);
  free(coeff_phi);
  free_dmatrix(matrix, 0, ncircles+1, 0, npoints+1);
  
/*
  if (0) {
    x = integral_over_sphere_new(level, ncircles, npoints, r, index);

    printf("integral new     = %.16e\n", x);
    printf("integral old     = %.16e\n", integral);
    printf("integral new-old = %.16e\n", x - integral);
  }
*/

  //printf("integral = %f\n",integral);
  return integral;
}




/* more efficient version doing sums locally

   an issue is that there will be points which by round-off error
   in the coordinates are attributed to more than one processor, i.e. we
   have to communicate who does what; this is unfortunate because this
   involves communication for each single point, while for integration
   in principle only the sum over all points per processor has to be
   communicated

   medium fast version:
   this processor does all the interpolation it can and saves the result
   then all processors communicate who did what, and only one of them
   is allowed to do the sum for the integral

   faster:
   - cache who did what for the next call, but this has to take into account
     that e.g. the boxes may move and that this routine will be called for
     various radii
   - work on lists of variables because the largest overhead is per point

   may want to extend this to variable order of interpolation

   BB 3/2007
*/
double integral_over_sphere_Simpson(tL *level, double x0, double y0, double z0,
                                    int ntheta, int nphi, double r, int index, int order)
{
  /* ntheta = number of theta values (without the poles) -> must be ODD  */
  /* nphi   = number of phi values                       -> must be EVEN */
  /* index  = index of the variable that we want to integrate            */
  
  int rank = bampi_rank();
  int i, j, n, l;
  double dtheta, dphi;
  double localsum = 0;
  double integral = 0;
  int flag;
  double value;
  double *coeff_theta = dmalloc(ntheta+2);
  double *coeff_phi   = dmalloc(nphi+2);
  int npoints = (ntheta+1)*(nphi+1);
  double *result = dmalloc(npoints);
  double *localp = dmalloc(npoints);
  double *globalp = dmalloc(npoints);
  double x, y, z, rp,pp,tp;
  double theta, phi;
  double c, cost, sint;
  double symfactor = 1;
  
  /* checks */
  if (ntheta%2 == 0)
    errorexit("Call to integral over sphere with even number of circles");

  if ((nphi+1)%2 == 0)
    errorexit("Call to integral over sphere with odd number of points");

  if (!(Getv("grid", "box")))
      errorexit("Integral over sphere only tested in Box grid");

  /* handle symmetries */
  if (Getv("grid", "bitant")) {
      /* z=0 reflection symmetry */
      dtheta = PI/(2*(ntheta+1));
      dphi = 2.0*PI/nphi;
      symfactor = 2;
  } 
  else if (Getv("grid", "rotant")) {
      /* x=y=0 inversion symmetry */
      dtheta = PI/(ntheta+1);
      dphi = PI/nphi;
      symfactor = 2;
  } 
  else if (Getv("grid", "quadrant")) {
      /* z=0 reflection and x=y=0 inversion symmetry */
      dtheta = PI/(2*(ntheta+1));
      dphi = PI/nphi;
      symfactor = 4;
  } 
  else if (Getv("grid", "octant")) {
      /* x=0, y=0, and z=0 reflection symmetry */
      dtheta = PI/(2*(ntheta+1));
      dphi = PI/(2*nphi);
      symfactor = 8;
  } 
  else if (Getv("grid", "qreflect")) {
      /* z=0 and y=0 reflection symmetry, without rotation */
      dtheta = PI/(2*(ntheta+1));
      dphi = PI/nphi;
      symfactor = 4;
  } 
  else {
    /* full grid */
    dtheta = PI/(ntheta+1);
    dphi = 2.0*PI/nphi;
    symfactor = 1;
  }
 
  /* Theta coefficients */
  for (i=1; i<=ntheta-2; i=i+2) {
    coeff_theta[i] = 2;
    coeff_theta[i+1] = 1;
  }
  coeff_theta[ntheta] = 2;
  coeff_theta[ntheta+1] = 0.5;

  /* Phi coefficients */
  coeff_phi[1] = 0.5;
  for (i=2; i<=nphi; i=i+2) {
    coeff_phi[i] = 2;
    coeff_phi[i+1] = 1;
  }
  coeff_phi[nphi+1] = 0.5;


  /* new parallelization */
  
  /* interpolate locally for all points for which this is possible */
  n = 0;
  for (i = 1; i <= ntheta+1; i++) {
    theta = i*dtheta;
    sint = sin(theta);
    cost = cos(theta);
    for (j = 1; j <= nphi+1; j++) {
      phi = (j-1)*dphi;
      c = 4*dtheta*dphi/9 * coeff_theta[i]*coeff_phi[j] * sint;
      x = x0 + r*sint*cos(phi);
      y = y0 + r*sint*sin(phi);
      z = z0 + r*cost;

      if (level->shells) {
        l = find_shellsbox_from_xyz(x,y,z);
        
        convert_box_to_shells(l, x,y,z, &rp,&pp,&tp);
        
        if (pp>=-PI/4. || pp<=PI/4. || tp>=-PI/4. || tp<=PI/4.) {
          flag = interpolate_xyz_localinbox_minimal(level->box[l], 
            rp,pp,tp, 1,&index,&value, order,LAGRANGE);
        } else {
          flag = 0;
        }
        
      } else {
        /* interpolate for local points */
        flag = interpolate_xyz_local_minimal(level, x, y, z, 1, 
					   &index, &value, order,LAGRANGE);
      }

      /* save result */
      if (flag) {
	result[n] = c * value;
	localp[n] = rank;
      } else 
	localp[n] = -1;
      
      n++;
    }
  }

  /* determine the processor with the highest rank for each point */
  bampi_allreduce_max_vector(localp, globalp, n); 

  /* sum all the results that this processor is responsible for
     and check at the same time whether any points where left over
  */
  localsum = 0;
  for (i = 0; i < n; i++) {
    if (globalp[i] == rank)
      localsum += result[i];
    else if (globalp[i] == -1)
      return 0;
  }
      
  /* the integral is the sum over all local sums */
  bampi_allreduce_sum_vector(&localsum, &integral, 1);

  /* multiply result by appropriate symmetry factor to get
     result for integration over entire sphere
  */
  integral *= symfactor;

  /* free vars and the dmatrix */
  free(coeff_theta);
  free(coeff_phi);
  free(result);
  free(localp);
  free(globalp);
  
  return integral;
}



/* more accurate version */
double integral_over_sphere_Simpson_Gauss_Chebyshev(tL *level, double x0,double y0,double z0,
                                int ntheta, int nphi, double r, int index)
{
  /* ntheta = number of theta values (without the poles) -> must be ODD  */
  /* nphi   = number of phi values                       -> must be EVEN */
  /* index  = index of the variable that we want to integrate            */
  
  int rank = bampi_rank();
  int order = Geti("order_RP");
  int i, j, n, l;
  double dtheta, dphi;
  double localsum = 0;
  double integral = 0;
  int flag;
  double value;
  double *coeff_phi   = dmalloc(nphi+2);
  int npoints = (ntheta+1)*(nphi+1);
  double *result = dmalloc(npoints);
  double *localp = dmalloc(npoints);
  double *globalp = dmalloc(npoints);
  double x, y, z, rp,pp,tp;
  double theta, phi;
  double w,g,Rn,xi;
  double symfactor = 1;
  
  timer_start(0, "integral_over_sphere_new");

  /* checks */
  if (ntheta%2 == 0)
    errorexit("Call to integral over sphere with even number of circles");

  if ((nphi+1)%2 == 0)
    errorexit("Call to integral over sphere with odd number of points");

  if (!(Getv("grid", "box")))
    errorexit("Integral over sphere only tested in Box grid");

  /* handle symmetries */
  if ((Getv("grid", "bitant")) || (Getv("grid", "rotant")) || (Getv("grid", "quadrant")) || 
        (Getv("grid", "octant")) || (Getv("grid", "qreflect")))
    errorexit("is not implemented yet");
  
  
  
  /* full grid */
  dphi = 2.0*PI/nphi;

  /* Phi coefficients for Simpson rule */
  coeff_phi[0] = 0.5;
  for (i=1; i<=nphi; i=i+2) {
    coeff_phi[i] = 2;
    coeff_phi[i+1] = 1;
  }
  coeff_phi[nphi] = 0.5;

  /* we need different interpolation when using shells */
  int shells = (level->l==level->grid->lmin && Getv("grid","shells"));

  
  /* now go through points, first phi */
  n = 0;
  for (i=0; i<=nphi; i++) {
    phi = i*dphi;
    
    for (j=1; j<=ntheta; j++) {
      
      xi    = cos( (2.*j-1.)*PI/(2.*ntheta) );
      theta = PI/2.*xi + PI/2.;
      
      /* find cartesian coordinates */
      x = x0 + r*sin(theta)*cos(phi);
      y = y0 + r*sin(theta)*sin(phi);
      z = z0 + r*cos(theta);
      
      /* interpolate the value */
      if (shells) {
        /* coordinated transformation -> this seems stupid/unefficiant I know */
        l = find_shellsbox_from_xyz(x,y,z);
        convert_box_to_shells(l, x,y,z, &rp,&pp,&tp);
        /* interpolate only if we are not at the shells ghostzones */
        if (pp>=-PI/4. || pp<=PI/4. || tp>=-PI/4. || tp<=PI/4.) {
          flag = interpolate_xyz_localinbox_minimal(level->box[l], rp,pp,tp, 1,&index,&value, order,LAGRANGE);
        } else {
          flag = 0;
        }
      } else {
        /* interpolate for local points */
        flag = interpolate_xyz_local_minimal(level, x,y,z, 1, &index, &value, order,LAGRANGE);
      }
      
      /* save result */
      if (flag) {
        
        w  = PI/ntheta;
        g  = PI/2. * sqrt(1.-xi*xi) * value *sin(theta);
        Rn = PI/(pow(2.,ntheta-1.)*fact(2.*ntheta)) * pow(g,2.*ntheta);
        
        result[n] = ( 2./3.*dphi*coeff_phi[i] ) * ( w*g + Rn );
        localp[n] = rank;
      } else {
        localp[n] = -1;
      }
      
      n++;
    }
  }

  /* determine the processor with the highest rank for each point */
  bampi_allreduce_max_vector(localp, globalp, n); 

  /* sum all the results that this processor is responsible for
  and check at the same time whether any points where left over
  */
  localsum = 0;
  for (i = 0; i < n; i++) {
    if (globalp[i] == rank)
      localsum += result[i];
    else if (globalp[i] == -1)
      return 0;
  }
      
  /* the integral is the sum over all local sums */
  bampi_allreduce_sum_vector(&localsum, &integral, 1);

  /* multiply result by appropriate symmetry factor to get
  result for integration over entire sphere
  */
  //integral *= symfactor;

  /* free vars and the dmatrix */
  free(coeff_phi);
  free(result);
  free(localp);
  free(globalp);
  
  return integral;
}


/* pick version */
double integral_over_sphere(tL *level, double x0,double y0,double z0,
                            int ntheta, int nphi,  double r, int index, int order)
{

  double v;
  timer_start(0, "integral_over_sphere");


  /* old version */
  if (Getv("integrate_over_sphere","old"))
    v = integral_over_sphere_old(level, ntheta, nphi, r, index);

  /* actual working version */
  else if (Getv("integrate_over_sphere","Simpson"))
    v = integral_over_sphere_Simpson(level, x0,y0,z0, ntheta, nphi, r, index, order);

  /* version with Simpson in phi direction, Gauss-Chebyshev in theta direction */
  else if (Getv("integrate_over_sphere","GaussCheb"))
    v = integral_over_sphere_Simpson_Gauss_Chebyshev(level, x0,y0,z0, ntheta,nphi, r, index);

  else
    errorexit("this integration method does not exist yet");


  timer_stop(0, "integral_over_sphere");
  if (0) printf("integral = %f\n",v);
  return v;
}


















