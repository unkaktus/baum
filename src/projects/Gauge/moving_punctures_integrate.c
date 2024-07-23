/* moving_punctures_integrate.c */
/* Bernd Bruegmann, 12/2005 */
/* Jose Gonzalez, 6/2006 */
/* Wolfgang Tichy, 10/2006 */
/* Bernd Bruegmann, 11/2006: 
  clean up logic, allow punctures on different levels */
/* mth 2012 simplified */
/* dtim 2013 add switch for trackmethods */

/* comments on integration:
   forward Euler seems to work and is quite accurate
   leapfrog goes unstable rather quickly with wild grid oscillations
   use ICN as often done for evolution, although we may want rk4 etc
   RK4 probably pointless because we don't have 4-th order shift

   high order time integration is complicated because we don't have
   the shift at intermediate times
   ICN with interpolate shift should be second order, but what about RK4?
   evolve puncture as scalar participant of standard RK evolutions?
*/ 

#include "bam.h"
#include "Gauge.h"




/* integrate moving_puncture */
int track_moving_puncture_integrate(tL* level, int np, double *pold, double *pnew)
{
  tG *grid = level->grid;
  double x0, y0, z0, x1, y1, z1;
  double betax0, betay0, betaz0, betax1, betay1, betaz1;
  double dt = level->dt;
  int i, ipick, j, n;
  int pr = 0;

  for (i=0; i<3; i++)
    pnew[i] = pold[i];
  
  /* we integrate only in the finest level, or in the level amr_lmax2 to avoid problems*/
     if (grid->lmaxpunc[np] != level->l) return 0;

   /* nothing to do for iteration 0 */
     if (level->iteration == 0) return 0;
  
  
  /* previous coordinates */
  x1 = pold[0];
  y1 = pold[1];
  z1 = pold[2];

  if (grid->half[0] && dless(x1,0) ) return 0;
  if (grid->half[1] && dless(y1,0) ) return 0;
  if (grid->half[2] && dless(z1,0) ) return 0;

  /* get shift from previous time step
     we run in ANALYZE, and the shift was adjusted in POST_OUTPUT 
     after ANALYZE of the previous time step
  */
  betax1 = interpolate_xyz_scalar(level, x1, y1, z1, Ind("betax_p"), Geti("order_RP"),LAGRANGE);
  betay1 = interpolate_xyz_scalar(level, x1, y1, z1, Ind("betay_p"), Geti("order_RP"),LAGRANGE);
  betaz1 = interpolate_xyz_scalar(level, x1, y1, z1, Ind("betaz_p"), Geti("order_RP"),LAGRANGE);

  /* 3 step ICN using averaged right hand sides
     n steps would be cheap, should be tried
  */
  x0 = x1;
  y0 = y1;
  z0 = z1;
  for (j = 0; j < 3; j++) {
    betax0 = interpolate_xyz_scalar(level, x0, y0, z0, Ind("betax"), Geti("order_RP"),LAGRANGE);
    betay0 = interpolate_xyz_scalar(level, x0, y0, z0, Ind("betay"), Geti("order_RP"),LAGRANGE);
    betaz0 = interpolate_xyz_scalar(level, x0, y0, z0, Ind("betaz"), Geti("order_RP"),LAGRANGE);
    x0 = x1 - dt * 0.5 * (betax0 + betax1);
    y0 = y1 - dt * 0.5 * (betay0 + betay1);
    z0 = z1 - dt * 0.5 * (betaz0 + betaz1);
  }

  /* save the result */
  pnew[0] = x0;
  pnew[1] = y0;
  pnew[2] = z0;
  if ( (pintoxyplane == 2) || (pintoxyplane== 3))  
                pnew[2] = pold[2];
  if (0) {
    printf("track i=%d l=%d\n", np, level->l);
    printf("  %.9f ->  %.9f\n", pold[0],pnew[0]);
    printf("  %.9f ->  %.9f\n", pold[1],pnew[1]);
    printf("  %.9f ->  %.9f\n", pold[2],pnew[2]);
  }


   if (Getv("grid", "quadrant")) {
     if (dless(pnew[1],-(level->dy))){
        if(np==0)  pnew[0] = -1.0*Getd("moving_puncture2_x");  
        else       pnew[0] = -1.0*Getd("moving_puncture1_x"); 

        if(np==0)  pnew[1] = -1.0*Getd("moving_puncture2_y");  
        else       pnew[1] = -1.0*Getd("moving_puncture1_y"); 
     }
    }

  return 1;
}




/* find maximum by fitting 1D polynomial */
double findextremum(double vm, double v0, double vp, int max)
{
    /*
  2nd order function to track extremum
  f = ax^2+bx+c   || x = [-1, 0, 1] , f = [vm, v0, vp]
  a = 0.5*(vp+vm)-v0;
  b = 0.5*(vp-vm);
  c = v0;
  f' == 0         =>  x = -b/(2.*a);
  delx = x*dx
    */
    
  double v;
    
  if (!finite(vm) || !finite(v0) || !finite(vp)) 
    return 0.;
    
  if (fabs(0.5*(vp+vm)-v0) <= 1e-10) {
    if      ((vp>v0) && (vp>vm)) v = 1.;
    else if ((vm>v0) && (vm>vp)) v =-1.;
    else                         v = 0.;
    v *= (max)?(1.):(-1.);
  } else {
    v = (vm-vp)/(4.*( 0.5*(vp+vm)-v0 ));
    v = (v> 1.)? 1.:v;
    v = (v<-1.)?-1.:v;
  }
    
  return v;
}

/* find maximum inside 3x3x3 level
   -> this is a verry bad solution, 
   but enough for this purpose
*/
int track_moving_puncture_extremum(tL *level, int np, double *pold, double *pnew )
{
  double v0,vm,vp, v;
  double dx,dy,dz;
  int var = Ind(Gets("moving_puncture_track_var"));
  int max = Getv("moving_puncture_track_mode","max");
  double minimal_move = Getd("moving_puncture_track_minmove")*level->dx;
  double x0, y0, z0, x1, y1, z1;
  int pr = 0;

  /* previous coordinates */
  x0 = pold[0];
  y0 = pold[1];
  z0 = pold[2];

  // value at the old puncture
  v0 = interpolate_xyz_scalar(level, x0,y0,z0, var, Geti("order_RP"),LAGRANGE);
  
  // look for rho around the 0 position and make a BAD new guess
  // x direction
  vm = interpolate_xyz_scalar(level, x0-level->dx,y0,z0, var, Geti("order_RP"),LAGRANGE);
  vp = interpolate_xyz_scalar(level, x0+level->dx,y0,z0, var, Geti("order_RP"),LAGRANGE);
  if (pr) printf("   %2.2e   %2.2e   %2.2e\n",vm,v0,vp);
  dx = findextremum(vm,v0,vp, max) * level->dx;
  
  // y direction
  vm = interpolate_xyz_scalar(level, x0,y0-level->dy,z0, var, Geti("order_RP"),LAGRANGE);
  vp = interpolate_xyz_scalar(level, x0,y0+level->dy,z0, var, Geti("order_RP"),LAGRANGE);
  if (pr) printf("   %2.2e   %2.2e   %2.2e\n",vm,v0,vp);
  dy = findextremum(vm,v0,vp, max) * level->dy;
  
  // z direction (like pthe puncture case)
  vm = interpolate_xyz_scalar(level, x0,y0,z0-level->dz, var, Geti("order_RP"),LAGRANGE);
  vp = interpolate_xyz_scalar(level, x0,y0,z0+level->dz, var, Geti("order_RP"),LAGRANGE);
  if (pr) printf("   %2.2e   %2.2e   %2.2e\n",vm,v0,vp);
  dz = findextremum(vm,v0,vp, max) * level->dz;

  if (pr) printf("  matter move:  %e %e %e\n",dx,dy,dz);
  
  v = interpolate_xyz_scalar(level, x0+dx,y0+dy,z0+dz, var, Geti("order_RP"),LAGRANGE);
  if (!finite(v) || !finite(v0) || 
       ((max) && (v<=v0)) ||
       ((!max) && (v>=v0)) ||
       (sqrt(dx*dx + dy*dy + dz*dz) < minimal_move))
    dx = dy = dz = 0.;
  
  x1 = x0 + dx;
  y1 = y0 + dy;
  z1 = z0 + dz;
  if ((pintoxyplane == 1) || (pintoxyplane == 3)) z1 = z0 ;
 
  /* save the result */
  pnew[0] = x1;
  pnew[1] = y1;
  pnew[2] = z1;
  
  if (0) {
    printf("track i=%d l=%d\n", np, level->l);
    printf("  %.9f ->  %.9f\n", pold[0],pnew[0]);
    printf("  %.9f ->  %.9f\n", pold[1],pnew[1]);
    printf("  %.9f ->  %.9f\n", pold[2],pnew[2]);
  }

  if (Getv("grid", "quadrant")) { 
    if (dless(pnew[1],-(level->dy))){
       if(np==0)  pnew[0] = -1.0*Getd("moving_puncture2_x");  
       else       pnew[0] = -1.0*Getd("moving_puncture1_x"); 

       if(np==0)  pnew[1] = -1.0*Getd("moving_puncture2_y");  
       else       pnew[1] = -1.0*Getd("moving_puncture1_y"); 

    }
  }

  return 1;
}





void compute_moving_puncture_distance(tL *level, double *p1, double *p2, double *length, double *clength)
{  
  double x0 = p1[0];
  double y0 = p1[1];
  double z0 = p1[2];
  double x1 = p2[0];
  double y1 = p2[1];
  double z1 = p2[2];
  double X1[4]={0,x0,y0,z0};
  double X2[4]={0,x1,y1,z1};
  
  
  if (Getv("physics","AHmod")) {
    if (Getv("compute_moving_puncture_distance","AHsimple")) {
      
      double l  = sqrt( (x1-x0)*(x1-x0) + (y1-y0)*(y1-y0) + (z1-z0)*(z1-z0));
      double r1 = Getd("ahf1_r");
      double r2 = Getd("ahf2_r");
      
      if (r1>=0.) {
        X1[1] = x0 + r1/l*(x1-x0);
        X1[2] = y0 + r1/l*(y1-y0);
        X1[3] = z0 + r1/l*(z1-z0);
      }
      if (r2>=0.) {
        X2[1] = x1 - r2/l*(x1-x0);
        X2[2] = y1 - r2/l*(y1-y0);
        X2[3] = z1 - r2/l*(z1-z0);
      }
            
     /* if (Getv("AHmod_verbose","yes")) {
        printf("compute minimized distance:\n");
        printf("  puncture1:  %2.4e %2.4e %2.4e   -> start at   %2.4e %2.4e %2.4e\n",x0,y0,z0,X1[1],X1[2],X1[3]);
        printf("  puncture2:  %2.4e %2.4e %2.4e   -> start at   %2.4e %2.4e %2.4e\n",x1,y1,z1,X2[1],X2[2],X2[3]);
      }*/
    } 
    else if (Getv("compute_moving_puncture_distance_method","AHsimpleNSBH")) {
      
      double l  = sqrt( (x1-x0)*(x1-x0) + (y1-y0)*(y1-y0) + (z1-z0)*(z1-z0));
      double r2 = Getd("ahf0_r");
      
      if (r2>=0.) {
        X2[1] = x1 - r2/l*(x1-x0);
        X2[2] = y1 - r2/l*(y1-y0);
        X2[3] = z1 - r2/l*(z1-z0);
      }
            
      if (Getv("AHmod_verbose","yes")) {
        printf("compute minimized distance:\n");
        printf("  puncture1:  %2.4e %2.4e %2.4e   -> start at   %2.4e %2.4e %2.4e\n",x0,y0,z0,X1[1],X1[2],X1[3]);
        printf("  puncture2:  %2.4e %2.4e %2.4e   -> start at   %2.4e %2.4e %2.4e\n",x1,y1,z1,X2[1],X2[2],X2[3]);
      }
    } else  if (Getv("compute_moving_puncture_distance","AH") && Getv("AHmod_output_lm","yes")) {

    }
  }
  
  
  // length is proper-distance and clength is coord-distance
  *length  = compute_path_integral(level->grid, X1,X2, Gets("compute_moving_puncture_distance"));
  *clength = sqrt( (X2[1]-X1[1])*(X2[1]-X1[1]) + 
                   (X2[2]-X1[2])*(X2[2]-X1[2]) + 
                   (X2[3]-X1[3])*(X2[3]-X1[3]) );


  double morigin = Getd("moving_puncture_mir_orgin");
  if (morigin > 0 ) {
    if (*length < morigin){ 
      trackmethod = 3;
     }
  }

  double forigin = Getd("moving_puncture_finboxfix");
  if (forigin > 0 ) {
    if (*length < forigin){ 
      trackmethod = 4;
     }
  }

}




















