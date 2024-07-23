/* corotation.c */
/* Bernd Bruegmann, 10/03 , Wolfgang Tichy 12/2003 */

#include "bam.h"
#include "Gauge.h"

#define EXC_BOUNDARY 2

double *commotion_x, *commotion_y;
int commotion_n = 0, commotion_nmax = 0;

/* global vars containing lapse at bh on x-, y-, and z-axis */
double alpha_bh_x1;
double alpha_bh_y1;
double alpha_bh_z1;
double alpha_bh_x2;
double alpha_bh_y2;
double alpha_bh_z2;




/* attenuation of psi^-psipower for punctures */
double punc_att(double psi, double psipower, int exp_to_one)
{
  double p = pow(psi, -psipower);
  
  if(exp_to_one)
  {
    double a = Attenuation01(p, 20, 0.88);
    return (1.0-a) * p + a;
  }
  else
    return p;
}



/* attenuation function for radial direction 
   r0 = 4:  p=1, c=0.5 to 0.25,  p=2, c=1
   this attenuation is applied in addition to psi^-psipower
   make sure that att(r0) = 1
*/
double rdot_att(double c, double p, double r0, double x, double y, 
                int exp_to_zero)
{
  double att, r2, r20;

  if (p == 0) return 1;  // backward compatibility 

  r2  = c*c + (x*x + y*y)/(r0*r0);
  r20 = c*c + 1;

  att = pow(r2/r20, -p/2);
  att /= r0;

  if(exp_to_zero)
    return Attenuation01(att, 40, 0.035) * att;
  else
    return att;
}




/* initialize corotation shift 
   called after initial data has been computed, so we assume that the 
   shift has already been enabled and initialized

   recall that the shift has opposite sign compared to velocities 
*/
int initialize_betarot(tL *level)
{
  double *x = Ptr(level, "x");
  double *y = Ptr(level, "y");
  double *z = Ptr(level, "z");
  double *psi = Ptr(level, "psi");
  double *betax = Ptr(level, "betax");
  double *betay = Ptr(level, "betay");
  double *betaz = Ptr(level, "betaz");
  double *betarotx = PtrEnable(level, "betarotx");
  double *betaroty = PtrEnable(level, "betaroty");
  double *betarotz = PtrEnable(level, "betarotz");
  double rdot  = Getd("corotation_rdot");
  double omega = Getd("corotation_omega");
  double scale = Getd("corotation_scale");
  double psipower = Getd("corotation_psipower") * Getv("physics", "punctures");
  double p, q, r;
  int a_rational = Getv("corotation_attenuate", "rational");
  int a_exponential = Getv("corotation_attenuate", "exponential");
  double a = Getd("corotation_attenuate_a");
  double b = Getd("corotation_attenuate_b");
  double c = Getd("corotation_attenuate_c");
  double omega_c = Getd("corotation_omega_c");
  double omega_p = Getd("corotation_omega_p");
  double rdot_c = Getd("corotation_rdot_c");
  double rdot_p = Getd("corotation_rdot_p");
  int exponFarAtt = !Getv("corotation_exponFarAtt", "no");
  double x1 = Getd("bhx1");
  double y1 = Getd("bhy1");
  double z1 = Getd("bhz1");
  double x2 = Getd("bhx2");
  double y2 = Getd("bhy2");
  double z2 = Getd("bhz2");
  int i;

  /* the omega parameter can be set from initial data */
  omega *= scale;  /* <- arbitrarily scale omaga */
  printf("Gauge: Initial corotation shift with:\n");
  printf("       omega = %e * %e  = %e\n", 
	 Getd("corotation_omega"), Getd("corotation_scale"), omega);

  /* the rdot parameter in principal, too */
  printf("       rdot = %e\n", rdot);


  /* for all points including boundary */
  forallpoints(level, i) {

    /* compute and save the corotation shift 
       eventually should allow arbitrary axis
    */
    p = rdot_att(omega_c, omega_p, y1, x[i], y[i], exponFarAtt);
    betarotx[i] = - omega * y[i] * p;
    betaroty[i] =   omega * x[i] * p;
    betarotz[i] = 0;

    /* radial correction */
    q = rdot_att(rdot_c, rdot_p, y1, x[i], y[i], exponFarAtt);
    betarotx[i] += - rdot * x[i] * q;
    betaroty[i] += - rdot * y[i] * q;
    betarotz[i] += 0;

    /* for puncture data, the shift has to vanish at the punctures */
    p = 1;

    if(Getv("corotation_puncture_attenuation", "with_psi"))
    {
      p *= punc_att(psi[i], psipower, exponFarAtt);
    }
    if(Getv("corotation_puncture_attenuation", "linear"))
    {
      double D  = sqrt( pow((x2-x1),2) + pow((y2-y1),2) + pow((z2-z1),2) );
      r = sqrt(x[i]*x[i] + y[i]*y[i] + z[i]*z[i]);
      p = 1 - ( (1 + pow(r/(0.5*D), 4))/(1 + pow(r/(0.5*D), 6)) );
    }

    /* we may want to attenuate the shift away near infinity */
    if (a_rational) {
      r = sqrt(x[i]*x[i] + y[i]*y[i] + z[i]*z[i]);
      p /= 1 + b*pow(r, c);
    }

    if (a_exponential) {
      r = sqrt(x[i]*x[i] + y[i]*y[i] + z[i]*z[i]);
      
      /* ... */
    }
    

    /* adjust corotation shift */
    if (p != 1) {
      betarotx[i] *= p;
      betaroty[i] *= p;
      betarotz[i] *= p;
    }
    
    /* adjust shift that was set with initial data */
    betax[i] += betarotx[i];
    betay[i] += betaroty[i];
    betaz[i] += betarotz[i];
  }

  return 0;
}




/* add or remove the omega shift from beta */
void add_betarot(tVarList *vl, double sign)
{
  int ibetax = vl->index[18];  // temporary hack
  tL *level = vl->level;
  double *betax = level->v[ibetax+0];
  double *betay = level->v[ibetax+1];
  double *betaz = level->v[ibetax+2];
  double *x = Ptr(level, "x");
  double *y = Ptr(level, "y");
  double *z = Ptr(level, "z");
  double *betarotx = Ptr(level, "betarotx");
  double *betaroty = Ptr(level, "betaroty");
  double *betarotz = Ptr(level, "betarotz");
  double rdot  = Getd("corotation_rdot");
  double omega = Getd("corotation_omega") * Getd("corotation_scale");
  double omega_c = Getd("corotation_omega_c");
  double omega_p = Getd("corotation_omega_p");
  double rdot_c = Getd("corotation_rdot_c");
  double rdot_p = Getd("corotation_rdot_p");
  int exponFarAtt = !Getv("corotation_exponFarAtt", "no");
  double y1 = Getd("bhy1");
  double p, q;
  int i;
  
  if (!strstr(VarName(ibetax), "betax"))
    errorexit("add_betarot: needs fix for choices other than ICN/BSSN");

  if (!Getv("corotation_attenuate", "no"))
    sign = 2*sign/fabs(sign);

  if (sign == 1) {
    forallpoints(level, i) {
      p = rdot_att(omega_c, omega_p, y1, x[i], y[i], exponFarAtt);
      q = rdot_att(rdot_c, rdot_p, y1, x[i], y[i], exponFarAtt);
      betax[i] += - omega * y[i] * p - rdot * x[i] * q;
      betay[i] +=   omega * x[i] * p - rdot * y[i] * q;
    }
  } 
  
  if (sign == -1) {
    forallpoints(level, i) {
      p = rdot_att(omega_c, omega_p, y1, x[i], y[i], exponFarAtt);
      q = rdot_att(rdot_c, rdot_p, y1, x[i], y[i], exponFarAtt);
      betax[i] -= - omega * y[i] * p - rdot * x[i] * q;
      betay[i] -=   omega * x[i] * p - rdot * y[i] * q;
    }
  } 
  
  if (sign == 2) {
    forallpoints(level, i) {
      betax[i] += betarotx[i];
      betay[i] += betaroty[i];
      betaz[i] += betarotz[i];
    }
  } 
  
  if (sign == -2) {
    forallpoints(level, i) {
      betax[i] -= betarotx[i];
      betay[i] -= betaroty[i];
      betaz[i] -= betarotz[i];
    }
  } 
}




/* monitor symmetry near black hole
   experimental first version
   optimize: search for puncture only once, not five times per output
*/
int monitor_bh_asymmetry(tL *level)
{
  tG *g = level->grid;

  /* only do it on finest level, and only when time aligned with coarsest */
  if (level->l == g->lmax && 
      dequal(level->time, g->level[g->lmin]->time)) {
    int pr = 0;
    double d, fcc, fmc, fpc, fcm, fcp;
    double min[2];
    double dx = level->dx;
    double dy = level->dy;
    double dz = level->dz;
    int i;
    char filename[1000];
    
    /* variable to look at, should become list of variables in a parameter */
    char *vname = "alpha";
    int vi = Ind(vname);
    
    /* puncture coordinates, should look at both punctures if not quadrant */
    double x = Getd("bhx1");
    double y = Getd("bhy1");
    double z = Getd("bhz1");

    /* asymmetry for excision: compute difference between right and left side
       (would also work for puncture)
       assumes there is only one excision region (say, binary in quadrant)
       z-reflection works out automatically for scalar
       we also assume that the excision region is a centered sphere or cube
    */
    if (Getv("boundary", "excision")) {
      double *xp = Ptr(level, "x");
      double *yp = Ptr(level, "y");
      double *zp = Ptr(level, "z");
      double *mask = Ptr(level, "excisionmask");
      double *v = level->v[vi];

      /* look in x and y direction */
      if (0) {
	fmc = fpc = fcm = fcp = 0;
	forallpoints(level, i) 
	  if (mask[i] == EXC_BOUNDARY)
	    if (dless(fabs(zp[i]-z), dz)) {
	      if (dless(fabs(xp[i]-x), dx)) {
		if (yp[i]-y > 0) fcp += v[i];
		else             fcm += v[i];
	      }
	      if (dless(fabs(yp[i]-y), dy)) {
		if (xp[i]-x > 0) fpc += v[i];
		else             fmc += v[i];
	      }
	    }
	min[0] = -(fpc - fmc)/2;  // could normalize with /(fpc + fmc);
	min[1] = -(fcp - fcm)/2;
      }

      /* find center of mass */
      if (1) {
	double sumv = 0, sumvx = 0, sumvy = 0;

	forallpoints(level, i) 
	  if (mask[i] == EXC_BOUNDARY && dless(fabs(zp[i]-z), dz)) {
	    sumvx += v[i]*(xp[i]-x);
	    sumvy += v[i]*(yp[i]-y);
	    sumv  += v[i];
	  }
	if (dequal(sumv, 0)) sumv = 1;
	min[0] = -sumvx/sumv;
	min[1] = -sumvy/sumv;
      }
    }

    /* asymmetry for punctures: minimum of parabola through puncture */
    if (!Getv("boundary", "excision") && !Getv("Gauge", "moving_puncture")) {

      /* value at puncture */
      fcc = interpolate_xyz_scalar(level, x, y, z, vi, pr);
      
      /* tangential direction, assumed to be x direction */
      fmc = interpolate_xyz_scalar(level, x-dx, y, z, vi, pr);
      fpc = interpolate_xyz_scalar(level, x+dx, y, z, vi, pr);
      
      /* radial direction, assumed to be y direction */
      fcm = interpolate_xyz_scalar(level, x, y-dy, z, vi, pr);
      fcp = interpolate_xyz_scalar(level, x, y+dy, z, vi, pr);
      
      /* location of minimum of fitted parabola */
      d = fmc - 2*fcc + fpc;
      if (dequal(d,0)) min[0] = 0; else
	min[0] = dx/2 * (fmc - fpc)/d;
      d = fcm - 2*fcc + fcp;
      if (dequal(d,0)) min[1] = 0; else
	min[1] = dy/2 * (fcm - fcp)/d;

      /* info */
      if (pr) {
	printf("minimum of fitted parabola in x is  %e, %e\n", x+min[0], 
	       interpolate_xyz_scalar(level, x+min[0], y, z, vi, pr));
	printf("minimum of fitted parabola in y is  %e, %e\n", y+min[1],
	       interpolate_xyz_scalar(level, x, y+min[1], z, vi, pr));
      }   
    }

    /* moving puncture */
    if (Getv("Gauge", "moving_puncture")) {
      min[0] = Getd("moving_puncture_x");
      min[1] = Getd("moving_puncture_y");
      vname = "puncture";
      if (0) printf("mp: got %f %f\n", min[0], min[1]);
    }
  
    /* output */
    for (i = 0; i < 2; i++) {
      snprintf(filename, 1000, "%s_asymmetry_%s", vname, (i ? "y":"x"));
      write_scalar(level, filename, min[i]);
    }
    if (Getv("corotation_monitor", "adjust")) {
      write_parameter_double(level, "corotation_scale");
      write_parameter_double(level, "corotation_rdot");
    }

    /* save asymmetries */
    Setd("commotion_x", min[0]);
    Setd("commotion_y", min[1]);
    if (commotion_n == commotion_nmax) {
      commotion_nmax += 256;
      commotion_x = realloc(commotion_x, commotion_nmax * sizeof(double));
      commotion_y = realloc(commotion_y, commotion_nmax * sizeof(double));
      if (!commotion_x || !commotion_y) 
	errorexit("monitor_bh_asymmetry: out of memory");
    }
    commotion_x[commotion_n] = min[0];
    commotion_y[commotion_n] = min[1];
    commotion_n += 1;
    if (0 && commotion_n >= 4) {
      double dt = g->level[g->lmin]->dt;
      double x0, x1, x2, x3, y0, y1, y2, y3, ddx, ddy;
      i = commotion_n - 1;
      x0 = commotion_x[i];
      x1 = commotion_x[i-1];
      x2 = commotion_x[i-2];
      x3 = commotion_x[i-3];
      y0 = commotion_y[i];
      y1 = commotion_y[i-1];
      y2 = commotion_y[i-2];
      y3 = commotion_y[i-3];
      ddx = (2*x0 - 5*x1 + 4*x2 - x3)/(dt*dt);
      ddy = (2*y0 - 5*y1 + 4*y2 - y3)/(dt*dt);
      write_scalar(level, "commotion_ddx", ddx);
      write_scalar(level, "commotion_ddy", ddy);
      ddx = (3*x0 - 4*x1 + x2)/(2*dt);
      ddy = (3*y0 - 4*y1 + y2)/(2*dt);
      write_scalar(level, "commotion_dx3", ddx);
      write_scalar(level, "commotion_dy3", ddy);
      ddx = (11*x0 - 18*x1 + 9*x2 - 2*x3)/(6*dt);
      ddy = (11*y0 - 18*y1 + 9*y2 - 2*y3)/(6*dt);
      write_scalar(level, "commotion_dx4", ddx);
      write_scalar(level, "commotion_dy4", ddy);
    }
  }

  return 0;
}




/* adjust shift on all levels */
void adjust_beta(tG *g, double domega, double drdot)
{
  double psipower = Getd("corotation_psipower") * 
                    Getv("physics", "punctures");
  double omega_c = Getd("corotation_omega_c");
  double omega_p = Getd("corotation_omega_p");
  double rdot_c = Getd("corotation_rdot_c");
  double rdot_p = Getd("corotation_rdot_p");
  double y1 = Getd("bhy1");
  int exponFarAtt = !Getv("corotation_exponFarAtt", "no");
  int adjust_betarot = !Getv("corotation_adjust_betarot", "no");
  double o, p, q;
  int i, l;

  for (l = g->lmin; l <= g->lmax; l++) {
    tL *level = g->level[l];
    double *betax = Ptr(level, "betax");
    double *betay = Ptr(level, "betay");
    double *betarotx = Ptr(level, "betarotx");
    double *betaroty = Ptr(level, "betaroty");
    double dbetx, dbety;
    double *x = Ptr(level, "x");
    double *y = Ptr(level, "y");
    double *psi = Ptr(level, "psi");

    if (0) printf("adjusting shift on level %d\n", level->l);

    forallpoints(level, i)
    {
      o = punc_att(psi[i], psipower, exponFarAtt);
      p = rdot_att(omega_c, omega_p, y1, x[i], y[i], exponFarAtt);
      q = rdot_att(rdot_c, rdot_p, y1, x[i], y[i], exponFarAtt);
      dbetx = o * (- domega * y[i] * p - drdot * x[i] * q);
      dbety = o * (  domega * x[i] * p - drdot * y[i] * q);
      betax[i] += dbetx;
      betay[i] += dbety;
      if(adjust_betarot)
      {
        betarotx[i] += dbetx;
        betaroty[i] += dbety;
      }
    }
  } 
}




/* monitor lapse value near BH  */
int monitor_bh_lapse(tL *level)
{
  tG *g = level->grid;

  /* only do it on finest level, and only when time aligned with coarsest */
  if(level->l == g->lmax && 
     dequal(level->time, g->level[g->lmin]->time))
  {
    int pr = 1;
    double dx = level->dx;
    double dy = level->dy;
    double dz = level->dz;
    int i;
    char filename[1000];
    double avrvx, avrvy, avrvz;
    int    numvx, numvy, numvz;
    
    char *vname = "alpha";
    int vi = Ind(vname);
    
    /* puncture coordinates, should look at both punctures if not quadrant */
    double x1 = Getd("bhx1");
    double y1 = Getd("bhy1");
    double z1 = Getd("bhz1");

    /* use a little patch of size patchsize to average the lapse 
       on the excision surface along each axis */
    avrvx=0; avrvy=0; avrvz=0;
    numvx=0; numvy=0; numvz=0;
    if (Getv("boundary", "excision"))
    {
      double patchsize=2.0;
      double *xp = Ptr(level, "x");
      double *yp = Ptr(level, "y");
      double *zp = Ptr(level, "z");
      double *mask = Ptr(level, "excisionmask");
      double *v = level->v[vi];

      forallpoints(level, i) 
	if(mask[i] == EXC_BOUNDARY)
	{
	  if(dless(fabs(yp[i]-y1), patchsize*dy) && 
	     dless(fabs(zp[i]-z1), patchsize*dz))
	  {
	    avrvx += v[i];
	    numvx++;
	    // printf("x: x=%f  y=%f  z=%f\n",xp[i], yp[i], zp[i]);
	  }
	  if(dless(fabs(xp[i]-x1), patchsize*dx) &&
	     dless(fabs(zp[i]-z1), patchsize*dz))
	  {
	    avrvy += v[i];
	    numvy++;
	    // printf("y: x=%f  y=%f  z=%f\n",xp[i], yp[i], zp[i]);
	  }
	  if(dless(fabs(yp[i]-y1), patchsize*dy) &&
	     dless(fabs(xp[i]-x1), patchsize*dx))
	  {
	    avrvz += v[i];
	    numvz++;
	    // printf("z: x=%f  y=%f  z=%f\n",xp[i], yp[i], zp[i]);
	  }
	}
	avrvx /= numvx;
	avrvy /= numvy;
	avrvz /= numvz;
    }

    /* output */
    snprintf(filename, 1000, "%s_bh_%s", vname, "x");
    write_scalar(level, filename, avrvx);
    snprintf(filename, 1000, "%s_bh_%s", vname, "y");
    write_scalar(level, filename, avrvy);
    snprintf(filename, 1000, "%s_bh_%s", vname, "z");
    write_scalar(level, filename, avrvz);

    /* save lapse at BH */
    alpha_bh_x1 = avrvx;
    alpha_bh_y1 = avrvy;
    alpha_bh_z1 = avrvz;
    
    alpha_bh_x2 = avrvx;
    alpha_bh_y2 = avrvy;
    alpha_bh_z2 = avrvz;
  }

  return 0;
}
