/* commotion.c */
/* Bernd Bruegmann, 11/03 , Wolfgang Tichy 12/2003 */

#include "bam.h"
#include "Gauge.h"

tVarList *vlu_c_old = 0, *vlu_p_old = 0;
extern double *commotion_x, *commotion_y;
extern int commotion_n, commotion_nmax;



/***************************************************************************/
/* adjust commotion parameters manually 
*/
int commotion_manual(tL *level)
{
  tG *g = level->grid;
  static int firstcall = 1, n;
  static double *t = 0, *s = 0, *r = 0;
  double domega, drdot, dscale;
  int i;
  int pr = 1;

  /* initialize */
  if (firstcall) {
    char *rlist = Gets("commotion_manual_r");
    char *slist = Gets("commotion_manual_s");
    char *tlist = Gets("commotion_manual_t");
    char *a;

    for (n = 0; NextEntry(tlist); n++);
    for (i = 0; NextEntry(rlist); i++);
    if (i != n) errorexit("commotion_manual: need equal length lists");
    for (i = 0; NextEntry(slist); i++);
    if (i != n) errorexit("commotion_manual: need equal length lists");
    r = dmalloc(n);
    s = dmalloc(n);
    t = dmalloc(n);

    for (i = 0; a = NextEntry(rlist); i++) sscanf(a, "%lf", r+i);
    for (i = 0; a = NextEntry(slist); i++) sscanf(a, "%lf", s+i);
    for (i = 0; a = NextEntry(tlist); i++) sscanf(a, "%lf", t+i);
    if (pr) {
      printf("Commotion parameters:\n");
      for (i = 0; i < n; i++)
	printf("  t = %8.6f  s = %8.6f  r = %8.6f\n", t[i], s[i], r[i]);
    }
    firstcall = 0; 
  }
    
  /* check whether there is anything to do */
  if (!Getv("corotation_monitor", "manual")) return 0;
  if (level->l != g->lmax ||
      !dequal(level->time, g->level[g->lmin]->time)) return 0;

  /* only do something if we have reached one of the given times */
  for (i = 0; i < n; i++)
    if (dequal(level->time, t[i])) break;
  if (i == n) return 0;

  /* compute corrections as determined by entry i */
  if (pr) printf("cr = %8.6f -> %8.6f,   ", Getd("corotation_scale"), s[i]);
  if (pr) printf("cm = %8.6f -> %8.6f\n",   Getd("corotation_rdot"),  r[i]);
  drdot  = r[i] - Getd("corotation_rdot");
  dscale = s[i] - Getd("corotation_scale");
  domega = dscale * Getd("corotation_omega");

  /* adjust parameters */
  Setd("corotation_rdot",  r[i]);
  Setd("corotation_scale", s[i]);
  
  /* adjust shift on all levels */
  adjust_beta(g, domega, drdot);

  return 0;
}




/***************************************************************************/
/* try something simple:
     damping in direction opposite of drift + 
     driving in direction of puncture
   assumes that finest is time aligned with coarsest because
   otherwise coarse level would evolve with different shift than fine
*/
int commotion_instantaneous(tL *level)
{
  tG *g = level->grid;
  tL *toplevel = g->level[g->lmin];
  double t  = toplevel->time;
  double dt = toplevel->dt;
  double t0, t1, dt0;
  double rdot, omega, scale, scale_max;
  double drdot, domega, dscale;
  double x0, x1, x2, x3, y0, y1, y2, y3;
  double xoffset, yoffset;
  double dx, dxp, dxp0, dy, dyp, dyp0;
  double vx, vx0, vy, vy0;
  double dvx, dvy;
  double c, d;
  int i0, i1, i2;
  int i, j;
  static int scale_max_flag = 0;
  int pr = 0;

  /* check whether there is anything to do */
  if (!Getv("corotation_monitor", "adjust")) return 0;
  if (level->l != g->lmax ||
      !dequal(level->time, g->level[g->lmin]->time)) return 0;

  /* initialize t1 */
  t1  = Getd("corotation_adjust_t1");
  
  /* if t>=t1 give up on correcting */
  if( (t>=t1) && (t1>0.0) )
  {
    char str[1000];
    char before[1000];
    char *after;
    char *at;
    
    snprintf(str, 1000, "%s", Gets("corotation_monitor") );
    at=strstr(str, "adjust");
    if(at==NULL) 
      errorexit("commotion.c: why is there no adjust in corotation_monitor?");
    after=at+6; /* adjust has 6 chars */
    at[0]=0;
    snprintf(before, 1000, "%s", str);
    snprintf(str, 1000, "%s %s", before, after );
    Sets("corotation_monitor", str);
    printf("commotion_instantaneous: stopping correction  t=%g  >=  t1=%g\n", 
           t, t1 );
    printf("setting  corotation_monitor = %s\n",
           Gets("corotation_monitor") );
    
    return 0;
  }

  /* initialize */
  t0  = Getd("corotation_adjust_t0");
  dt0 = Getd("corotation_adjust_dt"); 
  if (dless(dt0, dt)) {
    dt0 = dt;
    Setd("corotation_adjust_dt", dt);
  }
  if (dless(t0, dt0)) {
    t0 = dt0;
    Setd("corotation_adjust_t0", t0);
  }

  /* wait until t0+dt0 before doing first correction */
  if (dless(t, t0+dt0)) return 0;

  /* wait for multiple of dt0 */
  i = (t-t0)/dt + 0.5;
  j = dt0/dt + 0.5;
  if (i % j) return 0;
  if (pr) printf("adjust commotion: t = %f\n", level->time);
  if (pr) printf("dt0 = %f,  dt = %f,  j = %d, i = %d\n", dt0, dt, j, i);

  /* offset if wanted */
  xoffset = Getd("commotion_offset_x");
  yoffset = Getd("commotion_offset_y");

  /**********************************************************************/
  /* new */

  /* use up to last three points */
  /* previously: use first, middle and final point of previous interval
     because data may have been too noisy or too local */
  /* better: fit parabola? */
  i2 = commotion_n - 1;
  i1 = i2 - 1;
  i0 = i1 - 1;
  x0 = commotion_x[i0] - xoffset;
  x1 = commotion_x[i1] - xoffset;
  x2 = commotion_x[i2] - xoffset;
  y0 = commotion_y[i0] - yoffset;
  y1 = commotion_y[i1] - yoffset;
  y2 = commotion_y[i2] - yoffset;
  
  /* compute tangent to drift curve */
  if (Getv("commotion_DriftVelocity", "FromPreviousCorrections") ||
      Getv("commotion_DriftVelocity", "2ndorder")) 
  {
    /* one-sided second order derivative */
    vx = (3*x2 - 4*x1 + x0)/(2*dt);
    vy = (3*y2 - 4*y1 + y0)/(2*dt);
    if (pr) printf("x: %f %f %f  %f\n", x0, x1, x2, vx);
    if (pr) printf("y: %f %f %f  %f\n", y0, y1, y2, vy);
  }

  else if (Getv("commotion_DriftVelocity", "Instantaneous") ||
	   Getv("commotion_DriftVelocity", "1storder"))
  {
    /* one-sided first order derivative */
    vx = (x2 - x1)/dt;
    vy = (y2 - y1)/dt;
  }
  else
    errorexit("commotion_instantaneous: commotion_DriftVelocity is set wrong");
  
  /* compute direction toward place where we want the puncture to be */
  dxp0 = -(x2 - Getd("bhx1"));
  dyp0 = -(y2 - Getd("bhy1"));

  /* compensate within the next interval of dt0 */
  vx0 = dxp0/dt0;
  vy0 = dyp0/dt0;

  /* correction = damping in direction opposite of drift + 
                  driving in direction of puncture location
     corresponds to backward Euler:
       beta^n = beta^(n-1) + dt0 * (- c v^n - d (x^n - x0)) 
  */
  c = Getd("commotion_damping_x");
  d = Getd("commotion_driving_x");
  //d = pow(2*PI*d, 2);
  dvx = - c * vx + d * vx0;
  dvx *= dt0;
  c = Getd("commotion_damping_y");
  d = Getd("commotion_driving_y");
  //d = pow(2*PI*d, 2);
  dvy = - c * vy + d * vy0;
  dvy *= dt0;

  /* info */
  if (pr) {
    printf("vx = %f,  vx0 = %f,  dvx = %f\n", vx, vx0, dvx);
    printf("vy = %f,  vy0 = %f,  dvy = %f\n", vy, vy0, dvy);
  }


  /*--------------------------------------------------------------------*/
  /* old */
  if (Getv("corotation_monitor", "old2003")) {

    /* compute tangent to drift curve */
    if(Getv("commotion_DriftVelocity", "FromPreviousCorrections"))
    {
      /*  better: fit parabola  */
      i2 = commotion_n - 1;
      i1 = i2 - 1;
      i0 = i1 - 1;
      /* BB 12/05: something broke, perhaps dt0 used to be 2*dt0? */
      if (Getv("corotation_monitor", "old2003")) {
	i0 = i2 - j;
	  i1 = (i2+i0)/2;
      }
      if (pr) printf("i0 = %d, i1 = %d, i2 = %d\n", i0, i1, i2);
      x0 = commotion_x[i0] - xoffset;
      x1 = commotion_x[i1] - xoffset;
      x2 = commotion_x[i2] - xoffset;
      y0 = commotion_y[i0] - yoffset;
      y1 = commotion_y[i1] - yoffset;
      y2 = commotion_y[i2] - yoffset;
      dx = (3*x2 - 4*x1 + x0)/dt0;
      dy = (3*y2 - 4*y1 + y0)/dt0;
      if (!Getv("corotation_monitor", "old2003")) {
	dx /= 2;
	dy /= 2;
      }
      if (pr) printf("x: %f %f %f  %f\n", x0, x1, x2, dx);
      if (pr) printf("y: %f %f %f  %f\n", y0, y1, y2, dy);
    }
    else if(Getv("commotion_DriftVelocity", "Instantaneous"))
    {
      i2 = commotion_n - 1;
      i1 = i2 - 1;
      
      x1 = commotion_x[i1] - xoffset;
      x2 = commotion_x[i2] - xoffset;
      y1 = commotion_y[i1] - yoffset;
      y2 = commotion_y[i2] - yoffset;
      
      dx = (x2 - x1)/dt;
      dy = (y2 - y1)/dt;
    }
    else
      errorexit(
	"commotion_instantaneous: commotion_DriftVelocity is set wrong");
    
    /* compute direction toward puncture */
    dxp = -x2;
    dyp = -y2;  
    /* BB 12/05: this looks like a bug, or at least an unintentional feature */
    if (!Getv("corotation_monitor", "old2003")) {
      dyp += Getd("bhy1");
    }
    if (pr) printf("dx = %f,  dxp = %f\n", dx, dxp);
    if (pr) printf("dy = %f,  dyp = %f\n", dy, dyp);
    
    /* correction = damping in direction opposite of drift + 
       driving in direction of puncture
    */
    c = Getd("commotion_damping_x");
    d = Getd("commotion_driving_x");
    dvx = - c * dx + d * dxp;
    dvx *= dt0;
    c = Getd("commotion_damping_y");
    d = Getd("commotion_driving_y");
    dvy = - c * dy + d * dyp;
    dvy *= dt0;
  }
  /* end old */
  /**********************************************************************/


  /* set new corotation parameters */

  /* angular scale */
  y0 = commotion_y[commotion_n - 1];
  omega = Getd("corotation_omega");
  scale = Getd("corotation_scale");
  scale_max = Getd("corotation_adjust_scale_max");
  domega = 0;
  dscale = 0;
  if (y0 >= Getd("corotation_adjust_scale_y_min")) {
    if (scale < scale_max && !scale_max_flag) {
      domega = dvx/y0;
      dscale = domega/omega;
      scale += dscale;
      Setd("corotation_scale", scale);
      // stop adjusting once scale_max was reached
      if (scale > scale_max)
	scale_max_flag = 1;
    }
  }
  if (pr) printf("domega/omega0 = %e -> omega/omega0  %e\n", dscale, scale);

  /* radial rdot */
  rdot = Getd("corotation_rdot");
  drdot = dvy;
  rdot += drdot;
  Setd("corotation_rdot", rdot);
  if (pr) printf("drdot = %e -> rdot  %e\n", dvy, rdot);

  /* adjust shift on all levels */
  adjust_beta(g, domega, drdot);

  /* scale radial shift on all levels */
  if( Getd("corotation_shift_scale") != 0.0 ) 
    scale_radial_shift(g);

  /* increase/decrease radial shift on all levels */
  if( Getd("corotation_shift_additionFactor") != 0.0 ) 
    addto_radial_shift(g);

  /* adjust radial shift on all levels */
  if( Getd("corotation_shift_radialAdjustment") != 0.0 ) 
    adjust_radial_shift(g);
  
  return 0;
}




/***************************************************************************/
/* try something that should avoid adjustment induced oscillations
   restart evolutions at earlier time until happy with outcome
*/


/* save state information 
   initializes storage first time through
*/
void save_state(tG *g)
{
  tVarList *vlu_c = 0, *vlu_p = 0;
  tL *level;
  int l;

  evolve_vlretrieve(&vlu_c, &vlu_p, 0);

  if (!vlu_c_old) vlu_c_old = AddDuplicate(vlu_c, "_old");
  if (!vlu_p_old) vlu_p_old = AddDuplicate(vlu_p, "_old");

  for (l = g->lmin; l <= g->lmax; l++) {
    level = g->level[l];
    vlenablelevel(level, vlu_c_old);
    vlenablelevel(level, vlu_p_old);
    vlcopylevel(level, vlu_c_old, vlu_c);
    vlcopylevel(level, vlu_p_old, vlu_p);
  }
}




/* restore state information */
void restore_state(tG *g)
{
  tVarList *vlu_c = 0, *vlu_p = 0;
  int l;

  evolve_vlretrieve(&vlu_c, &vlu_p, 0);
  for (l = g->lmin; l <= g->lmax; l++) {
    vlcopylevel(g->level[l], vlu_c, vlu_c_old);
    vlcopylevel(g->level[l], vlu_p, vlu_p_old);
  }
}




/* reset state vector to earlier time */
void reset_grid_time(tG *g)
{
  tL *level;
  int l;

  /* overwrite current state vector with old information */
  restore_state(g);

  /* pretend nothing happened */
  for (l = g->lmin; l <= g->lmax; l++) {
    level = g->level[l];
    level->iteration = 0;
    level->time = level->iteration * level->dt;
  }
}




/* reiterate */
int commotion_reiterate(tL *level)
{
  tG *g = level->grid;
  tL *toplevel = g->level[g->lmin];
  double t  = toplevel->time;
  double dt = toplevel->dt;
  double x0, y0;
  double t0, dt0;
  int i, j;

  /* only do it when all levels are aligned */
  if (level->l != g->lmax || !dequal(level->time, t)) return 0;

  /* get parameters */
  x0  = Getd("commotion_x");
  y0  = Getd("commotion_y");
  t0  = Getd("corotation_reiterate_t0");
  dt0 = Getd("corotation_reiterate_dt"); 

  /* initialize */
  if (dless(dt0, dt)) {
    dt0 = dt;
    Setd("corotation_reiterate_dt", dt);
  }
  if (dless(t0, dt0)) {
    t0 = dt0;
    Setd("corotation_reiterate_t0", t0);
  }

  /* at t0 we save the state for the first time */
  if (dequal(t, t0))
    save_state(g);

  /* wait until t0+dt0 before doing first correction */
  if (dless(t, t0+dt0)) return 0;

  /* wait for multiple of dt0 */
  i = (t-t0)/dt + 0.5;
  j = dt0/dt + 0.5;
  if (i % j) return 0;
  printf("commotion: t=%f, x0=%f, y0=%f\n", level->time, x0, y0);
  
  /* unfinished: here we decide what to do */

  reset_grid_time(g);

  return 0;
}




