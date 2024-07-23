/* robin.c */
/* Bernd Bruegmann 1/03 */

/* Robin boundary condition for 1/r fall-off
   assume
     u(x,y,z) = u(r) = uinfinity + c/r
   which implies
     du/dx = (uinfinity - u)/r nx,  nx = x/r
*/

#include "bam.h"
#include "boundary.h"

#define EXC_SPHERE 0
#define EXC_CUBE   1

#define EXC_INTERIOR 0
#define EXC_EXTERIOR 1
#define EXC_BOUNDARY 2




/* find direction that is closest to a certain normal 
   here: average over available directions to estimate normal
   BB: why not average data along all directions? or some directions?
   fills variable with index to point in best direction
*/ 
void find_robin_normal(tL *level)
{
  double *dir = PtrEnable(level, "robinindex");
  tVarList *vlboundary = VLPtrEnable1(level, "robinflag");
  double *noghostboundary = VLPtr(vlboundary, 0);
  int ia[3][3][3];
  double vx, vy, vz, s, max, dotproduct;
  int i, normalindex;
  int pr = 0;

  if (pr) 
    printf("Computing normal for Robin boundary, level->l = %d\n", level->l);

  /* copy boundary flag into temporary variable so that we can synch (!) 
     we will loop over all phybound points as defined in standard boundary
     flag, which omits ghosts, but when we check whether a direction is 
     available we check against this "no ghost"/synchronized/symmetrized flag
  */
  forallpoints(level, i)
    noghostboundary[i] = boundaryflag(level, i);
  set_boundary_symmetry(level, vlboundary);
  bampi_vlsynchronize(vlboundary);


  /* estimate normal by averaging direction to all unmasked neighbors */
  forall27flag(level, PHYBOUND) {
    fillindexarray27(ia);
    vx = vy = vz = 0;

    /* vector in normal direction */
    if (0) {
      int i, j, k;
      for (i = 0; i < 3; i++)
      for (j = 0; j < 3; j++)
      for (k = 0; k < 3; k++) {
	if (noghostboundary[ia[i][j][k]] != PHYBOUND) {
	  vx += i-1;
	  vy += j-1;
	  vz += k-1;
	}
      }
    }

    /* BB: doesn't the above overemphasize diagonal directions a little bit?
       let's normalize each direction before taking the average
       note that the center point is automatically excluded, hence s != 0
    */
    if (1) {
      int i, j, k;
      for (i = 0; i < 3; i++)
      for (j = 0; j < 3; j++)
      for (k = 0; k < 3; k++) {
	if (noghostboundary[ia[i][j][k]] != PHYBOUND) {
	  s = sqrt((double)((i-1)*(i-1)+(j-1)*(j-1)+(k-1)*(k-1)));
	  vx += (i-1)/s;
	  vy += (j-1)/s;
	  vz += (k-1)/s;
	  if (pr) printf("%d %d %d\n", i, j, k);
	}
      }
      if (pr) printf("vxyz %2d %2d %2d  ", (int) vx, (int) vy, (int) vz);
    }

    /* find a grid direction that comes closest to the normal direction */
    /* this will be asymmetric walking around sphere, and it will
       dependent on round-off if several directions give the same 
       dot product
    */
    if (1) {
      int i, j, k;
      normalindex = ccc;
      max = 0;
      for (i = 0; i < 3; i++)
      for (j = 0; j < 3; j++)
      for (k = 0; k < 3; k++) {
	if (noghostboundary[ia[i][j][k]] != PHYBOUND) {
	  s = sqrt((double)((i-1)*(i-1)+(j-1)*(j-1)+(k-1)*(k-1)));
	  dotproduct = (vx*(i-1) + vy*(j-1) + vz*(k-1))/s;
	  
	  if (dotproduct > max) {
	    max = dotproduct;
	    normalindex = ia[i][j][k];
	  }
	}
      }
    }

    /* sanity check */
    if (normalindex == ccc)
      errorexit("Robin boundary did not find a good normal direction!\n");
    if (boundaryflag(level, normalindex) != 0 &&
	boundaryflag(level, normalindex) != SYMBOUND) 
      printf("Warning: Robin boundary picked a bad direction!\n");

    if (pr) {
      double *x = Ptr(level, "x");
      double *y = Ptr(level, "y");
      double *z = Ptr(level, "z");
      
      printf("%6.3f %6.3f %6.3f  %6.3f %6.3f %6.3f  %6.3f %6.3f %6.3f\n",
	     x[ccc], y[ccc], z[ccc], vx, vy, vz, 
	     x[normalindex], y[normalindex], z[normalindex]);
    }

    /* save index */
    dir[ccc] = normalindex;

  } endfor;
  
  /* clean up */
  VLDisableFree(vlboundary);
}




/* set Robin boundary for one variable */
void set_boundary_robin(tL *level, tVarList *ul) 
{
  double *xp = Ptr(level, "x");
  double *yp = Ptr(level, "y");
  double *zp = Ptr(level, "z");
  double *index = Ptr(level, "robinindex");
  double *u, ufalloff, uinfinity;
  int i, j, ui;
  double x0, y0, z0, x1, y1, z1, x, y, z, vx, vy, vz, c;
  double k1,r1;

  /* set Robin boundary for scalar variable */
  if (ul->n == 1) {
    ui = ul->index[0];

    /* catch those variables that are passed in but don't have
       proper boundary flags set 
    */
    if (0) 
      printf("level %d, Robin for scalar %s: off %e far %e\n", 
	     level->l, VarName(ui), VarFallOff(ui), VarFarLimit(ui));

    ufalloff = VarFallOff(ui);

    uinfinity = VarFarLimit(ui);
    u = level->v[ui];

    /* do not handle extreme cases */
    if (ufalloff < -5)
      return;

    /* if ufalloff = 0, use a dirichlet b.c., with the boundary value k1 */
    else if (ufalloff == 0) {
      k1 = Getd("robin_1oR_constant");
      forallpoints(level, i)
      {
	if (boundaryflag(level, i) == PHYBOUND)
        {
          j = index[i];
	  x1 = xp[i]; y1 = yp[i]; z1 = zp[i];
	  x0 = xp[j]; y0 = yp[j]; z0 = zp[j];
	  
	  x = 0.5*(x1+x0);
	  y = 0.5*(y1+y0);
	  z = 0.5*(z1+z0);
	  
	  r1 = sqrt(x*x + y*y + z*z);

          u[i] = -u[j] + 2*k1/r1;
	}
      }
    }

    else {
      forallpoints(level, i)
      {
	if (boundaryflag(level, i) == PHYBOUND)
        {
	  j = index[i];
	
	  x1 = xp[i]; y1 = yp[i]; z1 = zp[i];
	  x0 = xp[j]; y0 = yp[j]; z0 = zp[j];
	  
	  x = 0.5*(x1+x0);
	  y = 0.5*(y1+y0);
	  z = 0.5*(z1+z0);
	  
	  vx = x1-x0;
	  vy = y1-y0;
	  vz = z1-z0;
	  
	  c = 0.5 * ufalloff * (vx*x + vy*y + vz*z)/(x*x + y*y + z*z);
	  
	  /* 1-point Robin formula */
	  u[i] = 2*c*uinfinity + u[j]*(1-c)/(1+c);
	}
      }
    }
  }

  /* set Robin boundary for vector variable */
  /* for now treat each component as scalar */
  else {
    tVarList *vl = vlalloc(level);
    vlpushone(vl, 0);

    for (i = 0; i < ul->n; i++) {
      ui = ul->index[i];
      vl->index[0] = ui;

      if (0) {
	printf("Robin for vector component %s: off %e far %e\n", 
	       VarName(ui), VarFallOff(ui), VarFarLimit(ui));
      }

      /* recursive call for scalar */
      set_boundary_robin(level, vl);       
    }

    vlfree(vl);

#if 0
    /* candidate for most ugly hack of the year */
    for (i = 0; i < ul->n; i++) {
      ui = ul->index[i];
      u = level->v[ui];
      forallpoints(level, j) {
	if (level->node[j].boundary == PHYBOUND) {
	  u[j] = 0;
	}
      }
    }
    ui = ul->index[0];
    if (strcmp(VarName(ui), "betax") == 0 ||
	0 && strcmp(VarName(ui), "betax_iterative_ph") == 0 ||
	0 && strcmp(VarName(ui), "betax_iterative_sh") == 0) {
      double *ux = level->v[ul->index[0]];
      double *uy = level->v[ul->index[1]];
      double *uz = level->v[ul->index[2]];
      forallpoints(level, j) {
	if (level->node[j].boundary == PHYBOUND) {
	  void BY_Wofxyz(double, double, double, double *, double *, double *);
	  BY_Wofxyz(xp[j], yp[j], zp[j], ux+j, uy+j, uz+j); 
	  ux[j] *= 2;
	  uy[j] *= 2;
	  uz[j] *= 2;
	}
      }
    }
#endif

  }



  //prvare(level, VarName(ul->index[0]));
  bampi_vlsynchronize(ul);
  //prdivider(0);
  //prvare(level, VarName(ul->index[0]));
}




/* set scalar Robin boundary for all components of VarList 
   for all the points in tMlPointList *Mlplist              */
void set_boundary_ScalarRobin(tL *level, tMlPointList *Mlplist, tVarList *ul) 
{
  double *xp = Ptr(level, "x");
  double *yp = Ptr(level, "y");
  double *zp = Ptr(level, "z");
  double *index = Ptr(level, "robinindex");
  double *u, ufalloff, uinfinity;
  int i, j, ui;
  double x0, y0, z0, x1, y1, z1, x, y, z, vx, vy, vz, c;

  /* set Robin boundary for scalar variable */
  if (ul->n == 1)
  {
    ui = ul->index[0];

    /* catch those variables that are passed in but don't have
       proper boundary flags set 
    */
    if (0) printf("Robin for scalar %s: off %e far %e\n", 
		  VarName(ui), VarFallOff(ui), VarFarLimit(ui));

    ufalloff = VarFallOff(ui);
    if (ufalloff == 0) return;

    uinfinity = VarFarLimit(ui);
    u = level->v[ui];

    forMlPointList(level, Mlplist, i)
    {
      if (boundaryflag(level, i) == PHYBOUND)
      {
	j = index[i];
	
	x1 = xp[i]; y1 = yp[i]; z1 = zp[i];
	x0 = xp[j]; y0 = yp[j]; z0 = zp[j];
	
	x = 0.5*(x1+x0);
	y = 0.5*(y1+y0);
	z = 0.5*(z1+z0);
	
	vx = x1-x0;
	vy = y1-y0;
	vz = z1-z0;
	
	c = 0.5 * ufalloff * (vx*x + vy*y + vz*z)/(x*x + y*y + z*z);
	
	/* 1-point Robin formula */
	u[i] = 2*c*uinfinity + u[j]*(1-c)/(1+c);
      }
    } endfor;
  }

  /* set Robin boundary for vector variable */
  /* for now treat each component as scalar */
  else
  {
    tVarList *vl = vlalloc(level);
    vlpushone(vl, 0);

    for (i = 0; i < ul->n; i++)
    {
      ui = ul->index[i];
      vl->index[0] = ui;

      if (0) {
	printf("Robin for vector component %s: off %e far %e\n", 
	       VarName(ui), VarFallOff(ui), VarFarLimit(ui));
      }

      /* recursive call for scalar */
      set_boundary_ScalarRobin(level, Mlplist, vl);       
    }
    vlfree(vl);
  }

  bampi_vlsynchronize(ul);
}
