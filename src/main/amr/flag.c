/* flag.c */
/* Bernd Bruegmann, 12/99, 10/03 */

/* set various flags 
   flag system should be redone (BB 10/03):
     - bit flags attached to nodes?
     - status flag should not be special
     - keep variables like flagregrid for complicated/temporary work? 
     - excision on different levels?
*/ 

#include "bam.h"
#include "amr.h"

#define PR 0




/* set boundary flags: so far only correct for top level cube! 
   see setboundaryflags() in regrid.c for general version
*/
void set_boundary_flags(tL *level)
{
  int i, j, k;

  /* check direct neighbors:
     if there is no direct neighbor in any one direction, this point is at
     the outer boundary
  */
  forallpoints_ijk(level) {
    if (i == 0 || i == box->m - 1 ||
	j == 0 || j == box->n - 1 ||
	k == 0 || k == box->o - 1)
      level->boundary[ijk] = PHYBOUND;
  } endfor_ijk;

  
  
  /* flag symmetry boundary, this will be overridden for ghost zones */
  if (!level->shells) {
    set_boundary_flags_symmetry(level);
    set_boundary_flags_periodic(level);
  }
  
  /* call bampi to mark all ghost zones */
  bampi_set_boundary_flags(level);

  /* set one flag for each face of each box
     tested for reflection boundaries
     FIX ME: this may fail for periodic, tubes, planar, ...
  */
  forallboxes(level) {

    if ((Getv("grid", "octant") || Getv("grid", "quadrant") ||
	Getv("grid", "bitant") || Getv("grid", "rotant") || 
        Getv("grid", "qreflect")) && !level->shells) {

      /* annoying chicken and egg problem that I don't want to fix now
	 we may be called before the grid structure knows its symmetries ...
      */
      int half[3];
      half[0] = Getv("grid_half", "x");
      half[1] = Getv("grid_half", "y");
      half[2] = Getv("grid_half", "z");

      for (i = 0; i < 6; i++)              box->bflag[i] = PHYBOUND;
      for (i = 0; i < 6; i += 2) 
	if (half[i/2] && box->bbox[i] < 0) box->bflag[i] = SYMBOUND;  
      for (i = 0; i < 6; i++)
	if (box->com->nbrank[i] >= 0)      box->bflag[i] = GHOBOUND;
    } else {
      int imax = box->m-1;
      int jmax = box->n-1;
      int kmax = box->o-1;
      int imid = imax/2 + 1;
      int jmid = jmax/2 + 1;
      int kmid = kmax/2 + 1;
      box->bflag[0] = level->boundary[ijkofbox(box, 0,    jmid, kmid)];
      box->bflag[1] = level->boundary[ijkofbox(box, imax, jmid, kmid)];
      box->bflag[2] = level->boundary[ijkofbox(box, imid, 0,    kmid)];
      box->bflag[3] = level->boundary[ijkofbox(box, imid, jmax, kmid)];
      box->bflag[4] = level->boundary[ijkofbox(box, imid, jmid, 0   )];
      box->bflag[5] = level->boundary[ijkofbox(box, imid, jmid, kmax)];
    }

    if (0) {
      printbbox(level, box->bbox, 0);
      printf("set boundary flag for box:  %2d %2d  %2d %2d  %2d %2d\n",
	     box->bflag[0], box->bflag[1], box->bflag[2], 
	     box->bflag[3], box->bflag[4], box->bflag[5]);
    }
  } endforboxes;
}








/* set flags needed for restriction and prolongation */
/*
   roughly speaking:
   Pflag != 0: obtain these fine points by prolongation
   Rflag != 0: obtain these coarse points by restriction

   Pflag counts the points from the boundary
   Rflag >= 1: there are refinement points
         >= 2: needed for multigrid
         >= 3: needed for evolution

   nbuffer = 3 (e.g. 3 evolution steps, stencil radius 1, note 3 is odd)
   course, Rflag:    3   3   3   3   2   1   0   0
   fine,   Pflag:   0 0 0 0 0 0 0 0 0 3 2 1

   nbuffer = 6 (e.g. 3 evolution steps, stencil radius 2)
   course, Rflag:    3   3   3   2   2   1   0   0
   fine,   Pflag:   0 0 0 0 0 0 6 5 4 3 2 1
*/
void setRPflags(tL *lc, tL *lf)
{
  if (PR) printf("set RP flags %d %d\n", lc->l,lf->l);
  
  int nbuffer = Geti("amr_nbuffer");
  double *p = Ptr(lf, "flagprolong");
  double *r = Ptr(lc, "flagrestrict");
  int *half = lf->grid->half;
  int i, ii, jj, kk, n, b;
  int ia[3][3][3];
  double x, y, z, *bbox;


  /* first, on fine introduce as many buffer nodes as needed for evolution 
     these are the nodes obtained by prolongation
  */

  /* initialize based on boundary
     outer boundary is currently just 1 thick
     ignore symmetry and ghost boundaries which may be 2 thick
  */
  for (i = 0; i < lf->nnodes; i++)
    p[i] = (boundaryflag(lf, i) == PHYBOUND ||
	    boundaryflag(lf, i) == REFBOUND) ? 1 : 0;

  /* set P flag to count from outer boundary 
     use global bounding box to make it work with parallelization
  */
  forallpoints_ijk(lf) {
    if (boundaryflag(lf, ijk) == GHOBOUND || boundaryflag(lf, ijk) == SYMBOUND)
      continue;
    ii = jj = kk = box->npoints;

    /* global bbox */
    bbox = box->bbox;
    
    x = box->com->bbox[0] + i*lf->dx;
    y = box->com->bbox[2] + j*lf->dy;
    z = box->com->bbox[4] + k*lf->dz;

    if (!half[0] || bbox[0] > 0) {
      n = 1.5 + (x - bbox[0])/lf->dx;
      if (n <= nbuffer) ii = n;
    }
    if (1) {
      n = 1.5 - (x - bbox[1])/lf->dx;
      if (n <= nbuffer) ii = n;
    }
      
    if (!half[1] || bbox[2] > 0) {
      n = 1.5 + (y - bbox[2])/lf->dy;
      if (n <= nbuffer) jj = n;
    }
    if (1) {
      n = 1.5 - (y - bbox[3])/lf->dy;
      if (n <= nbuffer) jj = n;
    }

    if (!half[2] || bbox[4] > 0) {
      n = 1.5 + (z - bbox[4])/lf->dz;
      if (n <= nbuffer) kk = n;
    }
    if (1) {
      n = 1.5 - (z - bbox[5])/lf->dz;
      if (n <= nbuffer) kk = n;
    }

#if 0
    /* awaiting cleanup when BOX/BOXES are made one
       this version works for BOXES, but fails for BOX */
    if (!half[0] || box->bbox[0] > 0) {
      n = 1.5 + (xofbox(box, i) - box->bbox[0])/box->dx;
      if (n <= nbuffer) ii = n;
    }
    if (1) {
      n = 1.5 - (xofbox(box, i) - box->bbox[1])/box->dx;
      if (n <= nbuffer) ii = n;
    }
      
    if (!half[1] || box->bbox[2] > 0) {
      n = 1.5 + (yofbox(box, j) - box->bbox[2])/box->dy;
      if (n <= nbuffer) jj = n;
    }
    if (1) {
      n = 1.5 - (yofbox(box, j) - box->bbox[3])/box->dy;
      if (n <= nbuffer) jj = n;
    }

    if (!half[2] || box->bbox[4] > 0) {
      n = 1.5 + (zofbox(box, k) - box->bbox[4])/box->dz;
      if (n <= nbuffer) kk = n;
    }
    if (1) {
      n = 1.5 - (zofbox(box, k) - box->bbox[5])/box->dz;
      if (n <= nbuffer) kk = n;
    }
#endif

    if (ii > jj) ii = jj;
    if (ii > kk) ii = kk;
    if (ii < box->npoints)
      p[ijk] = ii;
    else
      p[ijk] = 0;
  } endfor_ijk;

  if (0) prvar01(lf, "flagprolong");


  /* second, on coarse determine nodes that are obtained by restriction
     needs odd nbuffer for now
     R = 0: no restrict, these points are ghosts for interpolation stencil
     R = 3,2,1: inner restrict region with buffer zone
  */
  // new 5/4/04, 12/14/05 BB

  /* initialize to 0 */
  forallpoints(lc, i)
    r[i] = 0;

  /* if P flag is 0, set to 3 */
  forallpoints_ijk(lf) {
    b = boundaryflag(lf, ijk);
    if (b != SYMBOUND && b != GHOBOUND && p[ijk] == 0)
      r[box_parent_ijk(box->pr, box, i, j, k)] = 3;
  } endfor_ijk;

  /* if P flag is >= 3, set to 2 
     (do it here so that it is not overwritten for nbuffer odd) */
  /* if P flag is 2, set to 1 */
  forallpoints_ijk(lf) {
    b = boundaryflag(lf, ijk);
    if (b != SYMBOUND && b != GHOBOUND) { 
      if (p[ijk] >= 3)       r[box_parent_ijk(box->pr, box, i, j, k)] = 2; 
      else if (p[ijk] == 2) r[box_parent_ijk(box->pr, box, i, j, k)] = 1;
    }
  } endfor_ijk;

  if (0) prvar01(lc, "flagrestrict");
}




/*************************************************************************/
/* flag points for refinement
   implement simple choices here
*/

/* helper function:
   turn off regrid flags that are too close to various boundaries
   just do nghosts = 2 or for BOX, nghosts = even 
   (odd would be problematic)
   for symmetry boundaries and ghost zones:
   set nghosts/2 flags so that nghosts points are created in refinement
   unsetting nghost/2 flags near outer or refinement boundaries
   guards against refinement too close to boundary
*/
void flagregridboundary(tL *level)
{
  double *f = Ptr(level, "flagregrid");
  int nghosts = Geti("bampi_nghosts");
  int n = nghosts/2;

  if (nghosts % 2 != 0)
    errorexit("flagregridboundary: requires nghosts = even");

  forallpoints_ijk(level) {
    if (f[ijk]) {
      if (i < n || j < n || k < n || imax-i < n || jmax-j < n || kmax-k < n)
	f[ijk] = 0;
    }
  } endfor_ijk;

}




/* flag nested boxes or spheres of constant number of points */
int flagsimplegrid(tL *level) 
{
  double *f = Ptr(level, "flagregrid");
  double *x = Ptr(level, "x");
  double *y = Ptr(level, "y");
  double *z = Ptr(level, "z");
  double r, r0 = level->dx * Geti("nx")/4; 
  int i;

  if (Geti("ny") != Geti("nx") || Geti("nz") != Geti("nx") ||
      Geti("ny") != Geti("nz"))
    errorexit("use nxyz for flagsimplegrid in flag.c");

  /* no refinement of finest level, results in all flags zero */
  if (level->l >= Geti("amr_lmax")) r0 = -1;

  /* nested boxes */
  if (Getv("amr_fmr", "nestedboxes")) {
    r0 += level->dx/8;
    for (i = 0; i < level->nnodes; i++) {
      f[i] = (x[i] < r0 && x[i] > -r0 &&
	      y[i] < r0 && y[i] > -r0 &&
	      z[i] < r0 && z[i] > -r0) ? 1 : 0;
    }
  }
  
  /* nested spheres */
  else if (Getv("amr_fmr", "nestedspheres")) {
    for (i = 0; i < level->nnodes; i++) {
      r = sqrt(x[i]*x[i] + y[i]*y[i] + z[i]*z[i]);
      f[i] = (r < r0) ? 1 : 0;
    }
  }

  /* graded boxes */
  else if (Getv("amr_fmr", "gradedboxes")) {
    flaggradedboxes(level);
    return 0;
  }

  /* unknown */
  else
    errorexits("flagsimpleregrid: amr_fmr = %s not known", Gets("amr_fmr"));


  /* adjust flags at boundaries */
  flagregridboundary(level);

  return 0;
}




/* flag centered boxes whose size varies gradually between two values
  
   given nx for the coarsest and nxfinest for the finest level, 
   increase the number of points by a constant factor per refinement

   example: nx = 64, nxfinest = 256, lmax = 3
            set constant factor to (256/64)^(1/3) = 1.59
        =>  64, 102, 161, 256

   alternatively, specify amr_fmr_rx_min etc which implies nxfinest
   note that different symmetry boundaries are handled implicitly
*/
int flaggradedboxes(tL *level) 
{
  double *f = Ptr(level, "flagregrid");
  double *x = Ptr(level, "x");
  double *y = Ptr(level, "y");
  double *z = Ptr(level, "z");
  double rx, ry, rz;
  double p, r;
  int ncube = Geti("order_RP")/2;
  int i, n;

  /* return if no more refinements wanted (leaves all flags zero) */
  if (level->l >= Geti("amr_lmax")) return 0;

  /* one can specify the smallest size either in points or coordinates
     transform the latter into the former
  */
  p = pow(2.0, Getd("amr_lmax") - level->l + 1.0);
  if ((r = Getd("amr_fmr_rx_min")) > 0)
    Seti("amr_fmr_nx_finest", p*r/level->dx + 2.5);
  if ((r = Getd("amr_fmr_ry_min")) > 0)
    Seti("amr_fmr_ny_finest", p*r/level->dy + 2.5);
  if ((r = Getd("amr_fmr_rz_min")) > 0)
    Seti("amr_fmr_nz_finest", p*r/level->dz + 2.5);

  /* the default is to keep the number of points constant */
  rx = level->dx * Geti("nx") / 4;
  ry = level->dy * Geti("ny") / 4;
  rz = level->dz * Geti("nz") / 4;

  /* but if the number of points on the final level changes,
     adjust with a scale factor for a gradual transition to that size
  */
  p = (1.0+level->l)/Getd("amr_lmax");
  if ((n = Geti("amr_fmr_nx_finest")) > 0) rx *= pow(n/Getd("nx"), p);
  if ((n = Geti("amr_fmr_ny_finest")) > 0) ry *= pow(n/Getd("ny"), p);
  if ((n = Geti("amr_fmr_nz_finest")) > 0) rz *= pow(n/Getd("nz"), p);

  /* sanity check */
  if (rx > level->bbox[1] - (ncube+0.5)*level->dx ||
      ry > level->bbox[3] - (ncube+0.5)*level->dy ||
      rz > level->bbox[5] - (ncube+0.5)*level->dz)
    errorexit("flaggradedboxes: boxes are getting too large!");
        
  /* guard against round-off */
  rx += level->dx/8;
  ry += level->dy/8;
  rz += level->dz/8;

  /* set flags */
  for (i = 0; i < level->nnodes; i++)
    if (x[i] < rx && x[i] > -rx &&
	y[i] < ry && y[i] > -ry &&
	z[i] < rz && z[i] > -rz)
      f[i] = 1;

  /* adjust flags at boundaries */
  flagregridboundary(level);

  return 0;
}




/* flag points for refinement centered around some point */
int FlagRegrid_BoxesAroundCenter(tL *level) 
{
  double *f = Ptr(level, "flagregrid");
  double x, *px = Ptr(level, "x");
  double y, *py = Ptr(level, "y");
  double z, *pz = Ptr(level, "z");
  double r, r0, rx, ry; 
  double ox, oy, oz;
  int i, ncube;

  /* simple nesting strategy: boxes/spheres of equal size at each level */
  r0 = level->dx * Geti("nx")/4;

  /* initialize flags to zero and return if nothing else to do */
  forallpoints(level, i) 
    f[i] = 0;
  if (level->l >= Geti("amr_lmax")) 
    return 0;

  /* get center of refinement region */
  ox = Getd("amr_fmr_centerx");
  oy = Getd("amr_fmr_centery");
  oz = Getd("amr_fmr_centerz");

  /* refinements need buffer to boundary for 4^3 interpolation box */
  ncube = Geti("order_RP")/2;
  if (dless(level->bbox[3], r0 + oy + (ncube-0.5)*level->dx))
    errorexit("FlagRegrid_BoxesAroundCenter: outer boundary is too close!");

  /* nested boxes */
  if (Getv("amr_fmr", "nestedboxes"))
  {
    double rx_min = Getd("amr_fmr_rx_min");
    double ry_min = Getd("amr_fmr_ry_min");

    /* guard against round-off ?! */
    r = r0 + level->dx/100; 
    /*      + fmod(r0-ox, level->dx/10)/100
            + fmod(r0-oy, level->dy/10)/100
	      + fmod(r0-oz, level->dz/10)/100; */

    rx = r; /* width of box (in x-direction) = rx+rx */
    ry = r; /* length of box (in y-direction) = ry+r */

    if( rx <= rx_min + 2.0*level->dx ) 
      rx = rx_min + 2.0*level->dx + level->dx/100.0;
    if( ry <= ry_min + 2.0*level->dy )  
      ry = ry_min + 2.0*level->dy + level->dy/100.0;

    printf("r = %e   rx = %e   ry = %e\n", r, rx, ry);

    /* set flags */
    forallpoints(level, i)
    {
      x = px[i] - ox;
      y = py[i] - oy;
      z = pz[i] - oz;
      if(x < rx && x > -rx  &&  y < ry && y > -ry  &&  z < r && z > -r)
        f[i] = 1;
    }
  }
  
  /* nested spheres */
  else if (Getv("amr_fmr", "nestedspheres")) {
    forallpoints(level, i) {
      x = px[i] - ox;
      y = py[i] - oy;
      z = pz[i] - oz;
      r = sqrt(x*x + y*y + z*z);
      if (r < r0)  f[i] = 1;
    }
  }

  /* unknown */
  else
    errorexits("FlagRegrid_BoxesAroundCenter:: amr_fmr = %s not known",
                Gets("amr_fmr"));
  

  /* adjust flags at boundaries */
  flagregridboundary(level);

  if (0) prvar01(level, "flagregrid");
  return 0;
}
