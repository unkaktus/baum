/* symmetry.c */
/* Bernd Bruegmann 11/02 */

/* Boundary conditions determined by reflection symmetry and rotations

   Possible reflections are about the x=0, y=0, z=0 planes
     e.g. bitant (z=0), octant (all three planes)

   Assumes that the grid either is aligned with the plane (so that
   points fall into the plane) or that the symmetry plane is staggered
   (so that it lies half way between grid points). In these cases
   the data is obtained by a simple copy operation of data that is 
      2*x/ds points
   away.

   Assumes that all data is available locally on this processor, in particular
   we do not check whether those points actually exist.
   (This breaks down for inversion through an axis.)
*/

#include "bam.h"
#include "boundary.h"




/* set reflection boundary for one plane */
void set_boundary_reflection_one(tL *level, tVarList *varlist, 
				 double *xp, double dx, int dir) 
{
  double *var;
  int nvars = varlist->n;
  int *ivar = varlist->index; 
  int nv;
  double x, g = -dx/4, h = dx/2;
  int i, ni, s, steps;
  double *sym;

  /* read sign for reflection 
     should probably store this somewhere else, inline, export, ...
  */
  sym = malloc(sizeof(double)*varlist->n);
  for (nv = 0; nv < nvars; nv++) {
    sym[nv] = VarSymmetry(ivar[nv], dir/2);
    if (0) printf("reflection for %s, dir %d, sym %3.0lf\n",
		  VarName(ivar[nv]), dir, sym[nv]);
  }

  for (nv = 0; nv < nvars; nv++) {
    double *var = level->v[ivar[nv]];
    int n, nneg, nmax, dn;

    /* for all boxes */
    forallboxes(level) {
      dn = box->id[dir/2];
      
      /* only the processor at the lower end should apply symmetry
	 there might be cases where the ghosts of the second processor
         reach into the symmetry ghosts, so don't rely on negative lower
         bound
      */
      if (box->com->nbrank[dir-1] >= 0) continue;

      /* number of points with negative coordinate */
      nneg = (int) (1 - box->bbox[dir-1]/dx);  // assumes staggering
      
      /* index of reflected point at lower end */
      nmax = 2*nneg -1;
      if (0) printf("symmetry: xmin %f, nneg %d nmax %d\n",
		    box->bbox[dir-1], nneg, nmax);

      /* for all points in the given plane */
      forplane_boxijk(box, dir) {
	
	/* for all points along orthogonal direction */
	for (n = 0; n < nneg; n++)
	  
	  /* copy value with proper sign */
	  var[ijk + n*dn] = sym[nv] * var[ijk + (nmax-n)*dn];
	
      } endfor_ijk;
    } endforboxes;
  }

  free(sym);
}





/* set reflection boundary for list of variables */
void set_boundary_symmetry(tL *level, tVarList *varlist) 
{
  double *xp = level->v[Ind("x")];
  double *yp = level->v[Ind("y")];
  double *zp = level->v[Ind("z")];
  
  if (level->shells) {
    set_boundary_shells_symmetry(level, varlist);
    return;
  }

  /* not needed except if there is ghost overreach, which can happen
     for nghost = 4 introduced with BOX */
  if (0) bampi_vlsynchronize(varlist);

  /* octant means three reflections */
  if (Getv("grid", "octant")) {
    /* since we first copy all in x direction, then y, then z, we have
       in the end also the correct data for the edges and corners */
    set_boundary_reflection_one(level, varlist, xp, level->dx, 1);
    set_boundary_reflection_one(level, varlist, yp, level->dy, 3);
    set_boundary_reflection_one(level, varlist, zp, level->dz, 5);
  }

  /* bitant and quadrant require the reflection about z = 0 */
  if (Getv("grid", "bitant") || Getv("grid", "quadrant") ||
      Getv("grid", "qreflect"))
    set_boundary_reflection_one(level, varlist, zp, level->dz, 5);

  if (Getv("grid", "qreflect"))
    set_boundary_reflection_one(level, varlist, yp, level->dy, 3);
    
  /* rotant and quadrant perform an inversion about the z-axis */
  if (Getv("grid", "rotant") || Getv("grid", "quadrant")) 
    set_boundary_inversion(level, varlist);

  /* tube boundary, has to be set before periodic */
  set_boundary_tube(level, varlist);

  /* plane boundary, has to be set before periodic */
  set_boundary_plane(level, varlist);

  /* setting periodic boundaries, note that periodic and reflection symmetry
     should not be used together for now */
  set_boundary_periodic(level, varlist);
}




/* flag symmetry boundary
   prevents application of physical boundary condition
   this will be overridden for ghost zones 
*/
void set_boundary_flags_symmetry(tL *level)
{
  double *xp = level->v[Ind("x")];
  double *yp = level->v[Ind("y")];
  double *zp = level->v[Ind("z")];
  double gx = -level->dx/4;
  double gy = -level->dy/4;
  double gz = -level->dz/4;
  double hx = -gx;
  double hy = -gy;
  double hz = -gz;
  int i;

  if (level->shells)
    return;
  
  if (Getv("grid", "octant")) {
    forallpoints(level, i) {
      if (xp[i] < gx || yp[i] < gy || zp[i] < gz)
	boundaryflag(level, i) = SYMBOUND;
    }
  }
  
  if (Getv("grid", "bitant") || Getv("grid", "quadrant") || 
      Getv("grid", "qreflect")) {
    forallpoints(level, i) {
      if (zp[i] < gz)
	boundaryflag(level, i) = SYMBOUND;
    }
  }

  if (Getv("grid", "rotant") || Getv("grid", "quadrant") ||
      Getv("grid", "qreflect")) {
    forallpoints(level, i) {
      if (yp[i] < gy)
	boundaryflag(level, i) = SYMBOUND;
    }
  }

  if (0) printlevel(level);
}














