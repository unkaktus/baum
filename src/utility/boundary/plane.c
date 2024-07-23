/* periodic.c */
/* Nina Jansen 04/03 */
/* $Id: plane.c 1204 2003-10-31 13:20:05Z bruegman $ tube
boundary. Not parallelized (but I don't think there is a need to
parallelize since all neighbours exist on the same processor, and
that's all we're using. */

#include "bam.h"
#include "boundary.h"

void set_boundary_plane(tL *level, tVarList *varlist) 
{
  double *var;
  int nvars = varlist->n;
  int *ivar = varlist->index; 
  int nv,n,copyneigh,flag;
  double *xp = level->v[Ind("x")];
  double *yp = level->v[Ind("y")];
  double *zp = level->v[Ind("z")];
  double dz = level->dz; 
  double xmin = level->bbox[0];
  double xmax = level->bbox[1];
  double ymin = level->bbox[2];
  double ymax = level->bbox[3];
  double zmin = level->bbox[4];
  double zmax = level->bbox[5];
  
  if (!Getv("grid", "plane")) return;

  forall27flag(level,PHYBOUND) {
    if (((zp[ccc] < zmin + 0.1*dz) && (zp[ccc] > zmin - 0.1*dz))) {
      copyneigh = ccp;
    }
    else if (((zp[ccc] < zmax + 0.1*dz) && (zp[ccc] > zmax - 0.1*dz))) {
      copyneigh = ccm;
    }
    else {
      errorexit("something is bungled: grid is a plane, "
                "but a phybound exists which does not have "
                "z=zmax or zmin");
    }
    for (nv = 0; nv < nvars; nv++) {
      var = level->v[ivar[nv]];
      var[ccc] = var[copyneigh];
    }
  } endfor;
  
}



