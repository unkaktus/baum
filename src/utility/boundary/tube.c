/* periodic.c */
/* Nina Jansen 04/03 */
/* $Id: tube.c 1204 2003-10-31 13:20:05Z bruegman $ tube
boundary. Not parallelized (but I don't think there is a need to
parallelize since all neighbours exist on the same processor, and
that's all we're using. */

#include "bam.h"
#include "boundary.h"

#undef dequaleps
#undef dless
#undef dequal

/* snap effect for grid coordinates, assumes d = level->dx etc */
#define dequaleps (0.1*d)
#define dless(a,b) ((a)<(b)-dequaleps)
#define dlequal(a,b) ((a)<(b)+dequaleps)
#define dequal(a,b) (dlequal(a,b)&&dlequal(b,a))

void set_boundary_tube(tL *level, tVarList *varlist) 
{
  double *var;
  int nvars = varlist->n;
  int *ivar = varlist->index; 
  int nv,n,copyneigh,flag;
  double *xp = level->v[Ind("x")];
  double *yp = level->v[Ind("y")];
  double *zp = level->v[Ind("z")];
  double d; 
  double xmin = level->bbox[0];
  double xmax = level->bbox[1];
  double ymin = level->bbox[2];
  double ymax = level->bbox[3];
  double zmin = level->bbox[4];
  double zmax = level->bbox[5];


  if (!Getv("grid", "tubex") && !Getv("grid", "tubez")) return;

  if (Getv("grid", "tubex")) {
    d = level->dx;

    forall27flag(level,PHYBOUND) {
      
      if (dequal(yp[ccc],ymin)&&dequal(zp[ccc],zmin)) {
	copyneigh = cpp;
      } 
      else if (dequal(yp[ccc],ymax)&&dequal(zp[ccc],zmin)) {
	copyneigh = cmp;
      }
      else if (dequal(yp[ccc],ymin)&&dequal(zp[ccc],zmax)) {
	copyneigh = cpm;
      }
      else if (dequal(yp[ccc],ymax)&&dequal(zp[ccc],zmax)) {
	copyneigh = cmm;
      } 
      else if (dequal(yp[ccc],ymin)) {
	copyneigh = cpc;
      }
      else if (dequal(yp[ccc],ymax)) {
	copyneigh = cmc;
      }
      else if (dequal(zp[ccc],zmin)) {
	copyneigh = ccp;
      }
      else if (dequal(zp[ccc],zmax)) {
	copyneigh = ccm;
      }
      else {
	errorexit("something is bungled: grid is a tube, "
                 "but a phybound exists which does not have "
                 "y=ymax or ymin or z=zmax or zmin");
      }
      for (nv = 0; nv < nvars; nv++) {
	var = level->v[ivar[nv]];
	var[ccc] = var[copyneigh];
      }
    } endfor;
  }
  
  if (Getv("grid", "tubez")) {
    d = level->dz;

    forall27flag(level,PHYBOUND) {
      if (dequal(xp[ccc],xmin)&&dequal(yp[ccc],ymin)) {
	copyneigh = ppc;
      } 
      else if (dequal(xp[ccc],xmax)&&dequal(yp[ccc],ymin)) {
	copyneigh = mpc;
      }
      else if (dequal(xp[ccc],xmin)&&dequal(yp[ccc],ymax)) {
	copyneigh = pmc;
      }
      else if (dequal(xp[ccc],xmax)&&dequal(yp[ccc],ymax)) {
	copyneigh = mmc;
      } 
      else if (dequal(xp[ccc],xmin)) {
	copyneigh = pcc;
      }
      else if (dequal(xp[ccc],xmax)) {
	copyneigh = mcc;
      }
      else if (dequal(yp[ccc],ymin)) {
	copyneigh = cpc;
      }
      else if (dequal(yp[ccc],ymax)) {
	copyneigh = cmc;
      }
      else {
	errorexit("something is bungled: grid is a tubez, " 
                 "but a phybound exists which does not have " 
                 "x=xmax or xmin or y=ymax or ymin");
      }
      for (nv = 0; nv < nvars; nv++) {
	var = level->v[ivar[nv]];
	var[ccc] = var[copyneigh];
      }
    } endfor;
  }
}



