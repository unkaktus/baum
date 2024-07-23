/* periodic.c */
/* Nina Jansen 04/03, parallelized BB 06/03 */
/* $Id: periodic.c 1205 2003-10-31 17:59:10Z bruegman $ */
/* Periodic boundaries. 
   You can elect periodic boundaries in x,y and z directions.
*/

#include "bam.h"
#include "boundary.h"


/* for now, this only works on unigrid, or if outer boundary cells
   always have the same size. This can be fixed by setting d to the
   smallest cell-size in the grid such that e.g. xp[ccc] < xmin -
   0.1*dmin. The point is that C has problem when determining when
   doubles are equal, therefor the tolerances must be introduced. The
   function dequal used elsewhere in bam which is not standard c does
   not always work either. */

/* BB: it doesn't? dequal is supposed to work in this context 
   I'll define it again for experimentation
*/
#undef dequaleps
#undef dless
#undef dequal

/* snap effect for grid coordinates, assumes d = level->dx, dy, dz */
#define dequaleps (1e-9)
#define dless(a,b) ((a)<(b)-dequaleps)
#define dlequal(a,b) ((a)<(b)+dequaleps)
#define dequal(a,b) (dlequal(a,b)&&dlequal(b,a))




/* set periodic boundary condition 
   special: O(N) operations by looping over each face separately,
            since then the points are ordered and one search sweep suffices

   note: not parallelized
         number of ghostpoints is hardwired to 2, but the new version
           is easy to generalize
	 ignores other boundary flags for now

   BUG: wrong values along diagonal
*/
void old2_set_boundary_periodic(tL *level, tVarList *varlist) 
{
  double *var;
  int nv, nvars = varlist->n;
  int *ivar = varlist->index; 
  double *xp = level->v[Ind("x")];
  double *yp = level->v[Ind("y")];
  double *zp = level->v[Ind("z")];
  double dx = level->dx; 
  double dy = level->dy; 
  double dz = level->dz; 
  double xmin = level->bbox[0] + 1.5*dx;
  double xmax = level->bbox[1] - 1.5*dx;
  double Dx   = xmax - xmin;
  double ymin = level->bbox[2] + 1.5*dy;
  double ymax = level->bbox[3] - 1.5*dy;
  double Dy   = ymax - ymin;
  double zmin = level->bbox[4] + 1.5*dz;
  double zmax = level->bbox[5] - 1.5*dz;
  double Dz   = zmax - zmin;
  double xc, yc, zc;
  int i, j;
  int n = level->nnodes;

  if (Getv("grid", "periodicx")) {
    for (i = j = 0; i < n; i++) {
      if (xp[i] < xmin) {
	xc = xp[i] + Dx;
	while (!dequal(xc,xp[j]) && ++j < n);
	for (nv = 0; nv < nvars; nv++) {
	  var = level->v[ivar[nv]];
	  var[i] = var[j];
	}
      }
    }
    for (i = j = 0; i < n; i++) {
      if (xp[i] > xmax) {
	xc = xp[i] - Dx;
	while (!dequal(xc,xp[j]) && ++j < n);
	for (nv = 0; nv < nvars; nv++) {
	  var = level->v[ivar[nv]];
	  var[i] = var[j];
	}
      }
    }
  }

  if (Getv("grid", "periodicy")) {
    for (i = j = 0; i < n; i++) {
      if (yp[i] < ymin) {
	yc = yp[i] + Dy;
	while (!dequal(yc,yp[j]) && ++j < n);
	for (nv = 0; nv < nvars; nv++) {
	  var = level->v[ivar[nv]];
	  var[i] = var[j];
	}
      }
    }
    for (i = j = 0; i < n; i++) {
      if (yp[i] > ymax) {
	yc = yp[i] - Dy;
	while (!dequal(yc,yp[j]) && ++j < n);
	for (nv = 0; nv < nvars; nv++) {
	  var = level->v[ivar[nv]];
	  var[i] = var[j];
	}
      }
    }
  }

  if (Getv("grid", "periodicz")) {
    for (i = j = 0; i < n; i++) {
      if (zp[i] < zmin) {
	zc = zp[i] + Dz;
	while (!dequal(zc,zp[j]) && ++j < n);
	for (nv = 0; nv < nvars; nv++) {
	  var = level->v[ivar[nv]];
	  var[i] = var[j];
	}
      }
    }
    for (i = j = 0; i < n; i++) {
      if (zp[i] > zmax) {
	zc = zp[i] - Dz;
	while (!dequal(zc,zp[j]) && ++j < n);
	for (nv = 0; nv < nvars; nv++) {
	  var = level->v[ivar[nv]];
	  var[i] = var[j];
	}
      }
    }
  }

}





/* set periodic boundary condition */
/* bad: O(N^2) for N = number of points, say N = 32^3 */
void old_set_boundary_periodic(tL *level, tVarList *varlist) 
{
  double *var;
  int nvars = varlist->n;
  int *ivar = varlist->index; 
  int *pbound;
  double *coords;

  int nv,n,perneigh,flag,ibound,i,j,interp;
  double *xp = level->v[Ind("x")];
  double *yp = level->v[Ind("y")];
  double *zp = level->v[Ind("z")];
  double dx = level->dx; 
  double dy = level->dy; 
  double dz = level->dz; 
  double xmin = level->bbox[0];
  double xmax = level->bbox[1];
  double ymin = level->bbox[2];
  double ymax = level->bbox[3];
  double zmin = level->bbox[4];
  double zmax = level->bbox[5];
  double xc, yc, zc;

  int nbound;

  int index[1];
  int nbuf, nbuffer;
  double *buf = 0, *buffer = 0;

  nbound = 0;
  interp = 0;

  forallpoints(level, i) {
    if (boundaryflag(level,i) == PERBOUND) nbound++;
  }

  pbound = imalloc(nbound);
  coords = malloc(sizeof(double) * 3 * nbound);

  if (!coords) 
    errorexit("Out of memory in periodic.");

  ibound = 0;

  forboundary13(level, PERBOUND) 
    {

      if (Getv("grid","periodicx")) {
	if ( ((xp[ccc] < xmin + 0.1*dx) && (xp[ccc] > xmin - 0.1*dx))) {
	  xc = xmax - 3*dx;
	} else if ( (xp[ccc] < xmax + 0.1*dx) && (xp[ccc] > xmax - 0.1*dx)) {
	  xc = xmin + 3*dx;
	}
	else if ((xp[ccc] > xmin + 0.9*dx) && (xp[ccc] < xmin + 1.1*dx)) {
	  xc = xmax-2*dx;
	}
	else if ((xp[ccc] > xmax - 1.1*dx) && (xp[ccc] < xmax - 0.9*dx)) {
	  xc = xmin + 2*dx;
	}
	else {
	  xc = xp[ccc];
	}
      } 
      else {
	xc = xp[ccc];
      }	

      if (Getv("grid","periodicy")) {
	if ( ((yp[ccc] < ymin + 0.1*dy) && (yp[ccc] > ymin - 0.1*dy))) {
	  yc = ymax - 3*dy;
	} else if ( (yp[ccc] < ymax + 0.1*dy) && (yp[ccc] > ymax - 0.1*dy)) {
	  yc = ymin + 3*dy;
	}
	else if ((yp[ccc] > ymin + 0.9*dy) && (yp[ccc] < ymin + 1.1*dy)) {
	  yc = ymax-2*dy;
	}
	else if ((yp[ccc] > ymax - 1.1*dy) && (yp[ccc] < ymax - 0.9*dy)) {
	  yc = ymin + 2*dy;
	}
	else {
	  yc = yp[ccc];
	}
      }
      else {
	yc = yp[ccc];
      }

      if (Getv("grid","periodicz")) {
	if ( ((zp[ccc] < zmin + 0.1*dz) && (zp[ccc] > zmin - 0.1*dz))) {
	  zc = zmax - 3*dz;
	} else if ( (zp[ccc] < zmax + 0.1*dz) && (zp[ccc] > zmax - 0.1*dz)) {
	  zc = zmin + 3*dz;
	}
	else if ((zp[ccc] > zmin + 0.9*dz) && (zp[ccc] < zmin + 1.1*dz)) {
	  zc = zmax-2*dz;
	}
	else if ((zp[ccc] > zmax - 1.1*dz) && (zp[ccc] < zmax - 0.9*dz)) {
	  zc = zmin + 2*dz;
	}
	else {
	  zc = zp[ccc];
	}
      }
      else {
	zc = zp[ccc];
      }
      
      pbound[ibound] = ccc;
      coords[ibound*3] = xc;
      coords[ibound*3+1] = yc;
      coords[ibound*3+2] = zc;
      
      ibound++;
    }endfor;
    
  /*  for (ibound = 0; ibound < nbound; ibound++)
    printf("%d %d %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\n", 
	   ibound, pbound[ibound],xp[pbound[ibound]],yp[pbound[ibound]],
	   zp[pbound[ibound]],
	   coords[ibound*3],coords[ibound*3+1],coords[ibound*3+2]);
  */
  /* allocate buffer for data on processor 0 */
  if (processor0) {
    nbuffer = nbound;
    buffer = malloc(sizeof(double) * nbuffer);
    if (!buffer)
      errorexit("Out of memory in periodic.");
  }
  
  /* for each variable */
  for (j = 0; j < nvars; j++) { 
    index[0] = ivar[j];
    /* fill buffer with data at given list of coordinates, allocates buf */

    buffernonorderedpoints(level, nbound, coords, 1, index, 
			   &nbuf, &buf, interp);    
    bampi_combinepointbuffers(level, nbuf, buf, buffer, 1); 
    
    if (processor0) {
      var = level->v[ivar[j]];     
      for (ibound = 0; ibound < nbound; ibound++) {
	if (insidelocalbbox(level,xp[pbound[ibound]],yp[pbound[ibound]],
			    zp[pbound[ibound]])) {
	  var[pbound[ibound]] = buffer[ibound];
	}
      }
    }
    bampi_barrier();
  }
  /* clean up */
  free(buffer);
  free(buf);
  free(coords); 
}




/* flag periodic boundary
   prevents application of physical boundary condition
   this will be overridden for ghost zones 
*/
void set_boundary_flags_periodic(tL *level)
{
  double *xp = level->v[Ind("x")];
  double *yp = level->v[Ind("y")];
  double *zp = level->v[Ind("z")];
  double dx = level->dx; 
  double dy = level->dy; 
  double dz = level->dz; 
  double xmin = level->bbox[0];
  double xmax = level->bbox[1];
  double ymin = level->bbox[2];
  double ymax = level->bbox[3];
  double zmin = level->bbox[4];
  double zmax = level->bbox[5];

  int i;

  if (Getv("grid", "periodicx")) {
    forallpoints(level, i) {
      if ( ((xp[i] < xmin + 0.1*dx) && (xp[i] > xmin - 0.1*dx)) ||
	   ((xp[i] < xmax + 0.1*dx) && (xp[i] > xmax - 0.1*dx)) ||
	   ((xp[i] > xmin + 0.9*dx) && (xp[i] < xmin + 1.1*dx)) ||
	   ((xp[i] > xmax - 1.1*dx) && (xp[i] < xmax - 0.9*dx)) ){
	boundaryflag(level, i) = PERBOUND;
      }
    }
  }
  
  if (Getv("grid", "periodicy")) {
    forallpoints(level, i) {
      if ( ((yp[i] < ymin + 0.1*dy) && (yp[i] > ymin - 0.1*dy)) ||
	   ((yp[i] < ymax + 0.1*dy) && (yp[i] > ymax - 0.1*dy)) ||
	   ((yp[i] > ymin + 0.9*dy) && (yp[i] < ymin + 1.1*dy)) ||
	   ((yp[i] > ymax - 1.1*dy) && (yp[i] < ymax - 0.9*dy)) ){
	boundaryflag(level, i) = PERBOUND;
      }
    }
  }
  
  if (Getv("grid", "periodicz")) {
    forallpoints(level, i) {
      if ( ((zp[i] < zmin + 0.1*dz) && (zp[i] > zmin - 0.1*dz)) ||
	   ((zp[i] < zmax + 0.1*dz) && (zp[i] > zmax - 0.1*dz)) ||
	   ((zp[i] > zmin + 0.9*dz) && (zp[i] < zmin + 1.1*dz)) ||
	   ((zp[i] > zmax - 1.1*dz) && (zp[i] < zmax - 0.9*dz)) ){
	boundaryflag(level, i) = PERBOUND;
      }
    }
  }  
}   





/***************************************************************************/
/* parallel version */

/* fill buffer with data at points specified by a bounding box */
void boxfillbuffer(tL *level, tVarList *vl, double *bbox, 
		   double *buffer, int n)
{
  double p[3];
  double *xp = Ptr(level, "x");
  double *yp = Ptr(level, "y");
  double *zp = Ptr(level, "z");
  int i, j = 0;
  int nv, nvars = vl->n;
  int *ivar = vl->index; 

  forallpoints(level, i) {
    p[0] = xp[i];
    p[1] = yp[i];
    p[2] = zp[i];
    if (insidebbox(bbox, p)) {
      for (nv = 0; nv < nvars; nv++) {
	buffer[j++] = level->v[ivar[nv]][i];
	if (j > n) errorexit("boxfillbuffer(): buffer size mismatch"); 
      }
    }
  }
}





/* read data from buffer at points specified by a bounding box */
void boxreadbuffer(tL *level, tVarList *vl, double *bbox, 
		   double *buffer, int n)
{
  double p[3];
  double *xp = Ptr(level, "x");
  double *yp = Ptr(level, "y");
  double *zp = Ptr(level, "z");
  int i, j = 0;
  int nv, nvars = vl->n;
  int *ivar = vl->index; 

  forallpoints(level, i) {
    p[0] = xp[i];
    p[1] = yp[i];
    p[2] = zp[i];
    if (insidebbox(bbox, p)) {
      for (nv = 0; nv < nvars; nv++) {
        level->v[ivar[nv]][i] = buffer[j++];
	if (j > n) errorexit("boxreadbuffer(): buffer size mismatch"); 
      }
    }
  }
}





/* set periodic boundary in one direction */
void set_boundary_periodic_1d(tL *level, tVarList *vl, int rank, int rank2, 
			      double *boxsend, double *boxrecv, 
			      int npoints)
{
  int n = npoints * vl->n;
  double *bufsend = malloc(sizeof(double) * n * 2);
  double *bufrecv = bufsend + n;

  if (!bufsend) errorexit("set_boundary_periodic_1d(): out of memory");
  if (0) printf("swap: n %d, rank %d, rank2 %d  ", n, rank, rank2);
  
  boxfillbuffer(level, vl, boxsend, bufsend, n);

  bampi_isendrecv(bufsend, bufrecv, n, rank2);

  boxreadbuffer(level, vl, boxrecv, bufrecv, n);

  free(bufsend);
}




/* special version if there is only one proc in this direction */
void set_boundary_periodic_1d_oneproc(tL *level, tVarList *vl, 
			      double *boxsendlo, double *boxrecvlo, 
			      double *boxsendhi, double *boxrecvhi, 
			      int npoints)
{
  int n = npoints * vl->n;
  double *buffer = malloc(sizeof(double) * n);

  if (!buffer)
    errorexit("set_boundary_periodic_1d_oneproc(): out of memory");

  boxfillbuffer(level, vl, boxsendlo, buffer, n);
  boxreadbuffer(level, vl, boxrecvhi, buffer, n);

  boxfillbuffer(level, vl, boxsendhi, buffer, n);
  boxreadbuffer(level, vl, boxrecvlo, buffer, n);

  free(buffer);
}





/* set periodic boundary condition */
void set_boundary_periodic(tL *level, tVarList *varlist) 
{
  int pr = 0;
  char periodic[16], *dirxyz[3] = {"x", "y", "z"};
  double h, o; 
  int c, d, e, i, n, r, s;
  double *bbox = level->bbox;
  double boxsendlo[6], boxrecvlo[6];
  double boxsendhi[6], boxrecvhi[6];
  int rank = level->com->myrank;
  int rank2;
  int *ibbox = level->ibbox;
  int npoints, nallpoints = (ibbox[1]+1)*(ibbox[3]+1)*(ibbox[5]+1);

  if (0) {
    old_set_boundary_periodic(level, varlist);
    return;
  }

  /* for all directions:  c = 0,1,2  d = 0,2,4  e = 1,3,5 */
  for (c = 0; c < 3; c++) {
    d = 2*c;
    e = d+1;
    if (c == 0) h = level->dx;
    if (c == 1) h = level->dy;
    if (c == 2) h = level->dz;
    sprintf(periodic, "periodic%s", dirxyz[c]);

    /* if direction is periodic */
    if (Getv("grid", periodic)) {

      /* check */
      if (Geti("bampi_nghosts") < 2)
	errorexit("set_boundary_periodic: need larger number");

      /* number of buffer points */
      npoints = Geti("bampi_nghosts") * nallpoints/(ibbox[e]+1);

      /* initialize bbox to physical bbox plus h/2 in each direction */ 
      for (i = 0; i < 6; i++) 
	boxsendlo[i] = boxrecvlo[i] = boxsendhi[i] = boxrecvhi[i] = 
	  bbox[i] - h/2 + (i%2)*h;

      /* key step: figure out the four bboxes for periodic boundary */
      o = bbox[d] - h/2;
      boxrecvlo[d] = o;
      boxrecvlo[e] = o + 2*h;
      boxsendlo[d] = o + 2*h;
      boxsendlo[e] = o + 4*h;
      o = bbox[e] + h/2;
      boxrecvhi[e] = o;
      boxrecvhi[d] = o - 2*h;
      boxsendhi[e] = o - 2*h;
      boxsendhi[d] = o - 4*h;
      
      /* determine rank information for parallelization */
      n = level->com->sizexyz[c];
      r = level->com->myrankxyz[c];
      s = n - r - 1;
      if (c == 0) 
	rank2 = rank + s - r;
      else if (c == 1)
	rank2 = rank + level->com->sizexyz[0] * (s - r);
      else
	rank2 = rank + level->com->sizexyz[0]*level->com->sizexyz[1] * (s - r);
      
      /* debug info */
      if (pr && r == 0) {
	printf("lo: n %d  r %d  s %d  ranks %d,%d\n", n, r, s, rank, rank2);  
	printbbox(level, boxrecvlo, 0);
	printbbox(level, boxsendlo, 0);
      }
      if (pr && s == 0) {
	printf("hi: n %d  r %d  s %d  ranks %d,%d\n", n, r, s, rank, rank2);  
	printbbox(level, boxrecvhi, 0);
	printbbox(level, boxsendhi, 0);
      }
      
      /* if there is only one processor in this direction */
      if (n == 1) {
	set_boundary_periodic_1d_oneproc(level, varlist, 
					 boxsendlo, boxrecvlo,
					 boxsendhi, boxrecvhi, npoints);
      } 

      /* else if this processor is at low end */
      else if (r == 0) {
	set_boundary_periodic_1d(level, varlist, rank, rank2, 
				 boxsendlo, boxrecvlo, npoints);
      }
      
      /* else if this processor is at high end */
      else if (s == 0) {
	set_boundary_periodic_1d(level, varlist, rank, rank2, 
				 boxsendhi, boxrecvhi, npoints);
      }

      /* to be save, wait until every processor gets here */
      bampi_barrier();
    }
  }
}





