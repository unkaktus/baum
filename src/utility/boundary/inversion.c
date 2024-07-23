/* inversion.c */
/* Bernd Bruegmann 10/03 */

/* set symmetry for inversion through axis: "rotant" or "quadrant" 
   this is a non-local operation that may require data from another processor

   for now we fix the z-axis as axis of inversion, and 
           we fix the y=0 plane as plane of symmetry

   explicitly: 
     x -> -x
     y -> -y 
     z -> +z
*/

#include "bam.h"
#include "boundary.h"

void set_boundary_inversion_boxes(tL *level, tVarList *vl);




/* xy-inversion: invert through z-axis along y=0 plane */
void box_xy_inversion(tL *level, tVarList *vl, double *bbox, 
		      double *buffer, int nbuffer)
{
  int pr = 0;
  int ni, nj, nk;
  int i, j, k, m, n, vn;
  int invi, invj, invm, invn;
  double *sym, tmp;

  /* compute sign for inversion */
  sym = dmalloc(vl->n);
  for (vn = 0; vn < vl->n; vn++) {
    sym[vn] = VarSymmetry(vl->index[vn], 0) * VarSymmetry(vl->index[vn], 1);
    if (pr) printf("xy-inversion for %s, sym %3.0lf\n",
		   VarName(vl->index[vn]), sym[vn]);
  }

  /* determine integer size of buffer */
  ni = (bbox[1] - bbox[0])/level->dx + 0.5;
  nj = (bbox[3] - bbox[2])/level->dy + 0.5;
  nk = (bbox[5] - bbox[4])/level->dz + 0.5;

  /* for all points in buffer */
  for (k = 0; k < nk; k++)
  for (j = 0; j < nj; j++)
  for (i = 0; i < ni; i++) {
    n = i + ni*(j + nj*k);
    
    /* invert indices */
    invi = ni - i - 1;
    invj = nj - j - 1;
    invn = invi + ni*(invj + nj*k);

    /* work only on triangle, otherwise swap will happen twice */
    if (invn < n) continue;

    if (pr) printf("n    %2d  i    %d  j    %d  k %d\n", n, i, j, k);
    if (pr) printf("invn %2d  invi %d  invj %d  k %d\n", invn, invi, invj, k);
    
    /* swap data with appropriate sign change, also works on diagonal */
    m = n * vl->n;
    invm = invn * vl->n;
    for (vn = 0; vn < vl->n; vn++, m++, invm++) {
      tmp          = sym[vn] * buffer[m];
      buffer[m]    = sym[vn] * buffer[invm];
      buffer[invm] = tmp;
    }
  }

  /* clean up */
  free(sym);
}




/* set inversion boundary */
void set_boundary_inversion_twoproc(
			 tL *level, tVarList *vl, int rank, int rank2, 
			 double *boxsend, double *boxrecv, 
			 int npoints)
{
  int n = npoints * vl->n;
  double *bufsend = malloc(sizeof(double) * n * 2);
  double *bufrecv = bufsend + n;

  if (!bufsend) errorexit("set_boundary_inversion(): out of memory");
  if (0) printf("swap: n %d, rank %d, rank2 %d  ", n, rank, rank2);
  
  boxfillbuffer(level, vl, boxsend, bufsend, n);
  box_xy_inversion(level, vl, boxsend, bufsend, n);

  bampi_isendrecv(bufsend, bufrecv, n, rank2);
  boxreadbuffer(level, vl, boxrecv, bufrecv, n);

  free(bufsend);
}




/* special version if there is only one proc in this direction */
void set_boundary_inversion_oneproc(tL *level, tVarList *vl, 
				 double *boxsend, double *boxrecv, 
				 int npoints)
{
  int n = npoints * vl->n;
  double *buffer = malloc(sizeof(double) * n);

  if (!buffer)
    errorexit("set_boundary_inversion_oneproc(): out of memory");

  boxfillbuffer(level, vl, boxsend, buffer, n);
  box_xy_inversion(level, vl, boxsend, buffer, n);

  boxreadbuffer(level, vl, boxrecv, buffer, n);

  free(buffer);
}




/* set inversion boundary condition 
   uses box method from periodic.c, keeps some of the arbitrary direction code 
*/
void set_boundary_inversion(tL *level, tVarList *varlist) 
{
  int pr = 0;
  int nghosts = Geti("bampi_nghosts");
  double h, o; 
  int c, d, e, i, n, r, s;
  double *bbox = level->com->bbox;
  int   *ibbox = level->com->ibbox;
  double boxsend[6], boxrecv[6];
  int rank = level->com->myrank;
  int rank2;
  int npoints;

  /* if there is more than one box, we cannot use the faster
     version that relies on symmetric processor arrangement in a single box
  */
  if (level->nboxes > 1) {
    set_boundary_inversion_boxes(level, varlist);
    return;
  }

  /* do nothing if box does not overlap x-z plane
     this works for special two puncture fmr, but not in general 
  */
  if (level->bbox[2] > level->dy/4) return;

  /* check for errors */
  // if (nghosts != 2)
  //   errorexit("set_boundary_inversion: need nghosts = 2, easy to change");
  if (!dequal(level->bbox[1], -level->bbox[0])) {
    printbbox(level, level->bbox, level->ibbox);
    errorexit("set_boundary_inversion: need xmax = -xmin");
  }

  /* initialize some constants */
  h = level->dy;
  npoints = 2 * (ibbox[1]+1)*(ibbox[5]+1)*nghosts;
  c = 0; //   x direction for 0,1,2 indexing
  d = 2; // - y direction for 0,2,4 indexing
  e = 3; // + y direction

  /* initialize bbox to physical bbox plus h/2 in each direction */ 
  for (i = 0; i < 6; i++) 
    boxsend[i] = boxrecv[i] = bbox[i] - h/2 + (i%2)*h;

  /* key step: figure out the two bboxes for inversion boundary */
  o = bbox[d] - h/2;
  boxrecv[d] = o;
  boxrecv[e] = o + nghosts*h;
  boxsend[d] = o + nghosts*h;
  boxsend[e] = o + 2*nghosts*h;
      
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
  if (pr && level->com->myrankxyz[1] == 0) {
    printf("boxrecv/send: n %d  r %d  s %d  ranks %d,%d\n",
	   n, r, s, rank, rank2);  
    printbbox(level, boxrecv, 0);
    printbbox(level, boxsend, 0);
  }
  
  /* if there is only one processor in this direction */
  if (n == 1) {
    set_boundary_inversion_oneproc(level, varlist, 
				   boxsend, boxrecv, npoints);
  } 

  /* else if this processor owns symmetry points */
  else if (level->com->myrankxyz[1] == 0) {
    set_boundary_inversion_twoproc(level, varlist, rank, rank2, 
				   boxsend, boxrecv, npoints);
  }

  /* to be safe, wait until every processor gets here */
  bampi_barrier();
}




/* given buffer with points with rotant symmetry,
   replace it by enlarged buffer with symmetry points added
   called by output routines, for example to turn quadrant into bitant
   - only works for rectangular boxes, but this is all output does right now
*/
void add_rotant_points_to_buffer(tL *level, int vi,
  int *ptr_nbuffer, double **ptr_buffer,
  double *px0, double *py0, double *pz0, int *pnx, int *pny, int *pnz)
{
  tVarList *vl;
  int nx = level->ibbox[1] + 1;
  int ny = level->ibbox[3] + 1;
  int nz = level->ibbox[5] + 1;
  double bbox[6];
  int nbuffer = *ptr_nbuffer;
  double *buffer = *ptr_buffer;
  int i, j, k, m, n, newny, nhalf, nfull;
  double *work;

  /* check whether there is rotant boundary on this level */
  if (!Getv("grid", "rotant") && !Getv("grid", "quadrant")) return; 
  if (level->bbox[2] > 0) return;

  /* change grid info */
  *py0 = -level->bbox[3];
  newny = *pny = 2*(ny - 2);
  nhalf = nbuffer - 2 * nx * nz;
  nfull = 2*nhalf;

  /* store all points at y > 0 in a temporary buffer */
  work = dmalloc(nhalf);
  m = 0;
  for (k = 0; k < nz; k++)
  for (j = 2; j < ny; j++)
  for (i = 0; i < nx; i++)
    work[m++] = buffer[i + nx*(j + ny*k)];

  /* make room for additional points
     reallocation will happen only once if we are called for different vars
  */
  nbuffer = *ptr_nbuffer = nfull;
  buffer  = *ptr_buffer  = realloc(*ptr_buffer, nbuffer * sizeof(double));

  /* store all points at y > 0 in enlarged buffer */
  for (k = 0; k < nz; k++)
  for (j = 0; j < newny/2; j++)
  for (i = 0; i < nx; i++) 
    buffer[i + nx*(j + newny/2 + newny*k)] = work[i + nx*(j + newny/2*k)];

  /* invert the upper half to obtain lower half 
     the available function expects a bounding box and a variable list
  */
  vl = vlalloc(level);
  vlpushone(vl, vi);
  bbox[0] = level->bbox[0] - level->dx/2;
  bbox[1] = level->bbox[1] + level->dx/2;
  bbox[2] = 0;
  bbox[3] = level->bbox[3] + level->dy/2;
  bbox[4] = level->bbox[4] - level->dz/2;
  bbox[5] = level->bbox[5] + level->dz/2;
  box_xy_inversion(level, vl, bbox, work, nhalf); 

  /* store lower half */
  for (k = 0; k < nz; k++)
  for (j = 0; j < newny/2; j++)
  for (i = 0; i < nx; i++)
    buffer[i + nx*(j + newny*k)] = work[i + nx*(j + newny/2*k)];

  /* cleanup */
  free(work);
  vlfree(vl);
}




/* set inversion boundary condition 
   handles more than one box using bampi_getdata, cmp. move_overlap_parallel
   might have to be optimized 
*/
void set_boundary_inversion_boxes(tL *level, tVarList *vl) 
{
  double *x = Ptr(level, "x");
  double *y = Ptr(level, "y");
  double *z = Ptr(level, "z");
  double *coord;
  double *data;
  int *index;
  int i, j, k;
  int npoints;
  int nfound;
  int nvl = vl->n;
  double *sym;

  /* do nothing if boxes do not overlap y=0 plane */
  forallboxes(level) {
    if (box->bbox[2] > level->dy/4) return;
  } endforboxes;

  /* compute sign for inversion */
  sym = dmalloc(nvl);
  for (i = 0; i < nvl; i++)
    sym[i] = VarSymmetry(vl->index[i], 0) * VarSymmetry(vl->index[i], 1);

  /* determine maximal number of points 
     (there may be less for ghost parents) */
  npoints = 0;
  for (i = 0; i < level->nboxes; i++) 
    npoints += level->box[i]->m * level->box[i]->o;
  npoints *= Geti("bampi_nghosts");

  /* storage */
  coord = dmalloc(npoints * 3);
  data  = dmalloc(npoints * (nvl+1));
  index = imalloc(npoints);

  /* create list of all the points that are needed for the symmetry ghosts 
     3d loop, should be replaced by 2d box version
  */
  j = 0;
  forallpoints(level, i) {
    if (y[i] > 0) continue;
    if (j >= 3*npoints) errorexit("found too many symmetry ghosts");
    index[j/3]   = i;
    coord[j++] = - x[i];
    coord[j++] = - y[i];
    coord[j++] = + z[i];
  }
  npoints = j/3;

  /* debug */
  if (0) {
    double bbox[6];
    bbox[0] = bbox[2] = bbox[4] =  DBL_MAX;
    bbox[1] = bbox[3] = bbox[5] = -DBL_MAX;
    for (i = 0; i < j/3; i++) {
      if (coord[3*i+0] < bbox[0]) bbox[0] = coord[3*i+0];
      if (coord[3*i+0] > bbox[1]) bbox[1] = coord[3*i+0];
      if (coord[3*i+1] < bbox[2]) bbox[2] = coord[3*i+1];
      if (coord[3*i+1] > bbox[3]) bbox[3] = coord[3*i+1];
      if (coord[3*i+2] < bbox[4]) bbox[4] = coord[3*i+2];
      if (coord[3*i+2] > bbox[5]) bbox[5] = coord[3*i+2];
    }
    printf("inversion boxes l%d:\nlooking for ", level->l);
    printbbox(level, bbox, 0);
    k = 0;
    for (i = 0; i < level->nboxes; i++) {
      printf("in          ");
      printbbox(level, level->box[i]->bbox, 0);
      if (box_ainb(bbox, level->box[i]->bbox)) k = 1;
    }
    if (!k) { 
      printf("=> this won't work\n");
    }
  }


  /* get data */
  bampi_getdata(level, vl, npoints, coord, data);


  /* store data with correct sign
     if data was not found, this is an error, we need this for the sym!
     nfound is sometimes > 1, say 2 for 2 overlapping ghost parents, divide!
  */
  for (i = 0; i < npoints; i++) {
    k = index[i];
    nfound = data[(nvl+1)*i + nvl];
    if (nfound > 0) {
      for (j = 0; j < nvl; j++)
	level->v[vl->index[j]][k] = sym[j] * data[(nvl+1)*i + j] / nfound;
    } else {
      printf("did not find point %f %f %f\n",
	     coord[3*i], coord[3*i+1], coord[3*i+2]);
      errorexit("set_boundary_inversion_boxes: can't do without it");
    }
  }

  /* clean up */
  free(index);
  free(data);
  free(coord);
  free(sym);
  // should not be needed since we have asked getdata also for the ghosts
  // bampi_vlsynchronize(vl);
}










