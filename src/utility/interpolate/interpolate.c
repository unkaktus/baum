/* interpolate.c */
/* Bernd Bruegmann 3/03, Wolfgang Tichy 1/2004, 9/2010 */
/* Jose Gonzalez 2/06 */


#include "bam.h"
#include "interpolate.h"







/* print datacube */
void interpolate_printdatacube(tL *level,
                               int uui[MAXCUBESIZE][MAXCUBESIZE][MAXCUBESIZE], 
                               double uuu[MAXCUBESIZE][MAXCUBESIZE][MAXCUBESIZE], 
                               int fillcubesize)
{
  int i, j, k;
  
  for (k = 0; k < fillcubesize; k++)
  for (j = 0; j < fillcubesize; j++)
  for (i = 0; i < fillcubesize; i++) {
    printf("uui[%d][%d][%d] = %8d   ", i, j, k, uui[i][j][k]);
    printf("uuu[%d][%d][%d] = %22.15e\n", i, j, k, uuu[i][j][k]);
  }
}


/* collect data in local cube */
void interpolate_filldatacube(tL *level, int vi, 
                              int uui[MAXCUBESIZE][MAXCUBESIZE][MAXCUBESIZE], 
                              double uuu[MAXCUBESIZE][MAXCUBESIZE][MAXCUBESIZE], 
                              int fillcubesize)
{
  double *v = level->v[vi];
  int i, j, k;
  //printf("fillcubesize %d\n", fillcubesize);
  
  for (k = 0; k < fillcubesize; k++)
    for (j = 0; j < fillcubesize; j++)
      for (i = 0; i < fillcubesize; i++)
        uuu[i][j][k] = v[uui[i][j][k]];
}


/* find indices of local cube given lower corner 
   return 0 if cube does not fit into one of the boxes
   return 1 if successful
*/
int interpolate_fillindexcube(tL *level, int icorner, 
                              int uui[MAXCUBESIZE][MAXCUBESIZE][MAXCUBESIZE], 
                              int fillcubesize)
{
  tB *box;
  int i, j, k;
  const int a = fillcubesize;

  /* preliminary */
  /* which box are we in? */
  if (0) {
    if (level->nboxes > 2)
      errorexit("BOX: implement interpolate_fillindexcube");
    box = level->box[0];
    if (level->nboxes == 2 && icorner >= level->box[1]->noffset)
      box = level->box[1];
  }
  
  /* generic solution */
  for (i=0; i<level->nboxes; i++) {
    if (icorner<level->box[i]->noffset)
      break;
    box = level->box[i];
  }

  /* don't leave the box */
  i = iofbox(box, icorner);
  j = jofbox(box, icorner);
  k = kofbox(box, icorner);
  if (i < 0 || j < 0 || k < 0 ||
      i+a > box->m || j+a > box->n || k+a > box->o)
    return 0;

  /* fill index cube */
  for (k = 0; k < a; k++)
    for (j = 0; j < a; j++)
      for (i = 0; i < a; i++)
        uui[i][j][k] = icorner + i*box->di + j*box->dj + k*box->dk;
  return 1;
}


/* find indices of local cube given lower corner 
   return 0 if cube does not fit into one of the boxes
   return 1 if successful
*/
int interpolate_fillindexcube_box(tB *box, int ii, int jj, int kk, 
                                  int uui[MAXCUBESIZE][MAXCUBESIZE][MAXCUBESIZE], 
                                  int fillcubesize)
{
  int i, j, k;
  const int a = fillcubesize;

  /* don't leave the box */
  if (ii < 0 || jj < 0 || kk < 0 ||
      ii+a > box->m || jj+a > box->n || kk+a > box->o)
    return 0;

  /* fill index cube */
  for (k = 0; k < a; k++)
    for (j = 0; j < a; j++)
      for (i = 0; i < a; i++)
        uui[i][j][k] = ijkofbox(box, ii+i, jj+j, kk+k);
  return 1;
}


/* general function in order to interpolate with certain order and scheme at one point */
double interpolate_TriN(double x, double y, double z, 
                        double xmin, double ymin, double zmin,
                        double dx, double dy, double dz,
                        double uuu[MAXCUBESIZE][MAXCUBESIZE][MAXCUBESIZE], int N, int scheme)
{

  double v[MAXCUBESIZE][MAXCUBESIZE], w[MAXCUBESIZE];
  double c[MAXCUBESIZE];
  int i, j;
  double sum = PI;
  
  /* optimized lagrange interpolation, with minimal computaion of coefficiants */
  if (scheme==LAGRANGE) {
    
    coefficients_lagrange_N(N, z, zmin, dz, c);
    for (i = 0; i < N; i++)
    for (j = 0; j < N; j++) {
      v[i][j] = interpolate_lagrange_N(N, z, zmin, dz, c, uuu[i][j]);
    }
 
    coefficients_lagrange_N(N, y, ymin, dy, c);
    for (i = 0; i < N; i++) {
      w[i] = interpolate_lagrange_N(N, y, ymin, dy, c, v[i]);
    }

    coefficients_lagrange_N(N, x, xmin, dx, c);
    sum = interpolate_lagrange_N(N, x, xmin, dx, c, w); 
    
  /* WENO interpolation */
  } else if (scheme==WENO) {
    
    for (i = 0; i < N; i++)
    for (j = 0; j < N; j++) {
      v[i][j] = interpolate_WENO_N(N, z, zmin, dz, c, uuu[i][j]);
    }
  
    for (i = 0; i < N; i++) {
      w[i] = interpolate_WENO_N(N, y, ymin, dy, c, v[i]);
    }
  
    sum = interpolate_WENO_N(N, x, xmin, dx, c, w); 
  }
  
  return sum;
}

/* wrapper for TriN used by multigrid */
void interpolate_TriN_varlist(tL *level0, tL *level1,
                              tVarList *u0, tVarList *u1, tB *box0, tB *box1,
                              int i, int j, int k, int order)
{
  int a = order;
  int b = a/2 - 1;
  int flag;
  int ii, jj, kk, ijk, j0, j1, n;
  double dx = level0->dx;
  double dy = level0->dy;
  double dz = level0->dz;
  double Dx,Dy,Dz;
  double uuu[MAXCUBESIZE][MAXCUBESIZE][MAXCUBESIZE];
  int uui[MAXCUBESIZE][MAXCUBESIZE][MAXCUBESIZE];
  
  
  /* find node in box0 that snaps to node i,j,k in box1 
  note that R and P are symmetric in this regard
  */
  box_point_indices(box0, box1, &ii, &jj, &kk, i, j, k);
  
  /* find origin of interpolation cube */
  ii -= b;
  jj -= b;
  kk -= b;
  
  /* fill index cube */
  flag = interpolate_fillindexcube_box(box0, ii, jj, kk, uui,a);
  if (!flag) {
    printf("point %d %d %d,  box 0 %d  0 %d  0 %d\n", 
           ii, jj, kk, box0->imax, box0->jmax, box0->kmax);
    errorexit("interpolation cube doesn't fit into box");
  }
  
  
  /* offset between origin and interpolation point */
  Dx = xofbox(box1, i) - xofbox(box0, ii);
  Dy = yofbox(box1, j) - yofbox(box0, jj);
  Dz = zofbox(box1, k) - zofbox(box0, kk);
  
  /* linear index where we want to store the result */
  ijk = ijkofbox(box1, i, j, k);
  /* for all variables */
  for (n = 0; n < u0->n; n++) {
    j0 = u0->index[n];
    j1 = u1->index[n];
    if (level0->v[j0] == 0 || level1->v[j1] == 0) continue;
    
    /* fill data cube */
    interpolate_filldatacube(level0, j0, uui,uuu,a);
    
    /* interpolate */
    level1->v[j1][ijk] = 
        interpolate_TriN(Dx, Dy, Dz, 0., 0., 0., dx, dy, dz, 
                         uuu,order,LAGRANGE);
  }
  
}



















/* check whether this processor might own a point for interpolation 
   based on its bounding box

   unused, it turned out to be simpler to let everyone interpolate
   who finds a cube and then to average the (identical) result
*/
int insidelocalbbox(tL *level, double x, double y, double z)
{
  double dx = level->dx/2;
  double dy = level->dy/2;
  double dz = level->dz/2;
  double xmin = level->com->bbox[0] - dx;
  double ymin = level->com->bbox[2] - dy;
  double zmin = level->com->bbox[4] - dz;
  double xmax = level->com->bbox[1] + dx;
  double ymax = level->com->bbox[3] + dy;
  double zmax = level->com->bbox[5] + dz;

  if (dless(x,xmin) || dless(y,ymin) || dless(z,zmin))
    return 0;

  if (dless(x,xmax) && dless(y,ymax) && dless(z,zmax))
    return 1;

  return 0;
}



/* map point into available grid region by symmetry
   set flag for those coordinates where the mapping actually happened
*/
void map_xyz_withsym(tL* l, double *x, double *y, double *z,
		     int *fx, int *fy, int *fz)
{
  tG* g = l->grid;
  
  /* remember if we used the symmetry to transform the point */
  *fx = *fy = *fz = 0;

  /* we do not have to check anything if we are in full mode 
     leave x, y, z unchanged   */
  if(g->full) return;

  /* Do nothing if we are in the shells. If shells ever use syms we have to
     fix this!!! NOTE: in shells: *x,*y,*z are not Cartesian coords. */
  if(l->shells) return;

  /* use symmetries on x, y, z */
  if(g->bitant)
  {
    if(*z<0.0) { *fz = 1; *z = -*z; }
  }
  else if(g->rotant)
  {
    if(*y<0.0) { *fx = *fy = 1; *y = -*y;  *x = -*x; }
  }
  else if(g->quadrant)
  {
    if(*z<0.0) { *fz = 1; *z = -*z; }
    if(*y<0.0) { *fx = *fy = 1; *y = -*y;  *x = -*x; }
  }
  else if(g->qreflect)
  {
    if(*z<0.0) { *fz = 1; *z = -*z; }
    if(*y<0.0) { *fy = 1; *y = -*y; }
  }
  else if(g->octant)
  {
    if(*z<0.0) { *fz = 1; *z = -*z; }
    if(*y<0.0) { *fy = 1; *y = -*y; }
    if(*x<0.0) { *fx = 1; *x = -*x; }
  }
  else 
    errorexit("unknown grid mode");
}

/* check if point is within a bbox, use level lev to figure out syms */
int xyzinsidebbox_withsym(tL* lev, 
                          double *bbox, double x, double y, double z)
{
  double coord[3];
  int fx, fy, fz;

  /* map point by symmetry as needed */
  map_xyz_withsym(lev, &x, &y, &z, &fx, &fy, &fz);
        
  coord[0] = x;
  coord[1] = y;
  coord[2] = z;
  return insidebbox(bbox, coord);
}

/* check whether interpolation cube fits on this processor */
int check_interpolation_cube_local_withsym(tL *level, 
                                           double x, double y, double z, int order)
{
  double xc, yc, zc;
  double xm = level->bbox[0];
  double ym = level->bbox[2];
  double zm = level->bbox[4];
  double dx = level->dx;
  double dy = level->dy;
  double dz = level->dz;
  int f;              // not used here
  int m = order - 1;  // order 4 and 6:  3 and 5
  int n = m/2;        // order 4 and 6:  1 and 2

  /* use symmetry to remap points */
  map_xyz_withsym(level, &x, &y, &z, &f, &f, &f);

  /* determine origin of cube 
     using the global bbox should guard against certain round-off problems */
  xc = xm + dx * (floor((x - xm + dequaleps)/dx) - n);
  yc = ym + dy * (floor((y - ym + dequaleps)/dy) - n);
  zc = zm + dz * (floor((z - zm + dequaleps)/dz) - n);


  /* look for lower and upper corner both in one box */
  int i=0;
  forallboxes(level) {
    if (xyzinsidebbox(box->com->bbox, xc, yc, zc) &&
	xyzinsidebbox(box->com->bbox, xc+m*dx, yc+m*dy, zc+m*dz))
      
      /* interpolation cube does fit */
      return 1;

  } endforboxes;

  /* interpolation cube does not fit */
  return 0;
}



/* check whether interpolation cube fits on this processor */
int check_interpolation_cubeinbox_local(tB *Box,
                                        double x, double y, double z, int order)
{
 tL* level=Box->level;
  double xc, yc, zc;
  double xm = Box->bbox[0];
  double ym = Box->bbox[2];
  double zm = Box->bbox[4];
  double dx = level->dx;
  double dy = level->dy;
  double dz = level->dz;
  int m = order - 1;  // order 4 and 6:  3 and 5
  int n = m/2;        // order 4 and 6:  1 and 2

  /* determine origin of cube 
  using the global bbox should guard against certain round-off problems */
  xc = xm + dx * (floor((x - xm + dequaleps)/dx) - n);
  yc = ym + dy * (floor((y - ym + dequaleps)/dy) - n);
  zc = zm + dz * (floor((z - zm + dequaleps)/dz) - n);
  /*
  printf("  %2.2f %2.2f %2.2f    %2.2f %2.2f %2.2f   %2.2f %2.2f %2.2f\n", x,y,z,xc,yc,zc,xc+m*dx, yc+m*dy, zc+m*dz);
  printf("  %2.2f %2.2f %2.2f %2.2f %2.2f %2.2f\n", Box->bbox[0],Box->bbox[1],Box->bbox[2],Box->bbox[3],Box->bbox[4],Box->bbox[5]);
  printf("  %2.2f %2.2f %2.2f %2.2f %2.2f %2.2f\n\n", Box->com->bbox[0],Box->com->bbox[1],Box->com->bbox[2], Box->com->bbox[3],Box->com->bbox[4],Box->com->bbox[5]);
*/
  /* look for lower and upper corner both in one box */
  if (xyzinsidebbox(Box->com->bbox, xc, yc, zc) &&
      xyzinsidebbox(Box->com->bbox, xc+m*dx, yc+m*dy, zc+m*dz)) {
      /* interpolation cube does fit */
      return 1;
      
      }
  /* interpolation cube does not fit */
  return 0;
}



/* interpolate scalar to one given point */
double interpolate_xyz_scalar(tL *level, double x, double y, double z, 
			      int vi, int order, int scheme)
{
  double *xp = Ptr(level, "x");
  double *yp = Ptr(level, "y");
  double *zp = Ptr(level, "z");
  double xc, yc, zc;
  double xm = level->com->bbox[0];
  double ym = level->com->bbox[2];
  double zm = level->com->bbox[4];
  double dx = level->dx;
  double dy = level->dy;
  double dz = level->dz;
  double vxyz;
  int flag;
  int pr = 0;
  int m = order - 1;    // order 4 and 6:  3 and 5
  int n = m/2;          // order 4 and 6:  1 and 2
  double uuu[MAXCUBESIZE][MAXCUBESIZE][MAXCUBESIZE];
  int uui[MAXCUBESIZE][MAXCUBESIZE][MAXCUBESIZE];

  /* determine origin of cube */
  xc = xm + dx * (floor((x - xm + dequaleps)/dx) - n);
  yc = ym + dy * (floor((y - ym + dequaleps)/dy) - n);
  zc = zm + dz * (floor((z - zm + dequaleps)/dz) - n);
  if (pr) printf("%7.3f %7.3f %7.3f  falls into cube  "
		 "[%.3f,%.3f]x[%.3f,%.3f]x[%.3f,%.3f]\n", 
		 x, y, z, xc, xc+m*dx, yc, yc+m*dy, zc, zc+m*dz);

  /* look for this point on this processor */
  for (n = 0; n < level->nnodes; n++)
    if (dequal(xc, xp[n]) && dequal(yc, yp[n]) && dequal(zc, zp[n]))
      break;

  /* if point found */
  if (n != level->nnodes)
    flag = interpolate_fillindexcube(level, n, uui,order);

  /* we can't interpolate if the cube does not exist on this processor */
  if (n == level->nnodes || !flag) {
    flag = 0;
    vxyz = 0;
  }

  /* we are ready to interpolate */
  else {
    flag = 1;
    interpolate_filldatacube(level, vi, uui,uuu,order);
    vxyz = interpolate_TriN(x, y, z, xc, yc, zc, dx, dy, dz, 
                            uuu,order,scheme);
  }

  /* each processor wants this data 
     tricky bit is that more than one processor may have the data if
     the point falls into ghostzones 
  */
  if (bampi_size() > 1) {
    double local[2], global[2];
    
    local[0] = flag;
    local[1] = vxyz;
    bampi_allreduce_sum_vector(local, global, 2);
    
    if (pr) printf("Processor %d is%s interpolating.\n", 
		   bampi_rank(), (flag?"":" not"));
    if (global[0] == 0) {
      if (pr) printf("No processor has interpolated.\n");
      vxyz = 0;
    } else {
      if (global[0] == 1) {
	if (pr) printf("One processor has interpolated.\n");
      } else 
	if (pr) printf("More than one processor has interpolated.\n");
      vxyz = global[1]/global[0];
    }
  }

  /* print info */
  if (0) {
    printf("The value of %s at %7.3f %7.3f %7.3f   is   %e.\n", 
	   VarName(vi), x, y, z, vxyz);
  }

  /* return value */
  return vxyz;
}



/* Interpolate scalar to one given point, using grid symmetries.
   This is a wrapper for interpolate_xyz_scalar */
double interpolate_xyz_scalar_withSym(tL *level, double x, double y, double z,
                                      int vi, int order, int scheme)
{
  double val;

  /* use symmetries on x, y, z */
  if(level->grid->bitant)
  {
    if(z<0.0)   z = -z;
  }
  if(level->grid->quadrant)
  {
    if(z<0.0)   z = -z;
    if(y<0.0) { y = -y;  x = -x; }
  }
  if(level->grid->qreflect)
  {
    if(z<0.0)   z = -z;
    if(y<0.0)   y = -y;
  }
  if(level->grid->octant)
  {
    if(z<0.0)   z = -z;
    if(y<0.0)   y = -y;
    if(x<0.0)   x = -x;
  }

  val=interpolate_xyz_scalar(level, x, y, z, vi, order, scheme);
  return val;
}



/* interpolate list of variables to a given point, processor local 
   called by bufferorderedpoints in amr/points.c
*/
int interpolate_xyz_local(tL *level, double x, double y, double z, 
                          int nv, int *iv, double *vinterp, 
                          int start, int order, int scheme)
{
  double xc, yc, zc;
  double xm = level->com->bbox[0];
  double ym = level->com->bbox[2];
  double zm = level->com->bbox[4];
  double xmin = level->bbox[0];
  double xmax = level->bbox[1];
  double ymin = level->bbox[2];
  double ymax = level->bbox[3];
  double zmin = level->bbox[4];
  double zmax = level->bbox[5];
  double dx = level->dx;
  double dy = level->dy;
  double dz = level->dz;
  int flag, i;
  int m = order - 1;    // order 4 and 6:  3 and 5
  int n = m/2;          // order 4 and 6:  1 and 2
  double uuu[MAXCUBESIZE][MAXCUBESIZE][MAXCUBESIZE];
  int uui[MAXCUBESIZE][MAXCUBESIZE][MAXCUBESIZE];

  if (0) printf("calling interp for %d vars, %d %s to %d %s\n", 
		nv, iv[0], VarName(iv[0]), iv[nv-1], VarName(iv[nv-1]));

  /* determine origin of cube */
  xc = xm + dx * (floor((x - xm + dequaleps)/dx) - n);
  yc = ym + dy * (floor((y - ym + dequaleps)/dy) - n);
  zc = zm + dz * (floor((z - zm + dequaleps)/dz) - n);

  /* useful for rectangular top-level boxes: 
     move cube off-center at outer boundary */
  if (dless(xc,xmin)) xc = xmin;
  if (dless(yc,ymin)) yc = ymin;
  if (dless(zc,zmin)) zc = zmin;
  if (dless(xmax,xc+m*dx)) xc = xmax-m*dx;
  if (dless(ymax,yc+m*dy)) yc = ymax-m*dy;
  if (dless(zmax,zc+m*dz)) zc = zmax-m*dz;

  if (0) printf("%7.3f %7.3f %7.3f  falls into cube  "
		"[%.3f,%.3f]x[%.3f,%.3f]x[%.3f,%.3f]\n", 
		x, y, z, xc, xc+m*dx, yc, yc+m*dy, zc, zc+m*dz);

  /* look for the origin on this processor */
  i = find_one_point_box(level, xc, yc, zc);

  /* debug */
  if (1 && i != -1) {
    double *xp = Ptr(level, "x");
    double *yp = Ptr(level, "y");
    double *zp = Ptr(level, "z");
    if (!(dequal(xc, xp[i]) && dequal(yc, yp[i]) && dequal(zc, zp[i]))) {
      printf("find_one_point_box returned invalid point %d\n", i);
      printf("%f %f   %f %f   %f %f\n", xc, xp[i], yc, yp[i], zc, zp[i]);
      errorexit("this is bad");
    }
  }
  if (0 && i == -1)
    printf("find_one_point_box did not find %f %f %f\n", xc, yc, zc);

  /* point was not found on this processor */
  if (i == -1) return -1;

  /* if point found */
  flag = interpolate_fillindexcube(level, i, uui,order);

  /* we can't interpolate if the cube does not exist on this processor */
  if (!flag) return -1;

  /* we are ready to interpolate */
  for (n = 0; n < nv; n++) {
    interpolate_filldatacube(level, iv[n], uui,uuu,order);
    vinterp[n] = interpolate_TriN(x, y, z, xc, yc, zc, dx, dy, dz,
                                  uuu,order,scheme);
    if (0) printf("vinterp[%d] = %e\n", n, vinterp[n]);
  }

  return i;
}



/* interpolate list of variables to a given point, processor local 
   version with only minimal checks and without shift for top level
   return 1 if successful, 0 if not successful
*/
int interpolate_xyz_local_minimal(tL *level, double x, double y, double z, 
                                  int nv, int *iv, double *vinterp, 
                                  int order, int scheme)
{
  double xc, yc, zc;
  double xm = level->bbox[0];
  double ym = level->bbox[2];
  double zm = level->bbox[4];
  double dx = level->dx;
  double dy = level->dy;
  double dz = level->dz;
  int m = order - 1;  // order 4 and 6:  3 and 5
  int n = m/2;        // order 4 and 6:  1 and 2
  int i;
  double uuu[MAXCUBESIZE][MAXCUBESIZE][MAXCUBESIZE];
  int uui[MAXCUBESIZE][MAXCUBESIZE][MAXCUBESIZE];

  if (0) printf("calling interp for %d vars, %d %s to %d %s\n", 
		nv, iv[0], VarName(iv[0]), iv[nv-1], VarName(iv[nv-1]));

  /* determine origin of cube 
     using the global bbox should guard against certain round-off problems */
  xc = xm + dx * (floor((x - xm + dequaleps)/dx) - n);
  yc = ym + dy * (floor((y - ym + dequaleps)/dy) - n);
  zc = zm + dz * (floor((z - zm + dequaleps)/dz) - n);

  /* look for point (loops over all boxes) */
  i = find_one_point_box(level, xc, yc, zc);
  if (i == -1) return 0;
  
  /* fill index cube, may fail if upper corner of cube does not fit */
  i = interpolate_fillindexcube(level, i, uui,order);
  if (i == 0) return 0;

  /* interpolate each variable */
  for (i = 0; i < nv; i++) {
    interpolate_filldatacube(level, iv[i], uui,uuu,order);
    vinterp[i] = interpolate_TriN(x, y, z, xc, yc, zc, dx, dy, dz,
                                  uuu,order,scheme);
  }

  /* if we got here, it worked */
  return 1;
}



/* interpolate list of variables to a given point, processor local 
   version with only minimal checks and without shift for top level
   return 1 if successful, 0 if not successful
   use only for shells because there we have 6 times the same coordinates
*/
int interpolate_xyz_localinbox_minimal(tB *box, double x, double y, double z, 
                                  int nv, int *iv, double *vinterp, 
                                  int order, int scheme)
{
  tL* level=box->level;
  double xm = box->bbox[0];
  double ym = box->bbox[2];
  double zm = box->bbox[4];
  double dx = level->dx;
  double dy = level->dy;
  double dz = level->dz;
  int m = order - 1;  // order 4 and 6:  3 and 5
  int n = m/2;        // order 4 and 6:  1 and 2
  int l;
  double uuu[MAXCUBESIZE][MAXCUBESIZE][MAXCUBESIZE];
  int uui[MAXCUBESIZE][MAXCUBESIZE][MAXCUBESIZE];
  
  if (0) printf("calling interp for %d vars, %d %s to %d %s\n", 
      nv, iv[0], VarName(iv[0]), iv[nv-1], VarName(iv[nv-1]));

  /* determine origin of cube using the global bbox should guard against 
     certain round-off problems */
  double xc = xm + dx * (floor((x - xm + dequaleps)/dx) - n);
  double yc = ym + dy * (floor((y - ym + dequaleps)/dy) - n);
  double zc = zm + dz * (floor((z - zm + dequaleps)/dz) - n);
  
  /* fill index cube, may fail if upper corner of cube does not fit */
  int i = floor((xc - box->com->bbox[0] + dequaleps)/dx);
  int j = floor((yc - box->com->bbox[2] + dequaleps)/dy);
  int k = floor((zc - box->com->bbox[4] + dequaleps)/dz);
  if (!interpolate_fillindexcube_box(box, i,j,k, uui,order))
    return 0;

  /* interpolate each variable */
  for (l = 0; l < nv; l++) {
    interpolate_filldatacube(level, iv[l], uui,uuu,order);
    vinterp[l] = interpolate_TriN(x, y, z, xc, yc, zc, dx, dy, dz, 
                                  uuu,order,scheme);
  }

  /* if we got here, it worked */
  return 1;

}



/* interpolate list of variables to a given point, processor local 
   version with only minimal checks and without shift for top level
   return 1 if successful, 0 if not successful
*/
int interpolate_xyz_local_minimal_withsym(tL *level, double x, double y, double z, 
                                          int nv, int *iv, double *vinterp, 
                                          int order, int scheme)
{
  double val;
  int i, j, flag;
  int fx, fy, fz;

  /* map point by symmetry as needed */
  map_xyz_withsym(level, &x, &y, &z, &fx, &fy, &fz);

  /* call standard interpolator with transformed point */
  flag=interpolate_xyz_local_minimal(level, x, y, z, nv, iv, vinterp, order,scheme);

  if(flag){

    /* correct interpolate variable with proper sign */
    if (fx || fy || fz)
      for (i = 0; i < nv; i++) {

	val=1.0;

	if (fx) val *= VarSymmetry(iv[i], 0);
	if (fy) val *= VarSymmetry(iv[i], 1);
	if (fz) val *= VarSymmetry(iv[i], 2);

	vinterp[i]*=val;

	if (0) {
	  printf("i = %d  fx %d fy %d fz %d\n",i,fx,fy,fz);
	  for (j = 0; j < 3; j++) 
	    printf("sign for %s, direction %d:  sym %d\n",
		   VarName(iv[i]), j, VarSymmetry(iv[i], j));

	  printf("val = %f\n",val);

	}
      }
  }
  
  return flag;
}



/* interpolate list of variables to a given point, processor local 
   version with only minimal checks and without shift for top level
   return 1 if successful, 0 if not successful
*/
int interpolate_xyz_localinbox_minimal_withsym(tB *box, double x, double y, double z, 
                                               int nv, int *iv, double *vinterp, 
                                               int order, int scheme)
{
  double val;
  int i, j, flag;
  int fx, fy, fz;

  /* map point by symmetry as needed */
  map_xyz_withsym(box->level, &x, &y, &z, &fx, &fy, &fz);

  /* call standard interpolator with transformed point */
  flag=interpolate_xyz_localinbox_minimal(box, x, y, z, nv, iv, vinterp, order,scheme);

  if(flag){

    /* correct interpolate variable with proper sign */
    if (fx || fy || fz)
      for (i = 0; i < nv; i++) {

      val=1.0;

      if (fx) val *= VarSymmetry(iv[i], 0);
      if (fy) val *= VarSymmetry(iv[i], 1);
      if (fz) val *= VarSymmetry(iv[i], 2);

      vinterp[i]*=val;

      if (0) {
        printf("i = %d  fx %d fy %d fz %d\n",i,fx,fy,fz);
        for (j = 0; j < 3; j++) 
          printf("sign for %s, direction %d:  sym %d\n",
                 VarName(iv[i]), j, VarSymmetry(iv[i], j));

        printf("val = %f\n",val);

      }
      }
  }
  
  return flag;
}








/* interpolate scalar to one given point remembering the last point used */
double interpolate_xyz_scalar_remember(tL *level, double x, double y, double z,
				       int vi, int order, int scheme)
{
  double *xp = Ptr(level, "x");
  double *yp = Ptr(level, "y");
  double *zp = Ptr(level, "z");
  double xc, yc, zc;
  double xm = level->com->bbox[0];
  double ym = level->com->bbox[2];
  double zm = level->com->bbox[4];
  double dx = level->dx;
  double dy = level->dy;
  double dz = level->dz;
  double vxyz;
  int flag, n;
  int pr = 0;
  static int last_n = 0;
  int m = order - 1;    // order 4 and 6:  3 and 5
  int o = m/2;          // order 4 and 6:  1 and 2
  double uuu[MAXCUBESIZE][MAXCUBESIZE][MAXCUBESIZE];
  int uui[MAXCUBESIZE][MAXCUBESIZE][MAXCUBESIZE];

  //printf("last_n = %d\n",last_n);

  /* determine origin of cube */
  xc = xm + dx * (floor((x - xm + dequaleps)/dx) - o);
  yc = ym + dy * (floor((y - ym + dequaleps)/dy) - o);
  zc = zm + dz * (floor((z - zm + dequaleps)/dz) - o);

  if (pr) printf("%7.3f %7.3f %7.3f  falls into cube  "
		 "[%.3f,%.3f]x[%.3f,%.3f]x[%.3f,%.3f]\n", 
		 x, y, z, xc, xc+m*dx, yc, yc+m*dy, zc, zc+m*dz);

  /* look for this point on this processor */
  for (n = last_n; n < level->nnodes; n++)
    if (dequal(xc, xp[n]) && dequal(yc, yp[n]) && dequal(zc, zp[n]))
      break;

  if (n == level->nnodes){
    for (n = last_n; n >= 0; n--)
    if (dequal(xc, xp[n]) && dequal(yc, yp[n]) && dequal(zc, zp[n]))
      break;
    if (n==-1) n=level->nnodes;
  }

  /* if point found */
  if (n != level->nnodes)
    flag = interpolate_fillindexcube(level, n, uui,order);

  /* we can't interpolate if the cube does not exist on this processor */
  if (n == level->nnodes || !flag) {
    flag = 0;
    vxyz = 0;
    last_n = 0;
  }

  /* we are ready to interpolate */
  else {
    flag = 1;
    interpolate_filldatacube(level, vi,  uui,uuu,order);
    vxyz = interpolate_TriN(x, y, z, xc, yc, zc, dx, dy, dz, 
                            uuu,order,scheme);
    last_n = n;
  }


  /* each processor wants this data 
     tricky bit is that more than one processor may have the data if
     the point falls into ghostzones 
  */
  if (bampi_size() > 1) {
    double local[2], global[2];
    
    local[0] = flag;
    local[1] = vxyz;
    bampi_allreduce_sum_vector(local, global, 2);
    
    if (pr) printf("Processor %d is%s interpolating.\n", 
		   bampi_rank(), (flag?"":" not"));
    if (global[0] == 0) {
      if (pr) printf("No processor has interpolated.\n");
      vxyz = 0;
    } else {
      if (global[0] == 1) {
	if (pr) printf("One processor has interpolated.\n");
      } else 
	if (pr) printf("More than one processor has interpolated.\n");
      vxyz = global[1]/global[0];
    }
  }

  /* print info */
  if (0) {
    printf("The value of %s at %7.3f %7.3f %7.3f   is   %e.\n", 
	   VarName(vi), x, y, z, vxyz);
  }

  /* return value */
  return vxyz;
}




/* interpolate scalar to one given point using a box grid */
double interpolate_xyz_scalar_box(tL *level, double x, double y, double z, 
                                  int vi, int order, int scheme)
{
  double xc, yc, zc;
  double xmin = level->com->bbox[0];
  double xmax = level->com->bbox[1];
  double ymin = level->com->bbox[2];
  double ymax = level->com->bbox[3];
  double zmin = level->com->bbox[4];
  double zmax = level->com->bbox[5];
  double dx = level->dx;
  double dy = level->dy;
  double dz = level->dz;
  double vxyz;
  int nx = level->com->ibbox[1]+1;
  int ny = level->com->ibbox[3]+1;
  int nz = level->com->ibbox[5]+1;
  int flag, n,i,j,k;
  int pr = 0;
  int m = order - 1;    // order 4 and 6:  3 and 5
  int o = m/2;          // order 4 and 6:  1 and 2
  double uuu[MAXCUBESIZE][MAXCUBESIZE][MAXCUBESIZE];
  int uui[MAXCUBESIZE][MAXCUBESIZE][MAXCUBESIZE];
  tB *fbox;

  /* get min and max of the correct box !!! */
  fbox = box_containing_xyz(level, x, y, z);
  xmin = fbox->com->bbox[0];
  xmax = fbox->com->bbox[1];
  ymin = fbox->com->bbox[2];
  ymax = fbox->com->bbox[3];
  zmin = fbox->com->bbox[4];
  zmax = fbox->com->bbox[5];

  if(pr)
  {
    printf("x,y,z=%g,%g,%g\n", x,y,z);
    printf("[xmin,xmax]x[ymin,ymax]x[zmin,zmax]=[%g,%g]x[%g,%g]x[%g,%g]\n",
           xmin,xmax, ymin,ymax, zmin,zmax);
  }

  /* determine origin of cube */
  xc = xmin + dx * (floor((x - xmin + dequaleps)/dx) - o);
  yc = ymin + dy * (floor((y - ymin + dequaleps)/dy) - o);
  zc = zmin + dz * (floor((z - zmin + dequaleps)/dz) - o);
  if(pr) printf("m=%d o=%d  xc,yc,zc=%g,%g,%g\n", m,o, xc,yc,zc);

  if(pr) printf("%g %g %g  falls into cube  "
                "[%g,%g]x[%g,%g]x[%g,%g]\n", 
                x, y, z, xc, xc+m*dx, yc, yc+m*dy, zc, zc+m*dz);

  /* if ((xc>=xmin)&&(xc<=xmax-m*dx)&&(yc>=ymin)&&(yc<=ymax-m*dy)&&(zc>=zmin)&&(zc<=zmax-m*dz)){ */
  if( dgreatereq(xc,xmin) && dlesseq(xc,xmax-m*dx) && 
      dgreatereq(yc,ymin) && dlesseq(yc,ymax-m*dy) &&  
      dgreatereq(zc,zmin) && dlesseq(zc,zmax-m*dz) )
  {
    if(pr) printf("point found in proccessor %d\n",bampi_rank());

    i = (xc-xmin+dequaleps)/dx;
    j = (yc-ymin+dequaleps)/dy;
    k = (zc-zmin+dequaleps)/dz;
    /* n = k*nx*ny + j*nx + i;   does not work for second box */
    n = k*nx*ny + j*nx + i + fbox->noffset;
    
    interpolate_fillindexcube(level, n, uui,order);
    flag = 1;
    interpolate_filldatacube(level, vi, uui,uuu,order);
    vxyz = interpolate_TriN(x, y, z, xc, yc, zc, dx, dy, dz, 
                            uuu,order,scheme);
  }else{

    if(pr) printf("point NOT founded in proccessor %d\n",bampi_rank());
    flag = 0;
    vxyz = 0;

  }

  /* each processor wants this data 
     tricky bit is that more than one processor may have the data if
     the point falls into ghostzones 
  */
  if (bampi_size() > 1) {
    double local[2], global[2];
    
    local[0] = flag;
    local[1] = vxyz;

    bampi_allreduce_sum_vector(local, global, 2);
    if (pr) printf("Processor %d is%s interpolating.\n", 
		   bampi_rank(), (flag?"":" not"));
    if (global[0] == 0) {
      if(pr)
      { 
        printf("No processor has interpolated.\n");
        errorexit("interpolate_xyz_scalar_box failed!");
      }
      vxyz = 0;
    } else {
      if (global[0] == 1) {
	if (pr) printf("One processor has interpolated.\n");
      } else 
	if (pr) printf("More than one processor has interpolated.\n");
      vxyz = global[1]/global[0];
    }
  }
  else if(flag==0) 
    if(pr)
    { 
      printf("No processor has interpolated.\n");
      errorexit("interpolate_xyz_scalar_box failed!");
    }

  /* print info */
  if (0) {
    printf("The value of %s at %g %g %g   is   %g.\n", 
	   VarName(vi), x, y, z, vxyz);
  }

  /* return value */
  return vxyz;
}



















