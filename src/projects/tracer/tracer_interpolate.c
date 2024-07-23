#include "bam.h"
#include "tracer.h"

#define MAXCUBESIZE 12


/* collect data in local cube */
void tracer_interpolate_filldatacube(tL *level, int vi, 
                              int uui[MAXCUBESIZE][MAXCUBESIZE][MAXCUBESIZE], 
                              double uuu[MAXCUBESIZE][MAXCUBESIZE][MAXCUBESIZE], 
                              int fillcubesize)
{
  double *v = level->v[vi];
  int i, j, k;
  
  for (k = 0; k < fillcubesize; k++)
    for (j = 0; j < fillcubesize; j++)
      for (i = 0; i < fillcubesize; i++)
        uuu[i][j][k] = v[uui[i][j][k]];
}


/* find indices of local cube given lower corner 
   return 0 if cube does not fit into one of the boxes
   return 1 if successful
*/
int tracer_interpolate_fillindexcube(tL *level, int icorner, 
                              int uui[MAXCUBESIZE][MAXCUBESIZE][MAXCUBESIZE], 
                              int fillcubesize)
{
  tB *box;
  int i, j, k;
  const int a = fillcubesize;

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


/* general function in order to interpolate with certain order and scheme at one point */
double tracer_interpolate_TriN(double x, double y, double z, 
                        double xmin, double ymin, double zmin,
                        double dx, double dy, double dz,
                        double uuu[MAXCUBESIZE][MAXCUBESIZE][MAXCUBESIZE], int N)
{

  double v[MAXCUBESIZE][MAXCUBESIZE], w[MAXCUBESIZE];
  double c[MAXCUBESIZE];
  int i, j,k;
  double sum = PI;
  
  /* optimized lagrange interpolation, with minimal computaion of coefficiants */
  coefficients_lagrange_N(N, z, zmin, dz, c);
  for (i = 0; i < N; i++)
  for (j = 0; j < N; j++) {
    v[i][j]=0;
    for (k=0; k<N; k++)
      v[i][j] += c[k]*uuu[i][j][k];
  }
 
  coefficients_lagrange_N(N, y, ymin, dy, c);
  for (i = 0; i < N; i++) {
    w[i]=0;
    for (k=0; k<N; k++)
      w[i] += c[k]*v[i][k];
  }

  coefficients_lagrange_N(N, x, xmin, dx, c);
  sum=0;
  for (k=0; k<N; k++)
      sum += c[k]*w[k];
  
  return sum;
}

int tracer_interpolate_xyz_local(tL *level, double x, double y, double z,
                          int iv, int order, double *result)
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
  double vinterp;

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

  /* look for the origin on this processor */
  i = find_one_point_box(level, xc, yc, zc);

  /* point was not found on this processor */
  if (i == -1)
    return 0;

  /* if point found */
  flag = tracer_interpolate_fillindexcube(level, i, uui,order);

  /* we can't interpolate if the cube does not exist on this processor */
  if (!flag)
    return 0;

  /* we are ready to interpolate */
  tracer_interpolate_filldatacube(level, iv, uui,uuu,order);
  vinterp = tracer_interpolate_TriN(x, y, z, xc, yc, zc, dx, dy, dz,
                                uuu,order);

  *result = vinterp;

  return 1;
}

