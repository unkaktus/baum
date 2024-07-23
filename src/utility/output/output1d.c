/* output1d.c */
/* Bernd Bruegmann 12/02 */

#include "bam.h"
#include "output.h"

#define PR 0




/***************************************************************************/
/* 1d output */


/* make list of ordered points for a given straight line based on global 
   bounding box
*/
void makeorderedlist1d(tL *level, int snapflag, 
		       double dx, double dy, double dz, 
		       int *npoints, double **coords)
{
  double x0 = Getd("outputx0");
  double y0 = Getd("outputy0");
  double z0 = Getd("outputz0");
  double a, b, coord[3];
  double *bbox = level->bbox;
  int nmax = level->ibbox[1] + level->ibbox[3] + level->ibbox[5];
  int i, j, m, n;
  double dxtmp = dx;

  /* choose whether point location should snap to next larger grid point */
  if (snapflag) {
    double dx = level->dx;
    double dy = level->dy;
    double dz = level->dz;
    double xmin = level->bbox[0];
    double ymin = level->bbox[2];
    double zmin = level->bbox[4];

    x0 = xmin + floor((x0-xmin)/dx - 0.001)*dx + dx;
    y0 = ymin + floor((y0-ymin)/dy - 0.001)*dy + dy;
    z0 = zmin + floor((z0-zmin)/dz - 0.001)*dz + dz;
  }

  /* this is one of several boxes
     for now just handle special case of coordinate lines if box is off-center
  */
  if (!xyzinsidebbox(bbox, x0, y0, z0)) {
    *npoints = 0;
    *coords = 0;
    
    /* ignore diagonals if off-center */
    if (dx && dy && dz) return;
    
    /* consider x, y, z lines */
    if (dx) x0 = bbox[0] + (!snapflag)*dx/2;
    if (dy) y0 = bbox[2] + (!snapflag)*dy/2;
    if (dz) z0 = bbox[4] + (!snapflag)*dz/2;
    if (!xyzinsidebbox(bbox, x0, y0, z0)) return;
  }
 
  /* for negative slope */
  if (dx < 0) dxtmp = - dx;


  /* walk line backwards */
  for (m = 0; m > -nmax; m--) {
    if (dless(x0 + m*dxtmp, bbox[0])) break;
    if (dless(y0 + m*dy, bbox[2])) break;
    if (dless(z0 + m*dz, bbox[4])) break;
  }
  m++;

  /* walk line forwards */
  for (n = 0; n < nmax; n++) {
    if (dless(bbox[1], x0 + n*dxtmp)) break;
    if (dless(bbox[3], y0 + n*dy)) break;
    if (dless(bbox[5], z0 + n*dz)) break;
  }
  n--;

  if (0) printf("makeorderedlist1d: l%d xyz0 %6.3f %6.3f %6.3f "
		"dxyz %6.3f %6.3f %6.3f range %d %d\n", 
		level->l, x0, y0, z0, dx, dy, dz, m, n);

  /* now we know the number of points */
  *npoints = n-m+1;
  if (*npoints <= 0) {
    *npoints = 0;
    return;
  }
  *coords = malloc(sizeof(double) * 3 * (*npoints));

  /* fill in coordinates */
  for (i = m, j = 0; i <= n; i++) {
    
   (*coords)[j++] = x0 + i*dx;
   (*coords)[j++] = y0 + i*dy;
   (*coords)[j++] = z0 + i*dz;

  
  }
  
}

/* create list of points through puncture positions */
/* only works for 2 black holes */
/* assumes center of mass frame */
int makeconnectinglist1d(tL *level, int snapflag, 
                       double dx, double dy, double dz, 
                       int *npoints, double **coords, int line)
{
  double x0 = Getd("outputx0");
  double y0 = Getd("outputy0");
  double z0 = Getd("outputz0");

  double a, b, coord[3];
  double *bbox = level->bbox;
  int l = line;
  double slope = 0;
  double invxdiff;
  double dynew;
  double xp1 = level->grid->puncpos[0][0];
  double yp1 = level->grid->puncpos[0][1];
  double xp2, yp2;
  double ri, divx;
  int m, n, i, j;
  int nmax = level->ibbox[1] + level->ibbox[3] + level->ibbox[5];

   /* grid symmetries */
  if (Getv("grid", "rotant") || Getv("grid", "octant")){
     printf("\nError: rotant and octant symmetry not handled by makeconnectinglist1d!\n");
     return -1;
   } else {
     xp2 = level->grid->puncpos[1][0];
     yp2 = level->grid->puncpos[1][1];
   }

   /* if the punctures are too close (e.g. after merger), do nothing */
   if (dequal(xp1, xp2) && dequal(yp1, yp2)) return -1;

  /* if we are on y-axis, use coordinate axis as before*/
  if (dequal(xp1, 0)) {
    makeorderedlist1d(level, snapflag,  0, dy,  0, npoints, coords);  
    /* calculate coordinates correctly */
    l = 2;
  }
  /* x-axis */
  else if (dequal(yp1, 0)){
    makeorderedlist1d(level, snapflag, dx, 0, 0, npoints, coords);
    l = 1; 
  }
  else {
    /* slope and vertical offset*/
    invxdiff = 1/(xp2 - xp1);

    slope = (yp2 - yp1) * invxdiff;
    if (line == 7) slope = 1/slope;

    if (1) 
      printf("\nslope = %f, line = %d\n", slope, line);


  /* choose whether point location should snap to next larger grid point */
  if (snapflag) {
    double dx = level->dx;
    double dy = level->dy;
    double dz = level->dz;
    double xmin = level->bbox[0];
    double ymin = level->bbox[2];
    double zmin = level->bbox[4];

    x0 = xmin + floor((x0-xmin)/dx - 0.001)*dx + dx;
    y0 = ymin + floor((y0-ymin)/dy - 0.001)*dy + dy;
    z0 = zmin + floor((z0-zmin)/dz - 0.001)*dz + dz;
  }

  /* do some of the checks of makeorderedlist1d */
  /* this is one of several boxes
     for now just handle special case of coordinate lines if box is off-center
  */
  if (!xyzinsidebbox(bbox, x0, y0, z0)) {

    *npoints = 0;
    *coords = 0;
    
    /* ignore diagonals if off-center */
    if (dx && dy && dz) return -1;
    
    /* consider x, y, z lines */
    if (dx) x0 = bbox[0] + (!snapflag)*dx/2;
    if (dy) y0 = bbox[2] + (!snapflag)*dy/2;
    if (dz) z0 = bbox[4] + (!snapflag)*dz/2;
    if (!xyzinsidebbox(bbox, x0, y0, z0)) return -1;
  }

    /* want as many points as on x- and y-axes to make sure no points come and go during simulation
       we use dx or dy from level information to count the points and give the radius of the circles
       all the points should lie on circles contained in the box with spacing
       given by dx 	 */

  /* count how many points fit on maximal radius according to dx; want as many
	points on this line as on coordinate axes 
  */

  /* walk x-axis backwards */
  for (m = 0; m > -nmax; m--) {
    if (dless(x0 + m*(level->dx), bbox[0])) break;
  }
  m++;
  /* walk x-axis forwards */
  for (n = 0; n < nmax; n++) {
    if (dless(bbox[1], x0 + n*(level->dx))) break;
  }
  n--;
  
  dynew = slope * dx;

  /* go to work */
  /* assume no vertical offset! */
 

  /* now we know the number of points */
  *npoints = n-m+1;

  if (*npoints <= 0) {
    *npoints = 0;
    return -1;
  }
  *coords = malloc(sizeof(double) * 3 * (*npoints));

  /* fill in coordinates */
  /* x_i^2 = r_i^2/(1 + (dynew/dx)^2)
     y_i = x_i * slope */
  for (i = m, j = 0; i <= n; i++) {
    ri = i*dx;  
    divx = 1/sqrt(1 + dynew*dynew/(dx*dx));
    
    (*coords)[j++] = (x0 + ((ri)*(divx))); 
	if (0) printf("coordx = %f\n", (*coords)[j-1]);

    
    (*coords)[j++] = y0 + slope*((*coords)[j-1] - x0);
	if (0) printf("coordy = %f\n", (*coords)[j-1]);
    
    if (Getv("grid", "quadrant")){
      if (((*coords)[j-2])<0 && ((*coords)[j-1])<0) 
      {
        (*coords)[j-2] = - (*coords)[j-2];
        (*coords)[j-1] = - (*coords)[j-1];
      }
    }
    (*coords)[j++] = z0; 
	if (0) printf("coordz = %f\n", (*coords)[j-1]);
    }
  }
 
  return l;
}

int makeperpconnectinglist1d(tL *level, int snapflag, 
                             double dx, double dy, double dz, 
                             int *npoints, double **coords, int line)
{

  int l = line;
  double xp1 = level->grid->puncpos[0][0];
  double yp1 = level->grid->puncpos[0][1];

  /* if punctures are on y-axis, use coordinate axis x as before*/
  if (dequal(xp1, 0)) {
    makeorderedlist1d(level, snapflag, dx, 0, 0, npoints, coords);
    /* calculate coordinates correctly */
    l = 1;
  }
  /* x-axis */
  else if (dequal(yp1, 0)){
    makeorderedlist1d(level, snapflag, 0, dy, 0, npoints, coords);
    l = 2;
  } else {

    makeconnectinglist1d(level, snapflag, dy, dx, 0, npoints, coords, 7);
    l = 7;

  }


  return l;

}




/* write 1d 
   uses generic but slow method to generate list of points 
   that is then searched for in the level
*/
void write_level_1d(tL *level, int nv, int *iv, int line, char *suffix)
{
  if (PR) printf("write_level_1d  line %d\n",line);
  
  int pr = 0;
  int i, j;
  FILE *fp;
  char filename[1000];
  double dx = level->dx;
  double dy = level->dy;
  double dz = level->dz;
  double *coords;
  int npoints, coordoff;
  int nbuf, nbuffer;
  double b, *buf = 0, *buffer = 0;
  int interp = Getv("1doutinterpolate", "yes");

  if (nv <= 0) return;
  i = !interp;
  //printf("level %d, %s\n", level->l, VarName(*iv));


  /* make list of points for this 1d line based on global bounding box */
  if (line == 0) 
    makeorderedlist1d(level, i, dx, dy, dz, &npoints, &coords);
  if (line == 1)
    makeorderedlist1d(level, i, dx,  0,  0, &npoints, &coords);
  if (line == 2)
    makeorderedlist1d(level, i,  0, dy,  0, &npoints, &coords);
  if (line == 3)
    makeorderedlist1d(level, i,  0,  0, dz, &npoints, &coords);
 
  /* create diagonal in xy-plane and another one perpendicular to the former*/
  if (line == 4)
    makeorderedlist1d(level, i, dx, dy,  0, &npoints, &coords);
  if (line == 5)
    makeorderedlist1d(level, i, dx,  0, dz, &npoints, &coords);
  if (line == 6)
    makeorderedlist1d(level, i,  0, dy, dz, &npoints, &coords);
  
  if (line == 7)
    makeorderedlist1d(level, i, -dx, dy, 0, &npoints, &coords);

  

  if (pr) 
    for (i = 0; i < npoints; i++)
      printf("%d %6.3f %6.3f %6.3f\n", 
	     i, coords[3*i], coords[3*i+1], coords[3*i+2]);

  /* the line may not intersect the points on this level */
  if (npoints == 0) return;

  /* fill buffer with data at given list of coordinates, allocates buf */
  bufferorderedpoints(level, npoints, coords, nv, iv, &nbuf, &buf, interp, output_order,output_scheme);
  if (pr) {
    printf("npoints = %d, nbuf = %d\n", npoints, nbuf);
    for (i = 0; i < nbuf; i++) {
      printf(" %10.3e", buf[i]);
      if ((i+1) % (nv+1) == 0) printf("\n");
    }
  }

  /* allocate buffer for data on processor 0 */
  if (processor0) {
    nbuffer = npoints * nv;
    buffer = malloc(sizeof(double) * nbuffer);

    /* points that are not found are masked with a special value */
    for (i = 0; i < nbuffer; i++)
      buffer[i] = -1; // NODATA;
  }

  /* combine data from different processors */
  bampi_combinepointbuffers(level, nbuf, buf, buffer, nv);

  /* processor 0 does the writing */
  if (processor0) {
    
    /* compute coordinates for the diagonals and set offset */
    if (line == 0 || line == 4 || line == 5 || line == 6 || line == 7) {

      for (i = 0; i < 3*npoints; i += 3)
	coords[i] = sqrt(coords[i]  *coords[i] + 
			 coords[i+1]*coords[i+1] +
			 coords[i+2]*coords[i+2]) * (coords[i]<0?-1:1);
      coordoff = 0;
      
   } else
      coordoff = line-1;

    /* for each variable */
    for (j = 0; j < nv; j++) {

      /* filenames */
      snprintf(filename, 1000, "%s/%s.%s%d%s", 
               Gets("outdir_1d"), VarName(iv[j]), suffix, level->l, boxstring);
      
      /* open files */
      fp = fopen(filename, "a");
      if (!fp) errorexits("failed opening %s", filename);
    
      /* write info */
      fprintf(fp, "\"Time = %f\"\n", level->time);
      
      /* write */
      for (i = 0; i < npoints; i++) {
	b = buffer[i*nv+j];
	if (1 || b != NODATA)
	  fprintf(fp, "%22.15e  %22.15e\n", coords[3*i+coordoff], b);
      }

      /* close file */
      fprintf(fp,"\n");
      fclose(fp);
    }
  }

  /* clean up */
  if (processor0) free(buffer);
  free(buf);
  free(coords);

  /* wait, could be optimized across different writes */
  bampi_barrier();
}




/* write 1d of every level */
void write_grid_1d(tL *level, int nv, int *iv, int line, char *suffix)
{
  if (PR) printf("write_grid_1d  line %d\n",line);
  
  int pr = 0;
  int i, j;
  FILE *fp;
  char filename[1000];
  double dx = level->dx;
  double dy = level->dy;
  double dz = level->dz;
  double *coords;
  int npoints, coordoff;
  int nbuf, nbuffer;
  double b, *buf = 0, *buffer = 0;
  int interp = Getv("1doutinterpolate", "yes");

  if (nv <= 0) return;
  i = !interp;
  //printf("level %d, %s\n", level->l, VarName(*iv));


  /* make list of points for this 1d line based on global bounding box */
  if (line == 0) 
    makeorderedlist1d(level, i, dx, dy, dz, &npoints, &coords);
  if (line == 1)
    makeorderedlist1d(level, i, dx,  0,  0, &npoints, &coords);
  if (line == 2)
    makeorderedlist1d(level, i,  0, dy,  0, &npoints, &coords);
  if (line == 3)
    makeorderedlist1d(level, i,  0,  0, dz, &npoints, &coords);
 
  /* create diagonal in xy-plane and another one perpendicular to the former*/
  if (line == 4)
    makeorderedlist1d(level, i, dx, dy,  0, &npoints, &coords);
  if (line == 5)
    makeorderedlist1d(level, i, dx,  0, dz, &npoints, &coords);
  if (line == 6)
    makeorderedlist1d(level, i,  0, dy, dz, &npoints, &coords);
  
  if (line == 7)
    makeorderedlist1d(level, i, -dx, dy, 0, &npoints, &coords);

  

  if (pr) 
    for (i = 0; i < npoints; i++)
      printf("%d %6.3f %6.3f %6.3f\n", 
             i, coords[3*i], coords[3*i+1], coords[3*i+2]);

  /* the line may not intersect the points on this level */
  if (npoints == 0) return;

  /* fill buffer with data at given list of coordinates, allocates buf */
  bufferorderedpoints(level, npoints, coords, nv, iv, &nbuf, &buf, interp, output_order,output_scheme);
  if (pr) {
    printf("npoints = %d, nbuf = %d\n", npoints, nbuf);
    for (i = 0; i < nbuf; i++) {
      printf(" %10.3e", buf[i]);
      if ((i+1) % (nv+1) == 0) printf("\n");
    }
  }

  /* allocate buffer for data on processor 0 */
  if (processor0) {
    nbuffer = npoints * nv;
    buffer = malloc(sizeof(double) * nbuffer);

    /* points that are not found are masked with a special value */
    for (i = 0; i < nbuffer; i++)
      buffer[i] = -1; // NODATA;
  }

  /* combine data from different processors */
  bampi_combinepointbuffers(level, nbuf, buf, buffer, nv);

  /* processor 0 does the writing */
  if (processor0) {
    
    /* compute coordinates for the diagonals and set offset */
    if (line == 0 || line == 4 || line == 5 || line == 6 || line == 7) {

      for (i = 0; i < 3*npoints; i += 3)
        coords[i] = sqrt(coords[i]  *coords[i] + 
            coords[i+1]*coords[i+1] +
            coords[i+2]*coords[i+2]) * (coords[i]<0?-1:1);
      coordoff = 0;
      
    } else
      coordoff = line-1;

      /* for each variable */
      for (j = 0; j < nv; j++) {

        /* filenames */
        snprintf(filename, 1000, "%s/%s.%s%d%s", 
                 Gets("outdir_1d"), VarName(iv[j]), suffix, level->l, boxstring);
      
        /* open files */
        fp = fopen(filename, "a");
        if (!fp) errorexits("failed opening %s", filename);
    
        /* write info */
        fprintf(fp, "\"Time = %f\"\n", level->time);
      
        /* write */
        for (i = 0; i < npoints; i++) {
          b = buffer[i*nv+j];
          if (1 || b != NODATA)
            fprintf(fp, "%22.15e  %22.15e\n", coords[3*i+coordoff], b);
        }

        /* close file */
        fprintf(fp,"\n");
        fclose(fp);
      }
  }

  /* clean up */
  if (processor0) free(buffer);
  free(buf);
  free(coords);

  /* wait, could be optimized across different writes */
  bampi_barrier();
}











