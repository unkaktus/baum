/* 2d output on a sphere of coordinate radius r*/
/* Sascha Husa 3/06 */
/* Wolfgang Tichy 3/2010 */
/* Roman Gold 11/2010 */
/* mth 11/2011 -> shells support*/

#include "bam.h"
#include "output.h"

#define PR 0

/* sphere output: single proc version (for testing and debugging )*/ 
void write_sphere(tL *level, int ncircles, int npoints, double r, char *suffix, int index){
      
  /* ncircles = number of circles for the sphere (without the poles) -> must be ODD  */
  /* npoints  = number of points for each circle of the sphere       -> must be EVEN */
  /* r        = coordinate radius of sphere                                          */
  /* suffix   = a string used in the file name, e.g. r1/r3/r3                        */
  /* index    = index of the variable that we want to output                         */
  
  char name[201];
  FILE *fp;

  int i, j, l = 0;
  double dtheta, dphi;
  double x, y, z, theta, phi;
  int np_proc, np_rem;

  int np_total = (ncircles+1) * (npoints+1);
  int init_point, end_point;

  double **matrix = dmatrix(0, ncircles+1, 0, npoints+1);


  if (ncircles%2==0){
    errorexit("Call to integral over sphere with even number of circles");
  }
  if ((npoints+1)%2==0){
    errorexit("Call to integral over sphere with odd number of points");
  }

  if (!(Getv("grid", "box")))
    errorexit("Integral over sphere only tested in Box grid");

  if (Getv("grid", "bitant")) {
    /* z=0 reflection symmetry */
    dtheta = PI/(2*(ncircles+1));
    dphi = 2.0*PI/npoints;

  } else if (Getv("grid", "rotant")) {
      /* x=y=0 inversion symmetry */
      dtheta = PI/(ncircles+1);
      dphi = PI/npoints;

  } else if (Getv("grid", "quadrant")) {
    /* z=0 and y=0 reflection symmetry */
    dtheta = PI/(2*(ncircles+1));
    dphi = PI/npoints;

  } else if (Getv("grid", "octant")) {
    /* x=0, y=0, and z=0 reflection symmetry */
    dtheta = PI/(2*(ncircles+1));
    dphi = PI/(2*npoints);
  } else {
    /* Full grid */
    dtheta = PI/(ncircles+1);
    dphi = 2.0*PI/npoints;
  }

  for   (i=1; i <= ncircles + 1; i++){
    for (j=1; j <= npoints  + 1; j++){
	  
      theta =  i    * dtheta;
      phi   = (j-1) * dphi;

      x = r * sin(theta) * cos(phi);
      y = r * sin(theta) * sin(phi);
      z = r * cos(theta);

      errorexit("fixme: outdated call to interpolate");
      //matrix[i][j] = interpolate_xyz_scalar_box(level, x, y, z, index, 0);
    }
  }

  sprintf(name, "%s/%s.%sl%d", Gets("outdir_r"), VarName(index), suffix, level->l);
  fp = fopen(name, "a");

  if (level->iteration == 0) { // write header only for first time level

    fprintf(fp, "# gridfunction = %s\n", VarName(index));
    fprintf(fp, "# r = %14.6e\n", r);
    fprintf(fp, "# grid = {%d, %d}\n", ncircles, npoints);
    fprintf(fp, "# theta = {%14.6e, %14.6e}\n", dtheta, theta);
    fprintf(fp, "# phi   = {%14.6e, %14.6e}\n", 0.0e0,  phi);
    fprintf(fp, "# dtheta = %14.6e\n", dtheta);
    fprintf(fp, "# dphi   = %14.6e\n", dphi);
  }

  fprintf(fp, "# iteration = %d\n",     level->iteration);
  fprintf(fp, "# time      = %14.6e\n", level->time);

  for   (i=1; i <= ncircles + 1; i++){
    for (j=1; j <= npoints  + 1; j++){

      fprintf(fp, "%14.6e", matrix[i][j]);
    }
    fprintf(fp, "\n");
  }
  fprintf(fp, "\n");

  fclose(fp);
  free_dmatrix(matrix, 0, ncircles+1, 0, npoints+1);
}





/* 2d output for a sphere: parallel version */
void write_level_sphere0(tL *level, int numvars, int *indexarray, 
                         int ntheta, int nphi, double r, char *suffix)
{
  /* Function Arguments:

  level        ...  level struct
  numvars      ...  number of gridfunctions for output
  indexarray   ...  array containing indices for output grid functions
  ntheta       ...  number of theta intervals 
  nphi         ...  number of phi intervals
  r            ...  extraction radius
  suffix       ...  file suffix
  */

  int pr = 0;       // debugging switch: 0 is off, 1 is on
  int interp = 1;   // always allow interpolation in interpolation function

  static int ncall = 0;

  int nseries;      // 	   \_  timing info used for file header etc.
  double dtseries;  //     /

  int i, j, k, k2;  //  counters

  FILE *fp;
  char filename[1000];

  double dtheta, dphi, theta, phi;

  double *global_bbox = level->bbox;
  double *local_bbox  = level->com->bbox;

  double *coords, *sph_coords;
  double coord[3];

  int npoints, nlocalpoints;

  int nbuf, nbuffer;
  double *buf = 0, *buffer = 0;

  /* END OF DECLARATION SECTION */

  if (numvars <= 0) return;  // no variables requested for output 

  timer_start(0, "write_level_sphere0");
  ncall++;

  if (pr && level->iteration == 0) printf("extract at radius %e at level %d\n", r, level->l);

  if (r > global_bbox[1] && r > global_bbox[3] && r > global_bbox[5]) {
    if (pr && level->iteration == 0) printf("bbox does not fit at this level\n");
    return;
  }

  npoints = (ntheta+1) * (nphi+1); // total # of points

  nseries = timeforoutput_di_dt(level, Geti("2doutiter"), Getd("2douttime"));
  dtseries = (nseries == 1) ? 0 : level->iteration/(nseries-1)*level->dt;

  if (Getv("grid", "bitant")) {
    /* z=0 reflection symmetry */
    dtheta = PI/(2*(ntheta+1));
    dphi = 2.0*PI/nphi;

  } else if (Getv("grid", "rotant")) {
      /* x=y=0 inversion symmetry */
      dtheta = PI/(ntheta+1);
      dphi = PI/nphi;

  } else if (Getv("grid", "quadrant")) {
    /* z=0 and y=0 reflection symmetry */
    dtheta = PI/(2*(ntheta+1));
    dphi = PI/nphi;

  } else if (Getv("grid", "octant")) {
    /* x=0, y=0, and z=0 reflection symmetry */
    dtheta = PI/(2*(ntheta+1));
    dphi = PI/(2*nphi);
  } else {
    /* Full grid */
    dtheta = PI/(ntheta+1);
    dphi = 2.0*PI/nphi;
  }


  /* count points for this processor */

  if (pr) {
    nlocalpoints = 0;
    for   (i=1; i <= ntheta + 1; i++){
      for (j=1; j <= nphi   + 1; j++){

	theta =  i    * dtheta;
	phi   = (j-1) * dphi;

	coord[0] = r * sin(theta) * cos(phi);
	coord[1] = r * sin(theta) * sin(phi);
	coord[2] = r * cos(theta);

	if (insidebbox(local_bbox, coord)) {
	  nlocalpoints++;
	} 
      } 
    }
  } else {
    nlocalpoints = -1;
  }

  if (pr && level->iteration == 0 && nlocalpoints > -1) {
    printf("found %d points from %d total points at proc %d\n", nlocalpoints, npoints, bampi_rank());
  }

  /* allocate and fill in coordinates */
  coords     = (double *) malloc(sizeof(double) * 3 * npoints);
  sph_coords = (double *) malloc(sizeof(double) * 2 * npoints);

  k  = 0;
  k2 = 0;
  for   (i=1; i <= ntheta + 1; i++){
    for (j=1; j <= nphi   + 1; j++){

      theta =  i    * dtheta;
      phi   = (j-1) * dphi;

      coord[0] = r * sin(theta) * cos(phi);
      coord[1] = r * sin(theta) * sin(phi);
      coord[2] = r * cos(theta);

      sph_coords[k2++] = theta;
      sph_coords[k2++] = phi;    

      coords[k++] = coord[0];
      coords[k++] = coord[1];
      coords[k++] = coord[2];
    }
  }

  if (PR) // tell the world what coord vals we got
    for (i = 0; i < npoints; i++) {
      printf("%d %6.3f %6.3f %6.3f\n", 
	     i, coords[3*i], coords[3*i+1], coords[3*i+2]);
    }

  /* fill buffer with data at given list of coordinates, allocates buf */
  if (level->shells || 1) {
    //FIXME: this is a more general scheme, should be tested and should replace the other function
    bufferpoints_full(level, npoints, coords, numvars, indexarray, &nbuf, &buf, output_order,output_scheme);
  } else {
    bufferorderedpoints(level, npoints, coords, numvars, indexarray, &nbuf, &buf, interp, output_order,output_scheme);
  }
  
  if (PR) printf("nbuf computed from bufferorderedpoints = %d\n", nbuf); 

  /* allocate buffer for data on processor 0 */
  if (processor0) {
    nbuffer = npoints * numvars;
    buffer = malloc(sizeof(double) * nbuffer);

    /* points that are not found are masked with a special value */
    for (i = 0; i < nbuffer; i++)
      buffer[i] = -1; // NODATA;
  }

  /* combine data from different Processors */
  /* copy nbuf buf to buffer */
  bampi_combinepointbuffers(level, nbuf, buf, buffer, numvars);

  /* processor 0 does the writing */
  if (processor0) {
   
    /* output for each variable */
    for (j = 0; j < numvars; j++) {
  
      /* filenames */
      snprintf(filename, 1000, "%s/%s.%s.l%d",
 	       Gets("outdir_r"), VarName(indexarray[j]), suffix, level->l);
 
      /* open file */
      fp = fopen(filename, "a");
      if (!fp) errorexits("failed opening %s", filename);      
      /* write data in gnuplot column format */
     
      if (level->iteration == 0) { // write header only for first time level
	
	fprintf(fp, "# gridfunction = %s\n", VarName(indexarray[j]));
	fprintf(fp, "# r = %14.6e\n", r);
	fprintf(fp, "# grid = {%d, %d}\n", ntheta, nphi);
	fprintf(fp, "# theta = {%14.6e, %14.6e}\n", dtheta, theta);
	fprintf(fp, "# phi   = {%14.6e, %14.6e}\n", 0.0e0,  phi);
	fprintf(fp, "# dtheta = %14.6e\n", dtheta);
	fprintf(fp, "# dphi   = %14.6e\n", dphi);
      }

      fprintf(fp, "# iteration = %d\n",     level->iteration);
      fprintf(fp, "# time      = %14.6e\n", level->time);
     
      for (i = 0; i < npoints; i++) {

        if (Getv("2doutputr_type", "cartesian")) {
	  fprintf(fp, "%16.9e %16.9e %16.9e  %22.15e\n",
		  coords[3*i], coords[3*i+1], coords[3*i+2], buffer[i * numvars + j]); 
	}

        if (Getv("2doutputr_type", "spherical")) {
	  fprintf(fp, "%16.9e %16.9e %22.15e\n",
		  sph_coords[2*i], sph_coords[2*i+1], buffer[i * numvars + j]);
	}

        if ( !Getv("2doutputr_type", "cartesian") && !Getv("2doutputr_type", "spherical")) {
          fprintf(fp, "%22.15e ", buffer[i * numvars + j]);
        }

        if ((i+1)%(nphi+1) == 0) fprintf(fp,"\n");
      }
      fprintf(fp,"\n");
	
      /* close file */
      fclose(fp);
    }
  }

  /* clean up */
  if (processor0) free(buffer);
  free(buf);
  free(coords);  
  free(sph_coords);
  timer_stop(0, "write_level_sphere0");

  /* wait, could be optimized across different writes */
  bampi_barrier();
}





/* 2d output for a sphere */
/* WT: just like write_level_sphere0, but output on the points used in 
   sgrid's SphericalDF */
void write_level_sphere_SphericalDF(tL *level, int numvars, int *indexarray, 
                                    int ntheta_f, int nphi_f, double r, 
                                    char *suffix)
{
  /* Function Arguments:

  level        ...  level struct
  numvars      ...  number of gridfunctions for output
  indexarray   ...  array containing indices for output grid functions
  ntheta_f     ...  number of points between theta=0 and theta=PI on full grid
  nphi_f       ...  number of points between phi=0 and phi=2PI on full grid
  r            ...  extraction radius
  suffix       ...  file suffix
  */

  int pr = 1;       // debugging switch: 0 is off, 1 is on
  int interp = 1;   // always allow interpolation in interpolation function

  static int ncall = 0;

  int nseries;      // 	   \_  timing info used for file header etc.
  double dtseries;  //     /

  int i, j, k, k2;  //  counters

  FILE *fp;
  char filename[1000];
  char gridtype[1000];

  int ntheta, nphi;  /* number of points in theta and phi on actual grid with syms */
  double dtheta, dphi, theta, phi;

  double *global_bbox = level->bbox;
  double *local_bbox  = level->com->bbox;

  double *coords, *sph_coords;
  double coord[3];

  int npoints, nlocalpoints;

  int nbuf, nbuffer;
  double *buf = 0, *buffer = 0;

  /* END OF DECLARATION SECTION */

  if (numvars <= 0) return;  // no variables requested for output 

  timer_start(0, "write_level_sphere_SphericalDF");
  ncall++;

  if (level->iteration == 0) printf("  extract at radius %e at level %d\n", r, level->l);

  if (r > global_bbox[1] && r > global_bbox[3] && r > global_bbox[5]) {
    if (level->iteration == 0) printf("bbox does not fit at this level\n");
    return;
  }


  nseries = timeforoutput_di_dt(level, Geti("2doutiter"), Getd("2douttime"));
  dtseries = (nseries == 1) ? 0 : level->iteration/(nseries-1)*level->dt;

  /* set grid spacing */
  dtheta = PI/ntheta_f;
  dphi   = 2.0*PI/nphi_f;

  /* set ntheta, nphi for full grid */
  ntheta = ntheta_f;
  nphi   = nphi_f; 
  
  /* check syms */
  if (Getv("grid", "bitant"))
  {
    /* z=0 reflection symmetry */
    ntheta = (ntheta_f + 1)/2;
    sprintf(gridtype, "%s", "bitant");
  }
  else if (Getv("grid", "rotant"))
  {
    /* x=y=0 inversion symmetry */
    nphi = (nphi_f + 1)/2;
    sprintf(gridtype, "%s", "rotant");
  }
  else if (Getv("grid", "quadrant"))
  {
    /* z=0 and y=0 reflection symmetry */
    ntheta = (ntheta_f + 1)/2;
    nphi = (nphi_f + 1)/2;
    sprintf(gridtype, "%s", "quadrant");
  }
  else if (Getv("grid", "octant"))
  {
    /* x=0, y=0, and z=0 reflection symmetry */
    ntheta = (ntheta_f + 1)/2;
    nphi = ( (nphi_f + 1)/2 + 2)/2;
    sprintf(gridtype, "%s", "octant");
  }
  else
  {
    /* Full grid */
    ntheta = ntheta_f;
    nphi   = nphi_f; 
    sprintf(gridtype, "%s", "full");
  }

  npoints = ntheta * nphi; // total # of points

  /* count points for this processor */
  nlocalpoints = 0;
  for(j=0; j < nphi;   j++)
    for(i=0; i < ntheta; i++)
    {
      theta = i * dtheta + 0.5*dtheta;
      phi   = j * dphi;

      coord[0] = r * sin(theta) * cos(phi);
      coord[1] = r * sin(theta) * sin(phi);
      coord[2] = r * cos(theta);
      
      if (level->shells) {
        convert_box_to_shells(find_shellsbox_from_xyz(coord[0],coord[1],coord[2]), 
                              coord[0],coord[1],coord[2],
                              &coord[0],&coord[1],&coord[2] );
      }
      
      if (insidebbox(local_bbox, coord))  nlocalpoints++;
    } 

  if (level->iteration == 0 && nlocalpoints > -1)
    printf("  found %d points from %d total points at proc %d\n", nlocalpoints, npoints, bampi_rank());

  /* allocate and fill in coordinates */
  coords     = (double *) malloc(sizeof(double) * 3 * npoints);
  sph_coords = (double *) malloc(sizeof(double) * 2 * npoints);

  k  = 0;
  k2 = 0;
  for(j=0; j < nphi; j++)
    for(i=0; i < ntheta; i++)
    {
      theta = i * dtheta + 0.5*dtheta;
      phi   = j * dphi;

      coord[0] = r * sin(theta) * cos(phi);
      coord[1] = r * sin(theta) * sin(phi);
      coord[2] = r * cos(theta);

      sph_coords[k2++] = theta;
      sph_coords[k2++] = phi;    

      coords[k++] = coord[0];
      coords[k++] = coord[1];
      coords[k++] = coord[2];
    }


  if (PR) // tell the world what coord vals we got
    for (i = 0; i < npoints; i++) {
      printf("%d %6.3f %6.3f %6.3f\n", 
	     i, coords[3*i], coords[3*i+1], coords[3*i+2]);
    }

  /* fill buffer with data at given list of coordinates, allocates buf */
  if (level->shells || 1)
  {
    //FIXME: this is a more general scheme, should be tested and should replace the other function
    bufferpoints_full(level, npoints, coords, numvars, indexarray, &nbuf, &buf, output_order,output_scheme);
  } 
  else
  {
    bufferorderedpoints(level, npoints, coords, numvars, indexarray, &nbuf, &buf, interp, output_order,output_scheme);
  }

  if (PR) printf("nbuf computed from bufferorderedpoints = %d\n", nbuf); 

  /* allocate buffer for data on processor 0 */
  if (processor0) {
    nbuffer = npoints * numvars;
    buffer = malloc(sizeof(double) * nbuffer);

    /* points that are not found are masked with a special value */
    for (i = 0; i < nbuffer; i++)
      buffer[i] = -1; // NODATA;
  }

  /* combine data from different Processors */
  /* copy nbuf buf to buffer */
  bampi_combinepointbuffers(level, nbuf, buf, buffer, numvars);

  /* processor 0 does the writing */
  if (processor0)
  {
    /* output for each variable */
    for (j = 0; j < numvars; j++)
    {
      /* filenames */
      snprintf(filename, 1000, "%s/%s_%s.l%d",
 	       Gets("outdir_r"), VarName(indexarray[j]), suffix, level->l);
 
      /* open file */
      fp = fopen(filename, "a");
      if (!fp) errorexits("failed opening %s", filename);      
      /* write data in gnuplot column format */
     
      if (level->iteration == 0)
      { // write header only for first time level
	fprintf(fp, "# gridfunction = %s\n", VarName(indexarray[j]));
	fprintf(fp, "# r = %-.16e\n", r);
	fprintf(fp, "# SphericalDF:   n2 = %d  n3 = %d\n", 2*ntheta_f, nphi_f);
	fprintf(fp, "# gridtype = %s\n", gridtype);
	fprintf(fp, "# ntheta = %d\n", ntheta);
	fprintf(fp, "# nphi   = %d\n", nphi);
	fprintf(fp, "# theta = { %-.16e, %-.16e }\n", 0.5*dtheta, theta);
	fprintf(fp, "# phi   = { %-.16e, %-.16e }\n", 0.0, phi);
	fprintf(fp, "# dtheta = %-.16e\n", dtheta);
	fprintf(fp, "# dphi   = %-.16e\n", dphi);
      }

      fprintf(fp, "# iteration = %d\n",     level->iteration);
      fprintf(fp, "# time      = %-.16e\n", level->time);
     
      for (i = 0; i < npoints; i++)
      {
        if (Getv("2doutputr_type", "cartesian"))
	  fprintf(fp, "%16.9e %16.9e %16.9e  %22.15e\n",
		  coords[3*i], coords[3*i+1], coords[3*i+2], buffer[i * numvars + j]); 

        if (Getv("2doutputr_type", "spherical"))
	  fprintf(fp, "%16.9e %16.9e %22.15e\n",
		  sph_coords[2*i], sph_coords[2*i+1], buffer[i * numvars + j]);

        if ( !Getv("2doutputr_type", "cartesian") && !Getv("2doutputr_type", "spherical"))
          fprintf(fp, "%22.15e ", buffer[i * numvars + j]);

        if ((i+1)%(ntheta) == 0) fprintf(fp,"\n");
      }
      fprintf(fp,"\n");
	
      /* close file */
      fclose(fp);
    }
  }

  /* clean up */
  if (processor0) free(buffer);
  free(buf);
  free(coords);  
  free(sph_coords);
  timer_stop(0, "write_level_sphere_SphericalDF");

  /* wait, could be optimized across different writes */
  bampi_barrier();
}





/* 2d sphere-output adjusted to spectral integration (Gauss-Legendre quadrature and composite trapezoidal rule), parallel */
/* Roman Gold, Nov. 2010 */
/* Evaluate Legendre Polynomials. Get abscissas and weights for Gauss-Legendre-quadrature */
void gauleg(double x1, double x2, double roots[], double w[], int n)
{
  int i, j, m;
  double eps = 1.0E-14;// BAM has it's own epsilon... Modify?

  double p1, p2, p3, pp, xl, xm, z, z1;
  
  m = (n+1)/2;
  xm = 0.5*(x2+x1);
  xl = 0.5*(x2-x1);
  for(i=1; i<=m; i++)
  {
    z = cos(PI*((double)i-0.25)/((double)n+0.5));
    while(1)
    {
      p1 = 1.0;
      p2 = 0.0;
      for(j=1; j<=n; j++)
      {
        p3 = p2;
        p2 = p1;
        p1 = ((2.0*(double)j-1.0)*z*p2-((double)j-1.0)*p3)/
              (double)j;
      }
      pp = (double)n*(z*p1-p2)/(z*z-1.0);
      z1 = z;
      z = z1 - p1/pp;
      if(fabs(z-z1) <= eps) break;
    }
    roots[i] = xm - xl*z;
    roots[n+1-i] = xm + xl*z;
    w[i] = 2.0*xl/((1.0-z*z)*pp*pp);
    w[n+1-i] = w[i];
  }
} /* end gauleg */

/* copied from write_level_sphere(). The difference is that in order to perform a Gauss-Legendre quadrature we need sphere data on different non-equally distant points */
void write_level_sphere_spectral(tL *level, int numvars, int *indexarray, 
                         int ntheta, int nphi, double r, char *suffix)
{
  /* Function Arguments:

  level        ...  level struct
  numvars      ...  number of gridfunctions for output
  indexarray   ...  array containing indices for output grid functions
  ntheta       ...  number of theta intervals 
  nphi         ...  number of phi intervals
  r            ...  extraction radius
  suffix       ...  file suffix
  */

  int pr = 0;       // debugging switch: 0 is off, 1 is on
  int interp = 1;   // always allow interpolation in interpolation function

  static int ncall = 0;

  int nseries;      // 	   \_  timing info used for file header etc.
  double dtseries;  //     /

  int i, j, k, k2;  //  counters

  FILE *fp;
  char filename[1000];

  double dtheta, dphi, theta, phi;

  double *global_bbox = level->bbox;
  double *local_bbox  = level->com->bbox;

  double *coords, *sph_coords;
  double coord[3];

  int npoints, nlocalpoints;

  int nbuf, nbuffer;
  double *buf = 0, *buffer = 0;

  double roots[ntheta+1];// gridpoints for Gauss-Legendre quadrature
  double w[ntheta+1];// Gauss-Legendre weights

  /* END OF DECLARATION SECTION */

  if (numvars <= 0) return;  // no variables requested for output 

  timer_start(0, "write_level_sphere0");
  ncall++;

  if (level->iteration == 0) printf("  extract at radius %e at level %d\n", r, level->l);

  if (r > global_bbox[1] && r > global_bbox[3] && r > global_bbox[5]) {
    if (level->iteration == 0) printf("bbox does not fit at this level\n");
    return;
  }

  npoints = (ntheta) * (nphi);


  nseries = timeforoutput_di_dt(level, Geti("2doutiter"), Getd("2douttime"));
  dtseries = (nseries == 1) ? 0 : level->iteration/(nseries-1)*level->dt;


/***************** SYMMETRIES **************************/
if (Getv("grid", "bitant")) {
//    /* z=0 reflection symmetry */
    ntheta = ntheta / 2 ;
    dphi = 2.0*PI/nphi;

  } else if (Getv("grid", "rotant")) {
//      /* x=y=0 inversion symmetry */
      ntheta = ntheta;
      dphi = PI/nphi;

  } else if (Getv("grid", "quadrant")) {
//    /* z=0 and y=0 reflection symmetry */
    ntheta = ntheta / 2;
    dphi = PI/nphi;

  } else if (Getv("grid", "octant")) {
//    /* x=0, y=0, and z=0 reflection symmetry */
    ntheta = ntheta / 2;
    dphi = PI/(2*nphi);
  } else {
    /* Full grid */
    ntheta = ntheta;
    dphi = 2.0 * PI / nphi;
  }
/*******************************************************/

if (0/*ntheta==47*/){
/****** We could precompute gauss-legendre quadrature here ********/
  }
else gauleg(-1.,1., roots, w, ntheta);

  /* count points for this processor */

  if (pr) {
    nlocalpoints = 0;

    // handle symmetries: compute roots/w's with ntheta, loop with ntheta_loop depends on symmetries 
    for (i=1; i <= ntheta  ; i++){

      for (j=1; j <= nphi ; j++){

	theta = acos(roots[i]); // as computed by gauleg() see above or precomputed...
	phi = j*dphi;

	coord[0] = r * sin(theta) * cos(phi);
	coord[1] = r * sin(theta) * sin(phi);
	coord[2] = r * cos(theta);

	if (insidebbox(local_bbox, coord)) {
	  nlocalpoints++;
	} 
      } 
    }
  } else {
    nlocalpoints = -1;
 }

  if (level->iteration == 0 && nlocalpoints > -1) {
    printf("  found %d points from %d total points at proc %d\n", nlocalpoints, npoints, bampi_rank());
  }

  /* allocate and fill in coordinates */
  coords     = (double *) malloc(sizeof(double) * 3 * npoints);
  sph_coords = (double *) malloc(sizeof(double) * 2 * npoints);

  k  = 0;
  k2 = 0;

  for   (i=1; i <= ntheta ; i++){

    for (j=1; j <= nphi ; j++){


      theta = acos(roots[i]); // as computed by gauleg() see above or precomputed...

      phi = j * dphi;

      coord[0] = r * sin(theta) * cos(phi);
      coord[1] = r * sin(theta) * sin(phi);
      coord[2] = r * cos(theta);

      sph_coords[k2++] = theta;
      sph_coords[k2++] = phi;    

      coords[k++] = coord[0];
      coords[k++] = coord[1];
      coords[k++] = coord[2];
    }
  }

#ifdef WSDEBUG
  if (pr) // tell the world what coord vals we got
    for (i = 0; i < npoints; i++) {
      printf("%d %6.3f %6.3f %6.3f\n", 
	     i, coords[3*i], coords[3*i+1], coords[3*i+2]);
    }
#endif

  /* fill buffer with data at given list of coordinates, allocates buf */
  bufferorderedpoints(level, npoints, coords, numvars, indexarray, &nbuf, &buf, interp, output_order,output_scheme);

#ifdef WSDEBUG
  printf("nbuf computed from bufferorderedpoints = %d\n", nbuf); 
#endif

  /* allocate buffer for data on processor 0 */
  if (processor0) {
    nbuffer = npoints * numvars;
    buffer = malloc(sizeof(double) * nbuffer);

    /* points that are not found are masked with a special value */
    for (i = 0; i < nbuffer; i++)
      buffer[i] = -1; // NODATA;
  }

  /* combine data from different Processors */
  /* copy nbuf buf to buffer */
  bampi_combinepointbuffers(level, nbuf, buf, buffer, numvars);

  /* processor 0 does the writing */
  if (processor0) {
   
    /* output for each variable */
    for (j = 0; j < numvars; j++) {

      /* filenames */
      snprintf(filename, 1000, "%s/%s.%s.l%d",
 	       Gets("outdir_r"), VarName(indexarray[j]), suffix, level->l);

      /* open file */
      fp = fopen(filename, "a");
      if (!fp) errorexits("failed opening %s", filename);      
      /* write data in gnuplot column format */
     
      if (level->iteration == 0) { // write header only for first time level
	fprintf(fp, "# gridfunction = %s\n", VarName(indexarray[j]));
	fprintf(fp, "# r = %14.6e\n", r);
	fprintf(fp, "# grid = {%d, %d}\n", ntheta, nphi);
	fprintf(fp, "# theta = {%14.6e, %14.6e}\n", dtheta, theta);
	fprintf(fp, "# theta = {%14.6e, %14.6e}\n", acos(roots[1]), acos(roots[ntheta]));// +1???
	fprintf(fp, "# phi   = {%14.6e, %14.6e}\n", 0.0e0,  phi);
	fprintf(fp, "# dphi   = %14.6e\n", dphi);
      }
      
      fprintf(fp, "# iteration = %d\n",     level->iteration);
      fprintf(fp, "# time      = %14.6e\n", level->time);
      
      for (i = 0; i < npoints; i++) {


        if (Getv("2doutputr_type", "cartesian")) {
	  fprintf(fp, "%16.9e %16.9e %16.9e  %22.15e\n",
		  coords[3*i], coords[3*i+1], coords[3*i+2], buffer[i * numvars + j]); 
	}

        if (Getv("2doutputr_type", "spherical")) {
	  fprintf(fp, "%16.9e %16.9e %22.15e\n",
		  sph_coords[2*i], sph_coords[2*i+1], buffer[i * numvars + j]);
	}

        if ( !Getv("2doutputr_type", "cartesian") && !Getv("2doutputr_type", "spherical")) {
          fprintf(fp, "%22.15e ", buffer[i * numvars + j]);
	}

        if ((i+1)%(nphi) == 0) fprintf(fp,"\n");

      }
      fprintf(fp,"\n");
	
      /* close file */
      fclose(fp);
    }
  }

  /* clean up */

  if (processor0) free(buffer);
  free(buf);
  free(coords);  
  free(sph_coords);
  timer_stop(0, "write_level_sphere0");

  /* wait, could be optimized across different writes */
  bampi_barrier();
}









/* wrapper for write_level_sphere0 and write_level_sphere_SphericalDF
   decides in which format we write sphere data */
void write_level_sphere(tL *level, int numvars, int *indexarray, 
                        int ntheta, int nphi, double r, char *suffix)
{
  if(Getv("2doutputr_type", "SphericalDF"))
    write_level_sphere_SphericalDF(level, numvars, indexarray, ntheta, nphi, r, suffix);
  else if (Getv("2doutputr_type","GaussLegendre")){
    write_level_sphere_spectral(level, numvars, indexarray, ntheta, nphi, r, suffix);}
  else
    write_level_sphere0(level, numvars, indexarray, ntheta, nphi, r, suffix);
}

