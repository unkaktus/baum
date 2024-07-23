/* output3d.c */
/* Bernd Bruegmann 12/02 */

#include "bam.h"
#include "output.h"

#define PR 0



/* make list of ordered points for the global grid
   this function is used for compatibility with the current 1d output
   function setup, but 3d could be optimized differently
*/
void makeorderedlist3d(tL *level, int interp, int *npoints, double **coords)
{
  double dx = level->dx;
  double dy = level->dy;
  double dz = level->dz;
  double x0, y0, z0, x1, y1, z1;
  int nx, ny, nz;
  int nmax = (level->ibbox[1]+1) * (level->ibbox[3]+1) *( level->ibbox[5]+1);
  int i, j, k, n;

  /* preliminary: assume cubical box */
  *npoints = nmax;
  *coords = malloc(sizeof(double) * 3 * nmax);
  if (!*coords) 
    errorexit("Out of memory in output3d.");

  /* preliminary: need MPI for
     global number of points
     global coordinates
  */

  /* define default for global box */
  x0 = level->bbox[0];
  y0 = level->bbox[2];
  z0 = level->bbox[4];
  x1 = level->bbox[1];
  y1 = level->bbox[3];
  z1 = level->bbox[5];

  /* amr output supports interpolation to nodes */
  if (interp) {
    x0 += dx/2;
    y0 += dy/2;
    z0 += dz/2;
    if (Getv("grid", "octant"))   x0 = y0 = z0 = 0;
    if (Getv("grid", "quadrant") || Getv("grid", "qreflect")) y0 = z0 = 0;
    if (Getv("grid", "rotant"))   y0 = 0;
    if (Getv("grid", "bitant"))   z0 = 0;
  }

  /* fill in coordinates */
  nx = (x1-x0)/dx + .1;
  ny = (y1-y0)/dy + .1;
  nz = (z1-z0)/dz + .1;
  n = 0;
  for (k = 0; k <= nz; k++)
  for (j = 0; j <= ny; j++)
  for (i = 0; i <= nx; i++) {
    (*coords)[n++] = x0 + i*dx;
    (*coords)[n++] = y0 + j*dy;
    (*coords)[n++] = z0 + k*dz;
  }
}


/* check whether this is a vector or not and remeber until next call*/
int  ckeckifvector(double *buffer, double *vbuffer, int ivn, char *varname, int nbuffer)
{
  int i;
  static char varname_static[100];
    
  sprintf(varname,"%s",VarName(ivn));
  varname[strlen(varname)-1]='\0';
    //printf("%s  %s\n",varname,VarName(ivn));
    
  if (VarComponent(ivn)==0) {
    if (VarNComponents(ivn) != 3) return 0;
    sprintf(varname_static,"%s",varname);
  } else {
    if (strcmp(varname,varname_static)!=0) return 0;
  }
    
  for (i=0; i<nbuffer; i++)
    vbuffer[VarComponent(ivn)*nbuffer + i] = buffer[i];
    
  if (VarComponent(ivn)==2) return 1;
    
  return 0;
}


/* 3d output: volume */
void write_level_3d(tL *level, int nv, int *iv, int sampling, char *suffix)
{
  if (PR) printf("write_level_3d \n");
  
  int pr = 0;
  int text   = Getv("3dformat", "text");
  int binary = Getv("3dformat", "binary");
  int opendx = Getv("3dformat", "opendx");
  int vtk    = Getv("3dformat", "vtk");
  int xdmf   = Getv("3dformat", "xdmf");
  int flt    = Getv("3dformat", "float");
  int dbl    = Getv("3dformat", "double");
  int cells  = Getv("3dformat", "cells");
  int interp = Getv("3doutinterpolate", "yes");
  int nseries;
  double dtseries;
  static int ncall = 0;
  float xfloat;
  double xdouble;
  int sfloat = sizeof(float);
  int sdouble = sizeof(double);
  int i, j, k, m, n;
  FILE *fp;
  char filename[1000];
  int npoints;
  double *coords;
  int index[1];
  int nbuf, nbuffer;
  double *buf = 0, *buffer = 0;
  double *vbuffer = 0;

  if (nv <= 0) return;
  ncall++;
  nseries = timeforoutput_di_dt(level, Geti("3doutiter"), Getd("3douttime"));
  dtseries = (nseries == 1) ? 0 : level->iteration/(nseries-1)*level->dt;
  if (interp) 
    errorexit("interpolation is not supported for 3d output,"
	      " but easily could be");
  /* parameter defaults */
  if (!flt && !dbl) flt = 1;
  if (!text && !binary) binary = 1;

  /* make global list of ordered points */
  makeorderedlist3d(level, interp, &npoints, &coords);
  if (pr) 
    for (i = 0; i < npoints; i++)
      printf("%d %6.3f %6.3f %6.3f\n", 
	     i, coords[3*i], coords[3*i+1], coords[3*i+2]);

  /* allocate buffer for data on processor 0 */
  if (processor0) {
    nbuffer = npoints;
    buffer = dmalloc(nbuffer);
    vbuffer = dmalloc(3*nbuffer);
  }

  /* for each variable */
  for (n = 0; n < nv; n++) {
    index[0] = iv[n];

    /* fill buffer with data at given list of coordinates, allocates buf */
    bufferorderedpoints(level, npoints, coords, 1, index, &nbuf, &buf, interp, output_order,output_scheme);
    if (pr) {
      printf("npoints = %d, nbuf = %d\n", npoints, nbuf);
      for (i = 0; i < nbuf; i++) {
	printf(" %10.3e", buf[i]);
	if ((i+1) % (nv+1) == 0) printf("\n");
      }
    }
   
    /* combine data from different processors */
    bampi_combinepointbuffers(level, nbuf, buf, buffer, 1);

    /* processor 0 does the writing */
    if (processor0) {
      double x0 = level->bbox[0];
      double y0 = level->bbox[2];
      double z0 = level->bbox[4];
      double dx = level->dx;
      double dy = level->dy;
      double dz = level->dz;
      int nx = level->ibbox[1] + 1;
      int ny = level->ibbox[3] + 1;
      int nz = level->ibbox[5] + 1;
      int nbuffer = npoints;

      /* we may want to transform the data */
      // output_transform(level, nbuffer, coords, buffer, &tnbuffer, &tbuffer,
      //       &x0, &y0, &z0, &dx, &dy, &dz, &nx, &ny, &nz);

      /* we may want to add more points by symmetry */
      if (Getv("3dformat", "add_rotant_points"))
	add_rotant_points_to_buffer(level, iv[n], &nbuffer, &buffer,
				    &x0, &y0, &z0, &nx, &ny, &nz);
   
      /* do output into name.x0 */
      if (!vtk && !xdmf) {

	/* filename */
	snprintf(filename, 1000, "%s/%s.%s%d%s", 
		 Gets("outdir_3d"), VarName(iv[n]), suffix, level->l, boxstring);

	/* open file */
	fp = fopen(filename, (text)?"a":"ab");
	if (!fp) errorexits("failed opening %s", filename);
	
	/* primitive column format */
	if (!opendx && !vtk) {
	  fprintf(fp, "\"Time = %f\"\n", level->time);
	  m = 0;
	  for (k = 0; k < nz; k++)
	  for (j = 0; j < ny; j++)
	  for (i = 0; i < nx; i++)
	    fprintf(fp, "%16.9e %16.9e %16.9e  %22.15e\n", 
		    x0 + i*dx, y0 + j*dy, z0 + k*dz, buffer[m++]);
	  fprintf(fp,"\n");
	}

	/* OpenDX */
	if (opendx) {
	  if (text) {
	    if (dbl) 
	      for (i = 0; i < nbuffer; i++)
		fprintf(fp, "%22.15e\n", buffer[i]);
	    else
	      for (i = 0; i < nbuffer; i++)
		fprintf(fp, "%16.9e\n", buffer[i]);
	  }
	  if (binary) {
	    if (dbl)
	      for (i = 0; i < nbuffer; i++) {
		xdouble = buffer[i];
		fwrite(&xdouble, sdouble, 1, fp);
	      }
	    else
	      for (i = 0; i < nbuffer; i++) {
		xfloat = buffer[i];
		fwrite(&xfloat, sfloat, 1, fp);
	      }
	  }
	}

	/* close file */
	fclose(fp);
      }



      /* OpenDX: general array format descriptor file
	 - since there doesn't appear to be a way to have a variable
	   length time series, we rewrite the descriptor file 
	   every time we appended a new time slice of data
         - the number of steps has to be given in the beginning of the file
           in the header: is it possible to change a fixed width number there
	   without rewriting/copying all of the 3d data?
      */
      if (opendx) {

	/* open descriptor file */
	strcat(filename, ".dx");
	fp = fopen(filename, "w");
	if (!fp) errorexits("failed opening %s", filename);

	/* write header */
	fprintf(fp, "file = %s.%s%d%s\n", 
		VarName(iv[n]), suffix, level->l, boxstring);
	fprintf(fp, "grid = %d x %d x %d\n", nx+cells, ny+cells, nz+cells);
	if (text)   fprintf(fp, "format = ascii\n");
	if (binary) fprintf(fp, "format = lsb binary\n");
	fprintf(fp, "interleaving = record\n");
	fprintf(fp, "majority = column\n");
	fprintf(fp, "series = %d, %f, %f\n", nseries, 0.0, dtseries);
	fprintf(fp, "field = %s\n", VarName(iv[n]));
	fprintf(fp, "structure = scalar\n");
	if (flt) fprintf(fp, "type = float\n");
	if (dbl) fprintf(fp, "type = double\n");
	if (cells) fprintf(fp, "dependency = connections\n");
	else       fprintf(fp, "dependency = positions\n");
	fprintf(fp, "positions = regular, regular, regular,");
	fprintf(fp, " %f, %f, %f, %f, %f, %f\n", x0, dx, y0, dy, z0, dz);
	fprintf(fp, "\nend\n");

	/* close file */
	fclose(fp);
      }


      /* VTK */
      /* output one file per time step in separate subdirectories */
      if (vtk) {

	/* open file (returns non-null file pointer) */
	fp = fopen_vtk(VarName(iv[n]), "3d",suffix, level->l, nseries-1);

	/* write header */
	fprintf(fp, "# vtk DataFile Version 2.0\n");
	fprintf(fp, "variable %s, level %d, time %16.9e\n", 
		VarName(iv[n]), level->l, level->time);
	fprintf(fp, binary ? "BINARY" : "ASCII\n");
	fprintf(fp, "\n");
	fprintf(fp, "DATASET STRUCTURED_POINTS\n");
	fprintf(fp, "DIMENSIONS %d %d %d\n", nx, ny, nz);
	fprintf(fp, "ORIGIN  %16.9e %16.9e %16.9e\n", x0, y0, z0);
	fprintf(fp, "SPACING %16.9e %16.9e %16.9e\n", dx, dy, dz);
	fprintf(fp, "\n");
	fprintf(fp, "POINT_DATA %d\n", nbuffer);
	fprintf(fp, "SCALARS scalars %s\n", dbl ? "double" : "float");
	fprintf(fp, "LOOKUP_TABLE default\n");

	/* write data */
	write_raw_vtk(fp, nbuffer, buffer, 1, 0, dbl, flt, text, binary);
	fclose(fp);
        
        /* add vector data to a new file, if the first component is write-active */
        char varname[1000];
        if ((Getv("3doutputall","yes")) && (ckeckifvector(buffer, vbuffer, iv[n], varname, nbuffer))) {
            
            /* open file (returns non-null file pointer) */
            fp = fopen_vtk(varname, "3d", suffix, level->l, nseries-1);

            /* write header */
            fprintf(fp, "# vtk DataFile Version 2.0\n");
            fprintf(fp, "variable %s, level %d, time %16.9e\n", 
                    varname, level->l, level->time);
            fprintf(fp, "ASCII\n");
            fprintf(fp, "\n");
            fprintf(fp, "DATASET STRUCTURED_POINTS\n");
            fprintf(fp, "DIMENSIONS %d %d %d\n", nx, ny, nz);
            fprintf(fp, "ORIGIN  %16.9e %16.9e %16.9e\n", x0, y0, z0);
            fprintf(fp, "SPACING %16.9e %16.9e %16.9e\n", dx, dy, dz);
            fprintf(fp, "\n");
            fprintf(fp, "POINT_DATA %d\n", nbuffer); 
            fprintf(fp, "VECTORS normal float\n");

            /* write data */
            write_raw_vec_vtk(fp, nbuffer,vbuffer, dbl, flt, text, binary);
            fclose(fp);
        }
      }

      if (xdmf) {
        double time = (nseries-1)*Getd("3douttime");

        /* write h5 */
        #ifdef XDMF
        if(GetvLax("xdmfgrid","regular")){
          write_xdmf3D(VarName(iv[n]), "3d",suffix, level->l, time,
                      nbuffer, buffer, dbl, flt, text, binary, nx, ny, nz, x0, y0, z0, dx, dy, dz);
          
          /* add vector data to a new file, if the first component is write-active */
          char varname[1000];
          if ((Getv("3doutputall","yes")) && (ckeckifvector(buffer, vbuffer, iv[n], varname, nbuffer))) {
            write_xdmf3D_vec(varname, "3d",suffix, level->l, time,
                      nbuffer, vbuffer, dbl, flt, text, binary, nx, ny, nz, x0, y0, z0, dx, dy, dz);
          }
        } else if(GetvLax("xdmfgrid","curvilinear")){
          write_xdmf3D_curvilinear(VarName(iv[n]), "3d",suffix, level->l, time,
                nbuffer, buffer, dbl, flt, text, binary, nx, ny, nz, x0, y0, z0, dx, dy, dz);

          /* add vector data to a new file, if the first component is write-active */
          char varname[1000];
          if ((Getv("3doutputall","yes")) && (ckeckifvector(buffer, vbuffer, iv[n], varname, nbuffer))) {
            write_xdmf3D_curvilinear_vec(varname, "3d",suffix, level->l, time,
                      nbuffer, vbuffer, dbl, flt, text, binary, nx, ny, nz, x0, y0, z0, dx, dy, dz);
          }
        } else if(GetvLax("xdmfgrid","rectilinear")){
          write_xdmf3D_rectilinear(VarName(iv[n]), "3d",suffix, level->l, time,
                nbuffer, buffer, dbl, flt, text, binary, nx, ny, nz, x0, y0, z0, dx, dy, dz);

          /* add vector data to a new file, if the first component is write-active */
          char varname[1000];
          if ((Getv("3doutputall","yes")) && (ckeckifvector(buffer, vbuffer, iv[n], varname, nbuffer))) {
            write_xdmf3D_rectilinear_vec(varname, "3d",suffix, level->l, time,
                      nbuffer, vbuffer, dbl, flt, text, binary, nx, ny, nz, x0, y0, z0, dx, dy, dz);
          }
        } else errorexit("XDMF grid format is not implemented!");
        #else 
        errorexit("XDMF library and flag must be installed to use xdmf output format!");
        #endif
      }
    }

    /* wait, could be optimized across different writes */
    bampi_barrier();
  }

  /* clean up */
  if (processor0) {
      free(buffer);
      free(vbuffer);
  }
  free(buf);
  free(coords);
}
