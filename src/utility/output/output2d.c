/* output2d.c */
/* Bernd Bruegmann 12/02 */

#include "bam.h"
#include "output.h"

#define PR 0


/* make list of ordered points for a given plane based on global bounding box
   this analogous to makeorderedlist1d, have a look at that one first
*/
void makeorderedlist2d(tL *level, int snapflag, 
		       double dx1, double dy1, double dz1, 
		       double dx2, double dy2, double dz2, 
		       double *ox, double *oy, double *oz, 
		       int *npoints1, int *npoints2,
		       int *npoints, double **coords)
{
  double x0 = Getd("outputx0");
  double y0 = Getd("outputy0");
  double z0 = Getd("outputz0");
  double coord[3];
  double dy;
  double *bbox = level->bbox;
  int nmax = level->ibbox[1] + level->ibbox[3] + level->ibbox[5];
  int i, j, k, m1, n1, m2, n2;

  *npoints = *npoints1 = *npoints2 = 0;
  *coords = 0;

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

  /* origin may not be in bbox, handle special amr/fmr cases here */
  if (!xyzinsidebbox(bbox, x0, y0, z0)) {

    /* x-y-plane */
    if ((dx1 && dy2 || dx2 && dy1) &&
        (!dless(bbox[5], z0) && !dless(z0, bbox[4]))) {
      x0 = bbox[0] + (!snapflag)*level->dx/2;
      y0 = bbox[2] + (!snapflag)*level->dy/2;
    }
    /* x-z-plane */
    else if ((dx1 && dz2 || dx2 && dz1) &&
        (!dless(bbox[3], y0) && !dless(y0, bbox[2]))) {
      x0 = bbox[0] + (!snapflag)*level->dx/2;
      z0 = bbox[4] + (!snapflag)*level->dz/2;
    }
    /* y-z-plane */
    else if ((dy1 && dz2 || dy2 && dz1) &&
        (!dless(bbox[1], z0) && !dless(z0, bbox[0]))) {
      y0 = bbox[2] + (!snapflag)*level->dy/2;
      z0 = bbox[4] + (!snapflag)*level->dz/2;
    }
    /* ignore diagonal planes if origin is not in box anyway */
    else
      return;
  }

  /* walk in direction one, backwards and forwards */
  for (m1 = 0; m1 > -nmax; m1--) {
    if (dless(x0 + m1*dx1, bbox[0])) break;
    if (dless(y0 + m1*dy1, bbox[2])) break;
    if (dless(z0 + m1*dz1, bbox[4])) break;
  }
  m1++;
  for (n1 = 0; n1 < nmax; n1++) {
    if (dless(bbox[1], x0 + n1*dx1)) break;
    if (dless(bbox[3], y0 + n1*dy1)) break;
    if (dless(bbox[5], z0 + n1*dz1)) break;
  }
  n1--;
  *npoints1 = n1 - m1 + 1;

  /* walk in direction two, backwards and forwards */
  for (m2 = 0; m2 > -nmax; m2--) {
    if (dless(x0 + m2*dx2, bbox[0])) break;
    if (dless(y0 + m2*dy2, bbox[2])) break;
    if (dless(z0 + m2*dz2, bbox[4])) break;
  }
  m2++;
  for (n2 = 0; n2 < nmax; n2++) {
    if (dless(bbox[1], x0 + n2*dx2)) break;
    if (dless(bbox[3], y0 + n2*dy2)) break;
    if (dless(bbox[5], z0 + n2*dz2)) break;
  }
  n2--;
  *npoints2 = n2 - m2 + 1;

  /* now we know the number of points */
  if (*npoints1 <= 0 || *npoints2 <= 0) return;
  *npoints = (*npoints1) * (*npoints2);
  *coords = malloc(sizeof(double) * 3 * (*npoints));

  /* fill in coordinates */
  k = 0;
  for (j = m2; j <= n2; j++)
  for (i = m1; i <= n1; i++) {
    (*coords)[k++] = x0 + i*dx1 + j*dx2;
    (*coords)[k++] = y0 + i*dy1 + j*dy2;
    (*coords)[k++] = z0 + i*dz1 + j*dz2;
  }

  /* fill in origin */
  *ox = x0 + m1*dx1 + m2*dx2;
  *oy = y0 + m1*dy1 + m2*dy2;
  *oz = z0 + m1*dz1 + m2*dz2;
}




/* 2d output: planes */
void write_level_2d(tL *level, int nv, int *iv, int plane, char *suffix)
{
  if (PR) printf("write_level_2d  plane %d\n",plane);
  
  int pr = 0;
  int text   = Getv("2dformat", "text");
  int binary = Getv("2dformat", "binary");
  int opendx = Getv("2dformat", "opendx");
  int vtk    = Getv("2dformat", "vtk");
  int xdmf   = Getv("2dformat", "xdmf");
  int flt    = Getv("2dformat", "float");
  int dbl    = Getv("2dformat", "double");
  int cells  = Getv("2dformat", "cells");
  int interp = Getv("2doutinterpolate", "yes");
  int nseries;
  double dtseries;
  static int ncall = 0;
  float xfloat;
  double xdouble;
  int sfloat = sizeof(float);
  int sdouble = sizeof(double);
  int i, j, k,l;
  FILE *fp;
  char filename[1000];
  char *format;
  double dx = level->dx;
  double dy = level->dy;
  double dz = level->dz;
  double exy = sqrt(dx*dx+dy*dy);
  double exz = sqrt(dx*dx+dz*dz);
  double eyz = sqrt(dy*dy+dz*dz);
  double d1, d2;
  double ox, oy, oz, o1, o2;
  double *coords;
  int npoints, npoints1, npoints2;
  int nbuf, nbuffer;
  double *buf = 0, *buffer = 0;

  if (nv <= 0) return;
  ncall++;
  nseries = timeforoutput_di_dt(level, Geti("2doutiter"), Getd("2douttime"));
  dtseries = (nseries == 1) ? 0 : level->iteration/(nseries-1)*level->dt;
  if (pr) {
    printf("%d. 2d output, plane %d: ", ncall, plane);
    for (i = 0; i < nv; i++)
      printf(" %s", VarName(iv[i]));
    printf("\n");
  }
  
  /* make list of points for this 2d plane based on global bounding box 
     - we pass three 3-tuples: the origin, direction 1, direction 2
     - note that we have to order the directions correctly so that the
       resulting point list is ordered (in particular plane "zd"!)
  */
  i = !interp;
  if (plane == 0)
    makeorderedlist2d(level, i, dx,0,0, 0,dy,0, 
		      &ox, &oy, &oz, &npoints1, &npoints2, &npoints, &coords);
  if (plane == 1)
    makeorderedlist2d(level, i, dx,0,0, 0,0,dz, 
		      &ox, &oy, &oz, &npoints1, &npoints2, &npoints, &coords);
  if (plane == 2)
    makeorderedlist2d(level, i, 0,dy,0, 0,0,dz, 
		      &ox, &oy, &oz, &npoints1, &npoints2, &npoints, &coords);
  if (plane == 3)
    makeorderedlist2d(level, i, dx,0,0, 0,dy,dz, 
		      &ox, &oy, &oz, &npoints1, &npoints2, &npoints, &coords);
  if (plane == 4)
    makeorderedlist2d(level, i, 0,dy,0, dx,0,dz, 
		      &ox, &oy, &oz, &npoints1, &npoints2, &npoints, &coords);
  if (plane == 5)
    makeorderedlist2d(level, i, dx,dy,0, 0,0,dz, 
		      &ox, &oy, &oz, &npoints1, &npoints2, &npoints, &coords);
  if (pr) 
    for (i = 0; i < npoints; i++)
      printf("%d(%d) %6.3f %6.3f %6.3f\n", 
	     i,npoints, coords[3*i], coords[3*i+1], coords[3*i+2]);

  /* the plane may not intersect the points on this level */
  if (npoints == 0) return;

  /* determine plane dependent grid spacing and origin 
     fix: assumes origin has negative coordinates 
  */
  if (plane == 0) {d1 = dx; d2 = dy;  o1 = ox; o2 = oy;}
  if (plane == 1) {d1 = dx; d2 = dz;  o1 = ox; o2 = oz;}
  if (plane == 2) {d1 = dy; d2 = dz;  o1 = oy; o2 = oz;}
  if (plane == 3) {d1 = dx; d2 = eyz; o1 = ox; o2 = -sqrt(oy*oy+oz*oz);}
  if (plane == 4) {d1 = dy; d2 = exz; o1 = oy; o2 = -sqrt(ox*ox+oz*oz);}
  if (plane == 5) {d1 = exy; d2 = dz; o1 = -sqrt(ox*ox+oy*oy); o2 = oz;}


  /* fill buffer with data at given list of coordinates, allocates buf */
  bufferorderedpoints(level, npoints, coords, nv, iv, &nbuf, &buf, interp, output_order,output_scheme);

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
    
    /* for each variable */
    for (j = 0; j < nv; j++) {

      /* do output into name.x0 */
      if (!vtk && !xdmf) {

	/* filenames */
	snprintf(filename, 1000, "%s/%s.%s%d%s", 
		 Gets("outdir_2d"), VarName(iv[j]), suffix, level->l, boxstring);
  
	/* open file */
	fp = fopen(filename, (text)?"a":"ab");
	if (!fp) errorexits("failed opening %s", filename);
	
	/* primitive column format */
	if (!opendx && !vtk && !xdmf) {
	  fprintf(fp, "\"Time = %f\"\n", level->time);
	  for (i = 0; i < npoints; i++)
	    fprintf(fp, "%16.9e %16.9e %16.9e  %22.15e\n", 
		    coords[3*i], coords[3*i+1], coords[3*i+2], buffer[i*nv+j]);
	  fprintf(fp,"\n");
	}
	
	/* OpenDX */
	if (opendx) {
	  if (text) {
	    if (dbl) 
	      for (i = 0; i < npoints; i++)
		fprintf(fp, "%22.15e\n", buffer[i*nv+j]);
	    else
	      for (i = 0; i < npoints; i++)
		fprintf(fp, "%16.9e\n", buffer[i*nv+j]);
	  }
	  if (binary) {
	    if (dbl)
	      for (i = 0; i < npoints; i++) {
		xdouble = buffer[i*nv+j];
		fwrite(&xdouble, sdouble, 1, fp);
	      }
	    else
	      for (i = 0; i < npoints; i++) {
		xfloat = buffer[i*nv+j];
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
		VarName(iv[j]), suffix, level->l, boxstring);
	fprintf(fp, "grid = %d x %d\n", npoints1+cells, npoints2+cells);
	if (text)   fprintf(fp, "format = ascii\n");
	if (binary) fprintf(fp, "format = lsb binary\n");
	fprintf(fp, "interleaving = record\n");
	fprintf(fp, "majority = column\n");
	fprintf(fp, "series = %d, %f, %f\n", nseries, 0.0, dtseries);
	fprintf(fp, "field = %s\n", VarName(iv[j]));
	fprintf(fp, "structure = scalar\n");
	if (flt) fprintf(fp, "type = float\n");
	if (dbl) fprintf(fp, "type = double\n");
	if (cells) fprintf(fp, "dependency = connections\n");
	else       fprintf(fp, "dependency = positions\n");
	fprintf(fp, "positions = regular, regular,");
	fprintf(fp, " %f, %f, %f, %f\n", o1, d1, o2, d2);
	fprintf(fp, "\nend\n");

	/* close file */
	fclose(fp);
      }


      /* VTK */
      /* output one file per time step in separate subdirectories */
      /* THIS OUTPUT IS NOT WELL DESIGNED -> xy is ok, but xz is plotted in xy plane
         so be aware of it if you want to look at this data */
      if (vtk) {

	/* open file (returns non-null file pointer) */
	fp = fopen_vtk(VarName(iv[j]),"2d", suffix, level->l, nseries-1);

	/* write header */
	fprintf(fp, "# vtk DataFile Version 2.0\n");
	fprintf(fp, "variable %s, level %d, time %16.9e\n", 
		VarName(iv[j]), level->l, level->time);
	fprintf(fp, binary ? "BINARY" : "ASCII\n");
	fprintf(fp, "\n");
	fprintf(fp, "DATASET STRUCTURED_POINTS\n");
	fprintf(fp, "DIMENSIONS %d %d %d\n", npoints1, npoints2, 1);
	fprintf(fp, "ORIGIN  %16.9e %16.9e %16.9e\n", o1, o2, 0.0);
	fprintf(fp, "SPACING %16.9e %16.9e %16.9e\n", d1, d2, (d1+d2)/2);
	fprintf(fp, "\n");
	fprintf(fp, "POINT_DATA %d\n", npoints);
	fprintf(fp, "SCALARS scalars %s\n", dbl ? "double" : "float");
	fprintf(fp, "LOOKUP_TABLE default\n");

	/* write data */
	write_raw_vtk(fp, npoints, buffer, nv, j, dbl, flt, text, binary);
	fclose(fp);
        
        
        /* this is for vector-output, it does not have so many output-options like the normal vtk output 
           does not work for diagonal output (xd,yd,zd)... this structure does not allow this easily 
           if you use yz or xz output the data is then always in the xy-plane ... lazynes*/
        if (strcmp(suffix,"xy")==0) { k=0; l=1; }
        else if (strcmp(suffix,"xz")==0) { k=0; l=2; }
        else if (strcmp(suffix,"yz")==0) { k=1; l=2; }
        else k=-1;
        if ((Getv("2doutputall","yes")) && (VarComponent(iv[j])==0) && (VarNComponents(iv[j]) == 3) && (k!=-1)) {
            
            /* open file (returns non-null file pointer) */
            char varname[1000];
            snprintf(varname,1000,"%s",VarName(iv[j]));
            varname[strlen(varname)-1]='\0';
            fp = fopen_vtk(varname,"2d", suffix, level->l, nseries-1);

            /* write header */
            fprintf(fp, "# vtk DataFile Version 2.0\n");
            fprintf(fp, "variable %s, level %d, time %16.9e\n", 
                    varname, level->l, level->time);
            fprintf(fp, "ASCII\n");
            fprintf(fp, "\n");
            fprintf(fp, "DATASET STRUCTURED_POINTS\n");
            fprintf(fp, "DIMENSIONS %d %d %d\n", npoints1, npoints2, 1);
            fprintf(fp, "ORIGIN  %16.9e %16.9e %16.9e\n", o1, o2, 0.0);
            fprintf(fp, "SPACING %16.9e %16.9e %16.9e\n", d1, d2, (d1+d2)/2);
            fprintf(fp, "\n");
            fprintf(fp, "POINT_DATA %d\n", npoints); 
            fprintf(fp, "VECTORS normal float\n");

            /* write data */
            int i;
            for (i = 0; i < npoints; i++)
                fprintf(fp, "%16.9e %16.9e %16.9e\n", (float)(buffer[nv*i+j+k]),(float)(buffer[nv*i+j+l]),(float)(0));
            fclose(fp);
        }
      }

      if (xdmf) {
        /* write h5 */
        double time = (nseries-1)*Getd("2douttime");

        #ifdef XDMF
        if(GetvLax("xdmfgrid","regular")){
          write_xdmf2D(VarName(iv[j]),"2d", suffix, level->l, time,
                      npoints, buffer, nv, j, dbl, flt, text, binary,npoints1, npoints2, o1, o2, d1, d2);

          /* add vector data to a new file, if the first component is write-active */
          if (strcmp(suffix,"xy")==0) { k=0; l=1; }
          else if (strcmp(suffix,"xz")==0) { k=0; l=2; }
          else if (strcmp(suffix,"yz")==0) { k=1; l=2; }
          else k=-1;
          if ((Getv("2doutputall","yes")) && (VarComponent(iv[j])==0) && (VarNComponents(iv[j]) == 3) && (k!=-1)) {
              char varname[1000];
              snprintf(varname,1000,"%s",VarName(iv[j]));
              varname[strlen(varname)-1]='\0';
              
              write_xdmf2D_vec(varname,"2d", suffix, level->l, time,
                      npoints, buffer, nv, j, dbl, flt, text, binary,npoints1, npoints2, o1, o2, d1, d2,k,l);
          }
        } else if(GetvLax("xdmfgrid","curvilinear")){
          double *sliced_buffer = dmalloc(npoints);
          for (i = 0; i < npoints; i++){
            sliced_buffer[i] = buffer[nv*i+j];
          }


          if(plane == 0)
            write_xdmf3D_curvilinear(VarName(iv[j]), "2d", suffix, level->l, time,
                        npoints, sliced_buffer, dbl, flt, text, binary,
                        npoints1, npoints2, 1, ox, oy, oz, dx, dy, dz);
          else if(plane == 1)
            write_xdmf3D_curvilinear(VarName(iv[j]), "2d", suffix, level->l, time,
                        npoints, sliced_buffer, dbl, flt, text, binary,
                        npoints1, 1, npoints2, ox, oy, oz, dx, dy, dz);
          else if(plane == 2)
            write_xdmf3D_curvilinear(VarName(iv[j]), "2d", suffix, level->l, time,
                        npoints, sliced_buffer, dbl, flt, text, binary,
                        1, npoints1, npoints2, ox, oy, oz, dx, dy, dz);
          else
            errorexit("diagonal grid not supported with XDMF. "
                      "Need to use unstructured mesh\n");

          free(sliced_buffer); sliced_buffer = NULL;


          /* add vector data to a new file, if the first component is write-active */
          if (Getv("2doutputall","yes")  && VarComponent(iv[j])==0 && VarNComponents(iv[j]) == 3) {

            sliced_buffer = dmalloc(3*npoints);
            for(size_t k = 0; k<3; k++) {
              for (i = 0; i < npoints; i++){
                sliced_buffer[k*npoints+i] = buffer[nv*i+j+k];
              }
            }

            char varname[1000];
            sprintf(varname,"%s",VarName(iv[j]));
            // removing component indicator from varname
            varname[strlen(varname)-1]='\0';

            if(plane == 0)
              write_xdmf3D_curvilinear_vec(varname, "2d", suffix, level->l, time,
                          npoints, sliced_buffer, dbl, flt, text, binary,
                          npoints1, npoints2, 1, ox, oy, oz, dx, dy, dz);
            else if(plane == 1)
              write_xdmf3D_curvilinear_vec(varname, "2d", suffix, level->l, time,
                          npoints, sliced_buffer, dbl, flt, text, binary,
                          npoints1, 1, npoints2, ox, oy, oz, dx, dy, dz);
            else if(plane == 2)
              write_xdmf3D_curvilinear_vec(varname, "2d", suffix, level->l, time,
                          npoints, sliced_buffer, dbl, flt, text, binary,
                          1, npoints1, npoints2, ox, oy, oz, dx, dy, dz);
            else
              errorexit("diagonal grid not supported with XDMF. "
                        "Need to use unstructured mesh\n");

            free(sliced_buffer); sliced_buffer = NULL;
          } 
        } else if(GetvLax("xdmfgrid","rectilinear")){
          double *sliced_buffer = dmalloc(npoints);
          for (i = 0; i < npoints; i++){
            sliced_buffer[i] = buffer[nv*i+j];
          }


          if(plane == 0)
            write_xdmf3D_rectilinear(VarName(iv[j]), "2d", suffix, level->l, time,
                        npoints, sliced_buffer, dbl, flt, text, binary,
                        npoints1, npoints2, 1, ox, oy, oz, dx, dy, dz);
          else if(plane == 1)
            write_xdmf3D_rectilinear(VarName(iv[j]), "2d", suffix, level->l, time,
                        npoints, sliced_buffer, dbl, flt, text, binary,
                        npoints1, 1, npoints2, ox, oy, oz, dx, dy, dz);
          else if(plane == 2)
            write_xdmf3D_rectilinear(VarName(iv[j]), "2d", suffix, level->l, time,
                        npoints, sliced_buffer, dbl, flt, text, binary,
                        1, npoints1, npoints2, ox, oy, oz, dx, dy, dz);
          else
            errorexit("diagonal grid not supported with XDMF. "
                      "Need to use unstructured mesh\n");

          free(sliced_buffer); sliced_buffer = NULL;


          /* add vector data to a new file, if the first component is write-active */
          if (Getv("2doutputall","yes")  && VarComponent(iv[j])==0 && VarNComponents(iv[j]) == 3) {

            sliced_buffer = dmalloc(3*npoints);
            for(size_t k = 0; k<3; k++) {
              for (i = 0; i < npoints; i++){
                sliced_buffer[k*npoints+i] = buffer[nv*i+j+k];
              }
            }
  
            char varname[1000];
            sprintf(varname,"%s",VarName(iv[j]));
            // removing component indicator from varname
            varname[strlen(varname)-1]='\0';

            if(plane == 0)
              write_xdmf3D_rectilinear_vec(varname, "2d", suffix, level->l, time,
                          npoints, sliced_buffer, dbl, flt, text, binary,
                          npoints1, npoints2, 1, ox, oy, oz, dx, dy, dz);
            else if(plane == 1)
              write_xdmf3D_rectilinear_vec(varname, "2d", suffix, level->l, time,
                          npoints, sliced_buffer, dbl, flt, text, binary,
                          npoints1, 1, npoints2, ox, oy, oz, dx, dy, dz);
            else if(plane == 2)
              write_xdmf3D_rectilinear_vec(varname, "2d", suffix, level->l, time,
                          npoints, sliced_buffer, dbl, flt, text, binary,
                          1, npoints1, npoints2, ox, oy, oz, dx, dy, dz);
            else
              errorexit("diagonal grid not supported with XDMF. "
                        "Need to use unstructured mesh\n");

            free(sliced_buffer); sliced_buffer = NULL;
          }
        } else errorexit("XDMF grid format is not implemented!");
        #else 
        errorexit("XDMF library and flag must be installed to use xdmf output format!");
        #endif
      }
    }
  }

  /* clean up */
  if (processor0) free(buffer);
  free(buf);
  free(coords);

  /* wait, could be optimized across different writes */
  bampi_barrier();
}



