/* outputFMR.c */
/* Bernd Bruegmann 4/04 
   OpenDX format based on Nina Jansen's perl scripts
   see 3Ddump2dx.perl for cell centered, 
   http://opendx.npaci.edu/docs/html/pages/usrgu068.htm#HDREDF

   Obsolete/defunct (BB 2/06): was only weakly tested, 
     has not been maintained since moving to VTK
*/


#include "bam.h"
#include "output.h"

/* make complete variable index list for output, all levels
   - merge two lists (say 1doutput and 1doutputx)
   - expand single components to all components if allflag is nonzero 
   - ignore variables without storage
*/


void makeoutputlist_all(tG *g, char *out0, char *out1, int allflag, 
		    int *pnindex, int **pindex)
{
  int pr = 0;
  int nindex = 0;
  int *index = malloc(sizeof(int)*globalnvariables);
  char *s[2], *name;
  int nstring;
  int i, j, k, l, n;
  int lmax = g->lmax;
  int lmin = g->lmin;
  tL *lvl;

  if (pr) {
    printf("out0 = %s\n", out0);
    printf("out1 = %s\n", out1);
  }

  /* for all lists of variable names */
  s[0] = out0;
  s[1] = out1;

  for (nstring = 0; nstring < 2; nstring++) {

    /* for each name */
    while (name = NextEntry(s[nstring])) {

      /* find range of variable indices */
      i = IndLax(name);
      if (i < 0) continue;
      n = 1;
      if (allflag) {
	i = IndComponent0(i);
	n = VarNComponents(i);
      }

      /*for all levels. this assumes a fixed grid 
	and that the grid has already been allocated*/
      
      for (l = lmin; l <= lmax; l++) {
	lvl = g->level[l];
	/* for all variable indices */
	for (j = i; j < i+n; j++) {
	  /* ignore variable if it is constant in time and iteration > 0 */
	  if (VarConstantFlag(j) && lvl->iteration) continue;

	  /* check whether index is already in the list */
	  for (k = 0; k < nindex; k++)
	    if (index[k] == j) break;
	  /* add index if not already in list */
	  if (k == nindex) 
	    index[nindex++] = j;
	}
      }
    }
  }

  /* return list of indices */
  *pnindex = nindex;
  *pindex = index;
  
  if (pr) {
    printf("nindex = %d, vars = \n", nindex); 
    for (i = 0; i < nindex; i++)
      printf(" %s", VarName(index[i]));
    printf("\n");
  }
}


int write_FMR(tL *level)
{
  tG *g = level->grid;
  static int di;
  static double dt;
  static char *ou;
  static int all;
  char *name, s[20];
  int d, l;
  int check, tfo;
  int nindex, *index;
  int lmax = g->lmax;
  int lmin = g->lmin;

  tL *lvl;
  if (!Getv("3dformat", "amrunion")) return 0;
 
  /* cache */

  sprintf(s, "3doutiter");
  di = Geti(s);
  sprintf(s, "3douttime");
  dt = Getd(s);
  sprintf(s, "3doutput");
  ou = Gets(s);
  sprintf(s, "3doutputall");
  all = Getv(s, "yes");

  if (Getv("outputall", "yes")) all = 1;

  /* 3d output */

  makeoutputlist_all(g, ou, "", all, &nindex, &index);

  /* assumes that time for output for all levels is the same as 
     time for output of lmax */

  tfo = 1;

  for (l = lmin; l <= lmax; l++) {
    lvl = g->level[l];    
    
    if (timeforoutput_di_dt(lvl, di, dt)) {
      check = 1; 
    }
    else {
      check = 0;     
    }
    
    tfo = tfo*check;      
  }
  
  if (tfo) {
    write_3d_amrunion(g, nindex, index);
  }
   
  free(index);
  return 0;
}



/* make list of ordered points for the global grid
   this function is used for compatibility with the standard 3d output
*/
void makeorderedlist3dFMR(tL *level, int interp, int *npoints, double **coords)
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
    x1 -= dx/2;
    y1 -= dy/2;
    z1 -= dz/2;

    /* disregard outermost points if not on outermost grid */
    if (level->l > 0) {
      x0 += dx;
      y0 += dy;
      z0 += dz;
      x1 -= dx;
      y1 -= dy;
      z1 -= dz;
    }

    /* symmetries */
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

void write_3d_amrunion(tG *g, int nv, int *iv)
{

  int pr = 0;
  int nseries;
  int interp = Getv("3doutinterpolate", "yes");
  double dtseries;
  static int ncall = 0;
  static int firstcall = 1;
  static int outnpoints = 0;
  FILE *fp, *fp1, *fp2, *fp3;
  char filename[1000], filename1fo[1000], filename2fo[1000];
  char filename1[1000],filename2[1000],filename3[1000];
  int lmax = g->lmax;
  int lmin = g->lmin;
  int i,j,k,l,lev, ind, unique,n;
  int npoints, allnpoints, oldnpoints;
  double *coords, *allcoords, *allflags, *flagregrid, *outcoords;
  double **corners, **ucorners;
  int *alllevels, *outlevels, *index2, *index3;
  double dist;
  double dx = g->level[lmin]->dx;
  double dy = g->level[lmin]->dy;
  double dz = g->level[lmin]->dz;
  int index[1];
  int nbuf, nbuffer;
  double *buf = 0, *buffer = 0;
  double symsign;
  double *timestamp;

  tL *lvl, *level;

  if (nv <= 0) return;
  if (bampi_size() > 1)
    errorexit("cannot output 3d data as amrunion on multiple processors yet");

  ncall++;

  if (firstcall ) {
    firstcall = 0;
    if ( processor0) {
      timestamp = (double *) malloc(sizeof(double));
      snprintf(filename, 1000, "%s/lxyz.xyz", Gets("outdir"));
      fp = fopen(filename, "a");      
      npoints = 0;
      allnpoints = 0;
      allcoords = (double *) malloc(sizeof(double));
      alllevels = (int *) malloc(sizeof(int));
      allflags  = (double *) malloc(sizeof(double));
      
      if (!allcoords) 
	errorexit("Out of memory in output3d.");
      
      for (l = lmin; l <= lmax; l++) {
	lvl = g->level[l];    
	flagregrid = lvl->v[Ind("flagregrid")];
	oldnpoints = allnpoints;
	makeorderedlist3dFMR(lvl, interp, &npoints, &coords);
	allnpoints += npoints;
	allcoords = (double *) realloc(allcoords, sizeof(double)*3*allnpoints);
	if (!allcoords) 
	  errorexit("Out of memory in output3d.");
	alllevels = (int *) realloc(alllevels, sizeof(int)*allnpoints);
	if (!alllevels) 
	  errorexit("Out of memory in output3d.");
	allflags =  (double *) realloc(allflags,  sizeof(double)*allnpoints);
	if (!allflags) 
	  errorexit("Out of memory in output3d.");
	for (i = 0; i < 3*npoints; i++)
	  allcoords[i+3*oldnpoints] = coords[i];           
	for (i = 0; i < npoints; i++)
	  alllevels[i+oldnpoints] = l;     
	for (i = 0; i < npoints; i++)
	  allflags[i+oldnpoints] = flagregrid[i];     
      }

      /* making a list of all points to be output */
      
      outcoords = (double *) malloc(sizeof(double)*3*allnpoints);
      if (!outcoords) 
	errorexit("Out of memory in output3d.");
      outlevels = (int *) malloc(sizeof(int)*allnpoints);
      if (!outlevels) 
	errorexit("Out of memory in output3d.");
      
      
      for (i = 0; i < allnpoints; i++) {
	if (!allflags[i]) { 
	  outcoords[3*outnpoints] = allcoords[3*i];
	  outcoords[3*outnpoints + 1] = allcoords[3*i+1];
	  outcoords[3*outnpoints + 2] = allcoords[3*i+2];
	  outlevels[outnpoints] = alllevels[i];
	  outnpoints++;
	}
      }
      free(allcoords);
      free(alllevels);
      free(allflags);
      
      outcoords = (double *) realloc(outcoords, sizeof(double)*3*outnpoints);
      if (!outcoords) 
	errorexit("Out of memory in output3d.");
      outlevels = (int *) realloc(outlevels, sizeof(int)*outnpoints);
      if (!outlevels) 
	errorexit("Out of memory in output3d.");
 
      corners = (double **) malloc(sizeof(double *)*outnpoints*8);
      
      for (i = 0; i < outnpoints*8; i++) {
	corners[i] = (double *) malloc(sizeof(double) * 4);
      }
      
      
      for (i = 0; i < outnpoints; i++) { 
	dist = 0.5*dx/(pow(2,outlevels[i]));
	
	corners[i*8][0] = i*8;
	corners[i*8][1] = outcoords[3*i]   - dist;
	corners[i*8][2] = outcoords[3*i+1] - dist;
	corners[i*8][3] = outcoords[3*i+2] - dist;
	
	corners[i*8+1][0] = i*8+1;
	corners[i*8+1][1] = outcoords[3*i]   - dist;
	corners[i*8+1][2] = outcoords[3*i+1] - dist;
	corners[i*8+1][3] = outcoords[3*i+2] + dist;
	
	corners[i*8+2][0] = i*8+2;
	corners[i*8+2][1] = outcoords[3*i]   - dist;
	corners[i*8+2][2] = outcoords[3*i+1] + dist;
	corners[i*8+2][3] = outcoords[3*i+2] - dist;

	corners[i*8+3][0] = i*8+3;
	corners[i*8+3][1] = outcoords[3*i]   - dist;
	corners[i*8+3][2] = outcoords[3*i+1] + dist;
	corners[i*8+3][3] = outcoords[3*i+2] + dist;

	corners[i*8+4][0] = i*8+4;
	corners[i*8+4][1] = outcoords[3*i]   + dist;
	corners[i*8+4][2] = outcoords[3*i+1] - dist;
	corners[i*8+4][3] = outcoords[3*i+2] - dist;

	corners[i*8+5][0] = i*8+5;
	corners[i*8+5][1] = outcoords[3*i]   + dist;
	corners[i*8+5][2] = outcoords[3*i+1] - dist;
	corners[i*8+5][3] = outcoords[3*i+2] + dist;

	corners[i*8+6][0] = i*8+6;
	corners[i*8+6][1] = outcoords[3*i]   + dist;
	corners[i*8+6][2] = outcoords[3*i+1] + dist;
	corners[i*8+6][3] = outcoords[3*i+2] - dist;

	corners[i*8+7][0] = i*8+7;
	corners[i*8+7][1] = outcoords[3*i]   + dist;
	corners[i*8+7][2] = outcoords[3*i+1] + dist;
	corners[i*8+7][3] = outcoords[3*i+2] + dist;
      }
      
      qsort(corners,outnpoints*8,sizeof( corners[0]),(int (*)()) &CoordComp);

      index2 = (int *) malloc(sizeof(int)*8*outnpoints);

      /* making list of indecess */
      for (i = 0; i < outnpoints*8 ; i++) {
	ind = corners[i][0];
	index2[ind] = i;
      }

      /* making array with unique corners */

      index3 = (int *) malloc(sizeof(int)*8*outnpoints);

      ucorners = (double **) malloc(sizeof(double *)*outnpoints*8);
      
      for (i = 0; i < outnpoints*8; i++) {
	ucorners[i] = (double *) malloc(sizeof(double) * 3);
      }

      index3[0] = 0;
      unique = 0;
      
      for ( k = 0; k < 3; k++) {
	ucorners[0][k] = corners[0][k+1];
      }


      for (i = 1; i < 8*outnpoints ; i++) {
	if (!dequal(corners[i][1],corners[i-1][1]) || !dequal(corners[i][2],corners[i-1][2]) || !dequal(corners[i][3],corners[i-1][3])) {
	  unique++;
	  for ( k = 0; k < 3; k++) {
	    ucorners[unique][k] = corners[i][k+1];
	  }
	}
	ind = corners[i][0];
	index3[ind] = unique;
      }

      fprintf(fp,"#object1 contains the positions \n\n");
      fprintf(fp,"object 1 class array type float rank 1 shape 3 items %d data follows\n", unique+1);

      for (i = 0; i < unique+1 ; i++) {
	fprintf(fp, "%.10f %.10f %.10f\n",ucorners[i][0],ucorners[i][1],ucorners[i][2]);
      }

      fprintf(fp, "\n\n");
      fprintf(fp, "#object2 contains the connections \n\n");
      fprintf(fp, "object 2 class array type int rank 1 shape 8 items %d data follows\n", outnpoints);

      for (i = 0; i < outnpoints ; i++) {
	fprintf (fp, "%d %d %d %d %d %d %d %d\n", index3[8*i], index3[8*i+1], index3[8*i+2], index3[8*i+3], index3[8*i+4], index3[8*i+5], index3[8*i+6], index3[8*i+7]);
      }
	
      fprintf(fp,"attribute \"element type\" string \"cubes\"\n");
      fprintf(fp,"attribute \"ref\" string \"positions\"\n\n");

      for (i = 0; i < outnpoints*8 ; i++) {
	free(corners[i]);
      }

      free(corners);

      for (i = 0; i < outnpoints*8 ; i++) {
	free(ucorners[i]);
      }

      free(ucorners);    
    
      free(index2);
      free(index3);

      fclose(fp);      
    }

    // defunct: what is n?!
    n = 0;
    snprintf(filename3, 1000, "%s/%s.dx", Gets("outdir"), VarName(iv[n]));    
  }


  /* write file containing data */
  
  /* for each level */
  for (l = lmin; l <= lmax; l++) {
    level = g->level[l];
    makeorderedlist3dFMR(level, interp, &npoints, &coords);
    flagregrid = level->v[Ind("flagregrid")];
    if (processor0) {
      /* allocate buffer for data on processor 0 */      
      nbuffer = npoints;
      buffer = dmalloc(nbuffer);
    }
    /* for each variable */
    for (n = 0; n < nv; n++) {
      index[0] = iv[n];
      
      /* fill buffer with data at given list of coordinates, allocates buf */
      bufferorderedpoints(level, npoints, coords, 1, index, &nbuf, &buf, interp);
      
      /* combine data from different processors */
      bampi_combinepointbuffers(level, nbuf, buf, buffer, 1);
      
      /* processor 0 does the writing */
      if (processor0) {	
	snprintf(filename1fo, 1000, "%s_data",VarName(iv[n]));
	snprintf(filename2fo, 1000, "%s_fields",VarName(iv[n]));

	/* filename */
	snprintf(filename1, 1000, "%s/%s_data", 
		 Gets("outdir"), VarName(iv[n]));
	
	/* open file */
	fp1 = fopen(filename1, "a");
	if (!fp1) errorexits("failed opening %s", filename1);


	snprintf(filename2, 1000, "%s/%s_fields", 
		 Gets("outdir"), VarName(iv[n]));
	
	/* open file */
	fp2 = fopen(filename2, "a");
	if (!fp2) errorexits("failed opening %s", filename2);


	snprintf(filename3, 1000, "%s/%s.dx", 
		 Gets("outdir"), VarName(iv[n]));
	
	/* open file */
	fp3 = fopen(filename3, "a");
	if (!fp3) errorexits("failed opening %s", filename3);	

	if (level->l == 0) {
	  timestamp = (double *) realloc(timestamp, sizeof(double)*ncall);
	  timestamp[ncall -1] = level->time;
	  fprintf(fp1,  "object %d class array type float rank 0 items %d data follows\n", ncall, outnpoints);
	}

	/* write data */
	if (Getv("3dformat", "add_rotant_points")) {
	  symsign = VarSymmetry(iv[n], 0) * VarSymmetry(iv[n], 1);
	  for (i = 0; i < npoints; i++)
	    if ((!flagregrid || !flagregrid[i]) && !dless(coords[3*i+1], 0)) {
	      fprintf(fp1, "%.9e\n", buffer[i]);
	      if (!dequal(coords[3*i+1], 0))
		fprintf(fp1, "%.9e\n", buffer[i] * symsign);
	    }
	}   
	
	else {
	  for (i = 0; i < npoints; i++)
	    if (!flagregrid || !flagregrid[i])  
	      fprintf(fp1, "%.9e\n", buffer[i]);
	}


	if (l == lmax) fprintf(fp1, "attribute \"dep\" string \"connections\"\n\n");
	
	if (level->l == 0) {
	  fprintf(fp2, "object %d class field\n", ncall);
	  fprintf(fp2, "component \"positions\" value file \"lxyz.xyz\",1\n");
	  fprintf(fp2, "component \"connections\" value file \"lxyz.xyz\",2\n");
	  fprintf(fp2, "component \"data\" value file %s,%d\n\n",filename1fo,ncall);      
	      
	  if (ncall == 1) fprintf(fp3, "object \"series\" class series\n");

	  fprintf(fp3, "member %d value file %s,%d position %.9e\n",ncall-1,filename2fo,ncall,timestamp[ncall-1]);
	}

	
        /* close files */
	fclose(fp1);
	fclose(fp2);
	fclose(fp3);
      }
      /* wait, could be optimized across different writes */
      bampi_barrier();
    }

    /* clean up */
    if (processor0) free(buffer);
    free(buf);
    free(coords);
  }
}




int CoordComp( double **p1, double **p2 ){

  double x1,x2,y1,y2,z1,z2;

  x1 = p1[0][1];
  x2 = p2[0][1];
  y1 = p1[0][2];
  y2 = p2[0][2];
  z1 = p1[0][3];
  z2 = p2[0][3];


  if (x1 < x2) {
    return(-1);
  }
  else if (x1 > x2) {
    return(1);
  }
  else if (x1 == x2) {
    if (y1 < y2) {
      return(-1);
    }
    else if (y1 > y2) {
      return(1);
    }
    else if (y1 == y2) {
      if (z1 < z2) {
	return(-1); 
      }
      else if (z1 > z2) {
	return(1);
      }
      else return(0);
    }
  }
  return(0);
}
