/* outputshells.c */
/* mth 04/11 */

#include "bam.h"
#include "output.h"

#define PR 0








void write_level_shells1d(tL* level, int nv, int *iv, int l, char *suffix)
{
  tB* box = level->box[0];
  
  /* adjust this for more output ... and adjust everything below */
  if (!((l==1 && boxnr==0) ||
        (l==1 && boxnr==1) ||
        (l==2 && boxnr==2) ||
        (l==2 && boxnr==3) ||
        (l==3 && boxnr==4) ||
        (l==3 && boxnr==5) ||
        (l==0 && boxnr==0) ||
        (l==0 && boxnr==1))) return;
  
  
  FILE *fp;
  char filename[1000];
  int flag;
  int i,j,k, ii,jj,kk, n,npoints;
  double b;
  
  int mm = level->ibbox[1]-level->ibbox[0]+1;
  int nn = level->ibbox[3]-level->ibbox[2]+1;
  int oo = level->ibbox[5]-level->ibbox[4]+1; 
  double *rp = level->v[Ind("shells_r")];
  double *pp = level->v[Ind("shells_phi")];
  double *tp = level->v[Ind("shells_theta")];
  int nbuf   = Geti("amr_nbuffer");
  int nbuf_r = Geti("amr_nbuffer");
  int nbuf_p = Geti("bampi_nghosts");
  
  npoints = mm;
  double *coords = (double*) malloc (npoints*sizeof(double));
  double *buffer = (double*) malloc ((nv+1)*(npoints)*sizeof(double));
  double vinterp[nv];
  
  for (ii=0; ii<npoints; ii++) {
    coords[ii] = level->bbox[0] + ii*level->dx;
    if (boxnr%2==0) 
      coords[ii] = 2*level->bbox[0] + (npoints-1)*level->dx - coords[ii];
    
    flag = interpolate_xyz_localinbox_minimal(box, coords[ii],(l==0)*PI/4.,(l==0)*PI/4., nv,iv,vinterp, output_order, output_scheme);
    
    /* if flag interpolate */
    if (flag) {
      for (n=0; n<nv; n++)
        buffer[ii*(nv+1)+n] = vinterp[n];
      buffer[ii*(nv+1)+nv] = 1.;
    } else {
      for (n=0; n<=nv; n++)
        buffer[ii*(nv+1)+n] = 0.;
    }
  }
  
  /* sync buffer -> all interpolated values are then available on all processors */
  if (bampi_size() > 1) {
    double *global = (double*) malloc ((nv+1)*(npoints)*sizeof(double));
    bampi_allreduce_sum_vector(buffer, global, (nv+1)*(npoints));
    for (ii=0; ii<npoints; ii++) {
      if (global[ii*(nv+1)+nv]==0.) {
        global[ii*(nv+1)+n] = 0.;
        //errorexit("box->shell interpolate: this point is nowhere");
      } else {//OPTIMIZE
        for (n=0; n<nv; n++)
          global[ii*(nv+1)+n] /= global[ii*(nv+1)+nv];
      }
    }
    free(buffer);
    buffer = global;
  }
  
  
  if (processor0) {
    /* for each variable */
    for (n = 0; n < nv; n++) {
  
      /* filenames */
      snprintf(filename, 1000, "%s/%s.%s%d", 
              Gets("outdir_1d"), VarName(iv[n]), suffix, level->l);
        
      /* open files */
      fp = fopen(filename, "a");
      if (!fp) errorexits("failed opening %s", filename);
      
      /* write info */
      if (boxnr%2==0) 
        fprintf(fp, "\"Time = %f\"\n", level->time);
        
      /* write */
      for (ii = nbuf; ii < npoints-nbuf; ii++) {
        b = buffer[ii*(nv+1)+n];
        fprintf(fp, "%22.15e  %22.15e\n", (boxnr%2==0?-1.:1.)*compute_fisheye(coords[ii],0), b);
      }
  
      /* close file */
      if (boxnr%2!=0)
        fprintf(fp,"\n");
      fclose(fp);
    }
  }
 
  free(buffer);
  free(coords);
}

void write_level_shells2d(tL* level, int nv, int *iv, int l, char *suffix)
{
  /* adjust this for more output ... and adjust everything below */
  //if (!(l==0 || l==1 || l==2 )) return;
  
  
  FILE *fp;
  char filename[1000];
  int flag;
  int i,j,k, ii,jj,kk, iijj, n,npoints,b;
  double r,p,t, x,y,z;
  
  int mm = level->ibbox[1]-level->ibbox[0]+1;
  int nn = level->ibbox[3]-level->ibbox[2]+1;
  int oo = level->ibbox[5]-level->ibbox[4]+1; 
  double *rp = level->v[Ind("shells_r")];
  double *pp = level->v[Ind("shells_phi")];
  double *tp = level->v[Ind("shells_theta")];
  int order  = Geti("order_RP");
  int nbuf   = Geti("amr_nbuffer");
  int nbuf_r = Geti("amr_nbuffer");
  int nbuf_p = Geti("bampi_nghosts");
  
  int nr = mm - 2*nbuf_r;
  int np = 4* (nn-2*nbuf_p-1)+1;
  
  npoints = nr*np;
  double *coords = (double*) malloc (3*npoints*sizeof(double));
  double *buffer = (double*) malloc ((nv+1)*(npoints)*sizeof(double));
  double vinterp[nv];
  
  int text   = Getv("2dformat", "text");
  int binary = Getv("2dformat", "binary");
  int flt    = Getv("2dformat", "float");
  int dbl    = Getv("2dformat", "double");
  
  iijj = 0;
  for (jj=0; jj<np; jj++) {
    for (ii=0; ii<nr; ii++) {
      
      r = level->bbox[0] + (nbuf_r+ii)*level->dx;
      p = level->bbox[2] + (nbuf_p+jj)*level->dy;
      
      r = compute_fisheye(r,0);
      
      if (l==0) { //xy
        x = r*cos(p);
        y = r*sin(p);
        z = 0.;
        b = find_shellsbox_from_xyz(x,y,z);
      } else if (l==1) { //xz
        x = r*cos(p);
        y = 0.;
        z = r*sin(p);
        b = find_shellsbox_from_xyz(x,y,z);
      } else if (l==2) { //yz
        x = 0.;
        y = r*sin(p);
        z = r*cos(p);
        b = find_shellsbox_from_xyz(x,y,z);
      } else if (l==3) { //xd
        x = r*cos(p);
        y = r*sin(p)/sqrt(2);
        z = r*sin(p)/sqrt(2);
        b = find_shellsbox_from_xyz(x,y,z);
      } else if (l==4) { //yd
        x = r*sin(p)/sqrt(2);
        y = r*cos(p);
        z = r*sin(p)/sqrt(2);
        b = find_shellsbox_from_xyz(x,y,z);
      } else if (l==5) { //zd
        x = r*sin(p)/sqrt(2);
        y = r*sin(p)/sqrt(2);
        z = r*cos(p);
        b = find_shellsbox_from_xyz(x,y,z);
      }
      
      convert_box_to_shells( b, x,y,z, &r,&p,&t);
     
      flag = interpolate_xyz_localinbox_minimal( level->box[b],
              r,p,t, nv,iv,vinterp, output_order,output_scheme);
      
      coords[iijj+0*npoints] = x;
      coords[iijj+1*npoints] = y;
      coords[iijj+2*npoints] = z;
      
      if (flag) {
        for (n=0; n<nv; n++)
          buffer[iijj*(nv+1)+n] = vinterp[n];
        buffer[iijj*(nv+1)+nv] = 1.;
      } else {
        for (n=0; n<=nv; n++)
          buffer[iijj*(nv+1)+n] = 0.;
      }
      
      iijj++;
    }
  }
  
  /*
  yo();
  printf(" %d %d %d    %d\n", nv,(npoints),sizeof(double), ((nv+1)*(npoints))*sizeof(double));
  */
  
  /* sync buffer -> all interpolated values are then available on all processors */
  if (bampi_size() > 1) {
    double* global = (double*) malloc (((nv+1)*(npoints))*sizeof(double));
    bampi_allreduce_sum_vector(buffer, global, (nv+1)*(npoints));
    for (ii=0; ii<npoints; ii++) {
      if (global[ii*(nv+1)+nv]==0.) {
        global[ii*(nv+1)+n] = 0.;
        //errorexit("box->shell interpolate: this point is nowhere");
      } else {//OPTIMIZE
        for (n=0; n<nv; n++)
          global[ii*(nv+1)+n] /= global[ii*(nv+1)+nv];
      }
    }
    free(buffer);
    buffer = global;
  }
  
  int nseries = timeforoutput_di_dt(level, Geti("2doutiter"), Getd("2douttime"));
  if (processor0) {
    /* for each variable */
    for (n = 0; n < nv; n++) {
  
      /* open files */
      fp = fopen_vtk(VarName(iv[n]), "2d", suffix, level->l, nseries-1);
      
      fprintf(fp, "# vtk DataFile Version 2.0\n");
      fprintf(fp, "variable %s, level %d, time %16.9e\n", 
              VarName(iv[n]), level->l, level->time);
      fprintf(fp, binary ? "BINARY" : "ASCII\n");
      fprintf(fp, "\n");
      fprintf(fp, "DATASET STRUCTURED_GRID\n");
      fprintf(fp, "DIMENSIONS %d %d %d\n", nr, np, 1);
      fprintf(fp, "POINTS %d %s\n", npoints, dbl ? "double" : "float");

      /* write points */
      write_raw_vec_vtk(fp, npoints, coords, dbl,flt,text,binary);
      
      fprintf(fp, "POINT_DATA %d\n",npoints);
      fprintf(fp, "SCALARS scalars %s\n", dbl ? "double" : "float");
      fprintf(fp, "LOOKUP_TABLE default\n");

      /* write data */
      write_raw_vtk(fp, npoints, buffer, nv+1,n, dbl,flt,text,binary);
      
      /* close file */
      fprintf(fp,"\n");
      fclose(fp);
    }
  }
 
  free(buffer);
  free(coords);
}

void write_level_shells3d(tL* level, int nv, int *iv, char *suffix)
{
  int var_x = Ind("x");
  int var_y = Ind("y");
  int var_z = Ind("z");
  
  double x,y,z;
  FILE *fp;
  char filename[1000];
  
  int nseries = timeforoutput_di_dt(level, Geti("3doutiter"), Getd("3douttime"));
  
  int i,j,k, ii,jj,kk, iijjkk, n;
  int mm = level->ibbox[1]-level->ibbox[0]+1;
  int nn = level->ibbox[3]-level->ibbox[2]+1;
  int oo = level->ibbox[5]-level->ibbox[4]+1;
  int npoints = mm*nn*oo;
  tB* box = level->box[0];
  double *rp = level->v[Ind("shells_r")];
  double *pp = level->v[Ind("shells_phi")];
  double *tp = level->v[Ind("shells_theta")];
  double *buffer,*global;
  
  int text   = Getv("3dformat", "text");
  int binary = Getv("3dformat", "binary");
  int flt    = Getv("3dformat", "float");
  int dbl    = Getv("3dformat", "double");
  
  /* for each variable */
  for (n = 0; n < nv; n++) {
    
    if (PR) printf("  output shells (box%d):  %d(%d)  %s\n", level->box[0]->i,n,nv,VarName(iv[n]));
    buffer = (double*) dcalloc (6*npoints);
    
    /* go through all global points and add information if this
        point exist in this processor */
    forallpoints_boxijk(box) {
      
      ii = i + floor( (box->com->bbox[0]-box->bbox[0])/box->dx+.0001 );
      jj = j + floor( (box->com->bbox[2]-box->bbox[2])/box->dy+.0001 );
      kk = k + floor( (box->com->bbox[4]-box->bbox[4])/box->dz+.0001 ); 
      iijjkk = mm*nn*kk + mm*jj + ii;
      if (iijjkk>=npoints)errorexit("something wrong");
      
      convert_shells_to_box(boxnr, rp[ijk],pp[ijk],tp[ijk],
                            &buffer[0*npoints+iijjkk],
                            &buffer[1*npoints+iijjkk],
                            &buffer[2*npoints+iijjkk]);
      
      if (xyzinsidebbox(box->com->bbox, rp[ijk],pp[ijk],tp[ijk])) {
        buffer[3*npoints + iijjkk] = level->v[iv[n]][ijk];
        buffer[4*npoints + iijjkk] = 1.;
      } else {
        buffer[0*npoints + iijjkk] = 0.;
        buffer[1*npoints + iijjkk] = 0.;
        buffer[2*npoints + iijjkk] = 0.;
        buffer[3*npoints + iijjkk] = 0.;
        buffer[4*npoints + iijjkk] = 0.;
      }
      
    } endfor_ijk;
    
    
    /* sync buffer -> all interpolated values are then available on all processors */
    if (bampi_size() > 1) {
      global = (double*) malloc (6*npoints*sizeof(double));
      bampi_allreduce_sum_vector(buffer, global, 6*npoints);
      for (i=0; i<npoints; i++) {
        
        if (global[4*npoints+i]==0.) 
          errorexit("this point is nowhere");
        
        global[0*npoints+i] /= global[4*npoints+i];
        global[1*npoints+i] /= global[4*npoints+i];
        global[2*npoints+i] /= global[4*npoints+i];
        global[3*npoints+i] /= global[4*npoints+i];
       
      }
      free(buffer);
      buffer = global;
    }
  
  
    /* put the data into a vtk file */
    if (processor0) {
    
      /* open file (returns non-null file pointer) */
      fp = fopen_vtk(VarName(iv[n]), "3d", suffix, level->l, nseries-1);
      
      /* write header */
      fprintf(fp, "# vtk DataFile Version 2.0\n");
      fprintf(fp, "variable %s, level %d, time %16.9e\n", 
              VarName(iv[n]), level->l, level->time);
      fprintf(fp, binary ? "BINARY" : "ASCII\n");
      fprintf(fp, "\n");
      fprintf(fp, "DATASET STRUCTURED_GRID\n");
      fprintf(fp, "DIMENSIONS %d %d %d\n", mm, nn, oo);
      fprintf(fp, "POINTS %d %s\n", npoints, dbl ? "double" : "float");

      /* write points */
      write_raw_vec_vtk(fp, npoints, buffer, dbl,flt,text,binary);
  
      fprintf(fp, "POINT_DATA %d\n",npoints);
      fprintf(fp, "SCALARS scalars %s\n", dbl ? "double" : "float");
      fprintf(fp, "LOOKUP_TABLE default\n");

      /* write data */
      write_raw_vtk(fp, npoints, &buffer[3*npoints], 1,0, dbl,flt,text,binary);
      
      /* close file */
      fclose(fp);
     
    }
    
    free(buffer);
  }
}










void write_level_shells3d_surface(tL* level, int nv, int *iv)
{
  if (!Getv("grid","shells") || level->l!=level->grid->lmin)
    return;
  
  if (0) printf("write_level_shells2d_surface  on level %d \n", level->l);
  
  FILE *fp;
  char filename[1000];
  
  int i,j,k,ijk,  ii,jj,kk, iijjkk, jjkk;
  int m,n,o, var;
  int nbuf_phi = Geti("order_RP_shells")/2+2;
  int mm = level->ibbox[1]-level->ibbox[0]+1;
  int nn = level->ibbox[3]-level->ibbox[2]+1-2*nbuf_phi;
  int oo = level->ibbox[5]-level->ibbox[4]+1-2*nbuf_phi;
  int npoints = 6*nn*oo;
  
  double *buffer,*global;
  double r,phi,theta;
  
  /* for each variable */
  for (var = 0; var < nv; var++) {
    
    
    /* filenames */
    snprintf(filename, 1000, "%s/%s.rpt.t%d", Gets("outdir_3d"), VarName(iv[var]), level->l);
    
    if (processor0) {
      /* open files */
      fp = fopen(filename, "a");
      if (!fp) errorexits("failed opening %s", filename);
      
      /* header */
      if (level->iteration == 0) { 
        
        fprintf(fp, "# data in boxshells mode %s\n", VarName(iv[var]));
        fprintf(fp, "# N    = %d x %d x %d\n", mm,nn,oo);
        fprintf(fp, "# Ntot = %d\n", 6*mm*nn*oo);
        fprintf(fp, "# r    = %14.6e %14.6e\n", level->bbox[0],level->bbox[1]);
      }
      fprintf(fp, "# iteration = %d\n",     level->iteration);
      fprintf(fp, "# time      = %14.6e\n", level->time);
    }
    
    
    /* find data */
    for (ii=level->ibbox[0]; ii<=level->ibbox[1]; ii++) {
      
      buffer = (double*) dcalloc (2*npoints);
      
      jjkk = 0;
      forallboxes(level) {
        
        m = box->com->ibbox[1]-box->com->ibbox[0]+1;
        n = box->com->ibbox[3]-box->com->ibbox[2]+1;
        o = box->com->ibbox[5]-box->com->ibbox[4]+1;
        
        for (kk=level->ibbox[4]+nbuf_phi; kk<=level->ibbox[5]-nbuf_phi; kk++) {
          for (jj=level->ibbox[2]+nbuf_phi; jj<=level->ibbox[3]-nbuf_phi; jj++) {
           
            r     = level->bbox[0] + ii*level->dx;
            phi   = level->bbox[2] + jj*level->dy;
            theta = level->bbox[4] + kk*level->dz;
            
            i = floor( (r    -box->com->bbox[0])/box->dx+.0001 );
            j = floor( (phi  -box->com->bbox[2])/box->dy+.0001 );
            k = floor( (theta-box->com->bbox[4])/box->dz+.0001 );
            ijk = box->noffset + m*n*(k) + m*(j) + i;
            
            if (i>=0 && j>=0 && k>=0 && i<=m && j<=n && k<=o) {
              buffer[jjkk] = level->v[ iv[var] ][ijk];
              buffer[npoints + jjkk] = 1.;
            } else {
              buffer[npoints + jjkk] = 0.;
            }
            jjkk++;
          }
        }
        
      } endforboxes;
      
      
      /* sync buffer -> all interpolated values are then available on all processors */
      if (bampi_size() > 1) {
        global = (double*) malloc (2*npoints*sizeof(double));
        bampi_allreduce_sum_vector(buffer, global, 2*npoints);
        for (i=0; i<npoints; i++) {
          
          if (global[npoints+i]==0.) 
            errorexit("this point is nowhere");
        
          global[i] /= global[npoints+i];
          
        }
        free(buffer);
        buffer = global;
      }
      
      
      /* fill the file */
      if (processor0) {
        for (jjkk=0; jjkk<npoints; jjkk++)
          fprintf(fp, "%14.6e", buffer[jjkk]);
        
      }
      
      
      free(buffer);
    }
    
    if (processor0)
      fclose(fp);
    
  }
  
  
  
  
  
  
  
  
  
}









