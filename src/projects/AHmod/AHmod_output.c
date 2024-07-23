/* AHmod_output.c */
/* mth 12/07 */
/* NL 02/09 */

#include "bam.h"
#include "AHmod.h"
#include "AHmod_output.h"

#include <sys/types.h>
#include <sys/stat.h>


// output in outdir/ah.xyz_vtk/ah.xyzMOTSNUM_TIMESLICE.vtk
// in vtk format
void AHmod_output(const int ntheta, const double dtheta, 
    const int nphi, const double dphi, const double xc, const double yc, 
    const double zc, double **rr, const char *outdir, 
    const int number, const double time, const double AHmod_time)
{
  int i,j,k,l;
  double theta, phi, x,y,z;
  char filename[300], path[1000];
  
  int slice = (int)(time/AHmod_time);
   
  // %.6d  = width=6 automatically fill with 0 
  sprintf(filename, "%s/ah.xyz_vtk/ah.xyz%d_%.6d.vtk",outdir,number,slice);
  
  FILE *fp1 = fopen(filename, "wb");
  
  /*check if file exists; if not, create file*/
  if (!fp1) {
    sprintf(path, "%s/ah.xyz_vtk",outdir);
    if (number==0) mkdir(path,0777);
    fp1 = fopen(filename, "wb");
    if (!fp1) errorexits("failed opening %s", filename);
  }
  
  
  printf("\n  Output AHmod (vtk-format)\n");  
  
  /*output header*/
  fprintf(fp1, "# vtk DataFile Version 2.0\n");
  fprintf(fp1, "Apparent Horizon %d, time=%f\n",number,time);
  fprintf(fp1, "ASCII\n\n");
  fprintf(fp1, "DATASET UNSTRUCTURED_GRID\n");
  
  /*output locations of all points*/
  fprintf(fp1, "POINTS %d float\n", nphi*ntheta);
  for(i=0;i<ntheta;i++){
    theta = AHmod_theta(dtheta, i);
    for(j=0;j<nphi;j++){
      phi = AHmod_phi(dphi, j);
      x = xc + rr[i][j]*sin(theta)*cos(phi);
      y = yc + rr[i][j]*sin(theta)*sin(phi);
      z = zc + rr[i][j]*cos(theta);
      fprintf(fp1, " %14.8f  %14.8f  %14.8f\n",x,y,z);
    }
  }
  fprintf(fp1, "\n");
    
  /*output data*/
  /*cells (nphi*ntheta squares and 2 polygons with nphi edges)*/
  k=0;
  fprintf(fp1, "CELLS %d %d\n", nphi*(ntheta-1)+2, 5*nphi*(ntheta-1)+(nphi+1)*2);  
  for(i=0;i<ntheta-1;i++){
    l=k;
    for(j=0;j<nphi;j++){      
      if (j==nphi-1) fprintf(fp1, " 4 %d %d %d %d\n",k,l,l+nphi,k+nphi);
      else fprintf(fp1, " 4 %d %d %d %d\n",k,k+1,k+nphi+1,k+nphi);
      k++;
    }
  }
  fprintf(fp1, " %d",nphi);
  for(j=0;j<nphi;j++) fprintf(fp1, " %d", j);
  fprintf(fp1, "\n %d",nphi);
  for(j=0;j<nphi;j++) fprintf(fp1, " %d", j+nphi*(ntheta-1));
  fprintf(fp1, "\n\n");
  
  /*celltypes (7=square, 9=polygon)*/
  fprintf(fp1, "CELL_TYPES %d\n", nphi*(ntheta-1)+2);
  for(i=0;i<(ntheta-1)*nphi;i++)  fprintf(fp1, " 9\n");
  fprintf(fp1, " 7\n 7\n");
  fprintf(fp1, "\n");
  
  /*value for each point (e.g. color)*/
  fprintf(fp1, "POINT_DATA %d\n", nphi*ntheta); 
  fprintf(fp1, "SCALARS const float\n");
  fprintf(fp1, "LOOKUP_TABLE default\n");
  for (i=1; i<=nphi*ntheta; i++)  fprintf(fp1, " 0.0\n");
    
  fclose(fp1);
}


// directly output the coefficients of the spherical expansion
void AHmod_output_sphercoeff(const int LMAX1,
      const double xc, const double yc, const double zc,
      double *a0,  double **ac, double **as, const char *outdir,
      const int number, const double time, const double AHmod_time)
{
  int l,m;
  char filename[300], path[1000];
  const int slice = AHmod_ComputeSlice(time,AHmod_time);

  // %.6d   width=6 automatically fill with 0
  sprintf(filename, "%s/AHmod_lm/MOTS%d_%.6d",outdir,number,slice);

  FILE *fp1= fopen(filename, "wb");

  /*check if file exists; if not, create file*/
  if (!fp1) {
    sprintf(path, "%s/AHmod_lm",outdir);
    if (number==0) mkdir(path,0777);
    fp1 = fopen(filename, "wb");
    if (!fp1) errorexits("failed opening %s", filename);
  }

  printf("  Output AHmod (spherical coefficients)\n");

  /*output header*/
  fprintf(fp1, "# Apparent Horizon %d, time=%f \n",number,time);
  fprintf(fp1, "# LMAX=%d, x=%g, y=%g, z=%g, \n",LMAX1,xc, yc, zc);

  /*output a0, ac and as components*/
  for(l=0;l<=LMAX1;l++){
    fprintf(fp1,"%e\n", a0[l]);
    for(m=0;m<=l;m++){
      fprintf(fp1,"%e\n",ac[l][m]);
      fprintf(fp1,"%e\n",as[l][m]);
    }
  }

  fclose(fp1);
}
    

// Break the output of the xy-slice into two parts.
// The first simply appends the 3 x ntheta values
// (each has x and y coordinate of middle points and time)
// to one raw-file.
// The second transforms this data into a valid vtk file.
// That has the advantage that the raw file can simply be appended
// and the vtk file is always generated consistently
void AHmod_output_xyt_vtk(const int ntheta, const double dtheta, 
      const int nphi, const double dphi, 
      const double xc, const double yc, const double zc, 
      double **rr, const char *outdir, const int number, const double time)
{
  char rawname[300], vtkname[300], path[300];
  FILE *fp;
    
  sprintf(rawname, "%s/ah.xyt_vtk/ah.xyt%d.raw",outdir,number);
  sprintf(vtkname, "%s/ah.xyt_vtk/ah.xyt%d.vtk",outdir,number);
  
  // Check whether raw-file exists, 
  // if not then create directory, raw file and vtk file 
  if ( AHmod_FileExists(rawname) == 0 ) {
    // make directory
    sprintf(path, "%s/ah.xyt_vtk",outdir);
    mkdir(path,0777);
    
    // make raw file
    fp = fopen(rawname, "w");
    if (!fp) errorexits("failed opening %s", rawname);
    fclose(fp);
    
    // make vtk file
    fp = fopen(vtkname, "w");
    if (!fp) errorexits("failed opening %s", vtkname);
    fclose(fp);
  }
  
  AHmod_Appendxyt(rawname, ntheta, dtheta, nphi, dphi, xc, yc, zc, rr, time);
  
  AHmod_ConvertRaw2VTK(rawname, vtkname, number, nphi);
  
}

void AHmod_Appendxyt(char *rawname,
      const int ntheta, const double dtheta, const int nphi, const double dphi, 
      const double xc, const double yc, const double zc, 
      double **rr, const double time)
{
  FILE *fp = fopen(rawname, "a");
  if (!fp) {
    errorexits("failed opening %s", rawname);
  }
  
  // if ntheta is odd take the appropriate theta values directly
  // if ntheta is even take the mean of the points a bit below and above z=0 plane
  const int nthetahalf=ntheta/2; // will be rounded towards zero
  int i1,i2,p;
  double x,y,theta1,theta2,phi;
  
  if ((ntheta%2)==0) {
    // calculate the mean of two values
    i1=nthetahalf;
    i2=nthetahalf+1;
    theta1 = AHmod_theta(dtheta, i1);
    theta2 = AHmod_theta(dtheta, i2);
    
    for (p=0;p<nphi;p++) {
      phi = AHmod_phi(dphi, p);
      x = xc + 0.5*( rr[i1][p]*sin(theta1) + rr[i2][p]*sin(theta2) )*cos(phi);
      y = yc + 0.5*( rr[i1][p]*sin(theta1) + rr[i2][p]*sin(theta2) )*sin(phi);
      fprintf(fp, " %f %f %f\n",x,y,time);
    }
  } else{
    i1=nthetahalf+1;
    // theta1 = AHmod_theta(dtheta, i1);
    // in the x-y-plane theta = PI/2, sin(theta)=1
    for(p=0;p<nphi;p++){
      phi = AHmod_phi(dphi,p);
      x = xc + rr[i1][p]*1.0*cos(phi);
      y = yc + rr[i1][p]*1.0*sin(phi);
      fprintf(fp, " %f %f %f\n",x,y,time);
    }
  }
  fclose(fp);
}

// Convert the specific raw data of AHmod xyt to vtk format
// Description of vtk file:
// http://www.cacr.caltech.edu/~slombey/asci/vtk/vtk_formats.simple.html
void AHmod_ConvertRaw2VTK(char *rawname, char *vtkname,
      const int number, const int nphi)
{
  int i;
    
  // open for reading
  FILE *fin = fopen(rawname, "r"); 
  if (!fin) {
    errorexits("failed opening %s", rawname);
  }
  // open for writing, replace existing file
  FILE *fout = fopen(vtkname, "w"); 
  if (!fout) {
    errorexits("failed opening %s", vtkname);
  }
    
  // get # of lines = # of points 
  int ptnum=0;
  char line[1000];
  while ( fgets(line, 1000, fin) != NULL) 
    ptnum++;
  rewind(fin);
    
  // start writing the vtk file
  fprintf(fout, "# vtk DataFile Version 2.0\n");
  fprintf(fout, "Apparent Horizon %d xyt\n",number);
  fprintf(fout, "ASCII\n");
  fprintf(fout, "\n");
  fprintf(fout, "DATASET UNSTRUCTURED_GRID\n");
  fprintf(fout, "POINTS %d float\n", ptnum);

  // copy the points from rawfile to vtkfile, close rawfile
  for (i=0; i<ptnum; i++) {
    fgets(line, 1000, fin);
    fputs(line, fout);
  }
  fprintf(fout,"\n");
  fclose(fin);
  
  // first value: number of points
  // second value: total number: 5 *ptnum, structure (1), corners (4)
  const int numfaces=ptnum-nphi;
  fprintf(fout,"CELLS %d %d\n",numfaces,5*numfaces);
  for(i=0;i<numfaces;i++) {
    if (((i+1)%nphi)==0) 
      fprintf(fout, " 4 %d %d %d %d\n",i, i+1-nphi, i+1, i+nphi);
    else 
      fprintf(fout, " 4 %d %d %d %d\n",i, i+1, i+nphi+1, i+nphi);
  }
  fprintf(fout,"\n\n");
  fprintf(fout, "CELL_TYPES %d\n",numfaces);
  for(i=0;i<numfaces;i++) 
    fprintf(fout, " 9\n");
  fprintf(fout, "\n");    
  fprintf(fout, "POINT_DATA %d\n",ptnum); 
  fprintf(fout, "SCALARS const float\n");
  fprintf(fout, "LOOKUP_TABLE default\n");
  for(i=0;i<ptnum;i++)
    fprintf(fout, " 0.0\n");
  fprintf(fout, "\n");
    
  fclose(fout);
}

