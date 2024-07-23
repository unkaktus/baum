/* AHF_output.c */
/* mth 12/07 */

#include "bam.h"
#include "AHF.h"

#include <sys/types.h>
#include <sys/stat.h>


void AHF_output(int ntheta, double dtheta, int nphi, double dphi, double xc, double yc, double zc, double ***rr, char **outdir, int number, double time, double ahf_time)
{
  int i,j,k,l;
  double theta, phi, x,y,z;
  char name[208], str[10], command[1000];
  FILE *fp1;
  
  int slice = (int)(time/ahf_time);
  
  /*create file-name*/
  strcpy(name,(*outdir));
  strcat(name,"/ah.xyz_vtk/ah.xyz");
  sprintf(str,"%d",number);
  strcat(name,str);  
  if (slice<10) strcat(name,"_000");
  else if (slice<100) strcat(name,"_00");
  else if (slice<1000) strcat(name,"_0");
  else if (slice>=1000) strcat(name,"_");
  sprintf(str,"%d",slice);
  strcat(name,str);
  strcat(name,".vtk");
  
  fp1 = fopen(name, "wb");
  
  /*check if file exists; if not, create file*/
  if (!fp1) {
    strcpy(command,(*outdir));
    strcat(command,"/ah.xyz_vtk");
    if (number==0) mkdir(command,0777);
    fp1 = fopen(name, "wb");
    if (!fp1) errorexits("failed opening %s", name);
  }
  
  
  printf("\n\noutput AHF (vtk-format)\n");
  
  
  /*output header*/
  fprintf(fp1, "# vtk DataFile Version 2.0\n");
  fprintf(fp1, "Apparent Horizon %d, time=%f\n",number,time);
  fprintf(fp1, "ASCII\n");
  fprintf(fp1, "\n");
  fprintf(fp1, "DATASET UNSTRUCTURED_GRID\n");
  /*output locations of all points*/
  fprintf(fp1, "POINTS %d float\n", nphi*ntheta);
  for(i=0;i<ntheta;i++){
    theta = dtheta*(1.0/2.0 + i);
    for(j=0;j<nphi;j++){
      phi = dphi*(1.0/2.0 + j);
      x = xc + (*rr)[i][j]*sin(theta)*cos(phi);
      y = yc + (*rr)[i][j]*sin(theta)*sin(phi);
      z = zc + (*rr)[i][j]*cos(theta);
      fprintf(fp1, " %f  %f  %f\n",x,y,z);
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





void AHF_output_xyt(int ntheta, double dtheta, int nphi, double dphi, double xc, double yc, double zc, double ***rr, char **outdir, int number, double time, double ahf_time)
{
  int i,j, i1,i2,n,c1,c2;
  double theta,theta1,theta2, phi, x,y;
  char name[208], str[1000], command[1000];
  FILE *fp1,*fp2;
  
  
  strcpy(name,(*outdir));
  strcat(name,"/ah.xyz_vtk/ah.xyt");
  sprintf(str,"%d",number);
  strcat(name,str);
  strcat(name,".vtk");
  
      
  fp1 = fopen(name, "r");

  if (!fp1) {
    fp1 = fopen(name, "wb");
    if (!fp1) errorexits("failed opening %s", name);    
    fprintf(fp1, "# vtk DataFile Version 2.0\n");
    fprintf(fp1, "Apparent Horizon %d xyt\n",number);
    fprintf(fp1, "ASCII\n");
    fprintf(fp1, "\n");
    fprintf(fp1, "DATASET UNSTRUCTURED_GRID\n");
    fprintf(fp1, "POINTS %d float\n", nphi);
    if ((ntheta%2)==0) {
      i1=(int)((double)(ntheta)/2.)-1;
      theta1 = dtheta*(1.0/2.0 + i1);
      i2=(int)((double)(ntheta)/2.);
      theta2 = dtheta*(1.0/2.0 + i2);
      for(j=0;j<nphi;j++){
        phi = dphi*(1.0/2.0 + j);
        x = xc + 0.5*((*rr)[i1][j]*sin(theta1)*cos(phi) + (*rr)[i2][j]*sin(theta2)*cos(phi));
        y = yc + 0.5*((*rr)[i1][j]*sin(theta1)*sin(phi) + (*rr)[i2][j]*sin(theta2)*sin(phi));
        //z = zc + ((*rr)[i1][j]*cos(theta1) + (*rr)[i2][j]*cos(theta2));
        fprintf(fp1, " %f  %f  %f\n",x,y,time);
      }      
    }
    else{
      i=(int)((double)(ntheta-1)/2.);
      theta = dtheta*(1.0/2.0 + i);
      for(j=0;j<nphi;j++){
	phi = dphi*(1.0/2.0 + j);
        x = xc + (*rr)[i][j]*sin(theta)*cos(phi);
        y = yc + (*rr)[i][j]*sin(theta)*sin(phi);
        //z = zc + (*rr)[i][j]*cos(theta);
        fprintf(fp1, " %f  %f  %f\n",x,y,time);
      }
    }
    fprintf(fp1, "\n");
    fprintf(fp1, "CELLS 0 0\n");
    fprintf(fp1, "\n");
    fprintf(fp1, "CELL_TYPES 0\n");
    fprintf(fp1, "\n");
    fprintf(fp1, "POINT_DATA 0\n"); 
    fprintf(fp1, "SCALARS const float\n");
    fprintf(fp1, "LOOKUP_TABLE default\n");
    fprintf(fp1, "\n");
    fclose(fp1);
    
    
  }
  else{
    strcpy(name,(*outdir));
    strcat(name,"/ah.xyz_vtk/tmp.vtk");
    
    fp2 = fopen(name, "wb"); 
    
    for (i=1;i<=5;i++) {
      fgets(str, 1000, fp1);
      fprintf(fp2,"%s", str);
    }
    fscanf(fp1, "POINTS %d float\n", &n);
    fprintf(fp2,"POINTS %d float\n", n+ntheta);
    for (i=7;i<7+n;i++) {
      fgets(str, 1000, fp1);
      fprintf(fp2,"%s", str);
    }
    if ((ntheta%2)==0) {
      i1=(int)((double)(ntheta)/2.)-1;
      theta1 = dtheta*(1.0/2.0 + i1);
      i2=(int)((double)(ntheta)/2.);
      theta2 = dtheta*(1.0/2.0 + i2);
      for(j=0;j<nphi;j++){
        phi = dphi*(1.0/2.0 + j);
        x = xc + 0.5*((*rr)[i1][j]*sin(theta1)*cos(phi) + (*rr)[i2][j]*sin(theta2)*cos(phi));
        y = yc + 0.5*((*rr)[i1][j]*sin(theta1)*sin(phi) + (*rr)[i2][j]*sin(theta2)*sin(phi));
        //z = zc + (*rr)[i1][j]*cos(theta1) + (*rr)[i2][j]*cos(theta2);
        fprintf(fp2, " %f  %f  %f\n",x,y,time);
      }      
    }
    else{
      i=(int)((double)(ntheta-1)/2.);
      theta = dtheta*(1.0/2.0 + i);
      for(j=0;j<nphi;j++){
	phi = dphi*(1.0/2.0 + j);
        x = xc + (*rr)[i][j]*sin(theta)*cos(phi);
        y = yc + (*rr)[i][j]*sin(theta)*sin(phi);
        //z = zc + (*rr)[i][j]*cos(theta);
        fprintf(fp2, " %f  %f  %f\n",x,y,time);
      }
    } 
    
    fscanf(fp1, "\n");
    fprintf(fp2,"\n");
    fscanf(fp1, "CELLS %d %d\n",&c1,&c2);
    fprintf(fp2,"CELLS %d %d\n",c1+nphi,c2+5*nphi);//????
    for(j=0;j<nphi+c1;j++) 
      if (((j+1)%nphi)==0) fprintf(fp2, " 4 %d %d %d %d\n",j, j+1-nphi, j+1, j+nphi);
      else fprintf(fp2, " 4 %d %d %d %d\n",j, j+1, j+nphi+1, j+nphi);
    
    fprintf(fp2,"\n");
    fprintf(fp2, "CELL_TYPES %d\n",c1+nphi);
    for(j=0;j<c1+nphi;j++) fprintf(fp2, " 9\n");
    fprintf(fp2, "\n");
    fprintf(fp2, "POINT_DATA %d\n",c1+nphi); 
    fprintf(fp2, "SCALARS const float\n");
    fprintf(fp2, "LOOKUP_TABLE default\n");
    for(j=0;j<c1+nphi;j++) fprintf(fp2, " 0.0\n");
    fprintf(fp2, "\n");
    
    fclose(fp2);
    fclose(fp1);
    
    sprintf(command,"mv %s/ah.xyz_vtk/tmp.vtk %s/ah.xyz_vtk/ah.xyt%d.vtk",(*outdir),(*outdir),number);
    system(command);
  }
 
}

