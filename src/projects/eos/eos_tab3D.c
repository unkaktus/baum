/* eos_tab3D.h */
/* mth 06/12 */

#include "bam.h"
#include "eos.h"

// #ifdef HDF5
// #include "hdf5.h"
// #endif

#define PR 0




struct {
  int Nr,NT,NY;
  double *Logr,*T,*Y;
  double dLogr,dT,dY;
  double *eps,*p;
} tab2d;

struct {
  int Nr,NT,NY;
  double *Logr,*LogT,*Y;
  double dLogr,dLogT,dY;
  double *eps,*p;
} tab3d;












/*****************************************************************************/
/* load 3d tables */
void eos_load_tab3d_shen(char *fname, int *n1,int *n2, int *n3,
                          double **x1, double **x2, double **x3, 
                          double **v1, double **v2)
{
  int pr = 0;
  
  tU units;
  set_units(&units);
  
  int P;
  int N1,N2,N3;
  int i,j,k, ijk, nj;
  double tmp, v_p,v_nB,v_LrB,v_Y,v_E, e,rho,nB,mu,E,p,eps;
  char line[1024];
  
  printf("  read in eos 3d table file:\n    %s\n",fname);
  i=j=k=ijk=0;
  N1=N2=N3=1;
  
  if (!(*x1)) (*x1) = (double*) malloc (1*sizeof(double));
  if (!(*x2)) (*x2) = (double*) malloc (1*sizeof(double));
  if (!(*x3)) (*x3) = (double*) malloc (1*sizeof(double));
  if (!(*v1)) (*v1) = (double*) malloc (1*sizeof(double));
  if (!(*v2)) (*v2) = (double*) malloc (1*sizeof(double));
  
  /* do this for each processor */
  for (P=0; P<bampi_size(); P++) {
    
    printf("    read on proc%d/%d\n",P,bampi_size());
    if (bampi_rank()==P) {

      /* open file */
      FILE *ifp = fopen(fname, "r");
      if (!ifp) 
        errorexits(" could not open file %s", fname);
      
      /* read file */
      while (fgets(line, 1024, ifp)) {
        
        if (strncmp(line," cccccccccccc",7)==0) {
          
          j=k=0;
          if (i>=N1) {
            N1=i+1;
            *x1 = (double*) realloc (*x1,N1*sizeof(double));
          }
          i++;
          fgets(line, 1024, ifp);
          fgets(line, 1024, ifp);
          sscanf(line," %le  %le\n",&((*x1)[i-1]),&tmp);
          fgets(line, 1024, ifp);
          nj = 0;
          
        } else if (strlen(line)==3) {
          
          k=0;
          nj++;
          
        } else if (strlen(line)>10) {
          
          if (nj) {j++; nj=0;}
          
          /* memory management */
          if (j>=N2) {
            N2  = j+1;
            *x2 = (double*) realloc (*x2, N2*sizeof(double));
          }
          if (k>=N3) {
            N3  = k+1;
            *x3 = (double*) realloc (*x3, N3*sizeof(double));
          }
          *v1 = (double*) realloc (*v1, N1*N2*N3*sizeof(double));
          *v2 = (double*) realloc (*v2, N1*N2*N3*sizeof(double));
          
          /* read actual data */
          sscanf(line,"  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le",
                  &v_LrB,&v_nB,&v_Y, &tmp,&v_E,&tmp,&tmp,&tmp,&tmp,&tmp,&tmp,&tmp,&tmp,&v_p,&tmp,&tmp);
          
          rho = pow( 10,v_LrB ) / units.Mdens_cgs;
          p   = v_p / units.Energy_MeV*units.Volume_fm3;
          mu  = 931.494;
          eps = v_E/mu;
          
          (*x2)[j]   = v_Y;
          (*x3)[k]   = log10(rho);
          (*v1)[ijk] = p;
          (*v2)[ijk] = eps;
          
          ijk++;
          k++;
        }
        
      }
      
      /* close file */
      fclose (ifp);
    }
  }
  
  //printf("%d %d %d\n", N1,N2,N3);
  *n1 = N1;
  *n2 = N2;
  *n3 = N3;
  
}

// #ifdef HDF5
// void eos_load_tab3d_ott()
// {
//   int pr = 0;
  
//   int i;
//   if (GetsLax("eos_tab_file")==0) errorexit(" need eos_tab_file");
//   char *fname = Gets("eos_tab_file");
//   printf("  read in eos 3d table file:\n    %s\n",fname);
  
//   /* hdf5 handles */
//   hid_t file,dset,dcpl;
//   herr_t status;
//   H5D_layout_t layout;
  
//   /* open file */
//   file = H5Fopen (fname, H5F_ACC_RDONLY, H5P_DEFAULT);
  
//   /* read # points */
//   dset   = H5Dopen  (file, "pointsrho");
//   status = H5Dread  (dset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(NLogrho));
//   status = H5Dclose (dset);
//   dset   = H5Dopen  (file, "pointstemp");
//   status = H5Dread  (dset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(NLogT));
//   status = H5Dclose (dset);
//   dset   = H5Dopen  (file, "pointsye");
//   status = H5Dread  (dset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(NYp));
//   status = H5Dclose (dset);
//   //printf("%d %d %d\n",NLogrho,NLogT,NYp);
  
//   /* read grid */
//   double *Logrho = (double*) malloc (NLogrho*sizeof(double));
//   double *LogT   = (double*) malloc (NLogT*sizeof(double));
//   double *Yp     = (double*) malloc (NYp*sizeof(double));
//   dset   = H5Dopen  (file, "logrho");
//   status = H5Dread  (dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, Logrho );
//   status = H5Dclose (dset);
//   dset   = H5Dopen  (file, "logtemp");
//   status = H5Dread  (dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, LogT );
//   status = H5Dclose (dset);
//   dset   = H5Dopen  (file, "ye");
//   status = H5Dread  (dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, Yp );
//   status = H5Dclose (dset);
//   /*
//   for (i=0; i<NLogrho; i++)
//     printf("%d   %e\n", i, Logrho[i]);
//   for (i=0; i<NLogT; i++)
//     printf("%d   %e\n", i, LogT[i]);
//   for (i=0; i<NYp; i++)
//     printf("%d   %e\n", i, Yp[i]);
//   */
//   Logrho0 = Logrho[0];
//   DLogrho = Logrho[1]-Logrho[0];
//   LogT0   = LogT[0];
//   DLogT   = LogT[1]-LogT[0];
//   Yp0     = Yp[0];
//   DYp     = Yp[1]-Yp[0];
//   free(Logrho);
//   free(LogT);
//   free(Yp);
  
//   /* read data */
//   eostab_Logp = (double*) malloc (NLogrho*NLogT*NYp*sizeof(double));
//   eostab_Loge = (double*) malloc (NLogrho*NLogT*NYp*sizeof(double));
//   eostab_dedT = (double*) malloc (NLogrho*NLogT*NYp*sizeof(double));
//   eostab_dpdrho = (double*) malloc (NLogrho*NLogT*NYp*sizeof(double));
//   eostab_dpdeps = (double*) malloc (NLogrho*NLogT*NYp*sizeof(double));
//   dset   = H5Dopen  (file, "logpress");
//   status = H5Dread  (dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, eostab_Logp );
//   status = H5Dclose (dset);
//   dset   = H5Dopen  (file, "logenergy");
//   status = H5Dread  (dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, eostab_Loge );
//   status = H5Dclose (dset);
//   dset   = H5Dopen  (file, "dedt");
//   status = H5Dread  (dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, eostab_dedT );
//   status = H5Dclose (dset);
//   dset   = H5Dopen  (file, "dpderho");
//   status = H5Dread  (dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, eostab_dpdeps );
//   status = H5Dclose (dset);
//   dset   = H5Dopen  (file, "dpdrhoe");
//   status = H5Dread  (dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, eostab_dpdrho );
//   status = H5Dclose (dset);
  
//   /* compute epsl table */
//   eostab_epsl = (double*) malloc (NLogrho*NLogT*NYp*sizeof(double));
//   for (i=0; i<NLogrho*NLogT*NYp; i++)
//     eostab_epsl[i] = pow(10.,eostab_Loge[i]);
  
//   /* close file */
//   status = H5Fclose (file);
// }
// #else
void eos_load_tab3d_ott()
{
  errorexit("you need to compile with -DHDF5 -lhdf5");
}
// #endif

void eos_load_tab3d()
{
  int i,j,k, ijk2,ijk3, ijk1,ijk4;
  
  
  /* T=0 table */
  eos_load_tab3d_shen(Gets("eos_tab_fileT0"),
                      &(tab2d.NT),&(tab2d.NY),&(tab2d.Nr),
                      &(tab2d.T), &(tab2d.Y), &(tab2d.Logr),
                      &(tab2d.p),&(tab2d.eps));
  
  
  /* full table */
  eos_load_tab3d_shen(Gets("eos_tab_file"), 
                      &(tab3d.NT),  &(tab3d.NY), &(tab3d.Nr),
                      &(tab3d.LogT),&(tab3d.Y), &(tab3d.Logr),
                      &(tab3d.p),&(tab3d.eps));
  
  
  /* expand T=0 table */
  tab2d.NT  = tab2d.NT + 1;
  tab2d.T   = (double*) realloc (tab2d.T,tab2d.NT*sizeof(double));
  tab2d.eps = (double*) realloc (tab2d.eps,tab2d.NT*tab2d.NY*tab2d.Nr*sizeof(double));
  tab2d.p   = (double*) realloc (tab2d.p,  tab2d.NT*tab2d.NY*tab2d.Nr*sizeof(double));
  if (!tab2d.eps || !tab2d.p) errorexit("realloc does not work");
  
  tab2d.T[0] = 0;
  tab2d.T[1] = pow(10,tab3d.LogT[0]);
  for (j=0; j<tab2d.NY; j++) {
    for (k=0; k<tab2d.Nr; k++) {
      
      ijk2 = 1*tab2d.NY*tab2d.Nr + j*tab2d.Nr + k;
      ijk3 = 0*tab3d.NY*tab3d.Nr + j*tab3d.Nr + k;
      tab2d.eps[ijk2] = tab3d.eps[ijk3];
      tab2d.p[ijk2]   = tab3d.p[ijk3];
      
      /*
      ijk1 = 0*tab2d.NY*tab2d.Nr + j*tab2d.Nr + k;
      tab2d.eps[ijk2] = tab2d.eps[ijk1];
      tab2d.p[ijk2]   = tab2d.p[ijk1];
      */
      
      /*
      printf("%e %e | %e %e |    %e  %e  %e  %e\n",tab2d.Logr[j],tab2d.Y[k], tab3d.Logr[j],tab3d.Y[k], 
             tab2d.eps[ijk1], tab2d.eps[ijk2], tab3d.eps[ijk3],tab3d.eps[ijk4]);
      */
    }
  }
  
  
  /* set faster access */
  tab2d.dLogr = tab2d.Logr[1]-tab2d.Logr[0];
  tab2d.dT    = tab2d.T[1]-tab2d.T[0];
  tab2d.dY    = tab2d.Y[1]-tab2d.Y[0];
  
  tab3d.dLogr = tab3d.Logr[1]-tab3d.Logr[0];
  tab3d.dLogT = tab3d.LogT[1]-tab3d.LogT[0];
  tab3d.dY    = tab3d.Y[1]-tab3d.Y[0];
  
  //errorexit("better stop");
}




/*****************************************************************************/
/* 3d interpolation, 2nd order optimized, with derivatives */
int interp_eos_3d(double xp,double yp,double zp,
                  double x0,double y0,double z0,
                  double dx,double dy,double dz,
                  int nx, int ny, int nz,
                  double *tab, double *var,
                  double *dvdx,double *dvdy,double *dvdz)
{
  int pr = 0;
  
  if (0) {
    printf("  %+e | %+e | %+e     (%e)\n", x0,xp,x0+(nx-1)*dx,dx);
    printf("  %+e | %+e | %+e     (%e)\n", y0,yp,y0+(ny-1)*dy,dy);
    printf("  %+e | %+e | %+e     (%e)\n", z0,zp,z0+(nz-1)*dz,dz);
  }
  
  double v[8], a[8];
  int px = floor((xp-x0+dequaleps)/dx);
  int py = floor((yp-y0+dequaleps)/dy);
  int pz = floor((zp-z0+dequaleps)/dz);
  if (px<0 || px>=nx-1) {
    if (pr) errorexit("point outside table 1");
    return 1;
  }
  if (py<0 || py>=ny-1) {
    if (pr) errorexit("point outside table 2");
    return 2;
  }
  if (pz<0 || pz>=nz-1) {
    if (pr) errorexit("point outside table 3");
    return 3;
  }
  
  double x = (xp - (x0 + px*dx))/dx;
  double y = (yp - (y0 + py*dy))/dy;
  double z = (zp - (z0 + pz*dz))/dz;
  int Dx = ny*nz;
  int Dy = nz;
  int Dz = 1;
  int ijk = px*Dx + py*Dy + pz*Dz;
  
  if (0) {
    int i;
    for (i=0; i<nx-1; i++) {
      ijk = i*Dx + py*Dy + pz*Dz;
      printf(" %d   %e  %e  %e\n", i, x0+i*dx, tab[ijk], tab[ijk+Dx] - tab[ijk]);
    }
    exit(0);
  }
  
  v[0] = tab[ijk + 0*Dx + 0*Dy + 0*Dz];
  v[1] = tab[ijk + 0*Dx + 0*Dy + 1*Dz];
  v[2] = tab[ijk + 0*Dx + 1*Dy + 0*Dz];
  v[3] = tab[ijk + 0*Dx + 1*Dy + 1*Dz];
  v[4] = tab[ijk + 1*Dx + 0*Dy + 0*Dz];
  v[5] = tab[ijk + 1*Dx + 0*Dy + 1*Dz];
  v[6] = tab[ijk + 1*Dx + 1*Dy + 0*Dz];
  v[7] = tab[ijk + 1*Dx + 1*Dy + 1*Dz];
  
  a[7] = v[0];
  a[6] = v[1]-v[0];
  a[5] = v[2]-v[0];
  a[4] = v[4]-v[0];
  a[3] = v[0]-v[1]-v[2]+v[3];
  a[2] = v[5]-v[4]-v[1]+v[0];
  a[1] = v[6]-v[4]-v[2]+v[0];
  a[0] = v[7]-v[6]-v[5]+v[4]-v[3]+v[2]+v[1]-v[0];
  
  /* value */
  *var  = a[0]*x*y*z + a[1]*x*y + a[2]*x*z + a[3]*y*z + a[4]*x + a[5]*y + a[6]*z + a[7];
  
  /* derivatives */
  *dvdx = (a[0]  *y*z + a[1]  *y + a[2]  *z            + a[4]                           )/dx;
  *dvdy = (a[0]*x  *z + a[1]*x              + a[3]*  z          + a[5]                  )/dy;
  *dvdz = (a[0]*x*y              + a[2]*x   + a[3]*y                     + a[6]         )/dz;
  
  return 0;
}




/*****************************************************************************/
/* interpolate EoS Tab */
void eos_tab3d_findT(double *rho, double *epsl, double *Y, double *T)
{
  if (PR) printf("%e %e %e -> %e\n", *rho,*epsl,*Y,*T);
  
  int pr = 0;
  
  int step;
  int    steps    = 100;
  double acc      = 1e-10;
  double minLogT  = tab3d.LogT[0];
  double maxLogT  = tab3d.LogT[tab3d.NT-2];
  
  double Logr = log10( *rho );
  double LogT = tab3d.LogT[4];
  double eps0 = *epsl;
  double eps,depsdLogT,depsdT, nLogT, d2,d3;
  
  if (pr) printf("%e %e %e %e %e\n", minLogT,maxLogT,Logr,LogT, eps0);
  
  /* do interpolation once with T=Tmin_3d  */
  interp_eos_3d(LogT,*Y,Logr,
                tab3d.LogT[0],tab3d.Y[0],tab3d.Logr[0],
                tab3d.dLogT,tab3d.dY,tab3d.dLogr,
                tab3d.NT,tab3d.NY,tab3d.Nr,
                tab3d.eps,&eps, 
                &depsdLogT,&d2,&d3);
  
  if (pr) printf("  %d   Logrho=%+e  LogT=%+e  Yp=%+e   ->   eps=%+e eps0= %+e   DedLT=%+e  %e \n", 
      -1, Logr,LogT,*Y, eps,eps0, eps-eps0,depsdLogT);
  
  /* change T in a way that the guess for eps is equal to the needed epsl */
  step = 0;
  while (fabs(eps-eps0)>acc && step<steps) {
    
    /* find new temperature */
    nLogT       = LogT - (eps-eps0)/depsdLogT;
    
    /* use T=0 table??? */
    if (nLogT<minLogT) {
      LogT  = nLogT;
      *T    = (nLogT<1e-20)?1e-19:pow(10,nLogT);
      
      /* FIXME: he we should test if eps is too small for the table
         which is equivalent to an unphysical point */
      interp_eos_3d(*T,*Y,Logr,
                    tab2d.T[0],tab2d.Y[0],tab2d.Logr[0],
                    tab2d.dT,tab2d.dY,tab2d.dLogr,
                    tab2d.NT,tab2d.NY,tab2d.Nr,
                    tab2d.eps,&eps, 
                    &depsdT,&d2,&d3);
      
      depsdLogT = (*T)*depsdT;
      
    } else {
      
      nLogT       = (nLogT<minLogT)?minLogT:nLogT;
      nLogT       = (nLogT>maxLogT)?maxLogT:nLogT;
      LogT        = nLogT;
      *T          = pow(10,nLogT);
      
      /* find eps for this value */
      interp_eos_3d(LogT,*Y,Logr,
                    tab3d.LogT[0],tab3d.Y[0],tab3d.Logr[0],
                    tab3d.dLogT,tab3d.dY,tab3d.dLogr,
                    tab3d.NT,tab3d.NY,tab3d.Nr,
                    tab3d.eps,&eps, 
                    &depsdLogT,&d2,&d3);
    }
    
    if (pr) printf("   %d   Logrho=%+e  LogT=%+e  Yp=%+e   ->   eps=%+e Delta=%+e   DedLT=%+e  \n", 
        step, Logr,LogT,*Y, eps, eps-eps0,depsdLogT);
    step++;
  }
  
  if (pr) {
    printf("Endwerte:\n");
    printf("  rho=%e  T=%e   epsl=%e<>%e   Yp=%e  \n", *rho,*T,*epsl, eps ,*Y);
  }
  
  if (PR) printf("%e %e %e -> %e\n", *rho,*epsl,*Y,*T);
}

void eos_tab3d_findp(double *rho, double *T, double *Y, 
                     double *pres, double *epsl, 
                     double *dpdr_T, double *dedr_T, double *dpdT_r, double *dedT_r)
{

  double p,eps, dpdLogr,dedLogr, d1,d2,d3;
  double Logr = log10( *rho );
  double minT = pow(10,tab3d.LogT[0]);
  if (*T<minT) {
    
    interp_eos_3d(*T,*Y,Logr,
                  tab2d.T[0],tab2d.Y[0],tab2d.Logr[0],
                  tab2d.dT,tab2d.dY,tab2d.dLogr,
                  tab2d.NT,tab2d.NY,tab2d.Nr,
                  tab2d.eps,&eps, 
                  &d1,dedT_r,&dedLogr);
    
    interp_eos_3d(*T,*Y,Logr,
                  tab2d.T[0],tab2d.Y[0],tab2d.Logr[0],
                  tab2d.dT,tab2d.dY,tab2d.dLogr,
                  tab2d.NT,tab2d.NY,tab2d.Nr,
                  tab2d.p,&p, 
                  &d1,dpdT_r,&dpdLogr);
  } else {
    
    double LogT = log10(*T);
    
    interp_eos_3d(LogT,*Y,Logr,
                  tab3d.LogT[0],tab3d.Y[0],tab3d.Logr[0],
                  tab3d.dLogT,tab3d.dY,tab3d.dLogr,
                  tab3d.NT,tab3d.NY,tab3d.Nr,
                  tab3d.eps,&eps,
                  &d1,dedT_r,&dedLogr);
    
    interp_eos_3d(LogT,*Y,Logr,
                  tab3d.LogT[0],tab3d.Y[0],tab3d.Logr[0],
                  tab3d.dLogT,tab3d.dY,tab3d.dLogr,
                  tab3d.NT,tab3d.NY,tab3d.Nr,
                  tab3d.p,&p,
                  &d1,dpdT_r,&dpdLogr);
    
  }
  
  *dpdr_T = 1./(*rho*log(10.)) * dpdLogr;      /* dp/drho | const T,Y */
  *dedr_T = 1./(*rho*log(10.)) * dedLogr;      /* depsl/drho | const T,Y */
  
  *pres   = p;
  *epsl   = eps;
  
  //if (PR) printf("%e %e %e %e -> %e %e %e %e\n", *rho,*epsl,*Y,*T, *pres,*dpdrho,*dpdepsl,*cs2);
}

void eos_cs2_rTp(double rho, double epsl, double p, 
                 double dpdr_T, double dedr_T, double dpdT_r, double dedT_r,
                 double *cs2, double *dpdr_e, double *dpde_r)
{
  double h;
  
  *dpdr_e = dpdr_T - dpdT_r * dedr_T/dedT_r;
  *dpde_r = dpdT_r/dedT_r;
  
  h       = rho + rho*epsl + p;
  *cs2    = 1./h * ((*dpdr_e) + p/rho * (*dpde_r));
}

int eos_tab3d_shen(double *rho, double *epsl, double *Y, 
                   double *p, double *T, double *cs2,
                   double *dpdrho, double *dpdepsl)
{
  double eps, dpde_T,dpdr_T,dedT_r,dpdT_r;
  
  /* find T for given epsl */
  eos_tab3d_findT(rho,epsl,Y, T);
  if (!finite(*T))
    return 1;
  
  /* find p for given T,Y,rho and get the derivatives */
  eos_tab3d_findp(rho,T,Y, p,&eps, &dpde_T,&dpdr_T,&dpdT_r,&dpdT_r);
  if (!finite(*p) || !finite(dpde_T) || !finite(dpdr_T) || !finite(dpdT_r) || !finite(dpdT_r))
    return 1;
  
  /* note, the tables are in T dependence, cs2 is compute in epsl dependence */
  eos_cs2_rTp(*rho,*epsl,*p, dpde_T,dpdr_T,dpdT_r,dpdT_r, cs2,dpdepsl,dpdrho);
  if (!finite(*cs2))
    return 1;
  
  return 0;
}

int eos_tab3d_shen_T(double *rho, double *epsl, double *Y, 
                     double *p, double *T, double *cs2,
                     double *dpdrho, double *dpdepsl)
{
  double dpde_T,dpdr_T,dedT_r,dpdT_r;
  
  /* find p for given T,Y,rho and get the derivatives */
  eos_tab3d_findp(rho,T,Y, p,epsl, &dpde_T,&dpdr_T,&dpdT_r,&dpdT_r);
  if (!finite(*p) || !finite(dpde_T) || !finite(dpdr_T) || !finite(dpdT_r) || !finite(dpdT_r))
    return 1;
  
  /* note, the tables are in T dependence, cs2 is compute in epsl dependence */
  eos_cs2_rTp(*rho,*epsl,*p, dpde_T,dpdr_T,dpdT_r,dpdT_r, cs2,dpdepsl,dpdrho);
  if (!finite(*cs2))
    return 1;
  
  return 0;
}

int eos_tab3d_shen_Y(double *p, double *cs2, double *dpdrho, double *dpdepsl, double *rho, double *epsl)
{
  if (!finite(*rho) || !finite(*epsl))
    return 1;
  
  double eps, dpde_T,dpdr_T,dedT_r,dpdT_r;
  double T;
  double Y = EOS.Y;
  eos_tab3d_findT(rho,epsl,&Y, &T);
  if (!finite(T))
    return 1;
  
  /* find p for given T,Y,rho and get the derivatives */
  eos_tab3d_findp(rho,&T,&Y, p,&eps, &dpde_T,&dpdr_T,&dpdT_r,&dpdT_r);
  if (!finite(*p) || !finite(dpde_T) || !finite(dpdr_T) || !finite(dpdT_r) || !finite(dpdT_r))
    return 1;
  
  /* note, the tables are in T dependence, cs2 is compute in epsl dependence */
  eos_cs2_rTp(*rho,*epsl,*p, dpde_T,dpdr_T,dpdT_r,dpdT_r, cs2,dpdepsl,dpdrho);
  if (!finite(*cs2))
    return 1;
  
  if (PR) printf("%e %e %e -> %e -> %e %e %e %e\n\n", *rho,*epsl,Y,T, *p,*dpdrho,*dpdepsl,*cs2);
  return 0;
}

int eos_tab3d_shen_T0_Y(double *p, double *cs2, double *dpdrho, double *dpdepsl, double *rho, double *epsl)
{
  double dpdr_T,dedr_T,dpdT_r,dedT_r;
  double T = 0.0;
  double Y = EOS.Y;
  
  if (!finite(*rho) || !finite(*epsl))
    return 1;
  
  /* find p for given T,Y,rho and get the derivatives */
  eos_tab3d_findp(rho,&T,&Y, p,epsl, &dpdr_T,&dedr_T,&dpdT_r,&dedT_r);
  if (!finite(*p) || !finite(dpdr_T) || !finite(dpdr_T) || !finite(dpdT_r) || !finite(dedT_r))
    return 1;
  
  /* note, the tables are in T dependence, cs2 is compute in epsl dependence */
  eos_cs2_rTp(*rho,*epsl,*p, dpdr_T,dedr_T,dpdT_r,dedT_r, cs2,dpdrho,dpdepsl);
  if (!finite(*cs2))
    return 1;
  
  return 0;
}












/*****************************************************************************/
/* EoS test function */
int eos_tab3d_test(tL *level)
{
  double p    = 0.; 
  double rho  = (pow(10,tab2d.Logr[0])+pow(10,tab2d.Logr[1]) )/2.; 
  double epsl = 4.13e-2;
  double Yp   = 0.35;
  
  double T,c,d1,d2,lr;
  
  
  rho  = pow(10,tab2d.Logr[0]);
  eos_tab3d_shen_T0_Y(&p, &c,&d1,&d2, &rho,&epsl);
  
  rho  = pow(10,(0.5*tab2d.Logr[0]+0.5*tab2d.Logr[1])); 
  eos_tab3d_shen_T0_Y(&p, &c,&d1,&d2, &rho,&epsl);
  
  rho  = pow(10,tab2d.Logr[1]); 
  eos_tab3d_shen_T0_Y(&p, &c,&d1,&d2, &rho,&epsl);
  
  int i;
  for (i=0; i<=100; i++) {
    lr = tab2d.Logr[0]*i/100. + (100.-i)/100.*tab2d.Logr[3];
    rho = pow(10,lr);
    eos_tab3d_shen_T0_Y(&p, &c,&d1,&d2, &rho,&epsl);
    printf("%e %e %e %e\n", rho, lr, p, (p-100.*rho*rho)/p);
  }
  
  
  errorexit("");
}







