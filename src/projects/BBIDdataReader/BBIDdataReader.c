/* BBIDdataReader.c */
/* mth 1/11 */

#include "bam.h"
#include "BBIDdataReader.h"

#define PR 0

int useMATTER = 0;

int BBID_read_ps_parameters(tL *level)
{
   read_ps_parameters(level);

   return 0;
}


/* set TOV star using spherstar/dookie data */
int     chebycoe_Extremes( double *u, int n, int inv ) 
{
  
  int k, j, isignum, N=n-1;
  double fac, sum, PioN, *c;
  
  c=(double *) malloc ( n * sizeof(double) );  
  
  PioN=PI/N;
  if (inv == 0) {
    fac=2.0/N;
    isignum=1;
    for (j=0; j<n; j++) {
      sum=0.5*(u[0]+u[N]*isignum);
      for (k=1; k<N; k++)
        sum += u[k]*cos(PioN*j*k);
      c[j]=fac*sum*isignum; 
      isignum = -isignum;
    }
    c[N] = 0.5*c[N];
  }
  else {
    for (j=0; j<n; j++) {
      sum=-0.5*u[0];
      isignum=1;
      for (k=0; k<n; k++){
        sum += u[k]*cos(PioN*j*k)*isignum;
        isignum = -isignum;
      }
      c[j]=sum;
    }
  }
  for (j=0;j<n;j++)
    u[j]=c[j];

  free(c);

  return 1;
}

double  chebyeva( double a, double b, double *c, int m, double x ) 
{
  
  double d=0.0,dd=0.0,sv,y,y2;
  int j;
  
  if ((x-a)*(x-b) > 0.0){ 
    printf( "ERROR in [interp.c] Routine 'chebyeva' :  x not in range [%g,%g].\n",a,b );
    return 0;
  }
  y2=2.0*(y=(2.0*x-a-b)/(b-a));
  for (j=m-1;j>=1;j--) {
    sv=d;
    d=y2*d-dd+c[j];
    dd=sv;
  }
  return( y*d-dd+0.5*c[0] );
}

double  invg4(double gd[4][4], double gu[4][4])
{
    // this was done by mathematica by drag and drop
  double a = gd[0][0];
  double b = gd[0][1];
  double c = gd[0][2];
  double d = gd[0][3];
  double e = gd[1][1];
  double f = gd[1][2];
  double g = gd[1][3];
  double h = gd[2][2];
  double i = gd[2][3];
  double j = gd[3][3];
    
  double detg = -(pow(c,2)*pow(g,2)) + a*pow(g,2)*h + 
        pow(d,2)*(-pow(f,2) + e*h) + 2*b*c*g*i - 2*a*f*g*i - 
        pow(b,2)*pow(i,2) + a*e*pow(i,2) + 
        2*d*(c*f*g - b*g*h - c*e*i + b*f*i) + (pow(c,2)*e - 
            2*b*c*f + a*pow(f,2) + pow(b,2)*h - a*e*h)*j;
  double oodetg = 1./detg;
    
  gu[0][0] = oodetg*(pow(g,2)*h - 2*f*g*i + pow(f,2)*j + e*(pow(i,2) - h*j));
  gu[0][1] = oodetg*(-(d*g*h) + d*f*i + c*g*i - b*pow(i,2) - c*f*j + b*h*j);
  gu[0][2] = oodetg*(d*f*g - c*pow(g,2) - d*e*i + b*g*i + c*e*j - b*f*j);
  gu[0][3] = oodetg*(-(d*pow(f,2)) + c*f*g + d*e*h - b*g*h - c*e*i + b*f*i);
  gu[1][1] = oodetg*(pow(d,2)*h - 2*c*d*i + pow(c,2)*j + a*(pow(i,2) - h*j));
  gu[1][2] = oodetg*(-(pow(d,2)*f) + c*d*g + b*d*i - a*g*i - b*c*j + a*f*j);
  gu[1][3] = oodetg*(c*d*f - pow(c,2)*g - b*d*h + a*g*h + b*c*i - a*f*i);
  gu[2][2] = oodetg*(pow(d,2)*e - 2*b*d*g + pow(b,2)*j + a*(pow(g,2) - e*j));
  gu[2][3] = oodetg*(-(c*d*e) + b*d*f + b*c*g - a*f*g - pow(b,2)*i + a*e*i);
  gu[3][3] = oodetg*(pow(c,2)*e - 2*b*c*f + pow(b,2)*h + a*(pow(f,2) - e*h));
    
  gu[1][0] = gu[0][1];
  gu[2][0] = gu[0][2];
  gu[3][0] = gu[0][3];
  gu[2][1] = gu[1][2];
  gu[3][1] = gu[1][3];
  gu[3][2] = gu[2][3];
    
  if (0) {
    int o,p,q;
    double test;
    for (o=0; o<=3; o++) {
      for (p=0; p<=3; p++) {
        test = 0.;
        for (q=0; q<=3; q++)
          test += gd[o][q]*gu[q][p];
        printf("%e ",test);
      }
      printf("\n");
    }
    printf("\n");
  }
    
  return detg;
}

double  invg3(double gd[3][3], double gu[3][3])
{
    
  gu[0][0] = gd[1][1]*gd[2][2] - gd[1][2]*gd[1][2];
  gu[0][1] = gd[0][2]*gd[1][2] - gd[0][1]*gd[2][2];
  gu[0][2] = gd[0][1]*gd[1][2] - gd[0][2]*gd[1][1];
  gu[1][1] = gd[0][0]*gd[2][2] - gd[0][2]*gd[0][2];
  gu[1][2] = gd[0][1]*gd[0][2] - gd[0][0]*gd[1][2];
  gu[2][2] = gd[0][0]*gd[1][1] - gd[0][1]*gd[0][1];
    
  double detg = gd[0][0]*gu[0][0] + gd[0][1]*gu[0][1] + gd[0][2]*gu[0][2];
  double oodetg = 1./detg;
    
  gu[0][0] *= oodetg;
  gu[0][1] *= oodetg;
  gu[0][2] *= oodetg;
  gu[1][1] *= oodetg;
  gu[1][2] *= oodetg;
  gu[2][2] *= oodetg;
    
  gu[1][0] = gu[0][1];
  gu[2][0] = gu[0][2];
  gu[2][1] = gu[1][2];
    
  return detg;
}

void    Gamma44(double g[4][4], double dg[4][4][4], double Gamma[4][4][4])
{
  int o,p,q,r;
  double gi[4][4];
    
  invg4(g,gi);
    
  for (o=0; o<=3; o++)
    for (p=0; p<=3; p++)
      for (q=0; q<=3; q++) {
    Gamma[o][p][q] = 0.;
    for (r=0; r<=3; r++) {
      Gamma[o][p][q] += 0.5*gi[o][r]* (dg[q][p][r] + dg[p][q][r] - dg[r][p][q]);
    }
      }
    
      if (0) {
        for (o=0; o<=3; o++) {
          printf("Gamma^%d_munu:\n",o);
          for (p=0; p<=3; p++) {
            for (q=0; q<=3; q++)
              printf(" %e",Gamma[o][p][q]);
            printf("\n");
          }
          printf("\n");
        }
        printf("\n");
      }
}

void    Gamma34(double g[4][4], double dg[4][4][4], double Gamma[3][3][3])
{
  int o,p,q,r;
  double gi[4][4];
    
  invg4(g,gi);
    
  for (o=1; o<=3; o++)
    for (p=1; p<=3; p++)
      for (q=1; q<=3; q++) {
    Gamma[o][p][q] = 0.;
    for (r=1; r<=3; r++) {
      Gamma[o][p][q] += 0.5*(gi[o][r] - gi[0][o]*gi[0][r]/gi[0][0])*(dg[q][p][r] + dg[p][q][r] - dg[r][p][q]);
    }
      }
}

void    Gamma33(double g[3][3], double dg[3][3][3], double Gamma[3][3][3])
{
  int o,p,q,r;
  double gi[3][3];
    
  invg3(g,gi);
    
  for (o=0; o<3; o++)
    for (p=0; p<3; p++)
      for (q=0; q<3; q++) {
    Gamma[o][p][q] = 0.;
    for (r=0; r<3; r++) {
      Gamma[o][p][q] += 0.5*gi[o][r]*(dg[q][p][r] + dg[p][q][r] - dg[r][p][q]);
    }
      }
}

void    setLAMBDA(double LAMBDA[4][4], double LAMBDAi[4][4], double xix,double xiy,double xiz)
{

  // 1e-10 saves the day if no boost is set and you have to divide by 0 inside LAMBDA
  double xi = sqrt(xix*xix + xiy*xiy + xiz*xiz)+1e-13;
  double gb = 1./sqrt(1.-xi*xi);
    
  LAMBDA[0][0] = gb;
  LAMBDA[0][1] = LAMBDA[1][0] = gb*xix;
  LAMBDA[0][2] = LAMBDA[2][0] = gb*xiy;
  LAMBDA[0][3] = LAMBDA[3][0] = gb*xiz;
  LAMBDA[1][1] = (1.+(gb-1.)*(xix*xix)/(xi*xi));
  LAMBDA[2][2] = (1.+(gb-1.)*(xiy*xiy)/(xi*xi));
  LAMBDA[3][3] = (1.+(gb-1.)*(xiz*xiz)/(xi*xi));
  LAMBDA[1][2] = LAMBDA[2][1] = (gb-1.)*xix*xiy/(xi*xi);
  LAMBDA[1][3] = LAMBDA[3][1] = (gb-1.)*xix*xiz/(xi*xi);
  LAMBDA[2][3] = LAMBDA[3][2] = (gb-1.)*xiy*xiz/(xi*xi);
    
  invg4(LAMBDA,LAMBDAi);
}

int     interp_locate(double *x, int Nx, double xval)
{
  int ju,jm,jl;
  int ascnd;
  
  jl=-1;
  ju=Nx;
  
  if (xval <= x[0]) { 
    if (xval < x[0]) if (PR) printf("  pt to locate is outside (xval<xx).\n"); 
    return 0; 
  } else if (xval >= x[Nx-1]) { 
    if (xval > x[Nx-1])if (PR)  printf("  pt to locate is outside (xval>xx).\n"); 
    return Nx-1; 
  }
  
  ascnd = (x[Nx-1] >= x[0]);
  
  while (ju-jl > 1) {
    
    jm = (ju+jl) >> 1;
    
    if (xval >= x[jm] == ascnd)
      jl=jm;
    else
      ju=jm;
    
  }
  
  return jl;
}

void    interp_lag4(double *f, double *x, int Nx, double xv, 
                    double *fv_p, double *dfv_p, double *ddfv_p )
{
  /* Given the values in xv, it returns the interpolated values fv and 
  its 1st and 2nd derivatives dfv, ddfv of the fuction f(x) 
  Lagrangian 4 pts interpolation is used */
  
  if (Nx < 4) errorexit(" too few points for interpolation");

  int i = interp_locate(x,Nx,xv);
  
  if( i < 1 ){ 
    if (PR) printf(" too few points on the left => interpolation maybe be inaccurate! (v=%e)\n",xv);
    i = 1;
  }
  if( i > (Nx-3) ){ 
    if (1+PR) printf(" too few points on the right => interpolation maybe be inaccurate! (v=%e   -> %e %e)\n",xv, x[Nx-2],x[Nx-1]);
    i = Nx-3;
  }
  
  double ximo =  x[i-1];
  double xi   =  x[i];
  double xipo =  x[i+1]; 
  double xipt =  x[i+2]; 
  
  double C1   = (f[i] - f[i-1])/(xi - ximo);
  double C2   = (-f[i] + f[i+1])/(-xi + xipo);
  double C3   = (-f[i+1] + f[i+2])/(-xipo + xipt);
  double CC1  = (-C1 + C2)/(-ximo + xipo);
  double CC2  = (-C2 + C3)/(-xi + xipt);
  double CCC1 = (-CC1 + CC2)/(-ximo + xipt);
  
  *fv_p   = f[i-1] + (-ximo + xv)*(C1 + (-xi + xv)*(CC1 + CCC1*(-xipo + xv)));
  *dfv_p  = C1 - (CC1 - CCC1*(xi + xipo - 2.*xv))*(ximo - xv)
      + (-xi + xv)*(CC1 + CCC1*(-xipo + xv));
  *ddfv_p = 2.*(CC1 - CCC1*(xi + ximo + xipo - 3.*xv));
}














/* read a TOV star at a point, think about the boost and set new coordinates */
void    read_TOV(tL* level, char *file, 
                 double px,double py,double pz, 
                 double mx,double my,double mz) 
{ 
 
/*****************************************************************************/
/* THIS IS FROM DOOKIE */
  FILE *fid;
  char line[1024];
  char svalue[256];
  
  unsigned int ix, ixR;
  double M, R;
  double *rs, *lap, *Dlap, *psi4_, *Dpsi4;
  double *rho_, *pres, *epsl_;
  unsigned int dim;

  fid = fopen( file, "r" );  
  if( fid == NULL ) {
    errorexit("ERROR in [initialdata.c] Routine 'id_STAR' : 'fopen' failed.\n" );
  } 
  int len = strlen(Gets("initial_file"));
  if (strncmp( &(Gets("initial_file")[len-4]),".ini",4)!=0) 
    errorexit("use a TOV-initialdata-file!!! hint: use *.ini");

  /* read dim of INTERIOR, 
  This is the only one of the 2 domain interpolated on the grid */
  fgets ( line , 1024 , fid );
  sscanf( line, "#dim=%s", svalue);
  dim=atoi( svalue );
  
  /* dummy line */
  fgets ( line , 1024 , fid );
  sscanf( line, "#cols=%s", svalue);

  /* read mass ... */
  fgets ( line , 1024 , fid );
  sscanf( line, "M=%s", svalue);
  M=atof( svalue );

  /* ... and radius */
  fgets ( line , 1024 , fid );
  sscanf( line, "R=%s", svalue);
  R=atof( svalue );

  /* dummy line */
  fgets ( line , 1024 , fid );
  sscanf( line, "#%s", svalue);

  /* temp vars on spectral grid */
  rs    = (double *) malloc ( dim * sizeof(double) );
  lap   = (double *) malloc ( dim * sizeof(double) );
  Dlap  = (double *) malloc ( dim * sizeof(double) );
  psi4_ = (double *) malloc ( dim * sizeof(double) );
  Dpsi4 = (double *) malloc ( dim * sizeof(double) );
  rho_  = (double *) malloc ( dim * sizeof(double) );
  pres  = (double *) malloc ( dim * sizeof(double) );
  epsl_ = (double *) malloc ( dim * sizeof(double) );
  
  
  for( ix = 0; ix < dim; ix++ ) {
    fgets ( line , 1024 , fid );
    sscanf( line, "%lf %lf %lf %lf %lf %lf %lf %lf\n", 
            &rs[ix], &lap[ix], &Dlap[ix], &psi4_[ix], &Dpsi4[ix], &rho_[ix], &pres[ix], &epsl_[ix] );
  }

  fclose( fid );
  
  chebycoe_Extremes( lap,   dim, 0 );
  chebycoe_Extremes( Dlap,  dim, 0 );
  chebycoe_Extremes( psi4_, dim, 0 );
  chebycoe_Extremes( Dpsi4, dim, 0 );
  chebycoe_Extremes( rho_,  dim, 0 );
  chebycoe_Extremes( pres,  dim, 0 );
  chebycoe_Extremes( epsl_, dim, 0 );
/* until here *****************************************************************/

  
  // bam grid stuff
  double r,r1,r2, x,y,z;
  
  double *rho   = Ptr(level, "grhd_rho");
  double *epsl  = Ptr(level, "grhd_epsl");
  double *alpha = Ptr(level, "alpha");
  double *xp    = Ptr(level, "x");
  double *yp    = Ptr(level, "y");
  double *zp    = Ptr(level, "z");
  
  double *psi4    = Ptr(level, "BBID_psi4");
  double *dalphax = Ptr(level, "BBID_dalphax");
  double *dalphay = Ptr(level, "BBID_dalphay");
  double *dalphaz = Ptr(level, "BBID_dalphaz");
  double *dpsi4x  = Ptr(level, "BBID_dpsi4x");
  double *dpsi4y  = Ptr(level, "BBID_dpsi4y");
  double *dpsi4z  = Ptr(level, "BBID_dpsi4z");
  
  double LAMBDA[4][4],LAMBDAi[4][4];
  setLAMBDA(LAMBDA,LAMBDAi, mx,my,mz);

  
  // set TOV to grid
  forallpoints_ijk(level) {
    
    
    // set transformed coordinates
    x = LAMBDA[1][1]*(xp[ijk]-px) + LAMBDA[1][2]*(yp[ijk]-py) + LAMBDA[1][3]*(zp[ijk]-pz);
    y = LAMBDA[2][1]*(xp[ijk]-px) + LAMBDA[2][2]*(yp[ijk]-py) + LAMBDA[2][3]*(zp[ijk]-pz);
    z = LAMBDA[3][1]*(xp[ijk]-px) + LAMBDA[3][2]*(yp[ijk]-py) + LAMBDA[3][3]*(zp[ijk]-pz);
    
    r = sqrt(x*x+y*y+z*z);
    
    
    // set derivatives which are needed later
    if (r<R) {
      psi4[ijk]    = chebyeva( 0.0, 1.0, psi4_, dim, r/R );
      alpha[ijk]   = chebyeva( 0.0, 1.0, lap,   dim, r/R );
      dalphax[ijk] = chebyeva( 0.0, 1.0, Dlap,  dim, r/R )*x/r;
      dalphay[ijk] = chebyeva( 0.0, 1.0, Dlap,  dim, r/R )*y/r;
      dalphaz[ijk] = chebyeva( 0.0, 1.0, Dlap,  dim, r/R )*z/r;
      dpsi4x[ijk]  = chebyeva( 0.0, 1.0, Dpsi4, dim, r/R )*x/r;
      dpsi4y[ijk]  = chebyeva( 0.0, 1.0, Dpsi4, dim, r/R )*y/r;
      dpsi4z[ijk]  = chebyeva( 0.0, 1.0, Dpsi4, dim, r/R )*z/r;
      rho[ijk]     = chebyeva( 0.0, 1.0, rho_,  dim, r/R );
      epsl[ijk]    = chebyeva( 0.0, 1.0, epsl_, dim, r/R );
    } else {
      psi4[ijk]    = pow(1.0+0.5*M/r, 4.);
      alpha[ijk]   = (1.0-0.5*M/r)/(1.0+0.5*M/r);
      dalphax[ijk] = M/((0.5*M+r)*(0.5*M+r))*x/r;
      dalphay[ijk] = M/((0.5*M+r)*(0.5*M+r))*y/r;
      dalphaz[ijk] = M/((0.5*M+r)*(0.5*M+r))*z/r;
      dpsi4x[ijk]  = -2.*M*pow(1.+0.5*M/r,3)/(r*r)*x/r;
      dpsi4y[ijk]  = -2.*M*pow(1.+0.5*M/r,3)/(r*r)*y/r;
      dpsi4z[ijk]  = -2.*M*pow(1.+0.5*M/r,3)/(r*r)*z/r;
      rho[ijk]     = 0.;
      epsl[ijk]    = 0.;
    }
    
  } endfor_ijk;
  
  free(rs); 
  free(lap); 
  free(Dlap); 
  free(psi4_); 
  free(Dpsi4);
  free(rho_); 
  free(pres); 
  free(epsl_);
  
}

/* set TOV by computing it on the fly */
void    set_TOV(tL* level, int freemem,
                double px,double py,double pz, 
                double mx,double my,double mz,
                double hc,double rhoc,double R0)
{


  // 14.01.2013 Now we pass the follwoing pars
  /* double hc   = Getd("BBID_tov_hc"); */
  /* double rhoc = Getd("BBID_tov_rhoc"); */
  /* double R0   = Getd("BBID_tov_R0"); */

  int useh    = Getv("BBID_tov_integrate","h");
  static double *p_r,*p_m,*p_h,*p_rho,*p_pre,*p_phi,*p_riso;
  
  /* alloc memory the first time */
  static int npts;
/* static int firstcall = 1;*/
  if (/*firstcall &&*/ freemem==0) {
    npts    = Geti("BBID_tov_npts");
    
    /* solve TOV only once */
    if (useh) {
      
      p_r   = (double*) malloc (npts*sizeof(double));
      p_m   = (double*) malloc (npts*sizeof(double));
      p_h   = (double*) malloc (npts*sizeof(double));
      p_rho = (double*) malloc (npts*sizeof(double));
      p_pre = (double*) malloc (npts*sizeof(double));
      p_phi = (double*) malloc (npts*sizeof(double));
      p_riso= (double*) malloc (npts*sizeof(double));
      
      tov_h(hc, npts,
	    p_r, p_m, p_h, 
	    p_rho, p_pre, p_phi,
	    p_riso);
      
    } else {
      
      /* size of the arrey is not yet known -> is handled automatically 
	 prototypes are different ! */

      double Lamb;

      tov_r(rhoc, R0, &npts, 
	    &p_r, &p_m, &p_rho, 
	    &p_pre, &p_phi,
	    &p_riso);//, &Lamb);
      
      
    } 
  }
  
  /* free memory and go out */
  if (freemem) {
    printf("  free TOV star memory\n");
    if (p_r) {
      free(p_r);
      free(p_m);
      free(p_rho);
      free(p_pre);
      free(p_phi);
      free(p_riso);
    }
    if (p_h) 
      free(p_h);
    return ;
  }
  
  
  
  
  /* bam grid stuff */
  double M = p_m[npts-1];
  double R = p_riso[npts-1];
  
  
  double r,r1,r2, x,y,z, T;
  double phi,dphi,rsch,drsch,e,pres,dummy,h;
  
  double *rho   = Ptr(level, "grhd_rho");
  double *epsl  = Ptr(level, "grhd_epsl");
  double *vx    = Ptr(level, "grhd_vx");
  double *vy    = Ptr(level, "grhd_vy");
  double *vz    = Ptr(level, "grhd_vz");
  double p;
  double *alpha = Ptr(level, "alpha");
  double *xp    = Ptr(level, "x");
  double *yp    = Ptr(level, "y");
  double *zp    = Ptr(level, "z");
  
  double *psi4    = Ptr(level, "BBID_psi4");
  double *dalphax = Ptr(level, "BBID_dalphax");
  double *dalphay = Ptr(level, "BBID_dalphay");
  double *dalphaz = Ptr(level, "BBID_dalphaz");
  double *dpsi4x  = Ptr(level, "BBID_dpsi4x");
  double *dpsi4y  = Ptr(level, "BBID_dpsi4y");
  double *dpsi4z  = Ptr(level, "BBID_dpsi4z");
  
  double LAMBDA[4][4],LAMBDAi[4][4];
  setLAMBDA(LAMBDA,LAMBDAi, mx,my,mz);
  
  
  int l         = Geti("BBID_tov_perturb_l");
  int m         = Geti("BBID_tov_perturb_m");
  double n      = Getd("BBID_tov_perturb_n");
  double lambda = Getd("BBID_tov_perturb_lambda");
  double lambda2= Getd("BBID_tov_perturb_lambda2");
  int perturb   = 1*Getv("BBID_tov_perturb","p") +
                  2*Getv("BBID_tov_perturb","v") +
                  3*Getv("BBID_tov_perturb","vrot");
  
  if ((perturb) && (level->l==0)) {
    printf("Perturbartion parameter:\n");
    if (perturb==1) printf("  dp      = +lambda * (p+rho*(1+epsl)) * sin((n+1)*pi*r/(2R_t)) * Ylm\n");
    if (perturb==2) printf("  dv      = -lambda * sin((n+1)*pi*r/(2R_t))\n");
    if (perturb==3) printf("  dvrot   = +lambda2*r/R\n");
    printf("  Ylm     = Y%d%d\n",l,m);
    printf("  n       = %e\n",n);
    printf("  lambda  = %e\n",lambda);
    printf("  lambda2 = %e\n",lambda2);
  }
  
  // set TOV to grid
  forallpoints_ijk(level) {
    
    // set transformed coordinates
    x = LAMBDA[1][1]*(xp[ijk]-px) + LAMBDA[1][2]*(yp[ijk]-py) + LAMBDA[1][3]*(zp[ijk]-pz);
    y = LAMBDA[2][1]*(xp[ijk]-px) + LAMBDA[2][2]*(yp[ijk]-py) + LAMBDA[2][3]*(zp[ijk]-pz);
    z = LAMBDA[3][1]*(xp[ijk]-px) + LAMBDA[3][2]*(yp[ijk]-py) + LAMBDA[3][3]*(zp[ijk]-pz);
    
    r = sqrt(x*x+y*y+z*z);
    
    // set derivatives which are needed later
    if (r<R) {
      
      interp_lag4(p_r, p_riso, npts, r, 
                  &rsch,&drsch,  &dummy);
      interp_lag4(p_phi, p_riso, npts, r, 
                  &phi,&dphi,  &dummy);
     
      alpha[ijk]   = exp(phi);
      dalphax[ijk] = exp(phi)*dphi * x/r;
      dalphay[ijk] = exp(phi)*dphi * y/r;
      dalphaz[ijk] = exp(phi)*dphi * z/r;
      psi4[ijk]    = pow(rsch/r,2.);
      dpsi4x[ijk]  = 2.* ( (rsch*drsch)/(r*r) - (rsch*rsch)/(r*r*r) ) * x/r;
      dpsi4y[ijk]  = 2.* ( (rsch*drsch)/(r*r) - (rsch*rsch)/(r*r*r) ) * y/r;
      dpsi4z[ijk]  = 2.* ( (rsch*drsch)/(r*r) - (rsch*rsch)/(r*r*r) ) * z/r;
      
      if (useh) {
        interp_lag4(p_h, p_riso, npts, r, 
                    &h,  &dummy,&dummy);
        EOS.comp("H","","","pre","","",
                 h, &(p),&(rho[ijk]),&(epsl[ijk]));
      } else {
        interp_lag4(p_rho, p_riso, npts, r, 
                    &rho[ijk],  &dummy,&dummy);
        EOS.comp("r","","","pe","","",
                  rho[ijk], &(p),&(epsl[ijk]));
      }
      
    } else {
      alpha[ijk]   = (1.0-0.5*M/r)/(1.0+0.5*M/r);
      dalphax[ijk] = M/((0.5*M+r)*(0.5*M+r))*x/r;
      dalphay[ijk] = M/((0.5*M+r)*(0.5*M+r))*y/r;
      dalphaz[ijk] = M/((0.5*M+r)*(0.5*M+r))*z/r;
      psi4[ijk]    = pow(1.0+0.5*M/r, 4.);
      dpsi4x[ijk]  = -2.*M*pow(1.+0.5*M/r,3.)/(r*r)*x/r;
      dpsi4y[ijk]  = -2.*M*pow(1.+0.5*M/r,3.)/(r*r)*y/r;
      dpsi4z[ijk]  = -2.*M*pow(1.+0.5*M/r,3.)/(r*r)*z/r;
      rho[ijk]     = 0.;
      epsl[ijk]    = 0.;
    }
    
    /* do a perturbation */
    if (r<R && perturb) {
      
      double phi,theta, rxy;
      double Yr,Yi,Ylm,Hl0,dp,dv,dvr;
      
      rxy     = sqrt(x*x + y*y);
      phi     = acos(x/rxy);
      phi     = (y>=0.)? (phi) : (2.*PI-phi);
      theta   = acos(z/r);
      
      // using same convention as in http://de.arxiv.org/abs/0808.4002
      SphericalHarmonicY( &Yr,&Yi, l,m, phi,cos(theta));
      Ylm      = Yr;
      Hl0      = lambda * sin((n+1.)*PI*r/(2.*R));
      dp       = (p+rho[ijk]*(1.+epsl[ijk])) * Hl0 * Ylm;
      dv       = -Hl0;
      dvr      = lambda2*r/R;
      
      // set modii
      if (perturb == 1) {
        p       += dp;
        p        = DMAX(0., p);
        rho[ijk] = pow( (p/EOS.K) ,  1./EOS.GAMMA );
      } else if (perturb == 2) {
        vx[ijk] += x/r * dv;
        vy[ijk] += y/r * dv;
        vz[ijk] += z/r * dv;
      } else if (perturb == 3) {
        vx[ijk] += x/r * dv - sin(phi) * dvr;
        vy[ijk] += y/r * dv + cos(phi) * dvr;
        vz[ijk] += z/r * dv;
      } else 
        errorexit("this perturbation method is not implemented");
    }
    
    
  } endfor_ijk;

}

/* set a puncture */
void    set_PUNC(tL* level, double M,
                 double px,double py,double pz, 
                 double mx,double my,double mz) 
{ 
  double r,r1,r2, x,y,z;
  
  double *rho,*epsl;
  if (useMATTER) {
    rho   = Ptr(level, "grhd_rho");
    epsl  = Ptr(level, "grhd_epsl");
  }
  double *alpha = Ptr(level, "alpha");
  double *xp    = Ptr(level, "x");
  double *yp    = Ptr(level, "y");
  double *zp    = Ptr(level, "z");
  
  double *psi4    = Ptr(level, "BBID_psi4");
  double *dalphax = Ptr(level, "BBID_dalphax");
  double *dalphay = Ptr(level, "BBID_dalphay");
  double *dalphaz = Ptr(level, "BBID_dalphaz");
  double *dpsi4x  = Ptr(level, "BBID_dpsi4x");
  double *dpsi4y  = Ptr(level, "BBID_dpsi4y");
  double *dpsi4z  = Ptr(level, "BBID_dpsi4z");
  
  double LAMBDA[4][4],LAMBDAi[4][4];
  setLAMBDA(LAMBDA,LAMBDAi, mx,my,mz);

  
  // set TOV to grid
  forallpoints_ijk(level) {
    
    
    // set transformed coordinates
    x = LAMBDA[1][1]*(xp[ijk]-px) + LAMBDA[1][2]*(yp[ijk]-py) + LAMBDA[1][3]*(zp[ijk]-pz);
    y = LAMBDA[2][1]*(xp[ijk]-px) + LAMBDA[2][2]*(yp[ijk]-py) + LAMBDA[2][3]*(zp[ijk]-pz);
    z = LAMBDA[3][1]*(xp[ijk]-px) + LAMBDA[3][2]*(yp[ijk]-py) + LAMBDA[3][3]*(zp[ijk]-pz);
    r = sqrt(x*x+y*y+z*z);
    
   
    psi4[ijk]    = pow(1.0+0.5*M/r,4);
    dalphax[ijk] = M/((0.5*M+r)*(0.5*M+r))*x/r;
    dalphay[ijk] = M/((0.5*M+r)*(0.5*M+r))*y/r;
    dalphaz[ijk] = M/((0.5*M+r)*(0.5*M+r))*z/r;
    dpsi4x[ijk]  = -2.*M*pow(1.+0.5*M/r,3)/(r*r)*x/r;
    dpsi4y[ijk]  = -2.*M*pow(1.+0.5*M/r,3)/(r*r)*y/r;
    dpsi4z[ijk]  = -2.*M*pow(1.+0.5*M/r,3)/(r*r)*z/r;
    
    alpha[ijk]  = 1;(1.0-0.5*M/r)/(1.0+0.5*M/r);
    
    if (useMATTER) {
      rho[ijk]    = 0.;
      epsl[ijk]   = 0.;
    }
      
  } endfor_ijk;
  
  
  
}


/* set a BH */
void    set_BH(tL* level, double M,
                 double px,double py,double pz, 
                 double mx,double my,double mz,
                 double sx,double sy,double sz) 
{ 

   tG *g = level->grid;
   int ijk,l;
   double r,r1,r2, x,y,z, oor, oor2, oor3, N[4], P[4], S[4], eNS[4], NP, psi4, psim2;
   double *xp    = Ptr(level, "x");
   double *yp    = Ptr(level, "y");
   double *zp    = Ptr(level, "z");
   double *alpha = Ptr(level, "BBID_Talpha");
   double *betax = Ptr(level, "BBID_Tbetax");
   double *betay = Ptr(level, "BBID_Tbetay");
   double *betaz = Ptr(level, "BBID_Tbetaz");
   double *gxx   = Ptr(level, "BBID_Tgxx");
   double *gxy   = Ptr(level, "BBID_Tgxy");
   double *gxz   = Ptr(level, "BBID_Tgxz");
   double *gyy   = Ptr(level, "BBID_Tgyy");
   double *gyz   = Ptr(level, "BBID_Tgyz");
   double *gzz   = Ptr(level, "BBID_Tgzz");

   double *Kxx   = Ptr(level, "BBID_TKxx");
   double *Kxy   = Ptr(level, "BBID_TKxy");
   double *Kxz   = Ptr(level, "BBID_TKxz");
   double *Kyy   = Ptr(level, "BBID_TKyy");
   double *Kyz   = Ptr(level, "BBID_TKyz");
   double *Kzz   = Ptr(level, "BBID_TKzz");

   double *rho,*epsl;
   double M2,px2,py2,pz2;
   

    if (useMATTER) {
      rho   = Ptr(level, "BBID_Trho");
      epsl  = Ptr(level, "BBID_Tepsl");
    }

  double *u   = Ptr(level, "incl_punctures_u");



  if (Getv("BBID_solve_external_method", "multigrid")) {

      M2 = Getd("mass1");
     px2 = Getd("px1");
     py2 = Getd("py1");
     pz2 = Getd("pz1");

    Setd("mass1", M);
    Setd("mass2",0);
    //Setd("mass2",M2);
    Setd("px1", px);
    Setd("px2", px2);
    Setd("py1", py);
    Setd("py2", py2);
    Setd("pz1", pz);
    Setd("pz2", pz2);

    Setd("mx1", mx);
    Setd("mx2", 0.0);
    Setd("my1", my);
    Setd("my2", 0.0);
    Setd("mz1", mz);
    Setd("mz2", 0.0);

    Setd("sx1", sx);
    Setd("sx2", 0.0);
    Setd("sy1", sy);
    Setd("sy2", 0.0);
    Setd("sz1", sz);
    Setd("sz2", 0.0);

    P[1] = mx;
    P[2] = my;
    P[3] = mz; 
 
    S[1] = sx;
    S[2] = sy;
    S[3] = sz;

     Sets("punctures_solver", "multigrid");
      PunctureData(level);
   }

  else if (Getv("BBID_solve_external_method", "spectral")) {
   
     Sets("punctures_solver", "spectral"); 
     printf("  set puncture solver to spectral \n");
         if(level->l == 0)  PunctureData(level);

    px   = Getd("px2");
    py   = Getd("py2");
    pz   = Getd("pz2");

    P[1] = Getd("mx2");
    P[2] = Getd("my2");
    P[3] = Getd("mz2");
 
    S[1] = Getd("sx2");
    S[2] = Getd("sy2");
    S[3] = Getd("sz2"); 
       M = Getd("mass2"); 

   } 

  else errorexit("This external method is not implemented!");

   forallpoints(level,ijk) {
    
      x = (xp[ijk]-px);
      y = (yp[ijk]-py);
      z = (zp[ijk]-pz);
    
      r = sqrt(x*x+y*y+z*z);
      oor= 1.0/r;
      oor2=oor*oor;
      oor3=oor2*oor;
    
      N[1] = x*oor;
      N[2] = y*oor;
      N[3] = z*oor;   

      NP = N[1]*P[1] + N[2]*P[2] + N[3]*P[3];
    
      eNS[1] = N[2]*S[3] - N[3]*S[2];
      eNS[2] = N[3]*S[1] - N[1]*S[3];
      eNS[3] = N[1]*S[2] - N[2]*S[1];

      psi4  = pow(1.0 + 0.5*M*oor + u[ijk], 4);
      psim2 = pow(1.0 + 0.5*M*oor + u[ijk],-2);

     /* BY-extrinsic curvature */
        Kxx[ijk] += psim2*(
        1.5*oor2 * (N[1]*P[1] + N[1]*P[1] - (1.0 - N[1]*N[1])*NP) -
        3.0*oor3 * (eNS[1]*N[1] + eNS[1]*N[1]));
        Kxy[ijk] += psim2*(
        1.5*oor2 * (N[1]*P[2] + N[2]*P[1] + N[1]*N[2]*NP) -
        3.0*oor3 * (eNS[1]*N[2] + eNS[2]*N[1]));
        Kxz[ijk] += psim2*(
        1.5*oor2 * (N[1]*P[3] + N[3]*P[1] + N[1]*N[3]*NP) -
        3.0*oor3 * (eNS[1]*N[3] + eNS[3]*N[1]));
        Kyy[ijk] += psim2*(
        1.5*oor2 * (N[2]*P[2] + N[2]*P[2] - (1.0 - N[2]*N[2])*NP) -
        3.0*oor3 * (eNS[2]*N[2] + eNS[2]*N[2]));
        Kyz[ijk] += psim2*(
        1.5*oor2 * (N[2]*P[3] + N[3]*P[2] + N[2]*N[3]*NP) -
        3.0*oor3 * (eNS[2]*N[3] + eNS[3]*N[2]));
        Kzz[ijk] += psim2*(
        1.5*oor2 * (N[3]*P[3] + N[3]*P[3] - (1.0 - N[3]*N[3])*NP) -
        3.0*oor3 * (eNS[3]*N[3] + eNS[3]*N[3]));

      alpha[ijk]  += 1;
      betax[ijk]  += 0;
      betay[ijk]  += 0;
      betaz[ijk]  += 0;
   
      gxx[ijk] += psi4;
      gxy[ijk] += 0;
      gxz[ijk] += 0;
      gyy[ijk] += psi4;
      gyz[ijk] += 0;
      gzz[ijk] += psi4;

        if (useMATTER) {
           rho[ijk]    += 0.;
           epsl[ijk]   += 0.;
         }    

    } //end loop over points

 }

/* set perturbation */
void    set_perturbation(tL* level, double R,
                         double px,double py,double pz, 
                         double mx,double my,double mz)
{
  double r, x,y,z;
  double phi,theta, rxy;
  double Yr,Yi,Ylm,Hl0,dp,dv,dvr;
  
  double *rho   = Ptr(level, "grhd_rho");
  double *epsl  = Ptr(level, "grhd_epsl");
  double *vx    = Ptr(level, "grhd_vx");
  double *vy    = Ptr(level, "grhd_vy");
  double *vz    = Ptr(level, "grhd_vz");
  double p;
  
  double *xp    = Ptr(level, "x");
  double *yp    = Ptr(level, "y");
  double *zp    = Ptr(level, "z");
  
  double LAMBDA[4][4],LAMBDAi[4][4];
  setLAMBDA(LAMBDA,LAMBDAi, mx,my,mz);
  
  int l         = Geti("BBID_tov_perturb_l");
  int m         = Geti("BBID_tov_perturb_m");
  double n      = Getd("BBID_tov_perturb_n");
  double lambda = Getd("BBID_tov_perturb_lambda");
  double lambda2= Getd("BBID_tov_perturb_lambda2");
  int perturb   = 1*Getv("BBID_tov_perturb","p") +
                  2*Getv("BBID_tov_perturb","v") +
                  3*Getv("BBID_tov_perturb","vrot");
  
  if (!perturb) return;
  
  if (level->l==0) {
      printf("Perturbartion parameter:\n");
      if (perturb==1) printf("  dp      = +lambda * (p+rho*(1+epsl)) * sin((n+1)*pi*r/(2R_t)) * Ylm\n");
      if (perturb==2) printf("  dv      = -lambda * sin((n+1)*pi*r/(2R_t))\n");
      if (perturb==3) printf("  dvrot   = +lambda2*r/R\n");
      printf("  Ylm     = Y%d%d\n",l,m);
      printf("  n       = %e\n",n);
      printf("  lambda  = %e\n",lambda);
      printf("  lambda2 = %e\n",lambda2);
  }
  
  // set perturbation to grid
  forallpoints_ijk(level) {
    
    // set transformed coordinates
    x = LAMBDA[1][1]*(xp[ijk]-px) + LAMBDA[1][2]*(yp[ijk]-py) + LAMBDA[1][3]*(zp[ijk]-pz);
    y = LAMBDA[2][1]*(xp[ijk]-px) + LAMBDA[2][2]*(yp[ijk]-py) + LAMBDA[2][3]*(zp[ijk]-pz);
    z = LAMBDA[3][1]*(xp[ijk]-px) + LAMBDA[3][2]*(yp[ijk]-py) + LAMBDA[3][3]*(zp[ijk]-pz);
    
    r = sqrt(x*x+y*y+z*z);
    
    /* do a perturbation */
    if (r<R) {
      
      rxy     = sqrt(x*x + y*y);
      phi     = acos(x/rxy);
      phi     = (y>=0.)? (phi) : (2.*PI-phi);
      theta   = acos(z/r);
      
      // using same convention as in http://de.arxiv.org/abs/0808.4002
      SphericalHarmonicY( &Yr,&Yi, l,m, phi,cos(theta));
      Ylm      = Yr;
      Hl0      = lambda * sin((n+1.)*PI*r/(2.*R));
      dp       = (p+rho[ijk]*(1.+epsl[ijk])) * Hl0 * Ylm;
      dv       = -Hl0;
      dvr      = lambda2*r/R;
      
      // set modii
      if (perturb == 1) {
        p       += dp;
        p        = DMAX(0., p);
        rho[ijk] = pow( (p/EOS.K) ,  1./EOS.GAMMA );
      } else if (perturb == 2) {
        vx[ijk] += x/r * dv;
        vy[ijk] += y/r * dv;
        vz[ijk] += z/r * dv;
      } else if (perturb == 3) {
        vx[ijk] += x/r * dv - sin(phi) * dvr;
        vy[ijk] += y/r * dv + cos(phi) * dvr;
        vz[ijk] += z/r * dv;
      } else 
        errorexit("this perturbation method is not implemented");
    }
    
    
  } endfor_ijk;

}







/* do a general boost and add it to a temporary buffer */
void    boost_spacetime(tL *level, double mx,double my,double mz)
{
  if (mx!=0. || my!=0. || mz!=0.) 
    printf("boost spacetime star on level %d   (%f %f %f)\n",level->l, mx,my,mz);
  
  int o,p,q,r,s,t;
  double g[4][4],  u[4],  delg[4][4][4],  Gamma[4][4][4];  // BEFORE BOOST
  double gp[4][4], up[4], delgp[4][4][4], Gammap[4][4][4]; // AFTER  BOOST (p for primed)
  double gip[4][4];
  double W,v2;
    
  double LAMBDA[4][4],LAMBDAi[4][4];
  setLAMBDA(LAMBDA,LAMBDAi, mx,my,mz);
  
  double *alpha = Ptr(level, "alpha");
  double *psi4  = Ptr(level, "BBID_psi4");
  double *dalpha[4],*dpsi4[4];
  dalpha[0] = NULL;
  dalpha[1] = Ptr(level, "BBID_dalphax");
  dalpha[2] = Ptr(level, "BBID_dalphay");
  dalpha[3] = Ptr(level, "BBID_dalphaz");
  dpsi4[0]  = NULL;
  dpsi4[1]  = Ptr(level, "BBID_dpsi4x");
  dpsi4[2]  = Ptr(level, "BBID_dpsi4y");
  dpsi4[3]  = Ptr(level, "BBID_dpsi4z");
  double *betax = Ptr(level, "betax");
  double *betay = Ptr(level, "betay");
  double *betaz = Ptr(level, "betaz");
  double *gxx   = Ptr(level, "adm_gxx");
  double *gxy   = Ptr(level, "adm_gxy");
  double *gxz   = Ptr(level, "adm_gxz");
  double *gyy   = Ptr(level, "adm_gyy");
  double *gyz   = Ptr(level, "adm_gyz");
  double *gzz   = Ptr(level, "adm_gzz");
  double *rho,*epsl,*vx,*vy,*vz;
  if (useMATTER) {
    rho   = Ptr(level, "grhd_rho");
    epsl  = Ptr(level, "grhd_epsl");
    vx    = Ptr(level, "grhd_vx");
    vy    = Ptr(level, "grhd_vy");
    vz    = Ptr(level, "grhd_vz");
  }
  double *Kxx   = Ptr(level, "adm_Kxx");
  double *Kxy   = Ptr(level, "adm_Kxy");
  double *Kxz   = Ptr(level, "adm_Kxz");
  double *Kyy   = Ptr(level, "adm_Kyy");
  double *Kyz   = Ptr(level, "adm_Kyz");
  double *Kzz   = Ptr(level, "adm_Kzz");
  
  double *Talpha = Ptr(level, "BBID_Talpha");
  double *Tbetax = Ptr(level, "BBID_Tbetax");
  double *Tbetay = Ptr(level, "BBID_Tbetay");
  double *Tbetaz = Ptr(level, "BBID_Tbetaz");
  double *Tgxx   = Ptr(level, "BBID_Tgxx");
  double *Tgxy   = Ptr(level, "BBID_Tgxy");
  double *Tgxz   = Ptr(level, "BBID_Tgxz");
  double *Tgyy   = Ptr(level, "BBID_Tgyy");
  double *Tgyz   = Ptr(level, "BBID_Tgyz");
  double *Tgzz   = Ptr(level, "BBID_Tgzz");
  double *Trho   = Ptr(level, "BBID_Trho");
  double *Tepsl  = Ptr(level, "BBID_Tepsl");
  double *Tvx    = Ptr(level, "BBID_Tvx");
  double *Tvy    = Ptr(level, "BBID_Tvy");
  double *Tvz    = Ptr(level, "BBID_Tvz");
  double *TKxx   = Ptr(level, "BBID_TKxx");
  double *TKxy   = Ptr(level, "BBID_TKxy");
  double *TKxz   = Ptr(level, "BBID_TKxz");
  double *TKyy   = Ptr(level, "BBID_TKyy");
  double *TKyz   = Ptr(level, "BBID_TKyz");
  double *TKzz   = Ptr(level, "BBID_TKzz");
  
  
  // set boost for each grid point
  forallpoints_ijk(level) {
    
    // do we need to boost?
    // for analytical issues see src/matter/doc/matter.tex 
    if (mx!=0. || my!=0. || mz!=0.) {
    
      // set T_munu ... but better to use/transform only u^mu
      // rho p epsl are invariant under transformation
      if (useMATTER) {
        v2 = 0.;
        W  = 1.0/sqrt(1.0 - v2);
        u[0] = W/alpha[ijk];
        u[1] = u[2] = u[3] = 0.;
      }
      
      // set g_\mu\nu
      g[0][0] = -(alpha[ijk]*alpha[ijk]);
      for (o=1; o<=3; o++)
        g[0][o] = g[o][0] = 0.;
      for (o=1; o<=3; o++)
      for (p=1; p<=3; p++)
        g[o][p] = psi4[ijk] * (double)(o==p);
      
      // time derivative of g_\mu\nu
      for (o=0; o<=3; o++)
      for (p=0; p<=3; p++)
        delg[0][o][p] = 0.;
          // spacial derivative of g_\mu\nu
      for (o=1; o<=3; o++) {
              // lapse
        delg[o][0][0] = -2.*alpha[ijk] * dalpha[o][ijk];
              // shift
        for (p=1; p<=3; p++)
          delg[o][p][0] = delg[o][0][p] = 0.;
              // metric
        for (p=1; p<=3; p++)
        for (q=1; q<=3; q++)
          delg[o][p][q] = dpsi4[o][ijk] * (double)(q==p);
      }
      
      // constract the 4 Gamma
      Gamma44(g,delg, Gamma);
      
      
      
      // make coordinate transformation for u^\mu, g_\mu\nu,  Gamma^\sigma_\mu\nu
      for (o=0; o<=3; o++)
      for (p=0; p<=3; p++)
      for (q=0; q<=3; q++) {
        up[o] = 0.;
        gp[o][p] = 0.;
        Gammap[o][p][q] = 0.;
        for (r=0; r<=3; r++) {
          up[o] += LAMBDAi[r][o] * u[r];
          for (s=0; s<=3; s++) {
            gp[o][p] += LAMBDA[o][r] * LAMBDA[p][s] * g[r][s];
            for (t=0; t<=3; t++) {
              Gammap[o][p][q] += LAMBDAi[o][r] * LAMBDA[p][s] * LAMBDA[q][t] * Gamma[r][s][t];
              Gammap[o][p][q] += LAMBDAi[o][r] * 0.; // lucky coincidence ... 0 because Lambda is constant
            }
          }
        }
      }
      
      
      
      // compute inverse gp^{\mu\nu}
      invg4(gp,gip);
      
      // set transformed g_munu values to grid
      gxx[ijk]    = gp[1][1];
      gxy[ijk]    = gp[1][2];
      gxz[ijk]    = gp[1][3];
      gyy[ijk]    = gp[2][2];
      gyz[ijk]    = gp[2][3];
      gzz[ijk]    = gp[3][3];
      betax[ijk]  = -gip[0][1]/gip[0][0];
      betay[ijk]  = -gip[0][2]/gip[0][0];
      betaz[ijk]  = -gip[0][3]/gip[0][0];
      alpha[ijk]  = sqrt(-1./gip[0][0]);
      
      // set transformed u^mu values to grid
      if (useMATTER) {
        W        = up[0]*alpha[ijk];
        vx[ijk]  = up[1]/W + betax[ijk]/alpha[ijk];
        vy[ijk]  = up[2]/W + betay[ijk]/alpha[ijk];
        vz[ijk]  = up[3]/W + betaz[ijk]/alpha[ijk];
      }
      
      // set K_ij 
      Kxx[ijk] = -alpha[ijk] * Gammap[0][1][1];
      Kxy[ijk] = -alpha[ijk] * Gammap[0][1][2];
      Kxz[ijk] = -alpha[ijk] * Gammap[0][1][3];
      Kyy[ijk] = -alpha[ijk] * Gammap[0][2][2];
      Kyz[ijk] = -alpha[ijk] * Gammap[0][2][3];
      Kzz[ijk] = -alpha[ijk] * Gammap[0][3][3];
      
    } else {
      
      gxx[ijk]    = psi4[ijk];
      gxy[ijk]    = 0.;
      gxz[ijk]    = 0.;
      gyy[ijk]    = psi4[ijk];
      gyz[ijk]    = 0.;
      gzz[ijk]    = psi4[ijk];
      
      
      
    }
    
    
    
    // copy/add stuff
    Talpha[ijk] += alpha[ijk];
    Tbetax[ijk] += betax[ijk];
    Tbetay[ijk] += betay[ijk];
    Tbetaz[ijk] += betaz[ijk];
    Tgxx[ijk]   += gxx[ijk];
    Tgxy[ijk]   += gxy[ijk];
    Tgxz[ijk]   += gxz[ijk];
    Tgyy[ijk]   += gyy[ijk];
    Tgyz[ijk]   += gyz[ijk];
    Tgzz[ijk]   += gzz[ijk];
    TKxx[ijk]   += Kxx[ijk];
    TKxy[ijk]   += Kxy[ijk];
    TKxz[ijk]   += Kxz[ijk];
    TKyy[ijk]   += Kyy[ijk];
    TKyz[ijk]   += Kyz[ijk];
    TKzz[ijk]   += Kzz[ijk];
    if (useMATTER) {
      Trho[ijk]   += rho[ijk];
      Tepsl[ijk]  += epsl[ijk];
      Tvx[ijk]    += vx[ijk];
      Tvy[ijk]    += vy[ijk];
      Tvz[ijk]    += vz[ijk];
      CheckForNANandINF(3, vx[ijk],vy[ijk],vz[ijk]);
    }
    
    // check ... to be sure
    CheckForNANandINF(16, alpha[ijk],betax[ijk],betay[ijk],betaz[ijk],
                      gxx[ijk],gxy[ijk],gxz[ijk],gyy[ijk],gyz[ijk],gzz[ijk], 
                      Kxx[ijk],Kxy[ijk],Kxz[ijk],Kyy[ijk],Kyz[ijk],Kzz[ijk]);
    
  } endfor_ijk;
}

/* all is added ... do we need to substract something... yes: */
void    add_spacetime(tL* level, int subMinkowsky)
{
  
  double *psi   = Ptr(level, "adm_psi");
  
  double *alpha = Ptr(level, "alpha");
  double *betax = Ptr(level, "betax");
  double *betay = Ptr(level, "betay");
  double *betaz = Ptr(level, "betaz");
  double *gxx   = Ptr(level, "adm_gxx");
  double *gxy   = Ptr(level, "adm_gxy");
  double *gxz   = Ptr(level, "adm_gxz");
  double *gyy   = Ptr(level, "adm_gyy");
  double *gyz   = Ptr(level, "adm_gyz");
  double *gzz   = Ptr(level, "adm_gzz");
  double *Kxx   = Ptr(level, "adm_Kxx");
  double *Kxy   = Ptr(level, "adm_Kxy");
  double *Kxz   = Ptr(level, "adm_Kxz");
  double *Kyy   = Ptr(level, "adm_Kyy");
  double *Kyz   = Ptr(level, "adm_Kyz");
  double *Kzz   = Ptr(level, "adm_Kzz");
  
  double *rho,*epsl,*vx,*vy,*vz;
  if (useMATTER) {
    rho   = Ptr(level, "grhd_rho");
    epsl  = Ptr(level, "grhd_epsl");
    vx    = Ptr(level, "grhd_vx");
    vy    = Ptr(level, "grhd_vy");
    vz    = Ptr(level, "grhd_vz");
  }
  
  double *Talpha = Ptr(level, "BBID_Talpha");
  double *Tbetax = Ptr(level, "BBID_Tbetax");
  double *Tbetay = Ptr(level, "BBID_Tbetay");
  double *Tbetaz = Ptr(level, "BBID_Tbetaz");
  double *Tgxx   = Ptr(level, "BBID_Tgxx");
  double *Tgxy   = Ptr(level, "BBID_Tgxy");
  double *Tgxz   = Ptr(level, "BBID_Tgxz");
  double *Tgyy   = Ptr(level, "BBID_Tgyy");
  double *Tgyz   = Ptr(level, "BBID_Tgyz");
  double *Tgzz   = Ptr(level, "BBID_Tgzz");
  double *Trho   = Ptr(level, "BBID_Trho");
  double *Tepsl  = Ptr(level, "BBID_Tepsl");
  double *Tvx    = Ptr(level, "BBID_Tvx");
  double *Tvy    = Ptr(level, "BBID_Tvy");
  double *Tvz    = Ptr(level, "BBID_Tvz");
  double *TKxx   = Ptr(level, "BBID_TKxx");
  double *TKxy   = Ptr(level, "BBID_TKxy");
  double *TKxz   = Ptr(level, "BBID_TKxz");
  double *TKyy   = Ptr(level, "BBID_TKyy");
  double *TKyz   = Ptr(level, "BBID_TKyz");
  double *TKzz   = Ptr(level, "BBID_TKzz");
 //if(Getv("physics","punctures")){
 if(0){
 forallpoints_ijk(level) {
    
    psi[ijk]   = 1.;
    
    alpha[ijk] = Talpha[ijk] + alpha[ijk]; 
    betax[ijk] = Tbetax[ijk] + betax[ijk]; 
    betay[ijk] = Tbetay[ijk] + betay[ijk]; 
    betaz[ijk] = Tbetaz[ijk] + betaz[ijk]; 
    
    gxx[ijk]   = Tgxx[ijk] + gxx[ijk]; 
    gxy[ijk]   = Tgxy[ijk] + gxy[ijk]; 
    gxz[ijk]   = Tgxz[ijk] + gxz[ijk]; 
    gyy[ijk]   = Tgyy[ijk] + gyy[ijk]; 
    gyz[ijk]   = Tgyz[ijk] + gyz[ijk]; 
    gzz[ijk]   = Tgzz[ijk] + gzz[ijk]; 
    
    Kxx[ijk]   = TKxx[ijk] + Kxx[ijk]; 
    Kxy[ijk]   = TKxy[ijk] + Kxy[ijk]; 
    Kxz[ijk]   = TKxz[ijk] + Kxz[ijk]; 
    Kyy[ijk]   = TKyy[ijk] + Kyy[ijk]; 
    Kyz[ijk]   = TKyz[ijk] + Kyz[ijk]; 
    Kzz[ijk]   = TKzz[ijk] + Kzz[ijk]; 
    
    if (useMATTER) {
      rho[ijk]   = Trho[ijk] + rho[ijk];
      epsl[ijk]  = Tepsl[ijk] + epsl[ijk];
      vx[ijk]    = Tvx[ijk] + vx[ijk];
      vy[ijk]    = Tvy[ijk] + vy[ijk];
      vz[ijk]    = Tvz[ijk] + vz[ijk];
    }
    if (subMinkowsky) {
      alpha[ijk] -= 1.;
      gxx[ijk]   -= 1.;
      gyy[ijk]   -= 1.;
      gzz[ijk]   -= 1.;
    }
    
  } endfor_ijk;

 }
 
 else { 
  // set boost for each grid point
  forallpoints_ijk(level) {
    
    psi[ijk]   = 1.;
    
    alpha[ijk] = Talpha[ijk]; 
    betax[ijk] = Tbetax[ijk]; 
    betay[ijk] = Tbetay[ijk]; 
    betaz[ijk] = Tbetaz[ijk]; 
    
    gxx[ijk]   = Tgxx[ijk]; 
    gxy[ijk]   = Tgxy[ijk]; 
    gxz[ijk]   = Tgxz[ijk]; 
    gyy[ijk]   = Tgyy[ijk]; 
    gyz[ijk]   = Tgyz[ijk]; 
    gzz[ijk]   = Tgzz[ijk]; 
    
    Kxx[ijk]   = TKxx[ijk]; 
    Kxy[ijk]   = TKxy[ijk]; 
    Kxz[ijk]   = TKxz[ijk]; 
    Kyy[ijk]   = TKyy[ijk]; 
    Kyz[ijk]   = TKyz[ijk]; 
    Kzz[ijk]   = TKzz[ijk]; 
    
    if (useMATTER) {
      rho[ijk]   = Trho[ijk];
      epsl[ijk]  = Tepsl[ijk];
      vx[ijk]    = Tvx[ijk];
      vy[ijk]    = Tvy[ijk];
      vz[ijk]    = Tvz[ijk];
    }
    if (subMinkowsky) {
      alpha[ijk] -= 1.;
      gxx[ijk]   -= 1.;
      gyy[ijk]   -= 1.;
      gzz[ijk]   -= 1.;
    }
    
   } endfor_ijk;
  }
}





/* superpose two TOVs ... maybe this can be easily done with a puncture */
int     set_binary_boost(tL* level) 
{  
   tG *g = level->grid;
  if (Getv("BBID_modus","NS") || Getv("BBID_modus","NSNS") || Getv("BBID_modus","NSBH")) 
    useMATTER++;

  // for matter
  adm_Minkowski(level); 
  if (useMATTER) {
    enablevar(level, Ind("adm_rho"));
    enablevar(level, Ind("adm_Sx"));
    enablevar(level, Ind("adm_SSxx"));
    enablevar(level, Ind("adm_ST"));
    enablevar(level, Ind("grhd_epsl"));
    enablevar(level, Ind("grhd_rho"));
    enablevar(level, Ind("grhd_vx"));
  }
  
  // for the boost
  enablevar(level, Ind("BBID_psi4"));
  enablevar(level, Ind("BBID_dalphax"));
  enablevar(level, Ind("BBID_dpsi4x"));
  
  // set these variables to zero at the beginning
  tVarList* vl = vlalloc(level);
  vlpush(vl, Ind("BBID_psi4"));
  vlpush(vl, Ind("BBID_dalphax"));
  vlpush(vl, Ind("BBID_dpsi4x"));
  vlpush(vl, Ind("BBID_Talpha"));
  vlpush(vl, Ind("BBID_Tbetax"));
  vlpush(vl, Ind("BBID_Tgxx"));
  vlpush(vl, Ind("BBID_TKxx"));
  if (useMATTER) {
    vlpush(vl, Ind("BBID_Trho"));
    vlpush(vl, Ind("BBID_Tepsl"));
    vlpush(vl, Ind("BBID_Tvx"));
  }
  enablevarlist(vl);
  vlsetconstant(vl, 0.0);
  
  int sub = 0;
  

    /* first star ****************************************************************/
  if (Getv("BBID_modus","NS") || Getv("BBID_modus","NSNS") || Getv("BBID_modus","NSBH")) {
    printf("set TOV star 1 on level %d\n",level->l);
    
    if (ExistPar("initial_file")) {
      
      read_TOV(level, Gets("initial_file"), 
               Getd("px1"),Getd("py1"),Getd("pz1"),
               Getd("mx1"),Getd("my1"),Getd("mz1"));
    } else {
      
      set_TOV(level,0,
              Getd("px1"),Getd("py1"),Getd("pz1"),
              Getd("mx1"),Getd("my1"),Getd("mz1"),
              Getd("BBID_tov_hc1"),Getd("BBID_tov_rhoc1"),Getd("BBID_tov_R01"));
    }
    
  } else {
    printf("set BH 1 on level %d\n",level->l);
    
    set_PUNC(level,Getd("mass1"), 
             Getd("px1"),Getd("py1"),Getd("pz1"),
             Getd("mx1"),Getd("my1"),Getd("mz1"));
    
  }
  
  boost_spacetime(level, Getd("mx1"),Getd("my1"),Getd("mz1"));
  
  
   /* second star ***************************************************************/
  if (Geti("nobjects")>1 && Getd("mass2")>0.) {
    
    if (Getv("BBID_modus","NSNS") || Getv("BBID_modus","BHNS")) {
      printf("set TOV star 2 on level %d\n",level->l);
      
      if (ExistPar("initial_file")) {
        
        if (ExistPar("initial_file2"))
          Sets("initial_file",Gets("initial_file2"));
        
        read_TOV(level, Gets("initial_file"), 
                 Getd("px2"),Getd("py2"),Getd("pz2"),
                 Getd("mx2"),Getd("my2"),Getd("mz2"));
        
      } else {
        
        set_TOV(level, 0,
                Getd("px2"),Getd("py2"),Getd("pz2"),
                Getd("mx2"),Getd("my2"),Getd("mz2"),
              Getd("BBID_tov_hc2"),Getd("BBID_tov_rhoc2"),Getd("BBID_tov_R02"));
      }
           
        boost_spacetime(level, Getd("mx2"),Getd("my2"),Getd("mz2"));
    
     } else { if(Getv("BBID_solve_external","yes")) {
      printf("compute BH 2 on level %d\n with another initial data project\n",level->l);
 
    
  double M1[4], P1[4], S1[4], M2[4], P2[4], S2[4], MASS1, MASS2;

    MASS1    = Getd("mass1");
    MASS2    = Getd("mass2");

    /*save the parameter of the initial data*/
    P1[1] = Getd("px1");
    P1[2] = Getd("py1");
    P1[3] = Getd("pz1"); 
    P2[1] = Getd("px2");
    P2[2] = Getd("py2");
    P2[3] = Getd("pz2");  

    M1[1] = Getd("mx1");
    M1[2] = Getd("my1");
    M1[3] = Getd("mz1"); 
    M2[1] = Getd("mx2");
    M2[2] = Getd("my2");
    M2[3] = Getd("mz2"); 

    S1[1] = Getd("sx1");
    S1[2] = Getd("sy1");
    S1[3] = Getd("sz1"); 
    S2[1] = Getd("sx2");
    S2[2] = Getd("sy2");
    S2[3] = Getd("sz2");

    /*We compute the puncture parameter with the existing puncture project*/
            set_BH(level,Getd("mass2"),
               Getd("px2"),Getd("py2"),Getd("pz2"),
               -Getd("mx2"),-Getd("my2"),-Getd("mz2"),
               Getd("sx2"),Getd("sy2"),Getd("sz2")); 

  /*set the parameter again to the correct values */

   Setd("mass1", MASS1);
   Setd("mass2", MASS2);
   Setd("px1", P1[1]);
   Setd("px2", P2[1]);
   Setd("py1", P1[2]);
   Setd("py2", P2[2]);
   Setd("pz1", P1[3]);
   Setd("pz2", P2[3]);
   Setd("mx1", M1[1]);
   Setd("mx2", M2[1]);
   Setd("my1", M1[2]);
   Setd("my2", M2[2]);
   Setd("mz1", M1[3]);
   Setd("mz2", M2[3]);
   Setd("sx1", S1[1]);
   Setd("sx2", S2[1]);
   Setd("sy1", S1[2]);
   Setd("sy2", S2[2]);
   Setd("sz1", S1[3]);
   Setd("sz2", S2[3]); 

           
     
    } else {  
          printf("set BH 2 on level %d\n  ",level->l); 
           set_PUNC(level,Getd("mass2"), 
             Getd("px2"),Getd("py2"),Getd("pz2"),
             Getd("mx2"),Getd("my2"),Getd("mz2")); 
       boost_spacetime(level, Getd("mx2"),Getd("my2"),Getd("mz2"));}
    }
       
      sub++;
    
  }
 
  
  /* add both ******************************************************************/
  add_spacetime(level, sub );
  
  
  
  
  // disable vars
  //vlfree(vl);

  if (level->l == level->grid->lmax) {
  int lev;  
  int lmin = g->lmin;
  int lmax = g->lmax;
    set_TOV(level, 1,0,0,0,0,0,0, 0,0,0);
    
  /* if(Getv("BBID_modus","BHNS")) {
      for (lev = lmin; lev <= lmax; lev++) {
      level = g->level[lev];
      disablevar(level, Ind("incl_punctures_u"));
      } 
    }*/
  }
   return 0;
}
