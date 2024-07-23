/* AHmod.c */
/* Jose Gonzalez 03/07 (original AHF.c) */
/* Norbert Lages, 02/09 */
/*  more modularized, */
/*  separate calculation of all MOTSs possible */
/* Ivan Markin 05/24 (OpenMP parallelizaion)*/

#include "bam.h"
#include "AHmod.h"
#include "AHmod_output.h"

#define AHmod_MaxSurfaces 30
static int AHmod_found[AHmod_MaxSurfaces];


// uses global constant dequaleps and function dequal
// isn't there something similar already in BAM?
int AHmod_IsMultipleTime(const double time, const double steptime)
{
  const int dtime = (time+dequaleps)/steptime; // note implicit conversion
  return dequal(time-dtime*steptime, 0);
}

// isn't there something similar in BAM already?
int AHmod_FileExists(const char * filename)
{
  FILE *file;
  if ( (file = fopen(filename, "r")) )
  {
    fclose(file);
    return 1;
  }
  return 0;
}

int AHmod_ComputeSlice(const double currenttime, const double timeInterval)
{
  return (int)(currenttime/timeInterval);
//  const int slice = (int)(time/AHmod_time);
}


// creates a list of the factorial up to number maxn
// calling procedure must allocate memory
// use double because of larger numbers and 
//            prevent computing 1/fac[n] as integer
// inaccurate for around maxn>17,  mantissa can not hold exact value
void AHmod_CalculateFactorialList(const int maxn, double *fac)
{
  if (maxn<1)
    errorexit("AHmod_CalculateFactorial: need maxn>=1");

  int i;
  fac[0]=1;
  for (i=1; i<=maxn; i++)
    fac[i]=fac[i-1]*i;
}

// have only one place to compute dtheta and dphi
inline 
double AHmod_dtheta(const int ntheta)
{
  return PI/ntheta;
}

inline 
double AHmod_dphi(const int nphi)
{
  return 2.0*PI/nphi;
}

// have only one place to compute actual theta and phi values,
// can be useful if one changes between normal and staggered grid
inline 
double AHmod_theta(const double dtheta, const int i)
{
  return dtheta*(0.5 + i);
}

inline 
double AHmod_phi(const double dphi, const int j)
{
  return dphi*(0.5 + j);
}

void AHmod_CheckAndCreateResultVariables(const int nhorizons) 
{
  int i;
  char parname[100];
  
  for (i=0; i<nhorizons; i++) {
    sprintf(parname,"ahf%d_x",i);
    if (ExistPar(parname)) continue;
    
    // otherwise create all ahf variables with this number
    AddPar(parname, "0.0","approximate mid point x");
    sprintf(parname,"ahf%d_y",i);
    AddPar(parname, "0.0","approximate mid point y");
    sprintf(parname,"ahf%d_z",i);
    AddPar(parname, "0.0","approximate mid point z");  
    sprintf(parname,"ahf%d_m",i);
    AddPar(parname, "0.0","mass calculated from apparent horizon");
    sprintf(parname,"ahf%d_r",i);
    AddPar(parname, "-1.0","mean radius of apparent horizon");
  }
}


void AHmod_CheckAndCreateSurfaceVariables(const int nhorizons) 
{
  int i;
  char parname[50];
  
  for (i=0; i<nhorizons; i++) {
    sprintf(parname,"AHmod_surface%d_centerx",i);
    if (ExistPar(parname)) continue;
    
    // otherwise create all surface variables with this number
    AddPar(parname, "0.0","center x of star shaped surface");
    sprintf(parname,"AHmod_surface%d_centery",i);
    AddPar(parname, "0.0","center y of star shaped surface");
    sprintf(parname,"AHmod_surface%d_centerz",i);
    AddPar(parname, "0.0","center z of star shaped surface");
  }
  
  for (i=0; i<nhorizons; i++) {
    sprintf(parname,"AHmod_surface%d_WaitUntilClosePunctures",i);
    if (ExistPar(parname)) continue;
    
    AddPar(parname, "no","Don't search for common horizon when punctures are far away");
  }
  
}

// Compute Legendre P[l][m] for l>=m and its first and 
// second derivatives (in theta) for given fixed theta.
// The calling procedure has to allocate memory for P, dPdt, dPdtdt.
void AHmod_ComputeLegendre(const double theta, const int LMAX1,
          double **P, double **dPdtheta, double **dPdthetadtheta)
{
  int l,m;
  const double costheta = cos(theta);
  const double sintheta = sin(theta);
    
  const int sizefactorial = 2*LMAX1+2;
  double *fac = dvector(0, sizefactorial);  

  AHmod_CalculateFactorialList(2*LMAX1+1,fac);  // precompute factorial

  //initialize P, dPdtheta and dPdthtetadtheta to 0.0
  for(l=0;l<=LMAX1;l++){
    for(m=0;m<LMAX1;m++){
      P[l][m] = 0.0;
      dPdtheta[l][m] = 0.0;
      dPdthetadtheta[l][m] = 0.0;
    }     
  }

  //Compute the Legendre functions
  //diagonal terms
  for(l=0;l<=LMAX1;l++){
    P[l][l] = sqrt((2*l+1)*fac[2*l]/(4.0*PI))/(pow(2,l)*fac[l])*pow((-sintheta),l);
  }

  // compare AHF.c
  // the special treatment of P[1][0] etc. is not needed any more
  // the loop has special treatment for all P[l][l-1] where P[l-2][l-1] is not needed.
  P[1][0] = sqrt(3)*costheta*P[0][0];

  for(l=2;l<=LMAX1;l++){
    for(m=0;m<l-1;m++){
      P[l][m] = sqrt((2*l+1)/(l*l-m*m))*(sqrt(2*l-1)*costheta*P[l-1][m]
                         - sqrt(((l-1)*(l-1)-m*m)/(2*l-3))*P[l-2][m]);
    }
    // do [l][l-1] separately otherwise P[l-2][m] not defined
    P[l][l-1] = sqrt((2*l+1)/(l*l-m*m))*(sqrt(2*l-1)*costheta*P[l-1][m]);
  }



  //Compute first derivatives of the Legendre functions
  for(l=0;l<=LMAX1;l++){
    dPdtheta[l][l] = sqrt((2*l+1)*fac[2*l]/(4.0*PI))/(pow(2,l)*fac[l])*l*pow((-sintheta),l-1)*(-costheta);
  }
 
  dPdtheta[1][0] = sqrt(3)*(-sintheta*P[0][0]+costheta*dPdtheta[0][0]);

  for(l=2;l<=LMAX1;l++){
    for(m=0;m<l-1;m++){
      dPdtheta[l][m] = sqrt((2*l+1)/(l*l-m*m))*(sqrt(2*l-1)*(-sintheta*P[l-1][m]
                       + costheta*dPdtheta[l-1][m])-sqrt(((l-1)*(l-1)-m*m)/(2*l-3))*dPdtheta[l-2][m]);
    } 
    dPdtheta[l][l-1] = sqrt((2*l+1)/(l*l-m*m))*(sqrt(2*l-1)*(-sintheta*P[l-1][m]
                       + costheta*dPdtheta[l-1][m]) );
  }



  //Compute second derivatives of the Legendre functions
  for(l=0;l<=LMAX1;l++){
    dPdthetadtheta[l][l] = sqrt((2*l+1)*fac[2*l]/(4.0*PI))/(pow(2,l)*fac[l])*l
                         *((l-1)*pow(-sintheta,l-2)*costheta*costheta 
                         + pow(-sintheta,l-1)*sintheta);
  }

  dPdthetadtheta[1][0] = sqrt(3)*(-costheta*P[0][0]-2.0*sintheta*dPdtheta[0][0]
                       + costheta*dPdthetadtheta[0][0]);

  for(l=2;l<=LMAX1;l++){
    for(m=0;m<l-1;m++){
      dPdthetadtheta[l][m] = sqrt((2*l+1)/(l*l-m*m))*(sqrt(2*l-1)*(-cos(theta)*P[l-1][m]
                           - 2.0*sin(theta)*dPdtheta[l-1][m] + cos(theta)*dPdthetadtheta[l-1][m])
                           - sqrt(((l-1)*(l-1)-m*m)/(2*l-3))*dPdthetadtheta[l-2][m]);
    }
    // m will be l-1
    dPdthetadtheta[l][m] = sqrt((2*l+1)/(l*l-m*m))*(sqrt(2*l-1)*(-cos(theta)*P[l-1][m]
                           - 2.0*sin(theta)*dPdtheta[l-1][m] + cos(theta)*dPdthetadtheta[l-1][m]) );
  }
  
  free_dvector(fac,0,sizefactorial);
} // end Legendre




// Compute spherical harmonics for grid of size ntheta*nphi.
// Results can be used for all horizons.
void AHmod_CalculateSphericalHarmonics(const int ntheta, const int nphi,
       const int LMAX1,
       double **Y0, double **dY0dtheta,  double **dY0dthetadtheta,
       double **Yc,              double **Ys, 
       double **dYcdtheta,       double **dYsdtheta,
       double **dYcdphi,         double **dYsdphi, 
       double **dYcdthetadtheta, double **dYsdthetadtheta,
       double **dYcdthetadphi,   double **dYsdthetadphi,
       double **dYcdphidphi,     double **dYsdphidphi)
{
  const double dtheta=AHmod_dtheta(ntheta);
  const double dphi=AHmod_dphi(nphi);
  const double sqrt2=sqrt(2.0);
  double theta,phi;
  int i,j,l,m,i1,l1;
     
  // need Legendre polynomials
  //    P, dPdtheta, dPdthetadtheta are not needed later in AHmod.
  double **P,**dPdtheta,**dPdthetadtheta;
  P  = dmatrix(0,LMAX1,0,LMAX1);
  dPdtheta = dmatrix(0,LMAX1,0,LMAX1);
  dPdthetadtheta = dmatrix(0,LMAX1,0,LMAX1);  
  
  
  for(i=0;i<ntheta;i++){    
    theta = AHmod_theta(dtheta, i);
    
    // Compute Legendre polynomial, depends only on theta
    AHmod_ComputeLegendre (theta, LMAX1, P, dPdtheta, dPdthetadtheta);
    
    for(j=0;j<nphi;j++){	
      phi = AHmod_phi(dphi,j);
       
      i1=i*nphi+j; // grid point, combine i and j into one index
      
      //Compute the spherical harmonics

      for(l=0;l<=LMAX1;l++){
        Y0[i1][l] = P[l][0];
      }

      for(l=1;l<=LMAX1;l++){
        for(m=1;m<=l;m++){
          l1=l*(LMAX1+1)+m;
          Yc[i1][l1] = sqrt2*P[l][m]*cos(m*phi);
          Ys[i1][l1] = sqrt2*P[l][m]*sin(m*phi);
        }
      }

      //Compute first derivatives of the spherical harmonics

      for(l=0;l<=LMAX1;l++){
        dY0dtheta[i1][l] = dPdtheta[l][0];
      }

      for(l=1;l<=LMAX1;l++){
        for(m=1;m<=l;m++){
          l1=l*(LMAX1+1)+m;
          dYcdtheta[i1][l1] = sqrt2*dPdtheta[l][m]*cos(m*phi);
          dYsdtheta[i1][l1] = sqrt2*dPdtheta[l][m]*sin(m*phi);
          dYcdphi[i1][l1] = -sqrt2*P[l][m]*m*sin(m*phi);
          dYsdphi[i1][l1] =  sqrt2*P[l][m]*m*cos(m*phi);
        }
      }

      //Compute second derivatives of the spherical harmonics

      for(l=0;l<=LMAX1;l++){
        dY0dthetadtheta[i1][l] = dPdthetadtheta[l][0];
      }
        
      for(l=1;l<=LMAX1;l++){
        for(m=1;m<=l;m++){
          l1=l*(LMAX1+1)+m;
          dYcdthetadtheta[i1][l1] = sqrt2*dPdthetadtheta[l][m]*cos(m*phi);
          dYcdthetadphi[i1][l1] = -sqrt2*dPdtheta[l][m]*m*sin(m*phi);
          dYsdthetadtheta[i1][l1] = sqrt2*dPdthetadtheta[l][m]*sin(m*phi);
          dYsdthetadphi[i1][l1] = sqrt2*dPdtheta[l][m]*m*cos(m*phi);
          dYcdphidphi[i1][l1] = -sqrt2*P[l][m]*m*m*cos(m*phi);
          dYsdphidphi[i1][l1] = -sqrt2*P[l][m]*m*m*sin(m*phi);
        }
      }
    } // end phi loop
  } // end theta loop

  free_dmatrix(P,0,LMAX1,0,LMAX1);
  free_dmatrix(dPdtheta,0,LMAX1,0,LMAX1);
  free_dmatrix(dPdthetadtheta,0,LMAX1,0,LMAX1);  
} // end spherical harmonics


// extract the information about which surfaces to search for
// from AHmod_searchMTS
// the data is arranged as follows
// n1 ts1 te1 num11 num12 ... num1n1   n2 ts2 te2 num21 num22 ... num2n2  ...
// n_:  number of punctures that this surfaces should surround
// ts_: time to search this MTS (start)
// te_: time to search this MTS (end)
// num_i: the (id-)number of each puncture that should be surrounded by that surface
// Memory is allocated for ListPunctures
void AHmod_GetSearchData(struct AHmod_surface_struct *surfaces,
                         const int nhorizons)
{
  int i,j;
  if (strlen(Gets("AHmod_searchMTS"))>299)
    printf("Problem in AHmod_GetSearchData: String too long.\n");

  if (!strcmp(Gets("AHmod_searchMTS"), "0"))
    errorexit("Missing parameter AHmod_searchMTS.");

  char ahsearchMTS[300],*token;

  strcpy(ahsearchMTS,Gets("AHmod_searchMTS"));

  token=strtok(ahsearchMTS, " ");  // get first token out of ahsearchMTS
  
  for (i=0; i<nhorizons; i++) {
   
    if (token == 0)
      errorexit("Too few parameters in AHmod_searchMTS.");
 
    surfaces[i].NumOfPunctures = atoi(token);
    
    // next tokens are start and end time
    token=strtok(NULL, " ");
    surfaces[i].SearchTimeStart = atof(token);
    token=strtok(NULL, " "); 
    surfaces[i].SearchTimeEnd = atof(token);
    
    // allocate space to hold the enclosed punctures for each surface
    if ((surfaces[i].ListPunctures=(int *) malloc(surfaces[i].NumOfPunctures*sizeof(int)))==NULL)
      printf("in AHmod_getSearchData: wrong malloc\n");
      
    // next tokens are the punctures
    for (j=0;j<surfaces[i].NumOfPunctures;j++) {
      token=strtok(NULL, " ");
      surfaces[i].ListPunctures[j]= atoi(token);
    }    
    // get next token
    token=strtok(NULL, " ");
  }
}   

void AHmod_PrintSurfacesData(const int nhorizons, struct AHmod_surface_struct *surfaces)
{
  int i,j;
  
  printf ("\nAHmod_surfaces: ");
  for (i=0; i<nhorizons; i++) {
    printf("\n  Start: %6.1f  End: %6.1f ",surfaces[i].SearchTimeStart,
                                           surfaces[i].SearchTimeEnd);
    printf (" NumOfPunctures: %d  namely ",surfaces[i].NumOfPunctures);
    for (j=0; j<surfaces[i].NumOfPunctures; j++)
      printf("  %d",surfaces[i].ListPunctures[j]);
    }
  printf("\n");
}



void AHmod_ResetCoefficientsZero(double *a0, double **ac, double **as, 
                                 const int LMAX1) {
  int l,m;
  for(l=0;l<=LMAX1;l++){
    a0[l] = 0.0; 
    for(m=0;m<=l;m++){
      ac[l][m] = 0.0; 
      as[l][m] = 0.0; 
     }
  }
}

// have only one place to create the file name
void AHmod_GetFilenameLastSurface(const int number, const char *outdir, 
                                  char *filename)
{
  sprintf(filename, "%s/AHmod_surface%d",outdir,number);
} 

// no error checking for wrong puncture number!
void AHmod_GetPunctureLocation(const int number, double *xp, double *yp, double *zp)
{
  char puncStr[30];
  
  sprintf(puncStr, "moving_puncture%d_x",number);
  *xp = Getd(puncStr);
 
  sprintf(puncStr, "moving_puncture%d_y",number);
  *yp = Getd(puncStr);
  
  sprintf(puncStr, "moving_puncture%d_z",number);
  *zp = Getd(puncStr);
}

// note that numbering starts at 1
double AHmod_GetPunctureMass(const int number) 
{
  char puncStr[20];
  sprintf(puncStr, "mass%d",number);
  return Getd(puncStr);
}

// note that puncture numbering starts at 1, automatically add 1
// try to read mass of puncture
int AHmod_ExistPuncture(const int number) 
{
  char parname[30];
  sprintf(parname,"mass%d",number+1);
  return ExistPar(parname);
}


double AHmod_GetLargestDistance(const int numpunc, int *punclist) 
{
  int p1,p2;
  double x1,y1,z1,x2,y2,z2;
  double maxdist=0.0;

  for (p1=0; p1<numpunc-1; p1++) 
  {
    AHmod_GetPunctureLocation(punclist[p1]+1, &x1,&y1,&z1);
    for (p2=p1+1; p2<numpunc; p2++) {
      AHmod_GetPunctureLocation(punclist[p2]+1, &x2,&y2,&z2);
      maxdist = DMAX(maxdist, sqrt(pow(x2-x1,2)+ pow(y2-y1,2)+ pow(z2-z1,2)));
    }
  }
  return maxdist;
}


// add up masses of the given punctures
// +1 because of naming of punctures
double AHmod_sumPunctureMasses(const int numpunc, int *punclist)
{
  int j;
  double mass=0;
  for (j=0; j<numpunc; j++)
    mass += AHmod_GetPunctureMass(punclist[j]+1);

  return mass;
}


// If a previous file exists 
//   it is read and the initial surface is set to these values.
// Otherwise get an estimate by 
//   either taking some standard guess
//   or looking at the masses and taking an appropriate sphere

void AHmod_GetInitialGuess(const int number, const char *outdir,
           double *a0, double **ac, double **as, const int LMAX1,
           struct AHmod_surface_struct *surfaces) 
{
  char filename[300];
  AHmod_GetFilenameLastSurface(number, outdir, filename);

  AHmod_ResetCoefficientsZero(a0,ac,as,LMAX1);
    
  if (AHmod_FileExists(filename)) { 
    AHmod_GetInitialGuessFromFile(filename, a0, ac, as, LMAX1);
  
    // shrink or expand the initial guess a bit
    a0[0]*=Getd("AHmod_initial_guess_expand");
  }
  else
  {
    if (surfaces[number].NumOfPunctures==0)
      // when no masses from punctures are available use some guess as radius
      a0[0] = sqrt(4.0*PI) * Getd("AHmod_initial_radius");
    else
      AHmod_GetInitialGuessFromPunctures(a0, surfaces[number].NumOfPunctures,
                                             surfaces[number].ListPunctures);
  }
}


void AHmod_GetInitialGuessFromPunctures(double *a0, 
                                        const int numpunc, int *punclist)
{                                        
  double mass = AHmod_sumPunctureMasses(numpunc,punclist);
  double largedist = AHmod_GetLargestDistance(numpunc,punclist);
  
  // For single BH in isotropic coordinates: horizon radius=m/2
  // but make sure it can surround all punctures comfortably, i. e.
  // make radius a bit larger than half the distance between any of the punctures
  a0[0] = DMAX(0.5*mass, 1.4 * 0.5*largedist); 

  // sqrt(4 pi) is normalization    
  a0[0] *= sqrt(4.0*PI);
}


// read data from previous surface
// no check is done whether the data is corrupted
void AHmod_GetInitialGuessFromFile(char *filename,
           double *a0, double **ac, double **as, const int LMAX1) 
{
  int l,m;
  FILE *fahmod;
  fahmod=fopen(filename, "r");
  if (fahmod==NULL) {
    printf("Problem opening file in AHmod_GetInitialGuessFromFile.\n");
    return;
  }

  // Get initial guess from file
  for(l=0;l<=LMAX1;l++){
    fscanf(fahmod,"%le", &a0[l]);
    for(m=0;m<=l;m++){
      fscanf(fahmod,"%le",&ac[l][m]);
      fscanf(fahmod,"%le",&as[l][m]);
    }
  }
  fclose(fahmod);
}



// write data of the horizon to have a good guess in the next iteration
void AHmod_WriteSurfaceShape(const int number, const char *outdir,
           double *a0, double **ac, double **as, const int LMAX1) 
{
  int l,m;
  FILE *fahmod;
  char filename[300];
  AHmod_GetFilenameLastSurface(number, outdir, filename);
  fahmod=fopen(filename, "w");
  
  // Write surface to file
  for(l=0;l<=LMAX1;l++){
    fprintf(fahmod,"%e\n", a0[l]); 
    for(m=0;m<=l;m++){
      fprintf(fahmod,"%e\n",ac[l][m]);
      fprintf(fahmod,"%e\n",as[l][m]);
    }
  }
  fclose(fahmod);
}

// Remove file with initial guess from last search
// (e.g. when the last search failed)
void AHmod_RemoveLastSurface(const int number, const char *outdir) 
{
  char filename[300];
  AHmod_GetFilenameLastSurface(number, outdir, filename);

  if( !AHmod_FileExists(filename) ) return;
  
  if( remove(filename) != 0 )
    printf( "Error deleting file in AHmod_RemoveLastSurface.\n" );
}

void AHmod_GetParameterCentralPoint(const int number, double *xc, double *yc, double *zc)
{
  char parname[40];

  sprintf(parname,"AHmod_surface%d_centerx",number);  *xc=Getd(parname);
  sprintf(parname,"AHmod_surface%d_centery",number);  *yc=Getd(parname);
  sprintf(parname,"AHmod_surface%d_centerz",number);  *zc=Getd(parname);
}

// take weighted mean of puncture locations
void AHmod_GetWeightedMassCentralPoint(const int numpunc, int* punclist,
                                       double *xc, double *yc, double *zc) 
{  
  double x1,y1,z1,m1; // location and mass of one BH
  double summx=0.0;   // sum of m_i*x_i
  double summy=0.0;
  double summz=0.0;
  double summ=0.0;    // sum of m_i
  int i;
  
  for (i=0; i<numpunc; i++) {
    // add 1 because puncture name numbering starts at 1
    AHmod_GetPunctureLocation(punclist[i]+1, &x1, &y1, &z1);
    m1 = AHmod_GetPunctureMass(punclist[i]+1) ;
    summ  += m1;
    summx += m1*x1;
    summy += m1*y1;
    summz += m1*z1;
  }
  *xc = summx/summ;
  *yc = summy/summ;
  *zc = summz/summ;
}


// Write results in log file outdir/horizon_number
// one time step == one line
void  AHmod_WriteSurfaceResults(const int number, const char *outdir, 
   const int verbose, const double time, 
   const double xc, const double yc, const double zc, const double mass,
   const double Sx, const double Sy, const double Sz, const double S,
   const double coarea, const int levelnum)
{
  char filename[300];
  FILE *fp;
  
  sprintf(filename,"%s/horizon_%d",outdir, number);
  if (verbose) printf("  Horizon log at %s\n",filename);

  // if file does not exist, start new file and write header line
  if (!AHmod_FileExists(filename)) {
    fp = fopen(filename, "wb");
    if (!fp) errorexits("failed opening %s", filename);
    fprintf(fp, "#    time              x               y               z              mass             Sx              Sy              Sz             S           coord. area   lq\n");
  }
  else {            // if it does exist, reopen for append
    fp = fopen(filename, "ab");
    if (!fp) errorexits("failed opening %s", filename);
  }
  
  // write data
  fprintf(fp, "%14.6e  %14.6e  %14.6e  %14.6e  %14.6e  %14.6e  %14.6e  %14.6e  %14.6e  %14.6e  %d\n",
              time,xc,yc,zc,mass,Sx,Sy,Sz,S,coarea,levelnum);
  fclose(fp);
}

// fills an array of 6 points (18 values) that are the extent of
// a sphere around x0,y0,z0 with radius r
void AHmod_SphereEndPoints(double* pts, 
           const double x0, const double y0, const double z0, const double r)
{
  pts[0] = x0 - r;  pts[1] = y0;     pts[2] = z0;
  pts[3] = x0 + r;  pts[4] = y0;     pts[5] = z0;
  pts[6] = x0;      pts[7] = y0 - r; pts[8] = z0;
  pts[9] = x0;      pts[10]= y0 + r; pts[11]= z0;
  pts[12]= x0;      pts[13]= y0;     pts[14]= z0 - r;
  pts[15]= x0;      pts[16]= y0;     pts[17]= z0 + r;
}


//  in interpolate/interpolate.c
void AHmod_makeSymmetryPoint(tL *level, double *x, double *y, double *z)
{
  int fx,fy,fz;  
  map_xyz_withsym(level, x, y, z, &fx, &fy, &fz);
}


int AHmod_sphere_inside_level(tL *level, 
             const double xc, const double yc, const double zc, const double r)
{
  int j,k;  
  double pts[18];
  int inside;
  
  // fill with maximal and minimal points
  AHmod_SphereEndPoints(pts, xc, yc, zc, r);
  
  // use symmetry on each point (if necessary)
  for (k = 0; k < 6; k++) {
    AHmod_makeSymmetryPoint(level, &pts[k*3+0], &pts[k*3+1], &pts[k*3+2]);
  }
   
  // test each box separately
  for (j=0; j<level->nboxes; j++) {
    
    // every point has to be inside
    inside=1;
    for (k = 0; k<6; k++) {
      if ( !xyzinsidebbox(level->box[j]->bbox, pts[k*3+0], pts[k*3+1], pts[k*3+2]) ) {
        inside=0;
        break;
      }
    }
    if (inside==1) return 1;
    // else continue with other boxes...
  }
  
  return 0;
}


// checks whether a sphere with given mid point and radius fits
// into the given level, but not in the one below
int AHmod_IsOptimalLevel(tL* level, const double xc, const double yc, 
                         const double zc, const double meanradius)
{  
  const double r_max=Getd("AHmod_box_savety_factor")*meanradius;

  int incurrentlevel= AHmod_sphere_inside_level(level, xc, yc, zc, r_max);
  
  int inchildlevel=0;
  // if not already in max level, check the one below
  if (incurrentlevel && (level->l < level->grid->lmax) ) {
    tL *childlevel=level->grid->level[level->l+1];
    inchildlevel=AHmod_sphere_inside_level(childlevel, xc, yc, zc, r_max);
  }
  
  // output for testing
  // printf("  IsOptimalLevel: level=%d incurrentlevel=%d, inchildlevel=%d\n",
  //          level->l, incurrentlevel, inchildlevel);
  
  return (incurrentlevel && !(inchildlevel) );
}



//Compute the radius of the surface r=a_lm Y_lm
//loop over all surface points
void GetRadiiFromSphericalHarmonics(double **rr, const double dtheta,
  const double dphi, const int ntheta, const int nphi, const int LMAX1,
  double *a0, double **ac, double **as, double **Y0, double **Yc, double **Ys)
{
  double theta,phi;
  int i1,i,j,l,m;
  
  for(i=0;i<ntheta;i++){
    theta = AHmod_theta(dtheta, i);
    for(j=0;j<nphi;j++){      
      phi = AHmod_phi(dphi, j);
      
      i1=i*nphi+j;  // label points
      
      for(l=0;l<=LMAX1;l++){
        rr[i][j] += a0[l]*Y0[i1][l];
      }
      for(l=1;l<=LMAX1;l++){
        for(m=1;m<=l;m++){
          rr[i][j] += Yc[i1][l*(LMAX1+1)+m]*ac[l][m] + Ys[i1][l*(LMAX1+1)+m]*as[l][m];
        }
      }
      
    }
  }
}


//Compute derivatives of (x,y,z) with respect to theta and phi
// loop over all points of surface
// use 2nd order finite differences
void ComputeDerivativesXYZvsTP(double **rr, const double dtheta, const double dphi,
  const int ntheta, const int nphi, double *dxdt, double *dydt, double *dzdt,
  double *dxdp, double *dydp, double *dzdp)
{
  int i,j,i1;
  const int nphihalf = (int)(nphi/2);
  double theta,phi,rtp1,rtm1,rpp1,rpm1,drdt,drdp;
  double sintheta,costheta,sinphi,cosphi;
  
  for(i=0;i<ntheta;i++){
    theta = AHmod_theta(dtheta,i);
    sintheta=sin(theta);
    costheta=cos(theta);
    
    for(j=0;j<nphi;j++){
      phi = AHmod_phi(dphi,j);
      sinphi= sin(phi);
      cosphi= cos(phi);

      i1=i*nphi+j;

      // first decide where the neighbors are
      // read pp1 as phi plus 1
      //      tm1 as theta minus 1
      if(i==0){
        rtp1 = rr[1][j];
        if(j<nphihalf)  rtm1 = rr[0][j+nphihalf];
        if(j>=nphihalf) rtm1 = rr[0][j-nphihalf];
      } else if (i==ntheta-1) {
        if(j<nphihalf)  rtp1 = rr[ntheta-1][j+nphihalf];
        if(j>=nphihalf) rtp1 = rr[ntheta-1][j-nphihalf];
        rtm1 = rr[ntheta-2][j];
      } else {                    //normal case
        rtp1 = rr[i+1][j];
        rtm1 = rr[i-1][j];
      }
  
      if(j==0){
        rpp1 = rr[i][1];
        rpm1 = rr[i][nphi-1];
      }else if(j==nphi-1){
        rpp1 = rr[i][0];
        rpm1 = rr[i][nphi-2];
      }else{                    //normal case
        rpp1 = rr[i][j+1];
        rpm1 = rr[i][j-1];
      }
            
      drdt = (rtp1-rtm1)/(2.0*dtheta);
      drdp = (rpp1-rpm1)/(2.0*dphi);

      //Derivatives of (x,y,z) with respect to theta
      dxdt[i1] = (drdt*sintheta + rr[i][j]*costheta)*cosphi;
      dydt[i1] = (drdt*sintheta + rr[i][j]*costheta)*sinphi;
      dzdt[i1] =  drdt*costheta - rr[i][j]*sintheta;
            
      //Derivatives of (x,y,z) with respect to phi
      dxdp[i1] = (drdp*cosphi - rr[i][j]*sinphi)*sintheta;
      dydp[i1] = (drdp*sinphi + rr[i][j]*cosphi)*sintheta;
      dzdp[i1] =  drdp*costheta;
    }
  }
}


// Find new spectral components
// use fast flow
void ComputeFlowSpectralComponents(double *a0, double **ac, double **as, 
         const int LMAX1, const int ntheta, const int nphi, const int rank,
         const double dtheta, const double dphi, 
         double **Y0, double **Yc, double **Ys, double **rho, 
         double *globalp)
{            
  int i,j,l,m,i1,l1;
  double theta;
  const double alpha=1.0;
  const double beta=0.5;
  const double A = alpha/(LMAX1*(LMAX1+1))+beta;
  const double B = beta/alpha;

  double** int1_spec0 = dmatrix(0,ntheta,0,LMAX1);
  double*  spec0 = dvector(0,LMAX1);
  double*  spec0_l = dvector(0,LMAX1);
  double** int1_spec_c = dmatrix(0,ntheta,0,(LMAX1+1)*(LMAX1+1));
  double*  spec_c = dvector(0,(LMAX1+1)*(LMAX1+1));
  double*  spec_cl = dvector(0,(LMAX1+1)*(LMAX1+1));
  double** int1_spec_s = dmatrix(0,ntheta,0,(LMAX1+1)*(LMAX1+1));
  double*  spec_s = dvector(0,(LMAX1+1)*(LMAX1+1));
  double*  spec_sl = dvector(0,(LMAX1+1)*(LMAX1+1));

  // loop over all modes 0 .. LMAX1
  for(l=0;l<=LMAX1;l++){   
    
    // compute change of a0 component
    
    //sum all the results that this processor is responsible of
    for(i=0;i<ntheta;i++){
      int1_spec0[i][l]=0.0;
      for(j=0;j<nphi;j++){ 
        i1=i*nphi+j;
        if (globalp[i1] == rank){
          int1_spec0[i][l]+=dphi*Y0[i1][l]*rho[i][j];
        }          
      } 
    }

    spec0_l[l] = 0.0;
    for(i=0;i<ntheta;i++){
      theta = AHmod_theta(dtheta,i);
      spec0_l[l]+=dtheta*int1_spec0[i][l]*sin(theta);
    }

    //the integral is the sum over all local sums
    bampi_allreduce_sum_vector(spec0_l+l, spec0+l, 1);

    a0[l] -= A/(1.0+B*l*(l+1))*spec0[l];
   
    // printf("a0[%d] = %22.16e\n",l,a0[l]);

    // compute changes of ac and as components

    for(m=0;m<=l;m++){
    
      l1=l*(LMAX1+1)+m;
      
      for(i=0;i<ntheta;i++){ 
        int1_spec_c[i][l1]=0.0;
        int1_spec_s[i][l1]=0.0;
        for(j=0;j<nphi;j++){ 
          i1=i*nphi+j;
          if (globalp[i1] == rank){
            int1_spec_c[i][l1]+=dphi*Yc[i1][l1]*rho[i][j];
            int1_spec_s[i][l1]+=dphi*Ys[i1][l1]*rho[i][j];
          }
        } 
      }

      spec_cl[l1] = 0.0;
      spec_sl[l1] = 0.0;        
      for(i=0;i<ntheta;i++){
        theta = AHmod_theta(dtheta,i);
        spec_cl[l1]+=dtheta*int1_spec_c[i][l1]*sin(theta);
        spec_sl[l1]+=dtheta*int1_spec_s[i][l1]*sin(theta);
      }

      bampi_allreduce_sum_vector(spec_cl+l1, spec_c+l1, 1);
      bampi_allreduce_sum_vector(spec_sl+l1, spec_s+l1, 1);

      ac[l][m] -= A/(1.0+B*l*(l+1))*spec_c[l1];
      as[l][m] -= A/(1.0+B*l*(l+1))*spec_s[l1];
	
      // printf("ac[%d][%d] = %22.16e\n",l,m,ac[l][m]);
      // printf("as[%d][%d] = %22.16e\n",l,m,as[l][m]);
    }
  }
  
  free_dmatrix(int1_spec0,0,ntheta,0,LMAX1);
  free_dvector(spec0,0,LMAX1);
  free_dvector(spec0_l,0,LMAX1);
  free_dmatrix(int1_spec_c,0,ntheta,0,(LMAX1+1)*(LMAX1+1));
  free_dvector(spec_c,0,(LMAX1+1)*(LMAX1+1));
  free_dvector(spec_cl,0,(LMAX1+1)*(LMAX1+1));
  free_dmatrix(int1_spec_s,0,ntheta,0,(LMAX1+1)*(LMAX1+1));
  free_dvector(spec_s,0,(LMAX1+1)*(LMAX1+1));
  free_dvector(spec_sl,0,(LMAX1+1)*(LMAX1+1));
}

// Perform integrals over the surface
// Integrate the area element, the mean of the expansion and the spin components
void IntegrateOverSurface(double *v, const int ntheta, const int nphi, 
     const double dtheta, const double dphi, const int rank, double *globalp, double **deth,
     double **H, double **intSx, double **intSy, double **intSz, double **rr)                          
{
  int i,j,i1;
  double theta;
  double* int1_area = dvector(0,ntheta-1);
  double* int1_coarea = dvector(0,ntheta-1);
  double* int1_hrms = dvector(0,ntheta-1);
  double* int1_hmean = dvector(0,ntheta-1);
  double* integral_Sx = dvector(0,ntheta-1);
  double* integral_Sy = dvector(0,ntheta-1);
  double* integral_Sz = dvector(0,ntheta-1);
  double v_l[7];
  
  for(i=0;i<ntheta;i++){
    int1_area[i]=0.0;
    int1_coarea[i]=0.0;
    int1_hrms[i]=0.0;
    int1_hmean[i]=0.0;
    integral_Sx[i]=0.0;
    integral_Sy[i]=0.0;
    integral_Sz[i]=0.0;
  }

  // use only points in local processor.
  // The check that every point is in some processor is done before
  for(i=0;i<ntheta;i++){
    for(j=0;j<nphi;j++){
      i1=i*nphi+j; 
      if (globalp[i1] == rank){
        theta=AHmod_theta(dtheta,i);
        int1_area[i]   += dphi* sqrt(deth[i][j]);
        int1_coarea[i] += dphi* rr[i][j]*rr[i][j]*sin(theta);
        int1_hrms[i]   += dphi* H[i][j]*H[i][j]*sqrt(deth[i][j]);
        int1_hmean[i]  += dphi* H[i][j]*sqrt(deth[i][j]);
        integral_Sx[i] += dphi* intSx[i][j];
        integral_Sy[i] += dphi* intSy[i][j];
        integral_Sz[i] += dphi* intSz[i][j];
     }     
   } 
  }

  double area_l  = 0.0;
  double coarea_l= 0.0;
  double hrms_l  = 0.0;
  double hmean_l = 0.0;
  double Sx_l    = 0.0;
  double Sy_l    = 0.0;
  double Sz_l    = 0.0;
  
  for(i=0;i<ntheta;i++){
    area_l   += dtheta* int1_area[i];
    coarea_l += dtheta* int1_coarea[i];
    hrms_l   += dtheta* int1_hrms[i];
    hmean_l  += dtheta* int1_hmean[i];
    Sx_l     += dtheta* integral_Sx[i];
    Sy_l     += dtheta* integral_Sy[i];
    Sz_l     += dtheta* integral_Sz[i];
  }

  // the integral is the sum over all local sums
  // use vector to sum over all processors  
  v_l[0] = area_l; 
  v_l[1] = hmean_l; 
  v_l[2] = hrms_l; 
  v_l[3] = Sx_l; 
  v_l[4] = Sy_l; 
  v_l[5] = Sz_l;
  v_l[6] = coarea_l;
  
  bampi_allreduce_sum_vector(v_l, v, 7);
  // return value is v
  
  free_dvector(int1_area,0,ntheta-1);
  free_dvector(int1_coarea,0,ntheta-1);
  free_dvector(int1_hrms,0,ntheta-1);
  free_dvector(int1_hmean,0,ntheta-1);
  free_dvector(integral_Sx,0,ntheta-1);
  free_dvector(integral_Sy,0,ntheta-1);
  free_dvector(integral_Sz,0,ntheta-1);
}

// returns 1 when the maximal distance between all punctures is below threshold
// otherwise returns 0
int AHmod_ArePuncturesClose(const int numpunc, int *punclist)
{
  const double totalmass = AHmod_sumPunctureMasses(numpunc, punclist);
  const double maxdist = AHmod_GetLargestDistance(numpunc, punclist);
  const double mergedist = Getd("AHmod_merger_distance");
  
  return (maxdist < mergedist*totalmass) ? 1 : 0;  
}

int AHmod_WaitUntilClosePunctures(const int number)
{
  char parname[50];
  sprintf(parname, "AHmod_surface%d_WaitUntilClosePunctures", number);
  return Getv(parname, "yes") ? 1 : 0;
}

// Checks whether any surface uses the WaitUntilClosePunctures feature
int AHmod_WaitUntilAnyClosePunctures(const int numsurfaces) 
{  
  int i;
  
  for (i=0; i<numsurfaces; i++) {
    if ( AHmod_WaitUntilClosePunctures(i) ) return 1;
  }
  return 0;
}

void AHmod_set_global_values(const int N, const double t, const double x, const double y,
                             const double z, const double r, const double m)
{
  char str[20];

  sprintf(str,"ahf%d_x",N);   Setd(str, x);
  sprintf(str,"ahf%d_y",N);   Setd(str, y);
  sprintf(str,"ahf%d_z",N);   Setd(str, z);
  sprintf(str,"ahf%d_m",N);   Setd(str, m);
  sprintf(str,"ahf%d_r",N);   Setd(str, r);
}




int AHmod(tL *level) {

    timer_start(0, "AHmod");
    
    // do not search in shells -> would be a mess to change everything 
    if (level->shells)
      return 0;

    const int nhorizons = Geti("AHmod_nhorizons");
    const int ntheta    = Geti("AHmod_ntheta");
    const int nphi      = Geti("AHmod_nphi");	
    const int verbose   =!Getv("AHmod_verbose", "no");
    const int rank      = bampi_rank();

    const double AHmod_time = Getd("AHmod_time");
    const double time       = level->time;
    const char *outdir      = Gets("outdir");

           
    // check whether ahf_ variables exist; if not, create them
    AHmod_CheckAndCreateResultVariables(nhorizons);
    
    // check whether AHmod_surface...  variables exist; if not, create them
    AHmod_CheckAndCreateSurfaceVariables(nhorizons);
    
    //-----------------------------------------------------
    // check several ways to exit AHmod without calculation
    
    // check symmetry
    if (Getv("grid", "1D")) errorexit("AHmod: not implemented in symmetry: 1D");
    if (Getv("grid", "2D")) errorexit("AHmod: not implemented in symmetry: 2D");

    // check number of MTS that are searched for
    if (nhorizons<=0) {
      if(rank==0 && verbose) {
        printf("No surfaces to be searched for, AHmod_nhorizons= %i \n",nhorizons);
      }  
      return(0);
    }

    // check time is multiple of AHmod_time
    if (!AHmod_IsMultipleTime(time, AHmod_time)) {
      return(0);
    }
    
    // check nphi to be even, needed for boundaries at the poles
    if ((nphi+1)%2==0) errorexit("AHmod: odd number of points in phi direction");
    //if (ntheta%2==0) errorexit("AHmod: even number of circles");
    
    // merger feature needs punctures
    if ( !Getv("physics", "punctures") &&
         AHmod_WaitUntilAnyClosePunctures(nhorizons) ) {
      errorexit("AHmod_surface?_WaitUntilClosePuncture = yes needs punctures project!\n");
    }
    
    // the global array AHmod has fixed size
    if (nhorizons >= AHmod_MaxSurfaces) {
      printf("AHmod_MaxSurfaces is reached, nhorizons=%d\n",nhorizons);
      errorexit("Problems may occur for larger values in AHmod");
    }
    
    // end of checks
    // ---------------------------------------------------

    const int LMAX1      = Geti("AHmod_LMAX");
    const int flow_iter  = Geti("AHmod_flow_iter");

    const int ahintorder = Geti("AHmod_interpolation_order");
    // setting interpolation order to the value given by AH, restore in the end
    
    int number,l,m,l1,i1,cont,u0;

    double sum;
    double xc,yc,zc;
    double xp,yp,zp,rp,rhop;
    double Kxx,Kxy,Kxz,Kyy,Kyz,Kzz;
    double gxx,gxy,gxz,gyy,gyz,gzz;
    double dgxxdx,dgxydx,dgxzdx,dgyydx,dgyzdx,dgzzdx;
    double dgxxdy,dgxydy,dgxzdy,dgyydy,dgyzdy,dgzzdy;
    double dgxxdz,dgxydz,dgxzdz,dgyydz,dgyzdz,dgzzdz;
    double ginvdetxx,ginvdetxy,ginvdetxz,ginvdetyy,ginvdetyz,ginvdetzz;
    double ginvxx,ginvxy,ginvxz,ginvyx,ginvyy,ginvyz,ginvzx,ginvzy,ginvzz,detg;
    double drdx,drdy,drdz,dthetadx,dthetady,dthetadz,dphidx,dphidy,dphidz;
    double drdxdx,drdxdy,drdxdz,drdydy,drdydz,drdzdz;
    double dthetadxdx,dthetadxdy,dthetadxdz,dthetadydy,dthetadydz,dthetadzdz;
    double dphidxdx,dphidxdy,dphidxdz,dphidydy,dphidydz,dphidzdz;
    double dFdxdx,dFdxdy,dFdxdz,dFdydx,dFdydy,dFdydz,dFdzdx,dFdzdy,dFdzdz;
    double d2F,dFdadFdbKab,dFdadFdbFdadb,dFupdx,dFupdy,dFupdz,K;
    double nnFxx,nnFxy,nnFxz,nnFyy,nnFyz,nnFzz;
    double h11,h12,h22;
    double area,coarea,mass,mass_old,meanradius;
    double Rx,Ry,Rz,Sx,Sy,Sz,S;
    double phix_x,phix_y,phix_z;
    double phiy_x,phiy_y,phiy_z;
    double phiz_x,phiz_y,phiz_z;
    double hrms,hmean,v[7];
    
    double *a0 =  dvector(0,LMAX1);
    double **ac = dmatrix(0,LMAX1,0,LMAX1);
    double **as = dmatrix(0,LMAX1,0,LMAX1);
    double **rr = dmatrix(0,ntheta-1,0,nphi-1);
    double **Y0 = dmatrix(0,ntheta*nphi,0,LMAX1);
    double **Yc = dmatrix(0,ntheta*nphi,0,(LMAX1+1)*(LMAX1+1));
    double **Ys = dmatrix(0,ntheta*nphi,0,(LMAX1+1)*(LMAX1+1)); 
    double **dY0dtheta = dmatrix(0,ntheta*nphi,0,LMAX1);
    double **dYcdtheta = dmatrix(0,ntheta*nphi,0,(LMAX1+1)*(LMAX1+1));
    double **dYsdtheta = dmatrix(0,ntheta*nphi,0,(LMAX1+1)*(LMAX1+1));
    double **dYcdphi = dmatrix(0,ntheta*nphi,0,(LMAX1+1)*(LMAX1+1));
    double **dYsdphi = dmatrix(0,ntheta*nphi,0,(LMAX1+1)*(LMAX1+1));

    double **u=dmatrix(0,ntheta-1,0,nphi-1);
    double **H = dmatrix(0,ntheta-1,0,nphi-1);
    double **rho = dmatrix(0,ntheta-1,0,nphi-1);
    double **dFdx = dmatrix(0,ntheta-1,0,nphi-1);
    double **dFdy = dmatrix(0,ntheta-1,0,nphi-1);
    double **dFdz = dmatrix(0,ntheta-1,0,nphi-1);
    double **sigma = dmatrix(0,ntheta-1,0,nphi-1);
    double **deth = dmatrix(0,ntheta-1,0,nphi-1);
    double *localp = dvector(0,ntheta*nphi);
    double *globalp = dvector(0,ntheta*nphi);

    double **dY0dthetadtheta = dmatrix(0,ntheta*nphi,0,LMAX1);
    double **dYcdthetadtheta = dmatrix(0,ntheta*nphi,0,(LMAX1+1)*(LMAX1+1));
    double **dYcdthetadphi = dmatrix(0,ntheta*nphi,0,(LMAX1+1)*(LMAX1+1));
    double **dYcdphidphi = dmatrix(0,ntheta*nphi,0,(LMAX1+1)*(LMAX1+1));
    double **dYsdthetadtheta = dmatrix(0,ntheta*nphi,0,(LMAX1+1)*(LMAX1+1));
    double **dYsdthetadphi = dmatrix(0,ntheta*nphi,0,(LMAX1+1)*(LMAX1+1));
    double **dYsdphidphi = dmatrix(0,ntheta*nphi,0,(LMAX1+1)*(LMAX1+1));
    
    double *dxdt = dvector(0,ntheta*nphi);
    double *dxdp = dvector(0,ntheta*nphi);
    double *dydt = dvector(0,ntheta*nphi);
    double *dydp = dvector(0,ntheta*nphi);
    double *dzdt = dvector(0,ntheta*nphi);
    double *dzdp = dvector(0,ntheta*nphi);
    double **intSx = dmatrix(0,ntheta-1,0,nphi-1);
    double **intSy = dmatrix(0,ntheta-1,0,nphi-1);
    double **intSz = dmatrix(0,ntheta-1,0,nphi-1);

    // need no allocation, will point to existing arrays
    //double *gxxg,*gxyg,*gxzg,*gyyg,*gyzg,*gzzg;
    double *dgxxdxg,*dgxydxg,*dgxzdxg,*dgyydxg,*dgyzdxg,*dgzzdxg;
    double *dgxxdyg,*dgxydyg,*dgxzdyg,*dgyydyg,*dgyzdyg,*dgzzdyg;
    double *dgxxdzg,*dgxydzg,*dgxzdzg,*dgyydzg,*dgyzdzg,*dgzzdzg;
    tVarList *vl,*wl;
    
    const double dtheta = AHmod_dtheta(ntheta);
    const double dphi   = AHmod_dphi(nphi);

    
    // Compute spherical harmonics only once, 
    // the spherical grid is the same for all surfaces

    AHmod_CalculateSphericalHarmonics(ntheta, nphi, LMAX1, 
           Y0, dY0dtheta,  dY0dthetadtheta, 
           Yc, Ys,  dYcdtheta, dYsdtheta,  dYcdphi,  dYsdphi, 
           dYcdthetadtheta, dYsdthetadtheta, dYcdthetadphi, dYsdthetadphi,   
           dYcdphidphi,  dYsdphidphi);


    // each time step starts at finest level
    // when in finest level erase information about previously found surfaces
    if (level->l == level->grid->lmax) {
      for (int i=0; i<nhorizons; i++)
        AHmod_found[i] = 0;
    }    
    
    // has information about surfaces
    struct AHmod_surface_struct *surfaces;   
    surfaces = (struct AHmod_surface_struct*) 
      malloc(nhorizons*sizeof(struct AHmod_surface_struct));
    AHmod_GetSearchData(surfaces, nhorizons);
   
    // testing
    //AHmod_PrintSurfacesData(nhorizons, surfaces);
    
    //Start loop over number of surfaces
    for(number=0;number<nhorizons;number++){

      // -------------------------------------------------------
      // Don't proceed with this surface in the following cases
      
      // time has to be in interval
      if ((time<surfaces[number].SearchTimeStart) || 
          (time>surfaces[number].SearchTimeEnd)) {
        continue;  // go to next surface
      }

      if (Geti("AHmod_start_level") >= 0) {
        if (level->l < Geti("AHmod_start_level")) {
          continue;
        }
      }
      // if this surface was already found in a finer level, don't search again
      if (Getv("AHmod_always_search","no")) {
        if ( AHmod_found[number] ) {
          continue;  // go to next surface
        }
      }
      
      // don't start too early for common horizons
      if (  AHmod_WaitUntilClosePunctures(number) && 
            !AHmod_ArePuncturesClose(surfaces[number].NumOfPunctures,
                                     surfaces[number].ListPunctures) ) {
        printf("\nNo search for common horizon yet because punctures are not close enough");
        continue;  // go to next surface
      }
      
      // end of possible early breaks
      // --------------------------------------------------------
      
      
      // Initial guess for surface
      AHmod_GetInitialGuess(number, outdir, a0, ac, as, LMAX1, surfaces);
      meanradius = a0[0]/sqrt(4.0*PI);
      
      // Get central point
      if (surfaces[number].NumOfPunctures == 0)  {  // no puncture, use parameter values
        AHmod_GetParameterCentralPoint(number, &xc, &yc, &zc);
      } else {  // get from weighted puncture location
        AHmod_GetWeightedMassCentralPoint(surfaces[number].NumOfPunctures,
                                          surfaces[number].ListPunctures, &xc, &yc, &zc);
      }
      
      
      // Old method, check whether level is optimal for this guess      
      // this is not the best way because this information is re-created at all levels 
      // This should not interfere with AHmod_found because it will
      // search in only one level anyway.
      if (Getv("AHmod_UseOptimalLevel","yes")) {
        if ( !AHmod_IsOptimalLevel(level, xc, yc, zc, meanradius) ) {
          continue;
        }
      }
      
      if(verbose && rank==0) {
        printf("\nSearching for horizon %d in level %d at time %f\n",number,level->l,time);
        printf("  centered in:  x = %f  y = %f  z = %f,  r_mean = %f\n",xc,yc,zc,meanradius);
      }

      //Compute the derivatives of the metric in the level
      //and store them in the variables dgxxdxg,....
      vl = vlalloc(level);
      vlpush(vl,Ind("AHmod_dgdxxx"));
      enablevarlist(vl);
      
      // give names to certain pointer locations
      double *gxxg = Ptr(level, "adm_gxx");
      double *gxyg = Ptr(level, "adm_gxy");
      double *gxzg = Ptr(level, "adm_gxz");
      double *gyyg = Ptr(level, "adm_gyy");
      double *gyzg = Ptr(level, "adm_gyz");
      double *gzzg = Ptr(level, "adm_gzz");
      dgxxdxg = vldataptr(vl,0);
      dgxydxg = vldataptr(vl,1);
      dgxzdxg = vldataptr(vl,2);
      dgyydxg = vldataptr(vl,3);
      dgyzdxg = vldataptr(vl,4);
      dgzzdxg = vldataptr(vl,5);
      dgxxdyg = vldataptr(vl,6);
      dgxydyg = vldataptr(vl,7);
      dgxzdyg = vldataptr(vl,8);
      dgyydyg = vldataptr(vl,9);
      dgyzdyg = vldataptr(vl,10);
      dgzzdyg = vldataptr(vl,11);
      dgxxdzg = vldataptr(vl,12);
      dgxydzg = vldataptr(vl,13);
      dgxzdzg = vldataptr(vl,14);
      dgyydzg = vldataptr(vl,15);
      dgyzdzg = vldataptr(vl,16);
      dgzzdzg = vldataptr(vl,17);

      const double idx = 1.0/(2.0*level->dx);
      const double idy = 1.0/(2.0*level->dy);
      const double idz = 1.0/(2.0*level->dz);

      forinnerpoints_ijk(level){
        dgxxdxg[ijk] = idx*(gxxg[ijk+di] - gxxg[ijk-di]);
        dgxydxg[ijk] = idx*(gxyg[ijk+di] - gxyg[ijk-di]);
        dgxzdxg[ijk] = idx*(gxzg[ijk+di] - gxzg[ijk-di]);
        dgyydxg[ijk] = idx*(gyyg[ijk+di] - gyyg[ijk-di]);
        dgyzdxg[ijk] = idx*(gyzg[ijk+di] - gyzg[ijk-di]);
        dgzzdxg[ijk] = idx*(gzzg[ijk+di] - gzzg[ijk-di]);
        dgxxdyg[ijk] = idy*(gxxg[ijk+dj] - gxxg[ijk-dj]);
        dgxydyg[ijk] = idy*(gxyg[ijk+dj] - gxyg[ijk-dj]);
        dgxzdyg[ijk] = idy*(gxzg[ijk+dj] - gxzg[ijk-dj]);
        dgyydyg[ijk] = idy*(gyyg[ijk+dj] - gyyg[ijk-dj]);
        dgyzdyg[ijk] = idy*(gyzg[ijk+dj] - gyzg[ijk-dj]);
        dgzzdyg[ijk] = idy*(gzzg[ijk+dj] - gzzg[ijk-dj]);
        dgxxdzg[ijk] = idz*(gxxg[ijk+dk] - gxxg[ijk-dk]);
        dgxydzg[ijk] = idz*(gxyg[ijk+dk] - gxyg[ijk-dk]);
        dgxzdzg[ijk] = idz*(gxzg[ijk+dk] - gxzg[ijk-dk]);
        dgyydzg[ijk] = idz*(gyyg[ijk+dk] - gyyg[ijk-dk]);
        dgyzdzg[ijk] = idz*(gyzg[ijk+dk] - gyzg[ijk-dk]);
        dgzzdzg[ijk] = idz*(gzzg[ijk+dk] - gzzg[ijk-dk]);
      } endfor_ijk;

      
      wl = vlalloc(level);
      vlpush(wl,Ind("adm_gxx"));
      vlpush(wl,Ind("adm_Kxx"));
      vlpushvl(wl,vl);

      // reserve room for interpolation data
      double **interpolated_variables = (double**) malloc( (size_t)( ntheta*nphi * sizeof(double*) ) );
      for (int direction=0; direction < ntheta*nphi; direction++) {
        interpolated_variables[direction] = (double*) malloc( (size_t)( wl->n * sizeof(double) )) ;
      }

      // set boundaries
      set_boundary_symmetry(level, vl);
   
      // synchronize
      bampi_vlsynchronize(vl);

            
      mass=0.0;  // have old mass stored
      
      // ----------------------   Flow loop  ------------------------------ 
      // ------------------------------------------------------------------
      for(int k=0;k<flow_iter;k++){

        if (verbose && (rank==0)) {
          printf("Iteration %d\n",k);
        }

        //Initialize
        for(int i=0;i<ntheta;i++){
          for(int j=0;j<nphi;j++){
            rr[i][j] = 0.0; 
            u[i][j]  = 0.0; 
            H[i][j]  = 0.0; 
            rho[i][j] = 0.0;
            deth[i][j] = 0.0;
          }
        }
        
        //Compute the radius(heights) of the surface r=a_lm Y_lm
        GetRadiiFromSphericalHarmonics(rr, dtheta, dphi, ntheta, nphi, LMAX1,
                                       a0, ac, as, Y0, Yc, Ys);
        
        // Compute derivatives of (x,y,z) with respect to theta and phi
        ComputeDerivativesXYZvsTP(rr, dtheta, dphi,ntheta, nphi, 
                                  dxdt, dydt, dzdt, dxdp, dydp, dzdp);
        
        timer_start(0, "AHmod_local");

        //loop over each point of the surface
        bampi_openmp_parallel_for_collapse2
        for(int i=0;i<ntheta;i++){
          for(int j=0;j<nphi;j++){
            double theta = AHmod_theta(dtheta,i);
            double phi = AHmod_phi(dphi,j);
           
            int point_index=i*nphi+j;
            
            //global coordinates of the surface
            double x = xc + rr[i][j]*sin(theta)*cos(phi);
            double y = yc + rr[i][j]*sin(theta)*sin(phi);
            double z = zc + rr[i][j]*cos(theta);
             
            //Interpolate the extrinsic curvature, the 3 metric and the
            //derivatives of the 3 metric over the surface
            //The interpolation is only done on the points in this processor

            int flag = check_interpolation_cube_local_withsym(level,x,y,z,ahintorder);
            if(!flag){
              localp[point_index] = -1;
              continue;  // don't calculate on this processor
            }
            localp[point_index] = rank;

            interpolate_xyz_local_minimal_withsym(level,x,y,z,wl->n,wl->index,interpolated_variables[point_index],ahintorder,LAGRANGE);
          }
        }

        //loop over each point of the surface
        for(int i=0;i<ntheta;i++){
          double theta = AHmod_theta(dtheta,i);
          for(int j=0;j<nphi;j++){
            double phi = AHmod_phi(dphi,j);

            int i1 = i*nphi+j;

            //global coordinates of the surface
            double x = xc + rr[i][j]*sin(theta)*cos(phi);
            double y = yc + rr[i][j]*sin(theta)*sin(phi);
            double z = zc + rr[i][j]*cos(theta);

            // Skip non-local points
            if (localp[i1] == -1) {
              continue;
            }

            // Use the data interpolated in the loop before
            double *vinterp = interpolated_variables[i1];
            int nv = 0;

            gxx = vinterp[nv++];  
            gxy = vinterp[nv++];  
            gxz = vinterp[nv++];  
            gyy = vinterp[nv++];  
            gyz = vinterp[nv++];  
            gzz = vinterp[nv++];  

            Kxx = vinterp[nv++];  
            Kxy = vinterp[nv++];  
            Kxz = vinterp[nv++];  
            Kyy = vinterp[nv++];  
            Kyz = vinterp[nv++];  
            Kzz = vinterp[nv++];  

            dgxxdx = vinterp[nv++];  
            dgxydx = vinterp[nv++];  
            dgxzdx = vinterp[nv++];  
            dgyydx = vinterp[nv++];  
            dgyzdx = vinterp[nv++];  
            dgzzdx = vinterp[nv++];  
            dgxxdy = vinterp[nv++];  
            dgxydy = vinterp[nv++];  
            dgxzdy = vinterp[nv++];  
            dgyydy = vinterp[nv++];  
            dgyzdy = vinterp[nv++];  
            dgzzdy = vinterp[nv++];  
            dgxxdz = vinterp[nv++];  
            dgxydz = vinterp[nv++];  
            dgxzdz = vinterp[nv++];  
            dgyydz = vinterp[nv++];  
            dgyzdz = vinterp[nv++];  
            dgzzdz = vinterp[nv++];  

            //Inverse 3-metric.  

            ginvdetxx = gyy*gzz - gyz*gyz;
            ginvdetxy = gxz*gyz - gxy*gzz;
            ginvdetxz = gxy*gyz - gxz*gyy;
            ginvdetyy = gxx*gzz - gxz*gxz;
            ginvdetyz = gxy*gxz - gxx*gyz;
            ginvdetzz = gxx*gyy - gxy*gxy;
            
            detg = gxx*ginvdetxx + gxy*ginvdetxy + gxz*ginvdetxz;
           
            ginvxx = ginvdetxx/detg;
            ginvxy = ginvdetxy/detg;
            ginvxz = ginvdetxz/detg;
            ginvyx = ginvdetxy/detg;
            ginvyy = ginvdetyy/detg;
            ginvyz = ginvdetyz/detg;
            ginvzx = ginvdetxz/detg;
            ginvzy = ginvdetyz/detg;
            ginvzz = ginvdetzz/detg;

            //Trace of K

            K = ginvxx*Kxx + ginvxy*Kxy + ginvxz*Kxz
              + ginvyx*Kxy + ginvyy*Kyy + ginvyz*Kyz
              + ginvzx*Kxz + ginvzy*Kyz + ginvzz*Kzz;
            
            //local coordinates of the surface

            xp   = x-xc;
            yp   = y-yc;
            zp   = z-zc;
            rp   = sqrt(xp*xp + yp*yp + zp*zp);
            rhop = sqrt(xp*xp + yp*yp);
            
            //first derivatives of (r,theta,phi) with respect to (x,y,z)

            drdx = xp/rp;
            drdy = yp/rp;
            drdz = zp/rp;

            dthetadx = zp*xp/(rp*rp*rhop);
            dthetady = zp*yp/(rp*rp*rhop);
            dthetadz = -rhop/(rp*rp);
            
            dphidx = -yp/(rhop*rhop);
            dphidy = xp/(rhop*rhop);
            dphidz = 0.0;

            //second derivatives of (r,theta,phi) with respect to (x,y,z)

            drdxdx = 1.0/rp - xp*xp/(rp*rp*rp);
            drdxdy = - xp*yp/(rp*rp*rp);
            drdxdz = - xp*zp/(rp*rp*rp);
            drdydy = 1.0/rp - yp*yp/(rp*rp*rp);
            drdydz = - yp*zp/(rp*rp*rp);
            drdzdz = 1.0/rp - zp*zp/(rp*rp*rp);
          
            dthetadxdx = zp*(-2.0*xp*xp*xp*xp-xp*xp*yp*yp+yp*yp*yp*yp+zp*zp*yp*yp)/(rp*rp*rp*rp*rhop*rhop*rhop);
            dthetadxdy = - xp*yp*zp*(3.0*xp*xp+3.0*yp*yp+zp*zp)/(rp*rp*rp*rp*rhop*rhop*rhop);
            dthetadxdz = xp*(xp*xp+yp*yp-zp*zp)/(rp*rp*rp*rp*rhop);
            dthetadydy = zp*(-2.0*yp*yp*yp*yp-yp*yp*xp*xp+xp*xp*xp*xp+zp*zp*xp*xp)/(rp*rp*rp*rp*rhop*rhop*rhop);
            dthetadydz = yp*(xp*xp+yp*yp-zp*zp)/(rp*rp*rp*rp*rhop);
            dthetadzdz = 2.0*zp*rhop/(rp*rp*rp*rp);
            
            dphidxdx = 2.0*yp*xp/(rhop*rhop*rhop*rhop);  
            dphidxdy = (yp*yp-xp*xp)/(rhop*rhop*rhop*rhop);  
            dphidxdz = 0.0;  
            dphidydy = - 2.0*yp*xp/(rhop*rhop*rhop*rhop);  
            dphidydz = 0.0;  
            dphidzdz = 0.0;  
          
            //Compute first derivatives of F: dFdi  

            dFdx[i][j] = drdx;  
            dFdy[i][j] = drdy;  
            dFdz[i][j] = drdz;  

            for(l=0;l<=LMAX1;l++){  
              dFdx[i][j] -= a0[l]*dthetadx*dY0dtheta[i1][l];   
              dFdy[i][j] -= a0[l]*dthetady*dY0dtheta[i1][l];   
              dFdz[i][j] -= a0[l]*dthetadz*dY0dtheta[i1][l];   
            }  
              
            for(l=1;l<=LMAX1;l++){  
              for(m=1;m<=l;m++){  
                l1=l*(LMAX1+1)+m;  
                dFdx[i][j] -= ac[l][m]*(dthetadx*dYcdtheta[i1][l1] + dphidx*dYcdphi[i1][l1])   
                            + as[l][m]*(dthetadx*dYsdtheta[i1][l1] + dphidx*dYsdphi[i1][l1]);   
                dFdy[i][j] -= ac[l][m]*(dthetady*dYcdtheta[i1][l1] + dphidy*dYcdphi[i1][l1])   
                            + as[l][m]*(dthetady*dYsdtheta[i1][l1] + dphidy*dYsdphi[i1][l1]);   
                dFdz[i][j] -= ac[l][m]*(dthetadz*dYcdtheta[i1][l1] + dphidz*dYcdphi[i1][l1])   
                            + as[l][m]*(dthetadz*dYsdtheta[i1][l1] + dphidz*dYsdphi[i1][l1]);   
              }  
            }  

            //Compute second derivatives of F: dFdidj  

            dFdxdx = drdxdx;  
            dFdxdy = drdxdy;  
            dFdxdz = drdxdz;  
            dFdydy = drdydy;  
            dFdydz = drdydz;  
            dFdzdz = drdzdz;  
            
            for(l=0;l<=LMAX1;l++){  
              dFdxdx -= a0[l]*(dthetadxdx*dY0dtheta[i1][l] + dthetadx*dthetadx*dY0dthetadtheta[i1][l]);   
              dFdxdy -= a0[l]*(dthetadxdy*dY0dtheta[i1][l] + dthetadx*dthetady*dY0dthetadtheta[i1][l]);   
              dFdxdz -= a0[l]*(dthetadxdz*dY0dtheta[i1][l] + dthetadx*dthetadz*dY0dthetadtheta[i1][l]);   
              dFdydx -= a0[l]*(dthetadxdy*dY0dtheta[i1][l] + dthetady*dthetadx*dY0dthetadtheta[i1][l]);   
              dFdydy -= a0[l]*(dthetadydy*dY0dtheta[i1][l] + dthetady*dthetady*dY0dthetadtheta[i1][l]);   
              dFdydz -= a0[l]*(dthetadydz*dY0dtheta[i1][l] + dthetady*dthetadz*dY0dthetadtheta[i1][l]);  
              dFdzdx -= a0[l]*(dthetadxdz*dY0dtheta[i1][l] + dthetadz*dthetadx*dY0dthetadtheta[i1][l]);   
              dFdzdy -= a0[l]*(dthetadydz*dY0dtheta[i1][l] + dthetadz*dthetady*dY0dthetadtheta[i1][l]);    
              dFdzdz -= a0[l]*(dthetadzdz*dY0dtheta[i1][l] + dthetadz*dthetadz*dY0dthetadtheta[i1][l]);   
            }  
          
            for(l=1;l<=LMAX1;l++){  
              for(m=1;m<=l;m++){  
                l1=l*(LMAX1+1)+m;  
                dFdxdx -= ac[l][m]*(dthetadxdx*dYcdtheta[i1][l1] + dthetadx*(dthetadx*dYcdthetadtheta[i1][l1]   
                        + dphidx*dYcdthetadphi[i1][l1]) + dphidxdx*dYcdphi[i1][l1]   
                        + dphidx*(dthetadx*dYcdthetadphi[i1][l1] + dphidx*dYcdphidphi[i1][l1]))   
                        + as[l][m]*(dthetadxdx*dYsdtheta[i1][l1] + dthetadx*(dthetadx*dYsdthetadtheta[i1][l1]   
                        + dphidx*dYsdthetadphi[i1][l1]) + dphidxdx*dYsdphi[i1][l1]   
                        + dphidx*(dthetadx*dYsdthetadphi[i1][l1] + dphidx*dYsdphidphi[i1][l1]));  
                dFdxdy -= ac[l][m]*(dthetadxdy*dYcdtheta[i1][l1] + dthetadx*(dthetady*dYcdthetadtheta[i1][l1]   
                        + dphidy*dYcdthetadphi[i1][l1]) + dphidxdy*dYcdphi[i1][l1]   
                        + dphidx*(dthetady*dYcdthetadphi[i1][l1] + dphidy*dYcdphidphi[i1][l1]))   
                        + as[l][m]*(dthetadxdy*dYsdtheta[i1][l1] + dthetadx*(dthetady*dYsdthetadtheta[i1][l1]   
                        + dphidy*dYsdthetadphi[i1][l1]) + dphidxdy*dYsdphi[i1][l1]   
                        + dphidx*(dthetady*dYsdthetadphi[i1][l1] + dphidy*dYsdphidphi[i1][l1]));  
                dFdxdz -= ac[l][m]*(dthetadxdz*dYcdtheta[i1][l1] + dthetadx*(dthetadz*dYcdthetadtheta[i1][l1]   
                        + dphidz*dYcdthetadphi[i1][l1]) + dphidxdz*dYcdphi[i1][l1]   
                        + dphidx*(dthetadz*dYcdthetadphi[i1][l1] + dphidz*dYcdphidphi[i1][l1]))   
                        + as[l][m]*(dthetadxdz*dYsdtheta[i1][l1] + dthetadx*(dthetadz*dYsdthetadtheta[i1][l1]   
                        + dphidz*dYsdthetadphi[i1][l1]) + dphidxdz*dYsdphi[i1][l1]   
                        + dphidx*(dthetadz*dYsdthetadphi[i1][l1] + dphidz*dYsdphidphi[i1][l1]));  
                dFdydx -= ac[l][m]*(dthetadxdy*dYcdtheta[i1][l1] + dthetady*(dthetadx*dYcdthetadtheta[i1][l1]   
                        + dphidx*dYcdthetadphi[i1][l1]) + dphidxdy*dYcdphi[i1][l1]   
                        + dphidy*(dthetadx*dYcdthetadphi[i1][l1] + dphidx*dYcdphidphi[i1][l1]))   
                        + as[l][m]*(dthetadxdy*dYsdtheta[i1][l1] + dthetady*(dthetadx*dYsdthetadtheta[i1][l1]   
                        + dphidx*dYsdthetadphi[i1][l1]) + dphidxdy*dYsdphi[i1][l1]   
                        + dphidy*(dthetadx*dYsdthetadphi[i1][l1] + dphidx*dYsdphidphi[i1][l1]));  
                dFdydy -= ac[l][m]*(dthetadydy*dYcdtheta[i1][l1] + dthetady*(dthetady*dYcdthetadtheta[i1][l1]   
                        + dphidy*dYcdthetadphi[i1][l1]) + dphidydy*dYcdphi[i1][l1]   
                        + dphidy*(dthetady*dYcdthetadphi[i1][l1] + dphidy*dYcdphidphi[i1][l1]))   
                        + as[l][m]*(dthetadydy*dYsdtheta[i1][l1] + dthetady*(dthetady*dYsdthetadtheta[i1][l1]   
                        + dphidy*dYsdthetadphi[i1][l1]) + dphidydy*dYsdphi[i1][l1]   
                        + dphidy*(dthetady*dYsdthetadphi[i1][l1] + dphidy*dYsdphidphi[i1][l1]));  
                dFdydz -= ac[l][m]*(dthetadydz*dYcdtheta[i1][l1] + dthetady*(dthetadz*dYcdthetadtheta[i1][l1]   
                        + dphidz*dYcdthetadphi[i1][l1]) + dphidydz*dYcdphi[i1][l1]   
                        + dphidy*(dthetadz*dYcdthetadphi[i1][l1] + dphidz*dYcdphidphi[i1][l1]))   
                        + as[l][m]*(dthetadydz*dYsdtheta[i1][l1] + dthetady*(dthetadz*dYsdthetadtheta[i1][l1]   
                        + dphidz*dYsdthetadphi[i1][l1]) + dphidydz*dYsdphi[i1][l1]   
                        + dphidy*(dthetadz*dYsdthetadphi[i1][l1] + dphidz*dYsdphidphi[i1][l1]));  
                dFdzdx -= ac[l][m]*(dthetadxdz*dYcdtheta[i1][l1] + dthetadz*(dthetadx*dYcdthetadtheta[i1][l1]   
                        + dphidx*dYcdthetadphi[i1][l1]) + dphidxdz*dYcdphi[i1][l1]   
                        + dphidz*(dthetadx*dYcdthetadphi[i1][l1] + dphidx*dYcdphidphi[i1][l1]))   
                        + as[l][m]*(dthetadxdz*dYsdtheta[i1][l1] + dthetadz*(dthetadx*dYsdthetadtheta[i1][l1]   
                        + dphidx*dYsdthetadphi[i1][l1]) + dphidxdz*dYsdphi[i1][l1]   
                        + dphidz*(dthetadx*dYsdthetadphi[i1][l1] + dphidx*dYsdphidphi[i1][l1]));  
                dFdzdy -= ac[l][m]*(dthetadydz*dYcdtheta[i1][l1] + dthetadz*(dthetady*dYcdthetadtheta[i1][l1]   
                        + dphidy*dYcdthetadphi[i1][l1]) + dphidydz*dYcdphi[i1][l1]   
                        + dphidz*(dthetady*dYcdthetadphi[i1][l1] + dphidy*dYcdphidphi[i1][l1]))   
                        + as[l][m]*(dthetadydz*dYsdtheta[i1][l1] + dthetadz*(dthetady*dYsdthetadtheta[i1][l1]   
                        + dphidy*dYsdthetadphi[i1][l1]) + dphidydz*dYsdphi[i1][l1]   
                        + dphidz*(dthetady*dYsdthetadphi[i1][l1] + dphidy*dYsdphidphi[i1][l1]));  
                dFdzdz -= ac[l][m]*(dthetadzdz*dYcdtheta[i1][l1] + dthetadz*(dthetadz*dYcdthetadtheta[i1][l1]   
                        + dphidz*dYcdthetadphi[i1][l1]) + dphidzdz*dYcdphi[i1][l1]   
                        + dphidz*(dthetadz*dYcdthetadphi[i1][l1] + dphidz*dYcdphidphi[i1][l1]))   
                        + as[l][m]*(dthetadzdz*dYsdtheta[i1][l1] + dthetadz*(dthetadz*dYsdthetadtheta[i1][l1]   
                        + dphidz*dYsdthetadphi[i1][l1]) + dphidzdz*dYsdphi[i1][l1]   
                        + dphidz*(dthetadz*dYsdthetadphi[i1][l1] + dphidz*dYsdphidphi[i1][l1]));  
              }  
            }  

            //Compute dFdi with the index up  
              
            dFupdx = ginvxx*dFdx[i][j] + ginvxy*dFdy[i][j] + ginvxz*dFdz[i][j];  
            dFupdy = ginvyx*dFdx[i][j] + ginvyy*dFdy[i][j] + ginvyz*dFdz[i][j];  
            dFupdz = ginvzx*dFdx[i][j] + ginvzy*dFdy[i][j] + ginvzz*dFdz[i][j];  

            //Compute norm of dFdi  

            sum = dFupdx*dFdx[i][j] + dFupdy*dFdy[i][j] + dFupdz*dFdz[i][j];  
            if(sum>=0){  
              u[i][j] = sqrt(sum);  
              u0=0;  
            }else{  
              u[i][j] = 0.0;  
              u0=1;  
            }  

            //Compute nabla_a nabla_b F  
          
            nnFxx = dFdxdx - 0.5*(dFupdx*dgxxdx + dFupdy*(2.0*dgxydx-dgxxdy) + dFupdz*(2.0*dgxzdx-dgxxdz));  
            nnFxy = dFdxdy - 0.5*(dFupdx*dgxxdy + dFupdy*dgyydx + dFupdz*(dgyzdx+dgxzdy-dgxydz));  
            nnFxz = dFdxdz - 0.5*(dFupdx*dgxxdz + dFupdy*(dgyzdx + dgxydz - dgxzdy) + dFupdz*dgzzdx);  
            nnFyy = dFdydy - 0.5*(dFupdy*dgyydy + dFupdx*(2.0*dgxydy - dgyydx) + dFupdz*(2.0*dgyzdy-dgyydz));  
            nnFyz = dFdydz - 0.5*(dFupdy*dgyydz + dFupdx*(dgxzdy + dgxydz - dgyzdx) + dFupdz*dgzzdy);  
            nnFzz = dFdzdz - 0.5*(dFupdz*dgzzdz + dFupdx*(2.0*dgxzdz - dgzzdx) + dFupdy*(2.0*dgyzdz-dgzzdy));  

            //Compute d2F = g^{ab} nabla_a nabla_b F   

            d2F = ginvxx*nnFxx + ginvxy*nnFxy + ginvxz*nnFxz   
                + ginvyx*nnFxy + ginvyy*nnFyy + ginvyz*nnFyz  
                + ginvzx*nnFxz + ginvzy*nnFyz + ginvzz*nnFzz;   

            //Compute dFd^a dFd^b Kab  
              
            dFdadFdbKab = dFupdx*dFupdx*Kxx + 2.0*dFupdx*dFupdy*Kxy + 2.0*dFupdx*dFupdz*Kxz  
                        + dFupdy*dFupdy*Kyy + 2.0*dFupdy*dFupdz*Kyz + dFupdz*dFupdz*Kzz;  

            //Compute dFd^a dFd^b nabla_a nabla_b F  
              
            dFdadFdbFdadb = dFupdx*dFupdx*nnFxx + 2.0*dFupdx*dFupdy*nnFxy + 2.0*dFupdx*dFupdz*nnFxz   
                          + dFupdy*dFupdy*nnFyy + 2.0*dFupdy*dFupdz*nnFyz + dFupdz*dFupdz*nnFzz;  

            //Compute H  

            if(u0==0){  
              H[i][j] = d2F/u[i][j] + dFdadFdbKab/(u[i][j]*u[i][j]) - dFdadFdbFdadb/(u[i][j]*u[i][j]*u[i][j]) - K;  
            }else{  
              H[i][j] = 0.0;  
            }  


            //Compute sigma  
            
            sigma[i][j] = 1.0;  

            //Compute rho = H*u*sigma;  

            rho[i][j] = H[i][j]*u[i][j]*sigma[i][j];  

            //Induced metric on the horizon.  
            h11 = dxdt[i1]*(dxdt[i1]*gxx+dydt[i1]*gxy+dzdt[i1]*gxz)   
                + dydt[i1]*(dxdt[i1]*gxy+dydt[i1]*gyy+dzdt[i1]*gyz)   
                + dzdt[i1]*(dxdt[i1]*gxz+dydt[i1]*gyz+dzdt[i1]*gzz);  
            h12 = dxdt[i1]*(dxdp[i1]*gxx+dydp[i1]*gxy+dzdp[i1]*gxz)   
                + dydt[i1]*(dxdp[i1]*gxy+dydp[i1]*gyy+dzdp[i1]*gyz)   
                + dzdt[i1]*(dxdp[i1]*gxz+dydp[i1]*gyz+dzdp[i1]*gzz);  
            h22 = dxdp[i1]*(dxdp[i1]*gxx+dydp[i1]*gxy+dzdp[i1]*gxz)   
                + dydp[i1]*(dxdp[i1]*gxy+dydp[i1]*gyy+dzdp[i1]*gyz)   
                + dzdp[i1]*(dxdp[i1]*gxz+dydp[i1]*gyz+dzdp[i1]*gzz);  
            
            //Determinant of the induced metric.  
            deth[i][j] = h11*h22-h12*h12;  
            if(deth[i][j]<0) deth[i][j] = 0.0;  

            //Flat space coordinate rotational killing vectors.  

            phix_x =  0;  
            phix_y = -(z-zc);  
            phix_z =  (y-yc);  
            phiy_x =  (z-zc);  
            phiy_y =  0;  
            phiy_z = -(x-xc);  
            phiz_x = -(y-yc);  
            phiz_y =  (x-xc);  
            phiz_z =  0;  

            //Normal vector  

            Rx = dFupdx/u[i][j];
            Ry = dFupdy/u[i][j];
            Rz = dFupdz/u[i][j];

            //Integrands of S.

            intSx[i][j] = phix_x*Rx*Kxx + phix_x*Ry*Kxy + phix_x*Rz*Kxz
                        + phix_y*Rx*Kxy + phix_y*Ry*Kyy + phix_y*Rz*Kyz
                        + phix_z*Rx*Kxz + phix_z*Ry*Kyz + phix_z*Rz*Kzz;
            intSy[i][j] = phiy_x*Rx*Kxx + phiy_x*Ry*Kxy + phiy_x*Rz*Kxz
                        + phiy_y*Rx*Kxy + phiy_y*Ry*Kyy + phiy_y*Rz*Kyz
                        + phiy_z*Rx*Kxz + phiy_z*Ry*Kyz + phiy_z*Rz*Kzz;
            intSz[i][j] = phiz_x*Rx*Kxx + phiz_x*Ry*Kxy + phiz_x*Rz*Kxz
                        + phiz_y*Rx*Kxy + phiz_y*Ry*Kyy + phiz_y*Rz*Kyz
                        + phiz_z*Rx*Kxz + phiz_z*Ry*Kyz + phiz_z*Rz*Kzz;

            intSx[i][j] = intSx[i][j]*sqrt(deth[i][j]);
            intSy[i][j] = intSy[i][j]*sqrt(deth[i][j]);
            intSz[i][j] = intSz[i][j]*sqrt(deth[i][j]);
          }
        } //end of loop over points of sphere

        
        timer_stop(0, "AHmod_local");

        //determine the processor with the highest rank for each point
        bampi_allreduce_max_vector(localp, globalp, ntheta*nphi); 

        // check that each point is at least in one processor
        cont = 1;
        for(i1=0;i1<ntheta*nphi;i1++){
          if (globalp[i1] == -1) cont=0;
        }
        if (cont==0){
          if (rank==0) printf("AHmod failed:  point not found in the processors!\n");
          break;
        }

        //Integrals over the surface, return in v
        IntegrateOverSurface(v, ntheta, nphi, dtheta, dphi, rank, globalp, 
                             deth, H, intSx, intSy, intSz, rr);

        area  = v[0];
        hmean = v[1];
        hrms  = v[2];
        Sx    = v[3];
        Sy    = v[4];
        Sz    = v[5];
        coarea= v[6];
 
        mass_old = mass;
        mass = sqrt(area/(16.0*PI));
        meanradius = a0[0]/sqrt(4.0*PI);
        hrms/= area;

        Sx = Sx/(8.0*PI);
        Sy = Sy/(8.0*PI);
        Sz = Sz/(8.0*PI);
        S = sqrt(Sx*Sx + Sy*Sy + Sz*Sz);

        if (verbose && (rank==0)) {
          printf("  mass=%18.12e,  meanradius=%18.12e,  hmean=%18.12e\n",
                    mass, meanradius, hmean );
        }
        

        // check several conditions to stop iteration
        cont = 1;

        if(fabs(hmean)>Getd("AHmod_hmean_tol")){
          cont = 0;
          if ((rank==0) && verbose) 
            printf(" AHmod failed:  hmean > hmean_tol\n");
        }
       
        if(k==flow_iter-1){
          cont = 0;
          if ((rank==0) && verbose) 
            printf(" AHmod failed:  k = flow_iter\n");
        }
               
        if (meanradius<0.0) {
          cont = 0;
          if ((rank==0) && verbose) 
            printf(" AHmod failed:  surface negative height\n");
        }
        
        if (cont == 0)
        {
          // value -1  indicates surface not found
          AHmod_set_global_values(number, time,xc,yc,zc, -1.,-1.);
          
          // erase save of last surface -> no initial guess from file
          AHmod_RemoveLastSurface(number, outdir);
          break; //leave flow loop
        }
        
        // End flow when mass difference is small
        if( fabs(mass_old-mass) < Getd("AHmod_mass_tol") )  // MOTS found
        {
          AHmod_found[number]=1;
          AHmod_set_global_values(number, time, xc,yc,zc, meanradius, mass);
          
          if(rank==0) 
          {      // let only one processor do the output
            
            if (verbose){
              printf("\nHorizon found\n");
              printf("  mass = %22.16e\n",mass);
              printf("  radius = %22.16e\n",meanradius);
              printf("  hrms = %22.16e\n",hrms);
              printf("  hmean = %22.16e\n",hmean);
              printf("  Sx = %22.16e\n",Sx);
              printf("  Sy = %22.16e\n",Sy);
              printf("  Sz = %22.16e\n",Sz);
              printf("  S = %22.16e\n\n",S);
            }

            AHmod_WriteSurfaceResults(number, outdir, verbose,
                                      time,xc,yc,zc,mass,Sx,Sy,Sz,S,coarea,level->l);

            if (Getv("AHmod_uselast","yes"))
              AHmod_WriteSurfaceShape(number, outdir, a0, ac, as, LMAX1);
            
            if (Getv("AHmod_output", "yes"))
              AHmod_output(ntheta, dtheta, nphi, dphi, xc, yc, zc, rr, outdir, number, time, AHmod_time);
            
            if (Getv("AHmod_output_xyt", "yes"))
              AHmod_output_xyt_vtk(ntheta, dtheta, nphi, dphi, xc, yc, zc, rr, outdir, number, time);
            
            if (Getv("AHmod_output_lm", "yes"))
              AHmod_output_sphercoeff(LMAX1, xc, yc, zc, a0, ac, as, outdir, number, time, AHmod_time);
            
          }
          break; // end current flow loop
        }
        
        // find new spectral components
        ComputeFlowSpectralComponents(a0, ac, as, LMAX1, ntheta, nphi, 
                                      rank, dtheta, dphi, Y0, Yc, Ys, rho, globalp);
        
      }//end flow's loop
      
      /* Free interploated data */
      for (int direction=0; direction < ntheta*nphi; direction++) {
        free(interpolated_variables[direction]);
      }
      free(interpolated_variables);
      /* free variable lists for 3d data */
      disablevarlist(vl);
      vlfree(vl);
      vlfree(wl);

    }//end loop over number of horizons
    
    /* free */
    for (int i=0;i<nhorizons; i++) {
      free(surfaces[i].ListPunctures); } 
    free(surfaces);

    free_dvector(a0,0,LMAX1);
    free_dmatrix(ac,0,LMAX1,0,LMAX1);
    free_dmatrix(as,0,LMAX1,0,LMAX1);
    free_dmatrix(rr,0,ntheta-1,0,nphi-1);
    free_dmatrix(u,0,ntheta-1,0,nphi-1);
    free_dmatrix(H,0,ntheta-1,0,nphi-1);
    free_dmatrix(rho,0,ntheta-1,0,nphi-1);
    free_dmatrix(Y0,0,ntheta*nphi,0,LMAX1);
    free_dmatrix(Yc,0,ntheta*nphi,0,(LMAX1+1)*(LMAX1+1));
    free_dmatrix(Ys,0,ntheta*nphi,0,(LMAX1+1)*(LMAX1+1));
    free_dmatrix(dY0dtheta,0,ntheta*nphi,0,LMAX1);
    free_dmatrix(dYcdtheta,0,ntheta*nphi,0,(LMAX1+1)*(LMAX1+1));
    free_dmatrix(dYsdtheta,0,ntheta*nphi,0,(LMAX1+1)*(LMAX1+1));
    free_dmatrix(dYcdphi,0,ntheta*nphi,0,(LMAX1+1)*(LMAX1+1));
    free_dmatrix(dYsdphi,0,ntheta*nphi,0,(LMAX1+1)*(LMAX1+1));
    free_dmatrix(dY0dthetadtheta,0,ntheta*nphi,0,LMAX1);
    free_dmatrix(dYcdthetadtheta,0,ntheta*nphi,0,(LMAX1+1)*(LMAX1+1));
    free_dmatrix(dYcdthetadphi,0,ntheta*nphi,0,(LMAX1+1)*(LMAX1+1));
    free_dmatrix(dYcdphidphi,0,ntheta*nphi,0,(LMAX1+1)*(LMAX1+1));
    free_dmatrix(dYsdthetadtheta,0,ntheta*nphi,0,(LMAX1+1)*(LMAX1+1));
    free_dmatrix(dYsdthetadphi,0,ntheta*nphi,0,(LMAX1+1)*(LMAX1+1));
    free_dmatrix(dYsdphidphi,0,ntheta*nphi,0,(LMAX1+1)*(LMAX1+1));
    free_dmatrix(dFdx,0,ntheta-1,0,nphi-1);
    free_dmatrix(dFdy,0,ntheta-1,0,nphi-1);
    free_dmatrix(dFdz,0,ntheta-1,0,nphi-1);
    free_dmatrix(sigma,0,ntheta-1,0,nphi-1);
    free_dmatrix(deth,0,ntheta-1,0,nphi-1);
    free_dvector(localp,0,ntheta*nphi);
    free_dvector(globalp,0,ntheta*nphi);
    free_dvector(dxdt,0,ntheta*nphi);
    free_dvector(dxdp,0,ntheta*nphi);
    free_dvector(dydt,0,ntheta*nphi);
    free_dvector(dydp,0,ntheta*nphi);
    free_dvector(dzdt,0,ntheta*nphi);
    free_dvector(dzdp,0,ntheta*nphi);
    free_dmatrix(intSx,0,ntheta-1,0,nphi-1);
    free_dmatrix(intSy,0,ntheta-1,0,nphi-1);
    free_dmatrix(intSz,0,ntheta-1,0,nphi-1);

    timer_stop(0, "AHmod");

    return 0;
}



