/* AHmod.h */
/* Jose Gonzalez 03/07 */
/* Norbert Lages 02/09 */

int AHmod(tL *level);

struct AHmod_surface_struct {
  int NumOfPunctures;
  double SearchTimeStart;
  double SearchTimeEnd;
  int *ListPunctures;
};

int AHmod_IsMultipleTime(const double time, const double steptime);
                                 
int AHmod_FileExists(const char * filename);

int AHmod_ComputeSlice(const double currenttime, const double timeInverval);

void AHmod_CalculateFactorialList(const int maxn, double *fac);

double AHmod_dtheta(const int ntheta);

double AHmod_dphi(const int nphi);

double AHmod_theta(const double dtheta, const int i);

double AHmod_phi(const double dphi, const int j);

void AHmod_CheckAndCreateResultVariables(const int nhorizons);

void AHmod_ComputeLegendre(const double theta, const int LMAX1,
  double **P, double **dPdtheta, double **dPdthetadtheta);
  
void AHmod_CalculateSphericalHarmonics(const int ntheta, const int nphi,
       const int LMAX1,
       double **Y0, double **dY0dtheta,  double **dY0dthetadtheta,
       double **Yc,              double **Ys, 
       double **dYcdtheta,       double **dYsdtheta,
       double **dYcdphi,         double **dYsdphi, 
       double **dYcdthetadtheta, double **dYsdthetadtheta,
       double **dYcdthetadphi,   double **dYsdthetadphi,
       double **dYcdphidphi,     double **dYsdphidphi);

void AHmod_GetSearchData(struct AHmod_surface_struct *surfaces, 
                         const int nhorizons);

void AHmod_ResetCoefficientsZero(double *a0, double **ac, double **as, 
                                 const int LMAX1);

void AHmod_GetFilenameLastSurface(const int number, const char *outdir, 
                                  char *filename);

void AHmod_GetPunctureLocation(const int number, double *xp, double *yp, double *zp);
                                  
double AHmod_GetPunctureMass(const int number);
     
int AHmod_ExistPuncture(const int number);

double AHmod_GetLargestDistance(const int numpunc, int *punclist);

double AHmod_sumPunctureMasses(const int numpunc, int *punclist);

void AHmod_GetInitialGuess(const int number, const char *outdir,
           double *a0, double **ac, double **as, const int LMAX1,
           struct AHmod_surface_struct *surfaces);

void AHmod_GetInitialGuessFromPunctures(double *a0, 
                                        const int numpunc, int *punclist);
                                        
void AHmod_GetInitialGuessFromFile(char *filename,
           double *a0, double **ac, double **as, const int LMAX1);

void AHmod_WriteSurfaceShape(const int number, const char *outdir,
           double *a0, double **ac, double **as, const int LMAX1);

void AHmod_RemoveLastSurface(const int number, const char *outdir);
    
void AHmod_GetCentralPoint(const int numpunc, int* punclist, 
                           double *xc, double *yc, double *zc);

double AHmod_Getrmax(const double xc, const double yc, const double zc, tL *level);

void  AHmod_WriteSurfaceResults(const int number, const char *outdir, 
   const int verbose, const double time, 
   const double xc, const double yc, const double zc, const double mass, 
   const double Sx, const double Sy, const double Sz, const double S, 
   const double coarea, const int level);

int GetRadiusFromSphericalHarmonics(double **rr, const double dtheta,
  const double dphi, const int ntheta, const int nphi, const int LMAX1,
  const double r_max, 
  double *a0, double **ac, double **as, double **Y0, double **Yc, double **Ys);
  
void ComputeDerivativesXYZvsTP(double **rr, const double dtheta, const double dphi,
  const int ntheta, const int nphi, double *dxdt, double *dydt, double *dzdt,
  double *dxdp, double *dydp, double *dzdp);

void ComputeFlowSpectralComponents(double *a0, double **ac, double **as, 
         const int LMAX1, const int ntheta, const int nphi, const int rank,
         const double dtheta, const double dphi,  
         double **Y0, double **Yc, double **Ys, double **rho, double *globalp);

void IntegrateOverSurface(double *v, const int ntheta, const int nphi, 
     const double dtheta, const double dphi, const int rank, double *globalp, double **deth,
     double **H, double **intSx, double **intSy, double **intSz, double **rr);

int AHmod_ArePuncturesClose(const int numpunc, int *punclist);
     
int AHmod_WaitUntilClosePunctures(const int number);

void AHmod_set_global_values(const int N, const double t, const double x, const double y,
                             const double z, const double r, const double m);
