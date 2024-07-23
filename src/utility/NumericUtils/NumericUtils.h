/* NumericUtils.h */
/* (c) Wolfgang Tichy 6.6.2003 */
/* header file for NumericUtils functions */
/* JG 2.06 */

/* Functions */

/* for Newton-Raphson with line searches */
void newton_lnsrch(double x[], int n, int *check,
                   void (*vecfunc)(int, double [], double []),
               	   int MAXITS, double TOLF);
void fd_jacobian(int n, double x[], double fvec[], double **df,
  	void (*vecfunc)(int, double [], double []));

/* 1D and 3D integrals */
double integral(double (*func)(double), double a, double b, double s, int max);
double rombintegral(double (*func)(double), double a, double b,
                    double eps, int max);
double integral3D(double (*int_meth)(double (*f_int)(double), 
                                     double a, double b,
                                     double eps, int max),
                  double (*func)(double, double, double), 
                  double x_limit1, double x_limit2,
                  double (*y_limit1)(double x), double (*y_limit2)(double x),
                  double (*z_limit1)(double x,double y), 
                  double (*z_limit2)(double x,double y),
                  double sx, double sy, double sz, 
                  int maxx, int maxy, int maxz);

/* Newton-Raphson with bracketing */
int rtsafe_itsP(double *x0, 
		void (*funcd)(double x, double *f,double *df, void *par),
		double x1, double x2, void *par, int MAXIT, double xacc);

/* Attenuation functions */
double Attenuation01(double x, double s, double p);


/* Integration over a sphere */
double integral_over_sphere(tL *level, double x0, double y0, double z0,
			     int ntheta, int nphi, double r, int index, int order);

/* Integration over a cubical volume */
double integral_over_volume(tL *level, int Mintvol, 
	int Mintsurfx, int Mintsurfy, int Mintsurfz, int inner_surface);

/* Integration over a cubical volume */
double integral_over_volume_4(tL *level, int Mintvol, 
	int Mintsurfx, int Mintsurfy, int Mintsurfz, int inner_surface);

/* path_integral.c */
double compute_path_integral(tG *g, double x1[4],double x2[4], char *method);

/* fact.c */
double ffact(double n);
double fact(double n);

void SphericalHarmonicY(double *rY, double *iY, 
                        int l, int m, double phi, double costheta);

void SphericalHarmonicYprecomp(double *rY, double *iY,
                               int l, int m, double phi, double theta);

void spinweightedSphericalHarmonic(double *rY, double *iY, 
                        int l, int m, double phi, double costheta);

void spinweightedSphericalHarmonicprecomp(double *rY, double *iY,
                        int l, int m, double phi, double theta);
