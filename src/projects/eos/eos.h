/* eos.h */
/* mth 06/12 */


/* additions for routines */
/* hgieg 03/20 */

/* eos.c */
int eos_startup(tL *level) ;

double cs2_from_eos_pkc(double rho, double epsl, double p, double dpdepsl, double dpdrho);

int eos_comp(char*in,char*din,char*ddin, char*out,char*dout,char*ddout, ...);







/* eos_analytic.c */
double eos_cs2_rep(double rho, double epsl, double p, double dpdrho, double dpdepsl);

int eos_dummy();

int eos_dust(double *p, double *cs2, double *dpdrho, double *dpdepsl, double *rho, double *epsl);

int eos_poly(double *p, double *cs2, double *dpdrho, double *dpdepsl, double *rho, double *epsl);

int eos_poly_H(double *h, double *rho, double *p, double *epsl);

void eos_load_pwp();
int eos_pwp(double *p, double *cs2, double *dpdrho, double *dpdepsl, double *rho, double *epsl);
int eos_pwp_H(double *h, double *rho, double *p, double *epsl);
int eos_pwpHot(double *p, double *cs2, double *dpdrho, double *dpdepsl, double *rho, double *epsl);

int eos_pwp_p (double *p, double *cs2, double *drhodp, double *dedp, double *rho, double *e);

int eos_ideal(double *p, double *cs2, double *dpdrho, double *dpdepsl, double *rho, double *epsl);

void eos_load_pwp_p();

double EOS_Gamma(double rho, double p, double epsl, double cs2);
double EOS3D_Gamma(double rho, double Y, double T);

/* eos_analytic_WT.c */
double eos_cs2_rep_WT(double rho, double epsl, double p, 
                      double dpdrho, double dpdepsl);
int eos_pwpHot_WT(double *p, double *cs2, double *dpdrho, double *dpdepsl,
                  double *rho, double *epsl);
int eos_ideal_WT(double *p, double *cs2, double *dpdrho, double *dpdepsl,
                 double *rho, double *epsl) ;


/* 1D tables */
int locate1d( double *x, int Nx, double xval);
void interp1d_her4( double *f, double *df, double *x, int Nx, double xv, double *fv_p, double *dfv_p, double *ddfv_p );

void eos_load_tab1d();
int eos_cold_tab1d(double *p, double *cs2, double *dpdrho, double *dpdepsl, double *rho, double *epsl);



/* 3D tables */
void eos_load_tab3d();

int eos_tab3d_shen(double *rho, double *epsl, double *Yp, 
                   double *p, double *T, double *cs2,
                   double *dpdrho, double *dpdepsl);
int eos_tab3d_shen_T(double *rho, double *epsl, double *Y, 
                     double *p, double *T, double *cs2,
                     double *dpdrho, double *dpdepsl);
int eos_tab3d_shen_Y(double *p, double *cs2, double *dpdrho, double *dpdepsl, double *rho, double *epsl);
int eos_tab3d_shen_T0_Y(double *p, double *cs2, double *dpdrho, double *dpdepsl, double *rho, double *epsl);



int eos_tab3d_test(tL *level);


/* eos_tab3D_gieg.c*/


/************************ Eos Loaders ***********************************************************************/

void eos_load_tabT(char *fname, int *n1, int *n2, int *n3,
                        double **x1, double **x2, double **x3,
                        double **v1, double **v2, double **v3, double **v4,
                        double **v5, double **v6, double **v7, double **v8,
                        double **v9, double **v10, double **v11, double *mb);

void eos_load_tabT0(char *fname, int *n1, int *n2, int *n3,
                        double **x1, double **x2, double **x3,
                        double **v1, double **v2, double **v3, double **v4,
                        double **v5, double **v6, double **v7, double **v8,
                        double **v9, double *mb);

void find_mb_T0 (char *fname1, char *fname2, char *fname3, double *mb);

void find_mb (char *fname1, char *fname2, double *mb, double *mb1d);

void eos_load_derivs(char *fname, double **v1, double **v2);

void eos_load_M1_eq(char *fname, double **v1, double **v2);

void eos_load_complete();


/********************************* Interpolations *************************************************/

//Trilinear
int intp3d(double x, double y, double z, double *f, int nx, int ny, int nz, 
           double dx, double dy, double dz, double *f_t, double *x_t, double *y_t, double *z_t,
           double *dfdx, double *dfdy, double *dfdz);

//Cubic splines

void spline(double *x, double *y, int n, double **y2);
void splint(double *xa, double *ya, double *y2a, int n, double x, double *y, double *dydx);

/****************************** Validity Range ****************************************************************/

void eos_validity(double *rho, double *Y, double *rho_min, double *rho_max, 
                  double *Y_min, double *Y_max, double *eps_min, double *eps_max, double *Tmin, 
		  double *Tmax, double *mb, double *rhomin1D, double *mb1D);


int extend_validity(double *rho_hat, double *Y, double *eps_hat, double *T);


/******************************************* Temperature finder and functions computations **************************************************/

double eps_T(double rho, double Y, double T);

int compute_temp(double *T, double rho, double Y, double epsl, double Tmin, double Tmax, double epsmin, double epsmax);

int find_temp(double *T, double *rho, double *Y, double *epsl, double *Tmin, double *Tmax, double *epsmin, double *epsmax);

void cs2_tab(double rho, double T, double Y, double dpdr_T, double dedr_T,
             double dpdT_r, double dedT_r, double *cs2, double *dpdr_e, double *dpde_r);

void peps_tab (double *rho, double *T, double *Y, double *p, double *eps,
               double *dpdr_T, double *dedr_T, double *dpdT_r, double *dedT_r);

int eos_tab3d (double *rho, double *eps, double *Y, double *p,
               double *T, double *cs2, double *dpdrho, double *dpdeps);


int eos_tab3d_T(double *rho, double *eps, double *Y, double *p, double *T, double *cs2,
                double *dpdrho, double *dpdeps);

int eos_tab3d_T_D(double *rho, double *T, double *Y, double *chi, double *kappa);

int eos_tab3d_Y(double *p, double *cs2, double *dpdrho, double *dpdeps, double *rho, double *eps);

int eos_tab3d_T0_Y(double *p, double *cs2, double *dpdrho, double *dpdeps, double *rho, double *eps);

int eos_tab3d_micro (double *rho, double *T, double *Y, double *mun, double *mup, double *mue, double *A, double *Z, double *nh);

int peps_1d (double *p, double *cs2, double *dpdrho,  double *dpdeps, double *rho, double *eps);

int beta_1d(double *rho, double *Y, double *T);

int peps_1d_lin (double *p, double *cs2, double *dpdrho, double *dpdeps, double *rho, double *eps);

int beta_1d_lin (double *rho, double *Y, double *T);

int peps_1d_steffen (double *p, double *cs2, double *dpdrho, double *dpdeps, double *rho, double *eps);

int beta_1d_steffen(double *rho, double *Y, double *T);

int eos_tab3d_mu (double Logr, double LogT, double Y, double *mun, double *mup, double *mue);

int eos_tab3d_eq (double Logr, double LogT, double Y, double *Yeq, double *Teq);


/* eos_tab1D_cold.c */

void eos_read_tab1d(char *fname, int *n, double **x1, double **x2, double **x3, double **v1, double **v2, double *mb);

void eos_tab1d_validity (double *rho, double *Y, double *rho_min, double *rho_max, double *Y_min, double *Y_max, 
			 double *eps_min, double *eps_max, double *mb);

int eos_tab1d_extend (double *rho, double *Y, double *eps);

int eos_tab1d (double *p, double *cs2, double *dpdrho, double *dpdeps, double *rho, double *eps);

void eos_load_tab1d_cold();

int eos_tab1d_beta (double *rho, double *Y, double *T);

int eos_tab1d_extrapolate (double *p, double *cs2, double *dpdrho, double *dpdeps, double *rho, double *eps);

int steffen_init (double *xa, double *ya, int N, double **dya, double **a, double **b, double **c, double **d);

int steffen_intp(double *xa, double *ya, double *a, double *b, double *c, double *d, int N, double x, double *y, double *dydx);

int eos_tab1d_stf (double *p, double *cs2, double *dpdrho, double *dpdeps, double *rho, double *eps);

int eos_tab1d_stf_p (double *p, double *cs2, double *drhodp, double *dedp, double *rho, double *e);

int eos_tab1d_extend_stf (double *rho, double *Y, double *eps);

int eos_tab1d_beta_stf (double *rho, double *Y, double *T);

int eos_tab1d_p (double *p, double *cs2, double *drhodp, double *dedp, double *rho, double *e);

int eos_tab1d_lin (double *p, double *cs2, double *dpdrho,  double *dpdeps, double *rho, double *eps);

int eos_tab1d_beta_lin (double *rho, double *Y, double *T);

int eos_tab1d_lin_p (double *p, double *cs2, double *drhodp, double *dedp, double *rho, double *e);

int eos_tab1d_extend_lin (double *rho, double *Y, double *eps);


/* eos_tab1D_hot.c */

void eos_read_tab1d_hot(char *fname, int *n, double **x1, double **x2, double **x3,
                        double **v1, double **v2, double **v3, double **v4,
                        double **v5, double **v6, double **v7, double **v8,
                        double **v9, double *mb, int key);

void eos_load_tab1d_hot ();

void eos_tab1d_hot_validity (double *rho, double *Y, double *rho_min, double *rho_max, double *Y_min, double *Y_max, 
                             double *eps_min, double *eps_max, double *mb);

int eos_tab1d_hot_extend (double *rho, double *Y, double *eps);

int eos_tab1d_hot (double *p, double *cs2, double *dpdrho,  double *dpdeps, double *rho, double *eps);

int eos_tab1d_hot_micro (double *rho, double *T, double *Y, double *mun, double *mup, double *mue, double *A, double *Z, double *nh);

int eos_tab1d_hot_beta (double *rho, double *Y, double *T);

// NLS

void nls_table(double **Qe, double **Qa, double **Qx, double **Re, double **Ra, double **Rx,
               double **zetae, double **zetaa, double **zetax, double **eta_nue,
               double **ka_0_nue, double **ka_1_nue, double **ks_nue, double **ka_0_nua, double **ka_1_nua, double **ks_nua,
               double **ka_0_nux, double **ka_1_nux, double **ks_nux);

int nls_free(double rho, double T, double Y, double *Qe, double *Qa, double *Qx,
             double *Re, double *Ra, double *Rx);

void NLS_compute_sources (double dx, double dy, double dz, double *gxx, double *gyy, double *gzz,
                          double *rho, double *T, double *Y, double **Q_eff, double **R_eff, double *sQ,
                          double *sR, double **leak_tau, int keyx, int keyy, int keyz);

void get_EoS_vars (double *rho, double *T, double *Y, double *nb, 
                   double **nn, double **np, double **nh,
                   double *eta_np, double *eta_pn, double **eta_e, 
                   double **eta_nue, double **A, double **Z);

void compute_local_emission (double nb, double T, double eta_np, double eta_pn, double eta_e, 
                             double eta_nue, double **Q_loc, double **R_loc);

void compute_diff_emission (double nb, double T, double *nn, double *np, double *nh, double *A, double *Z, double *eta_nue, double *eta_e,
			    double *gxx, double *gyy, double *gzz, double dx, double dy, double dz, 
			    double **R_diff, double **Q_diff, double **leak_tau, int keyx, int keyy, int keyz);

double fermi_integral(int k, double eta);

int nls_diff (double rho, double T, double Y, double chi_e, double chi_a, double chi_x,
              double *Qdiff, double *Rdiff, double *leak_tau);

int nls_zeta (double rho, double T, double Y, double *zetae, double *zetaa, double *zetax);


int eos_m1_rates (int ind, double rho, double T, double Y, double *Q, double *R, double *ka_0,
             double *ka_1, double *ks);








