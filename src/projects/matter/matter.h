/* matter.h */


#define ALPHAMIN 1e-6
#define DETGMIN 1e-6




// hrsc.c
int hrsc_startup(tL* level);


// hrsc_rec1d.c
double MM2( double x, double y );
double MC2( double x, double y );

void rec1d_godunov( double *u,int i, double *ul,double *ur );
void rec1d_avg    ( double *u,int i, double *ul,double *ur );
void rec1d_lintvd ( double *u,int i, double *ul,double *ur );
void rec1d_ceno3  ( double *u,int i, double *ul,double *ur );
void rec1d_ppm4   ( double *u,int i, double *ul,double *ur );
void rec1d_ppm6   ( double *u,int i, double *ul,double *ur );
void rec1d_muscl  ( double *u,int i, double *ul,double *ur );
void rec1d_lag4   ( double *u,int i, double *ul,double *ur );
void rec1d_weno5  ( double *u,int i, double *ul,double *ur );
void rec1d_wenoz  ( double *u,int i, double *ul,double *ur );
void rec1d_weno5  ( double *u,int i, double *up,double *um );
void rec1d_lag6   ( double *u,int i, double *ul,double *ur );

double rec1d_p_godunov(double *u, int i);
double rec1d_p_avg    (double *u, int i);
double rec1d_p_lintvd (double *u, int i);
double rec1d_p_ceno3  (double *u, int i);
double rec1d_p_ppm4   (double *u, int i);
double rec1d_p_muscl  (double *u, int i);
double rec1d_p_lag4   (double *u, int i);
double rec1d_p_weno5  (double *u, int i);
double rec1d_p_wenoz  (double *u, int i);
double rec1d_p_mp5    (double *u, int i);
double rec1d_p_ppm6   (double *u, int i);
double rec1d_p_lag6   (double *u, int i);

double rec1d_m_godunov(double *u, int i);
double rec1d_m_avg    (double *u, int i);
double rec1d_m_lintvd (double *u, int i);
double rec1d_m_ceno3  (double *u, int i);
double rec1d_m_ppm4   (double *u, int i);
double rec1d_m_muscl  (double *u, int i);
double rec1d_m_lag4   (double *u, int i);
double rec1d_m_weno5  (double *u, int i);
double rec1d_m_wenoz  (double *u, int i);
double rec1d_m_mp5    (double *u, int i);
double rec1d_m_ppm6   (double *u, int i);
double rec1d_m_lag6   (double *u, int i);

// hrsc_rec1d_ppm
void PPM_MCslope ( const double *a, double *delta_a, const int i );
double PPM_MCslope_pt ( const double *a, const int i );
void PPM_monoton1D_pt ( const double *a, double *amins, double *aplus, const int i );
void PPM_interp1D( const double *a, double *amins, double *aplus, 
		           const int i, 
		           const double *delta_a );
void PPM_cdsteep1D ( const double *a, const double *delta_a, double *amins, double *aplus, 
		             const int i, 
		             const double epsilon1, const double eta1, const double eta2 );
void PPM_flat1D ( const double *a, double *amins, double *aplus, 
		          const int i, 
		          const double omega_flat, const double om_omega_flat );
void PPM_monoton1D ( const double *a, double *amins, double *aplus, const int i );


// hrsc_numflx1d.c
void hrsc_numflx1d_llf(double *uL,double *uR, double *fL,double *fR, double *lamL,double *lamR, double *f, int nf);
void hrsc_numflx1d_hll(double *uL,double *uR, double *fL,double *fR, double *lamL,double *lamR, double *f, int nf);
void hrsc_numflx1d_hll_grmhd(double *uL,double *uR, double *fL,double *fR, double *lamL,double *lamR, double *f, int nf);

// hrsc_HLLD
double HLL_solver_mag_tet(double *qL, double *qR, 
		                  double *fL, double *fR,
                          double *lamL, double *lamR,
                          double v_interface, double *fHLL, double *qHLL,
                          int nf);
void HLLD_solver(double *qL, double *qR, 
		         double *fL, double *fR, 
		         double *lamL, double *lamR, 
                 double v_interface, double *qHLL, double *fHLL, 
                 double Ptot, double mu, double hL, double hR,
		         double *f, double *q, int nf);

// matter.c
int matter_startup(tL *level);
int matter_init(tL *level);
void matter_c2p(tVarList *u_p, tL *level);
void matter_rhs_prepare(tVarList *unew, tVarList *upre, double c, tVarList *ucur);
void matter_rhs_set_null(tVarList *unew, tVarList *upre, double c, tVarList *ucur);
void matter_rhs(tVarList *unew, tVarList *upre, double c, tVarList *ucur);


// matter_fluxes.c
void matter_fluxes(tVarList *unew, tVarList *upre, double c, tVarList *ucur, tVarList *primVars, tVarList *otherVars);
void add_flx1d_rhs(int direction, tB* box, int dir1, int dir2,
                   tVarList *unew, tVarList *upre, double c, tVarList *ucur, tVarList *primVars, 
                   double **fluxs1d);
void compute_flux1d(int direction, tB* box, int dir1, int dir2,
                    tVarList *unew, tVarList *upre, double c, tVarList *ucur, tVarList *primVars, tVarList *otherVars,
                    double **fluxs1d);

// matter_fluxes_grmhd.c
void matter_fluxes_grmhd(tVarList *unew, tVarList *upre, double c, tVarList *ucur, tVarList *primVars, tVarList *otherVars);
void compute_flux1d_grmhd(int direction, tB* box, int dir1, int dir2,
                    tVarList *unew, tVarList *upre, double c, tVarList *ucur, tVarList *primVars, tVarList *otherVars,
                    double **fluxs1d);

// c_amr.c
int c_amr_mask(tL* level);
int c_amr_time(tL* level);
void c_amr_direction(double value, int *dir1, int *dir2, int *dir3);
void correct_varlist(tG *g, int lfine, int lcoarse, tVarList *uf, tVarList *uc, int nbuffer);
