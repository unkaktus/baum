/* grhd.h */


#define DETGMIN 1e-5
#define ALPHAMIN 1e-5



// ATM enum 
enum {
  GRHD_ATM_NO,
  GRHD_ATM_COLDSTATIC,
  GRHD_ATM_HYDROSTATIC,
  GRHD_ATM_VACUUM
};




// grhd.c 
int  grhd_startup(tL *level);
void grhd_init(tL *level);

// grhd_flx1d.c
void grhd_flx1d_fvpr		(double **g1d, double **w1d, double **q1d, double **o1d, double *mask1d, int nx,int ngx, int direction, double **flux1ds);
void grhd_flx1d_fvprb		(double **g1d, double **w1d, double **q1d, double **o1d, double *mask1d, int nx,int ngx, int direction, double **flux1ds);
void grhd_flx1d_ho		(double **g1d, double **w1d, double **q1d, double **o1d, double *mask1d, int nx,int ngx, int direction, double **flux1ds);
void grhd_flx1d_ho_highdens	(double **g1d, double **w1d, double **q1d, double **o1d, double *mask1d, int nx,int ngx, int direction, double **flux1ds);

void grhd_flx1d_fvprb_Vac	(double **g1d, double **w1d, double **q1d, double **o1d, double *mask1d, int nx,int ngx, int direction, double **flux1d);
void grhd_flx1d_fvpr_re_rWv_Vac	(double **g1d, double **w1d, double **q1d, double **o1d, double *mask1d, int nx,int ngx, int direction, double **flux1d);
void grhd_flx1d_ho_vac_highdens	(double **g1d, double **w1d, double **q1d, double **o1d, double *mask1d, int nx,int ngx, int direction, double **flux1ds);

void grhd_flx1d_fvpr_turb	(double **g1d, double **w1d, double **q1d, double **o1d, double *mask1d, int nx,int ngx, int direction, double **flux1d);
void grhd_flx1d_fvprb_turb	(double **g1d, double **w1d, double **q1d, double **o1d, double *mask1d, int nx,int ngx, int direction, double **flux1d);

void grhd_flx1d_fvpr_HLLC(double **g1d, double **w1d, double **q1d, double **o1d, double *mask1d, int nx,int ngx, int direction, double **flux1d);
void grhd_flx1d_fvprb_HLLC(double **g1d, double **w1d, double **q1d, double **o1d, double *mask1d, int nx,int ngx, int direction, double **flux1d);

void grhd_flx1d_fvpr_ppm(double **g1d, double **w1d, double **q1d, double **o1d, double *mask1d, int nx, int ngx, int direction, double **flux1d);
void grhd_flx1d_fvpr_ppm_HLLC(double **g1d, double **w1d, double **q1d, double **o1d, double *mask1d, int nx, int ngx, int direction, double **flux1d);

// grhd_ppm.c
void Rec1D_PPM_grhd( 
		      double *rho, double *epsl, double *vx, double *vy, double *vz, double *pres, 
		      double *rhomins, double *epslmins, double *vxmins, double *vymins, double *vzmins,
		      double *rhoplus, double *epslplus, double *vxplus, double *vyplus, double *vzplus,
		      int nx);

// grhd_sources.c
void grhd_a2p(tL *level);
void grhd_p2a(tL *level);
void grhd_sources	(tVarList *unew, tVarList *upre, double c, tVarList *ucur, tVarList *primVars, tVarList *otherVars, tVarList *matterADMVars);
void grhd_sources_turb	(tVarList *unew, tVarList *upre, double c, tVarList *ucur, tVarList *primVars, tVarList *otherVars, tVarList *matterADMVars);
void grhd_sources_Vac	(tVarList *unew, tVarList *upre, double c, tVarList *ucur, tVarList *primVars, tVarList *otherVars, tVarList *matterADMVars);

// grhd_p2c.c
void grhd_compute_v2_pt(double vx, double vy, double vz,
                        double gxx, double gxy, double gxz, double gyy, double gyz, double gzz,
                        double *vlowx, double *vlowy, double *vlowz, double *v2);		
void grhd_compute_detg_invg_pt(double g11, double g12, double g13, 
			       double g22, double g23, double g33,
			       double *det,
			       double *i11, double *i12, double *i13, 
			       double *i22, double *i23, double *i33);
void grhd_compute_q_pt(double gxx, double gxy, double gxz, double gyy, double gyz, double gzz, double detg, 
                       double rho, double epsl, double p, double vlowx, double vlowy, double vlowz, double Wlor,
                       double *D, double *T, double *Sx, double *Sy, double *Sz);
void grhd_compute_fx_pt(double alpha, double betax, double detg, 
                        double pres, double vx, 
                        double D, double T, double Sx, double Sy, double Sz,
                        double* fD, double* fT, double* fSx, double* fSy, double* fSz);
void grhd_addturb_fx_pt(double tTau11, double tTau12, double tTau13, double tTau22, double tTau23, double tTau33,
			double alpha, double detg, double gup11, double gup12, double gup13,
			double* fSx, double* fSy, double* fSz);

// grhd_c2p.c, grhd_c2p_v2.c, grhd_c2p_hroot.c 
void grhd_c2p_proot(tVarList *ucur, tVarList *primVars);
void grhd_c2p_rroot_ColdStaticAtm(tVarList *ucur, tVarList *primVars);
void grhd_c2p_proot_ColdStaticAtm(tVarList *ucur, tVarList *primVars);
void grhd_c2p_proot_ColdStaticAtm_hybrid(tVarList *ucur, tVarList *primVars);
void grhd_c2p_hroot(tVarList *ucur, tVarList *primVars);
void grhd_c2p_hroot_ColdStaticAtm(tVarList *ucur, tVarList *primVars);
void grhd_c2p_vanal(tVarList *ucur, tVarList *primVars);
int  grhd_make_c2p(tL *level);

/* grhd_c2p_Vac.c */
void grhd_c2p_proot_Vac(tVarList *ucur, tVarList *primVars);

// grhd_eig.c
void grhd_eig(double alpha, double betax, double vx, double v2, double gxxup, double cs2,
              double *lam0, double *lamm, double *lamp);
void grhd_eigenvals(double alpha, double betax, double vx, double v2,
                    double gxxup, double cs2,
	            double *lam0, double *lamm, double *lamp );
void grhd_eigLR(double gxx, double gxy, double gxz, double gyy, double gyz, double gzz, 
                double alpha, double betax, double gxxup, double detg, 
                double h, double W, double rho, double cs2, double kappa, 
                double vx, double vy, double vz, double vlowx, double vlowy, double vlowz,
                double lam0, double lamm, double lamp, 
                double *L, double *R);
void grhd_eigenvecLR_gen(double gxx, double gxy, double gxz,
                         double gyy, double gyz, double gzz,
                         double alpha, double betax, double gxxup, double detg,
                         double h, double W,
                         double rho, double cs2, double kappa,
                         double vx, double vy, double vz,
                         double vlowx, double vlowy, double vlowz,
                         double lam0, double lamm, double lamp,
                         double *L, double *R);

// grhd_atm.c
void grhd_atm_set_prim_pt(double *p, double *rho, double *epsl, double *vx, double *vy, double *vz, double *v2);
void grhd_atm_set_prerhs(tVarList *unew, tVarList *upre, double c, tVarList *ucur, tVarList *primVars);
void grhd_atm_set_postrhs(tVarList *unew, tVarList *upre, double c, tVarList *ucur, tVarList *primVars);
void grhd_atm_set_postrhs_const(tVarList *unew, tVarList *upre, double c, tVarList *ucur, tVarList *primVars);
void grhd_atm_set_u(tVarList *u, tVarList *primVars);
void grhd_atm_set_mask(tL *level);
void grhd_atm_set_levels(tL *level);
void grhd_atm_set_RHOMAX_global(tL *level);

/* grhd_Vac.c */
void grhd_Vac_set_mask(tVarList *u, tVarList *primVars);
void grhd_Vac_set_u(tVarList *u, tVarList *primVars);
void grhd_Vac_set_prerhs(tVarList *unew, tVarList *upre, double c,
                         tVarList *ucur, tVarList *primVars);
void grhd_Vac_set_levels(tL *level);

// grhd_excision.c
void grhd_excision_init(tL *level);
void grhd_set_excision_fixr(tVarList *upre,tVarList *ucur);
void grhd_set_excision_atm(tVarList *upre,tVarList *ucur);
void grhd_set_excision_vmax(tVarList *upre,tVarList *ucur);
void grhd_set_excision_cut(tVarList *upre,tVarList *ucur);
void grhd_set_excision_lin(tVarList *upre,tVarList *ucur);










