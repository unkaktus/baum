/* bam_eos.c */
/* mth 06/12 */




/*****************************************************************************/
/* EoS wrappers:                                                             */
/*  - all eos functions should return 0 if ok                                */
/*    and return 1 if there is a problem                                     */
/*  - the call from outside should go through EOS.comp(...), which           */
/*    should handle all calls given by the strings                           */
/*  - usage:  EOS.comp('str_vars_in', 'str_dvars_in', 'str_ddvars_in',       */
/*                     'str_vars_out','str_dvars_out','str_ddvars_out',      */
/*                     vars(consistent with number of used chars)            */
/*                    );                                                     */
/*****************************************************************************/


enum {
  DUST, 
  POLY,
  PWP,
  TAB1D,
  IDEAL,
  PWPHOT,
  TAB1DHOT,
  TAB2D,
  TAB3D
};


/* definition of project configuration... fast access */
typedef struct {

  /* internal stuff, should not be called from outside */
  int type, COLD;
  
  double K,GAMMA,GAMMAMO;
  
  int    PWP_N;
  double PWP_LG10P1, PWP_LG10R1, PWP_K0,PWP_GAMMA0;
  double *PWP_RHO,*PWP_K,*PWP_GAMMA,*PWP_a;
  double *PWP_p,*PWP_epsl,*PWP_eta;
  
  double Y, mb, rhomin1D, mb1D;

  double Ymin, Ymax, epsmin, epsmax, rhomin, rhomax, Tmin, Tmax, T;
  
  int (*use1D) (double *p, double *cs2, double *dpdrho, double *dpdepsl, double *rho, double *epsl);
  int (*use1DH)(double *H, double *rho, double *p, double *epsl);
  int (*use2D) (double *p, double *cs2, double *dpdrho, double *dpdepsl, double *rho, double *epsl);
  int (*use3D) (double *rho, double *epsl, double *Y, 
                double *p, double *T, double *cs2, double *dpdepsl, double *dpdrho);
  int (*use3DT)(double *rho, double *epsl, double *Y, 
                double *p, double *T, double *cs2, double *dpdepsl, double *dpdrho);
  int (*betaEoS) (double *rho, double *Y, double *T);
  int (*micro) (double *rho, double *Y, double *T, double *mun, double *mup, double *mue, double *A, double *Z, double *nh);
  int (*mu) (double Logr, double LogT, double Y, double *mun, double *mup, double *mue);  
  int (*extrap) (double *p, double *cs2, double *dpdrho, double *dpdepsl, double *rho, double *epsl);
  int (*use1Dp)  (double *p, double *cs2, double *drhodp, double *dedp, double *rho, double *e);
  int (*use3DT_D) (double *rho, double *T, double *Y, double *chi, double *kappa);

  /* validity range */

  void (*validity) ();
  int (*extend) ();

  int USETABT0;

  /* external functions to call */
  int (*comp)(char*,char*,char*, char*,char*,char*, ...);
  int (*nls_free)(double rho, double T, double Y, double *Qe, double *Qa, double *Qx,
		double *Re, double *Ra, double *Rx);
  int (*nls_diff)(double rho, double T, double Y, double chi_e, double chi_a, double chi_x,
              double *Qdiff, double *Rdiff, double *leak_tau);
  int (*zeta) ();

  int (*M1) ();
  int (*beta_eq) (double Logr, double LogT, double Y, double *Yeq, double *Teq);  

} tEOS;


/* infos about all specific function ...*/
extern tEOS EOS;









