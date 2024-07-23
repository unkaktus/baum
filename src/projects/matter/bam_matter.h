/* bam_matter.h */

/* definition of project configuration... fast access */
typedef struct
{
  // number of variables used
  int NVq;   // no conservatives (evolved) vars
  int NVf;   // no fluxes (often = NVq)
  int NVw;   // no primitives vars
  int NVadm; // no ADM vars and det3g (11+1)
  int NVg;   // no metric vars (11) alpha, beta^i, g_{ij}
  int NVnls; // no nls vars (9) Q_eff, R_eff, luminosities
  int NVchi;
  int NVo; // no other variables (depend on project, set in grhd_startup)

  // indexes of 1st var in different var lists
  int INDX_VAR_q;   // conservative var
  int INDX_VAR_w;   // primtives
  int INDX_VAR_o;   // other
  int INDX_VAR_adm; // ADM
  int INDX_VAR_gxx; // 3-metric
  int INDX_VAR_alp; // lapse
  int INDX_VAR_bet; // shift
  int INDX_VAR_mask;
  int INDX_VAR_detg;
  int INDX_VAR_camr; // conservative mesh refinement
  int INDX_VAR_nls;  // nls
  int INDX_VAR_chi;

  // hydro atmosphere
  double ATM_RHOATM;
  double ATM_FATM;

  // high-resolution-shock-capturing scheme parameters
  int HRSC_NGHOST;
  int USEATM, USEATM_PRERHS, USEATM_POSTRHS, USEMASK, USECAMR, CAMRACTIVE;
  int USEXCISION;

  char **q_names_list;   // list of cons vars names
  char **w_names_list;   // list of prim vars names
  char **a_names_list;   // list of ADM  vars names
  char **c_names_list;   // list of camr vars names
  char **nls_names_list; // list of NLS  vars names
  char **chi_names_list; // list of chi  vars names
  char **o_names_list;   // list of other vars names

  // function pointers
  void (*init)(); // fkt to initialize matter vars cons+prim
  void (*c2p)();  // ftk to compute con2prim
  void (*sources)();
  void (*set_excision)();
  void (*set_atm_prerhs)();
  void (*set_atm_postrhs)();
  void (*hrsc_flx1d)();
  int (*check)();
  void (*set_mask)();

  // function pointers for hrsc
  void (*rec1d)();
  double (*rec1dl)();
  double (*rec1dr)();
  void (*rec1dm)();
  double (*rec1dml)();
  double (*rec1dmr)();
  double (*hrsc_TVDlim)();
  void (*hrsc_riemsol1d)();
  // HG: setting pointers for rec_change case rho_rec/epsl_rec < 0.
  double (*rec1dsl)(); // called rec1dSAFEr/l
  double (*rec1dsr)();

  // AN: setting pointers for lower rec in case higher rec fails
  double (*rec1dl_low)();
  double (*rec1dr_low)();

  // switch so change from HO-flux to 2nd order LLF-flux
  double flx_LLF_HO_RHO;
  double flx_LLF_HO_ALPHA;

  // NLS key
  int USENLS;
  int ADDNLS;

  // Safe recontruction
  int USESR;

  // function pointer for nls
  void (*nls_init)();
  void (*chi_compute)();
  void (*chi_compute_init)();
  void (*chi_prolong)();
} tMATTER;

/* infos about all specific function ...*/
extern tMATTER MATTER;

tVarList *give_consVars_vl(tL *level);
tVarList *give_matterADMVars_vl(tL *level);
tVarList *give_primVars_vl(tL *level);
tVarList *give_camrVars_vl(tL *level);
tVarList *give_nlsVars_vl(tL *level);
tVarList *give_chiVars_vl(tL *level);
