/* bssn.h */
/* Bernd Bruegmann, 6/02 */
/* Wolfgang Tichy  4/2004 */

extern int METRIC_bssn_INDX_VAR;


/* bssn.c */
void bssn_setbound_G(tL *level, int igi, int iG);
void setbound_eta(tL *level, int ieta);
void synchAMR_eta(tL *level);
int  bssn_to_adm(tL *level);
void bssn_to_adm_vl(tVarList *unew, tVarList *upre, double c, tVarList *ucur);
int  bssn_algcon_wrapper(tL *level);
int  bssn_startup(tL *level);

/* bssn_eta.c */
void bssn_eta_init(tG* grid);
double bssn_eta_set(double x,double y,double z, double absdpsim2, double psim2);



/* m files */
void bssn_algcon(tVarList *u);
void bssn_init(tVarList* vl);
void bssn_boundary(tVarList *unew, tVarList *upre, double c, tVarList *ucur);
void bssn_rhs_movpunc_N(tVarList *unew, tVarList *upre, double c, tVarList *ucur);
void bssn_rhs_add_source(tVarList *unew, tVarList *upre, double c, tVarList *ucur);


