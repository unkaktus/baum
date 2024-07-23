/* Z4d.h */
/* Wolfgang Tichy  4/2004 */

extern int METRIC_z4_INDX_VAR;


/* z4.c */
void z4_setbound_G(tL *level, int igi, int iG);
int  z4_algcon_wrapper(tL *level);
int  z4_to_adm(tL *level);
void z4_to_adm_vl(tVarList *unew, tVarList *upre, double c, tVarList *ucur);
int  z4_startup(tL *level);

/* m files */
void z4_algcon(tVarList *u);
void z4_init(tVarList* vl);
void z4_boundary_box(tVarList *unew, tVarList *upre, double c, tVarList *ucur);
void z4_boundary_box_rad(tVarList *unew, tVarList *upre, double c, tVarList *ucur);
void z4_boundary_shell(tVarList *unew, tVarList *upre, double c, tVarList *ucur);
void z4_boundary(tVarList *unew, tVarList *upre, double c, tVarList *ucur);
void z4_rhs_movpunc_N(tVarList *unew, tVarList *upre, double c, tVarList *ucur);
void z4_rhs_add_source(tVarList *unew, tVarList *upre, double c, tVarList *ucur);




