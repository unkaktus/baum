/* adm.h */
/* Bernd Bruegmann, 6/02 */





/* adm.c */
int init_physical_objects(tL *level);
int computeadmconstraints(tL *level);
int adm_init_conftrans(tL *level);
int adm_undo_conftrans(tL *level);
int set_K_initial(tL *level);


/* m files */
void adm_constraints_N(tVarList *ucon,tVarList *uadm);
