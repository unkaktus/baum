/* evolve.h */
/* Bernd Bruegmann, 6/02 */


extern tVarList *u_c, *u_p, *u_q, *u_r, *u_aux;

int evolve(tL *level);
int evolve_store_rhs(tL *level, int comp);
int evolve_test_startup(tL *level);
int evolve_test_analyze(tL *level);

void evolve_euler(tL *level, int comp);
void evolve_icn(tL *level, int comp);
void evolve_ngs(tL *level, int comp);
void evolve_rk(tL *level, char *name, int comp);
void evolve_rad_imex(tL *level);

int evolve_rk_norder(void);








