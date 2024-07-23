/* bam_evolve.h */
/* Bernd Bruegmann, 6/02 */

int evolve(tL *level);

tVarList *get_evolve_vlregister(tL *level);
void set_rhs_func_register(char *name, void (*f)(tVarList *, tVarList *, double, tVarList *));

void evolve_vlretrieve(
  tVarList **vlu_c, tVarList **vlu_p, tVarList **vlu_pp);
void evolve_vlretrieve_cpq(
  tVarList **vlu_c, tVarList **vlu_p, tVarList **vlu_q);
void evolve_vlretrieve_pn(
  tVarList **vlu_p2, tVarList **vlu_p3, tVarList **vlu_p4);

double get_dissipation_factor(tL *level);

typedef struct {

	int stage;

} tEVO;


extern tEVO EVO;




