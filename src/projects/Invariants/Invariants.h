/* Invariants.h */
/* Bernd Bruegmann, 12/2005 */
/* Jose Gonzalez 01.06 */




/* Invariants.c */
void AddVarModes();
int compute_curvature_invariants(tL *level);

void set_origin(double t, double *x0, double *y0, double *z0);
void rescale_psi4(tL *level);
void analytic_Psi4_kinnersley_teukolskywave(tL *level, char *outdir);


int check_mode_for_output(int l, int m);
void compute_modes(tL *level, int n0,int n, double *rlist, char *outdir);
void compute_energy(tL *level, int n0,int n, double *rlist, char *outdir);


void output_invariants(tL *level, int n0,int nr, double *rlist);




/**/
void curvature_invariants_N(tVarList *u);
void dphi_int2_psi4_N(tVarList *u);