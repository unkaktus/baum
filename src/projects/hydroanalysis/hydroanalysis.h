/* hydroanalysis.h */

int hydroanalysis(tL *level);

int modeproject_hydrovars(tL *level);
int Mbar_spheres_hydrovars(tL *level);

void compute_hydrovars(tVarList *u);
void postproc_hydrovars(tL *level);

void compute_ejecta_through_sphere(tL *level);
void compute_radiation_through_sphere(tL *level);
void compute_ejecta_angular(tL *level);
void hydroa_compute_detg_invg_pt(double g11, double g12, double g13, 
			       double g22, double g23, double g33,
			       double *det,
			       double *i11, double *i12, double *i13, 
			       double *i22, double *i23, double *i33);

void compute_magvars(tL *level);