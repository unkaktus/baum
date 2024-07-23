/* bam_adm.h */
/* Bernd Bruegmann, 6/02 */



/* adm_data.c */
int adm_Minkowski(tL *level);
int adm_puncture(tL *level);
int adm_gauge_wave(tL *level);


/* adm_utility.c */
double detg(double g11, double g12, double g13, 
	    double g22, double g23, double g33);
double invg(double g11, double g12, double g13, 
	    double g22, double g23, double g33,
	    double *i11, double *i12, double *i13, 
	    double *i22, double *i23, double *i33);
