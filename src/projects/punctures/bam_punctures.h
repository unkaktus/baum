/* bam_punctures.h */
/* Bernd Bruegmann, 6/02 , Wolfgang Tichy 6/2003*/

void bam_punctures();
int read_ps_parameters(tL *level);

/* functions available to other modules */
void ReadPunctureParameters(int *ptrKnonzero);
void BY_Wofxyz(double x, double y, double z,
               double *W1, double *W2, double *W3);
void BY_Kofxyz(double x, double y, double z,
               double *K11, double *K12, double *K13,
               double *K22, double *K23, double *K33);
double BY_KKofxyz(double x, double y, double z);

/* Functions to prepare moving punctures */
void Set_alpha_rtoN_atPunc(tL *level);
