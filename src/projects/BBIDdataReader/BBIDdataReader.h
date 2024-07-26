/* BBIDdataReader.h */
/* mth 11/11 */


/* BBIDdataReader.c */
void read_TOV(tL *level, char *file, 
              double bhx,double bhy,double bhz, 
              double bhpx,double bhpy,double bhpz);
void boost_spacetime(tL *level, double bhpx,double bhpy,double bhpz);
void add_spacetime(tL* level, int subMinkowsky);
void set_PUNC(tL* level, double M,
              double bhx,double bhy,double bhz, 
              double bhpx,double bhpy,double bhpz);
void set_BH(tL* level, double M,
              double bhx,double bhy,double bhz, 
              double bhpx,double bhpy,double bhpz,
              double bhsx,double bhsy,double bhsz);
int set_binary_boost(tL* level);

int BBID_read_ps_parameters(tL *level);

/* BBID_solve.c */
int BBID_solve(tL* level);


/* tov.c */
int tov_h(double hc, int npts, 
          double *p_r, double *p_m, double *p_h, 
          double *p_rho, double *p_pre,double *p_phi,
          double *p_riso);

int tov_r(double rhoc, double R, int *npts, 
          double **p_r, double **p_m,
          double **p_rho, double **p_pre,double **p_phi, 
          double **p_riso);//, double *Lamb);

double computeLambda (double M, double R, double G, double H); 

double computeRinv (double *r, double *m, int n);

int test_tov_MvsR(tL *level);

int tov_p(double pc, double R, int *npts, 
          double **p_r, double **p_m,
          double **p_rho, double **p_pre,double **p_phi, 
          double **p_riso, double *Lamb);


int test_tov_MvsR_p(tL *level);

