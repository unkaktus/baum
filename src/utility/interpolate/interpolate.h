/* interpolate.h */
/* Bernd Bruegmann, 3/03 */


#define MAXCUBESIZE 12



/* interpolate */
void interpolate_filldatacube(tL *level, int vi, 
                              int uui[MAXCUBESIZE][MAXCUBESIZE][MAXCUBESIZE], 
                              double uuu[MAXCUBESIZE][MAXCUBESIZE][MAXCUBESIZE], int fillcubesize);
int interpolate_fillindexcube(tL *level, int icorner, 
                              int uui[MAXCUBESIZE][MAXCUBESIZE][MAXCUBESIZE], int fillcubesize);
int interpolate_fillindexcube_box(tB *box, int ii, int jj, int kk, 
                                  int uui[MAXCUBESIZE][MAXCUBESIZE][MAXCUBESIZE], int fillcubesize);
double interpolate_TriN(double x, double y, double z, 
                        double xmin, double ymin, double zmin,
                        double dx, double dy, double dz,
                        double uuu[MAXCUBESIZE][MAXCUBESIZE][MAXCUBESIZE], int N, int scheme);
                        
int check_interpolation_cube_local_withsym(tL *level, double x, double y, double z, int order);
int check_interpolation_cubeinbox_local(tB *box, double x, double y, double z, int order);
void map_xyz_withsym(tL *l, double *x, double *y, double *z, int *fx, int *fy, int *fz);



/* interpolate_boxtobox.c */
void prolong_varlist_boxtobox_N(tL *lc, tL *lf, 
                                tVarList *uc, tVarList *uf, int nbuffer,
                                int order, int scheme);
void restrict_varlist_boxtobox_N(tL *lf, tL *lc, 
                                 tVarList *uf, tVarList *uc, int nbuffer,
                                 int order, int scheme);


/* lagrange.c */
void coefficients_lagrange_N(int n, double x, double xmin, double h, double *c);
double interpolate_lagrange_N(int N, double x, double x0, double h, double *c, double *u);



/* for matter vars */
/* WENO.c */
double interpolate_WENO_N(int N, double x, double xmin, double h, double *c, double *u);
double interpolate_WENOZ_N(int N, double x, double xmin, double h, double *c, double *u);

/* linear.c */
double interpolate_linear(int N, double x, double xmin, double h, double *c, double *u);
double interpolate_linear_lim(int N, double x, double xmin, double h, double *c, double *u);
double interpolate_linear_avg(int N, double x, double xmin, double h, double *c, double *u);


