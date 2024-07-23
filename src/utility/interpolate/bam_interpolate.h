/* bam_interpolate.h */
/* Bernd Bruegmann, 3/03, Wolfgang Tichy 1/2004 */
/* Jose Gonzalez 2/06 */


enum {
  LAGRANGE,
  WENO,
  WENOZ,  
  LIN,  
};





/* barycentric.c */
double barycentric_single(double x0, int n, int s, double *x, double *y);
void barycentric_omega(int n, int s, double *x, double *omega);
double barycentric(double x0, int n, int s, double *x, double *y,
                   double *omega);
double interpolate_tri_bar(const int order, double x, double y, double z,
                           int n1, int n2, int n3, 
                           double *x1, double *x2, double *x3, double *yp);




/* interpolate.c */
void interpolate_TriN_varlist(tL *level0, tL *level1,
                              tVarList *u0, tVarList *u1, tB *box0, tB *box1,
                              int i, int j, int k, int order);

int check_interpolation_cube_local_withsym(tL *level, 
                                           double x, double y, double z, int order);
int check_interpolation_cubeinbox_local(tB *Box,
                                        double x, double y, double z, int order);

double interpolate_xyz_scalar(tL *level, double x, double y, double z, 
                              int vi, int order, int scheme);
double interpolate_xyz_scalar_withSym(tL *level, double x, double y, double z,
                                      int vi, int order, int scheme);

int interpolate_xyz_local(tL *level, double x, double y, double z, 
                          int nv, int *iv, double *vinterp, 
                          int start, int order, int scheme);
int interpolate_xyz_local_minimal(tL *level, double x, double y, double z, 
                                  int nv, int *iv, double *vinterp, 
                                  int order, int scheme);
int interpolate_xyz_localinbox_minimal(tB *box, double x, double y, double z, 
                                       int nv, int *iv, double *vinterp, 
                                       int order, int scheme);
int interpolate_xyz_local_minimal_withsym(tL *level, double x, double y, double z, 
                                          int nv, int *iv, double *vinterp, 
                                          int order, int scheme);
int interpolate_xyz_localinbox_minimal_withsym(tB *box, double x, double y, double z, 
                                               int nv, int *iv, double *vinterp, 
                                               int order, int scheme);
void map_xyz_withsym(tL* l, double *x, double *y, double *z,
		     int *fx, int *fy, int *fz);
int xyzinsidebbox_withsym(tL* lev, 
                          double *bbox, double x, double y, double z);


/* interpolate_boxtobox.c */
void restrict_varlist(tG *g, int lfine, int lcoarse,
                      tVarList *uf, tVarList *uc, int nbuffer);
void prolong_varlist(tG *g, int lcoarse, int lfine, 
                     tVarList *uc, tVarList *uf, int nbuffer);
