/* boundary.h */
/* Bernd Bruegmann 01/00 */



/* mask.c */
void dirichlet3(tL *level, int nth, int b, int i0, int i1, int i2, 
		double **var, int nvars, double *x, double *y, double *z);

double lapoly3(double a, double b, double c, double ua, double ub, double uc);
double dintersectcube(double x0, double y0, double z0, 
		      double x1, double y1, double z1, 
		      double xc, double yc, double zc, double r);
double dintersectsphere(double x0, double y0, double z0, 
			double x1, double y1, double z1, 
			double xc, double yc, double zc, double r);

/* periodic.c */
void boxfillbuffer(tL *level, tVarList *vl, double *bbox, 
		   double *buffer, int n);
void boxreadbuffer(tL *level, tVarList *vl, double *bbox, 
		   double *buffer, int n);

/* symmetry.c */
void set_boundary_reflection_one(tL *level, tVarList *varlist, 
				 double *xp, double dx, int dir);




/* extrapolate.c background.c */
void extrapolate_layer(tVarList *ucur);
#define ifnewPHYBOUND(i,j,k, box, off)\
    int Di,Dj,Dk;\
    const int imax = box->m-1;\
    const int jmax = box->n-1;\
    const int kmax = box->o-1;\
    if ((i<0+off) || (j<0+off) || (k<0+off) || (i>imax-off) || (j>jmax-off) || (k>kmax-off)) continue;\
    if ((i!=0+off) && (j!=0+off) && (k!=0+off) && (i!=imax-off) && (j!=jmax-off) && (k!=kmax-off)) continue;\
    if (i==0+off)          Di = 1;\
    else if (i==imax-off)  Di =-1;\
    else                   Di = 0;\
    if (j==0+off)          Dj = 1;\
    else if (j==jmax-off)  Dj =-1;\
    else                   Dj = 0;\
    if (k==0+off)          Dk = 1;\
    else if (k==kmax-off)  Dk =-1;\
    else                   Dk = 0;\
    int d0 = ((Di>0) && (box->bflag[0]==PHYBOUND))?1:0;\
    int d1 = ((Di<0) && (box->bflag[1]==PHYBOUND))?1:0;\
    int d2 = ((Dj>0) && (box->bflag[2]==PHYBOUND))?1:0;\
    int d3 = ((Dj<0) && (box->bflag[3]==PHYBOUND))?1:0;\
    int d4 = ((Dk>0) && (box->bflag[4]==PHYBOUND))?1:0;\
    int d5 = ((Dk<0) && (box->bflag[5]==PHYBOUND))?1:0;\
    if (d0+d1+d2+d3+d4+d5==0) continue;
double extrapolate_N(int N, double *v);
void set_boundary_extrapolate_shells(tL *level, tVarList *vlu);


/* normals.c */
int give_cube_normal (int i, int j, int k, 
		      int imax, int jmax, int kmax, int o,
		      double *shat1, double *shat2, double *shat3);






