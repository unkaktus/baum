/* multigrid.h */
/* Bernd Bruegmann 11/02 */

/* for transition from old bam to new bam */
#define vb multigrid_vb
#define stack multigrid_stack
#define uh multigrid_uh
#define vh multigrid_vh
#define wh multigrid_wh
#define fh multigrid_fh
#define uH multigrid_uH
#define vH multigrid_vH
#define wH multigrid_wH
#define fH multigrid_fH

#define add multigrid_add
#define copy multigrid_copy
#define subtract multigrid_subtract

extern int vb;
extern tG *stack;
extern tVarList *uh, *vh, *wh, *fh;
extern tVarList *uH, *vH, *wH, *fH;
extern double tolres2, tolresi, factau2, factaui;


/* exactsolve.c */
void exactsolve(int l);

/* fasmg.c */
void coarsetofine(int lH, int lh);
int finetocoarse(int lh, int lH); 

/* makestack.c */
void makestack(tL *level);
void removestack(tL *level);
void initializestack(tL *level, tVarList *fine);

/* multigrid.c */
void setmgind(tL *g, int *u, int *v, int *w, int *f);
void lgs(tVarList *u, tVarList *v);
void lop(tVarList *u, tVarList *v);

/* restrictprolong.c */
void mg_restrict(tVarList *f, tVarList *c, int order);
void mg_prolong(tVarList *c, tVarList *f, int order);

/* print.c */
void prprim(char *s, tVarList *vl);
void prprime(char *s, tVarList *vl);

/* relax.c */
void rb_gauss_seidel(tVarList *u, tVarList *v, tVarList *f);
void relax(int l, int nsweeps);

/* residual.c */
double findres(int p, tL *gh);
void restau(int lh, double *res2, double *resi, double *tau2, double *taui, 
	    int p);
void prresiduals(int lmax);
void prnormres(int lh); 

/* solve.c */
void fullVdriver(tL *level, double tolres2);
void solve(int l, int lmax, int npre, int npost, int wnmin, int wnmax);

/* utility.c */
tL *levelof(int l);
tL *setlevelh(int lh);
tL *setlevelH(int lH);
void setlevelHh(int lH, int lh);
void setbound(tVarList *u);
void setconst(tVarList *vlu, double c);
void add(tVarList *u, tVarList *v, tVarList *w);
void copyall(tVarList *u, tVarList *v);
void subtract(tVarList *u, tVarList *v, tVarList *w);
double l2norm(tVarList *u);
double linorm(tVarList *u);
void initializeredblack(tL *level);
void prolong_u(int lH, int lh);




