/* bam_multigrid.h */
/* Bernd Bruegmann, 11/02 , Wolfgang Tichy 3/2005 */

int multigrid(tL *l, tVarList *x, tVarList *b, tVarList *r, tVarList *c,
	      int itmax, double tol, double *normres,
	      void (*lop)(tL *, tVarList *, tVarList *), 
	      void (*precon)(tL *, tVarList *, tVarList *));

/* pointer which can be used together with par mg_setbound to call
   special boundary functions in multigrid                         */
extern void (*multigrid_setbound)(tL *, tVarList *);

/* helper */
void prvar01(tL *level, char *name);
void prvar(tL *level, char *name);
void prvarbasic(tL *level, char *name, int formattype, int slices);
void prvare(tL *level, char *name);
