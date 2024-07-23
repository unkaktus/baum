/* bam_iterative.h */
/* Bernd Bruegmann 01/00 */

/* should probably use prefix, like bamel_jacobi */

int cgs(tL *l, tVarList *x, tVarList *b, tVarList *r, tVarList *c, 
	int imax, double tol, double *res,
	void (*lop)(tL *, tVarList *, tVarList *), 
	void (*precon)(tL *, tVarList *, tVarList *));

int bicgstab(tL *l, tVarList *x, tVarList *b, tVarList *r, tVarList *c, 
	     int imax, double tol, double *res,
	     void (*lop)(tL *, tVarList *, tVarList *), 
	     void (*precon)(tL *, tVarList *, tVarList *));

int sor(tL *l, tVarList *x, tVarList *b, tVarList *r, tVarList *c, 
	int imax, double tol, double *res,
	void (*lop)(tL *, tVarList *, tVarList *), 
	void (*precon)(tL *, tVarList *, tVarList *));

int jacobi(tL *level, tVarList *x, tVarList *b, tVarList *r, tVarList *c, 
	   int itmax, double tol, double *normres,
	   void (*lop)(tL *, tVarList *, tVarList *), 
	   void (*precon)(tL *, tVarList *, tVarList *));

int gaussseidel(tL *l, tVarList *x, tVarList *b, tVarList *r, tVarList *coeffs,
		int itmax, double tol, double *normres,
		void (*lop)(tL *, tVarList *, tVarList *), 
		void (*gs)(tL *, tVarList *, tVarList *));

int Newton(tL *level, 
  void  (*Fu)(tL *l, tVarList *vl_Fu,  tVarList *vl_u),
  void (*Jdu)(tL *l, tVarList *vl_Jdu, tVarList *vl_du),
  tVarList *vlu, tVarList *vlFu, int itmax, double tol,double *normres, int pr,
  int (*linSlover)(tL *l, tVarList *vl_du, tVarList *vl_Fu, 
                    tVarList *vl_res, tVarList *vl_Coeffs,
                    int lin_itmax, double lin_tol, double *lin_normres,
	            void (*J_du)(tL *, tVarList *, tVarList *), 
	            void (*lin_precon)(tL *, tVarList *, tVarList *)),
  void (*linPrecon)(tL *l, tVarList *Hinv_v, tVarList *v),
  tVarList *vldu, tVarList *vlres, tVarList *vlCoeffs,
  int linSolv_itmax, double linSolv_tolFac );
