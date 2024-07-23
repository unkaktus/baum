/* LinSolve_withHYPRE.c */
/* Wolfgang Tichy 8/14/2003, Bernd Bruegmann 11/2003 */

#include "bam.h"
#include "punctures.h"


/* BB: make it possible to compile without HYPRE by defining dummy functions
   ifdef's are uggly but simple
   to avoid a mess, let's keep all the ifdef stuff in one file
   this macro is to be set in MyConfig by  DFLAGS += -DHYPRE
*/
#ifndef HYPRE

int Punc_hypreWrapper(tL *level, tVarList *x, tVarList *b, tVarList *r, 
                        tVarList *c, int itmax, double tol, double *res,
                        void (*Atimes)(tL *, tVarList *, tVarList *),
                        void (*precon)(tL *, tVarList *, tVarList *))
{
  return 0;
}

int Punc_hypreCleanup(tL *level)
{
  return 0;
}

#else
/*************************************************************************/


/* global point lists */
extern tPointList *Punc_InnerList;
extern tPointList *Punc_PhyBoundList;

/* matrix A for HYPRE */
void *Punc_AMatrix_voidptr;



/* clean up function */
int Punc_hypreCleanup(tL *level)
{
  /* free Matrix A for HYPRE */
  if(Punc_AMatrix_voidptr!=NULL)
  {
    hypre_IJMatrixDestroy(Punc_AMatrix_voidptr);
    Punc_AMatrix_voidptr=NULL;
  }

  return 0;
}




/*************************************************************************/
/* compute finite differencing matrix */

/* compute matrix via reverse engineering */
/* retired for now because slow for 20^3 or more points */
#if 0
void dummy() 
{
  /* point lists could/should be precomputed */
  PunctureSpecialInit(level);
  
  /* set vldu = 0 everywhere */  
  vlConstForPointList(Punc_InnerList, vldu, 0);
  vlConstForPointList(Punc_PhyBoundList, vldu, 0);
  
  /* for all points */
  forallpoints(level, i) {
    
    /* for all variables */
    for(jvar = 0; jvar < nvar; jvar++) {
      
      /* construct short point list for the points in the stencil */
      plist = AllocatePointList(level);
      addtoPointList_check_7(plist, level, i);
      
      /* put a single 1 into du at point i */
      du = VLPtr(vldu, jvar);
      column = (i+offset)*nvar+jvar;
      du[i] = 1;
      
      /* compute column = i*nvar+jvar of matrix, 
	 s.t. => vlr = A_{z, column}  , for  0 <= z < nvar*npoints  */
      Lpuncture_lin_PointList(plist, vlr, vldu);
      
      /* remove the 1 in du at point i */
      du[i] = 0;
      
      /* remove entries at boundary */ 
      vlConstForPointList(Punc_PhyBoundList, vldu, 0);
      
      /* Determine if vlr = A_{row, column} has an entry at 
	 row=n*nvar+kvar.                                         */
      forPointList(plist, n) {
	for (kvar = 0; kvar < nvar; kvar++) {
	      row = (n+offset)*nvar+kvar;
	      r = VLPtr(vlr, kvar);
	      if (r[n] != 0) 
		AddToSparseVector(Arow[row], column, r[n]);
	}
      } endfor;
      FreePointList(plist);
    }
  }
}
#endif




/* construct stencil for linearized elliptic operator:  lu = L(u) */
/* black hole puncture data
   L(u) = Laplace(delta) u - 7 a b (1 + a (1+u0))^-8) u 
*/
void Lpuncture_linear_stencil7(tL *level, int i, 
			       double *u0, double *a, double *b,
			       int *index, double *coefficient)
{
  double cx = 1/(level->dx*level->dx);
  double cy = 1/(level->dy*level->dy);
  double cz = 1/(level->dz*level->dz);
  double cc = -2*(cx + cy + cz);
  int n = 0;

  errorexit("Lpuncture_linear_stencil7: not ready for BOX ... see older BAM version");
}




/* set up matrix problem Ax = b */
int Punc_SetAxb(tL *level, tSparseVector **Arow, 
		       tVarList *vldu, tVarList *vlf, int offset)
{
  tG *g = level->grid;
  double *a  = Ptr(level, "punctures_a");
  double *b  = Ptr(level, "punctures_b");
  double *u0 = Ptr(level, "punctures_u");
  double *du = VLPtr(vldu, 0);
  double *f  = VLPtr(vlf, 0);
  int nvar = vldu->n;
  int i, n;        /* point indices */
  int jvar, kvar;  /* var indices   */
  int column;      /* column of matrix where an entry should be made */
  int row;         /* row of matrix where an entry should be made */
  int offset2;
  int boundary;
  int index[64];
  double coefficient[64];

  if (1) printf("SetMatrixRows level %d, %p, offset %d\n",
		level->l, level, offset);

  /* for all points */
  forallpoints(level, i) {
    boundary = boundaryflag(level, i);
    
    /* set du = 0 everywhere */  
    du[i] = 0;

    /* set finite differencing stencil for interior */
    if (boundary == NOTBOUND) {

      /* get stencil */
      Lpuncture_linear_stencil7(level, i, u0, a, b, index, coefficient);

      /* enter stencil into matrix */
      for (n = 0; n < 7; n++) {
	column = index[n] + offset;
	row    =    i     + offset;
	AddToSparseVector(Arow[row], column, coefficient[n]);
      }
    }
    
    /* all boundaries want a one on the diagonal */
    if (boundary) {
      for (jvar = 0; jvar < nvar; jvar++) {
	column = (i+offset)*nvar+jvar;
	row =    (i+offset)*nvar+jvar;
	AddToSparseVector(Arow[row], column, 1.0);
      }

      /* also, set right hand side to zero */
      f[i] = 0;
    }

    /* set physical boundary */
    if (boundary == PHYBOUND) {

      /* for now have one on diagonal so that x_i = b_i (set b_i, too!) */
      for (jvar = 0; jvar < nvar; jvar++) {
        /* ... */
      }
    }
    
    /* set matrix for refinement boundary */
    if (boundary == REFBOUND) {
      
      /* get interpolation stencil */
      prolong_stencil64(level, i, index, coefficient);

      /* enter interpolation stencil */
      offset2 = offset - g->level[level->l-1]->nnodes;
      for (jvar = 0; jvar < nvar; jvar++) {
	for (n = 0; n < 64; n++) {
	  column = (index[n] + offset2)*nvar + jvar;
	  row =    (i        + offset )*nvar + jvar;
	  AddToSparseVector(Arow[row], column, -coefficient[n]);   /* - (!) */
	}
      }
    }

    /* set matrix for restriction boundary */
    if (boundary == EXCBOUND) {

      /* get interpolation stencil */
      restrict_stencil64(level, i, index, coefficient);

      /* enter interpolation stencil */
      offset2 = offset + level->nnodes;
      for (jvar = 0; jvar < nvar; jvar++) {
	for (n = 0; n < 64; n++) {
	  column = (index[n] + offset2)*nvar + jvar;
	  row =    (i        + offset )*nvar + jvar;
	  AddToSparseVector(Arow[row], column, -coefficient[n]);   /* - (!) */
	}
      }
    }
  } /* end forallpoints */

  /* debug */
  if (0) 
    if (level->l == 1)
      forallpoints(level,i)
	if (i == 1071)
	for(jvar = 0; jvar < nvar; jvar++) {
	  row = (i+offset)*nvar+jvar;
	  printf("%d ", row);
	  prSparseVector(Arow[row]);
	}
  
  return 0;
}




/* determine and set HYPRE parameters */
int Punc_hypreParameters(tL *level, int itmax, double tol, double *res)
{
  int pr = 1;
  char *options;
  char *tolstart;
  char *aftertol;
  char *dummy1;

  /* determine HYPRE tolerance */
  options  = strdup( Gets("HYPRE_options") );
  dummy1   = (char *) calloc(strlen(options) + 20, sizeof(*dummy1));

  tolstart = strstr(options, "-tol");

  if( !Getv("punctures_use_HYPRE_tol", "yes") )
  {
    if(tolstart!=NULL)
    {
      sscanf(tolstart, "%*s%*s%s", dummy1);
      aftertol = strstr(tolstart, dummy1);

      /* put a 0 into options at the start of -tol */
      tolstart[0]=0;
      if(dummy1[0]!=0)
        sprintf(dummy1,"%s%s%s%e", options, aftertol, " -tol ", tol);
      else
        sprintf(dummy1,"%s%s%e", options, "-tol ", tol);
    }
    else
      sprintf(dummy1,"%s %s %e", options, "-tol", tol);

    Sets("HYPRE_options", dummy1);
  }
  else
    if(tolstart!=NULL)
    {
      sscanf(tolstart, "%*s%s", dummy1);
      tol=atof(dummy1);
    }

  free(options);
  free(dummy1);

  if(pr) printf("LinSolve_withHYPRE:\nHYPRE_options = %s\n", 
                Gets("HYPRE_options"));

  if(pr) printf("  res=%e  tol=%e\n", *res, tol);  
  return 0;
}




/*************************************************************************/
/* Wrapper for linear solve of A x = b with HYPRE, 
   returns solution in x.
   Note: c and precon are dummy pointers, which are not used 
         and thus can be 0 when calling hypreLinSolvWrapper .
         Also, itmax is ignored.                                    */
int Punc_hypreWrapper(tL *level, tVarList *x, tVarList *b, tVarList *r, 
                        tVarList *c, int itmax, double tol, double *res,
                        void (*Atimes)(tL *, tVarList *, tVarList *),
                        void (*precon)(tL *, tVarList *, tVarList *))
{
  tG *g = level->grid;
  int lmin = level->l;
  int lmax = g->lmax;
  int l;
  int *offset;
  int nvar    = x->n;
  int npoints;
  int nrows;
  tSparseVector **Arow;
  int row;
  int pr = 1;

  /* determine and set HYPRE parameters */
  Punc_hypreParameters(level, itmax, tol, res);

  /* compute offset for global index list at this level
     for now we take ALL nodes along, even though those with children are
     could be deleted
  */
  offset = calloc(lmax+2, sizeof(int));
  for (l = lmin; l <= lmax; l++)
    offset[l+1] = offset[l] + g->level[l]->nnodes;

  /* allocate an array of matrix rows */
  npoints = offset[lmax+1];
  nrows = nvar*npoints;
  Arow = calloc(nvar*npoints, sizeof(*Arow));
  for (row = 0; row < nrows; row++)
    Arow[row] = AllocateSparseVector();

  /* set matrix rows for all levels */
  if (pr) printf("calling Punc_SetMatrixRows\n");
  for (l = lmin; l <= lmax; l++) {
    level = g->level[l];
    x->level = level;
    b->level = level;
    Punc_SetAxb(level, Arow, x, b, offset[l]);
  }
  if (pr) printf("done    Punc_SetMatrixRows\n");

  /* set matrix A from matrix rows in Arow */
  Punc_AMatrix_voidptr = hypre_SetMatrix_FromLines(nrows, nrows, Arow);

  /* free Arow and its content */
  for (row = 0; row < nrows; row++)
    FreeSparseVector(Arow[row]);
  free(Arow);

  /* call linear Solver in HYPRE */  
  if (pr) printf("calling hypre_SolveLinEq_ForVarList\n");
  level = g->level[lmin];
  r->level = level;
  x->level = level;
  b->level = level;
  hypre_SolveLinEq_ForVarList(level, Punc_AMatrix_voidptr, x, b);
  
  /* compute linear residuum */
  *res = 0;
  for (l = lmax; l >= lmin; l--) {
    PunctureSpecialInit(level);
    x->level = level;
    r->level = level;
    b->level = level;
    Atimes(level, r, x);
    vlsubtract(r, r, b);
    *res += pow(bampi_allreduce_allnorm2(r), 2);
  }
  *res = sqrt(*res/(lmax-lmin+1));

  /* clean up */
  free(offset);
  return 1;
}


/*************************************************************************/
/* end ifdef HYPRE */
#endif
