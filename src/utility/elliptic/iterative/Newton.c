/* Newton.c */
/* Wolfgang Tichy 7/2003 */

#include "bam.h"
#include "iterative.h"


int Newton(tL *level, 
  void  (*Fu)(tL *l, tVarList *vl_Fu,  tVarList *vl_u),
  void (*Jdu)(tL *l, tVarList *vl_Jdu, tVarList *vl_du),
  tVarList *vlu, tVarList *vlFu, int itmax, double tol,double *normres, int pr,
  int (*linSolver)(tL *l, tVarList *vl_du, tVarList *vl_Fu, 
                    tVarList *vl_res, tVarList *vl_Coeffs,
                    int lin_itmax, double lin_tol, double *lin_normres,
	            void (*J_du)(tL *, tVarList *, tVarList *), 
	            void (*lin_precon)(tL *, tVarList *, tVarList *)),
  void (*linPrecon)(tL *l, tVarList *Hinv_v, tVarList *v),
  tVarList *vldu, tVarList *vlres, tVarList *vlCoeffs,
  int linSolv_itmax, double linSolv_tolFac )
{
  tG *g = level->grid;
  int i, j, inewton, l;
  double res, sum, lin_normres = 0;
  int lin_its = 0;
  
  //set_refinement_excision(level);

  /* solve with Newton-Raphson iterations: */
  if(pr) printf("Newton:  starting Newton iterations ...\n");
  for (inewton = 0; inewton <= itmax; inewton++)
  {
    /* compute vlFu = F(u) */
    sum = 0;
    for (l = g->lmax; l >= level->l; l--) {
      vlFu->level = g->level[l];
      vlu ->level = g->level[l];
      Fu(g->level[l], vlFu, vlu);
      res = bampi_allreduce_allnorm2(vlFu);
      if (pr) printf("Non-linear residual level %d:  %.4e\n", l, res);
      sum += res*res;
    }
    *normres = sqrt(sum/(g->lmax-level->l+1));

    if(pr)
    {
      printf("Newton step %d: linSolves=%d res=%.4e"
             "  Newton residual = %.4e\n",
             inewton, lin_its, lin_normres, *normres);
      fflush(stdout);
    }
    
    if (*normres <= tol) break;

    /* solve linear equation */
    lin_its=linSolver(level, vldu, vlFu, vlres, vlCoeffs, 
                      linSolv_itmax, (*normres)*linSolv_tolFac, &lin_normres,
	              Jdu, linPrecon);
    if(pr) printf("Newton: after linSolver: %e\n",lin_normres);

    /* do Newton step: u^{n+1} = u^{n} - du */
    for (l = level->l; l <= g->lmax; l++)
      for(j = 0; j < vlu->n; j++) 
      {
	double *u  = g->level[l]->v[vlu ->index[j]]; 
	double *du = g->level[l]->v[vldu->index[j]]; 

	forallpoints(g->level[l], i)
        {
	  u[i] -= du[i];  /* do Newton step: u^{n+1} = u^{n} - du */
	  du[i] = 0;      /* reset du to zero */
	}
      }
    
    /* sync vlu. sync is not needed if du is synced */
    /* bampi_vlsynchronize(vlu); */

    /* boundary conditions for vlu are set in the function
       Fu(level, vlFu, vlu);  supplied by the user.         */
  } 

  /* warn if we didn't converge */
  if (inewton > itmax)
  {
    printf("Newton warning: *** Too many Newton steps! *** \n");
    sum = 0;
    for (l = g->lmax; l >= level->l; l--) {
      vlFu->level = g->level[l];
      vlu ->level = g->level[l];
      Fu(g->level[l], vlFu, vlu);
      res = bampi_allreduce_allnorm2(vlFu);
      if (pr) printf("Non-linear residual level %d:  %.4e\n", l, res);
      sum += res*res;
    }
    *normres = sqrt(sum/(g->lmax-level->l+1));
    printf("Newton: Residual after %d Newton steps:"
           "  Newton residual = %e\n", inewton, *normres);
  }

  //unset_refinement_excision(level);

  return inewton;
}

