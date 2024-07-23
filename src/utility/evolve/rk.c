/* rk.c */
/* Bernd Bruegmann 10/03 */

/* First version in admc_rk.c for Cactus, Bernd Bruegmann 6/99.
   Now generalized with the more general setup found in Mathematica
   http://documents.wolfram.com/v5/Built-inFunctions/NumericalComputation/
     EquationSolving/AdvancedDocumentation/NDSolve.html
   See also the description for Matlab, 
   "Behind and Beyond the MATLAB ODE Suite"
   okumedia.cc.osaka-kyoiku.ac.jp/~ashino/pdf/2651.pdf
   www.site.uottawa.ca/~remi/odesuite.ps

   Issues:
   - which of the "modern" Runge-Kutta features matter for method of lines?
     time step control? implicity?
   - currently not memory optimized, let's experiment first
*/

#include "bam.h"
#include "evolve.h"


/* pointers to temporary variable lists */
#define NSTAGES 6
tVarList *u_k[NSTAGES];
#define PR 0

tEVO EVO;


/* various Butcher tables that provide the coefficients for Runge-Kutta */
typedef struct {
  int nstages;
  double a[NSTAGES][NSTAGES];
  double b[NSTAGES];
  double c[NSTAGES];
} tButcher;

/* RK1: Euler */
tButcher rk1 = {
  1,
  {{0}},
  {1}, 
  {1}
};

/* RK2: modified Euler (?) */
tButcher rk2 = {
  2,
  {{0}, {1}},
  {0.5, 0.5}, 
  {1, 1}
};

/* RK2: midpoint method */
tButcher rk2a = {
  2,
  {{0}, {0.5}},
  {0, 1  }, 
  {0, 0.5}
};


/* RK3: classic Heun, version 1, Mathematica's 3(2) */
tButcher rk3 = {
  3,
  {{0}, {0.5}, {-1, 2}},
  {1./6, 2./3, 1./6}, 
  {0,    1./2, 1   }
};

/* RK3: classic Heun, version 2 */
tButcher rk3a = {
  3,
  {{0}, {1./3}, {0, 2./3}},
  {1./4, 0   , 3./4}, 
  {0,    1./3, 2./3}
};

/* RK3: Matlab's ode23, Ralston */
tButcher rk3b = {
  3,
  {{0}, {1./2}, {0, 3./4}},
  {2./9, 1./3, 4./9}, 
  {0,    1./2, 3./4}
};

/* RK3: SSP */
tButcher rk3SSP = {
    3,
    {{0.}, {1.}, {0.25, 0.25}},
    {1./6., 1./6., 2./3.}, 
    {0.,    1., 0.5}
};


/* RK4: classic Runge-Kutta */
tButcher rk4 = {
  4,
  {{0}, {0.5}, {0, 0.5}, {0, 0, 1}},
  {1./6, 1./3, 1./3, 1./6}, 
  {0,    1./2, 1./2, 1   }
};
tButcher rk4g = {
  4,
  {{0}, {0.5}, {0, 0.5}, {0, 0, 1}},
  {1./6, 1./3, 1./3, 1./6}, 
  {0,    1./2, 1./2, 1   }
};


/* RK4: Mathematica's 4(3) */
tButcher rk4a = {
  4,
  {{0}, {2./5}, {-3./20, 3./4}, {19./44, -15./44, 10./11}},
  {11./72, 25./72, 25./72, 11./72}, 
  {2./5, 3./5, 1, 1}
};


/* RK4: SSP with 5 entries */
tButcher rk4SSP = {
    5,
    {{0.}, 
    {0.39175222700392}, 
    {0.21766909633821,0.36841059262959}, 
    {0.08269208670950,0.13995850206999,0.25189177424738},
    {0.06796628370320,0.11503469844438,0.20703489864929,0.54497475021237}},
    {0.14681187618661,0.24848290924556,0.10425883036650,0.27443890091960,0.22600748319395}, 
    {0.              ,0.39175222700392,0.58607968896779,0.47454236302687,0.93501063100924}
};





/* RK5: Dormand-Prince, ode45 in Matlab */
tButcher rk5 = {
  6,
  {{0}, 
   {1./5}, 
   {3./40, 9./40}, 
   {44./45, -56./15, 32./9},
   {19372./6561, -25360./2187, 64448./6561, -212./729},
   {9017./3168, -355./33, 46732./5247, 49./176, -5103./18656}
  },
  {35./384, 0, 500./1113, 125./192, -2187./6784, 11./84}, 
  {1./5, 3./10, 4./5, 8./9, 1, 1}
};





/* for debugging: print Butcher table */
void prbutchertable(tButcher *b)
{
  int i, j;

  printf("Runge-Kutta: Butcher table for %s\n", Gets("evolution_method_rk")); 
  for (i = 0; i < b->nstages; i++) {
    printf("%6.3f  ", b->c[i]);
    for (j = 0; j < b->nstages; j++) 
      printf("%6.3f ", b->a[i][j]);
    printf("\n");
  }
  printf("        ");
  for (j = 0; j < b->nstages; j++)
    printf("%6.3f ", b->b[j]);
  printf("\n");
}





/* generic explicit Runge-Kutta */
void evolve_rk_generic(tL *level, tButcher *rk, int comp) 
{
  double dt = level->dt;
  int i, j;

  /* initialize: u_p = u_c */
  if (PR) printf("start RK%d on level %d\n",rk->nstages,level->l);
  vlcopy(u_p, u_c);

  /* iterate over stages */
  for (i = 0; i < rk->nstages; i++) {
    if (PR) printf("      RK%d step %d\n",rk->nstages,i);
    if (i) vlcopy(u_c, u_p);
    EVO.stage = i;

    /* u_c = u_p + Sum aij * kj
       note that vladd ignores terms with zero coefficient */
    for (j = 0; j < i; j++)
      vladdto(u_c, dt*rk->a[i][j], u_k[j]); 
    /* u_ki = F(u_c) 
       step size 0.0 means right-hand side is returned in u_ki
       u_p is not needed but passed in as dummy
    */
    evolve_rhs(u_k[i], u_p, 0.0, u_c, comp);

    /* synchronize */
    bampi_vlsynchronize(u_k[i]);
  }

  /* compute final time level, u_c = u_p + Sum bj * kj */
  vlcopy(u_c, u_p);
  for (j = 0; j < rk->nstages; j++)
    vladdto(u_c, dt*rk->b[j], u_k[j]);

  if (PR) printf("stop  RK%d\n",rk->nstages);
}





/* rk4 written out for optimization (uses only 2 instead of 4 additional variables) */
void evolve_rk_rk4(tL *level, int comp)
{

  double dt = level->dt;
  tVarList *r = u_k[0], *v = u_k[1];

  if (PR) printf("start RK4 on level %d\n",level->l);
  vlcopy(u_p, u_c);                 // u_p = u_c

  if (PR) printf("      RK4 step 1\n");
  EVO.stage = 0;
  evolve_rhs(r, u_p, 0.0, u_c, comp);     // r    = f(u_p)
  bampi_vlsynchronize(r);
  vladdto(u_c, dt/6, r);            // u_c += dt/6 r
  
  if (PR) printf("      RK4 step 2\n");
  EVO.stage = 1;
  vladd(v, 1.0, u_p, dt/2, r);      // v    = u_p + dt/2 r
  evolve_rhs(r, u_p, 0.0, v, comp);       // r    = f(v)
  bampi_vlsynchronize(r);
  vladdto(u_c, dt/3, r);            // u_c += dt/3 r
  
  if (PR) printf("      RK4 step 3\n");
  EVO.stage = 2;
  vladd(v, 1.0, u_p, dt/2, r);      // v    = u_p + dt/2 r
  evolve_rhs(r, u_p, 0.0, v, comp);       // r    = f(v)
  bampi_vlsynchronize(r);
  vladdto(u_c, dt/3, r);            // u_c += dt/3 r
  
  if (PR) printf("      RK4 step 4\n");
  EVO.stage = 3;
  vladd(v, 1.0, u_p, dt  , r);      // v    = u_p + dt r
  evolve_rhs(r, u_p, 0.0, v, comp);       // r    = f(v)
  bampi_vlsynchronize(r);
  vladdto(u_c, dt/6, r);            // u_c += dt/6 r
  
  if (PR) printf("stop  RK4\n");
}







/* explicit Runge-Kutta */
void evolve_rk(tL *level, char *name, int comp) 
{
  static int firstcall = 1;
  tButcher *rk;
  char s[3] = "_r";
  int i;


  /* pick method */
  if      (strcmp(name,"rk5"     )==0)   rk = &rk5;
  else if (strcmp(name,"rk4SSP"  )==0)   rk = &rk4SSP;
  else if (strcmp(name,"rk4a"    )==0)   rk = &rk4a;
  else if (strcmp(name,"rk4g"    )==0)   rk = &rk4g;
  else if (strcmp(name,"rk4"     )==0)   rk = &rk4;
  else if (strcmp(name,"rk3SSP"  )==0)   rk = &rk3SSP;
  else if (strcmp(name,"rk3b"    )==0)   rk = &rk3b;
  else if (strcmp(name,"rk3a"    )==0)   rk = &rk3a;
  else if (strcmp(name,"rk3"     )==0)   rk = &rk3;
  else if (strcmp(name,"rk2a"    )==0)   rk = &rk2a;
  else if (strcmp(name,"rk2"     )==0)   rk = &rk2;
  else if (strcmp(name,"rk1"     )==0)   rk = &rk1;
  else if (strcmp(name,"midpoint")==0)   rk = &rk2a;
  else if (strcmp(name,"euler"   )==0)   rk = &rk1;
  else
    errorexits("unknown radiation_evolution_method_rk %s", Gets("rad_evolution_method_rk"));

  /* catch special cases
     some methods are implemented without tables for efficiency */
  if (rk == &rk4)
    rk->nstages = 2;

  /* add auxiliary fields to data base, no storage yet 
     incrementing second letter in s gives _r, _s, _t, ...
  */
  if (firstcall) {
    firstcall = 0;
    for (i = 0; i < rk->nstages; i++, s[1]++) 
      u_k[i] = AddDuplicate(u_c, s);

    if (0) prbutchertable(rk);
  }

  /* store current level in existing variable lists */
  for (i = 0; i < rk->nstages; i++)
    if (u_k[i]) u_k[i]->level = level;

  /* turn on memory (if varlist is non null) */
  for (i = 0; i < rk->nstages; i++)
    enablevarlist(u_k[i]);

  /* evolve */
  if (rk == &rk4)
    evolve_rk_rk4(level, comp);
  else
    evolve_rk_generic(level, rk, comp);

  /* turn off memory */
  for (i = 0; i < rk->nstages; i++)
    disablevarlist(u_k[i]);
}




/* return order of Runge-Kutta (needed for amr in-time interpolation) */
int evolve_rk_norder(void)
{
  static int firstcall = 1;
  static int order;

  if (firstcall) {
    firstcall = 0;

    /* set order */
    if      (Getv("evolution_method_rk", "rk5"))      order = 5;
    else if (Getv("evolution_method_rk", "rk4SSP"))   order = 4;
    else if (Getv("evolution_method_rk", "rk4a"))     order = 4;
    else if (Getv("evolution_method_rk", "rk4"))      order = 4;
    else if (Getv("evolution_method_rk", "rk3SSP"))   order = 3;
    else if (Getv("evolution_method_rk", "rk3b"))     order = 3;
    else if (Getv("evolution_method_rk", "rk3a"))     order = 3;
    else if (Getv("evolution_method_rk", "rk3"))      order = 3;
    else if (Getv("evolution_method_rk", "rk2a"))     order = 2;
    else if (Getv("evolution_method_rk", "rk2"))      order = 2;
    else if (Getv("evolution_method_rk", "rk1"))      order = 1;
    else if (Getv("evolution_method_rk", "midpoint")) order = 2;
    else if (Getv("evolution_method_rk", "euler"))    order = 1;
    else
      errorexits("unknown evolution_method_rk %s",
		 Gets("evolution_method_rk"));
  }
  
  return order;
}


/* iterative Crank-Nicholson */
void evolve_rad_imex(tL *level)
{
  double dt = level->dt;
  int i;

  EVO.stage = 0;

  /* store current level in existing variable lists */
  for (i = 0; i < 2; i++)
    if (u_k[i]) u_k[i]->level = level;

  /* turn on memory */
  for (i = 0; i < 2; i++)
    enablevarlist(u_k[i]);

  evolve_radiation(u_c, dt, u_k[0], u_k[1]);

  /* turn off memory */
  for (i = 0; i < 2; i++)
    disablevarlist(u_k[i]);

}
