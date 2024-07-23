/* storage.c */
/* Bernd Bruegmann, 12/99 */
/* experiment with most primitive AMR: tree of single nodes */ 

#include "bam.h"
#include "amr.h"

#define PR 0

int storage_verbose = 0;



/**************************************************************************/
/* basic memory management */



/* allocate refinement level */
tL *alloc_level(tG *g, int l, int nnodes) 
{
  tL *new;

  /* make new level structure */
  new = (tL *) calloc(1, sizeof(tL));
  if (!new) errorexit("alloc_level:  out of memory");
  new->l = l;
  new->shells = 0;

  /* put new level into grid */
  /* assumes caller protects existing levels */
  replace_level(g, new, l);
  
  /* allocate storage for nodes, set node->i and node->l */
  new->nnodes = nnodes;

  /* allocate storage for data pointers, they default to NULL */
  realloc_levelvariables(new, g->nvariables);

  /* allocate communication data */
  new->com = bampi_alloc_com();

  //  /* miscellaneous */
  // new->box = (tBox *) calloc(3, sizeof(tBox));
  // if (!new->box) errorexit("alloc_level:  out of memory");

  return new;
} 




/* several utility functions to enable and disable storage for variables */

/* count (currently nenabled is too large, apparently misses some frees) */
int nenabled = 0;
int nenabledmax = 0;

void printenabled(tL *level)
{
  int n = level->nvariables;
  int i, j;
  
  for (i = j = 0; i < n; i++) 
    if (level->v[i]) j++;

  nenabled = j;
  if (j > nenabledmax) nenabledmax = j;

  printf("%d variables, %d off, %d ON, %d ON max\n", n, n-j, j, nenabledmax);
}




/* enable one component of variable 
   storage is initialized to zero
*/
void enablevarcomp(const tL *level, const int i)
{
  if (!level->v[i]) {
    level->v[i] = (double*) malloc( (level->nnodes) * sizeof(double));
    if (!level->v[i]) {
      printf("enable storage for variable %d = %s: out of memory\n",
	     i, VarName(i));
      errorexit("");
    }

    bampi_openmp_parallel_for
    for (int j = 0; j < level->nnodes; j++) {
        level->v[i][j] = 0.0;
    }

    nenabled++;
    //if (nenabledmax < nenabled) nenabledmax = nenabled;

    if (storage_verbose) {
      printf("enabling  variable %2d = %-20s on l%d => ",
	     i, VarName(i), level->l);
      printf("%10d * %2d * %2d = %10d bytes  (%p)\n", 
	     level->nnodes, nenabled, (int)sizeof(double), 
	     level->nnodes * nenabled * (int)sizeof(double),level->v[i]);
    } 
  }
}




/* disable one component of variable */
void disablevarcomp(const tL *level, const int i)
{
  if (level->v[i]) {
    if (storage_verbose) {
      printf("disabling variable %2d = %-20s on l%d => ", i, VarName(i), level->l);
      printf("%10d * %2d * %2d = %10d bytes  (%p)\n", 
	     level->nnodes, nenabled, (int)sizeof(double), 
             level->nnodes * nenabled * (int)sizeof(double),level->v[i]);
    }
    free(level->v[i]);
    level->v[i] = 0;
    nenabled--;
  }
}




/* enable all components of a variable */
void enablevar(tL *level, int i) 
{
  int j, n = VarNComponents(i);

  for (j = 0; j < n; j++)
    enablevarcomp(level, i+j);
}




/* disable all components of a variable */
void disablevar(const tL *level, const int i)
{
  int j, n = VarNComponents(i);

  for (j = 0; j < n; j++)
    disablevarcomp(level, i+j);
}




/* enable variable list */
void enablevarlist(const tVarList *v)
{
  if (!v) {
    return;
  }
  int i;
    for (i = 0; i < v->n; i++) 
      enablevarcomp(v->level, v->index[i]);
  }




/* disable variable list */
void disablevarlist(const tVarList *v)
{
  if (!v) {
    return;
  }
  int i;

    for (i = 0; i < v->n; i++)
      disablevarcomp(v->level, v->index[i]);
}




/* enable all variables found on given level */
void enablesamevars(tL *level, tL *newlevel)
{
  int i;

  if (level->nvariables != newlevel->nvariables) {
    printf("nvariables: old %d, new %d\n", 
	   level->nvariables, newlevel->nvariables);
    errorexit("enablesamevars: need same number of variables");
  }

  for (i = 0; i < level->nvariables; i++)
    if (level->v[i]) 
      enablevarcomp(newlevel, i);
}


 

/* find one component of variable by name and return pointer */
double *Ptr(const tL *level, char *name)
{
  return level->v[Ind(name)];
}




/* enable one component of variable by name and return pointer */
double *PtrEnable(const tL *level, const char *name)
{
  int i = Ind(name);

  enablevarcomp(level, i);
  return level->v[i];
}




/* disable one component of variable by name */
void PtrDisable(const tL *level, char *name)
{
  disablevarcomp(level, Ind(name));
}




/* create room for more variables */
void realloc_levelvariables(tL *level, int nvariables)
{
  int i;

  if (PR) printf("realloc_levelvariables from %d to %d\n", 
		 level->nvariables, nvariables);
  //level->v = (double **) realloc(level->v, sizeof(double *) * nvariables);
  double **tmp = (double **) realloc(level->v, sizeof(double *) * nvariables);
  if (!tmp) {
    tmp = (double **) malloc(sizeof(double *) * nvariables);
    free(level->v);
  }
  level->v = tmp;
  for (i = level->nvariables; i < nvariables; i++)
    level->v[i] = 0;
  level->nvariables = nvariables;
  
  /* should actually check whether this is correct */
  level->grid->nvariables = nvariables;
}




/* free level, leaves grid untouched */
void free_level(tL *level) 
{
  int i;
  tL *prglobal, *prlocal;

  if (!level) return;

  /* free ghost parent levels (but not the regular parent) */
  if (level->prlocal && 
      level->prlocal != level->grid->level[level->l-1] &&
      (bampi_size() > 1 || level->nboxes > 1)
     ) {
    free_level(level->prlocal);
  }

  /* free communication data */
  if (level->nboxes > 1) {
    for (i = 0; i < level->nboxes; i++)
      bampi_free_com(level->box[i]->com);
  } else 
    bampi_free_com(level->com);

  /* free boxes */
  for (i = 0; i < level->nboxes; i++) {
    free(level->box[i]);
  }
  free(level->box);
  free(level->sbox);
  
  /* free variable storage */
  if (0) 
    printf("nvariables %d %d\n", level->nvariables, level->grid->nvariables);
  for (i = 0; i < level->nvariables; i++)
    disablevarcomp(level, i);
  free(level->v);

  /* free node storage */
  free(level->boundary);

  /* and finally free the level structure itself */
  free(level);
}




/* remove level and all sublevels: free levels and adjust grid */
void remove_levels(tL *level)
{
  tG *grid = level->grid;
  int newlmax = level->l - 1;
  int l;

  /* sanity checks */
  if (!level) return;
  if (level->l < 1) return;

  /* free the levels */
  for (l = level->l; l <= grid->lmax; l++) 
    free_level(grid->level[l]);

  /* set new number of levels */
  grid->lmax = newlmax;
  grid->nlevels = newlmax + 1;
  Seti("amr_lmax", newlmax);
}




/* allocate grid */
tG *alloc_grid(int lmin, int lmax, int nvariables)
{
  tG *new;

  new = (tG *) calloc(1, sizeof(tG));
  if (!new) errorexit("alloc_grid:  out of memory");

  new->lmin = lmin;
  new->lmax = lmax;
  new->ltop = lmin;
  new->lbo = (Getv("amr", "bo")) ? 0 : INT_MAX;
  new->nvariables = nvariables;
  new->nlevels = lmax - lmin + 1;
  new->level = (tL **) calloc(new->nlevels, sizeof(tL *));
  if (!new->level) errorexit("alloc_grid:  out of memory 2");

  storage_verbose = Getv("storage_verbose", "yes");
  return new;
}




/* free grid without freeing levels */
void free_grid_only(tG *g)
{
  if (!g) return;
  if (g->level) free(g->level);
  free(g);
}




/* free grid */
void free_grid(tG *g)
{
  int i;

  if (!g) return;
  for (i = g->nlevels-1; i >= 0; i--) {
    if (PR) printf("free level %d\n", i);
    free_level(g->level[i]);
  }
  free_grid_only(g);
}




/* CLEAN UP: replace level and insert level can be handled by one function */
/* replace level: 
   overwrite level in grid assuming caller protects existing levels!
*/
void replace_level(tG *g, tL *level, int l)
{
  if (!g) {
    printf("error: grid not defined in replace_level\n");
    return;
  }
  if (!level) {
    printf("error: level not defined in replace_level\n");
    return;
  }

  /* allowed? */
  if (l < g->lmin || l > g->lmax+1)
    errorexit("replace_level: index out of range");

  /* if level does not exist yet, make room for new level */
  if (l > g->lmax) {
    g->level = (tL **) realloc(g->level, sizeof(tL *) * (l+1));
    if (!g->level) errorexit("replace_level:  out of memory");

    g->lmax = l;
    g->nlevels = g->lmax - g->lmin + 1;
  }

  /* put level into grid */
  g->level[l] = level;
  level->grid = g;
  level->l = l;
}




/* insert level into grid at specified index
   unfinished: also supports negative levels, i.e.
     m = 0:           largest physical level
     m =  1,  2, ...: physical refinements
     m = -1, -2, ...: multigrid coarsenings
   with suitable changes if the largest physical level is some m != 0
*/
void insert_level(tG *g, tL *l, int m)
{
  int i;

  if (!g) return;
  if (!l) return;

  /* insert allowed? */
  if (m < g->lmin-1 || m > g->lmax+1)
    errorexit("insert_level: index out of range, no gaps allowed");

  /* set info in level that refers to the grid */
  l->grid = g;
  l->l = m;
  if (0 && l->nvariables != g->nvariables)
    //errorexit("insert_level: wrong number of variables");
    printf("warning in insert_level: wrong number of variables\n");

  /* make room for level in grid structure */
  /* currently we only support push_level(), see below, with g->lmin = 0 */
  if (g->lmin)
    errorexit("insert_level(): implement g->lmin != 0"); 
  g->nlevels++;
  g->level = (tL **) realloc(g->level, g->nlevels * sizeof(tL *));
  if (!g->level)
    errorexit("insert_level: out of memory");

  /* append */
  if (m == g->lmax+1) {
    g->lmax++;
    g->ltop++;
    g->level[g->lmax] = l;
  }
    
  /* insert at 0, relabelling all other levels */
  else if (m == 0) {
    g->lmax++;
    g->ltop++;
    for (i = g->nlevels-1; i >= 1; i--) {
      g->level[i] = g->level[i-1];
      g->level[i]->l += 1;
    }
    for (i = 0; i < g->nlevels; i++)
      if (g->level[i]->prlocal) 
	g->level[i]->prlocal->l = g->level[i]->l - 1;

    g->level[g->lmin] = l;
  }

  /* default */
  else
    errorexit("insert_level: illegal operation");
}




/* standard wrappers for insert_level 
   currently, push and prepend cannot be mixed
   actually, append and prepend are not finished in insert_level()
*/
void push_level(tG *g, tL *l)
{
  insert_level(g, l, 0);
}



/* remove top level */
void remove_top_level(tG *g)
{
  int i;

  if (!g) return;

  if (g->level[0]->l == 0) {
    free_level(g->level[0]);
    g->nlevels--;
    g->lmax--;
    g->ltop--;
    for (i = 0; i < g->nlevels; i++) {
      g->level[i] = g->level[i+1];
      g->level[i]->l -= 1;
    }
    for (i = 0; i < g->nlevels; i++)
      if (g->level[i]->prlocal) 
	g->level[i]->prlocal->l = g->level[i]->l - 1;
  }

  else 
    errorexit("remove_top_level: illegal operation");
} 

