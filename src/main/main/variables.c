/* variables.c */
/* Bernd Bruegmann 12/99 */

#include "bam.h"

typedef struct
{
  double iotime[4];
  int ioiter[4];
  int ioflag[4];
  int ioatall;
} tIO;

typedef struct tVAR
{
  char *name;
  char *tensorindices;
  char *description;
  int index;
  int ncomponents;
  int component;
  tIO *io;
  double farlimit;
  double falloff;
  double propspeed;
  int gauge;
  int sym[3];
  int constant;
} tVar;

tVar *vdb = 0;
int nvdb = 0;
int globalnvariables = 0;

tVarList **vlldb = 0;
int nvlldb = 0;

/* add variable to data base */
void AddVar(const char *name, char *tensorindices, const char *description)
{
  tVar *new;
  int i, j;
  char fullname[100];
  int nilist;
  char *ilist[NINDEXLIST];
  int sym[3 * NINDEXLIST];
  char *symsigns[3] = {"-", "0", "+"}, **ss = symsigns + 1;

  if (0)
    printf("AddVar: name %s, tensorindices %s\n", name, tensorindices);

  /* construnct list of tensor indices */
  tensorindexlist(tensorindices, &nilist, ilist, sym);

  /* for each tensor index */
  for (j = 0; j < nilist; j++)
  {

    /* construct name of variable */
    snprintf(fullname, 100, "%s%s", name, ilist[j]);
    if (0)
      printf("  variable  %s\n", fullname);
    else
    {
      printf("  variable  %-32s", fullname);
      printf("sym %s%s%s\n", ss[sym[3 * j]], ss[sym[3 * j + 1]], ss[sym[3 * j + 2]]);
    }

    /* make sure that this variable does not exist yet */
    for (i = 0; i < nvdb; i++)
      if (!strcmp(vdb[i].name, fullname))
        errorexits("AddVar: variable \"%s\" already exists\n", fullname);

    /* variable does not exist, so add a new element to data base */
    vdb = (tVar *)realloc(vdb, sizeof(tVar) * (nvdb + 1));
    new = &vdb[nvdb];

    /* initialize and fill in structure */
    memset(new, 0, sizeof(tVar));
    new->name = strdup(fullname);
    new->tensorindices = strdup(tensorindices);
    new->description = strdup(description);
    new->index = nvdb;
    new->ncomponents = nilist;
    new->component = j;
    new->io = NULL;
    new->sym[0] = sym[3 * j];
    new->sym[1] = sym[3 * j + 1];
    new->sym[2] = sym[3 * j + 2];

    nvdb++;
    globalnvariables = nvdb;

    free(ilist[j]);
  }
}

/* add constant variable to data base */
void AddConstantVar(const char *name, char *tensorindices,
                    const char *description)
{
  AddVar(name, tensorindices, description);
  VarNameSetConstantFlag(vdb[nvdb - 1].name);
}

/* free all variable storage */
void freeVars()
{
  int i;
  for (i = 0; i < nvdb; i++)
  {
    if (0)
      printf("delete %d %s\n", i, vdb[i].name);
    free(vdb[i].name);
    free(vdb[i].tensorindices);
    free(vdb[i].description);
    free(vdb[i].io);
  }
  if (0)
    printf("%d\n", nvdb);
  free(vdb);
}

/* return index of variable */
int Ind(const char *name)
{
  int i;

  for (i = 0; i < nvdb; i++)
    if (!strcmp(vdb[i].name, name))
    {
      if (0)
        printf("index(%s) = %d\n", name, vdb[i].index);
      return vdb[i].index;
    }
  errorexits("Ind: variable \"%s\" does not exist\n", name);
  return 0;
}

/* return index of variable or -1 if it was not found */
int IndLax(const char *name)
{
  int i;

  for (i = 0; i < nvdb; i++)
    if (!strcmp(vdb[i].name, name))
    {
      if (0)
        printf("index(%s) = %d\n", name, vdb[i].index);
      return vdb[i].index;
    }
  return -1;
}

int Indvl(const char *name, const tVarList *vl)
{
  int i;

  for (i = 0; i < vl->n; i++)
    if (strncmp(VarName(vl->index[i]), name, strlen(name)) == 0)
      return vl->index[i];
  errorexits("Ind: variable \"%s\" does not exist inside varlist\n", name);
  return 0;
}
/* return index of variable given pointer */
int IndFromPtr(const tL *level, const double *p)
{
  int i;

  for (i = 0; i < nvdb; i++)
    if (level->v[i] == p)
      return vdb[i].index;
  return -1;
}

/* return name given index */
char *VarName(const int i)
{
  if (i < 0 || i >= nvdb)
    errorexit("VarName: index out of range");

  return vdb[i].name;
}

/* return number of components */
int VarNComponents(const int i)
{
  if (i < 0 || i >= nvdb)
    errorexit("VarNComponents: index out of range");
  if (vdb[i].component != 0)
    errorexit("VarNComponents: you have to use index of zeroth component");

  return vdb[i].ncomponents;
}

/* return component */
int VarComponent(const int i)
{
  return vdb[i].component;
}

/* return index of component 0 */
int IndComponent0(const int i)
{
  return i - vdb[i].component;
}

/* return name of component 0 for a given name */
char *VarNameComponent0(const char *name)
{
  return VarName(IndComponent0(Ind(name)));
}

/* return string with tensor indices */
char *VarTensorIndices(const int i)
{
  return vdb[i].tensorindices;
}

/* set information on how variable behaves at Boundary*/
void VarNameSetBoundaryInfo(const char *name, const double farlimit,
                            const double falloff, const double propspeed)
{
  int i = Ind(name);

  vdb[i].farlimit = farlimit;
  vdb[i].falloff = falloff;
  vdb[i].propspeed = propspeed;
}

/* set information on how variable behaves at Boundary*/
void VarNameSetConstantFlag(const char *name)
{
  int i, i0 = IndComponent0(Ind(name));
  int n = VarNComponents(i0);

  for (i = 0; i < n; i++)
  {
    vdb[i + i0].constant = 1;
    if (0)
      printf("  setting %s constant\n", vdb[i + i0].name);
  }
}

/* return various boundary information */
double VarFallOff(const int i) { return vdb[i].falloff; }
double VarFarLimit(const int i) { return vdb[i].farlimit; }
double VarPropSpeed(const int i) { return vdb[i].propspeed; }
int VarGaugeFlag(const int i) { return vdb[i].gauge; }
int VarConstantFlag(const int i) { return vdb[i].constant; }
int VarSymmetry(const int i, const int dir) { return vdb[i].sym[dir]; }

/* set symmetry information */
void VarSetSymmetry(const int i, const int dir, const int sym)
{
  vdb[i].sym[dir] = sym;
}
void VarInvertSymmetry(const int i)
{
  vdb[i].sym[0] *= -1;
  vdb[i].sym[1] *= -1;
  vdb[i].sym[2] *= -1;
}

/************************************************************************/
/* utility functions for variable lists */

/* print variable list */
void prvarlist(const tVarList *v)
{
  int i, j;

  printf("VarList:  n = %d    v = %p\n", v->n, v);
  for (i = 0; i < v->n; i++)
  {
    j = v->index[i];
    if (v->level)
      printf(" %2d i=%3d level=%p l=%d v=%p %s\n",
             i, j, v->level, v->level->l, v->level->v[j],
             VarName(v->index[i]));
    else
      printf(" %2d i=%3d level=%p %s\n", i, j, v->level,
             VarName(v->index[i]));
  }
}

/* allocate an empty variable list */
tVarList *vlalloc(tL *level)
{
  timer_start(level, "vlalloc");
  tVarList *u;

  u = calloc(1, sizeof(tVarList));
  u->level = level;

  vlldb = (tVarList **)realloc(vlldb, sizeof(tVarList *) * (nvlldb + 1));
  vlldb[nvlldb] = u;
  nvlldb++;

  timer_stop(level, "vlalloc");
  return u;
}

/* free a variable list */
void vlfree(tVarList *u)
{
  int i, j;

  if (u)
  {

    for (i = nvlldb - 1; i >= 0; i--)
      if (u == vlldb[i])
        break;

    if (i >= 0)
    {
      for (j = i + 1; j < nvlldb; j++)
      {
        vlldb[j - 1] = vlldb[j];
      }
      nvlldb--;
    }
    else
      printf("Warning: trying to remove tVarList that is not in the list.\n");

    if (u->index)
      free(u->index);
    free(u);
  }
}

void freevll()
{
  while (nvlldb)
    vlfree(vlldb[0]);
  free(vlldb);
}

/* add a variable (one component) to a variable list,
   but only if it does not exist yet
*/
void vlpushone(tVarList *v, const int vi)
{
  v->n++;
  v->index = realloc(v->index, sizeof(int) * v->n);
  v->index[v->n - 1] = vi;
}

/* add a variable with all its components to a variable list */
void vlpush(tVarList *v, const int vi)
{
  int i, n = VarNComponents(vi);

  for (i = 0; i < n; i++)
    vlpushone(v, vi + i);
}

/* add a variable list to a variable list */
void vlpushvl(tVarList *v, const tVarList *u)
{
  int i;

  if (!v || !u)
    return;

  for (i = 0; i < u->n; i++)
    vlpushone(v, u->index[i]);
}

/* drop a variable (one component) from a variable list */
void vldropone(tVarList *v, const int vi)
{
  int i;

  for (i = 0; i < v->n; i++)
    if (v->index[i] == vi)
    {
      v->n -= 1;
      for (; i < v->n; i++)
        v->index[i] = v->index[i + 1];
      v->index = realloc(v->index, sizeof(int) * v->n);
      break;
    }
}

/* drop a variable with all its components from a variable list */
void vldrop(tVarList *v, const int vi)
{
  int i, n = VarNComponents(vi);

  for (i = 0; i < n; i++)
    vldropone(v, vi + i);
}

/* drop last n variables from a variable list */
void vldropn(tVarList *v, const int n)
{
  if (n <= 0)
    return;
  if (n >= v->n)
    v->n = 0;
  else
    v->n -= n;
}

/* return new variable list containing the enabled variables of a given list */
tVarList *vlclean(const tVarList *v)
{
  int i;
  tVarList *u = vlalloc(v->level);

  for (i = 0; i < v->n; i++)
    if (v->level->v[v->index[i]])
      vlpushone(u, v->index[i]);

  return u;
}

/* duplicate variable list */
tVarList *vlduplicate(const tVarList *v)
{
  int i;
  tVarList *u = vlalloc(v->level);

  for (i = 0; i < v->n; i++)
    if (VarComponent(v->index[i]) == 0)
      vlpush(u, v->index[i]);

  return u;
}

/* enable all variables in a variable list */
void vlenable(const tVarList *v)
{
  enablevarlist(v);
}

void vlenablelevel(tL *level, tVarList *v)
{
  v->level = level;
  enablevarlist(v);
}

/* disable all variables in a variable list */
void vldisable(const tVarList *v)
{
  disablevarlist(v);
}

/* create, enable, return pointer for a 1 variable VarList */
tVarList *VLPtrEnable1(tL *level, const char *varname)
{
  tVarList *vl = vlalloc(level);
  int i = Ind(varname);

  enablevar(level, i);
  vlpush(vl, i);
  return vl;
}

/* disable variables in a VarList and free VarList */
void VLDisableFree(tVarList *vl)
{
  int i;

  disablevarlist(vl);
  vlfree(vl);
}

/* add variables based on an existing variable list and a postfix
   note that we add each component as a scalar but fix it later because
   we want gxx_p, gxy_p, ...  and not gxx_pxx, gxx_pxy ...
*/
tVarList *AddDuplicate(const tVarList *vl, const char *postfix)
{
  char name[1000];
  int i, j;
  int nadded = 0;
  tVarList *newvl;
  tVar *var, *newvar;

  /* new variable list with same number of indices */
  newvl = vlduplicate(vl);

  /* for all scalar variables in list */
  for (i = 0; i < vl->n; i++)
  {

    /* construct new name */
    var = &vdb[vl->index[i]];
    snprintf(name, 1000, "%s%s", var->name, postfix);

    /* if variable already exists, don't add it again */
    /* note that we nevertheless return a corresponding variable list */
    if ((j = IndLax(name)) >= 0)
    {
      newvl->index[i] = j;
      continue;
    }

    /* add scalar variable with new name to variable database */
    AddVar(name, "", var->description);
    nadded++;

    /* get index of new variable and overwrite index in duplicate */
    newvl->index[i] = Ind(name);

    /* get pointer to old variable again since AddVar reallocates vdb */
    var = &vdb[vl->index[i]];

    /* set structure in variable data base */
    newvar = &vdb[newvl->index[i]];
    free(newvar->tensorindices);
    newvar->tensorindices = strdup(var->tensorindices);
    newvar->component = var->component;
    newvar->ncomponents = var->ncomponents;
    newvar->farlimit = var->farlimit;
    newvar->falloff = var->falloff;
    newvar->propspeed = var->propspeed;
    newvar->gauge = var->gauge;
    newvar->constant = var->constant;
    for (j = 0; j < 3; j++)
      newvar->sym[j] = var->sym[j];
  }

  /* create storage for as many variables as have been actually added
     do it on all levels so that nvariables remains the same on all levels
  */
  if (vl->level && nadded)
  {
    tG *g = vl->level->grid;
    int n = vl->level->nvariables + nadded;
    int l;

    for (l = g->lmin; l <= g->lmax; l++)
      realloc_levelvariables(g->level[l], n);
  }
  if (0)
    printf("nvdb is now %d\n", nvdb);

  return newvl;
}

/* add duplicate and enable variables */
tVarList *AddDuplicateEnable(const tVarList *vl, const char *postfix)
{
  tVarList *newvl;

  newvl = AddDuplicate(vl, postfix);
  enablevarlist(newvl);
  return newvl;
}

/* set: u = c */
void vlsetconstant(const tVarList *u, const double c)
{
  tL *level = u->level;
  timer_start(level, "vlsetconstant");
  int nnodes = level->nnodes;

  bampi_openmp_parallel_for_collapse2
  for (int n = 0; n < u->n; n++)
  {
    for (int i = 0; i < nnodes; i++)
    {
      double *pu = level->v[u->index[n]];
      pu[i] = c;
    }
  }
  timer_stop(level, "vlsetconstant");
}

/* copy: v = u */
void vlcopy(const tVarList *v, const tVarList *u)
{
  tL *level = v->level;
  timer_start(level, "vlcopy");
  int nnodes = level->nnodes;

  bampi_openmp_parallel_for_collapse2
  for (int n = 0; n < v->n; n++)
  {
    for (int i = 0; i < nnodes; i++)
    {
      double *pu = level->v[u->index[n]];
      double *pv = level->v[v->index[n]];
      pv[i] = pu[i];
    }
  }
  timer_stop(level, "vlcopy");
}

void vlcopylevel(tL *level, tVarList *v, tVarList *u)
{
  if (!level || !v || !u)
    return;
  v->level = u->level = level;
  vlcopy(v, u);
}

/* average: r=(a+b)/2 */
void vlaverage(const tVarList *r, const tVarList *a, const tVarList *b)
{
  tL *level = r->level;
  timer_start(level, "vlaverage");
  int nnodes = level->nnodes;

  bampi_openmp_parallel_for_collapse2
  for (int n = 0; n < r->n; n++)
  {
    for (int i = 0; i < nnodes; i++)
    {
      double *pr = level->v[r->index[n]];
      double *pa = level->v[a->index[n]];
      double *pb = level->v[b->index[n]];
      pr[i] = 0.5 * (pa[i] + pb[i]);
    }
  }
  timer_stop(level, "vlaverage");
}

/* subtract two var lists: r = a - b
   can be called as vlsubtract(r,a,b); or vlsubtract(a,a,b); */
void vlsubtract(const tVarList *r, const tVarList *a, const tVarList *b)
{
  tL *level = r->level;
  timer_start(level, "vlsubtract");
  int nnodes = level->nnodes;

  bampi_openmp_parallel_for_collapse2
  for (int n = 0; n < r->n; n++)
  {
    for (int i = 0; i < nnodes; i++)
    {
      double *pr = level->v[r->index[n]];
      double *pa = level->v[a->index[n]];
      double *pb = level->v[b->index[n]];
      pr[i] = pa[i] - pb[i];
    }
  }
  timer_stop(level, "vlsubtract");
}

/* linear combination of two var lists: r = ca*a + cb*b
   should change function name
   one function can catch several special cases like cb == 0 (unfinished)
   important: if coefficient is zero we guarantee that memory is not accessed
*/
void vladd(const tVarList *r, const double ca, const tVarList *a,
           const double cb, const tVarList *b)
{
  tL *level = r->level;
  timer_start(level, "vladd");
  int nnodes = level->nnodes;

  if (ca == 0 && cb == 0)
  {
    bampi_openmp_parallel_for_collapse2
    for (int n = 0; n < r->n; n++)
    {
      for (int i = 0; i < nnodes; i++)
      {
        double *pr = level->v[r->index[n]];
        pr[i] = 0;
      }
    }
    timer_stop(level, "vladd");
    return;
  }

  if (ca == 0)
  {
    bampi_openmp_parallel_for_collapse2
    for (int n = 0; n < r->n; n++)
    {
      for (int i = 0; i < nnodes; i++)
      {
        double *pr = level->v[r->index[n]];
        double *pb = level->v[b->index[n]];
        pr[i] = cb * pb[i];
      }
    }
    timer_stop(level, "vladd");
    return;
  }

  if (cb == 0)
  {
    bampi_openmp_parallel_for_collapse2
    for (int n = 0; n < r->n; n++)
    {
      for (int i = 0; i < nnodes; i++)
      {
        double *pr = level->v[r->index[n]];
        double *pa = level->v[a->index[n]];
        pr[i] = ca * pa[i];
      }
    }
    timer_stop(level, "vladd");
    return;
  }

  bampi_openmp_parallel_for_collapse2
  for (int n = 0; n < r->n; n++)
  {
    for (int i = 0; i < nnodes; i++)
    {
      double *pr = level->v[r->index[n]];
      double *pa = level->v[a->index[n]];
      double *pb = level->v[b->index[n]];
      pr[i] = ca * pa[i] + cb * pb[i];
    }
  }
  timer_stop(level, "vladd");
}

/* add second var list to first: r += ca*a
   special treatment for ca = 1 and ca = -1
*/
void vladdto(const tVarList *r, const double ca, const tVarList *a)
{
  tL *level = r->level;
  timer_start(level, "vladdto");
  int nnodes = level->nnodes;

  if (ca == 0)
    return;

  bampi_openmp_parallel_for_collapse2
  for (int n = 0; n < r->n; n++)
  {
    for (int i = 0; i < nnodes; i++)
    {
      double *pr = level->v[r->index[n]];
      double *pa = level->v[a->index[n]];
      pr[i] += ca * pa[i];
    }
  }
  timer_stop(level, "vladdto");
}

/* swap content of variables in two var lists
   actually swaps pointers, make sure that this works in your case
*/
void vlswap(const tVarList *u, const tVarList *v)
{
  tL *level = v->level;

  for (int n = 0; n < v->n; n++)
  {
    double *pu = level->v[u->index[n]];
    double *pv = level->v[v->index[n]];
    level->v[u->index[n]] = pv;
    level->v[v->index[n]] = pu;
  }
}