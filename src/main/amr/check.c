/* CheckIfFinite.c */
/* Wolfgang Tichy 3/2003 */
/* check if ham contains NANs or INFs */

#include "bam.h"
#include "amr.h"

#define PR 0


/* check if var with index ivar is finite 
   and return the number of NANs found      */
int CheckIfFinite_VarIndexOne(tL* level, int ivar)
{
  double *var;
  double num=0.0;
  int ccc, ccc_old;
  double *x; 
  double *y;  
  double *z;  
  int messageflag=0;

  /* printf("Checking for INF or NAN in %s index %d\n", VarName(ivar), ivar); */
  /* var = level->v[ivar]; 
     printf("level=%p var=%p\n", level, var); */

  /* The stuff inside the brackets below is only needed because 
     errorexit(""); does not work properly!                          */
  /* alternative method that does not find out position of NAN
     incidentally, it makes sure every processor knows about it */
  {
    double *sum;
    tVarList *vl = vlalloc(level);

    enablevarcomp(level, ivar);  // This will enable the var!!
    vlpushone(vl, ivar);

    bampi_allreduce_sum(vl, &sum);

    if(!finite(*sum)) num=0.1;

    free(sum);
    vlfree(vl);
    /* Currently (5/30/2016) num is always 0, even if there are NANs.
       So skip next line and go on below and check everywhere */
    //if(num==0.0) return 0;
  }

  /* printf("Checking for INF or NAN in %s\n", VarName(ivar)); */
  var = level->v[ivar]; 
  
  forallpoints(level,ccc)
  {
    if(var==NULL) 
    {
      x=level->v[Ind("x")];
      y=level->v[Ind("y")];
      z=level->v[Ind("z")];
      printf("pointer to %s is NULL at ccc=%d:  x=%f y=%f z=%f\n",
	     VarName(ivar),ccc,x[ccc],y[ccc],z[ccc]); 
      continue;
    }
      
    if( !finite(var[ccc]) ) 
    {
      if(messageflag==0)
      {
	x=level->v[Ind("x")];
	y=level->v[Ind("y")];
	z=level->v[Ind("z")];
	
	printf("NAN/INF: %s=%g at ccc=%d:  x=%g y=%g z=%g\n", 
	       VarName(ivar), var[ccc], ccc, x[ccc], y[ccc], z[ccc]);
      }
      ccc_old=ccc;
      messageflag++;
      num++;
    }
    else
    {
      if(messageflag>1)
        printf("NAN/INF: %s=%g til    %d:  x=%g y=%g z=%g\n",
               VarName(ivar), var[ccc_old], ccc_old,
               x[ccc_old], y[ccc_old], z[ccc_old]);
                     
      messageflag=0;
    }
  }
  if(num>0)
  {
    printf("%d NAN/INFs were detected on level %d at time %g, iteration %d.\n",
           ((int) num), level->l, level->time, level->iteration);
    fflush(stdout);

    if(num<1) return -1;
  }

  return ((int) num);
}

/* check if var with name varname is finite */
int CheckIfFinite(tL* level, char *varname)
{
  int ret=0, ivar=Ind(varname);
  ret = CheckIfFinite_VarIndexOne(level, ivar);
  return ret;
}

/* check a var with index i with all its components 
   and return the number of NANs found               */
int CheckIfFinite_VarIndex(tL* level, int i)
{
  int j, ret=0, n = VarNComponents(i);
  
  for (j=0; j<n; j++)  ret = ret + CheckIfFinite_VarIndexOne(level, i+j);
  return ret;
}

/* check a var with index i with all its components 
   and return the number of NANs found               */
int CheckIfFinite_vl(tVarList *vl)
{
  int i, ret=0;
  
  for(i=0; i < vl->n; i++)
    ret = ret + CheckIfFinite_VarIndexOne(vl->level, vl->index[i]);
  if(ret) printf("Found NAN/INF in varlist %p!\n", vl);
  return ret;
}



/* BB 5/2007: 
   alternative version of CheckIfFinite
   insert calls to this function for debugging
   note that we have to communicate our result between processors
*/
int vlfinite(tVarList *vl, char *message, int pr)
{
  int flag = 0;
  int i;
  double *sum;

  if (pr) printf("%s, checking for NANs, l %d, t %.3f   (%p)\n",
	 message, vl->level->l, vl->level->time, vl); 

  bampi_allreduce_sum(vl, &sum);

  for (i = 0; i < vl->n; i++)
    if (!finite(sum[i])) {
      flag = 1;
      printf("found NAN in %s\n", VarName(vl->index[i]));
    }

  free(sum);

  return flag;
}

/* wrapper for single variable */
int varfinite(tL *level, char *name, char *message, int pr)
{
  int flag;
  tVarList *vl = vlalloc(level);

  vlpush(vl, Ind(name));
  flag = vlfinite(vl, message,pr);
  vlfree(vl);
  return flag;
}

int allfinite(tL *level, int pr)
{
    int i,flag;
    tVarList *vl = vlalloc(level);
    
    for (i = 0; i < level->grid->nvariables; i++) {
        if (pr) printf("%s   %d %d\n",VarName(i),i,level->grid->nvariables);
        if (level->v[i]) vlpush(vl, i);
        i += VarNComponents(i)-1;
    }
    flag = vlfinite(vl,"all",pr);
    vlfree(vl);
    return flag;
}





/* check for NANs, alternative version, to be registered 
   currently registered in POST_EVOLVE so that NANs are not output
   register in PRE_EVOLVE so that the previous output step has completed
*/
int ExitIfNAN(tL* level)
{
  if (Getv("ExitIfNAN","no")) return 0;
  
  int v      = Getv("ExitIfNAN","verbose");
  int all    = Getv("ExitIfNAN","all");
  char *vars = Gets("ExitIfNAN_vars");
  if (v) printf("Check for NAN's on level %d (%s component)\n",level->l,all?"all":"first");
  
  char *var;
  
  int flag = 0;
  while (var = NextEntry(vars)) {
    tVarList* vl = vlalloc(level);
    if (all)
      vlpush(vl,Ind(var));
    else
      vlpushone(vl,Ind(var));
    
   flag += vlfinite(vl, var, v); 
   vlfree(vl);
  }
  
  
  if (flag) {
    printf("ExitIfNAN: something is not finite after evolve!\n"  
        "  level %d, iteration %d\n"
        "  crashtime %f, outdir %s\n",
        level->l, level->iteration+1, 
        (level->iteration+1) * level->dt, Gets("outdir"));
    errorexit("Too bad.");
  }
  
  return 0;
}



/* check several variables if finite, can be 
   used in math files after rhs computation */
int CheckForNANandINF(int N, ...)
{
  int found = 0;
  int i;
  double value;
    
  va_list ptr;
  va_start(ptr,N);
    
  for (i=0; i<N; i++) {
    value = va_arg(ptr,double);
    if (!finite(value)) {
      if (found==0) printf("found NAN/INF in var:  ");
      printf("%d:%g ", i+1, value);
      found++;
    }
  }
  if (found) {printf(" of %d\n",N);}
    
  va_end(ptr);
    
  return found;
}






