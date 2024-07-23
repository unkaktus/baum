/* output.c */
/* Bernd Bruegmann 12/99 */
/* dtim 07/14 */

#include "bam.h"
#include "output.h"

char boxstring[2] = {0, 0};
int boxnr = 0;

int output_order = 4;
int output_scheme = LAGRANGE;





/* decide whether it is time for output 
   return 0 if not
   return number of output, 1 for t=0, 2 for next, ...
   tricky bit is that we have two criteria, di and dt, and dt could mean 
   also 'output after this time is exceeded' if dt and the timestep don't match
*/
int timeforoutput_di_dt(tL *level, int di, double dt) 
{
  /* time for output based on number of iterations */
  if (di > 0 && level->iteration % di == 0) 
    return level->iteration/di + 1;

  /* time for output based on time interval, assumes t >= 0 */
  double dt_min = level->grid->level[ level->grid->lmax ]->dt;
  if (dt > 0) {
    int i = (level->time + 0.1*dt_min)/dt;
    if (dequal(level->time-i*dt, 0))
      return i + 1;
  }
    
  /* not time for output */
  return 0;
}

int timeforoutput_dt(tL *level, double dt)
{

  /* time for output based on time interval, assumes t >= 0 */
  double dt_min = level->grid->level[ level->grid->lmax ]->dt;
  if (dt > 0) {
    int i = (level->time + 0.1*dt_min)/dt;
    if (dequal(level->time-i*dt, 0))
      return i + 1;
  }

  /* not time for output */
  return 0;
}


/* decide whether it is time for output in any format */
int timeforoutput_any(tL *level)
{
  return 
    timeforoutput_di_dt(level, Geti("0doutiter"), Getd("0douttime")) ||
    timeforoutput_di_dt(level, Geti("1doutiter"), Getd("1douttime")) ||
    timeforoutput_di_dt(level, Geti("2doutiter"), Getd("2douttime")) ||
    timeforoutput_di_dt(level, Geti("3doutiter"), Getd("3douttime"));
}




/* decide whether it is time for output for a particular variable 
   the variable has to be present in any one of the output parameters, and
   the iteration count has to be right
*/
#define NOUTPUT 4
int timeforoutput_index(tL *level, int index)
{
    
  static int firstcall = 1;
  static int di[NOUTPUT];
  static double dt[NOUTPUT];
  static char output[NOUTPUT][30];
  char *name, s[30];
  int i, d, n;

  if (firstcall) {
    
    firstcall = 0;

    /* cache */
    n=0;
    for (d = 0; d <= 3; d++,n++) {
      sprintf(s, "%ddoutiter", d);
      di[n] = Geti(s);
      sprintf(s, "%ddouttime", d);
      dt[n] = Getd(s);
      sprintf(output[n], "%ddoutput", d);
    }
  }

  name = VarName(index);

  for (n=0; n<NOUTPUT; n++) {
    
    if (Getv(output[n], name)) { 
      if (timeforoutput_di_dt(level, di[n], dt[n])) 
	return 1;
    }
  }
  
  return 0;
}




/* decide whether it is time for output for any member of a variable list */
int timeforoutput(tL *level, tVarList *vl)
{
  int i;
  
  for (i = 0; i < vl->n; i++) {
    if (timeforoutput_index(level, vl->index[i])) return 1;
  }
  
  return 0;
}




/* make complete variable index list for output
   - expand single components to all components if allflag is nonzero 
   - ignore variables without storage
*/
void makeoutputlist(tL *level, char *string, int allflag, 
		    int *pnindex, int **pindex)
{
  int pr = 0;
  int nindex = 0;
  int *index =  malloc (globalnvariables*sizeof(int));
  char *name;
  int i, j, k, n;

  /* for each name */
  while (name = NextEntry(string)) {
    
    /* find range of variable indices */
    i = IndLax(name);
    if (i < 0) continue;
    n = 1;
    if (allflag) {
      i = IndComponent0(i);
      n = VarNComponents(i);
    }
    if (allflag==2) {
      i = 0;
      n = globalnvariables;
    }

    /* for all variable indices */
    for (j = i; j < i+n; j++) {
    
      /* ignore variable if storage is not enabled */
      if (!level->v[j]) continue;

      /* ignore variable if it is constant in time and iteration > 0 */
      if (VarConstantFlag(j) && level->iteration) continue;

      /* check whether index is already in the list */
      for (k = 0; k < nindex; k++)
        if (index[k] == j) break;

      /* add index if not already in list */
      if (k == nindex) 
        index[nindex++] = j;
    }
  }

  /* return list of indices */
  *pnindex = nindex;
  *pindex = index;

  if (pr) {
    printf("nindex = %d, vars = \n", nindex); 
    for (i = 0; i < nindex; i++)
      printf(" %s", VarName(index[i]));
    printf("\n");
  }
}










#define N0D 11
#define N1D 8
#define N2D 6
/* master function to do all the writing */
int write_level(tL *level)
{
  int timer = timer_start(0, "write_level");
  char *name, *string, str[10000];
  int d, l, lmax,vi;
  int nindex, *index;
  double value,x,y,z;
  
  static int firstcall = 1;
  static int di[4];
  static double dt[4];
  static char *ou[4];
  static int all[4];
  static int ou0d[N0D],ou1d[N1D],ou2d[N2D];
  static char *string0d[N0D] = {"max", "min", "norm","norminf","integral","normmask","integralmask", "integralouter", "pmaxmin", "punc","pt"};
  static char *string1d[N1D] = {"d", "x", "y", "z", "xy","xz","yz", "pxy"};
  static char *string2d[N2D] = {"xy", "xz", "yz", "xd", "yd", "zd"};

  if (Getv("outputinterpolatescheme","WENO")) output_scheme = WENO;



  if (firstcall) {
    firstcall = 0;

    /* cache */
    for (d = 0; d <= 3; d++) {
      sprintf(str, "%ddoutiter", d);
      di[d] = Geti(str);
      sprintf(str, "%ddouttime", d);
      dt[d] = Getd(str);
      sprintf(str, "%ddoutput", d);
      ou[d] = Gets(str);
      sprintf(str, "%ddoutputall", d);
      all[d] = Getv(str, "yes");
    }
    for (d = 0; d < N0D; d++) {
      ou0d[d] = 0;
      string = Gets("0doutputmode");
      while (name = NextEntry(string))
        if (strcmp(name,string0d[d])==0) ou0d[d] = 1;
    }
    for (d = 0; d < N1D; d++) {
      ou1d[d] = 0;
      string = Gets("1doutputmode");
      while (name = NextEntry(string))
        if (strcmp(name,string1d[d])==0) ou1d[d] = 1;
    }
    for (d = 0; d < N2D; d++) {
      ou2d[d] = 0;
      string = Gets("2doutputmode");
      while (name = NextEntry(string))
        if (strcmp(name,string2d[d])==0) ou2d[d] = 1;
    }
  }



  /* loop over all boxes by creating temporary one box level */
  tL *level0 = level;
  int i;
  for (i = 0; i < level0->nboxes; i++) {
    level = one_box_level(level0, i);
    if (level0->nboxes == 1) boxstring[0] = 0;
    else boxstring[0] = 'a' + i;
    boxnr = i;

    if (0) {
      printf("l%d output loop %d\n", level->l, i);
      printbbox(level, level->bbox, level->ibbox);
      printbbox(level, level->com->bbox, level->com->ibbox);
    }
    
    
    
    /* 0d output (norm)*/
    d = 0;
    if (timeforoutput_di_dt(level, di[d], dt[d]) &&
        level->l >= Geti("0doutlmin") && level->l <= Geti("0doutlmax")) {
      
      makeoutputlist(level, ou[d], all[d], &nindex, &index);
      for (l = 0; l < N0D; l++) {
        if (ou0d[l]) {
          if (l<8) {
            write_level_0d(level, nindex, index, l, string0d[l]);
          } else {
            write_finestlevel_point(level, nindex, index, l-8, string0d[l]);
          }
        }
      }
      free(index);
      
    }
    
    
    
    /* 1d output */
    d = 1;
    if (timeforoutput_di_dt(level, di[d], dt[d]) &&
        level->l >= Geti("1doutlmin") && level->l <= Geti("1doutlmax")) {

      makeoutputlist(level, ou[d], all[d], &nindex, &index);
      for (l = 0; l < N1D; l++) {
        if (ou1d[l]) {
          if (level->shells) {
            boxstring[0] = 0;
            write_level_shells1d(level, nindex, index, l, string1d[l]);
          } else
            write_level_1d(level, nindex, index, l, string1d[l]);
        }
      }
      free(index);

    }
    
    
    
    /* 2d output */
    d = 2;
    if (timeforoutput_di_dt(level, di[d], dt[d]) &&
        level->l >= Geti("2doutlmin") && level->l <= Geti("2doutlmax")) {
      
      makeoutputlist(level, ou[d], all[d], &nindex, &index);
      for (l = 0; l < N2D; l++) {
        if (ou2d[l]) {
          if (level->shells) {
            // the function outputs data from all shells -> call it only once
            if (i==0) {
              boxstring[0] = 0;
              write_level_shells2d(level0, nindex, index, l, string2d[l]);
            }
          } else
            write_level_2d(level, nindex, index, l, string2d[l]);
        }
      }
      free(index);
    }
    
    
    
    /* 3d output */
    d = 3;
    if (timeforoutput_di_dt(level, di[d], dt[d]) &&
        level->l >= Geti("3doutlmin") && level->l <= Geti("3doutlmax")) {
      makeoutputlist(level, ou[d], all[d], &nindex, &index);
    
      if (level->shells) {
        
        /* 3-D output of the spherical shells in level 0 */
        boxstring[0] = 'a' + i;
        write_level_shells3d(level, nindex, index, "xyz");
        
      } else {
      
        /* 3-D output can be done in several formats - e.g. opendx and HDF5 */
        if (Getv("3dformat", "opendx") ||
            Getv("3dformat", "vtk") ||
            Getv("3dformat", "xdmf") ||
            Getv("3dformat", "text"))
          write_level_3d(level, nindex, index, 1, "xyz");
        
      }
      
      /* amrunion output no longer done here */
      free(index);
    }

  }

  
  timer_stop(0, "write_level");
  
  return 0;
}

void init_output()
{
  char str[1000];
  static int firstcall = 1;
  
  if (firstcall) {
    firstcall = 0;
    
    /* generate ?d-output folders */
    if (processor0) {
      sprintf(str,"%s/output_0d",Gets("outdir"));
      AddPar("outdir_0d",str,"");
      system_mkdir(str);
      sprintf(str,"%s/output_1d",Gets("outdir"));
      AddPar("outdir_1d",str,"");
      system_mkdir(str);
      sprintf(str,"%s/output_2d",Gets("outdir"));
      AddPar("outdir_2d",str,"");
      system_mkdir(str);
      sprintf(str,"%s/output_3d",Gets("outdir"));
      AddPar("outdir_3d",str,"");
      system_mkdir(str);
      sprintf(str,"%s/output_r",Gets("outdir"));
      AddPar("outdir_r",str,"");
      system_mkdir(str);
    }
  }
  
}








