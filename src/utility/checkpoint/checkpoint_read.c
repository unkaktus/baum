/* checkpoint_read.c */
/* Bernd Bruegmann, 01/2006, Wolfgang Tichy 01/2006 */

/* read checkpoint files 
*/

#include "bam.h"
#include "checkpoint.h"



/* read checkpoint file 
   - every processor reads its own file
*/
void checkpoint_read(tG *g)
{
  FILE *fp;
  int i, l, np;
  int oldread=0;
  char str[10000];

  /* if there is not a complete set of checkpoint files,
     don't do anything: assume that this run starts from the beginnning
  */
  fp = checkpoint_openfiles("_previous");
  if (!fp) return;

  /* for all processes: read the pars and iteration numbers */
  for (i = 0; i < bampi_size(); i++)
  {
    /* read if it is my turn */
    if (i == bampi_rank())
      checkpoint_read_ParsAndIterations_local(g, fp);

    /* everyone please wait here */
    bampi_barrier();
  }

  //printgrid(g);
  
  /* amr does use its own parameters, which have to be set */
  if (ExistPar("moving_puncture1_x")) {
    for (np=1; np<=Geti("nobjects"); np++) {
      for (i=0; i<3; i++) {
        sprintf(str,"moving_puncture%d_%c",np,(i==0)?'x':(i==1)?'y':'z');
        g->puncpos[np-1][i] =  Getd(str);
      }
    }
  }
  
  /* move all boxes according to what the new pars say */
  printf("moving boxes where they should be according to the parameters:\n");
  for(l = g->lmax; l >= g->lmin; l--)
    move_box(g, l);

  bampi_barrier();

  if ( (Getv("physics", "matter")) && GetvLax("conservative_amr","yes")) {
    for(l = g->lmax; l >= g->lmin; l--)
      CAMR.mask(g->level[l]);
  }

  //printgrid(g);
  
  /* for all processes: read the vars */
  for (i = 0; i < bampi_size(); i++)
  {
    /* read if it is my turn */
    if (i == bampi_rank())
      checkpoint_read_Vars_local(g, fp);

    /* everyone please wait here */
    bampi_barrier();
  }

  /* clean up */
  fclose(fp);
}




/* utility function: strip trailing newline */
void striptrailingnewline(char *s)
{
  int n = strlen(s) - 1;
  if (n >= 0 && s[n] == '\n')
    s[n] = 0;
}




/* read pars, iteration numbers and times from checkpoint file on 
   one processor */
void checkpoint_read_ParsAndIterations_local(tG *g, FILE *fp)
{
  char name[10000], value[10000], str[10000], *currentvalue;
  int i,j,ret;
  int lmin,lmax;
  int numberofpars, parsread;

  if (1) printf("reading parameters and iterations from checkpoint file:\n");

  /* read one line text header, check with "head -1 checkpoint.0" */
  if(fgotonext(fp, "$BEGIN_checkpoint_header:") == EOF)
    errorexit("$BEGIN_checkpoint_header: is missing!");
  fgetparameter(fp, "$lmin", value);  lmin = atoi(value);
  fgetparameter(fp, "$lmax", value);  lmax = atoi(value);
  fgetparameter(fp, "$time", value);
  printf("found $lmin = %d , $lmax = %d , $time = %s\n",
         lmin, lmax, value);

  /* WT: read iteration and time on each level */
  for(i = lmin; i <= lmax; i++)
  {
    tL *level_i = g->level[i];

    fgetparameter(fp, "level", value);
    j = atof(value);
    if(j != i)  errorexit("wrong level number");

    fgetparameter(fp, "iteration", value);
    level_i->iteration = atof(value);
  
    fgetparameter(fp, "time", value);
    level_i->time = atof(value);

    /* the time in the header is just for human inspection
       use same definition as in main.c, advance_level()   */
    level_i->time = level_i->iteration * level_i->dt;

    printf("  on level %d, iteration=%d, time=%g\n",
           i, level_i->iteration, level_i->time);
  }

  /* read parameter database
     some parameters are used like global variables and change during evolution
  */
  if(fgotonext(fp, "$BEGIN_parameter_database:") == EOF)
    errorexit("$BEGIN_parameter_database: is missing!");

  fgetparameter(fp, "$number_of_parameters", value);
  numberofpars = atoi(value);
  printf("  found %d parameters\n", numberofpars);

  parsread = 0;
  while( (ret=fscanline(fp, str)) != EOF )
  {
    while(ret==0)  ret = fscanline(fp, str);
    if( ret==EOF || strcmp(str,"$END_parameter_database")==0 ) break;

    if(0) printf("  str: %s\n", str);
    extrstr_before_after_EQ(str, name, value);
 
    /* remove spaces at the end of name and the beginning of value */
    strcpy(str, name);
    sscanf(str, "%s", name);
    strcpy(str, value);
    j = strspn(str, " ");  /* strcpy(value, str+j); */
    if(j>0) strcpy(value, str+1); /* <--we remove at most one space */
    if(0) printf("  found parameter %s = ", name);
    if(0) printf("%s\n", value);
    currentvalue = GetsLax(name);
    if(currentvalue==0)
    {
      printf("skipping %s = %s\n", name, value);
      continue;
    }
    if(strcmp(value, currentvalue) != 0)
    {
      /* handle list of parameters that are allowed to keep a new value
         (should be made actual list)
      */
      if(strcmp(name, "checkpoint") &&
	 strcmp(name, "finaltime")  &&
	 strcmp(name, "iterations")  && 
	 strcmp(name, "amr_lmax")
	 )
      {
	if (1) printf("overwriting %s = %s -> %s\n",
		      name, currentvalue, value);
	Sets(name, value);
      }
    }
    parsread++;
  }
  printf("  read %d parameters\n", parsread);
  
  /* delete all levels that we don't need any more 
  (after merger, some could have been removed) */
  if (0) {
    if(g->lmax > lmax) {
      printf("removing level %d and all sublevels\n", lmax+1);
      remove_levels(g->level[lmax + 1]);
    }
  } else {
    // better not delete, only ignore
    if(g->lmax > lmax) {
      printf("hiding level %d and all sublevels\n", lmax+1);
      g->lmax = lmax;
    }
  }
}



/* read vars from checkpoint file on one processor 
   - assumes identical executable!
   - should add checks for grid size, variable names etc.
*/
void checkpoint_read_Vars_local(tG *g, FILE *fp)
{
  tL *level;
  int ilax, i, j, l, nnodes_on_l;
  int lold = -11111;
  int n, nv;
  long long int nall = 0;
  int pr = 0;
  char varname[10000], value[10000];

  if (1) printf("reading variables from checkpoint file:\n");
  if(fgotonext(fp, "$BEGIN_variables:") == EOF)
    errorexit("$BEGIN_variables: is missing!");  
  
  l = lold;
  /* for all saved variables */
  while(fgetparameter(fp, "$variable", varname) != EOF)
  {
    ilax = IndLax(varname);
    if(ilax>=0)
      i = ilax;
    else
    {
      printf("variable %s in checkpoint file does not exist!\n", varname);
      i = 0;
    }

    fgetparameter(fp, "level", value);
    l = atoi(value);
    level = g->level[l];
    if(l != lold)
    {
      lold = l;
      n = nv = 0;
    }

    fgetparameter(fp, "level->nnodes", value);
    nnodes_on_l  = atoi(value);
    if(nnodes_on_l != level->nnodes) {
      printf("  %d  %d  var: %s    level: %d \n",nnodes_on_l,level->nnodes,varname,l);
      errorexit("level->nnodes in checkpoint file is not what it should be!");
    }
    
    /* read the trailing newline behind the level->nnodes par */
    fscanline(fp, value);

    /* read that variable */
    if(ilax>=0 && level->v[i] != 0)
    {
      /* read */
      fread(level->v[i], sizeof(double), nnodes_on_l, fp);

      nv++;
      n += nnodes_on_l;
      nall += nnodes_on_l;
      if(pr) printf("read level %d, nnodes=%d, variable=%s, index=%d,\n", 
		     l, nnodes_on_l, VarName(i), i);
    }
    else
    {
      //errorexit("we should overread the next nnodes_on_l doubles in fp");
      printf("skipping the %d double numbers of %s (on level %d) ...\n",
             nnodes_on_l, varname, l);
      for(j=1; j<=(int)(nnodes_on_l*sizeof(double)/sizeof(char)); j++) fgetc(fp);
    }
  }
  /* info */
  if (1) printf("read levels %d to %d, %d variables, %lld doubles, %lf GB\n",
		g->lmin, g->lmax, nv, nall, 
		(double)(nall * sizeof(double)) * 1e-9 );
}



