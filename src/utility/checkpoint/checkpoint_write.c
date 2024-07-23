/* checkpoint_write.c */
/* Bernd Bruegmann, 01/2006, Wolfgang Tichy 01/2006 */

/* write checkpoint files 
*/

#include "bam.h"
#include "checkpoint.h"

#define PR 1


/* write checkpoint file 
   - every processor writes its own file
   - assumes that every processor can write its own file
   - assumes that there is one global file system and that it is better
     to have the processors write one after another rather than hitting
     the disks with simultaneous requests from each processor
     (make parameter to choose this behaviour)
*/
void checkpoint_write(tG *g)
{
  if (PR) printf("start checkpoint_write\n");
  
  char suffix[2];
  // sprintf(suffix,""); // gives compiler warnings
  suffix[0] = '\0';
  char *filename = checkpoint_filename(suffix);
  char *filename_new;
  FILE *fp;
  int i;
  time_t t;
  
  /* name of backup file */
  if (PR) printf("checkpoint_write copy name\n");
  filename_new = cmalloc(strlen(filename) + 20);
  sprintf(filename_new, "%s_new", filename);

  /* for all processes */
  for(i = 0; i < bampi_size(); i++)
  {
    
    if (PR) {
      time(&t);
      printf("checkpoint_write cpu %04d          %s",i,ctime(&t));
    }
    
    /* if it is my turn */
    if(i == bampi_rank())
    {
      /* open file for writing */
      fp = fopen(filename_new, "wb");
      if (!fp) errorexits("failed opening %s", filename_new);

      /* now we are ready to go! */
      checkpoint_write_local(g, fp);

      /* done */
      fclose(fp);
    }
    
    if (PR) {
      time(&t);
      printf("checkpoint_write wait for others   %s",ctime(&t));
    }
    
    /* everyone please wait here */
    bampi_barrier();
  }

  /* Move filename_new to filename, but only after we are sure
     all has been written.                                       */
  fp = fopen(filename_new, "rb");
  if(fp)
  {
    fclose(fp);
    system_move(filename_new, filename); 
  }
  
  /* clean up */
  free(filename_new);
  free(filename);
}




/* write checkpoint file on one processor 
   - should add checks for grid size etc.
*/
void checkpoint_write_local(tG *g, FILE *fp)
{
  if (PR) printf("start checkpoint_write_local\n");
  
  tL *level = g->level[g->lmin];
  int i, j, l;
  int n, nv;
  long long int nall = 0;
  int pr = 0;
  int number_of_variables;

  /* write one line text header, check with "head -1 checkpoint.0" */
  fprintf(fp, "$BEGIN_checkpoint_header:  ");
  fprintf(fp, "$lmin = %d , $lmax = %d , ", g->lmin, g->lmax);
  fprintf(fp, "$time = %.3f\n", level->time);
  fflush(fp);

  for(i = g->lmin; i <= g->lmax; i++)
  {
    tL *level_i = g->level[i];
    fprintf(fp, "level %d : iteration = %d , time = %g\n",
            i, level_i->iteration, level_i->time);
  }
  fprintf(fp, "$END_checkpoint_header\n\n");

  /* write parameter database
     some parameters are used like global variables and change during evolution
  */
  fprintf(fp, "$BEGIN_parameter_database:  ");
  fprintf(fp, "$number_of_parameters = %d\n", GetnParameters());
  for (i = 0; i < GetnParameters(); i++)
  {
    fprintf(fp, "%s = ", GetnameInd(i));
    fprintf(fp, "%s\n", GetsInd(i));
  }
  fprintf(fp, "$END_parameter_database\n\n");


  fprintf(fp, "$BEGIN_variables:\n");
  /* for all levels */
  for (l = g->lmin; l <= g->lmax; l++)
  {
    level = g->level[l];
    n = nv = 0;

    /* test for nan's */
    if (1) {
      if (PR) printf("test for nan's\n");
      for (i = 0; i < g->nvariables; i++) {
        if (level->v[i]) {
          if (PR) printf("$variable = %s : level = %d , level->nnodes = %d\n", VarName(i), l, level->nnodes);
          for (j = 0; j < level->nnodes; j++) {
            if (isnan(level->v[i][j])) printf("there is a nan inside the vars\n");
          }
        }
      }
      if (PR) printf("test for nan's finished\n");
    }

    /* write all variables */
    if (Getv("checkpoint_variables", "all"))
    {
      /* for all variables with storage */
      for (i = 0; i < g->nvariables; i++)
      {
	if (level->v[i] != 0)
	{
	  /* write */
	  fprintf(fp, "$variable = %s : level = %d , level->nnodes = %d\n",
	          VarName(i), l, level->nnodes);
	  fwrite(level->v[i], sizeof(double), level->nnodes, fp);
          fprintf(fp, "\n");

	  n += level->nnodes;
	  nv++;
	  if (pr) printf("wrote level %d, nnodes=%d, variable=%s\n", 
			 l, level->nnodes, VarName(i));
	}
      }
    }

    /* write only those variables that are needed
       - needs to be tested, made more flexible!
    */
    else
    {
      tVarList *vl = checkpoint_varlist(level);
      
      /* for all variables in list with storage */
      for (j = 0; j < vl->n; j++)
      {
	i = vl->index[j];
	if (level->v[i] != 0)
	{
	  /* write */
	  fprintf(fp, "$variable = %s : level = %d , level->nnodes = %d\n",
	          VarName(i), l, level->nnodes);
	  fwrite(level->v[i], sizeof(double), level->nnodes, fp);
          fprintf(fp, "\n");

	  n += level->nnodes;
	  nv++;
          if (pr) printf("wrote level %d, nnodes=%d, variable=%s\n",
			 l, level->nnodes, VarName(i));
	}
      }
      vlfree(vl);
    }
    /* info */
    if (1) printf("wrote l %d, %d variables, total of %d points\n", l, nv, n);
    nall += n;
  }
  fprintf(fp, "$END_variables\n\n");

  /* info */
  if (1) printf("wrote levels %d to %d, %d variables, %lld doubles, %lf GB\n",
		g->lmin, g->lmax, nv, nall,
		(double)(nall * sizeof(double)) * 1e-9 );
}


