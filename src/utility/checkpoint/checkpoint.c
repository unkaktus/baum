/* checkpoint.c */
/* Bernd Bruegmann 01/2006, Wolfgang Tichy 11/2006 */

/* write and read checkpoint files 
*/

#include "bam.h"
#include "checkpoint.h"

#define PR 0


/* main entry into checkpointing */
void checkpoint(tG *g, int l)
{
  tL *level = g->level[l];
  int di, skip_checkpoint_i = 1, skip_checkpoint_t = 1, skip_checkpoint_T = 1;
  double uptime, uptime_since_checkpoint, dt, dT, dt_hours, dt_hours_quit;
  static double last_checkpoint_uptime = 0.;
  int pr = 1;

  uptime = bampi_uptime()/3600.;
  
  if (PR) printf("called checkpoint(.) at %f hours after MPI startup\n\n", uptime);
  bampi_barrier();

  /* wait for time alignment of all levels */
  if (l != g->lmin) return;

  /* check whether we are on */
  if (Getv("checkpoint", "no")) return;
  if (pr) printf("checkpoint called at top level iteration %d\n",
		 level->iteration);

  /* we may not be able to */
  if (Geti("order_timeinterpolation") > 3)
    errorexit("checkpoint: implement order_timeinterpolation > 3");
  // fixme: moving punctures


  /* if we are in the restart phase */
  if (Getv("checkpoint", "restart")) {

    /* here we could do restarts from initial data */ 
    /* if (level->iteration == 0) 1; */
    
    /* clobber data after one full evolution step */
    if (level->iteration == 1) {
      printf("checkpoint restart: read files\n");
      checkpoint_read(g);
      checkpoint_copy_output(g->level[0]->time);
      Sets("checkpoint", "yes");

      /* run for that much longer */
      dT = Getd("checkpoint_DeltaT");
      if (dT > 0) {
	Setd("finaltime", level->time + dT);
	Seti("iterations", ((int) ((level->time + dT)/level->dt + 0.5)));
	printf("new time is %.3f, adding %.3f to the clock, run until %.3f\n",
	       level->time, dT, level->time + dT);
	printf("embarking on iterations %d to %d\n",
	       level->iteration + 1, Geti("iterations")); 
      }
    }

    /* return since we are still starting up and 
       don't want new checkpoint files
    */
    return;
  }

  /* if we get here, we are not in the restart phase */

  /* determine how often we want to checkpoint */
  di            = Geti("checkpoint_di");
  dt            = Getd("checkpoint_dt");
  dt_hours      = Getd("checkpoint_dt_hours");
  dt_hours_quit = Getd("checkpoint_dt_hours_quit");
  if (pr) printf("checkpoint vars:  %d %e %e %e\n", di,dt,dt_hours,dt_hours_quit);
  if (di<1) errorexit("di has to be positive");
  

  if (di==0) errorexit("checkpoint_di = 0 is not allowed");
  if (pr) printf("di = %d\n", di);

  
  if (pr) printf("checkpoint test for step number\n");
  if (dt <= 0. && dt_hours <= 0. && dt_hours_quit <= 0.)  {        // checkpoint by checkpoint_di
     skip_checkpoint_i = (level->iteration % di);
     printf("checkpoint by iteration number, di = %d\n", di);
  }

  
  if (pr) printf("checkpoint test for step time\n");
  if (dt > 0) {
    di = (int)(dt/level->dt + dequaleps);                      // checkpoint by checkpoint_dt
    if (di<=1) errorexit("NO, increase checkpoint_dt");
    skip_checkpoint_t = level->iteration % di;
    if (!skip_checkpoint_t) 
      printf("checkpoint by physical time, level->time = %f, di=%d\n", level->time, di);
  }

  
  if (pr) printf("checkpoint test for walltime\n");
  uptime_since_checkpoint = uptime - last_checkpoint_uptime;   // checkpoint by walltime
  if ((dt_hours      > 0.  &&  dt_hours      <= uptime_since_checkpoint)   ||
      (dt_hours_quit > 0.  &&  dt_hours_quit <= uptime)) {

    skip_checkpoint_T = 0; // checkpoint now
    last_checkpoint_uptime = uptime;
    if (!skip_checkpoint_T) 
      printf("checkpoint by walltime, dt = %f hours\n", uptime_since_checkpoint);
  }


  /* check whether it is time to write a checkpoint */
  if (pr) printf("checkpoint test to write\n");
  if (skip_checkpoint_i && skip_checkpoint_t && skip_checkpoint_T) {
    if (pr) printf("not time to write checkpoint %f hours into the run at t = %f\n", uptime, level->time);
    return;
  }

  /* for now do not write checkpoints after initial data */
  if (level->iteration == 0) {
    return;
  }

  /* write checkpoint */
  if (pr) printf("checkpoint write\n");  
  checkpoint_write(g);
  printf("It took %f minutes to checkpoint %f hours into the run.\n", bampi_uptime()/60. - 60. * uptime, uptime);
  

  bampi_barrier();

  if (0 < dt_hours_quit && dt_hours_quit <= uptime) {

    if (processor0) {
      if (Getv("checkpoint_backup", "yes")){
        printf("backup checkpoint files \n");
        char checkdir[1000];
        sprintf(checkdir, "%s/checkpoint.*", Gets("outdir"));
        system3("cp",checkdir,Gets("checkpoint_backup_dir"));
      }
  
      if (Getv("checkpoint_resub", "yes")){    
        printf("resub simulation\n");
        system2("",Gets("checkpoint_resub_command"));   
      }
    }
  
    printf("Now kill bam before the queuing system does!\n");
    bampi_finalize(0, 0);
    exit(0);
  }

}




/* return pointer to filename, has to be freed by caller 
   this is actually quite tricky since machine dependent
   - assume for now that each processor can write its own output!
   - it may be more efficient to collect all data onto processor 0
     and then write one file
*/
char *checkpoint_filename_from_outdir_suffix(char *outdir, char *suffix)
{
  char *dir = cmalloc(strlen(outdir) + strlen(suffix) + 1);
  char *filename = cmalloc(1000);
  char formatstring[100];

  /* name of directory */
  sprintf(dir, "%s%s", outdir, suffix);
    
  /* variable length format for different number of processors
     (same as for stdout.01 etc.)
  */
  snprintf(formatstring, 100, "%%s/checkpoint.%%0%dd",
	   (int) log10(bampi_size())+1);
  snprintf(filename, 1000, formatstring, dir, bampi_rank());

  if (1) printf("  checkpoint filename = %s\n", filename);

  free(dir);
  return filename;
}

/* return pointer to filename, has to be freed by caller,
   -get outdir from par "outdir" */
char *checkpoint_filename(char *suffix)
{
  return checkpoint_filename_from_outdir_suffix(Gets("outdir"), suffix);
}


/* open checkpoint file for this processor 
   return 0 if not
   return file pointer to opened file
   this file has to be closed by caller
*/
FILE *checkpoint_openfiles(char *suffix)
{
  char *filename = checkpoint_filename(suffix);
  FILE *fp;
  int nfp;

  fp = fopen(filename, "rb");
  nfp = fp ? 1 : 0;
  nfp = bampi_allreduce_sum_int(nfp);

  if (nfp == 0) {
    if (1) printf(
      "  there are no checkpoint files to read yet\n");
    return 0;
  }

  if (nfp < bampi_size()) {
    if (fp) fclose(fp);
    errorexit("there are some, but not all of the expected checkpoint files");
  }

  free(filename);
  return fp;
}




/* check whether there is a checkpoint file for each processor */
int checkpoint_checkforfiles(char *suffix)
{
  FILE *fp = checkpoint_openfiles(suffix);

  if (fp) {
    fclose(fp);
    if (1) printf("  found checkpoint files, going into restart mode\n");
    return 1;
  }

  if (1) printf("  no checkpoint files, no restart\n");
  return 0;
}




/* copy output from previous run into current location so that all
   additional output can continue/append where we left off
*/
void checkpoint_copy_output(double wtime)
{
  int i;
  char s[1000], f[100], svn1[1000],svn2[1000];
  char *suffix   = "_previous";
  char *suffix2  = "_tmp";
  char *current  = Gets("outdir");
  char *previous = cmalloc(strlen(current) + strlen(suffix) + 1);
  char *tmp      = cmalloc(strlen(current) + strlen(suffix2) + 1);


  /* name of previous/tmp directory */
  sprintf(previous, "%s%s", current, suffix);
  sprintf(tmp, "%s%s", current, suffix2);
  
  /* checkpoint file name */
  sprintf(f, "%%s/stdout.%%0%dd", (int) log10(bampi_size())+1);
  
  /* redirect stdout temporary, so that we can move safely -> create tmp-stdout-files */
  if (processor0) {
    system_mkdir(tmp); 
    // create all files 
    for (i=0; i<bampi_size(); i++) {
      sprintf(s, f, tmp, i);
      FILE *file = fopen(s, "w");
      fclose(file);
    }
  }
  
  /* check files */
  if (processor0) {
    for (i=0; i<bampi_size(); i++) {
      FILE *file = fopen(s, "r");
      if (!file) errorexit("something is wrong with the tmp-folder");
      fclose(file);
    }
  }
  bampi_barrier();
  
  /* now redirect */
  if (!processor0) {
    sprintf(s, f, tmp, bampi_rank());
    freopen(s, "a", stdout);
    freopen(s, "a", stderr);
  }
  bampi_barrier();
  
  /* backup all restarts in svn_info */
  if (processor0) {
    FILE *file1;
    FILE *file2;
    sprintf(svn1, "%s/svn_info", current);
    sprintf(svn2, "%s/svn_info", previous);
    file1 = fopen(svn1, "r");
    file2 = fopen(svn2, "a+");
    if (!file1) printf("checkpoint_copy_output: cannot open %s\n", svn1);
    if (!file2) printf("checkpoint_copy_output: cannot open %s\n", svn2);

    if (file1!=NULL && file2!=NULL) {
      time_t t;
      time(&t);
      fprintf(file2, "\nRestart  %s  T = %2.6f\n",ctime(&t),wtime);
      int i=0;
      while (fgets(s,sizeof s,file1) != NULL) {
        if (i>1) fprintf(file2,"%s\n", s);
        i++;
      }
    }

    if (file1) fclose(file1);
    if (file2) fclose(file2);
  }
  bampi_barrier();
  
  
  
  

  /* remove current output directory */
  if (processor0) { 
    system_rmrf(current);
  }
  bampi_barrier();
  
  /* check if current is gone and previous still exist */
  if (processor0) {
    DIR *dir = opendir(current);
    if (dir) errorexit("current DIR is not removed");
    dir = opendir(previous);
    if (!dir) errorexit("previous DIR does not exist");
    closedir(dir);
  }
  bampi_barrier();
  
  
  
  /* move or copy previous output */
  if (processor0) {
    if (Getv("checkpoint_previous", "no"))
      system_move(previous, current);
    else
      system_cpr(previous, current);
  }
  bampi_barrier();
  
  /* check if current is gone and previous still exist */
  if (processor0) {
    DIR *dir = opendir(current);
    if (!dir) errorexit("current DIR does not exist");
    closedir(dir);
    dir = opendir(previous);
    if (dir) errorexit("previous DIR is not moved");
  }
  bampi_barrier();
  
  
  
  /* redirect stdout and stderr for multiprocessor jobs */
  if (!processor0) {
    sprintf(s, f, current, bampi_rank());
    freopen(s, "a", stdout);
    freopen(s, "a", stderr);
  }
  bampi_barrier();
  
  /* remove temp folder */
  if (processor0) {
    system_rmrf(tmp);
  }
  bampi_barrier();
  
  
  /* clean up */
  free(previous);
  free(tmp);
}




/* create list of variables that are really needed for checkpoint
   caller has to free storage
   - needs to be tested, made more flexible!
   - could reduce storage by factor 3 for AMR: 
     AMR currently stores three time levels, but we could restart
     the evolution using just one time level as new initial data
*/
tVarList *checkpoint_varlist(tL *level) 
{
  tVarList *u, *up, *upp;
  tVarList *vl = vlalloc(level);
  int i;
  
  /* print all variables */
  if (0) printvariables(level);

  /* get all variables that are maintained for evolution
     note that the ADM variables are derived quantities e.g. for BSSN
  */
  evolve_vlretrieve(&u, &up, &upp);
  vlpushvl(vl, u);
  vlpushvl(vl, up);
  vlpushvl(vl, upp);

  /* add all variables with storage contained in
     checkpoint_additional_variables */
  for (i = 0; i < level->grid->nvariables; i++)
  {
    if (level->v[i] != 0)
      if( Getv("checkpoint_additional_variables", VarName(i)) )
      {
        vldrop(vl, i); /* <--drop it in case it is there already */
        vlpush(vl, i); /* <--now add it again */
      }
  }

  /* there are other variables that we may need/want
     - x, y, z: available from startup
     - psi and derivs: should be set during initial data (!)
     - K_initial: ?
     some of these could be added to the par checkpoint_additional_variables
     by other modules...
  */

  /* print variables that are considered to be essential */
  if (0) prvarlist(vl);

  return vl;
}

/* init some things so that checkpoint works properly */
int checkpoint_init(tL *level)
{
  bampi_uptime(); /* make sure t0 in bampi_uptime is set early */

  return 0;
}
