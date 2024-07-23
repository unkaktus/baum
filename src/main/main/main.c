/* main.c */
/* Bernd Bruegmann, 12/99 */

#include "bam.h"
#include "main.h"

#define PR 1

/**************************************************************************/
/* main */
int main(int argc, char **argv)
{
  tG *g;

  bampi_initialize(&argc, &argv);
#ifdef OPENMP
  bampi_openmp_initialize(&argc, &argv, PR);
#endif

  read_command_line(argc, argv);
  parse_parameter_file(Gets("parameterfile"));
  initialize_libraries();
  timer_start(0, "main");

  while (iterate_parameters())
  {
    init_output();
    g = make_grid(1);
    prdivider(0);
    initialize_grid(g);
    evolve_grid(g);
    finalize_grid(g);
  }

  timer_stop(0, "main");
  finalize_libraries();
  free_stuff(g);

  bampi_finalize(argc, argv);
  return 0;
}

/* read command line */
void read_command_line(int argc, char **argv)
{
  int i;
  int abort = 0;
  FILE *fp;

  if (PR)
    for (i = 0; i < argc; i++)
      printf("argv[%d] = %s\n", i, argv[i]);

  if (argc != 2)
  {
    if (processor0)
    {
      printf("Welcome to bam. BAM or poof, that is the question.\n");
      printf("Usage:    bam name.par\n");
    }
    bampi_finalize(argc, argv);
    exit(0);
  }

  if (processor0)
  {
    prdivider(0);
    printf("Welcome to bam.\n");
    prdivider(0);
  }

  /* got two parameters */
  if (argc == 2)
  {
    char *parfile = (char *)calloc(sizeof(char), strlen(argv[1]) + 40);
    char *outdir = (char *)calloc(sizeof(char), strlen(argv[1]) + 40);
    char *outdirp = (char *)calloc(sizeof(char), strlen(argv[1]) + 40);
    char *outdirb = (char *)calloc(sizeof(char), strlen(argv[1]) + 40);

    /* determine name of parameter file and output directory */
    strcpy(parfile, argv[1]);
    if (!strstr(parfile, ".par"))
      strcat(parfile, ".par");
    strcpy(outdir, parfile);
    *strstr(outdir, ".par") = '\0';

    /* make output directory, save current one */
    strcpy(outdirp, outdir);
    strcat(outdirp, "_previous");
    strcpy(outdirb, outdir);
    strcat(outdirb, "_backup");

    /* check for checkpoint files in previous directory before overwriting it */
    if (processor0)
    {
      char *prev = checkpoint_filename_from_outdir_suffix(outdir, "_previous");
      char *curr = checkpoint_filename_from_outdir_suffix(outdir, "");
      FILE *fpprev = fopen(prev, "r");
      FILE *fpcurr = fopen(curr, "r");

      if (1)
        printf("Safety first: checking data before overwriting previous directory.\n");
      // printf("prev: %s, curr: %s\n", prev, curr);

      if (fpprev && !fpcurr)
      {
        /* print error message here because we haven't redirected stdout, yet, but don't
        want message for every proc */
        printf("Error: %s exists, while %s does not.\n", prev, curr);
        printf("Error: This could result in the loss of the checkpoint!\n");
        abort = 1;
      }
      if (fpprev)
        fclose(fpprev);
      if (fpcurr)
        fclose(fpcurr);
      free(prev);
      free(curr);
    }
    /* more error is above in print statement */
    if (bampi_or(abort))
      errorexit("You might loose data by overwriting the previous directory!");

    /* move old folder */
    if (processor0)
    {
      /* check if a shell is available to execute commands later */
      /* NOTE: many system_... functions do not need a shell */
      if (system(NULL) == 0)
        printf("WARNING: system(NULL)=0 => cannot execute shell commands!\n");

      /* creating a backup of previous and delete old backup*/
      /* move old output dir to previous */
      if (system_isdir(outdirb))
        system_rmrf(outdirb);
      if (system_isdir(outdirp))
        system_move(outdirp, outdirb);
      if (system_isdir(outdir))
        system_move(outdir, outdirp);

      /* make new dir and copy parfile */
      system_mkdir(outdir);
      system_cp(parfile, outdir);
    }
    bampi_barrier();

    /* test for folder location */
    if (PR)
    {
      for (i = 0; i < bampi_size(); i++)
      {
        if (bampi_rank() == i)
        {
          if (i == 0)
          {
            prdivider(0);
            printf("LOCATION: I am currently at\n");
          }
          system("pwd");
        }
        bampi_barrier();
      }
      if (processor0)
        prdivider(0);
    }

    /* redirect stdout and stderr for multiprocessor jobs
       all output is collected in outdir/stdout.001 etc
    */
    if (!processor0)
    {
      char s[1000], f[100];
      sprintf(f, "%%s/stdout.%%0%dd", (int)log10(bampi_size()) + 1);
      sprintf(s, f, outdir, bampi_rank());
      for (i = 100; i >= 0; i--)
      {
        FILE *fp = fopen(s, "r");
        if (!fp)
          break;
        if (i == 0)
        {
          printf("%s\n", s);
          system("pwd");
          errorexit("it seems that the old folder is not yet moved");
        }
        fclose(fp);
        sleep(100);
      }

      freopen(s, "w", stdout);
      freopen(s, "w", stderr);
      printf("generating %s\n", s);
    }

    /* first parameter initializes parameter data base */
    printf("Adding command line parameters\n");
    AddPar("outdir", outdir, "output directory");
    AddPar("parameterfile", parfile,
           "name of parameter file given on command line");
    // AddPar("trace_memory", "no", "enable memory tracing");

    free(parfile);
    free(outdir);
    free(outdirp);
  }

  /* more initialization */
}

/* initialize libraries
   the automatically generated file calls the initializers for each module
*/
void initialize_libraries(void)
{
  prdivider(0);
  printf("Initializing libraries\n");

#include "bam_automatic_initialize.c"
}

void finalize_libraries(void)
{
  prdivider(0);
  printf("Finalizing libraries\n");

#include "bam_automatic_finalize.c"
}

/* initialize grid */
void initialize_grid(tG *g)
{
  int timer = timer_start(0, "initialize_grid");
  int l;

  /* if you have to add infos not related to the grid */
  for (l = g->lmin; l <= g->lmax; l++)
    RunFun(PRE_PRE_PRE_INITIALDATA, g->level[l]);

  for (l = g->lmin; l <= g->lmax; l++)
    RunFun(PRE_PRE_INITIALDATA, g->level[l]);

  prdivider(0);

  /* find out where objects are located and adjust the grid parameters */
  init_physical_objects(g);

  prdivider(0);
  printf("Initializing grid\n");

  /* check if symmetries are allowed */
  test_grid_symmetry(g);

  /* introduce refinements as needed */
  for (l = g->lmin; l <= g->lmax; l++)
    regrid(g, l + 1, 0);

  /* check box setup */
  test_box_consistency(g, 0);

  /* add shells if needed */
  add_grid_shells(g);

  /* if you need to initialize things before writing initialdata */
  for (l = g->lmin; l <= g->lmax; l++)
    RunFun(PRE_INITIALDATA, g->level[l]);

  /* compute initial data */
  for (l = g->lmin; l <= g->lmax; l++)
    RunFun(INITIALDATA_SET, g->level[l]);

  /* some things need to be done right after initial data to obtail actual ADM data*/
  for (l = g->lmin; l <= g->lmax; l++)
    RunFun(INITIALDATA_FINISH, g->level[l]);

  /* some things need to be done after initial data -> GEOMETRY stuff*/
  for (l = g->lmin; l <= g->lmax; l++)
    RunFun(POST_INITIALDATA, g->level[l]);

  /* some things need to be done after initial data -> MATTER stuff*/
  for (l = g->lmin; l <= g->lmax; l++)
    RunFun(POST_POST_INITIALDATA, g->level[l]);

  /* initial data is just another new time slice */
  for (l = g->lmin; l <= g->lmax; l++)
    RunFun(POST_EVOLVE, g->level[l]);

  /* initial data complete */
  prdivider(0);
  printf("Done with initialization\n");
  printf("%2d %5d %7.3f\n", 0, g->level[0]->iteration, g->level[0]->time);

  /* analyze initial data */
  for (l = g->lmin; l <= g->lmax; l++)
    analyze_level(g, l);

  /* checkpoint, just in case we need it here already */
  for (l = g->lmin; l <= g->lmax; l++)
    checkpoint(g, l);

  /* memory tracing */
  prmemory("after initializing grid");

  /* timing */
  timer_stop(0, "initialize_grid");
  timer_print(g);
}

/* evolve grid */
void evolve_grid(tG *g)
{
  tL *tl = g->level[0];
  int l;
  int iterationmax = Geti("iterations");
  double timemax = Getd("finaltime");

  prdivider(0);

  if (timemax > 0)
    iterationmax = timemax / tl->dt + 0.5;

  if (iterationmax > 0)
  {
    printf("Evolving grid for %d top level iterations to time %.3f\n",
           iterationmax, iterationmax * tl->dt);

    if (Getv("stdout_verbose", "outtime"))
    {
      int it = 10000;
      int i;
      double time;
      double dt_min = g->level[g->lmax]->dt;
      printf("  Times:      dt_min=%2.16e \n", g->level[g->lmax]->dt);

      time = it * Getd("0douttime");
      i = (time + 0.1 * dt_min) / dt_min;
      if (time > 0)
      {
        printf("  outTimes:   dt_0d =%2.16e  -> Tout(%10d)=%2.16e\n", Getd("0douttime"), it, time);
        printf("                                             -> Tevo(%10d)=%2.16e\n", i, i * dt_min);
      }

      time = it * Getd("1douttime");
      i = (int)((time + 0.1 * dt_min) / dt_min);
      if (time > 0)
      {
        printf("  outTimes:   dt_1d =%2.16e  -> Tout(%10d)=%2.16e\n", Getd("1douttime"), it, time);
        printf("                                             -> Tevo(%10d)=%2.16e\n", i, i * dt_min);
      }

      time = it * Getd("2douttime");
      i = (time + 0.1 * dt_min) / dt_min;
      if (time > 0)
      {
        printf("  outTimes:   dt_2d =%2.16e  -> Tout(%10d)=%2.16e\n", Getd("2douttime"), it, time);
        printf("                                             -> Tevo(%10d)=%2.16e\n", i, i * dt_min);
      }

      time = it * Getd("3douttime");
      i = (time + 0.1 * dt_min) / dt_min;
      if (time > 0)
      {
        printf("  outTimes:   dt_3d =%2.16e  -> Tout(%10d)=%2.16e\n", Getd("3douttime"), it, time);
        printf("                                             -> Tevo(%10d)=%2.16e\n", i, i * dt_min);
      }

      char str[1000];
      sprintf(str, " %%%ds | \n", 64 + (int)(fabs(log10(dequaleps))));
      printf(str, " ");
      printf("  if you want to output after %d steps, up to this line the numbers HAVE to coincide\n", it);
      printf("  within accuracy=%e   -> %d digets\n", dequaleps, (int)(fabs(log10(dequaleps))));
    }
  }
  if (iterationmax <= 0)
    return;

  /* outer most evolution loop */
  while (tl->iteration < iterationmax)
  {

    /* evolve grid trough all levels */
    advance_levelandsublevels(g, 0);

    /* update since this may change during evolution, say when checkpointing */
    timemax = Getd("finaltime");
    iterationmax = (timemax > 0) ? timemax / tl->dt + 0.5 : Geti("iterations");

    /* check memory use */
    if (0)
      printvariables(tl);

    /* timing */
    timer_print(g);
  }
}

/* time step: recursive advance of all levels l to lmax
   the default currently is to refine in space but have uniform time step
   Berger-Oliger with refinement in space and time is not fully supported yet
*/
void advance_levelandsublevels(tG *g, int l)
{
  int pr = 0;

  if (pr)
    printf("entering alas %d\n", l);

  /* advance the given level */
  advance_level(g, l);

  /* advance all sub levels */
  if (l < g->lmax)
  {
    advance_levelandsublevels(g, l + 1);
    if (l >= g->lbo)
      advance_levelandsublevels(g, l + 1);
  }

  /* now that all sublevels have caught up with the given level,  communicate with parent level */
  if (l > g->lmin)
    restrict_prolong_evolve(g, l - 1, l);

  /* analyze */
  timer_start(0, "analyze");
  analyze_level(g, l);
  timer_stop(0, "analyze");

  /* move boxes */
  timer_start(0, "move_box");
  if (1)
    prmemory_adv("before move_box", g);
  if (move_box(g, l + 1))
  {

    /* do something like cons2prim after movebox to avoid nans */
    int ll;
    for (ll = l + 1; ll <= g->lmax; ll++)
      RunFun(POST_MOVEBOX, g->level[ll]);

    /* test box setting to see whether there is something wrong */
    test_box_consistency(g, l);
  }
  if (1)
    prmemory_adv("after move_box", g);
  timer_stop(0, "move_box");

  /* checkpoint */
  timer_start(0, "checkpoint");
  checkpoint(g, l);
  timer_start(0, "checkpoint");

  /* info */
  if (l == 0)
    prmemory("after advance");
  if (pr)
    printf("leaving alas %d\n", l);
}

/* advance level l for one time step */
void advance_level(tG *g, int l)
{
  tL *level = g->level[l];
  int i;

  /* pre evolve */
  RunFun(PRE_EVOLVE, level);

  /* evolve */
  if (0)
    printf("l%d evolve\n", l);
  RunFun(EVOLVE, level);

  /* set boundaries */
  if (0)
    printf("l%d setrefbound\n", l);
  if (0 && l)
    prvar(level, "g");
  // setrefbound(level);
  if (0 && l)
    prvar(level, "g");

  /* post evolve */
  RunFun(POST_EVOLVE, level);

  /* evolution step complete */
  level->iteration++;
  level->time = level->iteration * level->dt;

  /* print info */
  if (g->lmax == 0)
  {
    if (timeforoutput_any(level))
      printf("%10d %13.3f\n", level->iteration, level->time);
  }
  else
  {
    for (i = 0; i < l; i++)
      printf("  ");
    printf("%2d %10d %13.3f", l, level->iteration, level->time);

    if (Getv("storage_verbose", "mem"))
    {
      for (i = l; i < g->lmax; i++)
        printf("  ");
#ifdef __APPLE__
      printf("storage_verbose = mem: mallinfo not available on Apple\n");
#else
// mallinfo is deprecated, add -DUSE_MALLINFO2 to the compile flags to use
// the new mallinfo2
#ifdef USE_MALLINFO2
      struct mallinfo2 mem = mallinfo2();
      printf("    %15zu %15zu", mem.uordblks, mem.fordblks);
#else
      struct mallinfo mem = mallinfo();
      printf("    %15d %15d", mem.uordblks, mem.fordblks);
#endif
#endif
    }

    if (Getv("stdout_verbose", "dtdxdiss") && g->level[0]->iteration == 1)
      printf("    %13.6f %13.6f     %13.10f", level->dt, level->dx, get_dissipation_factor(level));

    printf("\n");
  }
  fflush(stdout);
}

/* analyze and output */
void analyze_level(tG *g, int l)
{
  int lmin = l;

  /* wait until this level is the coarsest level of its time */
  if (l > g->lmin && !dless(g->level[l]->time, g->level[l - 1]->time))
    return;

  /* if there is something to do before analyze
   => till now only for setting AHF-found flag */
  for (l = g->lmax; l >= lmin; l--)
    RunFun(PRE_ANALYZE, g->level[l]);

  /* now analyze this and all sublevels, may do its own output */
  for (l = g->lmax; l >= lmin; l--)
    RunFun(ANALYZE, g->level[l]);

  /* some functions, which use results known after ANALYZE, e.g. hydroanalysis*/
  for (l = lmin; l <= g->lmax; l++)
    RunFun(POST_ANALYZE, g->level[l]);

  /* output for permanent variables */
  for (l = lmin; l <= g->lmax; l++)
    RunFun(OUTPUT, g->level[l]);

  /* post output */
  for (l = lmin; l <= g->lmax; l++)
    RunFun(POST_OUTPUT, g->level[l]);
}

/* finalize grid */
void finalize_grid(tG *g)
{
  prdivider(0);
  free_grid(g);
  bampi_barrier();
}

void free_stuff(tG *g)
{

  freeParameters();
  freeFuns();
  freeVars();
  freeTimer();
  freevll();

  prmemory_active();
}
