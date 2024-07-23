/* profile.c */
/* Bernd Bruegmann 4/2006 */

/* Profiling MPI programs is not straightforward:
   - there is special commercial software (VTune etc.) which
     gives professional results but is complicated to use, not necessarily
     robust, machine dependent
   - the simplest general mechanism is based on gcc -pg and gprof, which
     for MPI programs however requires recompilation of MPI and linking
     to appropriate libraries since otherwise time spent in MPI routine is
     dropped
   - to be tried: TAU, gprof with recompiling MPI, VTune

   Here we do a first trivial experiment with gcc's
   "instrument-functions" option. This generates a call whenever a
   function is entered and exited. The following adds some minimal
   functionality to keep track of what happens. Accuracy is as
   expected atrocious but the results are not completely useless. Main
   limitation is the accuracy of the timer, MPI_Wtime, second is the
   overhead of instrumentation for short functions that we don't know
   how to subtract.

   The method is:
   - on function entry store time stamp
   - on function exit, add the difference between the exit and entry times
     to a running total for this function
   (Recursive functions have incorrect timing, which could be fixed by
   checking whether a function is already on the stack.)

   Note that the time for each function includes the time for all functions
   that it calls. This includes time for MPI functions, which is what we want.

   To compile with profiling use
     gcc -DPROFILING -finstrument-functions

   Run as normal.
   Each processor writes a outdir/profile.n.
*/


#ifdef PROFILING

#include "bam.h"
#include "bampi.h"




/* store information at index derived from function address */
/* lazy storage */
#define NCOUNT 5000000
int count_enter[NCOUNT];
int count_exit[NCOUNT];
char *symbol[NCOUNT];
double timespent[NCOUNT];

typedef unsigned long int tAddress;
tAddress min_address = 0, max_address = 0, base_address = 0;
#define indofadd(address) \
(((tAddress)(address)) - base_address + NCOUNT/2)

/* function call stack */
#define NSTACK 10000
void *stack[NSTACK];
double timestamp[NSTACK];
int nstack = -1;



/* protect the profiling function against profiling themselves */
void stack_push(void *p)
  __attribute__ ((no_instrument_function));
void stack_pop(void *p) 
  __attribute__ ((no_instrument_function));
tAddress profile_func_index(void *p)
  __attribute__ ((no_instrument_function));
void __cyg_profile_func_enter(void *this_fn, void *call_site) 
  __attribute__ ((no_instrument_function));
void __cyg_profile_func_exit(void *this_fn, void *call_site)
  __attribute__ ((no_instrument_function));
int profile_func_symbols(char *exename)
  __attribute__ ((no_instrument_function));
void bampi_profile(char *exename)
  __attribute__ ((no_instrument_function));




/* gcc instrumentation function */
void __cyg_profile_func_enter(void *this_fn, void *call_site) 
{
  tAddress i = profile_func_index(this_fn);

  count_enter[i]++;

  stack_push(this_fn);
  
  if (0) 
    printf("entering %p from %p at %f\n", this_fn, call_site, MPI_Wtime());
}




/* gcc instrumentation function */
void __cyg_profile_func_exit(void *this_fn, void *call_site)
{
  tAddress i = profile_func_index(this_fn);

  stack_pop(this_fn);

  count_exit[i]++;

  if (0) 
    printf("leaving %p from %p at %f\n", this_fn, call_site, MPI_Wtime());
}




/* push function address on stack, record timestamp */
void stack_push(void *p) 
{
  nstack++;
  if (nstack >= NSTACK) {
    printf("out of stack space\n");
    exit(1);
  }
  stack[nstack] = p;
  timestamp[nstack] = MPI_Wtime();
  if (0) printf("push %p onto %d at %f\n", p, nstack, timestamp[nstack]);
}




/* pop function address from stack, add up time spent in function */
void stack_pop(void *p) 
{
  double t = MPI_Wtime();
  double dt;
  int j;

  if (stack[nstack] != p) {
    printf("stack pop failed: %p, nstack %d\n", p, nstack);
    for (j = 0; j < nstack; j++)
      printf("%p\n", stack[j]);
    exit(1);
  }
  if (nstack < 0) {
    printf("no more entries in stack\n");
    exit(1);
  }
 
  dt = t-timestamp[nstack];
  timespent[indofadd(p)] += dt;
  nstack--;

  if (0) printf("pop  %p from %d at %f, dt = %f\n", p, nstack, t, dt);
}




/* return function address as integer, keep track of range */
tAddress profile_func_index(void *p)
{
  tAddress address = (tAddress) p;
  tAddress i;

  if (base_address == 0) 
    base_address = min_address = max_address = address;
  if (min_address > address) min_address = address;
  if (max_address < address) max_address = address;

  i = indofadd(address);
  if (0) printf("got %p %ld, returning %ld\n", p, address, i-NCOUNT/2);
  if (i < 0 || i >= NCOUNT) i = 0;
  return i;
}




/* perform address to symbol mapping
   relies on the unix commands nm and grep and on specific format!
   symbol table is written to file exename.symbols
*/
int profile_func_symbols(char *exename)
{
  char s[1000], *t;
  FILE *fp;
  int rank;
  tAddress address;
  int i;

  /* initialize */
  for (address = min_address; address <= max_address; address++)
    symbol[indofadd(address)] = "(no name for function)";

  /* write symbol file */
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (!rank) {
    sprintf(s, "nm -td %s | grep \" T \" > %s.symbols", exename, exename);
    printf("Executing command %s\n", s);
    system(s);
  }
  MPI_Barrier(MPI_COMM_WORLD);

  /* open file */
  sprintf(s, "%s.symbols", exename);
  fp = fopen(s, "r");
  if (!fp) {
    printf("could not open %s\n", s);
  }
  else {
    /* read file line by line */
    while (fgets(s, 1000, fp)) {
      address = atol(s);
      t = strstr(s, " T ") + 3;
      i = strlen(t)-1;
      if (t[i] == '\n') t[i] = 0;
      if (0) {
	printf("line: %s", s);
	printf("address %d\n", address);
	printf("symbol  %s\n", t);
      }
      i = indofadd(address);
      if (count_enter[i] || count_exit[i])
	symbol[i] = strdup(t);
    }
    fclose(fp);
  }
}




/* write results to stdout
   change to exename.profile.n?
*/
void bampi_profile(char *exename)
{
  char *outdir = Gets("outdir");
  char com[2000], s[1000], f[100];
  FILE *fp;
  double totaltime, totaltimeall, *timespentall;
  tAddress address, i, imin, iall;
  int *countall;
  int j, rank, size;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  if (0) printf("%ld %ld\n",
	       min_address-base_address, max_address-base_address);
  if (!exename) return;

  /* determine name of output file */
  snprintf(f, 100, "%%s/profile.%%0%dd", (int) log10(size)+1);
  snprintf(s, 1000, f, outdir, rank);
  printf("Writing profiling information to %s\n", s);

  /* the code is about to exit, but in our case there are still e.g.
     the calls to bampi_finalize and main on the stack 
     we fake popping those ourselves to get the timing
  */ 
  j = nstack;
  for (; nstack >= 0;)
    stack_pop(stack[nstack]); 
  nstack = j;

  /* get symbol names */
  profile_func_symbols(exename);

  /* get total time */
  for (address = min_address; address <= max_address; address++) {
    i = indofadd(address);
    if (strcmp(symbol[i], "main") == 0) 
      totaltime = timespent[i];
  }


  /* open file */
  fp = fopen(s, "w");
  if (!fp) {
    printf("could not open %s for writing\n", s);
    MPI_Finalize();
    exit(1);
  }

  /* print result in a table */
  for (address = min_address; address <= max_address; address++) {
    i = indofadd(address);
    if (count_enter[i]) {
      fprintf(fp, "%7.1f%% %9.2fs %12d    %s\n", 
	      100*timespent[i]/totaltime, timespent[i], 
	      count_enter[i], symbol[i]);
    }
  }
 
  /* debug, every enter needs an exit */
  if (0) {
    printf("counters: address, number of exits\n");
    for (address = min_address; address <= max_address; address++) {
      i = indofadd(address);
      if (count_exit[i]) 
	printf("%p %d\n", (void *) address, count_exit[i]);
    }
  }
  fclose(fp);

  /* sort */
  snprintf(com, 2000, "sort -n %s > %s~; mv %s~ %s", s, s, s, s);
  printf("Executing command %s\n", com);
  system(com);


  /* do it again for averaged total */

  /* average times across all processes */
  imin = indofadd(min_address);
  iall = indofadd(max_address) - imin + 1;
  countall = (int *) calloc(iall, sizeof(int));
  MPI_Allreduce(count_enter+imin, countall, iall, 
		MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  timespentall = (double *) calloc(iall, sizeof(double));
  MPI_Allreduce(timespent+imin, timespentall, iall, 
		MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  totaltimeall = 0;
  MPI_Allreduce(&totaltime, &totaltimeall, 1, 
		MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  /* process 0 writes result */
  if (!rank) { 

    /* open file */
    snprintf(s, 1000, "%s/profile.n", outdir);
    fp = fopen(s, "w");
    if (!fp) {
      printf("could not open %s for writing\n", s);
      MPI_Finalize();
      exit(1);
    }
    
    /* print result in a table */
    for (address = min_address; address <= max_address; address++) {
      i = indofadd(address);
      if (count_enter[i]) {
	fprintf(fp, "%7.1f%% %9.2fs %12d/%d    %s\n", 
		100*timespentall[i-imin]/totaltimeall,
		timespentall[i-imin]/size, 
		countall[i-imin], size, symbol[i]);
      }
    }
    fclose(fp);

    /* sort */
    snprintf(com, 2000, "sort -n %s > %s~; mv %s~ %s", s, s, s, s);
    printf("Executing command %s\n", com);
    system(com);
  }

  free(countall);
  free(timespentall);
}



/* #ifdef PROFILING */
#endif

#ifndef PROFILING
void bampi_profile(char *exename)
{
}
#endif

