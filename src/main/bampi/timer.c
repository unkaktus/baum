/* timer.c */
/* Bernd Bruegmann 3/2007 */

/* Since profiling turned out to be rather involved (see profile.c),
   in particular since profiled code runs very slowly,
   we provide here some simple timing information that can be turned
   on for standard runs.
*/

#include "bam.h"
#include "bampi.h"

#define PR 0


/* time bin data base */
typedef struct {
  char *name;
  int l;
  double start;
  double stop;
  double time;
  int n;
} tTimer;

static tTimer *tdb = 0;
int ntimers = 0;

static int timer_on = -1;
static int timer_barrier = 0;




/* return pointer to entry for name
   if name does not exist yet, allocate
*/
tTimer *timer_get(tL *level, char *name)
{
  int l = level ? level->l : -1;
  int i;

  for (i = 0; i < ntimers; i++)
    if (l == tdb[i].l && !strcmp(tdb[i].name, name))
      return tdb+i;
  tdb = realloc(tdb, (++ntimers) * sizeof(tTimer));
  tdb[i].name  = strdup(name);
  tdb[i].l     = l;
  tdb[i].start = -1;
  tdb[i].time  = 0;
  tdb[i].n     = 0;

  return tdb+i;
}




/* start timer */
int timer_start(tL *level, char *name)
{
  if (timer_on) {
    tTimer *t;

    if (timer_on == -1) {
      timer_on = Getv("bampi_timer_on", "yes");
      if (!timer_on) return 0;
      timer_barrier = Getv("bampi_timer_barrier", "yes");
    }

    t = timer_get(level, name);
    if (timer_barrier) bampi_barrier();
    t->start = bampi_time();

    if (0) printf("starting timer %s, %.15e\n", t->name, t->start);
  }
  return 1;
}




/* stop timer */
int timer_stop(tL *level, char *name)
{
  if (timer_on) {
    tTimer *t = timer_get(level, name);

    if (t->start < 0) return 0;
    if (timer_barrier) bampi_barrier();
    t->n += 1;
    t->time += bampi_time() - t->start;
    t->start = -1;
    
    if (0) printf("stopping timer %s, %e\n", t->name, t->time);
  }  
  return 1;
}




/* timer comparison for sorting */
int timer_compare(const void *a, const void *b)
{
  tTimer *t = (tTimer *) a;
  tTimer *u = (tTimer *) b;
  int i;

  if (t->l == -1 && u->l >= 0) return  1;
  if (t->l >= 0 && u->l == -1) return -1;

  if (t->l == -1 && u->l == -1) {
    if (t->time < u->time) return -1;
    if (t->time > u->time) return  1;
    return 0;
  }

  i = strcmp(t->name, u->name);
  if (!i) {
    if (t->l < u->l) return -1;
    if (t->l > u->l) return  1;
  }

  return i;
}




/* print timer info */
void timer_print(tG *g)
{
  FILE *fp;
  char *outdir = Gets("outdir");
  char f[100], s[1000];
  tTimer *t;
  double p, total;
  int i;

  if (!timer_on) return;

  /* open file */
  snprintf(f, 100, "%%s/timer.%%0%dd", (int) log10(bampi_size())+1);
  snprintf(s, 1000, f, outdir, bampi_rank());
  fp = fopen(s, "a");
  if (!fp) errorexits("could not open %s", s);
  fprintf(fp, "\nTimers after top level iteration %d, time %.3f\n",
	  g->level[0]->iteration, g->level[0]->time);
  
  /* update timer for main() */
  timer_stop(0, "main");
  timer_start(0, "main");
  t = timer_get(0, "main");
  total = t->time;
  t->n = 1;

  /* add up timers for different levels */
  for (i = 0; i < ntimers; i++) {
    if (tdb[i].l >= 0) {
      t = timer_get(0, tdb[i].name);
      t->time += tdb[i].time;
      t->n    += tdb[i].n;
    }
  }

  /* sort */
  if (1) qsort(tdb, ntimers, sizeof(tTimer), timer_compare);

  /* print */
  for (i = 0; i < ntimers; i++) {
    t = tdb+i;
    p = t->time/total * 100;

    if (t->l < 0) 
      fprintf(fp, "%-24s %5.1f  %9.2f %7d\n",
	      t->name, p, t->time, t->n);
 
    if (t->l >= 0)
      fprintf(fp, "%2d %-21s        %9.2f %7d\n",
	      t->l, t->name, t->time, t->n);
  }

  /* reset level specific timers */
  for (i = 0; i < ntimers; i++) {
    if (Getv("bampi_timer_reset_after_every_iteration", "yes") ||  (tdb[i].l >= 0)) {
      tdb[i].time = 0;
      tdb[i].n    = 0;
    }
  }

  /* done */
  fclose(fp);
}



void freeTimer()
{
  int i;
  for (i=0; i<ntimers; i++) {
    if (0) printf("delete %d %s\n", i,tdb[i].name);
    free(tdb[i].name);
  }
  if (0) printf("%d\n", ntimers);
  free(tdb);
}




/*
void timing_analyze_level(tG *g, int l)
{
  double dt = bampi_delta_time();

#undef analyze_level
  analyze_level(g, l);

  dt = bampi_delta_time();

  printf("analyze_level  dt = %e\n", dt);
}
*/

