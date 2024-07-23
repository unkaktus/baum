/* alpha_center.c */
/* Bernd Bruegmann, 12/2005, 3/2007 */

#include "bam.h"
#include "Gauge.h"


/* write value of lapse between punctures
   (compare BrillWave/lapse_origin.c)
   a lapse of less than 0.3 indicates roughly the merger of the 
   apparent horizons 
*/ 
int write_lapse_between_punctures(tL *level)
{
  double x = 0, y = 0, z = 0;
  FILE *fp;
  char *outdir = Gets("outdir");
  char *name = "alpha_center";
  int n = strlen(outdir) + strlen(name) + 200;
  char filename[10000];
  double alpha0;
  
  /* how often do we want to output anyway? */
  if (!timeforoutput_di_dt(level, Geti("0doutiter"), Getd("0douttime")) &&
      !timeforoutput_di_dt(level, Geti("1doutiter"), Getd("1douttime")))
    return 0;

  /* figure out what point we want to watch
     for equal masses it is half way between the punctures
     for quadrant and octant it is the origin
     unfinished ...
  */
  if (!Getv("grid", "quadrant") && !Getv("grid", "octant"))
    return 0;

  /* in principle it should be helpful to only output the lapse
     on the finest level that covers the point
  */
  {
    tG *g = level->grid;
    int flag, l;
    
    for (l = g->lmin; l <= g->lmax; l++) {
      flag = 0;
      forallboxes(g->level[l]) {
	if (xyzinsidebbox(box->bbox, x, y, z)) flag = 1;
      } endforboxes;
      if (!flag) break;
    }
    l = l-1; // this is the finest level covering point
    l = bampi_allreduce_max_int(l);
    if (level->l != l) 
      return 0;
  }

  /* get lapse */
  alpha0 = interpolate_xyz_scalar(level, x, y, z, Ind("alpha"), 0);

  /* all processors have to call interpolator, but only one does the writing */
  if (bampi_rank() != 0) return 0;

  /* filename: we used to output into different files per level,
     but since the level determination turned out to be reliable, use just one file
  */
  if (0) 
    snprintf(filename, n, "%s/%s.t%d", outdir, name, level->l);
  else
    snprintf(filename, n, "%s/%s.tl", outdir, name);

  /* write */
  fp = fopen(filename, "a");
  if (!fp) errorexits("failed opening %s", filename);

  fprintf(fp, "%22.15e %22.15e %d\n", level->time, alpha0, level->l);
  fclose(fp);


  /* determine and write merger time
     rule of thumb for merger of apparent horizons:
       lapse between black holes has dropped below ca. 0.3
     use parameters for compatibility with checkpointing
  */
  if (1) {
    double a1 = Getd("alpha_center_value");
    double a2 = alpha0;
    double am = Getd("alpha_center_merger");
    
    if (a1 > am && a2 <= am) {
      double t1 = Getd("alpha_center_time");
      double t2 = level->time;
      double tm;

      /* linear interpolation to determine time for alpha = alpha_merger */
      tm = t1 + (t2 - t1) * (am - a1) / (a2 - a1);

      /* write, the file contains the latest up/down crossing of threshold */
      snprintf(filename, n, "%s/%s_merger.tl", outdir, name);
      fp = fopen(filename, "w");
      if (!fp) errorexits("failed opening %s", filename);
      fprintf(fp, "%9.3f  %6.3f\n", tm, am);
      fclose (fp);
    }
    
    Setd("alpha_center_value", alpha0);
    Setd("alpha_center_time",  level->time);
  }


  /* done */
  return 0;
}



