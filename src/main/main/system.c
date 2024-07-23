/* system.c */
/* Bernd Bruegmann, 4/2006 */
/* Ivan Markin, 5/2024 */

/* monitor system status
   On Linux, it relies on /proc directory
   On macOS, it relies on system calls
*/

#include "bam.h"
#include "main.h"

// This function works all Unixes even without procfs
void read_loadavg(double *loadavg1, double *loadavg5, double *loadavg15) {
  double loadavg[3];
  if (!getloadavg(loadavg, 3)) {
    return;
  }
  *loadavg1 = loadavg[0];
  *loadavg5 = loadavg[1];
  *loadavg15 = loadavg[2];
}


#ifdef __linux__
  /* return content of file in string
    for now: read first N characters only
  */
  char *read_entire_file(char *name)
  {
    int pr = Getv("system_monitor_verbose", "yes");
    FILE *fp = fopen(name, "r");
    char *s;
    int n;

    if (!fp) {
      if (pr) printf("system monitor: could not open %s\n", name);
      return 0;
    }

    s = bcalloc(10001);
    n = fread(s, 1, 10000, fp);
    s[n] = 0;
    fclose(fp);

    if (pr) {
      printf("system monitor: read %d bytes from %s, got\n", n, name);
      printf(">>>>>>>\n%s<<<<<<<\n", s);
    }
    return s;
  }

  /* read one value given a keyword assuming lines like
      keyword:  multi word string
  */
  char *read_keyword_value(char *s, char *keyword, double *value)
  {
    char *t, *u, c;
    char *svalue;

    t = strstr(s, keyword);
    if (!t) {
      *value = 0;
      return 0;
    }
    for (t = t + strlen(keyword); *t == ' ' || *t == '\t'; t++);
    for (u = t; *u != 0 && *u != '\n' && *u != '\r'; u++);
    c = *u;
    *u = 0;
    svalue = strdup(t);
    *value = atof(svalue);
    *u = c;

    return svalue;
  }


  void read_meminfo(double *VmSize, double *VmRSS) {
    char *s = read_entire_file("/proc/self/status");
    if (s) {
      char * sVmSize = read_keyword_value(s, "VmSize:", VmSize);
      char * sVmRSS  = read_keyword_value(s, "VmRSS:", VmRSS);
      free(sVmSize);
      free(sVmRSS);
      free(s);
    }
  }

#endif


#ifdef __APPLE__

  #include <mach/mach.h>

  void read_meminfo(double *VmSize, double *VmRSS) {
    kern_return_t ret;
    mach_task_basic_info_data_t info;
    mach_msg_type_number_t count = MACH_TASK_BASIC_INFO_COUNT;

    ret = task_info(mach_task_self(), MACH_TASK_BASIC_INFO, (task_info_t)&info, &count);
    if (ret != KERN_SUCCESS || count != MACH_TASK_BASIC_INFO_COUNT)
    {
        return;
    }

    // Return sizes in kB
    *VmSize = (double) info.virtual_size/1000;
    *VmRSS = (double) info.resident_size/1000;
  }

#endif

/* monitor system status and performance */
int system_monitor(tL *level)
{
  static int firstcall = 1;
  static double time0 = 0, time1 = 0, phystime = 0;
  double time2, dtime0, dtime1, dphystime;
  int pr = 1;
  FILE *fp;
  char *s;
  char out[1000], form[100];
  int i, n;
  double loadavg1 = 0, loadavg5 = 0, loadavg15 = 0;
  double VmSize = 0, VmRSS = 0;
  double v[10], vall[10], vmin[10], vmax[10];

  i = Geti("system_monitor_every");
  if (i <= 0 || 
      level->l != level->grid->lmin ||
      level->iteration % i != 0)
    return 0;

  /* do this once */
  if (firstcall) {
    /* time of first call */
    time0 = bampi_time();
    phystime = level->time;
  }


  /* system load */
  read_loadavg(&loadavg1, &loadavg5, &loadavg15);

  /* process status */
  read_meminfo(&VmSize, &VmRSS);

  /* get the sum or average, minimal, maximal value */
  n = 0;
  v[n++] = loadavg5;
  v[n++] = VmSize/1000000;
  v[n++] = VmRSS /1000000;
  bampi_allreduce_sum_vector(v, vall, n);
  bampi_allreduce_max_vector(v, vmax, n);
  bampi_allreduce_min_vector(v, vmin, n);
  vall[0] /= bampi_size();
  vall[3] /= bampi_size();

  /* timing information */
  time2 = bampi_time();
  dtime1 = time2 - time1;
  dtime0 = time2 - time0;
  time1 = time2;
  dphystime = level->time - phystime;
  phystime = level->time;

  /* print information */
  if (processor0) {
    snprintf(out, 1000, "%s/system.log", Gets("outdir"));
    fp = fopen(out, "a");
    if (fp) {
      if (firstcall) 
	fprintf(fp, " PhysTime WallTime Phys/Wall RecentPh/Wa   "
		"  LoadAvg          VmSize           VmRSS\n");
      
      fprintf(fp, "%7.1fM %7.1fh %7.1fM/h %7.1fM/h",
	      level->time, dtime0/3600,
	      level->time/(dtime0/3600+1e-12),
	      dphystime/(dtime1/3600+1e-12));
      
      fprintf(fp, "   %3.1f (%3.1f %3.1f)", vall[0], vmin[0], vmax[0]);
      fprintf(fp, "   %3.1fg (%3.1f %3.1f)", vall[1], vmin[1], vmax[1]);
      fprintf(fp, "   %3.1fg (%3.1f %3.1f)", vall[2], vmin[2], vmax[2]);

      fprintf(fp, "\n");
      fclose(fp);
    }
  }
  
  if (Getv("stdout_verbose","speed")) {
    printf("speed:   %7.1fh %7.1fM/h\n",level->time/(dtime0/3600+1e-12), dphystime/(dtime1/3600+1e-12));
  }

  /* done */
  
  firstcall = 0;
  return 0;
}
