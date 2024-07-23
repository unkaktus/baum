/* bam_profile.h */
/* Bernd Bruegmann 4/2006 */

/* Exclude certain functions from profiling with gcc -finstrument-functions
   Needs more testing.
   Did we get the order of #includes right? How to attributes really work?

   Some functions do little work, are called millions of times, and the
   instrumentation overhead distorts the timing of the calling function.
   Such functions can be excluded here.

   To exclude entire files and projects, compile selectively with and without
   -finstrument-functions. For example, compile bam with instrumentation,
   then touch/update a file or project and do a make with instrumentation
   turned on in MyConfig. (Has to be tested again with the latest profile.c.)

   Another reason to exclude a function from profiling is that the 
   current profile.c does not handle recursive functions.
*/


#ifdef PROFILING

/* recursive functions */
void advance_levelandsublevels(tG *g, int l)
  __attribute__ ((no_instrument_function));

/* light-weight functions that are called so often that they cloud the issue */
int xyzinsidebbox(double *bbox, double x, double y, double z)
  __attribute__ ((no_instrument_function));

double interpolate_tripolynomial(double, double, double, double, double,
  double, double, double, double)
  __attribute__ ((no_instrument_function));

double interpolate_tri6(double, double, double, double, double,
  double, double, double, double)
  __attribute__ ((no_instrument_function));

void coefficients_lagrange6(double x, double xmin, double h,
			    double *k0, double *k1, double *k2, 
			    double *k3, double *k4, double *k5)
  __attribute__ ((no_instrument_function));

void interpolate_filldatacube(tL *level, int vi)
  __attribute__ ((no_instrument_function));

int find_one_point_box(tL *level, double x, double y, double z)
  __attribute__ ((no_instrument_function));

int insidebbox(double *bbox, double *coord)
  __attribute__ ((no_instrument_function));

#endif
