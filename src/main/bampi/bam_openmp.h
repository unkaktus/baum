/* bam_openmp.h */
/* mth 4/2012 */






#ifdef OPENMP

  #define bampi_openmp_start \
  _Pragma("omp parallel")\
  { 
  
  #define bampi_openmp_stop \
  _Pragma("omp barrier") \
  }

  #define bampi_openmp_loop\
  _Pragma("omp for schedule(static)")

  #define bampi_openmp_loop_collapse2\
  _Pragma("omp for collapse(2) schedule(static)")

  #define bampi_openmp_loop_collapse3\
  _Pragma("omp for collapse(3) schedule(static)")


  #define bampi_openmp_parallel_for\
  _Pragma("omp parallel for schedule(static)")

  #define bampi_openmp_parallel_for_collapse2\
  _Pragma("omp parallel for collapse(2) schedule(static)")

  #define bampi_openmp_parallel_for_collapse3\
  _Pragma("omp parallel for collapse(3) schedule(static)")


  #define bampi_openmp_barrier \
  _Pragma("omp barrier")


#else


  #define bampi_openmp_start \
  {

  #define bampi_openmp_stop \
  }

  #define bampi_openmp_loop \

  #define bampi_openmp_loop_collapse2
  #define bampi_openmp_loop_collapse3

  #define bampi_openmp_parallel_for
  #define bampi_openmp_parallel_for_collapse2
  #define bampi_openmp_parallel_for_collapse3

  #define bampi_openmp_barrier \


#endif












