/* openmp.c */
/* mth 2/12 */

#include "bam.h"
#include "bampi.h"

#define PR 0


#ifdef OPENMP

/* detailed info about openmp */
void bampi_openmp_pr_info()
{

  int nt, mnt, tn;

  nt = omp_get_num_threads();
  mnt = omp_get_max_threads();
  tn = omp_get_thread_num();
  printf("CPU %d of %d, Thread %d of %d\n"
         "omp_get_num_threads()=%d omp_get_max_threads()=%d "
         "omp_get_thread_num()=%d\n",
         bampi_rank(), bampi_size(), bampi_openmp_rank(), bampi_openmp_size(),
         nt, mnt, tn);
  fflush (stdout);
  system("echo OMP_NUM_THREADS=$OMP_NUM_THREADS");
  fflush (stdout);
}


void bampi_openmp_initialize(int *pargc, char ***pargv, int pr) 
{
  int i,j;
  if (processor0 && pr)
    prdivider(0);
  fflush (stdout); 
  bampi_barrier();
  
  omp_set_num_threads(1);
  omp_set_dynamic(0);
  for (i = 0; i < *pargc; i++) {
    if (strcmp((*pargv)[i],"-nt")==0) {
      if (i+1<*pargc) {
        omp_set_num_threads(atoi((*pargv)[i+1]));
        
        for (j=i; j<*pargc-2; j++)
          (*pargv)[j] = (*pargv)[j+2];
        
        *pargc -= 2;
      } else {
        errorexit("wrong openmp syntax, use -nt NUMBER");
      }
    }
  }
  
  
  for (i=0; i<bampi_size(); i++) {
    if (bampi_rank()==i && pr)
    {
      printf("Initializing OPENMP on cpu %d of %d with %d Threads \n",
             bampi_rank(),bampi_size(),bampi_openmp_size());
      bampi_openmp_pr_info();
    }
    fflush (stdout); 
    bampi_barrier();
  }
}

int bampi_openmp_rank() 
{
  return omp_get_thread_num();
}

int bampi_openmp_size() 
{
  if (omp_get_num_threads()>1) return omp_get_num_threads();
  return omp_get_max_threads();
}

#endif
