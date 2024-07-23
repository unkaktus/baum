/* bampi.h */
/* Bernd Bruegmann, 5/02 */

#include "mpi.h"

#ifdef OPENMP
#include "omp.h"
#endif



/* debug.c */
void prmpiresult(char *file, int line, int result, char *comment);
#define prmpiresult(r,c) prmpiresult(__FILE__, __LINE__, (r), (c))


#if 0
#define MPI_Allreduce debug_MPI_Allreduce
#define MPI_Bcast debug_MPI_Bcast
#define MPI_Irecv debug_MPI_Irecv
#define MPI_Isend debug_MPI_Isend
#define MPI_Recv debug_MPI_Recv
#define MPI_Reduce debug_MPI_Reduce
#define MPI_Send debug_MPI_Send
#endif

#if 0
#define MPI_Recv debug_MPI_Recv
#define MPI_Send debug_MPI_Send
#endif


/* profile.c */
void bampi_profile(char *exename);
