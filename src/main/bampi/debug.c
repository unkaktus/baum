/* debug.c */
/* Bernd Bruegmann 2/2006 */

/* note that we do not include any bam headers */
#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>

#define PR 0




/* print string for return value of MPI calls */
void prmpiresult(char *file, int line, int result, char *comment)
{
  int ns;
  char *s;

  if (!PR && result == MPI_SUCCESS) return;

  /* write same string on LAM and MPICH if there is no error */
  if (result == MPI_SUCCESS) {
    printf("%s: MPI SUCCESS ", comment);
  }

  /* write error string */
  else {
    s = calloc(sizeof(char), MPI_MAX_ERROR_STRING);
    MPI_Error_string(result, s, &ns);
    printf("%s: MPI ERROR, %s ", comment, s);
    free(s);
  }

  printf("(%s, line %d)\n", file, line);
  fflush(stdout);
}
#define prmpiresult(r,c) prmpiresult(__FILE__, __LINE__, (r), (c))




/* debug: wrap all MPI calls and compute checksums 

   there are facilities for LAM-MPI and MPICH, but they do not seem
   to be uniform across machines, and we want to do cross machine comparisons

   the surprising thing is that I didn't feel compelled to attempt such
   debugging until 2/2006, when for the first time running tests on single
   processor machines with say 8 processes was not able to reproduce the 
   issue on our linux cluster ...
*/
void debug_checksum(void *buf, int count, MPI_Datatype dtype, char *comment)
{
  int i;
  double sum = 0;
  static int n = 0;

  if (dtype == MPI_DOUBLE) {
    double *b = buf;
    for (i = 0; i < count; i++)
      sum += b[i];
  }

  if (dtype == MPI_INT) {
    int *b = buf;
    for (i = 0; i < count; i++)
      sum += b[i];
  }

  if (dtype == MPI_CHAR) {
    char *b = buf;
    for (i = 0; i < count; i++)
      sum += b[i];
  }

  printf("checksum %d %6d %23.16e  %s\n", n++, count, sum, comment);
  fflush(stdout);
}




int debug_MPI_Allreduce(void *sbuf, void *rbuf, int count,
		  MPI_Datatype dtype, MPI_Op op, MPI_Comm comm)
{
  int result = MPI_Allreduce(sbuf, rbuf, count, dtype, op, comm);
  char comment[1000];
  static int n = 0;

  snprintf(comment, 1000, "MPI_Allreduce %d", n++);
  debug_checksum(sbuf, count, dtype, comment);
  debug_checksum(rbuf, count, dtype, comment);
  prmpiresult(result, comment);
  return result;
}




int debug_MPI_Bcast(void *buff, int count, MPI_Datatype datatype,
	      int root, MPI_Comm comm)
{
  int result = MPI_Bcast(buff, count, datatype, root, comm);
  char comment[1000];
  static int n = 0;

  snprintf(comment, 1000, "MPI_Bcast %d, root %d", n++, root);
  debug_checksum(buff, count, datatype, comment);
  prmpiresult(result, comment);
  return result;
}




int debug_MPI_Irecv(void *buf, int count, MPI_Datatype dtype,
		    int src, int tag, MPI_Comm comm,
		    MPI_Request *req)
{
  int result = MPI_Irecv(buf, count, dtype, src, tag, comm, req);
  char comment[1000];
  static int n = 0;

  snprintf(comment, 1000, "MPI_Irecv %d, src %d, tag %d", n++, src, tag);
  // no sense checking since the actual receive probably has not yet happened
  debug_checksum(buf, 0, dtype, comment);
  prmpiresult(result, comment);
  return result;
}




int debug_MPI_Isend(void *buf, int count, MPI_Datatype dtype,
		    int dest, int tag, MPI_Comm comm,
		    MPI_Request *req)
{
  int result = MPI_Isend(buf, count, dtype, dest, tag, comm, req);
  char comment[1000];
  static int n = 0;

  snprintf(comment, 1000, "MPI_Isend %d, dst %d, tag %d", n++, dest, tag);
  debug_checksum(buf, count, dtype, comment);
  prmpiresult(result, comment);
  return result;
}




int debug_MPI_Recv(void *buf, int count, MPI_Datatype dtype,
		   int src, int tag, MPI_Comm comm, MPI_Status *stat)
{
  int result = MPI_Recv(buf, count, dtype, src, tag, comm, stat);
  char comment[1000];
  static int n = 0;

  snprintf(comment, 1000, "MPI_Recv %d, src %d, tag %d", n++, src, tag);
  debug_checksum(buf, count, dtype, comment);
  prmpiresult(result, comment);
  return result;
}




int debug_MPI_Reduce(void *sbuf, void *rbuf, int count, 
		     MPI_Datatype dtype, MPI_Op op, int root,
		     MPI_Comm comm)
{
  int result = MPI_Reduce(sbuf, rbuf, count, dtype, op, root, comm);
  char comment[1000];
  int rank;
  static int n = 0;

  if (dtype == MPI_INT)
    snprintf(comment, 1000, "MPI_Reduce INT %d, root %d", n++, root);
  else if (dtype == MPI_DOUBLE)
    snprintf(comment, 1000, "MPI_Reduce DBL %d, root %d", n++, root);
  else
    snprintf(comment, 1000, "MPI_Reduce xxx %d, root %d", n++, root);

  debug_checksum(sbuf, count, dtype, comment);

  MPI_Comm_rank(comm, &rank);
  if (root == rank) 
    debug_checksum(rbuf, count, dtype, comment);

  prmpiresult(result, comment);

  //if (n >= 5) errorexit("testing");
  return result;
}




int debug_MPI_Send(void *buf, int count, MPI_Datatype dtype,
		    int dest, int tag, MPI_Comm comm)
{
  int result = MPI_Send(buf, count, dtype, dest, tag, comm);
  char comment[1000];
  static int n = 0;

  snprintf(comment, 1000, "MPI_Send %d, dst %d, tag %d", n++, dest, tag);
  debug_checksum(buf, count, dtype, comment);
  prmpiresult(result, comment);
  return result;
}










