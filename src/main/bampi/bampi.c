/* bampi.c */
/* Bernd Bruegmann 5/02 */

#include "bam.h"
#include "bampi.h"

#define PR 0




/* standard MPI start-up */
void bampi_initialize(int *pargc, char ***pargv) 
{
  MPI_Init(pargc, pargv);
}




/* standard MPI shut-down */
void bampi_finalize(int argc, char **argv)
{
  if (argv) bampi_profile(argv[0]);

  printf("Thank you for running b a m.\n");
  if (0) printf("We know you had a pleasant run because you got this far.\n");
  fclose(stderr);
  fclose(stdout);
  MPI_Finalize();
}




/* standard MPI abort */
void bampi_abort(void)
{
  MPI_Abort(MPI_COMM_WORLD, 1);
}




/* rank of our process, often used to single out processor 0 for output */
/* for convenience there is macro processor0 in bam_bampi.h for a rank 0 check
   because we do not want to use global variables to avoid lib dependencies
*/
int bampi_rank(void) 
{
  int rank;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  return rank;
}




/* number of processors we are running on */
int bampi_size(void) 
{
  int size;

  MPI_Comm_size(MPI_COMM_WORLD, &size);
  return size;
}




/* timing */
/* return MPI time */
double bampi_time(void)
{
  return MPI_Wtime();
}


/* return MPI time since the first call to this function
*/
double bampi_uptime(void)
{
  static double t0 = -1;
  double t, dt;

  t = MPI_Wtime();
  if (t0 < 0) t0 = t;
  dt = t - t0;
/*
  printf("initial MPI wall time = %f\n", t0);
  printf("current MPI wall time = %f\n", t); 
  printf("elapsed MPI wall time = %f\n", dt); 
*/
  return dt;
}




/* return MPI time since the last call to this function
   call repeatedly to time pieces of code
*/
double bampi_delta_time(void)
{
  static double t0 = -1;
  double dt;

  if (t0 < 0) t0 = MPI_Wtime();
  dt = MPI_Wtime() - t0;
  t0 += dt;
  return dt;
}




/* barrier */
void bampi_barrier(void)
{
  MPI_Barrier(MPI_COMM_WORLD);
}




/* broadcast */
void bampi_bcast_int(void *buffer, int count, int root)
{
  int result = MPI_Bcast(buffer, count, MPI_INT, root, MPI_COMM_WORLD);
  prmpiresult(result, "Bcast int");
}
void bampi_bcast_double(void *buffer, int count, int root)
{
  int result = MPI_Bcast(buffer, count, MPI_DOUBLE, root, MPI_COMM_WORLD);
  prmpiresult(result, "Bcast double");
}
void bampi_bcast_char(void *buffer, int count, int root)
{
  int result = MPI_Bcast(buffer, count, MPI_CHAR, root, MPI_COMM_WORLD);
  prmpiresult(result, "Bcast char");
}


void bampi_alltoall(void *bufsend, void *bufrecv, int count)
{
  int result = MPI_Allgather(bufsend, count, MPI_DOUBLE, 
                           bufrecv, count, MPI_DOUBLE, MPI_COMM_WORLD);
 
  prmpiresult(result, "AllToAll");
}







/* parallel or-operator
   allows code like
     if (something == happened) lets_all_react();
   good: all processes learn about it
   bad:  this is an all to all communication with overhead
   for exit upon error, use bampi_abort()
*/
int bampi_or(int flag)
{
  if (flag != 0) flag = 1;
  return bampi_allreduce_max_int(flag);
}




/* swap buffer content between processors, 
   e.g. for periodic and rotant boundaries

   does not work for unknown reason
*/
void bampi_swap_buffer(void *buffer, int count, int srcrank, int dstrank)
{
  MPI_Status status;
  int sendtag = 0*srcrank;
  int recvtag = 0*dstrank;
  int result;

  printf("dst %d  src %d  ", dstrank, srcrank); 

  if (0)
  result = MPI_Sendrecv_replace(buffer, count, MPI_DOUBLE, dstrank,
				    sendtag, srcrank, recvtag,
				    MPI_COMM_WORLD, &status); 

  if (0) prmpiresult(result, "Sendrecv_replace");
}




/* exchange buffers with non-blocking send and receive
*/
void bampi_isendrecv(void *bufsend, void *bufrecv, int count, int rank)
{
  MPI_Status status[2];
  MPI_Request request[2];
  int result;

  if (0) {printf("myrank %d  hisrank %d  ", bampi_rank(), rank);}

  result = MPI_Irecv(bufrecv, count, MPI_DOUBLE, rank, 0,
		     MPI_COMM_WORLD, &request[1] );
  prmpiresult(result, "Irecv");

  result = MPI_Isend(bufsend, count, MPI_DOUBLE, rank, 0,
		     MPI_COMM_WORLD, &request[0] );
  prmpiresult(result, "Isend");

  result = MPI_Waitall(2, request, status);
  prmpiresult(result, "Waitall");
}




/* allocate one communication buffer */
tComBuf *bampi_alloc_combuf(void)
{
  tComBuf *b = calloc(1, sizeof(tComBuf));

  if (!b) errorexit("bampi_alloc_combuf: out of memory");
  return b;
}




/* free one communication buffer */
void bampi_free_combuf(tComBuf *b)
{
  if (!b) return;
  if (b->i) free(b->i);
  if (b->buffer) free(b->buffer);
  free(b);
}




/* allocate list of communication buffers */
tComBuf **bampi_alloc_combufs(int n)
{
  int i;
  tComBuf **b = malloc(n * sizeof(tComBuf *));

  if (!b) errorexit("bampi_alloc_combufs: out of memory");
  for (i = 0; i < n; i++)
    b[i] = bampi_alloc_combuf();
  return b;
}




/* free list of communication buffers */
void bampi_free_combufs(int n, tComBuf **b)
{
  int i;

  if (!b) return;
  for (i = 0; i < n; i++)
    bampi_free_combuf(b[i]);
  free(b);
}




/* allocate combuf list */
void bampi_alloc_combuflist(tL *level)
{
  int i;

  if (!level) return;
  if (!level->com) return;
  if (level->com->buflist) {
    if (0) errorexit("com->buflist should not be overwritten");
    if (0) printf("freeing com->buflist to avoid clobbering\n");
    bampi_free_combuflist(level);
  }
  else
    if (0) printf("NOT freeing com->buflist to avoid clobbering\n");
    

  level->com->buflist = pmalloc(5);
  for (i = 0; i < 5; i++)
    level->com->buflist[i] = bampi_alloc_combufs(bampi_size());
}




/* disable all lists of communication buffers (free and set null) */
void bampi_free_combuflist(tL *level)
{
  if (!level) return;
  if (!level->com) return;
  bampi_free_combuflist_com(level->com);
  level->com->buflist = 0;
}




/* free lists of communication buffers */
void bampi_free_combuflist_com(tCom *c)
{
  int i;
  tComBuf ***buflist;
  
  if (!c) return;
  buflist = c->buflist;
  if (!buflist) return;

  for (i = 0; i < 5; i++)
    bampi_free_combufs(bampi_size(), buflist[i]);
  free(buflist);
  c->buflist = 0;
}




/* allocate tCom * structure
   called from alloc_level() in amr/grid.c */
tCom *bampi_alloc_com(void)
{
  tCom *c;
  tComBuf *b;
  int i, s;

  c = (tCom *) calloc(sizeof(tCom), 1);
  c->myrank = bampi_rank();
  for (i = 0; i < 6; i++)
    c->nbrank[i] = -1;

  s = sizeof(tComBuf);
  b = (tComBuf *) calloc(s, 12);
  for (i = 0; i < 6; i++) {
    c->send[i] = b + 2*i;
    c->recv[i] = b + 2*i+1;
  }

  return c;
}




/* free communication data */
void bampi_free_com(tCom *c)
{
  int i;

  if (!c) return;

  bampi_free_combuflist_com(c);
  
  for (i = 0; i < 6; i++) {
    free(c->send[i]->i);
    free(c->send[i]->buffer);
    free(c->recv[i]->i);
    free(c->recv[i]->buffer);
  }
 
  free(c->send[0]);

  free(c);
}




/* set boundary flags for ghosts */
void bampi_set_boundary_flags(tL *l)
{
  int d, i, n, *ghostindex;

  for (d = 0; d < 6; d++) 
    if (l->com->nbrank[d] >= 0) {

      n = l->com->recv[d]->ni;
      ghostindex = l->com->recv[d]->i;

      if (PR) printf("d %d  n %d  ghostindex %p\n", d, n, ghostindex);
 
      for (i = 0; i < n; i++) 
	boundaryflag(l, ghostindex[i]) = GHOBOUND;
    }
}


