/* synchronize.c */
/* Bernd Bruegmann 5/02 */

/* Our synchronization strategy is based on the example 
   for Jacobi iterations with nonblocking recieves and sends found at
   http://www-unix.mcs.anl.gov/mpi/tutorial/mpiexmpl/contents.html
*/

#include "bam.h"
#include "bampi.h"

#define PR 0


/* fill a communication buffer with data from variables */
void fillbuffer(tL *l, tComBuf *cb, int nv, int *iv)
{
  int ni = cb->ni;
  int *i = cb->i;
  double *b = cb->buffer;
  double **v = l->v;

  bampi_openmp_parallel_for_collapse2
  for (int n = 0; n < nv; n++) {
    for (int j = 0; j < ni; j++) {
      int k = n * ni + j;
      double *u = v[iv[n]];
      b[k] = u[i[j]];
    }
  }
}

/* fill all buffers for each of the 6 directions with data from variables */
void fillbuffers(tL *l, int nv, int *iv)
{
  int i;
  
  for (i = 0; i < 6; i++) 
    if (l->com->nbrank[i] >= 0)
      fillbuffer(l, l->com->send[i], nv, iv);
}




/* copy data from a communication buffer into variables */
void readbuffer(tL *l, tComBuf *cb, int nv, int *iv)
{
  int ni = cb->ni;
  int *i = cb->i;
  double *b = cb->buffer;
  double **v = l->v;

  bampi_openmp_parallel_for_collapse2
  for (int n = 0; n < nv; n++) {
    for (int j = 0; j < ni; j++) {
      int k = n * ni + j;
      double *u = v[iv[n]];
      u[i[j]] = b[k];
    }
  }
}

/* copy data from all buffers for each of the 6 directions into variables */
void readbuffers(tL *l, int nv, int *iv)
{
  int i;
  
  for (i = 0; i < 6; i++) 
    if (l->com->nbrank[i] >= 0)
      readbuffer(l, l->com->recv[i], nv, iv);
}




/* synchronize buffers */
/* there are many possibilities to experiment with
   - let's start with a simple non-blocking scheme
   - first post two receives for opposing faces, 
     then the two sends, then wait
   - this way edges and corners are propagated around correctly,
*/
void nonblockrecv(tCom *c, int d, MPI_Request *request, int *nrequests)
{
  tComBuf *b;
  int srcrank;
  int result;
  int tag;

  srcrank = c->nbrank[d];
  if (srcrank >= 0) {
    (*nrequests)++;
    b = c->recv[d];
    tag = d;
    
    if (PR) printf("%d: Irecv from %d, nbuffer=%d, tag=%d\n", 
		  srcrank, c->myrank, b->nbuffer, tag);
    
    result = MPI_Irecv(b->buffer, b->nbuffer, MPI_DOUBLE, srcrank, tag, 
		       MPI_COMM_WORLD, request);

    prmpiresult(result, "Irecv");
  }
}

void nonblocksend(tCom *c, int d, MPI_Request *request, int *nrequests)
{
  tComBuf *b;
  int dstrank;
  int result;
  int tag;
  int ns;
  char s[10000];

  dstrank = c->nbrank[d];
  if (dstrank >= 0) {
    (*nrequests)++;
    b = c->send[d];
    tag = (d/2)*2 + 1 - d%2;  /* inverse direction to match recv */
    
    if (PR) printf("%d: Isend from %d, nbuffer=%d, tag=%d\n", 
		  dstrank, c->myrank, b->nbuffer, tag);
    
    result = MPI_Isend(b->buffer, b->nbuffer, MPI_DOUBLE, dstrank, tag, 
		       MPI_COMM_WORLD, request);

    prmpiresult(result, "Isend");
  }
}




/* synchronize variables using the helpers above */
void bampi_synchronize_work(tL *l, int nv, int *iv)
{
  tCom *c = l->com;
  int nrequests;
  MPI_Request request[4];
  MPI_Status status[4];
  int result;
  int d;

  if (bampi_size() == 1) return;

  if (PR) {
    if (0) prdivider(0);
    printf("Synchronizing %d variables %s ... on l%d\n",
	   nv, VarName(iv[0]), l->l);
  }

  /* make sure we have all the storage we need for all directions */
  for (d = 0; d < 6; d++)
    if (c->send[d]->ni &&
	c->send[d]->nbuffer != nv * c->send[d]->ni) {
      free(c->send[d]->buffer);
      free(c->recv[d]->buffer);
      c->send[d]->nbuffer = nv * c->send[d]->ni;
      c->recv[d]->nbuffer = nv * c->recv[d]->ni;
      c->send[d]->buffer = dmalloc(c->send[d]->nbuffer);
      c->recv[d]->buffer = dmalloc(c->recv[d]->nbuffer);
    }

  /* for 3 pairs of opposing faces */
  for (d = 0; d < 6; d += 2) {

    /* this variable counts the requests since the number of requests
       varies with the location of the processor in the processor grid
       (fewer requests at boundary)
    */ 
    nrequests = 0;

    /* receive */
    nonblockrecv(c, d,   &request[nrequests], &nrequests);
    nonblockrecv(c, d+1, &request[nrequests], &nrequests);

    /* fill send buffers */ 
    timer_start(l, "sync_fillbuffers");
    fillbuffer(l, c->send[d],   nv, iv);
    fillbuffer(l, c->send[d+1], nv, iv);
    timer_stop(l, "sync_fillbuffers");

    /* send */
    nonblocksend(c, d,   &request[nrequests], &nrequests);
    nonblocksend(c, d+1, &request[nrequests], &nrequests);

    /* wait until those two faces are dealt with */
    timer_start(l, "sync_waitall");
    result = MPI_Waitall(nrequests, request, status);
    prmpiresult(result, "Waitall");
    timer_stop(l, "sync_waitall");

    /* read receive buffers */
    timer_start(l, "sync_readbuffers");
    readbuffer(l, c->recv[d],   nv, iv);
    readbuffer(l, c->recv[d+1], nv, iv);
    timer_stop(l, "sync_readbuffers");
  }

  if (PR) {
    printf("Synchronizing finished on l%d\n", l->l);
    prdivider(0);
  }
}




/* wrapper for efficiency:
   depending on network speed, it sometimes is faster to collect
   all communication into one buffer (for high latency), or to send
   many small messages (for low latency to overlap with copy)
*/
void bampi_synchronize_variables(tL *l, int nv, int *iv)
{
  int i;
  timer_start(l, "synchronize");

  if (Getv("bampi_lowlatency", "no")) {
    bampi_synchronize_work(l, nv, iv);
  } else {
    for (i = 0; i < nv; i++)
      bampi_synchronize_work(l, 1, iv+i);
  }

  timer_stop(l, "synchronize");
}




/* synchronize all components of a variable given by its first component */
void bampi_synchronize(tL *level, int vi)
{
  tVarList *vl = vlalloc(level);

  vlpush(vl, vi);
  bampi_vlsynchronize(vl);
  vlfree(vl);
}




/* synchronize a variable list */
void bampi_vlsynchronize(tVarList *vl)
{
  if (0) printf("vlsynchronize %d on l%d, %s to %s\n",
		vl->n, vl->level->l,
		VarName(vl->index[0]), VarName(vl->index[vl->n-1]));
  if (0) prvarlist(vl);


  int nbox;
  tL *obl;

  for (nbox = 0; nbox < vl->level->nboxes; nbox++) {
    obl = one_box_level(vl->level, nbox);
    bampi_synchronize_variables(obl, vl->n, vl->index);
  }
  
  sync_shells(vl->level, vl->n, vl->index);
}

