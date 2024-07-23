/* synchronize_fmr.c */
/* Bernd Bruegmann 5/04 */

/* For FMR, with as of now unknown overhead for AMR, the following
   strategy for parallelization is implemented:
   - store local parent buffer that includes all points for P and R operations 
   - P and R only refer to local parent buffer
   - synchronize with global parent as needed
   Concretely,
     level->prlevel points to local parent buffer (when running in parallel)
     grid->level[level->l-1] points to local part of global parent

   Optimize:
   - only restrict as many points as are needed for evolution, but
     offer full restrict for output and multigrid
   - include intra-level ghostzones
*/

#include "bam.h"
#include "bampi.h"




/* fill a communication buffer with data from a variable list */
/* clean up together with buffer routines in synchronize.c */
void vlfillbuffer(tL *level, tComBuf *cb, tVarList *vl)
{
  int *i = cb->i;
  double *b = cb->buffer;

  bampi_openmp_parallel_for_collapse2
  for (int vi = 0; vi < vl->n; vi++) {
    for (int j = 0; j < cb->ni; j++) {
      double *v = level->v[vl->index[vi]];
      int k = vi * (cb->ni) + j;
      b[k] = v[i[j]];
    }
  }
}




/* read data from a communication buffer into a variable list */
void vlreadbuffer(tL *level, tComBuf *cb, tVarList *vl)
{
  int *i = cb->i;
  double *b = cb->buffer;

  bampi_openmp_parallel_for_collapse2
  for (int vi = 0; vi < vl->n; vi++) {
    for (int j = 0; j < cb->ni; j++) {
      double *v = level->v[vl->index[vi]];
      int k = vi * (cb->ni) + j;
      v[i[j]] = b[k];
    }
  }
}




/* print info about com buffer */
void printcombuf(tComBuf *b, char *name, int verbose)
{
  int i, j, nj;
  
  /* long */
  if (verbose) {
    printf("combuf %s\n", name);
    if (!b) return;

    for (i = 0; i < b->ni; i++) {
      printf("%4d ", b->i[i]);
      if (b->buffer) {
	nj = b->nbuffer/b->ni;
	for (j = 0; j < nj; j++)
	  printf("%6.3f ", b->buffer[j+i*nj]);
	printf("\n");
      }
    }
    printf("\n");
  }

  /* short */
  else {
    printf("combuf %s ", name);
    if (!b) {printf("\n"); return;}
    printf("ni=%d\n", b->ni);
  }
}




/* print info about list of combufs */
void printcombufs(tComBuf **b, char *name, int verbose)
{
  int i, n = bampi_size();
  char s[1000];

  for (i = 0; i < n; i++) {
    snprintf(s, 1000, "p%d %s[%d]", bampi_rank(), name, i);
    printcombuf(b[i], s, verbose);
  }
}




/****************************************************************************/
/* precompute com bufs that are needed for synchronization of data */

/* initialize local to global given point list in buffer b */
void init_buflist(tL *prlocal,  tComBuf **blocal,
		  tL *prglobal, tComBuf **bglobal, tComBuf *b)
{
  double *coordbuf;
  int *ibuf, *iwho;
  int i, j, k, n, p, q;
  int me = bampi_rank();
  int size = bampi_size();

  /* each processor gets its turn to ask who has the points in b */
  for (p = 0; p < size; p++) {

    /* broadcast list of points for processor p */
    n = b->ni;

    bampi_bcast_int(&n, 1, p); 
    coordbuf = dmalloc(3*n);
    if (p == me)
      for (i = 0; i < 3*n; i++) coordbuf[i] = b->buffer[i]; 

    bampi_bcast_double(p==me ? b->buffer : coordbuf, 3*n, p);

    // shouldn't that be:
    // bampi_bcast_double(coordbuf, 3*n, p);
    if (0) printf("me%d, p%d: %d points\n", me, p, n);
    if (0) printlevel(prglobal);
    if (0) printlevel(prlocal);

    /* each processor, p included, looks for the points on global level */
    ibuf = imalloc(n);
    if (bglobal[p]->i) errorexit("clobbering 1");
    bglobal[p]->i = imalloc(n);

    for (i = k = 0; i < n; i++) {

      j = find_one_point_box(prglobal, 
			     coordbuf[3*i], coordbuf[3*i+1], coordbuf[3*i+2]);

      /* if the point is found on this processor */
      if (j >= 0 && 
	  boundaryflag(prglobal, j) != GHOBOUND &&
	  boundaryflag(prglobal, j) != SYMBOUND) {
	bglobal[p]->i[k++] = j;    /* remember where I found it */
	ibuf[i] = me;              /* remember who found it */
      } else
	ibuf[i] = 0;
    }

    /* this processor now has a list of indices for points sent by p */
    if (0) printf("k %d, n %d\n", k, n);
    bglobal[p]->ni = k;
    if (bglobal[p]->ni) {
      int *tmp = (int*) realloc(bglobal[p]->i, k*sizeof(int));
      if (!tmp) {
        tmp = (int*) malloc (k*sizeof(int));
        int n = (k<n)?k:n;
        for (i = 0; i < n; i++)
          tmp[i] = bglobal[p]->i[i];
        free(bglobal[p]->i);
      }
      bglobal[p]->i = tmp;
    } else {
      free(bglobal[p]->i);
      bglobal[p]->i = imalloc( k*sizeof(int));
    }
    
    /* let p know who found what */
    /* should check whether more than one processor found point ... */
    iwho = imalloc(n);
    MPI_Reduce(ibuf, iwho, n, MPI_INT, MPI_SUM, p, MPI_COMM_WORLD);
    if (0 && p == me) {
      for (i = 0; i < n; i++)
	printf("%5d %d %d\n", i, ibuf[i], iwho[i]);
    }

    /* if this is me (and there is more than one of us),
       store which points were found where 
    */
    if (p == me && blocal != bglobal) {

      /* count number of points per processor */
      for (i = 0; i < n; i++) {
	blocal[iwho[i]]->ni++;
	if (0) printf("iwho[%d] = %d\n", i, iwho[i]);
      }

      /* allocate memory for indices */
      for (q = 0; q < size; q++) {
	if (0) printf("q%d  ni %d\n", q, blocal[q]->ni);
	if (blocal[q]->i) printf("clobbering 2\n");
	blocal[q]->i = imalloc(blocal[q]->ni);
	blocal[q]->ni = 0;
      }

      /* store indices */
      for (i = 0; i < n; i++) {
	q = iwho[i];
	blocal[q]->i[blocal[q]->ni++] = b->i[i]; 
      }
    }

    /* clean up */
    free(iwho);
    free(ibuf);
    free(coordbuf);
  }

  if (0) {
    printcombufs(blocal, "blocal", 0);
    printcombufs(bglobal, "bglobal", 0);
  }
  
  /* check whether each processor sends the same number of
     point which are received by all the other */
  if (1) {
    double send[size], recv[size], s1,s2;
    int q;
    
    s1 = 0.;
    for (q = 0; q < size; q++) {
        s1 += (double)(blocal[q]->ni);
        send[q] = 0.;
    }
    
    for (p = 0; p < size; p++) {
        send[0] = (double)(bglobal[p]->ni);
        bampi_alltoall(send,recv,1);
        s2 = 0.;
        for (q = 0; q < size; q++) s2 += recv[q];
        
        if ((p==me) && (s1 != s2)) {
            printf("WARNING:  at least one processor is sending/receiving more than it receives/sends\n");
            printf("          this can lead to MPI problems -> maybe one box does not fit in an other box\n");
            printf(" %d %d %1.0f %1.0f    %1.0f %1.0f %1.0f %1.0f -> %1.0f %1.0f %1.0f %1.0f\n",p,me,s1,s2, send[0],send[1],send[2],send[3],recv[0],recv[1],recv[2],recv[3]);
        }
    }
    
  }
}




/* construct point list in com buffer b based on flags */
void init_localpointlist(tL *prlocal, tComBuf *b,
			 int flag0, int flag1, int flag2, int flag3)
{
  double *R = Ptr(prlocal, "flagrestrict");
  double *x = Ptr(prlocal, "x");
  double *y = Ptr(prlocal, "y");
  double *z = Ptr(prlocal, "z");
  int i, k, n;

  /* debug */
  if (0) {
    double coord[3] = {2.0, 2.0, 18.0};
    i = find_one_point(prlocal, coord);
    if (i >= 0) R[i] = 9;
  }

  /* initialize */
  if (b->i) free(b->i);
  if (b->buffer) free(b->buffer);
  b->ni = b->nbuffer = 0;

  /* count points that match one of the flags */
  forallpoints(prlocal, i)
    if (boundaryflag(prlocal, i) != SYMBOUND &&
	boundaryflag(prlocal, i) != GHOBOUND
	|| 0)
      if (R[i] == flag0 || R[i] == flag1 || R[i] == flag2 || R[i] == flag3) 
	b->ni++;

  if (0) if (prlocal->l == 2) prvar01(prlocal, "flagrestrict");
  if (0) printf("R %d %d %d %d: %d flags\n", 
		flag0, flag1, flag2, flag3, b->ni);

  /* allocate storage */
  b->i = imalloc(b->ni);
  b->nbuffer = 3*b->ni;
  b->buffer = dmalloc(b->nbuffer);

  /* store indices and coordinates of points */
  k = 0;
  forallpoints(prlocal, i)
    if (boundaryflag(prlocal, i) != SYMBOUND &&
	boundaryflag(prlocal, i) != GHOBOUND
	|| 0)
      if (R[i] == flag0 || R[i] == flag1 || R[i] == flag2 || R[i] == flag3) {
	b->i[k] = i;
	b->buffer[3*k]   = x[i];
	b->buffer[3*k+1] = y[i];
	b->buffer[3*k+2] = z[i];
	k++;
      }

  /* info */
  if (0) {
    char s[1000];
    snprintf(s, 1000, "buffer %d %d %d %d", flag0, flag1, flag2, flag3);
    printcombuf(b, s, 0);
  }
}




/* initialize communication layout */
void bampi_syncparent_init(tL *level)
{
  tL *pr = level->grid->level[level->l-1];
  tL *prlocal = level->prlocal;
  tComBuf *b;

  /* allocate communication buffers for send and recv 
     each processor creates list of buffers for each other processor 
     there are several types of point lists 
  */
  bampi_alloc_combuflist(pr);
  if (pr != prlocal) 
    bampi_alloc_combuflist(prlocal);

  /* temporary buffer for list creation */
  b = bampi_alloc_combuf();


  /* create com buffers for various point lists */

  /* type 0: send local R=3 */
  init_localpointlist(prlocal, b, 3, 3, 3, 3);
  init_buflist(prlocal, prlocal->com->buflist[0], pr, pr->com->buflist[0], b);

  /* type 1: get local R=0,1,2 */
  init_localpointlist(prlocal, b, 0, 1, 2, 2);
  init_buflist(prlocal, prlocal->com->buflist[1], pr, pr->com->buflist[1], b);

  /* type 2: send local R=2,3 */
  init_localpointlist(prlocal, b, 2, 3, 3, 3);
  init_buflist(prlocal, prlocal->com->buflist[2], pr, pr->com->buflist[2], b);

  /* type 3: get local R=0,1,2,3 */
  init_localpointlist(prlocal, b, 0, 1, 2, 3);
  init_buflist(prlocal, prlocal->com->buflist[3], pr, pr->com->buflist[3], b);

  /* type 4: send local R=1,2,3 */
  init_localpointlist(prlocal, b, 1, 2, 3, 3);
  init_buflist(prlocal, prlocal->com->buflist[4], pr, pr->com->buflist[4], b);

  /* debug */
  if (0) {
    init_localpointlist(prlocal, b, 9, 9, 9, 9);
    init_buflist(prlocal, prlocal->com->buflist[3], pr, pr->com->buflist[3],b);
  }

  /* clean up */
  bampi_free_combuf(b);
}




/****************************************************************************/
/* given all the necessary com bufs, perform actual synchronization of data */

/* general purpose routine that is wrapped below for the different types */
void syncparent_type(int type, tL *prsend, tL *prrecv, tVarList *vl)
{
  int pr = 0;
  tComBuf *b, **bsend, **brecv;
  int nv = vl->n;
  int ni, p;
  int me = bampi_rank();
  int size = bampi_size();
  int result, tag;
  MPI_Status *status;
  MPI_Request *request;
  int nrequests = 0;

  /* check */
  if (bampi_size() == 1 && !Getv("amr", "newfmr")) return;
  if (prsend == 0 || prrecv == 0) return;
  bsend = prsend->com->buflist[type];
  brecv = prrecv->com->buflist[type];
  timer_start(0, "syncparent");
  
  /* check number of buffers */
  if (0) {
      printf("cpu %d of %d\n",me,size);
      for (p = 0; p < size; p++) 
          printf("   %d:   ni_send=%d ni_recv=%d\n",p,bsend[p]->ni,brecv[p]->ni);
  }

  /* allocate memory 
     shouldn't matter performance-wise too much, but we do not deallocate bufs
  */
  request = calloc(2*size, sizeof(MPI_Request));
  status  = calloc(2*size, sizeof(MPI_Status));
  for (p = 0; p < size; p++) {
    ni = bsend[p]->ni;
    if (nv * ni > bsend[p]->nbuffer) {
      free(bsend[p]->buffer);
      bsend[p]->buffer = dmalloc(nv * ni);
      bsend[p]->nbuffer = nv * ni;
    }
    ni = brecv[p]->ni;
    if (nv * ni > brecv[p]->nbuffer) {
      free(brecv[p]->buffer);
      brecv[p]->buffer = dmalloc(nv * ni);
      brecv[p]->nbuffer = nv * ni;
    }
  }

  /* post receives first so that send doesn't need extra buffers */ 
  for (p = 0; p < size; p++) {
    b = brecv[p];
    if (b->ni <= 0) continue;

    tag = me*size + p;
    result = MPI_Irecv(b->buffer, b->nbuffer, MPI_DOUBLE, p, tag,
                       MPI_COMM_WORLD, request + nrequests);
    nrequests++;

    if (pr) printf("%d: Irecv from %d, nbuffer=%d, tag=%d\n",
		  me, p, b->nbuffer, tag);
    prmpiresult(result, "Irecv");
  }

  /* fill buffers and post sends */
  for (p = 0; p < size; p++) {
    b = bsend[p];
    if (b->ni <= 0) continue;

    vlfillbuffer(prsend, b, vl);
 
    tag = p*size + me;
    result = MPI_Isend(b->buffer, b->nbuffer, MPI_DOUBLE, p, tag,
                       MPI_COMM_WORLD, request + nrequests);
    nrequests++;

    if (pr) printf("%d: Isend  to  %d, nbuffer=%d, tag=%d\n",
                  me, p, b->nbuffer, tag);
    prmpiresult(result, "Isend");
  }

  /* wait until all requests have been satisfied
     should use MPI_Probe so that copying can start as soon as the individual 
     messages arrive
  */
  result = MPI_Waitall(nrequests, request, status);
  prmpiresult(result, "Waitall");
  

  /* read buffers */
  for (p = 0; p < size; p++)
    vlreadbuffer(prrecv, brecv[p], vl);

  /* done with inter-level communications */
  free(status);
  free(request);
  timer_stop(0, "syncparent");

  /* this fixes symmetries:
       set_boundary_symmetry(prrecv, vl);
     but since it is not needed in all cases we do it case by case in wrappers
  */

  /* do intra-level synchronization */
  vl->level = prrecv;
  bampi_vlsynchronize(vl);
}




/* send R=3 from local to global */
void bampi_syncparent_send3(tL *level, tVarList *vl)
{
  tL *pr = level->grid->level[level->l-1];
  tL *prlocal = level->prlocal;
  
  if (!prlocal) return;
  syncparent_type(0, prlocal, pr, vl);
  set_boundary_symmetry(pr, vl);
}




/* send R=0,1,2 from global to local */
void bampi_syncparent_recv012(tL *level, tVarList *vl)
{
  tL *pr = level->grid->level[level->l-1];
  tL *prlocal = level->prlocal;

    if (!prlocal) return;
    syncparent_type(1, pr, prlocal, vl);
    set_boundary_symmetry(prlocal, vl);
}




/* send R=2,3 from local to global */
void bampi_syncparent_send23(tL *level, tVarList *vl)
{
  tL *pr = level->grid->level[level->l-1];
  tL *prlocal = level->prlocal;
    
  if (!prlocal) return;
  syncparent_type(2, prlocal, pr, vl);
}




/* send R=0,1,2,3 from global to local */
void bampi_syncparent_recv0123(tL *level, tVarList *vl)
{
  tL *pr = level->grid->level[level->l-1];
  tL *prlocal = level->prlocal;

  if (!prlocal) return;
  syncparent_type(3, pr, prlocal, vl);
  set_boundary_symmetry(prlocal, vl);
}

/* send R=0,1,2,3 from local to global */
void bampi_syncparent_send0123(tL *level, tVarList *vl)
{
  tL *pr = level->grid->level[level->l-1];
  tL *prlocal = level->prlocal;

  if (!prlocal) return;
  syncparent_type(3, prlocal, pr, vl);
  set_boundary_symmetry(pr, vl);
}


/* send R=1,2,3 from local to global */
void bampi_syncparent_send123(tL *level, tVarList *vl)
{
  tL *pr = level->grid->level[level->l-1];
  tL *prlocal = level->prlocal;

  if (!prlocal) return;
  syncparent_type(4, prlocal, pr, vl);
  set_boundary_symmetry(pr, vl);
}
