/* tracer_synchronize.c */
/* A Neuweiler 3/23 */

#include "bam.h"
#include "tracer.h"
#include "mpi.h"

/* debug.c */
void prmpiresult(char *file, int line, int result, char *comment);
#define prmpiresult(r,c) prmpiresult(__FILE__, __LINE__, (r), (c))

/* send and receive particles */
void send_receive_particles(tParticleList *plist)
{
    int nproc = bampi_size();
    int me = bampi_rank();

    tComParticle *cp = plist->comp;

    MPI_Request request[2*nproc];
    MPI_Status status[2*nproc];
    int result;
    int d;

    /* make sure to have empty storage for communication*/
    for (d = 0; d < nproc; d++)
        cp->recv[d]->ni      = 0;

    /* first exchange how many particles are send/received*/
    int nrequests=0;
    for (d = 0; d < nproc; d++)
    {

        /* receive */
        recv_nparticles(cp, d, &request[nrequests], &nrequests);

        /* send */
        send_nparticles(cp, d, &request[nrequests], &nrequests);

    }

    /* wait until everything is send and received */
    result = MPI_Waitall(nrequests, request, status);
    prmpiresult(result, "Waitall");
 
    /* ensure we have enough space to receive */
    for (d = 0; d < nproc; d++)
        cp->recv[d]->buffer  = dcalloc(cp->recv[d]->ni*6);
 
    /* now let's send particles to correct processors*/
    nrequests = 0;
    for (d = 0; d < nproc; d++)
    {

        /* receive */
        nonblockrecv_particles(cp, d, &request[nrequests], &nrequests);

        /* send */
        nonblocksend_particles(cp, d, &request[nrequests], &nrequests);

    }
  
    /* wait until everything is send and received */
    result = MPI_Waitall(nrequests, request, status);
    prmpiresult(result, "Waitall");
   
}

void nonblockrecv_particles(tComParticle *c, int srcrank, MPI_Request *request, int *nrequests)
{
    tComBuf *b;
    int result;
    int tag;

    if (srcrank >= 0)
    {
        (*nrequests)++;
        b = c->recv[srcrank];
        tag = 2;

        result = MPI_Irecv(b->buffer, b->ni*6, MPI_DOUBLE, srcrank, tag,
                           MPI_COMM_WORLD, request);

        prmpiresult(result, "Irecv");
    }
}

void nonblocksend_particles(tComParticle *c, int dstrank, MPI_Request *request, int *nrequests)
{
    tComBuf *b;
    int result;
    int tag;

    if (dstrank >= 0)
    {
        (*nrequests)++;
        b = c->send[dstrank];
        tag = 2; 

        result = MPI_Isend(b->buffer, b->ni*6, MPI_DOUBLE, dstrank, tag,
                           MPI_COMM_WORLD, request);

        prmpiresult(result, "Isend");
    }
}

void recv_nparticles(tComParticle *c, int srcrank, MPI_Request *request, int *nrequests)
{
    tComBuf *b;
    int result;
    int tag;

    if (srcrank >= 0)
    {
        (*nrequests)++;
        b = c->recv[srcrank];
        tag = 1;

        result = MPI_Irecv(&b->ni, 1, MPI_INT, srcrank, tag,
                           MPI_COMM_WORLD, request);

        prmpiresult(result, "Irecv");
    }
}

void send_nparticles(tComParticle *c, int dstrank, MPI_Request *request, int *nrequests)
{
    tComBuf *b;
    int result;
    int tag;
    int ns;

    if (dstrank >= 0)
    {
        (*nrequests)++;
        b = c->send[dstrank];
        tag = 1;

        result = MPI_Isend(&b->ni, 1, MPI_INT, dstrank, tag,
                           MPI_COMM_WORLD, request);

        prmpiresult(result, "Isend");
    }
}


void combine_particles(int nlb, double *localbuffer, double *buffer)
{
  int root = 0;
  int source, tag;
  int nproc = bampi_size();
  int me    = bampi_rank();
  int i, iproc;
  MPI_Status status;
  int result;
  int nmessages;
  int ntb;
  int nrequest = 1;
  double *tempbuffer;
  
  /* if this processor is not the root for output, all it has to do is 
     send the local buffer to root
  */
  if (me != root) {
    tag = me+nproc;

    /* first send number of elements */
    result = MPI_Send(&nlb, 1, MPI_INT, root, me, MPI_COMM_WORLD);
    prmpiresult(result, "Send");

    /* then send data */
    result = MPI_Send(localbuffer, nlb, MPI_DOUBLE, root, tag, MPI_COMM_WORLD);
    prmpiresult(result, "Send");
  }

  /* the output processor has to receive and combine the data */  
  if (me == root) {

    /* copy local data without sending */
    for (i = 0; i < nlb;i++) buffer[i] = localbuffer[i];
    
    /* deal with the messages one by one */
    for (iproc = 1; iproc < nproc; iproc++) {  
        tag = iproc+nproc;

        result = MPI_Recv(&ntb, 1, MPI_INT, iproc, iproc, MPI_COMM_WORLD, &status);
        prmpiresult(result, "Recv");

        tempbuffer = dcalloc(ntb);
        result = MPI_Recv(tempbuffer, ntb, MPI_DOUBLE, iproc, tag, MPI_COMM_WORLD, &status);
        prmpiresult(result, "Recv");

        /* copy received data */
        for (i=0; i < ntb;i++)
            buffer[i+nlb] = tempbuffer[i];

        /* done */
        free(tempbuffer);
        nlb += ntb;
    }
  }

   /* wait until everything is send and received */
   bampi_barrier();
}

void combine_particles_int(int nlb, int *localbuffer, int *buffer)
{
  int root = 0;
  int source, tag;
  int nproc = bampi_size();
  int me    = bampi_rank();
  int i, iproc;
  MPI_Status status;
  int result;
  int nmessages;
  int ntb;
  int nrequest = 1;
  int *tempbuffer;
  
  /* if this processor is not the root for output, all it has to do is 
     send the local buffer to root
  */
  if (me != root) {
    tag = me+nproc;

    /* first send number of elements */
    result = MPI_Send(&nlb, 1, MPI_INT, root, me, MPI_COMM_WORLD);
    prmpiresult(result, "Send");

    /* then send data */
    result = MPI_Send(localbuffer, nlb, MPI_INT,    root, tag, MPI_COMM_WORLD);
    prmpiresult(result, "Send");
  }

  /* the output processor has to receive and combine the data */  
  if (me == root) {

    /* copy local data without sending */
    for (i = 0; i < nlb;i++) buffer[i] = localbuffer[i];
    
    /* deal with the messages one by one */
    for (iproc = 1; iproc < nproc; iproc++) {  
        tag = iproc+nproc;

        result = MPI_Recv(&ntb, 1, MPI_INT, iproc, iproc, MPI_COMM_WORLD, &status);
        prmpiresult(result, "Recv");

        tempbuffer = icalloc(ntb);
        result = MPI_Recv(tempbuffer, ntb, MPI_INT, iproc, tag, MPI_COMM_WORLD, &status);
        prmpiresult(result, "Recv");

        /* copy received data */
        for (i=0; i < ntb;i++)
            buffer[i+nlb] = tempbuffer[i];

        /* done */
        free(tempbuffer);
        nlb += ntb;
    }
  }

   /* wait until everything is send and received */
   bampi_barrier();
}


int boxes_tracer(tL *level){

    /* if not activated do nothing!*/
    if (!Getv("tracer","activated"))
        return 0;

    tParticleList *plist;
    plist = get_tracer_particles(level);
    tG *grid = level->grid;

    int nproc = bampi_size();
    int me = bampi_rank();
    int lmax = grid->lmax;

    tComParticle *cp = plist->comp;
    MPI_Request request[2*(nproc-1)];
    MPI_Status status[2*(nproc-1)];
    int result;
    int d,l,i;

    /* make sure to have empty storage */
    for (d = 0; d < nproc; d++)
        for (l = 0; l<=lmax;l++)
            for (i=0; i<6;i++)
                cp->mpiboxes[d][l*6 + i] = 0;

    /* fill in own boxes */
    for (l = 0; l<=lmax;l++)
        for (i=0; i<6;i++)
            cp->mpiboxes[me][l*6+i] = grid->level[l]->com->bboxown[i];

    /* communicate boxes of other processors*/
    int nrequests=0;
    for (d = 0; d < nproc; d++)
    {
        if(d!=me){
            /* receive */
            result = MPI_Irecv(cp->mpiboxes[d], 6*(lmax+1), MPI_DOUBLE, d, 1,
                            MPI_COMM_WORLD, &request[nrequests]);
            prmpiresult(result, "Irecv");
            nrequests++;

            /* send */
            result = MPI_Isend(cp->mpiboxes[me], 6*(lmax+1), MPI_DOUBLE, d, 1,
                            MPI_COMM_WORLD, &request[nrequests]);
            prmpiresult(result, "Isend");
            nrequests++;
        }
    }

    /* wait until everything is send and received */
    result = MPI_Waitall(nrequests, request, status);
    prmpiresult(result, "Waitall");

    /* print results*/
    // for (d = 0; d < nproc; d++)
    //     for (l = 0; l<=lmax;l++){
    //         printf("MPI sub-box proc %d on level %d: \n", d, l);
    //         printf("    ");
    //         for (i=0; i<6;i++)
    //             printf("%e, ", cp->mpiboxes[d][l*6+i]);
    //         printf("\n");

    //     }

    
    return 0;
}


void get_local_tracer_number(tL *level, double rho_th, int npart, int nproc, int *nstart, int *nlocal){

    if (nproc ==1){
        *nstart = 0;
        *nlocal = npart;
        return;
    }

    int me = bampi_rank();
    if (rho_th <=0){
        *nstart = npart / nproc * me;
        *nlocal = npart / nproc;
        return;
    }

    /* distribute particles, uniform 
    where we have matter with rho > rho_th*/
    double ratio, npoint_global, npoint_loc = 0;
    int i, iproc;

    // how many grid locally with rho > rho_th
    double *rho   = Ptr(level, "grhd_rho");
    forallpoints(level,i)
        if ((boundaryflag(level, i) == 0) && (rho[i]>rho_th))
            npoint_loc+=1.0;

    // how many grid globally with rho > rho_th
    MPI_Allreduce(&npoint_loc, &npoint_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    ratio = npoint_loc/npoint_global;
    *nlocal = (int) ((double) npart)*ratio;

    // also need to know where we start
    int nsend,nrecv;
    for(iproc = 0; iproc < nproc; iproc++){
        if (me<iproc) nsend = (int) ((double) npart)*ratio;
        else          nsend = 0;

        MPI_Reduce(&nsend, &nrecv, 1, MPI_INT, MPI_SUM, iproc, MPI_COMM_WORLD);
    }

    *nstart = nrecv;
}
