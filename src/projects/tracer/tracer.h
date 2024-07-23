/* tracer.h */
#include "mpi.h"

typedef struct tparticle
{ 
  int ID;         /*particle identity*/
  int l;          /*level particle lives*/
  double time;    /*particle time*/
  double pos[3];  /*particle position x, y, z*/
} tParticle;

/* communicating particles */
typedef struct tCOMPARTICLE{

  /* also buffers for MPI communication */
  tComBuf **send;
  tComBuf **recv;

  double **mpiboxes;

} tComParticle;

typedef struct tparticlelist
{
  tParticle **particles; /* points to list with n particles*/
  int n;                        /* number of particles*/

  tComParticle *comp;     /* communicating particles with other processors */
} tParticleList;

/* tracer.c*/
tParticleList *get_tracer_particles(tL *level);
int tracer(tL *level);
void tracer_init(tL *level);
void ParticlePush(tParticleList *pl, int pi, int l, double time, double px, double py, double pz);
void ParticleDrop(tParticleList *p, int pi);


/* tracer_synchronize.c*/
void send_receive_particles(tParticleList *plist);
void nonblockrecv_particles(tComParticle *c, int d, MPI_Request *request, int *nrequests);
void nonblocksend_particles(tComParticle *c, int d, MPI_Request *request, int *nrequests);
void recv_nparticles(tComParticle *c, int d, MPI_Request *request, int *nrequests);
void send_nparticles(tComParticle *c, int d, MPI_Request *request, int *nrequests);
void combine_particles(int nlb, double *localbuffer, double *buffer);
void combine_particles_int(int nlb, int *localbuffer, int *buffer);
int boxes_tracer(tL *level);
void get_local_tracer_number(tL *level, double rho_th, int npart, int nproc, int *nstart, int *nlocal);

/* tracer_evolve.c*/
void evolve_tracer(tL *level);
int find_com(tL *level, double xp, double yp, double zp);
int find_com_nb(tL *level, double xp, double yp, double zp);
void communicate_particles(tL *level);

/* tracer_interpolate.c*/
int tracer_interpolate_xyz_local(tL *level, double x, double y, double z, 
                          int iv, int order, double *result);

/* tracer_output.c*/
int output_tracer(tL *level);
void write_tracer(tL *level, int nv, int *iv);
void init_output_tracer(tL *level);

/* tracer_checkpoint.c*/
void tracer_checkpoint(tL *level);
void write_tracer_checkpoint(tL *level);
void read_tracer_checkpoint(tL *level);
