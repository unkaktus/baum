/* tracer.c */
/* A Neuweiler 12/22 */

#include "bam.h"
#include "tracer.h"
#include <stdlib.h>  

#define PR 0

/* store paricles locally */
tParticleList *p_t = 0;
tParticleList *get_tracer_particles(tL *level)
{
  int i;
  int lmax=level->grid->lmax;
  if (!p_t){
    p_t = calloc(1, sizeof(tParticleList));

    tComParticle *comp = calloc(1, sizeof(tComParticle));
    comp->send = bampi_alloc_combufs(bampi_size());
    comp->recv = bampi_alloc_combufs(bampi_size());
    comp->mpiboxes = (double**) calloc(bampi_size(), sizeof(double*));
    for (i=0; i<bampi_size(); i++)
      comp->mpiboxes[i] = (double*) calloc((lmax+1)*6, sizeof(double));

    p_t->comp   = comp;
  }

  return p_t;
}

int tracer(tL *level){

  int timer = timer_start(0, "tracer");
  static int firstcall = 1;

  /* first check if anything needs to be done!*/
  if (Getv("tracer","no")){
    timer_stop(0, "tracer");
    return 0;
  }

  else if (Getv("tracer","activated")){
    if (Getv("checkpoint", "yes") && firstcall){
      boxes_tracer(level);
      firstcall = 0;
    }

    evolve_tracer(level);
  }
  
  /* check if Tracer Particles need to be initialized */
  else if (Getv("tracer","yes")) {

    /* don't initialize in the restart phase! */
    if (Getv("checkpoint", "restart")){
      timer_stop(0, "tracer");
      return 0;
    }
    
    /* are we on correct level ?*/
    int tracer_level = Geti("tracer_init_level");
    if (level->l != tracer_level) {
      timer_stop(0, "tracer");
      return 0;
    }

    /* check for time */
    if (Getd("tracer_init_time") > level->time + 1e-15){
      timer_stop(0, "tracer");
      return 0;
    }

    /* check if distance between compact objects is close enough*/
    if(Getd("tracer_init_distance")>0){
      double dist;
      double  p[2][3];
      int np, i;

      /* define puncture positions */
      for (np = 0; np < 2; np++){
        for (i = 0; i < 3; i++) {
          p[np][i] = level->grid->puncpos[np][i];
        }
      }

      dist = sqrt( (p[1][0]-p[0][0])*(p[1][0]-p[0][0]) +
                   (p[1][1]-p[0][1])*(p[1][1]-p[0][1]) +
                   (p[1][2]-p[0][2])*(p[1][2]-p[0][2]) );
      
      if (dist>Getd("tracer_init_distance")){
        timer_stop(0, "tracer");
        return 0;
      }

      Setd("tracer_init_time",level->time);
    }

    /* if we came so far we should initialize tracer!*/
    printf("ACTIVATE TRACERS! \n");
    Sets("tracer","activated");
    tracer_init(level);
    firstcall = 0;

    /* also have to evolve them already !*/
    evolve_tracer(level);

  } else errorexit("Value of tracer-Parameter is not valid!");
  
  timer_stop(0, "tracer");
  return 0;
}


void tracer_init(tL *level){
  tParticleList *plist;
  plist = get_tracer_particles(level);
  int nproc = bampi_size();

  int id, n, nstart, npart = Geti("tracer_number");  
  int order = Geti("tracer_interpolation_order");  
  double rho_th = Getd("tracer_init_density_threshold");
  double xnew_p,ynew_p,znew_p,rho_p;
  int check;
  tParticle *p;

  if (Getv("tracer_init_distribution","uniform")){
    
    get_local_tracer_number(level, rho_th, npart, nproc, &nstart, &n);

    for (id=nstart; id<n+nstart; id++){ 

      xnew_p = (double) rand()/RAND_MAX*(level->com->bboxown[1]-level->com->bboxown[0])     +level->com->bboxown[0];
      ynew_p = (double) rand()/RAND_MAX*(level->com->bboxown[3]-level->com->bboxown[2])     +level->com->bboxown[2];
      znew_p = (double) rand()/RAND_MAX*(level->com->bboxown[5]-level->com->bboxown[4])     +level->com->bboxown[4];
      
      check = tracer_interpolate_xyz_local(level,xnew_p,ynew_p,znew_p,Ind("grhd_rho"),order,&rho_p);

      if(rho_p>rho_th && check==1)
        ParticlePush(plist,id,level->l, level->time,xnew_p,ynew_p,znew_p);
      else {
        /* look again ! */
        id--;
        continue;
      }

      if(PR){
        p = plist->particles[plist->n-1];
        printf("set Tracer Particle with ID = %d at x=%e, y=%e, z=%e and time=%e \n",id, p->pos[0],p->pos[1],p->pos[2],p->time);
      }
    }
  } else errorexit("Distribution for Tracer-Particles unknown! ");

  /* initialise output */
  int nindex, *index;
  char *ou = Gets("tracer_output");
  makeoutputlist(level, ou, 0, &nindex, &index);
  write_tracer(level, nindex, index);
  printf("  initialised output for tracer!\n");
  free(index);

  /* set mpi boxes for communication */
  boxes_tracer(level);

}


/* add particle to particle list*/
void ParticlePush(tParticleList *pl, int pi, int l, double time, double px, double py, double pz)
{ 
  tParticle *particle = calloc(1, sizeof(tParticle));

  if (PR) printf("Pushing Particle %d \n",pi);
  particle->ID = pi;
  particle->pos[0] = px;
  particle->pos[1] = py;
  particle->pos[2] = pz;
  particle->l      = l;
  particle->time   = time;

  pl->n++;
  pl->particles = realloc(pl->particles, sizeof(tParticle) * pl->n); 
  pl->particles[pl->n-1] = particle;
}


/* drop particle from a particle list*/
void ParticleDrop(tParticleList *pl, int pi)
{ 
  int i;
  tParticle *p;
  if (PR) printf("Dropping Particle %d \n",pi);
  for (i = 0; i < pl->n; i++){
    p = pl->particles[i];
    if (p->ID == pi) {
      free(p);

      for (; i < pl->n; i++) 
          pl->particles[i] = pl->particles[i+1];

      pl->n--;
      pl->particles = realloc(pl->particles, sizeof(tParticle) * pl->n); 
      break;
    }
  }
}
