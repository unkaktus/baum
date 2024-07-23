/* tracer_evolve.c */
/* A Neuweiler 3/23 */

#include "bam.h"
#include "tracer.h"

#define PR 0

void evolve_tracer(tL *level)
{ 
  tParticleList *plist = get_tracer_particles(level);

  int order= Geti("tracer_interpolation_order");
  int nparticles = plist->n;
  int ip, *pdrop = imalloc(nparticles);

  bampi_openmp_start

  tParticle *p;

  double xp, yp, zp, dt;
  double vx, vy, vz, betax, betay, betaz, alpha;
  int comi,check;

  /* go through all particles of local particle list*/
  /* first move it and if it leaves local domain, send it to neighbors */
  bampi_openmp_loop
  for (ip=0; ip<nparticles;ip++){

    p = plist->particles[ip];
    pdrop[ip] = 0;

    if (p->l != level->l) continue;
    if (p->time >= level->time + level->dt){
      printf("skip Tracer-Particle %d: level time = %e, particle time = %e \n",p->ID, level->time,p->time);
      continue;
    }

    xp = p->pos[0];
    yp = p->pos[1];
    zp = p->pos[2];
    // time step using simple euler
    dt = level->time + level->dt - p->time;

    check = tracer_interpolate_xyz_local(level, xp, yp, zp, Ind("grhd_vx"), order, &vx);
    check = tracer_interpolate_xyz_local(level, xp, yp, zp, Ind("grhd_vy"), order, &vy);
    check = tracer_interpolate_xyz_local(level, xp, yp, zp, Ind("grhd_vz"), order, &vz);

    check = tracer_interpolate_xyz_local(level, xp, yp, zp, Ind("betax"), order, &betax);
    check = tracer_interpolate_xyz_local(level, xp, yp, zp, Ind("betay"), order, &betay);
    check = tracer_interpolate_xyz_local(level, xp, yp, zp, Ind("betaz"), order, &betaz);

    check = tracer_interpolate_xyz_local(level, xp, yp, zp, Ind("alpha"), order, &alpha);

    /* check if interpolation worked, otherwise drop this particle! */
    if (!check){
      printf("Interpolation to Tracer Particle faild! Better drop this particle (ID %d)... \n", p->ID);
      // ParticleDrop(plist,p->ID);
      // ip--;
      pdrop[ip] = 1;
      continue;
    }

    /* use simple euler schme for particles */
    xp += dt*(alpha*vx - betax);
    yp += dt*(alpha*vy - betay);
    zp += dt*(alpha*vz - betaz);

    /* check if it is still on level*/
    if(xyzinsidebbox(level->bbox, xp, yp, zp)){

      /* update particle */
      p->pos[0]=xp;
      p->pos[1]=yp;
      p->pos[2]=zp;
      p->time +=dt;

    } else {
      
      /* move to parent level!*/
      for ( ; p->l >= 0; p->l--){
        if (PR)
          printf("Tracer-Particle %d left level %d! \n", p->ID, p->l);

        /* check if particle lives on this level,
           otherwise continue the loop */
        if(xyzinsidebbox(level->grid->level[p->l]->bbox, xp, yp, zp))
              break;
      }
      
      /* check if Tracer Particle moved of simulation domain !*/
      if ((p->l<0) || (!xyzinsidebbox(level->grid->level[0]->bbox, xp, yp, zp))){
        if (PR) printf("Tracer-Particle %d left simulation domain! \n", p->ID);
        // ParticleDrop(plist,p->ID);
        // ip--;
        pdrop[ip] = 1;
      }

    }
  }
  bampi_openmp_stop

  // drop particles
  tParticle *part_d; 
  for (ip=0; ip<nparticles;ip++){
    if (pdrop[ip] > 0.5) {
        part_d = plist->particles[ip];
        ParticleDrop(plist,part_d->ID);
        ip--;
        nparticles--;
    }
  }

  free(pdrop);

  /* after evoling, might have left MPI sub-box */
  if (bampi_size()>1)
    communicate_particles(level);

}

void communicate_particles(tL *level){

  tParticleList *plist = get_tracer_particles(level);
  tComParticle *comp = plist->comp;
  tParticle *p; 
  tL *level_p;

  int nprocs = bampi_size();
  double xp, yp, zp;
  int nparticles = plist->n;
  int i, ip, comi, d;

  int *np = icalloc(nprocs);
  double *buffer[nprocs];

  for (d = 0; d < nprocs; d++){
    np[d]=0;
    buffer[d]= dcalloc(6*nparticles);
  }

  for (ip=0; ip<nparticles;ip++){
    p = plist->particles[ip];
    level_p = level->grid->level[p->l];
    xp = p->pos[0];
    yp = p->pos[1];
    zp = p->pos[2];
    
    if(!insideownership(level_p, xp, yp, zp)){

      if (PR){
        printf("Tracer-Particle %d left MPI sub-box! \n", p->ID);
        printf("  location of Tracer-Particle: x=%e,y=%e,z=%e \n",xp,yp,yp);
        printf("  size of MPI sub-box: x=[%e,%e],y=[%e,%e],z=[%e,%e] \n",level_p->com->bboxown[0],level_p->com->bboxown[1],
                                                                         level_p->com->bboxown[2],level_p->com->bboxown[3],
                                                                         level_p->com->bboxown[4],level_p->com->bboxown[5]);
      }

      /* check to which processor we have to send it */
      // comi = -1;
      /* if still on same level must be neighbour (easier and more save)*/
      // if (p->l == level->l)
      //   comi = find_com_nb(level_p, p->pos[0], p->pos[1], p->pos[2]);
      // /* otherwise check globally */
      // if (comi<0)
      comi = find_com(level_p, p->pos[0], p->pos[1], p->pos[2]);

      /* sanity check */
      if (comi >= bampi_size() || comi<0 ){
        printf("Couldn't find MPI sub-box for Tracer Particle! Better drop this particle (ID %d)... \n", p->ID);
        ParticleDrop(plist,p->ID);
        ip--;
        nparticles--;
        continue;
      }
      if (comi!=bampi_rank() && comi < bampi_size() && comi>=0 ){
        if (PR)
          printf("Tracer-Particle %d is sent to MPI rank %d! \n", p->ID, comi);

        /* fill buffer */
        buffer[comi][np[comi]*6]   = (double) p->ID;
        buffer[comi][np[comi]*6+1] = (double) p->l;
        buffer[comi][np[comi]*6+2] =          p->time;
        buffer[comi][np[comi]*6+3] =          xp;
        buffer[comi][np[comi]*6+4] =          yp;
        buffer[comi][np[comi]*6+5] =          zp;
        np[comi]++;
      }
    }
  }

  /* make sure to have empty storage
     and fill communication buffers*/
  for (d = 0; d < nprocs; d++){
    comp->send[d]->ni      = np[d];
    comp->send[d]->buffer  = dcalloc(comp->send[d]->ni*6);
    for (i=0;i<6*np[d];i++)
      comp->send[d]->buffer[i]=buffer[d][i];
  }

  /* communicate particles between MPI processes */
  send_receive_particles(plist);

  /* read received particles */
  for (d = 0; d < nprocs; d++)
    for (ip=0;ip<comp->recv[d]->ni;ip++)
        ParticlePush(plist, (int) comp->recv[d]->buffer[ip*6],
                            (int) comp->recv[d]->buffer[ip*6+1],
                                  comp->recv[d]->buffer[ip*6+2],
                                  comp->recv[d]->buffer[ip*6+3],
                                  comp->recv[d]->buffer[ip*6+4],
                                  comp->recv[d]->buffer[ip*6+5]);
  

  /* remove sent particles */
  for (d = 0; d < nprocs; d++)
    for (ip=0;ip<comp->send[d]->ni;ip++)
      ParticleDrop(plist, comp->send[d]->buffer[ip*6]);

  /* everything is sent and received, let's empty communication storage again*/
  for (d = 0; d < nprocs; d++){
    free(comp->send[d]->buffer);
    free(comp->recv[d]->buffer);
    comp->send[d]->ni      = 0;
    comp->recv[d]->ni      = 0;
    free(buffer[d]);
  }
  free(np);

}

int find_com(tL *level, double xp, double yp, double zp){
  int rank;
  int nproc = bampi_size();
  int l=level->l;

  tParticleList *plist = get_tracer_particles(level);
  tComParticle *comp = plist->comp;
  double **mpiboxes = comp->mpiboxes;
  double xmin, xmax, ymin, ymax, zmin, zmax;

  /* go through mpi boxes and find point */
  for (rank=0; rank<nproc; rank++){
    xmin = mpiboxes[rank][l*6];
    xmax = mpiboxes[rank][l*6+1];
    ymin = mpiboxes[rank][l*6+2];
    ymax = mpiboxes[rank][l*6+3];
    zmin = mpiboxes[rank][l*6+4];
    zmax = mpiboxes[rank][l*6+5];

    if ((xmin < xp && xp < xmax) &&
        (ymin < yp && yp < ymax) && 
        (zmin < zp && zp < zmax)  )
        break;
  }

  /* check if result is reasonable */
  if (rank>=nproc)
      return -1;
  
  /* return the rank of processor that covers xp, yp, zp */
  return rank;
}

int find_com_nb(tL *level, double xp, double yp, double zp){
  int rank;
  int xrank, yrank, zrank;

  /* the processors are arranged on a Cartesian grid */
  /* find out properties of the grid */
  int xsize = level->com->sizexyz[0];
  int ysize = level->com->sizexyz[1];
  int zsize = level->com->sizexyz[2];

  /* ownership */
  double xmin = level->com->bboxown[0];
  double xmax = level->com->bboxown[1];
  double ymin = level->com->bboxown[2];
  double ymax = level->com->bboxown[3];
  double zmin = level->com->bboxown[4];
  double zmax = level->com->bboxown[5];

  /* find x, y, z rank in processor grid */
  if (xp < xmin)
    xrank = level->com->myrankxyz[0];
  else if (xp > xmax)
    xrank = level->com->myrankxyz[0];
  else
    xrank = level->com->myrankxyz[0];

  if (yp < ymin)
    yrank = level->com->myrankxyz[1];
  else if (yp > ymax)
    yrank = level->com->myrankxyz[1];
  else
    yrank = level->com->myrankxyz[1];

  if (zp < zmin)
    zrank = level->com->myrankxyz[2];
  else if (zp > zmax)
    zrank = level->com->myrankxyz[2];
  else
    zrank = level->com->myrankxyz[2];

  /* check if result is reasonable */
  if ((xrank < 0 || level->com->sizexyz[0] <= xrank ) || 
      (yrank < 0 || level->com->sizexyz[1] <= yrank ) || 
      (zrank < 0 || level->com->sizexyz[2] <= zrank ) )
      return -1;
  
  /* return the rank of processor that covers xp, yp, zp */
  return xrank + xsize*yrank + xsize*ysize*zrank;
}

