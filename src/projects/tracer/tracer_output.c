/* tracer_output.c */
/* A Neuweiler 7/23 */

#include "bam.h"
#include "tracer.h"

#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "assert.h"

#include "hdf5.h"

int output_tracer(tL *level)
{
    /* only for level 0*/
    if (level->l!=0)
      return 0;
    
    int timer = timer_start(0, "tracer_output");

    int nindex, *index;
    int iteration;
    static int firstcall = 1;
    static double dt;
    static double inittime;
    static char *ou;

    if (Getv("tracer","activated")){
      if (firstcall) {
          firstcall = 0;
          dt = Getd("tracer_outtime");
          ou = Gets("tracer_output");
          inittime = Getd("tracer_init_time");
      }

      if (timefor_tracer_output(level, dt, inittime, &iteration)){
          makeoutputlist(level, ou, 0, &nindex, &index);
          write_tracer(level, nindex, index);
          free(index);
      }
    }

    timer_stop(0, "tracer_output");

    /* we also want to check if we have to write a checkpoint */
    timer = timer_start(0, "tracer_checkpoint");
    tracer_checkpoint(level);
    timer_stop(0, "tracer_checkpoint");

    return 0;
}

int timefor_tracer_output(tL *level, double dt, double inittime, int *iteration) 
{

  /* time for output based on time interval, assumes t >= 0 */
  double dt_min = level->grid->level[ level->grid->lmax ]->dt;
  int i;
  if (dt > 0) {
    i    = (level->time - inittime + 0.1*dt_min)/dt;
    if (dequal(level->time - inittime -i*dt, 0)){
      *iteration = i;
      return 1;
    }
  }
    
  /* not time for output */
  return 0;
}

void write_tracer(tL *level, int nv, int *iv)
{
  char dir[1000],filename[1000],time_id[1000];
  int order = Geti("tracer_interpolation_order");
  int nparticles = Geti("tracer_number");
  tG *grid = level->grid;
  FILE *fp, *fp2;

  /* get file name*/
  sprintf(dir,"%s/tracer",Gets("outdir"));
  snprintf(filename, 1000, "%s/tracer.h5", dir);

  /* for particles*/
  int ip, i, check;
  tParticle *p;
  tParticleList *plist= get_tracer_particles(level);
  tL *p_level;

  /* get output and store in buffer */
  int     nbuf,     nbuffer;
  double *buf,      *buffer, *var_buffer;
  double *xbuf,*ybuf,*zbuf;
  double *xbuffer,*ybuffer,*zbuffer;
  double *timebuf,  *timebuffer;
  int *pibuf, *pibuffer;
  double xp, yp, zp;

  nbuf = plist->n;
  buf   = dcalloc(nv*nbuf);
  pibuf = icalloc(nbuf);
  timebuf   = dcalloc(nbuf);
  xbuf = dcalloc(nbuf);
  ybuf = dcalloc(nbuf);
  zbuf = dcalloc(nbuf);
  
  for (ip=0; ip<plist->n;ip++){
        p = plist->particles[ip];
        xp=p->pos[0];
        yp=p->pos[1];
        zp=p->pos[2];
        p_level = grid->level[p->l];

        pibuf[ip] = p->ID;
        timebuf[ip]   = p->time;
        xbuf[ip] = xp;
        ybuf[ip] = yp;
        zbuf[ip] = zp;
        for (i=0; i<nv;i++)
            check = tracer_interpolate_xyz_local(p_level, xp, yp, zp, iv[i], order, &buf[i+ip*nv]);
  }

  /* get number of particles */
  MPI_Allreduce(&nbuf, &nbuffer, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  /* allocate buffer for data on processor 0 */
  if (processor0) {
    buffer  = dcalloc(nv*nbuffer);
    pibuffer  = icalloc(nbuffer);
    timebuffer   = dcalloc(nbuffer);
    xbuffer = dcalloc(nbuffer);
    ybuffer = dcalloc(nbuffer);
    zbuffer = dcalloc(nbuffer);
    var_buffer = dcalloc(nbuffer);
  }

  /* combine data from different processors */
  combine_particles(nv*nbuf, buf, buffer);
  combine_particles(nbuf, timebuf,   timebuffer);
  combine_particles(nbuf, xbuf, xbuffer);
  combine_particles(nbuf, ybuf, ybuffer);
  combine_particles(nbuf, zbuf, zbuffer);
  combine_particles_int(nbuf, pibuf, pibuffer);

  /* processor 0 does the writing */
  if(processor0){ 
    /* directory */
    fp = fopen(dir, "r");
    if (!fp)
      mkdir(dir, 0777);
    else
      fclose(fp);

    /* h5 stuff*/
    hid_t  file_id, t_id, dataset_id, dataspace_id;
    herr_t status;
    hsize_t dim[1] = {nbuffer};
    hsize_t nobj;

    
    /* Open file */
    fp2 = fopen(filename, "r");
    if (!fp2)
      file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    else
      file_id = H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT);

    /* Initilize / Write output */
    H5Gget_num_objs(file_id, &nobj);
    snprintf(time_id, 1000, "%06llu", nobj+1 );
    t_id = H5Gcreate(file_id, time_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    dataspace_id = H5Screate_simple(1, dim, NULL);
    dataset_id = H5Dcreate(t_id, "Particle-ID", H5T_NATIVE_INT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, pibuffer);
    status = H5Dclose(dataset_id);
    status = H5Sclose(dataspace_id);

    dataspace_id = H5Screate_simple(1, dim, NULL);
    dataset_id = H5Dcreate(t_id, "time", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, timebuffer);
    status = H5Dclose(dataset_id);
    status = H5Sclose(dataspace_id);

    dataspace_id = H5Screate_simple(1, dim, NULL);
    dataset_id = H5Dcreate(t_id, "x", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, xbuffer);
    status = H5Dclose(dataset_id);
    status = H5Sclose(dataspace_id);

    dataspace_id = H5Screate_simple(1, dim, NULL);
    dataset_id = H5Dcreate(t_id, "y", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, ybuffer);
    status = H5Dclose(dataset_id);
    status = H5Sclose(dataspace_id);

    dataspace_id = H5Screate_simple(1, dim, NULL);
    dataset_id = H5Dcreate(t_id, "z", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, zbuffer);
    status = H5Dclose(dataset_id);
    status = H5Sclose(dataspace_id);

    /* Initialize output variables */
    for (i=0; i<nv;i++){
      for (ip=0; ip<nbuffer; ip++)
        var_buffer[ip]= buffer[i+ip*nv];

      dataspace_id = H5Screate_simple(1, dim, NULL);
      dataset_id = H5Dcreate(t_id, VarName(iv[i]), H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, var_buffer);
      status = H5Dclose(dataset_id);
      status = H5Sclose(dataspace_id);
    }

    // Close the group
    status = H5Gclose(t_id);

    // Close the file
    status = H5Fclose(file_id);

    free(pibuffer);
    free(var_buffer);free(buffer);
    free(xbuffer);free(ybuffer);free(zbuffer);
    free(timebuffer); 
    
  }

  free(pibuf);
  free(buf);
  free(timebuf);
  free(xbuf);free(ybuf);free(zbuf);

}
