/* tracer_checkpoint.c */
/* A Neuweiler 7/23 */

#include "bam.h"
#include "tracer.h"

#define PR 0

void tracer_checkpoint(tL *level)
{

  /* need same routines as in checkpoint to find out when we have to write checkpoint files for tracers */
  int di, skip_checkpoint_i = 1, skip_checkpoint_t = 1, skip_checkpoint_T = 1;
  double uptime, uptime_since_checkpoint, dt, dT, dt_hours, dt_hours_quit;
  static double last_checkpoint_uptime = 0.;

  uptime = bampi_uptime()/3600.;
  bampi_barrier();

  /* check whether we are on */
  if (Getv("checkpoint", "no")) return;

  /* check whether we are restarting */
  if (Getv("checkpoint", "restart")){
    if (level->iteration == 1)
      read_tracer_checkpoint(level);
    return;
  }

  if (!Getv("tracer","activated")) 
        return;

  /* determine how often we want to checkpoint */
  di            = Geti("checkpoint_di");
  dt            = Getd("checkpoint_dt");
  dt_hours      = Getd("checkpoint_dt_hours");
  dt_hours_quit = Getd("checkpoint_dt_hours_quit");
  
  if (dt <= 0. && dt_hours <= 0. && dt_hours_quit <= 0.) 
    skip_checkpoint_i = (level->iteration % di);

  if (dt > 0) {
    di = (int)(dt/level->dt + dequaleps);
    skip_checkpoint_t = level->iteration % di;
  }

  uptime_since_checkpoint = uptime - last_checkpoint_uptime;
  if ((dt_hours      > 0.  &&  dt_hours      <= uptime_since_checkpoint)   ||
      (dt_hours_quit > 0.  &&  dt_hours_quit <= uptime)) {
    skip_checkpoint_T = 0;
    last_checkpoint_uptime = uptime;
  }

  /* check whether it is time to write a checkpoint */
  if (skip_checkpoint_i && skip_checkpoint_t && skip_checkpoint_T)
    return;

  /* write checkpoint for tracers! */
  write_tracer_checkpoint(level);
}

void write_tracer_checkpoint(tL *level)
{
  char dir[1000],filename[1000];
  char formatstring[1000];

  /* for particles*/
  int iproc, ip,i;
  tParticle *p;
  tParticleList *plist= get_tracer_particles(level);

  /* get file name*/
  sprintf(dir,"%s/tracer",Gets("outdir"));
  snprintf(formatstring, 100, "%%s/checkpoint_tracer.%%0%dd",
	   (int) log10(bampi_size())+1);
  snprintf(filename, 1000, formatstring, dir, bampi_rank());
  
  FILE *fp;

  /* for all processes */
  for(iproc = 0; iproc < bampi_size(); iproc++)
  {    
    /* if it is my turn */
    if(iproc == bampi_rank())
    {
      /* open file for writing */
      fp = fopen(filename, "wb");
      if (!fp) errorexits("failed opening %s", filename);

      /* now we are ready to go! */
      /* write one line text header, check with "head -1 checkpoint.0" */
      fprintf(fp, "$BEGIN_checkpoint_tracer_header:  ");
      fprintf(fp, "$time = %.3f\n", level->time);

      /* save info for each particle */
      for (ip=0; ip<plist->n;ip++)
      {
        p = plist->particles[ip];
        /* write */
	    fprintf(fp, "$ParticleID = %d : level = %d , time = %g , x = %g , y = %g , z = %g \n",
                    p->ID, p->l, p->time, p->pos[0], p->pos[1], p->pos[2]);
      }

      /* done */
      fclose(fp);
    }

    /* everyone please wait here */
    bampi_barrier();
  }
  
}


void read_tracer_checkpoint(tL *level)
{
  char dir[1000],filename[1000];
  char formatstring[1000];
  FILE *fp;

  /* for particles*/
  int iproc,i;
  tParticleList *plist= get_tracer_particles(level);
  tParticle *part;
  int id, l;
  double time, xp, yp, zp;
  char particle[10000], value[10000];

  /* get file name*/
  sprintf(dir,"%s_previous/tracer",Gets("outdir"));
  snprintf(formatstring, 100, "%%s/checkpoint_tracer.%%0%dd",
	   (int) log10(bampi_size())+1);
  snprintf(filename, 1000, formatstring, dir, bampi_rank());

  /* check if there is a complete set of checkpoint files */
  fp = fopen(filename, "rb");
  if (!fp) return;

  if (1) printf("reading Tracer-Particles from checkpoint file:\n");

  /* read one line text header, check with "head -1 checkpoint.0" */
  if(fgotonext(fp, "$BEGIN_checkpoint_tracer_header:") == EOF)
    errorexit("$BEGIN_checkpoint_tracer_header: is missing!");
  fgetparameter(fp, "$time", value);
  printf("   found $time = %s\n", value);

  /* for all processes: read the pars and iteration numbers */
  for (iproc = 0; iproc < bampi_size(); iproc++)
  {
    /* read if it is my turn */
    if (iproc == bampi_rank()){
      while(fgetparameter(fp, "$ParticleID", particle) != EOF) {

        id = atoi(particle);

        fgetparameter(fp, "level", value);
        l = atoi(value);

        fgetparameter(fp, "time", value);
        time  = atof(value);

        fgetparameter(fp, "x", value);
        xp  = atof(value);

        fgetparameter(fp, "y", value);
        yp  = atof(value);

        fgetparameter(fp, "z", value);
        zp  = atof(value);

        // safety, check if particle already exists
        for (i = 0; i < plist->n; i++){
          part = plist->particles[i];
            if (part->ID == id)
              break;
        }
        
        if (part->ID != id)
          ParticlePush(plist,id,l,time,xp,yp,zp);
      }
    }
  
    /* everyone please wait here */
    bampi_barrier();
  }

  /* clean up */
  fclose(fp);
}

