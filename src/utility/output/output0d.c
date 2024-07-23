/* output0d.c */
/* Bernd Bruegmann 12/02 */

#include "bam.h"
#include "output.h"

#define PR 0




/* output routine */
void write_scalar_N(char *name, double *value, int N, char *suffix)
{
  int i;
  char filename[1000];
  snprintf(filename, 1000,"%s/%s%s%s", Gets("outdir_0d"), name, suffix, boxstring);
  
  FILE *fp = fopen(filename, "a");
  if (!fp) errorexits("failed opening %s", filename);
  for (i=0; i<N; i++)
    fprintf(fp, "%22.15e  ",value[i]);
  fprintf(fp, "\n");
  fclose(fp);
}

/* 0d output: norms and extrema */
void write_level_0d(tL *level, int nv, int *iv, int type, char *suffix)
{
  if (PR) printf("write_level_0d  type %d  %s\n",type,suffix);
  
  int pr = 0;
  int i, j;
  tVarList *vl, *mask, *camrmask;
  double *result, v[2];
  char str[100];
  int imask = IndLax("matter_mask");
  int cmask = IndLax("camr_mask_A");

  if (nv <= 0) return;
  if (((type == 5) || (type==6)) && (imask==-1 || !level->v[imask])) return;
  if ((type == 7) &&  ((cmask==-1)|| (imask==-1) || (!level->v[cmask]) || (!level->v[imask]))) return;

  /* preliminary: should use variable lists everywhere in output */
  vl = vlalloc(level);
  if (imask)
    mask = vlalloc(level);
  if(cmask) 
    camrmask = vlalloc(level);

  for (j = 0; j < nv; j++)
    vlpushone(vl, iv[j]);
  if(imask)
    vlpushone(mask, imask);
  if(cmask)
    vlpushone(mask, cmask);

  /* get scalar */
  /* for output we should implement reduce instead of allreduce */
  if (type == 0) bampi_allreduce_max(vl, &result);
  if (type == 1) bampi_allreduce_min(vl, &result);
  if (type == 2) bampi_allreduce_norm(vl, &result);
  if (type == 3) bampi_allreduce_normInf(vl, &result);
  if (type == 4) bampi_allreduce_sum(vl, &result);
  if (type == 5) bampi_allreduce_norm2_mask(vl, &result, mask);
  if (type == 6) bampi_allreduce_sum_mask(vl, &result, mask);
  if (type == 7) bampi_allreduce_sum(vl, &result); //bampi_allreduce_sum_2masks(vl, &result, mask);

   /* correct norms*/
  /* otherwise */
  if(1){
    if ((type == 2) || (type == 5)){
      for (j = 0; j < nv; j++) result[j] *= sqrt(nv);
    }
  }

  if (level->shells) {
     /* integrals have to be changed since the domain is not cubic */
    if ((type==4) || (type==6) || (type==7)) {
      //FIXME: what is the correct answer??? is there a simple way, or
      //       do we have to change the bampi functions???
      for (j = 0; j < nv; j++) result[j] *= 0;
    }
    if ((type==2) || (type==3) || (type==5)) {
      //FIXME: dito
      for (j = 0; j < nv; j++) result[j] *= 0;
    }
  } else {
    /* if it is a integral multiply by dx^3 */
    if ((type==4) || (type==6) || (type==7)) {
      for (j = 0; j < nv; j++) result[j] *= pow(level->dx,3);
    }
  }

  /* processor 0 does the writing */
  if (processor0) {
    sprintf(str,"_%s.l%d",suffix,level->l);
    v[0] = level->time;
    
    /* for each variable */
    for (j = 0; j < nv; j++) {
      v[1] = result[j];
      write_scalar_N(VarName(iv[j]),v,2, str);
    }
  }

  /* wait, could be optimized across different writes */
  bampi_barrier();
  free(result);
  vlfree(vl);
  if(imask) vlfree(mask);
  if(cmask) vlfree(camrmask);

}

/* 0d output: point information */
void write_finestlevel_point(tL *level, int nv, int *iv, int type, char *suffix)
{
  if (level->l!=level->grid->lmax) return;
  if (PR) printf("write_finestlevel_point  type %d\n",type);
  
  int i,j;
  char str[100], name[100], *value;
  double v[5];
  double xmax, ymax, zmax, varmax;
  double xmin, ymin, zmin, varmin;
  
  // min/max
  if (type==0) {
    
    /* for each variable */
    for (j = 0; j < nv; j++) {
      sprintf(name,"%s", VarName(iv[j]));
      
      bampi_allreduce_maxminpos(level, name, 
                                &varmax, &xmax, &ymax, &zmax,
                                &varmin, &xmin, &ymin, &zmin);
      v[2] = xmax;
      v[3] = ymax;
      v[4] = zmax;
      v[0] = varmax;
      v[1] = level->time;
      sprintf(str,"_pmax.l%d",level->l);
      if (processor0)
        write_scalar_N(name,v,5, str);
      
      v[2] = xmin;
      v[3] = ymin;
      v[4] = zmin;
      v[1] = varmin;
      v[0] = level->time;
      sprintf(str,"_pmin.l%d",level->l);
      if (processor0)
        write_scalar_N(name,v,5, str);
    }
    
  }
  
  
  // amr punctures
  if (type==1) {
    
    //value = interpolate_xyz_scalar(level, x,y,z, vi,0);
  }
  
  // points
  if (type==2) {
    
    i=0;
    while (value = NextEntry( Gets("0doutputpoint"))) {
      v[2+i%3] = atof(value);
      i++;
      
      // if there are 3 values we have a point
      if (i%3==0) {
        for (j = 0; j < nv; j++) {
          v[1] = interpolate_xyz_scalar(level, v[2],v[3],v[4], iv[j], output_order,output_scheme);
          v[0] = level->time;
          sprintf(str,"_pt%d.l%d",i%3,level->l);
          if (processor0)
            write_scalar_N(VarName(iv[j]),v,5, str);
        }
      }
      
    }
    
    
    if (i%3!=0)
      errorexit("you need tripeds of coordinates at 0doutputpoint");
    
    //value = interpolate_xyz_scalar(level, x,y,z, vi,0);
  }
  
  
  
}




