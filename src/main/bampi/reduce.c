/* reduce.c */
/* Bernd Bruegmann 11/02 */

/* Perform reduction operations like computing the norm or the sum 
   of a variable over a level or grid.
   We sometimes want and sometimes do not want to include certain boundaries.
*/

#include "bam.h"
#include "bampi.h"

#define PR 0




/* experimental macros, keep if more useful than confusing */
#define int_pandv(vl,u,i)                 \
  int nvars  = vl->n;                     \
  int *ivar  = vl->index;                 \
  tL *level  = vl->level;                 \
  double **vars = level->v;               \
  double u;                               \
  int i, nv

#define forpandv(vl,u,i)                  \
  forallpoints(level, i) {                \
    if (boundaryflag(level, i) == 0) {	  \
      for (nv = 0; nv < nvars; nv++) {    \
	u = vars[ivar[nv]][i];        

#define endforpandv }}}

#define forpandvref(vl,u,i)                  \
  forallpoints(level, i) {                \
    if ( (boundaryflag(level, i) == 0 ) || (boundaryflag(level, i) < 2 )) {	  \
      for (nv = 0; nv < nvars; nv++) {    \
	u = vars[ivar[nv]][i];        


void bampi_allreduce(tVarList *vl, tVarList *wl, 
		     double **result, MPI_Op op, int typesum)
{
  int_pandv(vl,u,i);
  double *global, *local;
  int nresult = nvars;
  int n;

  /* space to store our result */
  if (op != MPI_MAX && op != MPI_MIN) nresult++;
  *result = calloc(sizeof(double), 2*nresult);
  global = *result;
  local  = *result + nresult;

  /* compute local result */
  if (op == MPI_MAX) {
    for (nv = 0; nv < nvars; nv++) 
      local[nv] = -DBL_MAX;
    if (typesum == 0) {
      forpandv(vl,u,i) {
	if (u > local[nv]) local[nv] = u;
      } endforpandv;
    } else {
      forpandv(vl,u,i) {
	u = fabs(u);
	if (u > local[nv]) local[nv] = u;
      } endforpandv;
    }
  }
  else if (op == MPI_MIN) {
    for (nv = 0; nv < nvars; nv++) 
      local[nv] = DBL_MAX;
    forpandv(vl,u,i) {
      if (u < local[nv]) local[nv] = u;
    } endforpandv;
  }
  else {
    n = 0;
    for (nv = 0; nv < nvars; nv++) 
      local[nv] = 0;
    if (typesum == 0) {
      forpandv(vl,u,i) {
	local[nv] += u;
	n++;
      } endforpandv;
    }
    if (typesum == 1) {
      forpandv(vl,u,i) {
	local[nv] += u * u;
	n++;
      } endforpandv;
    }
    if (typesum == 2) {
      forpandv(vl,u,i) {
	local[nv] += fabs(u);
	n++;
      } endforpandv;
    }
    if (typesum == 3) {
      forpandv(vl,u,i) {
	local[nv] += u * vars[wl->index[nv]][i];
	n++;
      } endforpandv;
    }
    if (typesum == 4) {
      forpandv(vl,u,i) {
        if (vars[wl->index[0]][i]<0.5)
          local[nv] += u * u;
        n++;
      } endforpandv;
    }
    if (typesum == 5) {
      forpandv(vl,u,i) {
        if (vars[wl->index[0]][i]<0.5)
          local[nv] += u;
        n++;
      } endforpandv;
    }
    if (typesum == 6) {
      forpandvref(vl,u,i) {
        if ((vars[wl->index[0]][i]<0.5) && (vars[wl->index[1]][i]<0.5))
            local[nv] += u;
        n++;
      } endforpandv;
    }
    local[nvars] = n;
  }

  /* reduce */
  MPI_Allreduce(local, global, nresult, MPI_DOUBLE, op, MPI_COMM_WORLD);
}




/* maximum */
void bampi_allreduce_max(tVarList *vl, double **result)
{
  bampi_allreduce(vl, 0, result, MPI_MAX, 0);
}



/* minimum */
void bampi_allreduce_min(tVarList *vl, double **result)
{
  bampi_allreduce(vl, 0, result, MPI_MIN, 0);
}



/* sum */
void bampi_allreduce_sum(tVarList *vl, double **result)
{
  bampi_allreduce(vl, 0, result, MPI_SUM, 0);
}

/* sum of all points with mask == 0 */
void bampi_allreduce_sum_mask(tVarList *vl, double **result, tVarList *mask)
{
  bampi_allreduce(vl, mask, result, MPI_SUM, 5);
}

void bampi_allreduce_sum_2masks(tVarList *vl, double **result, tVarList *mask)
{
  bampi_allreduce(vl, mask, result, MPI_SUM, 6);
}


/* dot product */
void bampi_allreduce_dot(tVarList *vl, tVarList *wl, double **result)
{
  bampi_allreduce(vl, wl, result, MPI_SUM, 3);
}



/* dot product for each pair of components, then sum of the results */
double bampi_allreduce_alldot(tVarList *vl, tVarList *wl)
{
  int i, n;
  double *result, sum = 0;

  bampi_allreduce_dot(vl, wl, &result);

  for (i = 0; i < vl->n; i++)
    sum += result[i];

  free(result);
  return sum;
}



/* l1 norm */
void bampi_allreduce_norm1(tVarList *vl, double **result)
{
  int i, n;

  bampi_allreduce(vl, 0, result, MPI_SUM, 2);
  n = (*result)[vl->n];
  for (i = 0; i < vl->n; i++) 
    (*result)[i] = (n) ? (*result)[i]/n : -1;
}



/* l2 norm */
void bampi_allreduce_norm2(tVarList *vl, double **result)
{
  int i, n;

  bampi_allreduce(vl, 0, result, MPI_SUM, 1);
  n = (*result)[vl->n];
  for (i = 0; i < vl->n; i++) 
    (*result)[i] = (n) ? sqrt((*result)[i]/n) : -1;
}


/* l2 norm of values with mask == 0 */
void bampi_allreduce_norm2_mask(tVarList *vl, double **result, tVarList *mask)
{
    int i, n;
    
    bampi_allreduce(vl, mask, result, MPI_SUM, 4);
    n = (*result)[vl->n];
    for (i = 0; i < vl->n; i++) 
        (*result)[i] = (n) ? sqrt((*result)[i]/n) : -1;
}


/* l2 norm of each variable, and then l2 norm of the results */
double bampi_allreduce_allnorm2(tVarList *vl)
{
  int i, n;
  double *result, sum = 0;

  bampi_allreduce_norm2(vl, &result);

  for (i = 0; i < vl->n; i++)
    sum += result[i] * result[i];

  free(result);
  return sqrt(sum/vl->n);
}


/* l2 norm alias */
void bampi_allreduce_norm(tVarList *vl, double **result) 
{
  bampi_allreduce_norm2(vl, result);
}





/* linfinite norm is maximum of the absolute value */
void bampi_allreduce_normInf(tVarList *vl, double **result)
{
  int i, n;

  bampi_allreduce(vl, 0, result, MPI_MAX, 1);
}


/* linfinite norm of each variable, and then linfinite norm of the results */
double bampi_allreduce_allnormInf(tVarList *vl)
{
  int i, n;
  double *result, sum = 0;
  double max = -DBL_MAX;
  
  bampi_allreduce_normInf(vl, &result);

  for (i = 0; i < vl->n; i++)
    if (result[i] > max) max = result[i];

  free(result);
  return max;
}






/***************************************************************************/
/* reduce vectors as opposed to variable lists */

void bampi_allreduce_sum_vector(void *local, void *global, int n)
{
  int result = MPI_Allreduce(local, global, n, 
			     MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  prmpiresult(result, "Allreduce sum_vector");
}

void bampi_allreduce_max_vector(void *local, void *global, int n)
{
  int result = MPI_Allreduce(local, global, n, 
			     MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  prmpiresult(result, "Allreduce max_vector");
}

void bampi_allreduce_min_vector(void *local, void *global, int n)
{
  int result = MPI_Allreduce(local, global, n, 
			     MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  prmpiresult(result, "Allreduce min_vector");
}


void bampi_allreduce_bbox(tL *level, double *bbox, int *ibbox) 
{
  int i;
  double local[6], global[6];

  for (i = 0; i < 3; i++) local[i+0] = - bbox[2*i+0];
  for (i = 0; i < 3; i++) local[i+3] =   bbox[2*i+1];
  
  bampi_allreduce_max_vector(local, global, 6);

  for (i = 0; i < 3; i++)  bbox[2*i+0] = - global[i+0];
  for (i = 0; i < 3; i++)  bbox[2*i+1] =   global[i+3];

  ibbox[1] = (bbox[1]-bbox[0])/level->dx + 0.5;
  ibbox[3] = (bbox[3]-bbox[2])/level->dy + 0.5;
  ibbox[5] = (bbox[5]-bbox[4])/level->dz + 0.5;
}




/***************************************************************************/
/* reduce scalars */

int bampi_allreduce_sum_int(int i)
{
  double local = i;
  double global;

  bampi_allreduce_sum_vector(&local, &global, 1);
  return (int) global;
}




int bampi_allreduce_max_int(int i)
{
  double local = i;
  double global;

  bampi_allreduce_max_vector(&local, &global, 1);
  return (int) global;
}



int bampi_allreduce_min_int(int i)
{
  double local = i;
  double global;

  bampi_allreduce_min_vector(&local, &global, 1);
  return (int) global;
}




/***************************************************************************/
/* reduce min/max with position */

void bampi_allgather_vector(int nlocal, void *local, void *global)
{
  int result = 
    MPI_Allgather(local,  nlocal, MPI_DOUBLE, 
		  global, nlocal, MPI_DOUBLE, MPI_COMM_WORLD);

  prmpiresult(result, "Allgather vector");
}




void bampi_allreduce_maxminpos(tL *level, char *name,
                double *pvarmax, double *pxmax, double *pymax, double *pzmax, 
                double *pvarmin, double *pxmin, double *pymin, double *pzmin)
{
  double *xp  = Ptr(level, "x");
  double *yp  = Ptr(level, "y");
  double *zp  = Ptr(level, "z");
  double *var = Ptr(level, name);
  double varmax = -DBL_MAX;
  double varmin =  DBL_MAX;
  int i, imin = 0, imax = 0;
  int nlocal, nglobal;
  double *local, *global;

  /* find local extremum */
  forinner1(level, i) {
    if (var[i] > varmax) {
      varmax = var[i];
      imax = i; 
    }
    if (var[i] < varmin) {
      varmin = var[i];
      imin = i; 
    }
  }

  /* create buffer with our local data in place */
  nlocal = 8;
  local = dmalloc(nlocal);
  i = 0;
  local[i++] = varmax;
  local[i++] = varmin;
  local[i++] = xp[imax]; 
  local[i++] = yp[imax]; 
  local[i++] = zp[imax];
  local[i++] = xp[imin]; 
  local[i++] = yp[imin]; 
  local[i++] = zp[imin];
  
  if (0) printf("local:  %f %f %f %f\n",
		local[1], local[5], local[6], local[7]);

  /* communicate with all our friends and share the result */
  nglobal = nlocal * bampi_size();
  global = dmalloc(nglobal);
  bampi_allgather_vector(nlocal, local, global);

  if (0) printf("global: %f %f %f %f\n", 
		global[1], global[5], global[6], global[7]);

  /* find the global winner */
  varmax = global[0];
  varmin = global[1];
  imin = imax = 0;
  for (i = 0; i < nglobal; i += nlocal) {
    if (global[i]   > varmax) {imax = i; varmax = global[i];}
    if (global[i+1] < varmin) {imin = i; varmin = global[i+1];}
  }
  
  /* store the result */
  if (pvarmax) {
    *pvarmax = varmax;
    i = imax + 2;
    *pxmax = global[i++];
    *pymax = global[i++];
    *pzmax = global[i++];
  }
  if (pvarmin) {
    *pvarmin = varmin;
    i = imin + 5;
    *pxmin = global[i++];
    *pymin = global[i++];
    *pzmin = global[i++];
  }

  /* done */
  free(local);
  free(global);
  if (0) printf("return: %f %f %f %f\n", *pvarmin, *pxmin, *pymin, *pzmin);
}




void bampi_allreduce_maxpos(tL *level, char *name,
		 double *pvarmax, double *pxmax, double *pymax, double *pzmax)
{
  bampi_allreduce_maxminpos(level, name, 
			    pvarmax, pxmax, pymax, pzmax, 0, 0, 0, 0);
}




void bampi_allreduce_minpos(tL *level, char *name,
                double *pvarmin, double *pxmin, double *pymin, double *pzmin)
{
  bampi_allreduce_maxminpos(level, name, 
			    0, 0, 0, 0, pvarmin, pxmin, pymin, pzmin);
}
