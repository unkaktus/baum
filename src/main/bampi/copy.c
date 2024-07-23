/* copy.c */
/* Bernd Bruegmann 4/04 */

/* copy data between processors for a given list of
   physical coordinates or for a pecomputed list of indices
*/ 

#include "bam.h"
#include "bampi.h"



/* get data given a point list
   inefficient in that points are searched for on all processors every time
   must be called on all processors
   returns data plus one additional entry on how often data was found
   sender does not send ghosts, so 
   - either point is not found 
   - or point is found exactly once
*/
void bampi_getdata(tL *level, tVarList *vl, 
		   int npoints, double *coord, double *data)
{
  int me = bampi_rank();
  int size = bampi_size();
  int *allnpoints = imalloc(size);
  int result;
  int maxnpoints = -1;
  int nvl = vl->n;
  int i, j, k, m, n, p;
  double *buf, *coordbuf, *databuf;
  int pr = 0;

  /* for the collective communications later on each processor has to
     know the number of points requested on each other processor
     optimize: use MPI_probe to determine size of incoming message
  */
  result = MPI_Allgather(&npoints, 1, MPI_INT,
			 allnpoints, 1, MPI_INT, MPI_COMM_WORLD);
  prmpiresult(result, "Allgather");
  if (pr) {
    printf("Allgather ");
    for (i = 0; i < size; i++) printf(" %d", allnpoints[i]);
    printf("\n");
  }

  /* allocate buffers of sufficient size */
  for (p = 0; p < size; p++)
    if (allnpoints[p] > maxnpoints) maxnpoints = allnpoints[p];
  coordbuf = dmalloc(3*maxnpoints);
  databuf = dmalloc((nvl+1)*maxnpoints);

  /* each processor gets its turn to ask for and to receive data */
  for (p = 0; p < size; p++) {
    n = allnpoints[p];
    if (pr) printf("p = %d\n", p);

    /* initialize data to zero */
    for (i = 0; i < (nvl+1)*n; i++) 
      databuf[i] = 0;

    /* broadcast list of points */
    buf = (p == me) ? coord : coordbuf;
    bampi_bcast_double(buf, 3*n, p); 

    /* look for these points and store data 
       store one extra entry to indicate whether a point was found or not
    */
    for (i = 0; i < n; i++) {

      j = find_one_point_box(level, buf[3*i], buf[3*i+1], buf[3*i+2]);
      
      if (j < 0) continue;

      if (boundaryflag(level, j) == GHOBOUND) continue;

      for (k = 0; k < nvl; k++)
	databuf[(nvl+1)*i + k] = level->v[vl->index[k]][j];

      databuf[(nvl+1)*i + nvl] = 1;
    }

    /* collect data on this processor by reduction with sum
       the extra data entry counts how often a point was found
    */
    result = MPI_Reduce(databuf, data, (nvl+1)*n,
			MPI_DOUBLE, MPI_SUM, p, MPI_COMM_WORLD);
    prmpiresult(result, "Reduce");
  }

  /* debug
     test the assumption that each point is found no more than once
  */
  if (0) {
    for (i = j = 0; i < npoints; i++)
      if (data[(nvl+1)*i + nvl] > 1) j++; 
    j = bampi_allreduce_sum_int(j);
    printf("global j %d\n", j);
    if (j) 
      errorexit("bampi_getdata: unexpected number of returned values");
  }

  /* cleanup */
  free(coordbuf);
  free(databuf);
  free(allnpoints);
}






#if 0

/* get data given a point list
   must be called on all processors
   returns data grouped by points:
     data[nv*i + j] = variable j at point i

   this version uses MPI_Reduce on global list
   unfinished
*/
void bampi_interpolate_0_reduce(tL *level, tVarList *vl, 
				int npoints, double *coord, double *data)
{
  int order = Geti("order_RP");
  int nv = vl->n;
  int nbuffer = nv * npoints;
  int root = 0;

  /* allocate buffers
     "The MPI standard states that the send and receive buffers in 
     MPI_Reduce, etc. must be distinct" 
  */
  localbuffer = dmalloc(nbuffer);
  reducbuffer = dcalloc(nbuffer);
  
  /* initialize */
  for (i = 0; i < nbuffer; i++)
    localbuffer[i] = DBL_MAX;

  /* for all points */
  for (i = 0; i < npoints; i++) {

    /* interpolate, does not change the buffer if interpolation not possible */
    interpolate_xyz_local_minimal(level,
	     coord[3*i], coord[3*i+1], coord[3*i+2],
	     vl->n, vl->index, localbuffer + nv*i, order);
  }

  /* collect data on root by reduction to minimum
     if no processor provided data, result is DBL_MAX
     if one processor provided data, result is min(DBL_MAX,data)=data
     if two processors provided data, it should be identical within round-off
       so the result is again min(DBL_MAX,data,data)=data
  */
  result = MPI_Reduce(localbuffer, reducbuffer, nbuffer,
		      MPI_DOUBLE, MPI_MIN, root, MPI_COMM_WORLD);

  /* copy data, only root may have storage
     here we could reorder the data, the interpolator returns data
     grouped by points
  */
  if (bampi_rank() == root) 
    for (i = 0; i < nbuffer; i++) {
      data[i] = reducbuffer[i];
      if (data[i] == DBL_MAX) 
	errorexit("bampi_interpolate_0 failed to get data for all points");
    }

  /* cleanup */
  free(localbuffer);
  free(reducbuffer);
}

#endif



