/* output.c */
/* Bernd Bruegmann 5/02 */

/* functions for the parallelization of utility/output
   in this case communication is special because
   - each processor knows how to compute the global point list
   - only processor 0 needs the data
*/

#include "bam.h"
#include "bampi.h"




/* helper function: copy point data into global buffer
   format of b:   index in point list, data variable 0, data variable 1, ...
   format of buffer: data variable 0, data variable 1, ...
*/
void copypointbuffer(int nb, double *b, double *buffer, int nvariables)
{
  int i, j, k;

  for (i = 0; i < nb;) {
    j = nvariables*b[i++];
    for (k = 0; k < nvariables; k++)
      buffer[j++] = b[i++];
  }
}




/* combine data from diffent point buffers into single buffer on processor 0
   note that the combined data is in the correct order because we
   also communicate the index of the point in the request list, 
   see copypointbuffer()
*/
void bampi_combinepointbuffers(tL *level, int nlb, double *localbuffer, 
			       double *buffer, int nvariables)
{
  int pr = 0;
  int root = 0;
  int source, tag;
  int size   = level->com->size;
  int myrank = level->com->myrank;
  int i, j, k;
  MPI_Status status;
  int result;
  int nmessages;
  int ntb;
  double *tempbuffer;
  
  /* if this processor is not the root for output, all it has to do is 
     send the local buffer to root
  */
  if (myrank != root) {
    tag = myrank;
    result = MPI_Send(localbuffer, nlb, MPI_DOUBLE, root, tag, MPI_COMM_WORLD);
    prmpiresult(result, "Send");
    if (pr) printf("Processor %d sends to %d, tag %d, count %d\n", 
		   myrank, root, tag, nlb);
  }
  
  /* the output processor has to receive and combine the data */  
  if (myrank == root) {

    /* copy local data without sending */
    copypointbuffer(nlb, localbuffer, buffer, nvariables);
    
    /* deal with the messages one by one */
    /* can and should work on messages in order that they arrive
       so that root does not have to wait unnecessarily */
    for (nmessages = size-1; nmessages > 0; nmessages--) {
      
      /* wait until there is an arbitrary message pending */ 
      MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
      source = status.MPI_SOURCE;
      tag    = status.MPI_TAG;

      /* this may be dangerous if we have missed some messages before */
      if (source != tag) {
	printf("Warning: bampi/ouput.c, Probe returned source %d != tag %d\n",
	       source, tag);
      };

      /* need temporary buffer for receive of size unkown to root */
      MPI_Get_count(&status, MPI_DOUBLE, &ntb);
      tempbuffer = malloc(sizeof(double) * ntb);
      if (pr) printf("Processor 0 receives source %d, tag %d, count %d\n", 
		     source, tag, ntb);

      /* receive */
      result = MPI_Recv(tempbuffer, ntb, MPI_DOUBLE, 
			source, tag, MPI_COMM_WORLD, &status);
      prmpiresult(result, "Recv");

      /* copy received data into output buffer */
      copypointbuffer(ntb, tempbuffer, buffer, nvariables);
    
      /* done */
      free(tempbuffer);
    }
  }
}




#if 0
/* unfinished */

/* collect data at root for a given list of points from the processors
   that owns them
   
   this version is non-blocking
   note that 
   - a processor may not have a point because we only look at the bounding box
     in this case some default empty data is returned
   - a processor may be asked to send data more than once whenever the the 
     of points intersects the box more than once
*/
void bampi_getpoints(tL *level, int npoints, double *coords, 
		     int iv, double *data)
{
  MPI_Request request[10000];  /* fix: allocate size */
  MPI_Status status[10000];
  int owner, root = 0;
  int i, istart;
  int tag = 0;
  int result;
  
  /* unfinished code that does not use Gatherv */
  owner = findowner(level, coords);
  for (i = istart = 0; i < npoints; i++) {
    if (owned(level, owner, coords+3*i)) continue;
    
    if (level->com->myrank == root)
      recvpointsfromowner(level, owner, tag, result, 
			  istart, i-istart, coords, iv, data); 
    else
      sendpointstoroot(level, root, tag, result,
		       istart, i-istart, coords, iv, data);
    tag++;
  }
  
  /* wait for communication
     number of requests equals tag 
  */
  result = MPI_Waitall(tag, request, status);
  prmpiresult(result, "Waitall in getpoints");
  
  /* needed? */
  bampi_barrier();
}




/* check whether a point falls within a local bounding box */
int owned(tL *level, int rank, double *coord) {
  int i;
  double *bb = level->com->bbox[rank];
  
  for (i = 0; i < 3; i++)
    if (dless(coord[i], bb[2*i]) ||
	dless(bb[2*i+1], coord[i]))
       return 0;
  return 1;
}


/* find owner of a point */
int findowner(tL *level, double *coord)
{



}

#endif





/* get interpolated data for a variable list given a point list
   returns data on processor 0
   must be called on all processors

   returns data for each of nv variables grouped by points:
     data[nv*i + j] = variable j at point i

   if data is found by several processors, it is assumed that the 
   data is identical to within round-off and one of them is picked at random
*/
void bampi_interpolate_0(tL *level, tVarList *vl, 
			 int npoints, double *coord, double *data)
{
  int order = Geti("order_RP");
  int nw = vl->n + 1;
  double *localbuffer;
  int i, j;
  int flag;

  /* allocate buffer large enough for entire data
     with one extra entry per point */
  localbuffer = dmalloc(nw * npoints);
  
  /* for all points */
  for (i = j = 0; i < npoints; i++) {

    /* interpolate
       does not touch the buffer if interpolation not possible */
    flag = interpolate_xyz_local_minimal(level,
	     coord[3*i], coord[3*i+1], coord[3*i+2],
	     vl->n, vl->index, localbuffer + nw*j + 1, order, LAGRANGE);
    
    /* if data available, save index in first place and increase counter */
    if (flag) localbuffer[nw*(j++)] = (double) i;
  }

  /* each processor has compiled buffer with the data it can provide 
     combine into complete list on root = 0
  */
  bampi_combinepointbuffers(level, nw*j, localbuffer, data, vl->n);

  /* cleanup */
  free(localbuffer);
}




/* test */
void test_bampi_interpolate_0(tL *level)
{    
  double r = 10; // (level->bbox[1] - level->bbox[0])/3;
  double theta, phi;
  double ntheta  = 4;
  double nphi    = 8;
  double npoints;
  double *coord;
  double *data;
  tVarList *vl;
  int i, j, k;

  /* create variable list, test with coordinates themselves */
  vl = vlalloc(level);
  vlpush(vl, Ind("x"));
  vlpush(vl, Ind("y"));
  vlpush(vl, Ind("z"));

  /* create list of points, test with points on sphere */
  npoints = ntheta * nphi;
  coord = dmalloc(3 * npoints);

  for (j = 0; j < ntheta; j++) {
    for (i = 0; i < nphi; i++) {
      theta = PI/ntheta * j;
      phi   = 2*PI/nphi * i;

      k = i + nphi * j;
      coord[3*k]   = r * sin(theta) * cos(phi);
      coord[3*k+1] = r * sin(theta) * sin(phi);
      coord[3*k+2] = r * cos(theta);
    }
  }
  
  /* on processor 0, create buffer for the data that we want to collect */
  if (bampi_rank() == 0)
    data = dcalloc(vl->n * npoints);
  else
    data = 0;

  /* interpolate to processor 0 
     all processors have to call this routine
  */
  bampi_interpolate_0(level, vl, npoints, coord, data);

  /* check result */
  if (bampi_rank() == 0) {
    printbbox(level, level->bbox, 0);
    printf("r = %7.3f\n", r);
    for (i = 0; i < npoints; i++) {
      printf("%5d %7.3f %7.3f %7.3f  %7.3f %7.3f %7.3f\n",
	     i, coord[3*i], coord[3*i+1], coord[3*i+2], 
	     data[vl->n*i], data[vl->n*i+1], data[vl->n*i+2]);
    }
  }

  /* done */
  free(vl);
  free(coord);
  free(data);
}
