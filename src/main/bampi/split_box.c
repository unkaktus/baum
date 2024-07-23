/* split_box.c */
/* Bernd Bruegmann 5/02 */

#include "bam.h"
#include "bampi.h"

#define PR 0



/* given the global extent of a rectangular box, return local box */

/* The method used here produces planar cuts across the whole level
   along the axes, which generalizes directly to other shapes like spheres.
   The result is a rectangular mesh such that
     size = xsize * ysize * zsize
     mlocal = m/xsize, nlocal = n/ysize, olocal = o/zsize
   where
     size is the number of processors
     xsize is the number of processors in the x direction
     mlocal is the local number of points in the x direction
*/


/* given the number of processors, make list of all possible factorizations 
   just do something simple and remember that n probably < 10000
*/
int factorization[10000][3];
int nfcts = 0;

void pushfactorization(int i, int j, int k)
{
  if (nfcts >= 10000)
    errorexit("Too many factors in processor number.");
  if (PR) 
    printf("%7d:   %4d  * %4d  * %4d  =  %d\n", nfcts, i, j, k, i*j*k);

  factorization[nfcts][0] = i;
  factorization[nfcts][1] = j;
  factorization[nfcts][2] = k;
  nfcts++;
}




/* factorize by trying to divide by all numbers smaller than the sqrt */
void factorize(int n)
{
  int intsqrt_m, intsqrt_n;
  int i, m, b, f1, f2;

  intsqrt_n = sqrt(n + 1e-9);     /* is sqrt exact for integer roots? */
  for (i = 1; i <= intsqrt_n; i++) {
    if (n % i) continue;
    f1 = i;
    m = n/f1;

    for (b = 0; b < 2; b++) {
      intsqrt_m = sqrt(m + 1e-9);
      for (f2 = 1; f2 <= intsqrt_m; f2++) {
	if (m % f2) continue;
	pushfactorization(f1, f2, m/f2);
	if (f2 != m/f2) 
	  pushfactorization(f1, m/f2, f2);
      }
      if (m == n/m) break;
      f1 = m;
      m = n/f1;
    }
  }
}



/* given number of procs and points, find the size and origin in 1d */
void split_box_1d(int procs, int points, 
		  int myproc, int *mypoints, int *myorigin)
{ 
  int n, rest;

  /* integer division gives number of points but ignores the left overs */
  n = points/procs;

  /* these points are left over */
  rest = points%procs;

  /* distribute the rest evenly over the first few processors */ 
  if (myproc < rest) 
    n++;
  
  /* these are my points */
  *mypoints = n;

  /* for the origin I have to add up all points on previous processors */
  n = (points/procs) * myproc;

  /* but again I have to correct for the left overs */
  if (myproc < rest)
    n += myproc;
  else
    n += rest;

  /* this is where my origin is */
  *myorigin = n;
}




/* given number of procs and points, find the size and origin in 1d 
   this version insists on symmetric point distribution which is needed
   for some symmetry conditions
*/
void split_box_1d_symmetric(int procs, int points, 
			    int myproc, int *mypoints, int *myorigin)
{ 
  int i, n, o, rest;
  int *add = imalloc(procs);

  /* integer division gives number of points but ignores the left overs */
  n = points/procs;

  /* these points are left over */
  rest = points%procs;

  /* distribute the rest symmetrically */ 
  for (i = 0; i < procs; i++) 
    add[i] = 0;
  if (rest % 2) { 
    if (procs % 2) {
      add[procs/2] = 1;
      rest--;
    }
    else
      errorexit("split_box_1d_symmetric:\n"
		"cannot distribute an odd number of points symmetrically\n"
		"over an even number of processors, change your parameters!");
  }
  for (i = 0; i < rest/2; i++)
    add[i] = add[procs-i-1] = 1;

  /* these are my points */
  *mypoints = n + add[myproc];

  /* for the origin I have to add up all points on previous processors */
  for (i = o = 0; i < myproc; i++)
    o += n + add[i];

  /* this is where my origin is */
  *myorigin = o;

  /* info */
  if (0) printf("procs %d, myproc %d, npoints %d, mypoints %d, myorigin %d\n", 
		procs, myproc, points, *mypoints, *myorigin);
  free(add);
}




/* split the box in all three directions */
void split_box(
  int xprocs, int xpoints, int xmyproc, int *xmypoints, int *xmyorigin,
  int yprocs, int ypoints, int ymyproc, int *ymypoints, int *ymyorigin,
  int zprocs, int zpoints, int zmyproc, int *zmypoints, int *zmyorigin)
{
  if (Getv("grid", "quadrant") || Getv("grid", "rotant"))
    split_box_1d_symmetric(xprocs, xpoints, xmyproc, xmypoints, xmyorigin);
  else
    split_box_1d(xprocs, xpoints, xmyproc, xmypoints, xmyorigin);

  split_box_1d(yprocs, ypoints, ymyproc, ymypoints, ymyorigin);
  split_box_1d(zprocs, zpoints, zmyproc, zmypoints, zmyorigin);
}




/* the processors are arranged on a Cartesian grid */
int rank(int xsize, int ysize, int zsize, int xrank, int yrank, int zrank)
{
  return xrank + xsize*yrank + xsize*ysize*zrank;
}

void xyzrank(int rank, int xsize, int ysize, int zsize, 
	      int *xrank, int *yrank, int *zrank)
{
  *zrank = rank/(xsize*ysize);
  *yrank = (rank - (*zrank)*xsize*ysize)/xsize;
  *xrank = rank - (*zrank)*xsize*ysize - (*yrank)*xsize;
}




/* read or compute processor distribution and split box */
void bampi_split_box(int m, int n, int o, int *xs, int *ys, int *zs, 
		     int *mlocal, int *nlocal, int *olocal,
		     int *ilocal, int *jlocal, int *klocal)
{
  int rank, size;
  int xrank, yrank, zrank;
  int xsize = Geti("bampi_xsize");
  int ysize = Geti("bampi_ysize");
  int zsize = Geti("bampi_zsize");
  int nghosts = Geti("bampi_nghosts");
  int pr = Getv("bampi_verbose", "yes");
  int requestedsize = xsize * ysize * zsize;
  long int v, vlocal, ng, nglocal;
  int i;
  int ml, nl, ol;
  long int ngmin = INT_MAX;
  int imin;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  
  /* for symmetries, remove symmetry ghosts but add them back later */
  if (Getv("bampi_symmpointsextra", "yes")) {
    if (Getv("grid_half", "x")) m -= nghosts;
    if (Getv("grid_half", "y")) n -= nghosts;
    if (Getv("grid_half", "z")) o -= nghosts;
  }

  /* if the processor layout has not been specified in the par file */
  if (requestedsize <= 0) {
    int rank;
    
    /* find all possible factorizations */
    if (PR) printf("Factorizing processor number %d\n", size);
		nfcts = 0;
    factorize(size);
    
    /* if we have shells test for special case (n x 1 x 1) 
       ...I think this is obsulete but how knows ... at least a speedup */
    if ( ((double)(m)/(double)(size)>n) && ((double)(m)/(double)(size)>o) ) {
      nfcts = 1;
      factorization[0][0] = size;
      factorization[0][1] = 1;
      factorization[0][2] = 1;
    }

    /* for all factorizations */
    for (i = 0; i < nfcts; i++) {
      xsize = factorization[i][0];
      ysize = factorization[i][1];
      zsize = factorization[i][2];
      if (PR) prdivider(0);
      if (PR) printf("%7d:   %4d  * %4d  * %4d\n", i, xsize, ysize, zsize);
      ng = 0;
      v = 0;

      /* for all processors */
      for (rank = 0; rank < size; rank++) {

	xyzrank(rank, xsize, ysize, zsize, &xrank, &yrank, &zrank);

	split_box(xsize, m, xrank, &ml, ilocal,
		  ysize, n, yrank, &nl, jlocal,
		  zsize, o, zrank, &ol, klocal);

	/* inner volume */
	vlocal = ml*nl*ol;

	/* increase box size for ghosts if a face is not at outer boundary */
	if (xrank > 0) ml += nghosts;
	if (yrank > 0) nl += nghosts;
	if (zrank > 0) ol += nghosts;
	if (xrank < xsize-1) ml += nghosts;
	if (yrank < ysize-1) nl += nghosts;
	if (zrank < zsize-1) ol += nghosts;

	/* each face contributes a certain number of ghost zones */
	/* note that the number of ghosts we are looking for is not
	   just the difference between old and new volume because
	   ghosts at the edges and corners will have to be communicated
           several times */
	/* add up number of ghosts */
	nglocal = 0;
	if (xrank > 0) nglocal += nl*ol;
	if (yrank > 0) nglocal += ml*ol;
	if (zrank > 0) nglocal += ml*nl;
	if (xrank < xsize-1) nglocal += nl*ol;
	if (yrank < ysize-1) nglocal += ml*ol;
	if (zrank < zsize-1) nglocal += ml*nl;
	nglocal *= nghosts;

	if (PR) printf("  proc %d (%d,%d,%d) owns %ld points and %ld ghosts\n", 
		       rank, xrank, yrank, zrank, vlocal, nglocal);
	v += vlocal;
	ng += nglocal;
      }
      if (PR) printf("=> %ld points and %ld ghosts total.\n", v, ng);

      /* merit function to determine best split 
	 should consider 
	 - total number of ghosts
	 - maximal volume per proc
         - cache issues since long lists in x direction are good for cache
           but bad for number of ghosts
      */
      /* remember configuration with minimal number of ghosts */
      if (ng < ngmin) {
	ngmin = ng;
	imin = i;
      }
      if (PR) printf("I remember to have %ld min ghosts from factorization %d. \n", ngmin,imin);
      if (PR) printf("Currently, I have %ld  ghosts in factorization %d. \n", ng,i);
    }

    /* best factorization found has index imin */
    i = imin;
    xsize = factorization[i][0];
    ysize = factorization[i][1];
    zsize = factorization[i][2];

    /* setting the size parameters here means that upon next call
       these same parameters will be used (no new factorization allowed)
    */
    if (Getv("bampi_sizes_constant", "yes")) {
      Seti("bampi_xsize", xsize);
      Seti("bampi_ysize", ysize);
      Seti("bampi_zsize", zsize);
    } else
      errorexit("bampi_sizes_constant=no is not fully implemented yet");
  }

  /* else a factorization has been requested by the parameters */ 
  else {
    /* is requested size allowed? */
    if (requestedsize != size) {
      printf(
        "The bampi_xsize ... parameters request %dx%dx%d = %d processors,\n", 
	 xsize, ysize, zsize, requestedsize);
      printf("but the program was started with %d processors.\n",
	     size);
      errorexit("");
    }
  }

  /* ok, we got a processor factorization */

  /* compute the coordinates of my processor in the processor grid */
  xyzrank(rank, xsize, ysize, zsize, &xrank, &yrank, &zrank);

  /* talk about it */
  if (pr || PR) {
    if (size > 1) 
      printf("%d processors are distributed as %dx%dx%d "
	     "with %d at %d,%d,%d.\n",
	     size, xsize, ysize, zsize, rank, xrank, yrank, zrank);
    else 
      printf("%d processor is distributed as %dx%dx%d "
	     "with %d at %d,%d,%d.\n",
	     size, xsize, ysize, zsize, rank, xrank, yrank, zrank);
    printf("There are %dx%dx%d = %d points.\n", m, n, o, m*n*o);
  }

  /* split the box */
  split_box(xsize, m, xrank, &ml, ilocal,
	    ysize, n, yrank, &nl, jlocal,
	    zsize, o, zrank, &ol, klocal);

  /* correct for ghosts */
  if (Getv("bampi_symmpointsextra", "yes")) {
    if (Getv("grid_half", "x")) {
      m += nghosts; if (xrank == 0) ml += nghosts; else *ilocal += nghosts;
    }
    if (Getv("grid_half", "y")) {
      n += nghosts; if (yrank == 0) nl += nghosts; else *jlocal += nghosts;
    }
    if (Getv("grid_half", "z")) {
      o += nghosts; if (zrank == 0) ol += nghosts; else *klocal += nghosts; 
    }
  }
  if (xrank > 0) {ml += nghosts; *ilocal -= nghosts;}
  if (yrank > 0) {nl += nghosts; *jlocal -= nghosts;}
  if (zrank > 0) {ol += nghosts; *klocal -= nghosts;}
  if (xrank < xsize-1) ml += nghosts;
  if (yrank < ysize-1) nl += nghosts;
  if (zrank < zsize-1) ol += nghosts;

  /* store result */
  *xs = xsize;
  *ys = ysize;
  *zs = zsize;
  *mlocal = ml;
  *nlocal = nl;
  *olocal = ol;

  /* diagnostics */
  if (PR) {
    int np, ml, nl, ol, il, jl, kl;

    printf("x-direction:  %d processors and %d points\n", xsize, m);
    for (np = 0; np < xsize; np++) {
      split_box(xsize, m, np, &ml, &il, 1, n, 1, &nl, &jl, 1, o, 1, &ol, &kl);
      printf("Processor %4d:  %4d points starting at %4d\n", np, ml, il); 
    }

    printf("y-direction:  %d processors and %d points\n", ysize, n);
    for (np = 0; np < ysize; np++) {
      split_box(1, m, 1, &ml, &il, ysize, n, np, &nl, &jl, 1, o, 1, &ol, &kl);
      printf("Processor %4d:  %4d points starting at %4d\n", np, nl, jl); 
    }

    printf("z-direction:  %d processors and %d points\n", zsize, o);
    for (np = 0; np < zsize; np++) {
      split_box(1, m, 1, &ml, &il, 1, n, 1, &nl, &jl, zsize, o, np, &ol, &kl);
      printf("Processor %4d:  %4d points starting at %4d\n", np, ol, kl); 
    }
  }
}




/* print communication structure */
void prcom(tCom *c) 
{
  int d, i;

  printf("prcom:  c = %p, ", c);
  printf("myrank = %d, size = %d, ", c->myrank, c->size);
  printf("sizexyz = %d %d %d\n", c->sizexyz[0], c->sizexyz[1], c->sizexyz[2]);

  if (0) 
    for (d = 0; d < 6; d++) {
      printf("%3d  send.ni=%d send.i=%p  recv.ni=%d recv.i=%p\n", c->nbrank[d],
	   c->send[d]->ni, c->send[d]->i,
	   c->recv[d]->ni, c->recv[d]->i);
      if (c->recv[d]->ni && c->recv[d]->i) {
	for (i = 0; i < c->recv[d]->ni; i++)
	  printf("%4d ", c->recv[d]->i[i]);
	printf("\n");
      }
    }
}




/* initialize communication structure */
void bampi_init_com(tL *level, 
		    int m, int n, int o, int xsize, int ysize, int zsize)
{
  if (PR)
    printf("  %d %d %d  %d %d %d\n",m,n,o, xsize,ysize,zsize);
  
  tCom *c = level->com;
  int ng = Geti("bampi_nghosts");
  int xrank, yrank, zrank;
  int bound[6], rmin[6], rmax[6], smin[6], smax[6], ir[6], is[6];
  int d, i, j, k, ijk;
  int index;
  int half[3];
  double ds[3];

  /* save pointer to com structure also in the first box */
  forallboxes(level) {
    box->com = c;
  } endforboxes;

  /* store my rank and size */
  c->size   = bampi_size();
  c->myrank = bampi_rank();
  xyzrank(c->myrank, xsize, ysize, zsize, &xrank, &yrank, &zrank);
  c->sizexyz[0] = xsize;
  c->sizexyz[1] = ysize;
  c->sizexyz[2] = zsize;
  c->myrankxyz[0] = xrank;
  c->myrankxyz[1] = yrank;
  c->myrankxyz[2] = zrank;
  if (PR)
    printf("  %d %d %d \n",xrank,yrank,zrank);

  /* store the rank of my neighbors in the processor grid 
     but if there is no neighbor in this direction, store a -1 */
  c->nbrank[0] = (xrank > 0) ? 
    rank(xsize, ysize, zsize, xrank-1, yrank, zrank) : -1;
  c->nbrank[2] = (yrank > 0) ? 
    rank(xsize, ysize, zsize, xrank, yrank-1, zrank) : -1;
  c->nbrank[4] = (zrank > 0) ? 
    rank(xsize, ysize, zsize, xrank, yrank, zrank-1) : -1;
  c->nbrank[1] = (xrank < xsize-1) ? 
    rank(xsize, ysize, zsize, xrank+1, yrank, zrank) : -1;
  c->nbrank[3] = (yrank < ysize-1) ? 
    rank(xsize, ysize, zsize, xrank, yrank+1, zrank) : -1;
  c->nbrank[5] = (zrank < zsize-1) ? 
    rank(xsize, ysize, zsize, xrank, yrank, zrank+1) : -1;


  /* set local bounding box */
  findbbox(level, c->bbox, c->ibbox);

  /* set local bounding box for the points we own
     ownership is defined by cell boundaries */
  ds[0] = level->dx;
  ds[1] = level->dy;
  ds[2] = level->dz;
  half[0] = Getv("grid_half", "x");
  half[1] = Getv("grid_half", "y");
  half[2] = Getv("grid_half", "z");
  for (i = 0; i < 6; i += 2) {
    j = (c->nbrank[i] < 0) ? 0 : ng;
    c->bboxown[i] = c->bbox[i] + (j - 0.5) * ds[i/2];
    if (c->nbrank[i] < 0 && half[i/2]) c->bboxown[i] = 0;
  }
  for (i = 1; i < 6; i += 2) {
    j = (c->nbrank[i] < 0) ? 0 : ng;
    c->bboxown[i] = c->bbox[i] - (j - 0.5) * ds[i/2];
  }


  /* anything else to do? */
  if (ng <= 0) {
    if (PR) prcom(c);
    return;
  }

  /* store the number of the nodes in the ghostzones in each direction */
  for (d = 0; d < 6; d++) if (c->nbrank[d] >= 0) {
    if (d < 2)      i = n * o * ng;
    else if (d < 4) i = m * o * ng;
    else            i = m * n * ng;
    
    c->send[d]->ni = c->recv[d]->ni = i;
  }

  /* allocate memory for the indices in all directions */
  for (d = 0; d < 6; d++) {
    if (c->send[d]->ni > 0) {
      c->send[d]->i = icalloc(c->send[d]->ni);
      c->recv[d]->i = icalloc(c->recv[d]->ni);
    }
  }

  /* store the indices for the ghosts
     this is real work, but later we can trivially loop over those indices
     when it is time to communicate
  */

  /* there are 12 different regions to communicate, 
     let's precompute the ranges
  */

  /* bound[d]: index of box boundary for each direction */
  for (d = 0; d < 6; d++)
    bound[d] = c->ibbox[d];

  /* for the receives, the ghostzone extends ng points inward from the boundary
     rmin[d]: minimal index for receive
     rmax[d]: maximal index for receive, where rmax - rmin = nghosts - 1
  */
  for (d = 0; d < 6; d++) {
    rmin[d] = bound[d] - (d%2)*(ng-1);
    rmax[d] = rmin[d] + (ng-1);
    if (PR && c->nbrank[d] >= 0) 
      printf("rmin %d, rmax %d\n", rmin[d], rmax[d]);
  }

  /* for the sends, we have to send ng points of those we own 
     at the lower bound (d=0,2,4), move ng points up
     at the upper bound (d=1,3,5), move ng points down
  */
  for (d = 0; d < 6; d++) {
    smin[d] = rmin[d] + (d%2 ? -ng : ng);
    smax[d] = rmax[d] + (d%2 ? -ng : ng);
    if (PR && c->nbrank[d] >= 0) 
      printf("smin %d, smax %d\n", smin[d], smax[d]);
  }

  /* now we could loop over each face separately, 
     but here is a lazy 3d version 
  */

  /* index ijk is the same as in make_toplevelbox, node[ijk].i = ijk */
  ijk = 0;

  /* we will push indices onto the index list as they happen to show up */
  for (d = 0; d < 6; d++) ir[d] = is[d] = 0;

  /* 3d loop over all points plus loop over all directions */
  for (k = 0; k < o; k++) { 
    for (j = 0; j < n; j++) {  
      for (i = 0; i < m; i++, ijk++) {
	for (d = 0; d < 6; d++) {

	  /* figure out which index belongs to which direction */
	  if (d < 2)      index = i;
	  else if (d < 4) index = j;
	  else            index = k;

	  /* if we are in the appropriate ghost zone, save the index */

	  /* note that we don't do anything if we are at the physical 
	     boundary, i.e. when recv[d]->i or send[d]->i are NULL;
	     this is a good thing because our index calculations 
	     are nonsense in this case anyway!
	  */

	  if (rmin[d] <= index && index <= rmax[d] && c->recv[d]->i) 
	    c->recv[d]->i[ir[d]++] = ijk;

	  if (smin[d] <= index && index <= smax[d] && c->send[d]->i) 
	    c->send[d]->i[is[d]++] = ijk;

	}
      }
    }
  }

  /* debug */
  if (PR) prcom(c);
}




/* check whether there is the correct minimal number of owned points per
   process
*/
void bampi_check_split_box(tL *level, tCom *com, int pr)
{
  int flag;
  int nghosts = Geti("bampi_nghosts");
  int i, j;
  int n, n0, n1;
  int size = com->size;
  int n01[9];
  int *alln01 = imalloc(9*size);
  int result;

  for (i = 0; i < 6; i += 2) {
    n = com->ibbox[i+1] + 1;
    n0 = n1 = 0;
    if (com->nbrank[i]   >= 0) n0 = nghosts;
    if (com->nbrank[i+1] >= 0) n1 = nghosts;
    if (com->nbrank[i] < 0)
      if (level->grid->half[i/2] && com->bbox[i] < 0) n0 = nghosts;
    n -= n0 + n1;
    n01[3*i/2+0] = n;
    n01[3*i/2+1] = n0;
    n01[3*i/2+2] = n1;
  }

  result = MPI_Allgather(n01,    9, MPI_INT,
			 alln01, 9, MPI_INT, MPI_COMM_WORLD);
  prmpiresult(result, "Allgather");
 
  flag = 0;
  for (i = 0; i < size; i++)
    for (j = 0; j < 3; j++)
      if (alln01[9*i+3*j] < nghosts) {flag = 1; break;}

  if (pr || flag) {
    for (i = 0; i < size; i++) {
      printf("p%2d  ", i);
      for (j = 0; j < 3; j++) { 
	n  = alln01[9*i+3*j+0];
	n0 = alln01[9*i+3*j+1];
	n1 = alln01[9*i+3*j+2];
	printf(" %3d = %d+%d+%d ", n0+n+n1, n0, n, n1);
      }
      printf("\n");
    }
  }

  if (flag) {
    printf("local number of points too small, has to be >= nghosts = %d\n",
	   nghosts);
    printf("WARNING: this can spoil parallelization, or lead to crashes\n");
    printf("         (inside multigrid it tends to be harmless)\n");
    printf("         enlarge nxyz or adjust number of processors\n");
    // errorexit("enlarge nxyz or adjust number of processors");
  }

  free(alln01);
}




void bampi_check_split(tL *level)
{
  int pr = Getv("bampi_verbose", "yes");

  if (pr) {
    printf("checking box split\n");
    if (level == level->grid->level[level->l])
      printf("level %d\n", level->l);
    else
      printf("ghost parent level %d\n", level->l);
  }

  forallboxes(level) {
    if (pr) printf("box %d\n", nbox);
    bampi_check_split_box(level, box->com, pr);
  } endforboxes;
}



















