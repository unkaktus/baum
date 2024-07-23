/* grid.c */
/* Bernd Bruegmann, 12/99 */
/* experiment with most primitive AMR: tree of single nodes */ 


#include "bam.h"
#include "amr.h"

#define PR 0




/* initialize grid 
   called in main()
   here we sort out the various parameter options before calling
   for example make_grid_box
*/
tG *make_grid(int pr)
{
  tG *grid;
  tL *level;
  int m = Geti("nx");
  int n = Geti("ny");
  int o = Geti("nz");
  int nghosts = Geti("bampi_nghosts");
  int half[3], symm[3];
  char *s;
  double rad; 
  double dx = Getd("dx");
  double dy = Getd("dy");
  double dz = Getd("dz");
  double dtfac = Getd("dtfac");
  double min, xmin, ymin, zmin, xmax, ymax, zmax;
  double x, y, z;
  int mlocal, nlocal, olocal;
  int ilocal, jlocal, klocal;
  int xprocs, yprocs, zprocs;
  int i;
  int c_order,i_order,ighosts;

  /* print info */
  if (pr) {
    prdivider(0);
    printf("Creating grid\n");
  }

  /* override for dxyz and nxyz parameter 
     do it only once since e.g. multigrid will call this routine several times
  */
  if (Getd("dxyz")) {
    dx = dy = dz = Getd("dxyz");
    Setd("dx", dx);
    Setd("dy", dx);
    Setd("dz", dx);
    Setd("dxyz", 0);
  } else {
    if (Getv("grid", "1d") || Getv("grid", "2d"))
      errorexit("Cartoon requires uniform grid spacing, could be relaxed");
  }

  if (Geti("nxyz")) {
    m = n = o = Geti("nxyz");
    Seti("nx", m);
    Seti("ny", m);
    Seti("nz", m);
    Seti("nxyz", 0);
  }

  if (GetdArrayN("amr_nxyz")) {
    m = n = o = GetdEntry("amr_nxyz", 0);
    Seti("nx", m);
    Seti("ny", m);
    Seti("nz", m);
    Seti("nxyz", 0);
  }

  /* provide a hook for special grid preparation */
  RunFun(PRE_GRID, 0);

  /* as a standard choice allow only grids with 4n+2 points in any direction
     start only for some minimal grid since for multigrid it may be different
  */
  if (Getv("grid_fournplustwo", "yes") && !Getv("grid","seminmax")) {
    if (m > 10) {
      m = (m/4)*4 + 2;
      Seti("nx", m);
    }
    if (n > 10) {
      n = (n/4)*4 + 2;
      Seti("ny", n);
    }
    if (o > 10) {
      o = (o/4)*4 + 2;
      Seti("nz", o);
    }
  }

  /*******************************************************/
  /* top level can be only a box */
  if (!Getv("grid", "box"))
    errorexits("initialize_grid: grid = \"%s\" unknown", Gets("grid"));

  /* these are the coordinates of the global box */
  xmin = -(m-1)*dx/2;  xmax = -xmin;
  ymin = -(n-1)*dy/2;  ymax = -ymin;
  zmin = -(o-1)*dz/2;  zmax = -zmin;
  printf("  Creating a rectangular box\n");


  /*******************************************************/
  /* grids for pseudospectral methods */
  if (Getv("grid", "spectral") ||
      GetsLax("coordtrans") && !Getv("coordtrans", "none")) {
    xmin = Getd("xmin");
    ymin = Getd("ymin");
    zmin = Getd("zmin");
    xmax = Getd("xmax");
    ymax = Getd("ymax");
    zmax = Getd("zmax");
    dx = Getd("dx");
    dy = Getd("dy");
    dz = Getd("dz");

    /* preliminary */
    if (Getd("ps_dtfac")) 
      Setd("dtfac", Getd("ps_dtfac"));
    else {
      double max = m; // max3(m,n,o);
      if (Getv("ps_basis_x", "Ce")) max--;
      Setd("dtfac", dtfac * 4/PI/max);
    }
  }
  
  /*******************************************************/

  if (Getv("grid", "neworigin")) {
    xmin = Getd("neworiginx")-dx*(m-1)/2.0;
    ymin = Getd("neworiginy")-dy*(n-1)/2.0;
    zmin = Getd("neworiginz")-dz*(o-1)/2.0;
    xmax = xmin + (m-1)*dx;
    ymax = ymin + (n-1)*dy;
    zmax = zmin + (o-1)*dz;
  }
  
 /*******************************************************/

  if (Getv("grid", "setminmax")) {
    xmin = Getd("xmin");
    ymin = Getd("ymin");
    zmin = Getd("zmin");
    xmax = Getd("xmax");
    ymax = Getd("ymax");
    zmax = Getd("zmax");   
    dx = (xmax-xmin)/(m-1);
    dy = (ymax-ymin)/(n-1);
    dz = (zmax-zmin)/(o-1);
    Setd("dx", dx); 
    Setd("dy", dy);
    Setd("dz", dz);
  }


   /*******************************************************/
  /* 1dtubex, for using a 3D code to do 1d computations */
  if (Getv("grid", "tubex")) {
    n = o = 3;
    ymin = -dy; zmin = -dz;
    ymax =  dy; zmax =  dz;
    if (pr) printf("Creating a 1d tube along x-axis\n");
  }
  
  if (Getv("grid", "tubez")) {
    m = n = 3;
    ymin = -dy; xmin = -dx;
    ymax =  dy; xmax =  dx;
    if (pr) printf("Creating a 1d tube along z-axis\n");
  }
  
  /* plane, for using a 3D code to do 2d computations */
  if (Getv("grid", "plane")) {
    o = 3;
    zmin = -dz;
    zmax =  dz;
    if (pr) printf("Creating a 2d xy plane\n");
  }


  /*******************************************************/
  /* if there are periodic boundaries extend the box such that the user
     specifies size of physical domain  */
  if (Getv("grid", "periodicx")) {
    m = m + 2*Geti("bampi_nghosts");
    xmin = xmin - Geti("bampi_nghosts")*dx;
    xmax = xmax + Geti("bampi_nghosts")*dx;
    if (pr) printf("periodic grid in x direction selected, " 
		   "setting size to %dx%dx%d\n", m, n, o);
  }
  if (Getv("grid", "periodicy")) {
    n = n + 2*Geti("bampi_nghosts");
    ymin = ymin - Geti("bampi_nghosts")*dy;
    ymax = ymax + Geti("bampi_nghosts")*dy;
    if (pr) printf("periodic grid in y direction selected, " 
		   "setting size to %dx%dx%d\n", m, n, o);
  }
  if (Getv("grid", "periodicz")) {
    o = o + 2*Geti("bampi_nghosts");
    zmin = zmin - Geti("bampi_nghosts")*dz;
    zmax = zmax + Geti("bampi_nghosts")*dz;
    if (pr) printf("periodic grid in z direction selected, "
		   "setting size to %dx%dx%d\n", m, n, o);
  }
  

  /*******************************************************/
  /* record grid extent and planar symmetries
     these flags are initialized to zero
     if a grid is reflection symmetric, only half the space is stored
     while inversion symmetry about the z-axis (rotant and quadrant) implies
     halving in the y-direction, but no y-reflection symmetry 
  */
  for (i = 0; i < 3; i++) half[i] = symm[i] = 0;

  if (Getv("grid", "octant"))   half[0] = half[1] = half[2] = 1;
  if (Getv("grid", "bitant"))                       half[2] = 1;
  if (Getv("grid", "rotant"))             half[1]           = 1;
  if (Getv("grid", "quadrant"))           half[1] = half[2] = 1;
  if (Getv("grid", "qreflect"))           half[1] = half[2] = 1;
  
  if (Getv("grid", "octant"))   symm[0] = symm[1] = symm[2] = 1;
  if (Getv("grid", "bitant"))                       symm[2] = 1;
  if (Getv("grid", "rotant"))                                 1;
  if (Getv("grid", "quadrant"))                     symm[2] = 1;
  if (Getv("grid", "qreflect"))           symm[1] = symm[2] = 1;

  /* this will be stored in grid, but at start up we need this before
     grid is available 
  */
  if (half[0]) Appends("grid_half", "x");
  if (half[1]) Appends("grid_half", "y");
  if (half[2]) Appends("grid_half", "z");

  /* reduce actual extent of box */
  if (half[0]) {
    m = m/2 + nghosts;
    xmin = dx * (0.5 - nghosts);
  }
  if (half[1]) {
    n = n/2 + nghosts;
    ymin = dy * (0.5 - nghosts);
  }
  if (half[2]) {
    o = o/2 + nghosts;
    zmin = dz * (0.5 - nghosts);
  }
  
  /* info */
  s = 0;
  if (Getv("grid", "bitant"))   s = "bitant";
  if (Getv("grid", "rotant"))   s = "rotant";
  if (Getv("grid", "quadrant")) s = "quadrant";
  if (Getv("grid", "qreflect")) s = "qreflect";
  if (Getv("grid", "octant"))   s = "octant";
  if (s) printf("  %s selected, setting size to %dx%dx%d\n", s, m, n, o);
  

  /*******************************************************/
  /* we are done sorting out all those grid parameters! */

  /* call this function to set up top level */
  /* should also pass xmax etc and check xmax = (m-1)*dx + xmin */
  grid = make_grid_box(0, m, n, o, xmin, ymin, zmin, dx, dy, dz, pr);

  /* remember symmetry */
  grid->bitant=grid->rotant=grid->quadrant=grid->octant=grid->qreflect=grid->full = 0;
  if(Getv("grid", "bitant"))        grid->bitant = 1;
  else if(Getv("grid", "rotant"))   grid->rotant = 1;
  else if(Getv("grid", "quadrant")) grid->quadrant = 1;
  else if(Getv("grid", "octant"))   grid->octant = 1;
  else if(Getv("grid", "qreflect")) grid->qreflect = 1;
  else                              grid->full = 1;
  
  
  /* define the shape of the amr, fixed max number of 10 */
  char str[100];
  grid->npunc = Geti("amr_npunctures");
  grid->lmaxpunc = (int*) malloc(10*sizeof(int));
  grid->puncpos  = (double**) malloc(10*sizeof(double*));
  printf("  With %d amr punctures\n", grid->npunc);  

  for (i=0; i<10; i++) { 
    grid->lmaxpunc[i] = Geti("amr_lmax");    
   if(Geti("amr_lmax2") > 0)
    grid->lmaxpunc[0] = Geti("amr_lmax2");
    
    //sprintf(str,"amr_lmax%d",i+1);
    //if (i!=1 && ExistPar(str))
    //  grid->lmaxpunc[i] = Geti(str);
    
    grid->puncpos[i] = (double*) malloc(10*sizeof(double));
    grid->puncpos[i][0] = 0.;
    grid->puncpos[i][1] = 0.;
    grid->puncpos[i][2] = 0.;
  }
  
  
  /* save symmetry information */
  for (i = 0; i < 3; i++) { 
    grid->half[i] = half[i];
    grid->symmetric[i] = symm[i];
  }

  /* sanity check */
  bampi_check_split(grid->level[0]);

  /* return pointer to newly created grid */
  return grid;
}


/* check whether symmetries are allowed with the given amr setup */
void test_grid_symmetry(tG* g)
{
  int i;
  
  if (Getv("grid", "bitant")) {
    for (i=0; i<g->npunc; i++) {
      if (g->puncpos[i][2]!=0.)
        errorexit("  bitant: amr puncture is not in xy plane");
    }
  }
  if (Getv("grid", "rotant")) {
    if (g->npunc>2) 
      printf("  rotant: you use more than two amr puncture?!?\n");
    if (g->puncpos[0][2]!=0. || g->puncpos[1][2]!=0.)
      errorexit("  rotant: amr puncture is not in  xy plane");
    if (g->puncpos[0][0]!=-g->puncpos[1][0] ||
        g->puncpos[0][2]!=-g->puncpos[1][2])
      errorexit("  rotant: amr punctures are not in rotant position");
  }
/* dtim: allow in bam14 also quadrant simulations with proper Gauge-project
  if (Getv("grid", "quadrant")) {
    if (g->npunc>1) 
      printf("  quadrant: you use more than one amr puncture?!?\n");
    for (i=0; i<g->npunc; i++) {
      if (g->puncpos[i][0]!=0. || g->puncpos[i][1]!=0. || g->puncpos[i][2]!=0.)
        errorexit("  quadrant: amr puncture is not in center");
    }
  }*/
  if (Getv("grid", "qreflect")) {
    errorexit("  qreflect: is this supported???"); 
  }
  if (Getv("grid", "octant")) {
    if (g->npunc>1) 
      printf("  octant: you use more than one amr puncture?!?\n");
    for (i=0; i<g->npunc; i++) {
      if (g->puncpos[i][0]!=0. || g->puncpos[i][1]!=0. || g->puncpos[i][2]!=0.)
        errorexit("  rotant: amr puncture is not in center");
    }
  }
}



/* wrapper of make_grid_box_local which discards the grid but keeps the level
   this is slightly backward, should change make_grid_box_local 
*/
tL *make_level_box_local(int l, int m, int n, int o,
			 double xmin, double ymin, double zmin, 
			 double dx, double dy, double dz)
{
  tG *grid = make_grid_box_local(l, m, n, o, xmin, ymin, zmin, dx, dy, dz);
  tL *newlevel = grid->level[0];

  free_grid_only(grid);
  return newlevel;
}




/* make grid that consists of a single box
   calls bampi to distribute this box across processors
   calls make_grid_box_local for local boxes
   the dimensions of the box are passed as arguments as opposed to
   make_grid, which reads the parameter database before calling this function
*/
tG *make_grid_box(int l, int m, int n, int o, double xmin, double ymin, double zmin, 
		  double dx, double dy, double dz, int pr)
{
  tG *grid;
  tL *level;
  int mlocal, nlocal, olocal;
  int ilocal, jlocal, klocal;
  int xprocs, yprocs, zprocs;
  double xmax = xmin + (m-1)*dx;
  double ymax = ymin + (n-1)*dy;
  double zmax = zmin + (o-1)*dz;
  int i;

  /* these are the box parameters */
  if (pr) { 
    printf("  nx = %3d, lx = %8.3f, xmin = %8.3f, xmax = %8.3f, dx = %2.16e\n",
	   m, xmax-xmin, xmin, xmax, dx); 
    printf("  ny = %3d, ly = %8.3f, ymin = %8.3f, ymax = %8.3f, dy = %2.16e\n",
	   n, ymax-ymin, ymin, ymax, dy); 
    printf("  nz = %3d, lz = %8.3f, zmin = %8.3f, zmax = %8.3f, dz = %2.16e\n",
	   o, zmax-zmin, zmin, zmax, dz); 
  }

  /* call bampi to find local size and coordinates */
  mlocal = m; 
  nlocal = n; 
  olocal = o;
  ilocal = jlocal = klocal = 0;
  xprocs = yprocs = zprocs = 1;
  bampi_split_box(m, n, o, 
		  &xprocs, &yprocs, &zprocs, 
		  &mlocal, &nlocal, &olocal,
		  &ilocal, &jlocal, &klocal);
  
  /* each processor creates its own top level grid */
  grid = make_grid_box_local(l, mlocal, nlocal, olocal, 
			  xmin+ilocal*dx, ymin+jlocal*dy, zmin+klocal*dz,
			  dx, dy, dz);
  
  /* now that we have a level, save global box information */
  level = grid->level[0];
  level->bbox[0] = xmin;
  level->bbox[2] = ymin;
  level->bbox[4] = zmin; 
  level->bbox[1] = xmax;
  level->bbox[3] = ymax;
  level->bbox[5] = zmax; 
  level->ibbox[0] = 0;
  level->ibbox[2] = 0;
  level->ibbox[4] = 0; 
  level->ibbox[1] = m-1;
  level->ibbox[3] = n-1;
  level->ibbox[5] = o-1;
  level->l = l;

  for (i = 0; i < 6; i++) {
      level->box[0]->bbox[i] = level->bbox[i];
      level->box[0]->ibbox[i] = level->ibbox[i];
  }


  /* call bampi to fill in communication structure of level */
  bampi_init_com(level, mlocal, nlocal, olocal, 
		 xprocs, yprocs, zprocs);

  /* boundary flags */
  set_boundary_flags(level);
  
  /* if we are doing amr, we need flagregrid */
  if (Geti("amr_lmax") > 0) {
    enablevar(grid->level[0], Ind("flagregrid"));
    enablevar(grid->level[0], Ind("flagrestrict"));
    enablevar(grid->level[0], Ind("flagprolong"));
  }

  /* return pointer to newly created grid */
  return grid;
}




/* wrapper of make_grid_box which discards the grid but keeps the level
   this is slightly backward, should change make_grid_box_local 
*/
tL *make_level_box(int l, int m, int n, int o, double xmin, double ymin, double zmin, 
		   double dx, double dy, double dz, int pr)
{
  tG *grid = make_grid_box(l, m, n, o, xmin, ymin, zmin, dx, dy, dz, pr);
  tL *newlevel = grid->level[0];

  free_grid_only(grid);
  return newlevel;
}




/* make level based on parameter data base
   calls make_grid to make stand alone level without children or parents
   for simplicity, we temporarily change the parameters
   should change this to call the new function make_grid_box
*/
tL *make_level(double dx, double dy, double dz, int m, int n, int o, int pr)
{
  tG *grid;
  tL *newlevel;

  /* save original parameters of the given level */
  int mm = Geti("nx");
  int nn = Geti("ny");
  int oo = Geti("nz");
  double ddx = Getd("dx");
  double ddy = Getd("dy");
  double ddz = Getd("dz");

  /* set new parameters */
  Seti("nx", m);
  Seti("ny", n);
  Seti("nz", o);
  Setd("dx", dx);
  Setd("dy", dy);
  Setd("dz", dz);

  /* make a new distributed level */
  grid = make_grid(pr);
  newlevel = grid->level[0];
  free_grid_only(grid);
  
  /* restore original parameters */
  Seti("nx", mm);
  Seti("ny", nn);
  Seti("nz", oo);
  Setd("dx", ddx);
  Setd("dy", ddy);
  Setd("dz", ddz);

  /* return new level */
  return newlevel;
}




/* make level given parameters with domain decomposition based on coarse 
   level: this ensures that each coarse node has 8 children or none as if
   we had done local refinement
*/
tL *make_level_box_coarsealigned(int l, int m, int n, int o, 
				 double xmin, double ymin, double zmin, 
				 double dx, double dy, double dz, int pr)
{
  tG *grid;
  tL *level;
  int mlocal, nlocal, olocal;
  int ilocal, jlocal, klocal;
  int xprocs, yprocs, zprocs;
  double xmax = xmin + (m-1)*dx;
  double ymax = ymin + (n-1)*dy;
  double zmax = zmin + (o-1)*dz;
  int nghosts;

  /* these are the global parameters of the refined box */
  if (pr) { 
    printf("  nx = %3d, lx = %8.3f, xmin = %8.3f, xmax = %8.3f, dx = %.6f\n",
	   m, xmax-xmin, xmin, xmax, dx); 
    printf("  ny = %3d, ly = %8.3f, ymin = %8.3f, ymax = %8.3f, dy = %.6f\n",
	   n, ymax-ymin, ymin, ymax, dy); 
    printf("  nz = %3d, lz = %8.3f, zmin = %8.3f, zmax = %8.3f, dz = %.6f\n",
	   o, zmax-zmin, zmin, zmax, dz); 
  }

  /* initialize coarse grid parameters for 1 processor */
  mlocal = m/2; 
  nlocal = n/2; 
  olocal = o/2;
  ilocal = jlocal = klocal = 0;
  xprocs = yprocs = zprocs = 1;

  /* call bampi to find local size of coarse box */
  nghosts = Geti("bampi_nghosts");
  Seti("bampi_nghosts", nghosts/2);
  bampi_split_box(m/2, n/2, o/2, 
		  &xprocs, &yprocs, &zprocs, 
		  &mlocal, &nlocal, &olocal,
		  &ilocal, &jlocal, &klocal);
  Seti("bampi_nghosts", nghosts);

  /* set parameters of fine box */
  ilocal *= 2;
  jlocal *= 2;
  klocal *= 2;
  mlocal *= 2;
  nlocal *= 2;
  olocal *= 2;

  /* each processor creates its own level */
  grid = make_grid_box_local(l, mlocal, nlocal, olocal, 
			     xmin+ilocal*dx, ymin+jlocal*dy, zmin+klocal*dz,
			     dx, dy, dz);
  level = grid->level[0];
  free_grid_only(grid);
  
  /* now that we have a level, save global box information */
  level->bbox[0] = xmin;
  level->bbox[2] = ymin;
  level->bbox[4] = zmin; 
  level->bbox[1] = xmax;
  level->bbox[3] = ymax;
  level->bbox[5] = zmax; 
  level->ibbox[0] = 0;
  level->ibbox[2] = 0;
  level->ibbox[4] = 0; 
  level->ibbox[1] = m-1;
  level->ibbox[3] = n-1;
  level->ibbox[5] = o-1; 

  /* call bampi to fill in communication structure of level */
  bampi_init_com(level, mlocal, nlocal, olocal, 
		 xprocs, yprocs, zprocs);
  
  
  /* boundary flags */
  set_boundary_flags(level);
  
  /* if we are doing amr, we need flagregrid */
  if (Geti("amr_lmax") > 0) {
    enablevar(level, Ind("flagregrid"));
    enablevar(level, Ind("flagrestrict"));
    enablevar(level, Ind("flagprolong"));
  }

  /* return pointer to newly created level */
  return level;
}



