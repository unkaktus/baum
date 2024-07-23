/* bam_main.c */
/* Bernd Bruegmann 12/99 */

#include "bam.h"
#include "main.h"



void bam_main() 
{
  printf("Adding main\n");

  /* functions */
  AddFun(PRE_INITIALDATA,  system_monitor, "write system status");
  AddFun(POST_OUTPUT,      system_monitor, "write system status");

  /* variables */

  /* parameters */
  AddPar("physics", "", "what problem to solve");

  AddPar("system_monitor_every", "1",
	 "write system status every n top level iterations [0 for off]");
  AddPar("system_monitor_verbose", "no",
	 "print entire system files to figure out format");

  AddPar("grid", "box", "how to specify grid layout [box, shells]");

  AddPar("nx"  , "5", "points in the x direction for rectangular grids");
  AddPar("ny"  , "5", "points in the y direction for rectangular grids");
  AddPar("nz"  , "5", "points in the z direction for rectangular grids");
  AddPar("nxyz", "0", "points in all three directions for cubical grid");

  AddPar("dx"  , "0.1", "grid spacing in the x direction");
  AddPar("dy"  , "0.1", "grid spacing in the y direction");
  AddPar("dz"  , "0.1", "grid spacing in the z direction");
  AddPar("dxyz", "0.0", "uniform gridspacing");

  AddPar("xmin"  , "0", "range of x");
  AddPar("xmax"  , "0", "range of x");
  AddPar("ymin"  , "0", "range of y");
  AddPar("ymax"  , "0", "range of y");
  AddPar("zmin"  , "0", "range of z");
  AddPar("zmax"  , "0", "range of z");

  AddPar("dtfac", "0.25", "dt = dtfac * ds");
  AddPar("iterations", "0", "number of toplevel iterations");
  AddPar("finaltime", "0", "iterate until toplevel reaches this time");
  AddPar("iterate_parameters", "no", "whether to iterate certain parameters");

  /* order n of approximation in finite differencing, O(h^n) 
     work in progress: test, create good defaults, move to proper place
  */
  AddPar("order_centered", "2", 
	 "finite differencing order of centered stencils [2,4]");
  AddPar("order_advection", "2",
	 "finite differencing order of advection terms [2,4]");
  AddPar("order_RP", "4",
	 "order of interpolation for restriction/prolongation [4,6]");
  AddPar("order_RP_shells","6",  
         "order shell-box interpolation [6,8]");
  AddPar("order_timeinterpolation", "3",
	 "order of interpolation in time for amr [3,4,5,6]");
  AddPar("mg_order", "2",
	 "finite differencing order [2]");
 
  AddPar("order_dissipation", "0",
	 "order of derivatives for dissipation, 0 to turn off [4,6]");
  AddPar("dissipation_factor", "0",
	 "coefficient of dissipation terms, will be divided by 2^order [0]");
  AddPar("dissipation_factor_level0", "-1",
         "coefficient of dissipation terms for level 0 will be divided by 2^order [0], for -1 same value as for other levels is used.");
  
  /* infos about physical objects */
  printf("Adding physical information of objects\n");
  int i=1;
  while( ExistPari("mass",i) ) {
    AddPari("mass", i, "0.", "initial mass of object i");
    AddPari("px",   i, "0.", "initial position in x direction of object i");
    AddPari("py",   i, "0.", "initial position in y direction of object i");
    AddPari("pz",   i, "0.", "initial position in z direction of object i");
    AddPari("mx",   i, "0.", "initial momentum in x direction of object i");
    AddPari("my",   i, "0.", "initial momentum in y direction of object i");
    AddPari("mz",   i, "0.", "initial momentum in z direction of object i");
    AddPari("sx",   i, "0.", "initial spin in x direction of object i");
    AddPari("sy",   i, "0.", "initial spin in y direction of object i");
    AddPari("sz",   i, "0.", "initial spin in z direction of object i");
    i++;
  }
  char str[100];
  sprintf(str,"%d",i-1);
  AddPar("nobjects", str, "Number of BH/NS ...");
  
}





void bam_main_final() 
{
  printf("Removing main\n");
  
}



