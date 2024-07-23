/* bam_io.c */
/* Bernd Bruegmann 12/99 */

#include "bam.h"
#include "output.h"



void bam_output() 
{
  printf("Adding output\n");

  /* functions */
  AddFun(OUTPUT,   write_level, "write data to files");
  
  /* variables */

  /* parameters */
  AddPar("0doutiter", "-1", "when to output based on iterations");
  AddPar("1doutiter", "-1", "when to output based on iterations");
  AddPar("2doutiter", "-1", "when to output based on iterations");
  AddPar("3doutiter", "-1", "when to output based on iterations");

  AddPar("0douttime", "-1", "when to output based on time");
  AddPar("1douttime", "-1", "when to output based on time");
  AddPar("2douttime", "-1", "when to output based on time");
  AddPar("3douttime", "-1", "when to output based on time");

  AddPar("0doutlmin", "0", "range of amr levels to output");
  AddPar("1doutlmin", "0", "range of amr levels to output");
  AddPar("2doutlmin", "0", "range of amr levels to output");
  AddPar("3doutlmin", "0", "range of amr levels to output");

  AddPar("0doutlmax", "10000", "range of amr levels to output");
  AddPar("1doutlmax", "10000", "range of amr levels to output");
  AddPar("2doutlmax", "10000", "range of amr levels to output");
  AddPar("3doutlmax", "10000", "range of amr levels to output");
  
  AddPar("outputinterpolatescheme", "lagrange", "lagrange,WENO");
  AddPar("1doutinterpolate", "no", 
	 "whether to interpolate onto standard coord");
  AddPar("2doutinterpolate", "no", 
	 "whether to interpolate onto standard coord");
  AddPar("3doutinterpolate", "no", 
	 "whether to interpolate onto standard coord");

  AddPar("0doutput", "", "variables to output at point or norm");
  AddPar("1doutput", "", "variables to output along axes");
  AddPar("2doutput", "", "variables to output on coordinate planes");
  AddPar("3doutput", "", "variables to output in 3d");

  AddPar("0doutputall", "no", "whether to output all components");
  AddPar("1doutputall", "no", "whether to output all components");
  AddPar("2doutputall", "no", "whether to output all components");
  AddPar("3doutputall", "no", "whether to output all components");
  
  AddPar("0doutputmode", "max min norm norminf integral normmask integralmask integralouter pmaxmin punc pt",    "");
  AddPar("1doutputmode", "x y z d", "");
  AddPar("2doutputmode", "xy",      "");

  
  AddPar("2doutputr", "",  "variables to output on spheres r=const.");
  AddPar("2doutputr_type", "matrix SphericalDF", 
    "output format on spheres r=const: 1 of cartesian, spherical, matrix, "
    "[cartesian, spherical, matrix] SphericalDF");

  AddPar("2dformat", "vtk binary float",
	 "format for 2d output (opendx,vtk,xdmf,text,binary,float,double)"); 

  AddPar("3dformat", "vtk binary float", "format for 3d output "
	 "(opendx,vtk,xdmf,amrunion,text,binary,float,double,"
	 "add_rotant_points)"); 

  if(Getv("2dformat","xdmf") || Getv("3dformat","xdmf")){
    AddPar("xdmfgrid", "regular", "grid format for xdmf output (regular, curvilinear, rectilinear)"); 
  }
  
  AddPar("outputx0", "0", "origin for output");
  AddPar("outputy0", "0", "origin for output");
  AddPar("outputz0", "0", "origin for output");
  
  AddPar("0doutputpoint", "0.0 0.0 0.0", "output at a point (or several)");
  
}
