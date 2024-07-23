/* bam_boundary.c */
/* Bernd Bruegmann 01/00 */

#include "bam.h"
#include "boundary.h"

void (*set_boundary_elliptic_ptr)(tL *level, tVarList *u);

void bam_boundary() 
{
  printf("Adding boundary\n");

  /* parameters */
  AddPar("boundary_verbose", "yes", "talk about it");
  AddPar("boundary", "", "boundary condition");
 

  /* this is not generallay implemented, only for background a shells radiative */
  AddPar("boundary_order_extrapolate","6","");
  AddPar("boundary_N_extrapolate",    "3","");


  /* radiative boundary */
  if (Getv("boundary", "radiative"))
  {
    AddPar("boundary_radpower", "0", "exponent of non-wave term");
    AddPar("boundary_radconstant", "", 
	   "keep these variables constant for radiative boundary");
    AddPar("boundary_rad_shiftadvection", "no",
           "yes: use advection along beta^i - v*n^i, "
           "no: use advection only along -v*n^i");
  }


  /* robin boundary */
  AddPar("robin_1oR_constant", "0", "1 over r fall-off constant for dirichlet b.c.");
  if (Getv("boundary", "robin")) {
    AddVar("robinindex", "", "index for robin boundary");
    AddVar("robinflag",  "", "temporary flag for robin boundary");
  }

  
  /* excision */
  if (Getv("boundary", "excision")) {

    /* determine when excision mask is set 
       we could add more possibilities for excision_initmask, e.g. the 
       number of the timestep when set_mask_bhdata should be called  
    */
    AddPar("excision_initmask", "PRE_INITIALDATA",
           "when to set excision mask for the first time "
	   "[PRE_INITIALDATA, POST_INITIALDATA]");
    if (Getv("excision_initmask", "POST_INITIALDATA"))
      AddFun(POST_INITIALDATA, set_boundary_flags_excision,                  
             "regridding and thus setting excision flags");
    else
      AddFun(PRE_INITIALDATA, set_boundary_flags_excision,
      	     "regridding and thus setting excision flags");    

  
    AddVar("excisionmask", "", "flags for excision");
    AddVar("excisionindex", "", "index for excision");

    AddPar("excision_type", "lego", "type of excision [lego]");
    AddPar("excision_shape", "sphere", 
	   "shape of excision region [sphere, cube]");
    AddPar("excision_shrink", "0.9", 
	   "shrink by this factor to create buffer zone");
    AddPar("excision_mask_normalfinder", "simple", 
           "how to find the index of the point considered normal, "
           "i.e. the point from which we copy [simple, WT]");
    AddPar("excision_InteriorValue", "VarFarLimit", 
           "Value to be written into points marked EXC_INTERIOR "
           "[any number, VarFarLimit]");
    AddPar("excision_OnlyOnFinestLevel", "yes", 
           "whether excision is used only on the finest level");

    AddPar("excision_set", "bh", 
	   "should excision be set from bh parameter, manually or "
	   "where the lapse is negative? [bh,manual,lapseneg]");
    AddPar("excision_minLapse", "0.0", 
           "value of lapse, below which lapseneg will excise");
 
   /* these parameters are unused, used for manual exision */
    AddPar("excision_r1", "0", "excision mask radius 1");
    AddPar("excision_x1", "0", "excision mask x of center 1");
    AddPar("excision_y1", "0", "excision mask y of center 1");
    AddPar("excision_z1", "0", "excision mask z of center 1");
    AddPar("excision_r2", "0", "excision mask radius 2");
    AddPar("excision_x2", "0", "excision mask x of center 2");
    AddPar("excision_y2", "0", "excision mask y of center 2");
    AddPar("excision_z2", "0", "excision mask z of center 2");
  }

  /* general boundary masks, used only in elliptic solver experiments */
  if (Getv("boundary", "mask")) {
    AddVar("boundarymask", "", "boundary flag (for testing)");
    AddVar("boundarysource", "", "boundary source");
    AddPar("boundary_r1", "0", "boundary mask radius 1 (< 0 for exterior)");
    AddPar("boundary_x1", "0", "boundary mask x of center 1");
    AddPar("boundary_y1", "0", "boundary mask y of center 1");
    AddPar("boundary_z1", "0", "boundary mask z of center 1");
  }

}
