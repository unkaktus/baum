/* boundary.c */
/* Bernd Bruegmann, Wolfgang Tichy 1/03 */

#include "bam.h"
#include "boundary.h"


/* functions provided elsewhere */
void (*Gauge_add_betarot)(tVarList *vl, double sign) = 0;
void (*boundary_fkt_ptr)(tVarList *unew, tVarList *upre, double c, tVarList *ucur) = 0;


/* set boundary pointer */
/* 18.02.2015 at the moment 'boundary_fkt_ptr' is called in set_boundary()
   only with shells and boundary = background */
void boundary_register(void (*f)(tVarList *unew, tVarList *upre, double c, tVarList *ucur)) 
{
  boundary_fkt_ptr = f;
}


/* apply all active boundary conditions after a standard evolution step 
   18.02.2015 this is a new version, hopefully cleaner 
   include an option for background+boxes */
void set_boundary(tVarList *unew, tVarList *upre, double c, tVarList *ucur)
{
  tL *level = unew->level;
  timer_start(level, "set_boundary");
  int i, j;
  static double harmonic_f = 1;
  double v0, var0, varfalloff;

  
  /* ****************************************************************
     radiative boundary condition */
  if ( Getv("boundary", "radiative") ) {
    
    /* there're several 'radiative' bcs : */
    
    /* special radiative boundary condition for shells */
    if (level->shells) {

      int use_centred = Getv("boundary", "centered");

      /* for all variables */
      for (j = 0; j < unew->n; j++) {
        i    = unew->index[j];
        v0   = VarPropSpeed(i);
        var0 = VarFarLimit(i);
  
        /* keep this variable constant, just copy it */
        if (Getv("boundary_radconstant", VarName(i))) v0 = 0;
  
	/* do it */
        if (use_centred) {
          set_boundary_radiative_shells_centered(level, 
						 unew->index[j], upre->index[j], c, ucur->index[j], var0, v0);
        } else {
          set_boundary_radiative_shells(level, 
					unew->index[j], upre->index[j], c, ucur->index[j], var0, v0);
        }
	
      }

      if (use_centred)
        set_boundary_extrapolate_shells(level, unew);
    }
    

    /* radiative boundary condition */
    else if ( Getv("boundary", "norotation") ) {

      /* if (GetvLax("Gauge", "corotation")) add_betarot(ucur, -1); */
      /* if (Gauge_add_betarot) Gauge_add_betarot(ucur, -1); */
  
      /* monopole fall off rate */
      // not implemented yet !
      varfalloff = 0;

      /* for all variables */
      for (j = 0; j < unew->n; j++) {
        i    = unew->index[j];
        v0   = VarPropSpeed(i);
        var0 = VarFarLimit(i);
	
        /* keep this variable constant, just copy it */
        if (Getv("boundary_radconstant", VarName(i))) v0 = 0.;
	
        /* do it */
        set_boundary_radiative(level, 
			       unew->index[j], upre->index[j], c, ucur->index[j], var0, v0);
      }

      /* if (GetvLax("Gauge", "corotation")) add_betarot(ucur, 1); */
      if (Gauge_add_betarot) Gauge_add_betarot(ucur, -1);
    }


    /* radiative boundary condition, with rigid rotation */
    else if ( !Getv("boundary", "norotation") ) {
      
      int ibetax_cur=0;
      int ibetax_far=0;
      int rot=0;
      int betaAdv=Getv("boundary_rad_shiftadvection", "yes");
      
      /* search for shift in ucur */
      for (j = 0; j < ucur->n; j++)
	{
	  char name[8];
	  
	  i = ucur->index[j];
	  snprintf(name, 6, "%s", VarName(i));
	  name[5]=0;
	  if( strcmp(name, "betax")==0)
	    ibetax_cur=i;
	}
      if (GetvLax("Gauge", "corotation"))  
	{
	  rot=1;
	  ibetax_far=Ind("betarotx");
	}
      
      /* for all variables */
      for (j = 0; j < unew->n; j++)
	{
	  i    = unew->index[j];
	  v0   = VarPropSpeed(i);
	  var0 = VarFarLimit(i);
	  
	  /* keep this variable constant, just copy it */
	  if (Getv("boundary_radconstant", VarName(i))) v0 = 0;
	  
	  set_boundary_rot_rad(level,
			       unew->index[j], upre->index[j], c, ucur->index[j], 
			       var0, v0, ibetax_cur, ibetax_far, rot, betaAdv);
      }
      
    }
    
  }
  
  
  /* ****************************************************************
     extrapolate boundary condition */
  if ( Getv("boundary", "extrapolate") )
  {

    bampi_vlsynchronize(unew);
    for (j = 0; j < unew->n; j++) {
      set_boundary_extrapolate(level, unew->index[j]);
    }

  }
  
  
  /* ****************************************************************
     excision boundary */
  if ( Getv("boundary", "excision") ) {

    bampi_vlsynchronize(unew);
    set_boundary_symmetry(level, unew); 
    bampi_vlsynchronize(unew);

    for (j = 0; j < unew->n; j++)
      set_boundary_excision(level, unew->index[j], upre->index[j]);
  }
  

  /* ****************************************************************
     with background values */
  if ( Getv("boundary", "background")  && (level->l==0) ) {
    
    /* special radiative boundary condition for shells */
    if (level->shells) {

      /* sync ... to be sure */
      bampi_vlsynchronize(unew);
      
      /* compute background boundary condition 
	 -> e.g. function z4_boundary_shell() from Z4c project */
      boundary_fkt_ptr(unew,upre, c, ucur);
      
      /* set extrapolated values -> we compute boundary with centered stencils 
	 ... dookie style */ 
      set_boundary_extrapolate_shells(level, unew);
      
      /* sync new values */
      bampi_vlsynchronize(unew); 

    } else {

      /* sync ... to be sure */
      bampi_vlsynchronize(unew);

      /* set extrapolated field values */
      for (j = 0; j < ucur->n; j++) {
	set_boundary_extrapolate(level, ucur->index[j]);
      }

      /* compute background boundary condition 
	 -> e.g. function z4_boundary_box() from Z4c project */
      boundary_fkt_ptr(unew,upre, c, ucur);
      
      /* set extrapolated rhs values */
      for (j = 0; j < unew->n; j++) {
      	set_boundary_extrapolate(level, unew->index[j]);
      }
      
      /* sync new values */
      bampi_vlsynchronize(unew); 

    }

  }
  
  
  /* symmetry boundary 
     has to be after radiative and excision
     can be before or after synchronization, but before is more efficient  */
  set_boundary_symmetry(level, unew); 
  
  timer_stop(level, "set_boundary");
}


/* apply boundary routines for elliptic problems */
void set_boundary_elliptic(tL *level, tVarList *u)
{
  static int flag = -1;

  if (flag == -1) flag = Getv("grid", "spectral");

  if (flag) 
    set_boundary_elliptic_ptr(level, u);
  else {
    set_boundary_robin(level, u);
    set_boundary_symmetry(level, u);
  }
}


/* apply boundary routines for elliptic problems in 1d ignoring
   the cartoon boundary (assumes special 1d elliptic operator 
*/
void set_boundary_elliptic_1d(tL *level, tVarList *u)
{
  double *zp = Ptr(level, "z");

  set_boundary_robin(level, u);
  set_boundary_reflection_one(level, u, zp, level->dz, 5);
}








