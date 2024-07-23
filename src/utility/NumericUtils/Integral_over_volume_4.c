/* integral_over_volume_4.c */
/* Pedro Marronetti 11/06 */

#include "bam.h"
#include "../../main/bampi/bampi.h"

double integral_over_volume_4(tL *level, int inner_surface, int Mintvol, 
	int Mintsurfx, int Mintsurfy, int Mintsurfz){

/* DESCRIPTION

	This subr. calculates the volume integral of the field indexed by Mintvol.
  It also allows the calculation of a "inner surface" integral, to accommodate for
  applications of Gauss theorem, like in the case of the calculation of the ADM
  mass and angular momentum. This surface integral is performed over the sides of
  a cube at one of the maximum refinement levels, with Mintsurfx, Mintsurfy and Mintsurfz
  the integrands corresponding to the x=const, y=const and z=const sides 
  respectively.
  
  	The subr. has three operational modes:
	
  1) Full volume integral: No inner surfaces are considered.
  	 inner_surface = 0
	 
  2) Volume integral + Inner surface integral: Volume integral outside inner cube
     plus surface integral on sides of the inner cube.
  	inner_surface = 1
	 
  3) Volume integral excluding inner cube: Volume integral outside inner cube
     only. This option might be useful for Hamiltonian Constr. convergence tests.
  	inner_surface = 2

  4) Like 2) but using the second, third, fourth highest refinement level.
  	inner_surface = 3, 4, 5
  
  The volume integral is done using an extended trapezoidal algorithm with the
  following characteristics:
	1) The error goes like O(N^(-4)).
	2) The inner boundary coincides with the grid points of the level (closed
	   formula).
	3) The outer boundary is set at half point between two grid points (open
	   formula). This is done to match the volume boundary of the next coarser
	   level. 

  The inner surface integral uses an extended version of the trapezoidal rule 
  with fourth order convergence (Num. Rec. in Fortran eq. 4.1.14).
  
  Future Upgrades: 
  There is a lot of room for improving the efficiency of this function, for 
  instance, by removing the many branches inside the loops.
  
*/


    tG *grid = level->grid;
    tL *lmax = grid->level[grid->lmax];
    tL *lc  = grid->level[grid->lmax-1];
    
    int ii, j, jj, jmin, jmax;
    double global[2], local[2];
    double sum, sumi, dV, dVi, symm_factor=1.0, Tiny = 1.e-10;
    double xmin_i[2][3], xmax_i[2][3], ymin_i[2][3];
    double ymax_i[2][3], zmin_i[2][3], zmax_i[2][3];
    double xmin_o[2][3], xmax_o[2][3], ymin_o[2][3];
    double ymax_o[2][3], zmin_o[2][3], zmax_o[2][3];
    double xmin_s[2][3], xmax_s[2][3], ymin_s[2][3];
    double ymax_s[2][3], zmin_s[2][3], zmax_s[2][3];
    int var_x = Ind("x"), var_y = Ind("y"), var_z = Ind("z");
    int flgprlng = Ind("flagprolong");

    double dx, dy, dz, dxy, dxz, dyz, sgn;

    int number_special_boxes, jc, number_special_boxes_lc; 

    double c[3], o[3];

    double f1 = 3./8., f2 = 7./6., f3 = 23./24.;
    double f11 = f1*f1, f12 = f1*f2, f13 = f1*f3;
    double f22 = f2*f2, f23 = f2*f3, f33 = f3*f3;

    double testvar;

    /* The xmin_s, etc define the inner surface where the surface integral */
    /* will be evaluated. The volume inside the inner surface is excluded  */
    /* from the volume integral and is located in the maximum refinement   */
    /* level. The parameter "offset" controls how far in into this level   */
    /* we set the surface.						   */

    double offset = 1.5, offset2 = offset + 0.5;


    /* Initialization */
    
    jmin = grid->lmin; jmax = grid->lmax;
    
    if (inner_surface > 2 ) {		/* Define Inner Cube as the 2nd, 3rd, 4th */
    lmax = grid->level[grid->lmax-inner_surface+2];	/* highest refinement level 	*/
    lc   = grid->level[grid->lmax-inner_surface+1];
    jmax = grid->lmax - inner_surface+2;
    }
    
    local[0] = 0.0;  local[1] = 0.0;
    
    dx = lmax->dx, dy = lmax->dy, dz = lmax->dz;

    c[0] = 3./8.,   c[1] = 7./6., c[2] = 23./24.;
    o[0] = 13./12., o[1] = 7./8., o[2] = 25./24.;

    /* Set the boundaries of the inner surface */

    if (inner_surface > 0) {

      for (ii = 0; ii < 2; ii++) {
      
        
      if ((lc->nsboxes==0) && (ii==0)) {
        xmin_s[ii][0] = lmax->bbox[0]+offset2*dx;
	xmax_s[ii][0] = lmax->bbox[1]-offset2*dx;
        ymin_s[ii][0] = lmax->bbox[2]+offset2*dy;
	ymax_s[ii][0] = lmax->bbox[3]-offset2*dy;
        zmin_s[ii][0] = lmax->bbox[4]+offset2*dz;
	zmax_s[ii][0] = lmax->bbox[5]-offset2*dz;
      } else if ((lc->nsboxes==0) && (ii==1)) {
        xmin_s[ii][0] = 0.; xmax_s[ii][0] = 0.;
        ymin_s[ii][0] = 0.; ymax_s[ii][0] = 0.;
        zmin_s[ii][0] = 0.; zmax_s[ii][0] = 0.;
      } else {
        xmin_s[ii][0] = lc->sbox[ii].bbox[0]+offset*dx;
	xmax_s[ii][0] = lc->sbox[ii].bbox[1]-offset*dx;
        ymin_s[ii][0] = lc->sbox[ii].bbox[2]+offset*dy;
	ymax_s[ii][0] = lc->sbox[ii].bbox[3]-offset*dy;
        zmin_s[ii][0] = lc->sbox[ii].bbox[4]+offset*dz;
	zmax_s[ii][0] = lc->sbox[ii].bbox[5]-offset*dz;
      }
      
        xmin_s[ii][1] = xmin_s[ii][0]+dx; xmin_s[ii][2] = xmin_s[ii][1]+dx;
	xmax_s[ii][1] = xmax_s[ii][0]-dx; xmax_s[ii][2] = xmax_s[ii][1]-dx;
        ymin_s[ii][1] = ymin_s[ii][0]+dy; ymin_s[ii][2] = ymin_s[ii][1]+dy;
	ymax_s[ii][1] = ymax_s[ii][0]-dy; ymax_s[ii][2] = ymax_s[ii][1]-dy;
        zmin_s[ii][1] = zmin_s[ii][0]+dz; zmin_s[ii][2] = zmin_s[ii][1]+dz;
	zmax_s[ii][1] = zmax_s[ii][0]-dz; zmax_s[ii][2] = zmax_s[ii][1]-dz;
	
	if (xmin_s[ii][0] > xmax_s[ii][0]) {
	  xmin_s[ii][0] = 0.0; xmin_s[ii][1] = 0.0; xmin_s[ii][2] = 0.0;
	  xmax_s[ii][0] = 0.0; xmax_s[ii][1] = 0.0; xmax_s[ii][2] = 0.0;
	}

	if (ymin_s[ii][0] > ymax_s[ii][0]) {
	  ymin_s[ii][0] = 0.0; ymin_s[ii][1] = 0.0; ymin_s[ii][2] = 0.0;
	  ymax_s[ii][0] = 0.0; ymax_s[ii][1] = 0.0; ymax_s[ii][2] = 0.0;
	}

	if (zmin_s[ii][0] > zmax_s[ii][0]) {
	  zmin_s[ii][0] = 0.0; zmin_s[ii][1] = 0.0; zmin_s[ii][2] = 0.0;
	  zmax_s[ii][0] = 0.0; zmax_s[ii][1] = 0.0; zmax_s[ii][2] = 0.0;
	}
      }

    } else {
    
      for (ii = 0; ii < 2; ii++) {
        xmin_s[ii][0] = 0.0; xmax_s[ii][0] = 0.0;
        ymin_s[ii][0] = 0.0; ymax_s[ii][0] = 0.0;
        zmin_s[ii][0] = 0.0; zmax_s[ii][0] = 0.0;
        xmin_s[ii][1] = 0.0; xmax_s[ii][1] = 0.0;
        ymin_s[ii][1] = 0.0; ymax_s[ii][1] = 0.0;
        zmin_s[ii][1] = 0.0; zmax_s[ii][1] = 0.0;
        xmin_s[ii][2] = 0.0; xmax_s[ii][2] = 0.0;
        ymin_s[ii][2] = 0.0; ymax_s[ii][2] = 0.0;
        zmin_s[ii][2] = 0.0; zmax_s[ii][2] = 0.0;
      }
      
    }
    
    /* Correct them for symmetry */
    
    if (Getv("grid", "bitant")) {
      symm_factor = 2.0; 
      for (ii = 0; ii < 2; ii++) { 
      
	if ((zmin_s[ii][0]-offset*dz) <= 0.0) {
	  zmin_s[ii][0] = 0.0; zmin_s[ii][1] = 0.0; zmin_s[ii][2] = 0.0; }
	 
	  if (zmax_s[ii][0] <= 0.0) {
            zmax_s[ii][0] = 0.0; zmax_s[ii][1] = 0.0; zmax_s[ii][2] = 0.0; }
          
      }
      
    } else if (Getv("grid", "quadrant")) {
      symm_factor = 4.0;
      for (ii = 0; ii < 2; ii++) {
        if ((ymin_s[ii][0]-offset*dy) <= 0.0) {
          ymin_s[ii][0] = 0.0; ymin_s[ii][1] = 0.0; ymin_s[ii][2] = 0.0; }
	  
	if (ymax_s[ii][0] <= 0.0) {
          ymax_s[ii][0] = 0.0; ymax_s[ii][1] = 0.0; ymax_s[ii][2] = 0.0; }
	
	if ((zmin_s[ii][0]-offset*dz) <= 0.0) {
	  zmin_s[ii][0] = 0.0; zmin_s[ii][1] = 0.0; zmin_s[ii][2] = 0.0; }
	 
	if (zmax_s[ii][0] <= 0.0) {
          zmax_s[ii][0] = 0.0; zmax_s[ii][1] = 0.0; zmax_s[ii][2] = 0.0; }
      }
    
    } else if (Getv("grid", "octant")) {
      symm_factor = 8.0;
      for (ii = 0; ii < 2; ii++) { 
        if ((xmin_s[ii][0]-offset*dx) <= 0.0) {
          xmin_s[ii][0] = 0.0; xmin_s[ii][1] = 0.0; xmin_s[ii][2] = 0.0; }
	
	if (xmax_s[ii][0] <= 0.0) {
          xmax_s[ii][0] = 0.0; xmax_s[ii][1] = 0.0; xmax_s[ii][2] = 0.0; } 
	  
        if ((ymin_s[ii][0]-offset*dy) <= 0.0) {
          ymin_s[ii][0] = 0.0; ymin_s[ii][1] = 0.0; ymin_s[ii][2] = 0.0; }
	
	if (ymax_s[ii][0] <= 0.0) {
          ymax_s[ii][0] = 0.0; ymax_s[ii][1] = 0.0; ymax_s[ii][2] = 0.0; }
	
	if ((zmin_s[ii][0]-offset*dz) <= 0.0) {
	  zmin_s[ii][0] = 0.0; zmin_s[ii][1] = 0.0; zmin_s[ii][2] = 0.0; }
	 
	if (zmax_s[ii][0] <= 0.0) {
          zmax_s[ii][0] = 0.0; zmax_s[ii][1] = 0.0; zmax_s[ii][2] = 0.0; } 
      }
    }
    
    /* First, we take care of the volume integral */
    
    for(j = jmin; j < jmax; j++) { 
     
        tL *l = grid->level[j];
        tL *lf = grid->level[j+1];
	dx = l->dx, dy = l->dy, dz = l->dz;
        number_special_boxes    = l->nsboxes;
	jc = ( j == 0 ) ? 0 : j-1;
        number_special_boxes_lc = grid->level[jc]->nsboxes;
        sum = 0.0, sumi =0.0;
	
	/* If we are not at the inner cube level, we need to exclude the	*/
	/* area covered by the next finer level from the integration 		*/
	/* domain. Note that the fine level box is global.			*/
	
	if (j < jmax ) {
	
          for (ii = 0; ii < 2; ii++) {
            
            if(lc->nsboxes==0) {
              xmin_i[ii][0] = lf->bbox[0]+1.25*dx, xmax_i[ii][0] = lf->bbox[1]-1.25*dx;
              ymin_i[ii][0] = lf->bbox[2]+1.25*dy, ymax_i[ii][0] = lf->bbox[3]-1.25*dy;
              zmin_i[ii][0] = lf->bbox[4]+1.25*dz, zmax_i[ii][0] = lf->bbox[5]-1.25*dz;
            } else {
              xmin_i[ii][0] = l->sbox[ii].bbox[0]+dx, xmax_i[ii][0] = l->sbox[ii].bbox[1]-dx;
              ymin_i[ii][0] = l->sbox[ii].bbox[2]+dy, ymax_i[ii][0] = l->sbox[ii].bbox[3]-dy;
              zmin_i[ii][0] = l->sbox[ii].bbox[4]+dz, zmax_i[ii][0] = l->sbox[ii].bbox[5]-dz;
            }

          }
	} else {
	
          for (ii = 0; ii < 2; ii++) {
	    xmin_i[ii][0] = xmin_s[ii][0], xmax_i[ii][0] = xmax_s[ii][0];
	    ymin_i[ii][0] = ymin_s[ii][0], ymax_i[ii][0] = ymax_s[ii][0];
	    zmin_i[ii][0] = zmin_s[ii][0], zmax_i[ii][0] = zmax_s[ii][0];
	  }

	}
	
        for (ii = 0; ii < 2; ii++) {

	  xmin_i[ii][1] = xmin_i[ii][0]+dx, xmin_i[ii][2] = xmin_i[ii][1]+dx;
	  xmax_i[ii][1] = xmax_i[ii][0]-dx, xmax_i[ii][2] = xmax_i[ii][1]-dx;
	  ymin_i[ii][1] = ymin_i[ii][0]+dy, ymin_i[ii][2] = ymin_i[ii][1]+dy;
	  ymax_i[ii][1] = ymax_i[ii][0]-dy, ymax_i[ii][2] = ymax_i[ii][1]-dy;
	  zmin_i[ii][1] = zmin_i[ii][0]+dz, zmin_i[ii][2] = zmin_i[ii][1]+dz;
	  zmax_i[ii][1] = zmax_i[ii][0]-dz, zmax_i[ii][2] = zmax_i[ii][1]-dz;

          xmin_o[ii][0] = l->bbox[0]+3.*dx, xmax_o[ii][0] = l->bbox[1]-3.*dx;
          ymin_o[ii][0] = l->bbox[2]+3.*dy, ymax_o[ii][0] = l->bbox[3]-3.*dy;
          zmin_o[ii][0] = l->bbox[4]+3.*dz, zmax_o[ii][0] = l->bbox[5]-3.*dz;
	  
	  xmin_o[ii][1] = xmin_o[ii][0]+dx, xmin_o[ii][2] = xmin_o[ii][1]+dx;
	  xmax_o[ii][1] = xmax_o[ii][0]-dx, xmax_o[ii][2] = xmax_o[ii][1]-dx;
	  ymin_o[ii][1] = ymin_o[ii][0]+dy, ymin_o[ii][2] = ymin_o[ii][1]+dy;
	  ymax_o[ii][1] = ymax_o[ii][0]-dy, ymax_o[ii][2] = ymax_o[ii][1]-dy;
	  zmin_o[ii][1] = zmin_o[ii][0]+dz, zmin_o[ii][2] = zmin_o[ii][1]+dz;
	  zmax_o[ii][1] = zmax_o[ii][0]-dz, zmax_o[ii][2] = zmax_o[ii][1]-dz;

	  if (xmin_i[ii][0] > xmax_i[ii][0]) {
	    xmin_i[ii][0] = 0.0; xmin_i[ii][1] = 0.0; xmin_i[ii][2] = 0.0;
	    xmax_i[ii][0] = 0.0; xmax_i[ii][1] = 0.0; xmax_i[ii][2] = 0.0;
	  }
	
	  if (ymin_i[ii][0] > ymax_i[ii][0]) {
	    ymin_i[ii][0] = 0.0; ymin_i[ii][1] = 0.0; ymin_i[ii][2] = 0.0;
	    ymax_i[ii][0] = 0.0; ymax_i[ii][1] = 0.0; ymax_i[ii][2] = 0.0;
	  }
	
	  if (zmin_i[ii][0] > zmax_i[ii][0]) {
	    zmin_i[ii][0] = 0.0; zmin_i[ii][1] = 0.0; zmin_i[ii][2] = 0.0;
	    zmax_i[ii][0] = 0.0; zmax_i[ii][1] = 0.0; zmax_i[ii][2] = 0.0;
	  }
	}

        /* Correct them for symmetry */
    
        if (j < jmax) {
        if (Getv("grid", "bitant")) {
          for (ii = 0; ii < 2; ii++) { 
	    if ((zmin_i[ii][0]-dz) <= 0.0) {
	      zmin_i[ii][0] = 0.0; zmin_i[ii][1] = 0.0; zmin_i[ii][2] = 0.0; }
	 
	    if (zmax_i[ii][0] <= 0.0) {
              zmax_i[ii][0] = 0.0; zmax_i[ii][1] = 0.0; zmax_i[ii][2] = 0.0; }

	    if (l->bbox[4] <= 0.0) {
	      zmin_o[ii][0] = 0.0; zmin_o[ii][1] = 0.0; zmin_o[ii][2] = 0.0; }
	 
	    if (zmax_o[ii][0] <= 0.0) {
              zmax_o[ii][0] = 0.0; zmax_o[ii][1] = 0.0; zmax_o[ii][2] = 0.0; }
        }
      
        } else if (Getv("grid", "quadrant")) {
          for (ii = 0; ii < 2; ii++) {
            if ((ymin_i[ii][0]-dy) <= 0.0) {
              ymin_i[ii][0] = 0.0; ymin_i[ii][1] = 0.0; ymin_i[ii][2] = 0.0; }
	  
	    if (ymax_i[ii][0] <= 0.0) {
              ymax_i[ii][0] = 0.0; ymax_i[ii][1] = 0.0; ymax_i[ii][2] = 0.0; }
	
            if (l->bbox[2] <= 0.0) {
              ymin_o[ii][0] = 0.0; ymin_o[ii][1] = 0.0; ymin_o[ii][2] = 0.0; }
	  
	    if (ymax_o[ii][0] <= 0.0) {
              ymax_o[ii][0] = 0.0; ymax_o[ii][1] = 0.0; ymax_o[ii][2] = 0.0; }
	
	    if ((zmin_i[ii][0]-dz) <= 0.0) {
	      zmin_i[ii][0] = 0.0; zmin_i[ii][1] = 0.0; zmin_i[ii][2] = 0.0; }
	 
	    if (zmax_i[ii][0] <= 0.0) {
              zmax_i[ii][0] = 0.0; zmax_i[ii][1] = 0.0; zmax_i[ii][2] = 0.0; }

	    if (l->bbox[4] <= 0.0) {
	      zmin_o[ii][0] = 0.0; zmin_o[ii][1] = 0.0; zmin_o[ii][2] = 0.0; }
	 
	    if (zmax_o[ii][0] <= 0.0) {
              zmax_o[ii][0] = 0.0; zmax_o[ii][1] = 0.0; zmax_o[ii][2] = 0.0; }
        }
    
        } else if (Getv("grid", "octant")) {
          for (ii = 0; ii < 2; ii++) { 
            if ((xmin_i[ii][0]-dx) <= 0.0) {
              xmin_i[ii][0] = 0.0; xmin_i[ii][1] = 0.0; xmin_i[ii][2] = 0.0; }
	
	    if (xmax_i[ii][0] <= 0.0) {
              xmax_i[ii][0] = 0.0; xmax_i[ii][1] = 0.0; xmax_i[ii][2] = 0.0; } 

            if (l->bbox[0] <= 0.0) {
              xmin_o[ii][0] = 0.0; xmin_o[ii][1] = 0.0; xmin_o[ii][2] = 0.0; }
	
	    if (xmax_o[ii][0] <= 0.0) {
              xmax_o[ii][0] = 0.0; xmax_o[ii][1] = 0.0; xmax_o[ii][2] = 0.0; } 
	  
            if ((ymin_i[ii][0]-dy) <= 0.0) {
              ymin_i[ii][0] = 0.0; ymin_i[ii][1] = 0.0; ymin_i[ii][2] = 0.0; }
	  
	    if (ymax_i[ii][0] <= 0.0) {
              ymax_i[ii][0] = 0.0; ymax_i[ii][1] = 0.0; ymax_i[ii][2] = 0.0; }
	
            if (l->bbox[2] <= 0.0) {
              ymin_o[ii][0] = 0.0; ymin_o[ii][1] = 0.0; ymin_o[ii][2] = 0.0; }
	  
	    if (ymax_o[ii][0] <= 0.0) {
              ymax_o[ii][0] = 0.0; ymax_o[ii][1] = 0.0; ymax_o[ii][2] = 0.0; }
	
	    if ((zmin_i[ii][0]-dz) <= 0.0) {
	      zmin_i[ii][0] = 0.0; zmin_i[ii][1] = 0.0; zmin_i[ii][2] = 0.0; }
	 
	    if (zmax_i[ii][0] <= 0.0) {
              zmax_i[ii][0] = 0.0; zmax_i[ii][1] = 0.0; zmax_i[ii][2] = 0.0; }

	    if (l->bbox[4] <= 0.0) {
	      zmin_o[ii][0] = 0.0; zmin_o[ii][1] = 0.0; zmin_o[ii][2] = 0.0; }
	 
	    if (zmax_o[ii][0] <= 0.0) {
              zmax_o[ii][0] = 0.0; zmax_o[ii][1] = 0.0; zmax_o[ii][2] = 0.0; }
          }
        }}
    
	/* Now we calculate the volume integral. We do two sums: one 	*/
	/* over the whole level (sum) and one over the part that 	*/
	/* overlaps with the next higher res level (sumi). The volume	*/
	/* integral is computed as sum-sumi. The variables dV and dVi	*/
	/* scale the cell volume according to whether the point ijk is	*/
	/* inside the box, on a surface, on a line or a is corner point.*/

        forinnerpoints_ijk(l) {

	  double u   = l->v[Mintvol][ijk];
	  double xx  = l->v[var_x][ijk];
          double yy  = l->v[var_y][ijk];
          double zz  = l->v[var_z][ijk];
	  double fp  = l->v[flgprlng][ijk];
	  
          dV = 1.0, dVi = 1.0;

	  if ((fp != 1) && (xmin_o[0][0] < ( xx+Tiny ) && ( xx-Tiny ) < xmax_o[0][0] &&
	                    ymin_o[0][0] < ( yy+Tiny ) && ( yy-Tiny ) < ymax_o[0][0] &&
	                    zmin_o[0][0] < ( zz+Tiny ) && ( zz-Tiny ) < zmax_o[0][0] ) ) {
	    
            for (jj = 0; jj < 3; jj++) { 

	      if (fabs(xmin_o[0][jj] - xx) < Tiny || fabs(xmax_o[0][jj] - xx) < Tiny ) dV *= o[jj];
	      if (fabs(ymin_o[0][jj] - yy) < Tiny || fabs(ymax_o[0][jj] - yy) < Tiny ) dV *= o[jj];
	      if (fabs(zmin_o[0][jj] - zz) < Tiny || fabs(zmax_o[0][jj] - zz) < Tiny ) dV *= o[jj];
   	      }

            sum += u * dV;

          }

	  if ((fp != 1) && (xmin_i[0][0] < ( xx+Tiny ) && ( xx-Tiny ) < xmax_i[0][0] &&
	                    ymin_i[0][0] < ( yy+Tiny ) && ( yy-Tiny ) < ymax_i[0][0] &&
	                    zmin_i[0][0] < ( zz+Tiny ) && ( zz-Tiny ) < zmax_i[0][0] ) ) {
	    
            for (jj = 0; jj < 3; jj++) { 

	      if (fabs(xmin_i[0][jj] - xx) < Tiny || fabs(xmax_i[0][jj] - xx) < Tiny ) dVi *= c[jj];
	      if (fabs(ymin_i[0][jj] - yy) < Tiny || fabs(ymax_i[0][jj] - yy) < Tiny ) dVi *= c[jj];
	      if (fabs(zmin_i[0][jj] - zz) < Tiny || fabs(zmax_i[0][jj] - zz) < Tiny ) dVi *= c[jj];
   	      }

            sumi += u * dVi;

          }

	  if ( number_special_boxes > 1 || number_special_boxes_lc > 1) {
	    
	    if (xmin_o[1][0] != xmin_o[0][0] || xmax_o[1][0] != xmax_o[0][0] 
             || ymin_o[1][0] != ymin_o[0][0] || ymax_o[1][0] != ymax_o[0][0] ) {

	      if ((fp != 1) && (xmin_o[1][0] < ( xx+Tiny ) && ( xx-Tiny ) < xmax_o[1][0] &&
	                        ymin_o[1][0] < ( yy+Tiny ) && ( yy-Tiny ) < ymax_o[1][0] &&
	                        zmin_o[1][0] < ( zz+Tiny ) && ( zz-Tiny ) < zmax_o[1][0] ) ) {
	    
                for (jj = 0; jj < 3; jj++) { 

	          if (fabs(xmin_o[1][jj] - xx) < Tiny || fabs(xmax_o[1][jj] - xx) < Tiny ) dV *= o[jj];
	          if (fabs(ymin_o[1][jj] - yy) < Tiny || fabs(ymax_o[1][jj] - yy) < Tiny ) dV *= o[jj];
	          if (fabs(zmin_o[1][jj] - zz) < Tiny || fabs(zmax_o[1][jj] - zz) < Tiny ) dV *= o[jj];
   	        }

                sum += u * dV;
              }
            }

	    if ((fp != 1) && (xmin_i[1][0] < ( xx+Tiny ) && ( xx-Tiny ) < xmax_i[1][0] &&
	                      ymin_i[1][0] < ( yy+Tiny ) && ( yy-Tiny ) < ymax_i[1][0] &&
	                      zmin_i[1][0] < ( zz+Tiny ) && ( zz-Tiny ) < zmax_i[1][0] ) ) {
	    
            for (jj = 0; jj < 3; jj++) { 

	      if (fabs(xmin_i[1][jj] - xx) < Tiny || fabs(xmax_i[1][jj] - xx) < Tiny ) dVi *= c[jj];
	      if (fabs(ymin_i[1][jj] - yy) < Tiny || fabs(ymax_i[1][jj] - yy) < Tiny ) dVi *= c[jj];
	      if (fabs(zmin_i[1][jj] - zz) < Tiny || fabs(zmax_i[1][jj] - zz) < Tiny ) dVi *= c[jj];
   	      }

              sumi += u * dVi;

   	    }
          }

      } endfor_ijk;
      
	local[0] += (sum - sumi) * dx * dy * dz; 
    
    } 
  
  
    /* Second, we take care of the surface integral */
    /* Note that we do the integral in each one of the special boxes */
    
    sum = 0.0;

    /* do we need the next line ? */
    dx = lmax->dx, dy = lmax->dy, dz = lmax->dz;
    dxy = dx * dy, dxz = dx * dz, dyz = dy * dz;

    if (inner_surface == 1 || inner_surface > 2) {
      
      forinnerpoints_ijk(lmax) {
      
        double ux  = lmax->v[Mintsurfx][ijk];
        double uy  = lmax->v[Mintsurfy][ijk];
        double uz  = lmax->v[Mintsurfz][ijk];
        double xx  = lmax->v[var_x][ijk];
        double yy  = lmax->v[var_y][ijk];
        double zz  = lmax->v[var_z][ijk];
	
	for (ii = 0; ii < 2; ii++) {
	 
        if ( fabs(xx - xmin_s[ii][0]) < Tiny || fabs(xx - xmax_s[ii][0]) < Tiny ) {
	  
	  if (fabs(xx - xmin_s[ii][0]) < Tiny) sgn = -1.0;
	  if (fabs(xx - xmax_s[ii][0]) < Tiny) sgn =  1.0;
	  	   
          if ( (ymin_s[ii][2] + Tiny) < yy && (yy + Tiny) < ymax_s[ii][2]  
	    && (zmin_s[ii][2] + Tiny) < zz && (zz + Tiny) < zmax_s[ii][2] )  {
	    sum += sgn * ux * dyz; }
	   
          if ( (fabs(ymin_s[ii][0] - yy) < Tiny || fabs(yy - ymax_s[ii][0]) < Tiny )  
	     && (zmin_s[ii][2] + Tiny) < zz  && (zz + Tiny) < zmax_s[ii][2] )  {
	     sum += sgn * f1 * ux * dyz; }
	  
	  if ( (fabs(ymin_s[ii][1] - yy) < Tiny || fabs(yy - ymax_s[ii][1]) < Tiny )  
	     && (zmin_s[ii][2] + Tiny) < zz  && (zz + Tiny) < zmax_s[ii][2] )  {
	     sum += sgn * f2 * ux * dyz; }
	  
	  if ( (fabs(ymin_s[ii][2] - yy) < Tiny || fabs(yy - ymax_s[ii][2]) < Tiny )  
	     && (zmin_s[ii][2] + Tiny) < zz  && (zz + Tiny) < zmax_s[ii][2] )  {
	     sum += sgn * f3 * ux * dyz; }
	   
          if ( (ymin_s[ii][2] + Tiny) < yy  && (yy + Tiny) < ymax_s[ii][2]  
	   && (fabs(zmin_s[ii][0] - zz) < Tiny || fabs(zz - zmax_s[ii][0]) < Tiny) )  {
	   sum += sgn * f1 * ux * dyz; }
	  
	  if ( (ymin_s[ii][2] + Tiny) < yy && (yy + Tiny) < ymax_s[ii][2]  
	   && (fabs(zmin_s[ii][1] - zz) < Tiny || fabs(zz - zmax_s[ii][1]) < Tiny) )  {
	   sum += sgn * f2 * ux * dyz; }
	  
	  if ( (ymin_s[ii][2] + Tiny) < yy  && (yy + Tiny) < ymax_s[ii][2]  
	   && (fabs(zmin_s[ii][2] - zz) < Tiny || fabs(zz - zmax_s[ii][2]) < Tiny) )  {
	   sum += sgn * f3 * ux * dyz; }
          
	  if ( (fabs(ymin_s[ii][0] - yy) < Tiny || fabs(yy - ymax_s[ii][0]) < Tiny) 
	    && (fabs(zmin_s[ii][0] - zz) < Tiny || fabs(zz - zmax_s[ii][0]) < Tiny) )  {
	    sum += sgn * f11 * ux * dyz; }
	  
	  if ( (fabs(ymin_s[ii][1] - yy) < Tiny || fabs(yy - ymax_s[ii][1]) < Tiny) 
	    && (fabs(zmin_s[ii][0] - zz) < Tiny || fabs(zz - zmax_s[ii][0]) < Tiny) )  {
	    sum += sgn * f12 * ux * dyz; }
	  
	  if ( (fabs(ymin_s[ii][2] - yy) < Tiny || fabs(yy - ymax_s[ii][2]) < Tiny) 
	    && (fabs(zmin_s[ii][0] - zz) < Tiny || fabs(zz - zmax_s[ii][0]) < Tiny) )  {
	    sum += sgn * f13 * ux * dyz; }
	  
	  if ( (fabs(ymin_s[ii][0] - yy) < Tiny || fabs(yy - ymax_s[ii][0]) < Tiny) 
	    && (fabs(zmin_s[ii][1] - zz) < Tiny || fabs(zz - zmax_s[ii][1]) < Tiny) )  {
	    sum += sgn * f12 * ux * dyz; }
	  
	  if ( (fabs(ymin_s[ii][1] - yy) < Tiny || fabs(yy - ymax_s[ii][1]) < Tiny) 
	    && (fabs(zmin_s[ii][1] - zz) < Tiny || fabs(zz - zmax_s[ii][1]) < Tiny) )  {
	    sum += sgn * f22 * ux * dyz; }
	  
	  if ( (fabs(ymin_s[ii][2] - yy) < Tiny || fabs(yy - ymax_s[ii][2]) < Tiny) 
	    && (fabs(zmin_s[ii][1] - zz) < Tiny || fabs(zz - zmax_s[ii][1]) < Tiny) )  {
	    sum += sgn * f23 * ux * dyz; }
	  
	  if ( (fabs(ymin_s[ii][0] - yy) < Tiny || fabs(yy - ymax_s[ii][0]) < Tiny) 
	    && (fabs(zmin_s[ii][2] - zz) < Tiny || fabs(zz - zmax_s[ii][2]) < Tiny) )  {
	    sum += sgn * f13 * ux * dyz; }
	  
	  if ( (fabs(ymin_s[ii][1] - yy) < Tiny || fabs(yy - ymax_s[ii][1]) < Tiny) 
	    && (fabs(zmin_s[ii][2] - zz) < Tiny || fabs(zz - zmax_s[ii][2]) < Tiny) )  {
	    sum += sgn * f23 * ux * dyz; }
	  
	  if ( (fabs(ymin_s[ii][2] - yy) < Tiny || fabs(yy - ymax_s[ii][2]) < Tiny) 
	    && (fabs(zmin_s[ii][2] - zz) < Tiny || fabs(zz - zmax_s[ii][2]) < Tiny) )  {
	    sum += sgn * f33 * ux * dyz; }

        } 
	 
	/*-----------------------*/ 
	   
        if ( fabs(yy - ymin_s[ii][0]) < Tiny || fabs(yy - ymax_s[ii][0]) < Tiny) {
	   
	  if (fabs(yy - ymin_s[ii][0]) < Tiny) sgn = -1.0;
	  if (fabs(yy - ymax_s[ii][0]) < Tiny) sgn =  1.0;
	  
          if ( (xmin_s[ii][2] + Tiny) < xx && (xx + Tiny) < xmax_s[ii][2]  
	    && (zmin_s[ii][2] + Tiny) < zz && (zz + Tiny) < zmax_s[ii][2] )  {
	    sum += sgn * uy * dxz; }
	   
          if ( (fabs(xmin_s[ii][0] - xx) < Tiny || fabs(xx - xmax_s[ii][0]) < Tiny)  
	     && (zmin_s[ii][2] + Tiny) < zz  && (zz + Tiny) < zmax_s[ii][2] )  {
	     sum += sgn * f1 * uy * dxz; }
	  
	  if ( (fabs(xmin_s[ii][1] - xx) < Tiny || fabs(xx - xmax_s[ii][1]) < Tiny)  
	     && (zmin_s[ii][2] + Tiny) < zz  && (zz + Tiny) < zmax_s[ii][2] )  {
	     sum += sgn * f2 * uy * dxz; }
	  
	  if ( (fabs(xmin_s[ii][2] - xx) < Tiny || fabs(xx - xmax_s[ii][2]) < Tiny)  
	     && (zmin_s[ii][2] + Tiny) < zz  && (zz + Tiny) < zmax_s[ii][2] )  {
	     sum += sgn * f3 * uy * dxz; }
	   
          if ( (xmin_s[ii][2] + Tiny) < xx  && (xx + Tiny) < xmax_s[ii][2]  
	   && (fabs(zmin_s[ii][0] - zz) < Tiny || fabs(zz - zmax_s[ii][0]) < Tiny) )  {
	   sum += sgn * f1 * uy * dxz; }
	  
	  if ( (xmin_s[ii][2] + Tiny) < xx && (xx + Tiny) < xmax_s[ii][2]  
	   && (fabs(zmin_s[ii][1] - zz) < Tiny || fabs(zz - zmax_s[ii][1]) < Tiny) )  {
	   sum += sgn * f2 * uy * dxz; }
	  
	  if ( (xmin_s[ii][2] + Tiny) < xx  && (xx + Tiny) < xmax_s[ii][2]  
	   && (fabs(zmin_s[ii][2] - zz) < Tiny || fabs(zz - zmax_s[ii][2]) < Tiny) )  {
	   sum += sgn * f3 * uy * dxz; }
          
	  if ( (fabs(xmin_s[ii][0] - xx) < Tiny || fabs(xx - xmax_s[ii][0]) < Tiny) 
	    && (fabs(zmin_s[ii][0] - zz) < Tiny || fabs(zz - zmax_s[ii][0]) < Tiny) )  {
	    sum += sgn * f11 * uy * dxz; }
	  
	  if ( (fabs(xmin_s[ii][1] - xx) < Tiny || fabs(xx - xmax_s[ii][1]) < Tiny) 
	    && (fabs(zmin_s[ii][0] - zz) < Tiny || fabs(zz - zmax_s[ii][0]) < Tiny) )  {
	    sum += sgn * f12 * uy * dxz; }
	  
	  if ( (fabs(xmin_s[ii][2] - xx) < Tiny || fabs(xx - xmax_s[ii][2]) < Tiny) 
	    && (fabs(zmin_s[ii][0] - zz) < Tiny || fabs(zz - zmax_s[ii][0]) < Tiny) )  {
	    sum += sgn * f13 * uy * dxz; }
	  
	  if ( (fabs(xmin_s[ii][0] - xx) < Tiny || fabs(xx - xmax_s[ii][0]) < Tiny) 
	    && (fabs(zmin_s[ii][1] - zz) < Tiny || fabs(zz - zmax_s[ii][1]) < Tiny) )  {
	    sum += sgn * f12 * uy * dxz; }
	  
	  if ( (fabs(xmin_s[ii][1] - xx) < Tiny || fabs(xx - xmax_s[ii][1]) < Tiny) 
	    && (fabs(zmin_s[ii][1] - zz) < Tiny || fabs(zz - zmax_s[ii][1]) < Tiny) )  {
	    sum += sgn * f22 * uy * dxz; }
	  
	  if ( (fabs(xmin_s[ii][2] - xx) < Tiny || fabs(xx - xmax_s[ii][2]) < Tiny) 
	    && (fabs(zmin_s[ii][1] - zz) < Tiny || fabs(zz - zmax_s[ii][1]) < Tiny) )  {
	    sum += sgn * f23 * uy * dxz; }
	  
	  if ( (fabs(xmin_s[ii][0] - xx) < Tiny || fabs(xx - xmax_s[ii][0]) < Tiny) 
	    && (fabs(zmin_s[ii][2] - zz) < Tiny || fabs(zz - zmax_s[ii][2]) < Tiny) )  {
	    sum += sgn * f13 * uy * dxz; }
	  
	  if ( (fabs(xmin_s[ii][1] - xx) < Tiny || fabs(xx - xmax_s[ii][1]) < Tiny) 
	    && (fabs(zmin_s[ii][2] - zz) < Tiny || fabs(zz - zmax_s[ii][2]) < Tiny) )  {
	    sum += sgn * f23 * uy * dxz; }
	  
	  if ( (fabs(xmin_s[ii][2] - xx) < Tiny || fabs(xx - xmax_s[ii][2]) < Tiny) 
	    && (fabs(zmin_s[ii][2] - zz) < Tiny || fabs(zz - zmax_s[ii][2]) < Tiny) )  {
	    sum += sgn * f33 * uy * dxz; }

        } 
	   
	/*-----------------------*/ 
		   
        if ( fabs(zz - zmin_s[ii][0]) < Tiny || fabs(zz - zmax_s[ii][0]) < Tiny) {
	  
	  if (fabs(zz - zmin_s[ii][0]) < Tiny) sgn = -1.0;
	  if (fabs(zz - zmax_s[ii][0]) < Tiny) sgn =  1.0;
	  
          if ( (ymin_s[ii][2] + Tiny) < yy && (yy + Tiny) < ymax_s[ii][2]  
	    && (xmin_s[ii][2] + Tiny) < xx && (xx + Tiny) < xmax_s[ii][2] )  {
	    sum += sgn * uz * dxy; }
	   
          if ( (fabs(ymin_s[ii][0] - yy) < Tiny || fabs(yy - ymax_s[ii][0]) < Tiny)  
	     && (xmin_s[ii][2] + Tiny) < xx  && (xx + Tiny) < xmax_s[ii][2] )  {
	     sum += sgn * f1 * uz * dxy; }
	  
	  if ( (fabs(ymin_s[ii][1] - yy) < Tiny || fabs(yy - ymax_s[ii][1]) < Tiny)  
	     && (xmin_s[ii][2] + Tiny) < xx  && (xx + Tiny) < xmax_s[ii][2] )  {
	     sum += sgn * f2 * uz * dxy; }
	  
	  if ( (fabs(ymin_s[ii][2] - yy) < Tiny || fabs(yy - ymax_s[ii][2]) < Tiny)  
	     && (xmin_s[ii][2] + Tiny) < xx  && (xx + Tiny) < xmax_s[ii][2] )  {
	     sum += sgn * f3 * uz * dxy; }
	   
          if ( (ymin_s[ii][2] + Tiny) < yy  && (yy + Tiny) < ymax_s[ii][2]  
	   && (fabs(xmin_s[ii][0] - xx) < Tiny || fabs(xx - xmax_s[ii][0]) < Tiny) )  {
	   sum += sgn * f1 * uz * dxy; }
	  
	  if ( (ymin_s[ii][2] + Tiny) < yy && (yy + Tiny) < ymax_s[ii][2]  
	   && (fabs(xmin_s[ii][1] - xx) < Tiny || fabs(xx - xmax_s[ii][1]) < Tiny) )  {
	   sum += sgn * f2 * uz * dxy; }
	  
	  if ( (ymin_s[ii][2] + Tiny) < yy  && (yy + Tiny) < ymax_s[ii][2]  
	   && (fabs(xmin_s[ii][2] - xx) < Tiny || fabs(xx - xmax_s[ii][2]) < Tiny) )  {
	   sum += sgn * f3 * uz * dxy; }
          
	  if ( (fabs(ymin_s[ii][0] - yy) < Tiny || fabs(yy - ymax_s[ii][0]) < Tiny) 
	    && (fabs(xmin_s[ii][0] - xx) < Tiny || fabs(xx - xmax_s[ii][0]) < Tiny) )  {
	    sum += sgn * f11 * uz * dxy; }
	  
	  if ( (fabs(ymin_s[ii][1] - yy) < Tiny || fabs(yy - ymax_s[ii][1]) < Tiny) 
	    && (fabs(xmin_s[ii][0] - xx) < Tiny || fabs(xx - xmax_s[ii][0]) < Tiny) )  {
	    sum += sgn * f12 * uz * dxy; }
	  
	  if ( (fabs(ymin_s[ii][2] - yy) < Tiny || fabs(yy - ymax_s[ii][2]) < Tiny) 
	    && (fabs(xmin_s[ii][0] - xx) < Tiny || fabs(xx - xmax_s[ii][0]) < Tiny) )  {
	    sum += sgn * f13 * uz * dxy; }
	  
	  if ( (fabs(ymin_s[ii][0] - yy) < Tiny || fabs(yy - ymax_s[ii][0]) < Tiny) 
	    && (fabs(xmin_s[ii][1] - xx) < Tiny || fabs(xx - xmax_s[ii][1]) < Tiny) )  {
	    sum += sgn * f12 * uz * dxy; }
	  
	  if ( (fabs(ymin_s[ii][1] - yy) < Tiny || fabs(yy - ymax_s[ii][1]) < Tiny) 
	    && (fabs(xmin_s[ii][1] - xx) < Tiny || fabs(xx - xmax_s[ii][1]) < Tiny) )  {
	    sum += sgn * f22 * uz * dxy; }
	  
	  if ( (fabs(ymin_s[ii][2] - yy) < Tiny || fabs(yy - ymax_s[ii][2]) < Tiny) 
	    && (fabs(xmin_s[ii][1] - xx) < Tiny || fabs(xx - xmax_s[ii][1]) < Tiny) )  {
	    sum += sgn * f23 * uz * dxy; }
	  
	  if ( (fabs(ymin_s[ii][0] - yy) < Tiny || fabs(yy - ymax_s[ii][0]) < Tiny) 
	    && (fabs(xmin_s[ii][2] - xx) < Tiny || fabs(xx - xmax_s[ii][2]) < Tiny) )  {
	    sum += sgn * f13 * uz * dxy; }
	  
	  if ( (fabs(ymin_s[ii][1] - yy) < Tiny || fabs(yy - ymax_s[ii][1]) < Tiny) 
	    && (fabs(xmin_s[ii][2] - xx) < Tiny || fabs(xx - xmax_s[ii][2]) < Tiny) )  {
	    sum += sgn * f23 * uz * dxy; }
	  
	  if ( (fabs(ymin_s[ii][2] - yy) < Tiny || fabs(yy - ymax_s[ii][2]) < Tiny) 
	    && (fabs(xmin_s[ii][2] - xx) < Tiny || fabs(xx - xmax_s[ii][2]) < Tiny) )  {
	    sum += sgn * f33 * uz * dxy; }
	    
        } 

	}
		
      } endfor_ijk;
      
      local[0] += sum;
      local[1]  = sum;
      
    }
	
    local[0] *= symm_factor;	/* symmetry factor */
    local[1] *= symm_factor;
     
    /* reduce */
  
    MPI_Allreduce(local, global, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  
    return global[0];

}


