/* excision.c */
/* Bernd Bruegmann 1/03, Wolfgang Tichy 12/2003 */

/* Excision boundary condition
   These routines are based in part on the "Simple Excision" idea of 
     http://xxx.lanl.gov/abs/gr-qc/0008067
   The key idea is to "copy the time derivative" of variables to obtain
   boundary values as opposed to copying or extrapolating the variabes 
   themselves.

   Originally, a cubical excision region was used, but here we 
   generalize to (quite) arbitrary excision masks.
   For arbitrary surfaces one has to define a direction for the
   copy/extrapolation. There are several methods, some of which
   use the center of the black holes (Laguna/Maya et al, Yo et al). 
   Algorithms without explicit reference to a center include 
   Miguel Alcubierre's LegoExcision in Cactus, developed with
   some input from Deirdre Shoemaker and Eric Schnetter. Schnetter's 
   more complicated method to define a normal does not appear to be in 
   use in Cactus. Another early method without reference to a center
   is used in BB 1996, http://xxx.lanl.gov/abs/gr-qc/9608050.

   Initially we plan to implement LegoExcision and BB 1996 + copy of 
   time derivative. 
*/

#include "bam.h"
#include "boundary.h"

#define EXC_SPHERE 0
#define EXC_CUBE   1

#define EXC_INTERIOR 0
#define EXC_EXTERIOR 1
#define EXC_BOUNDARY 2




/* find direction that is closest to a certain normal 
   here: average over available directions to estimate normal
   BB: why not average data along all directions? or some directions?
   fills variable with index to point in best direction
*/ 
void find_mask_normal_BB(tL *level)
{
  double *m = Ptr(level, "excisionmask");
  double *dir = Ptr(level, "excisionindex");
  int ia[3][3][3];
  int i, j, k;
  double vx, vy, vz, s, max, dotproduct;
  int normalindex;

  /* estimate normal by averaging direction to all unmasked neighbors */
  forall27flag(level, EXCBOUND) {
    if (m[ccc] == EXC_BOUNDARY) {  // should be redundant
      fillindexarray27(ia);

      vx = vy = vz = 0;

      /* vector in normal direction as in MA's LegoExcision */
      if (1) {
        for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
        for (k = 0; k < 3; k++) {
	  if (m[ia[i][j][k]] == EXC_EXTERIOR) {
 	    vx += i-1;
	    vy += j-1;
	    vz += k-1;
	  }
        }
      }

      /* BB: doesn't the above overemphasize diagonal directions a little bit?
	 let's normalize each direction before taking the average
	 note that the center point is automatically excluded, hence s != 0
      */
      if (0) {
        for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
        for (k = 0; k < 3; k++) {
	  if (m[ia[i][j][k]] == EXC_EXTERIOR) {
	    s = sqrt((double)((i-1)*(i-1)+(j-1)*(j-1)+(k-1)*(k-1)));
 	    vx += (i-1)/s;
	    vy += (j-1)/s;
	    vz += (k-1)/s;
	  }
        }
      }

      /* find a grid direction that comes closest to the normal direction */
      /* this will be asymmetric walking around sphere, and it will
         dependent on round-off if several directions give the same 
         dot product
      */
      if (1) {
	normalindex = ccc;
	max = 0;
        for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
        for (k = 0; k < 3; k++) {
	  if (m[ia[i][j][k]] == EXC_EXTERIOR) {
	    s = sqrt((double)((i-1)*(i-1)+(j-1)*(j-1)+(k-1)*(k-1)));
	    dotproduct = (vx*(i-1) + vy*(j-1) + vz*(k-1))/s;

	    if (dotproduct > max) {
	      max = dotproduct;
	      normalindex = ia[i][j][k];
	    }
	  }
        }
      }

      /* sanity check */
      if (normalindex == ccc)
	errorexit("Excision did not find a good normal direction!\n");
      if (m[normalindex] != EXC_EXTERIOR)
	printf("Warning: Excision picked a bad direction!\n");

      if (0) {
	double *x = Ptr(level, "x");
	double *y = Ptr(level, "y");
	double *z = Ptr(level, "z");

	printf("%6.3f %6.3f %6.3f  %6.3f %6.3f %6.3f  %6.3f %6.3f %6.3f\n",
	       x[ccc], y[ccc], z[ccc], vx, vy, vz, 
	       x[normalindex], y[normalindex], z[normalindex]);
      }

      /* save index */
      dir[ccc] = normalindex;
    }
  } endfor;

  bampi_synchronize(level, Ind("excisionindex"));
}




/* WT: Fill in an array with the 6 Next Nearest neighbors, and ccc otherwise. 
   This is needed in the function below to easily access those 
   6  Next Nearest neighbors.
   Maybe this define could go into bam_amr_loops.h */
#define fillindexarray6NextNearest(indexarray) \
          indexarray[1][1][1] = ccc; \
          indexarray[0][1][1] = Mcc; \
	  indexarray[2][1][1] = Pcc; \
	  indexarray[1][0][1] = cMc; \
	  indexarray[1][2][1] = cPc; \
	  indexarray[1][1][0] = ccM; \
	  indexarray[1][1][2] = ccP; \
	  indexarray[0][0][1] = ccc; \
	  indexarray[0][1][0] = ccc; \
	  indexarray[1][0][0] = ccc; \
	  indexarray[2][2][1] = ccc; \
	  indexarray[2][1][2] = ccc; \
	  indexarray[1][2][2] = ccc; \
	  indexarray[0][2][1] = ccc; \
	  indexarray[0][1][2] = ccc; \
	  indexarray[1][0][2] = ccc; \
	  indexarray[2][0][1] = ccc; \
	  indexarray[2][1][0] = ccc; \
	  indexarray[1][2][0] = ccc; \
	  indexarray[0][0][0] = ccc; \
	  indexarray[0][0][2] = ccc; \
	  indexarray[0][2][0] = ccc; \
	  indexarray[0][2][2] = ccc; \
	  indexarray[2][0][0] = ccc; \
	  indexarray[2][0][2] = ccc; \
	  indexarray[2][2][0] = ccc; \
	  indexarray[2][2][2] = ccc
/* WT find index of point from which we want to copy the time derivative */
void find_mask_normal_WT(tL *level)
{
  double *m = Ptr(level, "excisionmask");
  double *dir = Ptr(level, "excisionindex");
  int ia[3][3][3];
  int iann[3][3][3];
  int i, j, k;
  double vx, vy, vz, s, max, dotproduct;
  int normalindex;

  /* estimate normal by averaging direction to all unmasked neighbors */
  forall27flag(level, EXCBOUND)
  {
    if(m[ccc] == EXC_BOUNDARY) // should be redundant 
    { 
      /* declare and set Mcc, Pcc,  cMc, cPc,  ccM, ccP */
      int_nbinds_faces2;
      set_nbinds_check2_6(nodes);

      /* fill in arrays ia and iann */
      fillindexarray27(ia);
      fillindexarray6NextNearest(iann);
         
      /* vector in normal direction */
      vx = vy = vz = 0;

      /* average v to find a vector in the normal direction, 
         here I include the next nearest neighbors 
         Mcc, Pcc,  cMc, cPc,  ccM, ccP                   */
      if(1) 
      {
        /* check nearest neigbors */
        for(i = 0; i < 3; i++)
        for(j = 0; j < 3; j++)
        for(k = 0; k < 3; k++) 
        {
	  if(m[ia[i][j][k]] == EXC_EXTERIOR) 
	  {
	    s = ((double)((i-1)*(i-1)+(j-1)*(j-1)+(k-1)*(k-1)));

 	    vx += (i-1)/s;
	    vy += (j-1)/s;
	    vz += (k-1)/s;
	  }
        }
        /* check the 6 next nearest neigbors */
        for(i = 0; i < 3; i++)
        for(j = 0; j < 3; j++)
        for(k = 0; k < 3; k++)
        {
          if( (m[iann[i][j][k]] == EXC_EXTERIOR) && (iann[i][j][k] != ccc) )
          {
            vx += (i-1);
            vy += (j-1);
            vz += (k-1);
          }
        }		
      }

      /* find a grid direction that comes closest to the normal direction */
      /* this will be asymmetric walking around sphere, and it will
         dependent on round-off if several directions give the same 
         dot product
      */
      if(1) 
      {
	normalindex = ccc;
	max = 0;
	
	/* check nearest neigbors */
        for(i = 0; i < 3; i++)
        for(j = 0; j < 3; j++)
        for(k = 0; k < 3; k++)
        {
	  if (m[ia[i][j][k]] == EXC_EXTERIOR)
	  {
	    s = sqrt((double)((i-1)*(i-1)+(j-1)*(j-1)+(k-1)*(k-1)));
	    dotproduct = (vx*(i-1) + vy*(j-1) + vz*(k-1))/s;

	    if (dotproduct > max)
	    {
	      max = dotproduct;
	      normalindex = ia[i][j][k];
	    }
	  }
        }
        /* check the 6 next nearest neigbors */
        for(i = 0; i < 3; i++)
        for(j = 0; j < 3; j++)
        for(k = 0; k < 3; k++)
        {
          if( (m[iann[i][j][k]] == EXC_EXTERIOR) && (iann[i][j][k] != ccc) )
          {
            s = ((double)((i-1)*(i-1)+(j-1)*(j-1)+(k-1)*(k-1)));
            /* here I put a 0.9 in front to de-emphasize the
               6 next nearest neigbors somewhat */
            dotproduct = 0.9*(vx*(i-1) + vy*(j-1) + vz*(k-1))/s;

            if (dotproduct > max)
	    {
	      max = dotproduct;
	      normalindex = iann[i][j][k];
	    }
          }
        }
      }

      /* sanity check */
      if (normalindex == ccc)
	errorexit("Excision did not find a good normal direction!\n");
      if (m[normalindex] != EXC_EXTERIOR)
	printf("Warning: Excision picked a bad direction!\n");

      if (0) {
	double *x = Ptr(level, "x");
	double *y = Ptr(level, "y");
	double *z = Ptr(level, "z");

	printf("%6.3f %6.3f %6.3f  %6.3f %6.3f %6.3f  %6.3f %6.3f %6.3f\n",
	       x[ccc], y[ccc], z[ccc], vx, vy, vz, 
	       x[normalindex], y[normalindex], z[normalindex]);
      }

      /* save index */
      dir[ccc] = normalindex;
    }
  } endfor;

  bampi_synchronize(level, Ind("excisionindex"));
}




/* wrapper */
void find_mask_normal(tL *level)
{
 if(Getv("excision_mask_normalfinder", "simple")) 
   find_mask_normal_BB(level);
 else
   find_mask_normal_WT(level);
}




/* find and mark points at the boundary of a mask 
   after calling this function the meaning of the flags is:
     1: point is not masked
     0: point is masked and has no unmasked neighbors
     2: point is masked and has unmasked neighbors 
*/
void find_mask_boundary(tL *level)
{
  double *m = Ptr(level, "excisionmask");
  int i, j, k;
  int ia[3][3][3];

  forinner27(level) {
    if (m[ccc] == EXC_INTERIOR) {
      /* reset mask at mask boundary */
      fillindexarray27(ia);
      for (i = 0; i < 3; i++)
      for (j = 0; j < 3; j++)
      for (k = 0; k < 3; k++)
	if (m[ia[i][j][k]] == EXC_EXTERIOR) {
	  m[ccc] = EXC_BOUNDARY;
	  i = j = k = 3;
	}

      /* set the boundary flag */
      if (boundaryflag(level, ccc) == 0) {
	if (m[ccc] == EXC_BOUNDARY)
	  boundaryflag(level, ccc) = EXCBOUND;
	else
	  boundaryflag(level, ccc) = EXCINTER;
      }
    }
  } endfor;

  bampi_synchronize(level, Ind("excisionmask"));
}




/* set flags in entire excision mask to exterior */
void set_entire_mask_to_EXTERIOR(tL *level)
{
  int i;
  double *mask = level->v[Ind("excisionmask")];

  /* set flags in mask */
  for (i = 0; i < level->nnodes; i++)
    mask[i] = EXC_EXTERIOR;
}




/* set flags in excision mask for a sphere or a cube */
void set_mask_inside_cubeorsphere(tL *level, int type, 
			          double r0, double x0, double y0, double z0)
{
  int i;
  double *x = level->v[Ind("x")];
  double *y = level->v[Ind("y")];
  double *z = level->v[Ind("z")];
  double *mask = level->v[Ind("excisionmask")]; 
  double dx, dy, dz, r2;

  if (r0 == 0) return;

  /* adjust r0 to avoid irregular shapes when r0 is multiple of grid spacing */
  r0 -= level->grid->level[0]->dx/1e9; 
  r2 = r0 * r0;

  /* set flags in mask */
  for (i = 0; i < level->nnodes; i++) {
    dx = x[i] - x0;
    dy = y[i] - y0;
    dz = z[i] - z0;

    if ((type == EXC_SPHERE && (dx*dx + dy*dy + dz*dz) < r2) ||
	(type == EXC_CUBE && fabs(dx) < r0 && fabs(dy) < r0 && fabs(dz) < r0))
      mask[i] = EXC_INTERIOR;
  }
}




/* set flags in excision mask in points with negative lapse */
void set_lapseneg(tL *level)
{
  int i;
  double *mask = level->v[Ind("excisionmask")]; 
  double *alpha = level->v[Ind("alpha")];
  double minLapse = Getd("excision_minLapse");
  
  /* set flags in mask */
  for (i = 0; i < level->nnodes; i++) {
    if (dless(alpha[i], minLapse)) {
      mask[i] = EXC_INTERIOR;
    }
    else {
      mask[i] = EXC_EXTERIOR;
    }
  }
}




/* set excision mask based on static black hole info */
int set_mask_bhdata(tL *level)
{
  int i, type;
  double r[2], x[2], y[2], z[2];
  double rah;
  double *mask; 
  double ds = level->dx;

  if (!Getv("boundary", "excision") ||
      !(Getv("physics", "KerrSchild") || Getv("physics", "RealisticBBH") ||
	Getv("physics", "punctures") ||  Getv("physics", "Schwarzschild"))) {
    return 0;
  }
  type = Getv("excision_shape", "sphere") ? EXC_SPHERE : EXC_CUBE;

  /* rah = radius of apparent horizon in M */
  if( Getv("physics", "punctures") )
    rah = 0.5;
  else if( Getv("physics", "KerrSchild") )
    rah = 2.0;
  else
    rah = 1.0;

  /* shrink the size to create a buffer zone */
  rah *= Getd("excision_shrink");

  /* for a cube, rah refers to the distance to the face, so shrink
     to move corners into sphere */
  if (type == EXC_CUBE) 
    rah /= sqrt(3.0);

  /* get black hole parameters */
  r[0] = rah * Getd("bhmass1");
  x[0] = Getd("bhx1");
  y[0] = Getd("bhy1");
  z[0] = Getd("bhz1");
  r[1] = rah * Getd("bhmass2");
  x[1] = Getd("bhx2");
  y[1] = Getd("bhy2");
  z[1] = Getd("bhz2");

  if (Getv("boundary", "excision") && Getv("grid", "1d")) {
    if (!dequal(r[0],0.0) && r[0] < 2*ds ) {
      errorexit("You have chosen cartoon and excision, but you need two points"
                " inside the excision region for this to work. "
		"Change your resolution\n");
    }
    if (!dequal(r[1],0.0) && r[1] < 2*ds ) {
      errorexit("You have chosen cartoon and excision, but you need two points"
                " inside the excision region for this to work. "
		"Change your resolution\n");
    }
  }

  /* set mask for each hole */
  for (i = 0; i < 2; i++)
    set_mask_inside_cubeorsphere(level, type, r[i], x[i], y[i], z[i]);

  /* mark boundary of mask */
  find_mask_boundary(level);

  /* find index of point in the direction closest to normal */
  find_mask_normal(level);

  /* print info */
  if (0) prvar01(level, "excisionmask");
  if (1 && r[0] > 0) {
    printf("Turning on Simple Lego Excision with %s:\n",
	   !type ? "spheres" : "cubes");
    for (i = 0; i < 2; i++)
      if (r[i] > 0) 
	printf("  r%d = %6.3f, center at %6.3f %6.3f %6.3f\n", 
	       i+1, r[i], x[i], y[i], z[i]);
  }
  return 0;
}




/* set mask manually */
int set_mask_manual(tL *level)
{
  int i, type;
  double r[2], x[2], y[2], z[2];
  double rah;
  double *mask; 
  double ds = level->dx;
  
  if (!Getv("boundary", "excision") ||
      !Getv("excision_set", "manual") )
    return 0;
  
  type = Getv("excision_shape", "sphere") ? EXC_SPHERE : EXC_CUBE;
  
  /* get black hole parameters */
  r[0] = Getd("excision_r0");
  x[0] = Getd("excision_x0");
  y[0] = Getd("excision_y0");
  z[0] = Getd("excision_z0"); 

  r[1] = Getd("excision_r1");
  x[1] = Getd("excision_x1");
  y[1] = Getd("excision_y1");
  z[1] = Getd("excision_z1"); 

  if (Getv("boundary", "excision") && Getv("grid", "1d")) {
    if (!dequal(r[0],0.0) && r[0] < 2*ds ) {
      errorexit("You have chosen cartoon and excision, but you need two points"
                " inside the excision region for this to work. "
		"Change your resolution\n");
    }
    if (!dequal(r[1],0.0) && r[1] < 2*ds ) {
      errorexit("You have chosen cartoon and excision, but you need two points"
                " inside the excision region for this to work. "
		"Change your resolution\n");
    }
  }

  /* set mask for each hole */
  for (i = 0; i < 2; i++)
    set_mask_inside_cubeorsphere(level, type, r[i], x[i], y[i], z[i]);
  
  /* mark boundary of mask */
  find_mask_boundary(level);
  
  /* find index of point in the direction closest to normal */
  find_mask_normal(level);

  /* print info */
  if (0) prvar01(level, "excisionmask");
  if (1 && r[0] > 0) {
    printf("Turning on Simple Lego Excision with %s:\n",
	   !type ? "spheres" : "cubes");
    for (i = 0; i < 2; i++)
      if (r[i] > 0) 
	printf("  r%d = %6.3f, center at %6.3f %6.3f %6.3f\n", 
	       i+1, r[i], x[i], y[i], z[i]);
  }
  return 0;
}




/* set excision mask based on static black hole info */
int set_mask_lapseneg(tL *level)
{
  int i, type;
  double *mask; 

  if (!Getv("boundary", "excision") ||!Getv("excision_set", "lapseneg") ||
      !(Getv("physics", "KerrSchild") || Getv("physics", "punctures") ||  
	Getv("physics", "Schwarzschild"))) {
    return 0;
  }

  set_lapseneg(level);

  /* mark boundary of mask */
  find_mask_boundary(level);

  /* find index of point in the direction closest to normal */
  find_mask_normal(level);

  /* print info */
  if (1) {
    printf("Turning on excision for regions with negative lapse");
  }
  return 0;
}




/* set excision boundary for one variable */
void set_boundary_excision(tL *level, int unew, int upre) 
{
  double *varpre = level->v[upre];
  double *varnew = level->v[unew];
  double *mask = Ptr(level, "excisionmask");
  double *index = Ptr(level, "excisionindex");
  double var0;
  int i, j;

  if(Getv("excision_InteriorValue", "VarFarLimit"))
    var0 = VarFarLimit(unew);
  else
    var0 = Getd("excision_InteriorValue");

  /* for all points, but could store list of excision points instead */
  forallpoints(level, i) {

    /* set points in boundary */
    if (mask[i] == EXC_BOUNDARY) {
      j = index[i];

      /* perhaps the key trick of simple excision: copy the time derivative
         i.e. evolve each point on the excision boundary by adding 
         the rate of change at some point along the outgoing normal
      */
      varnew[i] = varpre[i] + varnew[j] - varpre[j];
    }

    /* overwrite interior since that may have been evolved for simplicity */
    else if (1 && mask[i] == EXC_INTERIOR)
      varnew[i] = var0;
  }

  /* note: calling routine has to synchronize! */
}




/* copy from index point (ouside Excision boundary) 
   onto Excision boundary for varlist vlu               */
void copy_onto_excision_boundary(tL *level, tVarList *vlu)
{
  double *varnew;
  double *mask =  Ptr(level, "excisionmask");
  double *index = Ptr(level, "excisionindex");
  int i, j, ind;

  /* for all points, but could store list of excision points instead */
  forallpoints(level, i)
  {
    if(mask[i] == EXC_BOUNDARY)
      for(j = 0; j < vlu->n; j++)
      {
        varnew = VLPtr(vlu, j);
        /* copy var from outside at point index[i] */
        ind = index[i];
        varnew[i] = varnew[ind];
      }
    else if(1 && mask[i] == EXC_INTERIOR)
      for(j = 0; j < vlu->n; j++)
      {
        varnew = VLPtr(vlu, j);
        varnew[i] = 0.0;
      }
  }
  /* note: calling routine has to synchronize! */
}




/* set excision flags */
int set_boundary_flags_excision(tL *level) 
{
  /* why do you call me, anyway? */
  if (!Getv("boundary", "excision")) return 0;

  /* enable storage for the mask, which is initialized to zero */
  enablevar(level, Ind("excisionmask"));
  enablevar(level, Ind("excisionindex"));
  set_entire_mask_to_EXTERIOR(level);

  /* excise only on finest level? (to avoid problems with too coarse grids) */
  if( Getv("excision_OnlyOnFinestLevel", "yes") )
    if (level->l < level->grid->lmax)
      return 0;
    
  printf("Setting excision flags and mask on level %d:\n", level->l);
  
  /* pick method */
  if (Getv("boundary", "excision") && Getv("excision_set","bh")) {
    printf ("Setting bh mask\n");
    set_mask_bhdata(level);
  }
  else if 
    (Getv("boundary", "excision") && Getv("excision_set","manual")) {
    set_mask_manual(level);
  }
  else if 
    (Getv("boundary", "excision") && Getv("excision_set","lapseneg")) {
    set_mask_lapseneg(level);
  }
  return 0;
}




/* FIX */
/* function to register when we want to set excision flags */
int regrid_for_excision(tL *level) {
  regrid(level->grid, level->l,0);
  return 0;
}
