/* box.c */
/* Bernd Bruegmann, 2/2006 */

#include "bam.h"
#include "amr.h"

#define PR 0



/* allocate and initialize rectangular box
*/
tB *alloc_box(int m, int n, int o, double x0, double y0, double z0, 
	      double dx, double dy, double dz)
{
  tB *b;
  int var_x = Ind("x");
  int var_y = Ind("y");
  int var_z = Ind("z");
  int i, j, k, ijk, di, dj, dk;
  double dmin, dmax;

  /* print input */
  if (0) 
    printf("make_box_local: m=%d n=%d o=%d  %6.3f %6.3f %6.3f  h=%.3f\n",
	   m, n, o, x0, y0, z0, dx);

  /* allocate box */
  b = calloc(1, sizeof(tB));
  if (!b) errorexit("make_box_local: out of memory");

  /* set available parameters */
  b->bbox[0] = x0;
  b->bbox[1] = x0 + (m-1)*dx;
  b->bbox[2] = y0;
  b->bbox[3] = y0 + (n-1)*dy;
  b->bbox[4] = z0;
  b->bbox[5] = z0 + (o-1)*dz;

  b->x0 = x0;
  b->y0 = y0;
  b->z0 = z0;

  b->dd[0] = b->dx = dx;
  b->dd[1] = b->dy = dy;
  b->dd[2] = b->dz = dz;

  b->m = m;
  b->n = n;
  b->o = o;
  b->ibbox[1] = b->imax = m-1;
  b->ibbox[3] = b->jmax = n-1;
  b->ibbox[5] = b->kmax = o-1;

  b->npoints = m*n*o;

  b->id[0] = b->di = 1;
  b->id[1] = b->dj = m;
  b->id[2] = b->dk = m*n;

  return b;
}




/* add box into existing level */
void level_add_box(tL *level, tB *box)
{
  int pr = 0;

  box->i = level->nboxes;
  box->l = level->l;
  box->level = level;

  box->noffset = level->npoints;
  level->npoints += box->npoints;
  level->nnodes = level->npoints;

  level->nboxes++;
  level->box = realloc(level->box, level->nboxes * sizeof(tB *));
  level->box[level->nboxes - 1] = box;

  if (level->nboxes > 1) {
    if (pr) printf("adding box: joined bbox becomes\n");
    if (pr) printbbox(level, level->bbox, level->ibbox);
    box_join(level->bbox, box->bbox);
    level->ibbox[1] = (level->bbox[1] - level->bbox[0])/level->dx + 0.5;
    level->ibbox[3] = (level->bbox[3] - level->bbox[2])/level->dy + 0.5;
    level->ibbox[5] = (level->bbox[5] - level->bbox[4])/level->dz + 0.5;
    if (pr) printbbox(level, level->bbox, level->ibbox);
  }
}




/* add one box level into existing level, add data
   currently only boundary and x,y,z data are copied
*/
void level_add_box_withdata(tL *level, tL *level0)
{
  int nsofar, ntotal;
  int i, j;

  nsofar = level->npoints;
  level_add_box(level, level0->box[0]);

  ntotal = level->npoints;
  //level->boundary = realloc(level->boundary, ntotal*sizeof(int));
  int *tmp = (int *) realloc(level->boundary, ntotal*sizeof(int));
  if (!tmp) {
    tmp = (int *) malloc(ntotal*sizeof(int));
    for (i = 0; i < nsofar; i++)
      tmp[i] = level->boundary[i];
    free(level->boundary);
  }
  level->boundary = tmp;
  
  
  
  if (!level->boundary) printf("boundary of current level is not defined\n");
  if (!level0->boundary) printf("boundary of additional box is not defined\n");
  
  for (i = nsofar; i < ntotal; i++) 
    level->boundary[i] = level0->boundary[i-nsofar];
  
  if (Ind("x") != 0)
    errorexit("x, y, z should be first in list of variables");

  for (j = 0; j < level->nvariables; j++) {
    if (level->v[j]) {
      double *tmp = (double *) realloc(level->v[j], ntotal*sizeof(double));
      if (!tmp) {
        if (!level->v[j]) printf("realloc failed\n");
        tmp = (double *) malloc(ntotal*sizeof(double));
        for (i = 0; i < nsofar; i++)
          tmp[i] = level->v[j][i];
        free(level->v[j]);
      }
      level->v[j] = tmp;
      if (!level->v[j]) printf("realloc failed totally\n");
    }
    if (j < 6) {
      if (!level->v[j])
        printf(" x y z is not enabled during level creation?!?!  %d\n",j);
      for (i = nsofar; i < ntotal; i++)
        level->v[j][i] = level0->v[j][i-nsofar];
    }
  }
}




/* make top level that defines physical outer boundary 
   here: rectangular box
   called by each processor for its own local box 
*/
tG *make_grid_box_local(int l, int m, int n, int o, double x0, double y0, double z0, 
			double dx, double dy, double dz)
{
  tG *g;
  tL *level;
  tB *box;
  int var_x = Ind("x");
  int var_y = Ind("y");
  int var_z = Ind("z");
  int i, j, k, ijk, di, dj, dk;
  double dmin, dmax;

  /* print info */
  if (PR) printf("\n+++++++++++++ new box based amr ++++++++++++++\n\n");
  if (PR) 
    printf("make_grid_box_local: m=%d n=%d o=%d  %6.3f %6.3f %6.3f  h=%.3f\n",
	   m, n, o, x0, y0, z0, dx);

  /* make grid structure with lmin=0, lmax=0, nvariables */
  g = alloc_grid(0, 0, globalnvariables);
  if (0) printgrid(g);

  /* make level structure at level=0 with so many nodes and add it to g */
  level = alloc_level(g, 0, m*n*o);
  level->dx = dx;
  level->dy = dy;
  level->dz = dz;
  level->boundary = icalloc(level->nnodes);
  level->l  = l;

  /* add one box */
  box = alloc_box(m, n, o, x0, y0, z0, dx, dy, dz);
  level_add_box(level, box);

  /* time step */
  dmax = max3(dx, dy, dz);
  dmin = min3((m>1)?dx:dmax, (n>1)?dy:dmax, (o>1)?dz:dmax);
  level->dt = dmin * Getd("dtfac");
  if (PR) printlevel(level);

  /* enable storage of the coordinate variables */
  enablevar(level, var_x);
  enablevar(level, var_y);
  enablevar(level, var_z);
  
  /* fill in coordinates of points/nodes */
  forallpoints_ijk(level) {
    /* coordinates */
    level->v[var_x][ijk] = box->x0 + i*level->dx;
    level->v[var_y][ijk] = box->y0 + j*level->dy;
    level->v[var_z][ijk] = box->z0 + k*level->dz;
  } endfor_ijk;

  /* done! */
  return g;
}




/* find parent box pointers */
void level_setboxpr(tL *parent, tL *level)
{
  int i, j;

  for (j = 0; j < level->nboxes; j++) {
    for (i = 0; i < parent->nboxes; i++) {
      if (0) printbbox(level,   level->box[j]->com->bbox, 0);
      if (0) printbbox(parent, parent->box[i]->com->bbox, 0);
      if (box_ainb(level->box[j]->com->bbox, parent->box[i]->com->bbox)) {
	level->box[j]->pr = parent->box[i];
	break;
      }
    }
    if (!level->box[j]->pr) errorexit("no parent for child");
  }
}




/* find parent pointers
   given level l, we set parent box pointers on level l and l+1 if available
   unfinished elsewhere: don't do l+1
*/
void level_setboxprch(tL *level)
{
  int i, j;
  tL *parent = level->prlocal;
  tL *child = (level->l < level->grid->lmax) ? 
              level->grid->level[level->l + 1] : 0;

  if (!level)  errorexit("level_setprch: need level");
  if (!parent) errorexit("level_setprch: need local parent");

  if (0) {
    printf("level_setboxprch: l%d\n", level->l);
    printf(" level %p, nboxes %d\n", level,   level->nboxes);
    printf("parent %p, nboxes %d\n", parent, parent->nboxes);
    if (0) printf(" child %p, nboxes %d\n", child,  child ? child->nboxes : 0);
  }

  level_setboxpr(parent, level);
  if (0) if (child) level_setboxpr(level, child);
}




/* parent bounding box */
void box_parent_bbox(double *p, double *f, double dx, double dy, double dz,
		     int *bflag, int nbuffer)
{
} 




/* get index in parent box 
   fx = fx0 + fi*fdx;
   ci = (fx - cx0)/cdx;
   recall how points are staggered:
   this is an integer i-1/4 or i+1/4, snap to i by [ci+0.5]
*/
int box_parent_index(double cx0, double cdx, double fx0, double fdx, int fi)
{
  return (int)(0.5 + (fx0 + fi*fdx - cx0)/cdx);
}




/* get index in child box 
   as for parent, except that 2i+1/2 snaps to 2i
   as a curiosity note that with "+ 3/8" these two calculations are identical
*/
int box_child_index(double fx0, double fdx, double cx0, double cdx, int ci)
{
  return (int)((cx0 + ci*cdx - fx0)/fdx);
}



/* get index for different coordinates, snap down */
int box_point_index(double X0, double dX, double x0, double dx, int i)
{
  return (int)((x0 + i*dx - X0)/dX);
} 




/* get indices for different coordinates, snap down */
void box_point_indices(tB *c, tB *f, int *pci, int *pcj, int *pck,
			int fi, int fj, int fk)
{
  *pci = box_point_index(c->x0, c->dx, f->x0, f->dx, fi);
  *pcj = box_point_index(c->y0, c->dy, f->y0, f->dy, fj);
  *pck = box_point_index(c->z0, c->dz, f->z0, f->dz, fk);
}




/* get indices in parent box */
void box_parent_indices(tB *c, tB *f,
			int *pci, int *pcj, int *pck,
			int fi, int fj, int fk)
{
  *pci = box_parent_index(c->x0, c->dx, f->x0, f->dx, fi);
  *pcj = box_parent_index(c->y0, c->dy, f->y0, f->dy, fj);
  *pck = box_parent_index(c->z0, c->dz, f->z0, f->dz, fk);
}




/* get linear index in parent box */
int box_parent_ijk(tB *c, tB *f, int fi, int fj, int fk)
{
  return ijkofbox(c,
		  box_parent_index(c->x0, c->dx, f->x0, f->dx, fi),
		  box_parent_index(c->y0, c->dy, f->y0, f->dy, fj),
		  box_parent_index(c->z0, c->dz, f->z0, f->dz, fk));
}




/* get linear index in child box */
int box_child_ijk(tB *f, tB *c, int ci, int cj, int ck)
{
  return ijkofbox(f,
		  box_child_index(f->x0, f->dx, c->x0, c->dx, ci),
		  box_child_index(f->y0, f->dy, c->y0, c->dy, cj),
		  box_child_index(f->z0, f->dz, c->z0, c->dz, ck));
}




/* get box containing a given point */
tB *box_containing_xyz(tL *level, double x, double y, double z) 
{
  if (xyzinsidebbox(level->box[0]->bbox, x, y, z)) 
    return level->box[0];
  if (xyzinsidebbox(level->box[1]->bbox, x, y, z)) 
    return level->box[1];

  printf("%6.3f %6.3f %6.3f\n", x, y, z);
  printbbox(level, level->box[0]->bbox, 0);
  printbbox(level, level->box[1]->bbox, 0);
  errorexit("box containing xyz not found");
  return 0;
}




/* get box containing a point given by indices in another box */
tB *box_containing_ijk(tL *level, tB *b, int i, int j, int k)
{
  return box_containing_xyz(level, xofbox(b, i), yofbox(b, j), zofbox(b, k)); 
}




/* find local bbox minus parallelization and symmetry ghosts */
void box_findbbox_notghost(tB *box, double *bbox, int *ibbox)
{
  int nghosts = Geti("bampi_nghosts");
  int b, i;

  for (i = 0; i < 6; i++) {
    bbox[i] = box->com->bbox[i];
    ibbox[i] = box->com->ibbox[i];
  }

  i = 0;
  b = box->bflag[i];
  if (b == GHOBOUND || b == SYMBOUND) {
    bbox[i] += nghosts * box->dx;
    ibbox[i+1] -= nghosts;
  }
  b = box->bflag[++i];
  if (b == GHOBOUND || b == SYMBOUND) {
    bbox[i] -= nghosts * box->dx;
    ibbox[i] -= nghosts;
  }
  b = box->bflag[++i];
  if (b == GHOBOUND || b == SYMBOUND) {
    bbox[i] += nghosts * box->dy;
    ibbox[i+1] -= nghosts;
  }
  b = box->bflag[++i];
  if (b == GHOBOUND || b == SYMBOUND) {
    bbox[i] -= nghosts * box->dy;
    ibbox[i] -= nghosts;
  }
  b = box->bflag[++i];
  if (b == GHOBOUND || b == SYMBOUND) {
    bbox[i] += nghosts * box->dz;
    ibbox[i+1] -= nghosts;
  }
  b = box->bflag[++i];
  if (b == GHOBOUND || b == SYMBOUND) {
    bbox[i] -= nghosts * box->dz;
    ibbox[i] -= nghosts;
  }
}





/* make level based on global nesting function */
tL *make_level_nested_boxes(tL *level, int pr)
{
  int nboxes;
  tSBox *sbox = bcalloc(10 * sizeof(tSBox));
  tL *newlevel;

  /* determine bounding boxes based on global nesting function */
  find_nested_bboxes(level, &nboxes, sbox);

  /* make new level */
  newlevel = make_level_bboxes(level, nboxes, sbox, pr);

  /* save bbox information on coarse level */
  if (1) {
    free(level->sbox);
    level->sbox = sbox;
    level->nsboxes = nboxes;
  } else // for debugging 
    free(sbox);

  /* return new level */
  return newlevel;
}




/* make level based on given boxes */
tL *make_level_bboxes(tL *level, int nboxes, tSBox *box, int pr)
{
  tG *g = level->grid;
  int l = level->l;
  tL *newlevel, *newlevel0, *prlevel, *prlevel0;
  int i, j, nbox;
  int m, n, o;
  double dx, dy, dz, xmin, ymin, zmin;
  double *flagregrid;
  double *bbox;
  int *ibbox;

  /* for all boxes */
  for (nbox = 0; nbox < nboxes; nbox++) {
    double *bbox = box[nbox].bbox;
    int *ibbox = box[nbox].ibbox;

    /* turn parameters into those of the refined level
       flags correctly implement symmetries, see flagregridboundary()
    */ 
    m = 2 * (ibbox[1] + 1);
    n = 2 * (ibbox[3] + 1);
    o = 2 * (ibbox[5] + 1);
    dx = level->dx/2;
    dy = level->dy/2;
    dz = level->dz/2;
    xmin = bbox[0] - dx/2;
    ymin = bbox[2] - dy/2;
    zmin = bbox[4] - dz/2;

    /* create new level consisting of the box we want */
    newlevel0 = 
      make_level_box_coarsealigned(l+1, m, n, o, xmin, ymin, zmin, dx, dy, dz, pr);
    newlevel0->l = level->l+1;
    newlevel0->grid = level->grid;

    if (!newlevel0->boundary) 
      printf("boundary does not exist");
    //printf("after make (%d):   %p\n",nbox,newlevel0->boundary);
    if (!newlevel0->v[0])
      printf("  x y z is not enabled during level creation?!?!\n");

    
    /* save global bounding box (the local bbox is already stored in com) */
    for (i = 0; i < 6; i++) {
      newlevel0->box[0]->bbox[i] = newlevel0->bbox[i];
      newlevel0->box[0]->ibbox[i] = newlevel0->ibbox[i];
    }

    /* if this is the first box, insert new level into grid */
    if (nbox == 0) {
      newlevel = newlevel0;
      replace_level(g, newlevel, l+1);
      //printf("after replace :   %p %p   %p %p %p\n",newlevel->boundary,newlevel0->boundary, newlevel,newlevel0,g->level[l+1]);
      if (!g->level[l+1]->v[0])
        printf("  x y z is not enabled during level creation?!?!\n");
    }

    /* if this isn't the first box, add box to existing level, discard level */
    else {
      //printf("before add :   %p %p\n",newlevel->boundary,newlevel0->boundary);
      if (!newlevel->v[0])
        printf("  x y z is not enabled during level creation?!?!\n");
      level_add_box_withdata(newlevel, newlevel0);
      newlevel0->box[0] = 0;
      newlevel0->com = 0;
      free_level(newlevel0);
    }
  }

  /* the new level now has list of boxes stored */
  if (pr) {
    printf("local newlevel:\n");
    for (i = 0; i < newlevel->nboxes; i++)
      printbbox(newlevel,
		newlevel->box[i]->com->bbox, newlevel->box[i]->com->ibbox);
  }

  /* points flagged as physical boundary become the refinement boundary */
  forallpoints(newlevel, i) {
    if (boundaryflag(newlevel, i) == PHYBOUND) 
      boundaryflag(newlevel, i) = REFBOUND;
  }
  forallboxes(newlevel) {
    for (i = 0; i < 6; i++)
      if (box->bflag[i] == PHYBOUND) 
	box->bflag[i] = REFBOUND; 
  } endforboxes;

  for (i = 0; i < g->nvariables; i++)
    if (level->v[i])
      enablevarcomp(newlevel, i);
  
  /* sort out parent/child relationships 
     there is some reverse engineering because at this point we work
     with global nesting function which does not know parents
     currently there can be 1 or 2 parents with 1 or 2 children
  */

  /* if we are on a single processor with single box, use existing parent */
  if (bampi_size() == 1 && newlevel->nboxes == 1) {
    prlevel = newlevel->prlocal = level;
  }

  /* if there is more than one processor or more than one box */
  if (bampi_size() > 1 || newlevel->nboxes > 1) {
    int nnn;

    forallboxes(newlevel) {
      double bbox[6];
      int ibbox[6];

      /* find local bounding box of new level minus ghost points */
      box_findbbox_notghost(box, bbox, ibbox);
      if (pr) printf("local newlevel minus ghosts:\n");
      if (pr) printbbox(newlevel, bbox, ibbox);
    
      /* determine grid parameters including buffer zone */
      /* preliminary: 
	 - assume that we need nghosts both for P and symmetries
	 which should be safe but we should get away with fewer points
	 - alternatively, try order_RP/2 points which should be exact
	 for non-symmetry boundaries, but some symmetry boundaries
	 might actually insist on nghost points
      */
      if (1) nnn = Geti("bampi_nghosts");
      else   nnn = Geti("order_RP") / 2;   // currently fails
      dx = level->dx;
      dy = level->dy;
      dz = level->dz;
      m = (ibbox[1]+1)/2 + 2*nnn;
      n = (ibbox[3]+1)/2 + 2*nnn;
      o = (ibbox[5]+1)/2 + 2*nnn;
      xmin = bbox[0] - (nnn - 0.25)*dx;
      ymin = bbox[2] - (nnn - 0.25)*dy;
      zmin = bbox[4] - (nnn - 0.25)*dz;
    
      /* create local parent level */
      prlevel0 = 
	make_level_box_local(l, m, n, o, xmin, ymin, zmin, dx, dy, dz);
      prlevel0->l = level->l;
      prlevel0->grid = level->grid;

      /* set bounding box for local parent level */
      findbbox(prlevel0, prlevel0->bbox, prlevel0->ibbox);
      if (pr) printf("local parent plus buffer:\n");
      if (pr) printbbox(prlevel0, prlevel0->bbox, prlevel0->ibbox);
      
      bampi_allreduce_bbox(prlevel0, prlevel0->bbox, prlevel0->ibbox);
      if (pr) printf("global local parent:\n");
      if (pr) printbbox(prlevel0, prlevel0->bbox, prlevel0->ibbox);
      for (i = 0; i < 6; i++) {
	prlevel0->box[0]->bbox[i] = prlevel0->bbox[i];
	prlevel0->box[0]->ibbox[i] = prlevel0->ibbox[i];
      }

      /* call bampi to fill in communication structure of local parent level */
      bampi_init_com(prlevel0, m, n, o, 
		     level->com->sizexyz[0], 
		     level->com->sizexyz[1], 
		     level->com->sizexyz[2]);

      /* local parent needs same variables enabled as global parent */
      for (i = 0; i < g->nvariables; i++)
	if (level->v[i])
	  enablevarcomp(prlevel0, i);
      
      /* set boundary flags */
      set_boundary_flags(prlevel0);
      

      /* if this is the first box */
      if (nbox == 0) {
	prlevel = prlevel0;
      }

      /* if this isn't the first box */
      else {
        level_add_box_withdata(prlevel, prlevel0);
	prlevel0->box[0] = 0;
	prlevel0->com = 0;
	free_level(prlevel0);
      }

    } endforboxes;

    /* set parent pointer */
    newlevel->prlocal = prlevel;
  }

  /* sort out box parent pointers */
  level_setboxprch(newlevel);

  /* set flags for restriction and prolongation */
  setRPflags(prlevel, newlevel);

  /* initialize flagregrid to make this the bottom level */
  flagregrid = PtrEnable(newlevel, "flagregrid");
  forallpoints(newlevel, i) flagregrid[i] = 0;

  /* initialize parent synchronization */
  bampi_syncparent_init(newlevel);

  /* sanity check */
  bampi_check_split(newlevel);
  bampi_check_split(newlevel->prlocal);

  /* return new level */
  return newlevel;
}




/* given a level, return pointer to level containing just one box
   helps using routines that expect one box only
   main customer: src/utility/output
   one copy only!
   don't ever free this level!
*/
tL *one_box_level(tL *level, int nbox)
{
  int i;
  static tL *obl = 0;
  static tB *obb = 0;
  static tB **pobb = 0;
  static double **pobv = 0;

  /* one time initialization of static storage */
  /* here we have to hide theses pointers to not show up inside the memory tracker 
     mth: this function is a sort of lazy function ... but it is needed */
  if (!obl) {
    obl = bcalloc(sizeof(tL));
    hidepointer(obl);
  }
  if (!obb) {
    obb = bcalloc(sizeof(tB));
    hidepointer(obb);
  }
  if (!pobb) {
    pobb = bcalloc(3*sizeof(tB *));
    hidepointer(pobb);
  }
  if (!pobv) {
    pobv = bcalloc(100*level->nvariables*sizeof(double *)); // fix me
    hidepointer(pobv);
  }
  
  /* copy entire level data structure, overwrites all previous information
     how robust is this? 
  */
  memcpy(obl, level, sizeof(tL));

  /* reset data if there is more than one box */
  if (level->nboxes > 1) {
    obl->box = pobb;
    obl->box[0] = obb;
    obl->v = pobv;

    /* copy box data structure */
    memcpy(obb, level->box[nbox], sizeof(tB));

    /* adjust data pointers to point to beginning of data for this box */
    obl->boundary = level->boundary + obb->noffset;
    for (i = 0; i < level->nvariables; i++) {
      if (level->v[i]) 
        obl->v[i] = level->v[i] + obb->noffset;    
      else 
        obl->v[i] = NULL;
    }

    /* set parameters for single box */
    obb->level = obl;
    obb->i = 0;
    obb->noffset = 0;
    obl->nboxes = 1;
    obl->npoints = obl->nnodes = obb->npoints;
    for (i = 0; i < 6; i++) obl->bbox[i] = obb->bbox[i];
    for (i = 0; i < 6; i++) obl->ibbox[i] = obb->ibbox[i];
    obl->com = obb->com;
  }

  /* return pointer to one box level */
  return obl;
}






/* function for boundaryaway(N) 
   gives the number of points to the next boundary back
*/
int boundaryaway_(int N, int i,int j,int k, int imax,int jmax,int kmax)
{
    int n;
    for (n = 0; n < N; n++) {
        if (boundaryNaway(n)) 
            return n;
    }
    return -1;
}

int test_box_in_box(tB *bf, tB *bc, int pr)
{
    int v = 0;
    int i,n;
    double bboxc[6];
    int ibboxc[6];
    
    if (bc->level->l!=bf->level->l-1) errorexit("watch your boxes");
    if (pr) printf("test if box %d(%d) in level %d   is inside   box %d(%d) in level %d\n",bf->i,bf->level->nboxes,bf->l,bc->i,bc->level->nboxes,bc->l);
    
    
    // test full box
    v += box_ainb(bf->bbox,bc->bbox);
    if (pr) {
        printf("  big box:     ");
        printbbox(bc->level,bc->bbox,bc->ibbox);
        printf("  small box:   ");
        printbbox(bf->level,bf->bbox,bf->ibbox);
        printf("  -> %d   (full box)\n",v);
    }
    
    
    // test local boxes
    /* I hope I have to use RestrictProlongOrder, It seem to be 
       the right one because you only have to care about interpolation
       not abount bampi-buffer or amr-buffer(which has to bigger than order_RP)
    */
    int gh = Geti("order_RP");
    for (i=0; i<6; i++) {
        // set boundary-cases -> I hope it is correct 
        // to substract points only at the coarser box
        // and then compare
        if      (bc->bflag[i] == SYMBOUND) n = 0;
        else if (bc->bflag[i] == GHOBOUND) n = 0;
        else if (bc->bflag[i] == PHYBOUND) n = ((i%2==0)?(1):(-1)); //physical boundary has only one point
        else                               n = ((i%2==0)?(1):(-1)) * gh+1;
        
        bboxc[i]  = bc->bbox[i]  + n * bc->dd[i/2];
        ibboxc[i] = bc->ibbox[i] + n;
    }
    v += box_ainb(bf->bbox,bboxc);
    if (pr) {
        printf("  big box:     ");
        printbbox(bc->level,bboxc,ibboxc);
        printf("  small box:   ");
        printbbox(bf->level,bf->bbox,bf->ibbox);
        printf("  -> %d   (full box without ghosts)\n",v);
        if (v!=2) printf("  ==> THIS BOX DOES NOT FIT\n");
        printf("\n");
    }
    
    
    return (v==2)?1:0;
}

int test_box_touch(tB *b1, tB *b2, int pr)
{
    if (b1->level->l!=b2->level->l) errorexit("watch your boxes");
    
    return box_overlap(b1->bbox,b2->bbox);
}

int test_box_consistency(tG* g, int level)
{

  if (g->level[level]->shells) return 1;

    int pr = 0;
    int l,b1,b2,sum;
    int problem = 0;
    if (pr) printf("test_consistency starting at level %d\n", level);
  
    // go through all finer levels and test if child is inside parent
    for (l=level+1; l<=g->lmax; l++) {
        for (b1=0; b1<g->level[l]->nboxes; b1++) {
            sum = 0;
            for (b2=0; b2<g->level[l-1]->nboxes; b2++) {
                sum += test_box_in_box( g->level[l]->box[b1], g->level[l-1]->box[b2],pr);
            }
            if (sum!=1) problem++;
        }
    }
    if (problem) errorexit("AMR-boxsetup fails\n  At least one box is not inside the parent one, mostly l1 is bigger l0,\n  since l0 does not have amr_buffers.\n  It can help to increase nxyz OR decrease amr_move_nxyz or you can decrease amr_nbuffer.\n  Set pr=1 inside test_box_consistency() to see a detailed box size printmessage\n");
    
    // go through all finer levels and test if boxes touch
    for (l=level; l<=g->lmax; l++) {
        for (b1=0; b1<g->level[l]->nboxes; b1++) {
            for (b2=b1+1; b2<g->level[l]->nboxes; b2++) {
                problem += test_box_touch( g->level[l]->box[b1], g->level[l]->box[b2],pr);
            }
        }
    }
    if (problem) errorexit("at leat one box touches/intersects one other box");

    return problem;
}

int test_box_consistency_level(tL* lc, tL* lf)
{
    int pr = 0;
    int b1,b2,sum;
    int problem = 0;
    if (pr) printf("test_consistency at level %p\n", lf);
    
     // go through finer level and test if child is inside parent
    for (b1=0; b1<lf->nboxes; b1++) {
        sum = 0;
        for (b2=0; b2<lc->nboxes; b2++) {
            sum += test_box_in_box( lf->box[b1], lc->box[b2],pr);
        }
        if (sum!=1) problem++;
    }
    if (problem) errorexit("\n  at least one box is not inside the parent one\n  it can help to increase nxyz OR decrease amr_move_nxyz\n  if not set pr=1 inside test_box_consistency()\n  to see a detailed box size printmessage\n");
    
    // go through finer leves and test if boxes touch
    for (b1=0; b1<lf->nboxes; b1++) {
        for (b2=b1+1; b2<lf->nboxes; b2++) {
            problem += test_box_touch(lf->box[b1], lf->box[b2],pr);
        }
    }
    if (problem) errorexit("at leat one box touches/intersects one other box");
    
    return problem;
}






