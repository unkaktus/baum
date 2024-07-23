/* move.c */
/* Bernd Bruegmann 12/2005 */

#include "bam.h"
#include "amr.h"

#define PR 0


/* still experimental and mostly unoptimized code, 12/2005
   some cases work, but not completely tested, 3/2006
   added rotant, now works in all known cases, 2/2007
*/






/* check whether two boxes are equal */
int box_equal(double *a, double *b)
{
  int i; 

  for (i = 0; i < 6; i++)
    if (!dequal(a[i], b[i])) return 0;
  return 1;
}




/* check whether two intervals overlap */
int interval_overlap(double *a, double *b)
{
  if (dless(a[1], b[0]) || dless(b[1], a[0]))
    return 0;
  return 1;
}

/* check whether two boxes overlap */
int box_overlap(double *a, double *b)
{
  if (interval_overlap(a+0, b+0) &&
      interval_overlap(a+2, b+2) &&
      interval_overlap(a+4, b+4))
    return 1;
  return 0;
}




/* check whether box a is contained in box b */
int box_ainb(double *a, double *b)
{
  if (xyzinsidebbox(b, a[0], a[2], a[4]) &&
      xyzinsidebbox(b, a[1], a[3], a[5]))
    return 1;
  return 0;
}




/* find separation between two intervals, defined by shortest distance */
int interval_separation(double h, double *a, double *b)
{
  double o = a[0];
  double round(double x);
  int a0 = round((a[0]-o)/h);
  int a1 = round((a[1]-o)/h);
  int b0 = round((b[0]-o)/h);
  int b1 = round((b[1]-o)/h);

  if (0) printf("%f  %f %f  %f %f  %d %d  %d %d\n", 
		h, a[0], a[1], b[0], b[1], a0, a1, b0, b1);

  if (b0 >= a1) return b0 - a1;
  if (a0 >= b1) return a0 - b1;
  return 0;
}

/* find separation between two boxes, defined by shortest distance */
int box_separation(double *h, double *a, double *b)
{
  int i = interval_separation(h[0], a+0, b+0);
  int j = interval_separation(h[1], a+2, b+2);
  int k = interval_separation(h[2], a+4, b+4);

  if (0) printf("i %d  j %d  k %d\n", i, j, k);

  if (j > i) i = j;
  if (k > i) i = k;
  return i;
}




/* given two boxes, return joined bounding box in first argument */
void box_join(double *a, double *b)
{
  if (a[0] > b[0]) a[0] = b[0];
  if (a[1] < b[1]) a[1] = b[1];
  if (a[2] > b[2]) a[2] = b[2];
  if (a[3] < b[3]) a[3] = b[3];
  if (a[4] > b[4]) a[4] = b[4];
  if (a[5] < b[5]) a[5] = b[5];
}




/* arrange nested boxes around two centers (two punctures) 
   overlapping boxes generate one big box
   a one point buffer zone around each box on coarse guarantees
   that fine boxes use coarse points for interpolation that are
   not owned by another fine box (may not be strictly necessary)
*/
void find_nested_bboxes(tL *level, int *pnboxes, tSBox *box)
{
  int nghosts = Geti("bampi_nghosts");
  int np = level->grid->npunc;
  double p[10][3];
  int n[3];
  double h[3];
  double bbox[10][6];
  double minymax;
  int flag[2];
  int nboxes;
  int i, j, k, l;
  int overlap, separation;
  int pr = 0;

  
  /* only move boxes with l > lcube, use fmr cubes for l <= lcube 
  note that we work with the coarse level and determine child boxes 
  we enforce one centered cube by putting all puncture points on the origin
  */
  for (j = 0; j < np; j++)
  for (i = 0; i < 3; i++) {
    if (level->l+1 <= Geti("amr_move_lcube"))
      p[j][i] = 0.;
    else 
      p[j][i] = level->grid->puncpos[j][i];
  } 
  
  
  if (pr) {
    for (j = 0; j < np; j++) 
      printf("point %d:  %9.6f %9.6f %9.6f\n", j+1, p[j][0], p[j][1], p[j][2]);
  }
  
  
  /* for testing: force motion */
  if (np<=2) {
    i = Geti("amr_move_force");
    if (i == 0) {
      p[0][0] -= level->time;
    }
    else if (i == 1) {
      /* single puncture, move box on circle, to be used with move.par */
      double r = 1.;
      double o = Getd("amr_move_force_o");
      p[0][0] = -r * 0.9*sin(o * level->time);
      p[0][1] =  r * 0.8*cos(o * level->time);
      printf("setting x %f, y %f\n", p[0][0], p[0][1]);
    }  
    else if (i == 2) {
      /* two punctures, move boxes on straight lines */
      p[0][1] -= level->time;
      p[1][1] += level->time;
      printf("setting x %f, y %f\n", p[0][0], p[0][1]);
    }
    else if (i == 3) {
      /* two punctures, move two boxes on circle */
      double r = Getd("amr_move_force_r"); 
      double o = Getd("amr_move_force_o");
      p[0][0] =  r * 0.9*sin(PI+o * level->time);
      p[0][1] =  r *     cos(PI+o * level->time);
      p[1][0] =  r * 0.9*sin(o * level->time);
      p[1][1] =  r *     cos(o * level->time);
      printf("setting p1:  x %f, y %f\n", p[0][0], p[0][1]);
      printf("setting p2:  x %f, y %f\n", p[1][0], p[1][1]);
    }
    else if (i == 4) {
      /* single puncture moves */
      double o = Getd("amr_move_force_o");
      p[0][0] = o*level->time;
      p[0][1] = 0;
      printf("setting x %f, y %f\n", p[0][0], p[0][1]);
    }
  }


  /* get second point by symmetry if we are running with symmetry */
  /* preliminary: we don't cover all possible arrangements here! */
  if (np<=2) {
    if (Getv("grid", "quadrant") || Getv("grid", "rotant")) {
      p[1][0] = - p[0][0];
      p[1][1] = - p[0][1];
      p[1][2] =   p[0][2];
    }
    else if (Getv("grid", "qreflect")) {
      p[1][0] = p[0][0];
      p[1][1] = p[0][1];
      p[1][2] = p[0][2];
    }
    else if (Getv("grid", "octant")) {
      p[1][0] = - p[0][0];
      p[1][1] = - p[0][1];
      p[1][2] = - p[0][2];
    }
    /* info */
    if (pr) {
      if (1) printf("find_nested_bboxes: level %d\n", level->l);
      if (1 || level->l == 0)
        for (j = 0; j < np; j++) 
          printf("puncture %d:  %9.6f %9.6f %9.6f\n",
                j, p[j][0], p[j][1], p[j][2]);
    }
  }

  if (Geti("amr_move_force")>=0) {
    for (j = 0; j < np; j++)
      for (i = 0; i < 3; i++) 
        level->grid->puncpos[j][i] = p[j][i];
  }
  
  
  

  /*** determine box size ***/

  /* size of basic top/coarsest level box */
  n[0] = Geti("nx"); 
  n[1] = Geti("nx"); 
/* HERE IS A BUG,   n[1] = Geti("ny"); gives error for more than one level????*/ 
  n[2] = Geti("nz");
  if (n[0]%2 || n[1]%2 || n[2]%2)
    errorexit("nested boxes need even nx, ny, nz");

  /* adjust basic size of moving box */
  if (level->l+1 > Geti("amr_move_lcube")) {
    j = Geti("amr_move_nxyz");
    if (j) n[0] = n[1] = n[2] = j;
  }

  /* if available, use list of box sizes for each level */
  i = GetdArrayN("amr_nxyz");
  if (i) {
    j = (i-1 > level->l+1) ? level->l+1 : i-1;
    n[0] = n[1] = n[2] = GetdEntry("amr_nxyz", j);
    if (n[0]%4 || n[1]%4 || n[2]%4)
      errorexit("box sizes in amr_nxyz have to be multiples of 4");
  }

  /* for large buffer zones it is clearer to add buffer as extra points */
  if (Getv("amr_buffer_addextra", "yes"))
    for (i = 0; i < 3; i++) n[i] += 2 * Geti("amr_nbuffer");

  /* basic factor two nesting, rounding up to multiple of two
     (which happens to be what flagregrid does)
  */
  for (i = 0; i < 3; i++) 
    n[i] = n[i]/2 + (n[i]/2)%2;


  /* grid spacing */
  h[0] = level->dx;
  h[1] = level->dy;
  h[2] = level->dz;

  /* determine boxes around each point */
  for (j = 0; j < np; j++) 
  for (i = 0; i < 3; i++) {

    /* snap to cell boundary
       note that (int) is symmetric wrt origin, which is essential for syms */
    p[j][i] = ((int) (p[j][i]/h[i])) * h[i];
    bbox[j][2*i+0] = p[j][i] - (n[i]-1)*h[i]/2;
    bbox[j][2*i+1] = p[j][i] + (n[i]-1)*h[i]/2;
  }
  if (pr) {
    for (j = 0; j < np; j++) 
      printf("point %d:  %9.6f %9.6f %9.6f\n", j+1, p[j][0], p[j][1], p[j][2]);
    printf("h[0] = %f, n[0] = %d\n", h[0], n[0]);
    for (j = 0; j < np; j++) 
      printbbox(level, &bbox[j][0], 0);
  }

  
  

  /*** given individual boxes, join boxes as needed ***/

  /* we always start with two boxes */
  nboxes = np;

  /* do these boxes overlap? */
  for (i=0; i<nboxes-1; i++) {
    for (j=i+1; j<nboxes; j++) {
      overlap = box_overlap(&bbox[i][0], &bbox[j][0]);
      separation = box_separation(h, &bbox[i][0], &bbox[j][0]); 
      if (0) printf("box %d %d (%d) overlap = %d, separation = %d\n", i,j,nboxes,overlap, separation);
      
      /* if overlap, create one box containing both */
      /* if no overlap, keep the two boxes */
      if (overlap) {
        box_join(&bbox[i][0], &bbox[j][0]);
        for (k=j; k<nboxes-1; k++) {
          for (l=0; l<6; l++)
            bbox[k][l] = bbox[k+1][l];
        }
        nboxes--;
        j=i;
      }
      
      if (0)
        for (k=0; k<nboxes; k++)
          printbbox(level, &bbox[k][0], 0);
    }
  }
  
  

  
  
  
  /* handle special quadrant configurations (would be similar for rotant) */
  if (Getv("grid", "quadrant") && nboxes == 2 && np<=2) {

    if (Getv("grid", "onequadrantbox"))
      /* here comes a neat idea to avoid two boxes for quadrant
	 there only are two boxes visible in quadrant when the boxes cross the
	 y=0 plane, so detect that case and temporarily use one box
	 memory peaks at the bar configuration value, but speed will be better
      */
      /* check for y=0 overlap, note that we have to check only one box */
      if (dless(bbox[0][2], nghosts*level->dy) && 
	  dless(-nghosts*level->dy, bbox[0][3])) {
	box_join(&bbox[0][0], &bbox[1][0]);
	nboxes = 1;
      }
  }

  /* adjust boxes for symmetries
     this used to be handled automatically and elegantly by the flag functions
  */
  for (j = 0; j < nboxes; j++) {
    int *half = level->grid->half;
    flag[j] = 0;
    if (half[0]) {
      if (bbox[j][0] < 0) bbox[j][0] = (0.5 - nghosts/2) * level->dx;
      if (bbox[j][1] < 0) flag[j] = 1;
    }
    if (half[1]) {
      if (bbox[j][2] < 0) bbox[j][2] = (0.5 - nghosts/2) * level->dy;
      if (bbox[j][3] < 0) flag[j] = 1;
    }
    if (half[2]) {
      if (bbox[j][4] < 0) bbox[j][4] = (0.5 - nghosts/2) * level->dz; 
      if (bbox[j][5] < 0) flag[j] = 1;
    }
  }

  /* handle quadrant and rotant, two boxes */ 
  if (nboxes == 2 && np<=2 &&
      (Getv("grid", "quadrant") || Getv("grid", "rotant")) 
      && !Getv("grid", "onequadrantbox")) {

    /* pick box that reaches larger y values (implies ymax > 0) */
    k = (bbox[0][3] > bbox[1][3]) ? 0 : 1;
    j = 1 - k;

    /* if ymin is less than zero, but also if it is too close to zero
       for ghost parent level to have ymin > 0, set ymin to symmetry range */
    if (bbox[k][2] < nghosts*level->dy) {
      bbox[k][2] = (0.5 - nghosts/2) * level->dy;
      bbox[j][2] = + bbox[k][2];

      /* make sure there are enough points in y-direction */
      if (Getv("bampi_sizes_constant", "yes"))
	minymax = 2 * Geti("bampi_ysize") * nghosts;
      else
	minymax = 2 * nghosts;
      minymax = 0.5*(minymax - 1) * level->dy;
      if (bbox[j][3] < minymax) bbox[j][3] = minymax;
      if (dless(bbox[k][3], bbox[j][3]))
	errorexit("two box rotant/quadrant: both boxes are too small");
    }

    /* make sure the boxes are symmetric in x-direction */
    bbox[j][0] = - bbox[k][1];
    bbox[j][1] = - bbox[k][0];

    /* did we create a new overlap? */
    if (box_overlap(&bbox[0][0], &bbox[1][0])) {
      box_join(&bbox[0][0], &bbox[1][0]);
      bbox[0][2] = (0.5 - nghosts/2) * level->dy;
      nboxes = 1;
    }

    /* set flag for removal by symmetry */
    flag[0] = flag[1] = 0;
    if (bbox[0][3] < 0) flag[0] = 1;
    if (nboxes == 2 && bbox[1][3] < 0) flag[1] = 1;
  }

  /* remove boxes which are outside a symmetry domain */
  for (k = j = 0; j < nboxes; j++) {
    if (!flag[j]) {
      for (i = 0; i < 6; i++) 
	bbox[k][i] = bbox[j][i];
      k++;
    }
  }
  if (0) printf("nboxes was %d, now setting to %d\n", nboxes, k);
  nboxes = k;


 /* different number of refinements around different centers 
     lmax = largest level number allowed
     lmax2 = largest level number allowed for the 'larger' puncture, means puncture 0 
  */
  if (nboxes == 2 && 
      Geti("amr_lmax2") > 0 &&
      Geti("amr_lmax2") <= level->l) {

   //  one of the boxes has to be removed 
     printf("removing one box: l=%d, lmax=%d, lmax2=%d\n", 
		  level->l+1, Geti("amr_lmax"), Geti("amr_lmax2"));
    
    // we will always remove box 0, puncture 1 will have lmax
      for (i = 0; i < 6; i++) 
	bbox[0][i] = bbox[1][i];
    nboxes = 1;
  } 

  if (nboxes > 2) {
     char str[100];
     for (i = 1; i < np; i++) {
      if (level->grid->lmaxpunc[i] != level->grid->lmaxpunc[0])
       errorexit("differnt number of levels for more than 2 puncture not implemented");
     }
   } 


  /*** final result ***/

  /* store boxes */
  *pnboxes = nboxes;
  for (j = 0; j < nboxes; j++) {
    for (i = 0; i < 6; i++) 
      box[j].bbox[i] = bbox[j][i];
    for (i = 0; i < 6; i+=2)
      box[j].ibbox[i] = 0;
    for (i = 1; i < 6; i+=2)
      box[j].ibbox[i] = (bbox[j][i] - bbox[j][i-1])/h[i/2] + 0.5;
  }

  /* info */
  if (pr) {
    for (j = 0; j < np; j++) 
      printf("point %d:  %9.6f %9.6f %9.6f\n", j+1, p[j][0], p[j][1], p[j][2]);
    for (j = 0; j < nboxes; j++) 
      printbbox(level, box[j].bbox, box[j].ibbox);
  }
}




/* compare points */
int compare_xyz(double x, double y, double z, double a, double b, double c)
{
  if (dless(z,c)) return -1;
  if (dless(c,z)) return  1;  
  if (dless(y,b)) return -1;
  if (dless(b,y)) return  1;
  if (dless(x,a)) return -1;
  if (dless(a,x)) return  1;
  return 0;
}




/* given two boxes, return local overlap as index ranges (ibbox) */
#define iofx(x,x0,dx) floor(((x) - (x0) + dequaleps) / (dx));

int box_overlap_range(int local, tB *boxa, tB *boxb, int *i, int *j)
{
  double *a, *b;
  int *ia, *ib;
  double o[6];
  int k;
  int pr = 0;

  if (local) {
    a = boxa->com->bbox;
    b = boxb->com->bbox;
    ia = boxa->com->ibbox;
    ib = boxb->com->ibbox;
  } else {
    a = boxa->bbox;
    b = boxb->bbox;
    ia = boxa->ibbox;
    ib = boxb->ibbox;
  }

  if (pr) {
    printf("box_overlap_range\n");
    printbbox(boxa->level, a, ia);
    printbbox(boxb->level, b, ib);
  }

  /* initialize empty index ranges */
  i[0] = i[2] = i[4] = j[0] = j[2] = j[4] =  0;
  i[1] = i[3] = i[5] = j[1] = j[3] = j[5] = -1;

  /* find min/max overlap */
  o[0] = max2(a[0], b[0]);
  o[1] = min2(a[1], b[1]);
  o[2] = max2(a[2], b[2]);
  o[3] = min2(a[3], b[3]);
  o[4] = max2(a[4], b[4]);
  o[5] = min2(a[5], b[5]);

  /* return if overlap is empty */
  if (dless(o[1], o[0]) || dless(o[3], o[2]) || dless(o[5], o[4])) {
    if (pr) printf("-> no overlap\n");
    return 0;
  }

  /* convert floating point overlap into index ranges */
  for (k = 0; k < 6; k++) {
    i[k] = iofx(o[k], a[(k/2)*2], boxa->dd[k/2]);
    j[k] = iofx(o[k], b[(k/2)*2], boxb->dd[k/2]);
  }

  if (pr) {
    printf("-> overlap\n");
    printbbox(boxa->level, o, 0);
    printf("%2d %2d  %2d %2d  %2d %2d\n", i[0], i[1], i[2], i[3], i[4], i[5]);
    printf("%2d %2d  %2d %2d  %2d %2d\n", j[0], j[1], j[2], j[3], j[4], j[5]);
  }

  return 1;
}




/* non-parallel move overlap
   should loop over points come first?
*/
void move_overlap_nonparallel(tL *lold, tL *lnew, tVarList *all)
{
  int c, i, iall, iv, j, k;
  double *flag = Ptr(lnew, "flagprolong");
  double *xold = Ptr(lold, "x");
  double *yold = Ptr(lold, "y");
  double *zold = Ptr(lold, "z");
  double *xnew = Ptr(lnew, "x");
  double *ynew = Ptr(lnew, "y");
  double *znew = Ptr(lnew, "z");
  double x, y, z;
  double *old, *new;
  int pr = Getv("amr_verbose","all")+Getv("amr_verbose","move");
  
  /* initialize flag */
  forallpoints(lnew, i) flag[i] = 1;

  /* pairwise box copy */
  int iold[6], inew[6];
  int io, jo, ko, no;
  int in, jn, kn, nn;
  if (pr) printf("move_overlap_nonparallel\n");

  if (0) {
    printf("old global:\n");
    for (j = 0; j < lold->nboxes; j++) 
      printbbox(lold, lold->box[j]->bbox, lold->box[j]->ibbox);
    printf("new global:\n");
    for (j = 0; j < lnew->nboxes; j++) 
      printbbox(lnew, lnew->box[j]->bbox, lnew->box[j]->ibbox);
    printf("\nold local:\n");
    for (j = 0; j < lold->nboxes; j++) 
      printbbox(lold, lold->box[j]->com->bbox, lold->box[j]->com->ibbox);
    printf("new local:\n");
    for (j = 0; j < lnew->nboxes; j++) 
      printbbox(lnew, lnew->box[j]->com->bbox, lnew->box[j]->com->ibbox);
  }      

  /* consider all pairs of old and new boxes */ 
  for (i = 0; i < lold->nboxes; i++)  
  for (j = 0; j < lnew->nboxes; j++) {

    /* find overlap in the form of index ranges */
    k = box_overlap_range(1, lold->box[i], lnew->box[j], iold, inew); 
    
    /* nothing to do if there is no overlap */
    if (!k) continue;

    /* for all variables */
    for (iv = 0; iv < all->n; iv++) {
      iall = all->index[iv];
      old = lold->v[iall];
      new = lnew->v[iall];
      if (!old || !new) continue;

      /* copy overlap */
      for (ko = iold[4], kn = inew[4]; ko <= iold[5]; ko++, kn++)
      for (jo = iold[2], jn = inew[2]; jo <= iold[3]; jo++, jn++)
      for (io = iold[0], in = inew[0]; io <= iold[1]; io++, in++) {
        no = ijkofbox(lold->box[i], io, jo, ko); // not the fastest version
        nn = ijkofbox(lnew->box[j], in, jn, kn);

        /* copy inner or copy all */
        const int copyall = 1;
        if (copyall || boundaryflag(lnew, nn) == 0) {
          flag[nn] = 0;
          new[nn] = old[no];
        }
      }
    }
  }


  /* synchronization takes care of certain points flagged in the interior */
  /* -> assume that sync is called elsewhere */
}




/* parallel move overlap */
void move_overlap_parallel(tL *lold, tL *lnew, tVarList *vl)
{
  double *flag = Ptr(lnew, "flagprolong");
  double *xnew = Ptr(lnew, "x");
  double *ynew = Ptr(lnew, "y");
  double *znew = Ptr(lnew, "z");
  double *coord = dmalloc(3*lnew->nnodes);
  double *data;
  int i, j, nfound;
  int nvl = vl->n;
  int npoints = 0;
  int pr = Getv("amr_verbose","all")+Getv("amr_verbose","move");
  
  if (pr) printf("move_overlap_parallel\n");
  
  /* create list of all points  */
  forallpoints(lnew, i) {
    coord[3*i] = xnew[i];
    coord[3*i+1] = ynew[i];
    coord[3*i+2] = znew[i];
    npoints++;
  }

  /* get data */
  data = dmalloc(npoints * (nvl+1));
  timer_start(0, "move_box_getd");
  bampi_getdata(lold, vl, npoints, coord, data);
  timer_stop (0, "move_box_getd");


  /* if data was found, store data, else set flag
      stay away from ghosts and symmetry points, they are handled elsewhere
  */
  /* use same loop for getdata above ... */
  forallpoints(lnew, i) { 
    if (boundaryflag(lnew, i) == SYMBOUND ||
        boundaryflag(lnew, i) == GHOBOUND)
      continue;
    nfound = data[(nvl+1)*i + nvl];
    if (nfound > 0) {
      for (j = 0; j < nvl; j++)
        lnew->v[vl->index[j]][i] = data[(nvl+1)*i + j] / nfound;
    }
  }

    /* added to fix parallisation probem 
     2023 - > to be watched! */
    if (1) bampi_vlsynchronize(vl);

  
  /* clean up */
  free(data);
  free(coord);
}




/* move overlap from old to new boxes 
   requires interpolation into new points from coarse level
   plus a copy from available fine data
   - most frequent case, and the case that dominates run-time, is that
     the boxes only change slightly:
     - e.g. a box does not change in size but moves by 2 points
     - e.g. a box changes size by 4 points because its internal boxes move by 2
     --> requires only a local copy if bampi's sync buffer is large enough (!)
     --> punctures usually do not move faster than 1 point per time step
   - sometimes, but only rarely over the course of a run, 
     the box change is larger than the sync buffers
     - e.g. when two boxes merge to one, or one splits into two
     --> requires general parallel copy operation (slow)
*/
void move_overlap(tL *lc, tL *lold, tL *lnew, tVarList *vlc, tVarList *vlnew)
{
  int iold[6], inew[6];
  double *flagnewsave;
  int i, j, k, max;
  int local = 0;
  int pr = Getv("amr_verbose","all")+Getv("amr_verbose","move");
    
  /* figure out whether only a local copy is required */
  /* could be refined for more special cases, but the main goal is
     to detect the most frequent case in which the boxes change only slightly 
  */
  /* consider only equal number of boxes as candidates for local copies */
  if (lold->nboxes == lnew->nboxes) {

    /* find the maximum number of points that any of the six faces of any 
       of the boxes moves */
    max = -1;
    for (i = 0; i < lnew->nboxes; i++) {

      /* find index ranges for overlap of global box */
      k = box_overlap_range(0, lold->box[i], lnew->box[i], iold, inew); 
      /* should also have a local box version, which needs a communication
         but would be completely safe */

      if (!k) 
	max = 1000; // if the boxes don't overlap, give up
      else
	for (k = 0; k < 6; k++)
	  if ((j = abs(inew[k]-iold[k])) > max) max = j; 
      // if (max > 4) break;
    }

    /* typically we have 4 ghosts, and even if opposite sides move by 4 each
       the point distribution across processors should be fine
       check: this criterion may not be sufficient in general! */
    if (max <= Geti("bampi_nghosts")) 
      local = 1;
  }

  /* general method that uses parallel copy (slow) */
  if (!local) {

    if (pr) printf("move overlap: global\n");
    
    /* prolong from coarse into the new level
    -1: prolong everywhere using polynomial interpolation
        using flag[] is not working due to the prolong algorithm
    */
    timer_start(0, "move_box_prol");
    prolong_varlist(lc->grid, lc->l, lnew->l, vlc, vlnew, -1); 
    timer_stop(0, "move_box_prol");
    
    /* copy data that was already available on the old level 
    uses general purpose parallelized copy 
    */
    timer_start(0, "move_box_copy");
    move_overlap_parallel(lold, lnew, vlnew);
    timer_stop(0, "move_box_copy"); 

  /* special method that uses local copy (fast) */
  } else {

    /* copy data that was already available on the old level 
       flag those points for which data was not available
       interpolate with special -3 flag which does a syncparent all but
       interpolates only into flagged points
    */
    if (pr) printf("move overlap: local\n");
    
     /* do not overwrite the prolong informations, use a new
    temporary field instead */
    j = Ind("flagprolong");
    flagnewsave = lnew->v[j];
    lnew->v[j] = dcalloc(lnew->npoints);
    
    timer_start(0, "move_box_copy");
    move_overlap_nonparallel(lold, lnew, vlnew);
    timer_stop(0, "move_box_copy");
    
    timer_start(0, "move_box_prol");
    prolong_varlist(lc->grid, lc->l, lnew->l, vlc, vlnew, -3);
    timer_stop(0, "move_box_prol");
    
    free(lnew->v[j]);
    lnew->v[j] = flagnewsave;
  }
  
}




/* move a refinement level
   also grows/shrinks boxes, i.e. it does a regrid
*/
int move_box(tG *g, int l) 
{ int lint = l;
  tL *lold, *lnew, *lc;
  tVarList *vlc, *vlnew;
  tVarList *u, *up, *upp, *up2, *up3, *up4;
  int i, j;
  int pr = Getv("amr_verbose","all")+Getv("amr_verbose","move");

  /* check whether we are on */
  if (!Getv("amr", "move")) return 0;

  /* never move top level */
  if (l <= g->lmin) return 0;

  /* never move level1 when we have shells */
  if (l == g->lmin+1 && Getv("grid","shells")) return 0;
  
  /* the level up for movement may not exist if move_box called on l+1 after 
     advance of finest level l
  */
  if (l > g->lmax) return 0;
  if (0) printf("l = %d\n", l);

  /* wait until level l-1 is the coarsest level of its time 
     in particular, l-1 to lmax are time aligned but l-2 (if it exists) is not
     we are called after advancing l-1 and all sublevels, so they are aligned
  */
  if (Getv("grid","shells")) {
    /* in case of shells we need one more level which does not move */
    if (l-3 >= g->lmin && !dless(g->level[l-1]->time, g->level[l-2]->time))
      return 0;
  } else {
    if (l-2 >= g->lmin && !dless(g->level[l-1]->time, g->level[l-2]->time))
      return 0;
  }
  
  /* timer */
  timer_start(0, "move_box");


  /* now move all boxes that need it from coarsest to finest
     this way the existing grid generation takes care of the ch/pr pointers
     (as opposed to moving a grid in the middle while resetting all pointers)
  */
  int moved = 0;
  for (; l <= g->lmax; l++) {
    tL *lc = g->level[l-1];
    double *f = Ptr(lc, "flagregrid");

    double *x = Ptr(lc, "x");
    double *y = Ptr(lc, "y");
    double *z = Ptr(lc, "z");
    tSBox *box = bcalloc(3 * sizeof(tSBox));
    int nboxes;
    
    if (0) printf("move_box l=%d\n", l);
    
    /* determine bounding boxes based on global nesting function */
    find_nested_bboxes(lc, &nboxes, box);
    if (0) {
      printf("lc->sbox %p\n", lc->sbox);
      if (lc->sbox) {
	printf("lc->nsboxes %d\n", lc->nsboxes);
	for (j = 0; j < lc->nsboxes; j++) 
	  printbbox(lc, lc->sbox[j].bbox, lc->sbox[j].ibbox);
      }
    }

    /* if the previous boxes are equal to the new ones, continue */
    if (lc->sbox && lc->nsboxes == nboxes) {
      for (j = 0; j < nboxes; j++) 
	if (!box_equal(lc->sbox[j].bbox, box[j].bbox)) break;
      if (j == nboxes) {
	free(box);
	continue;
      }
    }
    
    /* now we know that we move something */
    moved = 1;

    /* info */
    if (pr) {
      if (lc->sbox && lc->nsboxes) {
	printf("\nmoving box(es) on level %d from\n", l);
	for (j = 0; j < lc->nsboxes; j++) 
	  printbbox(lc, lc->sbox[j].bbox, lc->sbox[j].ibbox);
	printf("to\n");
      } else {
	printf("\nsetting box(es) on level %d to\n", l);
      }
      for (j = 0; j < nboxes; j++) 
	printbbox(lc, box[j].bbox, box[j].ibbox);
    }

    /* new boxes, store */
    free(lc->sbox);
    lc->sbox = box;
    lc->nsboxes = nboxes;
    
    /* save pointer to old level */
    lold = g->level[l];

    /* free ghost parent now so that there is no confusion later 
       tricky bit is that the grid information may have changed one level up
    */
    if (lold->prlocal) {
      if (bampi_size() > 1 || lold->nboxes > 1)
	free_level(lold->prlocal);

      /* either prlocal was freed because it was a ghost, or 
	 prlocal is pointing to a non-ghost level that no longer exists
	 set to zero so that free(lold) below doesn't touch it
      */
      lold->prlocal = 0;
    }

    /* create new level based on bounding boxes */
    lnew = make_level_bboxes(lc, nboxes, box, 0);

    /* the new level starts at time of old level */
    lnew->time = lold->time;
    lnew->iteration = lold->iteration;
    lnew->dt = lold->dt;
    
    /* create list of all variables that have to be transfered */
    if (0) printvariables(lnew);
    vlc = vlalloc(lc);
    evolve_vlretrieve(&u, &up, &upp);
    vlpushvl(vlc, u);
    vlpushvl(vlc, up);
    vlpushvl(vlc, upp);
    evolve_vlretrieve_pn(&up2, &up3, &up4);
    vlpushvl(vlc, up2);
    vlpushvl(vlc, up3);
    vlpushvl(vlc, up4);
    if (0) prvarlist(vlc);
    vlnew = vlalloc(lnew);
    vlpushvl(vlnew, vlc);

    /* test new level if it is REALLY consistent */
    test_box_consistency_level(lc,lnew);
    
    /* move overlap
       copy where data is available, interpolate from coarse where needed
    */
    move_overlap(lc, lold, lnew, vlc, vlnew);

    /* remove old level */

    /* by design the new level does not yet have the information
       on the current boxes, i.e. level->sbox = 0
       this triggers a "move box" where the points stay in place
       -> reinitializes level (is this really needed, say for BOX mode?)
       could transfer box to new level before deleting it
       
       lnew->sbox = lold->sbox;        // or = box
       lnew->nsboxes = lold->nsboxes;  // or = nboxes
       lold->sbox = 0;            
       lold->nboxes = 0; // so that free_level(lold) doesn't remove it
    */ 

    /* now remove old level */
    free_level(lold);
    
    /* clean up */
    vlfree(vlc);
    vlfree(vlnew);

    /* hack for camr, one could also try to use c_amr_mask in POST_MOVEBOX
    if ((Getv("conservative_amr", "yes")) && (Getv("physics", "matter"))) { 
        if (MATTER.CAMRACTIVE == 1){ 
         if (PR) printf("recompute camr mask \n"); 
         c_amr_mask(g->level[l]) ; }
     }*/

  }
  timer_stop(0, "move_box");

  return moved;
}
