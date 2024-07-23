/* points.c */
/* Bernd Bruegmann 5/02 */

/* Basic problem with a node based amr system is to locate points.
   Provide some basic point location routines here.
*/

#include "bam.h"
#include "amr.h"




/* find bounding box for the points in this level 
   also returns index range for box if ibbox is not NULL
*/
void findbbox(tL *level, double *bbox, int *ibbox) 
{
  findbbox_flag(level, bbox, ibbox, 0);
}




/* find bounding box for the points in this level that have a given flag set 
   consider all points if flag is NULL
   also returns index range for box if ibbox is not NULL
*/
void findbbox_flag(tL *level, double *bbox, int *ibbox, double *flag) 
{
  int i, n = level->nnodes;
  double *x  = level->v[Ind("x")];
  double *y  = level->v[Ind("y")];
  double *z  = level->v[Ind("z")];
  
  bbox[0] = bbox[2] = bbox[4] =  DBL_MAX;
  bbox[1] = bbox[3] = bbox[5] = -DBL_MAX;

  for (i = 0; i < n; i++) {
    if (flag && !flag[i]) continue;
    if (x[i] < bbox[0]) bbox[0] = x[i];
    if (x[i] > bbox[1]) bbox[1] = x[i];
    if (y[i] < bbox[2]) bbox[2] = y[i];
    if (y[i] > bbox[3]) bbox[3] = y[i];
    if (z[i] < bbox[4]) bbox[4] = z[i];
    if (z[i] > bbox[5]) bbox[5] = z[i];
  }

  if (ibbox) {
    ibbox[0] = 0;
    ibbox[1] = (bbox[1] - bbox[0])/level->dx + 0.5;
    ibbox[2] = 0;
    ibbox[3] = (bbox[3] - bbox[2])/level->dy + 0.5;
    ibbox[4] = 0;
    ibbox[5] = (bbox[5] - bbox[4])/level->dz + 0.5;
  }  

  if (0) printbbox(level, bbox, ibbox);
}




/* find bounding box for the points in this level that are not ghosts */
void findbbox_notghost(tL *level, double *bbox, int *ibbox)
{
  int b, i, n = level->nnodes;
  double *x  = level->v[Ind("x")];
  double *y  = level->v[Ind("y")];
  double *z  = level->v[Ind("z")];
  
  bbox[0] = bbox[2] = bbox[4] =  DBL_MAX;
  bbox[1] = bbox[3] = bbox[5] = -DBL_MAX;

  for (i = 0; i < n; i++) {
    b = boundaryflag(level, i);
    if (b == GHOBOUND || b == SYMBOUND) continue;
    if (x[i] < bbox[0]) bbox[0] = x[i];
    if (x[i] > bbox[1]) bbox[1] = x[i];
    if (y[i] < bbox[2]) bbox[2] = y[i];
    if (y[i] > bbox[3]) bbox[3] = y[i];
    if (z[i] < bbox[4]) bbox[4] = z[i];
    if (z[i] > bbox[5]) bbox[5] = z[i];
  }

  if (ibbox) {
    ibbox[0] = 0;
    ibbox[1] = (bbox[1] - bbox[0])/level->dx + 0.5;
    ibbox[2] = 0;
    ibbox[3] = (bbox[3] - bbox[2])/level->dy + 0.5;
    ibbox[4] = 0;
    ibbox[5] = (bbox[5] - bbox[4])/level->dz + 0.5;
  }  
}




/* check whether a point is inside the bounding box */
int insidebbox(double *bbox, double *coord)
{
  if (dless(coord[0],bbox[0]) ||
      dless(coord[1],bbox[2]) ||
      dless(coord[2],bbox[4]) ||
      dless(bbox[1],coord[0]) ||
      dless(bbox[3],coord[1]) ||
      dless(bbox[5],coord[2]))
    return 0;
  return 1;
}

int xyzinsidebbox(double *bbox, double x, double y, double z)
{
  double coord[3];
  coord[0] = x;
  coord[1] = y;
  coord[2] = z;
  return insidebbox(bbox, coord);
}

int xyzinsidelevel(tL* level, double x, double y, double z)
{
  int inside = 0;
  
  forallboxes(level) {
    inside += xyzinsidebbox(box->bbox, x,y,z);
  } endforboxes;
  
  if (inside>1 && !level->shells) {
    errorexit("point is in multiple boxes?!?!");
  } 
  return inside;
}

/* check wether this ijk point is also inside finer level*/
int ijkinsidefinerlevel(tB *box, int ijk)
{
    tL* level = box->level;
    if (level->l==level->grid->lmax) return 0;
    
    // bbox of the finer level
    double *bbox = level->grid->level[level->l+1]->bbox;
    double h     = level->grid->level[level->l+1]->dx;
    double off   = h * Getd("amr_nbuffer");
    
    double x = Ptr(level,"x")[ijk];
    double y = Ptr(level,"y")[ijk];
    double z = Ptr(level,"z")[ijk];
     
    int inside = 1;
    
    if (dless(x,bbox[0]+off) ||
        dless(y,bbox[2]+off) ||
        dless(z,bbox[4]+off) ||
        dless(bbox[1]-off,x) ||
        dless(bbox[3]-off,y) ||
        dless(bbox[5]-off,z))
        inside = 0;
    
    if (level->boundary[ijk]) inside = -1;
    
    return inside;
}





/* check whether a point is owned by this processor */
int insideownership(tL *level, double x, double y, double z)
{
  forallboxes(level) {

    if (xyzinsidebbox(box->com->bboxown, x, y, z))
      return 1;

  } endforboxes;

  return 0;
}




/* check whether a sphere is inside a bounding box 
   the caller has to set the bbox correctly, say in the case of symmetries
*/
int sphere_inside_bbox(double *bbox, 
		       double x0, double y0, double z0, double r)
{
  double a[6];

  a[0] = x0 - r;
  a[1] = x0 + r;
  a[2] = y0 - r;
  a[3] = y0 + r;
  a[4] = z0 - r;
  a[5] = z0 + r;

  return box_ainb(a, bbox);
}




/* check whether a sphere is inside a level 
   return 0 if there is more than one box
   return 0 if there is one box that does not cover origin
   takes care of symmetries
*/
int sphere_inside_level(tL *level, 
			double x0, double y0, double z0, double r)
{
  double b[6];
  int *half = level->grid->half;
  int i;

  if (level->shells) {
    double rmin = level->bbox[0] + Geti("order_RP") * level->dx;
    double rmax = level->bbox[1] - Geti("order_RP") * level->dx;
    r -= Geti("order_RP") * level->dx;
    double r1 = r - sqrt(x0*x0+y0*y0+z0*z0);
    double r2 = r + sqrt(x0*x0+y0*y0+z0*z0);
    if (r1<=0. || r1<rmin || r2<rmin || r1>rmax || r2>rmax )
      return 0;
    else
      return 1;
  }

  if (level->nboxes > 1) return 0;
  if (!xyzinsidebbox(level->bbox, x0, y0, z0)) return 0;

  for (i = 0; i < 6; i++) b[i] = level->bbox[i];
  for (i = 0; i < 3; i++)
    if (half[i]) b[2*i] = - b[2*i+1];

  return sphere_inside_bbox(b, x0, y0, z0, r);
}




/* check whether a sphere centered on origin is inside a one box level
   including the points for centered interpolation
*/
int centered_sphere_safely_inside_level(tL *level, double r)
{
  r += Geti("order_RP") * level->dx;
  return sphere_inside_level(level, 0, 0, 0, r);
}

int sphere_safely_inside_level(tL *level, 
			       double x0, double y0, double z0, double r)
{
  r += Geti("order_RP") * level->dx;
  return sphere_inside_level(level, x0, y0, y0, r);
}




/* return buffer with data that is found on this processor for the given 
   list of points 

   caller has to free buffer

   assumes ordered points
   ordered points makes this an O(N) operation

   if points are not ordered but also not completely random, then
   one could reverse linear search direction when the sequence of points
   changes order

   format of buf:   index in point list, data variable 0, data variable 1, ...
   (!) the key trick is to store the index, too, so that we have no problem
   combining data from different processors
*/
void bufferorderedpoints(tL *level, int npoints, double *coords, 
			 int nv, int *iv, int *nbuf, double **buf, 
			 int interpolate, int order, int scheme)
{
  double *bbox = level->com->bbox;
  double *vinterp = calloc(nv, sizeof(double));
  int n = level->nnodes;
  int i, ii, j, k;
  double *buffer;
  int nb;
  double xp, yp, zp, xn, yn, zn; 
  double *x  = level->v[Ind("x")];
  double *y  = level->v[Ind("y")];
  double *z  = level->v[Ind("z")];
  
  /* if not done already, allocate buffer of sufficient size */
  if (!*buf) {
    nb = 0;
    for (j = 0; j < npoints; j++)
      if (insidebbox(bbox, coords+3*j)) nb++;
    buffer = calloc(sizeof(double), nb*(nv+1));
  } else
    buffer = *buf;
  nb = 0;

  /* for all points inside the box */
  for (i = j = 0; i < n && j < npoints; j++) {
    if (!insidebbox(bbox, coords+3*j)) continue;
    xp = coords[3*j];
    yp = coords[3*j+1];
    zp = coords[3*j+2];

    /* if we are not allowed to interpolate */
    if (!interpolate) {

      /* try to find point as is in node list */
      for (ii = i; ii < n; ii++) {
	xn = x[ii];
	yn = y[ii];
	zn = z[ii];
	if (dequal(xp,xn) && dequal(yp,yn) && dequal(zp,zn)) break;
      }

      /* if point found */
      if (ii < n) {
	
	/* advance i */
	i = ii;
	
	/* don't collect ghost points */
	if (boundaryflag(level, i) != GHOBOUND) {
	  
	  /* save index into point list */
	  buffer[nb++] = (double) j;
	  
	  /* save data, loop should be inverted */
	  for (k = 0; k < nv; k++) 
	    buffer[nb++] = level->v[iv[k]][i];
	}
      }
    }

    /* if we are asked to interpolate */
    else {
      
      ii = interpolate_xyz_local(level, xp, yp, zp, nv, iv, vinterp, i, order,scheme);
      
      /* if interpolation was possible */
      if (ii != -1) {
	
	/* advance i */
	i = ii;
	
	/* save index into point list */
	buffer[nb++] = (double) j;
	
	/* save data */
	for (k = 0; k < nv; k++) 
	  buffer[nb++] = vinterp[k];
      }
    }

    /* if point not found or interpolation not possible, 
       do nothing, in particular
       do not advance i because if there is another point
       to be found, it will be for some index > i
    */
  }

  /* return result */
  free(vinterp);
  if (!*buf) {
    *nbuf = nb;
    *buf = buffer;
  }
}






/* return buffer with data that is found on this processor for the given 
   list of points 

   caller has to free buffer

   does not assume ordered points

   this is for periodic boundaries, and is a temporary fix because I
   don't know how to sort a list of coordinates.

*/
void buffernonorderedpoints(tL *level, int npoints, double *coords, 
			    int nv, int *iv, int *nbuf, double **buf, 
			    int interpolate)
{
  double *bbox= level->com->bbox;
  double *vinterp = calloc(sizeof(double), nv);
  int n = level->nnodes;
  int i, ii, j, k;
  double *buffer;
  int nb;
  double xp, yp, zp, xn, yn, zn; 
  double *x  = level->v[Ind("x")];
  double *y  = level->v[Ind("y")];
  double *z  = level->v[Ind("z")];
  /* d should be the d of the smallest level if there are more than 1*/
  double dx = level->dx; 
  double dy = level->dy; 
  double dz = level->dz; 

  /* get upper limit on required buffer size */
  nb = 0;
  for (j = 0; j < npoints; j++) {
    if (insidebbox(bbox, coords+3*j)) nb++;
  }
  buffer = calloc(sizeof(double), nb*(nv+1));

  nb = 0;

  /* for all points inside the box */
  for (j = 0; j < npoints; j++) {
    if (!insidebbox(bbox, coords+3*j)) continue;
    xp = coords[3*j];
    yp = coords[3*j+1];
    zp = coords[3*j+2];
    /* if we are not allowed to interpolate */
    if (!interpolate) {
      /* try to find point as is in node list */
      for (ii = 0; ii < n; ii++) {
	xn = x[ii];
	yn = y[ii];
	zn = z[ii];
	if (((xn > xp-0.1*dx)&&(xn < xp+0.1*dx))&&
	    ((yn > yp-0.1*dy)&&(yn < yp+0.1*dy))&&
	    ((zn > zp-0.1*dz)&&(zn < zp+0.1*dz))) break;
      }

      /* if point found */
      if (ii < n) {
	/* don't collect ghost points */
	if (boundaryflag(level, ii) != GHOBOUND) {
	  /* save index into point list */
	  buffer[nb++] = (double) j;	  
	  /* save data, loop should be inverted */
	  for (k = 0; k < nv; k++) 
	    buffer[nb++] = level->v[iv[k]][ii];
	}
      }
    }
    /* if we are asked to interpolate */
    else {
      errorexit("sorry, can't interpolate non-ordered points");
    }

    /* if point not found or interpolation not possible, 
       do nothing, in particular
       do not advance i because if there is another point
       to be found, it will be for some index > i
    */
  }

  /* return result */
  free(vinterp);
  *nbuf = nb;
  *buf = buffer;
}





/* return buffer with data that is found on this processor for the given 
   list of points for box and shells ... slow
*/
void bufferpoints_full(tL *level, int npoints, double *coords, 
                       int nv, int *iv, int *nbuf, double **buf, 
                       int order, int scheme)
{
  int i, j;
  int nb = 0;
  double xp, yp, zp;
  double *vinterp = calloc(nv, sizeof(double));
  double *buffer = calloc(sizeof(double), npoints*(nv+1));
  
  
  /* for all points in the coords array */
  for (i = 0; i < npoints; i++)
  {
    forallboxes(level)
    {
      /* check if point i is in this box, if not continue with next box */
      if(level->shells)
      {
        if(box->i != 
           find_shellsbox_from_xyz(coords[3*i],coords[3*i+1],coords[3*i+2]))
          continue;
      }
      else
      {
        if(!xyzinsidebbox_withsym(level, box->bbox, 
                                  coords[3*i],coords[3*i+1],coords[3*i+2]))
          continue;
      }

      /* store point, convert to spherical coordinates if we are in shells */
      if(level->shells)
      {
        convert_box_to_shells(box->i, 
          coords[3*i],coords[3*i+1],coords[3*i+2], &xp,&yp,&zp);
          /* -> xp,yp,zp are now the shell coordinates */
      }
      else
      {
        xp = coords[3*i];
        yp = coords[3*i+1];
        zp = coords[3*i+2];
      }
      
      /* interpolate only if we have the point */
      if(xyzinsidebbox_withsym(level, box->com->bbox, xp,yp,zp))
        if(interpolate_xyz_localinbox_minimal_withsym(box, xp,yp,zp, nv, iv, vinterp, order, scheme))
        {
          //printf("%d %d %d    %e\n",nb,i,npoints, vinterp[0]);

          /* save index into point list */
          buffer[nb++] = (double) i;

          /* save data */
          for (j = 0; j < nv; j++) 
            buffer[nb++] = vinterp[j];
        
          break;
        }

    } endforboxes;
  }
  
  free(vinterp);
  
  if (*buf)
    free(*buf);
  
  *buf  = buffer;
  *nbuf = nb;
}






/* find one single point
   somewhat better than searching through all points
   for ordered list of points, use the routines above for better efficiency
   possible improvements:
   - use parent oct tree for log(N) search
   - walk in all three directions (may want to require convex level)
*/
#define compareandcheck(level,x,y,z,xp,yp,zp,i) \
  (dequal(x,xp[i]) && dequal(y,yp[i]) && dequal(z,zp[i]) && \
   (boundaryflag(level,i) == NOTBOUND ||	\
    boundaryflag(level,i) == PHYBOUND || \
    boundaryflag(level,i) == EXCBOUND))

int find_one_point(tL *level, double *coord)
{
  static int i0 = 0;
  static tL *level0 = 0;
  static double *xp = 0, *yp = 0, *zp = 0;
  double x, y, z;
  int i;

  if (0) printf("find_one_point %p %p\n", level0, level);

  /* reset */
  if (level == 0) {
    level0 = level;
    return -1;
  }
  
  /* reinitialize */
  if (level != level0) {
    level0 = level;
    i0 = 0;
    xp = Ptr(level, "x");
    yp = Ptr(level, "y");
    zp = Ptr(level, "z");
  }

  /* check bounding box */
  if (!insidebbox(level->com->bbox, coord))
    return -1;

  /* coordinates of point that we are looking for */
  x = coord[0];
  y = coord[1];
  z = coord[2];

  /* search forward exploiting that coordinates are ordered */
  if (dless(zp[i0],z) || dless(yp[i0],y) || dless(xp[i0],x)) {
    for (i = i0+1; i < level->nnodes; i++)
      if (compareandcheck(level, x, y, z, xp, yp, zp, i))
	return (i0 = i);
  } 

  /* search backward, include current point */
  else {
    for (i = i0; i >= 0; i--)
      if (compareandcheck(level, x, y, z, xp, yp, zp, i))
	return (i0 = i);
  }

  /* point not found */
  return -1;
}





/* find one single point
   assumes box (could be convex)
   -> walk into 3d direction of point
*/
/* fix: done like this there is unwanted/uncontrolled snap to grid 
        effect if the point is not part of the grid
*/ 
int find_one_point_box(tL *level, double x, double y, double z)
{
  forallboxes(level) {

    if (xyzinsidebbox(box->com->bbox, x, y, z)) {

      int i = floor((x - box->x0)/box->dx + 0.5);
      int j = floor((y - box->y0)/box->dy + 0.5);
      int k = floor((z - box->z0)/box->dz + 0.5);

      return ijkofbox(box, i, j, k);
    }

  } endforboxes;

  return -1;
}





/* ******************************************************************** */
/* WT: Here are some utility routines to use point lists of the         */
/*     type tPointList *, defined in bam_amr.h                          */
/* ******************************************************************** */

/* Allocate memory for a PointList */
tPointList *AllocatePointList(tL *level)
{
 tPointList *PL;
 
 PL=calloc(1, sizeof(*PL) );

 PL->level=level;
 PL->npoints=0;
 PL->point=NULL;

 // printf("AllocatePointList: PL=%p\n",PL); 
 return PL;
}




/* add one point to a PointList */
void AddToPointList(tL *level, tPointList *PL, int newpoint)
{
 void *ret;
 
 if(level!=PL->level) 
   errorexit("AddToPointList: It is forbidden to add points from different "
             "levels to one PointList");

 ret=realloc(PL->point, (sizeof( *(PL->point) ))*(PL->npoints+1) );
 if(ret==NULL) return;
 
 PL->point=ret;
 PL->point[PL->npoints]=newpoint;
 PL->npoints++;
}




/* free a PointList */
void FreePointList(tPointList *PL)
{
 if(PL!=NULL)
 {
   free(PL->point);
   free(PL);
 }
}




/* print a PointList */
void prPointList(tPointList *PL)
{
  int i;
  
  printf("PointList=%p  PointList->level=%p  PointList->npoints=%d\n",
         PL, PL->level, PL->npoints);

  printf("PointList->point = ");
  for(i=0; i<PL->npoints; i++)
    printf("%d ",PL->point[i]);
  printf("\n");
}




/* ******************************************************************** */
/* WT: Here are some utility routines to use Multi level point          */
/*     lists of the type tMlPointList *, defined in bam_amr.h             */
/* ******************************************************************** */

/* Allocate memory for a MlPointList */
tMlPointList *AllocateMlPointList(void)
{
 tMlPointList *MlPL;
 
 MlPL=calloc(1, sizeof(*MlPL) );

 MlPL->nlevels=0;
 
 MlPL->level=NULL;
 MlPL->PointList=NULL;
 MlPL->SlPointListResult=NULL;
 
 return MlPL;
}


/* find the single point list MlPL->PointList[j] in MlPL, which contains the 
   points of level level */
tPointList *SlPointList(tL *level, tMlPointList *MlPL)
{
 int j;
 
 for(j=0; j<MlPL->nlevels; j++)  if(MlPL->level[j] == level) break;
 // printf("j=%d\n",j);
 if(j>=MlPL->nlevels)  MlPL->SlPointListResult = NULL;
 else                  MlPL->SlPointListResult = MlPL->PointList[j];
 
 return MlPL->SlPointListResult;
}

/* add a level without any points to a MlPL */
void AddLevelToMlPointList(tL *level, tMlPointList *MlPL)
{
 void *ret;
 tPointList *PL;

 PL = SlPointList(level, MlPL);
 if( PL == NULL )
 {
   ret=realloc(MlPL->level, (sizeof( *(MlPL->level) ))*(MlPL->nlevels+1) );
   if(ret==NULL) return;
   MlPL->level=ret;
   MlPL->level[MlPL->nlevels]=level;
   
   ret=realloc(MlPL->PointList, 
               (sizeof( *(MlPL->PointList) ))*(MlPL->nlevels+1) );
   if(ret==NULL) return;
   MlPL->PointList=ret;
   MlPL->PointList[MlPL->nlevels]=AllocatePointList(level);
   MlPL->nlevels++;
 }
}


/* add one point to a MlPL */
void AddToMlPointList(tL *level, tMlPointList *MlPL, int newpoint)
{
 tPointList *PL;

 AddLevelToMlPointList(level, MlPL);
 PL = SlPointList(level, MlPL);
 AddToPointList(level, PL, newpoint);
}


/* free a MlPointList */
void FreeMlPointList(tMlPointList *MlPL)
{
 int j;
 if(MlPL!=NULL)
 {
   for(j=0; j<MlPL->nlevels; j++)
     FreePointList(MlPL->PointList[j]);
   free(MlPL->level);
   free(MlPL->PointList);
   free(MlPL);
 }
}


/* print a MlPointList */
void prMlPointList(tMlPointList *MlPL)
{
  int i;
  
  printf("MlPointList=%p MlPointList->nlevels=%d "
         "MlPL->SlPointListResult=%p\n",
         MlPL, MlPL->nlevels, MlPL->SlPointListResult);

  for(i=0; i<MlPL->nlevels; i++)
    prPointList(MlPL->PointList[i]);
}

