/* shells.c */
/* mth 04-06/11 */
/* adding spherical shells to 0 level  -> 6 boxes 
   prolong restrict is handled in a special (non-bam-generic) way
   added possibility to use symmetry for boxes and no symmetry for shells
*/ 


#include "bam.h"
#include "amr.h"

#define PR 0
#define debug 1

int nbuf_x = 6;
int nbuf_r = 6;
int nbuf_phi = 4;

double fisheye_s;
double fisheye_R;
double fisheye_e;




/* coordinate transo from cartesian to spherical shells */ 
void convert_box_to_shells(int b, 
                           double x, double y, double z,
                           double *r, double *phi, double *theta)
{
  // see Thornburg paper
  // I use a different convention ... britty arbitrary 
  *r = sqrt(x*x + y*y + z*z);
  
  if (b==0 || b==1) {
    *phi   = atan(y/x);
    *theta = atan(z/x);
  } else if (b==2 || b==3) {
    *phi   = atan(x/y);
    *theta = atan(z/y);
  } else {
    *phi   = atan(x/z); 
    *theta = atan(y/z);
  }
}

/* coordinate transo from  spherical shells to cartesian */ 
void convert_shells_to_box(int b, 
                           double r, double phi, double theta,
                           double *x, double *y, double *z)
{
  // put everything in mathematica -> tata  
  // attention, you have to flip nominator/denominar sometimes
  double oos = 1./sqrt( 1./(cos(phi)*cos(phi)) + tan(theta)*tan(theta) );
  oos *= (b==0 || b==2 || b==4)?(-1.):(1.);
  
  if (b==0 || b==1) {
    *x   = r*oos;
    *y   = r*oos*tan(phi);
    *z   = r*oos*tan(theta);
  } else if (b==2 || b==3) {
    *x   = r*oos*tan(phi);
    *y   = r*oos;
    *z   = r*oos*tan(theta);
  } else {
    *x   = r*oos*tan(phi);
    *y   = r*oos*tan(theta);
    *z   = r*oos;
  } 
}

/* determine in which of the 6 shells the point is */
int find_shellsbox_from_xyz(double x, double y, double z)
{
  if (fabs(x) > fabs(y)) {
    if (fabs(x) > fabs(z))
      return (x<0) ? 0 : 1;
    else
      return (z<0) ? 4 : 5;
  } else {
    if (fabs(y) > fabs(z))
      return (y<0) ? 2 : 3;
    else
      return (z<0) ? 4 : 5;
  }
}

/* fisheye coordinates */
double logcosh(double x)
{
  
  if (x>500.) return 1.;
  else if (x<-500.) return -1.;
  else return log(cosh(x));
}

double compute_fisheye(double r, int del)
{
  double s = fisheye_s;
  double e = fisheye_e;
  double R = fisheye_R;
  
  /* look at jacobian.m, in order to change the transformation */
#if 0
  // this is the transformation from dennis and christinan
  // -> bad for the hamiltonconstraint at shell-box interface
  double a = (e+R*(R+sqrt(e+R*R)*s))/(e+R*(R+sqrt(e+R*R)));
  double b = (s-a)*sqrt(e);
  if (del==0) {
  double f  = a*(r-R) + b*sqrt(1. + (r-R)*(r-R)/E);
  double f0 = a*(0-R) + b*sqrt(1. + (0-R)*(0-R)/E);
  return f-f0;
} else if (del==1) {
  double df =  a + b*(r-R)/(E*sqrt(1. + (r-R)*(r-R)/E));
  return df;
}
#else
  // mth transfo
  // better for the transition phase, could bring numerical problems
  double a = (s-1.)*0.5*(1.+exp(-2.*R/e));
  double b = (s-a);
  if (del==0) {
    double f  = b*r + a*e*logcosh((r-R)/e);
    double f0 = b*0 + a*e*logcosh((0-R)/e);
    if (!finite(f-f0)) errorexit("prob with fisheye coordinates");
    return f-f0;
  } else if (del==1) {
    // this is the model !!!
    double df = a*tanh((r-R)/e) + b;
    if (!finite(df)) errorexit("prob with fisheye coordinates");
    return df;
  }
#endif
  return 0;
}










/* this is similar to make_grid_box/make_level_box */
tL *make_level_shellsbox(int l, int m, int n, int o, double xmin, double ymin, double zmin, 
                         double dx, double dy, double dz)
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
  level->shells = 1;
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
  
  /* return pointer to newly created level */
  free_grid_only(grid);
  return level;
}

/* add 6 boxes to grid
   -> add level 0 and move all other
   -> add 6 boxes to 0 level
   -> due it in such a way that this spherical coordinates 
      can be handeled like normal cartesian ones
   -> do not care about mpi... this SHOULD work
   -> set boundary again (physical and ref...), this would be wrong 
*/
void add_grid_shells(tG *grid)
{
  if (!Getv("grid", "shells")) 
    return;
  printf("Adding spherical shells\n");
  
  /* this is for shell-box interp (RP), after full (RK)time step */
  nbuf_r    = Geti("amr_nbuffer");
  nbuf_x    = Geti("amr_nbuffer");
  /* this is for shell-shell interp (sync), after mpi sync = after every sub-(RK)-step */
  nbuf_phi  = Geti("bampi_nghosts");
  /* number of total/global points in radial and angular direction */
  int m     = Geti("amr_shells_nr")+2*nbuf_phi;
  int n     = Geti("amr_shells_nphi")+2*nbuf_phi;
  
  fisheye_s = Getd("amr_shells_stretch");
  fisheye_R = Getd("amr_shells_r0");
  fisheye_e = Getd("amr_shells_eps");

  /* nxyz is set to 0 in grid.c -> use nx */
  if ( Geti("amr_shells_nr")==-1 ) {
    m = Geti("nx");
    if (PR) printf("set amr_shells_nr to %d\n",m);
    Seti("amr_shells_nr",m);
  }
  if ( Geti("amr_shells_nphi")==-1 ) {
    n = Geti("nx");
    if (PR) printf("set amr_shells_nphi to %d\n",n);
    Seti("amr_shells_nphi",m);
  }
  if (n%2==1 || n==0) 
    //errorexit("use even number of angle points");
    printf("You are using an odd number of angle points -> THIS IS NOT TESTED\n");
  if (grid->lbo<1)
    errorexit("you have to use berger oliger at the shell level!!!\n  -> set  amr_bo_dxmax <= dxyz/2 (there are also other ways)");
  if (Geti("amr_move_lcube")<0)
    errorexit("you have to increase amr_move_lcube to at least 1, since level1 have to be cubic for shells");
  int o = n; //important!!! do not change
  
  
  
  
  if (PR) printf("  Deleting level 0 -> no box on level 0\n"); 
  
  free_level(grid->level[grid->lmin]);
  
  
  
  if (PR) printf("  Adding 6 boxes to the new level 0\n");
  
  double r0   = grid->level[1]->bbox[1];
  double dx   = grid->level[1]->dx;
  double dy   = PI/2./(n-1-2*nbuf_phi);
  double dz   = PI/2./(o-1-2*nbuf_phi);
  double xmin = r0 - (2.*nbuf_r+0.5*Getv("amr_shells_stagger","no"))*dx;   //FIXME: use a better aproach for this number... THINK!!! don't only try
  double ymin = -PI/4. - nbuf_phi*dy;
  double zmin = -PI/4. - nbuf_phi*dz;
  
  /* prepare to determine a NEW parallelization of the shells */  
  Seti("bampi_xsize",bampi_size());
  Seti("bampi_ysize",1);
  Seti("bampi_zsize",1);
  
  int nbox,i,l;
  tL* newlevel,*newlevel0;
  for (nbox = 0; nbox < 6; nbox++) {

    /* create new level consisting of the box we want */
    newlevel0 = 
        make_level_shellsbox(0, m, n, o, xmin, ymin, zmin, dx, dy, dz);
    newlevel0->l = grid->lmin;
    newlevel0->grid = grid;

    /* save global bounding box (the local bbox is already stored in com) */
    for (i = 0; i < 6; i++) {
      newlevel0->box[0]->bbox[i] = newlevel0->bbox[i];
      newlevel0->box[0]->ibbox[i] = newlevel0->ibbox[i];
    }

    /* if this is the first box, insert new level into grid */
    if (nbox == 0) {
      newlevel = newlevel0;
      replace_level(grid, newlevel, 0);
    }

    /* if this isn't the first box, add box to existing level, discard level */
    else {
      level_add_box_withdata(newlevel, newlevel0);
      newlevel0->box[0] = 0;
      newlevel0->com = 0;
      free_level(newlevel0);
    }
  }
  newlevel->dt = grid->level[grid->lmin+1]->dt;
  newlevel->shells = 1;
  
  /* set parallelization to 0, so bam will recomute setup for later levels */
  Seti("bampi_xsize",0);
  Seti("bampi_ysize",0);
  Seti("bampi_zsize",0);
  
  
  if (PR) printf("  Set coordinates on level 0\n");
  
  int var_x = Ind("x");
  int var_y = Ind("y");
  int var_z = Ind("z");
  int var_r = Ind("shells_r");
  int var_R = Ind("shells_R");
  int var_p = Ind("shells_phi");
  int var_t = Ind("shells_theta");
  int var_rg= Ind("flagregrid");
  int var_rs= Ind("flagrestrict");
  int var_pl= Ind("flagprolong");
  enablevar(newlevel, var_x);
  enablevar(newlevel, var_y);
  enablevar(newlevel, var_z);
  enablevar(newlevel, var_r);
  enablevar(newlevel, var_R);
  enablevar(newlevel, var_p);
  enablevar(newlevel, var_t);

  /* fill in coordinates of points/nodes */
  double r,R,phi,theta,x,y,z;
  forallpoints_ijk(newlevel) {

    /* coordinates with transformation IMPORTANT and 
    THESE should be used everywhere */
    r     = box->com->bbox[0] + i*dx;
    phi   = box->com->bbox[2] + j*dy;
    theta = box->com->bbox[4] + k*dz;
    
    // fisheycoordinates
    R     = compute_fisheye(r,0);
     
    /* coordinate transformation */
    convert_shells_to_box(nbox, R,phi,theta, &x,&y,&z);
    
    /* store coordinates */
    box->level->v[var_x][ijk] = x;
    box->level->v[var_y][ijk] = y;
    box->level->v[var_z][ijk] = z;
    box->level->v[var_r][ijk] = r;
    box->level->v[var_R][ijk] = R;
    box->level->v[var_p][ijk] = phi;
    box->level->v[var_t][ijk] = theta;
    
    /* no restrict prolong flags are needed -> consistently set to 0 */
    box->level->v[var_rg][ijk] = 0;
    box->level->v[var_rs][ijk] = 0;
    box->level->v[var_pl][ijk] = 0;
    
  } endfor_ijk;
  
  
  
  if (PR) printf("  Set boundaries on level 0\n");
  
  /* locate physical boundary and set mpi boundary*/
  forallpoints_ijk(newlevel) {
    
    r = newlevel->v[var_r][ijk];
    
    /* normaly the bampibuffers in angular direction should
       be handles as such, since shell-shell is done within synchronisation, 
       till now, they are like amr buffers (for RP) and cost more 
       RHS computation */
    if (0) {
      if (i == 0 || i == box->m - 1 ||
          j == 0 || j == box->n - 1 ||
          k == 0 || k == box->o - 1)
        newlevel->boundary[ijk] = REFBOUND;
    } else {
      if (i == 0 || i > box->m - 1   ||
          j < nbuf_phi || j > box->n - 1 - nbuf_phi  ||
          k < nbuf_phi || k > box->o - 1 - nbuf_phi )
        newlevel->boundary[ijk] = REFBOUND;
    }
    
    if (dequal(box->bbox[1],r) || dless(box->bbox[1],r))
      newlevel->boundary[ijk] = PHYBOUND; 
    
  } endfor_ijk;
  
  
  forallboxes(newlevel) {
    
    box->bflag[0] = REFBOUND;
    box->bflag[1] = PHYBOUND;
    box->bflag[2] = REFBOUND;
    box->bflag[3] = REFBOUND;
    box->bflag[4] = REFBOUND;
    box->bflag[5] = REFBOUND;
    
    if (0) {
      printf("set boundary flag for box:  %2d %2d  %2d %2d  %2d %2d\n",
             box->bflag[0], box->bflag[1], box->bflag[2], 
             box->bflag[3], box->bflag[4], box->bflag[5]);
    }
    
  } endforboxes;
  
  
  /* decouple shells from boxes  TESTVERSION*/
  if ((Getv("amr_shells_RP","no"))) {
    forallpoints_ijk(newlevel) {
      if (dless( Ptr(newlevel,"shells_r")[ijk] , newlevel->bbox[0]))
        newlevel->boundary[ijk] = PHYBOUND;
    } endfor_ijk;
  }
  
  
  
  /* show the new grid setup */
  bampi_barrier();
  if (1) {
    printf("\nShells are added, here the new gridsetup:\n");
    if (0) {
      for (l=grid->lmin; l<=grid->lmax; l++) 
        printlevel(grid->level[l]);
    }
    printgrid(grid);
    
    /* here we define the outer boundary for shells, this SHOULD be 
       similar to the location used in other files */
    printf("\nOuter boundary (shells) should be at:\n");
    printf("  Rb = %f\n",newlevel->bbox[1]-Geti("boundary_N_extrapolate")*newlevel->dx);
  }
  prdivider(0);
}









/* find points (ijk or xyz) in one shell box at one of the surfaces */
inline int find_shell_local_ijk(tB* box, int ii,int jj,int kk, int a, int *ijk, double *pos)
{
  // a==0: phi = [-pi/4 .. +pi/4], theta = -pi/4   
  // a==1: phi = +pi/4,            theta = [-pi/4 .. +pi/4]
  // a==2: phi = [-pi/4 .. +pi/4], theta = +pi/4
  // a==3: phi = -pi/4,            theta = [-pi/4 .. +pi/4]
  // kk is moving, jj is constant
  tL* level=box->level;
  double dr,dp,dt, r,p,t;
  
  if (a==0) {
    // phi = [-pi/4 .. pi/4], theta = -pi/4   
    // -> phi == kk, theta == jj
    dr = level->dx;
    dp = level->dz;
    dt = level->dy;
    
    r = box->bbox[0] + ii*dr;
    p = box->bbox[2] + kk*dp;
    t = box->bbox[4] + jj*dt;
  } else if (a==1) {
    dr = level->dx;
    dp = level->dy;
    dt = level->dz;
    
    r = box->bbox[0] + ii*dr;
    p = box->bbox[3] - jj*dp;
    t = box->bbox[4] + kk*dt;
  } else if (a==2) {
    dr = level->dx;
    dp = level->dz;
    dt = level->dy;
    
    r = box->bbox[0] + ii*dr;
    p = box->bbox[3] - kk*dp;
    t = box->bbox[5] - jj*dt;
  } else if (a==3) {
    dr = level->dx;
    dp = level->dy;
    dt = level->dz;
    
    r = box->bbox[0] + ii*dr;
    p = box->bbox[2] + jj*dp;
    t = box->bbox[5] - kk*dt;
  }
  
  int i = floor((r - box->com->bbox[0])/dr + 0.001);
  int j = floor((p - box->com->bbox[2])/dp + 0.001);
  int k = floor((t - box->com->bbox[4])/dt + 0.001);
  
  if (i<box->com->ibbox[0] || i>box->com->ibbox[1] || 
      j<box->com->ibbox[2] || j>box->com->ibbox[3] || 
      k<box->com->ibbox[4] || k>box->com->ibbox[5] )
    return 0;
  
  int m = box->com->ibbox[1]-box->com->ibbox[0]+1;
  int n = box->com->ibbox[3]-box->com->ibbox[2]+1;
  int o = box->com->ibbox[5]-box->com->ibbox[4]+1;
  *ijk = box->noffset + m*n*(k) + m*(j) + i;
  
  pos[0] = r;
  pos[1] = p;
  pos[2] = t;
  
  return 1;
}

// force instantiation of inline function (in case we compile with -g)
int find_shell_local_ijk(tB* box, int ii,int jj,int kk, int a, int *ijk, double *pos);

/* special interpolation of two arrays */
inline void interpolate_1D(double *x1,double *d1, double *x2,double *d2, int N,int order,int gh,int layer, int nv,int*iv)
{
  /* both fields have N pts
  we have to interpolate on x2 
  and fill x1
  x1 has to filled at [Ngh .. N-Ngh-1] because of the edges (make a picture and you will see)
  x2 has to be used only inside -> [order/2 .. N-order/2-1]
  
  x2 is always bigger than x1
  */
  
  int i1,i2, j;
  int m = order - 1;  // order 4 and 6:  3 and 5
  int n = m/2;        // order 4 and 6:  1 and 2
  double c[order];
  double h = x2[1]-x2[0]; //normal angular direction ... h = dphi = dtheta
  int nn;
  
  if (0) {
    printf("%d\n",gh);
    for (j=0; j<N; j++) 
      printf("  %d    x1=%e   x2=%e   \n",j,-x1[j],x2[j]);
  }
  
  i2 = gh;
  /* go through all x1 points */
  for (i1=0; i1<N; i1++) {
  
    // do not interpolate in bufferzones
    if ((i1<layer) || (i1>=N-layer)) {
      d1[nv*N+i1] = 0.; 
      continue;
    }
    
    // set flag because here we can interpolate
    d1[nv*N+i1] = 1.;
    
    // find location for the interpolation stencil
    // i2 = [layer .. N-layer-1]
    // FIXME: OPTIMIZE !!
    //i2=gh;
    if (x1[0]<x1[1]) {
      while (x2[i2+n+1]<x1[i1] && i2<N-gh-m-1)
        i2++;
    } else {
      while (x2[i2+n+1]<-x1[i1] && i2<N-gh-m-1)
        i2++;
    }
    
    // test stencil
    if (0) {
      printf("%d (%d)   ",i1,N);
      for (j=i2; j<i2+order; j++)
        printf("  %2.4e",x2[j]);
      printf("    <-> %2.4e\n",x1[i1]);
    }
    
    // now interpolate and compute coefficiants
    if (x1[0]<x1[1])
      coefficients_lagrange_N(order,  x1[i1], x2[i2],h, c);
    else 
      coefficients_lagrange_N(order, -x1[i1], x2[i2],h, c);
    
    if (0) {
      printf(" %d   %e %e %e",order,x1[i1],x2[i2],h);
      for (j=0; j<order; j++)
        printf(" %e",c[j]);
      printf("\n\n");
    }
    
    // sum values for each variable
    for (nn=0; nn<nv; nn++) {
      
      d1[nn*N+i1] = 0.;
      for (j=0; j<order; j++)
        d1[nn*N+i1] += c[j] * d2[nn*N+i2+j];
      
      // test interpolation
      if (0) {
        printf(" %e %e\n", x1[i1],d1[i1]);
        for (j=0; j<order; j++)
          printf("      %e %e\n",x2[i2+j],d2[i2+j]);
        printf("\n\n");
      }
      
    }
  }
  
}

// force instantiation of inline function (in case we compile with -g)
void interpolate_1D(double *x1,double *d1, double *x2,double *d2, int N,int order,int gh,int layer, int nv,int*iv);

/* take ghostpoints from box1 and store in global buffer
   find location in box2
   find stencil in box2
   interpolate with values of box2
   sync buffer
   put buffer values back to box1
*/
void sync_shells_surface_interpolate(tB* box1, tB* box2, int a1, int a2, double f, int nv, int *iv, int order)
{
  int i,j, ii,jj,kk, ijk, n;
  int gh    = nbuf_phi;
  int flag,nr;
  double *buffer,*global,*vinterp;
  double pos[3];
  
  /* assume phi and theta are same size */
  int pts = box1->ibbox[5]-box1->ibbox[4]+1;
  int npoints = (box1->ibbox[1]-box1->ibbox[0]+1)*pts*gh;
  double x1[pts],x2[gh][pts], data2[pts],h2;
  double dx = box1->dy;
  
  /* test if these surfaces fit together, does only work for bam on one processor */
  if (0 && bampi_size()==1) {
    double x,y,z;
    
    printf("box%d\n",box1->i);
    
    find_shell_local_ijk(box1, 0,gh,pts-1-gh, a1, &ijk, pos);
    convert_shells_to_box(box1->i, pos[0],pos[1],pos[2], &x,&y,&z);
    printf("  %+2.4f %+2.4f %+2.4f    %+2.4f %+2.4f %+2.4f\n",pos[0],pos[1],pos[2],x,y,z);
    find_shell_local_ijk(box1, box1->ibbox[1]-box1->ibbox[0],gh,pts-1-gh, a1, &ijk, pos);
    convert_shells_to_box(box1->i, pos[0],pos[1],pos[2], &x,&y,&z);
    printf("  %+2.4f %+2.4f %+2.4f    %+2.4f %+2.4f %+2.4f\n",pos[0],pos[1],pos[2],x,y,z);
    find_shell_local_ijk(box1, box1->ibbox[1]-box1->ibbox[0],gh,gh, a1, &ijk, pos);
    convert_shells_to_box(box1->i, pos[0],pos[1],pos[2], &x,&y,&z);
    printf("  %+2.4f %+2.4f %+2.4f    %+2.4f %+2.4f %+2.4f\n",pos[0],pos[1],pos[2],x,y,z);
    find_shell_local_ijk(box1, 0,gh,gh, a1, &ijk, pos);
    convert_shells_to_box(box1->i, pos[0],pos[1],pos[2], &x,&y,&z);
    printf("  %+2.4f %+2.4f %+2.4f    %+2.4f %+2.4f %+2.4f\n",pos[0],pos[1],pos[2],x,y,z);
    
    printf("box%d\n",box2->i);
    find_shell_local_ijk(box2, 0,gh,pts-1-gh, a2, &ijk, pos);
    convert_shells_to_box(box2->i, pos[0],pos[1],pos[2], &x,&y,&z);
    printf("  %+2.4f %+2.4f %+2.4f    %+2.4f %+2.4f %+2.4f\n",pos[0],pos[1],pos[2],x,y,z);
    find_shell_local_ijk(box2, box2->ibbox[1]-box2->ibbox[0],gh,pts-1-gh, a2, &ijk, pos);
    convert_shells_to_box(box2->i, pos[0],pos[1],pos[2], &x,&y,&z);
    printf("  %+2.4f %+2.4f %+2.4f    %+2.4f %+2.4f %+2.4f\n",pos[0],pos[1],pos[2],x,y,z);
    find_shell_local_ijk(box2, box2->ibbox[1]-box2->ibbox[0],gh,gh, a2, &ijk, pos);
    convert_shells_to_box(box2->i, pos[0],pos[1],pos[2], &x,&y,&z);
    printf("  %+2.4f %+2.4f %+2.4f    %+2.4f %+2.4f %+2.4f\n",pos[0],pos[1],pos[2],x,y,z);
    find_shell_local_ijk(box2, 0,gh,gh, a2, &ijk, pos);
    convert_shells_to_box(box2->i, pos[0],pos[1],pos[2], &x,&y,&z);
    printf("  %+2.4f %+2.4f %+2.4f    %+2.4f %+2.4f %+2.4f\n",pos[0],pos[1],pos[2],x,y,z);
   
    printf("\n\n");
  }
  
  /* compare both 1D coordinate systems */
  //FIXME: do a static array!!!
  for (i=0; i<pts; i++) {
    x1[i]  = box1->bbox[4] +i*dx;
    for (j=0; j<gh; j++) 
      x2[j][i]  = atan( tan(x1[i]) / tan( ( box1->bbox[4] +j*dx) ) );
  }
  
  /* go through all expected global poins */
  for (ii=box1->ibbox[0]; ii<=box1->ibbox[1];  ii++)
  for (jj=0; jj<gh;  jj++) {
    
    /* do memory managment */
    buffer  = (double*) dcalloc ((nv+1)*pts);
    vinterp = (double*) dcalloc ((nv+1)*pts);
    
    //FIXME test if this is necessary, maybe by thinking abount reusing array instead of free+alloc
    for (kk=0; kk<pts; kk++)
      for (n=0; n<=nv; n++)
        buffer[n*pts + kk] = vinterp[n*pts + kk] = 0.;
    
    /* find data in box2 and store it */
    for (kk=0; kk<pts; kk++) {
      flag = find_shell_local_ijk(box2, ii,2*gh-jj,kk, a2, &ijk, pos);
      
      if (flag) {
        for (n=0; n<nv; n++) 
          buffer[n*pts+kk] =  box2->level->v[iv[n]][ijk];
        buffer[nv*pts+kk] = 1.;
      } else {
        for (n=0; n<=nv; n++)
          buffer[n*pts+kk] = 0.;
      }
    }
    
    /* sync points -> 1D data-stream is on all processors */
    if (bampi_size() > 1 ) {
      global = (double*) malloc ((nv+1)*pts*sizeof(double));
      bampi_allreduce_sum_vector(buffer, global, (nv+1)*pts);
      for (i=0; i<pts; i++) {
        if (global[nv*pts+i]==0.) {
          errorexit("shell->shell interpolate: this point is nowhere");
        } else { //OPTIMIZE
          for (n=0; n<nv; n++)
            global[n*pts+i] /= global[nv*pts+i];
        }
      }
      free(buffer);
      buffer = global;
    }
    
    /* test it for one variable */
    if (0) {
      for (kk=0; kk<pts; kk++) {
        printf("%d   %e %e    |",kk,x2[jj][kk],buffer[kk]);
        
        flag = find_shell_local_ijk(box1, ii,jj,kk, a1, &i, pos);
        if (flag)
          printf("  %e %e",x1[kk],box1->level->v[iv[0]][i]);
        printf("\n");
      }
    }
    
    /* interpolate 1D data-stream locally, attention there are bufferpoints */
    /* note: due to symmetry I switch the fields -> now the field is equdistant */
    interpolate_1D(x2[jj],vinterp, x1,buffer, pts,order,gh,jj, nv,iv);
    
    /* copy interpolated data to field */
    for (kk=0; kk<pts; kk++) {
      
      flag = find_shell_local_ijk(box1, ii,jj,kk, a1, &ijk, pos);
      
      /* attention we have to say whether the angles are co or conter rotating */
      if (flag && vinterp[nv*pts+kk]==1. )
        for (n=0; n<nv; n++) {
         
          // version 1
          //box1->level->v[iv[n]][ijk] = vinterp[ n*pts+ (f==1) ? (kk) : (pts-1-kk) ];
          
          // version 2 -> Doreen's bugfix 
          box1->level->v[iv[n]][ijk] = vinterp[ n*pts+ ((f==1) ? (kk) : (pts-1-kk)) ];
          
          // version 3
          /*nr = n*pts;
          if (f==1)
            nr += kk;
          else
            nr += pts-1-kk;
          
          box1->level->v[iv[n]][ijk] = vinterp[ nr ];
          */
          
        }
     
    }

    
    free(buffer);
    free(vinterp);
    
  }
  
}

/* efficient version if we have optimal paralellization (n x 1 x 1) */
void sync_shells_surface_interpolate_opt(tB* box1, tB* box2, int a1, int a2, double f, int nv, int *iv, int order)
{
  
  //bampi_openmp_start

  int i,j, ii,jj,kk, ii0, ijk, n;
  int gh    = nbuf_phi;
  int flag,nr;
  double pos[3], h2;
  
  /* assume phi and theta are same size */
  int pts = box1->com->ibbox[5]-box1->com->ibbox[4]+1;
  int npoints = (box1->com->ibbox[1]-box1->com->ibbox[0]+1)*pts*gh;
  //double x1[pts],x2[gh][pts], data2[pts],h2;
  double dx = box1->dy;
  
  /* do memory managment */
  double *buffer  = (double*) dcalloc ((nv+1)*pts);
  double *vinterp = (double*) dcalloc ((nv+1)*pts);
  double *x1      = (double*) dcalloc (pts);
  double *data    = (double*) dcalloc (pts);
  double **x2     = (double**) malloc (gh*sizeof(double*));
  for (i=0; i<gh;  i++)
    x2[i] = (double*) dcalloc (pts);
  
  
  /* compare both 1D coordinate systems */
  for (i=0; i<pts; i++) {
    x1[i]  = box1->com->bbox[4] +i*dx;
    for (j=0; j<gh; j++) 
      x2[j][i]  = atan( tan(x1[i]) / tan( ( box1->com->bbox[4] +j*dx) ) );
  }
  
  /* go through all expected global poins */
  //bampi_openmp_loop
  for (ii=box1->com->ibbox[0]; ii<=box1->com->ibbox[1];  ii++) {
    for (jj=0; jj<gh;  jj++) {
      
      /* set buffers to 0 */
      for (kk=0; kk<pts; kk++)
        for (n=0; n<=nv; n++)
          buffer[n*pts + kk] = vinterp[n*pts + kk] = 0.;
      
      /* set kk0 because we have an offset due to radial parallelization */
      ii0 = (int)( floor ((box1->com->bbox[0]-box1->bbox[0])/box1->dx +.0001 ) );
      
      /* find data in box2 and store it */
      for (kk=0; kk<pts; kk++) {
        flag = find_shell_local_ijk(box2, ii0+ii,2*gh-jj,kk, a2, &ijk, pos);
        
        if (flag) {
          for (n=0; n<nv; n++) 
            buffer[n*pts+kk] =  box2->level->v[iv[n]][ijk];
          buffer[nv*pts+kk] = 1.;
        } else {
          for (n=0; n<=nv; n++)
            buffer[n*pts+kk] = 0.;
        }
      }
      
      /* interpolate 1D data-stream locally, attention there are bufferpoints */
      /* note: due to symmetry I switch the fields -> now the field is equdistant */
      interpolate_1D(x2[jj],vinterp, x1,buffer, pts,order,gh,jj, nv,iv);
      
      /* copy interpolated data to field */
      for (kk=0; kk<pts; kk++) {
        flag = find_shell_local_ijk(box1, ii0+ii,jj,kk, a1, &ijk, pos);
        
        /* attention we have to say whether the angles are co or conter rotating */
        if (flag && vinterp[nv*pts+kk]==1.) {
          for (n=0; n<nv; n++)
            box1->level->v[iv[n]][ijk] = vinterp[ n*pts+ ((f==1) ? (kk) : (pts-1-kk)) ];
        }
      }
      
    }
  }
  
  free(buffer);
  free(vinterp);
  free(x1);
  free(data);
  for (i=0; i<gh; i++)
    free(x2[i]);
  free(x2);
  
  //bampi_openmp_stop
}

/* find out if we can use optimized version (only one proc in phi and theta direction)
   => every processor can sync with itself, not mpi needed */
int shells_optimize(tL* level)
{
  int optimize = Getv("amr_shells_optimize","yes");
  forallboxes(level) {
    if ((box->com->sizexyz[1]!=1) || (box->com->sizexyz[2]!=1))
      optimize = 0;
  } endforboxes;
  
  return optimize;
}

/* sync shells with themselfes
   this is slow and stupid but we do not need it very fast
   this can not be handeld by the standart bampi stuff ...
   at least I do not see how (mth 04/11)
*/
void sync_shells(tL* level0, int nv, int *iv)
{
  if (level0->l!=level0->grid->lmin)
    return;
  if (!Getv("grid","shells"))
    return;
  if (!nv)
    return;
  
  if (PR) printf("interpolate between shells\n");
  timer_start(level0, "sync_shells");
  
  int n,i,j;
  int order = Geti("order_RP_shells");
  
  if (PR) {
    printf("Variables to sync (shells <-> shells):\n");
    for (n=0; n<nv; n++)
      printf("  %s\n", VarName(iv[n]));
  }
  
  
  /* these numbers (nr 3,4,5) are adjusted by eye... 
     I took visit and tried until all patches fit together
     this took very long -> don't change them */
  if (shells_optimize(level0)) {
    if (PR) 
      printf("  optimized version (no mpi communication is needed)\n");
    
    //sync box0 ... -x
    sync_shells_surface_interpolate_opt(level0->box[0],level0->box[2],1,1, 1., nv,iv, order); // -y
    sync_shells_surface_interpolate_opt(level0->box[0],level0->box[3],3,3,-1., nv,iv, order); // +y
    sync_shells_surface_interpolate_opt(level0->box[0],level0->box[4],2,1,-1., nv,iv, order); // -z
    sync_shells_surface_interpolate_opt(level0->box[0],level0->box[5],0,3, 1., nv,iv, order); // +z
    
    //sync box1 ... +x
    sync_shells_surface_interpolate_opt(level0->box[1],level0->box[2],3,3,-1., nv,iv, order); // -y
    sync_shells_surface_interpolate_opt(level0->box[1],level0->box[3],1,1, 1., nv,iv, order); // +y
    sync_shells_surface_interpolate_opt(level0->box[1],level0->box[4],0,3, 1., nv,iv, order); // -z
    sync_shells_surface_interpolate_opt(level0->box[1],level0->box[5],2,1,-1., nv,iv, order); // +z
  
    //sync box2 ... -y
    sync_shells_surface_interpolate_opt(level0->box[2],level0->box[0],1,1, 1., nv,iv, order); // -y
    sync_shells_surface_interpolate_opt(level0->box[2],level0->box[1],3,3,-1., nv,iv, order); // +y
    sync_shells_surface_interpolate_opt(level0->box[2],level0->box[4],2,2, 1., nv,iv, order); // -z
    sync_shells_surface_interpolate_opt(level0->box[2],level0->box[5],0,0,-1., nv,iv, order); // +z
    
    //sync box3 ... +y
    sync_shells_surface_interpolate_opt(level0->box[3],level0->box[0],3,3,-1., nv,iv, order); // -y
    sync_shells_surface_interpolate_opt(level0->box[3],level0->box[1],1,1, 1., nv,iv, order); // +y
    sync_shells_surface_interpolate_opt(level0->box[3],level0->box[4],0,0,-1., nv,iv, order); // -z
    sync_shells_surface_interpolate_opt(level0->box[3],level0->box[5],2,2, 1., nv,iv, order); // +z
    
    //sync box4 ... -z
    sync_shells_surface_interpolate_opt(level0->box[4],level0->box[0],1,2,-1., nv,iv, order); // -y
    sync_shells_surface_interpolate_opt(level0->box[4],level0->box[1],3,0, 1., nv,iv, order); // +y
    sync_shells_surface_interpolate_opt(level0->box[4],level0->box[2],2,2, 1., nv,iv, order); // -z
    sync_shells_surface_interpolate_opt(level0->box[4],level0->box[3],0,0,-1., nv,iv, order); // +z
    
    //sync box5 ... +z
    sync_shells_surface_interpolate_opt(level0->box[5],level0->box[0],3,0, 1., nv,iv, order); // -y
    sync_shells_surface_interpolate_opt(level0->box[5],level0->box[1],1,2,-1., nv,iv, order); // +y
    sync_shells_surface_interpolate_opt(level0->box[5],level0->box[2],0,0,-1., nv,iv, order); // -z
    sync_shells_surface_interpolate_opt(level0->box[5],level0->box[3],2,2, 1., nv,iv, order); // +z
    
    //now sync 
    /* normaly  bampi_vlsynchronize(u); is correct, hewever it would call itself ... 
       this is bad, therefore we do it like inside bampi_vlsynchronize(): */
    int nbox;
    tL *obl;
    for (nbox = 0; nbox < 6; nbox++) {
      obl = one_box_level(level0, nbox);
      bampi_synchronize_variables(obl, nv, iv);
    }
    
    
    
  } else {
    //sync box0 ... -x
    sync_shells_surface_interpolate(level0->box[0],level0->box[2],1,1, 1., nv,iv, order); // -y
    sync_shells_surface_interpolate(level0->box[0],level0->box[3],3,3,-1., nv,iv, order); // +y
    sync_shells_surface_interpolate(level0->box[0],level0->box[4],2,1,-1., nv,iv, order); // -z
    sync_shells_surface_interpolate(level0->box[0],level0->box[5],0,3, 1., nv,iv, order); // +z
    
    //sync box1 ... +x
    sync_shells_surface_interpolate(level0->box[1],level0->box[2],3,3,-1., nv,iv, order); // -y
    sync_shells_surface_interpolate(level0->box[1],level0->box[3],1,1, 1., nv,iv, order); // +y
    sync_shells_surface_interpolate(level0->box[1],level0->box[4],0,3, 1., nv,iv, order); // -z
    sync_shells_surface_interpolate(level0->box[1],level0->box[5],2,1,-1., nv,iv, order); // +z
  
    //sync box2 ... -y
    sync_shells_surface_interpolate(level0->box[2],level0->box[0],1,1, 1., nv,iv, order); // -y
    sync_shells_surface_interpolate(level0->box[2],level0->box[1],3,3,-1., nv,iv, order); // +y
    sync_shells_surface_interpolate(level0->box[2],level0->box[4],2,2, 1., nv,iv, order); // -z
    sync_shells_surface_interpolate(level0->box[2],level0->box[5],0,0,-1., nv,iv, order); // +z
    
    //sync box3 ... +y
    sync_shells_surface_interpolate(level0->box[3],level0->box[0],3,3,-1., nv,iv, order); // -y
    sync_shells_surface_interpolate(level0->box[3],level0->box[1],1,1, 1., nv,iv, order); // +y
    sync_shells_surface_interpolate(level0->box[3],level0->box[4],0,0,-1., nv,iv, order); // -z
    sync_shells_surface_interpolate(level0->box[3],level0->box[5],2,2, 1., nv,iv, order); // +z
    
    //sync box4 ... -z
    sync_shells_surface_interpolate(level0->box[4],level0->box[0],1,2,-1., nv,iv, order); // -y
    sync_shells_surface_interpolate(level0->box[4],level0->box[1],3,0, 1., nv,iv, order); // +y
    sync_shells_surface_interpolate(level0->box[4],level0->box[2],2,2, 1., nv,iv, order); // -z
    sync_shells_surface_interpolate(level0->box[4],level0->box[3],0,0,-1., nv,iv, order); // +z
    
    //sync box5 ... +z
    sync_shells_surface_interpolate(level0->box[5],level0->box[0],3,0, 1., nv,iv, order); // -y
    sync_shells_surface_interpolate(level0->box[5],level0->box[1],1,2,-1., nv,iv, order); // +y
    sync_shells_surface_interpolate(level0->box[5],level0->box[2],0,0,-1., nv,iv, order); // -z
    sync_shells_surface_interpolate(level0->box[5],level0->box[3],2,2, 1., nv,iv, order); // +z
    
  }
  
  if (PR) printf("interpolate between shells ... finished\n");
  bampi_barrier();
  timer_stop(level0, "sync_shells");
}



/* interpolation function in one direction */
void interpolate_box_to_shells(tL* level0, tL* level1, int nv, int *iv)
{
  
  /* this was a first "guess" which worked, I do not know why */
  int mm = (int)( (sqrt(3)*level1->bbox[1] - sqrt(3)*nbuf_x*level1->dx - level0->bbox[0])/(level0->dx) +1.5);
  int nn = level0->ibbox[3]-level0->ibbox[2]+1-2*nbuf_phi;
  int oo = level0->ibbox[5]-level0->ibbox[4]+1-2*nbuf_phi;
  int npoints = mm*nn*oo;
  double *buffer = (double*) dcalloc ((nv+1)*6*npoints);
  
  
  
  /* parallelize ONLY the interpolation, not the bampi-sync nor copying */
  bampi_openmp_start
  
  int order = Geti("order_RP_shells");
  
  int flag;
  int ii,jj,kk, iijjkk,n;
  double r,phi,theta,x,y,z, f;
  double vinterp[nv];

  /* find all points inside the shell */
  forallboxes(level0) {
    
    /* go through all points in one shell */ 
    bampi_openmp_loop
    for (kk=0; kk<oo; kk++) {
      for (jj=0; jj<nn; jj++) {
        for (ii=0; ii<mm; ii++) {
          
          iijjkk = nbox*npoints + kk*nn*mm+jj*mm+ii;
      
          r     = level0->bbox[0] + ii*level0->dx;
          phi   = level0->bbox[2] + (jj+nbuf_phi)*level0->dy;
          theta = level0->bbox[4] + (kk+nbuf_phi)*level0->dz;
          
          // fisheycoordinates
          r     = compute_fisheye(r,0);
          
          convert_shells_to_box(nbox, r,phi,theta, &x,&y,&z);
          
          /* test if xyz is inside prozessor */
          flag = check_interpolation_cube_local_withsym(level1,x,y,z, order);
          
          /* test if we have enough points to interpolate */
          if (flag) {
            if (nbox==0 && dless(x,level1->bbox[0]+(nbuf_x-1)*level1->dx)) flag = 0;
            if (nbox==1 && dless(level1->bbox[1]-(nbuf_x-1)*level1->dx,x)) flag = 0;
            if (nbox==2 && dless(y,level1->bbox[2]+(nbuf_x-1)*level1->dy)) flag = 0;
            if (nbox==3 && dless(level1->bbox[3]-(nbuf_x-1)*level1->dy,y)) flag = 0;
            if (nbox==4 && dless(z,level1->bbox[4]+(nbuf_x-1)*level1->dz)) flag = 0;
            if (nbox==5 && dless(level1->bbox[5]-(nbuf_x-1)*level1->dz,z)) flag = 0;
          }
          
          /* if flag interpolate with symmetry */
          if (flag)
            flag = interpolate_xyz_local_minimal_withsym(level1, x,y,z, nv,iv, vinterp, order,LAGRANGE);
          
          /* if flag, store value */
          if (flag) {
            for (n=0; n<nv; n++)
              buffer[n*6*npoints+iijjkk] = vinterp[n];
            buffer[nv*6*npoints+iijjkk] = 1.;
          } else {
            for (n=0; n<=nv; n++)
              buffer[n*6*npoints+iijjkk] = 0.;
          }
          
        }
      }
    }
    
  } endforboxes;
  
  bampi_openmp_stop
  
  
  
  /* sync buffer -> all interpolated values are then available on all processors */
  int i,n;
  if (bampi_size() > 1) {
    double *global = (double*) malloc ((nv+1)*6*npoints*sizeof(double));
    bampi_allreduce_sum_vector(buffer, global, (nv+1)*6*npoints);
    for (i=0; i<6*npoints; i++) {
      if (global[nv*6*npoints+i]==0) {
        /* ATTENTION: this can happen because due to the round shape of the shells
        there are points inside the shell radial  buffer which are not inside the 
        box -> we want to interpolate at more points */
        //errorexit("shell->box interpolate: this point is nowhere");
      } else {//OPTIMIZE
        for (n=0; n<nv; n++)
          global[n*6*npoints+i] /= global[nv*6*npoints+i];
      }
    }
    free(buffer);
    buffer = global;
  }
  
  
  
  /* put the values at the correct local position */
  int ii,jj,kk,iijjkk;
  forallboxes(level0) {
    forallpoints_boxijk(box) {
    
      ii = i + floor( (box->com->bbox[0]-box->bbox[0])/box->dx+.0001 );
      jj = j + floor( (box->com->bbox[2]-box->bbox[2])/box->dy+.0001 )-nbuf_phi;
      kk = k + floor( (box->com->bbox[4]-box->bbox[4])/box->dz+.0001 )-nbuf_phi; 
      
      /* set value at correct location */
      if (ii>=0 && jj>=0 && kk>=0 && ii<mm && jj<nn && kk<oo) {
        iijjkk = nbox*npoints + mm*nn*kk + mm*jj + ii;
        if (iijjkk>=6*npoints)
          errorexit("something wrong");
        
        if (buffer[nv*6*npoints+iijjkk]!=0) {
          for (n=0; n<nv; n++)
            level0->v[iv[n]][ijk] = buffer[n*6*npoints+iijjkk];
        }
      }
    
    } endfor_ijk;
    
  } endforboxes;
  
  /* free buffer, clear for next usage */
  free(buffer);
}

/* interpolation function in other direction */
void interpolate_shells_to_box(tL* level0, tL* level1, int nv, int *iv)
{
  /* set number of points */
  int i,mml[6],mmr[6],nnl[6],nnr[6],ool[6],oor[6];  
  for (i=0; i<6; i++) {
    mml[i] = nnl[i] = ool[i] = 0;
    mmr[i] = level1->ibbox[1]-level1->ibbox[0]+1; 
    nnr[i] = level1->ibbox[3]-level1->ibbox[2]+1;
    oor[i] = level1->ibbox[5]-level1->ibbox[4]+1;
  }
  mmr[0] = mml[0]+nbuf_x;
  mml[1] = mmr[1]-nbuf_x;
  nnr[2] = nnl[2]+nbuf_x;
  nnl[3] = nnr[3]-nbuf_x;
  oor[4] = ool[4]+nbuf_x;
  ool[5] = oor[5]-nbuf_x;
 
  int npoints_off[6];
  int npoints = 0;
  for (i=0; i<6; i++) {
    npoints_off[i] = npoints;
    npoints += (mmr[i]-mml[i])*(nnr[i]-nnl[i])*(oor[i]-ool[i]);
    if (0) printf("  alloc %d  +=  %d %d %d\n",npoints, (mmr[i]-mml[i]),(nnr[i]-nnl[i]),(oor[i]-ool[i]));
  }
  double *buffer = (double*) dcalloc ((nv+1)*npoints);
  
  
  
  
  /* parallelize ONLY the interpolation, not the bampi-sync nor copying */
  bampi_openmp_start
  
  int order = Geti("order_RP_shells");
  
  int flag;
  int ii,jj,kk, iijjkk, n,l;
  double r,phi,theta,x,y,z, f;
  double vinterp[nv];
  
  /* find all points inside the box */
  forallboxes(level1) {
    
    /* go through all 6 sides */
    for (l=0; l<6; l++) {
      
      /* go through global points */
      bampi_openmp_loop
      for (kk=ool[l]; kk<oor[l]; kk++) {
        for (jj=nnl[l]; jj<nnr[l]; jj++) { 
          for (ii=mml[l]; ii<mmr[l]; ii++) {
            
            iijjkk = npoints_off[l] +
                     (kk-ool[l])*(nnr[l]-nnl[l])*(mmr[l]-mml[l]) + 
                     (jj-nnl[l])*(mmr[l]-mml[l]) + (ii-mml[l]);
        
            /* handle symmetries */
            if (level1->grid->bitant && l==4) {
              buffer[nv*npoints+iijjkk] = 1.;
            } else {
              x = level1->bbox[0] + ii*level1->dx;
              y = level1->bbox[2] + jj*level1->dy;
              z = level1->bbox[4] + kk*level1->dz;
              
              // fisheycoordinates
              f  = compute_fisheye(sqrt(x*x+y*y+z*z),0) / sqrt(x*x+y*y+z*z);
              x *= f;
              y *= f;
              z *= f;
          
              /* dirty trick, use only the shell-box which is at the expected side
              -> some points could go wrong ... lets see */
              convert_box_to_shells(l, x,y,z, &r,&phi,&theta);
          
              /* test if xyz is inside prozessor */
              flag = interpolate_xyz_localinbox_minimal(level0->box[l], 
                  r,phi,theta, nv,iv,vinterp, order,LAGRANGE);
          
              /* if flag interpolate */
              if (flag) {
                for (n=0; n<nv; n++)
                  buffer[n*npoints+iijjkk] = vinterp[n];
                buffer[nv*npoints+iijjkk] = 1.;
              } else {
                for (n=0; n<=nv; n++)
                  buffer[n*npoints+iijjkk] = 0.;
              }
            }
            
          }
        }
      }
    }
    
  } endforboxes;
  
  bampi_openmp_stop
  
  
  
  
  /* sync buffer -> all interpolated values are then available on all processors */
  int l,n;
  if (bampi_size() > 1) {
    double *global = (double*) malloc ((nv+1)*npoints*sizeof(double));
    bampi_allreduce_sum_vector(buffer, global, (nv+1)*npoints);
    for (i=0; i<npoints; i++) {
      if (global[nv*npoints+i]==0) {
        errorexit("box->shell interpolate: this point is nowhere");
        //printf("box->shell interpolate: this point is nowhere  %d  %d  %e\n",i,npoints,(double)(npoints)/(double)(i));
      } else {//OPTIMIZE
        for (n=0; n<nv; n++)
          global[n*npoints+i] /= global[nv*npoints+i];
      }
    }
    free(buffer);
    buffer = global;
  }
  
  
  /* put the values at the correct local position */
  int ii,jj,kk,iijjkk,flag;
  double x,y,z;
  forallboxes(level1) {
    
    /* go through all 6 sides */
    for (l=0; l<6; l++) {
      
      /* handle symmetries */
      if (level1->grid->bitant && l==4) continue;
      
      /* go through global points */
      forallpoints_boxijk(box) {

        ii = i + floor( (box->com->bbox[0]-box->bbox[0])/box->dx+.0001 );
        jj = j + floor( (box->com->bbox[2]-box->bbox[2])/box->dy+.0001 );
        kk = k + floor( (box->com->bbox[4]-box->bbox[4])/box->dz+.0001 );
        
        if (ii>=mml[l] && jj>=nnl[l] && kk>=ool[l] && ii<mmr[l] && jj<nnr[l] && kk<oor[l]) {
          
          iijjkk  = npoints_off[l] + 
              (mmr[l]-mml[l])*(nnr[l]-nnl[l])*(kk-ool[l]) + 
              (mmr[l]-mml[l])*(jj-nnl[l]) + 
              (ii-mml[l]);
          
          if (iijjkk>=npoints) 
            errorexit("something is wrong");
          
          /* test if this point depends only to the interior of this shell */
          flag = 0;
          if (buffer[nv*npoints+iijjkk]>0.) {
            x = fabs(box->bbox[0] + ii*box->dx);
            y = fabs(box->bbox[2] + jj*box->dy);
            z = fabs(box->bbox[4] + kk*box->dz);
            if      (l<=1) flag = (x>=y-0.1*box->dy && x>=z-0.1*box->dz);
            else if (l<=3) flag = (y>=x-0.1*box->dx && y>=z-0.1*box->dz);
            else if (l<=5) flag = (z>=x-0.1*box->dx && z>=y-0.1*box->dy);
          }
          
          if (flag) {
            for (n=0; n<nv; n++)
              level1->v[iv[n]][ijk] = buffer[n*npoints+iijjkk];
          }
        }
      } endfor_ijk;
    }
    
  } endforboxes;
  
  /* free buffer again, clear for next usage */
  free(buffer);
  
}



/* sync shells with the box in the next level
   this is slow and stupid but we do not need it very fast
   this can not be handeld by the standart bampi stuff ...
   at least I do not see how (mth 04/11)
*/
void restrict_prolong_shells(tG* g, int lfine, int lcoarse, tVarList *uf, tVarList *uc)
{
  /* decouple shells from boxes  DEBUGVERSION */
  if (Getv("amr_shells_RP","no"))
    return ;
  
  /* shell[b0..5,l0] <-> box[b0,l1] */
  if (PR) 
    printf("prolong restrict sync  (special sync between shells in level 0 and box in level 1)\n");
  
  if (lfine!=g->lmin+1 || lcoarse!=g->lmin)
    errorexit("there is something wrong with sync_shells");

  
  /* find all variables to sync */
  if (!uc || !uc->n || !uf || !uf->n)
    return;
  
  tL* level0 = g->level[lcoarse];
  tL* level1 = g->level[lfine];
  int nv  = uc->n; 
  int *iv = uc->index;
  if (PR) {
    int n;
    printf("Variables to RP (shells <-> box):\n");
    for (n=0; n<nv; n++)
      printf("  %s\n", VarName(iv[n]));
  }
  
  
  
  /* ************************************************************************* */
  /* Interpolate from BOX to SHELL                                             */
  /* ************************************************************************* */
  if (PR) printf("interpolate from box to shells\n");
  interpolate_box_to_shells(level0, level1, nv, iv);
  bampi_barrier();
  
  
  
  
  
  
  /* ************************************************************************* */
  /* Interpolate between SHELLS in 1D                                          */
  /* ************************************************************************* */
  if (PR) printf("do actual shells patches sync\n");
  sync_shells(level0, nv, iv);
  
  
  
  
  
  
  /* ************************************************************************* */
  /* Set symmetries in SHELLS                                                  */
  /* ************************************************************************* */
  if (PR) printf("do symmetries in shells\n");
  set_boundary_symmetry(level0, uc);
  
  
  
  
  
  
  /* ************************************************************************* */
  /* Interpolate from SHELLS to BOX                                            */
  /* ************************************************************************* */
  if (PR) printf("interpolate from shells to box\n");
  interpolate_shells_to_box(level0, level1, nv, iv);
  bampi_barrier();
  
  
  
  
  
  
  /* ************************************************************************* */
  /* Set symmetries in BOX                                                     */
  /* ************************************************************************* */
  if (PR) printf("do symmetries in box\n");
  set_boundary_symmetry(level1, uf);
  
  
  
  
  
  /* better way for all procs */
  if (PR) printf("prolong restrict sync ... finished\n");
  bampi_barrier();

}









/* apply symmetry for shells ... this is stupid but can
   reduce/cure grid modes 
   this should be ONLY called by set_boundary_symmetry();
*/
void set_boundary_shells_symmetry(tL *level, tVarList *varlist)
{
 /* do we want to introduce artificial symmetry */
  if (Getv("amr_shells_symmetry","no") || !level->shells)
    return;
  
  /* if we use one of the symmetries, go on, if not symm is
     used, do nothing */
  if (! (Getv("grid", "bitant") || Getv("grid", "rotant") || 
         Getv("grid", "quadrant") || Getv("grid", "qreflect") || 
         Getv("grid", "octant")) ) 
    return;

  if (!shells_optimize(level)) 
    errorexit("does not work with non-optimized shell version... yet?!?");
  
  double *var, sym;
  int nvars = varlist->n;
  int *ivar = varlist->index; 
  int nv, iijjkk, d;
  
  double *zp = level->v[Ind("z")];
  
  
  if (Getv("grid", "bitant")) {
    
    if (PR) printf("set bitant symmetries in shells\n");
    
    for (nv = 0; nv < nvars; nv++) {
      
      sym = VarSymmetry(ivar[nv], 2);
      var = level->v[ivar[nv]];
      if (PR) printf("  set var: %s for b0-3     sym: %2.2f\n", VarName(ivar[nv]), sym);
    
      // x+/- and y+/- ... copy values from +z to -z inside each box
      for (d=1; d<4; d++) {
        
        forallpoints_boxijk(level->box[d]) {
         
          if (zp[ijk]<0.) {
            iijjkk   = level->box[d]->noffset + (kmax-k)*dk + (j)*dj + (i)*di ;
            var[ijk] = sym * var[iijjkk];
            
            if (debug)
            if (!(dequal(Ptr(level,"x")[ijk], Ptr(level,"x")[iijjkk]) && 
                  dequal(Ptr(level,"y")[ijk], Ptr(level,"y")[iijjkk]) && 
                  dequal(Ptr(level,"z")[ijk],-Ptr(level,"z")[iijjkk])))
              printf("xyz(b%d): %e %e %e <--> %e %e %e\n", d,
                     Ptr(level,"x")[ijk],Ptr(level,"y")[ijk],Ptr(level,"z")[ijk],
                     Ptr(level,"x")[iijjkk],Ptr(level,"y")[iijjkk],Ptr(level,"z")[iijjkk]);
          }
        } endfor_ijk;
      }
      
      if (PR) printf("  set var: %s for b4\n", VarName(ivar[nv]));
      
      // z direction (box 4 and 5) copy full box
      forallpoints_boxijk(level->box[4]) {
      
        iijjkk = level->box[5]->noffset +
                 level->box[5]->m*level->box[5]->n * (kmax-k) +
                 level->box[5]->m* (jmax-j) + (i);
        var[ijk] = sym * var[iijjkk];
        
        if (debug) {
          if (zp[ijk]>0.) 
            printf("write to z>0 ???\n");
          if (!(dequal(Ptr(level,"x")[ijk], Ptr(level,"x")[iijjkk]) && 
                dequal(Ptr(level,"y")[ijk], Ptr(level,"y")[iijjkk]) && 
                dequal(Ptr(level,"z")[ijk],-Ptr(level,"z")[iijjkk])))
            printf("xyz(b%d): %e %e %e <--> %e %e %e\n", 4,
                   Ptr(level,"x")[ijk],Ptr(level,"y")[ijk],Ptr(level,"z")[ijk],
                   Ptr(level,"x")[iijjkk],Ptr(level,"y")[iijjkk],Ptr(level,"z")[iijjkk]);
        }
        
      } endfor_ijk;
    
      
    }
    
  } else if (Getv("grid", "rotant")) {
    
    errorexit("rotant is not implemented!");
    
  } else {
    errorexit("shell-symmetry does not exist yet, implement me!");
  }
  
}





