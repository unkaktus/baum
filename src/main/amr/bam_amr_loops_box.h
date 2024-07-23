/* bam_amr_loops_box.h */
/* Bernd Bruegmann, 2/2006 */

/* Loops are performed by macros so that the user has to know very little
   about the implementation details.

   Basic considerations are
   - looping over different point sets, e.g. inner points, boundary points, ...
   - accessing neighbors of points


   Do we want to keep this for boxes?
   --
   Accessing neighbors is simplified and also made relatively efficient by
   providing a standard set of indices by which to access the nodes in 
   a stencil:
     ccc is the index of the center node 'i,j,k'
     mcc is the index of the i-1,j,k node
     pcc is the index of the i+1,j,k node
     Mcc is the index of the i-2,j,k node
     Pcc is the index of the i+2,j,k node
     m3cc is the index of the i-3,j,k node
     p3cc is the index of the i+3,j,k node
     nnn is ordered as ijk and xyz

   For loops without stencil the index variable is an argument for the macro.
*/





/****************************************************************************/
/* loops that are specific to boxes */

/* note: these loops define several variables implicitly! (for convenience)
         make sure there are no conflicts in the loop
	 e.g. since i,j,k are loop variables, don't use i for something else
*/


#define forallboxes(level) {				\
    int nbox; 						\
    for (nbox = 0; nbox < (level)->nboxes; nbox++) {    \
      tB *box = (level)->box[nbox];			\

#define endforboxes }}



/****************************************************************************/
/* NORMAL STUFF */

#define forallpoints_ijk(level) {			\
    const int *boundary = (level)->boundary;		\
    int nbox, ijk, i, j, k;                             \
    for (nbox = 0; nbox < (level)->nboxes; nbox++) {    \
      tB *box = (level)->box[nbox];                     \
      const int imax = box->m-1;			\
      const int jmax = box->n-1;			\
      const int kmax = box->o-1;			\
      const int di = box->di;				\
      const int dj = box->dj;				\
      const int dk = box->dk;				\
      ijk = box->noffset;  \
      for (k = 0; k <= kmax; k++)			\
      for (j = 0; j <= jmax; j++)                       \
      for (i = 0; i <= imax; i++, ijk++)		\

#define forinnerpoints_ijk(level)\
    forallpoints_ijk(level) if (!boundary[ijk])


#define forallpoints_boxijk(box) {{     \
    int ijk, i, j, k;         \
    const int imax = (box)->m-1;      \
    const int jmax = (box)->n-1;      \
    const int kmax = (box)->o-1;      \
    const int di = (box)->di;       \
    const int dj = (box)->dj;       \
    const int dk = (box)->dk;       \
    ijk = (box)->noffset;       \
    for (k = 0; k <= kmax; k++)       \
    for (j = 0; j <= jmax; j++)                         \
    for (i = 0; i <= imax; i++, ijk++)           \


#define forinnerpoints_boxijk(level,box)                \
    const int *boundary = (level)->boundary;    \
    forallpoints_boxijk(box) if (!boundary[ijk])      \


#define forplane_boxijk(forplanebox,direction) {      \
    const int imax = ((direction) == 1) ? 0 : (forplanebox)->m-1; \
    const int jmax = ((direction) == 3) ? 0 : (forplanebox)->n-1; \
    const int kmax = ((direction) == 5) ? 0 : (forplanebox)->o-1; \
    const int ijkoffset = (forplanebox)->noffset;     \
    const int di = (forplanebox)->di;         \
    const int dj = (forplanebox)->dj;         \
    const int dk = (forplanebox)->dk;         \
    int ijk, i, j, k;             \
    for (k = 0; k <= kmax; k++) \
    for (j = 0; j <= jmax; j++) \
    for (i = 0; i <= imax; i++) {                 \
      ijk = ijkoffset + i*di + j*dj + k*dk;       \

#define endfor_ijk }}






/****************************************************************************/
/* OPENMP STUFF */

#define forallpoints_ijk_openmp(level) {                  \
    const int *boundary = (level)->boundary;              \
    int nbox;                                             \
    for (nbox = 0; nbox < (level)->nboxes; nbox++) {      \
      tB *box = (level)->box[nbox];                       \
      const int imax = box->m-1;                          \
      const int jmax = box->n-1;                          \
      const int kmax = box->o-1;                          \
      const int di = box->di;                             \
      const int dj = box->dj;                             \
      const int dk = box->dk;                             \
      int ijk_0;                                          \
      bampi_openmp_loop                                   \
      for (ijk_0 = 0; ijk_0 < box->npoints; ijk_0++) {    \
        int k = ijk_0 / dk;                               \
        int j = (ijk_0 - dk*k) / dj;                      \
        int i = ijk_0 - dk*k - dj*j;                      \
        int ijk = ijk_0 + box->noffset;                   \

//      printf(" %d/%d     %d %d %d    %d %d %d     %d %d\n",  bampi_openmp_rank(), bampi_openmp_size(), box->m,box->n,box->o, i,j,k, box->noffset, ijk);

#define forinnerpoints_ijk_openmp(level)\
  forallpoints_ijk_openmp(level) if (!boundary[ijk])


#define forallpoints_ijk_openmp_collapse(level) {             \
  const int *boundary = (level)->boundary;         \
  int nbox, i,j,k, ijk;                            \
  for (nbox = 0; nbox < (level)->nboxes; nbox++) { \
    tB *box = (level)->box[nbox];                  \
    const int imax = box->m-1;  \
    const int jmax = box->n-1;  \
    const int kmax = box->o-1;  \
    const int di = box->di;     \
    const int dj = box->dj;     \
    const int dk = box->dk;     \
    bampi_openmp_loop_collapse3  \
    for (k = 0; k <= kmax; k++) \
    for (j = 0; j <= jmax; j++) \
    for (i = 0; i <= imax; i++) \
    {                           \
    ijk = box->noffset + i + j*box->m + k*box->m*box->n; \
      
//      printf(" %d/%d     %d %d %d    %d %d %d     %d %d\n",  bampi_openmp_rank(), bampi_openmp_size(), box->m,box->n,box->o, i,j,k, box->noffset, ijk);

#define forinnerpoints_ijk_openmp_collapse(level)\
  forallpoints_ijk_openmp_collapse(level) if (!boundary[ijk])


#define endfor_ijk_openmp }}}


#define forallpoints_openmp(level,i)  \
  int i;                              \
  bampi_openmp_loop                   \
  for (i = 0; i < level->npoints; i++)


/* End of OPENMP STUFF */
/****************************************************************************/ 




/* linear index, not the fastest version for tight loops */
#define ijkofbox(box,i,j,k) \
  ((box)->noffset + (i)*(box)->di + (j)*(box)->dj + (k)*(box)->dk)


/* 3d indices from linear index */
#define iofbox(box,ijk) ( ((ijk) - (box)->noffset) % (box)->dj)
#define jofbox(box,ijk) ((((ijk) - (box)->noffset) % (box)->dk) / (box)->dj)
#define kofbox(box,ijk) ( ((ijk) - (box)->noffset) / (box)->dk)


/* get coordinate without referring to coordinate variables */
#define xofbox(box,i) ((box)->x0 + (i)*(box)->dx)
#define yofbox(box,j) ((box)->y0 + (j)*(box)->dy)
#define zofbox(box,k) ((box)->z0 + (k)*(box)->dz)



/* check how far away the boundary is
   this is supposed to be used to determine the radius for centered stencils
   check the boundary flag beforehand (or use forinner) to avoid sym/ghosts
*/
#define boundary2ormore\
  (i >= 2 && j >= 2 && k >= 2 && i <= imax-2 && j <= jmax-2 && k <= kmax-2)
#define boundary3ormore\
  (i >= 3 && j >= 3 && k >= 3 && i <= imax-3 && j <= jmax-3 && k <= kmax-3)
#define boundary4ormore\
  (i >= 4 && j >= 4 && k >= 4 && i <= imax-4 && j <= jmax-4 && k <= kmax-4)

#define boundary1away\
  (i == 1 || j == 1 || k == 1 || i == imax-1 || j == jmax-1 || k == kmax-1)
#define boundary2away\
  (i == 2 || j == 2 || k == 2 || i == imax-2 || j == jmax-2 || k == kmax-2)

#define boundaryNaway(N)\
  (i == (N) || j == (N) || k == (N) ||\
   i == imax-(N) || j == jmax-(N) || k == kmax-(N))
#define boundaryNormore(N)\
  (i >= (N) && j >= (N) && k >= (N) &&\
   i <= imax-(N) && j <= jmax-(N) && k <= kmax-(N))

#define boundaryaway(N)\
  (boundaryaway_(N, i,j,k, imax,jmax,kmax))


/* is this useful?
#define MIN2(a,b) (((a) < (b)) ? (a) : (b))
#define MIN3(a,b,c) (((a) < (b)) ? MIN2(a,c) : MIN2(b,c))
*/



/****************************************************************************/
/* replacements for non-box loops */

/* for all points without stencil:
   for all is simple because this is how we store them 
   the name of the index is passed to the macro 
*/
#define forallpoints(level,i) \
  for (i = 0; i < level->npoints; i++)

/* for all inner points */
#define forinner1(level,i) \
  for (i = 0; i < level->npoints; i++) if (!level->boundary[i])

/* for all inner points including ghosts */
#define forallinner(level,i) \
  for (i = 0; i < level->npoints; i++) \
    if ((!level->boundary[i]) || level->boundary[i] == GHOBOUND)




/****************************************************************************/
/* declaration of neighbor indices */

#define int_nbinds_faces   int mcc, pcc, cmc, cpc, ccm, ccp
#define int_nbinds_faces2  int Mcc, Pcc, cMc, cPc, ccM, ccP
#define int_nbinds_faces3  int m3cc, p3cc, cm3c, cp3c, ccm3, ccp3
#define int_nbinds_corners int mmm, mmp, mpm, mpp, pmm, pmp, ppm, ppp
#define int_nbinds_edges \
  int mmc, mcm, cmm, ppc, pcp, cpp, mpc, mcp, cmp, pmc, pcm, cpm

#define int_nbinds_7  int ccc;       int_nbinds_faces
#define int_nbinds_13 int_nbinds_7;  int_nbinds_faces2
#define int_nbinds_19 int_nbinds_7;  int_nbinds_edges
#define int_nbinds_25 int_nbinds_19; int_nbinds_faces2
#define int_nbinds_27 int_nbinds_19; int_nbinds_corners




/****************************************************************************/
/* set neighbor indices 
   recall that the neighbors are ordered as
     0   1   2   3   4   5
   -di +di -dj +dj -dk +dk

   key issues are that 
   - neighbors may not exist
     in this case we introduce a special rule: set index to ccc
     later we can compare say Mcc == ccc to find out whether neighbor 
     existed without repeating pointer stuff
   - for special loops we actually guarantee by construction that certain
     neighbors exist, for example for inner points all nodes in the 
     27 node cube exist
*/


/* 7 point stencil assuming all nodes exist (inner points) */
#define set_nbinds_7(node) \
  ccc = ijk;		\
  mcc = ijk - di;	\
  pcc = ijk + di;	\
  cmc = ijk - dj;	\
  cpc = ijk + dj;	\
  ccm = ijk - dk;	\
  ccp = ijk + dk;	\

/* 19 point stencil assuming all nodes exist (inner points) */
/* should use previously defined nbinds_7 ... */
#define set_nbinds_19(node)			\
  set_nbinds_7(node);				\
  cmm = ijk - dj - dk;				\
  cpm = ijk + dj - dk;				\
  cmp = ijk - dj + dk;     			\
  cpp = ijk + dj + dk;     			\
  mmc = ijk - di - dj;     			\
  mpc = ijk - di + dj;     			\
  mcm = ijk - di - dk;     			\
  mcp = ijk - di + dk;     			\
  pmc = ijk + di - dj;     			\
  ppc = ijk + di + dj;     			\
  pcm = ijk + di - dk;     			\
  pcp = ijk + di + dk;     			\

/* 27 point cube assuming all nodes exist (inner points) */
#define set_nbinds_27(node)			\
  set_nbinds_19(node);				\
  mmm = ijk - di - dj - dk;       		\
  pmm = ijk + di - dj - dk;       		\
  mpm = ijk - di + dj - dk;       		\
  mmp = ijk - di - dj + dk;       		\
  mpp = ijk - di + dj + dk;       		\
  pmp = ijk + di - dj + dk;       		\
  ppm = ijk + di + dj - dk;       		\
  ppp = ijk + di + dj + dk;       		\


/* 7 point stencil where some nodes may not exist */
#define set_nbinds_check_7(node)			\
  ccc = ijk;  					\
  mcc = (i > 0) ? ijk-di : ijk;			\
  cmc = (j > 0) ? ijk-dj : ijk;			\
  ccm = (k > 0) ? ijk-dk : ijk;			\
  pcc = (i+1 < box->m) ? ijk+di : ijk;		\
  cpc = (j+1 < box->n) ? ijk+dj : ijk;		\
  ccp = (k+1 < box->o) ? ijk+dk : ijk;		\

/* 12 points for 12 edges where some nodes may not exist */
#define set_nbinds_check_12(nodes) \
mmc = (i>0 && j>0) ? ijk-di-dj : ijk;\
mcm = (i>0 && k>0) ? ijk-di-dk : ijk;\
cmm = (j>0 && k>0) ? ijk-dj-dk : ijk;\
mpc = (i>0 && j+1<box->n) ? ijk-di+dj : ijk;\
mcp = (i>0 && k+1<box->o) ? ijk-di+dk : ijk;\
cmp = (j>0 && k+1<box->o) ? ijk-dj+dk : ijk;\
pmc = (j>0 && i+1<box->m) ? ijk+di-dj : ijk;\
pcm = (k>0 && i+1<box->m) ? ijk+di-dk : ijk;\
cpm = (k>0 && j+1<box->m) ? ijk+dj-dk : ijk;\
ppc = (i+1<box->m && j+1<box->n) ? ijk+di+dj : ijk;\
pcp = (i+1<box->m && k+1<box->o) ? ijk+di+dk : ijk;\
cpp = (j+1<box->n && k+1<box->o) ? ijk+dj+dk : ijk;\

/* 8 points for 8 corners where some nodes may not exist */
#define set_nbinds_check_8(nodes) \
mmm = (mmc != ccc && k>0       ) ? mmc-dk : ccc;\
mmp = (mmc != ccc && k+1<box->o) ? mmc+dk : ccc;\
mpm = (mpc != ccc && k>0       ) ? mpc-dk : ccc;\
mpp = (mpc != ccc && k+1<box->o) ? mpc+dk : ccc;\
pmm = (pmc != ccc && k>0       ) ? pmc-dk : ccc;\
pmp = (pmc != ccc && k+1<box->o) ? pmc+dk : ccc;\
ppm = (ppc != ccc && k>0       ) ? ppc-dk : ccc;\
ppp = (ppc != ccc && k+1<box->o) ? ppc+dk : ccc;\


/* 6 points for step 2 faces where some nodes may not exist */
#define set_nbinds_check2_6(nodes) \
Mcc = (i>1       ) ? mcc-di : ccc;\
Pcc = (i+2<box->m) ? pcc+di : ccc;\
cMc = (j>1       ) ? cmc-dj : ccc;\
cPc = (j+2<box->n) ? cpc+dj : ccc;\
ccM = (k>1       ) ? ccm-dk : ccc;\
ccP = (k+2<box->o) ? ccp+dk : ccc;\


/* 6 points for step 3 faces where some nodes may not exist */
#define set_nbinds_check3_6(nodes) \
m3cc = (i>2       ) ? Mcc-di : ccc;\
p3cc = (i+3<box->m) ? Pcc+di : ccc;\
cm3c = (j>2       ) ? cMc-dj : ccc;\
cp3c = (j+3<box->n) ? cPc+dj : ccc;\
ccm3 = (k>2       ) ? ccM-dk : ccc;\
ccp3 = (k+3<box->n) ? ccP+dk : ccc;\


/* 13 point stencil where some points may not exist
   note:
   if ccm == ccc, then in set_nbinds_check2_6 we are just repeating 
   what was done in set_nbinds_check_6 with the result that ccM == ccc, too!
*/
#define set_nbinds_check_13(n,nodes) \
set_nbinds_check_7(n);\
set_nbinds_check2_6(nodes);\

/* 19 point stencil where some points may not exits */
#define set_nbinds_check_19(n,nodes) \
set_nbinds_check_7(n);\
set_nbinds_check_12(nodes);\

/* 25 point stencil where the step 2 points may not exist */
#define set_nbinds_check2_25(n,nodes) \
set_nbinds_19(n);\
set_nbinds_check2_6(nodes);\

/* 27 point stencil where some points may not exits */
#define set_nbinds_check_27(n,nodes) \
set_nbinds_check_7(n);\
set_nbinds_check_12(nodes);\
set_nbinds_check_8(nodes);\





/****************************************************************************/
/* special indices and weights for advection stencils (one-sided derivative)
 
   The numerical values of the weights should be defined elsewhere,
   but this is our first example for such a stencil so we leave it here 
   for now.

   For example, v^x del_x u 
   becomes      v^x_i (3 u_i - 4 u_(i-1) + u_(i-2)) / (2h)
   if v^x is negative (upwinding from the left).

   Note that for us del_t u = + v^x del_x u which advects to the left for
   v^x positive.

   Here we provide special neighbor indices maa, Maa, ... and corresponding
   weights wmaa, wMaa, ..., ('a' as in advection)
   while the actual calculation is done elsewhere.
*/

/* set weights and indices for one given direction */
#define set_advection_onedir(vx,waaa,wmaa,wMaa,maa,Maa,Mcc,mcc,pcc,Pcc) \
wwww = vx;                                                              \
if (wwww < 0 && Mcc != ccc) {                                           \
    waaa += 3*wwww; wmaa = -4*wwww; wMaa =  wwww; maa = mcc; Maa = Mcc; \
}                                                                       \
else if (wwww > 0 && Pcc != ccc) {                                      \
    waaa -= 3*wwww; wmaa =  4*wwww; wMaa = -wwww; maa = pcc; Maa = Pcc; \
}                                                                       \
else {                                                                  \
                    wmaa =   -wwww; wMaa =  wwww; maa = mcc; Maa = pcc; \
} 



/* set weights and indices for all three directions */
#define set_advection(vx,vy,vz)                                         \
waaa = 0;                                                               \
set_advection_onedir(vx,waaa,wmaa,wMaa,maa,Maa,Mcc,mcc,pcc,Pcc);        \
set_advection_onedir(vy,waaa,wama,waMa,ama,aMa,cMc,cmc,cpc,cPc);        \
set_advection_onedir(vz,waaa,waam,waaM,aam,aaM,ccM,ccm,ccp,ccP);        \


/* debugging */
#define debugadvection \
{ double vx = oo2dx*beta1[ccc];\
  double vy = oo2dy*beta2[ccc];\
  double vz = oo2dz*beta3[ccc];\
printf("%10.3e %10.3e %10.3e  %10.3e  %10.3e %10.3e  %10.3e %10.3e  %10.3e %10.3e\n",\
vx, vy, vz, waaa, wmaa/vx, wMaa/vx, wama/vy, waMa/vy, waam/vz, waaM/vz);} \


/* declaration of weights and indices for advection */
#define int_advection                                                   \
double wwww, waaa, wmaa, wama, waam, wMaa, waMa, waaM;                  \
int maa, ama, aam, Maa, aMa, aaM                                        \




/****************************************************************************/
/****************************************************************************/
/* loop over nodes with indices set for various stencils 
   could be made even shorter but perhaps less clear by allowing
   an argument for the stencil being passed, say forinner(27,level),
   with a preprocessor #if case distinction
*/


/* local shorthand for just this file */
#define forloopstart(level)                             \
  int nbox, ijk, i, j, k, di, dj, dk;			\
  int *boundary = level->boundary;			\
  int node, nodes; /* dummies */			\
  for (nbox = 0; nbox < (level)->nboxes; nbox++) {	\
    tB *box = (level)->box[nbox];			\
    ijk = box->noffset;					\
    di = box->di;					\
    dj = box->dj;					\
    dk = box->dk;					\
    const int imax = (box)->m-1;			\
    const int jmax = (box)->n-1;			\
    const int kmax = (box)->o-1;			\
    for (k = 0; k < box->o; k++)			\
    for (j = 0; j < box->n; j++)                        \
    for (i = 0; i < box->m; i++, ijk++)


/* for inner points with neighbors that are guaranteed to exist */
#define forinner7(level) \
{\
  int_nbinds_7;\
  forloopstart(level)\
  if (!boundary[ijk]) {\
      set_nbinds_7(node);\

  
#define forinner19(level) {\
  int_nbinds_19;\
  forloopstart(level)\
  if (!boundary[ijk]) {\
    set_nbinds_19(node);\

#define forinner27(level) {\
  int_nbinds_27;\
  forloopstart(level)\
  if (!boundary[ijk]) {\
      set_nbinds_27(node);

#define forinner27nostatus(level) {\
  int_nbinds_27;\
  forloopstart(level)\
  if (!boundary[ijk]) {\
      set_nbinds_27(node);

#define forinnerpoints forinner27
  

/* for all points with special flag checking for one way existence */
#define forall27flag(level,flag) {\
  int_nbinds_27;\
  forloopstart(level)\
  if (boundary[ijk] == flag) {\
    set_nbinds_check_27(node,nodes);\



/* for inner points including e.g. excision */
#define forinner19plus(level,flag) {\
  int_nbinds_19;\
  forloopstart(level)\
  if ((!boundary[ijk]) || boundary[ijk] == flag) {\
    set_nbinds_19(node);\

/* for inner points with additional 2 step advection stencil */
#define forinner25(level,vx,vy,vz) {\
  int_nbinds_25;\
  int_advection;\
  forloopstart(level)\
  if (!boundary[ijk]) {\
    set_nbinds_check2_25(node,nodes);\
    set_advection(vx,vy,vz);\

/* for boundary points with 13 point stencil */
#define forboundary13(level,flag) {\
  int_nbinds_13;\
  forloopstart(level)\
  if (boundary[ijk] == flag) {\
    set_nbinds_check_13(node,nodes);\


/* every loop must be terminated */
#define endfor\
    }\
  }\
}\

#define endforinner endfor  



/**************************************************************************/
/* point lists */
typedef struct tPL
{
  tL  *level;   /* level to which points belong */
  int npoints;  /* number of points in list */
  int *point;   /* array containing indices of all the points in the list */ 
} tPointList;


/* loop over all points i in list pointlist (end this loop with endfor;) */
#define forPointList(pointlist, i) { {\
  int forPointList_i;\
  for( forPointList_i=0; forPointList_i<(pointlist)->npoints; forPointList_i++) {\
    i=(pointlist)->point[forPointList_i];

/* ***********************************************************************
NOTE: the macros
 forPointList19, forPointList27, 
 addtoPointList_check_7, addtoPointList_check_19, addtoPointList_check_27
rely on 
 set_nbinds_19, set_nbinds_27,
 set_nbinds_check_7, set_nbinds_check_19, set_nbinds_check_27
which in turn only work if
 ijk, i,j,k, di,dj,dk, box
are defined and set correctly! Below I define declare_i_j_k_di_dj_dk_box
and compute_i_j_k_di_dj_dk_box to do just that. But they are very
inefficient...
************************************************************************** */
/* define and set i,j,k, di,dj,dk, box */
#define declare_i_j_k_di_dj_dk_box \
  int i,j,k, di,dj,dk;	\
  tB *box;		\
  int ni, nj, nk;	\
  int nbox

#define compute_i_j_k_di_dj_dk_box(level) \
  for(nbox = 0, box=(level)->box[0]; nbox < (level)->nboxes; nbox++, box=(level)->box[nbox])\
    if(ijk >= box->noffset) break;				\
  ni = box->m;							\
  nj = box->n;							\
  nk = box->o;							\
  di = box->di;							\
  dj = box->dj;							\
  dk = box->dk;							\
  k = ijk/(ni*nj);   /* ijk = i + ni*j + ni*nj*k */		\
  j = (ijk-ni*nj*k)/ni;						\
  i = ijk-ni*nj*k-ni*j


/* loop over all points i in list pointlist and set the 19 points
   ccc, mcc, ... cpm  at each point i (end this loop with endfor;)   */
#define forPointList19(pointlist) {\
  int_nbinds_19;\
  int ijk;\
  int forPointList_i;\
  for( forPointList_i=0; forPointList_i<(pointlist)->npoints; forPointList_i++) {\
    ijk = (pointlist)->point[forPointList_i];\
    if (!(pointlist)->level->boundary[ijk]) {\
      declare_i_j_k_di_dj_dk_box;\
      compute_i_j_k_di_dj_dk_box((pointlist)->level);\
      set_nbinds_19(node);  


/* loop over all points i in list pointlist and set the 27 points
   ccc, mcc, ... ppp  at each point i (end this loop with endfor;)   */
#define forPointList27(pointlist) {\
  int_nbinds_27;\
  int ijk;\
  int forPointList_i;\
  for( forPointList_i=0; forPointList_i<(pointlist)->npoints; forPointList_i++) {\
    ijk = (pointlist)->point[forPointList_i];\
    if (!(pointlist)->level->boundary[ijk]) {\
      declare_i_j_k_di_dj_dk_box;\
      compute_i_j_k_di_dj_dk_box((pointlist)->level);\
      set_nbinds_27(node);  


/* add the 7 points ccc, mcc, ... ccp at point i to pointlist */
#define addtoPointList_check_7(pointlist, level, i) {\
  int_nbinds_7;\
  int ijk = i; /* tN *node = (level)->node + (i); */ \
  declare_i_j_k_di_dj_dk_box;\
  compute_i_j_k_di_dj_dk_box((level));\
  set_nbinds_check_7(node);\
  AddToPointList((level), (pointlist), ccc); \
  if(mcc!=ccc) AddToPointList((level), (pointlist), mcc); \
  if(pcc!=ccc) AddToPointList((level), (pointlist), pcc); \
  if(cmc!=ccc) AddToPointList((level), (pointlist), cmc); \
  if(cpc!=ccc) AddToPointList((level), (pointlist), cpc); \
  if(ccm!=ccc) AddToPointList((level), (pointlist), ccm); \
  if(ccp!=ccc) AddToPointList((level), (pointlist), ccp);   }


/* add the 19 points ccc, mcc, ... cpm at point i to pointlist */
#define addtoPointList_check_19(pointlist, level, i) {\
  int_nbinds_19;\
  int ijk = i; /* tN *node = (level)->node + (i); */	\
  int nodes;   /* tN *nodes = level->node; */		\
  declare_i_j_k_di_dj_dk_box;\
  compute_i_j_k_di_dj_dk_box((level));\
  set_nbinds_check_19(node,nodes);\
  AddToPointList((level), (pointlist), ccc); \
  if(mcc!=ccc) AddToPointList((level), (pointlist), mcc); \
  if(pcc!=ccc) AddToPointList((level), (pointlist), pcc); \
  if(cmc!=ccc) AddToPointList((level), (pointlist), cmc); \
  if(cpc!=ccc) AddToPointList((level), (pointlist), cpc); \
  if(ccm!=ccc) AddToPointList((level), (pointlist), ccm); \
  if(ccp!=ccc) AddToPointList((level), (pointlist), ccp); \
  if(mmc!=ccc) AddToPointList((level), (pointlist), mmc); \
  if(mcm!=ccc) AddToPointList((level), (pointlist), mcm); \
  if(cmm!=ccc) AddToPointList((level), (pointlist), cmm); \
  if(ppc!=ccc) AddToPointList((level), (pointlist), ppc); \
  if(pcp!=ccc) AddToPointList((level), (pointlist), pcp); \
  if(cpp!=ccc) AddToPointList((level), (pointlist), cpp); \
  if(mpc!=ccc) AddToPointList((level), (pointlist), mpc); \
  if(mcp!=ccc) AddToPointList((level), (pointlist), mcp); \
  if(cmp!=ccc) AddToPointList((level), (pointlist), cmp); \
  if(pmc!=ccc) AddToPointList((level), (pointlist), pmc); \
  if(pcm!=ccc) AddToPointList((level), (pointlist), pcm); \
  if(cpm!=ccc) AddToPointList((level), (pointlist), cpm);   }


/* add the 27 points ccc, mcc, ... ppp at point i to pointlist */
#define addtoPointList_check_27(pointlist, level, i) {\
  int_nbinds_27;\
  int ijk = i; /* tN *node = (level)->node + (i); */	\
  int nodes;   /* tN *nodes = level->node; */		\
  declare_i_j_k_di_dj_dk_box;\
  compute_i_j_k_di_dj_dk_box((level));\
  set_nbinds_check_27(node,nodes);\
  AddToPointList((level), (pointlist), ccc); \
  if(mcc!=ccc) AddToPointList((level), (pointlist), mcc); \
  if(pcc!=ccc) AddToPointList((level), (pointlist), pcc); \
  if(cmc!=ccc) AddToPointList((level), (pointlist), cmc); \
  if(cpc!=ccc) AddToPointList((level), (pointlist), cpc); \
  if(ccm!=ccc) AddToPointList((level), (pointlist), ccm); \
  if(ccp!=ccc) AddToPointList((level), (pointlist), ccp); \
  if(mmc!=ccc) AddToPointList((level), (pointlist), mmc); \
  if(mcm!=ccc) AddToPointList((level), (pointlist), mcm); \
  if(cmm!=ccc) AddToPointList((level), (pointlist), cmm); \
  if(ppc!=ccc) AddToPointList((level), (pointlist), ppc); \
  if(pcp!=ccc) AddToPointList((level), (pointlist), pcp); \
  if(cpp!=ccc) AddToPointList((level), (pointlist), cpp); \
  if(mpc!=ccc) AddToPointList((level), (pointlist), mpc); \
  if(mcp!=ccc) AddToPointList((level), (pointlist), mcp); \
  if(cmp!=ccc) AddToPointList((level), (pointlist), cmp); \
  if(pmc!=ccc) AddToPointList((level), (pointlist), pmc); \
  if(pcm!=ccc) AddToPointList((level), (pointlist), pcm); \
  if(cpm!=ccc) AddToPointList((level), (pointlist), cpm); \
  if(mmm!=ccc) AddToPointList((level), (pointlist), mmm); \
  if(mmp!=ccc) AddToPointList((level), (pointlist), mmp); \
  if(mpm!=ccc) AddToPointList((level), (pointlist), mpm); \
  if(mpp!=ccc) AddToPointList((level), (pointlist), mpp); \
  if(pmm!=ccc) AddToPointList((level), (pointlist), pmm); \
  if(pmp!=ccc) AddToPointList((level), (pointlist), pmp); \
  if(ppm!=ccc) AddToPointList((level), (pointlist), ppm); \
  if(ppp!=ccc) AddToPointList((level), (pointlist), ppp);   }





/**************************************************************************/
/* Multi level(Ml) point lists */
typedef struct tMlPL
{
  int	     nlevels;  	     /* number of levels in list */
  tL  	     **level;	     /* array containing all the levels in the list */
  tPointList **PointList;    /* array containing all the single level
				PointLists in the list */
  tPointList *SlPointListResult; /* temp pointer to a single level PointList,
				    needed to store result of func SlPointList */
} tMlPointList;


/* loop over all points i in list Mlpointlist on level Level
   (end this loop with endfor;)                              */
#define forMlPointList(Level,Mlpointlist, i) \
  SlPointList((Level),(Mlpointlist)); \
  if( (Mlpointlist)->SlPointListResult != NULL) \
    forPointList( ((Mlpointlist)->SlPointListResult), (i))
  

/* loop over all points i in list Mlpointlist on level Level
   and set the 19 points
   ccc, mcc, ... cpm  at each point i (end this loop with endfor;)   */
#define forMlPointList19(Level,Mlpointlist) \
  SlPointList((Level),(Mlpointlist)); \
  if( (Mlpointlist)->SlPointListResult != NULL) \
    forPointList19( ((Mlpointlist)->SlPointListResult) )


/* loop over all points i in list Mlpointlist on level Level
   and set the 27 points
   ccc, mcc, ... ppp  at each point i (end this loop with endfor;)   */
#define forMlPointList27(Level,Mlpointlist) \
  SlPointList((Level),(Mlpointlist)); \
  if( (Mlpointlist)->SlPointListResult != NULL) \
    forPointList27( ((Mlpointlist)->SlPointListResult) )


/* add the 7 points ccc, mcc, ... ccp at point i on level Level
   to Mlpointlist */
#define addtoMlPointList_check_7(Level, Mlpointlist, i) \
  AddLevelToMlPointList((Level),(Mlpointlist));\
  addtoPointList_check_7(SlPointList((Level),(Mlpointlist)),(Level),(i))


/* add the 19 points ccc, mcc, ... cpm at point i on level Level
   to Mlpointlist */
#define addtoMlPointList_check_19(Level, Mlpointlist, i) \
  AddLevelToMlPointList((Level),(Mlpointlist));\
  addtoPointList_check_19(SlPointList((Level),(Mlpointlist)),(Level),(i))


/* add the 27 points ccc, mcc, ... ppp at point i on level Level
   to pointlist */
#define addtoMlPointList_check_27(Level, Mlpointlist, i) \
  AddLevelToMlPointList((Level),(Mlpointlist));\
  addtoPointList_check_27(SlPointList((Level),(Mlpointlist)),(Level),(i))

        



/****************************************************************************/
/* fill in an array with the available indices 
   useful for loops over 3 dimensions
*/
#define fillindexarray27(indexarray) \
          indexarray[1][1][1] = ccc; \
          indexarray[0][1][1] = mcc; \
	  indexarray[2][1][1] = pcc; \
	  indexarray[1][0][1] = cmc; \
	  indexarray[1][2][1] = cpc; \
	  indexarray[1][1][0] = ccm; \
	  indexarray[1][1][2] = ccp; \
	  indexarray[0][0][1] = mmc; \
	  indexarray[0][1][0] = mcm; \
	  indexarray[1][0][0] = cmm; \
	  indexarray[2][2][1] = ppc; \
	  indexarray[2][1][2] = pcp; \
	  indexarray[1][2][2] = cpp; \
	  indexarray[0][2][1] = mpc; \
	  indexarray[0][1][2] = mcp; \
	  indexarray[1][0][2] = cmp; \
	  indexarray[2][0][1] = pmc; \
	  indexarray[2][1][0] = pcm; \
	  indexarray[1][2][0] = cpm; \
	  indexarray[0][0][0] = mmm; \
	  indexarray[0][0][2] = mmp; \
	  indexarray[0][2][0] = mpm; \
	  indexarray[0][2][2] = mpp; \
	  indexarray[2][0][0] = pmm; \
	  indexarray[2][0][2] = pmp; \
	  indexarray[2][2][0] = ppm; \
	  indexarray[2][2][2] = ppp

