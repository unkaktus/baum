/* bam_amr.h */
/* Bernd Bruegmann, 12/99, 2/2006 */

#pragma once

/*************************************************************************/
/* basic structures */

/* a box */
typedef struct tBOX {
  struct tBOX *pr;            /* parent if available */
  struct tLEVEL *level;       /* level */
  int l;                      /* number of level */
  int i;                      /* index of box in list of boxes of level */
  double dx, dy, dz;          /* grid spacing */
  double dd[3];               /* grid spacing for each direction */
  int di, dj, dk;             /* stride to increment linear index */
  int id[3];                  /* increments for each direction */
  double bbox[6];             /* global bounding box xmin, xmax, ... */
  int ibbox[6];               /* global bounding box imin, imax, ... */

                              /* local information */
  struct tCOM *com;           /* communication info for parallelization */
  int npoints;                /* total number of points */
  int noffset;                /* offset into linear storage of level */
  double x0, y0, z0;          /* origin of coordinates, x0 = xmin, ... */
  int m, n, o;                /* number of points in the three directions */
  int imax, jmax, kmax;       /* number of points, imax=m-1, ... */
  int bflag[6];               /* boundary flag for each face */
} tB;

typedef tB tBox;




/* access boundary flags, which are stored either in tN or tL
   note: this macro is also used as lvalue (!)
*/
#define boundaryflag(level,i) level->boundary[i]





/* a level */
typedef struct tLEVEL {
  struct tGRID *grid;         /* pointer to parent grid */
  struct tLEVEL *prlocal;     /* pointer to ghost parent level */
  int l;                      /* index of this level */
  int shells;                 /* flag whether if it is a cartesian box 
                                 or shell coordinates */

  int nboxes;                 /* number of pointers to boxes */
  tB **box;                   /* list of pointers to boxes */

  int npoints;                /* number of points = npoints = nnodes */
  int nnodes;                 /* number of points = npoints = nnodes */

  int *boundary;              /* list of boundary flags */

  int nvariables;             /* number of pointers to variables */
  double **v;                 /* list of data pointers to variables */
  int iteration;              /* number of iterations */
  double time;                /* physical time */
  double dx, dy, dz;          /* spatial grid spacing */
  double dt;                  /* time step */
  double bbox[6];             /* global bounding box */
  int ibbox[6];               /* global bounding box in index range */
  struct tSBOX *sbox;         /* special boxes */
  int nsboxes;                /* number of special boxes */
  struct tCOM *com;           /* communication info for parallelization
                                 this structure is defined in a header file 
                                 that is included in bam.h */
} tL;




/* several levels make up a numerical grid */
typedef struct tGRID {
  tL **level;                 /* list of pointers to levels */
  int nlevels;                /* number of levels */
  int lmin;                   /* index of coarsest level */
  int lmax;                   /* index of finest level */
  int ltop;                   /* index of coarsest non-multigrid level */
  int lbo;                    /* index of coarsest Berger-Oliger level */
  int nvariables;             /* number of variables */
  int half[3];                /* flags for x/y/z halves */
  int symmetric[3];           /* flags for x/y/z reflection symmetries */
  
  /* symmetrie flags, ... no more statics!!! */
  int full,bitant,rotant,quadrant,octant,qreflect;
  
  int npunc;                  /* number of punctures we use */
  int *lmaxpunc;              /* for differnt numbers of levels */
  double **puncpos;           /* actual position of the puncture */
} tG;


//Conservative Adaptive Mesh refinement routines
typedef struct tCAMR {
  void (*correct_varlist) ();
  void (*mask) ();
} tCAMR;
extern tCAMR CAMR;


/**************************************************************************/



/* various values of boundary flag */
#define REFBOUND -1
#define NOTBOUND  0
#define PHYBOUND  1
#define GHOBOUND  2
#define SYMBOUND  3
#define EXCBOUND  4
#define EXCINTER  5
#define PERBOUND  6
#define OUTBOUND  7


/* variable lists 
   there should be a better place for it */
typedef struct tVARLIST {
  struct tLEVEL *level;
  int n;
  int *index;
} tVarList;


/* sometimes we want to keep track of several boxes within a level */
typedef struct tSBOX {
  double bbox[6];
  int ibbox[6];
} tSBox;




/**************************************************************************/
/* loops */
#include "bam_amr_loops_box.h"

/**************************************************************************/
/* advection */
#include "bam_amr_advection.h"


/**************************************************************************/
/* functions */

/* box.c */
int box_parent_index(double cx0, double cdx, double fx0, double fdx, int fi);
int box_child_index(double fx0, double fdx, double cx0, double cdx, int ci);
int box_point_index(double X0, double dX, double x0, double dx, int i);
void box_parent_indices(tB *c, tB *f, int *pci, int *pcj, int *pck,
			int fi, int fj, int fk);
void box_point_indices(tB *c, tB *f, int *pci, int *pcj, int *pck,
		       int fi, int fj, int fk);
int box_parent_ijk(tB *c, tB *f, int fi, int fj, int fk);
int box_child_ijk(tB *f, tB *c, int ci, int cj, int ck);
tB *box_containing_xyz(tL *level, double x, double y, double z); 
tB *box_containing_ijk(tL *level, tB *b, int i, int j, int k);
tL *one_box_level(tL *level, int nbox);
int boundaryaway_(int N, int i,int j,int k, int imax,int jmax,int kmax);
int test_box_consistency(tG* g, int level);
int test_box_consistency_level(tL* lc, tL* lf);

/* check.c */
int ExitIfNAN(tL* level);
int CheckForNANandINF(int N, ...);
int CheckIfFinite_VarIndexOne(tL* level, int ivar);
int CheckIfFinite(tL* level, char *varname);
int CheckIfFinite_VarIndex(tL* level, int i);
int CheckIfFinite_vl(tVarList *vl);

/* flag.c */
int flagsimplegrid(tL *level);
int flaggradedboxes(tL *level);
int FlagRegrid_BoxesAroundCenter(tL *level);
void flagregridboundary(tL *level);
void set_refinement_excision(tL *level);
void unset_refinement_excision(tL *level);

/* grid.c */
tG *make_grid(int pr);
tL *make_level(double dx, double dy, double dz, int m, int n, int o, int pr);
void test_grid_symmetry(tG* g);

/* inject.c */
int parentaligned(tL *level);
void enableparent(tL *lf, tL **plc);
void restrict_prolong_varlist(tG *g, int lcoarse, int lfine, 
			      tVarList *uc, tVarList *uf, int nbuffer);
void restrict_prolong_evolve(tG *g, int lcoarse, int lfine);
void restrict_prolong_grid(tG *g, tVarList *u);
void set_boundary_refinement_cf(tG *g, int lcoarse, int lfine, tVarList *u);
void set_boundary_refinement(tL *level, tVarList *u); 

/* move.c */
int box_overlap(double *a, double *b);
int box_ainb(double *a, double *b);
int move_box(tG *g, int l);

/* points.c */
void findbbox(tL *level, double *bbox, int *ibbox);
void findbbox_flag(tL *level, double *bbox, int *ibbox, double *flag);
void findbbox_notghost(tL *level, double *bbox, int *ibbox);
void bufferorderedpoints(tL *level, int npoints, double *coords, 
			 int nv, int *iv, int *nbuf, double **buf,
			 int interpolate, int order, int scheme);
void buffernonorderedpoints(tL *level, int npoints, double *coords, 
			    int nv, int *iv, int *nbuf, double **buf, 
			    int interpolate);
void bufferpoints_full(tL *level, int npoints, double *coords, 
                       int nv, int *iv, int *nbuf, double **buf, 
                       int order, int scheme);
int find_one_point(tL *level, double *coord);
int find_one_point_box(tL *level, double x, double y, double z);
int insidebbox(double *bbox, double *coord);
int xyzinsidebbox(double *bbox, double x, double y, double z);
int xyzinsidelevel(tL* level, double x, double y, double z);
int ijkinsidefinerlevel(tB *box, int ijk);
int insideownership(tL *level, double x, double y, double z);
int sphere_inside_bbox(double *bbox, 
		       double x0, double y0, double z0, double r);
int sphere_inside_level(tL *level, 
			double x0, double y0, double z0, double r);
int sphere_safely_inside_level(tL *level, 
			       double x0, double y0, double z0, double r);
int centered_sphere_safely_inside_level(tL *level, double r);

/* Single level point lists */
tPointList *AllocatePointList(tL *level);
void AddToPointList(tL *level, tPointList *PL, int newpoint);
void FreePointList(tPointList *PL);
void prPointList(tPointList *PL);

/* Multi level point lists */
tMlPointList *AllocateMlPointList(void);
tPointList *SlPointList(tL *level, tMlPointList *MlPL);
void AddLevelToMlPointList(tL *level, tMlPointList *MlPL);
void AddToMlPointList(tL *level, tMlPointList *MlPL, int newpoint);
void FreeMlPointList(tMlPointList *MlPL);
void prMlPointList(tMlPointList *MlPL);

/* print.c */
void printgrid(tG *g);
void printlevel(tL *l);
void printbbox(tL *level, double *bbox, int *ibbox);
void printboundary(tL *level);
void printvariables(tL *level);
void printvarinfo(tL *level, int i);

/* regrid.c */
void regrid(tG *g, int l, int secondcall);

/* shells.c */
void add_grid_shells(tG *grid);
void convert_box_to_shells(int b, 
                           double x, double y, double z,
                           double *r, double *phi, double *theta);
void convert_shells_to_box(int b, 
                           double r, double phi, double theta,
                           double *x, double *y, double *z);
int find_shellsbox_from_xyz(double x, double y, double z);
double compute_fisheye(double x, int del);
void set_boundary_shells_symmetry(tL *level, tVarList *varlist);
void sync_shells(tL* level0, int nv, int *iv);

/* storage.c */
void printenabled(tL *level);
void enablevarcomp(const tL *level, const int i);
void disablevarcomp(const tL *level, const int i);
void enablevar(tL *level, int i);
void disablevar(const tL *level, const int i);
void enablevarlist(const tVarList *v);
void disablevarlist(const tVarList *v);
void enablesamevars(tL *level, tL *newlevel);
double *Ptr(const tL *level, char *name);
double *PtrEnable(const tL *level, const char *name);
void PtrDisable(const tL *level, char *name);

void realloc_levelvariables(tL *level, int nvariables);
void free_level(tL *level);
void remove_levels(tL *level);
tG *alloc_grid(int lmin, int lmax, int nvariables);
void free_grid_only(tG *g);
void free_grid(tG *g);

void replace_level(tG *g, tL *level, int l);
void insert_level(tG *g, tL *l, int m);
void append_level(tG *g, tL *l);
void prepend_level(tG *g, tL *l);
void push_level(tG *g, tL *l);
void remove_top_level(tG *g);




