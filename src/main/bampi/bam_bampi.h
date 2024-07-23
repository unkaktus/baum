/* bam_bampi.h */
/* Bernd Bruegmann, 5/02 */

#include "bam_profile.h"
#include "bam_openmp.h"


/* buffer for parallelization */
typedef struct {
  int ni;            /* number of indices */
  int *i;            /* list of indices */
  int nbuffer;       /* number of data entries, can be a multiple of ni */
  double *buffer;    /* data */
} tComBuf;

/* communication data */
typedef struct tCOM {
  int size;          /* number of processors */
  int myrank;        /* my rank, which in MPI means my processor number */
  int nbrank[6];     /* ranks of my neighbors in cartesian processor grid */
  int sizexyz[3];    /* number of processors in xyz direction */
  int myrankxyz[3];  /* my rank in xyz direction */
  double bbox[6];    /* local bounding box */
  int ibbox[6];      /* local bounding box in index range */
  double bboxown[6]; /* local bounding box for ownership */
  tComBuf *send[6];  /* buffers for sending for each direction */
  tComBuf *recv[6];  /* buffers for receiving for each direction */
  tComBuf ***buflist;/* buffers for FMR parent communication */
} tCom;



/* for convenience here is a macro to check for processor 0
   we do not want to use global variables to avoid lib dependencies */
#define processor0 (!bampi_rank())



/* bampi.c */
void bampi_initialize(int *pargc, char ***pargv);
void bampi_finalize(int argc, char **argv);
void bampi_abort(void);
int bampi_rank(void);
int bampi_size(void);
double bampi_time(void);
double bampi_uptime(void);
double bampi_delta_time(void);
void bampi_barrier(void);
void bampi_bcast_int(void *buffer, int count, int root);
void bampi_bcast_char(void *buffer, int count, int root);
void bampi_bcast_double(void *buffer, int count, int root);
void bampi_alltoall(void *bufsend, void *bufrecv, int count);


int bampi_or(int flag);
void bampi_swap_buffer(void *buffer, int count, int srcrank, int dstrank);
void bampi_isendrecv(void *bufsend, void *bufrecv, int count, int rank);
tComBuf *bampi_alloc_combuf(void);
void bampi_free_combuf(tComBuf *b);
tComBuf **bampi_alloc_combufs(int n);
void bampi_free_combufs(int n, tComBuf **b);
void bampi_alloc_combuflist(tL *level);
void bampi_free_combuflist(tL *level);
void bampi_free_combuflist_com(tCom *c);
tCom *bampi_alloc_com(void);
void bampi_free_com(tCom *c);
void bampi_set_boundary_flags(tL *l);

/* copy.c */
void bampi_getdata(tL *level, tVarList *vl, 
		   int npoints, double *coords, double *data);
void test_bampi_getdata(tL *level);

/* output.c */
void bampi_combinepointbuffers(tL *level, int nbuf, double *buf, 
			       double *buffer, int nvariables);

/* openmp.c */
void bampi_openmp_initialize(int *pargc, char ***pargv, int pr);
int bampi_openmp_rank();
int bampi_openmp_size();
void bampi_openmp_test();

/* reduce.c */
void bampi_allreduce_max(tVarList *vl, double **result);
void bampi_allreduce_min(tVarList *vl, double **result);
void bampi_allreduce_sum(tVarList *vl, double **result);
void bampi_allreduce_sum_2masks(tVarList *vl, double **result, tVarList *mask);
void bampi_allreduce_sum_mask(tVarList *vl, double **result, tVarList *mask);
void bampi_allreduce_norm1(tVarList *vl, double **result);
void bampi_allreduce_norm2(tVarList *vl, double **result);
void bampi_allreduce_norm2_mask(tVarList *vl, double **result, tVarList *mask);
void bampi_allreduce_norm(tVarList *vl, double **result); 
void bampi_allreduce_normInf(tVarList *vl, double **result);
double bampi_allreduce_allnormInf(tVarList *vl);
double bampi_allreduce_allnorm2(tVarList *vl);
void bampi_allreduce_dot(tVarList *vl, tVarList *wl, double **result);
double bampi_allreduce_alldot(tVarList *vl, tVarList *wl);
void bampi_allreduce_sum_vector(void *local, void *global, int n);
void bampi_allreduce_max_vector(void *local, void *global, int n);
void bampi_allreduce_min_vector(void *local, void *global, int n);
void bampi_allreduce_bbox(tL *level, double *bbox, int *ibbox);
int bampi_allreduce_sum_int(int i);
int bampi_allreduce_max_int(int i);
int bampi_allreduce_min_int(int i);
void bampi_allgather_vector(int nlocal, void *local, void *global);
void bampi_allreduce_maxminpos(tL *level, char *name,
                double *pvarmax, double *pxmax, double *pymax, double *pzmax, 
                double *pvarmin, double *pxmin, double *pymin, double *pzmin);
void bampi_allreduce_maxpos(tL *level, char *name,
		double *pvarmax, double *pxmax, double *pymax, double *pzmax);
void bampi_allreduce_minpos(tL *level, char *name,
                double *pvarmin, double *pxmin, double *pymin, double *pzmin);

/* split_box.c */
void bampi_split_box(int m, int n, int o, int *xs, int *ys, int *zs, 
		     int *mlocal, int *nlocal, int *olocal,
		     int *ilocal, int *jlocal, int *klocal);
void bampi_init_com(tL* level, 
		    int m, int n, int o, int xsize, int ysize, int zsize);
void bampi_check_split(tL *level);

/* synchronize.c */
void bampi_synchronize(tL *level, int vi);
void bampi_synchronize_variables(tL *l, int nv, int *iv);
void bampi_vlsynchronize(tVarList *vl);

/* synchronize_fmr.c */
void bampi_syncparent_init(tL *level);
void bampi_syncparent_send3(tL *level, tVarList *vl);
void bampi_syncparent_recv012(tL *level, tVarList *vl);
void bampi_syncparent_send23(tL *level, tVarList *vl);
void bampi_syncparent_recv0123(tL *level, tVarList *vl);
void bampi_syncparent_send0123(tL *level, tVarList *vl);
void bampi_syncparent_send123(tL *level, tVarList *vl);

/* timer.c */
int timer_start(tL *level, char *name);
int timer_stop(tL *level, char *name);
void timer_print(tG *g);
void freeTimer();












