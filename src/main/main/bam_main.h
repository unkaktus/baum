/* bam_main.h */
/* Bernd Bruegmann, 12/99 */

/* memory tracing */
// #define TRACEMEMORY

#include <stdlib.h>
#include <stdio.h>

/* constants */
#ifdef PI
#undef PI
#endif
#define PI  3.1415926535897932
#define PIh 1.5707963267948966
#define PIq 0.78539816339744831

/* snap effect for grid coordinates */
#define dequaleps 1e-10
#define dless(a,b) ((a)<(b)-dequaleps)
#define dequal(a,b) (!(dless(a,b)||dless(b,a)))

/* approx. <= and >= with dequaleps tolerance */
#define dlesseq(a,b) ( (a)<(b)+dequaleps )
#define dgreatereq(a,b) ( (a)>(b)-dequaleps )
#define dlessgreater(a,b,c) (dless(a,b)&&dless(b,c))
#define dlessgreatereq(a,b,c) (dlesseq(a,b)&&dlesseq(b,c))

/* parameters.c */
void AddPar(const char *name, const char *value, const char *description);
void AddPari(const char *name, const int i, const char *value, const char *description);
void AppPar(const char *name, const char *value);
int ExistPar(const char *name);
int ExistPari(const char *name, const int i);
void Sets(const char *name, const char *value);
void Seti(const char *name, const int i);
void Setd(const char *name, const double d);
char *Gets(const char *name);
char *GetsLax(const char *name);
int Geti(const char *name);
int GetiLax(const char *name);
double Getd(const char *name);
double GetdLax(const char *name);
int Getv(const char *name, const char *value);
int GetvLax(const char *name, const char *value);
void SetdArray(const char *name, const int n, const double *a);
double *GetdArray(const char *name, int *n);
int GetdArrayN(const char *name);
double GetdEntry(const char *name, const int n);
char *NextEntry(const char *list);
void Appends(const char *name, const char *value);
char *GetsInd(const int i);
char *GetnameInd(const int i);
int GetnParameters();

/* skeleton.c */
enum {
  PRE_GRID,
  PRE_PRE_PRE_INITIALDATA,
  PRE_PRE_INITIALDATA,
  PRE_INITIALDATA,
  INITIALDATA_SET,
  INITIALDATA_FINISH,
  POST_INITIALDATA,
  POST_POST_INITIALDATA,
  PRE_EVOLVE,
  EVOLVE,
  POST_EVOLVE,
  POST_MOVEBOX,
  PRE_ANALYZE,
  ANALYZE,
  POST_ANALYZE,
  OUTPUT,
  POST_OUTPUT,
  NFUNCTIONS
};
void AddFun(int step, int (*f)(tL *), char *name);
void RunFun(int step, tL *level);

/* system.c */
int system_monitor(tL *level);

/* tensors.c */
#define NINDEXLIST 100
void tensorindexlist(char *tensorindices, int *nilist, char **ilist, int *sym);

/* utilities.c */
extern tVarList *vl_analyze;

#ifdef TRACEMEMORY
#define malloc bam_malloc
#define calloc bam_calloc
#define realloc bam_realloc
#define free bam_free
#ifdef strdup
#undef strdup
#endif
#define strdup bam_strdup
#endif
void prmemory_on(void);
void prmemory_off(void);
void prmemory(char *s);
void prmemory_adv(char *s, tG *g);
void prmemory_active(void);
void hidepointer(void *ptr);
void *bam_malloc(const size_t size);
void *bam_calloc(const size_t nmemb, const size_t size);
void *bam_realloc(void *ptr, const size_t size);
void bam_free(void *ptr);
char *bam_strdup(const char *cptr);

void errorexit(const char *file, const int line, const char *s);
void errorexits(const char *file, const int line, const char *s, const char *t);
void errorexiti(const char *file, const int line, const char *s, const int i);
#define errorexit(s) errorexit(__FILE__, __LINE__, (s))
#define errorexits(s,t) errorexits(__FILE__, __LINE__, (s), (t))
#define errorexiti(s,i) errorexiti(__FILE__, __LINE__, (s), (i))

void yo(), yo1(), yo2(), yo3(), yo4(), yo5(), yo6(), yo7(), yo8(), yo9(), yob();
void prdivider(const int n);
int aprintf(char *target, const char *format, ...);
double min2(const double x, const double y);
double min3(const double x, const double y, const double z);
double max2(const double x, const double y);
double max3(const double x, const double y, const double z);
int system2(const char *s1, const char *s2);
int system3(const char *s1, const char *s2, const char *s3);
int construct_argv(char *str, char ***argv);
int system_emu(const char *command);
int lock_curr_til_EOF(FILE *out);
int unlock_curr_til_EOF(FILE *out);
int system_isfile(const char *s1);
int system_isdir(const char *s1);
int system_mkdir(const char *s);
int system_cp(const char *s1, const char *s2);
int system_cpr(const char *s1, const char *s2);
int system_move(const char *s1, const char *s2);
int system_rmrf(const char *s);
int system_chmod(const char *s, const char *m);

double *dmalloc(const int n);
double **ddmalloc(const int n1, const int n2);
void ddfree(double **p, const int n);
double *dcalloc(const int n);
int *imalloc(const int n);
int *icalloc(const int n);
char *cmalloc(const int n);
void *pmalloc(const int n);
void *bcalloc(const int n);

/* variables.c */
#define PtrFromInd(l,i)     ((l)->v[(i)])
#define VLPtr(vl,i)         ((i)<(vl)->n?(vl)->level->v[(vl)->index[i]]:0)
#define vldataptr(vl,i)     ((i)<(vl)->n?(vl)->level->v[(vl)->index[i]]:0)
#define vlldataptr(vl,l,i)  ((l)->v[(vl)->index[i]])
extern int globalnvariables;
int Ind(const char *name);
int IndLax(const char *name);
int Indvl(const char *name, const tVarList *vl);
void AddVar(const char *name, char *indices, const char *description);
void AddConstantVar(const char *name, char *tensorindices,
		    const char *description);
tVarList *AddDuplicate(const tVarList *vl, const char *postfix);
tVarList *AddDuplicateEnable(const tVarList *vl, const char *postfix);

char *VarName(const int i);
int VarComponent(const int i);
int VarNComponents(const int i);
int IndComponent0(const int i);
char *VarNameComponent0(const char *name);
char *VarTensorIndices(const int i);
void VarNameSetBoundaryInfo(const char *name, const double farlimit,
			    const double falloff, const double propspeed);
double VarFallOff(const int i);
double VarFarLimit(const int i);
double VarPropSpeed(const int i);
int VarGaugeFlag(const int i);
void VarNameSetConstantFlag(const char *name);
int VarConstantFlag(const int i);
int VarSymmetry(const int i, const int dir);
void VarSetSymmetry(const int i, const int dir, const int sym);
void VarInvertSymmetry(const int i);

void prvarlist(const tVarList *v);
tVarList *vlalloc(tL *level);
void vlenable(const tVarList *v);
void vldisable(const tVarList *v);
void vlenablelevel(tL *level, tVarList *v);
void vlfree(tVarList *u);
void freevll();
void vlpushone(tVarList *v, const int vi);
void vlpush(tVarList *v, const int vi);
void vlpushvl(tVarList *v, const tVarList *u);
void vldropone(tVarList *v, const int vi);
void vldrop(tVarList *v, const int vi);
void vldropn(tVarList *v, const int n);
tVarList *vlclean(const tVarList *v);
tVarList *vlduplicate(const tVarList *v);
void vlsetconstant(const tVarList *u, const double c);
void vlcopy(const tVarList *v, const tVarList *u);
void vlcopylevel(tL *level, tVarList *v, tVarList *u);
void vlaverage(const tVarList *r, const tVarList *a, const tVarList *b);
void vlsubtract(const tVarList *r, const tVarList *a, const tVarList *b);
void vladd(const tVarList *r, const double ca, const tVarList *a,
	   const double cb, const tVarList *b);
void vladdto(const tVarList *r, const double ca, const tVarList *a);
void vlswap(const tVarList *u, const tVarList *v);


tVarList *VLPtrEnable1(tL *level, const char *varname);
void VLDisableFree(tVarList *vl);

void SetSymmetry(int i, int dir, int sym);

















