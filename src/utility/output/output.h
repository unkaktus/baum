/* output.h */

#include <sys/stat.h>
#include <sys/types.h>

extern char boxstring[2];
extern int boxnr;

extern int output_order;
extern int output_scheme;


/* output0d.c */
void write_level_0d(tL *level, int nv, int *iv, int type, char *suffix);
void write_finestlevel_point(tL *level, int nv, int *iv, int type, char *suffix);

/* output1d.c */
void write_level_1d(tL *level, int nv, int *iv, int line, char *suffix);

/* output2d.c */
void write_level_2d(tL *level, int nv, int *iv, int plane, char *suffix);

/* output3d.c */
void write_level_3d(tL *level, int nv, int *iv, int sampling, char *suffix);

/* outputShells.c */
void write_level_shells1d(tL* level,int nv, int *iv, int l, char *suffix);
void write_level_shells2d(tL* level,int nv, int *iv, int l, char *suffix);
void write_level_shells3d(tL* level,int nv, int *iv, char *suffix);

/* outputVTK.c */
FILE *fopen_vtk(char *varname, char *dirsuffix, char *suffix, int l, int n);
void write_raw_vtk(FILE *fp, int n, double *buffer, int nv, int j,
                   int dbl, int flt, int text, int binary);
void write_raw_vec_vtk(FILE *fp, int n, double *buffer,
                       int dbl, int flt, int text, int binary);

/* outputXDMF.c */
void write_xdmf3D(char *varname, char *dirsuffix, char *suffix, int level, double t, 
        int n, double *buffer, int dbl, int flt, int text, int binary, 
        int nx, int ny, int nz, double x0, double y0, double z0, double dx, double dy, double dz);
void write_xdmf3D_vec(char *varname, char *dirsuffix, char *suffix, int level, double t,
        int n, double *buffer, int dbl, int flt, int text, int binary,
        int nx, int ny, int nz, double x0, double y0, double z0, double dx, double dy, double dz);
void write_xdmf3D_rectilinear(char *varname, char *dirsuffix, char *suffix, int level, double t,
        int n, double *buffer, int dbl, int flt, int text, int binary,
        int nx, int ny, int nz, double x0, double y0, double z0, double dx, double dy, double dz);
void write_xdmf3D_rectilinear_vec(char *varname, char *dirsuffix, char *suffix, int level, double t,
        int n, double *buffer, int dbl, int flt, int text, int binary,
        int nx, int ny, int nz, double x0, double y0, double z0, double dx, double dy, double dz);
void write_xdmf3D_curvilinear(char *varname, char *dirsuffix, char *suffix, int level, double t,
        int n, double *buffer, int dbl, int flt, int text, int binary,
        int nx, int ny, int nz, double x0, double y0, double z0, double dx, double dy, double dz);
void write_xdmf3D_curvilinear_vec(char *varname, char *dirsuffix, char *suffix, int level, double t,
        int n, double *buffer, int dbl, int flt, int text, int binary,
        int nx, int ny, int nz, double x0, double y0, double z0, double dx, double dy, double dz);
void write_xdmf2D(char *varname, char *dirsuffix, char *suffix, int level, double t, 
        int n, double *buffer, int nv, int j, int dbl, int flt, int text, int binary, 
        int nx, int ny, double x0, double y0, double dx, double dy);
void write_xdmf2D_vec(char *varname, char *dirsuffix, char *suffix, int level, double t, 
        int n, double *buffer, int nv, int j, int dbl, int flt, int text, int binary, 
        int nx, int ny, double x0, double y0, double dx, double dy, int k, int l);

