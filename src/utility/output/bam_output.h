/* bam_output.h */
/* Bernd Bruegmann 12/99 */

    
/* this value marks points in the output for which no data is available */
#define NODATA (1.234e+56)




/* output.c */
int timeforoutput(tL *level, tVarList *vl);
int timeforoutput_index(tL *level, int index);
int timeforoutput_any(tL *level);
int timeforoutput_di_dt(tL *level, int di, double dt);
int timeforoutput_dt(tL *level, double dt);
void makeoutputlist(tL *level, char *out0, int allflag,
                    int *pnindex, int **pindex);


int write_level(tL *level);
void init_output();


void write_sphere(tL *level, int ncircles, int npoints, double r, char* suffix, int index);
void write_level_sphere(tL *level, int numvars, int *indexarray, 
                        int ncircles, int npoints, double r, char *suffix);


/* outputShells.c */
void write_level_shells3d_surface(tL* level, int nv, int *iv);



