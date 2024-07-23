/* checkpoint.h */
/* Bernd Bruegmann 01/2006, Wolfgang Tichy 11/2006 */


/* checkpoint.c */
char *checkpoint_filename(char *suffix);
FILE *checkpoint_openfiles(char *suffix);
tVarList *checkpoint_varlist(tL *level);
void checkpoint_copy_output(double time);

void checkpoint_read(tG *g);
void checkpoint_read_ParsAndIterations_local(tG *g, FILE *fp);
void checkpoint_read_Vars_local(tG *g, FILE *fp);

void checkpoint_write(tG *g);
void checkpoint_write_local(tG *g, FILE *fp);

int checkpoint_init(tL *level);
