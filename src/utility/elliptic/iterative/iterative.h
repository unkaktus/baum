/* iterative.h */
/* Bernd Bruegmann 01/00 */


int poisson(tL *level);
void DPflatlinear(tL *level, tVarList *v, tVarList *u);
void Lflatlinear(tL *level, tVarList *lu, tVarList *u);


#define forvector(i) for (i = 0; i < nvector; i++)
extern int nvector;
extern int DPflag;

double dot(tVarList *u, tVarList *v);
double norm2(tVarList *u);

void minus(tVarList *vlw, tVarList *vlu, tVarList *vlv);
