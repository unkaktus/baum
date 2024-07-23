/* poisson.c */
/* Bernd Bruegmann 01/00 */

/* test elliptic solvers with Poisson equation 
   adapted from:
*/

/* thorn_BAM_Elliptic: Multigrid Elliptic Solver for Cactus
   (C) 1998 Bernd Bruegmann 
*/
/* poisson.c */
/* Bernd Bruegmann, 4/97 */

#include "bam.h"
#include "iterative.h"


#define EPOL 0
int DPflag = 0;



int poisson(tL *level) 
{
  errorexit("BOX: implement poisson, see older BAM version");
  return 0;
}



  
/* diagonal preconditioning:   v = u/Lii */
void DPflatlinear(tL *level, tVarList *v, tVarList *u)
{

}




/* apply elliptic operator:  lu = Laplace^flat u */
void Lflatlinear(tL *level, tVarList *vllu, tVarList *vlu)
{

  errorexit("BOX: implement Lflatlinear, see older BAM version!");

}


