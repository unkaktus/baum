/* physics.c */
/* mth, 04/12 */

#include "bam.h"
#include "main.h"






/* initialize Object data */
void init_physical_objects(tG* grid) 
{

  printf("Initializing Objects (given by parfile, used by Initialdata project)\n");
  int pr = 1;
  int i,j;
  char str[100];
  
  for (i=1; i<=Geti("nobjects"); i++) {
    sprintf(str,"mass%d",i);
    if (Getd(str)==0) break;
    
    if (pr) {
      printf("  Objects%d   = \n",i);
      sprintf(str,"mass%d",i);
      printf("    mass      = %+2.4e\n",Getd(str));
      printf("    position  = ");
      for (j=0; j<3; j++) {
        sprintf(str,"p%c%d",(j==0)?'x':((j==1)?'y':'z'),i);
        printf(" %+2.4e",Getd(str));
      }
      printf("\n");
      printf("    momentum  = ");
      for (j=0; j<3; j++) {
        sprintf(str,"m%c%d",(j==0)?'x':((j==1)?'y':'z'),i);
        printf(" %+2.4e",Getd(str));
      }
      printf("\n");
      printf("    spin      = ");
      for (j=0; j<3; j++) {
        sprintf(str,"s%c%d",(j==0)?'x':((j==1)?'y':'z'),i);
        printf(" %+2.4e",Getd(str));
      }
      printf("\n");
    }
    
    // set values to global structure which will see the amr
    for (j=0; j<3; j++) {
      sprintf(str,"p%c%d",(j==0)?'x':((j==1)?'y':'z'),i);
      grid->puncpos[i-1][j] = Getd(str);
    }
  }
  Seti("amr_npunctures", Geti("nobjects"));
  
  if (Geti("amr_npunctures")==0) {
    if (Geti("amr_lmax") == 0)
      Seti("amr_npunctures",1);
    else if (Getv("physics", "srhdtestID"))
      Seti("amr_npunctures",1);
    else
      errorexit("  There is no BH or NS on the grid -> I better stop, look at the parfile");
  }
  
  prdivider(0);
  grid->npunc = Geti("amr_npunctures");
  printf("Adjusting grid\n");
  printf("  Change to %d amr puncture(s) at:\n",grid->npunc);
  for (i=0; i<grid->npunc; i++) {
    printf("  p%d:   x=%+2.4f y=%+2.4f z=%+2.4f   lmax=%d\n",i+1,
           grid->puncpos[i][0],grid->puncpos[i][1],grid->puncpos[i][2],
           grid->lmaxpunc[i]);
  }
 
  
}






