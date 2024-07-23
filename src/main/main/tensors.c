/* tensors.c */
/* Bernd Bruegmann 6/02, 11/02 */

#include "bam.h"
#include <ctype.h>



/* helper function: 
   decode tensor index string
   return list of indices into variable list 
   return sign under reflections for symmetry boundaries
   note that we treat 3d indices i,j,k,... and 4d indices a,b,c,...

   should be made automatic, but for now this is simpler 
*/
void tensorindexlist(char *t, int *nilist, char **ilist, int *sym)
{
  /* name of coordinates, could be made variable */
  char *coord[3]  = {"x", "y", "z"};
  char *coord4[4] = {"t", "x", "y", "z"};
  int i, j, k, l;
  int n = 0;
  char *tensorindices = strdup(t);

  /* convert local copy to lower case since we ignore co/contra-variance */
  for (t = tensorindices; *t; t++)
    *t = tolower(*t);

  /* initialize symmetries */
  for (i = 0; i < 3*NINDEXLIST; i++) 
    sym[i] = 1;


  /* now treat each case separately */

  /* scalar */
  if (strcmp(tensorindices, "") == 0) {
    ilist[n] = calloc(sizeof(char), 8);
    // sprintf(ilist[n++], ""); // gives gcc warning
    // *ilist[n++] = 0;         // not necessary since calloc() initializes
    n++;
  }
  
  /* 3d indices */
  if (strcmp(tensorindices, "i") == 0) {
    for (i = 0; i < 3; i++) {
      sym[3*n+i] *= -1;
      ilist[n] = calloc(sizeof(char), 8);
      sprintf(ilist[n++], "%s", coord[i]);
    }
  }

  if (strcmp(tensorindices, "ij") == 0) {
    for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++) {
      sym[3*n+i] *= -1; 
      sym[3*n+j] *= -1;
      ilist[n] = calloc(sizeof(char), 8);
      sprintf(ilist[n++], "%s%s", coord[i], coord[j]);
    }
  }

  if(strcmp(tensorindices, "ij+ji") == 0 ||
     strcmp(tensorindices, "(ij)" ) == 0) {
    for (i = 0; i < 3; i++)
    for (j = i; j < 3; j++) {
      sym[3*n+i] *= -1; 
      sym[3*n+j] *= -1;
      ilist[n] = calloc(sizeof(char), 8);
      sprintf(ilist[n++], "%s%s", coord[i], coord[j]);
    }
  }
  
  if (strcmp(tensorindices, "ij-ji") == 0) {
    for (i = 0; i < 3; i++)
    for (j = i+1; j < 3; j++) {
      sym[3*n+i] *= -1; 
      sym[3*n+j] *= -1;
      ilist[n] = calloc(sizeof(char), 8);
      sprintf(ilist[n++], "%s%s", coord[i], coord[j]);
    }
  }
  
  if (strcmp(tensorindices, "ijk") == 0) {
    for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++) 
    for (k = 0; k < 3; k++) {
      sym[3*n+i] *= -1; 
      sym[3*n+j] *= -1;
      sym[3*n+k] *= -1;
      ilist[n] = calloc(sizeof(char), 8);
      sprintf(ilist[n++], "%s%s%s", coord[i], coord[j], coord[k]);
    }
  }
  
  if (strcmp(tensorindices, "ijk+ikj") == 0) {
    for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++) 
    for (k = j; k < 3; k++) {
      sym[3*n+i] *= -1; 
      sym[3*n+j] *= -1;
      sym[3*n+k] *= -1;
      ilist[n] = calloc(sizeof(char), 8);
      sprintf(ilist[n++], "%s%s%s", coord[i], coord[j], coord[k]);
    }
  }
  
  if (strcmp(tensorindices, "ijk+jik") == 0) {
    for (i = 0; i < 3; i++)
    for (j = i; j < 3; j++) 
    for (k = 0; k < 3; k++) {
      sym[3*n+i] *= -1; 
      sym[3*n+j] *= -1;
      sym[3*n+k] *= -1;
      ilist[n] = calloc(sizeof(char), 8);
      sprintf(ilist[n++], "%s%s%s", coord[i], coord[j], coord[k]);
    }
  }
  
  if (strcmp(tensorindices, "ijk-jik") == 0) {
    for (i = 0; i < 3; i++)
    for (j = i+1; j < 3; j++) 
    for (k = 0; k < 3; k++) {
      sym[3*n+i] *= -1; 
      sym[3*n+j] *= -1;
      sym[3*n+k] *= -1;
      ilist[n] = calloc(sizeof(char), 8);
      sprintf(ilist[n++], "%s%s%s", coord[i], coord[j], coord[k]);
    }
  }
  
  if (strcmp(tensorindices, "ijkl") == 0) {
    for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++) 
    for (k = 0; k < 3; k++)
    for (l = 0; l < 3; l++)
    {
      sym[3*n+i] *= -1; 
      sym[3*n+j] *= -1;
      sym[3*n+k] *= -1;
      sym[3*n+l] *= -1;
      ilist[n] = calloc(sizeof(char), 8);
      sprintf(ilist[n++], "%s%s%s%s", coord[i], coord[j], coord[k], coord[l]);
    }
  }

  if(strcmp(tensorindices, "ijkl+ijlk+jikl+jilk)") == 0 ||
     strcmp(tensorindices, "(ij)(kl)") == 0) {
    for (i = 0; i < 3; i++)
    for (j = i; j < 3; j++) 
    for (k = 0; k < 3; k++)
    for (l = k; l < 3; l++)
    {
      sym[3*n+i] *= -1; 
      sym[3*n+j] *= -1;
      sym[3*n+k] *= -1;
      sym[3*n+l] *= -1;
      ilist[n] = calloc(sizeof(char), 8);
      sprintf(ilist[n++], "%s%s%s%s", coord[i], coord[j], coord[k], coord[l]);
    }
  }

  /* 4d indices */
  if (strcmp(tensorindices, "a") == 0) {
    for (i = 0; i <= 3; i++) {
      if (i > 0) sym[3*n+i-1] *= -1;
      ilist[n] = calloc(sizeof(char), 8);
      sprintf(ilist[n++], "%s", coord4[i]);
    }
  }
  
  if (strcmp(tensorindices, "ab") == 0) {
    for (i = 0; i <= 3; i++)
    for (j = 0; j <= 3; j++) {
      if (i > 0) sym[3*n+i-1] *= -1; 
      if (j > 0) sym[3*n+j-1] *= -1;
      ilist[n] = calloc(sizeof(char), 8);
      sprintf(ilist[n++], "%s%s", coord4[i], coord4[j]);
    }
  }

  if (strcmp(tensorindices, "ai") == 0) {
    for (i = 0; i <= 3; i++)
    for (j = 0; j <  3; j++) {
      if (i > 0) sym[3*n+i-1] *= -1; 
      sym[3*n+j] *= -1;
      ilist[n] = calloc(sizeof(char), 8);
      sprintf(ilist[n++], "%s%s", coord4[i], coord[j]);
    }
  }

  if (strcmp(tensorindices, "ab+ba") == 0) {
    for (i = 0; i <= 3; i++)
    for (j = i; j <= 3; j++) {
      if (i > 0) sym[3*n+i-1] *= -1; 
      if (j > 0) sym[3*n+j-1] *= -1;
      ilist[n] = calloc(sizeof(char), 8);
      sprintf(ilist[n++], "%s%s", coord4[i], coord4[j]);
    }
  }
  
  if (strcmp(tensorindices, "abc+acb") == 0) {
    for (i = 0; i <= 3; i++)
    for (j = 0; j <= 3; j++) 
    for (k = j; k <= 3; k++) {
      if (i > 0) sym[3*n+i-1] *= -1; 
      if (j > 0) sym[3*n+j-1] *= -1;
      if (k > 0) sym[3*n+k-1] *= -1;
      ilist[n] = calloc(sizeof(char), 8);
      sprintf(ilist[n++], "%s%s%s", coord4[i], coord4[j], coord4[k]);
    }
  }

  if(strcmp(tensorindices, "aij+aji") == 0 ||
     strcmp(tensorindices, "a(ij)") == 0) {
    for (i = 0; i <= 3; i++)
    for (j = 0; j <  3; j++) 
    for (k = j; k <  3; k++) {
      if (i > 0) sym[3*n+i-1] *= -1; 
      sym[3*n+j] *= -1;
      sym[3*n+k] *= -1;
      ilist[n] = calloc(sizeof(char), 8);
      sprintf(ilist[n++], "%s%s%s", coord4[i], coord[j], coord[k]);
    }
  }

  if(strcmp(tensorindices, "abc+bac") == 0 ||
     strcmp(tensorindices, "(ab)c")   == 0) {
    for (i = 0; i <= 3; i++)
    for (j = i; j <= 3; j++) 
    for (k = 0; k <= 3; k++) {
      if (i > 0) sym[3*n+i-1] *= -1; 
      if (j > 0) sym[3*n+j-1] *= -1;
      if (k > 0) sym[3*n+k-1] *= -1;
      ilist[n] = calloc(sizeof(char), 8);
      sprintf(ilist[n++], "%s%s%s", coord4[i], coord4[j], coord4[k]);
    }
  }

  if(strcmp(tensorindices, "abi+bai")   == 0 ||
     strcmp(tensorindices, "(ab)i")   == 0) {
    for (i = 0; i <= 3; i++)
    for (j = i; j <= 3; j++) 
    for (k = 0; k <  3; k++) {
      if (i > 0) sym[3*n+i-1] *= -1; 
      if (j > 0) sym[3*n+j-1] *= -1;
      sym[3*n+k] *= -1;
      ilist[n] = calloc(sizeof(char), 8);
      sprintf(ilist[n++], "%s%s%s", coord4[i], coord4[j], coord[k]);
    }
  }

  if(strcmp(tensorindices, "abij+abji+baij+baji)") == 0 ||
     strcmp(tensorindices, "(ab)(ij)") == 0) {
    for (i = 0; i <= 3; i++)
    for (j = i; j <= 3; j++) 
    for (k = 0; k < 3; k++)
    for (l = k; l < 3; l++)
    {
      if (i > 0) sym[3*n+i] *= -1; 
      if (j > 0) sym[3*n+j] *= -1;
      sym[3*n+k] *= -1;
      sym[3*n+l] *= -1;
      ilist[n] = calloc(sizeof(char), 8);
      sprintf(ilist[n++], "%s%s%s%s", coord4[i], coord4[j], coord[k], coord[l]);
    }
  }

  /* special case: 3d and 4d axial vectors */
  if (strcmp(tensorindices, "i-axial") == 0) {
    for (i = 0; i < 3; i++) {
      for (j = 0; j < 3; j++)
        if (j != i) sym[3*n+j] *= -1;
      ilist[n] = calloc(sizeof(char), 8);
      sprintf(ilist[n++], "%s", coord[i]);
    }
  }

  if (strcmp(tensorindices, "a-axial") == 0) {
    for (i = 0; i <= 3; i++) {
      if (i > 0) 
        for (j = 1; j < 3; j++)
          if (j != i) sym[3*n+j-1] *= -1;
      ilist[n] = calloc(sizeof(char), 8);
      sprintf(ilist[n++], "%s", coord4[i]);
    }
  }
  
  /* special case: 3d vectors all positive*/
  if (strcmp(tensorindices, "i-positive") == 0) {
    for (i = 0; i < 3; i++) {
      for (j = 0; j < 3; j++)
        if (j != i) sym[3*n+j] *= +1;
      ilist[n] = calloc(sizeof(char), 8);
      sprintf(ilist[n++], "%s", coord[i]);
    }
  }

  /* special case: scalar all negative */
  if (strcmp(tensorindices, "pseudo") == 0) {
    ilist[n] = calloc(sizeof(char), 8);
    for (j = 0; j < 3; j++)
      sym[3*n+j] *= -1;
    n++;
  }

  /* error */
  if (n == 0) {
    printf("Error in index string %s.\n", tensorindices);
    printf("Legal combinations besides the empty string are\n");
    printf("i, ij, ij+ji, ij-ji, ijk, ijk+ikj, ijk+jik\n");
    printf("a, ab, ab+ba, abc+acb\n");
    printf("(ab)c, (ab)i, (ab)(ij)\n");
    printf("i-axial, a-axial, pseudo\n");
    printf("Anything else can be easily added to main/tensors.c.\n");
    errorexit("");
  }

  *nilist = n;
  free(tensorindices);
}




