/* print.c */
/*  */

#include "bam.h"
#include "multigrid.h"




void prvarbasic1(tL *level, char *name, int formattype, int slices)
{
  int n, nn = level->nnodes;
  int i, j, k, max[3];
  double *bbox = level->com->bbox;
  double u, *a, zz;
  double *v = Ptr(level, name);
  double *x = Ptr(level, "x");
  double *y = Ptr(level, "y");
  double *z = Ptr(level, "z");
  char *format, *empty;

  for (i = 0; i < 3; i++)
    max[i] = level->com->ibbox[2*i+1] + 1;
  
  if (formattype == 0) {
    format = "%6.3f";
    empty  = "......";
  } else if (formattype == 1) {
    format = " %2.0f";
    empty  = "  ";
  } else if (formattype == 2) {
    format = "%10.3e";
    empty  = "          ";
  } else 
    errorexit("prvarbasic: unknown format");

    n = max[0]*max[1]*max[2];
    a = (double *) malloc(sizeof(double) * n);
    for (i = 0; i < n; i++)
      a[i] = DBL_MIN;

    for (n = 0; n < nn; n++) {
      i = (x[n] - bbox[0])/level->dx + 0.5;
      j = (y[n] - bbox[2])/level->dy + 0.5;
      k = (z[n] - bbox[4])/level->dz + 0.5;
      a[i + max[0] * (j + max[1] * k)] = v[n];
    }

    for (k = 0; k < max[2]; k++) {
      zz = bbox[4] + k*level->dz;
    
      /* enable if only one plane is to be output */
      if (1) 
        if (dless(zz, 0) || !dless(zz, level->dz)) continue;

      printf("%s l%d  k %d  z %.3f\n", name, level->l, k, zz);

      for (j = max[1]-1; j >= 0; j--) {
        for (i = 0; i < max[0]; i++) {
          u = a[i + max[0] * (j + max[1] * k)];
          if (u == DBL_MIN) 
            printf("%s", empty);
          else 
            printf(format, u); 
        }
        printf("\n");
      }
    }

    free(a);
}

void prvarbasic(tL *level, char *name, int formattype, int slices)
{
  forallboxes(level) {
    if (level->nboxes > 1) printf("box %d of %d\n", nbox, level->nboxes);
    prvarbasic1(one_box_level(level, nbox), name, formattype, slices);
  } endforboxes;
}




void prvar01(tL *level, char *name) {
  prvarbasic(level, name, 1, 0);
}

void prvar(tL *level, char *name) {
  prvarbasic(level, name, 0, 0);
}

void prvare(tL *level, char *name) {
  prvarbasic(level, name, 2, 0);
}







void prprim(char *s, tVarList *vl)
{
  int i;

  for (i = 0; i < vl->n; i++) {
    prdivider(0);
    printf("%s\n", s);
    prvar(vl->level, VarName(vl->index[i]));
  }
}

void prprime(char *s, tVarList *vl)
{
  int i;

  for (i = 0; i < vl->n; i++) {
    prdivider(0);
    printf("%s\n", s);
    prvare(vl->level, VarName(vl->index[i]));
  }
}