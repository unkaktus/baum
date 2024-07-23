/* skeleton.c */
/* Bernd Bruegmann 12/99 */

/* better name? fps function pointer skeleton */

#include "bam.h"
#include "main.h"

#define PR 0


typedef struct tTODO {
  struct tTODO *next;
  int (*f)(tL *);
  char *name;
} tTodo;

tTodo *fps[NFUNCTIONS];




void AddFun(int step, int (*f)(tL *), char *name)
{
  tTodo *t;

  if (1) printf("  function  %-30s  =  (ENUM %d)\n", name,step);
  
  if (!fps[step]) fps[step] = (tTodo *) calloc(sizeof(tTodo), 1);
  
  for (t = fps[step]; t->next; t = t->next);
  t->next = (tTodo *) calloc(sizeof(tTodo), 1);
  t->f = f;
  t->name = (char *) calloc(sizeof(char), strlen(name)+1);
  strcpy(t->name, name);
  
}




void RunFun(int step, tL *level) 
{
  tTodo *t;
  
  if (!fps[step]) return;

  for (t = fps[step]; t->next; t = t->next) {
    if (PR)
      printf("======== fps l%d %s ========\n", (level)?level->l:0, t->name);
    if (0) prmemory(t->name);
    (*(t->f))(level);
    if (0) prmemory(t->name);
  }
}




void freeFuns()
{
  int i,j;
  tTodo *t;
  
  for (i=0; i<NFUNCTIONS; i++) {
    if (fps[i]) {
      do {
        j = 0;
        for (t = fps[i]; t->next->next; t = t->next) j++;
        if (0) printf("freeFun:  %d  %d  %p %p %p %p\n",i,j,fps[i],t,t->name,t->next);
        free(t->next);
        t->next = 0;
        free(t->name);
      } while (t!=fps[i]);
    }
    free(fps[i]);
  }
}




