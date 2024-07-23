/* parameters.c */
/* Bernd Bruegmann 12/99 */

/* basic parameter file syntax:
   parameter_name = parameter_value parameter_value ...

   special characters:
   #: this and the rest of the line is ignored
   ": counts as white space

   before interpretation the file content is converted to the form
   par=val ... val par=val ... val ...
*/

#include <ctype.h>
#include "bam.h"

/* parameter data base */
typedef struct
{
  char *name;
  char *value;
  char *description;
} tParameter;

tParameter *pdb;
int npdb, npdbmax = 1000;

tParameter *findparameter(const char *name, const int fatal);
void setparameter(const char *name, const char *value);
void makeparameter(const char *name, const char *value, const char *description);
void printparameters(void);
void translatevalue(char **value);

/* parse given parameter file */
void parse_parameter_file(const char *parfile)
{
  FILE *fp;
  int c, i, j;
  int nbuffer;
  char *buffer;
  char *par, *val;
  int lpar, lval;

  if (1)
    printf("Reading parameter file \"%s\"\n", parfile);

  /* read file into memory, add one space at end and beginning */
  if (processor0)
  {
    fp = fopen(parfile, "r");
    if (!fp)
    {
      printf("Could not open parameter file \"%s\"\n", parfile);
      errorexit("");
    }
    buffer = 0;
    for (i = nbuffer = 0;; i++)
    {
      if (i >= nbuffer - 2)
      {
        if (nbuffer > 1000000)
          errorexit("Sanity forbids parameter files bigger than 1MB");
        buffer = (char *)realloc(buffer, sizeof(char) * (nbuffer += 1000));
        if (!buffer)
          errorexit("Out of memory while reading parameter file.");
      }
      if (i == 0)
        buffer[i++] = ' ';
      if ((c = fgetc(fp)) == EOF)
        break;
      buffer[i] = c;
    }
    fclose(fp);
    buffer[i++] = ' ';
    buffer[i] = '\0';
    nbuffer = strlen(buffer);
  }
  bampi_bcast_int(&nbuffer, 1, 0);
  if (!processor0)
  {
    buffer = malloc(sizeof(char) * (nbuffer + 1));
  }
  bampi_bcast_char(buffer, nbuffer + 1, 0);
  if (0)
    printf("%s", buffer);

  /* white out comments and quotes */
  for (i = 0; i < nbuffer; i++)
  {
    if (buffer[i] == '#')
      while (i < nbuffer && buffer[i] != '\n')
        buffer[i++] = ' ';
    if (buffer[i] == '"')
      buffer[i] = ' ';
  }

  /* collapse all white space to single space */
  for (i = j = 1; i < nbuffer; i++)
  {
    if (!isspace(buffer[i]))
      buffer[j++] = buffer[i];
    else if (!isspace(buffer[j - 1]))
      buffer[j++] = ' ';
  }
  buffer[j] = '\0';
  nbuffer = strlen(buffer);
  if (0)
    printf("|%s|\n", buffer);

  /* and remove spaces around = */
  for (i = j = 1; i < nbuffer; i++)
  {
    if (buffer[i] != ' ' || buffer[i - 1] != '=' && buffer[i + 1] != '=')
      buffer[j++] = buffer[i];
  }
  buffer[j - 1] = '\0';
  nbuffer = strlen(buffer);
  if (0)
    printf("|%s|\n", buffer);

  /* now buffer is
     | par=val ... val par=val ... val|
  */

  /* split parameter names and values by inserting zeros */
  for (i = 1; i < nbuffer; i++)
  {
    if (buffer[i] == '=')
    {
      buffer[i] = '\0';
      for (j = i - 1; j > 0 && buffer[j - 1] != ' '; j--)
        ;
      buffer[j - 1] = '\0';
    }
  }

  /* loop over all parameter/value pairs */
  for (i = 1; i < nbuffer; i += lpar + lval + 2)
  {
    par = buffer + i;
    lpar = strlen(par);
    val = par + lpar + 1;
    lval = strlen(val);
    if (0)
      printf("%s = |%s|\n", par, val);

    if (!findparameter(par, 0))
      makeparameter(par, val, "not in libs");
    else
      setparameter(par, val);
  }

  /* print parameters */
  if (0)
  {
    printf("after reading the parameterfile:\n");
    printparameters();
  }
  free(buffer);

  /* for debugging, this may be needed right away so let's do it here */
  AddPar("stdout_flush", "no",
         "for debugging only: turn off buffering for stdout and stderr [no]");
  if (Getv("stdout_flush", "yes"))
  {
    setvbuf(stdout, 0, _IONBF, 0);
    setvbuf(stderr, 0, _IONBF, 0);
  }

  /* redirection of stdout on process 0 is determined by parameter file,
     so this is the first place we can do that
  */
  AddPar("stdout_redirect", "no",
         "redirect stdout from process 0 [no]");
  if (Getv("stdout_redirect", "yes"))
  {
    if (bampi_rank() == 0)
    {
      char s[1000], f[100];
      sprintf(f, "%%s/stdout.%%0%dd", (int)log10(bampi_size()) + 1);
      sprintf(s, f, Gets("outdir"), bampi_rank());
      freopen(s, "w", stdout);
      freopen(s, "w", stderr);
      prdivider(0);
      printf("Welcome to bam.\n");
      prdivider(0);
      printf("Parameterfile = %s\n", Gets("parameterfile"));
    }
  }

  AddPar("stdout_verbose", "dtdxdiss outtime ",
         "[no,dtdxdiss,outtime ]");
}

/***************************************************************************/
/* parameter data base */

/* make new parameter in parameter data base, merge if already there */
void makeparameter(const char *name, const char *value, const char *description)
{
  static int firstcall = 1;
  tParameter *p;

  if (0)
    printf("Makp %s = %s,  %s\n", name, value, description);

  if (firstcall)
  {
    firstcall = 0;
    pdb = (tParameter *)calloc(sizeof(tParameter), npdbmax);
    if (!pdb)
      errorexit("makeparameter: out of memory");
    npdb = 0;
  }

  p = findparameter(name, 0);

  if (!p)
  {
    p = &pdb[npdb++];
    p->name = (char *)calloc(sizeof(char), strlen(name) + 1);
    p->value = (char *)calloc(sizeof(char), strlen(value) + 1);
    strcpy(p->name, name);
    strcpy(p->value, value);
    translatevalue(&p->value);
  }
  else
  {
    free(p->description);
  }
  p->description = (char *)calloc(sizeof(char), strlen(description) + 1);
  strcpy(p->description, description);

  if (npdb >= npdbmax)
    errorexit("makeparameter: lazy coding, no more space for new parameters");

  if (0)
    printparameters();
}

/* give memory free */
void freeParameters()
{
  int i;
  for (i = 0; i < npdb; i++)
  {
    if (0)
      printf("delete %d %s\n", i, pdb[i].name);
    free(pdb[i].name);
    free(pdb[i].value);
    free(pdb[i].description);
  }
  if (0)
    printf("%d\n", npdb);
  free(pdb);
}

/* set parameter */
void setparameter(const char *name, const char *value)
{
  tParameter *p = findparameter(name, 1);

  free(p->value);
  p->value = strdup(value);
  translatevalue(&p->value);
}

/* find parameter */
tParameter *findparameter(const char *name, const int fatal)
{
  int i;

  if (!name)
  {
    errorexit("findparameter: called without parameter name");
  }

  for (i = 0; i < npdb; i++)
    if (!strcmp(pdb[i].name, name))
      return &pdb[i];

  if (fatal)
  {
    printf("Could not find parameter \"%s\"\n", name);
    errorexit("this one is required!");
  }
  return 0;
}

/* translate parameter value
   should be much more elaborate, say, implement simple arithmetic
*/
void translatevalue(char **value)
{
  double x = 0;

  if (strcmp(*value, "pi") == 0)
    x = PI;
  if (strcmp(*value, "-pi") == 0)
    x = -PI;
  if (strcmp(*value, "pi/2") == 0)
    x = PI / 2;
  if (strcmp(*value, "-pi/2") == 0)
    x = -PI / 2;
  if (strcmp(*value, "2*pi") == 0)
    x = 2 * PI;
  if (strcmp(*value, "-2*pi") == 0)
    x = -2 * PI;

  if (x)
  {
    char newvalue[100];
    sprintf(newvalue, "%.18e", x);
    free(*value);
    *value = strdup(newvalue);
  }
}

/* print parameters */
void printparameters(void)
{
  int i;

  for (i = 0; i < npdb; i++)
    printf("pdb[%3d]:  %16s = %-16s,  %s\n",
           i, pdb[i].name, pdb[i].value, pdb[i].description);
}

/***************************************************************************/
/* functions for external calls */

/* creation functions */
void AddPar(const char *name, const char *value, const char *description)
{
  makeparameter(name, value, description);
  printf("  parameter %-30s  =  %s\n", name, Gets(name));
}

void AddPari(const char *name, const int i, const char *value, const char *description)
{
  char new[100];
  sprintf(new, "%s%d", name, i);
  makeparameter(new, value, description);
  printf("  parameter %-30s  =  %s\n", new, Gets(new));
}

void AppPar(const char *name, const char *value)
{
  if (!*value)
    return;

  char *v;
  if (findparameter(name, 0) == 0)
    makeparameter(name, value, "");
  else
  {
    while ((v = NextEntry(value)))
    {
      Appends(name, v);
    }
  }
  printf("  parameter %-30s  =  %s\n", name, Gets(name));
}

int ExistPar(const char *name)
{
  tParameter *p = findparameter(name, 0);
  return p ? 1 : 0;
}

int ExistPari(const char *name, int i)
{
  char new[100];
  sprintf(new, "%s%d", name, i);
  tParameter *p = findparameter(new, 0);
  return p ? 1 : 0;
}

/* assignment functions */
void Sets(const char *name, const char *value)
{
  setparameter(name, value);
}

void Seti(const char *name, const int i)
{
  char value[100];
  sprintf(value, "%d", i);
  setparameter(name, value);
}

void Setd(const char *name, const double d)
{
  char value[100];
  sprintf(value, "%.20e", d);
  setparameter(name, value);
}

void Appends(const char *name, const char *value)
{
  if (Getv(name, value))
    return;

  {
    char *oldvalue = Gets(name);
    char *newvalue = cmalloc(strlen(oldvalue) + strlen(value) + 2);
    sprintf(newvalue, "%s %s", oldvalue, value);
    while (newvalue[0] == ' ')
      sprintf(newvalue, "%s", &(newvalue[1]));
    setparameter(name, newvalue);
    free(newvalue);
  }
}

/* in order to evaluate analytic expression inside the parfile */
typedef struct tDB
{
  struct tDB *next;
  struct tDB *prev;
  struct tDB *db;
  int b, o;
  double v;
} tdb;

void EvalBracket(tdb *first, tdb *last)
{
  int b;
  tdb *ptr = first;
  tdb *f = NULL;
  tdb *l = NULL;
  tdb *tmp;
  tdb *ptr1;
  int pr = 0;

  /* look for a bracket */
  while (1)
  {

    // locate an "("
    if (ptr->o == 1 && !f)
      f = ptr;

    // locate the corresponding ")"
    if (ptr->o == 1 && f != ptr)
    {
      if (f->b == ptr->b)
      {
        l = ptr;
        if (f->next == l)
          errorexit("something is wrong with the brackets");

        /* find inner bracket */
        EvalBracket(f->next, l->prev);

        ptr = f->next;

        /* delete brackets around this computed value */
        if (f->o == 1)
        {
          f->prev->next = f->next;
          f->next->prev = f->prev;
          if (f == first)
            first = first->next;
          free(f);
        }
        if (l->o == 1)
        {
          if (l->next)
          {
            l->next->prev = l->prev;
            l->prev->next = l->next;
          }
          else
          {
            l->prev->next = NULL;
          }
          if (l == last)
            last = last->prev;
          free(l);
        }

        f = NULL;
      }
    }

    /* move ptr */
    if (ptr == last)
      break;
    ptr = ptr->next;
  }

  /* speed up */
  if (last == first)
    return;

  /* test data */
  if (pr)
  {
    ptr = first->db;
    while (ptr->next)
    {
      ptr = ptr->next;
      printf("  %d  %d   %e   %p\n", ptr->o, ptr->b, ptr->v, ptr);
    }
    printf("\n");
  }

  /* now we are in the innermost bracket */
  int i, del;
  for (i = 8; i > 1; i--)
  {

    ptr = first;
    while (first != last)
    {

      // do the simple calculation with 2 numbers
      if (i == 8 && ptr == first && ptr->o == 2)
      {
        if (!ptr->next || ptr->next->o != 0)
          errorexit("wrong syntax8");
        ptr->v = +ptr->next->v;
        ptr->o = 0;
        del = 1;
      }
      else if (i == 8 && ptr->prev && ptr->prev->o == 6 && ptr->o == 2)
      {
        if (!ptr->next || ptr->next->o != 0)
          errorexit("wrong syntax8.2");
        ptr->v = +ptr->next->v;
        ptr->o = 0;
        del = 1;
      }
      else if (i == 7 && ptr == first && ptr->o == 3)
      {
        if (!ptr->next || ptr->next->o != 0)
          errorexit("wrong syntax7");
        ptr->v = -ptr->next->v;
        ptr->o = 0;
        del = 1;
      }
      else if (i == 7 && ptr->prev && ptr->prev->o == 6 && ptr->o == 3)
      {
        if (!ptr->next || ptr->next->o != 0)
          errorexit("wrong syntax7.2");
        ptr->v = -ptr->next->v;
        ptr->o = 0;
        del = 1;
      }
      else if (i == 6 && ptr->next && ptr->next->o == 6)
      {
        if (!ptr->next->next || ptr->o != 0 || ptr->next->next->o != 0)
          errorexit("wrong syntax6");
        ptr->v = pow(ptr->v, ptr->next->next->v);
        del = 2;
      }
      else if (i == 5 && ptr->next && ptr->next->o == 5)
      {
        if (!ptr->next->next || ptr->o != 0 || ptr->next->next->o != 0)
          errorexit("wrong syntax5");
        ptr->v = ptr->v / ptr->next->next->v;
        del = 2;
      }
      else if (i == 4 && ptr->next && ptr->next->o == 4)
      {
        if (!ptr->next->next || ptr->o != 0 || ptr->next->next->o != 0)
          errorexit("wrong syntax4");
        ptr->v = ptr->v * ptr->next->next->v;
        del = 2;
      }
      else if (i == 3 && ptr->next && ptr->next->o == 3)
      {
        if (!ptr->next->next || ptr->o != 0 || ptr->next->next->o != 0)
          errorexit("wrong syntax3");
        ptr->v = ptr->v - ptr->next->next->v;
        del = 2;
      }
      else if (i == 2 && ptr->next && ptr->next->o == 2)
      {
        if (!ptr->next->next || ptr->o != 0 || ptr->next->next->o != 0)
          errorexit("wrong syntax2");
        ptr->v = ptr->v + ptr->next->next->v;
        del = 2;
      }
      else
      {
        del = 0;
      }

      // delete the operator and the 2nd value
      if (del == 2)
      {

        if (ptr->next->next->next)
        {
          tmp = ptr->next->next->next;
          if (ptr->next->next == last)
            last = ptr;
          free(ptr->next->next);
          free(ptr->next);
          ptr->next = tmp;
          ptr->next->prev = ptr;
        }
        else
        {
          last = ptr;
          free(ptr->next->next);
          free(ptr->next);
          ptr->next = NULL;
        }
        ptr = first;
      }
      else if (del == 1)
      {

        if (ptr->next->next)
        {
          tmp = ptr->next->next;
          if (ptr->next == last)
            last = ptr;
          free(ptr->next);
          ptr->next = tmp;
          ptr->next->prev = ptr;
        }
        else
        {
          last = ptr;
          free(ptr->next);
          ptr->next = NULL;
        }
        ptr = first;
      }
      else
      {

        // next ptr
        if (ptr == last)
          break;
        ptr = ptr->next;
      }

      /* test data */
      if (pr && ptr == first)
      {
        ptr1 = first->db;
        while (ptr1->next)
        {
          ptr1 = ptr1->next;
          printf("  %d  %d   %e   %p\n", ptr1->o, ptr1->b, ptr1->v, ptr1);
        }
        printf("\n");
      }
    }

    /* final test data */
    if (pr)
    {
      ptr = first->db;
      while (ptr->next)
      {
        ptr = ptr->next;
        printf("  %d  %d   %e   %p\n", ptr->o, ptr->b, ptr->v, ptr);
      }
      printf("\n");
    }
  }

  /* test calculation */
  if (pr)
  {
    ptr = first;
    while (ptr->next)
    {
      ptr = ptr->next;
      printf("  %d  %d   %e   %p\n", ptr->o, ptr->b, ptr->v, ptr);
    }
  }

  /* final check */
  if (first != last)
    errorexit("something is wrong with the analysis parser");
}

double EvalString(const char *name, const char *str)
{

  if (strlen(str) >= 0)
  {
    if (!strchr(str, '+') && !strchr(str, '-') && !strchr(str, '*') &&
        !strchr(str, '/') && !strchr(str, '^'))
      return atof(str);
  }

  int i, j;
  int b = 0;
  int ns = 0;
  char s[40];
  int pr = 0;

  tdb *db = malloc(sizeof(tdb));
  tdb *ptr = db;
  db->next = NULL;
  db->prev = NULL;

  /* parse expression */
  for (i = 0; i < strlen(str); i++)
  {

    if (str[i] == ' ')
      continue;

    if ((str[i] == '.') || (str[i] >= '0' && str[i] <= '9'))
    {

      // read in number ... Note, there could be an e
      while ((str[i] == '.') ||
             (str[i] >= '0' && str[i] <= '9') ||
             (str[i] == 'e' && i + 1 < strlen(str) && str[i + 1] >= '0' && str[i + 1] <= '9') ||
             (str[i] == 'e' && i + 1 < strlen(str) && (str[i + 1] == '-' || str[i + 1] == '+')) ||
             (ns > 0 && str[i - 1] == 'e' && (str[i] == '-' || str[i] == '+')))
        s[ns++] = str[i++];
      i--;
    }
    else
    {

      // create an element if we now the number
      if (ns != 0)
      {
        s[ns] = '\0';
        ptr->next = (tdb *)malloc(sizeof(tdb));
        ptr->next->prev = ptr;
        ptr = ptr->next;
        ptr->next = NULL;
        ptr->o = 0;
        ptr->b = 0;
        ptr->v = atof(s);
        ptr->db = db;
        ns = 0;
      }

      // create an element for an operator
      ptr->next = (tdb *)malloc(sizeof(tdb));
      ptr->next->prev = ptr;
      ptr = ptr->next;
      ptr->next = NULL;
      ptr->db = db;

      // set the values for the operator
      ptr->v = 0;
      ptr->b = 0;
      if (str[i] == '(')
      {
        ptr->b = ++b;
        ptr->o = 1;
      }
      else if (str[i] == ')')
      {
        ptr->b = b--;
        ptr->o = 1;
      }
      else if (str[i] == '+')
      {
        ptr->o = 2;
      }
      else if (str[i] == '-')
      {
        ptr->o = 3;
      }
      else if (str[i] == '*')
      {
        ptr->o = 4;
      }
      else if (str[i] == '/')
      {
        ptr->o = 5;
      }
      else if (str[i] == '^')
      {
        ptr->o = 6;
      }
      else
      {
        printf("EXPR:  %s\n", str);
        errorexit("wrong operator used");
      }
    }
  }

  if (ns != 0)
  {
    s[ns] = '\0';
    ptr->next = (tdb *)malloc(sizeof(tdb));
    ptr->next->prev = ptr;
    ptr = ptr->next;
    ptr->next = NULL;
    ptr->o = 0;
    ptr->b = 0;
    ptr->v = atof(s);
    ptr->db = db;
    ns = 0;
  }

  if (pr)
  {
    printf("-------------------------------------------\n");
    printf("-->%s<-- = -->%s<--\n", name, str);
    ptr = db;
    while (ptr->next)
    {
      ptr = ptr->next;
      printf("  %d  %d   %e   %p\n", ptr->o, ptr->b, ptr->v, ptr);
    }
    printf("\n");
  }

  /* evaluate the stored expression */
  EvalBracket(db->next, ptr);
  double v = db->next->v;
  free(db);

  if (pr)
  {
    printf("result = %e\n", v);
    printf("-------------------------------------------\n");
  }

  return v;
}

/* query functions */
char *Gets(const char *name)
{
  tParameter *p = findparameter(name, 1);
  return p->value;
}

char *GetsLax(const char *name)
{
  tParameter *p = findparameter(name, 0);
  return p ? p->value : 0;
}

int Geti(const char *name)
{
#if 0
  tParameter *p = findparameter(name, 1);
  return atoi(p->value);
#else
  tParameter *p = findparameter(name, 1);
  return (int)(EvalString(name, p->value));
#endif
}

int GetiLax(const char *name)
{
  tParameter *p = findparameter(name, 0);
  return p ? atoi(p->value) : 0;
}

double Getd(const char *name)
{
  /* switch between version for parsing math expressions */
#if 0
  tParameter *p = findparameter(name, 1);
  return atof(p->value);
#else
  tParameter *p = findparameter(name, 1);
  return EvalString(name, p->value);
#endif
}

double GetdLax(const char *name)
{
  tParameter *p = findparameter(name, 0);
  return p ? atof(p->value) : 0.;
}

/* "get value" returns 1 if value is in the list of values and 0 else
   (not equivalent to value being a substring of string parameter)
*/
int GetvFlag(const char *name, const char *value, const int fatal)
{
  tParameter *p = findparameter(name, fatal);
  char *s, c;
  int startok, endok;

  if (!p)
    return 0;
  s = strstr(p->value, value);
  if (!s)
    return 0;
  startok = (s == p->value || *(s - 1) == ' ');
  c = s[strlen(value)];
  endok = (c == 0 || c == ' ');
  return startok && endok ? 1 : 0;
}

int GetvLax(const char *name, const char *value)
{
  return GetvFlag(name, value, 0);
}

int Getv(const char *name, const char *value)
{
  return GetvFlag(name, value, 1);
}

/* weakly tested */
/* arrays, could be optimized if needed */
double *GetdArray(const char *name, int *n)
{
  char *s = Gets(name), *t;
  double *a;
  int i;

  for (i = 0; NextEntry(s); i++)
    ; /* count */
  a = dmalloc(i);
  for (i = 0; (t = NextEntry(s)); i++)
    a[i] = atof(t);

  *n = i;
  return a;
}

int GetdArrayN(const char *name)
{
  char *s = Gets(name);
  int i;

  for (i = 0; NextEntry(s); i++)
    ;

  return i;
}

double GetdEntry(const char *name, const int n)
{
  int na;
  double *a = GetdArray(name, &na);
  double x;

  if (n < 0 || n >= na)
  {
    printf("no entry number %d in parameter %s\n", n, name);
    errorexit("GetdEntry: please check your parameters");
  }

  x = a[n];
  free(a);
  return x;
}

void SetdArray(const char *name, const int n, const double *a)
{
  tParameter *p = findparameter(name, 1);
  char *t;
  int i;

  if (n < 1 || a == 0)
    return;

  free(p->value);
  t = p->value = cmalloc(25 * n);
  *t = 0;
  for (i = 0;;)
  {
    t += strlen(t);
    sprintf(t, "%.16e", a[i]);
    if (++i == n)
      break;
    strcat(t, " ");
  }

  if (0)
    printf("SetdArray: %s = %s\n", name, p->value);
}

/* "for each"
   return next string delimited by SINGLE space (as in parameter data base)
   0 if end of list, then restart
   list = 0 or *list = "" restarts, too
*/
/* how universal is strdup? strsep? */
char *NextEntry(const char *list)
{
  static char *copyoflist = 0;
  static int i = 0, l = 0;
  static int firstcall = 1;
  char *s;

  if (!list || !*list || i && i == l)
  {
    i = l = 0;
    return 0;
  }

  if (!i)
  {
    l = strlen(list);
    copyoflist = (char *)realloc(copyoflist, sizeof(char) * (l + 1));
    if (firstcall)
    {
      hidepointer(copyoflist);
      firstcall = 0;
    }
    strcpy(copyoflist, list);
  }
  else
    i++;

  s = copyoflist + i;
  for (; i < l && copyoflist[i] != ' '; i++)
    ;
  if (copyoflist[i] == ' ')
    copyoflist[i] = 0;

  return s;
}

/* direct index access to parameter data base
   should not be needed except for, say, dumping the database for checkpointing
*/
char *GetsInd(const int i)
{
  if (0 <= i && i < npdb)
    return pdb[i].value;
  return 0;
}

char *GetnameInd(const int i)
{
  if (0 <= i && i < npdb)
    return pdb[i].name;
  return 0;
}

int GetnParameters()
{
  return npdb;
}