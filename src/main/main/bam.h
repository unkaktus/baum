/* bam.h */
/* Bernd Bruegmann, 12/99, Wolfgang Tichy 7/2015 */



#include <dirent.h>
#include <fcntl.h>
#include <float.h>
#include <limits.h>
#ifndef __APPLE__
// malloc.h is non-standard and not avaialable for Apple devices
#include <malloc.h>
#endif
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <time.h>
#ifndef __APPLE__
// not available on Mac OS, but I think this is not used anywhere, so could also be removed
//#include <sys/sendfile.h>
#endif
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/wait.h>   /* for waitpid */
#include <unistd.h>



#ifdef REDUCEORDERTO2
#define REDUCEORDERTO4
#endif
#ifdef REDUCEORDERTO4
#define REDUCEORDERTO6
#endif
#ifdef REDUCEORDERTO6
#define REDUCEORDERTO8
#endif


//mth: NO SUPPORT FOR NOT USING THESE TWO FLAGS ANY LONGER
#ifndef BOX
#define BOX
#endif

#ifndef BOXES
#define BOXES
#endif




#include "bam_automatic_include.h"
