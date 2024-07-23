/* bam_checkpoint.h */
/* Bernd Bruegmann 01/2006 */

#include "wolfio.h"

void checkpoint(tG *g, int l);
int checkpoint_checkforfiles(char *suffix);
char *checkpoint_filename_from_outdir_suffix(char *outdir, char *suffix);
