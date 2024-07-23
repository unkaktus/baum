#include <stdbool.h>
#include "helpers.h"
#include "bam.h"
#include "fields.h"
#include "binary.h"
#include "fukaccia.h"

void set_variables_from_fields(tL *level, Fields fields, bool hydro_enabled);

// Hooks
int fuka_reader_pre_grid(tL *level);
int fuka_reader_set_fields(tL *level);
int fuka_reader_finish(tL *level);
