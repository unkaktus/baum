#include <stdbool.h>
#include "fukaccia.h"

double *allocate_double(int n);
Fields allocate_fields(int n_points);

void write_variable(FILE *file, char *varname, double *data, int length);
void read_variable(FILE *file, char *varname, double *data, int length);

void save_fields(tL *level, Fields fields, bool hydro_enabled);
Fields load_fields(tL *level, bool hydro_enabled);

void interpolated_data_mark_ready();
int interpolated_data_is_ready();