#pragma once
#include "fukaccia.h"

typedef struct SpatialVector
{
    double *x, *y, *z;
} SpatialVector;

typedef struct SpatialRank2Tensor
{
    double *xx, *xy, *xz;
    double *yy, *yz;
    double *zz;

} SpatialRank2Tensor;

void allocate_variable(tL *level, const char *varname);
void allocate_variables(tL *level, char **var_list, size_t n_vars);

double *get_scalar(tL *level, const char *varname);

SpatialVector *get_spatial_vector(tL *level, const char *varname);
SpatialRank2Tensor *get_spatial_rank2_tensor(tL *level, const char *varname);

int det(SpatialRank2Tensor *g, int i);

Grid get_grid(tL *level);
