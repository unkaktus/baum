#include "bam.h"
#include "fields.h"
#include "helpers.h"

void allocate_variable(tL *level, const char *varname)
{
    int i = Ind(varname);
    enablevar(level, i);
}

void allocate_variables(tL *level, char **var_list, size_t n_vars)
{
    for (int i = 0; i < n_vars; i++)
    {
        char *varname = var_list[i];
        allocate_variable(level, varname);
    }
}

double *get_scalar(tL *level, const char *varname)
{
    int i = Ind(varname);
    if (i == 0)
    {
        log_printf("index variable %s is not found", varname);
    }
    double *pointer = level->v[i];
    if (pointer == NULL)
    {
        log_printf("index variable %s is found (%d), but the pointer is NULL", varname, i);
    }
    return pointer;
}

SpatialVector *get_spatial_vector(tL *level, const char *varname)
{
    int i = Ind(varname);
    SpatialVector *vector = malloc(sizeof(*vector));
    vector->x = level->v[i + 0];
    vector->y = level->v[i + 1];
    vector->z = level->v[i + 2];
    return vector;
}

SpatialRank2Tensor *get_spatial_rank2_tensor(tL *level, const char *varname)
{
    int i = Ind(varname);
    SpatialRank2Tensor *tensor = malloc(sizeof(*tensor));
    tensor->xx = level->v[i + 0];
    tensor->xy = level->v[i + 1];
    tensor->xz = level->v[i + 2];
    tensor->yy = level->v[i + 3];
    tensor->yz = level->v[i + 4];
    tensor->zz = level->v[i + 5];

    return tensor;
}

int det(SpatialRank2Tensor *g, int i)
{
    double det_g = (2. * g->xy[i] * g->xz[i] * g->yz[i] +
                    g->xx[i] * g->yy[i] * g->zz[i] -
                    g->zz[i] * g->xy[i] * g->xy[i] -
                    g->yy[i] * g->xz[i] * g->xz[i] -
                    g->xx[i] * g->yz[i] * g->yz[i]);
    return det_g;
}

// Retrieve grid coordinates from the level
Grid get_grid(tL *level)
{
    Grid grid = {
        .x = level->v[Ind("x")],
        .y = level->v[Ind("y")],
        .z = level->v[Ind("z")],
        .n_points = level->npoints,
    };

    return grid;
}
