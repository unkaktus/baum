#include "helpers.h"
#include "fuka_reader.h"
#include "levelfile.h"
#include <omp.h>

void set_variables_from_fields(tL *level, Fields fields, bool hydro_enabled)
{
    for (int i = 0; i < level->npoints; i++)
    {
        double *alpha = get_scalar(level, "alpha");
        alpha[i] = fields.alpha[i];

        SpatialVector *beta = get_spatial_vector(level, "betax");
        beta->x[i] = fields.beta_x[i];
        beta->y[i] = fields.beta_y[i];
        beta->z[i] = fields.beta_z[i];

        SpatialRank2Tensor *gamma = get_spatial_rank2_tensor(level, "adm_gxx");
        gamma->xx[i] = fields.gamma_xx[i];
        gamma->xy[i] = fields.gamma_xy[i];
        gamma->xz[i] = fields.gamma_xz[i];
        gamma->yy[i] = fields.gamma_yy[i];
        gamma->yz[i] = fields.gamma_yz[i];
        gamma->zz[i] = fields.gamma_zz[i];

        double det_gamma = det(gamma, i);
        if (det_gamma <= 0)
        {
            log_printf("det_gamma = %lf <=0, stopping.", det_gamma);
            graceful_shutdown(1);
        }

        SpatialRank2Tensor *K = get_spatial_rank2_tensor(level, "adm_Kxx");
        K->xx[i] = fields.K_xx[i];
        K->xy[i] = fields.K_xy[i];
        K->xz[i] = fields.K_xz[i];
        K->yy[i] = fields.K_yy[i];
        K->yz[i] = fields.K_yz[i];
        K->zz[i] = fields.K_zz[i];

        // Set matter fields
        if (hydro_enabled)
        {
            double *grhd_rho = get_scalar(level, "grhd_rho");
            grhd_rho[i] = fields.rho[i];

            double *grhd_p = get_scalar(level, "grhd_p");
            grhd_p[i] = fields.pressure[i];

            double *grhd_epsl = get_scalar(level, "grhd_epsl");
            grhd_epsl[i] = fields.epsilon[i];

            SpatialVector *grhd_v = get_spatial_vector(level, "grhd_vx");
            grhd_v->x[i] = fields.v_x[i];
            grhd_v->y[i] = fields.v_y[i];
            grhd_v->z[i] = fields.v_z[i];
        }
    }
}

int fuka_reader_pre_grid(tL *level)
{
    char *binary_type_str = Gets("fuka_reader_binary_type");
    if (strcmp(binary_type_str, "") == 0)
    {
        log_printf("\"fuka_reader_binary_type\" is not set");
        graceful_shutdown(1);
    }
    log_printf("this is a %s binary", binary_type_str);

    char *info_filename = Gets("fuka_reader_info_filename");
    if (strcmp(info_filename, "") == 0)
    {
        log_printf("\"fuka_reader_info_filename\" is not set");
        graceful_shutdown(1);
    }
    log_printf("info_filename: %s", info_filename);

    BinaryType binary_type = binary_type_from_string(binary_type_str);
    if (binary_type < 0)
    {
        log_printf("unknown binary type: %s", binary_type_str);
        graceful_shutdown(1);
    }
    BinaryInfo binary_info = read_binary_info(binary_type, info_filename);

    set_binary_parameters(&binary_info);

    return 0;
}

int fuka_reader_set_fields(tL *level)
{
    log_printf("level %d: started initial data filling routine", level->l);

    BinaryType binary_type = binary_type_from_string(Gets("fuka_reader_binary_type"));

    bool hydro_enabled = false;
    if (binary_type == BNS || binary_type == BHNS)
    {
        hydro_enabled = true;
    }

    char *vacuum_variables[] = {
        "adm_psi",         // Conformal factor
        "adm_dpsiopsix",   // Conformal factor derivative
        "adm_ddpsiopsixx", // Conformal factor second derivative
        "alpha",           // Lapse
        "betax",           // Shift
        "adm_gxx",         // Induced spatial metric
        "adm_Kxx",         // Extrinsic curvature
    };

    char *hydro_variables[] = {
        "grhd_rho",  // Primitive density
        "grhd_epsl", // Specific energy density
        "grhd_p",    // Primitive pressure
        "grhd_vx",   // Velocity
        "adm_rho",   // ADM density
    };

    // Allocate memory for all the variables
    allocate_variables(level, vacuum_variables,
                       sizeof(vacuum_variables) / sizeof(char *));
    if (hydro_enabled)
    {
        log_printf("matter is present, allocating hydro variables");
        allocate_variables(level, hydro_variables,
                           sizeof(hydro_variables) / sizeof(char *));
    }

    // Set the stage for the ID
    adm_Minkowski(level);

    log_printf("level %d: grid has %d points", level->l, get_grid(level).n_points);

    Grid grid = get_grid(level);

    Fields fields;

    if (interpolated_data_is_ready())
    {
        log_printf("level %d: there is already interpolated data present, loading it", level->l);
        fields = load_fields(level, hydro_enabled);
    }
    else
    {
        log_printf("level %d: interpolating ID onto the level grid", level->l);

        FUKAInterpolateRequest req = {
            .binary_type = binary_type,
            .info_filename = Gets("fuka_reader_info_filename"),
            .grid = &grid,
            .interpolation_offset = Getd("fuka_reader_interpolation_offset"),
            .interpolation_order = Geti("fuka_reader_interpolation_order"),
            .relative_dr_spacing = Getd("fuka_reader_relative_dr_spacing"),
        };

        fields = fukaccia_interpolate(&req, omp_get_max_threads());

        log_printf("level %d: saving the interpolated data", level->l);
        save_fields(level, fields, hydro_enabled);
    }

    // Fill the fields with the FUKA ID
    log_printf("setting variables from fields");
    set_variables_from_fields(level, fields, hydro_enabled);
    fukaccia_finalize();

    return 0;
}

int fuka_reader_finish(tL *level)
{
    bampi_barrier();

    // Mark the saved interpolated data as ready
    if (processor0)
    {
        log_printf("marking the interpolated data as ready");
        interpolated_data_mark_ready();
    }

    return 0;
}