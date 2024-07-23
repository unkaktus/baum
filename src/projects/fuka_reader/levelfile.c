#define _GNU_SOURCE
#include <string.h>
#include <errno.h>

#include "bam.h"
#include "levelfile.h"
#include "helpers.h"
#include "fukaccia.h"

double *allocate_double(int n)
{
    double *pointer = (double *)calloc(n, sizeof(double));
    if (pointer == NULL)
    {
        printf("cannot allocate memory\n");
        exit(1);
    }
    return pointer;
}

Fields allocate_fields(int n_points)
{
    Fields fields = {
        .alpha = allocate_double(n_points),
        .beta_x = allocate_double(n_points),
        .beta_y = allocate_double(n_points),
        .beta_z = allocate_double(n_points),
        .gamma_xx = allocate_double(n_points),
        .gamma_xy = allocate_double(n_points),
        .gamma_xz = allocate_double(n_points),
        .gamma_yy = allocate_double(n_points),
        .gamma_yz = allocate_double(n_points),
        .gamma_zz = allocate_double(n_points),
        .K_xx = allocate_double(n_points),
        .K_xy = allocate_double(n_points),
        .K_xz = allocate_double(n_points),
        .K_yy = allocate_double(n_points),
        .K_yz = allocate_double(n_points),
        .K_zz = allocate_double(n_points),
        // Hydro
        .rho = allocate_double(n_points),
        .epsilon = allocate_double(n_points),
        .pressure = allocate_double(n_points),
        .v_x = allocate_double(n_points),
        .v_y = allocate_double(n_points),
        .v_z = allocate_double(n_points),
    };
    return fields;
}

void write_variable(FILE *file, char *varname, double *data, int length)
{
    fprintf(file, "$variable = %s : length = %d\n", varname, length);
    fwrite(data, sizeof(double), length, file);
    fprintf(file, "\n");
}

void read_variable(FILE *file, char *varname, double *data, int length)
{
    char varname_file[10000];
    int status = fgetparameter(file, "$variable", varname_file);
    if (status == EOF)
    {
        log_printf("EOF reached");
        graceful_shutdown(1);
    }

    if (strcmp(varname_file, varname) != 0)
    {
        log_printf("the name of variable in the file does match the requested one");
        graceful_shutdown(1);
    }

    char length_str[10000];
    fgetparameter(file, "length", length_str);
    int length_file = atoi(length_str);
    if (length_file != length)
    {
        log_printf("number of points is different: n_points=%d, level->npoints=%d var: %s\n", length_file, length, varname);
        graceful_shutdown(1);
    }

    // Read the trailing newline
    char buffer[10000];
    fscanline(file, buffer);

    // Read the variable data
    fread(data, sizeof(double), length, file);
}

void save_fields(tL *level, Fields fields, bool hydro_enabled)
{
    // Create a separate directory for the interpolated data
    char *output_dir;
    asprintf(&output_dir, "%s/fuka-interpolated-id", Gets("outdir"));

    if (!system_isdir(output_dir))
    {
        system_mkdir(output_dir);
    }

    // Open the level file
    char *level_filename;
    asprintf(&level_filename, "%s/level%d.%d", output_dir, level->l, bampi_rank());

    FILE *level_file = fopen(level_filename, "wb");
    if (level_file == NULL)
    {
        log_printf("cannot open file %s: %s", level_filename, strerror(errno));
        graceful_shutdown(1);
    }

    // Write the header
    fprintf(level_file, "$BEGIN_variables:\n");

    // Write all the variables
    write_variable(level_file, "alpha", fields.alpha, level->npoints);
    write_variable(level_file, "beta_x", fields.beta_x, level->npoints);
    write_variable(level_file, "beta_y", fields.beta_y, level->npoints);
    write_variable(level_file, "beta_z", fields.beta_z, level->npoints);
    write_variable(level_file, "gamma_xx", fields.gamma_xx, level->npoints);
    write_variable(level_file, "gamma_xy", fields.gamma_xy, level->npoints);
    write_variable(level_file, "gamma_xz", fields.gamma_xz, level->npoints);
    write_variable(level_file, "gamma_yy", fields.gamma_yy, level->npoints);
    write_variable(level_file, "gamma_yz", fields.gamma_yz, level->npoints);
    write_variable(level_file, "gamma_zz", fields.gamma_zz, level->npoints);
    write_variable(level_file, "K_xx", fields.K_xx, level->npoints);
    write_variable(level_file, "K_xy", fields.K_xy, level->npoints);
    write_variable(level_file, "K_xz", fields.K_xz, level->npoints);
    write_variable(level_file, "K_yy", fields.K_yy, level->npoints);
    write_variable(level_file, "K_yz", fields.K_yz, level->npoints);
    write_variable(level_file, "K_zz", fields.K_zz, level->npoints);

    if (hydro_enabled)
    {
        write_variable(level_file, "rho", fields.rho, level->npoints);
        write_variable(level_file, "epsilon", fields.epsilon, level->npoints);
        write_variable(level_file, "pressure", fields.pressure, level->npoints);
        write_variable(level_file, "v_x", fields.v_x, level->npoints);
        write_variable(level_file, "v_y", fields.v_y, level->npoints);
        write_variable(level_file, "v_z", fields.v_z, level->npoints);
    }

    fclose(level_file);

    free(output_dir);
    free(level_filename);
}

Fields load_fields(tL *level, bool hydro_enabled)
{
    Fields fields = allocate_fields(level->npoints);

    char *output_dir;
    asprintf(&output_dir, "%s_previous/fuka-interpolated-id", Gets("outdir"));

    // Open the level file
    char *level_filename;
    asprintf(&level_filename, "%s/level%d.%d", output_dir, level->l, bampi_rank());

    FILE *level_file = fopen(level_filename, "rb");
    if (level_file == NULL)
    {
        log_printf("cannot open file %s: %s", level_filename, strerror(errno));
        graceful_shutdown(1);
    }

    // Find the variable section
    if (fgotonext(level_file, "$BEGIN_variables:") == EOF)
    {
        log_printf("$BEGIN_variables: is missing");
        graceful_shutdown(1);
    }

    // Read all the variables
    read_variable(level_file, "alpha", fields.alpha, level->npoints);
    read_variable(level_file, "beta_x", fields.beta_x, level->npoints);
    read_variable(level_file, "beta_y", fields.beta_y, level->npoints);
    read_variable(level_file, "beta_z", fields.beta_z, level->npoints);
    read_variable(level_file, "gamma_xx", fields.gamma_xx, level->npoints);
    read_variable(level_file, "gamma_xy", fields.gamma_xy, level->npoints);
    read_variable(level_file, "gamma_xz", fields.gamma_xz, level->npoints);
    read_variable(level_file, "gamma_yy", fields.gamma_yy, level->npoints);
    read_variable(level_file, "gamma_yz", fields.gamma_yz, level->npoints);
    read_variable(level_file, "gamma_zz", fields.gamma_zz, level->npoints);
    read_variable(level_file, "K_xx", fields.K_xx, level->npoints);
    read_variable(level_file, "K_xy", fields.K_xy, level->npoints);
    read_variable(level_file, "K_xz", fields.K_xz, level->npoints);
    read_variable(level_file, "K_yy", fields.K_yy, level->npoints);
    read_variable(level_file, "K_yz", fields.K_yz, level->npoints);
    read_variable(level_file, "K_zz", fields.K_zz, level->npoints);

    if (hydro_enabled)
    {
        read_variable(level_file, "rho", fields.rho, level->npoints);
        read_variable(level_file, "epsilon", fields.epsilon, level->npoints);
        read_variable(level_file, "pressure", fields.pressure, level->npoints);
        read_variable(level_file, "v_x", fields.v_x, level->npoints);
        read_variable(level_file, "v_y", fields.v_y, level->npoints);
        read_variable(level_file, "v_z", fields.v_z, level->npoints);
    }

    fclose(level_file);
    return fields;
}

void interpolated_data_mark_ready()
{
    char *readiness_filename;
    asprintf(&readiness_filename, "%s/fuka-interpolated-id/ready", Gets("outdir"));
    touch(readiness_filename);
    free(readiness_filename);
}

int interpolated_data_is_ready()
{
    char *readiness_filename;
    asprintf(&readiness_filename, "%s_previous/fuka-interpolated-id/ready", Gets("outdir"));
    int result = system_isfile(readiness_filename);
    free(readiness_filename);
    return result;
}